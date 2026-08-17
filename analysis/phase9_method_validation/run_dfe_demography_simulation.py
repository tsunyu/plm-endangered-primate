#!/usr/bin/env python3
"""Forward DFE simulations under inferred demography for transient drift load."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

from dfe_demography_common import (
    BASE,
    BIN_LABELS,
    PARAM_FILE,
    Demography,
    bin_frequencies,
    demography_table,
    load_demography_mle,
    load_observed_bins,
    simulate_replicate,
    summarize_replicate,
    write_metadata,
)

OUT_ROOT = BASE / "output/method_validation/dfe_simulation"
SEED = 20260710

SCENARIOS = {
    "transient_drift": {
        "description": "Inferred fastsimcoal2 piecewise Ne with deleterious gamma-DFE",
        "demography": "inferred",
        "s_mean": 0.0015,
        "s_shape": 0.25,
        "h": 0.10,
        "seed_standing": True,
    },
    "lof_stronger": {
        "description": "Inferred demography with stronger LoF-like selection",
        "demography": "inferred",
        "s_mean": 0.008,
        "s_shape": 0.25,
        "h": 0.10,
        "seed_standing": True,
    },
    "constant_current_ne": {
        "description": "Constant Ne=150 with same DFE (recent-small-population null)",
        "demography": "constant_current",
        "s_mean": 0.0015,
        "s_shape": 0.25,
        "h": 0.10,
        "seed_standing": True,
    },
    "constant_ancestral_ne": {
        "description": "Constant Ne=NANC with same DFE (no bottleneck null)",
        "demography": "constant_ancestral",
        "s_mean": 0.0015,
        "s_shape": 0.25,
        "h": 0.10,
        "seed_standing": True,
    },
    "neutral_drift": {
        "description": "Inferred demography with neutral alleles (s=0)",
        "demography": "inferred",
        "s_mean": 0.0,
        "s_shape": 1.0,
        "h": 0.25,
    },
}


def build_demography(kind: str, mle: Demography) -> Demography:
    if kind == "inferred":
        return mle
    if kind == "constant_current":
        return Demography(
            nanc=mle.ncur,
            nbot=mle.ncur,
            nrecover=mle.ncur,
            ncur=mle.ncur,
            tbot_old=mle.tbot_old,
            trecovery_old=mle.trecovery_old,
            trecent=mle.trecent,
            burn_in=mle.burn_in,
        )
    if kind == "constant_ancestral":
        return Demography(
            nanc=mle.nanc,
            nbot=mle.nanc,
            nrecover=mle.nanc,
            ncur=mle.nanc,
            tbot_old=mle.tbot_old,
            trecovery_old=mle.trecovery_old,
            trecent=mle.trecent,
            burn_in=mle.burn_in,
        )
    raise ValueError(f"Unknown demography kind: {kind}")


def simulate_neutral_replicate(demography: Demography, n_loci: int, rng: np.random.Generator) -> np.ndarray:
    p = rng.uniform(1e-4, 0.05, size=n_loci)
    for gen in range(demography.total_generations):
        ne = demography.ne_at(gen)
        counts = rng.binomial(max(int(round(2.0 * ne)), 2), p)
        p = counts / (2.0 * ne)
        if np.all((p <= 0.0) | (p >= 1.0)):
            break
    return p.astype(float)


def run_scenario(
    scenario: str,
    config: dict,
    mle: Demography,
    n_replicates: int,
    n_loci: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    demography = build_demography(config["demography"], mle)
    rows = []
    for rep in range(n_replicates):
        if config["s_mean"] <= 0:
            freqs = simulate_neutral_replicate(demography, n_loci, rng)
        else:
            freqs = simulate_replicate(
                demography,
                n_loci=n_loci,
                s_mean=config["s_mean"],
                s_shape=config["s_shape"],
                h=config["h"],
                rng=rng,
                seed_standing=bool(config.get("seed_standing", False)),
            )
        summary = summarize_replicate(freqs)
        summary.update(
            {
                "scenario": scenario,
                "replicate": rep,
                "n_loci": n_loci,
            }
        )
        rows.append(summary)
    return pd.DataFrame(rows)


def compare_to_observed(
    scenario_summary: pd.DataFrame,
    observed: pd.DataFrame,
    variant_class: str,
) -> pd.DataFrame:
    obs = observed[observed["variant_class"] == variant_class]
    obs_fracs = obs[obs["bin"].isin(BIN_LABELS)].set_index("bin")["fraction"]
    sim = scenario_summary[
        [f"frac_{label}" for label in BIN_LABELS]
    ].mean()
    sim_fracs = pd.Series({label: sim[f"frac_{label}"] for label in BIN_LABELS})
    table = pd.DataFrame(
        {
            "bin": BIN_LABELS,
            "observed_fraction": obs_fracs.values,
            "simulated_fraction": sim_fracs.values,
        }
    )
    obs_counts = (obs_fracs * 1000).round().astype(int).values
    sim_counts = (sim_fracs * 1000).round().astype(int).values
    nonzero = (obs_counts + sim_counts) > 0
    if (
        nonzero.sum() >= 2
        and obs_counts[nonzero].sum() > 0
        and sim_counts[nonzero].sum() > 0
    ):
        chi2, p_value, _, _ = chi2_contingency([obs_counts[nonzero], sim_counts[nonzero]])
    else:
        chi2, p_value = np.nan, np.nan
    table["variant_class"] = variant_class
    table["chi2"] = chi2
    table["chi2_pvalue"] = p_value
    return table


def plot_results(
    observed: pd.DataFrame,
    scenario_summary: pd.DataFrame,
    comparisons: pd.DataFrame,
    output: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    categories = ["Rare\n(<0.05)", "Low\n(0.05-0.20)", "Intermediate\n(0.20-0.50)", "Common\n(0.50-0.95)", "Fixed\n(≥0.95)"]

    ax = axes[0, 0]
    for variant_class, color in (("LoF", "#D55E00"), ("Deleterious Missense", "#0072B2")):
        obs = observed[observed["variant_class"] == variant_class]
        fracs = obs[obs["bin"].isin(BIN_LABELS)].set_index("bin").loc[BIN_LABELS, "fraction"].values * 100
        ax.plot(categories, fracs, marker="o", linewidth=2, label=f"Observed {variant_class}")
    ax.set_ylabel("Percentage of variants (%)")
    ax.set_title("Observed deleterious allele-frequency spectra")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    for scenario, group in scenario_summary.groupby("scenario"):
        means = group[[f"frac_{label}" for label in BIN_LABELS]].mean() * 100
        ax.plot(
            categories,
            [means[f"frac_{label}"] for label in BIN_LABELS],
            marker="o",
            linewidth=1.5,
            alpha=0.85,
            label=scenario,
        )
    ax.set_ylabel("Percentage of variants (%)")
    ax.set_title("Simulated spectra across scenarios")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="best")

    ax = axes[1, 0]
    dem = demography_table(load_demography_mle())
    ax.step(dem["generation"], dem["ne"], where="post", color="#009E73", linewidth=2)
    ax.set_yscale("log")
    ax.set_xlabel("Generation")
    ax.set_ylabel("Diploid Ne")
    ax.set_title("Inferred piecewise Ne trajectory")
    ax.grid(alpha=0.3)

    ax = axes[1, 1]
    focus = comparisons[comparisons["variant_class"] == "Deleterious Missense"]
    x = np.arange(len(BIN_LABELS))
    width = 0.35
    ax.bar(x - width / 2, focus["observed_fraction"] * 100, width, label="Observed missense", color="#0072B2", alpha=0.8)
    ax.bar(x + width / 2, focus["simulated_fraction"] * 100, width, label="Transient drift sim.", color="#E69F00", alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=8)
    ax.set_ylabel("Percentage (%)")
    ax.set_title("Observed vs simulated missense spectrum")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")

    plt.tight_layout()
    fig.savefig(output / "dfe_simulation_summary.png", dpi=300, bbox_inches="tight")
    plt.close()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUT_ROOT)
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--loci", type=int, default=1000)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    mle = load_demography_mle()
    observed = load_observed_bins()

    all_rows = []
    for scenario, config in SCENARIOS.items():
        frame = run_scenario(scenario, config, mle, args.replicates, args.loci, rng)
        frame.to_csv(args.output / f"{scenario}_replicates.csv", index=False)
        all_rows.append(frame)
    scenario_summary = pd.concat(all_rows, ignore_index=True)
    scenario_summary.to_csv(args.output / "scenario_replicate_summary.csv", index=False)

    comparisons = []
    for variant_class in ("LoF", "Deleterious Missense"):
        sim_class = "lof_stronger" if variant_class == "LoF" else "transient_drift"
        comp = compare_to_observed(
            scenario_summary[scenario_summary["scenario"] == sim_class],
            observed,
            variant_class,
        )
        comp["scenario"] = sim_class
        comparisons.append(comp)
    comparison_df = pd.concat(comparisons, ignore_index=True)
    comparison_df.to_csv(args.output / "observed_vs_simulated_bins.csv", index=False)

    plot_results(observed, scenario_summary, comparison_df, args.output)

    write_metadata(
        args.output / "simulation_metadata.json",
        {
            "engine": "wright-fisher",
            "slim_crosscheck": str(args.output / "slim_crosscheck_summary.csv"),
            "seed": SEED,
            "replicates": args.replicates,
            "loci_per_replicate": args.loci,
            "demography_mle": mle.__dict__,
            "scenarios": SCENARIOS,
            "interpretation": (
                "Forward simulations test whether the inferred bottleneck-with-recent-contraction "
                "history can elevate weakly deleterious alleles to intermediate frequencies without fixation."
            ),
        },
    )

    summary = scenario_summary.groupby("scenario").agg(
        mean_daf=("mean_daf", "mean"),
        frac_intermediate=("frac_intermediate", "mean"),
        frac_common=("frac_common", "mean"),
        frac_fixed=("frac_fixed", "mean"),
    )
    summary.to_csv(args.output / "scenario_means.csv")
    (args.output / "SUMMARY.md").write_text(
        "\n".join(
            [
                "# DFE / demography forward simulation",
                "",
                f"- Replicates per scenario: {args.replicates}",
                f"- Loci per replicate: {args.loci}",
                f"- Inferred demography source: `{PARAM_FILE}`",
                "",
                "## Scenario means",
                "",
                "```csv",
                summary.to_csv(),
                "```",
                "",
                "## Comparison to observed bins",
                "",
                "```csv",
                comparison_df.to_csv(index=False),
                "```",
            ]
        ),
        encoding="utf-8",
    )
    print(f"Wrote simulation outputs to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
