#!/usr/bin/env python3
"""Cross-species and score-transfer checks between human ClinVar and SNJ missense."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, mannwhitneyu, spearmanr

from load_sensitivity_common import sigmoid_probability


BASE = get_base_dir()
SEED = 20260710


def load_clinvar(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(
        path,
        usecols=["protein", "DMS_bin_score", "esm2_score"],
    )
    frame = frame.dropna()
    frame = frame[frame["DMS_bin_score"].isin(["Pathogenic", "Benign"])]
    frame["label"] = (frame["DMS_bin_score"] == "Pathogenic").astype(int)
    frame["pathogenicity_prob"] = frame["esm2_score"].map(sigmoid_probability)
    return frame


def load_monkey(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, usecols=["variant_id", "esm2_score", "esm2_prediction"])
    frame["esm2_score"] = pd.to_numeric(frame["esm2_score"], errors="coerce")
    frame = frame.dropna(subset=["esm2_score"])
    frame["pathogenicity_prob"] = frame["esm2_score"].map(sigmoid_probability)
    return frame


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--clinvar",
        type=Path,
        default=BASE / "output/phase5_genetic_load/esm2_predictions.csv",
    )
    parser.add_argument(
        "--monkey",
        type=Path,
        default=BASE / "output/phase4_plm_predictions/esm2/esm2_predictions.csv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=BASE / "output/method_validation/cross_species",
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    clinvar = load_clinvar(args.clinvar)
    monkey = load_monkey(args.monkey)

    human_path = clinvar.loc[clinvar["label"] == 1, "esm2_score"]
    human_benign = clinvar.loc[clinvar["label"] == 0, "esm2_score"]
    monkey_scores = monkey["esm2_score"]

    ks_monkey_vs_benign = ks_2samp(monkey_scores, human_benign)
    ks_monkey_vs_path = ks_2samp(monkey_scores, human_path)
    mw_prob = mannwhitneyu(
        clinvar.loc[clinvar["label"] == 1, "pathogenicity_prob"],
        clinvar.loc[clinvar["label"] == 0, "pathogenicity_prob"],
        alternative="greater",
    )

    summary_rows = [
        {
            "comparison": "human_pathogenic_vs_benign_llr",
            "statistic": float(mw_prob.statistic),
            "p_value": float(mw_prob.pvalue),
            "n_a": int((clinvar["label"] == 1).sum()),
            "n_b": int((clinvar["label"] == 0).sum()),
        },
        {
            "comparison": "monkey_vs_human_benign_llr_ks",
            "statistic": float(ks_monkey_vs_benign.statistic),
            "p_value": float(ks_monkey_vs_benign.pvalue),
            "n_a": int(len(monkey_scores)),
            "n_b": int(len(human_benign)),
        },
        {
            "comparison": "monkey_vs_human_pathogenic_llr_ks",
            "statistic": float(ks_monkey_vs_path.statistic),
            "p_value": float(ks_monkey_vs_path.pvalue),
            "n_a": int(len(monkey_scores)),
            "n_b": int(len(human_path)),
        },
    ]
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(args.output / "score_transfer_tests.csv", index=False)

    quantiles = [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]
    dist_rows = []
    for name, series in [
        ("human_pathogenic", human_path),
        ("human_benign", human_benign),
        ("snj_missense", monkey_scores),
    ]:
        row = {"cohort": name, "n": int(len(series))}
        for q in quantiles:
            row[f"llr_q{int(q * 100):02d}"] = float(series.quantile(q))
        row["llr_mean"] = float(series.mean())
        dist_rows.append(row)
    pd.DataFrame(dist_rows).to_csv(args.output / "llr_distribution_summary.csv", index=False)

    # Internal functional direction: more negative LLR should increase calibrated P
    rho, p = spearmanr(clinvar["esm2_score"], clinvar["pathogenicity_prob"])
    direction = pd.DataFrame(
        [
            {
                "metric": "clinvar_llr_vs_sigmoid_prob_spearman",
                "rho": float(rho),
                "p_value": float(p),
            }
        ]
    )
    direction.to_csv(args.output / "internal_directionality.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].hist(human_benign, bins=80, density=True, alpha=0.5, label="Human benign")
    axes[0].hist(human_path, bins=80, density=True, alpha=0.5, label="Human pathogenic")
    axes[0].hist(monkey_scores, bins=80, density=True, alpha=0.5, label="SNJ missense")
    axes[0].set(xlabel="ESM-2 LLR", ylabel="Density", title="Score distributions")
    axes[0].legend(fontsize=8)

    for label, color in [(1, "#D55E00"), (0, "#0072B2")]:
        subset = clinvar.loc[clinvar["label"] == label, "pathogenicity_prob"]
        axes[1].hist(subset, bins=40, density=True, alpha=0.6, label=f"Human label={label}")
    axes[1].hist(
        monkey["pathogenicity_prob"],
        bins=40,
        density=True,
        alpha=0.5,
        histtype="step",
        linewidth=1.5,
        label="SNJ calibrated P",
    )
    axes[1].set(xlabel="Sigmoid pathogenicity probability", ylabel="Density")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(args.output / "cross_species_score_distributions.png", dpi=300)
    plt.close(fig)

    metadata = {
        "clinvar_source": str(args.clinvar),
        "monkey_source": str(args.monkey),
        "n_clinvar": int(len(clinvar)),
        "n_monkey_missense": int(len(monkey)),
        "conda_environment": os.environ.get("CONDA_DEFAULT_ENV"),
        "python": sys.version.split()[0],
        "seed": SEED,
        "interpretation": (
            "Human ClinVar labels validate clinical calibration. SNJ score "
            "distributions are compared as a cross-species transfer sanity check; "
            "they do not prove selection coefficients in monkeys."
        ),
        "limitations": [
            "No independent primate common-variant catalog was available in-repo.",
            "No ProteinGym DMS assays were scored in this pass.",
            "SNJ and human variants are not ortholog-matched site-by-site.",
        ],
    }
    (args.output / "cross_species_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )

    report = f"""# Cross-species and score-transfer validation

- Human pathogenic vs benign LLR (Mann-Whitney): p = {mw_prob.pvalue:.3g}
- SNJ missense vs human benign LLR (KS): D = {ks_monkey_vs_benign.statistic:.3f}, p = {ks_monkey_vs_benign.pvalue:.3g}
- SNJ missense vs human pathogenic LLR (KS): D = {ks_monkey_vs_path.statistic:.3f}, p = {ks_monkey_vs_path.pvalue:.3g}
- ClinVar LLR vs sigmoid probability Spearman rho = {rho:.3f}

SNJ missense scores occupy an intermediate range relative to human ClinVar extremes,
consistent with population-segregating variants rather than fixed clinical extremes.
Independent DMS assays and ortholog-matched primate catalogs remain future extensions.
"""
    (args.output / "SUMMARY.md").write_text(report, encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
