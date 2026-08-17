#!/usr/bin/env python3
"""Run genetic-load sensitivity analyses from one annotated-VCF scan."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from load_sensitivity_common import (
    VariantMatrix,
    association_rows,
    calculate_load,
    load_predictions,
    nearest_maf_callrate_match,
    scan_vcf_once,
    scenario_concordance,
    summarize_matching,
    write_metadata,
)


BASE = get_base_dir()
REFERENCE_SCENARIO = "main_all_scored_continuous_h0.25_lof0.95"
SEED = 20260710


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Single-scan individual genetic-load sensitivity analysis"
    )
    parser.add_argument(
        "--vcf",
        type=Path,
        default=BASE
        / "output/phase2_annotation/snpeff_annotation/annotated_variants.vcf.gz",
    )
    parser.add_argument(
        "--predictions",
        type=Path,
        default=BASE
        / "output/phase4_plm_predictions/ensemble/ensemble_predictions.csv",
    )
    parser.add_argument(
        "--froh",
        type=Path,
        default=BASE
        / "output/phase3a_population_genomics/roh_analysis/roh_summary_per_individual.csv",
    )
    parser.add_argument(
        "--phase5-load",
        type=Path,
        default=BASE
        / "output/phase5_genetic_load/individual_load/individual_genetic_load.csv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=BASE / "output/method_validation/load_sensitivity",
    )
    parser.add_argument("--negative-control-replicates", type=int, default=100)
    parser.add_argument("--seed", type=int, default=SEED)
    return parser.parse_args()


def burden_scenarios(matrix: VariantMatrix, samples: list[str]) -> pd.DataFrame:
    outputs: list[pd.DataFrame] = []

    def add(name: str, **kwargs) -> None:
        outputs.append(calculate_load(matrix, samples, scenario=name, **kwargs))

    add(REFERENCE_SCENARIO, h=0.25, lof_weight=0.95)

    for h in (0.0, 0.1, 0.25, 0.5):
        add(f"dominance_h{h:g}", h=h, lof_weight=0.95)
    for lof_weight in (0.5, 0.75, 0.95, 1.0):
        add(
            f"lof_weight_{lof_weight:g}",
            h=0.25,
            lof_weight=lof_weight,
        )

    for percent in (5, 10, 15):
        threshold = float(np.quantile(matrix.probabilities, 1.0 - percent / 100.0))
        add(
            f"top_{percent}pct_pathogenicity",
            h=0.25,
            lof_weight=0.95,
            mask=matrix.probabilities >= threshold,
        )
    for threshold in (0.5, 0.8):
        add(
            f"pathogenicity_p_gt_{threshold:g}",
            h=0.25,
            lof_weight=0.95,
            mask=matrix.probabilities > threshold,
        )

    add(
        "lof_only",
        h=0.25,
        lof_weight=0.95,
        mask=matrix.kinds == "lof",
    )
    add(
        "missense_only_all_scored",
        h=0.25,
        mask=matrix.kinds == "missense",
    )
    add(
        "count_carrier_dominance_h0.25",
        h=0.25,
        weights=np.ones(len(matrix.variant_ids)),
    )
    add(
        "count_alt_alleles",
        weights=np.ones(len(matrix.variant_ids)),
        alt_allele_count=True,
    )
    return pd.concat(outputs, ignore_index=True)


def negative_controls(
    burden: VariantMatrix,
    synonymous: VariantMatrix,
    samples: list[str],
    replicates: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    load_frames = []
    matching_rows = []
    target_weights = burden.probabilities.copy()
    target_weights[burden.kinds == "lof"] = 0.95
    for replicate in range(1, replicates + 1):
        target_indices, control_indices = nearest_maf_callrate_match(
            burden, synonymous, rng
        )
        if len(target_indices) == 0:
            raise RuntimeError("No MAF/call-rate matched synonymous controls found")
        control_matrix = synonymous.subset(control_indices)
        load_frames.append(
            calculate_load(
                control_matrix,
                samples,
                scenario=f"negative_synonymous_maf_matched_{replicate:03d}",
                h=0.25,
                weights=target_weights[target_indices],
            )
        )
        matching_rows.append(
            {
                "replicate": replicate,
                **summarize_matching(
                    burden,
                    synonymous,
                    target_indices,
                    control_indices,
                ),
            }
        )
    return pd.concat(load_frames, ignore_index=True), pd.DataFrame(matching_rows)


def phase5_concordance(
    loads: pd.DataFrame, phase5_path: Path
) -> pd.DataFrame:
    from scipy.stats import pearsonr, spearmanr

    if not phase5_path.exists():
        return pd.DataFrame()
    old = pd.read_csv(phase5_path)
    main = loads.loc[
        loads["scenario"] == REFERENCE_SCENARIO,
        ["IID", "raw_load", "load_per_1000_callable"],
    ]
    merged = main.merge(
        old[["IID", "Realized_Load", "Total_Genetic_Load"]],
        on="IID",
        how="inner",
    )
    rows = []
    for new_metric in ("raw_load", "load_per_1000_callable"):
        for old_metric in ("Realized_Load", "Total_Genetic_Load"):
            pearson = pearsonr(merged[new_metric], merged[old_metric])
            spearman = spearmanr(merged[new_metric], merged[old_metric])
            rows.append(
                {
                    "new_metric": new_metric,
                    "phase5_metric": old_metric,
                    "n": len(merged),
                    "pearson_r": pearson.statistic,
                    "pearson_p": pearson.pvalue,
                    "spearman_rho": spearman.statistic,
                    "spearman_p": spearman.pvalue,
                }
            )
    return pd.DataFrame(rows)


def make_plots(
    scenario_loads: pd.DataFrame,
    concordance: pd.DataFrame,
    associations: pd.DataFrame,
    output: Path,
) -> None:
    non_null = concordance[
        ~concordance["scenario"].str.startswith("negative_")
    ].sort_values("spearman_rho")
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.barh(non_null["scenario"], non_null["spearman_rho"], color="#0072B2")
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set(
        xlabel=f"Spearman rank correlation with {REFERENCE_SCENARIO}",
        ylabel="Sensitivity scenario",
        xlim=(-1, 1),
    )
    fig.tight_layout()
    fig.savefig(output / "scenario_rank_concordance.png", dpi=300)
    plt.close(fig)

    pivot = scenario_loads[
        ~scenario_loads["scenario"].str.startswith("negative_")
    ].pivot(index="IID", columns="scenario", values="z_load")
    fig, ax = plt.subplots(figsize=(12, 10))
    sns.heatmap(pivot.corr(method="spearman"), cmap="vlag", center=0, ax=ax)
    ax.set_title("Spearman concordance of standardized load scenarios")
    fig.tight_layout()
    fig.savefig(output / "scenario_correlation_heatmap.png", dpi=300)
    plt.close(fig)

    froh = associations[associations["external_metric"] == "F_ROH"].copy()
    null = froh[froh["scenario"].str.startswith("negative_")]["spearman_rho"]
    observed = froh.loc[
        froh["scenario"] == REFERENCE_SCENARIO, "spearman_rho"
    ]
    if len(null) and len(observed):
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.hist(null, bins=20, color="#999999", alpha=0.8)
        ax.axvline(
            observed.iloc[0],
            color="#D55E00",
            linewidth=2,
            label=f"Observed = {observed.iloc[0]:.3f}",
        )
        ax.set(
            xlabel="Spearman correlation with F_ROH",
            ylabel="Matched-control replicates",
            title="MAF/call-rate matched synonymous negative control",
        )
        ax.legend()
        fig.tight_layout()
        fig.savefig(output / "negative_control_froh.png", dpi=300)
        plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    print(f"Loading predictions: {args.predictions}", flush=True)
    predictions = load_predictions(args.predictions)
    print(f"Single VCF scan: {args.vcf}", flush=True)
    samples, burden, synonymous, qc = scan_vcf_once(args.vcf, predictions)
    if not samples or len(burden.variant_ids) == 0:
        raise RuntimeError("No samples or scored burden variants were obtained")
    print(
        f"Retained {len(burden.variant_ids)} burden and "
        f"{len(synonymous.variant_ids)} synonymous variants",
        flush=True,
    )

    scenario_loads = burden_scenarios(burden, samples)
    null_loads, matching = negative_controls(
        burden,
        synonymous,
        samples,
        args.negative_control_replicates,
        args.seed,
    )
    all_loads = pd.concat([scenario_loads, null_loads], ignore_index=True)
    all_loads.to_csv(args.output / "individual_load_sensitivity.csv", index=False)
    matching.to_csv(args.output / "negative_control_matching_qc.csv", index=False)

    covariates = pd.read_csv(args.froh, usecols=["IID", "F_ROH"])
    associations = association_rows(all_loads, covariates)
    associations.to_csv(args.output / "froh_associations.csv", index=False)
    concordance = scenario_concordance(all_loads, REFERENCE_SCENARIO)
    concordance.to_csv(args.output / "scenario_concordance.csv", index=False)

    scenario_summary = (
        all_loads.groupby("scenario", sort=False)
        .agg(
            n_individuals=("IID", "nunique"),
            n_variants=("n_variants", "first"),
            mean_raw_load=("raw_load", "mean"),
            sd_raw_load=("raw_load", "std"),
            mean_callable_loci=("callable_loci", "mean"),
            mean_load_per_1000=("load_per_1000_callable", "mean"),
            sd_load_per_1000=("load_per_1000_callable", "std"),
        )
        .reset_index()
    )
    scenario_summary.to_csv(args.output / "scenario_summary.csv", index=False)
    old_concordance = phase5_concordance(scenario_loads, args.phase5_load)
    old_concordance.to_csv(args.output / "phase5_concordance.csv", index=False)

    observed_row = associations[
        associations["scenario"] == REFERENCE_SCENARIO
    ].iloc[0]
    null_rho = associations[
        associations["scenario"].str.startswith("negative_")
    ]["spearman_rho"].to_numpy()
    observed_rho = float(observed_row["spearman_rho"])
    empirical_p = float(
        (1 + np.sum(np.abs(null_rho) >= abs(observed_rho)))
        / (1 + len(null_rho))
    )
    metadata = {
        "inputs": {
            "vcf": str(args.vcf),
            "predictions": str(args.predictions),
            "froh": str(args.froh),
        },
        "seed": args.seed,
        "n_samples": len(samples),
        "prediction_records": len(predictions),
        "vcf_single_pass": True,
        "vcf_qc": qc,
        "definitions": {
            "main": (
                "All VCF LoF plus every missense variant having a finite ESM-2 "
                "score; continuous sigmoid weights, h=0.25, LoF weight=0.95."
            ),
            "realized_load": "sum(weight*h) for heterozygotes plus sum(weight) for alternate homozygotes",
            "standardized_load": (
                "raw load divided by per-individual callable scenario loci, "
                "scaled to 1,000 loci, then z-standardized across individuals"
            ),
            "top_percent": "upper x% of pathogenicity probabilities across retained LoF and scored missense loci",
            "negative_control": (
                "Synonymous loci matched without replacement to target loci by "
                "MAF (caliper 0.02) and call rate (caliper 0.05); target weights "
                "are transferred to matched synonymous genotypes."
            ),
        },
        "negative_control_replicates": args.negative_control_replicates,
        "negative_control_match_qc_mean": matching.mean(numeric_only=True).to_dict(),
        "main_froh_spearman_rho": observed_rho,
        "negative_control_two_sided_empirical_p": empirical_p,
    }
    write_metadata(args.output / "analysis_metadata.json", metadata)
    make_plots(all_loads, concordance, associations, args.output)

    non_null_concordance = concordance[
        ~concordance["scenario"].str.startswith("negative_")
    ]
    minimum = non_null_concordance.loc[
        non_null_concordance["spearman_rho"].idxmin()
    ]
    summary = f"""# Load sensitivity analysis

- Samples: {len(samples)}
- Single-pass VCF records: {qc['vcf_records']:,}
- Retained burden loci: {len(burden.variant_ids):,} ({qc['lof']:,} LoF; {qc['scored_missense']:,} scored missense)
- Synonymous control pool: {len(synonymous.variant_ids):,}
- Main standardized load vs F_ROH: Spearman rho = {observed_rho:.3f}, p = {observed_row['spearman_p']:.3g}
- MAF/call-rate matched synonymous null: {args.negative_control_replicates} replicates; two-sided empirical p = {empirical_p:.3g}
- Lowest rank agreement among requested non-null sensitivities: {minimum['scenario']} (rho = {minimum['spearman_rho']:.3f})
- Mean matched fraction: {matching['match_fraction'].mean():.3f}; mean absolute MAF difference: {matching['mean_absolute_maf_difference'].mean():.5f}

Loads are reported raw, per 1,000 callable loci, and as across-individual z-scores.
See `analysis_metadata.json` for exact scenario and negative-control definitions.
"""
    (args.output / "SUMMARY.md").write_text(summary, encoding="utf-8")
    print(json.dumps(metadata, indent=2), flush=True)


if __name__ == "__main__":
    main()
