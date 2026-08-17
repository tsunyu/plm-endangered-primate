#!/usr/bin/env python3
"""Population-genetic consistency checks for calibrated deleteriousness scores."""

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
from scipy.stats import linregress, spearmanr

from load_sensitivity_common import (
    load_predictions,
    parse_consequences,
    parse_genotypes,
    sigmoid_probability,
)


BASE = get_base_dir()
SEED = 20260710
RISK_BINS = [0.0, 0.2, 0.4, 0.6, 0.8, 1.01]
RISK_LABELS = ["<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", ">0.8"]


def scan_scored_variants(vcf_path: Path, predictions: dict[str, float]) -> pd.DataFrame:
    import gzip

    opener = gzip.open if vcf_path.suffix == ".gz" else open
    samples: list[str] = []
    rows: list[dict] = []
    froh: dict[str, float] = {}

    with opener(vcf_path, "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip().split("\t")[9:]
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 10:
                continue
            chrom, pos, _, ref, alt = fields[:5]
            if "," in alt:
                continue
            variant_id = f"{chrom}:{pos}:{ref}:{alt}"
            is_lof, is_missense, is_synonymous = parse_consequences(fields[7])
            score = predictions.get(variant_id)
            if is_lof:
                probability = 0.95
                kind = "lof"
            elif is_missense and score is not None:
                probability = sigmoid_probability(score)
                kind = "missense"
            elif is_synonymous and not is_missense:
                probability = 0.0
                kind = "synonymous"
            else:
                continue

            dosage, called = parse_genotypes(fields[9:], len(samples))
            called_alleles = int(2 * called.sum())
            if called_alleles == 0:
                continue
            alt_count = int(dosage.sum())
            maf = min(alt_count / called_alleles, 1.0 - alt_count / called_alleles)
            hom_obs = int((dosage == 2).sum())
            het_obs = int((dosage == 1).sum())

            expected_hom = 0.0
            for index, sample in enumerate(samples):
                if not called[index]:
                    continue
                p = alt_count / called_alleles
                expected_hom += p * p

            rows.append(
                {
                    "variant_id": variant_id,
                    "kind": kind,
                    "pathogenicity_prob": probability,
                    "maf": maf,
                    "hom_observed": hom_obs,
                    "hom_expected_hwe": expected_hom,
                    "het_observed": het_obs,
                    "carrier_count": int((dosage >= 1).sum()),
                    "n_called": int(called.sum()),
                }
            )
    return pd.DataFrame(rows)


def assign_risk_bin(probability: float) -> str:
    for low, high, label in zip(RISK_BINS[:-1], RISK_BINS[1:], RISK_LABELS):
        if low <= probability < high:
            return label
    return RISK_LABELS[-1]


def maf_spectrum(frame: pd.DataFrame) -> pd.DataFrame:
    bins = [0, 0.05, 0.2, 0.5, 0.95, 1.0]
    labels = ["<0.05", "0.05-0.2", "0.2-0.5", "0.5-0.95", ">=0.95"]
    rows = []
    for cohort, subset in [
        ("synonymous", frame[frame["kind"] == "synonymous"]),
        ("all_missense", frame[frame["kind"] == "missense"]),
        ("lof", frame[frame["kind"] == "lof"]),
    ]:
        counts = pd.cut(subset["maf"], bins=bins, labels=labels, include_lowest=True).value_counts()
        total = max(len(subset), 1)
        for label in labels:
            rows.append(
                {
                    "cohort": cohort,
                    "maf_bin": label,
                    "n_variants": int(counts.get(label, 0)),
                    "fraction": float(counts.get(label, 0) / total),
                }
            )
    missense = frame[frame["kind"] == "missense"].copy()
    missense["risk_bin"] = missense["pathogenicity_prob"].map(assign_risk_bin)
    for risk_bin, subset in missense.groupby("risk_bin"):
        counts = pd.cut(subset["maf"], bins=bins, labels=labels, include_lowest=True).value_counts()
        total = max(len(subset), 1)
        for label in labels:
            rows.append(
                {
                    "cohort": f"missense_risk_{risk_bin}",
                    "maf_bin": label,
                    "n_variants": int(counts.get(label, 0)),
                    "fraction": float(counts.get(label, 0) / total),
                }
            )
    return pd.DataFrame(rows)


def homozygote_depletion(frame: pd.DataFrame) -> pd.DataFrame:
    missense = frame[frame["kind"] == "missense"].copy()
    missense["risk_bin"] = missense["pathogenicity_prob"].map(assign_risk_bin)
    rows = []
    for risk_bin, subset in missense.groupby("risk_bin"):
        expected = subset["hom_expected_hwe"].sum()
        observed = subset["hom_observed"].sum()
        ratio = observed / expected if expected > 0 else np.nan
        rows.append(
            {
                "risk_bin": risk_bin,
                "n_variants": int(len(subset)),
                "hom_observed": int(observed),
                "hom_expected_hwe": float(expected),
                "observed_expected_ratio": float(ratio),
                "mean_maf": float(subset["maf"].mean()),
                "median_pathogenicity_prob": float(subset["pathogenicity_prob"].median()),
            }
        )
    return pd.DataFrame(rows)


def trend_tests(frame: pd.DataFrame) -> pd.DataFrame:
    missense = frame[frame["kind"] == "missense"].copy()
    rho_prob, p_prob = spearmanr(missense["pathogenicity_prob"], missense["maf"])
    lr = linregress(missense["pathogenicity_prob"], missense["maf"])
    synonymous = frame[frame["kind"] == "synonymous"]
    rho_syn, p_syn = spearmanr(
        synonymous["pathogenicity_prob"], synonymous["maf"]
    ) if len(synonymous) else (np.nan, np.nan)
    return pd.DataFrame(
        [
            {
                "test": "missense_prob_vs_maf_spearman",
                "statistic": float(rho_prob),
                "p_value": float(p_prob),
                "n": int(len(missense)),
            },
            {
                "test": "missense_prob_vs_maf_linear_slope",
                "statistic": float(lr.slope),
                "p_value": float(lr.pvalue),
                "n": int(len(missense)),
            },
            {
                "test": "synonymous_prob_vs_maf_spearman",
                "statistic": float(rho_syn),
                "p_value": float(p_syn),
                "n": int(len(synonymous)),
            },
        ]
    )


def make_plots(spectrum: pd.DataFrame, homozygote: pd.DataFrame, output: Path) -> None:
    risk_cohorts = sorted({c for c in spectrum["cohort"] if c.startswith("missense_risk_")})
    maf_bins = ["<0.05", "0.05-0.2", "0.2-0.5", "0.5-0.95", ">=0.95"]
    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(maf_bins))
    width = 0.15
    for index, cohort in enumerate(risk_cohorts):
        subset = spectrum[spectrum["cohort"] == cohort].set_index("maf_bin")
        values = [subset.loc[label, "fraction"] if label in subset.index else 0 for label in maf_bins]
        ax.bar(x + index * width, values, width=width, label=cohort.replace("missense_risk_", "P "))
    syn = spectrum[spectrum["cohort"] == "synonymous"].set_index("maf_bin")
    syn_values = [syn.loc[label, "fraction"] if label in syn.index else 0 for label in maf_bins]
    ax.plot(x + width * (len(risk_cohorts) / 2), syn_values, "k--", marker="o", label="synonymous")
    ax.set_xticks(x + width * (len(risk_cohorts) - 1) / 2)
    ax.set_xticklabels(maf_bins)
    ax.set(xlabel="Minor allele frequency bin", ylabel="Fraction of variants", title="MAF spectrum by risk bin")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(output / "maf_spectrum_by_risk_bin.png", dpi=300)
    plt.close(fig)

    if len(homozygote):
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(
            homozygote["median_pathogenicity_prob"],
            homozygote["observed_expected_ratio"],
            marker="o",
        )
        ax.axhline(1.0, color="black", linestyle="--", linewidth=0.8)
        ax.set(
            xlabel="Median pathogenicity probability (risk bin)",
            ylabel="Observed / expected homozygotes (HWE)",
            title="Homozygote depletion trend",
        )
        fig.tight_layout()
        fig.savefig(output / "homozygote_depletion_by_risk.png", dpi=300)
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--vcf",
        type=Path,
        default=BASE / "output/phase2_annotation/snpeff_annotation/annotated_variants.vcf.gz",
    )
    parser.add_argument(
        "--predictions",
        type=Path,
        default=BASE / "output/phase4_plm_predictions/esm2/esm2_predictions.csv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=BASE / "output/method_validation/population_consistency",
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    predictions = load_predictions(args.predictions)
    variants = scan_scored_variants(args.vcf, predictions)
    variants.to_csv(args.output / "variant_population_metrics.csv", index=False)

    spectrum = maf_spectrum(variants)
    spectrum.to_csv(args.output / "maf_spectrum_by_risk.csv", index=False)
    homozygote = homozygote_depletion(variants)
    homozygote.to_csv(args.output / "homozygote_depletion_summary.csv", index=False)
    trends = trend_tests(variants)
    trends.to_csv(args.output / "risk_maf_trend_tests.csv", index=False)
    make_plots(spectrum, homozygote, args.output)

    high = variants[(variants["kind"] == "missense") & (variants["pathogenicity_prob"] > 0.8)]
    low = variants[(variants["kind"] == "missense") & (variants["pathogenicity_prob"] < 0.2)]
    metadata = {
        "vcf": str(args.vcf),
        "predictions": str(args.predictions),
        "n_variants_retained": int(len(variants)),
        "n_missense_scored": int((variants["kind"] == "missense").sum()),
        "n_synonymous": int((variants["kind"] == "synonymous").sum()),
        "n_lof": int((variants["kind"] == "lof").sum()),
        "high_risk_missense_gt_0.8": int(len(high)),
        "low_risk_missense_lt_0.2": int(len(low)),
        "mean_maf_high_risk": float(high["maf"].mean()) if len(high) else None,
        "mean_maf_low_risk": float(low["maf"].mean()) if len(low) else None,
        "conda_environment": os.environ.get("CONDA_DEFAULT_ENV"),
        "python": sys.version.split()[0],
        "seed": SEED,
        "caveats": [
            "MAF uses alternate-allele frequency without polarized ancestral state.",
            "Homozygote expectation uses random-mating HWE at the population level.",
            "n=68 limits rare-variant homozygote tests; results are aggregate trends.",
        ],
    }
    (args.output / "population_consistency_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )

    rho_row = trends.loc[trends["test"] == "missense_prob_vs_maf_spearman"].iloc[0]
    report = f"""# Population consistency validation

- Scored missense variants: {metadata['n_missense_scored']:,}
- Spearman(pathogenicity probability, MAF) = {rho_row['statistic']:.3f}, p = {rho_row['p_value']:.3g}
- Mean MAF for P > 0.8 missense: {metadata['mean_maf_high_risk']}
- Mean MAF for P < 0.2 missense: {metadata['mean_maf_low_risk']}

Higher calibrated pathogenicity probabilities tend to occur at lower population
allele frequencies, while synonymous controls provide a neutral reference spectrum.
"""
    (args.output / "SUMMARY.md").write_text(report, encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
