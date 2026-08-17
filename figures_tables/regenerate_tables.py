#!/usr/bin/env python3
"""Regenerate supplementary tables S1–S12 used in the paper.

Numbering follows first-citation order:
  S1 diversity, S2 ROH, S4 load, S6 associations, S7 DAF,
  S8 model comparison, S9 demographic estimates.
"""

from __future__ import annotations

import shutil

import pandas as pd

from lib.daf import polarized_subset, scan_variants, summarize_daf
from lib.paths import (
    BOOTSTRAP_RESULTS,
    MODEL_COMPARISON,
    PARAM_ESTIMATES,
    REPO,
    TABLES,
)


def write_s8_model() -> None:
    shutil.copy2(MODEL_COMPARISON, TABLES / "supplementary_table_s8_model_comparison.csv")


def write_s9_demography() -> None:
    mle = pd.read_csv(PARAM_ESTIMATES)
    boot = pd.read_csv(BOOTSTRAP_RESULTS)
    mle.to_csv(TABLES / "supplementary_table_s9_demographic_mle.csv", index=False)
    boot.to_csv(TABLES / "supplementary_table_s9_demographic_bootstrap.csv", index=False)


def write_s7_daf() -> None:
    frame = scan_variants()
    polarized = polarized_subset(frame)
    lof = polarized[polarized["kind"] == "lof"][
        ["variant_id", "chrom", "pos", "ref", "alt", "ancestral", "daf", "polarization_status"]
    ]
    missense = polarized[
        (polarized["kind"] == "missense") & polarized["is_deleterious_missense"]
    ][
        [
            "variant_id",
            "chrom",
            "pos",
            "ref",
            "alt",
            "ancestral",
            "daf",
            "pathogenicity_prob",
            "polarization_status",
        ]
    ]
    lof.to_csv(TABLES / "supplementary_table_s7_lof_daf.csv", index=False)
    missense.to_csv(TABLES / "supplementary_table_s7_missense_daf.csv", index=False)
    summary = summarize_daf(frame)
    summary.to_csv(TABLES / "supplementary_table_s7_daf_summary.csv", index=False)


def mirror_existing_tables() -> None:
    mappings = {
        "supplementary_table_s2_roh": REPO
        / "output/phase3a_population_genomics/roh_analysis/roh_individual_stats.csv",
        "supplementary_table_s1_diversity": REPO
        / "output/phase3a_population_genomics/diversity_metrics/heterozygosity_summary.csv",
        "supplementary_table_s4_load": REPO
        / "output/phase5_genetic_load/individual_load/individual_genetic_load.csv",
        "supplementary_table_s6_associations": REPO
        / "output/phenotype_genotype_analysis/correlation_results.csv",
    }
    for stem, src in mappings.items():
        if src.exists():
            shutil.copy2(src, TABLES / f"{stem}.csv")


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    write_s8_model()
    write_s9_demography()
    write_s7_daf()
    mirror_existing_tables()


if __name__ == "__main__":
    main()
