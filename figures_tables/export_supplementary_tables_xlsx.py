#!/usr/bin/env python3
"""Merge supplementary tables S1–S12 into one Excel workbook (one sheet per table)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from lib.paths import REPO, TABLES

OUT = TABLES / "supplementary_tables.xlsx"


def read_csv(path: Path, **kwargs) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    sep = "\t" if path.suffix.lower() == ".tsv" else kwargs.pop("sep", ",")
    # Detect tab-separated CSVs that use a .csv suffix
    sample = path.read_text(encoding="utf-8", errors="replace")[:4096]
    if path.suffix.lower() == ".csv" and sample.count("\t") > sample.count(","):
        sep = "\t"
    return pd.read_csv(path, sep=sep, **kwargs)


def write_sections(writer: pd.ExcelWriter, sheet: str, sections: list[tuple[str, pd.DataFrame]]) -> None:
    """Write one or more titled blocks into a single sheet."""
    start = 0
    for title, frame in sections:
        if frame is None or frame.empty:
            continue
        pd.DataFrame([[title]]).to_excel(
            writer, sheet_name=sheet, startrow=start, header=False, index=False
        )
        frame.to_excel(writer, sheet_name=sheet, startrow=start + 1, index=False)
        start += len(frame) + 3


def table_s5_validation() -> list[tuple[str, pd.DataFrame]]:
    a = pd.DataFrame(
        [
            [
                "External calibration",
                "Held-out ClinVar evaluation: 13,147 variants from 505 proteins",
                "ROC-AUC 0.875; PR-AUC 0.895; balanced accuracy 0.813",
                "The fixed sigmoid discriminates known human clinical classes in proteins excluded from fitting.",
            ],
            [
                "Calibration reproducibility",
                "Five-fold cross-validation with proteins separated among folds",
                "ROC-AUC 0.875 ± 0.010 (mean ± s.d.)",
                "Discrimination is stable across protein partitions.",
            ],
            [
                "Cross-species transfer",
                "SNJ LLR distribution relative to ClinVar",
                "KS D = 0.072 versus benign; 0.650 versus pathogenic",
                "SNJ variants occupy an intermediate range.",
            ],
            [
                "Derived allele-frequency trend",
                "Probability versus polarized SNJ DAF",
                "Spearman rho = 0.020; P = 0.045 (n = 9,701)",
                "Weak trend after ancestral-state polarization.",
            ],
            [
                "Derived-homozygote depletion",
                "Observed/expected derived-homozygote ratio by probability bin",
                "Stronger depletion at higher predicted pathogenicity",
                "Consistent with purifying selection against predicted deleterious alleles.",
            ],
            [
                "Load sensitivity",
                "Rank concordance across alternative models",
                "Minimum non-null rho = 0.611",
                "Rankings are not determined by one parameter setting.",
            ],
            [
                "Synonymous control",
                "100 allele-frequency/call-rate-matched replicates",
                "Empirical two-sided P = 0.0099",
                "Main association is unlikely to arise from allele frequency, call rate or variant count alone.",
            ],
        ],
        columns=["Validation layer", "Analysis", "Result", "Interpretation"],
    )
    b = pd.DataFrame(
        [
            [
                "GRM mixed model",
                "Total load per s.d., conditional on sex and F_ROH",
                "β = 0.716; s.e. = 0.274",
                0.0088,
                "Higher load is associated with higher CHS after relatedness correction.",
            ],
            [
                "GRM mixed model",
                "F_ROH per s.d., conditional on sex and load",
                "β = −0.080; s.e. = 0.264",
                0.761,
                "Autozygosity adds no independent CHS association.",
            ],
            [
                "Freedman–Lane",
                "Load conditional on F_ROH; 4,999 permutations",
                "partial F = 6.454",
                0.0152,
                "Supports an incremental load association.",
            ],
            [
                "Freedman–Lane",
                "F_ROH conditional on load; 4,999 permutations",
                "partial F = 0.087",
                0.782,
                "Does not support incremental F_ROH.",
            ],
            [
                "Grouped cross-validation",
                "Covariates + load; five folds, n = 68",
                "RMSE = 2.295",
                None,
                "Lowest RMSE among the four prespecified models.",
            ],
        ],
        columns=["Analysis", "Model or comparison", "Estimate/performance", "P value", "Interpretation"],
    )
    return [
        ("a | External calibration, population-genetic consistency and sensitivity", a),
        ("b | Relatedness-aware morbidity criterion", b),
    ]


def table_s10_evidence() -> pd.DataFrame:
    return pd.DataFrame(
        [
            [
                "Morbidity criterion",
                "Load, especially heterozygous/potential load, predicts CHS after kinship correction; F_ROH does not",
                "Health tracks shared segregating burden rather than individual autozygosity",
                "Predicts stronger F_ROH and homozygous-load effects",
                "Fixation would reduce among-individual segregating-burden variation",
            ],
            [
                "ROH architecture",
                "Short ROH dominate; no ROH >5 Mb",
                "Autozygosity is historical rather than recent close-kin mating",
                "Predicts long young tracts",
                "Largely neutral",
            ],
            [
                "Deleterious AFS",
                "Intermediate-frequency enrichment; no fixed PLM-flagged missense variants",
                "Alleles are elevated but remain segregating",
                "Predicts a stronger rare-recessive/long-ROH pattern",
                "Predicts more fixed deleterious alleles",
            ],
            [
                "Demographic inference",
                "Bottleneck and recent contraction; reduced-size episodes are shorter than 4Ne",
                "History can shift frequencies before expected fixation",
                "Only partly explains the inferred timing",
                "Inconsistent with a long constant-small population",
            ],
            [
                "Forward simulation",
                "Inferred Ne(t) yields near-zero fixation with nonzero intermediate alleles; constant Ne = 150 yields higher fixation and 0% intermediate",
                "Provides qualitative mechanistic plausibility",
                "Constant recent-inbreeding framing does not reproduce the spectrum",
                "Constant-small null differs qualitatively",
            ],
        ],
        columns=[
            "Evidence stream",
            "Key observation",
            "Support for transient drift load",
            "Limitation of recent-inbreeding explanation",
            "Limitation of fixation-only explanation",
        ],
    )


def main() -> None:
    """Write sheets in manuscript first-citation order (S1–S12)."""
    TABLES.mkdir(parents=True, exist_ok=True)

    s1 = read_csv(TABLES / "supplementary_table_s1_diversity.csv")
    s2 = read_csv(TABLES / "supplementary_table_s2_roh.csv")

    s3_pred = read_csv(REPO / "output/phase4_plm_predictions/esm2/esm2_predictions.csv")
    s3 = s3_pred.drop(
        columns=[c for c in s3_pred.columns if c.lower() in {"wt_sequence", "sequence"}],
        errors="ignore",
    )

    s4 = read_csv(TABLES / "supplementary_table_s4_load.csv")

    s6_corr = read_csv(TABLES / "supplementary_table_s6_associations.csv")
    s6_cc_path = REPO / "output/phenotype_genotype_analysis/case_control_analysis.csv"
    s6_sections = [("Rank correlations", s6_corr)]
    if s6_cc_path.exists():
        s6_sections.append(("Case-control analysis", read_csv(s6_cc_path)))

    s7_sum = read_csv(TABLES / "supplementary_table_s7_daf_summary.csv")
    s7_lof = read_csv(TABLES / "supplementary_table_s7_lof_daf.csv")
    s7_mis = read_csv(TABLES / "supplementary_table_s7_missense_daf.csv")

    s8 = read_csv(TABLES / "supplementary_table_s8_model_comparison.csv")
    s9_mle = read_csv(TABLES / "supplementary_table_s9_demographic_mle.csv")
    s9_boot = read_csv(TABLES / "supplementary_table_s9_demographic_bootstrap.csv")

    # S11: LD-region annotation (former S12); GWAS top-hit sheet removed
    s11_sections: list[tuple[str, pd.DataFrame]] = []
    s11_ann = REPO / "output/phenotype_genotype_analysis/gwas/gwas_ld_region_roh_missense_annotation.csv"
    s11_var = REPO / "output/phenotype_genotype_analysis/gwas/gwas_ld_region_missense_variants.csv"
    if s11_ann.exists():
        s11_sections.append(("LD-region ROH and missense annotation", read_csv(s11_ann)))
    if s11_var.exists():
        s11_sections.append(("LD-region linked missense variants", read_csv(s11_var)))

    # S12: functional annotation (former S13)
    s12_sections: list[tuple[str, pd.DataFrame]] = []
    for title, path in [
        (
            "Variant counts by impact",
            REPO / "output/phase2_annotation/functional_annotation/v2_variant_counts_by_impact.tsv",
        ),
        (
            "Variant counts by effect",
            REPO / "output/phase2_annotation/functional_annotation/v2_variant_counts_by_effect.tsv",
        ),
        (
            "LoF variants",
            REPO / "output/phase2_annotation/snpeff_annotation/lof_variants.csv",
        ),
        (
            "High-impact variants",
            REPO / "output/phase2_annotation/snpeff_annotation/high_impact_variants.csv",
        ),
    ]:
        if path.exists():
            s12_sections.append((title, read_csv(path)))

    with pd.ExcelWriter(OUT, engine="openpyxl") as writer:
        s1.to_excel(writer, sheet_name="Table S1", index=False)
        s2.to_excel(writer, sheet_name="Table S2", index=False)
        s3.to_excel(writer, sheet_name="Table S3", index=False)
        s4.to_excel(writer, sheet_name="Table S4", index=False)
        write_sections(writer, "Table S5", table_s5_validation())
        write_sections(writer, "Table S6", s6_sections)
        write_sections(
            writer,
            "Table S7",
            [
                ("DAF class summary", s7_sum),
                ("LoF polarized DAF (variant-level)", s7_lof),
                ("Top-decile missense polarized DAF (variant-level)", s7_mis),
            ],
        )
        s8.to_excel(writer, sheet_name="Table S8", index=False)
        write_sections(
            writer,
            "Table S9",
            [
                ("Maximum-likelihood estimates and bootstrap CIs", s9_mle),
                ("Parametric bootstrap replicates", s9_boot),
            ],
        )
        table_s10_evidence().to_excel(writer, sheet_name="Table S10", index=False)
        if s11_sections:
            write_sections(writer, "Table S11", s11_sections)
        else:
            pd.DataFrame({"note": ["S11 source files not found"]}).to_excel(
                writer, sheet_name="Table S11", index=False
            )
        if s12_sections:
            write_sections(writer, "Table S12", s12_sections)
        else:
            pd.DataFrame({"note": ["S12 source files not found"]}).to_excel(
                writer, sheet_name="Table S12", index=False
            )

    print(f"Wrote {OUT}")
    xl = pd.ExcelFile(OUT)
    print("Sheets:", xl.sheet_names)
    for name in xl.sheet_names:
        df = pd.read_excel(OUT, sheet_name=name, header=None)
        print(f"  {name}: {df.shape[0]} rows x {df.shape[1]} cols")


if __name__ == "__main__":
    main()
