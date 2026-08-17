#!/usr/bin/env python3
"""Redraw supplementary figures S1–S16 with unified publication style."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

from lib.bootstrap import bootstrap_cohens_d_ci, bootstrap_spearman_ci
from lib.daf import (
    RISK_LABELS,
    assign_risk_bin,
    homozygote_depletion_by_bin,
    polarized_subset,
    probability_daf_trend,
    scan_variants,
)
from lib.export import save_figure
from lib.paths import (
    CASE_CONTROL,
    CORRELATION,
    DFE_REPLICATES,
    ESM2_PREDICTIONS,
    FITNESS_CV,
    FITNESS_FIXED,
    FITNESS_PERM,
    HET_SUMMARY,
    INDIVIDUAL_LOAD,
    MERGED_PHENO,
    MODEL_COMPARISON,
    NUC_DIV,
    PCA_RESULTS,
    REPO,
    ROH_SUMMARY,
    TAJIMAS_D,
)
from lib.style import COLORS, PALETTE, add_panel_label, apply_style, despine

GWAS_DIR = REPO / "output/phenotype_genotype_analysis/gwas"
LOAD_SENS = REPO / "output/method_validation/load_sensitivity"
CLINVAR = REPO / "output/method_validation/clinvar"
CROSS = REPO / "output/method_validation/cross_species"
SFS_DIR = REPO / "output/phase3b_fastsimcoal2/sfs"
# Study-wide Bonferroni threshold (339,508 LD-pruned SNPs; Methods)
GWAS_BONFERRONI_P = 1.12e-8


def save_supp(fig, stem: str) -> None:
    save_figure(fig, stem, supplementary=True)


def plot_s1() -> None:
    pca = pd.read_csv(PCA_RESULTS)
    taj = pd.read_csv(TAJIMAS_D) if TAJIMAS_D.exists() else None
    nuc = pd.read_csv(NUC_DIV) if NUC_DIV.exists() else None
    fig, axes = plt.subplots(1, 3, figsize=(7.09, 2.4))
    axes[0].scatter(pca["PC1"], pca["PC2"], s=12, color=COLORS["blue"], alpha=0.8)
    axes[0].set_xlabel("PC1")
    axes[0].set_ylabel("PC2")
    pi_col = None
    if nuc is not None:
        for candidate in ("PI", "pi", "Pi"):
            if candidate in nuc.columns:
                pi_col = candidate
                break
    if pi_col is not None:
        axes[1].hist(
            nuc[pi_col].dropna(),
            bins=40,
            color=COLORS["orange"],
            edgecolor="white",
            linewidth=0.2,
        )
        axes[1].set_xlabel(r"Nucleotide diversity $\pi$ (100-kb windows)")
        axes[1].set_ylabel("Windows")
    else:
        axes[1].text(0.5, 0.5, "Nucleotide diversity\nunavailable", ha="center", va="center")
    if taj is not None and "Tajima_D" in taj.columns:
        axes[2].hist(taj["Tajima_D"].dropna(), bins=30, color=COLORS["teal"], edgecolor="white", linewidth=0.2)
        axes[2].set_xlabel("Tajima's D")
        axes[2].set_ylabel("Windows")
    else:
        pi_path = REPO / "output/phase3a_population_genomics/diversity_metrics/tajimas_d.Tajima.D"
        if pi_path.exists():
            tdf = pd.read_csv(pi_path, sep=r"\s+")
            axes[2].hist(tdf.iloc[:, -1], bins=30, color=COLORS["teal"], edgecolor="white", linewidth=0.2)
            axes[2].set_xlabel("Tajima's D")
    for ax, label in zip(axes, "abc"):
        add_panel_label(ax, label)
        despine(ax)
    fig.subplots_adjust(wspace=0.35)
    save_supp(fig, "supplementary_figure_s1_population_diversity")


def plot_s2() -> None:
    func = REPO / "output/phase2_annotation/functional_annotation"
    impact = pd.read_csv(func / "v2_variant_counts_by_impact.tsv", sep="\t")
    effect = pd.read_csv(func / "v2_variant_counts_by_effect.tsv", sep="\t")
    hi = pd.read_csv(REPO / "output/phase2_annotation/snpeff_annotation/high_impact_variants.csv")
    fig, axes = plt.subplots(1, 3, figsize=(7.09, 2.7), gridspec_kw={"width_ratios": [0.95, 1.2, 1.2]})
    order = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
    impact = impact.set_index("Impact").reindex(order).dropna(how="all").reset_index()
    axes[0].bar(np.arange(len(impact)), impact["Count"], color=COLORS["blue"], width=0.7)
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Annotated consequences (log)")
    axes[0].set_xticks(np.arange(len(impact)))
    axes[0].set_xticklabels(impact["Impact"], rotation=25, ha="right", fontsize=3.6)

    def _short_effect(label: str) -> str:
        s = str(label).replace("_", " ").replace(" & ", " + ").replace("&", " + ")
        s = s.replace(" variant", "")
        if len(s) <= 13:
            return s
        parts = s.replace(" + ", " +|").split()
        # Prefer break after '+'
        if "|" in s.replace(" + ", " +|"):
            left, right = s.split(" + ", 1)
            return f"{left} +\n{right}"
        mid = (len(parts) + 1) // 2
        return " ".join(parts[:mid]) + "\n" + " ".join(parts[mid:])

    def _barh_labeled(ax, labels, values, color: str, xlabel: str) -> None:
        y = np.arange(len(values))
        vals = np.asarray(values, dtype=float)
        ax.barh(y, vals, color=color, height=0.7)
        ax.set_yticks(y)
        ax.set_yticklabels([])
        xmax = float(vals.max()) if len(vals) else 1.0
        ax.set_xlim(0, xmax * 1.75)
        for yi, lab, val in zip(y, labels, vals):
            ax.text(
                val + 0.03 * xmax,
                yi,
                lab,
                ha="left",
                va="center",
                fontsize=3.2,
                color=COLORS["black"],
                clip_on=True,
            )
        ax.set_xlabel(xlabel)

    top_eff = effect.sort_values("Count", ascending=False).head(8).iloc[::-1]
    _barh_labeled(
        axes[1],
        [_short_effect(e) for e in top_eff["Effect"]],
        (top_eff["Count"] / 1e6).to_numpy(),
        COLORS["orange"],
        r"Count ($\times 10^6$)",
    )
    hi_counts = hi["effect"].value_counts().head(6).iloc[::-1]
    _barh_labeled(
        axes[2],
        [_short_effect(e) for e in hi_counts.index],
        hi_counts.to_numpy(dtype=float),
        COLORS["teal"],
        "High-impact variants",
    )
    for ax, label in zip(axes, "abc"):
        add_panel_label(ax, label)
        despine(ax)
    fig.subplots_adjust(left=0.09, right=0.99, wspace=0.38, bottom=0.18, top=0.92)
    save_supp(fig, "supplementary_figure_s2_variant_annotation")


def plot_s3() -> None:
    roh = pd.read_csv(ROH_SUMMARY).sort_values("Num_ROH")
    fig, axes = plt.subplots(1, 2, figsize=(7.09, 2.5))
    axes[0].bar(np.arange(len(roh)), roh["Num_ROH"], color=COLORS["blue"], width=1.0)
    axes[0].set_xlabel("Individuals (ordered)")
    axes[0].set_ylabel("ROH count")
    axes[1].scatter(roh["Num_ROH"], roh["F_ROH"], s=14, color=COLORS["orange"], alpha=0.8)
    axes[1].set_xlabel("ROH count")
    axes[1].set_ylabel(r"$F_{ROH}$")
    for ax, label in zip(axes, "ab"):
        add_panel_label(ax, label)
        despine(ax)
    fig.subplots_adjust(wspace=0.35)
    save_supp(fig, "supplementary_figure_s3_roh_diagnostics")


def load_unfolded_sfs() -> pd.DataFrame:
    """Parse fastsimcoal2 .obs SFS (skip monomorphic bin 0)."""
    obs = SFS_DIR / "SNJ_DAFpop0.obs"
    lines = [ln for ln in obs.read_text().splitlines() if ln.strip()]
    header = lines[1].lstrip("\t").split("\t")
    counts = [float(x) for x in lines[2].lstrip("\t").split("\t")]
    frame = pd.DataFrame({"derived_count": np.arange(len(counts)), "snps": counts})
    return frame[frame["derived_count"] > 0].copy()


def plot_s4() -> None:
    """SFS, absolute AIC and bootstrap NCUR. ΔAIC omitted (shown in main Fig. 6a)."""
    comp = pd.read_csv(MODEL_COMPARISON).sort_values("AIC")
    fig, axes = plt.subplots(1, 3, figsize=(7.09, 2.6))
    sfs = load_unfolded_sfs()
    axes[0].bar(sfs["derived_count"], sfs["snps"], color=COLORS["blue"], width=1.0)
    axes[0].set_xlabel("Derived allele count")
    axes[0].set_ylabel("SNPs")
    axes[0].set_yscale("log")
    labels = [m.replace("_", "\n") for m in comp["Model"]]
    axes[1].barh(np.arange(len(comp)), comp["AIC"], color=COLORS["orange"])
    axes[1].set_yticks(np.arange(len(comp)))
    axes[1].set_yticklabels(labels, fontsize=4)
    axes[1].set_xlabel("AIC")
    axes[1].invert_yaxis()
    boot = pd.read_csv(
        REPO / "output/phase3b_fastsimcoal2/bootstrap/bottleneck_recent_contraction/bootstrap_results.csv"
    )
    ncur = [float(str(x).split()[0]) for x in boot["parameters"]]
    axes[2].hist(ncur, bins=20, color=COLORS["purple"], edgecolor="white", linewidth=0.2)
    axes[2].set_xlabel(r"Present $N_e$ (bootstrap)")
    axes[2].set_ylabel("Replicates")
    for ax, label in zip(axes, "abc"):
        add_panel_label(ax, label)
        despine(ax)
    fig.subplots_adjust(wspace=0.45)
    save_supp(fig, "supplementary_figure_s4_sfs_demography")


def plot_s5() -> None:
    esm = pd.read_csv(ESM2_PREDICTIONS, usecols=["esm2_score", "esm2_percentile"]).dropna()
    fig, axes = plt.subplots(1, 2, figsize=(7.09, 2.4))
    axes[0].hist(esm["esm2_percentile"], bins=30, color=COLORS["blue"], edgecolor="white", linewidth=0.2)
    axes[0].set_xlabel("ESM-2 percentile")
    axes[0].set_ylabel("Variants")
    axes[1].hist(esm["esm2_score"], bins=40, color=COLORS["orange"], edgecolor="white", linewidth=0.2)
    axes[1].axvline(0, color=COLORS["grey"], linestyle="--", linewidth=0.6)
    axes[1].set_xlabel("ESM-2 LLR")
    axes[1].set_ylabel("Variants")
    for ax, label in zip(axes, "ab"):
        add_panel_label(ax, label)
        despine(ax)
    fig.subplots_adjust(wspace=0.35)
    save_supp(fig, "supplementary_figure_s5_esm2_diagnostics")


def plot_s6() -> None:
    load = pd.read_csv(INDIVIDUAL_LOAD)
    pheno = pd.read_csv(MERGED_PHENO)
    merged = pheno.merge(load, on="IID", how="left", suffixes=("", "_load"))
    fig, axes = plt.subplots(2, 2, figsize=(7.09, 4.5))
    for ax, col, label in zip(
        axes.flat,
        ["Total_Genetic_Load", "Potential_Load", "Realized_Load", "Hom_Realized_Load"],
        "abcd",
    ):
        if col in merged.columns:
            ax.hist(merged[col].dropna(), bins=18, color=COLORS["blue"], edgecolor="white", linewidth=0.2)
            ax.set_xlabel(col.replace("_", " "))
        add_panel_label(ax, label)
        despine(ax)
    fig.subplots_adjust(hspace=0.45, wspace=0.35)
    save_supp(fig, "supplementary_figure_s6_genetic_load")


def plot_s7() -> None:
    pheno = pd.read_csv(MERGED_PHENO)
    fig, axes = plt.subplots(2, 2, figsize=(7.09, 4.5))
    axes[0, 0].scatter(pheno["F_ROH"], pheno["Load_in_ROH"], s=12, color=COLORS["blue"], alpha=0.7)
    axes[0, 0].set_xlabel(r"$F_{ROH}$")
    axes[0, 0].set_ylabel("Load in ROH")
    axes[0, 1].scatter(pheno["F_ROH_SHORT"], pheno["Load_in_Short_ROH"], s=12, color=COLORS["orange"], alpha=0.7)
    axes[0, 1].set_xlabel(r"$F_{ROH}$ short")
    axes[1, 0].scatter(pheno["Hom_Realized_Load"], pheno["Het_Realized_Load"], s=12, color=COLORS["teal"], alpha=0.7)
    axes[1, 0].set_xlabel("Homozygous realized load")
    axes[1, 0].set_ylabel("Heterozygous realized load")
    axes[1, 1].scatter(pheno["Total_Genetic_Load"], pheno["Potential_Load"], s=12, color=COLORS["purple"], alpha=0.7)
    axes[1, 1].set_xlabel("Total genetic load")
    axes[1, 1].set_ylabel("Potential load")
    for ax, label in zip(axes.flat, "abcd"):
        add_panel_label(ax, label)
        despine(ax)
    fig.subplots_adjust(hspace=0.45, wspace=0.35)
    save_supp(fig, "supplementary_figure_s7_roh_load")


def plot_s8() -> None:
    pheno = pd.read_csv(MERGED_PHENO)
    corr = pd.read_csv(CORRELATION)
    cc = pd.read_csv(CASE_CONTROL)
    metrics = ["Total_Genetic_Load", "Potential_Load", "Hom_Realized_Load", "F_ROH", "Het_Realized_Load", "Realized_Load"]
    fig = plt.figure(figsize=(7.09, 5.5))
    gs = GridSpec(2, 2, figure=fig)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])
    case_vals = [pheno[pheno["Has_Disease"] == 1][m].mean() for m in metrics if m in pheno.columns]
    ctrl_vals = [pheno[pheno["Has_Disease"] == 0][m].mean() for m in metrics if m in pheno.columns]
    x = np.arange(len(case_vals))
    ax_a.bar(x - 0.2, ctrl_vals, width=0.4, label="Unaffected", color=COLORS["sky"])
    ax_a.bar(x + 0.2, case_vals, width=0.4, label="Cases", color=COLORS["vermillion"])
    ax_a.set_xticks(x)
    ax_a.set_xticklabels([m.replace("_", "\n") for m in metrics if m in pheno.columns], fontsize=4, rotation=0)
    ax_a.legend(frameon=False, fontsize=4)
    chs = corr[corr["Phenotype_Variable"] == "CHS"].sort_values("Spearman_rho")
    ax_b.barh(chs["Genomic_Variable"], chs["Spearman_rho"], color=COLORS["blue"])
    ax_b.set_xlabel("Spearman rho with CHS")
    y = np.arange(len(cc))[::-1]
    ax_c.barh(y, cc["Cohens_d"], color=[COLORS["blue"] if d > 0 else COLORS["orange"] for d in cc["Cohens_d"]])
    ax_c.set_yticks(y)
    ax_c.set_yticklabels(cc["Variable"], fontsize=4)
    ax_c.set_xlabel("Cohen's d")
    for ax, label in zip([ax_a, ax_b, ax_c], "abc"):
        add_panel_label(ax, label)
        despine(ax)
    fig.subplots_adjust(hspace=0.45, wspace=0.4)
    save_supp(fig, "supplementary_figure_s8_phenotype_genotype")


def read_gwas(path: Path, kind: Optional[str] = None) -> pd.DataFrame:
    """Load CHR/BP/P from PLINK, GCTA-MLMA, or GEMMA association outputs."""
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n")
    raw_cols = header.split()
    upper_map = {c.upper().lstrip("#"): c for c in raw_cols}
    keys = set(upper_map)

    if kind is None:
        if path.suffix == ".mlma" or {"CHR", "BP", "P"} <= keys or (
            "CHR" in keys and "BP" in keys and "P" in keys
        ):
            # GCTA .mlma uses Chr/bp/p
            if "BP" in keys and "P" in keys:
                kind = "gcta_mlma"
            elif "PS" in keys and ("P_WALD" in keys or "P_LRT" in keys):
                kind = "gemma"
            else:
                kind = "plink"
        elif "PS" in keys and ("P_WALD" in keys or "P_LRT" in keys):
            kind = "gemma"
        else:
            kind = "plink"

    if kind == "gcta_mlma":
        chrom = upper_map.get("CHR")
        pos = upper_map.get("BP")
        pcol = upper_map.get("P")
    elif kind == "gemma":
        chrom = upper_map.get("CHR")
        pos = upper_map.get("PS")
        pcol = upper_map.get("P_WALD") or upper_map.get("P_LRT")
    else:
        chrom = upper_map.get("CHROM") or upper_map.get("CHR")
        pos = upper_map.get("POS") or upper_map.get("BP")
        pcol = upper_map.get("P") or upper_map.get("PVAL")

    usecols = [c for c in (chrom, pos, pcol) if c is not None]
    df = pd.read_csv(path, sep=r"\s+", usecols=usecols, low_memory=False)
    df.columns = [c.upper().lstrip("#") for c in df.columns]
    rename = {}
    if "CHROM" in df.columns:
        rename["CHROM"] = "CHR"
    if "PS" in df.columns:
        rename["PS"] = "BP"
    if "POS" in df.columns:
        rename["POS"] = "BP"
    if "P_WALD" in df.columns:
        rename["P_WALD"] = "P"
    elif "P_LRT" in df.columns:
        rename["P_LRT"] = "P"
    df = df.rename(columns=rename)
    df = df.dropna(subset=["P", "CHR"]).copy()
    df = df[np.isfinite(df["P"]) & (df["P"] > 0)]
    df["CHR"] = df["CHR"].astype(str).str.replace("NC_044549.1", "1", regex=False)
    if "BP" not in df.columns and "POS" in df.columns:
        df = df.rename(columns={"POS": "BP"})
    return df[["CHR", "BP", "P"]] if "BP" in df.columns else df


def manhattan_qq(
    ax_m,
    ax_q,
    gwas_path: Path,
    title: str,
    kind: Optional[str] = None,
) -> None:
    if not gwas_path.exists():
        ax_m.text(0.5, 0.5, f"Missing\n{gwas_path.name}", ha="center", va="center")
        ax_q.axis("off")
        return
    df = read_gwas(gwas_path, kind=kind)
    df["logp"] = -np.log10(df["P"].clip(lower=1e-300))

    def _chrom_key(x: str):
        try:
            return (0, int(float(x)))
        except ValueError:
            return (1, str(x))

    chrom_order = sorted(df["CHR"].unique(), key=_chrom_key)
    cursor = 0.0
    xticks, xticklabels = [], []
    rng = np.random.default_rng(1)
    for chrom in chrom_order:
        sub = df[df["CHR"] == chrom]
        xpos = sub["BP"].to_numpy(dtype=float) if "BP" in sub.columns else np.arange(len(sub), dtype=float)
        logp = sub["logp"].to_numpy(dtype=float)
        keep = logp >= 3.0
        n_bg = int((~keep).sum())
        if n_bg > 40_000:
            bg_idx = np.flatnonzero(~keep)
            chosen = rng.choice(bg_idx, size=40_000, replace=False)
            mask = keep.copy()
            mask[chosen] = True
        else:
            mask = np.ones(len(sub), dtype=bool)
        span = float(np.nanmax(xpos)) if len(xpos) else 0.0
        xticks.append(cursor + span / 2.0)
        xticklabels.append(str(int(float(chrom))) if str(chrom).replace(".", "", 1).isdigit() else str(chrom))
        ax_m.scatter(
            xpos[mask] + cursor,
            logp[mask],
            s=3,
            alpha=0.45,
            rasterized=True,
        )
        cursor += span + max(span * 0.02, 1.0)
    thr = -np.log10(GWAS_BONFERRONI_P)
    ax_m.axhline(thr, color=COLORS["vermillion"], linestyle="--", linewidth=0.7, zorder=3)
    if len(chrom_order) > 1:
        step = max(1, len(chrom_order) // 12)
        ax_m.set_xticks(xticks[::step])
        ax_m.set_xticklabels(xticklabels[::step], fontsize=3.5)
        ax_m.set_xlabel("Chromosome")
    else:
        ax_m.set_xlabel("Position")
    ax_m.set_ylabel(r"$-\log_{10}(P)$")
    ax_m.set_title(title, fontsize=5, loc="left")
    obs = np.sort(df["P"].clip(lower=1e-300).to_numpy())
    exp = np.linspace(1 / len(obs), 1, len(obs))
    if len(obs) > 200_000:
        keep = np.unique(
            np.concatenate(
                [
                    np.linspace(0, len(obs) - 1, 100_000, dtype=int),
                    np.arange(0, min(5000, len(obs))),
                ]
            )
        )
        obs_q = obs[keep]
        exp_q = exp[keep]
    else:
        obs_q, exp_q = obs, exp
    ax_q.scatter(-np.log10(exp_q), -np.log10(obs_q), s=3, color=COLORS["blue"], alpha=0.45, rasterized=True)
    ax_q.plot([0, max(-np.log10(exp))], [0, max(-np.log10(exp))], color=COLORS["grey"], linewidth=0.6)
    ax_q.set_xlabel("Expected")
    ax_q.set_ylabel("Observed")


def plot_s15() -> None:
    # Trait-specific GWAS: GCTA-MLMA (CHS) + GEMMA LMM (binary traits)
    traits = [
        (GWAS_DIR / "mlma_chs.mlma", "CHS", "gcta_mlma"),
        (
            GWAS_DIR / "output" / "gemma_has_eye_disease.assoc.txt",
            "Eye disease",
            "gemma",
        ),
        (
            GWAS_DIR / "output" / "gemma_has_finger_joint_abnormality.assoc.txt",
            "Finger joint",
            "gemma",
        ),
    ]
    fig, axes = plt.subplots(3, 2, figsize=(7.09, 6.5))
    for row, (path, title, kind) in enumerate(traits):
        manhattan_qq(axes[row, 0], axes[row, 1], path, title, kind=kind)
        add_panel_label(axes[row, 0], "abcdef"[row * 2])
        add_panel_label(axes[row, 1], "abcdef"[row * 2 + 1])
        for ax in axes[row]:
            despine(ax)
    fig.subplots_adjust(hspace=0.55, wspace=0.35)
    save_supp(fig, "supplementary_figure_s15_trait_gwas")


def _enrichment_dir() -> Path:
    folder = REPO / "works" / "enrichment"
    required = folder / "gprofiler_top10pct_missense_plot_terms.csv"
    if not required.is_file():
        raise SystemExit(
            "Supp. Fig. S16 reads g:Profiler tables from the analysis root, "
            f"not from this code repository:\n  {required}\n"
            "Place gene lists and FDR tables under $PLM_BASE_DIR/works/enrichment/."
        )
    return folder


SOURCE_LABEL = {
    "GO:MF": "GO MF",
    "GO:BP": "GO BP",
    "GO:CC": "GO CC",
    "KEGG": "KEGG",
    "HP": "HPO",
}
SOURCE_COLOR = {
    "GO:MF": "#4C72B0",
    "GO:BP": "#55A868",
    "GO:CC": "#C44E52",
    "KEGG": "#8172B3",
    "HP": "#CCB974",
}
S16_BAR_H = 0.72


def _plot_s16_panel_a(ax: plt.Axes, top_plot: pd.DataFrame) -> None:
    d = top_plot.sort_values("negative_log10_of_adjusted_p_value", ascending=True).reset_index(
        drop=True
    )
    y = np.arange(len(d))
    ax.barh(
        y,
        d["negative_log10_of_adjusted_p_value"],
        color=[SOURCE_COLOR.get(s, "#777") for s in d["source"]],
        height=S16_BAR_H,
    )
    labels = []
    for s, name in zip(d["source"], d["term_name"]):
        lab = f"{SOURCE_LABEL.get(s, s)}  {name}"
        labels.append(lab if len(lab) <= 56 else lab[:53] + "…")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    xmax = float(d["negative_log10_of_adjusted_p_value"].max())
    for yi, n, x in zip(y, d["intersection_size"], d["negative_log10_of_adjusted_p_value"]):
        ax.text(x + 0.04 * xmax, yi, f"n={int(n)}", va="center", fontsize=7, color="#333")
    ax.set_xlim(0, xmax * 1.28)
    ax.set_ylim(-0.5, len(d) - 0.5)
    ax.axvline(-np.log10(0.05), color="#888", ls="--", lw=0.8, zorder=0)
    ax.set_title("a  Top-decile missense (P ≥ 0.612)", loc="left", fontsize=10, fontweight="bold", pad=6)
    ax.set_xlabel(r"−log$_{10}$(FDR)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _plot_s16_panel_b(ax: plt.Axes, lof_plot: pd.DataFrame, n_slots: int) -> None:
    """Same y-span / bar height as panel a so a single hit is not stretched."""
    hit = lof_plot.sort_values("adjusted_p_value").iloc[0]
    hit_val = float(hit["negative_log10_of_adjusted_p_value"])
    y_hit, y_hpo = n_slots - 1, n_slots - 2

    ax.axvline(-np.log10(0.05), color="#888", ls="--", lw=0.8, zorder=0)
    ax.barh(y_hit, hit_val, color=SOURCE_COLOR["GO:MF"], height=S16_BAR_H)
    ax.text(
        hit_val + 0.05,
        y_hit,
        f"n={int(hit['intersection_size'])}",
        va="center",
        fontsize=7,
        color="#333",
    )
    ax.text(0.08, y_hpo, "—", va="center", ha="left", fontsize=11, color="#666", fontweight="bold")
    ax.set_yticks([y_hpo, y_hit])
    ax.set_yticklabels(
        ["HPO  no FDR-significant terms", "GO MF  inhibitory MHC I receptor activity"],
        fontsize=8,
    )
    ax.set_ylim(-0.5, n_slots - 0.5)
    ax.set_xlim(0, max(2.2, hit_val * 1.55))
    ax.set_title("b  LoF variants", loc="left", fontsize=10, fontweight="bold", pad=6)
    ax.set_xlabel(r"−log$_{10}$(FDR)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_s16() -> None:
    enrichment = _enrichment_dir()
    top_plot = pd.read_csv(enrichment / "gprofiler_top10pct_missense_plot_terms.csv")
    lof_plot = pd.read_csv(enrichment / "gprofiler_lof_plot_terms.csv")
    if lof_plot.empty:
        lof_plot = pd.read_csv(enrichment / "gprofiler_lof_informative_terms.csv")

    fig = plt.figure(figsize=(10.8, 5.5))
    gs = fig.add_gridspec(
        1,
        2,
        width_ratios=[1.7, 0.95],
        left=0.30,
        right=0.98,
        top=0.86,
        bottom=0.11,
        wspace=0.58,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    _plot_s16_panel_a(ax_a, top_plot)
    _plot_s16_panel_b(ax_b, lof_plot, n_slots=len(top_plot))

    handles = [
        Patch(color=SOURCE_COLOR[k], label=SOURCE_LABEL[k])
        for k in ["GO:MF", "GO:BP", "GO:CC", "KEGG", "HP"]
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=5,
        frameon=False,
        bbox_to_anchor=(0.55, 0.92),
        columnspacing=1.2,
        handletextpad=0.5,
    )
    save_supp(fig, "supplementary_figure_s16_burden_gene_enrichment")


def plot_s10() -> None:
    summary = pd.read_csv(CROSS / "llr_distribution_summary.csv")
    fig, ax = plt.subplots(figsize=(3.5, 2.5))
    cohorts = summary["cohort"].tolist()
    means = summary["llr_mean"].tolist()
    colors = [COLORS["blue"], COLORS["sky"], COLORS["orange"]]
    ax.bar(cohorts, means, color=colors[: len(cohorts)])
    ax.set_ylabel("Mean ESM-2 LLR")
    ax.set_xticklabels(cohorts, rotation=20, ha="right", fontsize=4)
    despine(ax)
    save_supp(fig, "supplementary_figure_s10_cross_species")


def plot_s11() -> None:
    frame = scan_variants()
    trend = probability_daf_trend(frame)
    depletion = homozygote_depletion_by_bin(frame)
    polarized = polarized_subset(frame)
    missense = polarized[(polarized["kind"] == "missense") & polarized["pathogenicity_prob"].notna()].copy()
    missense["risk_bin"] = missense["pathogenicity_prob"].map(assign_risk_bin)
    fig, axes = plt.subplots(1, 2, figsize=(7.09, 2.8))
    for label, sub in missense.groupby("risk_bin"):
        if label in RISK_LABELS:
            axes[0].hist(sub["daf"], bins=20, alpha=0.5, label=label)
    axes[0].set_xlabel("Polarized DAF")
    axes[0].set_ylabel("Variants")
    axes[0].legend(frameon=False, fontsize=4)
    axes[1].bar(np.arange(len(depletion)), depletion["oe_ratio"], color=COLORS["teal"])
    axes[1].axhline(1, color=COLORS["grey"], linestyle="--", linewidth=0.6)
    axes[1].set_xticks(np.arange(len(depletion)))
    axes[1].set_xticklabels(depletion["risk_bin"])
    axes[1].set_ylabel("Obs/exp derived homozygotes")
    axes[0].text(0.03, 0.97, f"rho={trend['spearman_rho']:.3f}, P={trend['p_value']:.2e}", transform=axes[0].transAxes, va="top", fontsize=5)
    for ax, label in zip(axes, "ab"):
        add_panel_label(ax, label)
        despine(ax)
    fig.subplots_adjust(wspace=0.35)
    save_supp(fig, "supplementary_figure_s11_polarized_consistency")


def plot_s12() -> None:
    conc = pd.read_csv(LOAD_SENS / "scenario_concordance.csv")
    froh = pd.read_csv(LOAD_SENS / "froh_associations.csv")
    ind = pd.read_csv(LOAD_SENS / "individual_load_sensitivity.csv")
    main = "main_all_scored_continuous_h0.25_lof0.95"
    rank = (
        conc[(conc["reference"] == main) & (~conc["scenario"].str.startswith("negative_"))]
        .drop_duplicates("scenario")
        .sort_values("spearman_rho")
    )
    keep = [c for c in ind["scenario"].unique() if not str(c).startswith("negative_")]
    wide = ind[ind["scenario"].isin(keep)].pivot_table(
        index="IID", columns="scenario", values="z_load", aggfunc="first"
    )
    # Prefer a compact set of primary sensitivity scenarios for the heatmap
    preferred = [
        main,
        "dominance_h0",
        "dominance_h0.1",
        "dominance_h0.5",
        "lof_weight_0.5",
        "lof_weight_1",
        "top_10pct_pathogenicity",
        "lof_only",
        "missense_only_all_scored",
    ]
    cols = [c for c in preferred if c in wide.columns]
    corr = wide[cols].corr(method="pearson")
    syn = froh[froh["scenario"].str.startswith("negative_synonymous")]
    main_rho = froh.loc[froh["scenario"] == main, "spearman_rho"].iloc[0]

    fig, axes = plt.subplots(1, 3, figsize=(7.09, 2.6))
    short = (
        rank["scenario"]
        .str.replace(main, "main", regex=False)
        .str.replace("dominance_", "", regex=False)
        .str.replace("lof_weight_", "LoF w=", regex=False)
        .str.replace("pathogenicity_p_gt_", "P>", regex=False)
        .str.replace("top_", "top ", regex=False)
        .str.replace("pct_pathogenicity", "% P", regex=False)
        .str.replace("missense_only_all_scored", "missense only", regex=False)
        .str.replace("count_carrier_dominance_h0.25", "count carrier", regex=False)
        .str.replace("count_alt_alleles", "count ALT alleles", regex=False)
        .str.replace("lof_only", "LoF only", regex=False)
        .str.replace("_", " ")
    )
    axes[0].barh(short, rank["spearman_rho"], color=COLORS["blue"])
    axes[0].set_xlabel(r"Spearman $\rho$ vs main load")
    axes[0].set_xlim(0, 1.05)
    im = axes[1].imshow(corr.to_numpy(), cmap="viridis", vmin=0.5, vmax=1.0, aspect="auto")
    axes[1].set_xticks(np.arange(len(cols)))
    axes[1].set_yticks(np.arange(len(cols)))
    tick = [c.replace(main, "main").replace("_", "\n") for c in cols]
    axes[1].set_xticklabels(tick, fontsize=3.5, rotation=90)
    axes[1].set_yticklabels(tick, fontsize=3.5)
    fig.colorbar(im, ax=axes[1], fraction=0.046, pad=0.04)
    axes[2].hist(syn["spearman_rho"], bins=20, color=COLORS["teal"], edgecolor="white", linewidth=0.2)
    axes[2].axvline(main_rho, color=COLORS["vermillion"], linestyle="--", linewidth=0.8, label="Main load")
    axes[2].set_xlabel(r"Synonymous-control Spearman $\rho$ with $F_{ROH}$")
    axes[2].set_ylabel("Replicates")
    axes[2].legend(frameon=False, fontsize=4)
    for ax, label in zip(axes, "abc"):
        add_panel_label(ax, label)
        despine(ax)
    fig.subplots_adjust(wspace=0.45)
    save_supp(fig, "supplementary_figure_s12_load_sensitivity")


def plot_s13() -> None:
    fixed = pd.read_csv(FITNESS_FIXED)
    perm = pd.read_csv(FITNESS_PERM)
    cv = pd.read_csv(FITNESS_CV)
    fig, axes = plt.subplots(1, 3, figsize=(7.09, 2.5))
    both = fixed[fixed["model"] == "covariates_both"]
    axes[0].barh(both["term"], both["beta"], color=COLORS["blue"])
    axes[0].set_xlabel("Beta")
    axes[1].barh(perm["test"], perm["p_freedman_lane"], color=COLORS["orange"])
    axes[1].set_xlabel("Permutation P")
    axes[2].bar(cv["model"], cv["RMSE"], color=COLORS["teal"])
    axes[2].set_ylabel("Grouped-CV RMSE")
    for ax, label in zip(axes, "abc"):
        add_panel_label(ax, label)
        despine(ax)
    fig.subplots_adjust(wspace=0.45)
    save_supp(fig, "supplementary_figure_s13_morbidity_validation")


def _ne_schedule_from_metadata() -> tuple[np.ndarray, np.ndarray]:
    import json

    meta = json.loads((REPO / "output/method_validation/dfe_simulation/simulation_metadata.json").read_text())
    dem = meta["demography_mle"]
    gen_time = 10.0
    # Piecewise history used by the transient-drift scenario (generations ago → Ne)
    edges = [
        (0.0, dem["trecent"], dem["ncur"]),
        (dem["trecent"], dem["trecovery_old"], dem["nrecover"]),
        (dem["trecovery_old"], dem["tbot_old"], dem["nbot"]),
        (dem["tbot_old"], dem["tbot_old"] + dem.get("burn_in", 2000), dem["nanc"]),
    ]
    times, nes = [], []
    for start, end, ne in edges:
        times.extend([start * gen_time / 1000.0, end * gen_time / 1000.0])
        nes.extend([ne, ne])
    return np.asarray(times), np.asarray(nes)


def plot_s14() -> None:
    reps = pd.read_csv(DFE_REPLICATES)
    bins = pd.read_csv(REPO / "output/method_validation/dfe_simulation/observed_vs_simulated_bins.csv")
    fig, axes = plt.subplots(2, 2, figsize=(7.09, 4.5))
    for scenario, color in zip(
        ["transient_drift", "constant_current_ne", "neutral_drift"],
        [COLORS["blue"], COLORS["orange"], COLORS["teal"]],
    ):
        sub = reps[reps["scenario"] == scenario]
        axes[0, 0].hist(sub["frac_fixed"], bins=20, alpha=0.45, label=scenario.replace("_", " "), color=color)
    axes[0, 0].set_xlabel("Fixed fraction")
    axes[0, 0].set_ylabel("Replicates")
    axes[0, 0].legend(frameon=False, fontsize=4)
    td = bins[bins["scenario"] == "transient_drift"]
    order = ["rare", "low", "intermediate", "common", "fixed"]
    td = td.set_index("bin").reindex(order).reset_index()
    x = np.arange(len(order))
    axes[0, 1].bar(x - 0.18, td["observed_fraction"], width=0.36, color=COLORS["blue"], label="Observed")
    axes[0, 1].bar(x + 0.18, td["simulated_fraction"], width=0.36, color=COLORS["orange"], label="Simulated")
    axes[0, 1].set_xticks(x)
    axes[0, 1].set_xticklabels(order, fontsize=4, rotation=20)
    axes[0, 1].set_ylabel("Fraction")
    axes[0, 1].legend(frameon=False, fontsize=4)
    t_kya, ne = _ne_schedule_from_metadata()
    axes[1, 0].step(t_kya, ne, where="post", color=COLORS["blue"], linewidth=0.9)
    axes[1, 0].set_xscale("log")
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_xlabel("Time before present (kya, log)")
    axes[1, 0].set_ylabel(r"$N_e$ (log)")
    axes[1, 0].invert_xaxis()
    obs = reps[reps["scenario"] == "transient_drift"]
    axes[1, 1].scatter(obs["frac_intermediate"], obs["frac_fixed"], s=8, alpha=0.55, color=COLORS["blue"])
    axes[1, 1].set_xlabel("Intermediate fraction")
    axes[1, 1].set_ylabel("Fixed fraction")
    for ax, label in zip(axes.flat, "abcd"):
        add_panel_label(ax, label)
        despine(ax)
    fig.subplots_adjust(hspace=0.45, wspace=0.35)
    save_supp(fig, "supplementary_figure_s14_dfe_simulation")


def plot_s9() -> None:
    metrics = pd.read_csv(CLINVAR / "clinvar_validation_metrics.csv")
    gkf = pd.read_csv(CLINVAR / "clinvar_groupkfold_summary.csv")
    row = metrics[metrics["model"] == "fixed_sigmoid_protein_held_out"].iloc[0]
    fig, axes = plt.subplots(1, 2, figsize=(7.09, 2.5))
    bars = ["ROC-AUC", "PR-AUC", "Balanced accuracy"]
    vals = [row["roc_auc"], row["pr_auc"], row["balanced_accuracy_at_0.5"]]
    axes[0].bar(bars, vals, color=[COLORS["blue"], COLORS["orange"], COLORS["teal"]])
    axes[0].set_ylim(0, 1)
    axes[0].set_ylabel("Metric value")
    axes[0].set_title("Held-out ClinVar (n = 13,147)", fontsize=5, loc="left")
    if "roc_auc_mean" in gkf.columns:
        axes[1].bar(
            ["ROC-AUC"],
            [gkf["roc_auc_mean"].iloc[0]],
            yerr=[gkf["roc_auc_std"].iloc[0]],
            color=COLORS["blue"],
            capsize=3,
        )
        axes[1].set_ylim(0, 1)
        axes[1].set_title("Five-fold protein-group CV", fontsize=5, loc="left")
    for ax, label in zip(axes, "ab"):
        add_panel_label(ax, label)
        despine(ax)
    fig.subplots_adjust(wspace=0.35)
    save_supp(fig, "supplementary_figure_s9_clinvar_validation")


def main() -> None:
    apply_style()
    plot_s1()
    plot_s2()
    plot_s3()
    plot_s4()
    plot_s5()
    plot_s6()
    plot_s7()
    plot_s8()
    plot_s9()
    plot_s10()
    plot_s11()
    plot_s12()
    plot_s13()
    plot_s14()
    plot_s15()
    plot_s16()


if __name__ == "__main__":
    main()
