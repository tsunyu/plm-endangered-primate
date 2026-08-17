#!/usr/bin/env python3
"""Figure 4 — morbidity tracks genetic load rather than recent inbreeding."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lib.bootstrap import bootstrap_cohens_d_ci, bootstrap_spearman_ci
from lib.export import save_figure
from lib.paths import CORRELATION, FITNESS_FIXED, FITNESS_PERM, MERGED_PHENO, RESULTS
from lib.style import COLORS, add_panel_label, despine

# Shared metric order keeps panels a and c visually aligned.
FOREST_METRICS = [
    ("Total_Deleterious", "Total deleterious"),
    ("Total_Genetic_Load", "Total genetic load"),
    ("Realized_Load", "Total realized load"),
    ("Het_Realized_Load", "Het realized load"),
    ("Potential_Load", "Potential load"),
    ("Hom_Realized_Load", "Hom realized load"),
    ("F_ROH", r"$F_{ROH}$"),
]

JOINT_MODEL_TERMS = [
    ("Total_Genetic_Load_z", "Total genetic load"),
    ("F_ROH_z", r"$F_{ROH}$"),
]

PERMUTATION_TESTS = [
    ("load_marginal", "Load (marginal)", COLORS["blue"]),
    ("load_conditional_on_F_ROH", r"Load | $F_{ROH}$", COLORS["blue"]),
    ("F_ROH_marginal", r"$F_{ROH}$ (marginal)", COLORS["orange"]),
    ("F_ROH_conditional_on_load", r"$F_{ROH}$ | load", COLORS["orange"]),
]

FOREST_KW = dict(
    markersize=4.5,
    capsize=2.5,
    elinewidth=0.7,
    markeredgewidth=0,
    zorder=3,
)


def load_chs_correlations() -> pd.DataFrame:
    corr = pd.read_csv(CORRELATION)
    return corr[corr["Phenotype_Variable"] == "CHS"].set_index("Genomic_Variable")


def _forest_y_positions(n: int) -> np.ndarray:
    """First metric at top (highest y), last metric at bottom."""
    return np.arange(n - 1, -1, -1)


def _style_forest_axis(ax, y_positions: np.ndarray, show_ylabels: bool) -> None:
    ax.axvline(0, color=COLORS["grey"], linewidth=0.6, zorder=1)
    ax.set_yticks(y_positions)
    if show_ylabels:
        ax.set_yticklabels([label for _, label in FOREST_METRICS], fontsize=5)
    else:
        ax.set_yticklabels([])
        ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.6, len(FOREST_METRICS) - 0.4)
    despine(ax)


def _plot_forest_points(ax, values: list[tuple[float, float, float]], y_positions: np.ndarray) -> None:
    for y, (est, lo, hi) in zip(y_positions, values):
        color = COLORS["blue"] if est > 0 else COLORS["orange"]
        ax.errorbar(
            est,
            y,
            xerr=[[est - lo], [hi - est]],
            fmt="o",
            color=color,
            markerfacecolor=color,
            **FOREST_KW,
        )


def plot_chs_forest(ax, pheno: pd.DataFrame, reported: pd.DataFrame, *, show_ylabels: bool = True) -> pd.DataFrame:
    rows = []
    values = []
    y_positions = _forest_y_positions(len(FOREST_METRICS))
    for col, label in FOREST_METRICS:
        rho, lo, hi, _ = bootstrap_spearman_ci(pheno, col, "CHS")
        fdr = reported.loc[col, "Spearman_p_adj"] if col in reported.index else np.nan
        sig = bool(fdr < 0.05) if pd.notna(fdr) else False
        rows.append(
            {
                "metric": label,
                "column": col,
                "rho": rho,
                "ci_lower": lo,
                "ci_upper": hi,
                "p_value": reported.loc[col, "Spearman_p"] if col in reported.index else np.nan,
                "fdr": fdr,
                "fdr_significant": sig,
            }
        )
        values.append((rho, lo, hi))
    _plot_forest_points(ax, values, y_positions)
    _style_forest_axis(ax, y_positions, show_ylabels)
    ax.set_xlabel("Spearman ρ with CHS")
    ax.set_xlim(-0.35, 0.72)
    return pd.DataFrame(rows)


def plot_joint_model(ax) -> None:
    fixed = pd.read_csv(FITNESS_FIXED)
    subset = fixed[
        (fixed["model"] == "covariates_both") & fixed["term"].isin([term for term, _ in JOINT_MODEL_TERMS])
    ].set_index("term")
    y_positions = _forest_y_positions(len(JOINT_MODEL_TERMS))
    values = []
    labels = []
    for term, label in JOINT_MODEL_TERMS:
        row = subset.loc[term]
        values.append((row["beta"], row["beta"] - 1.96 * row["SE"], row["beta"] + 1.96 * row["SE"]))
        labels.append(label)
    _plot_forest_points(ax, values, y_positions)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, fontsize=5)
    ax.set_xlabel(r"Standardized $\beta$ (joint GRM model)")
    ax.set_xlim(-0.45, 1.05)
    ax.set_ylim(-0.45, 1.45)
    ax.axvline(0, color=COLORS["grey"], linewidth=0.6, zorder=1)
    despine(ax)


def plot_case_control_forest(ax, pheno: pd.DataFrame, *, show_ylabels: bool = True) -> pd.DataFrame:
    rows = []
    values = []
    y_positions = _forest_y_positions(len(FOREST_METRICS))
    for col, label in FOREST_METRICS:
        d, lo, hi, p_value, _ = bootstrap_cohens_d_ci(pheno, col)
        rows.append({"metric": label, "column": col, "cohens_d": d, "ci_lower": lo, "ci_upper": hi, "p_value": p_value})
        values.append((d, lo, hi))
    _plot_forest_points(ax, values, y_positions)
    _style_forest_axis(ax, y_positions, show_ylabels)
    ax.set_xlabel("Cohen's d (cases vs controls)")
    ax.set_xlim(-0.55, 1.35)
    return pd.DataFrame(rows)


def plot_freedman_lane(ax) -> pd.DataFrame:
    perm = pd.read_csv(FITNESS_PERM).set_index("test")
    rows = []
    y_positions = _forest_y_positions(len(PERMUTATION_TESTS))
    neg_log_p = []
    for idx, (test_id, label, color) in enumerate(PERMUTATION_TESTS):
        p_value = float(perm.loc[test_id, "p_freedman_lane"])
        value = -np.log10(max(p_value, 1e-4))
        neg_log_p.append(value)
        rows.append(
            {
                "test": test_id,
                "label": label,
                "p_freedman_lane": p_value,
                "neg_log10_p": value,
                "permutations": int(perm.loc[test_id, "permutations"]),
            }
        )
        ax.plot(value, y_positions[idx], "o", color=color, markersize=4.5, zorder=3)
    threshold = -np.log10(0.05)
    ax.axvline(threshold, color=COLORS["grey"], linewidth=0.6, linestyle="--", zorder=1)
    ax.set_yticks(y_positions)
    ax.set_yticklabels([label for _, label, _ in PERMUTATION_TESTS], fontsize=5)
    ax.set_xlabel(r"$-\log_{10}$ Freedman–Lane $P$")
    ax.set_xlim(-0.05, max(neg_log_p) * 1.25 + 0.15)
    ax.set_ylim(-0.45, len(PERMUTATION_TESTS) - 0.55)
    despine(ax)
    return pd.DataFrame(rows)


def main() -> None:
    pheno = pd.read_csv(MERGED_PHENO)
    reported = load_chs_correlations()
    fig = plt.figure(figsize=(7.09, 6.2))
    gs = fig.add_gridspec(2, 2, hspace=0.48, wspace=0.44)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    chs_df = plot_chs_forest(ax_a, pheno, reported, show_ylabels=True)
    plot_joint_model(ax_b)
    cc_df = plot_case_control_forest(ax_c, pheno, show_ylabels=True)
    perm_df = plot_freedman_lane(ax_d)
    fig.subplots_adjust(left=0.26, right=0.98, top=0.96, bottom=0.10)
    add_panel_label(ax_a, "a", x=-0.34)
    add_panel_label(ax_b, "b", x=-0.24)
    add_panel_label(ax_c, "c", x=-0.34)
    add_panel_label(ax_d, "d", x=-0.24)
    chs_df.to_csv(RESULTS / "figure4_chs_forest.csv", index=False)
    cc_df.to_csv(RESULTS / "figure4_case_control_forest.csv", index=False)
    perm_df.to_csv(RESULTS / "figure4_freedman_lane.csv", index=False)
    save_figure(fig, "figure_4_morbidity_load", tight=False)


if __name__ == "__main__":
    main()
