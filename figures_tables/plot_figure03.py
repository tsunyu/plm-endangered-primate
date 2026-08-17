#!/usr/bin/env python3
"""Figure 3 — ESM-2 landscape and polarized population-genetic depletion."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from lib.daf import homozygote_depletion_by_bin, probability_daf_trend, scan_variants, summarize_daf
from lib.export import save_figure
from lib.paths import ESM2_PREDICTIONS, RESULTS, SIGMOID_K, SIGMOID_X0
from lib.style import COLORS, add_panel_label, despine, figure


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-np.clip(SIGMOID_K * (x - SIGMOID_X0), -700, 700)))


def plot_llr_distribution(ax, inset_ax) -> None:
    esm = pd.read_csv(ESM2_PREDICTIONS, usecols=["esm2_score"])
    scores = esm["esm2_score"].dropna()
    probs = sigmoid(scores.values)
    ax.hist(scores, bins=50, color=COLORS["blue"], edgecolor="white", linewidth=0.2)
    neg_pct = (scores < 0).mean() * 100
    ax.set_xlabel("ESM-2 log-likelihood ratio")
    ax.set_ylabel("Missense variants")
    ax.text(
        0.97,
        0.95,
        f"n = {len(scores):,}\n{neg_pct:.1f}% negative LLR",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5,
    )
    inset_ax.hist(probs, bins=40, color=COLORS["orange"], edgecolor="white", linewidth=0.2)
    inset_ax.set_xlabel("Calibrated probability", fontsize=4.5)
    inset_ax.set_ylabel("Count", fontsize=4.5)
    inset_ax.tick_params(labelsize=4)
    despine(ax)
    despine(inset_ax)


def plot_probability_vs_daf(ax, frame: pd.DataFrame, trend: dict) -> None:
    missense = frame[
        (frame["kind"] == "missense")
        & (frame["polarization_status"] == "polarized")
        & frame["pathogenicity_prob"].notna()
    ]
    x = missense["pathogenicity_prob"].values
    y = missense["daf"].values
    hb = ax.hist2d(
        x,
        y,
        bins=(20, 20),
        range=[[0, 1], [0, 1]],
        cmap="Blues",
        cmin=1,
    )
    cbar = plt.colorbar(hb[3], ax=ax, fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=4)
    cbar.set_label("Variant count", fontsize=5)
    ax.set_xlabel("Calibrated pathogenicity probability")
    ax.set_ylabel("Polarized derived allele frequency")
    ax.text(
        0.03,
        0.97,
        f"n = {trend['n']:,}\nSpearman rho = {trend['spearman_rho']:.3f}\nP = {trend['p_value']:.2e}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5,
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.85),
    )
    despine(ax)


def plot_homozygote_depletion(ax, depletion: pd.DataFrame) -> None:
    x = np.arange(len(depletion))
    ax.bar(x, depletion["oe_ratio"], color=COLORS["teal"], width=0.65)
    ax.axhline(1.0, color=COLORS["grey"], linestyle="--", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(depletion["risk_bin"], rotation=0)
    ax.set_xlabel("Pathogenicity-probability bin")
    ax.set_ylabel("O/E derived homozygotes")
    despine(ax)


def main() -> None:
    frame = scan_variants()
    summarize_daf(frame)
    trend = probability_daf_trend(frame)
    depletion = homozygote_depletion_by_bin(frame)
    depletion.to_csv(RESULTS / "figure3_homozygote_depletion.csv", index=False)
    pd.DataFrame([trend]).to_csv(RESULTS / "figure3_daf_trend.csv", index=False)

    fig = figure(width_mm=180, height_mm=120)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.1, 1], hspace=0.45, wspace=0.42)
    ax_a = fig.add_subplot(gs[0, :])
    inset = inset_axes(ax_a, width="28%", height="55%", loc="upper right", borderpad=1)
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    plot_llr_distribution(ax_a, inset)
    plot_probability_vs_daf(ax_b, frame, trend)
    plot_homozygote_depletion(ax_c, depletion)
    add_panel_label(ax_a, "a")
    add_panel_label(ax_b, "b", x=-0.14)
    add_panel_label(ax_c, "c", x=-0.22)
    save_figure(fig, "figure_3_plm_population_genetics")


if __name__ == "__main__":
    main()
