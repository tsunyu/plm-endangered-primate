#!/usr/bin/env python3
"""Figure 2 — cohort, heterozygosity, F_ROH and ROH architecture."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from lib.bootstrap import bootstrap_cohens_d_ci
from lib.export import save_figure
from lib.paths import HET_SUMMARY, MERGED_PHENO, ROH_SUMMARY
from lib.style import COLORS, add_panel_label, despine, figure


def plot_cohort_composition(ax) -> None:
    pheno = pd.read_csv(MERGED_PHENO)
    pheno["Sex"] = pheno["Sex"].fillna("Unknown")
    pheno["disease_label"] = np.where(pheno["Has_Disease"] == 1, "Disease case", "Unaffected")
    counts = (
        pheno.groupby(["Sex", "disease_label"], observed=True)
        .size()
        .reset_index(name="count")
    )
    sex_order = ["Female", "Male"]
    disease_order = ["Unaffected", "Disease case"]
    bottom = {sex: 0 for sex in sex_order}
    colors = {"Unaffected": COLORS["sky"], "Disease case": COLORS["vermillion"]}
    for disease in disease_order:
        values = []
        for sex in sex_order:
            row = counts[(counts["Sex"] == sex) & (counts["disease_label"] == disease)]
            values.append(int(row["count"].iloc[0]) if not row.empty else 0)
        ax.bar(sex_order, values, bottom=[bottom[s] for s in sex_order], color=colors[disease], width=0.55)
        for sex, val in zip(sex_order, values):
            if val:
                ax.text(sex, bottom[sex] + val / 2, str(val), ha="center", va="center", fontsize=5)
            bottom[sex] += val
    ax.set_ylabel("Individuals (n = 68)")
    ax.set_xlabel("")
    ax.legend(
        handles=[Patch(facecolor=colors[d], label=d) for d in disease_order],
        frameon=False,
        loc="upper right",
    )
    despine(ax)


def plot_heterozygosity(ax) -> None:
    het = pd.read_csv(HET_SUMMARY)
    values = het["OBS_HET"]
    ax.hist(values, bins=18, color=COLORS["blue"], edgecolor="white", linewidth=0.3)
    mean_val = values.mean()
    ax.axvline(mean_val, color=COLORS["vermillion"], linestyle="--", linewidth=0.8)
    ax.text(
        0.97,
        0.95,
        f"mean = {mean_val:.3f}\ns.d. = {values.std():.3f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5,
    )
    ax.set_xlabel("Observed heterozygosity")
    ax.set_ylabel("Count")
    despine(ax)


def plot_froh(ax) -> None:
    roh = pd.read_csv(ROH_SUMMARY)
    values = roh["F_ROH"]
    ax.hist(values, bins=18, color=COLORS["teal"], edgecolor="white", linewidth=0.3)
    mean_val = values.mean()
    ax.axvline(mean_val, color=COLORS["vermillion"], linestyle="--", linewidth=0.8, label="cohort mean")
    ax.set_xlabel(r"$F_{ROH}$")
    ax.set_ylabel("Count")
    ax.legend(frameon=False, loc="upper right")
    despine(ax)


def plot_roh_classes(ax) -> None:
    roh = pd.read_csv(ROH_SUMMARY).sort_values("Num_ROH")
    short = roh["Short_ROH_Count"]
    medium = roh["Long_ROH_Count"]
    long = roh["VeryLong_ROH_Count"]
    x = np.arange(len(roh))
    ax.bar(x, short, color=COLORS["blue"], width=1.0, label="Short (0.1–1 Mb)")
    ax.bar(x, medium, bottom=short, color=COLORS["orange"], width=1.0, label="Medium (1–5 Mb)")
    ax.bar(x, long, bottom=short + medium, color=COLORS["purple"], width=1.0, label="Long (>5 Mb)")
    ax.set_xlabel("Individuals (ranked by ROH count)")
    ax.set_ylabel("ROH segments")
    ax.set_xticks([])
    ax.legend(frameon=False, ncol=1, loc="upper left", fontsize=4.5)
    despine(ax)


def main() -> None:
    fig, axes = plt.subplots(2, 2, figsize=(7.09, 5.5))
    plot_cohort_composition(axes[0, 0])
    plot_heterozygosity(axes[0, 1])
    plot_froh(axes[1, 0])
    plot_roh_classes(axes[1, 1])
    for ax, label in zip(axes.flat, "abcd"):
        add_panel_label(ax, label)
    fig.subplots_adjust(wspace=0.35, hspace=0.45)
    save_figure(fig, "figure_2_cohort_diversity_roh")


if __name__ == "__main__":
    main()
