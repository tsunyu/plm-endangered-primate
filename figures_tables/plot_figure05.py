#!/usr/bin/env python3
"""Figure 5 — polarized deleterious derived allele-frequency spectra."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lib.daf import (
    DAF_CLASS_LABELS,
    annotate_deleterious_missense,
    daf_class_counts,
    polarized_subset,
    scan_variants,
    summarize_daf,
)
from lib.export import save_figure
from lib.paths import RESULTS
from lib.style import COLORS, despine

CATEGORY_STYLE = {
    "LoF": {"color": COLORS["blue"], "label": "Loss-of-function"},
    "Missense": {"color": COLORS["orange"], "label": "Top 10% missense"},
}


def _class_bar_data(subset: pd.DataFrame) -> pd.DataFrame:
    classes = daf_class_counts(subset["daf"].dropna())
    classes["percent"] = classes["fraction"] * 100
    return classes


def plot_grouped_daf_bars(ax, lof: pd.DataFrame, missense: pd.DataFrame, threshold: float) -> pd.DataFrame:
    lof_classes = _class_bar_data(lof)
    missense_classes = _class_bar_data(missense)
    x = np.arange(len(DAF_CLASS_LABELS))
    width = 0.36
    rows = []
    for offset, key, classes in [
        (-width / 2, "LoF", lof_classes),
        (width / 2, "Missense", missense_classes),
    ]:
        style = CATEGORY_STYLE[key]
        ax.bar(
            x + offset,
            classes["percent"],
            width,
            color=style["color"],
            label=style["label"],
            error_kw={"linewidth": 0.6},
        )
        for idx, row in classes.iterrows():
            ax.text(
                x[idx] + offset,
                row["percent"] + 1.2,
                f"{int(row['count'])}",
                ha="center",
                va="bottom",
                fontsize=4.5,
            )
            rows.append(
                {
                    "category": key,
                    "class": row["class"],
                    "count": int(row["count"]),
                    "percent": float(row["percent"]),
                    "n_total": int(classes["count"].sum()),
                }
            )
    ax.set_xticks(x)
    ax.set_xticklabels(DAF_CLASS_LABELS, rotation=30, ha="right")
    ax.set_xlim(-0.55, len(DAF_CLASS_LABELS) - 0.25)
    ax.set_ylabel("Variants (%)")
    ax.set_ylim(0, max(lof_classes["percent"].max(), missense_classes["percent"].max()) + 8)
    ax.legend(frameon=False, fontsize=5, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.02))
    ax.text(
        0.02,
        0.98,
        f"LoF n = {len(lof):,}\nMissense n = {len(missense):,}\nP ≥ {threshold:.3f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5,
    )
    despine(ax)
    return pd.DataFrame(rows)


def main() -> None:
    frame = scan_variants()
    frame, threshold = annotate_deleterious_missense(frame)
    summarize_daf(frame)
    polarized = polarized_subset(frame)
    lof = polarized[polarized["kind"] == "lof"]
    missense = polarized[
        (polarized["kind"] == "missense") & polarized["is_deleterious_missense"]
    ]

    fig, ax = plt.subplots(figsize=(7.09, 3.15))
    bar_df = plot_grouped_daf_bars(ax, lof, missense, threshold)
    ax.set_xlabel("Polarized derived allele frequency class", labelpad=4)
    bar_df.to_csv(RESULTS / "figure5_daf_grouped_bars.csv", index=False)
    pd.DataFrame(
        [
            {
                "threshold_quantile": 0.90,
                "calibrated_p_threshold": threshold,
                "n_lof_polarized": len(lof),
                "n_missense_polarized": len(missense),
                "n_scored_missense_all": int(
                    ((frame["kind"] == "missense") & frame["pathogenicity_prob"].notna()).sum()
                ),
            }
        ]
    ).to_csv(RESULTS / "figure5_deleterious_missense_threshold.csv", index=False)
    fig.subplots_adjust(left=0.11, right=0.94, bottom=0.34, top=0.80)
    save_figure(fig, "figure_5_polarized_daf_spectra", tight=False)


if __name__ == "__main__":
    main()
