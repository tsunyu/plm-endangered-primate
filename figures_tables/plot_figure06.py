#!/usr/bin/env python3
"""Figure 6 — demographic contraction, fixation timescales and forward simulation."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, ScalarFormatter

from lib.export import save_figure
from lib.paths import (
    BOOTSTRAP_RESULTS,
    DFE_REPLICATES,
    MODEL_COMPARISON,
    PARAM_ESTIMATES,
    RESULTS,
)
from lib.style import COLORS, add_panel_label, apply_style, despine

PARAM_ORDER = [
    "NCUR",
    "NRECOVER",
    "NBOT",
    "NANC",
    "TRECENT",
    "T_REC",
    "T_BOT",
    "TRECOVERY_OLD",
    "TBOT_OLD",
]
GENERATION_TIME = 10

MODEL_LABELS = {
    "bottleneck_recent_contraction": "Bottleneck +\nrecent contraction",
    "complex_multi_event": "Complex multi-event",
    "bottleneck_continuous_decline": "Continuous decline",
    "two_consecutive_bottlenecks": "Two bottlenecks",
    "single_bottleneck": "Single bottleneck",
    "constant_ne": "Constant Ne",
}


def load_params() -> dict[str, float]:
    df = pd.read_csv(PARAM_ESTIMATES)
    return dict(zip(df["Parameter"], df["Estimate"]))


def parse_bootstrap() -> pd.DataFrame:
    raw = pd.read_csv(BOOTSTRAP_RESULTS)
    rows = []
    for line in raw["parameters"]:
        values = [float(v) for v in str(line).split()]
        row = dict(zip(PARAM_ORDER, values[: len(PARAM_ORDER)]))
        rows.append(row)
    return pd.DataFrame(rows)


def ne_trajectory(params: dict[str, float], max_time_gen: int = 20000) -> tuple[np.ndarray, np.ndarray]:
    time_gen = np.linspace(0, max_time_gen, 1000)
    ne_values = np.zeros_like(time_gen)
    trecent = params["TRECENT"]
    trecovery_old = params["TRECOVERY_OLD"]
    tbot_old = params["TBOT_OLD"]
    ne_values[time_gen <= trecent] = params["NCUR"]
    ne_values[(time_gen > trecent) & (time_gen <= trecovery_old)] = params["NRECOVER"]
    ne_values[(time_gen > trecovery_old) & (time_gen <= tbot_old)] = params["NBOT"]
    ne_values[time_gen > tbot_old] = params["NANC"]
    return time_gen * GENERATION_TIME, ne_values


def fixation_ratios(params: dict[str, float]) -> tuple[float, float]:
    ncur = params["NCUR"]
    nbot = params["NBOT"]
    trecent = params["TRECENT"]
    trecovery = params["TRECOVERY_OLD"]
    tbot = params["TBOT_OLD"]
    recent = trecent / (4 * ncur)
    ancient = (tbot - trecovery) / (4 * nbot)
    return ancient, recent


def bootstrap_fixation(boot: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    valid = boot[boot["T_BOT"] > boot["T_REC"]].copy()
    ancient = (valid["T_BOT"] - valid["T_REC"]) / (4 * valid["NBOT"])
    recent = valid["TRECENT"] / (4 * valid["NCUR"])
    return ancient, recent


def plot_aic(ax) -> None:
    comp = pd.read_csv(MODEL_COMPARISON).sort_values("DeltaAIC")
    labels = [MODEL_LABELS.get(m, m.replace("_", " ")) for m in comp["Model"]]
    deltas = comp["DeltaAIC"].to_numpy()
    colors = [COLORS["blue"] if d == 0 else COLORS["grey"] for d in deltas]
    y = np.arange(len(comp))
    ax.barh(y, deltas, color=colors, height=0.62)
    ax.scatter(
        [0],
        [0],
        marker="o",
        s=12,
        facecolor=COLORS["blue"],
        edgecolor="white",
        linewidth=0.4,
        zorder=4,
        clip_on=False,
    )
    ax.text(220, 0, "0 (best)", va="center", ha="left", fontsize=4.5, color=COLORS["blue"])
    for yi, delta in zip(y[1:], deltas[1:]):
        label = f"{delta / 1e3:.2f}k" if delta < 1e6 else f"{delta / 1e6:.2f}M"
        ax.text(delta * 1.08, yi, label, va="center", ha="left", fontsize=4.2)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xscale("symlog", linthresh=1000, linscale=1.2)
    ax.set_xticks([0, 1e3, 1e4, 1e5, 1e6])
    ax.set_xticklabels(["0", "10³", "10⁴", "10⁵", "10⁶"])
    ax.set_xlabel(r"$\Delta$AIC (symmetric log scale)")
    ax.set_xlim(0, 4e6)
    ax.invert_yaxis()
    despine(ax)


def plot_ne_trajectory(ax, params: dict[str, float]) -> None:
    time_y, ne = ne_trajectory(params)
    time_kya = np.maximum(time_y / 1000.0, 0.05)
    ax.step(time_kya, ne, where="post", color=COLORS["blue"], linewidth=0.9, zorder=3)
    ax.fill_between(
        time_kya,
        ne,
        ne.min() * 0.85,
        step="post",
        color=COLORS["sky"],
        alpha=0.15,
        zorder=1,
    )
    trecent_kya = params["TRECENT"] * GENERATION_TIME / 1000.0
    ax.axvline(trecent_kya, color=COLORS["vermillion"], linestyle=":", linewidth=0.6, alpha=0.85)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(max(0.05, float(time_kya[time_kya > 0].min())), float(time_kya.max()))
    ax.set_ylim(float(ne[ne > 0].min()) * 0.85, float(ne.max()) * 1.25)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.set_xlabel("Time before present (kya)")
    ax.set_ylabel(r"Diploid $N_e$")
    ax.text(
        0.03,
        0.97,
        f"Present $N_e$ = {params['NCUR']:.0f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5,
    )
    despine(ax)


def plot_fixation_timescales(ax, params: dict[str, float], boot: pd.DataFrame) -> None:
    ancient, recent = fixation_ratios(params)
    boot_ancient, boot_recent = bootstrap_fixation(boot)
    means = [boot_ancient.median(), boot_recent.median()]
    lowers = [np.percentile(boot_ancient, 2.5), np.percentile(boot_recent, 2.5)]
    uppers = [np.percentile(boot_ancient, 97.5), np.percentile(boot_recent, 97.5)]
    labels = ["Ancient", "Recent"]
    x = np.arange(2)
    ax.bar(x, means, color=[COLORS["purple"], COLORS["orange"]], width=0.55)
    ax.errorbar(
        x,
        means,
        yerr=[np.array(means) - np.array(lowers), np.array(uppers) - np.array(means)],
        fmt="none",
        ecolor=COLORS["black"],
        capsize=2,
        linewidth=0.6,
    )
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Episode duration / 4$N_e$")
    ax.set_ylim(0, max(uppers) * 1.2)
    despine(ax)
    pd.DataFrame(
        {
            "episode": ["Ancient bottleneck", "Recent contraction"],
            "point_estimate": [ancient, recent],
            "bootstrap_median": means,
            "ci_lower": lowers,
            "ci_upper": uppers,
        }
    ).to_csv(RESULTS / "figure6_fixation_timescales.csv", index=False)


def _format_fraction_label(value: float) -> str:
    if value < 1e-6:
        return "0"
    if value < 0.001:
        return f"{value * 100:.3f}%"
    return f"{value * 100:.2f}%"


def plot_simulation_summary(ax) -> None:
    """Compare simulated fixed vs intermediate fractions under two demographies.

    Empirical spectrum is not plotted here: its intermediate class (~35%) is an
    order of magnitude larger than the simulated fractions, so joint axes or
    insets obscure the contrast this panel is meant to show.
    """
    reps = pd.read_csv(DFE_REPLICATES)
    focus = reps[reps["scenario"].isin(["transient_drift", "constant_current_ne"])]
    metrics = focus.groupby("scenario")[["frac_fixed", "frac_intermediate"]].agg(["mean", "std"])
    scenarios = ["transient_drift", "constant_current_ne"]
    labels = ["Transient drift", "Constant\n$N_e$ = 150"]
    width = 0.34
    x = np.arange(2)
    fixed = [metrics.loc[s, ("frac_fixed", "mean")] for s in scenarios]
    fixed_err = [metrics.loc[s, ("frac_fixed", "std")] for s in scenarios]
    inter = [metrics.loc[s, ("frac_intermediate", "mean")] for s in scenarios]
    inter_err = [metrics.loc[s, ("frac_intermediate", "std")] for s in scenarios]

    ax.bar(
        x - width / 2,
        fixed,
        width,
        yerr=fixed_err,
        label="Fixed",
        color=COLORS["blue"],
        capsize=2,
        error_kw={"linewidth": 0.6},
    )
    ax.bar(
        x + width / 2,
        inter,
        width,
        yerr=inter_err,
        label="Intermediate",
        color=COLORS["orange"],
        capsize=2,
        error_kw={"linewidth": 0.6},
    )

    # Near-zero bars are invisible on a linear scale; mark them explicitly.
    visibility_floor = 0.002
    ymax = max(
        max(f + e for f, e in zip(fixed, fixed_err)),
        max(i + e for i, e in zip(inter, inter_err)),
        visibility_floor,
    ) * 1.45
    for xi, value, err, color, offset in [
        *[(x[i], fixed[i], fixed_err[i], COLORS["blue"], -width / 2) for i in range(2)],
        *[(x[i], inter[i], inter_err[i], COLORS["orange"], width / 2) for i in range(2)],
    ]:
        xpos = xi + offset
        label = _format_fraction_label(value)
        if value < visibility_floor:
            ax.scatter(
                [xpos],
                [0],
                marker="v",
                s=14,
                color=color,
                edgecolor="white",
                linewidth=0.3,
                zorder=4,
                clip_on=False,
            )
            ax.text(
                xpos,
                ymax * 0.08,
                label,
                ha="center",
                va="bottom",
                fontsize=4.5,
                color=color,
            )
        else:
            ax.text(
                xpos,
                value + err + ymax * 0.03,
                label,
                ha="center",
                va="bottom",
                fontsize=4.5,
                color=COLORS["black"],
            )

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Simulated class fraction")
    ax.set_ylim(0, ymax)
    ax.legend(frameon=False, fontsize=5, loc="upper left", handlelength=1.0, handletextpad=0.4)
    despine(ax)


def main() -> None:
    apply_style()
    params = load_params()
    boot = parse_bootstrap()
    fig, axes = plt.subplots(2, 2, figsize=(7.09, 5.8))
    plot_aic(axes[0, 0])
    plot_ne_trajectory(axes[0, 1], params)
    plot_fixation_timescales(axes[1, 0], params, boot)
    plot_simulation_summary(axes[1, 1])
    for ax, label in zip(axes.flat, "abcd"):
        col = 0 if label in {"a", "c"} else 1
        add_panel_label(ax, label, x=-0.20 if col == 0 else -0.14, y=1.04)
    fig.subplots_adjust(left=0.26, right=0.98, bottom=0.16, top=0.95, wspace=0.38, hspace=0.48)
    save_figure(fig, "figure_6_demography_simulation", tight=False)


if __name__ == "__main__":
    main()
