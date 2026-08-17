"""Nature-style matplotlib theme (colorblind-friendly, 5–7 pt Arial)."""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt

# Okabe–Ito palette without red/green pairing
COLORS = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "teal": "#009E73",
    "yellow": "#F0E442",
    "grey": "#999999",
    "black": "#333333",
    "light_grey": "#E6E6E6",
}

PALETTE = [
    COLORS["blue"],
    COLORS["orange"],
    COLORS["sky"],
    COLORS["vermillion"],
    COLORS["purple"],
    COLORS["teal"],
]

MM_PER_INCH = 25.4
MAX_WIDTH_MM = 180
MAX_WIDTH_IN = MAX_WIDTH_MM / MM_PER_INCH


def apply_style() -> None:
    """Apply publication defaults once per process."""
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 6,
            "axes.titlesize": 7,
            "axes.labelsize": 6,
            "xtick.labelsize": 5,
            "ytick.labelsize": 5,
            "legend.fontsize": 5,
            "axes.linewidth": 0.4,
            "axes.edgecolor": COLORS["black"],
            "axes.facecolor": "white",
            "figure.facecolor": "white",
            "axes.grid": False,
            "xtick.major.width": 0.4,
            "ytick.major.width": 0.4,
            "xtick.major.size": 2,
            "ytick.major.size": 2,
            "lines.linewidth": 0.8,
            "patch.linewidth": 0.4,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 300,
            "savefig.facecolor": "white",
        }
    )


def figure(width_mm: float = 180, height_mm: float | None = None, nrows: int = 1, ncols: int = 1):
    apply_style()
    if height_mm is None:
        height_mm = width_mm * 0.75 * nrows / max(ncols, 1)
    width_in = min(width_mm, MAX_WIDTH_MM) / MM_PER_INCH
    height_in = height_mm / MM_PER_INCH
    return plt.figure(figsize=(width_in, height_in))


def add_panel_label(ax, label: str, x: float = -0.12, y: float = 1.08) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=7,
        fontweight="bold",
        va="top",
        ha="left",
    )


def despine(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
