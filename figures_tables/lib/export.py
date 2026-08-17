"""Export figures as editable PDF and RGB preview PNG under works/figures only."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt

from .paths import FIGURES, SUPP_FIGURES


def save_figure(
    fig: plt.Figure,
    stem: str,
    *,
    supplementary: bool = False,
    tight: bool = True,
    pad_inches: float = 0.02,
) -> dict[str, Path]:
    """Save figure under works/figures (or works/figures/supplementary)."""
    base = SUPP_FIGURES if supplementary else FIGURES
    outputs: dict[str, Path] = {}
    kwargs = {"bbox_inches": "tight", "pad_inches": pad_inches} if tight else {}
    for ext in ("pdf", "png"):
        path = base / f"{stem}.{ext}"
        fig.savefig(path, format=ext, dpi=300, **kwargs)
        outputs[ext] = path
    plt.close(fig)
    return outputs
