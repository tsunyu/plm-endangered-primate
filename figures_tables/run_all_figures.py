#!/usr/bin/env python3
"""Redraw publication figures and tables (except Figure 1)."""

from __future__ import annotations

import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from lib.style import apply_style


def run(module: str) -> None:
    print(f"\n=== {module} ===")
    start = time.time()
    mod = __import__(module)
    if hasattr(mod, "main"):
        mod.main()
    print(f"Completed in {time.time() - start:.1f}s")


def main() -> None:
    apply_style()
    modules = [
        "plot_figure02",
        "plot_figure03",
        "plot_figure04",
        "plot_figure05",
        "plot_figure06",
        "plot_supplementary",
        "regenerate_tables",
    ]
    for module in modules:
        run(module)
    print("\nAll figures and tables complete. Figure 1 remains manual.")


if __name__ == "__main__":
    main()
