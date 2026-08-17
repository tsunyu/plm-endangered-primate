"""Resolve the analysis root (directory containing data/ and output/)."""

from __future__ import annotations

import os
from pathlib import Path


def get_base_dir() -> Path:
    env = os.environ.get("PLM_BASE_DIR", "").strip()
    if env:
        return Path(env).expanduser().resolve()

    env_file = Path(__file__).resolve().parent / "base_dir.env"
    if env_file.is_file():
        for line in env_file.read_text(encoding="utf-8").splitlines():
            stripped = line.strip()
            if stripped.startswith("export "):
                stripped = stripped[len("export ") :].strip()
            if stripped.startswith("PLM_BASE_DIR="):
                value = stripped.split("=", 1)[1].strip().strip("'").strip('"')
                if value:
                    return Path(value).expanduser().resolve()

    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "data").is_dir() and (parent / "output").is_dir():
            return parent

    raise SystemExit(
        "Could not find the analysis root (a directory containing data/ and output/).\n"
        "Run:  bash configure_base_dir.sh /path/to/analysis_root\n"
        "  or: export PLM_BASE_DIR=/path/to/analysis_root"
    )


def analysis_package_dir() -> Path:
    """Directory that contains utils.py (this package's analysis/ folder)."""
    return Path(__file__).resolve().parent
