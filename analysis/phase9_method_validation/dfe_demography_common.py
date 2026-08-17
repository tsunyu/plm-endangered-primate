#!/usr/bin/env python3
"""Shared utilities for transient-drift-load forward simulations."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir

import numpy as np
import pandas as pd

BASE = get_base_dir()
PARAM_FILE = BASE / "output/phase3b_fastsimcoal2/parameter_estimates.csv"
OBS_STATS_FILE = BASE / "output/allele_frequency_spectrum/frequency_statistics.csv"

BIN_EDGES = [0.0, 0.05, 0.20, 0.50, 0.95, 1.01]
BIN_LABELS = ["rare", "low", "intermediate", "common", "fixed"]


@dataclass(frozen=True)
class Demography:
    nanc: float
    nbot: float
    nrecover: float
    ncur: float
    tbot_old: int
    trecovery_old: int
    trecent: int
    burn_in: int = 2000

    @property
    def total_generations(self) -> int:
        return self.burn_in + self.tbot_old

    def ne_at(self, generation: int) -> float:
        if generation < self.burn_in:
            return self.nanc
        since_bottle = generation - self.burn_in
        if since_bottle < self.tbot_old - self.trecovery_old:
            return self.nbot
        if since_bottle < self.tbot_old - self.trecent:
            return self.nrecover
        return self.ncur


def load_demography_mle() -> Demography:
    params = pd.read_csv(PARAM_FILE).set_index("Parameter")["Estimate"]
    return Demography(
        nanc=float(params["NANC"]),
        nbot=float(params["NBOT"]),
        nrecover=float(params["NRECOVER"]),
        ncur=float(params["NCUR"]),
        tbot_old=int(round(float(params["TBOT_OLD"]))),
        trecovery_old=int(round(float(params["TRECOVERY_OLD"]))),
        trecent=int(round(float(params["TRECENT"]))),
    )


def load_observed_bins() -> pd.DataFrame:
    stats = pd.read_csv(OBS_STATS_FILE, index_col=0)
    rows = []
    for label in ("LoF", "Deleterious Missense"):
        row = stats.loc[label]
        total = row["total"]
        for bin_label in BIN_LABELS:
            rows.append(
                {
                    "variant_class": label,
                    "bin": bin_label,
                    "count": int(row[bin_label]),
                    "fraction": row[bin_label] / total,
                }
            )
        rows.append(
            {
                "variant_class": label,
                "bin": "mean_daf",
                "count": np.nan,
                "fraction": row["mean_daf"],
            }
        )
    return pd.DataFrame(rows)


def bin_frequencies(freqs: np.ndarray) -> dict[str, float]:
    freqs = np.asarray(freqs, dtype=float)
    freqs = freqs[np.isfinite(freqs)]
    if freqs.size == 0:
        return {label: 0.0 for label in BIN_LABELS}
    counts = np.histogram(freqs, bins=BIN_EDGES)[0]
    total = counts.sum()
    if total == 0:
        return {label: 0.0 for label in BIN_LABELS}
    return {
        label: count / total
        for label, count in zip(BIN_LABELS, counts, strict=True)
    }


def mean_fitness(p: float, s: float, h: float) -> float:
    return 1.0 - 2.0 * p * (1.0 - p) * h * s - (p**2) * s


def post_selection_frequency(p: float, s: float, h: float) -> float:
    w_aa = 1.0 - s
    w_aA = 1.0 - h * s
    numerator = (p**2) * w_aa + p * (1.0 - p) * w_aA
    denominator = mean_fitness(p, s, h)
    if denominator <= 0:
        return 0.0
    return float(np.clip(numerator / denominator, 0.0, 1.0))


def wright_fisher_step(p: float, ne: float, s: float, h: float, rng: np.random.Generator) -> float:
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0
    p_sel = post_selection_frequency(p, s, h)
    count = rng.binomial(int(round(2.0 * ne)), p_sel)
    return count / (2.0 * ne)


def simulate_locus(
    demography: Demography,
    s: float,
    h: float,
    rng: np.random.Generator,
    init_p: float | None = None,
) -> float:
    if init_p is None:
        eq = demography.nanc
        init_p = min(max(rng.gamma(2.0, 0.5) * 1e-4 / max(s, 1e-6), 1e-6), 0.05)
        _ = eq
    p = init_p
    for gen in range(demography.total_generations):
        ne = demography.ne_at(gen)
        p = wright_fisher_step(p, ne, s, h, rng)
        if p <= 0.0 or p >= 1.0:
            break
    return p


def post_selection_frequency_vec(p: np.ndarray, s: np.ndarray, h: float) -> np.ndarray:
    w_aa = 1.0 - s
    w_aA = 1.0 - h * s
    numerator = (p**2) * w_aa + p * (1.0 - p) * w_aA
    denominator = 1.0 - 2.0 * p * (1.0 - p) * h * s - (p**2) * s
    return np.clip(np.divide(numerator, denominator, out=np.zeros_like(p), where=denominator > 0), 0.0, 1.0)


def simulate_replicate(
    demography: Demography,
    n_loci: int,
    s_mean: float,
    s_shape: float,
    h: float,
    rng: np.random.Generator,
    seed_standing: bool = False,
) -> np.ndarray:
    scales = s_mean / max(s_shape, 1e-8)
    s_values = np.clip(rng.gamma(s_shape, scales, size=n_loci), 1e-6, 0.05)
    p = np.clip(rng.gamma(2.0, 0.5, size=n_loci) * 5e-5 / np.maximum(s_values, 1e-8), 1e-6, 0.02)
    if seed_standing:
        standing = rng.uniform(0.0, 1.0, size=n_loci) < 0.35
        p[standing] = rng.uniform(0.05, 0.45, size=int(standing.sum()))
    active = np.ones(n_loci, dtype=bool)
    for gen in range(demography.total_generations):
        if not active.any():
            break
        ne = demography.ne_at(gen)
        p_sel = post_selection_frequency_vec(p, s_values, h)
        counts = rng.binomial(max(int(round(2.0 * ne)), 2), p_sel)
        p = counts / (2.0 * ne)
        active = active & (p > 0.0) & (p < 1.0)
    return p


def summarize_replicate(freqs: np.ndarray) -> dict[str, float]:
    bins = bin_frequencies(freqs)
    summary = {f"frac_{key}": value for key, value in bins.items()}
    summary["mean_daf"] = float(np.mean(freqs)) if freqs.size else np.nan
    summary["frac_polymorphic"] = float(np.mean((freqs > 0.0) & (freqs < 1.0)))
    summary["frac_fixed"] = float(np.mean(freqs >= 0.95))
    return summary


def demography_table(demography: Demography) -> pd.DataFrame:
    rows = []
    for gen in range(0, demography.total_generations + 1, 25):
        rows.append({"generation": gen, "ne": demography.ne_at(min(gen, demography.total_generations - 1))})
    rows.append({"generation": demography.total_generations, "ne": demography.ncur})
    return pd.DataFrame(rows).drop_duplicates(subset=["generation"])


def write_metadata(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
