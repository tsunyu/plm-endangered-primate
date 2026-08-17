"""Bootstrap confidence intervals for manuscript forest plots."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

from .paths import BOOTSTRAP_SEED, N_BOOTSTRAP


def bootstrap_spearman_ci(
    frame: pd.DataFrame, x_col: str, y_col: str, n_boot: int = N_BOOTSTRAP, seed: int = BOOTSTRAP_SEED
) -> tuple[float, float, float, float]:
    rng = np.random.default_rng(seed)
    observed = spearmanr(frame[x_col], frame[y_col], nan_policy="omit").statistic
    boot = []
    n = len(frame)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        sample = frame.iloc[idx]
        boot.append(spearmanr(sample[x_col], sample[y_col], nan_policy="omit").statistic)
    lower, upper = np.nanpercentile(boot, [2.5, 97.5])
    return float(observed), float(lower), float(upper), float(np.nanstd(boot))


def bootstrap_cohens_d_ci(
    frame: pd.DataFrame,
    value_col: str,
    group_col: str = "Has_Disease",
    n_boot: int = N_BOOTSTRAP,
    seed: int = BOOTSTRAP_SEED,
) -> tuple[float, float, float, float, float]:
    cases = frame[frame[group_col] == 1][value_col]
    controls = frame[frame[group_col] == 0][value_col]
    pooled_sd = np.sqrt(
        ((len(cases) - 1) * cases.std(ddof=1) ** 2 + (len(controls) - 1) * controls.std(ddof=1) ** 2)
        / (len(cases) + len(controls) - 2)
    )
    observed = (cases.mean() - controls.mean()) / pooled_sd if pooled_sd > 0 else np.nan
    _, p_value = mannwhitneyu(cases, controls, alternative="two-sided")
    rng = np.random.default_rng(seed)
    boot = []
    n = len(frame)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        sample = frame.iloc[idx]
        c = sample[sample[group_col] == 1][value_col]
        h = sample[sample[group_col] == 0][value_col]
        sd = np.sqrt(
            ((len(c) - 1) * c.std(ddof=1) ** 2 + (len(h) - 1) * h.std(ddof=1) ** 2)
            / (len(c) + len(h) - 2)
        )
        boot.append((c.mean() - h.mean()) / sd if sd > 0 else np.nan)
    lower, upper = np.nanpercentile(boot, [2.5, 97.5])
    return float(observed), float(lower), float(upper), float(p_value), float(np.nanstd(boot))
