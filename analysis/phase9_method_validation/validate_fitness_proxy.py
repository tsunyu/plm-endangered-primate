#!/usr/bin/env python3
"""Validate genomic load as a morbidity proxy while accounting for relatedness.

The primary outcome is the Composite Health Score (CHS), a field-recorded
morbidity composite rather than survival or reproductive fitness. Age is
deliberately excluded from primary models because it is unavailable for most
animals and is recorded mainly for cases. The binary disease analysis is
exploratory and is reported as discrimination, not as a prevalence model.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize_scalar
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    mean_absolute_error,
    mean_squared_error,
    r2_score,
    roc_auc_score,
)
from sklearn.model_selection import GroupKFold


BASE = get_base_dir()
SEED = 20260710
RELATIONSHIP_THRESHOLD = 0.125
MODEL_TERMS = {
    "covariates": [],
    "covariates_froh": ["F_ROH_z"],
    "covariates_load": ["Total_Genetic_Load_z"],
    "covariates_both": ["F_ROH_z", "Total_Genetic_Load_z"],
}


@dataclass
class LMMFit:
    beta: np.ndarray
    se: np.ndarray
    pvalue: np.ndarray
    fitted: np.ndarray
    residual: np.ndarray
    loglik: float
    aic: float
    bic: float
    delta: float
    h2: float
    sigma2_g: float
    covariance: np.ndarray
    v: np.ndarray


def read_ids(path: Path) -> list[str]:
    frame = pd.read_csv(path, sep=r"\s+", header=None, dtype=str)
    if frame.shape[1] < 2:
        raise ValueError(f"Expected FID and IID in {path}")
    ids = frame.iloc[:, 1].astype(str).tolist()
    if len(ids) != len(set(ids)):
        raise ValueError(f"Duplicate IIDs in {path}")
    return ids


def load_gcta_grm(prefix: Path) -> tuple[np.ndarray, list[str], dict]:
    id_path = Path(f"{prefix}.grm.id")
    bin_path = Path(f"{prefix}.grm.bin")
    n_path = Path(f"{prefix}.grm.N.bin")
    if not all(path.exists() for path in (id_path, bin_path, n_path)):
        raise FileNotFoundError("Complete GCTA GRM triplet is not present")
    ids = read_ids(id_path)
    n = len(ids)
    packed = np.fromfile(bin_path, dtype=np.float32)
    expected = n * (n + 1) // 2
    if packed.size != expected:
        raise ValueError(f"GCTA GRM has {packed.size} values; expected {expected}")
    matrix = np.zeros((n, n), dtype=float)
    lower = np.tril_indices(n)
    matrix[lower] = packed
    matrix = matrix + matrix.T - np.diag(np.diag(matrix))
    return matrix, ids, {
        "source_type": "GCTA_GRM",
        "matrix_path": str(bin_path),
        "id_path": str(id_path),
        "n_path": str(n_path),
    }


def load_gemma_kinship(path: Path, fam_path: Path) -> tuple[np.ndarray, list[str], dict]:
    matrix = np.loadtxt(path)
    ids = read_ids(fam_path)
    if matrix.shape != (len(ids), len(ids)):
        raise ValueError("GEMMA kinship dimensions do not match FAM IDs")
    return matrix, ids, {
        "source_type": "GEMMA_centered_K",
        "matrix_path": str(path),
        "id_path": str(fam_path),
    }


def choose_kinship(args: argparse.Namespace) -> tuple[np.ndarray, list[str], dict]:
    try:
        return load_gcta_grm(args.gcta_prefix)
    except FileNotFoundError:
        return load_gemma_kinship(args.gemma_kinship, args.gemma_fam)


def regularize_psd(matrix: np.ndarray) -> tuple[np.ndarray, dict]:
    matrix = (matrix + matrix.T) / 2
    values, vectors = np.linalg.eigh(matrix)
    floor = max(1e-8, float(values.max()) * 1e-8)
    corrected = np.maximum(values, floor)
    output = (vectors * corrected) @ vectors.T
    output = (output + output.T) / 2
    return output, {
        "minimum_eigenvalue_raw": float(values.min()),
        "maximum_eigenvalue_raw": float(values.max()),
        "eigenvalue_floor": float(floor),
        "eigenvalues_adjusted": int(np.sum(values < floor)),
        "maximum_absolute_adjustment": float(np.max(np.abs(output - matrix))),
    }


def align_data(
    phenotype_path: Path, kinship: np.ndarray, kinship_ids: list[str]
) -> tuple[pd.DataFrame, np.ndarray, pd.DataFrame]:
    frame = pd.read_csv(phenotype_path, dtype={"IID": str})
    required = {
        "IID",
        "Sex",
        "CHS",
        "Has_Disease",
        "Age",
        "F_ROH",
        "Total_Genetic_Load",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Missing phenotype columns: {sorted(missing)}")
    if frame["IID"].duplicated().any():
        raise ValueError("Duplicate IIDs in phenotype table")

    index = {iid: i for i, iid in enumerate(kinship_ids)}
    phenotype_ids = set(frame["IID"])
    common = [iid for iid in kinship_ids if iid in phenotype_ids]
    if not common:
        raise ValueError("No overlapping IIDs between phenotype data and GRM")
    frame_indexed = frame.set_index("IID")
    aligned = frame_indexed.loc[common].reset_index()
    matrix_index = [index[iid] for iid in common]
    aligned_kinship = kinship[np.ix_(matrix_index, matrix_index)]

    complete = aligned[
        ["Sex", "CHS", "Has_Disease", "F_ROH", "Total_Genetic_Load"]
    ].notna().all(axis=1)
    status = []
    for iid in sorted(set(kinship_ids).union(phenotype_ids)):
        in_grm = iid in index
        in_pheno = iid in phenotype_ids
        used = iid in set(aligned.loc[complete, "IID"])
        reason = "included" if used else (
            "missing_from_grm" if not in_grm else
            "missing_from_phenotype" if not in_pheno else
            "missing_primary_variable"
        )
        status.append(
            {"IID": iid, "in_grm": in_grm, "in_phenotype": in_pheno, "used": used, "status": reason}
        )
    keep = np.flatnonzero(complete.to_numpy())
    return (
        aligned.iloc[keep].reset_index(drop=True),
        aligned_kinship[np.ix_(keep, keep)],
        pd.DataFrame(status),
    )


def prepare_data(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    frame = frame.copy()
    sex_levels = sorted(frame["Sex"].astype(str).unique())
    if not set(sex_levels).issubset({"Female", "Male"}):
        raise ValueError(f"Unexpected Sex values: {sex_levels}")
    frame["Sex_Male"] = (frame["Sex"] == "Male").astype(float)
    scaling = {}
    for column in ("F_ROH", "Total_Genetic_Load"):
        mean = float(frame[column].mean())
        sd = float(frame[column].std(ddof=1))
        if not np.isfinite(sd) or sd == 0:
            raise ValueError(f"{column} has zero or invalid variance")
        frame[f"{column}_z"] = (frame[column] - mean) / sd
        scaling[column] = {"mean": mean, "sd": sd}
    return frame, scaling


def design_matrix(frame: pd.DataFrame, terms: list[str]) -> tuple[np.ndarray, list[str]]:
    names = ["Intercept", "Sex_Male", *terms]
    columns = [np.ones(len(frame)), frame["Sex_Male"].to_numpy(dtype=float)]
    columns.extend(frame[term].to_numpy(dtype=float) for term in terms)
    matrix = np.column_stack(columns)
    if np.linalg.matrix_rank(matrix) < matrix.shape[1]:
        raise ValueError(f"Rank-deficient design: {names}")
    return matrix, names


def fit_lmm(y: np.ndarray, x: np.ndarray, kinship: np.ndarray) -> LMMFit:
    """Fit y = Xb + u + e by profile ML and return the equivalent GLS fit."""
    n, p = x.shape
    eigenvalues, eigenvectors = np.linalg.eigh((kinship + kinship.T) / 2)
    eigenvalues = np.maximum(eigenvalues, 0)
    yt = eigenvectors.T @ y
    xt = eigenvectors.T @ x

    def evaluate(log_delta: float, full: bool = False):
        delta = math.exp(log_delta)
        variance = eigenvalues + delta
        weights = 1.0 / variance
        xtwx = xt.T @ (weights[:, None] * xt)
        inverse = np.linalg.pinv(xtwx)
        beta = inverse @ (xt.T @ (weights * yt))
        residual_t = yt - xt @ beta
        quadratic = float(np.sum(weights * residual_t**2))
        sigma2 = max(quadratic / n, np.finfo(float).tiny)
        loglik = -0.5 * (
            n * (math.log(2 * math.pi) + 1 + math.log(sigma2))
            + float(np.log(variance).sum())
        )
        if full:
            return delta, variance, beta, inverse, sigma2, loglik
        return -loglik

    optimum = minimize_scalar(evaluate, bounds=(-12, 12), method="bounded")
    delta, variance, beta, inverse, sigma2_g, loglik = evaluate(optimum.x, full=True)
    covariance = sigma2_g * inverse
    se = np.sqrt(np.maximum(np.diag(covariance), 0))
    zscore = np.divide(beta, se, out=np.zeros_like(beta), where=se > 0)
    pvalue = 2 * stats.norm.sf(np.abs(zscore))
    fitted = x @ beta
    residual = y - fitted
    v = sigma2_g * (
        kinship + delta * np.eye(n)
    )
    parameter_count = p + 2
    return LMMFit(
        beta=beta,
        se=se,
        pvalue=pvalue,
        fitted=fitted,
        residual=residual,
        loglik=float(loglik),
        aic=float(-2 * loglik + 2 * parameter_count),
        bic=float(-2 * loglik + math.log(n) * parameter_count),
        delta=float(delta),
        h2=float(1 / (1 + delta)),
        sigma2_g=float(sigma2_g),
        covariance=covariance,
        v=v,
    )


def connected_components(kinship: np.ndarray, threshold: float) -> np.ndarray:
    adjacency = kinship >= threshold
    np.fill_diagonal(adjacency, False)
    groups = np.full(len(kinship), -1, dtype=int)
    group = 0
    for start in range(len(kinship)):
        if groups[start] >= 0:
            continue
        stack = [start]
        groups[start] = group
        while stack:
            current = stack.pop()
            for neighbor in np.flatnonzero(adjacency[current]):
                if groups[neighbor] < 0:
                    groups[neighbor] = group
                    stack.append(int(neighbor))
        group += 1
    return groups


def fit_primary_models(
    frame: pd.DataFrame, kinship: np.ndarray
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, LMMFit]]:
    y = frame["CHS"].to_numpy(dtype=float)
    effect_rows = []
    fit_rows = []
    fits = {}
    for model_name, terms in MODEL_TERMS.items():
        x, names = design_matrix(frame, terms)
        fit = fit_lmm(y, x, kinship)
        fits[model_name] = fit
        fit_rows.append(
            {
                "outcome": "CHS",
                "model": model_name,
                "n": len(y),
                "fixed_effect_parameters": x.shape[1],
                "log_likelihood_ml": fit.loglik,
                "AIC": fit.aic,
                "BIC": fit.bic,
                "delta_residual_to_genetic": fit.delta,
                "h2_profile_ml": fit.h2,
            }
        )
        for term, beta, se, pvalue in zip(names, fit.beta, fit.se, fit.pvalue):
            effect_rows.append(
                {
                    "outcome": "CHS",
                    "model": model_name,
                    "term": term,
                    "beta": beta,
                    "SE": se,
                    "z": beta / se if se else np.nan,
                    "p_wald": pvalue,
                }
            )
    return pd.DataFrame(effect_rows), pd.DataFrame(fit_rows), fits


def likelihood_comparisons(fits: dict[str, LMMFit]) -> pd.DataFrame:
    comparisons = [
        ("covariates", "covariates_froh", 1),
        ("covariates", "covariates_load", 1),
        ("covariates", "covariates_both", 2),
        ("covariates_froh", "covariates_both", 1),
        ("covariates_load", "covariates_both", 1),
    ]
    rows = []
    for reduced, full, df in comparisons:
        statistic = max(0.0, 2 * (fits[full].loglik - fits[reduced].loglik))
        rows.append(
            {
                "outcome": "CHS",
                "reduced_model": reduced,
                "full_model": full,
                "df": df,
                "likelihood_ratio": statistic,
                "p_chi2": stats.chi2.sf(statistic, df),
                "delta_AIC_full_minus_reduced": fits[full].aic - fits[reduced].aic,
                "delta_BIC_full_minus_reduced": fits[full].bic - fits[reduced].bic,
            }
        )
    return pd.DataFrame(rows)


def grouped_cv(
    frame: pd.DataFrame,
    kinship: np.ndarray,
    groups: np.ndarray,
    folds: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    unique_groups = np.unique(groups)
    n_splits = min(folds, len(unique_groups))
    if n_splits < 2:
        raise ValueError("Need at least two relatedness clusters for grouped CV")
    splitter = GroupKFold(n_splits=n_splits)
    continuous_rows = []
    binary_rows = []
    prediction_rows = []
    for model_name, terms in MODEL_TERMS.items():
        prediction = np.full(len(frame), np.nan)
        for fold, (train, test) in enumerate(
            splitter.split(frame, groups=groups), start=1
        ):
            x_train, _ = design_matrix(frame.iloc[train], terms)
            x_test, _ = design_matrix(frame.iloc[test], terms)
            fit = fit_lmm(
                frame.iloc[train]["CHS"].to_numpy(dtype=float),
                x_train,
                kinship[np.ix_(train, train)],
            )
            prediction[test] = x_test @ fit.beta
            for index in test:
                prediction_rows.append(
                    {
                        "IID": frame.iloc[index]["IID"],
                        "model": model_name,
                        "fold": fold,
                        "relatedness_cluster": int(groups[index]),
                        "CHS_observed": frame.iloc[index]["CHS"],
                        "CHS_prediction": prediction[index],
                        "Has_Disease": int(frame.iloc[index]["Has_Disease"]),
                    }
                )
        observed = frame["CHS"].to_numpy(dtype=float)
        disease = frame["Has_Disease"].to_numpy(dtype=int)
        rho, rho_p = stats.spearmanr(observed, prediction)
        continuous_rows.append(
            {
                "model": model_name,
                "n": len(frame),
                "folds": n_splits,
                "R2": r2_score(observed, prediction),
                "RMSE": math.sqrt(mean_squared_error(observed, prediction)),
                "MAE": mean_absolute_error(observed, prediction),
                "spearman_rho": rho,
                "spearman_p_descriptive": rho_p,
            }
        )
        binary_prediction = (prediction >= 0.5).astype(int)
        binary_rows.append(
            {
                "model": model_name,
                "n": len(frame),
                "cases": int(disease.sum()),
                "folds": n_splits,
                "roc_auc": roc_auc_score(disease, prediction),
                "pr_auc": average_precision_score(disease, prediction),
                "brier_on_clipped_CHS_score": brier_score_loss(
                    disease, np.clip(prediction, 0, 1)
                ),
                "balanced_accuracy_threshold_CHS_0.5": balanced_accuracy_score(
                    disease, binary_prediction
                ),
            }
        )
    return (
        pd.DataFrame(continuous_rows),
        pd.DataFrame(binary_rows),
        pd.DataFrame(prediction_rows),
    )


def whiten_matrix(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    values, vectors = np.linalg.eigh((matrix + matrix.T) / 2)
    floor = max(float(values.max()) * 1e-10, 1e-12)
    values = np.maximum(values, floor)
    root = (vectors * np.sqrt(values)) @ vectors.T
    inverse_root = (vectors * (1 / np.sqrt(values))) @ vectors.T
    return root, inverse_root


def partial_f(y: np.ndarray, reduced: np.ndarray, full: np.ndarray) -> float:
    residual_reduced = y - reduced @ np.linalg.lstsq(reduced, y, rcond=None)[0]
    residual_full = y - full @ np.linalg.lstsq(full, y, rcond=None)[0]
    rss_reduced = float(residual_reduced @ residual_reduced)
    rss_full = float(residual_full @ residual_full)
    q = full.shape[1] - reduced.shape[1]
    denominator_df = len(y) - full.shape[1]
    return max(0.0, ((rss_reduced - rss_full) / q) / (rss_full / denominator_df))


def freedman_lane_test(
    y: np.ndarray,
    reduced_x: np.ndarray,
    full_x: np.ndarray,
    kinship: np.ndarray,
    iterations: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    null_fit = fit_lmm(y, reduced_x, kinship)
    _, inverse_root = whiten_matrix(null_fit.v)
    yw = inverse_root @ y
    reduced_w = inverse_root @ reduced_x
    full_w = inverse_root @ full_x
    reduced_beta = np.linalg.lstsq(reduced_w, yw, rcond=None)[0]
    fitted = reduced_w @ reduced_beta
    residual = yw - fitted
    observed = partial_f(yw, reduced_w, full_w)
    exceed = 0
    for _ in range(iterations):
        permuted_y = fitted + rng.permutation(residual)
        exceed += partial_f(permuted_y, reduced_w, full_w) >= observed
    return observed, (exceed + 1) / (iterations + 1)


def permutation_tests(
    frame: pd.DataFrame, kinship: np.ndarray, iterations: int
) -> pd.DataFrame:
    y = frame["CHS"].to_numpy(dtype=float)
    tests = [
        ("F_ROH_marginal", [], ["F_ROH_z"]),
        ("load_marginal", [], ["Total_Genetic_Load_z"]),
        ("F_ROH_conditional_on_load", ["Total_Genetic_Load_z"], ["Total_Genetic_Load_z", "F_ROH_z"]),
        ("load_conditional_on_F_ROH", ["F_ROH_z"], ["F_ROH_z", "Total_Genetic_Load_z"]),
        ("F_ROH_and_load_joint", [], ["F_ROH_z", "Total_Genetic_Load_z"]),
    ]
    rng = np.random.default_rng(SEED)
    rows = []
    for name, reduced_terms, full_terms in tests:
        reduced, _ = design_matrix(frame, reduced_terms)
        full, _ = design_matrix(frame, full_terms)
        statistic, pvalue = freedman_lane_test(
            y, reduced, full, kinship, iterations, rng
        )
        rows.append(
            {
                "outcome": "CHS",
                "test": name,
                "added_df": full.shape[1] - reduced.shape[1],
                "partial_F_whitened": statistic,
                "permutations": iterations,
                "p_freedman_lane": pvalue,
                "variance_model": "reduced-model profile-ML LMM",
            }
        )
    return pd.DataFrame(rows)


def save_plot(
    frame: pd.DataFrame, predictions: pd.DataFrame, output: Path
) -> None:
    selected = predictions[predictions["model"].isin(["covariates", "covariates_both"])]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for model, color in [("covariates", "#777777"), ("covariates_both", "#0072B2")]:
        data = selected[selected["model"] == model]
        axes[0].scatter(
            data["CHS_observed"], data["CHS_prediction"], alpha=0.75, label=model, color=color
        )
    axes[0].axline((0, 0), slope=1, linestyle="--", color="black", linewidth=0.8)
    axes[0].set(xlabel="Observed CHS (morbidity composite)", ylabel="Grouped-CV predicted CHS", title="Out-of-cluster morbidity prediction")
    axes[0].legend(frameon=False)
    axes[1].scatter(
        frame["F_ROH_z"],
        frame["CHS"],
        c=frame["Total_Genetic_Load_z"],
        cmap="viridis",
        edgecolor="white",
        linewidth=0.4,
    )
    axes[1].set(xlabel="Standardized F_ROH", ylabel="CHS (morbidity composite)", title="Load (color), inbreeding, and morbidity")
    fig.tight_layout()
    fig.savefig(output / "fitness_proxy_validation.png", dpi=300)
    plt.close(fig)


def write_report(
    output: Path,
    frame: pd.DataFrame,
    source: dict,
    psd: dict,
    groups: np.ndarray,
    effects: pd.DataFrame,
    cv: pd.DataFrame,
    binary: pd.DataFrame,
    permutation: pd.DataFrame,
    age_summary: dict,
) -> None:
    both = effects[
        (effects["model"] == "covariates_both")
        & effects["term"].isin(["F_ROH_z", "Total_Genetic_Load_z"])
    ]
    lines = [
        "# Morbidity-proxy method validation",
        "",
        f"- Primary sample: {len(frame)} IID-aligned individuals.",
        f"- Relatedness matrix: `{source['source_type']}` from `{source['matrix_path']}`.",
        f"- PSD correction: {psd['eigenvalues_adjusted']} eigenvalues adjusted; raw minimum {psd['minimum_eigenvalue_raw']:.6g}.",
        f"- Relatedness clusters: {len(np.unique(groups))}, using GCTA relationship >= {RELATIONSHIP_THRESHOLD}.",
        f"- Age observed: {age_summary['observed']}/{age_summary['total']} ({100 * age_summary['fraction_observed']:.1f}%). Age was excluded from all primary models.",
        "",
        "## Primary conditional fixed effects",
        "",
    ]
    for row in both.itertuples():
        lines.append(
            f"- {row.term}: beta={row.beta:.4g} per SD, SE={row.SE:.4g}, Wald p={row.p_wald:.4g}."
        )
    best = cv.sort_values("RMSE").iloc[0]
    lines.extend(
        [
            "",
            "## Relatedness-cluster grouped cross-validation",
            "",
            f"- Lowest RMSE model: {best['model']} (RMSE={best['RMSE']:.4g}, R2={best['R2']:.4g}, Spearman rho={best['spearman_rho']:.4g}).",
            "",
            "## Whitened-residual permutation",
            "",
        ]
    )
    for row in permutation.itertuples():
        lines.append(
            f"- {row.test}: partial F={row.partial_F_whitened:.4g}, Freedman-Lane p={row.p_freedman_lane:.4g} ({row.permutations} permutations)."
        )
    best_binary = binary.sort_values("roc_auc", ascending=False).iloc[0]
    lines.extend(
        [
            "",
            "## Exploratory binary outcome",
            "",
            f"- Highest grouped-CV ROC AUC: {best_binary['model']} ({best_binary['roc_auc']:.3f}); PR AUC={best_binary['pr_auc']:.3f}.",
            "- Binary metrics reuse the continuous CHS LMM score. They are exploratory discrimination metrics, not a calibrated logistic mixed model.",
            "",
            "## Limitations",
            "",
            "- CHS is a morbidity proxy, not direct lifetime reproductive fitness or survival.",
            "- The cohort is small (68 animals), disease ascertainment may be incomplete, and cluster-grouped CV has high sampling variance.",
            "- Age is mostly missing and structurally concentrated among recorded cases; complete-case age adjustment would induce selection bias and is not used.",
            "- Sex is the only primary covariate. Unmeasured environment, observation effort, and shared husbandry/social conditions may confound associations.",
            "- Profile-ML variance components and Wald/LRT approximations can be unstable at this sample size; whitened Freedman-Lane results are the preferred robustness check.",
            "- Relatedness clusters are threshold-defined connected components and do not establish pedigrees.",
            "- Binary Brier and threshold metrics use clipped CHS scores and should not be interpreted as calibrated disease probabilities.",
        ]
    )
    (output / "validation_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--phenotype",
        type=Path,
        default=BASE / "output/phenotype_genotype_analysis/merged_phenotype_genotype.csv",
    )
    parser.add_argument(
        "--gcta-prefix",
        type=Path,
        default=BASE / "output/phenotype_genotype_analysis/gwas/grm",
    )
    parser.add_argument(
        "--gemma-kinship",
        type=Path,
        default=BASE / "output/phenotype_genotype_analysis/gwas/output/gemma_grm.cXX.txt",
    )
    parser.add_argument(
        "--gemma-fam",
        type=Path,
        default=BASE / "output/phenotype_genotype_analysis/gwas/tmp_gemma_base.fam",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=BASE / "output/method_validation/fitness",
    )
    parser.add_argument("--permutations", type=int, default=4999)
    parser.add_argument("--folds", type=int, default=5)
    args = parser.parse_args()
    if args.permutations < 99:
        raise ValueError("--permutations must be at least 99")
    args.output.mkdir(parents=True, exist_ok=True)

    raw_kinship, kinship_ids, source = choose_kinship(args)
    kinship, psd = regularize_psd(raw_kinship)
    frame, kinship, alignment = align_data(args.phenotype, kinship, kinship_ids)
    frame, scaling = prepare_data(frame)
    groups = connected_components(kinship, RELATIONSHIP_THRESHOLD)
    cluster_frame = pd.DataFrame(
        {
            "IID": frame["IID"],
            "relatedness_cluster": groups,
            "cluster_size": pd.Series(groups).map(pd.Series(groups).value_counts()).to_numpy(),
        }
    )

    effects, model_fits, fits = fit_primary_models(frame, kinship)
    comparisons = likelihood_comparisons(fits)
    cv, binary, predictions = grouped_cv(frame, kinship, groups, args.folds)
    permutation = permutation_tests(frame, kinship, args.permutations)

    age_observed = int(frame["Age"].notna().sum())
    age_summary = {
        "total": int(len(frame)),
        "observed": age_observed,
        "missing": int(len(frame) - age_observed),
        "fraction_observed": float(age_observed / len(frame)),
        "included_in_primary_models": False,
        "reason": "High and outcome-dependent missingness; recorded primarily for cases.",
    }
    metadata = {
        "seed": SEED,
        "phenotype_source": str(args.phenotype),
        "kinship": {**source, **psd},
        "n_aligned_complete": int(len(frame)),
        "scaling": scaling,
        "age": age_summary,
        "primary_outcome": "CHS",
        "primary_covariates": ["Sex_Male"],
        "load_proxy": "Total_Genetic_Load",
        "relationship_cluster_threshold": RELATIONSHIP_THRESHOLD,
        "n_relationship_clusters": int(len(np.unique(groups))),
        "cluster_sizes": sorted(
            [int(x) for x in pd.Series(groups).value_counts().tolist()], reverse=True
        ),
        "permutations": args.permutations,
        "limitations": [
            "CHS measures observed morbidity rather than direct reproductive fitness or survival.",
            "Small cohort and sparse outcomes limit precision and external validation.",
            "Age was excluded because missingness is extensive and outcome-dependent.",
            "Binary results are exploratory discrimination metrics from continuous-score models.",
            "Relatedness clusters depend on a pre-specified GRM threshold.",
        ],
    }

    alignment.to_csv(args.output / "iid_alignment.csv", index=False)
    cluster_frame.to_csv(args.output / "relatedness_clusters.csv", index=False)
    effects.to_csv(args.output / "fixed_effect_tests.csv", index=False)
    model_fits.to_csv(args.output / "model_fit_statistics.csv", index=False)
    comparisons.to_csv(args.output / "model_comparisons.csv", index=False)
    cv.to_csv(args.output / "grouped_cv_metrics.csv", index=False)
    predictions.to_csv(args.output / "grouped_cv_predictions.csv", index=False)
    binary.to_csv(args.output / "binary_exploratory_metrics.csv", index=False)
    permutation.to_csv(args.output / "freedman_lane_tests.csv", index=False)
    (args.output / "validation_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    save_plot(frame, predictions, args.output)
    write_report(
        args.output,
        frame,
        source,
        psd,
        groups,
        effects,
        cv,
        binary,
        permutation,
        age_summary,
    )
    print(f"Fitness-proxy validation complete: {args.output}")


if __name__ == "__main__":
    main()
