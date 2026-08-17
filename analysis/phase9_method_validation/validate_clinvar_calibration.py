#!/usr/bin/env python3
"""Validate the ClinVar-to-ESM-2 sigmoid without variant-level leakage.

The available table does not contain ClinVar review status or gene symbols.
Protein accession is therefore used as the conservative grouping unit so that
no protein contributes variants to both training and test partitions.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn
from sklearn.calibration import calibration_curve
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    log_loss,
    matthews_corrcoef,
    roc_auc_score,
)
from sklearn.model_selection import GroupShuffleSplit, train_test_split


FIXED_K = 0.5287
FIXED_X0 = -6.8920
SEED = 20260710


def sigmoid_probability(score: np.ndarray, k: float = FIXED_K, x0: float = FIXED_X0) -> np.ndarray:
    exponent = np.clip(k * (score - x0), -700, 700)
    return 1.0 / (1.0 + np.exp(exponent))


def metric_row(name: str, y: np.ndarray, probability: np.ndarray) -> dict[str, float | str | int]:
    prediction = (probability >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, prediction, labels=[0, 1]).ravel()
    return {
        "model": name,
        "n": int(len(y)),
        "prevalence": float(y.mean()),
        "roc_auc": float(roc_auc_score(y, probability)),
        "pr_auc": float(average_precision_score(y, probability)),
        "mcc_at_0.5": float(matthews_corrcoef(y, prediction)),
        "balanced_accuracy_at_0.5": float(balanced_accuracy_score(y, prediction)),
        "sensitivity_at_0.5": float(tp / (tp + fn)) if tp + fn else np.nan,
        "specificity_at_0.5": float(tn / (tn + fp)) if tn + fp else np.nan,
        "brier": float(brier_score_loss(y, probability)),
        "log_loss": float(log_loss(y, np.clip(probability, 1e-12, 1 - 1e-12))),
    }


def calibration_fit(y: np.ndarray, probability: np.ndarray) -> tuple[float, float]:
    p = np.clip(probability, 1e-8, 1 - 1e-8)
    logit = np.log(p / (1 - p)).reshape(-1, 1)
    model = LogisticRegression(C=1e6, solver="lbfgs")
    model.fit(logit, y)
    return float(model.intercept_[0]), float(model.coef_[0, 0])


def protein_macro_auc(frame: pd.DataFrame, probability: np.ndarray) -> tuple[float, int]:
    work = frame.assign(probability=probability)
    aucs = []
    for _, group in work.groupby("protein_group"):
        if group["label"].nunique() == 2:
            aucs.append(roc_auc_score(group["label"], group["probability"]))
    return (float(np.mean(aucs)) if aucs else np.nan, len(aucs))


def clustered_bootstrap(
    frame: pd.DataFrame, probability: np.ndarray, iterations: int
) -> dict[str, tuple[float, float]]:
    rng = np.random.default_rng(SEED)
    labels = frame["label"].to_numpy()
    proteins = frame["protein_group"].to_numpy()
    unique_proteins = pd.unique(proteins)
    indices_by_protein = {
        protein: np.flatnonzero(proteins == protein) for protein in unique_proteins
    }
    values: dict[str, list[float]] = {"roc_auc": [], "pr_auc": [], "brier": []}
    for _ in range(iterations):
        sampled = rng.choice(unique_proteins, size=len(unique_proteins), replace=True)
        boot_index = np.concatenate([indices_by_protein[protein] for protein in sampled])
        boot_label = labels[boot_index]
        boot_probability = probability[boot_index]
        if np.unique(boot_label).size < 2:
            continue
        values["roc_auc"].append(roc_auc_score(boot_label, boot_probability))
        values["pr_auc"].append(average_precision_score(boot_label, boot_probability))
        values["brier"].append(brier_score_loss(boot_label, boot_probability))
    return {
        key: (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)))
        for key, vals in values.items()
        if vals
    }


def load_data(path: Path) -> pd.DataFrame:
    usecols = ["protein", "mutant", "DMS_bin_score", "esm2_score"]
    frame = pd.read_csv(path, usecols=usecols)
    frame = frame.dropna(subset=usecols).copy()
    frame = frame[frame["DMS_bin_score"].isin(["Pathogenic", "Benign"])].copy()
    frame["protein_group"] = frame["protein"].str.replace(r"\.\d+$", "", regex=True)
    frame["variant_key"] = frame["protein_group"] + ":" + frame["mutant"].astype(str)
    frame["label"] = (frame["DMS_bin_score"] == "Pathogenic").astype(int)
    frame["esm2_score"] = pd.to_numeric(frame["esm2_score"], errors="coerce")
    frame = frame.dropna(subset=["esm2_score"])

    conflicts = frame.groupby("variant_key")["label"].nunique()
    conflicting_keys = set(conflicts[conflicts > 1].index)
    frame = frame[~frame["variant_key"].isin(conflicting_keys)]
    frame = frame.drop_duplicates("variant_key", keep="first").reset_index(drop=True)
    frame.attrs["conflicting_variants_removed"] = len(conflicting_keys)
    return frame


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=get_base_dir() / "output/phase5_genetic_load/esm2_predictions.csv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=get_base_dir() / "output/method_validation/clinvar",
    )
    parser.add_argument("--bootstrap", type=int, default=1000)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    frame = load_data(args.input)
    splitter = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=SEED)
    train_index, test_index = next(
        splitter.split(frame, frame["label"], groups=frame["protein_group"])
    )
    train = frame.iloc[train_index].copy()
    test = frame.iloc[test_index].copy()

    model = LogisticRegression(C=1e6, solver="lbfgs", max_iter=2000)
    model.fit(train[["esm2_score"]], train["label"])
    refit_probability = model.predict_proba(test[["esm2_score"]])[:, 1]
    fixed_probability = sigmoid_probability(test["esm2_score"].to_numpy())

    # Naive split is retained only to quantify optimism from variant-level splitting.
    naive_train, naive_test = train_test_split(
        frame, test_size=0.2, random_state=SEED, stratify=frame["label"]
    )
    naive_model = LogisticRegression(C=1e6, solver="lbfgs", max_iter=2000)
    naive_model.fit(naive_train[["esm2_score"]], naive_train["label"])
    naive_probability = naive_model.predict_proba(naive_test[["esm2_score"]])[:, 1]

    metrics = [
        metric_row("fixed_sigmoid_protein_held_out", test["label"].to_numpy(), fixed_probability),
        metric_row("refit_logistic_protein_held_out", test["label"].to_numpy(), refit_probability),
        metric_row(
            "refit_logistic_naive_variant_split",
            naive_test["label"].to_numpy(),
            naive_probability,
        ),
    ]
    for row, data, probability in [
        (metrics[0], test, fixed_probability),
        (metrics[1], test, refit_probability),
        (metrics[2], naive_test, naive_probability),
    ]:
        intercept, slope = calibration_fit(data["label"].to_numpy(), probability)
        macro_auc, evaluable = protein_macro_auc(data, probability)
        row["calibration_intercept"] = intercept
        row["calibration_slope"] = slope
        row["protein_macro_auc"] = macro_auc
        row["proteins_with_both_labels"] = evaluable

    pd.DataFrame(metrics).to_csv(args.output / "clinvar_validation_metrics.csv", index=False)

    ci = clustered_bootstrap(test, fixed_probability, args.bootstrap)
    fitted_k = float(-model.coef_[0, 0])
    fitted_x0 = float(-model.intercept_[0] / model.coef_[0, 0])
    metadata = {
        "source": str(args.input),
        "runtime": {
            "conda_environment": os.environ.get("CONDA_DEFAULT_ENV"),
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scikit_learn": sklearn.__version__,
            "random_seed": SEED,
        },
        "n_after_qc": int(len(frame)),
        "n_proteins": int(frame["protein_group"].nunique()),
        "n_train": int(len(train)),
        "n_test": int(len(test)),
        "train_proteins": int(train["protein_group"].nunique()),
        "test_proteins": int(test["protein_group"].nunique()),
        "protein_overlap": int(
            len(set(train["protein_group"]).intersection(set(test["protein_group"])))
        ),
        "conflicting_variants_removed": int(frame.attrs["conflicting_variants_removed"]),
        "fixed_parameters": {"k": FIXED_K, "x0": FIXED_X0},
        "refit_parameters": {
            "sklearn_coef_for_score": float(model.coef_[0, 0]),
            "sklearn_intercept": float(model.intercept_[0]),
            "equivalent_k": fitted_k,
            "equivalent_x0": fitted_x0,
        },
        "fixed_sigmoid_cluster_bootstrap_95ci": ci,
        "known_limitations": [
            "No ClinVar review status is present in the supplied table.",
            "No gene symbol is present; RefSeq protein accession is the held-out grouping unit.",
            "ClinVar discrimination validates clinical-label prediction, not a selection coefficient.",
        ],
    }
    (args.output / "clinvar_validation_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].hist(
        frame.loc[frame["label"] == 0, "esm2_score"],
        bins=60,
        density=True,
        alpha=0.6,
        label="Benign",
    )
    axes[0].hist(
        frame.loc[frame["label"] == 1, "esm2_score"],
        bins=60,
        density=True,
        alpha=0.6,
        label="Pathogenic",
    )
    axes[0].set(xlabel="ESM-2 LLR", ylabel="Density", title="ClinVar score distributions")
    axes[0].legend()

    for name, probability in [
        ("Fixed sigmoid", fixed_probability),
        ("Refit logistic", refit_probability),
    ]:
        observed, predicted = calibration_curve(
            test["label"], probability, n_bins=10, strategy="quantile"
        )
        axes[1].plot(predicted, observed, marker="o", label=name)
    axes[1].plot([0, 1], [0, 1], linestyle="--", color="black", label="Perfect calibration")
    axes[1].set(
        xlabel="Mean predicted pathogenicity probability",
        ylabel="Observed pathogenic fraction",
        title="Protein-held-out calibration",
    )
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(args.output / "clinvar_validation.png", dpi=300)
    plt.close(fig)

    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
