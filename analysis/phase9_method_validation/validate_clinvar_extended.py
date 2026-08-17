#!/usr/bin/env python3
"""Extended ClinVar calibration: GroupKFold and protein-aware calibration."""

from __future__ import annotations

import argparse
import json
import os
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
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import GroupKFold

from validate_clinvar_calibration import (
    FIXED_K,
    FIXED_X0,
    SEED,
    load_data,
    sigmoid_probability,
)


def global_calibration(train: pd.DataFrame, test: pd.DataFrame) -> np.ndarray:
    model = LogisticRegression(C=1e6, solver="lbfgs", max_iter=2000)
    model.fit(train[["esm2_score"]], train["label"])
    return model.predict_proba(test[["esm2_score"]])[:, 1]


def protein_z_calibration(train: pd.DataFrame, test: pd.DataFrame) -> np.ndarray:
    """Per-protein z-score of LLR within train, then global logistic on z."""
    stats = (
        train.groupby("protein_group")["esm2_score"]
        .agg(["mean", "std"])
        .rename(columns={"mean": "mu", "std": "sigma"})
    )
    stats["sigma"] = stats["sigma"].replace(0, np.nan)

    def transform(frame: pd.DataFrame) -> pd.DataFrame:
        merged = frame.merge(stats, left_on="protein_group", right_index=True, how="left")
        global_mu = frame["esm2_score"].mean()
        global_sigma = frame["esm2_score"].std(ddof=0) or 1.0
        merged["mu"] = merged["mu"].fillna(global_mu)
        merged["sigma"] = merged["sigma"].fillna(global_sigma)
        merged["z_score"] = (merged["esm2_score"] - merged["mu"]) / merged["sigma"]
        return merged

    train_z = transform(train)
    test_z = transform(test)
    model = LogisticRegression(C=1e6, solver="lbfgs", max_iter=2000)
    model.fit(train_z[["z_score"]], train_z["label"])
    return model.predict_proba(test_z[["z_score"]])[:, 1]


def calibration_slope(y: np.ndarray, p: np.ndarray) -> tuple[float, float]:
    p = np.clip(p, 1e-8, 1 - 1e-8)
    logit = np.log(p / (1 - p)).reshape(-1, 1)
    model = LogisticRegression(C=1e6, solver="lbfgs")
    model.fit(logit, y)
    return float(model.intercept_[0]), float(model.coef_[0, 0])


def evaluate_fold(
    train: pd.DataFrame, test: pd.DataFrame, name: str, probability: np.ndarray
) -> dict:
    y = test["label"].to_numpy()
    intercept, slope = calibration_slope(y, probability)
    return {
        "model": name,
        "n_test": int(len(test)),
        "roc_auc": float(roc_auc_score(y, probability)),
        "pr_auc": float(average_precision_score(y, probability)),
        "brier": float(brier_score_loss(y, probability)),
        "calibration_intercept": intercept,
        "calibration_slope": slope,
    }


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
    parser.add_argument("--folds", type=int, default=5)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    frame = load_data(args.input)
    eligible = (
        frame.groupby("protein_group")["label"]
        .agg(n_total="size", n_classes="nunique")
        .query("n_classes == 2 and n_total >= 5")
        .index
    )
    frame = frame[frame["protein_group"].isin(eligible)].reset_index(drop=True)

    gkf = GroupKFold(n_splits=args.folds)
    fold_rows: list[dict] = []
    for fold, (train_index, test_index) in enumerate(
        gkf.split(frame, frame["label"], groups=frame["protein_group"])
    ):
        train = frame.iloc[train_index]
        test = frame.iloc[test_index]
        fixed = sigmoid_probability(test["esm2_score"].to_numpy())
        global_prob = global_calibration(train, test)
        protein_prob = protein_z_calibration(train, test)
        for name, probability in [
            ("fixed_sigmoid", fixed),
            ("global_logistic", global_prob),
            ("protein_z_logistic", protein_prob),
        ]:
            row = evaluate_fold(train, test, name, probability)
            row["fold"] = fold
            fold_rows.append(row)

    fold_df = pd.DataFrame(fold_rows)
    fold_df.to_csv(args.output / "clinvar_groupkfold_metrics.csv", index=False)

    summary = (
        fold_df.groupby("model")
        .agg(
            roc_auc_mean=("roc_auc", "mean"),
            roc_auc_std=("roc_auc", "std"),
            pr_auc_mean=("pr_auc", "mean"),
            brier_mean=("brier", "mean"),
            calibration_slope_mean=("calibration_slope", "mean"),
        )
        .reset_index()
    )
    summary.to_csv(args.output / "clinvar_groupkfold_summary.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for model, group in fold_df.groupby("model"):
        axes[0].errorbar(
            [model],
            [group["roc_auc"].mean()],
            yerr=[group["roc_auc"].std()],
            fmt="o",
            capsize=4,
            label=model,
        )
    axes[0].set_ylabel("ROC-AUC")
    axes[0].set_title(f"{args.folds}-fold protein GroupKFold")
    axes[0].tick_params(axis="x", rotation=20)

    pivot = fold_df.pivot(index="fold", columns="model", values="brier")
    pivot.plot(kind="bar", ax=axes[1], rot=0)
    axes[1].set_ylabel("Brier score")
    axes[1].set_title("Calibration error by fold")
    fig.tight_layout()
    fig.savefig(args.output / "clinvar_groupkfold_comparison.png", dpi=300)
    plt.close(fig)

    metadata = {
        "n_variants": int(len(frame)),
        "n_proteins": int(frame["protein_group"].nunique()),
        "folds": args.folds,
        "fixed_parameters": {"k": FIXED_K, "x0": FIXED_X0},
        "conda_environment": os.environ.get("CONDA_DEFAULT_ENV"),
        "python": sys.version.split()[0],
        "scikit_learn": sklearn.__version__,
        "seed": SEED,
        "note": (
            "Protein-aware calibration z-scores LLR within training proteins; "
            "gene symbols and ClinVar review status remain unavailable."
        ),
    }
    (args.output / "clinvar_groupkfold_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
