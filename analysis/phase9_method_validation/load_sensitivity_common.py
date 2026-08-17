#!/usr/bin/env python3
"""Shared primitives for single-pass genotype load sensitivity analysis."""

from __future__ import annotations

import gzip
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


SIGMOID_K = 0.5287
SIGMOID_X0 = -6.8920


@dataclass
class VariantMatrix:
    variant_ids: np.ndarray
    kinds: np.ndarray
    probabilities: np.ndarray
    genotypes: np.ndarray
    called: np.ndarray
    maf: np.ndarray
    call_rate: np.ndarray

    def subset(self, mask: np.ndarray) -> "VariantMatrix":
        return VariantMatrix(
            self.variant_ids[mask],
            self.kinds[mask],
            self.probabilities[mask],
            self.genotypes[mask],
            self.called[mask],
            self.maf[mask],
            self.call_rate[mask],
        )


def sigmoid_probability(score: float) -> float:
    exponent = np.clip(SIGMOID_K * (score - SIGMOID_X0), -700, 700)
    return float(1.0 / (1.0 + np.exp(exponent)))


def load_predictions(path: Path) -> dict[str, float]:
    frame = pd.read_csv(path, usecols=["variant_id", "esm2_score"])
    frame["esm2_score"] = pd.to_numeric(frame["esm2_score"], errors="coerce")
    frame = frame.dropna(subset=["variant_id", "esm2_score"])
    return dict(zip(frame["variant_id"], frame["esm2_score"]))


def parse_consequences(info: str) -> tuple[bool, bool, bool]:
    effects: list[str] = []
    impacts: list[str] = []
    for item in info.split(";"):
        if not item.startswith("ANN="):
            continue
        for annotation in item[4:].split(","):
            fields = annotation.split("|")
            if len(fields) >= 3:
                effects.append(fields[1].lower())
                impacts.append(fields[2].upper())
    is_lof = "HIGH" in impacts
    is_missense = any("missense_variant" in effect for effect in effects)
    is_synonymous = any("synonymous_variant" in effect for effect in effects)
    return is_lof, is_missense, is_synonymous


def parse_genotypes(fields: list[str], n_samples: int) -> tuple[np.ndarray, np.ndarray]:
    dosage = np.zeros(n_samples, dtype=np.int8)
    called = np.zeros(n_samples, dtype=bool)
    for index, sample_field in enumerate(fields[:n_samples]):
        genotype = sample_field.split(":", 1)[0].replace("|", "/")
        alleles = genotype.split("/")
        if len(alleles) != 2 or "." in alleles:
            continue
        try:
            values = [int(allele) for allele in alleles]
        except ValueError:
            continue
        if any(value not in (0, 1) for value in values):
            continue
        called[index] = True
        dosage[index] = sum(values)
    return dosage, called


def _to_matrix(records: list[tuple], n_samples: int) -> VariantMatrix:
    if not records:
        return VariantMatrix(
            np.array([], dtype=str),
            np.array([], dtype=str),
            np.array([], dtype=float),
            np.empty((0, n_samples), dtype=np.int8),
            np.empty((0, n_samples), dtype=bool),
            np.array([], dtype=float),
            np.array([], dtype=float),
        )
    ids, kinds, probabilities, genotypes, called, mafs, call_rates = zip(*records)
    return VariantMatrix(
        np.asarray(ids),
        np.asarray(kinds),
        np.asarray(probabilities, dtype=float),
        np.stack(genotypes),
        np.stack(called),
        np.asarray(mafs, dtype=float),
        np.asarray(call_rates, dtype=float),
    )


def scan_vcf_once(
    vcf_path: Path,
    predictions: dict[str, float],
    lof_probability: float = 0.95,
) -> tuple[list[str], VariantMatrix, VariantMatrix, dict[str, int]]:
    """Scan a biallelic VCF once and retain burden and synonymous contributions."""
    opener = gzip.open if vcf_path.suffix == ".gz" else open
    samples: list[str] = []
    burden_records: list[tuple] = []
    synonymous_records: list[tuple] = []
    qc = {
        "vcf_records": 0,
        "multiallelic_skipped": 0,
        "scored_missense": 0,
        "lof": 0,
        "synonymous_controls": 0,
    }
    with opener(vcf_path, "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip().split("\t")[9:]
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 10 or not samples:
                continue
            qc["vcf_records"] += 1
            chrom, pos, _, ref, alt = fields[:5]
            if "," in alt:
                qc["multiallelic_skipped"] += 1
                continue
            variant_id = f"{chrom}:{pos}:{ref}:{alt}"
            is_lof, is_missense, is_synonymous = parse_consequences(fields[7])
            score = predictions.get(variant_id)
            kind = ""
            probability = np.nan
            if is_lof:
                kind = "lof"
                probability = lof_probability
            elif is_missense and score is not None:
                kind = "missense"
                probability = sigmoid_probability(score)
            elif is_synonymous and not is_missense:
                kind = "synonymous"
            else:
                continue

            dosage, called = parse_genotypes(fields[9:], len(samples))
            called_alleles = 2 * int(called.sum())
            if called_alleles == 0:
                continue
            alt_frequency = float(dosage.sum() / called_alleles)
            maf = min(alt_frequency, 1.0 - alt_frequency)
            call_rate = float(called.mean())
            record = (
                variant_id,
                kind,
                probability,
                dosage,
                called,
                maf,
                call_rate,
            )
            if kind == "synonymous":
                synonymous_records.append(record)
                qc["synonymous_controls"] += 1
            else:
                burden_records.append(record)
                qc["lof" if kind == "lof" else "scored_missense"] += 1
    return (
        samples,
        _to_matrix(burden_records, len(samples)),
        _to_matrix(synonymous_records, len(samples)),
        qc,
    )


def calculate_load(
    matrix: VariantMatrix,
    samples: list[str],
    *,
    scenario: str,
    h: float = 0.25,
    lof_weight: float = 0.95,
    mask: np.ndarray | None = None,
    weights: np.ndarray | None = None,
    alt_allele_count: bool = False,
) -> pd.DataFrame:
    if mask is None:
        mask = np.ones(len(matrix.variant_ids), dtype=bool)
    selected = matrix.subset(mask)
    if weights is None:
        weights = selected.probabilities.copy()
        weights[selected.kinds == "lof"] = lof_weight
    if len(weights) != len(selected.variant_ids):
        raise ValueError("Weight vector and selected variant count differ")

    if alt_allele_count:
        contribution = selected.genotypes.astype(float) * weights[:, None]
    else:
        expression = np.where(
            selected.genotypes == 2,
            1.0,
            np.where(selected.genotypes == 1, h, 0.0),
        )
        contribution = expression * weights[:, None]
    contribution[~selected.called] = 0.0
    raw = contribution.sum(axis=0)
    callable_loci = selected.called.sum(axis=0)
    per_1000 = np.divide(
        raw * 1000.0,
        callable_loci,
        out=np.full_like(raw, np.nan),
        where=callable_loci > 0,
    )
    sd = float(np.nanstd(per_1000, ddof=1))
    zscore = (
        (per_1000 - float(np.nanmean(per_1000))) / sd
        if sd > 0
        else np.zeros_like(per_1000)
    )
    return pd.DataFrame(
        {
            "IID": samples,
            "scenario": scenario,
            "raw_load": raw,
            "callable_loci": callable_loci,
            "load_per_1000_callable": per_1000,
            "z_load": zscore,
            "n_variants": len(selected.variant_ids),
        }
    )


def nearest_maf_callrate_match(
    target: VariantMatrix,
    controls: VariantMatrix,
    rng: np.random.Generator,
    maf_caliper: float = 0.02,
    call_rate_caliper: float = 0.05,
) -> tuple[np.ndarray, np.ndarray]:
    """Randomized coarsened-exact matching within calipers, without replacement."""
    maf_bin_width = 0.01
    call_bin_width = 0.02
    control_bins: dict[tuple[int, int], list[int]] = {}
    for index, (maf, call_rate) in enumerate(
        zip(controls.maf, controls.call_rate)
    ):
        key = (
            int(np.floor(maf / maf_bin_width)),
            int(np.floor(call_rate / call_bin_width)),
        )
        control_bins.setdefault(key, []).append(index)
    for indices in control_bins.values():
        rng.shuffle(indices)

    target_order = rng.permutation(len(target.variant_ids))
    matched_target: list[int] = []
    matched_control: list[int] = []
    for target_index in target_order:
        maf_bin = int(np.floor(target.maf[target_index] / maf_bin_width))
        call_bin = int(
            np.floor(target.call_rate[target_index] / call_bin_width)
        )
        keys = [
            (maf_bin + maf_offset, call_bin + call_offset)
            for maf_offset in range(-2, 3)
            for call_offset in range(-3, 4)
        ]
        rng.shuffle(keys)
        keys.sort(
            key=lambda key: abs(key[0] - maf_bin) + abs(key[1] - call_bin)
        )
        control_index = None
        for key in keys:
            candidates = control_bins.get(key, [])
            while candidates:
                candidate = candidates.pop()
                if (
                    abs(controls.maf[candidate] - target.maf[target_index])
                    <= maf_caliper
                    and abs(
                        controls.call_rate[candidate]
                        - target.call_rate[target_index]
                    )
                    <= call_rate_caliper
                ):
                    control_index = candidate
                    break
            if control_index is not None:
                break
        if control_index is None:
            continue
        matched_target.append(int(target_index))
        matched_control.append(control_index)
    order = np.argsort(matched_target)
    return np.asarray(matched_target)[order], np.asarray(matched_control)[order]


def association_rows(
    loads: pd.DataFrame, covariates: pd.DataFrame
) -> pd.DataFrame:
    from scipy.stats import pearsonr, spearmanr

    merged = loads.merge(covariates, on="IID", how="left")
    rows = []
    for scenario, group in merged.groupby("scenario", sort=False):
        for metric in ("F_ROH",):
            valid = group[["z_load", metric]].dropna()
            if len(valid) < 3 or valid[metric].nunique() < 2:
                continue
            pearson = pearsonr(valid["z_load"], valid[metric])
            spearman = spearmanr(valid["z_load"], valid[metric])
            rows.append(
                {
                    "scenario": scenario,
                    "external_metric": metric,
                    "n": len(valid),
                    "pearson_r": pearson.statistic,
                    "pearson_p": pearson.pvalue,
                    "spearman_rho": spearman.statistic,
                    "spearman_p": spearman.pvalue,
                }
            )
    return pd.DataFrame(rows)


def scenario_concordance(loads: pd.DataFrame, reference: str) -> pd.DataFrame:
    from scipy.stats import pearsonr, spearmanr

    pivot = loads.pivot(index="IID", columns="scenario", values="z_load")
    reference_values = pivot[reference]
    rows = []
    for scenario in pivot:
        valid = pd.concat([reference_values, pivot[scenario]], axis=1).dropna()
        valid.columns = ["reference", "scenario"]
        pearson = pearsonr(valid["reference"], valid["scenario"])
        spearman = spearmanr(valid["reference"], valid["scenario"])
        rows.append(
            {
                "reference": reference,
                "scenario": scenario,
                "n": len(valid),
                "pearson_r": pearson.statistic,
                "spearman_rho": spearman.statistic,
            }
        )
    return pd.DataFrame(rows)


def write_metadata(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


def summarize_matching(
    target: VariantMatrix, controls: VariantMatrix, target_indices: Iterable[int], control_indices: Iterable[int]
) -> dict[str, float | int]:
    target_indices = np.asarray(list(target_indices), dtype=int)
    control_indices = np.asarray(list(control_indices), dtype=int)
    return {
        "n_target": int(len(target.variant_ids)),
        "n_matched": int(len(target_indices)),
        "match_fraction": float(len(target_indices) / len(target.variant_ids)),
        "mean_absolute_maf_difference": float(
            np.mean(np.abs(target.maf[target_indices] - controls.maf[control_indices]))
        ),
        "max_absolute_maf_difference": float(
            np.max(np.abs(target.maf[target_indices] - controls.maf[control_indices]))
        ),
        "mean_absolute_call_rate_difference": float(
            np.mean(
                np.abs(
                    target.call_rate[target_indices]
                    - controls.call_rate[control_indices]
                )
            )
        ),
    }
