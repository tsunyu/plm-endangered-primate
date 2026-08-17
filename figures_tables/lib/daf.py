"""Polarized derived-allele-frequency utilities."""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .paths import (
    ANNOTATED_VCF,
    DAF_CACHE,
    DAF_SUMMARY,
    ESM2_PREDICTIONS,
    LOF_PATHOGENICITY,
    RESULTS,
    SIGMOID_K,
    SIGMOID_X0,
)

DELETERIOUS_MISSENSE_QUANTILE = 0.90

LOF_EFFECTS = {
    "stop_gained",
    "stop_lost",
    "start_lost",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
}

RISK_BINS = [0.0, 0.2, 0.4, 0.6, 0.8, 1.01]
RISK_LABELS = ["<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", ">0.8"]
DAF_CLASS_BINS = [0, 0.05, 0.20, 0.50, 0.95, 1.000001]
DAF_CLASS_LABELS = ["<0.05", "0.05-0.20", "0.20-0.50", "0.50-0.95", ">=0.95"]


def sigmoid_probability(score: float) -> float:
    exponent = np.clip(SIGMOID_K * (score - SIGMOID_X0), -700, 700)
    return float(1.0 / (1.0 + np.exp(exponent)))


def parse_ancestral(info: str) -> str | None:
    for field in info.split(";"):
        if field.startswith("AA="):
            allele = field.split("=", 1)[1]
            if allele not in {".", "N", "-"}:
                return allele
    return None


def parse_effect(info: str) -> tuple[str | None, str | None]:
    for field in info.split(";"):
        if not field.startswith("ANN="):
            continue
        parts = field[4:].split("|")
        if len(parts) >= 3:
            return parts[1].lower(), parts[2].upper()
    return None, None


def parse_genotypes(fields: list[str]) -> tuple[int, int, int, int, int]:
    """Return ref_count, alt_count, derived_hom, called_samples, n_called_alleles."""
    ref_count = alt_count = derived_hom = called_samples = 0
    for sample in fields:
        gt = sample.split(":", 1)[0].replace("|", "/")
        alleles = gt.split("/")
        if len(alleles) != 2 or "." in alleles:
            continue
        try:
            a1, a2 = (int(a) for a in alleles)
        except ValueError:
            continue
        called_samples += 1
        ref_count += (2 - a1 - a2)
        alt_count += (a1 + a2)
    n_called = ref_count + alt_count
    return ref_count, alt_count, derived_hom, called_samples, n_called


def polarize_counts(
    ref: str, alt: str, ancestral: str, ref_count: int, alt_count: int, genotypes: list[str]
) -> tuple[str | None, str | None, int, int, int, str]:
    if ancestral == ref:
        derived, risk = alt, "alt_derived"
    elif ancestral == alt:
        derived, risk = ref, "ref_derived"
    else:
        return None, None, 0, 0, 0, "ambiguous"

    derived_count = alt_count if derived == alt else ref_count
    derived_hom = 0
    derived_index = 1 if derived == alt else 0
    for sample in genotypes:
        gt = sample.split(":", 1)[0].replace("|", "/")
        alleles = gt.split("/")
        if len(alleles) != 2 or "." in alleles:
            continue
        try:
            a1, a2 = (int(a) for a in alleles)
        except ValueError:
            continue
        if a1 == derived_index and a2 == derived_index:
            derived_hom += 1
    return ancestral, derived, derived_count, derived_hom, "polarized"


def load_esm_lookup() -> dict[str, float]:
    esm = pd.read_csv(ESM2_PREDICTIONS, usecols=["variant_id", "esm2_score"])
    return {
        row.variant_id: sigmoid_probability(row.esm2_score)
        for row in esm.itertuples(index=False)
        if pd.notna(row.esm2_score)
    }


def deleterious_missense_threshold(
    frame: pd.DataFrame | None = None, quantile: float = DELETERIOUS_MISSENSE_QUANTILE
) -> float:
    """Top-decile threshold on ClinVar-calibrated P across all scored missense variants."""
    if frame is not None:
        probs = frame.loc[
            (frame["kind"] == "missense") & frame["pathogenicity_prob"].notna(),
            "pathogenicity_prob",
        ]
        if not probs.empty:
            return float(probs.quantile(quantile))
    lookup = load_esm_lookup()
    if not lookup:
        return float("nan")
    probs = pd.Series(list(lookup.values()))
    return float(probs.quantile(quantile))


def annotate_deleterious_missense(
    frame: pd.DataFrame, quantile: float = DELETERIOUS_MISSENSE_QUANTILE
) -> tuple[pd.DataFrame, float]:
    esm_lookup = load_esm_lookup()
    threshold = deleterious_missense_threshold(quantile=quantile)
    frame = frame.copy()
    missense_mask = frame["kind"] == "missense"
    frame.loc[missense_mask, "pathogenicity_prob"] = frame.loc[missense_mask, "variant_id"].map(
        esm_lookup
    )
    frame.loc[frame["kind"] == "lof", "pathogenicity_prob"] = LOF_PATHOGENICITY
    frame.loc[frame["kind"] == "synonymous", "pathogenicity_prob"] = 0.0
    frame["is_deleterious_missense"] = (
        missense_mask
        & frame["pathogenicity_prob"].notna()
        & (frame["pathogenicity_prob"] >= threshold)
    )
    return frame, threshold


def scan_variants(force: bool = False) -> pd.DataFrame:
    if DAF_CACHE.exists() and not force:
        frame = pd.read_csv(DAF_CACHE)
        frame, _ = annotate_deleterious_missense(frame)
        return frame

    esm_lookup = load_esm_lookup()
    rows: list[dict] = []

    with gzip.open(ANNOTATED_VCF, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, _, ref, alt = fields[:5]
            if "," in alt:
                continue
            info = fields[7]
            genotypes = fields[9:]
            ancestral = parse_ancestral(info)
            effect, impact = parse_effect(info)
            if effect is None:
                continue

            variant_id = f"{chrom}:{pos}:{ref}:{alt}"
            ref_count, alt_count, _, _, n_called = parse_genotypes(genotypes)
            if n_called == 0:
                continue

            status = "missing_ancestral"
            daf = np.nan
            derived_hom = 0
            polarization = "unresolved"
            if ancestral is not None:
                _, _, derived_count, derived_hom, polarization = polarize_counts(
                    ref, alt, ancestral, ref_count, alt_count, genotypes
                )
                if polarization == "polarized":
                    daf = derived_count / n_called
                    status = "polarized"
                else:
                    status = "ambiguous"

            is_lof = effect in LOF_EFFECTS or impact == "HIGH"
            is_missense = "missense_variant" in effect
            is_synonymous = "synonymous_variant" in effect
            if not (is_lof or is_missense or is_synonymous):
                continue
            probability = np.nan
            kind = "other"
            if is_lof:
                kind = "lof"
                probability = LOF_PATHOGENICITY
            elif is_missense and variant_id in esm_lookup:
                kind = "missense"
                probability = esm_lookup[variant_id]
            elif is_synonymous:
                kind = "synonymous"
                probability = 0.0

            rows.append(
                {
                    "variant_id": variant_id,
                    "chrom": chrom,
                    "pos": int(pos),
                    "ref": ref,
                    "alt": alt,
                    "ancestral": ancestral,
                    "kind": kind,
                    "is_deleterious_missense": False,
                    "pathogenicity_prob": probability,
                    "daf": daf,
                    "derived_hom_obs": derived_hom,
                    "n_called_samples": int(n_called / 2),
                    "n_called_alleles": n_called,
                    "polarization_status": status,
                }
            )

    frame = pd.DataFrame(rows)
    frame, threshold = annotate_deleterious_missense(frame)
    metadata = pd.DataFrame(
        [
            {
                "deleterious_missense_quantile": DELETERIOUS_MISSENSE_QUANTILE,
                "deleterious_missense_threshold": threshold,
                "n_scored_missense": int(
                    ((frame["kind"] == "missense") & frame["pathogenicity_prob"].notna()).sum()
                ),
                "n_deleterious_missense": int(frame["is_deleterious_missense"].sum()),
            }
        ]
    )
    DAF_CACHE.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(DAF_CACHE, index=False)
    metadata.to_csv(RESULTS / "figure5_deleterious_missense_threshold.csv", index=False)
    return frame


def polarized_subset(frame: pd.DataFrame) -> pd.DataFrame:
    return frame[frame["polarization_status"] == "polarized"].copy()


def assign_risk_bin(probability: float) -> str:
    for low, high, label in zip(RISK_BINS[:-1], RISK_BINS[1:], RISK_LABELS):
        if low <= probability < high:
            return label
    return RISK_LABELS[-1]


def daf_class_counts(daf_values: pd.Series) -> pd.DataFrame:
    counts = pd.cut(
        daf_values, bins=DAF_CLASS_BINS, labels=DAF_CLASS_LABELS, include_lowest=True
    ).value_counts()
    total = len(daf_values)
    return pd.DataFrame(
        {
            "class": DAF_CLASS_LABELS,
            "count": [int(counts.get(label, 0)) for label in DAF_CLASS_LABELS],
            "fraction": [counts.get(label, 0) / total if total else 0 for label in DAF_CLASS_LABELS],
        }
    )


def homozygote_depletion_by_bin(frame: pd.DataFrame) -> pd.DataFrame:
    missense = frame[(frame["kind"] == "missense") & frame["pathogenicity_prob"].notna()].copy()
    missense["risk_bin"] = missense["pathogenicity_prob"].map(assign_risk_bin)
    rows = []
    for label in RISK_LABELS:
        subset = missense[missense["risk_bin"] == label]
        if subset.empty:
            continue
        expected = 0.0
        observed = 0
        callable_total = 0
        for row in subset.itertuples():
            p = row.daf
            n = row.n_called_samples
            expected += n * (p**2)
            observed += row.derived_hom_obs
            callable_total += n
        ratio = observed / expected if expected > 0 else np.nan
        rows.append(
            {
                "risk_bin": label,
                "n_variants": len(subset),
                "callable_samples": callable_total,
                "observed_hom": observed,
                "expected_hom": expected,
                "oe_ratio": ratio,
            }
        )
    return pd.DataFrame(rows)


def summarize_daf(frame: pd.DataFrame) -> pd.DataFrame:
    polarized = polarized_subset(frame)
    summaries = []
    for category, subset in [
        ("lof", polarized[polarized["kind"] == "lof"]),
        (
            "deleterious_missense",
            polarized[
                (polarized["kind"] == "missense") & polarized["is_deleterious_missense"]
            ],
        ),
        (
            "scored_missense",
            polarized[(polarized["kind"] == "missense") & polarized["pathogenicity_prob"].notna()],
        ),
    ]:
        if subset.empty:
            continue
        class_df = daf_class_counts(subset["daf"])
        class_df["category"] = category
        class_df["n_variants"] = len(subset)
        summaries.append(class_df)
    summary = pd.concat(summaries, ignore_index=True)
    DAF_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(DAF_SUMMARY, index=False)
    return summary


def probability_daf_trend(frame: pd.DataFrame) -> dict:
    missense = polarized_subset(frame)
    missense = missense[
        (missense["kind"] == "missense") & missense["pathogenicity_prob"].notna()
    ]
    rho, p_value = spearmanr(missense["pathogenicity_prob"], missense["daf"])
    return {"n": len(missense), "spearman_rho": float(rho), "p_value": float(p_value)}
