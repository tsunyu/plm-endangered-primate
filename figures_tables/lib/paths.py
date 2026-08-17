"""Analysis-root paths for publication figure production.

Resolution order for the analysis root (directory with data/ and output/):
1. Environment variable PLM_BASE_DIR
2. analysis/base_dir.env (written by configure_base_dir.sh)
3. Nearest ancestor that contains both data/ and output/
"""

from __future__ import annotations

import os
from pathlib import Path


def _read_env_file(path: Path) -> str:
    if not path.is_file():
        return ""
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped.startswith("export "):
            stripped = stripped[len("export ") :].strip()
        if stripped.startswith("PLM_BASE_DIR="):
            return stripped.split("=", 1)[1].strip().strip("'").strip('"')
    return ""


def _detect_repo() -> Path:
    env = os.environ.get("PLM_BASE_DIR", "").strip()
    if env:
        return Path(env).expanduser().resolve()

    pkg_root = Path(__file__).resolve().parents[2]
    from_file = _read_env_file(pkg_root / "analysis" / "base_dir.env")
    if from_file:
        return Path(from_file).expanduser().resolve()

    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "data").is_dir() and (parent / "output").is_dir():
            return parent

    raise SystemExit(
        "Could not find the analysis root (a directory containing data/ and output/).\n"
        "Run:  bash configure_base_dir.sh /path/to/analysis_root\n"
        "  or: export PLM_BASE_DIR=/path/to/analysis_root"
    )


REPO = _detect_repo()
WORKS = REPO / "works"
FIGURES = WORKS / "figures"
SUPP_FIGURES = FIGURES / "supplementary"
TABLES = WORKS / "tables"
RESULTS = WORKS / "results"

ANNOTATED_VCF = REPO / "output/phase2_annotation/snpeff_annotation/annotated_variants.vcf.gz"
ANCESTRAL_VCF = REPO / "output/phase2_annotation/ancestral_states/variants_with_ancestral.vcf.gz"
ESM2_PREDICTIONS = REPO / "output/phase4_plm_predictions/esm2/esm2_predictions.csv"
ENSEMBLE_PREDICTIONS = REPO / "output/phase4_plm_predictions/ensemble/ensemble_predictions.csv"
MERGED_PHENO = REPO / "output/phenotype_genotype_analysis/merged_phenotype_genotype.csv"
HET_SUMMARY = REPO / "output/phase3a_population_genomics/diversity_metrics/heterozygosity_summary.csv"
ROH_SUMMARY = REPO / "output/phase3a_population_genomics/roh_analysis/roh_summary_per_individual.csv"
PCA_RESULTS = REPO / "output/phase3a_population_genomics/population_structure/pca_results.csv"
NUC_DIV = REPO / "output/phase3a_population_genomics/diversity_metrics/nucleotide_diversity_summary.csv"
TAJIMAS_D = REPO / "output/phase3a_population_genomics/diversity_metrics/tajimas_d_summary.csv"
CORRELATION = REPO / "output/phenotype_genotype_analysis/correlation_results.csv"
CASE_CONTROL = REPO / "output/phenotype_genotype_analysis/case_control_analysis.csv"
MODEL_COMPARISON = REPO / "output/phase3b_fastsimcoal2/model_comparison/model_comparison.csv"
PARAM_ESTIMATES = REPO / "output/phase3b_fastsimcoal2/parameter_estimates.csv"
BOOTSTRAP_RESULTS = (
    REPO / "output/phase3b_fastsimcoal2/bootstrap/bottleneck_recent_contraction/bootstrap_results.csv"
)
FITNESS_FIXED = REPO / "output/method_validation/fitness/fixed_effect_tests.csv"
FITNESS_PERM = REPO / "output/method_validation/fitness/freedman_lane_tests.csv"
FITNESS_CV = REPO / "output/method_validation/fitness/grouped_cv_metrics.csv"
DFE_REPLICATES = REPO / "output/method_validation/dfe_simulation/scenario_replicate_summary.csv"
DFE_MEANS = REPO / "output/method_validation/dfe_simulation/scenario_means.csv"
INDIVIDUAL_LOAD = REPO / "output/phase5_genetic_load/individual_load/individual_genetic_load.csv"
DAF_CACHE = RESULTS / "daf_variant_metrics.csv"
DAF_SUMMARY = RESULTS / "daf_frequency_statistics.csv"

SIGMOID_K = 0.5287
SIGMOID_X0 = -6.8920
LOF_PATHOGENICITY = 0.95
BOOTSTRAP_SEED = 20260710
N_BOOTSTRAP = 2000

for directory in (FIGURES, SUPP_FIGURES, TABLES, RESULTS):
    directory.mkdir(parents=True, exist_ok=True)
