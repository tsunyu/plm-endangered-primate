#!/usr/bin/env bash
# Documentary orchestrator for the paper analysis steps.
# Default mode lists steps only and does not execute analyses.
#
# Usage:
#   bash pipeline.sh --list
#   bash pipeline.sh --run          # execute in paper order (long-running)
#   bash pipeline.sh --run-figures  # figures/tables only

set -euo pipefail

PKG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS="${PKG_DIR}/analysis"
FIGURES="${PKG_DIR}/figures_tables"

if [[ -z "${PLM_BASE_DIR:-}" && -f "${ANALYSIS}/base_dir.env" ]]; then
  set -a
  # shellcheck disable=SC1091
  source "${ANALYSIS}/base_dir.env"
  set +a
fi

MODE="list"

usage() {
  sed -n '1,12p' "$0" | sed 's/^# \{0,1\}//'
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --list) MODE="list"; shift ;;
    --run) MODE="run"; shift ;;
    --run-figures) MODE="run-figures"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

STEPS=(
  "1|Annotation|${ANALYSIS}/phase2_annotation/phase2_step0_annotation_pipeline.sh|bash"
  "2|Population genomics|${ANALYSIS}/phase3a_population_genomics/phase3a_population_genomics_pipeline.sh|bash"
  "3|Demography (fastsimcoal2)|${ANALYSIS}/phase3b_fastsimcoal2/run_complete_pipeline.sh|bash"
  "4|PLM / ESM-2 scoring|${ANALYSIS}/phase4_plm_predictions/phase4_step0_plm_pipeline.sh|bash"
  "5|Genetic load|${ANALYSIS}/phase5_genetic_load/phase5_step0_genetic_load_pipeline.sh|bash"
  "6|Phenotype / GWAS|${ANALYSIS}/phase8_phenotype_genotype_analysis/phase8_step0_genotype_phenotype_analysis.py|python"
  "7|Method validation|${ANALYSIS}/phase9_method_validation/run_all_validation.sh|bash"
  "8|Publication figures|${FIGURES}/run_all_figures.py|python"
  "9|Publication tables|${FIGURES}/regenerate_tables.py|python"
  "10|Supplementary Excel|${FIGURES}/export_supplementary_tables_xlsx.py|python"
)

list_steps() {
  echo "PLM endangered-primate analysis pipeline"
  if [[ -n "${PLM_BASE_DIR:-}" ]]; then
    echo "PLM_BASE_DIR is set"
  else
    echo "PLM_BASE_DIR: (unset — run configure_base_dir.sh or export it)"
  fi
  echo ""
  echo "Paper: Protein language modelling reveals a latent drift load"
  echo "       underlying health risk in an endangered primate"
  echo ""
  echo "Figure 1: manual artwork (no script)"
  echo ""
  printf '%-4s %-28s %s\n' "Step" "Module" "Entry point"
  printf '%-4s %-28s %s\n' "----" "----------------------------" "----------"
  for item in "${STEPS[@]}"; do
    IFS='|' read -r num name path runner <<<"${item}"
    rel="${path#${PKG_DIR}/}"
    printf '%-4s %-28s %s\n' "${num}" "${name}" "${rel}"
  done
}

run_one() {
  local path="$1"
  local runner="$2"
  if [[ ! -e "${path}" ]]; then
    echo "ERROR: missing ${path}" >&2
    exit 1
  fi
  echo ""
  echo "=== Running: ${path} ==="
  case "${runner}" in
    bash) bash "${path}" ;;
    python) python3 "${path}" ;;
    *) echo "Unknown runner: ${runner}" >&2; exit 1 ;;
  esac
}

case "${MODE}" in
  list)
    list_steps
    echo ""
    echo "Default is --list only (no execution). Use --run or --run-figures to execute."
    ;;
  run)
    list_steps
    for item in "${STEPS[@]}"; do
      IFS='|' read -r num name path runner <<<"${item}"
      run_one "${path}" "${runner}"
    done
    ;;
  run-figures)
    for item in "${STEPS[@]}"; do
      IFS='|' read -r num name path runner <<<"${item}"
      if [[ "${num}" -ge 8 ]]; then
        run_one "${path}" "${runner}"
      fi
    done
    ;;
esac
