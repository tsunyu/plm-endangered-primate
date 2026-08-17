#!/usr/bin/env bash
set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../load_base_dir.sh"
OUT_ROOT="${BASE_DIR}/output/method_validation"

mkdir -p "${OUT_ROOT}"

echo "=== ClinVar protein-held-out validation ==="
python3 "${SCRIPT_DIR}/validate_clinvar_calibration.py" \
    2>&1 | tee "${OUT_ROOT}/clinvar/clinvar_validation.log"

echo "=== ClinVar GroupKFold / protein-aware calibration ==="
python3 "${SCRIPT_DIR}/validate_clinvar_extended.py" \
    2>&1 | tee "${OUT_ROOT}/clinvar/clinvar_groupkfold.log"

echo "=== Cross-species score transfer ==="
python3 "${SCRIPT_DIR}/validate_cross_species_functional.py" \
    2>&1 | tee "${OUT_ROOT}/cross_species/cross_species.log"

echo "=== Population consistency (MAF / homozygote depletion) ==="
python3 "${SCRIPT_DIR}/validate_population_consistency.py" \
    2>&1 | tee "${OUT_ROOT}/population_consistency/population_consistency.log"

echo "=== Load sensitivity ==="
bash "${SCRIPT_DIR}/run_load_sensitivity.sh"

echo "=== Fitness-proxy validation (GRM-adjusted) ==="
python3 "${SCRIPT_DIR}/validate_fitness_proxy.py" \
    2>&1 | tee "${OUT_ROOT}/fitness/fitness_validation.log"

echo "=== DFE / demography forward simulations ==="
bash "${SCRIPT_DIR}/run_dfe_simulation.sh"

echo "All method validation steps completed."
