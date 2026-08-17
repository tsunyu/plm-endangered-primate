#!/usr/bin/env bash
set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../load_base_dir.sh"
OUT_ROOT="${BASE_DIR}/output/method_validation/dfe_simulation"

mkdir -p "${OUT_ROOT}"

echo "=== Wright-Fisher DFE / demography forward simulations ==="
python3 "${SCRIPT_DIR}/run_dfe_demography_simulation.py" \
    --output "${OUT_ROOT}" \
    --replicates 50 \
    --loci 500 \
    2>&1 | tee "${OUT_ROOT}/dfe_simulation.log"

echo "=== SLiM 5.1 cross-check (subset, optional) ==="
if python3 "${SCRIPT_DIR}/run_slim_dfe_crosscheck.py" \
    --output "${OUT_ROOT}" \
    --replicates 2 \
    --loci 50 \
    2>&1 | tee "${OUT_ROOT}/slim_crosscheck.log"; then
  echo "SLiM cross-check succeeded"
else
  echo "SLiM cross-check skipped or failed; primary inference uses Wright-Fisher forward simulations." | tee -a "${OUT_ROOT}/slim_crosscheck.log"
fi

echo "DFE simulation completed."
