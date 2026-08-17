#!/usr/bin/env bash
set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../load_base_dir.sh"
OUTPUT_DIR="${BASE_DIR}/output/method_validation/load_sensitivity"

mkdir -p "${OUTPUT_DIR}"

python3 "${SCRIPT_DIR}/run_load_sensitivity.py" \
    --output "${OUTPUT_DIR}" \
    "$@" 2>&1 | tee "${OUTPUT_DIR}/load_sensitivity.log"
