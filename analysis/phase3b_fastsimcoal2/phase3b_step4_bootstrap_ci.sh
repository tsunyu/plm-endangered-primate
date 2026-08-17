#!/bin/bash
################################################################################
# Bootstrap Confidence Intervals for fastsimcoal2 Parameters
################################################################################
#
# This script generates confidence intervals for demographic parameters using
# parametric bootstrap. It simulates SFS from the best-fit model and re-estimates
# parameters to quantify uncertainty.
#
# Input:  Best-fit model parameters
# Output: Bootstrap parameter distributions and 95% confidence intervals
#
# Author: Demographic Analysis Pipeline
# Date: 2026-01-26
################################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=phase3b_bootstrap_common.sh
source "${SCRIPT_DIR}/phase3b_bootstrap_common.sh"

# Configuration
source "${SCRIPT_DIR}/../load_base_dir.sh"
MODEL_DIR="${BASE_DIR}/output/phase3b_fastsimcoal2/models"
OUTPUT_DIR="${BASE_DIR}/output/phase3b_fastsimcoal2/bootstrap"
N_BOOTSTRAP=100  # Number of bootstrap replicates
N_CORES=4
BOOT_N_LOCI=200000
BOOT_DNA_LENGTH=100

# User can specify model name as argument, otherwise auto-detect best model
BEST_MODEL="${1:-auto}"

# Log file
LOG_FILE="${OUTPUT_DIR}/bootstrap.log"
mkdir -p "${OUTPUT_DIR}"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "========================================================================"
echo "BOOTSTRAP CONFIDENCE INTERVALS FOR DEMOGRAPHIC PARAMETERS"
echo "========================================================================"
echo "Start time: $(date)"
echo ""

# Check if fastsimcoal2 is installed (prefer fsc28, fallback to fsc27/fsc26)
if ! command -v fsc28 &> /dev/null && ! command -v fsc27 &> /dev/null && ! command -v fsc26 &> /dev/null; then
    echo "ERROR: fastsimcoal2 not found in PATH"
    exit 1
fi

# Determine fastsimcoal2 binary to use
if command -v fsc28 &> /dev/null; then
    FSC_CMD="fsc28"
elif command -v fsc27 &> /dev/null; then
    FSC_CMD="fsc27"
elif command -v fsc26 &> /dev/null; then
    FSC_CMD="fsc26"
fi

echo "Using: $(which $FSC_CMD)"
echo "Bootstrap replicates: ${N_BOOTSTRAP}"
echo ""

################################################################################
# Determine best model (if not specified)
################################################################################

if [ "${BEST_MODEL}" = "auto" ]; then
    echo "Auto-detecting best model from AIC comparison..."
    
    COMPARISON_FILE="${BASE_DIR}/output/phase3b_fastsimcoal2/model_comparison/model_comparison.csv"
    
    if [ ! -f "${COMPARISON_FILE}" ]; then
        echo "ERROR: Model comparison file not found: ${COMPARISON_FILE}"
        echo "Please run phase3b_step3_model_comparison.py first"
        exit 1
    fi
    
    # Extract best model (first line after header, first column)
    BEST_MODEL=$(tail -n +2 "${COMPARISON_FILE}" | head -1 | cut -d',' -f1)
    echo "Best model: ${BEST_MODEL}"
else
    echo "Using specified model: ${BEST_MODEL}"
fi

echo ""

# Check if best model results exist
BEST_RUN_DIR="${MODEL_DIR}/${BEST_MODEL}/best_run"
if [ ! -d "${BEST_RUN_DIR}" ]; then
    echo "ERROR: Best run directory not found: ${BEST_RUN_DIR}"
    echo "Please run phase3b_step2_run_fastsimcoal2.sh first"
    exit 1
fi

echo "Best model directory: ${BEST_RUN_DIR}"
echo ""

################################################################################
# Extract best-fit parameters
################################################################################

echo "Extracting best-fit parameters..."

PAR_FILE="${BEST_RUN_DIR}/${BEST_MODEL}/${BEST_MODEL}.bestlhoods"
if [ ! -f "${PAR_FILE}" ]; then
    PAR_FILE="${BEST_RUN_DIR}/${BEST_MODEL}.bestlhoods"
fi
if [ ! -f "${PAR_FILE}" ]; then
    echo "ERROR: Parameter file not found: ${PAR_FILE}"
    exit 1
fi

# Get parameter values from last line
PARAMS=$(tail -1 "${PAR_FILE}")
echo "Best parameters: ${PARAMS}"
echo ""

################################################################################
# Parametric Bootstrap
################################################################################

echo "========================================================================"
echo "PARAMETRIC BOOTSTRAP"
echo "========================================================================"
echo ""
echo "Strategy:"
echo "  1. Simulate SFS from best-fit model (${N_BOOTSTRAP} replicates)"
echo "  2. Re-estimate parameters for each simulated SFS"
echo "  3. Calculate 95% CI from bootstrap distribution"
echo ""

# Create bootstrap directory for this model
BOOT_DIR="${OUTPUT_DIR}/${BEST_MODEL}"
mkdir -p "${BOOT_DIR}"

# Copy model files
TPL_FILE="${SCRIPT_DIR}/models/${BEST_MODEL}.tpl"
EST_FILE="${SCRIPT_DIR}/models/${BEST_MODEL}.est"

if [ ! -f "${TPL_FILE}" ] || [ ! -f "${EST_FILE}" ]; then
    echo "ERROR: Model template files not found"
    exit 1
fi

INIT_VALUES_FILE="${BEST_RUN_DIR}/${BEST_MODEL}.pv"
if [ ! -f "${INIT_VALUES_FILE}" ]; then
    INIT_VALUES_FILE=$(find "${BEST_RUN_DIR}" -name "${BEST_MODEL}.pv" -type f | head -1)
fi
if [ -n "${INIT_VALUES_FILE}" ] && [ -f "${INIT_VALUES_FILE}" ]; then
    echo "Using initial values: ${INIT_VALUES_FILE}"
else
    INIT_VALUES_FILE=""
    echo "No .pv initial-values file found (bootstrap re-estimation will start from priors)"
fi

# Create parameter file for simulation (.par file from .bestlhoods)
# This requires converting the best parameters to .par format
echo "Setting up parameter file for simulation..."

# Copy the best .par file if it exists (try multiple possible locations and names)
BEST_PAR="${BEST_RUN_DIR}/${BEST_MODEL}/${BEST_MODEL}.par"
if [ ! -f "${BEST_PAR}" ]; then
    BEST_PAR="${BEST_RUN_DIR}/${BEST_MODEL}/${BEST_MODEL}_maxL.par"
fi
if [ ! -f "${BEST_PAR}" ]; then
    BEST_PAR="${BEST_RUN_DIR}/${BEST_MODEL}_maxL.par"
fi
if [ ! -f "${BEST_PAR}" ]; then
    BEST_PAR="${BEST_RUN_DIR}/${BEST_MODEL}.par"
fi
if [ ! -f "${BEST_PAR}" ]; then
    FOUND_PAR=$(find "${BEST_RUN_DIR}" -name "${BEST_MODEL}*.par" -type f | head -1)
    if [ -n "${FOUND_PAR}" ]; then
        BEST_PAR="${FOUND_PAR}"
    fi
fi

if [ -f "${BEST_PAR}" ]; then
    cp "${BEST_PAR}" "${BOOT_DIR}/${BEST_MODEL}_maxL.par"
    prepare_bootstrap_sim_par \
        "${BOOT_DIR}/${BEST_MODEL}_maxL.par" \
        "${BOOT_DIR}/${BEST_MODEL}_boot.par" \
        "${BOOT_N_LOCI}" \
        "${BOOT_DNA_LENGTH}"
    echo "✓ Bootstrap simulation parameter file ready"
else
    echo "ERROR: Could not find .par file for bootstrap simulation."
    echo "Searched in: ${BEST_RUN_DIR}"
    exit 1
fi

echo ""

################################################################################
# Run bootstrap replicates
################################################################################

echo "Running ${N_BOOTSTRAP} bootstrap replicates..."
echo ""

# Initialize results file
echo "replicate,max_est_likelihood,parameters" > "${BOOT_DIR}/bootstrap_results.csv"

for rep in $(seq 1 ${N_BOOTSTRAP}); do
    echo "Bootstrap replicate ${rep}/${N_BOOTSTRAP}"
    
    REP_DIR="${BOOT_DIR}/rep_${rep}"
    mkdir -p "${REP_DIR}"
    
    # Step 1: Simulate SFS from best-fit parameters
    cd "${REP_DIR}"
    
    # Copy parameter file and model files
    boot_prefix="${BEST_MODEL}_boot"
    cp "${BOOT_DIR}/${boot_prefix}.par" "./${boot_prefix}.par"
    cp "${TPL_FILE}" "./${BEST_MODEL}.tpl"
    cp "${EST_FILE}" "./${BEST_MODEL}.est"

    init_args=()
    if [ -n "${INIT_VALUES_FILE}" ] && [ -f "${INIT_VALUES_FILE}" ]; then
        cp "${INIT_VALUES_FILE}" "./${BEST_MODEL}.pv"
        init_args=(--initvalues "${BEST_MODEL}.pv")
    fi

    ${FSC_CMD} -i "${boot_prefix}.par" -n 1 -j -d -s 0 -x -I -q > sim.log 2>&1

    SIM_SFS=$(find_simulated_sfs "${REP_DIR}" "${boot_prefix}" || true)
    if [ -n "${SIM_SFS}" ] && [ -f "${SIM_SFS}" ]; then
        cp "${SIM_SFS}" "./SNJ_DAFpop0.obs"
        ln -sf SNJ_DAFpop0.obs "./${BEST_MODEL}_jointDAFpop0_0.obs"
        ln -sf SNJ_DAFpop0.obs "./${BEST_MODEL}_DAFpop0.obs"

        ${FSC_CMD} -t ${BEST_MODEL}.tpl \
                   -e ${BEST_MODEL}.est \
                   -d \
                   -0 \
                   -n 100000 \
                   -L 40 \
                   -M \
                   -c ${N_CORES} \
                   -q \
                   "${init_args[@]}" \
                   > est.log 2>&1

        if [ -f "${BEST_MODEL}/${BEST_MODEL}.bestlhoods" ]; then
            BOOT_PARAMS=$(tail -1 "${BEST_MODEL}/${BEST_MODEL}.bestlhoods")
            BOOT_LHOOD=$(extract_max_est_lhood "${BEST_MODEL}/${BEST_MODEL}.bestlhoods")
            echo "${rep},${BOOT_LHOOD},\"${BOOT_PARAMS}\"" >> "${BOOT_DIR}/bootstrap_results.csv"
            echo "  MaxEstLhood: ${BOOT_LHOOD}"
        else
            echo "  WARNING: Estimation failed"
            echo "${rep},NA,NA" >> "${BOOT_DIR}/bootstrap_results.csv"
        fi
    else
        echo "  WARNING: Simulation failed"
        echo "${rep},NA,NA" >> "${BOOT_DIR}/bootstrap_results.csv"
    fi
    
    cd - > /dev/null
    
    # Clean up to save space (optional)
    # rm -rf "${REP_DIR}"
done

echo ""
echo "Bootstrap replicates complete"
echo ""

################################################################################
# Analyze bootstrap results
################################################################################

echo "========================================================================"
echo "BOOTSTRAP RESULTS ANALYSIS"
echo "========================================================================"
echo ""

# This analysis will be done by the Python script in step 5
# For now, just report completion

BOOT_RESULTS="${BOOT_DIR}/bootstrap_results.csv"
if [ -f "${BOOT_RESULTS}" ]; then
    N_SUCCESS=$(tail -n +2 "${BOOT_RESULTS}" | grep -v "NA" | wc -l)
    echo "Successful bootstrap replicates: ${N_SUCCESS}/${N_BOOTSTRAP}"
    echo ""
    echo "Bootstrap results saved to: ${BOOT_RESULTS}"
else
    echo "WARNING: No bootstrap results file generated"
fi

echo ""
echo "========================================================================"
echo "BOOTSTRAP COMPLETE"
echo "========================================================================"
echo "End time: $(date)"
echo ""

echo "Output directory: ${BOOT_DIR}"
echo ""

echo "Next steps:"
echo "  1. Run phase3b_step5_analyze_results.py to calculate confidence intervals"
echo "  2. Review bootstrap parameter distributions"
echo ""

exit 0
