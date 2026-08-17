#!/bin/bash
################################################################################
# Bootstrap CI (PARALLEL VERSION for 32-core CPU)
################################################################################
#
# Optimized parametric bootstrap using parallel execution
# Designed for 32-core CPU with 128GB RAM
#
# Author: Demographic Analysis Pipeline
# Date: 2026-02-12
################################################################################

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=phase3b_bootstrap_common.sh
source "${SCRIPT_DIR}/phase3b_bootstrap_common.sh"

# Configuration - OPTIMIZED FOR 32-CORE CPU
source "${SCRIPT_DIR}/../load_base_dir.sh"
MODEL_DIR="${BASE_DIR}/output/phase3b_fastsimcoal2/models"
OUTPUT_DIR="${BASE_DIR}/output/phase3b_fastsimcoal2/bootstrap"

N_CORES_PER_JOB=2      # Cores per fastsimcoal2 instance
MAX_PARALLEL_JOBS=16   # Maximum parallel jobs (32 cores / 2)
N_BOOTSTRAP=100        # Number of bootstrap replicates
BOOT_N_LOCI=200000     # Independent loci for parametric bootstrap simulation
BOOT_DNA_LENGTH=100    # bp per locus (fastsimcoal2 manual default)

# User can specify model name
BEST_MODEL="${1:-auto}"

# Log file
LOG_FILE="${OUTPUT_DIR}/bootstrap_parallel.log"
mkdir -p "${OUTPUT_DIR}"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "========================================================================"
echo "BOOTSTRAP CONFIDENCE INTERVALS (PARALLEL VERSION)"
echo "========================================================================"
echo "Start time: $(date)"
echo ""
echo "Hardware configuration:"
echo "  Cores per job: ${N_CORES_PER_JOB}"
echo "  Max parallel jobs: ${MAX_PARALLEL_JOBS}"
echo "  Bootstrap replicates: ${N_BOOTSTRAP}"
echo "  Expected speedup: ~8-10x vs serial version"
echo ""

# Check GNU parallel
if ! command -v parallel &> /dev/null; then
    echo "WARNING: GNU parallel not found. Performance will be suboptimal."
    echo "  Install with: sudo apt-get install parallel"
    USE_GNU_PARALLEL=false
else
    USE_GNU_PARALLEL=true
fi

# Check fastsimcoal2
if ! command -v fsc28 &> /dev/null && ! command -v fsc27 &> /dev/null && ! command -v fsc26 &> /dev/null; then
    echo "ERROR: fastsimcoal2 not found in PATH"
    exit 1
fi

if command -v fsc28 &> /dev/null; then
    FSC_CMD="fsc28"
elif command -v fsc27 &> /dev/null; then
    FSC_CMD="fsc27"
elif command -v fsc26 &> /dev/null; then
    FSC_CMD="fsc26"
fi

echo "Using: $(which $FSC_CMD)"
echo ""

################################################################################
# Determine best model
################################################################################

if [ "${BEST_MODEL}" = "auto" ]; then
    COMPARISON_FILE="${BASE_DIR}/output/phase3b_fastsimcoal2/model_comparison/model_comparison.csv"
    
    if [ ! -f "${COMPARISON_FILE}" ]; then
        echo "ERROR: Model comparison file not found"
        exit 1
    fi
    
    BEST_MODEL=$(tail -n +2 "${COMPARISON_FILE}" | head -1 | cut -d',' -f1)
    echo "Best model: ${BEST_MODEL}"
else
    echo "Using specified model: ${BEST_MODEL}"
fi

echo ""

# Check best run exists
BEST_RUN_DIR="${MODEL_DIR}/${BEST_MODEL}/best_run"
if [ ! -d "${BEST_RUN_DIR}" ]; then
    echo "ERROR: Best run directory not found: ${BEST_RUN_DIR}"
    exit 1
fi

################################################################################
# Setup parameter file
################################################################################

echo "Setting up parameter file for simulation..."

BOOT_DIR="${OUTPUT_DIR}/${BEST_MODEL}"
mkdir -p "${BOOT_DIR}"

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

# Find .par file (try multiple possible locations and names)
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
    else
        echo "ERROR: Could not find .par file for bootstrap"
        exit 1
    fi
fi

cp "${BEST_PAR}" "${BOOT_DIR}/${BEST_MODEL}_maxL.par"
prepare_bootstrap_sim_par \
    "${BOOT_DIR}/${BEST_MODEL}_maxL.par" \
    "${BOOT_DIR}/${BEST_MODEL}_boot.par" \
    "${BOOT_N_LOCI}" \
    "${BOOT_DNA_LENGTH}"
echo "✓ Bootstrap simulation parameter file ready"
echo ""

################################################################################
# Function to run single bootstrap replicate
################################################################################
run_bootstrap_replicate() {
    local rep=$1
    local boot_dir=$2
    local model_name=$3
    
    REP_DIR="${boot_dir}/rep_${rep}"
    mkdir -p "${REP_DIR}"
    cd "${REP_DIR}"
    
    local boot_prefix="${model_name}_boot"
    cp "${boot_dir}/${boot_prefix}.par" "./${boot_prefix}.par"
    cp "${TPL_FILE}" "./${model_name}.tpl"
    cp "${EST_FILE}" "./${model_name}.est"
    
    local init_args=()
    if [ -n "${INIT_VALUES_FILE}" ] && [ -f "${INIT_VALUES_FILE}" ]; then
        cp "${INIT_VALUES_FILE}" "./${model_name}.pv"
        init_args=(--initvalues "${model_name}.pv")
    fi
    
    # Simulate pseudo-observed SFS (official fsc28 parametric bootstrap workflow)
    ${FSC_CMD} -i "${boot_prefix}.par" -n 1 -j -d -s 0 -x -I -q > sim.log 2>&1
    
    SIM_SFS=$(find_simulated_sfs "${REP_DIR}" "${boot_prefix}" || true)
    if [ -n "${SIM_SFS}" ] && [ -f "${SIM_SFS}" ]; then
        cp "${SIM_SFS}" "./SNJ_DAFpop0.obs"
        ln -sf SNJ_DAFpop0.obs "./${model_name}_jointDAFpop0_0.obs"
        ln -sf SNJ_DAFpop0.obs "./${model_name}_DAFpop0.obs"
        
        ${FSC_CMD} -t ${model_name}.tpl \
                   -e ${model_name}.est \
                   -d \
                   -0 \
                   -n 100000 \
                   -L 40 \
                   -M \
                   -c ${N_CORES_PER_JOB} \
                   -q \
                   "${init_args[@]}" \
                   > est.log 2>&1
        
        if [ -f "${model_name}/${model_name}.bestlhoods" ]; then
            BOOT_PARAMS=$(tail -1 "${model_name}/${model_name}.bestlhoods")
            BOOT_LHOOD=$(extract_max_est_lhood "${model_name}/${model_name}.bestlhoods")
            echo "${rep},${BOOT_LHOOD},\"${BOOT_PARAMS}\"" > "${boot_dir}/rep_${rep}_bootstrap_result.csv"
            echo "[$(date +%H:%M:%S)] Bootstrap ${rep}: MaxEstLhood=${BOOT_LHOOD}"
        else
            echo "${rep},NA,NA" > "${boot_dir}/rep_${rep}_bootstrap_result.csv"
            echo "[$(date +%H:%M:%S)] Bootstrap ${rep}: Estimation failed"
        fi
    else
        echo "${rep},NA,NA" > "${boot_dir}/rep_${rep}_bootstrap_result.csv"
        echo "[$(date +%H:%M:%S)] Bootstrap ${rep}: Simulation failed"
    fi
    
    cd - > /dev/null
}

export -f run_bootstrap_replicate
export -f extract_max_est_lhood
export -f prepare_bootstrap_sim_par
export -f find_simulated_sfs
export BOOT_DIR TPL_FILE EST_FILE FSC_CMD N_CORES_PER_JOB INIT_VALUES_FILE

################################################################################
# Run bootstrap replicates in parallel
################################################################################

echo "========================================================================"
echo "RUNNING ${N_BOOTSTRAP} BOOTSTRAP REPLICATES IN PARALLEL"
echo "========================================================================"
echo ""

# Initialize results file
rm -f "${BOOT_DIR}"/rep_*_bootstrap_result.csv
echo "replicate,max_est_likelihood,parameters" > "${BOOT_DIR}/bootstrap_results.csv"

if [ "$USE_GNU_PARALLEL" = true ]; then
    # Use GNU parallel for optimal performance
    seq 1 ${N_BOOTSTRAP} | parallel -j ${MAX_PARALLEL_JOBS} \
        run_bootstrap_replicate {} ${BOOT_DIR} ${BEST_MODEL}
else
    # Fallback: run in batches with background jobs
    for batch_start in $(seq 1 ${MAX_PARALLEL_JOBS} ${N_BOOTSTRAP}); do
        batch_end=$((batch_start + MAX_PARALLEL_JOBS - 1))
        if [ ${batch_end} -gt ${N_BOOTSTRAP} ]; then
            batch_end=${N_BOOTSTRAP}
        fi
        
        echo "Batch: replicates ${batch_start}-${batch_end}"
        
        for rep in $(seq ${batch_start} ${batch_end}); do
            run_bootstrap_replicate ${rep} ${BOOT_DIR} ${BEST_MODEL} &
        done
        
        wait
    done
fi

# Merge per-replicate result files after parallel execution to avoid concurrent
# appends corrupting bootstrap_results.csv.
echo "replicate,max_est_likelihood,parameters" > "${BOOT_DIR}/bootstrap_results.csv"
for rep in $(seq 1 ${N_BOOTSTRAP}); do
    REP_RESULT="${BOOT_DIR}/rep_${rep}_bootstrap_result.csv"
    if [ -f "${REP_RESULT}" ]; then
        awk 'NF > 0 { print }' "${REP_RESULT}" >> "${BOOT_DIR}/bootstrap_results.csv"
    else
        echo "${rep},NA,NA" >> "${BOOT_DIR}/bootstrap_results.csv"
    fi
done

################################################################################
# Summary
################################################################################

echo ""
echo "========================================================================"
echo "BOOTSTRAP COMPLETE"
echo "========================================================================"
echo "End time: $(date)"
echo ""

BOOT_RESULTS="${BOOT_DIR}/bootstrap_results.csv"
if [ -f "${BOOT_RESULTS}" ]; then
    N_SUCCESS=$(tail -n +2 "${BOOT_RESULTS}" | grep -v "NA" | wc -l)
    echo "Successful bootstrap replicates: ${N_SUCCESS}/${N_BOOTSTRAP}"
    echo "Success rate: $((N_SUCCESS * 100 / N_BOOTSTRAP))%"
    echo ""
    echo "Bootstrap results: ${BOOT_RESULTS}"
fi

echo ""
echo "Next steps:"
echo "  1. Run phase3b_step5_analyze_results.py to calculate confidence intervals"
echo ""

exit 0
