#!/bin/bash
################################################################################
# Run fastsimcoal2 Parameter Estimation
################################################################################
#
# This script runs fastsimcoal2 to estimate demographic parameters for
# multiple models using the site frequency spectrum.
#
# Input:  SFS file (SNJ_DAFpop0.obs) and model templates (.tpl, .est)
# Output: Maximum likelihood estimates for each model
#
# Author: Demographic Analysis Pipeline
# Date: 2026-01-26
################################################################################

set -e

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../load_base_dir.sh"
MODEL_DIR="${SCRIPT_DIR}/models"
SFS_DIR="${BASE_DIR}/output/phase3b_fastsimcoal2/sfs"
OUTPUT_DIR="${BASE_DIR}/output/phase3b_fastsimcoal2/models"
N_CORES=4
N_RUNS=50  # Number of independent optimization runs per model

# Log file
LOG_FILE="${OUTPUT_DIR}/fastsimcoal2_runs.log"
mkdir -p "${OUTPUT_DIR}"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "========================================================================"
echo "FASTSIMCOAL2 DEMOGRAPHIC PARAMETER ESTIMATION"
echo "========================================================================"
echo "Start time: $(date)"
echo ""

# Check if fastsimcoal2 is installed (prefer fsc28, fallback to fsc27/fsc26)
if ! command -v fsc28 &> /dev/null && ! command -v fsc27 &> /dev/null && ! command -v fsc26 &> /dev/null; then
    echo "ERROR: fastsimcoal2 not found in PATH"
    echo "Please install fastsimcoal2 (fsc28, fsc27 or fsc26)"
    echo "See: official fastsimcoal2 website for the latest version"
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

echo "Using fastsimcoal2: $(which $FSC_CMD)"
echo ""

extract_max_est_lhood() {
    local bestlhoods_file=$1
    awk '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                key = tolower($i)
                if (key == "maxestlhood" || key ~ /^maxestlhood/) {
                    idx = i
                }
            }
            next
        }
        NF > 0 { last = $0 }
        END {
            if (last == "") {
                exit 1
            }
            n = split(last, fields, /[[:space:]]+/)
            if (idx > 0 && idx <= n) {
                print fields[idx]
            } else if (n >= 2) {
                print fields[n - 1]
            } else {
                exit 1
            }
        }
    ' "${bestlhoods_file}"
}

# Check SFS file exists
SFS_FILE="${SFS_DIR}/SNJ_DAFpop0.obs"
if [ ! -f "${SFS_FILE}" ]; then
    echo "ERROR: SFS file not found: ${SFS_FILE}"
    echo "Please run phase3b_step0_prepare_sfs.py first"
    exit 1
fi

echo "SFS file: ${SFS_FILE}"
echo "Number of optimization runs per model: ${N_RUNS}"
echo "Cores per run: ${N_CORES}"
echo ""

# List of models to run
MODELS=(
    "constant_ne"
    "single_bottleneck"
    "two_consecutive_bottlenecks"
    "bottleneck_continuous_decline"
    "bottleneck_recent_contraction"
    "complex_multi_event"
)

echo "Models to estimate:"
for model in "${MODELS[@]}"; do
    echo "  - ${model}"
done
echo ""

################################################################################
# Function to run fastsimcoal2 for a single model
################################################################################
run_model() {
    local model_name=$1
    local run_id=$2
    
    echo "  Run ${run_id}/${N_RUNS}: ${model_name}"
    
    # Create run-specific directory
    local run_dir="${OUTPUT_DIR}/${model_name}/run_${run_id}"
    mkdir -p "${run_dir}"
    
    # Copy model files and SFS to run directory
    cp "${MODEL_DIR}/${model_name}.tpl" "${run_dir}/${model_name}.tpl"
    cp "${MODEL_DIR}/${model_name}.est" "${run_dir}/${model_name}.est"
    cp "${SFS_FILE}" "${run_dir}/SNJ_DAFpop0.obs"
    
    # fsc28 with -d option may expect either <prefix>_jointDAFpop0_0.obs
    # or <prefix>_DAFpop0.obs depending on version/build details.
    ln -sf SNJ_DAFpop0.obs "${run_dir}/${model_name}_jointDAFpop0_0.obs"
    ln -sf SNJ_DAFpop0.obs "${run_dir}/${model_name}_DAFpop0.obs"
    
    # Run fastsimcoal2
    cd "${run_dir}"
    
    # Run with different random seeds for each run
    SEED=$((RANDOM * run_id))
    
    ${FSC_CMD} -t ${model_name}.tpl \
               -e ${model_name}.est \
               -d \
               -0 \
               -n 100000 \
               -L 40 \
               -M \
               -c ${N_CORES} \
               -q \
               -s ${SEED} \
               --numBatches 1 \
               > run_${run_id}.log 2>&1
    
    # Extract MaxEstLhood. MaxObsLhood is the final column in fsc output and
    # should not be used for choosing the best optimization run.
    if [ -f "${model_name}/${model_name}.bestlhoods" ]; then
        MAX_LHOOD=$(extract_max_est_lhood "${model_name}/${model_name}.bestlhoods")
        echo "    MaxEstLhood: ${MAX_LHOOD}"
        echo "${run_id},${MAX_LHOOD}" >> "../run_likelihoods.csv"
    else
        echo "    WARNING: No .bestlhoods file generated"
        echo "${run_id},NA" >> "../run_likelihoods.csv"
    fi
    
    cd - > /dev/null
}

################################################################################
# Run all models
################################################################################

for model in "${MODELS[@]}"; do
    echo ""
    echo "========================================================================"
    echo "MODEL: ${model}"
    echo "========================================================================"
    
    # Create model output directory
    MODEL_OUT="${OUTPUT_DIR}/${model}"
    mkdir -p "${MODEL_OUT}"
    
    # Initialize likelihood tracking file
    echo "run_id,max_est_likelihood" > "${MODEL_OUT}/run_likelihoods.csv"
    
    # Check if model files exist
    if [ ! -f "${MODEL_DIR}/${model}.tpl" ] || [ ! -f "${MODEL_DIR}/${model}.est" ]; then
        echo "ERROR: Model files not found for ${model}"
        echo "Please run phase3b_step1_create_model_templates.py first"
        continue
    fi
    
    # Run multiple independent optimizations
    echo "Running ${N_RUNS} independent optimizations..."
    echo ""
    
    for run in $(seq 1 ${N_RUNS}); do
        run_model "${model}" "${run}"
    done
    
    echo ""
    echo "Completed ${N_RUNS} runs for ${model}"
    
    # Find best run (highest likelihood)
    if [ -f "${MODEL_OUT}/run_likelihoods.csv" ]; then
        # Sort by MaxEstLhood (column 2), take highest
        BEST_RUN=$(tail -n +2 "${MODEL_OUT}/run_likelihoods.csv" | \
                   grep -v "NA" | \
                   sort -t',' -k2 -rn | \
                   head -1 | \
                   cut -d',' -f1)
        
        if [ -n "${BEST_RUN}" ]; then
            echo "Best run: ${BEST_RUN}"
            
            # Copy best run results to model directory
            BEST_DIR="${MODEL_OUT}/run_${BEST_RUN}/${model}"
            if [ -d "${BEST_DIR}" ]; then
                cp -r "${BEST_DIR}" "${MODEL_OUT}/best_run"
                echo "Best run files copied to ${MODEL_OUT}/best_run"
                
                # Display best parameters
                if [ -f "${MODEL_OUT}/best_run/${model}.bestlhoods" ]; then
                    echo ""
                    echo "Best likelihood and parameters:"
                    tail -1 "${MODEL_OUT}/best_run/${model}.bestlhoods"
                fi
            fi
        else
            echo "WARNING: Could not determine best run (all runs may have failed)"
        fi
    fi
    
    echo ""
done

################################################################################
# Summary
################################################################################

echo ""
echo "========================================================================"
echo "PARAMETER ESTIMATION COMPLETE"
echo "========================================================================"
echo "End time: $(date)"
echo ""

echo "Output directory: ${OUTPUT_DIR}"
echo ""

echo "Results summary:"
for model in "${MODELS[@]}"; do
    if [ -f "${OUTPUT_DIR}/${model}/best_run/${model}.bestlhoods" ]; then
        MAX_L=$(extract_max_est_lhood "${OUTPUT_DIR}/${model}/best_run/${model}.bestlhoods")
        echo "  ${model}: MaxEstLhood = ${MAX_L}"
    else
        echo "  ${model}: No results"
    fi
done
echo ""

echo "Next steps:"
echo "  1. Run phase3b_step3_model_comparison.py to compare models"
echo "  2. Run phase3b_step4_bootstrap_ci.sh for confidence intervals"
echo "  3. Run phase3b_step5_analyze_results.py to extract parameters"
echo ""

exit 0
