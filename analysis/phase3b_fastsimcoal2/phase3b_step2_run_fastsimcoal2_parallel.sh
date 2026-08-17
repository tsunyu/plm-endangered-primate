#!/bin/bash
################################################################################
# Run fastsimcoal2 Parameter Estimation (PARALLEL VERSION for 32-core CPU)
################################################################################
#
# This optimized version runs multiple models and runs in parallel using GNU parallel
# Designed for 32-core CPU with 128GB RAM
#
# Optimization strategy:
#   - Run multiple models simultaneously
#   - Run multiple optimization runs per model in parallel
#   - Each fastsimcoal2 instance uses 2 cores (allows 16 parallel jobs)
#
# Author: Demographic Analysis Pipeline
# Date: 2026-02-12
################################################################################

set -e

# Configuration - OPTIMIZED FOR 32-CORE CPU
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../load_base_dir.sh"
MODEL_DIR="${SCRIPT_DIR}/models"
SFS_DIR="${BASE_DIR}/output/phase3b_fastsimcoal2/sfs"
OUTPUT_DIR="${BASE_DIR}/output/phase3b_fastsimcoal2/models"

# Parallel configuration
N_CORES_PER_JOB=2      # Cores per fastsimcoal2 instance
MAX_PARALLEL_JOBS=16   # Maximum parallel jobs (32 cores / 2 cores per job)
N_RUNS=50              # Number of independent optimization runs per model

# Log file
LOG_FILE="${OUTPUT_DIR}/fastsimcoal2_runs_parallel.log"
mkdir -p "${OUTPUT_DIR}"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "========================================================================"
echo "FASTSIMCOAL2 PARAMETER ESTIMATION (PARALLEL VERSION)"
echo "========================================================================"
echo "Start time: $(date)"
echo ""
echo "Hardware configuration:"
echo "  Total CPU cores: 32"
echo "  Cores per job: ${N_CORES_PER_JOB}"
echo "  Max parallel jobs: ${MAX_PARALLEL_JOBS}"
echo "  Expected speedup: ~8-10x vs serial version"
echo ""

# Check if GNU parallel is installed
if ! command -v parallel &> /dev/null; then
    echo "WARNING: GNU parallel not found. Installing is recommended for best performance."
    echo "  Install with: sudo apt-get install parallel"
    echo ""
    echo "Falling back to basic parallel execution with background jobs..."
    USE_GNU_PARALLEL=false
else
    echo "Using GNU parallel: $(which parallel)"
    USE_GNU_PARALLEL=true
fi
echo ""

# Check if fastsimcoal2 is installed
if ! command -v fsc28 &> /dev/null && ! command -v fsc27 &> /dev/null && ! command -v fsc26 &> /dev/null; then
    echo "ERROR: fastsimcoal2 not found in PATH"
    exit 1
fi

# Determine fastsimcoal2 binary
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

# Check SFS file
SFS_FILE="${SFS_DIR}/SNJ_DAFpop0.obs"
if [ ! -f "${SFS_FILE}" ]; then
    echo "ERROR: SFS file not found: ${SFS_FILE}"
    exit 1
fi

echo "SFS file: ${SFS_FILE}"
echo "Number of optimization runs per model: ${N_RUNS}"
echo ""

# List of models
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
# Function to run a single fastsimcoal2 optimization
################################################################################
run_single_optimization() {
    local model_name=$1
    local run_id=$2
    local output_dir=$3
    
    # Create run-specific directory
    local run_dir="${output_dir}/${model_name}/run_${run_id}"
    mkdir -p "${run_dir}"
    
    # Copy model files and SFS
    cp "${MODEL_DIR}/${model_name}.tpl" "${run_dir}/${model_name}.tpl"
    cp "${MODEL_DIR}/${model_name}.est" "${run_dir}/${model_name}.est"
    cp "${SFS_FILE}" "${run_dir}/SNJ_DAFpop0.obs"
    
    # fsc28 compatibility: create expected filenames
    ln -sf SNJ_DAFpop0.obs "${run_dir}/${model_name}_jointDAFpop0_0.obs"
    ln -sf SNJ_DAFpop0.obs "${run_dir}/${model_name}_DAFpop0.obs"
    
    # Run fastsimcoal2
    cd "${run_dir}"
    
    # Unique seed for each run
    SEED=$((12345 + run_id * 1000))
    
    ${FSC_CMD} -t ${model_name}.tpl \
               -e ${model_name}.est \
               -d \
               -0 \
               -n 100000 \
               -L 40 \
               -M \
               -c ${N_CORES_PER_JOB} \
               -q \
               -s ${SEED} \
               --numBatches 1 \
               > run_${run_id}.log 2>&1
    
    # Extract MaxEstLhood. MaxObsLhood is the final column in fsc output and
    # is constant across runs for a fixed observed SFS.
    if [ -f "${model_name}/${model_name}.bestlhoods" ]; then
        MAX_LHOOD=$(extract_max_est_lhood "${model_name}/${model_name}.bestlhoods")
        echo "${run_id},${MAX_LHOOD}" > "../run_${run_id}_likelihood.csv"
        echo "[$(date +%H:%M:%S)] Model ${model_name} run ${run_id}: MaxEstLhood=${MAX_LHOOD}"
    else
        echo "${run_id},NA" > "../run_${run_id}_likelihood.csv"
        echo "[$(date +%H:%M:%S)] Model ${model_name} run ${run_id}: FAILED"
    fi
    
    cd - > /dev/null
}

export -f run_single_optimization
export -f extract_max_est_lhood
export FSC_CMD MODEL_DIR SFS_FILE OUTPUT_DIR N_CORES_PER_JOB

################################################################################
# Run all models in parallel
################################################################################

if [ "$USE_GNU_PARALLEL" = true ]; then
    echo "========================================================================"
    echo "RUNNING WITH GNU PARALLEL (OPTIMAL PERFORMANCE)"
    echo "========================================================================"
    echo ""
    
    # Create job list: model_name run_id
    JOB_LIST="${OUTPUT_DIR}/job_list.txt"
    > "${JOB_LIST}"
    
    for model in "${MODELS[@]}"; do
        # Initialize model directory
        MODEL_OUT="${OUTPUT_DIR}/${model}"
        mkdir -p "${MODEL_OUT}"
        rm -f "${MODEL_OUT}"/run_*_likelihood.csv
        echo "run_id,max_est_likelihood" > "${MODEL_OUT}/run_likelihoods.csv"
        
        # Add jobs to list
        for run in $(seq 1 ${N_RUNS}); do
            echo "${model} ${run}" >> "${JOB_LIST}"
        done
    done
    
    echo "Total jobs: $(wc -l < ${JOB_LIST})"
    echo "Running with ${MAX_PARALLEL_JOBS} parallel jobs..."
    echo ""
    
    # Run all jobs in parallel
    cat "${JOB_LIST}" | parallel -j ${MAX_PARALLEL_JOBS} --colsep ' ' \
        run_single_optimization {1} {2} ${OUTPUT_DIR}
    
else
    echo "========================================================================"
    echo "RUNNING WITH BACKGROUND JOBS (NO GNU PARALLEL)"
    echo "========================================================================"
    echo ""
    
    for model in "${MODELS[@]}"; do
        echo ""
        echo "========================================================================"
        echo "MODEL: ${model}"
        echo "========================================================================"
        
        MODEL_OUT="${OUTPUT_DIR}/${model}"
        mkdir -p "${MODEL_OUT}"
        rm -f "${MODEL_OUT}"/run_*_likelihood.csv
        echo "run_id,max_est_likelihood" > "${MODEL_OUT}/run_likelihoods.csv"
        
        if [ ! -f "${MODEL_DIR}/${model}.tpl" ] || [ ! -f "${MODEL_DIR}/${model}.est" ]; then
            echo "ERROR: Model files not found for ${model}"
            continue
        fi
        
        echo "Running ${N_RUNS} optimizations with up to ${MAX_PARALLEL_JOBS} parallel jobs..."
        
        # Run in batches to avoid overwhelming the system
        for batch_start in $(seq 1 ${MAX_PARALLEL_JOBS} ${N_RUNS}); do
            batch_end=$((batch_start + MAX_PARALLEL_JOBS - 1))
            if [ ${batch_end} -gt ${N_RUNS} ]; then
                batch_end=${N_RUNS}
            fi
            
            echo "  Batch: runs ${batch_start}-${batch_end}"
            
            # Launch batch in background
            for run in $(seq ${batch_start} ${batch_end}); do
                run_single_optimization "${model}" "${run}" "${OUTPUT_DIR}" &
            done
            
            # Wait for batch to complete
            wait
        done
        
        echo ""
        echo "Completed ${N_RUNS} runs for ${model}"
    done
fi

# Merge per-run likelihood files after parallel execution to avoid concurrent
# appends corrupting run_likelihoods.csv.
for model in "${MODELS[@]}"; do
    MODEL_OUT="${OUTPUT_DIR}/${model}"
    if [ -d "${MODEL_OUT}" ]; then
        echo "run_id,max_est_likelihood" > "${MODEL_OUT}/run_likelihoods.csv"
        for run in $(seq 1 ${N_RUNS}); do
            RUN_LHOOD_FILE="${MODEL_OUT}/run_${run}_likelihood.csv"
            if [ -f "${RUN_LHOOD_FILE}" ]; then
                awk 'NF > 0 { print }' "${RUN_LHOOD_FILE}" >> "${MODEL_OUT}/run_likelihoods.csv"
            fi
        done
    fi
done

################################################################################
# Post-processing: Find best run for each model
################################################################################

echo ""
echo "========================================================================"
echo "POST-PROCESSING: FINDING BEST RUNS"
echo "========================================================================"
echo ""

for model in "${MODELS[@]}"; do
    MODEL_OUT="${OUTPUT_DIR}/${model}"
    
    if [ -f "${MODEL_OUT}/run_likelihoods.csv" ]; then
        # Sort by MaxEstLhood, take highest
        BEST_RUN=$(tail -n +2 "${MODEL_OUT}/run_likelihoods.csv" | \
                   grep -v "NA" | \
                   sort -t',' -k2 -rn | \
                   head -1 | \
                   cut -d',' -f1)
        
        if [ -n "${BEST_RUN}" ]; then
            echo "Model ${model}: Best run = ${BEST_RUN}"
            
            # Copy best run results
            BEST_DIR="${MODEL_OUT}/run_${BEST_RUN}/${model}"
            if [ -d "${BEST_DIR}" ]; then
                rm -rf "${MODEL_OUT}/best_run"
                cp -r "${BEST_DIR}" "${MODEL_OUT}/best_run"
                
                # Display best parameters
                if [ -f "${MODEL_OUT}/best_run/${model}.bestlhoods" ]; then
                    MAX_L=$(extract_max_est_lhood "${MODEL_OUT}/best_run/${model}.bestlhoods")
                    echo "  MaxEstLhood: ${MAX_L}"
                fi
            fi
        else
            echo "Model ${model}: WARNING - No successful runs"
        fi
    fi
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
