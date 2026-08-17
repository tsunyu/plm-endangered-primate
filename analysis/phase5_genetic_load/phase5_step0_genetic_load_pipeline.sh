#!/bin/bash
#
# Phase 5: Genetic Load Calculation Pipeline
#
# Calculates:
# 1. Individual genetic load (LOF and deleterious missense)
# 2. Population-level load statistics
# 3. Distribution of Fitness Effects (DFE)
# 4. Functional enrichment analysis
#
# Usage: bash phase5_step0_genetic_load_pipeline.sh [OPTIONS]
# Example: bash phase5_step0_genetic_load_pipeline.sh --only-step individual_load
#
# Options:
#   --force                Clear step checkpoints and rerun all steps
#   --clear-checkpoints    Clear all checkpoints before starting
#   --start-from STEP      Start from specific step (individual_load, population_load, enrichment)
#   --only-step STEP       Run only the specified step and skip all others
#   --help, -h             Show this help message
#

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/load_base_dir.sh"
OUTPUT_DIR="${BASE_DIR}/output/phase5_genetic_load"

mkdir -p "${OUTPUT_DIR}"/{individual_load,population_load,dfe,enrichment}
mkdir -p "${OUTPUT_DIR}/individual_load/visualizations"

# Checkpoint management
CHECKPOINT_DIR="${OUTPUT_DIR}/.checkpoints"
mkdir -p "${CHECKPOINT_DIR}"

# CLI flags
FORCE=false
CLEAR_CHECKPOINTS=false
START_FROM=""
ONLY_STEP=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --force)
            FORCE=true
            shift
            ;;
        --clear-checkpoints)
            CLEAR_CHECKPOINTS=true
            shift
            ;;
        --start-from)
            START_FROM="$2"
            shift 2
            ;;
        --only-step)
            ONLY_STEP="$2"
            shift 2
            ;;
        --help|-h)
            grep "^#" "$0" | grep -v "#!/bin/bash" | sed 's/^# *//'
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Clear checkpoints if requested
if [ "$CLEAR_CHECKPOINTS" = true ] || [ "$FORCE" = true ]; then
    echo "Clearing all checkpoints..."
    rm -rf "${CHECKPOINT_DIR}"
    mkdir -p "${CHECKPOINT_DIR}"
fi

LOGFILE="${OUTPUT_DIR}/phase5_pipeline.log"
exec > >(tee -a "$LOGFILE") 2>&1

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

run_step() {
    local step_name="$1"
    local step_key="$2"
    local script_path="$3"
    local checkpoint_file="${CHECKPOINT_DIR}/${step_key}.done"

    log_message "====================================================================="
    log_message "Step: ${step_name}"
    log_message "====================================================================="

    # If only a single step is requested, skip all others
    if [ -n "$ONLY_STEP" ] && [ "$step_key" != "$ONLY_STEP" ]; then
        log_message "SKIPPING: Only running step ($ONLY_STEP)"
        return 0
    fi

    # Check if we should start from this step
    if [ -n "$START_FROM" ]; then
        if [ "$step_key" != "$START_FROM" ] && [ ! -f "$checkpoint_file" ]; then
            log_message "SKIPPING: Not yet reached start point (--start-from $START_FROM)"
            return 0
        fi
    fi

    # Skip if checkpoint exists and not forcing
    if [ -f "$checkpoint_file" ] && [ "$FORCE" = false ]; then
        local timestamp=$(cat "$checkpoint_file")
        log_message "CHECKPOINT: ${step_name} already completed at ${timestamp}"
        log_message "SKIPPING: Use --force to rerun"
        return 0
    fi

    log_message "RUNNING: ${step_name}"

    if [[ "$script_path" =~ \.py$ ]]; then
        python3 "$script_path"
    elif [[ "$script_path" =~ \.R$ ]]; then
        Rscript "$script_path"
    else
        bash "$script_path"
    fi

    local exit_code=$?
    if [ $exit_code -eq 0 ]; then
        date '+%Y-%m-%d %H:%M:%S' > "$checkpoint_file"
        log_message "SUCCESS: ${step_name} completed"
    else
        log_message "ERROR: ${step_name} failed with exit code $exit_code"
        exit 1
    fi
}

log_message "====================================================================="
log_message "PHASE 5: GENETIC LOAD CALCULATION"
log_message "====================================================================="
log_message ""

if [ "$FORCE" = true ]; then
    log_message "Force mode: All steps will be rerun"
fi
if [ -n "$START_FROM" ]; then
    log_message "Starting from step: $START_FROM"
fi

# Step 1: Calculate individual genetic load (includes visualizations)
run_step "Individual Genetic Load" "individual_load" "${SCRIPT_DIR}/phase5_step1_calculate_individual_load.py"

# Step 2: Population-level statistics (Python)
run_step "Population Load Analysis" "population_load" "${SCRIPT_DIR}/phase5_step2_population_load_analysis.py"



log_message ""
log_message "====================================================================="
log_message "PHASE 5 ANALYSIS COMPLETE"
log_message "====================================================================="
log_message "Output directory: ${OUTPUT_DIR}"
log_message ""
log_message "Checkpoints saved in: ${CHECKPOINT_DIR}"
log_message "To rerun all steps: bash phase5_step0_genetic_load_pipeline.sh --force"
log_message ""
log_message "Next: phenotype / GWAS (pipeline step 6) or method validation (step 7)"
log_message "====================================================================="

exit 0


