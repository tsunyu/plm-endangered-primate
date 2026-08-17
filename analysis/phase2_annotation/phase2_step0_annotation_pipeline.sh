#!/bin/bash
#
# Phase 2: Ancestral State and Variant Annotation Pipeline
#
# This script orchestrates:
# 1. Ancestral state inference from outgroups
# 2. SnpEff functional annotation
# 3. Additional functional summaries from SnpEff annotations
# 4. Integration of all annotations
#
# Usage: bash phase2_step0_annotation_pipeline.sh
#

set -euo pipefail

# ============================================================================
# SETUP
# ============================================================================

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Load memory configuration
source "${SCRIPT_DIR}/../memory_config.sh"

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/load_base_dir.sh"
DATA_DIR="${BASE_DIR}/data"
OUTPUT_DIR="${BASE_DIR}/output/phase2_annotation"
VCF="${DATA_DIR}/monkey_snp_sex_qc.vcf"

# Checkpointing
CHECKPOINT_DIR="${OUTPUT_DIR}/.checkpoints"
mkdir -p "${CHECKPOINT_DIR}"

# CLI flags
FORCE=false
CLEAR_CHECKPOINTS=false
RESUME_FROM=""

# Parse arguments
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
        --resume-from)
            RESUME_FROM="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

if [ "$CLEAR_CHECKPOINTS" = true ]; then
    rm -rf "${CHECKPOINT_DIR}"
    mkdir -p "${CHECKPOINT_DIR}"
fi

# Memory monitoring
MAX_MEMORY_GB=120  # Reserve 8GB for system
MEMORY_CHECK_INTERVAL=60  # Check memory every 60 seconds

mkdir -p "${OUTPUT_DIR}"/{ancestral_states,snpeff_annotation,functional_annotation,integrated}

LOGFILE="${OUTPUT_DIR}/phase2_pipeline.log"
exec > >(tee -a "$LOGFILE") 2>&1

# ============================================================================
# FUNCTIONS
# ============================================================================

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

check_memory() {
    local used_memory=$(free -g | awk 'NR==2{print $3}')
    if [ "$used_memory" -gt "$MAX_MEMORY_GB" ]; then
        log_message "WARNING: High memory usage detected: ${used_memory}GB"
        log_message "Consider reducing parallel processes or increasing swap"
        return 1
    fi
    return 0
}

run_step() {
    local step_name="$1"
    local script_path="$2"
    local step_key=$(echo "$step_name" | tr ' ' '_' | tr 'A-Z' 'a-z')
    local checkpoint_file="${CHECKPOINT_DIR}/${step_key}.done"
    
    log_message "====================================================================="
    log_message "Starting: ${step_name}"
    log_message "====================================================================="
    
    # Resume logic: skip until reaching the requested step
    if [ -n "$RESUME_FROM" ] && [ "${step_name}" != "$RESUME_FROM" ] && [ ! -f "$CHECKPOINT_FILE" ]; then
        # If resuming from a later step, and this step isn't done, skip
        log_message "RESUME: Skipping '${step_name}' until '${RESUME_FROM}'"
        return 0
    fi

    # Checkpoint: skip completed steps unless forced
    if [ -f "$checkpoint_file" ] && [ "$FORCE" = false ]; then
        log_message "CHECKPOINT: '${step_name}' already completed. Skipping. (use --force to rerun)"
        return 0
    fi

    # Check memory before starting
    if ! check_memory; then
        log_message "ERROR: Insufficient memory for ${step_name}"
        exit 1
    fi
    
    if [ -f "$script_path" ]; then
        # Start memory monitoring in background
        (
            while true; do
                sleep $MEMORY_CHECK_INTERVAL
                check_memory || break
            done
        ) &
        local monitor_pid=$!
        
        bash "$script_path"
        local exit_code=$?
        
        # Stop memory monitoring
        kill $monitor_pid 2>/dev/null || true
        
        if [ $exit_code -eq 0 ]; then
            log_message "SUCCESS: ${step_name} completed"
            # Write checkpoint
            date > "$checkpoint_file"
        else
            log_message "ERROR: ${step_name} failed"
            log_message "Cleaning up temporary files..."
            find "${OUTPUT_DIR}" -name "*.tmp" -delete 2>/dev/null || true
            exit 1
        fi
    else
        log_message "ERROR: Script not found: $script_path"
        exit 1
    fi
}

# ============================================================================
# ANALYSIS
# ============================================================================

log_message "====================================================================="
log_message "PHASE 2: ANCESTRAL STATE AND VARIANT ANNOTATION"
log_message "====================================================================="
log_message ""

# Step 1: Ancestral State Inference
run_step "Ancestral State Inference" "${SCRIPT_DIR}/phase2_step1_ancestral_states.sh"

# Step 2: SnpEff Annotation
run_step "SnpEff Annotation" "${SCRIPT_DIR}/phase2_step2_snpeff_annotation.sh"

# Step 3: Functional Annotation
run_step "Additional Functional Annotation" "${SCRIPT_DIR}/phase2_step3_functional_annotation.sh"

# ============================================================================
# SUMMARY
# ============================================================================

log_message ""
log_message "====================================================================="
log_message "PHASE 2 ANALYSIS COMPLETE"
log_message "====================================================================="
log_message "Output directory: ${OUTPUT_DIR}"
log_message ""
log_message "Proceed to Phase 4: Protein Language Model Predictions"
log_message "====================================================================="

exit 0


