#!/bin/bash
#
# Phase 4: Protein Language Model Predictions Pipeline
#
# Runs mutation effect predictions using:
# - ESM-2 (structure-aware)
#
# Note: Uses deduplicated variants (missense_variants_for_prediction_unique.csv)
#       to avoid redundant predictions on the same genomic position.
#       Missense variants are scored with ESM-2 only.
#
# Usage: bash phase4_step0_plm_pipeline.sh [OPTIONS]
#
# Options:
#   --force                Clear all checkpoints and force rerun all steps
#   --clear-checkpoints    Clear all checkpoints before starting
#   --start-from STEP      Start from specific step (sequence_prep, esm2, ensemble, visualization)
#   --only-step STEP       Run only the specified step and skip all others
#   --help, -h             Show this help message
#
# Examples:
#   bash phase4_step0_plm_pipeline.sh                    # Normal run with checkpointing
#   bash phase4_step0_plm_pipeline.sh --force            # Force rerun all steps
#   bash phase4_step0_plm_pipeline.sh --start-from esm2  # Start from ESM-2 predictions
#

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/load_base_dir.sh"
OUTPUT_DIR="${BASE_DIR}/output/phase4_plm_predictions"

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

mkdir -p "${OUTPUT_DIR}"/{esm2,ensemble}

LOGFILE="${OUTPUT_DIR}/phase4_pipeline.log"
exec > >(tee -a "$LOGFILE") 2>&1

# ============================================================================
# FUNCTIONS
# ============================================================================

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
    
    # Check checkpoint: skip if already completed
    if [ -f "$checkpoint_file" ] && [ "$FORCE" = false ]; then
        local timestamp=$(cat "$checkpoint_file")
        log_message "CHECKPOINT: ${step_name} already completed at ${timestamp}"
        log_message "SKIPPING: Use --force to rerun"
        return 0
    fi
    
    log_message "RUNNING: ${step_name}"
    
    if [ -f "$script_path" ]; then
        python3 "$script_path"
        local exit_code=$?
        
        if [ $exit_code -eq 0 ]; then
            # Mark checkpoint
            date '+%Y-%m-%d %H:%M:%S' > "$checkpoint_file"
            log_message "SUCCESS: ${step_name} completed"
        else
            log_message "ERROR: ${step_name} failed with exit code $exit_code"
            exit 1
        fi
    else
        log_message "ERROR: Script not found: $script_path"
        exit 1
    fi
}

run_bash_block() {
    local step_name="$1"
    local step_key="$2"
    local checkpoint_file="${CHECKPOINT_DIR}/${step_key}.done"
    
    log_message "====================================================================="
    log_message "Step: ${step_name}"
    log_message "====================================================================="
    
    # If only a single step is requested, skip all others
    if [ -n "$ONLY_STEP" ] && [ "$step_key" != "$ONLY_STEP" ]; then
        log_message "SKIPPING: Only running step ($ONLY_STEP)"
        return 1  # Return 1 to skip (don't execute block)
    fi

    # Check if we should start from this step
    if [ -n "$START_FROM" ]; then
        if [ "$step_key" != "$START_FROM" ] && [ ! -f "$checkpoint_file" ]; then
            log_message "SKIPPING: Not yet reached start point (--start-from $START_FROM)"
            return 1  # Return 1 to skip
        fi
    fi
    
    # Check checkpoint
    if [ -f "$checkpoint_file" ] && [ "$FORCE" = false ]; then
        local timestamp=$(cat "$checkpoint_file")
        log_message "CHECKPOINT: ${step_name} already completed at ${timestamp}"
        log_message "SKIPPING: Use --force to rerun"
        return 1  # Return 1 to skip
    fi
    
    log_message "RUNNING: ${step_name}"
    return 0  # Return 0 to execute the block
}

log_message "====================================================================="
log_message "PHASE 4: PROTEIN LANGUAGE MODEL PREDICTIONS"
log_message "====================================================================="
log_message ""

if [ -n "$START_FROM" ]; then
    log_message "Starting from step: $START_FROM"
fi
if [ "$FORCE" = true ]; then
    log_message "Force mode: All steps will be rerun"
fi
log_message ""

# Step 1: Prepare sequences and variants
run_step "Sequence Preparation" "sequence_prep" "${SCRIPT_DIR}/phase4_step1_prepare_sequences.py"

# Step 2: ESM-2 predictions (only model used)
run_step "ESM-2 Predictions" "esm2" "${SCRIPT_DIR}/phase4_step2_esm2_predictions.py"

# Step 3: Create compatible output for phase5 (directly using ESM-2 scores, no ensemble integration)
if run_bash_block "Create Phase 5 Compatible Output" "ensemble"; then
python3 - << 'EOF'
import pandas as pd
import sys
import os

BASE_DIR = os.environ["PLM_BASE_DIR"]
ESM2_OUTPUT = f"{BASE_DIR}/output/phase4_plm_predictions/esm2/esm2_predictions.csv"
ENSEMBLE_OUTPUT = f"{BASE_DIR}/output/phase4_plm_predictions/ensemble/ensemble_predictions.csv"

# Load ESM-2 predictions
df = pd.read_csv(ESM2_OUTPUT)

# Create compatible format for phase5 (directly using ESM-2 scores, no ensemble)
# Normalize ESM-2 score to 0-1 scale (lower ESM-2 score = more deleterious)
score_col = 'esm2_score'
# Filter out NaN values for normalization
valid_scores = df[score_col].dropna()
if len(valid_scores) == 0:
    print("ERROR: No valid ESM-2 scores found!")
    sys.exit(1)

min_val = valid_scores.min()
max_val = valid_scores.max()

# Handle case where all scores are the same
if max_val == min_val:
    print("WARNING: All ESM-2 scores are identical, using uniform ensemble scores")
    df['ensemble_score'] = 0.5
else:
    # Invert so lower (more negative) scores become higher normalized values
    df['ensemble_score'] = 1 - (df[score_col] - min_val) / (max_val - min_val)
    # Fill NaN values with median
    df['ensemble_score'] = df['ensemble_score'].fillna(df['ensemble_score'].median())

# Percentile ranking
df['ensemble_percentile'] = df['ensemble_score'].rank(pct=True)

# Classify as deleterious (top 10%)
threshold = df['ensemble_score'].quantile(0.90)
df['is_deleterious'] = df['ensemble_score'] >= threshold

# Ensure variant_id column exists (used by phase5)
if 'variant_id' not in df.columns:
    # Create variant_id from available columns if missing
    if all(col in df.columns for col in ['chrom', 'pos', 'ref', 'alt']):
        # Format: chrom:pos:ref:alt (same as phase4_step1_prepare_sequences.py)
        df['variant_id'] = df['chrom'].astype(str) + ':' + df['pos'].astype(str) + ':' + df['ref'] + ':' + df['alt']
        print("Created variant_id column from chrom, pos, ref, alt")
    else:
        print("WARNING: Cannot create variant_id column. Please check ESM-2 output format.")
        sys.exit(1)
else:
    print(f"Using existing variant_id column: {len(df['variant_id'].unique())} unique variants")

# Save as ensemble output for phase5
os.makedirs(os.path.dirname(ENSEMBLE_OUTPUT), exist_ok=True)
df.to_csv(ENSEMBLE_OUTPUT, index=False)
print(f"Created compatible output: {ENSEMBLE_OUTPUT}")
print(f"  Total variants: {len(df)}")
print(f"  Deleterious variants (top 10%): {df['is_deleterious'].sum()}")
EOF

    if [ $? -eq 0 ]; then
        # Mark checkpoint
        date '+%Y-%m-%d %H:%M:%S' > "${CHECKPOINT_DIR}/ensemble.done"
        log_message "SUCCESS: Compatible output created"
    else
        log_message "ERROR: Failed to create compatible output"
        exit 1
    fi
fi

# Step 4: Generate visualizations
run_step "Generate Visualizations" "visualization" "${SCRIPT_DIR}/phase4_step3_esm2_visualizations.py"

log_message ""
log_message "====================================================================="
log_message "PHASE 4 ANALYSIS COMPLETE"
log_message "====================================================================="
log_message "Output directory: ${OUTPUT_DIR}"
log_message ""
log_message "Key output files:"
log_message "  - esm2/esm2_predictions.csv              : ESM-2 raw predictions"
log_message "  - ensemble/ensemble_predictions.csv      : Phase 5 compatible format (based on ESM-2 only)"
log_message "  - esm2/visualizations/*.png              : Visualization plots (based on ESM-2 scores)"
log_message ""
log_message "Checkpoints saved in: ${CHECKPOINT_DIR}"
log_message "To rerun all steps: bash phase4_step0_plm_pipeline.sh --force"
log_message "To start from a specific step: bash phase4_step0_plm_pipeline.sh --start-from STEP"
log_message ""
log_message "Proceed to Phase 5: Genetic Load Calculation"
log_message "====================================================================="

exit 0


