#!/bin/bash
################################################################################
# Complete fastsimcoal2 Demographic Inference Pipeline
################################################################################
#
# This master script runs the complete demographic inference pipeline
# with error checking and progress reporting.
#
# Usage:
#   bash run_complete_pipeline.sh [OPTIONS]
#
# Options:
#   --skip-bootstrap    Skip bootstrap CI calculation (faster, no CI)
#   --skip-step0        Skip Step 0 (SFS preparation)
#   --skip-step0b       Skip Step 0b (monomorphic site estimation)
#   --test             Run with reduced settings for testing
#   --no-parallel      Use serial version instead of parallel (slower)
#   --serial           Same as --no-parallel
#
# Author: Demographic Analysis Pipeline
# Date: 2026-01-26
################################################################################

set -e  # Exit on error

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/load_base_dir.sh"

# Parse arguments
SKIP_BOOTSTRAP=false
SKIP_STEP0=false
SKIP_STEP0B=false
TEST_MODE=false
USE_PARALLEL=true  # Default: use parallel version for better performance

for arg in "$@"; do
    case $arg in
        --skip-bootstrap)
            SKIP_BOOTSTRAP=true
            shift
            ;;
        --skip-step0)
            SKIP_STEP0=true
            shift
            ;;
        --skip-step0b)
            SKIP_STEP0B=true
            shift
            ;;
        --test)
            TEST_MODE=true
            shift
            ;;
        --no-parallel)
            USE_PARALLEL=false
            shift
            ;;
        --serial)
            USE_PARALLEL=false
            shift
            ;;
        *)
            ;;
    esac
done

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${BASE_DIR}/output/phase3b_fastsimcoal2"
LOG_FILE="${OUTPUT_DIR}/pipeline.log"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Logging function
log_message() {
    local level=$1
    shift
    local message="$@"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    case $level in
        INFO)
            echo -e "${BLUE}[${timestamp}] ${message}${NC}" | tee -a "${LOG_FILE}"
            ;;
        SUCCESS)
            echo -e "${GREEN}[${timestamp}] ✓ ${message}${NC}" | tee -a "${LOG_FILE}"
            ;;
        WARNING)
            echo -e "${YELLOW}[${timestamp}] ⚠ ${message}${NC}" | tee -a "${LOG_FILE}"
            ;;
        ERROR)
            echo -e "${RED}[${timestamp}] ✗ ${message}${NC}" | tee -a "${LOG_FILE}"
            ;;
    esac
}

# Check prerequisites
check_prerequisites() {
    log_message INFO "Checking prerequisites..."
    
    # Check fastsimcoal2 (prefer fsc28, fallback to fsc27/fsc26)
    if ! command -v fsc28 &> /dev/null && ! command -v fsc27 &> /dev/null && ! command -v fsc26 &> /dev/null; then
        log_message ERROR "fastsimcoal2 not found in PATH"
        log_message ERROR "Please install fastsimcoal2 (fsc28, fsc27 or fsc26)"
        exit 1
    fi
    
    # Check Python
    if ! command -v python3 &> /dev/null; then
        log_message ERROR "python3 not found"
        exit 1
    fi
    
    # Check required Python packages
    python3 -c "import numpy, pandas, matplotlib" 2>/dev/null
    if [ $? -ne 0 ]; then
        log_message WARNING "Some Python packages may be missing"
        log_message WARNING "Install with: pip install numpy pandas matplotlib"
    fi
    
    # Check GNU parallel for optimal performance
    if [ "$USE_PARALLEL" = true ]; then
        if ! command -v parallel &> /dev/null; then
            log_message WARNING "GNU parallel not found (parallel mode will use fallback)"
            log_message WARNING "For best performance, install with: sudo apt-get install parallel"
        else
            log_message INFO "GNU parallel detected: optimal performance enabled"
        fi
    fi
    
    # Check input files
    if [ ! -f "${BASE_DIR}/data/hardfilted.snp.pass.autosomes.vcf.gz" ]; then
        log_message ERROR "Main VCF not found"
        exit 1
    fi
    
    if [ ! -f "${BASE_DIR}/output/phase2_annotation/ancestral_states/variants_with_ancestral.vcf.gz" ]; then
        log_message ERROR "Ancestral VCF not found"
        exit 1
    fi
    
    log_message SUCCESS "Prerequisites check passed"
}

# Step 0: Prepare SFS
run_step0() {
    log_message INFO "==================================================================="
    log_message INFO "STEP 0: Preparing Site Frequency Spectrum"
    log_message INFO "==================================================================="
    
    cd "${SCRIPT_DIR}"
    
    if python3 phase3b_step0_prepare_sfs.py; then
        log_message SUCCESS "Step 0 completed: SFS generated"
    else
        log_message ERROR "Step 0 failed"
        exit 1
    fi
}

# Step 0b: Estimate monomorphic sites
run_step0b() {
    log_message INFO "==================================================================="
    log_message INFO "STEP 0b: Estimating Monomorphic Sites (Recommended)"
    log_message INFO "==================================================================="
    
    cd "${SCRIPT_DIR}"
    
    if python3 phase3b_step0b_estimate_monomorphic_sites.py; then
        log_message SUCCESS "Step 0b completed: Monomorphic sites estimated and SFS updated"
    else
        log_message WARNING "Step 0b failed (non-fatal, continuing without monomorphic sites)"
    fi
}

# Step 1: Create model templates
run_step1() {
    log_message INFO "==================================================================="
    log_message INFO "STEP 1: Creating Demographic Model Templates"
    log_message INFO "==================================================================="
    
    cd "${SCRIPT_DIR}"
    
    if python3 phase3b_step1_create_model_templates.py; then
        log_message SUCCESS "Step 1 completed: Model templates created"
    else
        log_message ERROR "Step 1 failed"
        exit 1
    fi
}

# Step 2: Run fastsimcoal2
run_step2() {
    log_message INFO "==================================================================="
    log_message INFO "STEP 2: Running fastsimcoal2 Parameter Estimation"
    log_message INFO "==================================================================="
    
    # Determine which script to use
    local step2_script="phase3b_step2_run_fastsimcoal2.sh"
    if [ "$USE_PARALLEL" = true ]; then
        if [ -f "${SCRIPT_DIR}/phase3b_step2_run_fastsimcoal2_parallel.sh" ]; then
            step2_script="phase3b_step2_run_fastsimcoal2_parallel.sh"
            log_message INFO "Using PARALLEL version (8-10x speedup)"
            log_message INFO "  CPU cores: $(nproc)"
            log_message INFO "  Parallel jobs: 16"
            log_message INFO "  Cores per job: 2"
        else
            log_message WARNING "Parallel script not found, using serial version"
            step2_script="phase3b_step2_run_fastsimcoal2.sh"
        fi
    else
        log_message INFO "Using SERIAL version (--no-parallel flag set)"
    fi
    
    if [ "$TEST_MODE" = true ]; then
        log_message WARNING "Running in TEST MODE (fewer runs, faster)"
        export N_RUNS=10
    else
        log_message INFO "Running full analysis (50 runs per model)"
        if [ "$USE_PARALLEL" = true ] && [ "$step2_script" = "phase3b_step2_run_fastsimcoal2_parallel.sh" ]; then
            log_message INFO "Estimated time: ~4-5 hours (parallel)"
        else
            log_message INFO "Estimated time: ~30-40 hours (serial)"
        fi
    fi
    
    cd "${SCRIPT_DIR}"
    
    if bash "${step2_script}"; then
        log_message SUCCESS "Step 2 completed: Parameter estimation finished"
    else
        log_message ERROR "Step 2 failed"
        exit 1
    fi
}

# Step 3: Model comparison
run_step3() {
    log_message INFO "==================================================================="
    log_message INFO "STEP 3: Comparing Demographic Models"
    log_message INFO "==================================================================="
    
    cd "${SCRIPT_DIR}"
    
    if python3 phase3b_step3_model_comparison.py; then
        log_message SUCCESS "Step 3 completed: Model comparison finished"
        
        # Display best model
        if [ -f "${OUTPUT_DIR}/model_comparison/model_comparison.csv" ]; then
            BEST_MODEL=$(tail -n +2 "${OUTPUT_DIR}/model_comparison/model_comparison.csv" | head -1 | cut -d',' -f1)
            log_message INFO "Best model: ${BEST_MODEL}"
        fi
    else
        log_message ERROR "Step 3 failed"
        exit 1
    fi
}

# Step 4: Bootstrap CI
run_step4() {
    if [ "$SKIP_BOOTSTRAP" = true ]; then
        log_message WARNING "Skipping bootstrap (--skip-bootstrap flag set)"
        return 0
    fi
    
    log_message INFO "==================================================================="
    log_message INFO "STEP 4: Bootstrap Confidence Intervals"
    log_message INFO "==================================================================="
    
    # Determine which script to use
    local step4_script="phase3b_step4_bootstrap_ci.sh"
    if [ "$USE_PARALLEL" = true ]; then
        if [ -f "${SCRIPT_DIR}/phase3b_step4_bootstrap_ci_parallel.sh" ]; then
            step4_script="phase3b_step4_bootstrap_ci_parallel.sh"
            log_message INFO "Using PARALLEL version (8-10x speedup)"
        else
            log_message WARNING "Parallel bootstrap script not found, using serial version"
            step4_script="phase3b_step4_bootstrap_ci.sh"
        fi
    else
        log_message INFO "Using SERIAL version (--no-parallel flag set)"
    fi
    
    if [ "$TEST_MODE" = true ]; then
        log_message WARNING "Running in TEST MODE (fewer bootstrap replicates)"
        export N_BOOTSTRAP=20
    else
        log_message INFO "Running 100 bootstrap replicates"
        if [ "$USE_PARALLEL" = true ] && [ "$step4_script" = "phase3b_step4_bootstrap_ci_parallel.sh" ]; then
            log_message INFO "Estimated time: ~1-2 hours (parallel)"
        else
            log_message INFO "Estimated time: ~15-20 hours (serial)"
        fi
    fi
    
    cd "${SCRIPT_DIR}"
    
    if bash "${step4_script}"; then
        log_message SUCCESS "Step 4 completed: Bootstrap finished"
    else
        log_message WARNING "Step 4 failed (bootstrap errors are non-fatal)"
    fi
}

# Step 5: Analyze results
run_step5() {
    log_message INFO "==================================================================="
    log_message INFO "STEP 5: Analyzing Results"
    log_message INFO "==================================================================="
    
    cd "${SCRIPT_DIR}"
    
    if python3 phase3b_step5_analyze_results.py; then
        log_message SUCCESS "Step 5 completed: Results analyzed"
        
        # Display key results
        if [ -f "${OUTPUT_DIR}/parameter_estimates.txt" ]; then
            log_message INFO "Parameter estimates saved to: parameter_estimates.txt"
        fi
    else
        log_message ERROR "Step 5 failed"
        exit 1
    fi
}

# Step 6: Visualize
run_step6() {
    log_message INFO "==================================================================="
    log_message INFO "STEP 6: Creating Visualizations"
    log_message INFO "==================================================================="
    
    cd "${SCRIPT_DIR}"
    
    if python3 phase3b_step6_visualize_demographic.py; then
        log_message SUCCESS "Step 6 completed: Visualizations created"
        
        # List created plots
        if [ -d "${OUTPUT_DIR}/plots" ]; then
            log_message INFO "Plots saved to: ${OUTPUT_DIR}/plots/"
            ls -1 "${OUTPUT_DIR}/plots"/*.png 2>/dev/null | while read plot; do
                log_message INFO "  - $(basename ${plot})"
            done
        fi
    else
        log_message ERROR "Step 6 failed"
        exit 1
    fi
}

# Main execution
main() {
    echo ""
    log_message INFO "==================================================================="
    log_message INFO "fastsimcoal2 DEMOGRAPHIC INFERENCE PIPELINE"
    log_message INFO "==================================================================="
    echo ""
    
    START_TIME=$(date +%s)
    
    log_message INFO "Configuration:"
    log_message INFO "  Parallel mode: ${USE_PARALLEL}"
    log_message INFO "  Skip bootstrap: ${SKIP_BOOTSTRAP}"
    log_message INFO "  Skip Step 0: ${SKIP_STEP0}"
    log_message INFO "  Skip Step 0b: ${SKIP_STEP0B}"
    log_message INFO "  Test mode: ${TEST_MODE}"
    log_message INFO "  CPU cores available: $(nproc)"
    log_message INFO "  Output directory: ${OUTPUT_DIR}"
    log_message INFO "  Log file: ${LOG_FILE}"
    echo ""
    
    if [ "$USE_PARALLEL" = true ]; then
        log_message INFO "Performance optimization:"
        log_message INFO "  ✓ Parallel execution enabled (8-10x speedup)"
        log_message INFO "  ✓ Expected total time: ~6-8 hours"
        log_message INFO "  To use serial version: add --no-parallel flag"
    else
        log_message WARNING "Running in SERIAL mode"
        log_message WARNING "  Expected total time: ~50-60 hours"
        log_message WARNING "  For faster execution, remove --no-parallel flag"
    fi
    echo ""
    
    # Run checks
    check_prerequisites
    echo ""
    
    # Run pipeline steps
    if [ "$SKIP_STEP0" = true ]; then
        log_message WARNING "Skipping Step 0 (--skip-step0 flag set)"
    else
        run_step0
    fi
    echo ""
    
    if [ "$SKIP_STEP0B" = true ]; then
        log_message WARNING "Skipping Step 0b (--skip-step0b flag set)"
    else
        run_step0b
    fi
    echo ""
    
    run_step1
    echo ""
    
    run_step2
    echo ""
    
    run_step3
    echo ""
    
    run_step4
    echo ""
    
    run_step5
    echo ""
    
    run_step6
    echo ""
    
    # Calculate total time
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    HOURS=$((ELAPSED / 3600))
    MINUTES=$(( (ELAPSED % 3600) / 60 ))
    
    log_message INFO "==================================================================="
    log_message SUCCESS "PIPELINE COMPLETE!"
    log_message INFO "==================================================================="
    echo ""
    log_message INFO "Total time: ${HOURS}h ${MINUTES}m"
    echo ""
    log_message INFO "Results summary:"
    log_message INFO "  - Parameter estimates: ${OUTPUT_DIR}/parameter_estimates.txt"
    log_message INFO "  - Model comparison: ${OUTPUT_DIR}/model_comparison/"
    log_message INFO "  - Plots: ${OUTPUT_DIR}/plots/"
    log_message INFO "  - Full log: ${LOG_FILE}"
    echo ""
    log_message INFO "Next steps:"
    log_message INFO "  1. Review ${OUTPUT_DIR}/parameter_estimates.txt"
    log_message INFO "  2. Check plots in ${OUTPUT_DIR}/plots/"
    log_message INFO "  3. Verify Ne scaling with verify_ne_scaling.py"
    log_message INFO "  4. Use parameter estimates and plots with the paper figures"
    echo ""
}

# Run main function
main

exit 0
