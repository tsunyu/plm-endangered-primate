#!/bin/bash
#
# Phase 3a: Population Genomics Analysis Pipeline
# Golden Snub-Nosed Monkey (Rhinopithecus roxellana) - Shennongjia Population
#
# This script orchestrates all Phase 3a analyses:
# 1. ROH analysis
# 2. Genetic diversity metrics
# 3. Effective population size estimation
# 4. Population structure analysis
#
# Usage: 
#   bash phase3a_population_genomics_pipeline.sh                    # Run all steps
#   bash phase3a_population_genomics_pipeline.sh --resume-from 2     # Resume from step 2
#   bash phase3a_population_genomics_pipeline.sh --step 3            # Run only step 3
#   bash phase3a_population_genomics_pipeline.sh --help              # Show help
#
# Author: Inbreeding Analysis Pipeline
# Date: 2025

set -euo pipefail

# ============================================================================
# COMMAND LINE ARGUMENTS
# ============================================================================

# Default values
RESUME_FROM=0
RUN_STEP=0
SHOW_HELP=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --resume-from)
            RESUME_FROM="$2"
            shift 2
            ;;
        --step)
            RUN_STEP="$2"
            shift 2
            ;;
        --help|-h)
            SHOW_HELP=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Show help if requested
if [ "$SHOW_HELP" = true ]; then
    cat << EOF
Phase 3a: Population Genomics Analysis Pipeline

USAGE:
    bash phase3a_population_genomics_pipeline.sh [OPTIONS]

OPTIONS:
    --resume-from N    Resume pipeline from step N (1-4)
    --step N          Run only step N (1-4)
    --help, -h        Show this help message

EXAMPLES:
    # Run all steps
    bash phase3a_population_genomics_pipeline.sh
    
    # Resume from step 2 (diversity metrics)
    bash phase3a_population_genomics_pipeline.sh --resume-from 2
    
    # Run only step 3 (Ne estimation)
    bash phase3a_population_genomics_pipeline.sh --step 3

STEPS:
    1. ROH Analysis
    2. Genetic Diversity Metrics  
    3. Effective Population Size Estimation
    4. Population Structure Analysis

EOF
    exit 0
fi

# Validate arguments
if [ "$RUN_STEP" -gt 0 ] && [ "$RESUME_FROM" -gt 0 ]; then
    echo "ERROR: Cannot specify both --step and --resume-from"
    exit 1
fi

if [ "$RUN_STEP" -lt 0 ] || [ "$RUN_STEP" -gt 4 ]; then
    echo "ERROR: --step must be between 1 and 4"
    exit 1
fi

if [ "$RESUME_FROM" -lt 0 ] || [ "$RESUME_FROM" -gt 4 ]; then
    echo "ERROR: --resume-from must be between 1 and 4"
    exit 1
fi

# ============================================================================
# SETUP
# ============================================================================

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Load configuration
CONFIG_FILE="${SCRIPT_DIR}/../config.yaml"
if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Configuration file not found: $CONFIG_FILE"
    exit 1
fi

# Parse configuration (basic bash parsing)
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/load_base_dir.sh"
DATA_DIR="${BASE_DIR}/data"
OUTPUT_DIR="${BASE_DIR}/output/phase3a_population_genomics"
PLINK_PREFIX="${DATA_DIR}/monkey_snp_sex_qc"
VCF="${DATA_DIR}/monkey_snp_sex_qc.vcf"

# Create output directories
mkdir -p "${OUTPUT_DIR}"/{roh_analysis,diversity_metrics,ne_estimation,population_structure,plots}

# Checkpoint file for tracking completed steps
CHECKPOINT_FILE="${OUTPUT_DIR}/checkpoints.txt"

# Log file
LOGFILE="${OUTPUT_DIR}/phase3a_pipeline.log"
exec > >(tee -a "$LOGFILE") 2>&1

# ============================================================================
# FUNCTIONS
# ============================================================================

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

check_command() {
    if ! command -v "$1" &> /dev/null; then
        log_message "ERROR: Required command '$1' not found in PATH"
        log_message "Please install $1 and try again"
        exit 1
    fi
}

check_file() {
    if [ ! -f "$1" ]; then
        log_message "ERROR: Required file not found: $1"
        exit 1
    fi
}

# Checkpoint management functions
mark_step_completed() {
    local step_num="$1"
    local step_name="$2"
    echo "${step_num}:${step_name}:$(date '+%Y-%m-%d %H:%M:%S')" >> "$CHECKPOINT_FILE"
}

is_step_completed() {
    local step_num="$1"
    if [ -f "$CHECKPOINT_FILE" ]; then
        grep -q "^${step_num}:" "$CHECKPOINT_FILE"
    else
        return 1
    fi
}

get_completed_steps() {
    if [ -f "$CHECKPOINT_FILE" ]; then
        cut -d: -f1 "$CHECKPOINT_FILE" | sort -n
    fi
}

show_checkpoint_status() {
    log_message "Checkpoint Status:"
    if [ -f "$CHECKPOINT_FILE" ]; then
        while IFS=: read -r step_num step_name timestamp; do
            log_message "  Step ${step_num}: ${step_name} (completed at ${timestamp})"
        done < "$CHECKPOINT_FILE"
    else
        log_message "  No completed steps found"
    fi
    log_message ""
}

run_step() {
    local step_num="$1"
    local step_name="$2"
    local script_path="$3"
    
    # Check if this step should be skipped
    if [ "$RUN_STEP" -gt 0 ] && [ "$RUN_STEP" -ne "$step_num" ]; then
        log_message "Skipping step ${step_num}: ${step_name} (running only step ${RUN_STEP})"
        return 0
    fi
    
    # Check if this step should be skipped due to resume-from
    if [ "$RESUME_FROM" -gt "$step_num" ]; then
        log_message "Skipping step ${step_num}: ${step_name} (resuming from step ${RESUME_FROM})"
        return 0
    fi
    
    # Check if step is already completed
    if is_step_completed "$step_num"; then
        log_message "Step ${step_num}: ${step_name} already completed - skipping"
        return 0
    fi
    
    log_message "====================================================================="
    log_message "Starting: Step ${step_num} - ${step_name}"
    log_message "====================================================================="
    
    if [ -f "$script_path" ]; then
        bash "$script_path"
        if [ $? -eq 0 ]; then
            log_message "SUCCESS: Step ${step_num} - ${step_name} completed"
            mark_step_completed "$step_num" "$step_name"
        else
            log_message "ERROR: Step ${step_num} - ${step_name} failed"
            log_message "To resume from this step after fixing the issue, run:"
            log_message "  bash $0 --resume-from $step_num"
            exit 1
        fi
    else
        log_message "ERROR: Script not found: $script_path"
        exit 1
    fi
    
    log_message ""
}

# ============================================================================
# PRE-FLIGHT CHECKS
# ============================================================================

log_message "====================================================================="
log_message "PHASE 3A: POPULATION GENOMICS ANALYSIS"
log_message "====================================================================="
log_message "Project: Golden Snub-Nosed Monkey Inbreeding Analysis"
log_message "Population: Shennongjia"
log_message "Species: Rhinopithecus roxellana"
log_message ""

# Show run mode
if [ "$RUN_STEP" -gt 0 ]; then
    log_message "RUN MODE: Running only step ${RUN_STEP}"
elif [ "$RESUME_FROM" -gt 0 ]; then
    log_message "RUN MODE: Resuming from step ${RESUME_FROM}"
else
    log_message "RUN MODE: Running all steps"
fi

# Show checkpoint status
show_checkpoint_status

log_message "Checking prerequisites..."

# Check required commands
REQUIRED_COMMANDS=("plink" "bcftools" "vcftools" "python3" "Rscript")
for cmd in "${REQUIRED_COMMANDS[@]}"; do
    check_command "$cmd"
done
log_message "All required commands found"

# Check input files
check_file "${PLINK_PREFIX}.bed"
check_file "${PLINK_PREFIX}.bim"
check_file "${PLINK_PREFIX}.fam"
check_file "${VCF}"
log_message "All input files found"

# Display input file info
log_message ""
log_message "Input Files:"
log_message "  PLINK prefix: ${PLINK_PREFIX}"
log_message "  VCF file: ${VCF}"

# Count samples
NUM_SAMPLES=$(wc -l < "${PLINK_PREFIX}.fam")
NUM_VARIANTS=$(wc -l < "${PLINK_PREFIX}.bim")

log_message ""
log_message "Dataset Information:"
log_message "  Number of samples: ${NUM_SAMPLES}"
log_message "  Number of variants: ${NUM_VARIANTS}"
log_message ""

# Prepare VCF (compress and index once for all analyses)
if [ -f "${VCF}" ]; then
    if [[ ! "${VCF}" =~ \.gz$ ]]; then
        log_message "Compressing VCF file for downstream analyses..."
        if [ ! -f "${VCF}.gz" ]; then
            log_message "Running bgzip compression..."
            if bgzip -c "${VCF}" > "${VCF}.gz"; then
                log_message "VCF compression successful"
            else
                log_message "ERROR: VCF compression failed"
                exit 1
            fi
        fi
        VCF="${VCF}.gz"
    fi
    
    # Verify compression format
    if ! file "${VCF}" | grep -q "BGZF"; then
        log_message "ERROR: VCF file is not in BGZF format. Recompressing..."
        rm -f "${VCF}"
        if bgzip -c "${VCF%.gz}" > "${VCF}"; then
            log_message "VCF recompression successful"
        else
            log_message "ERROR: VCF recompression failed"
            exit 1
        fi
    fi
    
    if [ ! -f "${VCF}.tbi" ]; then
        log_message "Indexing VCF file..."
        if tabix -p vcf "${VCF}"; then
            log_message "VCF indexing successful"
        else
            log_message "ERROR: VCF indexing failed"
            exit 1
        fi
    fi
    
    log_message "VCF prepared: ${VCF}"
    log_message ""
fi

# ============================================================================
# ANALYSIS STEPS
# ============================================================================

# Step 1: ROH Analysis
run_step 1 "ROH Analysis" "${SCRIPT_DIR}/phase3a_roh_analysis.sh"

# Step 2: Genetic Diversity Metrics
run_step 2 "Genetic Diversity Metrics" "${SCRIPT_DIR}/phase3_diversity_metrics.sh"

# Step 3: Effective Population Size Estimation
run_step 3 "Ne Estimation" "${SCRIPT_DIR}/phase3a_ne_estimation.sh"

# Step 4: Population Structure Analysis
run_step 4 "Population Structure" "${SCRIPT_DIR}/phase3a_population_structure.sh"

# ============================================================================
# GENERATE SUMMARY VISUALIZATIONS
# ============================================================================

log_message "====================================================================="
log_message "Generating summary visualizations..."
log_message "====================================================================="

if [ -f "${SCRIPT_DIR}/phase3a_generate_plots.R" ]; then
    Rscript "${SCRIPT_DIR}/phase3a_generate_plots.R" \
        --output-dir "${OUTPUT_DIR}/plots"
    
    if [ $? -eq 0 ]; then
        log_message "Summary plots generated successfully"
    else
        log_message "WARNING: Some plots may have failed to generate"
    fi
else
    log_message "WARNING: Plotting script not found"
fi

# ============================================================================
# GENERATE SUMMARY REPORT
# ============================================================================

log_message ""
log_message "====================================================================="
log_message "Generating Phase 3a Summary Report"
log_message "====================================================================="

SUMMARY_FILE="${OUTPUT_DIR}/phase3a_summary_report.txt"

cat > "$SUMMARY_FILE" << EOF
========================================================================
PHASE 3A: POPULATION GENOMICS ANALYSIS - SUMMARY REPORT
========================================================================

Project: Golden Snub-Nosed Monkey Inbreeding Analysis
Population: Shennongjia
Species: Rhinopithecus roxellana
Analysis Date: $(date '+%Y-%m-%d %H:%M:%S')

------------------------------------------------------------------------
DATASET INFORMATION
------------------------------------------------------------------------
Number of samples: ${NUM_SAMPLES}
Number of variants: ${NUM_VARIANTS}

------------------------------------------------------------------------
ANALYSES COMPLETED
------------------------------------------------------------------------
1. Runs of Homozygosity (ROH) Analysis
   - PLINK-based ROH detection
   - BCFtools RoH validation
   - F_ROH calculation
   - ROH categorization (short/medium/long)
   
2. Genetic Diversity Metrics
   - Individual heterozygosity
   - Nucleotide diversity (π) in 100kb windows
   - Tajima's D statistics
   
3. Effective Population Size (Ne) Estimation
   - LD-based Ne
   - Site frequency spectrum (SFS)
   - SFS-based demographic models (fastsimcoal2)
   
4. Population Structure Analysis
   - Principal Component Analysis (PCA)
   - Kinship analysis (KING)
   - Identity-by-descent (IBD) segments

------------------------------------------------------------------------
OUTPUT DIRECTORIES
------------------------------------------------------------------------
Main output: ${OUTPUT_DIR}

Subdirectories:
  - roh_analysis/          ROH results and statistics
  - diversity_metrics/     Heterozygosity, π, Tajima's D
  - ne_estimation/         Ne estimates and demographic models
  - population_structure/  PCA, kinship, IBD results
  - plots/                 Summary visualizations

------------------------------------------------------------------------
KEY FILES
------------------------------------------------------------------------
EOF

# Add key results files if they exist
if [ -f "${OUTPUT_DIR}/roh_analysis/roh_summary_per_individual.csv" ]; then
    echo "ROH Summary: roh_analysis/roh_summary_per_individual.csv" >> "$SUMMARY_FILE"
fi

if [ -f "${OUTPUT_DIR}/diversity_metrics/heterozygosity_summary.csv" ]; then
    echo "Heterozygosity: diversity_metrics/heterozygosity_summary.csv" >> "$SUMMARY_FILE"
fi

if [ -f "${OUTPUT_DIR}/ne_estimation/ne_estimates_summary.txt" ]; then
    echo "Ne Estimates: ne_estimation/ne_estimates_summary.txt" >> "$SUMMARY_FILE"
fi

cat >> "$SUMMARY_FILE" << EOF

------------------------------------------------------------------------
NEXT STEPS
------------------------------------------------------------------------
1. Review summary statistics and visualizations
2. Validate results across methods (e.g., PLINK vs BCFtools ROH)
3. Proceed to demography (fastsimcoal2)
4. Use Ne estimates to cross-check phase 3b demography

For detailed results, see individual analysis output directories.

========================================================================
END OF PHASE 3A SUMMARY
========================================================================
EOF

log_message "Summary report written to: ${SUMMARY_FILE}"

# Display summary
cat "$SUMMARY_FILE"

# ============================================================================
# COMPLETION
# ============================================================================

log_message ""
log_message "====================================================================="
log_message "PHASE 3A ANALYSIS COMPLETE!"
log_message "====================================================================="
log_message "Total runtime: $SECONDS seconds"
log_message "Output directory: ${OUTPUT_DIR}"
log_message "Log file: ${LOGFILE}"
log_message "Checkpoint file: ${CHECKPOINT_FILE}"
log_message ""

# Show final checkpoint status
log_message "Final Checkpoint Status:"
show_checkpoint_status

log_message "Review the summary report and proceed to Phase 3 when ready."
log_message ""
log_message "To clean up checkpoints and start fresh, run:"
log_message "  rm ${CHECKPOINT_FILE}"
log_message "====================================================================="

exit 0


