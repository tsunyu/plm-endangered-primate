#!/bin/bash
#
# Phase 3: Create All Annotation Figures
#
# Unified script to generate all annotation visualization figures
#
# Usage: bash phase2_step4_create_figures.sh [options]
#
# Options:
#   --snpeff-only        Generate only SnpEff figures
#   --functional-only    Generate only functional annotation figures
#   --dpi DPI           Figure DPI (default: 300)
#   --format FORMAT     Figure format: png, pdf, svg (default: png)
#   --help              Show this help message
#

set -euo pipefail

# ============================================================================
# PARAMETERS
# ============================================================================

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/load_base_dir.sh"
SNPEFF_VCF="${BASE_DIR}/output/phase2_annotation/snpeff_annotation/annotated_variants.vcf.gz"
SNPEFF_OUTPUT="${BASE_DIR}/output/phase2_annotation/snpeff_annotation"
FUNCTIONAL_DIR="${BASE_DIR}/output/phase2_annotation/functional_annotation"

DPI=300
FORMAT="png"
SNPEFF_ONLY=false
FUNCTIONAL_ONLY=false

# ============================================================================
# PARSE ARGUMENTS
# ============================================================================

while [[ $# -gt 0 ]]; do
    case $1 in
        --snpeff-only)
            SNPEFF_ONLY=true
            shift
            ;;
        --functional-only)
            FUNCTIONAL_ONLY=true
            shift
            ;;
        --dpi)
            DPI="$2"
            shift 2
            ;;
        --format)
            FORMAT="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: $0 [options]"
            echo ""
            echo "Options:"
            echo "  --snpeff-only        Generate only SnpEff figures"
            echo "  --functional-only    Generate only functional annotation figures"
            echo "  --dpi DPI           Figure DPI (default: 300)"
            echo "  --format FORMAT     Figure format: png, pdf, svg (default: png)"
            echo "  --help              Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# ============================================================================
# FUNCTIONS
# ============================================================================

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# ============================================================================
# MAIN
# ============================================================================

log_message "====================================================================="
log_message "PHASE 3 ANNOTATION VISUALIZATION"
log_message "====================================================================="
log_message ""
log_message "Configuration:"
log_message "  DPI: $DPI"
log_message "  Format: $FORMAT"
log_message ""

# ============================================================================
# SNPEFF VISUALIZATION
# ============================================================================

if [ "$FUNCTIONAL_ONLY" = false ]; then
    log_message "====================================================================="
    log_message "SnpEff Annotation Visualization"
    log_message "====================================================================="
    log_message ""
    
    if [ ! -f "$SNPEFF_VCF" ]; then
        log_message "WARNING: SnpEff annotated VCF not found: $SNPEFF_VCF"
        log_message "Skipping SnpEff visualization"
    else
        log_message "Generating SnpEff annotation figures..."
        
        python3 "${SCRIPT_DIR}/phase2_step4.1_visualize_snpeff.py" \
            --input-vcf "$SNPEFF_VCF" \
            --output-dir "$SNPEFF_OUTPUT" \
            --dpi "$DPI" \
            --format "$FORMAT" \
            --top-effects 20 \
            --top-genes 30
        
        log_message "SnpEff visualization complete"
    fi
    
    log_message ""
fi

# ============================================================================
# FUNCTIONAL ANNOTATION VISUALIZATION
# ============================================================================

if [ "$SNPEFF_ONLY" = false ]; then
    log_message "====================================================================="
    log_message "Functional Annotation Visualization"
    log_message "====================================================================="
    log_message ""
    
    if [ ! -d "$FUNCTIONAL_DIR" ]; then
        log_message "WARNING: Functional annotation directory not found: $FUNCTIONAL_DIR"
        log_message "Skipping functional annotation visualization"
    else
        log_message "Generating functional annotation figures..."
        
        # Generate plots for all gene types
        python3 "${SCRIPT_DIR}/phase2_step4.2_visualize_functional.py" \
            --input-dir "$FUNCTIONAL_DIR" \
            --output-dir "$FUNCTIONAL_DIR" \
            --dpi "$DPI" \
            --format "$FORMAT" \
            --all-types \
            --comparison \
            --show-counts
        
        log_message "Functional annotation visualization complete"
    fi
    
    log_message ""
fi

# ============================================================================
# COMPLETION
# ============================================================================

log_message "====================================================================="
log_message "VISUALIZATION COMPLETE"
log_message "====================================================================="
log_message ""

if [ "$FUNCTIONAL_ONLY" = false ]; then
    log_message "SnpEff figures saved to: $SNPEFF_OUTPUT"
    if [ -d "$SNPEFF_OUTPUT" ]; then
        ls -lh "$SNPEFF_OUTPUT"/*."$FORMAT" 2>/dev/null | while read line; do
            log_message "  $(basename $(echo $line | awk '{print $9}'))"
        done
    fi
    log_message ""
fi

if [ "$SNPEFF_ONLY" = false ]; then
    log_message "Functional annotation figures saved to: $FUNCTIONAL_DIR"
    if [ -d "$FUNCTIONAL_DIR" ]; then
        ls -lh "$FUNCTIONAL_DIR"/*."$FORMAT" 2>/dev/null | while read line; do
            log_message "  $(basename $(echo $line | awk '{print $9}'))"
        done
    fi
    log_message ""
fi

log_message "To regenerate specific figures, use:"
log_message "  python3 ${SCRIPT_DIR}/phase2_step4.1_visualize_snpeff.py --help"
log_message "  python3 ${SCRIPT_DIR}/phase2_step4.2_visualize_functional.py --help"
log_message ""

exit 0

