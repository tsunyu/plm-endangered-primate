#!/bin/bash
#
# Filter VCF to Include Only Autosomes (Chromosomes 1-21)
# ========================================================
#
# fastsimcoal2 standard practice: Use AUTOSOMES ONLY
# This script excludes sex chromosomes from the analysis.
#
# Usage:
#   bash filter_autosomes_only.sh
#
# Input:  data/monkey_snp_sex_qc.vcf.gz (22 chromosomes)
# Output: data/monkey_snp_autosomes_only.vcf.gz (21 chromosomes)
#

set -euo pipefail

# Paths
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/load_base_dir.sh"
INPUT_VCF="${BASE_DIR}/data/monkey_snp_sex_qc.vcf.gz"
OUTPUT_VCF="${BASE_DIR}/data/monkey_snp_autosomes_only.vcf.gz"
TEMP_VCF="${BASE_DIR}/data/monkey_snp_autosomes_only.temp.vcf"

echo "================================================================================"
echo "FILTERING VCF TO AUTOSOMES ONLY"
echo "================================================================================"
echo ""
echo "Standard practice for fastsimcoal2: Use AUTOSOMES ONLY"
echo ""
echo "Input:  ${INPUT_VCF}"
echo "Output: ${OUTPUT_VCF}"
echo ""
echo "Chromosomes to KEEP: 1-21 (autosomes)"
echo "Chromosomes to EXCLUDE: 22 (sex chromosome)"
echo ""
echo "================================================================================"
echo ""

# Check if input exists
if [ ! -f "${INPUT_VCF}" ]; then
    echo "ERROR: Input VCF not found: ${INPUT_VCF}"
    exit 1
fi

# Check if bcftools is available
if ! command -v bcftools &> /dev/null; then
    echo "ERROR: bcftools is not installed"
    echo ""
    echo "Install with:"
    echo "  conda install -c bioconda bcftools"
    echo "  # or"
    echo "  apt-get install bcftools"
    exit 1
fi

echo "Step 1: Extracting autosomes (chromosomes 1-21)..."
echo ""

# Method 1: Using bcftools view (preferred)
bcftools view \
    --regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21 \
    "${INPUT_VCF}" \
    -O z \
    -o "${OUTPUT_VCF}"

echo "Step 2: Indexing output VCF..."
bcftools index "${OUTPUT_VCF}"

echo ""
echo "================================================================================"
echo "FILTERING COMPLETE"
echo "================================================================================"
echo ""

# Count variants
echo "Variant Counts:"
echo "---------------"
echo ""

TOTAL_BEFORE=$(zcat "${INPUT_VCF}" | grep -v "^#" | wc -l)
TOTAL_AFTER=$(zcat "${OUTPUT_VCF}" | grep -v "^#" | wc -l)
EXCLUDED=$((TOTAL_BEFORE - TOTAL_AFTER))

echo "Before filtering (22 chr): $(printf "%'d" ${TOTAL_BEFORE}) variants"
echo "After filtering (21 chr):  $(printf "%'d" ${TOTAL_AFTER}) variants"
echo "Excluded (chr 22):         $(printf "%'d" ${EXCLUDED}) variants"
echo ""

# Check chromosomes
echo "Chromosomes in filtered VCF:"
zcat "${OUTPUT_VCF}" | grep -v "^##" | cut -f1 | grep -v "^#CHROM" | sort -u | tr '\n' ' '
echo ""
echo ""

echo "================================================================================"
echo "NEXT STEPS"
echo "================================================================================"
echo ""
echo "1. Update configuration to use autosome-only VCF:"
echo "   Edit: analysis/phase3b_fastsimcoal2/phase3b_step0_prepare_sfs.py"
echo ""
echo "   Change:"
echo "     MAIN_VCF = BASE_DIR / \"data/monkey_snp_sex_qc.vcf.gz\""
echo "   To:"
echo "     MAIN_VCF = BASE_DIR / \"data/monkey_snp_autosomes_only.vcf.gz\""
echo ""
echo "2. Update reference genome size:"
echo "   Edit: analysis/phase3b_fastsimcoal2/phase3b_step0b_estimate_monomorphic_sites.py"
echo ""
echo "   The auto-detection will now detect 21 chromosomes and use:"
echo "     Reference: 2,948,446,826 bp (21 autosomes)"
echo ""
echo "3. Re-run the pipeline:"
echo "   cd analysis/phase3b_fastsimcoal2"
echo "   python3 phase3b_step0_prepare_sfs.py"
echo "   python3 phase3b_step0b_estimate_monomorphic_sites.py"
echo "   bash run_complete_pipeline.sh"
echo ""
echo "================================================================================"
echo ""
echo "✅ Filtering complete!"
echo ""
