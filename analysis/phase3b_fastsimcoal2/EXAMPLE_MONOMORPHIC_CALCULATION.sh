#!/bin/bash
################################################################################
# Example: Calculate Monomorphic Sites for fastsimcoal2
################################################################################
#
# This script demonstrates how to calculate monomorphic sites count
# using the known reference genome size.
#
# Reference genome: 2,948,446,826 bp (21 autosomes)
#
################################################################################

set -e

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/load_base_dir.sh"

echo "=============================================================================="
echo "MONOMORPHIC SITES CALCULATION - EXAMPLE"
echo "=============================================================================="
echo ""

# Step 1: Define reference genome size
REF_GENOME=2948446826
CALLABLE_PROP=0.85

echo "Step 1: Reference genome information"
echo "  Total genome size: $(printf "%'d" $REF_GENOME) bp"
echo "  Callable proportion: ${CALLABLE_PROP} (85%)"
echo ""

# Step 2: Calculate callable regions
CALLABLE=$(echo "$REF_GENOME * $CALLABLE_PROP" | bc | cut -d'.' -f1)
echo "Step 2: Estimated callable genome"
echo "  Callable regions: $(printf "%'d" $CALLABLE) bp"
echo ""

# Step 3: Count variants from VCF
echo "Step 3: Counting variants from VCF..."
VCF_FILE="${BASE_DIR}/data/monkey_snp_autosomes_only.vcf.gz"

if [ -f "$VCF_FILE" ]; then
    VARIANTS=$(zcat "$VCF_FILE" | grep -v "^#" | wc -l)
    echo "  Total variants: $(printf "%'d" $VARIANTS)"
    echo ""
    
    # Step 4: Calculate monomorphic sites
    MONOMORPHIC=$(($CALLABLE - $VARIANTS))
    echo "Step 4: Calculate monomorphic sites"
    echo "  Formula: Callable genome - Total variants"
    echo "  Monomorphic sites: $(printf "%'d" $MONOMORPHIC)"
    echo ""
    
    # Step 5: Show how to update SFS
    echo "Step 5: Update SFS file (optional)"
    echo "  Run the following command to update:"
    echo ""
    echo "  python3 -c \""
    echo "  from phase3b_step0b_estimate_monomorphic_sites import update_sfs_file"
    echo "  from pathlib import Path"
    echo "  sfs = Path('output/phase3b_fastsimcoal2/sfs/SNJ_DAFpop0.obs')"
    echo "  update_sfs_file(sfs, $MONOMORPHIC)"
    echo "  \""
    echo ""
else
    echo "  VCF file not found: $VCF_FILE"
    echo "  Using example variant count: 1,500,000"
    VARIANTS=1500000
    MONOMORPHIC=$(($CALLABLE - $VARIANTS))
    echo ""
    echo "Step 4: Example calculation"
    echo "  Variants: $(printf "%'d" $VARIANTS)"
    echo "  Monomorphic: $(printf "%'d" $MONOMORPHIC)"
    echo ""
fi

echo "=============================================================================="
echo "RECOMMENDED WORKFLOW"
echo "=============================================================================="
echo ""
echo "For actual analysis, use the automated tool:"
echo ""
echo "  cd analysis/phase3b_fastsimcoal2"
echo "  python3 phase3b_step0b_estimate_monomorphic_sites.py"
echo ""
echo "This will:"
echo "  - Use known genome size automatically"
echo "  - Count variants from your VCF"
echo "  - Calculate optimal monomorphic count"
echo "  - Offer to update SFS file"
echo ""

exit 0
