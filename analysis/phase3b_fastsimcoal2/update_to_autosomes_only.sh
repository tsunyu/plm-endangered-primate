#!/bin/bash
#
# Automatic Update: Switch to Autosomes-Only VCF
# ===============================================
#
# This script automatically updates all configuration to use
# the autosomes-only VCF (chromosomes 1-21).
#
# Usage:
#   bash update_to_autosomes_only.sh
#

set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/load_base_dir.sh"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "================================================================================"
echo "UPDATING CONFIGURATION TO USE AUTOSOMES ONLY"
echo "================================================================================"
echo ""
echo "This will update:"
echo "  1. phase3b_step0_prepare_sfs.py"
echo "  2. phase3b_step0b_estimate_monomorphic_sites.py"
echo "  3. check_vcf_chromosomes.py"
echo ""
echo "Changing from:"
echo "  data/monkey_snp_sex_qc.vcf.gz (22 chromosomes)"
echo ""
echo "To:"
echo "  data/monkey_snp_autosomes_only.vcf.gz (21 chromosomes)"
echo ""
echo "================================================================================"
echo ""

# Check if autosomes-only VCF exists
AUTOSOMES_VCF="${BASE_DIR}/data/monkey_snp_autosomes_only.vcf.gz"

if [ ! -f "${AUTOSOMES_VCF}" ]; then
    echo "⚠️  WARNING: Autosomes-only VCF not found!"
    echo ""
    echo "Please run the filtering script first:"
    echo "  bash filter_autosomes_only.sh"
    echo ""
    exit 1
fi

echo "✅ Found autosomes-only VCF: ${AUTOSOMES_VCF}"
echo ""

# Backup original files
echo "Creating backups..."
cp "${SCRIPT_DIR}/phase3b_step0_prepare_sfs.py" "${SCRIPT_DIR}/phase3b_step0_prepare_sfs.py.backup_with_sexchr"
cp "${SCRIPT_DIR}/phase3b_step0b_estimate_monomorphic_sites.py" "${SCRIPT_DIR}/phase3b_step0b_estimate_monomorphic_sites.py.backup_with_sexchr"
cp "${SCRIPT_DIR}/check_vcf_chromosomes.py" "${SCRIPT_DIR}/check_vcf_chromosomes.py.backup_with_sexchr"
echo "✅ Backups created (.backup_with_sexchr)"
echo ""

# Update phase3b_step0_prepare_sfs.py
echo "Updating phase3b_step0_prepare_sfs.py..."
sed -i 's|"data/monkey_snp_sex_qc.vcf.gz"|"data/monkey_snp_autosomes_only.vcf.gz"|g' \
    "${SCRIPT_DIR}/phase3b_step0_prepare_sfs.py"
echo "✅ Updated"

# Update phase3b_step0b_estimate_monomorphic_sites.py
echo "Updating phase3b_step0b_estimate_monomorphic_sites.py..."
sed -i 's|"data/monkey_snp_sex_qc.vcf.gz"|"data/monkey_snp_autosomes_only.vcf.gz"|g' \
    "${SCRIPT_DIR}/phase3b_step0b_estimate_monomorphic_sites.py"
echo "✅ Updated"

# Update check_vcf_chromosomes.py
echo "Updating check_vcf_chromosomes.py..."
sed -i 's|"data/monkey_snp_sex_qc.vcf.gz"|"data/monkey_snp_autosomes_only.vcf.gz"|g' \
    "${SCRIPT_DIR}/check_vcf_chromosomes.py"
echo "✅ Updated"

echo ""
echo "================================================================================"
echo "UPDATE COMPLETE"
echo "================================================================================"
echo ""
echo "Changes made:"
echo "  ✅ All scripts now use: data/monkey_snp_autosomes_only.vcf.gz"
echo "  ✅ Backups saved with .backup_with_sexchr extension"
echo ""
echo "To verify the update:"
echo "  cd ${SCRIPT_DIR}"
echo "  python3 check_vcf_chromosomes.py"
echo ""
echo "Expected output:"
echo "  - Total chromosomes: 21"
echo "  - Reference genome: 2,948,446,826 bp (21 autosomes)"
echo "  - Monomorphic sites: ~2,501,617,687"
echo ""
echo "To restore original configuration:"
echo "  mv phase3b_step0_prepare_sfs.py.backup_with_sexchr phase3b_step0_prepare_sfs.py"
echo "  mv phase3b_step0b_estimate_monomorphic_sites.py.backup_with_sexchr phase3b_step0b_estimate_monomorphic_sites.py"
echo "  mv check_vcf_chromosomes.py.backup_with_sexchr check_vcf_chromosomes.py"
echo ""
echo "================================================================================"
echo ""
echo "✅ Ready to run analysis with AUTOSOMES ONLY!"
echo ""
