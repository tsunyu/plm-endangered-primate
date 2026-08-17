#!/bin/bash
#
# Phase 3a: Population Structure Analysis
#
# Performs:
# 1. Principal Component Analysis (PCA)
# 2. Kinship analysis (KING)
# 3. Identity-by-descent (IBD) segment detection
#
# Usage: bash phase3a_population_structure.sh
#

set -euo pipefail

# ============================================================================
# SETUP
# ============================================================================

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/load_base_dir.sh"
DATA_DIR="${BASE_DIR}/data"
OUTPUT_DIR="${BASE_DIR}/output/phase3a_population_genomics/population_structure"
PLINK_PREFIX="${DATA_DIR}/monkey_snp_sex_qc"

mkdir -p "${OUTPUT_DIR}"

LOGFILE="${OUTPUT_DIR}/population_structure.log"
exec > >(tee -a "$LOGFILE") 2>&1

# ============================================================================
# PARAMETERS
# ============================================================================

NUM_PCS=10
KINSHIP_THRESHOLD=0.0884  # 1st/2nd degree relatives
THREADS=8

# ============================================================================
# FUNCTIONS
# ============================================================================

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# ============================================================================
# ANALYSIS
# ============================================================================

log_message "====================================================================="
log_message "POPULATION STRUCTURE ANALYSIS"
log_message "====================================================================="
log_message ""

# ============================================================================
# Step 1: Principal Component Analysis
# ============================================================================

log_message "Step 1: Running Principal Component Analysis..."

# Prune SNPs for LD
log_message "  Pruning SNPs for LD..."
plink \
    --bfile "${PLINK_PREFIX}" \
    --indep-pairwise 50 5 0.2 \
    --threads ${THREADS} \
    --out "${OUTPUT_DIR}/pca_pruned"

# Run PCA
log_message "  Running PCA..."
plink \
    --bfile "${PLINK_PREFIX}" \
    --extract "${OUTPUT_DIR}/pca_pruned.prune.in" \
    --pca ${NUM_PCS} \
    --threads ${THREADS} \
    --out "${OUTPUT_DIR}/pca"

if [ $? -eq 0 ]; then
    log_message "PCA completed successfully"
    
    # Analyze PCA results
    python3 << 'EOF'
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

output_dir = os.environ["PLM_BASE_DIR"] + "/output/phase3a_population_genomics/population_structure"

# Read eigenvalues
eigenval_file = f"{output_dir}/pca.eigenval"
eigenvals = np.loadtxt(eigenval_file)

# Calculate variance explained
var_explained = 100 * eigenvals / eigenvals.sum()

# Read eigenvectors (PCs)
eigenvec_file = f"{output_dir}/pca.eigenvec"
pca_df = pd.read_csv(eigenvec_file, sep=r'\s+', header=None)
pca_df.columns = ['FID', 'IID'] + [f'PC{i}' for i in range(1, len(pca_df.columns)-1)]

# Save with variance explained
pca_summary = pca_df.copy()
pca_summary.to_csv(f"{output_dir}/pca_results.csv", index=False)

print("\n" + "="*70)
print("PRINCIPAL COMPONENT ANALYSIS SUMMARY")
print("="*70)
print(f"\nVariance explained by top 10 PCs:")
for i, var in enumerate(var_explained[:10], 1):
    print(f"  PC{i}: {var:.2f}%")

print(f"\nCumulative variance (PC1-PC10): {var_explained[:10].sum():.2f}%")

# Plot variance explained
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.bar(range(1, 11), var_explained[:10], color='steelblue', alpha=0.7, edgecolor='black')
ax1.set_xlabel('Principal Component')
ax1.set_ylabel('Variance Explained (%)')
ax1.set_title('Variance Explained by PCs')
ax1.set_xticks(range(1, 11))

# PC1 vs PC2 plot
ax2.scatter(pca_df['PC1'], pca_df['PC2'], s=50, alpha=0.6, color='coral', edgecolors='black')
ax2.set_xlabel(f'PC1 ({var_explained[0]:.2f}%)')
ax2.set_ylabel(f'PC2 ({var_explained[1]:.2f}%)')
ax2.set_title('PCA: PC1 vs PC2')
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(f"{output_dir}/pca_plots.png", dpi=300, bbox_inches='tight')
plt.close()

print("\nPCA plots saved")
print("="*70)
EOF

else
    log_message "ERROR: PCA failed"
    exit 1
fi

log_message ""

# ============================================================================
# Step 2: KING Kinship Analysis
# ============================================================================

log_message "Step 2: Running KING kinship analysis..."

# Run KING
king \
    -b "${PLINK_PREFIX}.bed" \
    --kinship \
    --prefix "${OUTPUT_DIR}/king"

if [ $? -eq 0 ]; then
    log_message "KING kinship analysis completed"
    
    # Analyze kinship results
    python3 << 'EOF'
import os
import pandas as pd
import numpy as np

output_dir = os.environ["PLM_BASE_DIR"] + "/output/phase3a_population_genomics/population_structure"

# Read kinship file
kinship_file = f"{output_dir}/king.kin0"

try:
    kin_df = pd.read_csv(kinship_file, sep='\t')
    
    print("\n" + "="*70)
    print("KING KINSHIP ANALYSIS SUMMARY")
    print("="*70)
    print(f"\nTotal pairs analyzed: {len(kin_df)}")
    
    # Classify relationships
    # Kinship > 0.354: Duplicate/MZ twins
    # 0.177 < kinship ≤ 0.354: 1st-degree (parent-offspring, full siblings)
    # 0.0884 < kinship ≤ 0.177: 2nd-degree (half-siblings, grandparent-grandchild)
    # 0.0442 < kinship ≤ 0.0884: 3rd-degree (first cousins)
    
    duplicates = kin_df[kin_df['Kinship'] > 0.354]
    first_degree = kin_df[(kin_df['Kinship'] > 0.177) & (kin_df['Kinship'] <= 0.354)]
    second_degree = kin_df[(kin_df['Kinship'] > 0.0884) & (kin_df['Kinship'] <= 0.177)]
    third_degree = kin_df[(kin_df['Kinship'] > 0.0442) & (kin_df['Kinship'] <= 0.0884)]
    
    print(f"\nRelationship classification:")
    print(f"  Duplicates/MZ twins (kinship > 0.354): {len(duplicates)}")
    print(f"  1st-degree relatives (0.177-0.354): {len(first_degree)}")
    print(f"  2nd-degree relatives (0.0884-0.177): {len(second_degree)}")
    print(f"  3rd-degree relatives (0.0442-0.0884): {len(third_degree)}")
    
    # Close relatives to consider for removal (1st/2nd degree)
    close_relatives = kin_df[kin_df['Kinship'] > 0.0884]
    
    print(f"\nClose relatives (1st/2nd degree) to consider: {len(close_relatives)}")
    
    if len(close_relatives) > 0:
        print("\nTop 10 most related pairs:")
        for idx, row in close_relatives.nlargest(10, 'Kinship').iterrows():
            print(f"  {row['ID1']} - {row['ID2']}: kinship = {row['Kinship']:.4f}")
    
    # Save relatedness summary
    close_relatives.to_csv(f"{output_dir}/close_relatives.csv", index=False)
    
    print("\n" + "="*70)
    
except Exception as e:
    print(f"Error processing kinship file: {e}")
EOF

else
    log_message "WARNING: KING analysis failed"
fi

log_message ""

# ============================================================================
# Step 3: IBD Segment Detection
# ============================================================================

log_message "Step 3: Detecting identity-by-descent (IBD) segments..."

# Use PLINK for IBD detection
plink \
    --bfile "${PLINK_PREFIX}" \
    --genome \
    --min 0.05 \
    --threads ${THREADS} \
    --out "${OUTPUT_DIR}/ibd"

if [ $? -eq 0 ]; then
    log_message "IBD analysis completed"
    
    python3 << 'EOF'
import os
import pandas as pd
import numpy as np

output_dir = os.environ["PLM_BASE_DIR"] + "/output/phase3a_population_genomics/population_structure"

# Read IBD file
ibd_file = f"{output_dir}/ibd.genome"

try:
    ibd_df = pd.read_csv(ibd_file, sep=r'\s+')
    
    print("\n" + "="*70)
    print("IDENTITY-BY-DESCENT (IBD) SUMMARY")
    print("="*70)
    print(f"\nTotal pairs analyzed: {len(ibd_df)}")
    
    # PI_HAT is proportion of genome shared IBD
    print(f"\nProportion IBD (PI_HAT) statistics:")
    print(f"  Mean: {ibd_df['PI_HAT'].mean():.4f}")
    print(f"  Median: {ibd_df['PI_HAT'].median():.4f}")
    print(f"  Min: {ibd_df['PI_HAT'].min():.4f}")
    print(f"  Max: {ibd_df['PI_HAT'].max():.4f}")
    
    # Pairs with high IBD sharing
    high_ibd = ibd_df[ibd_df['PI_HAT'] > 0.1]
    print(f"\nPairs with >10% genome shared IBD: {len(high_ibd)}")
    
    if len(high_ibd) > 0:
        print("\nTop 10 pairs with highest IBD sharing:")
        for idx, row in high_ibd.nlargest(10, 'PI_HAT').iterrows():
            print(f"  {row['IID1']} - {row['IID2']}: PI_HAT = {row['PI_HAT']:.4f}")
    
    # Z scores for relatedness inference
    print(f"\nZ-score statistics:")
    print(f"  Z0 (0 alleles IBD): {ibd_df['Z0'].mean():.4f}")
    print(f"  Z1 (1 allele IBD): {ibd_df['Z1'].mean():.4f}")
    print(f"  Z2 (2 alleles IBD): {ibd_df['Z2'].mean():.4f}")
    
    # Save summary
    high_ibd.to_csv(f"{output_dir}/high_ibd_pairs.csv", index=False)
    
    print("="*70)
    
except Exception as e:
    print(f"Error processing IBD file: {e}")
EOF

else
    log_message "WARNING: IBD analysis failed"
fi

log_message ""

# ============================================================================
# COMPLETION
# ============================================================================

log_message "====================================================================="
log_message "POPULATION STRUCTURE ANALYSIS COMPLETE"
log_message "====================================================================="
log_message "Output directory: ${OUTPUT_DIR}"
log_message ""
log_message "Key output files:"
log_message "  - pca_results.csv            : PCA results with variance explained"
log_message "  - pca_plots.png              : PCA visualization"
log_message "  - king.kin0                  : Pairwise kinship coefficients"
log_message "  - close_relatives.csv        : Close relatives (1st/2nd degree)"
log_message "  - ibd.genome                 : IBD sharing estimates"
log_message "  - high_ibd_pairs.csv         : Pairs with high IBD sharing"
log_message "====================================================================="

exit 0


