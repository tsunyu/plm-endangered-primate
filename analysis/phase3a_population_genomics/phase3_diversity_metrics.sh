#!/bin/bash
#
# Phase 3a: Genetic Diversity Metrics
#
# Calculates:
# - Individual heterozygosity
# - Nucleotide diversity (π) in 100kb windows
# - Tajima's D
#
# Usage: bash phase3_diversity_metrics.sh
#

set -euo pipefail

# ============================================================================
# SETUP
# ============================================================================

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Paths
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/load_base_dir.sh"
DATA_DIR="${BASE_DIR}/data"
OUTPUT_DIR="${BASE_DIR}/output/phase3a_population_genomics/diversity_metrics"
PLINK_PREFIX="${DATA_DIR}/monkey_snp_sex_qc"
VCF="${DATA_DIR}/monkey_snp_sex_qc.vcf"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Log file
LOGFILE="${OUTPUT_DIR}/diversity_metrics.log"
exec > >(tee -a "$LOGFILE") 2>&1

# ============================================================================
# PARAMETERS
# ============================================================================

# Window parameters
WINDOW_SIZE=100000       # 100 kb
WINDOW_STEP=50000        # 50 kb step

# Threads
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
log_message "GENETIC DIVERSITY METRICS"
log_message "====================================================================="
log_message ""

# ============================================================================
# Step 1: Individual Heterozygosity
# ============================================================================

log_message "Step 1: Calculating individual heterozygosity..."

plink \
    --bfile "${PLINK_PREFIX}" \
    --het \
    --threads ${THREADS} \
    --out "${OUTPUT_DIR}/heterozygosity"

if [ $? -eq 0 ]; then
    log_message "Heterozygosity calculation completed"
    
    # Parse and summarize results
    python3 << 'EOF'
import os
import pandas as pd
import numpy as np

# Read heterozygosity file
het_file = os.environ["PLM_BASE_DIR"] + "/output/phase3a_population_genomics/diversity_metrics/heterozygosity.het"
het_df = pd.read_csv(het_file, sep=r'\s+')

# Calculate observed heterozygosity
het_df['OBS_HET'] = (het_df['N(NM)'] - het_df['O(HOM)']) / het_df['N(NM)']

# Calculate F coefficient (inbreeding coefficient from heterozygosity)
het_df['F'] = het_df['F']

# Sort by heterozygosity
het_df_sorted = het_df.sort_values('OBS_HET', ascending=False)

# Save results
output_file = os.environ["PLM_BASE_DIR"] + "/output/phase3a_population_genomics/diversity_metrics/heterozygosity_summary.csv"
het_df_sorted.to_csv(output_file, index=False)

# Summary statistics
print("\n" + "="*70)
print("HETEROZYGOSITY SUMMARY")
print("="*70)
print(f"Number of individuals: {len(het_df)}")
print(f"\nObserved Heterozygosity:")
print(f"  Mean: {het_df['OBS_HET'].mean():.6f}")
print(f"  SD: {het_df['OBS_HET'].std():.6f}")
print(f"  Min: {het_df['OBS_HET'].min():.6f}")
print(f"  Max: {het_df['OBS_HET'].max():.6f}")
print(f"  Median: {het_df['OBS_HET'].median():.6f}")
print(f"\nInbreeding Coefficient (F):")
print(f"  Mean: {het_df['F'].mean():.6f}")
print(f"  SD: {het_df['F'].std():.6f}")
print(f"  Min: {het_df['F'].min():.6f}")
print(f"  Max: {het_df['F'].max():.6f}")

# Top 5 most heterozygous
print(f"\nTop 5 most heterozygous individuals:")
for idx, row in het_df_sorted.head(5).iterrows():
    print(f"  {row['IID']}: {row['OBS_HET']:.6f}")

# Top 5 least heterozygous
print(f"\nTop 5 least heterozygous individuals:")
for idx, row in het_df_sorted.tail(5).iterrows():
    print(f"  {row['IID']}: {row['OBS_HET']:.6f}")

print("="*70)
EOF

else
    log_message "ERROR: Heterozygosity calculation failed"
    exit 1
fi

log_message ""

# ============================================================================
# Step 2: Nucleotide Diversity (π) in Windows
# ============================================================================

log_message "Step 2: Calculating nucleotide diversity (π) in ${WINDOW_SIZE} bp windows..."

# Prepare VCF if needed
if [[ ! "$VCF" =~ \.gz$ ]]; then
    log_message "  Compressing VCF..."
    bgzip -c "$VCF" > "${VCF}.gz"
    VCF="${VCF}.gz"
fi

if [ ! -f "${VCF}.tbi" ]; then
    log_message "  Indexing VCF..."
    tabix -p vcf "$VCF"
fi

# Calculate nucleotide diversity with VCFtools
vcftools \
    --gzvcf "$VCF" \
    --window-pi ${WINDOW_SIZE} \
    --window-pi-step ${WINDOW_STEP} \
    --out "${OUTPUT_DIR}/nucleotide_diversity"

if [ $? -eq 0 ]; then
    log_message "Nucleotide diversity calculation completed"
    
    # Summarize results
    python3 << 'EOF'
import os
import pandas as pd
import numpy as np

# Read pi file
pi_file = os.environ["PLM_BASE_DIR"] + "/output/phase3a_population_genomics/diversity_metrics/nucleotide_diversity.windowed.pi"
pi_df = pd.read_csv(pi_file, sep='\t')

# Remove windows with insufficient data
pi_df_filtered = pi_df[pi_df['N_VARIANTS'] > 0].copy()

print("\n" + "="*70)
print("NUCLEOTIDE DIVERSITY (π) SUMMARY")
print("="*70)
print(f"Total windows: {len(pi_df)}")
print(f"Windows with variants: {len(pi_df_filtered)}")
print(f"\nNucleotide Diversity (π):")
print(f"  Mean: {pi_df_filtered['PI'].mean():.8f}")
print(f"  SD: {pi_df_filtered['PI'].std():.8f}")
print(f"  Min: {pi_df_filtered['PI'].min():.8f}")
print(f"  Max: {pi_df_filtered['PI'].max():.8f}")
print(f"  Median: {pi_df_filtered['PI'].median():.8f}")

# Per-chromosome statistics
print("\nPer-chromosome π:")
for chrom in sorted(pi_df_filtered['CHROM'].unique()):
    chrom_data = pi_df_filtered[pi_df_filtered['CHROM'] == chrom]
    print(f"  Chr {chrom}: {chrom_data['PI'].mean():.8f} (n={len(chrom_data)} windows)")

# Save summary
output_file = os.environ["PLM_BASE_DIR"] + "/output/phase3a_population_genomics/diversity_metrics/nucleotide_diversity_summary.csv"
pi_df_filtered.to_csv(output_file, index=False)

print("="*70)
EOF

else
    log_message "ERROR: Nucleotide diversity calculation failed"
    exit 1
fi

log_message ""

# ============================================================================
# Step 3: Tajima's D
# ============================================================================

log_message "Step 3: Calculating Tajima's D in ${WINDOW_SIZE} bp windows..."

vcftools \
    --gzvcf "$VCF" \
    --TajimaD ${WINDOW_SIZE} \
    --out "${OUTPUT_DIR}/tajimas_d"

if [ $? -eq 0 ]; then
    log_message "Tajima's D calculation completed"
    
    # Summarize results
    python3 << 'EOF'
import os
import pandas as pd
import numpy as np

# Read Tajima's D file
td_file = os.environ["PLM_BASE_DIR"] + "/output/phase3a_population_genomics/diversity_metrics/tajimas_d.Tajima.D"
td_df = pd.read_csv(td_file, sep='\t')

# Remove windows with NaN
td_df_filtered = td_df[~td_df['TajimaD'].isna()].copy()

print("\n" + "="*70)
print("TAJIMA'S D SUMMARY")
print("="*70)
print(f"Total windows: {len(td_df)}")
print(f"Windows with valid Tajima's D: {len(td_df_filtered)}")
print(f"\nTajima's D:")
print(f"  Mean: {td_df_filtered['TajimaD'].mean():.6f}")
print(f"  SD: {td_df_filtered['TajimaD'].std():.6f}")
print(f"  Min: {td_df_filtered['TajimaD'].min():.6f}")
print(f"  Max: {td_df_filtered['TajimaD'].max():.6f}")
print(f"  Median: {td_df_filtered['TajimaD'].median():.6f}")

# Count windows in different ranges
negative = len(td_df_filtered[td_df_filtered['TajimaD'] < -2])
neutral = len(td_df_filtered[(td_df_filtered['TajimaD'] >= -2) & (td_df_filtered['TajimaD'] <= 2)])
positive = len(td_df_filtered[td_df_filtered['TajimaD'] > 2])

print(f"\nWindows by Tajima's D range:")
print(f"  Negative (D < -2): {negative} ({100*negative/len(td_df_filtered):.1f}%)")
print(f"  Neutral (-2 ≤ D ≤ 2): {neutral} ({100*neutral/len(td_df_filtered):.1f}%)")
print(f"  Positive (D > 2): {positive} ({100*positive/len(td_df_filtered):.1f}%)")

# Interpretation
mean_d = td_df_filtered['TajimaD'].mean()
print(f"\nInterpretation:")
if mean_d < -1:
    print("  Mean Tajima's D < -1: Suggests population expansion or positive selection")
elif mean_d > 1:
    print("  Mean Tajima's D > 1: Suggests population bottleneck or balancing selection")
else:
    print("  Mean Tajima's D near 0: Consistent with neutral evolution")

# Save summary
output_file = os.environ["PLM_BASE_DIR"] + "/output/phase3a_population_genomics/diversity_metrics/tajimas_d_summary.csv"
td_df_filtered.to_csv(output_file, index=False)

print("="*70)
EOF

else
    log_message "ERROR: Tajima's D calculation failed"
    exit 1
fi

log_message ""

# ============================================================================
# Step 4: Generate Visualizations
# ============================================================================

log_message "Step 4: Generating diversity visualizations..."

cat > "${OUTPUT_DIR}/plot_diversity.R" << 'EOF'
#!/usr/bin/env Rscript
# Diversity Metrics Visualization

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

source(file.path(Sys.getenv("PLM_ANALYSIS_DIR"), "utils.R"))

output_dir <- file.path(Sys.getenv("PLM_BASE_DIR"), "output/phase3a_population_genomics/diversity_metrics")

# Read data
het_df <- read.csv(file.path(output_dir, "heterozygosity_summary.csv"))
pi_df <- read.csv(file.path(output_dir, "nucleotide_diversity_summary.csv"))
td_df <- read.csv(file.path(output_dir, "tajimas_d_summary.csv"))

# 1. Heterozygosity distribution
p1 <- ggplot(het_df, aes(x = OBS_HET)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(aes(xintercept = mean(OBS_HET)), color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Distribution of Observed Heterozygosity",
    x = "Observed Heterozygosity",
    y = "Count"
  ) +
  theme_publication()

ggsave(file.path(output_dir, "heterozygosity_distribution.png"), p1, width = 8, height = 6, dpi = 300)

# 2. Nucleotide diversity genome-wide
p2 <- plot_manhattan(pi_df, "PI", output_file = file.path(output_dir, "pi_manhattan.png"))

# 3. Tajima's D genome-wide
p3 <- plot_manhattan(td_df, "TajimaD", output_file = file.path(output_dir, "tajima_d_manhattan.png"))

# 4. Tajima's D distribution
p4 <- ggplot(td_df, aes(x = TajimaD)) +
  geom_histogram(bins = 50, fill = "coral", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0, color = "black", linetype = "solid", linewidth = 1) +
  geom_vline(aes(xintercept = mean(TajimaD)), color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = c(-2, 2), color = "blue", linetype = "dotted", linewidth = 0.8) +
  labs(
    title = "Distribution of Tajima's D",
    x = "Tajima's D",
    y = "Count"
  ) +
  theme_publication()

ggsave(file.path(output_dir, "tajima_d_distribution.png"), p4, width = 8, height = 6, dpi = 300)

# 5. Combined summary plot
p5 <- ggplot(pi_df, aes(x = PI)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  labs(title = "Nucleotide Diversity (π)", x = "π", y = "Count") +
  theme_publication()

combined <- p1 + p5 + p4 + plot_layout(ncol = 2)
ggsave(file.path(output_dir, "diversity_summary.png"), combined, width = 12, height = 10, dpi = 300)

cat("Diversity visualizations generated successfully\n")
EOF

Rscript "${OUTPUT_DIR}/plot_diversity.R"

if [ $? -eq 0 ]; then
    log_message "Visualizations generated successfully"
else
    log_message "WARNING: Some visualizations may have failed"
fi

# ============================================================================
# COMPLETION
# ============================================================================

log_message ""
log_message "====================================================================="
log_message "GENETIC DIVERSITY METRICS ANALYSIS COMPLETE"
log_message "====================================================================="
log_message "Output directory: ${OUTPUT_DIR}"
log_message ""
log_message "Key output files:"
log_message "  - heterozygosity_summary.csv         : Individual heterozygosity"
log_message "  - nucleotide_diversity_summary.csv   : π in 100kb windows"
log_message "  - tajimas_d_summary.csv              : Tajima's D in 100kb windows"
log_message "  - heterozygosity_distribution.png    : Het distribution plot"
log_message "  - pi_manhattan.png                   : Genome-wide π"
log_message "  - tajima_d_manhattan.png             : Genome-wide Tajima's D"
log_message "====================================================================="

exit 0


