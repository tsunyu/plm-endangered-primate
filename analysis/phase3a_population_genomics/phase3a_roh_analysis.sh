#!/bin/bash
#
# Phase 3a: Runs of Homozygosity (ROH) Analysis
#
# Performs ROH detection using PLINK and BCFtools
# Calculates F_ROH and categorizes ROH by length
#
# Usage: bash phase3a_roh_analysis.sh
#

set -euo pipefail

# ============================================================================
# SETUP
# ============================================================================

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Paths
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/load_base_dir.sh"
DATA_DIR="${BASE_DIR}/data"
OUTPUT_DIR="${BASE_DIR}/output/phase3a_population_genomics/roh_analysis"
PLINK_PREFIX="${DATA_DIR}/monkey_snp_sex_qc"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Log file
LOGFILE="${OUTPUT_DIR}/roh_analysis.log"
exec > >(tee -a "$LOGFILE") 2>&1

# ============================================================================
# PARAMETERS (from proposal)
# ============================================================================

# ROH detection parameters
MIN_SNPS=30              # Minimum SNPs per ROH
MIN_LENGTH_KB=100        # Minimum ROH length (kb)
MAX_HET=1                # Maximum heterozygous calls allowed
MAX_MISSING=5            # Maximum missing calls allowed
MAX_GAP_KB=1000          # Maximum gap between consecutive SNPs (kb)
WINDOW_SNPS=50           # Sliding window size (SNPs)

# ROH length categories (bp)
SHORT_MIN=100000         # 100 kb
SHORT_MAX=1000000        # 1 Mb
LONG_MIN=1000000         # 1 Mb
LONG_MAX=5000000         # 5 Mb
VERYLONG_MIN=5000000     # 5 Mb

# Genome parameters
AUTOSOMAL_GENOME_LENGTH=2800000000  # 2.8 Gb

# Number of threads
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
log_message "ROH ANALYSIS"
log_message "====================================================================="
log_message ""
log_message "Parameters:"
log_message "  Minimum SNPs: ${MIN_SNPS}"
log_message "  Minimum length: ${MIN_LENGTH_KB} kb"
log_message "  Max heterozygous calls: ${MAX_HET}"
log_message "  Max missing calls: ${MAX_MISSING}"
log_message "  Max gap: ${MAX_GAP_KB} kb"
log_message "  Window size: ${WINDOW_SNPS} SNPs"
log_message ""

# ============================================================================
# Step 1: PLINK ROH Detection
# ============================================================================

log_message "Step 1: Running PLINK ROH detection..."

plink \
    --bfile "${PLINK_PREFIX}" \
    --homozyg \
    --homozyg-snp ${MIN_SNPS} \
    --homozyg-kb ${MIN_LENGTH_KB} \
    --homozyg-window-het ${MAX_HET} \
    --homozyg-window-missing ${MAX_MISSING} \
    --homozyg-gap ${MAX_GAP_KB} \
    --homozyg-window-snp ${WINDOW_SNPS} \
    --threads ${THREADS} \
    --out "${OUTPUT_DIR}/plink_roh"

if [ $? -eq 0 ]; then
    log_message "PLINK ROH detection completed"
    
    # Count ROH segments
    if [ -f "${OUTPUT_DIR}/plink_roh.hom" ]; then
        NUM_ROH=$(tail -n +2 "${OUTPUT_DIR}/plink_roh.hom" | wc -l)
        log_message "  Total ROH segments detected: ${NUM_ROH}"
    fi
else
    log_message "ERROR: PLINK ROH detection failed"
    exit 1
fi

log_message ""

# ============================================================================
# Step 2: BCFtools RoH Detection (HMM-based validation)
# ============================================================================

log_message "Step 2: Running BCFtools RoH detection (HMM-based)..."

# BCFtools RoH requires VCF input
VCF="${DATA_DIR}/monkey_snp_sex_qc.vcf"

if [ -f "$VCF" ]; then
    # Check if VCF is compressed, compress if needed
    if [[ ! "$VCF" =~ \.gz$ ]]; then
        log_message "  Compressing VCF file..."
        bgzip -c "$VCF" > "${VCF}.gz"
        VCF="${VCF}.gz"
    fi
    
    # Index VCF if not already indexed
    if [ ! -f "${VCF}.tbi" ]; then
        log_message "  Indexing VCF file..."
        tabix -p vcf "$VCF"
    fi
    
    # Run BCFtools RoH
    bcftools roh \
        --AF-dflt 0.4 \
        --skip-indels \
        -G 30 \
        -O r \
        "$VCF" \
        > "${OUTPUT_DIR}/bcftools_roh.txt"
    
    if [ $? -eq 0 ]; then
        log_message "BCFtools RoH detection completed"
        
        # Count ROH from BCFtools (RG lines)
        NUM_BCF_ROH=$(grep "^RG" "${OUTPUT_DIR}/bcftools_roh.txt" | wc -l)
        log_message "  Total ROH segments detected: ${NUM_BCF_ROH}"
    else
        log_message "WARNING: BCFtools RoH detection failed"
    fi
else
    log_message "WARNING: VCF file not found, skipping BCFtools RoH"
fi

log_message ""

# ============================================================================
# Step 3: Calculate ROH Statistics with Python
# ============================================================================

log_message "Step 3: Calculating ROH statistics..."

# Create Python script for ROH analysis
cat > "${OUTPUT_DIR}/analyze_roh.py" << 'EOF'
#!/usr/bin/env python3
"""
Analyze PLINK ROH output and calculate statistics
"""

import sys
import pandas as pd
import numpy as np

def classify_roh_length(length_kb):
    """Classify ROH by length"""
    length_bp = length_kb * 1000
    
    if length_bp < 1e6:
        return 'short'
    elif length_bp < 5e6:
        return 'long'
    else:
        return 'verylong'

def main():
    # Parameters
    genome_length = 2.8e9  # 2.8 Gb
    output_dir = sys.argv[1]
    
    # Read PLINK ROH output
    roh_file = f"{output_dir}/plink_roh.hom"
    
    try:
        roh_df = pd.read_csv(roh_file, sep=r'\s+')
    except:
        print("ERROR: Could not read PLINK ROH file")
        sys.exit(1)
    
    print(f"Loaded {len(roh_df)} ROH segments")
    
    # Add ROH category
    roh_df['ROH_Category'] = roh_df['KB'].apply(classify_roh_length)
    roh_df['Length_BP'] = roh_df['KB'] * 1000
    
    # Calculate per-individual statistics
    individual_stats = []
    
    for iid in roh_df['IID'].unique():
        ind_roh = roh_df[roh_df['IID'] == iid]
        
        # Total ROH length and count
        total_length = ind_roh['Length_BP'].sum()
        num_roh = len(ind_roh)
        
        # F_ROH
        f_roh = total_length / genome_length
        
        # By category
        short_count = len(ind_roh[ind_roh['ROH_Category'] == 'short'])
        long_count = len(ind_roh[ind_roh['ROH_Category'] == 'long'])
        verylong_count = len(ind_roh[ind_roh['ROH_Category'] == 'verylong'])
        
        short_length = ind_roh[ind_roh['ROH_Category'] == 'short']['Length_BP'].sum()
        long_length = ind_roh[ind_roh['ROH_Category'] == 'long']['Length_BP'].sum()
        verylong_length = ind_roh[ind_roh['ROH_Category'] == 'verylong']['Length_BP'].sum()
        
        # Mean and max ROH length
        mean_length = ind_roh['Length_BP'].mean()
        max_length = ind_roh['Length_BP'].max()
        
        individual_stats.append({
            'IID': iid,
            'FID': ind_roh['FID'].iloc[0],
            'Total_ROH_Length_BP': total_length,
            'Total_ROH_Length_MB': total_length / 1e6,
            'Num_ROH': num_roh,
            'F_ROH': f_roh,
            'Mean_ROH_Length_KB': mean_length / 1000,
            'Max_ROH_Length_KB': max_length / 1000,
            'Short_ROH_Count': short_count,
            'Long_ROH_Count': long_count,
            'VeryLong_ROH_Count': verylong_count,
            'Short_ROH_Length_MB': short_length / 1e6,
            'Long_ROH_Length_MB': long_length / 1e6,
            'VeryLong_ROH_Length_MB': verylong_length / 1e6,
            'Pct_Genome_in_ROH': f_roh * 100
        })
    
    # Create DataFrame
    stats_df = pd.DataFrame(individual_stats)
    
    # Sort by F_ROH
    stats_df = stats_df.sort_values('F_ROH', ascending=False)
    
    # Save individual statistics
    stats_df.to_csv(f"{output_dir}/roh_summary_per_individual.csv", index=False)
    print(f"Individual statistics saved to: {output_dir}/roh_summary_per_individual.csv")
    
    # Population-level summary
    pop_summary = {
        'Metric': [
            'Mean F_ROH',
            'SD F_ROH',
            'Min F_ROH',
            'Max F_ROH',
            'Mean Total ROH Length (Mb)',
            'Mean Number of ROH',
            'Mean Short ROH Count',
            'Mean Long ROH Count',
            'Mean VeryLong ROH Count',
            'Total ROH Segments (all individuals)',
            'Mean % Genome in ROH'
        ],
        'Value': [
            stats_df['F_ROH'].mean(),
            stats_df['F_ROH'].std(),
            stats_df['F_ROH'].min(),
            stats_df['F_ROH'].max(),
            stats_df['Total_ROH_Length_MB'].mean(),
            stats_df['Num_ROH'].mean(),
            stats_df['Short_ROH_Count'].mean(),
            stats_df['Long_ROH_Count'].mean(),
            stats_df['VeryLong_ROH_Count'].mean(),
            len(roh_df),
            stats_df['Pct_Genome_in_ROH'].mean()
        ]
    }
    
    pop_df = pd.DataFrame(pop_summary)
    pop_df.to_csv(f"{output_dir}/roh_population_summary.csv", index=False)
    print(f"Population summary saved to: {output_dir}/roh_population_summary.csv")
    
    # ROH length distribution by category
    category_dist = roh_df['ROH_Category'].value_counts().sort_index()
    category_dist.to_csv(f"{output_dir}/roh_category_distribution.csv", header=['Count'])
    print(f"Category distribution saved to: {output_dir}/roh_category_distribution.csv")
    
    # Print summary statistics
    print("\n" + "="*70)
    print("ROH ANALYSIS SUMMARY")
    print("="*70)
    print(f"\nNumber of individuals: {len(stats_df)}")
    print(f"Total ROH segments: {len(roh_df)}")
    print(f"\nPopulation-level statistics:")
    print(f"  Mean F_ROH: {stats_df['F_ROH'].mean():.6f} ± {stats_df['F_ROH'].std():.6f}")
    print(f"  Range: {stats_df['F_ROH'].min():.6f} - {stats_df['F_ROH'].max():.6f}")
    print(f"  Mean total ROH length: {stats_df['Total_ROH_Length_MB'].mean():.2f} Mb")
    print(f"  Mean number of ROH: {stats_df['Num_ROH'].mean():.1f}")
    print(f"\nROH by category:")
    print(f"  Short (<1Mb): {category_dist.get('short', 0)} segments")
    print(f"  Long (1-5Mb): {category_dist.get('long', 0)} segments")
    print(f"  Very Long (>5Mb): {category_dist.get('verylong', 0)} segments")
    
    # Individuals with highest F_ROH
    print(f"\nTop 5 individuals with highest F_ROH:")
    for idx, row in stats_df.head(5).iterrows():
        print(f"  {row['IID']}: F_ROH = {row['F_ROH']:.6f} ({row['Num_ROH']} ROH)")
    
    print("\n" + "="*70)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python analyze_roh.py <output_dir>")
        sys.exit(1)
    main()
EOF

# Run Python analysis
python3 "${OUTPUT_DIR}/analyze_roh.py" "${OUTPUT_DIR}"

log_message ""

# ============================================================================
# Step 4: Generate ROH Visualizations
# ============================================================================

log_message "Step 4: Generating ROH visualizations..."

# Create R script for visualization
cat > "${OUTPUT_DIR}/plot_roh.R" << 'EOF'
#!/usr/bin/env Rscript
# ROH Visualization Script

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# Load utility functions
source(file.path(Sys.getenv("PLM_ANALYSIS_DIR"), "utils.R"))

# Output directory
output_dir <- file.path(Sys.getenv("PLM_BASE_DIR"), "output/phase3a_population_genomics/roh_analysis")

# Read data
stats_df <- read.csv(file.path(output_dir, "roh_summary_per_individual.csv"))
roh_df <- read.table(file.path(output_dir, "plink_roh.hom"), header = TRUE)

# Add ROH category
roh_df$ROH_Category <- cut(
  roh_df$KB * 1000,
  breaks = c(0, 1e6, 5e6, Inf),
  labels = c("short", "long", "verylong")
)

# 1. F_ROH distribution
p1 <- ggplot(stats_df, aes(x = F_ROH)) +
  geom_histogram(bins = 20, fill = "coral", color = "black", alpha = 0.7) +
  geom_vline(aes(xintercept = mean(F_ROH)), color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Distribution of Genomic Inbreeding Coefficient",
    x = expression(F[ROH]),
    y = "Count"
  ) +
  theme_publication()

ggsave(file.path(output_dir, "froh_distribution.png"), p1, width = 8, height = 6, dpi = 300)

# 2. ROH length distribution by category
p2 <- ggplot(roh_df, aes(x = KB, fill = ROH_Category)) +
  geom_histogram(bins = 50, alpha = 0.7, color = "black") +
  scale_x_log10(labels = scales::comma) +
  scale_fill_manual(
    values = c("short" = "lightblue", "long" = "orange", "verylong" = "darkred"),
    labels = c("Short (<1Mb)", "Long (1-5Mb)", "Very Long (>5Mb)")
  ) +
  labs(
    title = "Distribution of ROH Lengths",
    x = "ROH Length (kb, log scale)",
    y = "Count",
    fill = "ROH Category"
  ) +
  theme_publication()

ggsave(file.path(output_dir, "roh_length_distribution.png"), p2, width = 10, height = 6, dpi = 300)

# 3. Number of ROH per individual
p3 <- ggplot(stats_df, aes(x = reorder(IID, Num_ROH), y = Num_ROH)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  labs(
    title = "Number of ROH per Individual",
    x = "Individual",
    y = "Number of ROH"
  ) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

ggsave(file.path(output_dir, "roh_count_per_individual.png"), p3, width = 12, height = 6, dpi = 300)

# 4. ROH category composition per individual
stats_long <- stats_df %>%
  select(IID, Short_ROH_Count, Long_ROH_Count, VeryLong_ROH_Count) %>%
  pivot_longer(cols = -IID, names_to = "Category", values_to = "Count")

stats_long$Category <- factor(
  stats_long$Category,
  levels = c("Short_ROH_Count", "Long_ROH_Count", "VeryLong_ROH_Count"),
  labels = c("Short", "Long", "Very Long")
)

p4 <- ggplot(stats_long, aes(x = reorder(IID, Count), y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Short" = "lightblue", "Long" = "orange", "Very Long" = "darkred")) +
  labs(
    title = "ROH Composition by Individual",
    x = "Individual",
    y = "Number of ROH",
    fill = "ROH Category"
  ) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

ggsave(file.path(output_dir, "roh_composition_per_individual.png"), p4, width = 12, height = 6, dpi = 300)

# 5. F_ROH vs Number of ROH
p5 <- ggplot(stats_df, aes(x = Num_ROH, y = F_ROH)) +
  geom_point(size = 3, alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
  labs(
    title = expression(paste(F[ROH], " vs Number of ROH")),
    x = "Number of ROH",
    y = expression(F[ROH])
  ) +
  theme_publication()

ggsave(file.path(output_dir, "froh_vs_numroh.png"), p5, width = 8, height = 6, dpi = 300)

cat("ROH visualizations generated successfully\n")
EOF

# Run R script
Rscript "${OUTPUT_DIR}/plot_roh.R"

if [ $? -eq 0 ]; then
    log_message "ROH visualizations generated successfully"
else
    log_message "WARNING: Some visualizations may have failed"
fi

# ============================================================================
# COMPLETION
# ============================================================================

log_message ""
log_message "====================================================================="
log_message "ROH ANALYSIS COMPLETE"
log_message "====================================================================="
log_message "Output directory: ${OUTPUT_DIR}"
log_message ""
log_message "Key output files:"
log_message "  - plink_roh.hom                      : PLINK ROH segments"
log_message "  - roh_summary_per_individual.csv     : Per-individual statistics"
log_message "  - roh_population_summary.csv         : Population-level summary"
log_message "  - froh_distribution.png              : F_ROH distribution plot"
log_message "  - roh_length_distribution.png        : ROH length distribution"
log_message "====================================================================="

exit 0


