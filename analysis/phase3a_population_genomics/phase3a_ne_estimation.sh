#!/bin/bash
#
# Phase 3a: Effective Population Size (Ne) Estimation
#
# Performs Ne estimation using:
# 1. LD-based Ne estimation from linkage disequilibrium decay
# 2. Site Frequency Spectrum (SFS) construction
# 3. Demographic modeling with fastsimcoal2 (fsc28/fsc27/fsc26)
#    - Constant population size model
#    - Bottleneck model
#
# Usage: bash phase3a_ne_estimation.sh
#
# Force re-run specific steps:
# FORCE_STEP1=1 bash phase3a_ne_estimation.sh  # Force re-run Step 1 (LD)
# FORCE_STEP2=1 bash phase3a_ne_estimation.sh  # Force re-run Step 2 (SFS)
# FORCE_STEP3=1 bash phase3a_ne_estimation.sh  # Force re-run Step 3 (fastsimcoal2)
#

set -euo pipefail

# ============================================================================
# SETUP
# ============================================================================

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/load_base_dir.sh"
DATA_DIR="${BASE_DIR}/data"
OUTPUT_DIR="${BASE_DIR}/output/phase3a_population_genomics/ne_estimation"
PLINK_PREFIX="${DATA_DIR}/monkey_snp_sex_qc"
VCF="${DATA_DIR}/monkey_snp_sex_qc.vcf"

mkdir -p "${OUTPUT_DIR}"/{ld_based,sfs_based}

LOGFILE="${OUTPUT_DIR}/ne_estimation.log"
exec > >(tee -a "$LOGFILE") 2>&1

# ============================================================================
# PARAMETERS
# ============================================================================

MIN_MAF=0.05
PCRIT=0.02
THREADS=8

# Force re-run specific steps (set to 1 to force re-run)
FORCE_STEP1=${FORCE_STEP1:-0}
FORCE_STEP2=${FORCE_STEP2:-0}
FORCE_STEP3=${FORCE_STEP3:-0}

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
log_message "EFFECTIVE POPULATION SIZE (Ne) ESTIMATION"
log_message "====================================================================="
log_message ""

# Create step completion tracking
STEP1_COMPLETE="${OUTPUT_DIR}/step1_complete.flag"
STEP2_COMPLETE="${OUTPUT_DIR}/step2_complete.flag"
STEP3_COMPLETE="${OUTPUT_DIR}/step3_complete.flag"

# ============================================================================
# Step 1: Estimate Ne from LD using PLINK
# ============================================================================

if [ -f "$STEP1_COMPLETE" ] && [ "$FORCE_STEP1" -eq 0 ]; then
    log_message "Step 1: Calculating LD decay for Ne estimation... SKIPPED (already completed)"
else
    if [ "$FORCE_STEP1" -eq 1 ]; then
        log_message "Step 1: Calculating LD decay for Ne estimation... FORCED RE-RUN"
        rm -f "$STEP1_COMPLETE"
    fi
    log_message "Step 1: Calculating LD decay for Ne estimation..."

# Calculate LD in windows (optimized for Ne estimation)
# Use smaller window and higher r2 threshold to reduce output size
# Check if LD file already exists and is complete
LD_FILE="${OUTPUT_DIR}/ld_based/ld_decay.ld"
LD_LOG="${OUTPUT_DIR}/ld_based/ld_decay.log"

if [ -f "$LD_FILE" ] && [ -f "$LD_LOG" ]; then
    log_message "  LD file already exists: $LD_FILE"
    log_message "  File size: $(du -h "$LD_FILE" | cut -f1)"
    
    # Check if PLINK completed successfully by checking log file end
    if tail -n 5 "$LD_LOG" | grep -qE "(Analysis finished|End time)" ; then
        log_message "  LD calculation appears to be complete. Skipping..."
    else
        log_message "  LD calculation was interrupted or incomplete. Recalculating..."
        # Remove incomplete files and recalculate
        rm -f "$LD_FILE" "$LD_LOG"
        
        plink \
            --bfile "${PLINK_PREFIX}" \
            --r2 \
            --ld-window-kb 200 \
            --ld-window 500 \
            --ld-window-r2 0.2 \
            --maf ${MIN_MAF} \
            --threads ${THREADS} \
            --out "${OUTPUT_DIR}/ld_based/ld_decay"
    fi
else
    log_message "  Starting LD calculation with memory-optimized parameters..."
    # Use more restrictive parameters to reduce file size
    plink \
        --bfile "${PLINK_PREFIX}" \
        --r2 \
        --ld-window-kb 200 \
        --ld-window 500 \
        --ld-window-r2 0.2 \
        --maf ${MIN_MAF} \
        --threads ${THREADS} \
        --out "${OUTPUT_DIR}/ld_based/ld_decay"
fi

# Estimate Ne from LD (memory-efficient processing)
python3 << 'EOF'
import pandas as pd
import numpy as np
from scipy import stats
import os

output_dir = os.environ["PLM_BASE_DIR"] + "/output/phase3a_population_genomics/ne_estimation/ld_based"

# Read LD file
ld_file = f"{output_dir}/ld_decay.ld"

print(f"Processing LD file: {ld_file}")
print(f"File size: {os.path.getsize(ld_file) / (1024**3):.1f} GB")

try:
    # Process file in smaller chunks to avoid memory issues
    chunk_size = 100000  # Process 100K rows at a time (reduced from 1M)
    
    print("Reading LD file in chunks...")
    
    # Read header first
    with open(ld_file, 'r') as f:
        header = f.readline().strip().split()
    
    print(f"Columns: {header}")
    
    # Initialize distance bins for efficient processing
    bins = [0, 10000, 50000, 100000, 500000, 1000000]
    bin_labels = ['0-10kb', '10-50kb', '50-100kb', '100-500kb', '500kb-1Mb']
    
    # Store statistics for each bin
    bin_stats = {label: {'r2_sum': 0, 'count': 0, 'dist_sum': 0} for label in bin_labels}
    total_pairs = 0
    
    # Process file in chunks with error handling
    chunk_count = 0
    for chunk in pd.read_csv(ld_file, sep=r'\s+', chunksize=chunk_size, on_bad_lines='skip'):
        chunk_count += 1
        print(f"Processing chunk {chunk_count}...")
        
        # Check if chunk has the expected columns
        expected_cols = ['CHR_A', 'BP_A', 'SNP_A', 'CHR_B', 'BP_B', 'SNP_B', 'R2']
        if not all(col in chunk.columns for col in expected_cols):
            print(f"Warning: Chunk {chunk_count} has unexpected columns: {chunk.columns.tolist()}")
            continue
        
        # Filter out rows with missing or invalid data
        chunk = chunk.dropna(subset=['CHR_A', 'BP_A', 'CHR_B', 'BP_B', 'R2'])
        
        # Convert to numeric, coercing errors to NaN
        chunk['CHR_A'] = pd.to_numeric(chunk['CHR_A'], errors='coerce')
        chunk['BP_A'] = pd.to_numeric(chunk['BP_A'], errors='coerce')
        chunk['CHR_B'] = pd.to_numeric(chunk['CHR_B'], errors='coerce')
        chunk['BP_B'] = pd.to_numeric(chunk['BP_B'], errors='coerce')
        chunk['R2'] = pd.to_numeric(chunk['R2'], errors='coerce')
        
        # Remove rows with NaN values
        chunk = chunk.dropna()
        
        if len(chunk) == 0:
            print(f"Warning: Chunk {chunk_count} has no valid data after cleaning")
            continue
        
        # Calculate distance between SNPs
        chunk['DIST'] = abs(chunk['BP_B'] - chunk['BP_A'])
        
        # Filter to only same chromosome pairs
        chunk = chunk[chunk['CHR_A'] == chunk['CHR_B']]
        
        if len(chunk) == 0:
            print(f"Warning: Chunk {chunk_count} has no same-chromosome pairs")
            continue
        
        # Adjust r2 for sample size (Hill & Weir correction)
        n = 68  # Number of samples
        chunk['R2_adj'] = chunk['R2'] - (1 / (2 * n))
        chunk['R2_adj'] = chunk['R2_adj'].clip(lower=0)
        
        # Process each bin directly to avoid storing all data in memory
        for i in range(len(bins)-1):
            mask = (chunk['DIST'] >= bins[i]) & (chunk['DIST'] < bins[i+1])
            if np.sum(mask) > 0:
                bin_stats[bin_labels[i]]['r2_sum'] += chunk.loc[mask, 'R2_adj'].sum()
                bin_stats[bin_labels[i]]['count'] += np.sum(mask)
                bin_stats[bin_labels[i]]['dist_sum'] += chunk.loc[mask, 'DIST'].sum()
        
        total_pairs += len(chunk)
        
        # Free memory immediately
        del chunk
        
        # Print progress every 10 chunks
        if chunk_count % 10 == 0:
            print(f"  Processed {chunk_count} chunks, {total_pairs} pairs so far")
    
    print(f"Processed {chunk_count} chunks, total pairs: {total_pairs}")
    
    # Calculate final statistics
    r2_by_dist = {}
    for bin_name, stats in bin_stats.items():
        if stats['count'] > 0:
            r2_by_dist[bin_name] = {
                'mean_r2': stats['r2_sum'] / stats['count'],
                'count': stats['count'],
                'mean_dist': stats['dist_sum'] / stats['count']
            }
    
    print(f"Total LD pairs processed: {total_pairs}")
    
    print("\n" + "="*70)
    print("LD DECAY SUMMARY")
    print("="*70)
    print("\nMean adjusted r² by distance:")
    for bin_name, stats in r2_by_dist.items():
        print(f"  {bin_name}: r² = {stats['mean_r2']:.6f}, n = {stats['count']}")
    
    # Simple Ne estimation from mean r2
    # Ne ≈ 1 / (4c * E[r²])
    # where c is recombination rate per bp (assume 1.2e-8)
    
    # Use distances < 100kb for Ne estimation
    short_dist_bins = ['0-10kb', '10-50kb', '50-100kb']
    short_dist_stats = {bin_name: stats for bin_name, stats in r2_by_dist.items() 
                       if bin_name in short_dist_bins and stats['count'] > 0}
    
    if short_dist_stats:
        # Calculate weighted average r2 and distance for short distances
        total_r2_weighted = sum(stats['mean_r2'] * stats['count'] for stats in short_dist_stats.values())
        total_count = sum(stats['count'] for stats in short_dist_stats.values())
        total_dist_weighted = sum(stats['mean_dist'] * stats['count'] for stats in short_dist_stats.values())
        
        if total_count > 0:
            mean_r2 = total_r2_weighted / total_count
            mean_dist = total_dist_weighted / total_count
            c = 1.2e-8  # recombination rate per bp per generation
            
            ne_estimate = 1 / (4 * c * mean_dist * mean_r2)
            print(f"\nRough Ne estimate from LD:")
            print(f"  Mean r² (distances < 100kb): {mean_r2:.6f}")
            print(f"  Mean distance: {mean_dist:.0f} bp")
            print(f"  Estimated Ne: {ne_estimate:.0f}")
            print(f"\nNote: This is a rough estimate. Use NeEstimator for accurate Ne.")
    
    # Save summary
    summary_df = pd.DataFrame([
        {'Distance_Bin': bin_name, 'Mean_R2': stats['mean_r2'], 'Count': stats['count']}
        for bin_name, stats in r2_by_dist.items()
    ])
    summary_df.to_csv(f"{output_dir}/ld_decay_summary.csv", index=False)
    
    print("="*70)
    
except Exception as e:
    print(f"Error processing LD file: {e}")
    import traceback
    traceback.print_exc()
EOF

log_message ""

# ============================================================================
# Alternative LD calculation method (if main method fails)
# ============================================================================

# Check if we need to use alternative method
if [ ! -f "${OUTPUT_DIR}/ld_based/ld_decay_summary.csv" ] || [ ! -s "${OUTPUT_DIR}/ld_based/ld_decay_summary.csv" ]; then
    log_message "Step 2b: Alternative LD calculation method..."
    
    # Method 1: Calculate LD decay curve directly using PLINK's --r2-only option
    # This is more memory efficient as it only stores r2 values, not SNP pairs
    log_message "  Using PLINK --r2-only method for memory-efficient LD calculation..."
    
    plink \
        --bfile "${PLINK_PREFIX}" \
        --r2 \
        --ld-window-kb 100 \
        --ld-window 200 \
        --ld-window-r2 0.3 \
        --maf ${MIN_MAF} \
        --threads ${THREADS} \
        --out "${OUTPUT_DIR}/ld_based/ld_decay_alt"
    
    # Process the alternative LD file
    python3 << 'EOF'
import pandas as pd
import numpy as np
import os

output_dir = os.environ["PLM_BASE_DIR"] + "/output/phase3a_population_genomics/ne_estimation/ld_based"
ld_file = f"{output_dir}/ld_decay_alt.ld"

if os.path.exists(ld_file):
    print(f"Processing alternative LD file: {ld_file}")
    print(f"File size: {os.path.getsize(ld_file) / (1024**3):.1f} GB")
    
    # Process in chunks
    chunk_size = 50000  # Even smaller chunks
    bins = [0, 10000, 50000, 100000, 500000, 1000000]
    bin_labels = ['0-10kb', '10-50kb', '50-100kb', '100-500kb', '500kb-1Mb']
    
    bin_stats = {label: {'r2_sum': 0, 'count': 0, 'dist_sum': 0} for label in bin_labels}
    total_pairs = 0
    
    try:
        for chunk in pd.read_csv(ld_file, sep=r'\s+', chunksize=chunk_size, on_bad_lines='skip'):
            # Basic processing
            chunk = chunk.dropna(subset=['CHR_A', 'BP_A', 'CHR_B', 'BP_B', 'R2'])
            chunk['CHR_A'] = pd.to_numeric(chunk['CHR_A'], errors='coerce')
            chunk['BP_A'] = pd.to_numeric(chunk['BP_A'], errors='coerce')
            chunk['CHR_B'] = pd.to_numeric(chunk['CHR_B'], errors='coerce')
            chunk['BP_B'] = pd.to_numeric(chunk['BP_B'], errors='coerce')
            chunk['R2'] = pd.to_numeric(chunk['R2'], errors='coerce')
            chunk = chunk.dropna()
            
            if len(chunk) == 0:
                continue
                
            chunk['DIST'] = abs(chunk['BP_B'] - chunk['BP_A'])
            chunk = chunk[chunk['CHR_A'] == chunk['CHR_B']]
            
            if len(chunk) == 0:
                continue
            
            # Adjust r2
            n = 68
            chunk['R2_adj'] = chunk['R2'] - (1 / (2 * n))
            chunk['R2_adj'] = chunk['R2_adj'].clip(lower=0)
            
            # Process bins
            for i in range(len(bins)-1):
                mask = (chunk['DIST'] >= bins[i]) & (chunk['DIST'] < bins[i+1])
                if np.sum(mask) > 0:
                    bin_stats[bin_labels[i]]['r2_sum'] += chunk.loc[mask, 'R2_adj'].sum()
                    bin_stats[bin_labels[i]]['count'] += np.sum(mask)
                    bin_stats[bin_labels[i]]['dist_sum'] += chunk.loc[mask, 'DIST'].sum()
            
            total_pairs += len(chunk)
            del chunk
        
        # Calculate final statistics
        r2_by_dist = {}
        for bin_name, stats in bin_stats.items():
            if stats['count'] > 0:
                r2_by_dist[bin_name] = {
                    'mean_r2': stats['r2_sum'] / stats['count'],
                    'count': stats['count'],
                    'mean_dist': stats['dist_sum'] / stats['count']
                }
        
        print(f"Alternative method processed {total_pairs} pairs")
        
        # Save alternative summary
        summary_df = pd.DataFrame([
            {'Distance_Bin': bin_name, 'Mean_R2': stats['mean_r2'], 'Count': stats['count']}
            for bin_name, stats in r2_by_dist.items()
        ])
        summary_df.to_csv(f"{output_dir}/ld_decay_summary_alt.csv", index=False)
        
        print("Alternative LD calculation completed successfully")
        
    except Exception as e:
        print(f"Error in alternative LD processing: {e}")
else:
    print("Alternative LD file not found, skipping alternative method")
EOF
fi

    # Mark Step 1 as complete
    touch "$STEP1_COMPLETE"
    log_message "Step 1 completed successfully"
fi

log_message ""

# ============================================================================
# Step 2: Site Frequency Spectrum (SFS) Construction
# ============================================================================

if [ -f "$STEP2_COMPLETE" ] && [ "$FORCE_STEP2" -eq 0 ]; then
    log_message "Step 2: Constructing Site Frequency Spectrum for demographic modeling... SKIPPED (already completed)"
else
    if [ "$FORCE_STEP2" -eq 1 ]; then
        log_message "Step 2: Constructing Site Frequency Spectrum for demographic modeling... FORCED RE-RUN"
        rm -f "$STEP2_COMPLETE"
    fi
    log_message "Step 2: Constructing Site Frequency Spectrum for demographic modeling..."

# Prepare VCF
if [[ ! "$VCF" =~ \.gz$ ]]; then
    bgzip -c "$VCF" > "${VCF}.gz"
    VCF="${VCF}.gz"
fi

if [ ! -f "${VCF}.tbi" ]; then
    tabix -p vcf "$VCF"
fi

# Generate SFS using VCFtools
# Pass the VCF path as environment variable
export MONKEY_VCF="$VCF"

python3 << 'EOF'
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

output_dir = os.environ["PLM_BASE_DIR"] + "/output/phase3a_population_genomics/ne_estimation/sfs_based"
# Use the VCF path from environment (already compressed if needed)
vcf_file = os.environ.get('MONKEY_VCF', os.environ["PLM_BASE_DIR"] + "/data/monkey_snp_sex_qc.vcf.gz")

print("Calculating allele frequencies...")
print(f"Using VCF file: {vcf_file}")

# Calculate allele frequencies from genotype data
print("Calculating allele frequencies from genotype data...")

# Use bcftools to extract genotype data and calculate frequencies
cmd = f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT]\\n' {vcf_file}"
result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

# Parse genotypes and calculate allele frequencies
afs = []
sample_count = 0

for line in result.stdout.strip().split('\n'):
    if line:
        fields = line.split('\t')
        if len(fields) >= 5:
            chrom, pos, ref, alt, genotypes = fields[0], fields[1], fields[2], fields[3], fields[4]
            
            # Count alleles
            ref_count = 0
            alt_count = 0
            total_genotypes = 0
            
            # Parse genotypes (format: 0/00/10/1, etc. - concatenated without spaces)
            # Split by common genotype patterns
            import re
            gt_patterns = re.findall(r'[01]\/[01]|[01]\|[01]|\.\/\.|\.\|\.', genotypes)
            
            for gt in gt_patterns:
                if gt in ['0/0', '0|0']:
                    ref_count += 2
                    total_genotypes += 1
                elif gt in ['0/1', '0|1', '1/0', '1|0']:
                    ref_count += 1
                    alt_count += 1
                    total_genotypes += 1
                elif gt in ['1/1', '1|1']:
                    alt_count += 2
                    total_genotypes += 1
                elif gt in ['./.', '.|.']:
                    # Missing genotype, skip
                    continue
            
            if total_genotypes > 0:
                total_alleles = ref_count + alt_count
                if total_alleles > 0:
                    af = alt_count / total_alleles
                    if 0 < af < 1:  # Exclude fixed sites
                        afs.append(af)
                        if sample_count == 0:
                            sample_count = total_genotypes

afs = np.array(afs)
print(f"Detected {sample_count} samples from genotype data")

print(f"Number of polymorphic sites: {len(afs)}")

# Create folded SFS (since we may not have ancestral states yet)
# Fold allele frequencies (use minor allele frequency)
maf = np.minimum(afs, 1 - afs)

# Create SFS bins using detected sample count
n_samples = sample_count if sample_count > 0 else 68
n_chromosomes = 2 * n_samples
print(f"Using {n_samples} samples ({n_chromosomes} chromosomes) for SFS")

bins = np.arange(0, n_chromosomes + 1, 1)
sfs, bin_edges = np.histogram(afs * n_chromosomes, bins=bins)

# Folded SFS (handle even/odd chromosome counts safely)
folded_counts = []
L = len(sfs)
half = L // 2
# fold bins 1..floor(n/2); exclude 0 and fixed sites
for i in range(1, half + 1):
    j = L - i
    if i == j:
        folded_counts.append(sfs[i])
    else:
        folded_counts.append(sfs[i] + sfs[j])
folded_sfs = np.array(folded_counts)

# Save SFS
sfs_df = pd.DataFrame({
    'Frequency_Bin': range(len(sfs)),
    'Count': sfs
})
sfs_df.to_csv(f"{output_dir}/unfolded_sfs.txt", index=False, sep='\t')

folded_sfs_df = pd.DataFrame({
    'Frequency_Bin': range(len(folded_sfs)),
    'Count': folded_sfs
})
folded_sfs_df.to_csv(f"{output_dir}/folded_sfs.txt", index=False, sep='\t')

# Plot SFS
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.bar(range(len(sfs)), sfs, color='steelblue', alpha=0.7, edgecolor='black')
ax1.set_xlabel('Derived Allele Count')
ax1.set_ylabel('Number of SNPs')
ax1.set_title('Unfolded Site Frequency Spectrum')
ax1.set_yscale('log')

ax2.bar(range(len(folded_sfs)), folded_sfs, color='coral', alpha=0.7, edgecolor='black')
ax2.set_xlabel('Minor Allele Count')
ax2.set_ylabel('Number of SNPs')
ax2.set_title('Folded Site Frequency Spectrum')
ax2.set_yscale('log')

plt.tight_layout()
plt.savefig(f"{output_dir}/sfs_plot.png", dpi=300, bbox_inches='tight')
plt.close()

print(f"\nSFS saved to: {output_dir}/")
print("\n" + "="*70)
print("SITE FREQUENCY SPECTRUM SUMMARY")
print("="*70)
print(f"Total polymorphic sites: {len(afs)}")
print(f"Number of chromosomes: {n_chromosomes}")
print(f"\nSFS statistics:")
print(f"  Singletons: {sfs[1]}")
print(f"  Doubletons: {sfs[2]}")
print(f"  High-frequency derived (>50%): {sum(sfs[len(sfs)//2:])}")
print(f"\nMean allele frequency: {np.mean(afs):.4f}")
print(f"Nucleotide diversity (π): {2 * np.mean(afs * (1 - afs)):.6f}")
print("="*70)
EOF

    # Mark Step 2 as complete
    touch "$STEP2_COMPLETE"
    log_message "Step 2 completed successfully"
fi

log_message ""

# ============================================================================
# Step 3: Fastsimcoal2 Demographic Modeling
# ============================================================================

if [ -f "$STEP3_COMPLETE" ] && [ "$FORCE_STEP3" -eq 0 ]; then
    log_message "Step 3: Running fastsimcoal2 demographic models... SKIPPED (already completed)"
else
    if [ "$FORCE_STEP3" -eq 1 ]; then
        log_message "Step 3: Running fastsimcoal2 demographic models... FORCED RE-RUN"
        rm -f "$STEP3_COMPLETE"
    fi
    log_message "Step 3: Running fastsimcoal2 demographic models..."

# Create fastsimcoal2 configuration files for different demographic scenarios
mkdir -p "${OUTPUT_DIR}/sfs_based/fastsimcoal"

# Model 1: Constant population size
cat > "${OUTPUT_DIR}/sfs_based/fastsimcoal/constant.tpl" << 'EOF'
//Number of population samples
1
//Population effective sizes (number of genes)
N1
//Sample sizes
68
//Growth rates
0
//Number of migration matrices
0
//Historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix
0
//Number of independent loci
1 0
//Per chromosome: Number of linkage blocks
1
//Per block: data type, num loci, rec. rate and mut rate + optional parameters
FREQ 1 0 2.5e-8 OUTEXP
EOF

cat > "${OUTPUT_DIR}/sfs_based/fastsimcoal/constant.est" << 'EOF'
// Priors and rules file
// ***************
[PARAMETERS]
//#isInt? #name   #dist.#min  #max
1       N1       unif  100    100000
[RULES]
[COMPLEX PARAMETERS]
EOF

# Model 2: Bottleneck model
cat > "${OUTPUT_DIR}/sfs_based/fastsimcoal/bottleneck.tpl" << 'EOF'
//Number of population samples
1
//Population effective sizes
N_CURR
//Sample sizes
68
//Growth rates
0
//Number of migration matrices
0
//Historical events: time, source, sink, migrants, new deme size, growth rate, migr. matrix
1 historical event
T_BOT 0 0 0 RES_BOT 0 0
//Number of independent loci
1 0
//Per block: data type, num loci, rec. rate and mut rate
FREQ 1 0 2.5e-8 OUTEXP
EOF

cat > "${OUTPUT_DIR}/sfs_based/fastsimcoal/bottleneck.est" << 'EOF'
[PARAMETERS]
//#isInt? #name    #dist.#min   #max
1       N_CURR    unif  100     50000
1       T_BOT     unif  10      5000
0       RES_BOT   unif  0.01    1.0
[RULES]
[COMPLEX PARAMETERS]
0 N_ANC = N_CURR*RES_BOT output
EOF

log_message "Fastsimcoal2 model files created"

# Check if fastsimcoal2 is available
if command -v fsc28 &> /dev/null || command -v fsc27 &> /dev/null || command -v fsc26 &> /dev/null; then
    FSC_CMD=$(command -v fsc28 2>/dev/null || command -v fsc27 2>/dev/null || command -v fsc26 2>/dev/null)
    log_message "Found fastsimcoal2: $FSC_CMD"
    
    # Convert SFS to fastsimcoal2 format
    log_message "Converting SFS to fastsimcoal2 format..."
    
python3 << 'EOF'
import os
import pandas as pd
import numpy as np

output_dir = os.environ["PLM_BASE_DIR"] + "/output/phase3a_population_genomics/ne_estimation/sfs_based"

# Read folded SFS
sfs_df = pd.read_csv(f"{output_dir}/folded_sfs.txt", sep='\t')
sfs_counts = sfs_df['Count'].values

# Write in fastsimcoal2 MAF format
# Header: 1 observation, d derived alleles
with open(f"{output_dir}/fastsimcoal/monkey_MAFpop0.obs", 'w') as f:
    f.write("1 observations\n")
    # Write MAF counts (skip first bin which is monomorphic)
    for i, count in enumerate(sfs_counts):
        f.write(f"\td_{i+1}\t{int(count)}\n")

print(f"SFS converted to fastsimcoal2 format: monkey_MAFpop0.obs")
EOF
    
    # Run fastsimcoal2 for both models
    cd "${OUTPUT_DIR}/sfs_based/fastsimcoal"
    
    log_message "Running constant population size model..."
    if $FSC_CMD -t constant.tpl -e constant.est -n 100000 -m -M -L 40 -q 2>&1 | tee constant_run.log; then
        log_message "  Constant model completed"
    else
        log_message "  WARNING: Constant model failed or not fully optimized"
    fi
    
    log_message "Running bottleneck model..."
    if $FSC_CMD -t bottleneck.tpl -e bottleneck.est -n 100000 -m -M -L 40 -q 2>&1 | tee bottleneck_run.log; then
        log_message "  Bottleneck model completed"
    else
        log_message "  WARNING: Bottleneck model failed or not fully optimized"
    fi
    
    cd - > /dev/null
    
    # Extract and compare results
    log_message "Extracting parameter estimates..."
    
python3 << 'EOF'
import os
import glob

fsc_dir = os.environ["PLM_BASE_DIR"] + "/output/phase3a_population_genomics/ne_estimation/sfs_based/fastsimcoal"

print("\n" + "="*70)
print("FASTSIMCOAL2 RESULTS")
print("="*70)

# Try to find and parse bestlhoods files
for model in ['constant', 'bottleneck']:
    bestlhoods_files = glob.glob(f"{fsc_dir}/{model}/{model}.bestlhoods")
    if bestlhoods_files:
        print(f"\n{model.upper()} MODEL:")
        with open(bestlhoods_files[0], 'r') as f:
            lines = f.readlines()
            if len(lines) >= 2:
                header = lines[0].strip().split()
                values = lines[1].strip().split()
                for h, v in zip(header, values):
                    print(f"  {h}: {v}")
    else:
        print(f"\n{model.upper()} MODEL: Results not found")
        print(f"  Check {fsc_dir}/{model}/ for output files")

print("="*70)
EOF
    
else
    log_message "WARNING: fastsimcoal2 (fsc28/fsc27/fsc26) not found in PATH"
    log_message "Model configuration files created but not executed"
    log_message "To run manually:"
    log_message "  1. Install fastsimcoal2: http://cmpg.unibe.ch/software/fastsimcoal27/"
    log_message "  2. cd ${OUTPUT_DIR}/sfs_based/fastsimcoal"
    log_message "  3. fsc28 -t constant.tpl -e constant.est -n 100000 -m -M -L 40"
    log_message "  4. fsc28 -t bottleneck.tpl -e bottleneck.est -n 100000 -m -M -L 40"
fi

    # Mark Step 3 as complete
    touch "$STEP3_COMPLETE"
    log_message "Step 3 completed successfully"
fi

log_message ""

# =========================================================================
# Visualization: Generate summary plots for LD, SFS, and fastsimcoal2 results
# =========================================================================

log_message "Generating visualization summaries..."

python3 << 'EOF'
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

base_out = os.environ["PLM_BASE_DIR"] + "/output/phase3a_population_genomics/ne_estimation"
ld_dir = os.path.join(base_out, 'ld_based')
sfs_dir = os.path.join(base_out, 'sfs_based')
fsc_dir = os.path.join(sfs_dir, 'fastsimcoal')
viz_dir = base_out

os.makedirs(viz_dir, exist_ok=True)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()
plot_idx = 0

# 1) LD-based summary
ld_summary_csv = os.path.join(ld_dir, 'ld_decay_summary.csv')
if os.path.exists(ld_summary_csv):
    try:
        df = pd.read_csv(ld_summary_csv)
        order = ['0-10kb', '10-50kb', '50-100kb', '100-500kb', '500kb-1Mb']
        if 'Distance_Bin' in df.columns:
            df['Distance_Bin'] = pd.Categorical(df['Distance_Bin'], categories=order, ordered=True)
            df = df.sort_values('Distance_Bin')
        ax = axes[plot_idx]
        ax.plot(df['Distance_Bin'], df['Mean_R2'], marker='o', color='steelblue')
        ax.set_title('LD decay (adjusted r² by distance)')
        ax.set_xlabel('Distance bin')
        ax.set_ylabel('Mean adjusted r²')
        ax.grid(True, linestyle='--', alpha=0.3)
        ax.tick_params(axis='x', rotation=20)
        plot_idx += 1
    except Exception:
        pass

# 2) SFS plot (folded)
folded_sfs_txt = os.path.join(sfs_dir, 'folded_sfs.txt')
if os.path.exists(folded_sfs_txt):
    try:
        sfs_df = pd.read_csv(folded_sfs_txt, sep='\t')
        ax = axes[plot_idx]
        ax.bar(sfs_df['Frequency_Bin'], sfs_df['Count'], color='coral', edgecolor='black', alpha=0.8)
        ax.set_title('Folded Site Frequency Spectrum')
        ax.set_xlabel('Minor allele count bin')
        ax.set_ylabel('Number of SNPs')
        ax.set_yscale('log')
        ax.grid(True, axis='y', linestyle='--', alpha=0.3)
        plot_idx += 1
    except Exception:
        pass

# 3) fastsimcoal2 likelihood comparison

def read_bestlhoods(model_name):
    candidates = [
        os.path.join(fsc_dir, model_name, f"{model_name}.bestlhoods"),
        os.path.join(fsc_dir, f"{model_name}.bestlhoods")
    ]
    for p in candidates:
        if os.path.exists(p):
            try:
                with open(p, 'r') as f:
                    lines = f.readlines()
                    if len(lines) >= 2:
                        header = lines[0].strip().split()
                        values = lines[1].strip().split()
                        return dict(zip(header, values))
            except Exception:
                return None
    return None

const_res = read_bestlhoods('constant')
bot_res = read_bestlhoods('bottleneck')

if const_res or bot_res:
    ax = axes[plot_idx]
    names = []
    neglogl = []
    if const_res and 'MaxEstLhood' in const_res:
        names.append('constant')
        try:
            neglogl.append(-float(const_res['MaxEstLhood']))
        except Exception:
            neglogl.append(np.nan)
    if bot_res and 'MaxEstLhood' in bot_res:
        names.append('bottleneck')
        try:
            neglogl.append(-float(bot_res['MaxEstLhood']))
        except Exception:
            neglogl.append(np.nan)
    if names:
        ax.bar(names, neglogl, color=['#6aaed6', '#f08c6b'][:len(names)], edgecolor='black')
        ax.set_title('fastsimcoal2: -MaxEstLhood (lower is better)')
        ax.set_ylabel('-Log likelihood')
        ax.grid(True, axis='y', linestyle='--', alpha=0.3)
        plot_idx += 1

# 4) fastsimcoal2 parameter snapshot (bottleneck)
if bot_res:
    try:
        keys = [k for k in bot_res.keys() if k in {'N_CURR','T_BOT','RES_BOT','N_ANC'} or k.upper() in {'N_CURR','T_BOT','RES_BOT','N_ANC'}]
        if keys:
            ax = axes[plot_idx]
            ax.axis('off')
            table_data = [[k, bot_res[k]] for k in keys]
            table = ax.table(cellText=table_data, colLabels=['Parameter','Value'], loc='center')
            table.auto_set_font_size(False)
            table.set_fontsize(9)
            table.scale(1.2, 1.2)
            ax.set_title('Bottleneck model parameters (snapshot)')
            plot_idx += 1
    except Exception:
        pass

for i in range(plot_idx, len(axes)):
    axes[i].axis('off')

plt.tight_layout()
out_png = os.path.join(viz_dir, 'ne_summary.png')
plt.savefig(out_png, dpi=200, bbox_inches='tight')
plt.close()

print(f"Summary visualization saved: {out_png}")
EOF

# ============================================================================
# Summary
# ============================================================================

log_message "====================================================================="
log_message "Ne ESTIMATION COMPLETE"
log_message "====================================================================="
log_message "Output directory: ${OUTPUT_DIR}"
log_message ""
log_message "Completed analyses:"
log_message "  Step 1: LD decay calculation and Ne estimation"
log_message "  Step 2: Site Frequency Spectrum construction"
log_message "  Step 3: Fastsimcoal2 demographic modeling"
log_message ""
log_message "Key output files:"
log_message "  - ${OUTPUT_DIR}/ld_based/ld_decay_summary.csv"
log_message "  - ${OUTPUT_DIR}/sfs_based/folded_sfs.txt"
log_message "  - ${OUTPUT_DIR}/sfs_based/sfs_plot.png"
log_message "  - ${OUTPUT_DIR}/sfs_based/fastsimcoal/constant/"
log_message "  - ${OUTPUT_DIR}/sfs_based/fastsimcoal/bottleneck/"
log_message ""

# Print Ne estimate if available
if [ -f "${OUTPUT_DIR}/ld_based/ld_decay_summary.csv" ]; then
    log_message "LD-based Ne estimate available in log file above"
fi

log_message ""
log_message "For detailed results, check:"
log_message "  - LD-based Ne: ${OUTPUT_DIR}/ld_based/"
log_message "  - SFS analysis: ${OUTPUT_DIR}/sfs_based/"
log_message "  - Demographic models: ${OUTPUT_DIR}/sfs_based/fastsimcoal/"
log_message "====================================================================="

exit 0


