#!/bin/bash
#
# Phase 3: SnpEff Functional Annotation
#
# Annotates variants with:
# - Gene regions and coding effects
# - Loss-of-function (LOF) variants
# - Impact categories (HIGH, MODERATE, LOW, MODIFIER)
#
# Prerequisites:
# - Custom SnpEff database must be pre-downloaded/built and available
# - Ensure the database name matches SNPEFF_DB below
#
# Usage: bash phase2_step2_snpeff_annotation.sh
#

set -euo pipefail

# ============================================================================
# SETUP
# ============================================================================

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/load_base_dir.sh"
DATA_DIR="${BASE_DIR}/data"
OUTPUT_DIR="${BASE_DIR}/output/phase2_annotation/snpeff_annotation"
VCF_IN="${BASE_DIR}/output/phase2_annotation/ancestral_states/variants_with_ancestral.vcf.gz"

# Fallback to original VCF if ancestral not available
if [ ! -f "$VCF_IN" ]; then
    VCF_IN="${DATA_DIR}/monkey_snp_sex_qc.vcf"
fi

mkdir -p "${OUTPUT_DIR}"

LOGFILE="${OUTPUT_DIR}/snpeff_annotation.log"
exec > >(tee -a "$LOGFILE") 2>&1

# ============================================================================
# PARAMETERS
# ============================================================================

# SnpEff database name (custom built for ASM756505v1)
SNPEFF_DB="Rhinopithecus_roxellana_ASM756505v1"

THREADS=8
MAX_MEMORY_GB=100  # Reserve memory for SnpEff
JAVA_OPTS="-Xmx80g -Xms4g"  # Java memory settings for SnpEff

# Resolve SnpEff command: SNPEFF_BIN, or java -jar $SNPEFF_JAR, or snpEff on PATH.
if [ -n "${SNPEFF_BIN:-}" ]; then
  :
elif [ -n "${SNPEFF_JAR:-}" ] && [ -f "${SNPEFF_JAR}" ]; then
  SNPEFF_BIN="java -jar ${SNPEFF_JAR}"
elif command -v snpEff >/dev/null 2>&1; then
  SNPEFF_BIN="snpEff"
else
  SNPEFF_BIN="snpEff"
fi

# Helper to execute SnpEff even if SNPEFF_BIN contains spaces
snpeff_cmd() {
    # Usage: snpeff_cmd <args...>
    # Example: snpeff_cmd ann -v DB input.vcf
    local quoted_args
    # shellcheck disable=SC2145
    quoted_args="$@"
    eval "${SNPEFF_BIN} ${quoted_args}"
}

# Ensure required external tools are available
if ! command -v bgzip >/dev/null 2>&1; then
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: 'bgzip' not found in PATH. Install htslib/tabix." >&2
  exit 1
fi
if ! command -v tabix >/dev/null 2>&1; then
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: 'tabix' not found in PATH. Install htslib/tabix." >&2
  exit 1
fi

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
log_message "SNPEFF FUNCTIONAL ANNOTATION"
log_message "====================================================================="
log_message ""

# ============================================================================
# Step 1: Check Pre-downloaded SnpEff Database
# ============================================================================

log_message "Step 1: Checking pre-downloaded SnpEff database..."

# Ensure SnpEff command works (supports SNPEFF_BIN with spaces)
# Note: SnpEff -h returns non-zero exit code, so we check if output contains "SnpEff version"
SNPEFF_CHECK=$(snpeff_cmd -h 2>&1 || true)
if ! echo "$SNPEFF_CHECK" | grep -q "SnpEff version"; then
    # Extra diagnostics
    if [[ "$SNPEFF_BIN" == *"-jar"* ]] || [[ -n "${SNPEFF_JAR:-}" ]]; then
        JAR_PATH="${SNPEFF_JAR:-}"
        if [ -z "$JAR_PATH" ]; then
            JAR_PATH=$(echo "$SNPEFF_BIN" | awk -F ' -jar ' '{print $2}')
        fi
        if [ -n "$JAR_PATH" ] && [ ! -r "$JAR_PATH" ]; then
            log_message "ERROR: SnpEff JAR not readable: $JAR_PATH"
        fi
        if ! command -v java >/dev/null 2>&1; then
            log_message "ERROR: 'java' not found in PATH. Install Java (e.g., apt install default-jre)."
        fi
    fi
    log_message "ERROR: SnpEff command not found or not executable: ${SNPEFF_BIN}"
    log_message "Hint: put snpEff on PATH, or set SNPEFF_JAR=/path/to/snpEff.jar"
    log_message "      or SNPEFF_BIN=\"java -jar /path/to/snpEff.jar\""
    exit 1
fi
log_message "SnpEff command verified: ${SNPEFF_BIN}"

log_message "Using SnpEff database: ${SNPEFF_DB}"

# Verify database is accessible via data.dir (more reliable for custom DBs)
DATA_DIR=$(snpeff_cmd config 2>/dev/null | grep -E '^data.dir' | awk -F'=' '{print $2}' | tr -d ' ' || true)
if [ -n "$DATA_DIR" ] && [ -d "$DATA_DIR/${SNPEFF_DB}" ]; then
    if [ -f "$DATA_DIR/${SNPEFF_DB}/snpEffectPredictor.bin" ] || [ -f "$DATA_DIR/${SNPEFF_DB}/genes.gtf" ]; then
        log_message "Database directory detected: $DATA_DIR/${SNPEFF_DB}"
    else
        log_message "WARNING: Database directory found but predictor files missing: $DATA_DIR/${SNPEFF_DB}"
    fi
else
    log_message "WARNING: Could not confirm database in data.dir (${DATA_DIR:-unknown})"
    log_message "Info: Official database list check (non-fatal):"
    snpeff_cmd databases | grep -F "${SNPEFF_DB}" || log_message "Database ${SNPEFF_DB} not listed (custom DBs often do not appear)"
fi

log_message ""

# ============================================================================
# Step 2: Annotate Variants with SnpEff
# ============================================================================

log_message "Step 2: Preparing VCF with correct chromosome names..."

# Convert chromosome names from simple numbers (1, 2, 3...) to RefSeq format (NC_044549.1, NC_044550.1...)
# This is required because SnpEff database uses RefSeq chromosome names
VCF_CONVERTED="${OUTPUT_DIR}/variants_converted_chrom.vcf.gz"

log_message "Converting chromosome names to RefSeq format..."

python3 - "$VCF_IN" "$VCF_CONVERTED" << 'EOPYTHON'
import gzip
import sys

# Chromosome name mapping: simple numbers -> RefSeq accessions
chrom_map = {
    '1': 'NC_044549.1', '2': 'NC_044550.1', '3': 'NC_044551.1',
    '4': 'NC_044552.1', '5': 'NC_044553.1', '6': 'NC_044554.1',
    '7': 'NC_044555.1', '8': 'NC_044556.1', '9': 'NC_044557.1',
    '10': 'NC_044558.1', '11': 'NC_044559.1', '12': 'NC_044560.1',
    '13': 'NC_044561.1', '14': 'NC_044562.1', '15': 'NC_044563.1',
    '16': 'NC_044564.1', '17': 'NC_044565.1', '18': 'NC_044566.1',
    '19': 'NC_044567.1', '20': 'NC_044568.1', '21': 'NC_044569.1',
    '22': 'NC_044570.1'
}

vcf_in = sys.argv[1]
vcf_out = sys.argv[2]

print(f"Converting chromosomes in {vcf_in}...", file=sys.stderr)

opener_in = gzip.open if vcf_in.endswith('.gz') else open
with opener_in(vcf_in, 'rt') as f_in, gzip.open(vcf_out, 'wt') as f_out:
    converted_count = 0
    unmapped_count = 0
    for line in f_in:
        if line.startswith('##contig='):
            # Update contig lines
            for old_chr, new_chr in chrom_map.items():
                if f'<ID={old_chr},' in line or f'<ID={old_chr}>' in line:
                    line = line.replace(f'ID={old_chr}', f'ID={new_chr}')
                    break
        elif line.startswith('#'):
            # Keep other header lines as-is
            pass
        else:
            # Convert chromosome names in variant lines
            fields = line.split('\t', 1)
            chrom = fields[0]
            if chrom in chrom_map:
                line = chrom_map[chrom] + '\t' + fields[1]
                converted_count += 1
            else:
                unmapped_count += 1
        
        f_out.write(line)

print(f"Converted {converted_count} variant chromosome names", file=sys.stderr)
if unmapped_count:
    print(f"WARNING: {unmapped_count} variant lines had unmapped chromosomes (left unchanged)", file=sys.stderr)
print(f"Output written to {vcf_out}", file=sys.stderr)
EOPYTHON

log_message "Chromosome conversion complete"
log_message ""

log_message "Step 3: Annotating variants with SnpEff..."

# Use converted VCF for SnpEff
VCF_IN="$VCF_CONVERTED"

# Run SnpEff annotation using custom-built database
log_message "Running SnpEff annotation with database: ${SNPEFF_DB}"
log_message "Note: Using custom database built from ASM756505v1 assembly"

export JAVA_OPTS="${JAVA_OPTS}"
ulimit -v $((MAX_MEMORY_GB * 1024 * 1024))  # Set virtual memory limit

snpeff_cmd ann \
    -v \
    -stats "${OUTPUT_DIR}/snpeff_stats.html" \
    ${SNPEFF_DB} \
    "$VCF_IN" | \
    bgzip > "${OUTPUT_DIR}/annotated_variants.vcf.gz"

if [ $? -eq 0 ]; then
    log_message "SnpEff annotation completed successfully"
    
    # Index output VCF
    tabix -p vcf "${OUTPUT_DIR}/annotated_variants.vcf.gz"
    
else
    log_message "ERROR: SnpEff annotation failed"
    exit 1
fi

log_message ""

# ============================================================================
# Step 3: Extract and Summarize Annotations
# ============================================================================

log_message "Step 3: Extracting and summarizing functional annotations..."

python3 << 'EOF'
import os
"""
Parse SnpEff annotations and generate summary
"""

import gzip
from collections import Counter, defaultdict

def parse_snpeff_vcf(vcf_file, output_dir):
    """Parse SnpEff annotated VCF and summarize"""
    
    print("\n" + "="*70)
    print("SNPEFF ANNOTATION SUMMARY")
    print("="*70)
    
    impact_counts = Counter()
    effect_counts = Counter()
    lof_counts = Counter()
    gene_impacts = defaultdict(list)
    
    lof_variants = []
    high_impact_variants = []
    
    with gzip.open(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            chrom, pos, vid, ref, alt, qual, filt, info = fields[:8]
            
            # Parse INFO field
            info_dict = {}
            for item in info.split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    info_dict[key] = value
            
            # Parse ANN field (SnpEff annotation)
            if 'ANN' in info_dict:
                annotations = info_dict['ANN'].split(',')
                
                for ann in annotations:
                    fields_ann = ann.split('|')
                    if len(fields_ann) >= 4:
                        effect = fields_ann[1]
                        impact = fields_ann[2]
                        gene = fields_ann[3]
                        
                        impact_counts[impact] += 1
                        effect_counts[effect] += 1
                        
                        if impact == 'HIGH':
                            high_impact_variants.append({
                                'chrom': chrom,
                                'pos': pos,
                                'ref': ref,
                                'alt': alt,
                                'gene': gene,
                                'effect': effect
                            })
                            if gene:
                                gene_impacts[gene].append(impact)
            
            # Check for LOF
            if 'LOF' in info_dict:
                lof_counts['total'] += 1
                lof_variants.append({
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'lof_info': info_dict['LOF']
                })
    
    # Print summary
    print(f"\nVariant Impact Distribution:")
    for impact in ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']:
        count = impact_counts[impact]
        pct = 100 * count / sum(impact_counts.values()) if impact_counts else 0
        print(f"  {impact}: {count} ({pct:.1f}%)")
    
    print(f"\nTop 10 Effect Types:")
    for effect, count in effect_counts.most_common(10):
        print(f"  {effect}: {count}")
    
    print(f"\nLoss-of-Function (LOF) Variants: {lof_counts['total']}")
    
    print(f"\nGenes with HIGH impact variants: {len(gene_impacts)}")
    
    # Save detailed results
    import csv
    
    # HIGH impact variants
    if high_impact_variants:
        with open(f"{output_dir}/high_impact_variants.csv", 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['chrom', 'pos', 'ref', 'alt', 'gene', 'effect'])
            writer.writeheader()
            writer.writerows(high_impact_variants)
        print(f"\nHIGH impact variants saved: {output_dir}/high_impact_variants.csv")
    
    # LOF variants
    if lof_variants:
        with open(f"{output_dir}/lof_variants.csv", 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['chrom', 'pos', 'ref', 'alt', 'lof_info'])
            writer.writeheader()
            writer.writerows(lof_variants)
        print(f"LOF variants saved: {output_dir}/lof_variants.csv")
    
    # Summary stats
    with open(f"{output_dir}/annotation_summary.txt", 'w') as f:
        f.write("SnpEff Annotation Summary\n")
        f.write("="*70 + "\n\n")
        f.write("Impact Distribution:\n")
        for impact, count in impact_counts.most_common():
            f.write(f"  {impact}: {count}\n")
        f.write("\nTop Effect Types:\n")
        for effect, count in effect_counts.most_common(20):
            f.write(f"  {effect}: {count}\n")
        f.write(f"\nLOF variants: {lof_counts['total']}\n")
        f.write(f"Genes with HIGH impact: {len(gene_impacts)}\n")
    
    print("\n" + "="*70)

    # Visualization
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from collections import OrderedDict

        # Impact distribution plot
        ordered_impacts = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']
        counts = [impact_counts[i] for i in ordered_impacts]
        plt.figure(figsize=(6,4))
        plt.bar(ordered_impacts, counts, color=["#d62728", "#ff7f0e", "#1f77b4", "#2ca02c"])
        plt.title('Variant Impact Distribution')
        plt.ylabel('Count')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/impact_distribution.png", dpi=150)
        plt.close()

        # Top 10 effects plot
        top_effects = effect_counts.most_common(10)
        if top_effects:
            labels = [e for e,_ in top_effects]
            values = [c for _,c in top_effects]
            plt.figure(figsize=(8,5))
            plt.barh(labels[::-1], values[::-1], color="#4e79a7")
            plt.xlabel('Count')
            plt.title('Top 10 Effect Types')
            plt.tight_layout()
            plt.savefig(f"{output_dir}/top_effect_types.png", dpi=150)
            plt.close()

        # HIGH impact genes count (bar)
        if gene_impacts:
            plt.figure(figsize=(5,3))
            plt.bar(["HIGH impact genes"], [len(gene_impacts)], color="#d62728")
            plt.ylabel('Number of genes')
            plt.tight_layout()
            plt.savefig(f"{output_dir}/high_impact_genes_count.png", dpi=150)
            plt.close()

        print(f"Visualization saved: {output_dir}/impact_distribution.png")
        if top_effects:
            print(f"Visualization saved: {output_dir}/top_effect_types.png")
        if gene_impacts:
            print(f"Visualization saved: {output_dir}/high_impact_genes_count.png")
    except Exception as e:
        print(f"WARNING: Plotting failed ({e}). Install matplotlib to enable plots.", file=sys.stderr)

# Run analysis
vcf_file = os.environ["PLM_BASE_DIR"] + "/output/phase2_annotation/snpeff_annotation/annotated_variants.vcf.gz"
output_dir = os.environ["PLM_BASE_DIR"] + "/output/phase2_annotation/snpeff_annotation"

parse_snpeff_vcf(vcf_file, output_dir)
EOF

# ============================================================================
# COMPLETION
# ============================================================================

log_message ""
log_message "====================================================================="
log_message "SNPEFF ANNOTATION COMPLETE"
log_message "====================================================================="
log_message "Output directory: ${OUTPUT_DIR}"
log_message ""
log_message "Key output files:"
log_message "  - annotated_variants.vcf.gz      : Functionally annotated VCF"
log_message "  - snpeff_stats.html              : SnpEff statistics report"
log_message "  - annotation_summary.txt         : Summary of annotations"
log_message "  - high_impact_variants.csv       : HIGH impact variants"
log_message "  - lof_variants.csv               : Loss-of-function variants"
log_message "====================================================================="

exit 0


