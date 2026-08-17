#!/bin/bash
#
# Phase 3: Ancestral State Inference
#
# Infers ancestral alleles using phylogenetic outgroups
# Uses parsimony approach with multiple outgroup species
# Uses minimap2 for fast genome-to-genome alignment
#
# Usage: bash phase2_step1_ancestral_states.sh [--clear-checkpoints]
#
# Options:
#   --clear-checkpoints  Clear all checkpoints and restart from beginning
#

set -euo pipefail

# ============================================================================
# COMMAND LINE ARGUMENTS
# ============================================================================

# Parse command line arguments
CLEAR_CHECKPOINTS=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --clear-checkpoints)
            CLEAR_CHECKPOINTS=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [--clear-checkpoints]"
            echo ""
            echo "Options:"
            echo "  --clear-checkpoints  Clear all checkpoints and restart from beginning"
            echo "  -h, --help          Show this help message"
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
# SETUP
# ============================================================================

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/load_base_dir.sh"
DATA_DIR="${BASE_DIR}/data"
OUTPUT_DIR="${BASE_DIR}/output/phase2_annotation/ancestral_states"
VCF="${DATA_DIR}/monkey_snp_sex_qc.vcf"

# Reference genomes (UPDATE THESE PATHS)
REF_GENOME="${BASE_DIR}/data/reference/Rhinopithecus_roxellana.fna"
OUTGROUP1="${BASE_DIR}/data/reference/Pygathrix_nemaeus.fna"      # Primary: Douc langur
OUTGROUP2="${BASE_DIR}/data/reference/Rhinopithecus_bieti.fna"    # Secondary: Black snub-nosed monkey
OUTGROUP3="${BASE_DIR}/data/reference/Macaca_mulatta.fna"         # Tertiary: Rhesus macaque

# Memory optimization parameters
MAX_MEMORY_GB=100  # Reserve 28GB for system
CHUNK_SIZE=1000000  # Process variants in chunks

mkdir -p "${OUTPUT_DIR}"

LOGFILE="${OUTPUT_DIR}/ancestral_inference.log"
exec > >(tee -a "$LOGFILE") 2>&1

# ============================================================================
# PARAMETERS
# ============================================================================

MIN_OUTGROUPS_AGREE=2  # Minimum outgroups with same allele for confidence
MIN_ALIGNMENT_QUALITY=20
THREADS=32
MEMORY_LIMIT="110G"  # Limit memory usage for alignment tools

# ============================================================================
# FUNCTIONS
# ============================================================================

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Checkpoint functions for resuming interrupted runs
checkpoint_file="${OUTPUT_DIR}/checkpoints.txt"

mark_checkpoint() {
    local step_name="$1"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "${step_name}:${timestamp}" >> "$checkpoint_file"
    log_message "Checkpoint marked: $step_name"
}

is_checkpoint_completed() {
    local step_name="$1"
    if [ -f "$checkpoint_file" ]; then
        grep -q "^${step_name}:" "$checkpoint_file"
        return $?
    fi
    return 1
}

check_output_files() {
    local step_name="$1"
    shift
    local files=("$@")
    
    for file in "${files[@]}"; do
        if [ ! -f "$file" ] || [ ! -s "$file" ]; then
            log_message "Output file missing or empty: $file"
            return 1
        fi
    done
    return 0
}

should_skip_step() {
    local step_name="$1"
    shift
    local output_files=("$@")
    
    if is_checkpoint_completed "$step_name" && check_output_files "$step_name" "${output_files[@]}"; then
        log_message "Step '$step_name' already completed - skipping"
        return 0
    else
        log_message "Step '$step_name' needs to run"
        return 1
    fi
}

clear_checkpoints() {
    if [ -f "$checkpoint_file" ]; then
        rm "$checkpoint_file"
        log_message "Checkpoints cleared"
    fi
}

# ============================================================================
# ANALYSIS
# ============================================================================

log_message "====================================================================="
log_message "ANCESTRAL STATE INFERENCE"
log_message "====================================================================="
log_message ""

# Clear checkpoints if requested
if [ "$CLEAR_CHECKPOINTS" = true ]; then
    clear_checkpoints
    log_message "All checkpoints cleared - starting fresh run"
    log_message ""
fi

# Check if reference genomes exist
if [ ! -f "$REF_GENOME" ]; then
    log_message "ERROR: Reference genome not found: $REF_GENOME"
    log_message "Please update the REF_GENOME path before running in production."
    exit 1
fi

for OG in "$OUTGROUP1" "$OUTGROUP2" "$OUTGROUP3"; do
    if [ ! -f "$OG" ]; then
        log_message "ERROR: Outgroup genome not found: $OG"
        log_message "Please update outgroup paths before running in production."
        exit 1
    fi
done

# ============================================================================
# Step 1: Genome Alignment to Outgroups
# ============================================================================

STEP1_NAME="genome_alignment"
STEP1_OUTPUTS=(
    "${OUTPUT_DIR}/alignment_outgroup1.paf"
    "${OUTPUT_DIR}/alignment_outgroup2.paf"
    "${OUTPUT_DIR}/alignment_outgroup3.paf"
)

if should_skip_step "$STEP1_NAME" "${STEP1_OUTPUTS[@]}"; then
    log_message "Step 1 already completed - skipping genome alignment"
else
    log_message "Step 1: Aligning reference genome to outgroup genomes using minimap2..."
    log_message "Note: minimap2 is much faster than LAST for whole-genome alignment"
    log_message ""

    # Using minimap2 for pairwise alignment (much faster than LAST)
    alignment_success=true

    for i in 1 2 3; do
        OUTGROUP_VAR="OUTGROUP${i}"
        OUTGROUP="${!OUTGROUP_VAR}"
        
        if [ ! -f "$OUTGROUP" ]; then
            log_message "WARNING: Outgroup $i genome not found: $OUTGROUP"
            continue
        fi
        
        log_message "Aligning to outgroup $i: $(basename $OUTGROUP)"
        
        # Run minimap2 alignment (with memory limits)
        log_message "  Running minimap2 alignment..."
        ulimit -v $((MAX_MEMORY_GB * 1024 * 1024))  # Set virtual memory limit
        
        # minimap2 parameters:
        # -x asm20: preset for genome-to-genome alignment
        # -t: number of threads
        # -N 50: retain up to 50 secondary alignments
        # -f 0.02: min fraction of matches
        # -c: output CIGAR in PAF (cg:Z: tag) for coordinate translation
        if minimap2 -x asm20 -t ${THREADS} -N 50 -f 0.02 -c \
           "$OUTGROUP" "$REF_GENOME" > "${OUTPUT_DIR}/alignment_outgroup${i}.paf"; then
            log_message "  Alignment complete for outgroup $i"
        else
            log_message "  ERROR: Alignment failed for outgroup $i"
            alignment_success=false
        fi
    done

    if [ "$alignment_success" = true ]; then
        mark_checkpoint "$STEP1_NAME"
        log_message "Step 1 completed successfully"
    else
        log_message "ERROR: Step 1 failed - some alignments did not complete"
        exit 1
    fi
fi

log_message ""

# ============================================================================
# Step 2: Extract Outgroup Alleles at Variant Positions
# ============================================================================

STEP2_NAME="ancestral_inference"
STEP2_OUTPUT="${OUTPUT_DIR}/variants_with_ancestral.vcf.gz"

if should_skip_step "$STEP2_NAME" "$STEP2_OUTPUT"; then
    log_message "Step 2 already completed - skipping ancestral inference"
else
    log_message "Step 2: Extracting outgroup alleles at variant positions..."

python3 << 'EOF'
"""
Extract outgroup alleles from PAF alignments and infer ancestral state
Production version: maps reference positions to outgroup coordinates using PAF CIGAR
and fetches bases from outgroup FASTAs.
"""

import sys
import os
import gzip
from collections import defaultdict, Counter

try:
    import pysam
    HAS_PYSAM = True
except Exception:
    HAS_PYSAM = False
import subprocess

def fetch_base(fasta_handle, chrom, pos1, strand='+'):
    """Fetch 1-based base from FASTA file handle.
    
    Args:
        fasta_handle: pysam.FastaFile object
        chrom: chromosome/scaffold name
        pos1: 1-based position
        strand: '+' or '-'
    
    Returns:
        Base at position (A, C, G, T, or N)
    """
    base = 'N'
    try:
        seq = fasta_handle.fetch(chrom, pos1-1, pos1)
        base = seq.upper() if seq else 'N'
    except Exception:
        base = 'N'
    
    if strand == '-':
        comp = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
        base = comp.get(base, 'N')
    return base

def parse_cigar(cigar):
    num = ''
    for ch in cigar:
        if ch.isdigit():
            num += ch
        else:
            yield int(num), ch
            num = ''

def build_alignment_index(paf_file):
    """Build index of alignments from PAF file.
    
    Index by query (reference) chromosome since we want to look up
    reference positions to find outgroup bases.
    """
    alignments_by_query = defaultdict(list)
    skipped_no_cigar = 0
    
    if not os.path.exists(paf_file):
        print(f"ERROR: PAF file not found: {paf_file}", file=sys.stderr)
        return alignments_by_query
    
    try:
        with open(paf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 12:
                    continue
                q_name = fields[0]  # Reference genome chromosome
                q_start = int(fields[2])
                q_end = int(fields[3])
                strand = fields[4]
                t_name = fields[5]  # Outgroup genome scaffold
                t_start = int(fields[7])
                t_end = int(fields[8])
                tags = fields[12:]
                cigar = None
                for tag in tags:
                    if tag.startswith('cg:Z:'):
                        cigar = tag.split(':',2)[2]
                        break
                if cigar is None:
                    skipped_no_cigar += 1
                    continue
                # Index by query (reference) chromosome for fast lookup
                alignments_by_query[q_name].append({
                    'q_start': q_start,
                    'q_end': q_end,
                    't_name': t_name,
                    't_start': t_start,
                    't_end': t_end,
                    'strand': strand,
                    'cigar': cigar
                })
        
        for chrom in alignments_by_query:
            alignments_by_query[chrom].sort(key=lambda a: a['q_start'])
        
        total_alns = sum(len(v) for v in alignments_by_query.values())
        print(f"  Loaded {total_alns} alignments from {paf_file}", file=sys.stderr)
        print(f"  Indexed by {len(alignments_by_query)} reference chromosomes", file=sys.stderr)
        if skipped_no_cigar > 0:
            print(f"  Skipped {skipped_no_cigar} alignments without CIGAR", file=sys.stderr)
            
    except Exception as e:
        print(f"ERROR parsing PAF file {paf_file}: {e}", file=sys.stderr)
    
    return alignments_by_query

def map_query_pos_to_target(aln, q_pos1):
    """Map 1-based query (reference) position to target (outgroup) position using PAF CIGAR.
    
    Args:
        aln: Alignment dict with q_start, q_end, t_name, t_start, t_end, strand, cigar
        q_pos1: 1-based position in query (reference) sequence
    
    Returns:
        Tuple of (t_name, t_pos1, strand) or None if position is in a gap
    """
    q = aln['q_start'] + 1  # convert to 1-based
    t = aln['t_start'] + 1
    
    for length, op in parse_cigar(aln['cigar']):
        if op in ('M', '=', 'X'):
            # Match/mismatch: both sequences advance
            if q_pos1 < q:
                return None
            if q_pos1 <= q + length - 1:
                offset = q_pos1 - q
                t_pos1 = t + offset
                return (aln['t_name'], t_pos1, aln['strand'])
            q += length
            t += length
        elif op == 'I':
            # Insertion to target: advances target only
            t += length
        elif op == 'D':
            # Deletion from target (gap in target): advances query only
            if q_pos1 < q:
                return None
            if q_pos1 <= q + length - 1:
                return None  # This query position is in a gap in the target
            q += length
        elif op in ('S', 'H'):
            # Soft/hard clips - don't advance coordinates
            pass
        else:
            # Other operations (N, P, etc) - skip
            continue
    return None

def infer_ancestral_allele(ref_allele, alt_alleles, outgroup_alleles, min_agree=2):
    if len(outgroup_alleles) == 0:
        return ref_allele, 0.0
    allele_counts = Counter()
    valid_set = set([ref_allele] + alt_alleles)
    for _, allele in outgroup_alleles:
        if allele in valid_set:
            allele_counts[allele] += 1
    if not allele_counts:
        return ref_allele, 0.0
    ancestral, count = allele_counts.most_common(1)[0]
    confidence = count / len(outgroup_alleles)
    if count < min_agree:
        confidence = confidence * 0.5
    return ancestral, confidence

def process_vcf(vcf_file, paf_files, outgroup_fastas, output_file, min_agree=2):
    print("\n" + "="*70)
    print("ANCESTRAL ALLELES FROM OUTGROUP ALIGNMENTS (PRODUCTION)")
    print("="*70)

    # Check if pysam is available
    if not HAS_PYSAM:
        print("\nERROR: pysam module is required for efficient FASTA access", file=sys.stderr)
        print("Please install: conda install -c bioconda pysam", file=sys.stderr)
        sys.exit(1)

    # Chromosome name mapping: VCF names -> RefSeq accessions
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
    print(f"\nChromosome mapping loaded: {len(chrom_map)} chromosomes", file=sys.stderr)

    # Build indices
    print("\nLoading PAF alignments...", file=sys.stderr)
    paf_indexes = [build_alignment_index(p) for p in paf_files]
    
    # Open FASTA files for outgroups (keep handles open for efficiency)
    print("\nOpening outgroup FASTA files...", file=sys.stderr)
    fasta_handles = []
    for fasta_path in outgroup_fastas:
        try:
            fh = pysam.FastaFile(fasta_path)
            fasta_handles.append(fh)
            print(f"  Opened: {os.path.basename(fasta_path)}", file=sys.stderr)
        except Exception as e:
            print(f"  ERROR opening {fasta_path}: {e}", file=sys.stderr)
            sys.exit(1)

    # Statistics counters
    stats = {
        'total_variants': 0,
        'with_outgroup_data': 0,
        'no_outgroup_data': 0,
        'high_confidence': 0,  # >= 0.8
        'medium_confidence': 0,  # 0.5-0.8
        'low_confidence': 0,  # < 0.5
        'outgroup_hits': [0] * len(paf_files),
        'chromosomes_seen': set(),
        'unmapped_chroms': set()
    }

    print("\nProcessing variants...", file=sys.stderr)
    opener = gzip.open if vcf_file.endswith('.gz') else open
    with opener(vcf_file, 'rt') as f_in, gzip.open(output_file, 'wt') as f_out:
        for line in f_in:
            if line.startswith('##'):
                f_out.write(line)
                continue
            if line.startswith('#CHROM'):
                f_out.write('##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele (inferred from outgroups)">\n')
                f_out.write('##INFO=<ID=AA_CONF,Number=1,Type=Float,Description="Ancestral Allele Confidence (0-1)">\n')
                f_out.write('##INFO=<ID=AA_METHOD,Number=1,Type=String,Description="Ancestral inference method">\n')
                f_out.write(line)
                continue

            stats['total_variants'] += 1
            
            fields = line.strip().split('\t')
            chrom, pos, vid, ref, alt, qual, filt, info = fields[:8]
            rest = fields[8:] if len(fields) > 8 else []
            pos1 = int(pos)
            alt_alleles = alt.split(',')

            # Map VCF chromosome name to RefSeq accession
            stats['chromosomes_seen'].add(chrom)
            ref_chrom = chrom_map.get(chrom, chrom)
            if chrom not in chrom_map:
                stats['unmapped_chroms'].add(chrom)

            outgroup_alleles = []
            for og_idx, paf_idx in enumerate(paf_indexes):
                alignments = paf_idx.get(ref_chrom, [])
                allele = None
                for aln in alignments:
                    # Check if variant position falls within this alignment
                    if pos1 < aln['q_start'] + 1 or pos1 > aln['q_end']:
                        continue
                    # Map reference position to outgroup position
                    result = map_query_pos_to_target(aln, pos1)
                    if result is None:
                        continue
                    t_name, t_pos1, strand = result
                    # Fetch base from outgroup genome using pre-opened file handle
                    allele = fetch_base(fasta_handles[og_idx], t_name, t_pos1, strand=strand)
                    break
                if allele is not None and allele != 'N':
                    outgroup_alleles.append((og_idx+1, allele))
                    stats['outgroup_hits'][og_idx] += 1

            if len(outgroup_alleles) > 0:
                stats['with_outgroup_data'] += 1
            else:
                stats['no_outgroup_data'] += 1

            ancestral, confidence = infer_ancestral_allele(ref, alt_alleles, outgroup_alleles, min_agree)
            
            # Classify confidence
            if confidence >= 0.8:
                stats['high_confidence'] += 1
            elif confidence >= 0.5:
                stats['medium_confidence'] += 1
            else:
                stats['low_confidence'] += 1
            
            info += f";AA={ancestral};AA_CONF={confidence:.3f};AA_METHOD=parsimony_minimap2"
            out_line = '\t'.join([chrom, pos, vid, ref, alt, qual, filt, info] + rest) + '\n'
            f_out.write(out_line)
            
            # Progress reporting
            if stats['total_variants'] % 10000 == 0:
                print(f"  Processed {stats['total_variants']:,} variants...", file=sys.stderr)
    
    # Print summary
    print("\n" + "="*70, file=sys.stderr)
    print("ANCESTRAL INFERENCE SUMMARY", file=sys.stderr)
    print("="*70, file=sys.stderr)
    print(f"\nTotal variants: {stats['total_variants']:,}", file=sys.stderr)
    print(f"  With outgroup data: {stats['with_outgroup_data']:,} ({100*stats['with_outgroup_data']/stats['total_variants']:.1f}%)", file=sys.stderr)
    print(f"  No outgroup data: {stats['no_outgroup_data']:,} ({100*stats['no_outgroup_data']/stats['total_variants']:.1f}%)", file=sys.stderr)
    print(f"\nConfidence distribution:", file=sys.stderr)
    print(f"  High (≥0.8): {stats['high_confidence']:,} ({100*stats['high_confidence']/stats['total_variants']:.1f}%)", file=sys.stderr)
    print(f"  Medium (0.5-0.8): {stats['medium_confidence']:,} ({100*stats['medium_confidence']/stats['total_variants']:.1f}%)", file=sys.stderr)
    print(f"  Low (<0.5): {stats['low_confidence']:,} ({100*stats['low_confidence']/stats['total_variants']:.1f}%)", file=sys.stderr)
    print(f"\nOutgroup coverage:", file=sys.stderr)
    for i, hits in enumerate(stats['outgroup_hits'], 1):
        print(f"  Outgroup {i}: {hits:,} variants ({100*hits/stats['total_variants']:.1f}%)", file=sys.stderr)
    print(f"\nChromosome mapping:", file=sys.stderr)
    print(f"  Chromosomes in VCF: {sorted(stats['chromosomes_seen'])}", file=sys.stderr)
    if stats['unmapped_chroms']:
        print(f"  WARNING: Unmapped chromosomes: {sorted(stats['unmapped_chroms'])}", file=sys.stderr)
    print(f"\nOutput written to: {output_file}", file=sys.stderr)
    print("="*70, file=sys.stderr)
    
    # Close FASTA file handles
    for fh in fasta_handles:
        fh.close()

vcf_file = os.environ["PLM_BASE_DIR"] + "/data/monkey_snp_sex_qc.vcf"
paf_files = [
    os.environ["PLM_BASE_DIR"] + "/output/phase2_annotation/ancestral_states/alignment_outgroup1.paf",
    os.environ["PLM_BASE_DIR"] + "/output/phase2_annotation/ancestral_states/alignment_outgroup2.paf",
    os.environ["PLM_BASE_DIR"] + "/output/phase2_annotation/ancestral_states/alignment_outgroup3.paf",
]
outgroup_fastas = [
    os.environ["PLM_BASE_DIR"] + "/data/reference/Pygathrix_nemaeus.fna",
    os.environ["PLM_BASE_DIR"] + "/data/reference/Rhinopithecus_bieti.fna",
    os.environ["PLM_BASE_DIR"] + "/data/reference/Macaca_mulatta.fna",
]
output_file = os.environ["PLM_BASE_DIR"] + "/output/phase2_annotation/ancestral_states/variants_with_ancestral.vcf.gz"

process_vcf(vcf_file, paf_files, outgroup_fastas, output_file, min_agree=2)
EOF

    mark_checkpoint "$STEP2_NAME"
    log_message "Step 2 completed successfully"
fi

log_message ""

# ============================================================================
# Step 3: Identify Derived Alleles
# ============================================================================

STEP3_NAME="derived_classification"
STEP3_OUTPUT="${OUTPUT_DIR}/ancestral_summary.txt"

if should_skip_step "$STEP3_NAME" "$STEP3_OUTPUT"; then
    log_message "Step 3 already completed - skipping derived allele classification"
else
    log_message "Step 3: Identifying derived alleles..."

python3 << 'EOF'
import os
"""
Identify and classify derived alleles
"""

import gzip
from collections import Counter

def classify_variants_by_ancestry(vcf_file, output_summary):
    """
    Classify variants by ancestral/derived status
    """
    
    print("\n" + "="*70)
    print("DERIVED ALLELE CLASSIFICATION")
    print("="*70)
    
    stats = Counter()
    derived_variants = []
    
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
            
            if 'AA' not in info_dict:
                stats['no_ancestral'] += 1
                continue
            
            ancestral = info_dict['AA']
            confidence = float(info_dict.get('AA_CONF', 0))
            
            # Determine derived allele(s)
            if ancestral == ref:
                # ALT is derived
                derived = alt
                stats['ref_ancestral'] += 1
            elif ancestral in alt.split(','):
                # REF is derived
                derived = ref
                stats['alt_ancestral'] += 1
            else:
                # Ambiguous
                stats['ambiguous'] += 1
                continue
            
            derived_variants.append({
                'chrom': chrom,
                'pos': pos,
                'ancestral': ancestral,
                'derived': derived,
                'confidence': confidence
            })
    
    print(f"\nAncestral state statistics:")
    print(f"  REF is ancestral: {stats['ref_ancestral']}")
    print(f"  ALT is ancestral: {stats['alt_ancestral']}")
    print(f"  Ambiguous: {stats['ambiguous']}")
    print(f"  No ancestral info: {stats['no_ancestral']}")
    
    # Save summary
    with open(output_summary, 'w') as f:
        f.write("Metric\tCount\n")
        for key, value in stats.items():
            f.write(f"{key}\t{value}\n")
    
    print(f"\nSummary saved to: {output_summary}")
    print("="*70)

vcf_file = os.environ["PLM_BASE_DIR"] + "/output/phase2_annotation/ancestral_states/variants_with_ancestral.vcf.gz"
output_summary = os.environ["PLM_BASE_DIR"] + "/output/phase2_annotation/ancestral_states/ancestral_summary.txt"

classify_variants_by_ancestry(vcf_file, output_summary)
EOF

    mark_checkpoint "$STEP3_NAME"
    log_message "Step 3 completed successfully"
fi

# ============================================================================
# COMPLETION
# ============================================================================

log_message ""
log_message "====================================================================="
log_message "ANCESTRAL STATE INFERENCE COMPLETE"
log_message "====================================================================="
log_message "Output directory: ${OUTPUT_DIR}"
log_message ""
log_message "Key output files:"
log_message "  - variants_with_ancestral.vcf.gz  : VCF with ancestral allele annotations"
log_message "  - ancestral_summary.txt           : Summary statistics"
log_message ""
log_message "Note: For production analysis, ensure outgroup genomes are available"
log_message "      and paths are correctly configured"
log_message "====================================================================="

exit 0


