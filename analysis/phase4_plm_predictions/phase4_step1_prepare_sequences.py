#!/usr/bin/env python3
"""
Phase 4: Prepare Protein Sequences and Missense Variants for PLM Predictions

Extracts protein sequences and generates wild-type/mutant pairs for all
missense variants from SnpEff annotated VCF.

Usage: python3 phase4_step1_prepare_sequences.py
"""

import sys
import os
import gzip
import csv
import logging
from collections import defaultdict
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils import load_config, setup_logger, extract_snpeff_annotation

# ============================================================================
# SETUP
# ============================================================================

BASE_DIR = get_base_dir()
OUTPUT_DIR = f"{BASE_DIR}/output/phase4_plm_predictions"
ANNOTATION_DIR = f"{BASE_DIR}/output/phase2_annotation/snpeff_annotation"

# Input files
ANNOTATED_VCF = f"{ANNOTATION_DIR}/annotated_variants.vcf.gz"
PROTEIN_FASTA = str(get_base_dir() / "data/reference/Rhinopithecus_roxellana_protein.faa")  # UPDATE THIS PATH
GFF_FILE = str(get_base_dir() / "data/reference/Rhinopithecus_roxellana.gff")

# Output files
OUTPUT_VARIANTS = f"{OUTPUT_DIR}/missense_variants_for_prediction.csv"
OUTPUT_VARIANTS_UNIQUE = f"{OUTPUT_DIR}/missense_variants_for_prediction_unique.csv"
OUTPUT_SEQUENCES = f"{OUTPUT_DIR}/protein_sequences.fasta"

# Setup logging
logger = setup_logger("phase4_prep", f"{OUTPUT_DIR}/sequence_preparation.log")

# ============================================================================
# FUNCTIONS
# ============================================================================

def load_protein_sequences(fasta_file):
    """
    Load protein sequences from FASTA file
    
    Returns:
        dict: {protein_id: sequence}
    """
    sequences = {}
    
    if not os.path.exists(fasta_file):
        logger.warning(f"Protein FASTA not found: {fasta_file}")
        return sequences
    
    current_id = None
    current_seq = []
    
    opener = gzip.open if fasta_file.endswith('.gz') else open
    
    with opener(fasta_file, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                # Take first token as primary ID (e.g., XP_..., NP_...)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    # Build alias map to improve matching (strip version suffix)
    augmented_sequences = {}
    for prot_id, seq in sequences.items():
        augmented_sequences[prot_id] = seq
        # Add versionless alias if applicable (e.g., XP_12345.1 -> XP_12345)
        if '.' in prot_id and prot_id.split('.')[-1].isdigit():
            versionless = prot_id.split('.')[0]
            augmented_sequences[versionless] = seq
    logger.info(f"Loaded {len(sequences)} protein sequences (augmented keys: {len(augmented_sequences)})")
    return augmented_sequences

def load_transcript_protein_map(gff_file):
    """
    Build a mapping from transcript RefSeq IDs (e.g., XM_..., NM_...) to protein RefSeq IDs (XP_..., NP_...)
    by scanning the GFF attributes where CDS features have 'Parent=rna-<XM_..>' and 'protein_id=XP_..'.
    
    Returns:
        dict: {transcript_id (with and without version): protein_id (with and without version)}
    """
    mapping = {}
    if not os.path.exists(gff_file):
        logger.warning(f"GFF file not found: {gff_file}")
        return mapping
    try:
        with open(gff_file, 'r') as gf:
            for line in gf:
                if not line or line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 9:
                    continue
                feature_type = parts[2]
                if feature_type != 'CDS':
                    continue
                attrs = parts[8]
                # Quick filters to avoid heavy parsing
                if 'Parent=rna-' not in attrs or 'protein_id=' not in attrs:
                    continue
                # Extract transcript and protein IDs from attributes
                # Example: Parent=rna-XM_030939677.1;...;protein_id=XP_030795537.1
                tid = None
                pid = None
                for field in attrs.split(';'):
                    if field.startswith('Parent=rna-'):
                        tid = field.split('Parent=rna-', 1)[1]
                    elif field.startswith('protein_id='):
                        pid = field.split('protein_id=', 1)[1]
                if not tid or not pid:
                    continue
                # Record mapping for versioned and versionless keys
                def add_pair(t, p):
                    mapping[t] = p
                    if '.' in t and t.split('.')[-1].isdigit():
                        mapping[t.split('.')[0]] = p
                    if '.' in p and p.split('.')[-1].isdigit():
                        mapping[t] = p
                        mapping[t.split('.')[0]] = p.split('.')[0]
                add_pair(tid, pid)
        logger.info(f"Loaded transcript→protein mappings from GFF: {len(mapping)} entries")
    except Exception as e:
        logger.warning(f"Failed to load transcript→protein map from GFF: {e}")
    return mapping

def parse_hgvs_protein(hgvs_p):
    """
    Parse HGVS protein notation to extract mutation info
    
    Args:
        hgvs_p: HGVS protein notation (e.g., "p.Lys41Glu")
    
    Returns:
        tuple: (position, wt_aa, mut_aa) or None if invalid
    """
    if not hgvs_p or hgvs_p == '' or hgvs_p == '.':
        return None
    
    # Remove 'p.' prefix
    if hgvs_p.startswith('p.'):
        hgvs_p = hgvs_p[2:]
    
    # Parse pattern like Lys41Glu or K41E
    import re
    
    # Three-letter code pattern
    pattern3 = r'([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})'
    match = re.match(pattern3, hgvs_p)
    
    if match:
        wt_aa3, pos, mut_aa3 = match.groups()
        
        # Convert three-letter to one-letter
        aa_map = {
            'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
            'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
            'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
            'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
        }
        
        wt_aa = aa_map.get(wt_aa3, wt_aa3)
        mut_aa = aa_map.get(mut_aa3, mut_aa3)
        
        return (int(pos), wt_aa, mut_aa)
    
    # One-letter code pattern
    pattern1 = r'([A-Z])(\d+)([A-Z])'
    match = re.match(pattern1, hgvs_p)
    
    if match:
        wt_aa, pos, mut_aa = match.groups()
        return (int(pos), wt_aa, mut_aa)
    
    return None

def extract_missense_variants(vcf_file):
    """
    Extract all missense variants from SnpEff annotated VCF
    
    Returns:
        list: List of variant dictionaries
    """
    logger.info("Extracting missense variants from VCF...")
    
    missense_variants = []
    
    with gzip.open(vcf_file, 'rt') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            chrom, pos, vid, ref, alt, qual, filt, info = fields[:8]
            
            # Parse INFO field for ANN
            info_dict = {}
            for item in info.split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    info_dict[key] = value
            
            if 'ANN' not in info_dict:
                continue
            
            # Parse all annotations
            annotations = info_dict['ANN'].split(',')
            
            for ann in annotations:
                ann_fields = ann.split('|')
                
                if len(ann_fields) < 11:
                    continue
                
                allele = ann_fields[0]
                effect = ann_fields[1]
                impact = ann_fields[2]
                gene_name = ann_fields[3]
                gene_id = ann_fields[4]
                feature_id = ann_fields[6]  # Transcript ID
                hgvs_p = ann_fields[10]     # HGVS protein
                
                # Filter for missense variants
                if 'missense' in effect.lower():
                    # Parse HGVS notation
                    mutation_info = parse_hgvs_protein(hgvs_p)
                    
                    if mutation_info:
                        pos_aa, wt_aa, mut_aa = mutation_info
                        
                        variant = {
                            'chrom': chrom,
                            'pos': pos,
                            'ref': ref,
                            'alt': alt,
                            'gene_name': gene_name,
                            'gene_id': gene_id,
                            'transcript_id': feature_id,
                            'protein_id': feature_id.replace('transcript:', 'protein:'),
                            'effect': effect,
                            'impact': impact,
                            'hgvs_p': hgvs_p,
                            'aa_pos': pos_aa,
                            'wt_aa': wt_aa,
                            'mut_aa': mut_aa,
                            'variant_id': f"{chrom}:{pos}:{ref}:{alt}"
                        }
                        
                        missense_variants.append(variant)
            
            if line_num % 10000 == 0:
                logger.info(f"  Processed {line_num} variants, found {len(missense_variants)} missense")
    
    logger.info(f"Total missense variants extracted: {len(missense_variants)}")
    return missense_variants

def generate_mutant_sequences(variants, protein_sequences):
    """
    Generate wild-type and mutant protein sequences for each variant
    
    Args:
        variants: List of variant dictionaries
        protein_sequences: Dictionary of protein sequences
    
    Returns:
        list: Variants with sequences added
    """
    logger.info("Generating wild-type and mutant sequences...")
    
    variants_with_seqs = []
    no_sequence_count = 0
    position_mismatch = 0
    
    # Load transcript→protein mapping from GFF
    transcript_to_protein = load_transcript_protein_map(GFF_FILE)
    logger.info(f"Loaded {len(transcript_to_protein)} transcript→protein ID mappings from GFF")

    def normalize_transcript_id(raw_tid):
        if not raw_tid:
            return raw_tid
        if raw_tid.startswith('transcript:'):
            return raw_tid.split(':', 1)[1]
        if raw_tid.startswith('rna:'):
            return raw_tid.split(':', 1)[1]
        if raw_tid.startswith('rna-'):
            return raw_tid.split('-', 1)[1]
        return raw_tid

    def candidate_ids_for_variant(v):
        ids = []
        tid = normalize_transcript_id(v['transcript_id'])
        pid = v['protein_id']
        gid = v['gene_id']
        gname = v['gene_name']
        
        def add(idv):
            if not idv or idv == '.':
                return
            ids.append(idv)
            # Also try versionless if it looks like RefSeq with version
            if '.' in idv and idv.split('.')[-1].isdigit():
                ids.append(idv.split('.')[0])
        
        # Provided IDs
        add(pid)
        add(tid)
        add(gid)
        add(gname)
        
        # Common RefSeq transcript->protein mappings
        if tid.startswith('XM_'):
            add('XP_' + tid[3:])
        if tid.startswith('NM_'):
            add('NP_' + tid[3:])
        # Sometimes feature_id may be prefixed (already normalized above)
        if pid.startswith('protein:'):
            add(pid.split(':', 1)[1])

        # Use GFF-derived mapping if available
        if tid in transcript_to_protein:
            add(transcript_to_protein[tid])
        # Also try versionless
        if '.' in tid and tid.split('.')[-1].isdigit():
            tv = tid.split('.')[0]
            if tv in transcript_to_protein:
                add(transcript_to_protein[tv])
        
        # De-duplicate while preserving order
        seen = set()
        ordered = []
        for x in ids:
            if x not in seen:
                seen.add(x)
                ordered.append(x)
        return ordered
    
    for var in variants:
        # Try to find sequence by various candidate IDs
        sequence = None
        matched_protein_id = None
        tid = normalize_transcript_id(var['transcript_id'])
        candidates = candidate_ids_for_variant(var)
        
        for seq_id in candidates:
            if seq_id in protein_sequences:
                sequence = protein_sequences[seq_id]
                matched_protein_id = seq_id
                break
        
        if not sequence:
            no_sequence_count += 1
            continue
        
        aa_pos = var['aa_pos']
        wt_aa = var['wt_aa']
        mut_aa = var['mut_aa']
        
        # Validate position
        if aa_pos < 1 or aa_pos > len(sequence):
            position_mismatch += 1
            continue
        
        # Check if wild-type matches
        actual_aa = sequence[aa_pos - 1]  # 0-indexed
        
        if actual_aa != wt_aa:
            # Allow some flexibility for non-standard amino acids or ambiguous calls
            if actual_aa not in ['X', '*', 'B', 'Z', 'J']:
                position_mismatch += 1
                continue
        
        # Generate mutant sequence
        mutant_sequence = sequence[:aa_pos-1] + mut_aa + sequence[aa_pos:]
        
        # Use the actual matched protein ID (e.g., XP_...) instead of transcript ID
        # Try to get proper protein ID from GFF mapping
        if tid in transcript_to_protein:
            proper_protein_id = transcript_to_protein[tid]
        elif matched_protein_id and matched_protein_id.startswith('XP_'):
            proper_protein_id = matched_protein_id
        elif matched_protein_id and matched_protein_id.startswith('NP_'):
            proper_protein_id = matched_protein_id
        else:
            proper_protein_id = matched_protein_id
        
        var['wt_sequence'] = sequence
        var['mut_sequence'] = mutant_sequence
        var['sequence_length'] = len(sequence)
        var['matched_protein_id'] = proper_protein_id  # Store the actual protein ID used
        
        variants_with_seqs.append(var)
    
    logger.info(f"Generated sequences for {len(variants_with_seqs)} variants")
    logger.info(f"  No sequence found: {no_sequence_count}")
    logger.info(f"  Position mismatch: {position_mismatch}")
    
    return variants_with_seqs

# ============================================================================
# MAIN
# ============================================================================

def main():
    logger.info("="*70)
    logger.info("PHASE 4: SEQUENCE PREPARATION FOR PLM PREDICTIONS")
    logger.info("="*70)
    logger.info("")
    
    # Step 1: Load protein sequences
    logger.info("Step 1: Loading protein sequences...")
    protein_sequences = load_protein_sequences(PROTEIN_FASTA)
    
    if not protein_sequences:
        logger.error(
            "No protein sequences loaded from %s. "
            "Provide the reference proteome under data/reference/ and retry.",
            PROTEIN_FASTA,
        )
        return 1

    logger.info("")
    
    # Step 2: Extract missense variants
    logger.info("Step 2: Extracting missense variants from annotated VCF...")
    missense_variants = extract_missense_variants(ANNOTATED_VCF)
    
    if not missense_variants:
        logger.error("No missense variants found!")
        return 1
    
    logger.info("")
    
    # Step 3: Generate sequences
    logger.info("Step 3: Generating wild-type and mutant sequences...")
    variants_with_seqs = generate_mutant_sequences(missense_variants, protein_sequences)
    
    if not variants_with_seqs:
        logger.error("No variants with valid sequences generated!")
        return 1
    
    logger.info("")
    
    # Step 4: Save results
    logger.info("Step 4: Saving results...")
    
    fieldnames = [
        'variant_id', 'chrom', 'pos', 'ref', 'alt',
        'gene_name', 'gene_id', 'transcript_id', 'protein_id', 'matched_protein_id',
        'effect', 'impact', 'hgvs_p',
        'aa_pos', 'wt_aa', 'mut_aa', 'sequence_length'
    ]
    
    # Save complete variant list with all isoforms
    with open(OUTPUT_VARIANTS, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        
        for var in variants_with_seqs:
            row = {k: var.get(k, '') for k in fieldnames}
            writer.writerow(row)
    
    logger.info(f"Complete variant list saved: {OUTPUT_VARIANTS}")
    logger.info(f"  Total entries (with isoforms): {len(variants_with_seqs)}")
    
    # Deduplicate variants: keep one entry per unique genomic position
    # Priority: prefer canonical/longest isoform
    unique_variants = {}
    for var in variants_with_seqs:
        var_id = var['variant_id']
        if var_id not in unique_variants:
            unique_variants[var_id] = var
        else:
            # Keep the one with longer protein sequence (usually more canonical)
            if var.get('sequence_length', 0) > unique_variants[var_id].get('sequence_length', 0):
                unique_variants[var_id] = var
    
    # Save deduplicated variant list
    with open(OUTPUT_VARIANTS_UNIQUE, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        
        for var in unique_variants.values():
            row = {k: var.get(k, '') for k in fieldnames}
            writer.writerow(row)
    
    logger.info(f"Deduplicated variant list saved: {OUTPUT_VARIANTS_UNIQUE}")
    logger.info(f"  Unique genomic positions: {len(unique_variants)}")
    
    # Save unique protein sequences using matched protein IDs
    unique_proteins = {}
    for var in variants_with_seqs:
        prot_id = var.get('matched_protein_id', var['protein_id'])
        if prot_id not in unique_proteins:
            unique_proteins[prot_id] = var['wt_sequence']
    
    with open(OUTPUT_SEQUENCES, 'w') as f:
        for prot_id, seq in unique_proteins.items():
            f.write(f">{prot_id}\n")
            # Write in 60-character lines
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')
    
    logger.info(f"Protein sequences saved: {OUTPUT_SEQUENCES}")
    logger.info(f"Unique proteins: {len(unique_proteins)}")
    
    # Summary statistics
    logger.info("")
    logger.info("="*70)
    logger.info("SUMMARY")
    logger.info("="*70)
    logger.info(f"Total variant entries (with isoforms): {len(variants_with_seqs)}")
    logger.info(f"Unique genomic variants: {len(unique_variants)}")
    logger.info(f"Unique proteins: {len(unique_proteins)}")
    logger.info(f"")
    logger.info(f"Output files:")
    logger.info(f"  - {OUTPUT_VARIANTS} (all isoforms)")
    logger.info(f"  - {OUTPUT_VARIANTS_UNIQUE} (deduplicated)")
    logger.info(f"  - {OUTPUT_SEQUENCES} (protein sequences)")
    
    # Impact distribution
    from collections import Counter
    impact_counts = Counter(var['impact'] for var in variants_with_seqs)
    logger.info(f"\nImpact distribution:")
    for impact, count in impact_counts.most_common():
        logger.info(f"  {impact}: {count}")
    
    # Sequence length statistics
    import statistics
    seq_lengths = [var['sequence_length'] for var in variants_with_seqs]
    logger.info(f"\nProtein sequence length statistics:")
    logger.info(f"  Mean: {statistics.mean(seq_lengths):.0f} aa")
    logger.info(f"  Median: {statistics.median(seq_lengths):.0f} aa")
    logger.info(f"  Min: {min(seq_lengths)} aa")
    logger.info(f"  Max: {max(seq_lengths)} aa")
    
    logger.info("")
    logger.info("Sequence preparation complete!")
    logger.info("Ready for PLM predictions")
    logger.info("="*70)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())


