#!/usr/bin/env python3
"""
Extract Gene and Protein ID/Name Mapping from GFF File

Parses GFF annotation file to extract comprehensive mappings between:
- Gene ID and Gene Name
- Transcript ID and Protein ID
- Gene → Transcript → Protein relationships
- Protein names and descriptions

Usage:
    python3 phase2_prep_extract_gene_protein_mapping.py [OPTIONS]

Options:
    --gff FILE              GFF file path (default: reference/Rhinopithecus_roxellana.gff)
    --output DIR            Output directory (default: output/gene_protein_mapping)
    --format FORMAT         Output format: csv, tsv, or both (default: both)
    --help, -h              Show this help message

Output files:
    gene_info.csv                    : Gene ID, Gene Name
    transcript_protein_mapping.csv   : Transcript ID → Protein ID
    gene_transcript_protein.csv      : Complete Gene → Transcript → Protein mapping
    protein_info.csv                 : Protein ID, Protein Name, Description
    complete_annotation.csv          : All information combined

Examples:
    # Extract with defaults
    python3 phase2_prep_extract_gene_protein_mapping.py
    
    # Specify custom GFF file
    python3 phase2_prep_extract_gene_protein_mapping.py --gff /path/to/annotation.gff
    
    # Output only CSV format
    python3 phase2_prep_extract_gene_protein_mapping.py --format csv
"""

import sys
import os
import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils import setup_logger

# ============================================================================
# SETUP
# ============================================================================

BASE_DIR = get_base_dir()
DEFAULT_GFF = f"{BASE_DIR}/data/reference/Rhinopithecus_roxellana.gff"
DEFAULT_OUTPUT_DIR = f"{BASE_DIR}/output/gene_protein_mapping"

# ============================================================================
# FUNCTIONS
# ============================================================================

def parse_gff_attributes(attr_string):
    """
    Parse GFF attribute string into dictionary
    
    Args:
        attr_string: Attribute string from GFF (column 9)
    
    Returns:
        dict: Parsed attributes
    """
    attrs = {}
    for item in attr_string.split(';'):
        item = item.strip()
        if '=' in item:
            key, value = item.split('=', 1)
            attrs[key] = value
    return attrs


def extract_gene_info(gff_file, logger):
    """
    Extract gene information (ID, Name, Description)
    
    Returns:
        dict: {gene_id: {'name': str, 'description': str, 'chrom': str}}
    """
    logger.info("Extracting gene information...")
    
    genes = {}
    opener = gzip.open if gff_file.endswith('.gz') else open
    
    with opener(gff_file, 'rt') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue
            
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            
            chrom, source, feature_type, start, end, score, strand, phase, attributes = parts
            
            # Extract gene features
            if feature_type == 'gene':
                attrs = parse_gff_attributes(attributes)
                
                gene_id = attrs.get('ID', '').replace('gene-', '')
                gene_name = attrs.get('Name', attrs.get('gene', ''))
                description = attrs.get('description', attrs.get('product', ''))
                
                if gene_id:
                    genes[gene_id] = {
                        'name': gene_name or gene_id,
                        'description': description,
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand
                    }
            
            if line_num % 100000 == 0:
                logger.info(f"  Processed {line_num:,} lines, found {len(genes):,} genes")
    
    logger.info(f"Extracted {len(genes):,} genes")
    return genes


def extract_transcript_protein_mapping(gff_file, logger):
    """
    Extract transcript → protein ID mapping from CDS features
    
    Returns:
        dict: {transcript_id: protein_id}
    """
    logger.info("Extracting transcript → protein mapping...")
    
    mapping = {}
    opener = gzip.open if gff_file.endswith('.gz') else open
    
    with opener(gff_file, 'rt') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue
            
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            
            feature_type = parts[2]
            
            # Extract from CDS features
            if feature_type == 'CDS':
                attrs = parse_gff_attributes(parts[8])
                
                # Get transcript ID from Parent
                parent = attrs.get('Parent', '')
                transcript_id = parent.replace('rna-', '').replace('transcript:', '')
                
                # Get protein ID
                protein_id = attrs.get('protein_id', '')
                
                if transcript_id and protein_id:
                    mapping[transcript_id] = protein_id
            
            if line_num % 100000 == 0:
                logger.info(f"  Processed {line_num:,} lines, found {len(mapping):,} mappings")
    
    logger.info(f"Extracted {len(mapping):,} transcript → protein mappings")
    return mapping


def extract_protein_info(gff_file, logger):
    """
    Extract protein information (ID, Name, Product/Description)
    
    Returns:
        dict: {protein_id: {'name': str, 'product': str}}
    """
    logger.info("Extracting protein information...")
    
    proteins = {}
    opener = gzip.open if gff_file.endswith('.gz') else open
    
    with opener(gff_file, 'rt') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue
            
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            
            feature_type = parts[2]
            
            # Extract from CDS features
            if feature_type == 'CDS':
                attrs = parse_gff_attributes(parts[8])
                
                protein_id = attrs.get('protein_id', '')
                protein_name = attrs.get('Name', protein_id)
                product = attrs.get('product', '')
                gene = attrs.get('gene', '')
                
                if protein_id:
                    # Keep first occurrence (most complete info usually)
                    if protein_id not in proteins:
                        proteins[protein_id] = {
                            'name': protein_name,
                            'product': product,
                            'gene': gene
                        }
            
            if line_num % 100000 == 0:
                logger.info(f"  Processed {line_num:,} lines, found {len(proteins):,} proteins")
    
    logger.info(f"Extracted {len(proteins):,} proteins")
    return proteins


def extract_gene_transcript_relationships(gff_file, logger):
    """
    Extract gene → transcript relationships
    
    Returns:
        dict: {gene_id: [transcript_ids]}
    """
    logger.info("Extracting gene → transcript relationships...")
    
    relationships = defaultdict(list)
    opener = gzip.open if gff_file.endswith('.gz') else open
    
    with opener(gff_file, 'rt') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue
            
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            
            feature_type = parts[2]
            
            # Extract from mRNA/transcript features
            if feature_type in ['mRNA', 'transcript']:
                attrs = parse_gff_attributes(parts[8])
                
                transcript_id = attrs.get('ID', '').replace('rna-', '').replace('transcript:', '')
                parent = attrs.get('Parent', '').replace('gene-', '')
                
                if transcript_id and parent:
                    relationships[parent].append(transcript_id)
            
            if line_num % 100000 == 0:
                logger.info(f"  Processed {line_num:,} lines")
    
    total_transcripts = sum(len(v) for v in relationships.values())
    logger.info(f"Extracted {len(relationships):,} genes with {total_transcripts:,} transcripts")
    return relationships


def save_gene_info(genes, output_dir, output_format, logger):
    """Save gene information table"""
    logger.info("Saving gene information...")
    
    formats = [output_format] if output_format != 'both' else ['csv', 'tsv']
    
    for fmt in formats:
        delimiter = ',' if fmt == 'csv' else '\t'
        output_file = f"{output_dir}/gene_info.{fmt}"
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter=delimiter)
            writer.writerow(['gene_id', 'gene_name', 'description', 'chromosome', 'start', 'end', 'strand'])
            
            for gene_id in sorted(genes.keys()):
                info = genes[gene_id]
                writer.writerow([
                    gene_id,
                    info['name'],
                    info['description'],
                    info['chrom'],
                    info['start'],
                    info['end'],
                    info['strand']
                ])
        
        logger.info(f"  Saved: {output_file}")


def save_transcript_protein_mapping(mapping, output_dir, output_format, logger):
    """Save transcript → protein mapping table"""
    logger.info("Saving transcript → protein mapping...")
    
    formats = [output_format] if output_format != 'both' else ['csv', 'tsv']
    
    for fmt in formats:
        delimiter = ',' if fmt == 'csv' else '\t'
        output_file = f"{output_dir}/transcript_protein_mapping.{fmt}"
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter=delimiter)
            writer.writerow(['transcript_id', 'protein_id'])
            
            for transcript_id in sorted(mapping.keys()):
                protein_id = mapping[transcript_id]
                writer.writerow([transcript_id, protein_id])
        
        logger.info(f"  Saved: {output_file}")


def save_protein_info(proteins, output_dir, output_format, logger):
    """Save protein information table"""
    logger.info("Saving protein information...")
    
    formats = [output_format] if output_format != 'both' else ['csv', 'tsv']
    
    for fmt in formats:
        delimiter = ',' if fmt == 'csv' else '\t'
        output_file = f"{output_dir}/protein_info.{fmt}"
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter=delimiter)
            writer.writerow(['protein_id', 'protein_name', 'product_description', 'gene_name'])
            
            for protein_id in sorted(proteins.keys()):
                info = proteins[protein_id]
                writer.writerow([
                    protein_id,
                    info['name'],
                    info['product'],
                    info['gene']
                ])
        
        logger.info(f"  Saved: {output_file}")


def save_complete_annotation(genes, transcript_protein, proteins, gene_transcripts, 
                            output_dir, output_format, logger):
    """Save complete gene → transcript → protein annotation table"""
    logger.info("Saving complete annotation table...")
    
    formats = [output_format] if output_format != 'both' else ['csv', 'tsv']
    
    for fmt in formats:
        delimiter = ',' if fmt == 'csv' else '\t'
        output_file = f"{output_dir}/gene_transcript_protein.{fmt}"
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter=delimiter)
            writer.writerow([
                'gene_id', 'gene_name', 'transcript_id', 'protein_id', 
                'protein_name', 'product_description', 'chromosome'
            ])
            
            for gene_id in sorted(genes.keys()):
                gene_info = genes[gene_id]
                transcripts = gene_transcripts.get(gene_id, [])
                
                if transcripts:
                    for transcript_id in transcripts:
                        protein_id = transcript_protein.get(transcript_id, '')
                        protein_info = proteins.get(protein_id, {'name': '', 'product': ''})
                        
                        writer.writerow([
                            gene_id,
                            gene_info['name'],
                            transcript_id,
                            protein_id,
                            protein_info.get('name', ''),
                            protein_info.get('product', ''),
                            gene_info['chrom']
                        ])
                else:
                    # Gene without transcripts
                    writer.writerow([
                        gene_id,
                        gene_info['name'],
                        '',
                        '',
                        '',
                        '',
                        gene_info['chrom']
                    ])
        
        logger.info(f"  Saved: {output_file}")


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Extract gene and protein ID/name mapping from GFF file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument('--gff', 
                       default=DEFAULT_GFF,
                       help=f'GFF file path (default: {DEFAULT_GFF})')
    
    parser.add_argument('--output', 
                       default=DEFAULT_OUTPUT_DIR,
                       help=f'Output directory (default: {DEFAULT_OUTPUT_DIR})')
    
    parser.add_argument('--format', 
                       choices=['csv', 'tsv', 'both'],
                       default='both',
                       help='Output format (default: both)')
    
    args = parser.parse_args()
    
    # Setup output directory and logging
    Path(args.output).mkdir(parents=True, exist_ok=True)
    logger = setup_logger("gene_protein_mapping", f"{args.output}/extraction.log")
    
    logger.info("="*70)
    logger.info("GENE AND PROTEIN ID/NAME MAPPING EXTRACTION")
    logger.info("="*70)
    logger.info(f"GFF file: {args.gff}")
    logger.info(f"Output directory: {args.output}")
    logger.info(f"Output format: {args.format}")
    logger.info("")
    
    # Check if GFF file exists
    if not os.path.exists(args.gff):
        logger.error(f"GFF file not found: {args.gff}")
        return 1
    
    # Extract information
    logger.info("Step 1: Extracting gene information...")
    genes = extract_gene_info(args.gff, logger)
    logger.info("")
    
    logger.info("Step 2: Extracting transcript → protein mapping...")
    transcript_protein = extract_transcript_protein_mapping(args.gff, logger)
    logger.info("")
    
    logger.info("Step 3: Extracting protein information...")
    proteins = extract_protein_info(args.gff, logger)
    logger.info("")
    
    logger.info("Step 4: Extracting gene → transcript relationships...")
    gene_transcripts = extract_gene_transcript_relationships(args.gff, logger)
    logger.info("")
    
    # Save results
    logger.info("Step 5: Saving results...")
    save_gene_info(genes, args.output, args.format, logger)
    save_transcript_protein_mapping(transcript_protein, args.output, args.format, logger)
    save_protein_info(proteins, args.output, args.format, logger)
    save_complete_annotation(genes, transcript_protein, proteins, gene_transcripts, 
                            args.output, args.format, logger)
    
    # Summary
    logger.info("")
    logger.info("="*70)
    logger.info("SUMMARY")
    logger.info("="*70)
    logger.info(f"Genes extracted: {len(genes):,}")
    logger.info(f"Transcript → Protein mappings: {len(transcript_protein):,}")
    logger.info(f"Unique proteins: {len(proteins):,}")
    logger.info(f"Genes with transcripts: {len(gene_transcripts):,}")
    
    total_transcripts = sum(len(v) for v in gene_transcripts.values())
    if gene_transcripts:
        avg_transcripts = total_transcripts / len(gene_transcripts)
        logger.info(f"Average transcripts per gene: {avg_transcripts:.2f}")
    
    logger.info("")
    logger.info("Output files:")
    formats = [args.format] if args.format != 'both' else ['csv', 'tsv']
    for fmt in formats:
        logger.info(f"  - gene_info.{fmt}")
        logger.info(f"  - transcript_protein_mapping.{fmt}")
        logger.info(f"  - protein_info.{fmt}")
        logger.info(f"  - gene_transcript_protein.{fmt}")
    
    logger.info("")
    logger.info("Extraction complete!")
    logger.info("="*70)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

