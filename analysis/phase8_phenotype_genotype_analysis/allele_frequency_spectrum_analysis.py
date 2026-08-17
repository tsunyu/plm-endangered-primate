#!/usr/bin/env python3
"""
Deleterious Allele Frequency Spectrum Analysis
===============================================

Analyzes the frequency distribution of deleterious alleles to distinguish
between inbreeding load (Hypothesis A) vs historical bottleneck (Hypothesis B).

Hypothesis A predictions:
- High-frequency deleterious alleles (20-80%) should dominate
- Few fixed deleterious alleles

Hypothesis B predictions:
- Low-frequency deleterious alleles should dominate
- Substantial proportion of fixed deleterious alleles (>10%)

Author: Analysis Pipeline
Date: 2025-12-10
"""

import os
import sys
import gzip
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir
from collections import defaultdict
import warnings

warnings.filterwarnings('ignore')

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils import setup_logger

# =============================================================================
# SETUP
# =============================================================================

BASE_DIR = get_base_dir()
OUTPUT_DIR = BASE_DIR / "output" / "allele_frequency_spectrum"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

VCF_FILE = BASE_DIR / "output/phase2_annotation/snpeff_annotation/annotated_variants.vcf.gz"
ESM2_FILE = BASE_DIR / "output/phase4_plm_predictions/esm2/esm2_predictions.csv"
ENSEMBLE_FILE = BASE_DIR / "output/phase4_plm_predictions/ensemble/ensemble_predictions.csv"
SIGMOID_K = 0.5287
SIGMOID_X0 = -6.8920
LOF_PATHOGENICITY = 0.95

logger = setup_logger("allele_frequency_spectrum", str(OUTPUT_DIR / "afs_analysis.log"))

# Plotting style
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3
plt.rcParams['grid.color'] = '#cccccc'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8

# Colorblind-friendly palette
CB_COLORS = {
    'blue': '#0072B2',
    'orange': '#E69F00',
    'vermillion': '#D55E00',
    'sky_blue': '#56B4E9',
    'green': '#009E73',
}

# =============================================================================
# VCF PARSING
# =============================================================================

def parse_vcf_line(line):
    """Parse a VCF line and extract key information."""
    fields = line.strip().split('\t')
    chrom = fields[0]
    pos = int(fields[1])
    ref = fields[3]
    alt = fields[4]
    info = fields[7]
    
    # Parse genotypes (starting from column 9)
    genotypes = fields[9:]
    
    return {
        'chrom': chrom,
        'pos': pos,
        'ref': ref,
        'alt': alt,
        'info': info,
        'genotypes': genotypes
    }


def extract_annotation(info_field):
    """Extract SnpEff annotation from INFO field."""
    annotations = {}
    
    for item in info_field.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            annotations[key] = value
    
    # Extract ANN field (SnpEff annotation)
    if 'ANN' in annotations:
        ann_fields = annotations['ANN'].split('|')
        if len(ann_fields) >= 4:
            annotations['effect'] = ann_fields[1]
            annotations['impact'] = ann_fields[2]
            annotations['gene'] = ann_fields[3]
    
    return annotations


def calculate_allele_frequency(genotypes):
    """
    Calculate derived allele frequency from genotypes.
    Assumes derived allele is ALT (will be refined with ancestral state).
    
    Genotypes: list of GT format strings (e.g., '0/0', '0/1', '1/1')
    """
    allele_counts = {'ref': 0, 'alt': 0, 'missing': 0}
    
    for gt in genotypes:
        # Extract GT field (first field in FORMAT)
        gt_field = gt.split(':')[0]
        
        if '/' in gt_field:
            alleles = gt_field.split('/')
        elif '|' in gt_field:
            alleles = gt_field.split('|')
        else:
            allele_counts['missing'] += 2
            continue
        
        for allele in alleles:
            if allele == '.':
                allele_counts['missing'] += 1
            elif allele == '0':
                allele_counts['ref'] += 1
            elif allele == '1':
                allele_counts['alt'] += 1
    
    total_alleles = allele_counts['ref'] + allele_counts['alt']
    
    if total_alleles == 0:
        return None
    
    daf = allele_counts['alt'] / total_alleles
    
    return daf


def is_lof_variant(effect):
    """Check if variant is Loss-of-Function."""
    lof_effects = [
        'stop_gained', 'stop_lost', 'start_lost', 
        'frameshift_variant', 'splice_acceptor_variant', 
        'splice_donor_variant'
    ]
    
    return any(lof in effect for lof in lof_effects)


def is_missense_variant(effect):
    """Check if variant is missense."""
    return 'missense_variant' in effect


def sigmoid_pathogenicity(esm2_score):
    exponent = SIGMOID_K * (esm2_score - SIGMOID_X0)
    return 1.0 / (1.0 + np.exp(exponent))


def load_esm2_lookup(esm2_file, ensemble_file):
    """Build variant_id -> pathogenicity probability lookup."""
    lookup = {}
    if esm2_file.exists():
        esm2_df = pd.read_csv(esm2_file, usecols=['variant_id', 'esm2_score'])
        for _, row in esm2_df.iterrows():
            score = row.get('esm2_score')
            if pd.notna(score):
                lookup[row['variant_id']] = sigmoid_pathogenicity(score)
    if ensemble_file.exists():
        ensemble_df = pd.read_csv(
            ensemble_file, usecols=['variant_id', 'is_deleterious']
        )
        deleterious_ids = set(
            ensemble_df.loc[ensemble_df['is_deleterious'], 'variant_id']
        )
    else:
        deleterious_ids = set()
    return lookup, deleterious_ids


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def analyze_frequency_spectrum():
    """
    Main function to analyze allele frequency spectrum of deleterious variants.
    """
    logger.info("=" * 70)
    logger.info("DELETERIOUS ALLELE FREQUENCY SPECTRUM ANALYSIS")
    logger.info("=" * 70)
    
    logger.info("\nLoading ESM-2 predictions...")
    esm2_lookup, deleterious_ids = load_esm2_lookup(ESM2_FILE, ENSEMBLE_FILE)
    logger.info(
        f"Loaded {len(esm2_lookup)} scored missense variants; "
        f"{len(deleterious_ids)} PLM-flagged deleterious IDs available"
    )
    
    # Initialize data structures
    lof_variants = []
    deleterious_missense = []
    neutral_missense = []
    
    logger.info(f"\nParsing VCF file: {VCF_FILE}")
    
    with gzip.open(VCF_FILE, 'rt') as f:
        for i, line in enumerate(f):
            if line.startswith('#'):
                # Store header to get sample names
                if line.startswith('#CHROM'):
                    header = line.strip().split('\t')
                    samples = header[9:]
                    n_samples = len(samples)
                    logger.info(f"Found {n_samples} samples")
                continue
            
            # Progress update
            if (i + 1) % 100000 == 0:
                logger.info(f"  Processed {i+1} variants...")
            
            # Parse variant
            variant = parse_vcf_line(line)
            annotations = extract_annotation(variant['info'])
            
            # Skip if no effect annotation
            if 'effect' not in annotations:
                continue
            
            effect = annotations['effect']
            
            # Calculate allele frequency
            daf = calculate_allele_frequency(variant['genotypes'])
            if daf is None:
                continue
            
            # Classify and store
            variant_info = {
                'chrom': variant['chrom'],
                'pos': variant['pos'],
                'ref': variant['ref'],
                'alt': variant['alt'],
                'daf': daf,
                'gene': annotations.get('gene', 'Unknown'),
                'effect': effect,
                'impact': annotations.get('impact', 'Unknown')
            }
            
            variant_id = f"{variant['chrom']}:{variant['pos']}:{variant['ref']}:{variant['alt']}"
            if is_lof_variant(effect):
                variant_info['pathogenicity_prob'] = LOF_PATHOGENICITY
                lof_variants.append(variant_info)
            elif is_missense_variant(effect):
                if variant_id not in deleterious_ids:
                    continue
                score = esm2_lookup.get(variant_id)
                if score is None:
                    continue
                variant_info['pathogenicity_prob'] = score
                deleterious_missense.append(variant_info)
    
    logger.info(f"\nTotal variants processed: {i+1}")
    logger.info(f"LoF variants: {len(lof_variants)}")
    logger.info(f"Deleterious missense variants: {len(deleterious_missense)}")
    
    # Convert to DataFrames
    lof_df = pd.DataFrame(lof_variants)
    missense_df = pd.DataFrame(deleterious_missense)
    
    # Save raw data
    if len(lof_df) > 0:
        lof_df.to_csv(OUTPUT_DIR / "lof_frequency_spectrum.csv", index=False)
    if len(missense_df) > 0:
        missense_df.to_csv(OUTPUT_DIR / "missense_frequency_spectrum.csv", index=False)
    
    return lof_df, missense_df


def calculate_frequency_statistics(lof_df, missense_df):
    """Calculate frequency spectrum statistics."""
    logger.info("\n" + "=" * 70)
    logger.info("FREQUENCY SPECTRUM STATISTICS")
    logger.info("=" * 70)
    
    results = {}
    
    for name, df in [('LoF', lof_df), ('Deleterious Missense', missense_df)]:
        if len(df) == 0:
            logger.warning(f"No {name} variants found")
            continue
        
        logger.info(f"\n{name} Variants (n={len(df)}):")
        logger.info("-" * 50)
        
        # Calculate frequency bins
        rare = (df['daf'] < 0.05).sum()
        low = ((df['daf'] >= 0.05) & (df['daf'] < 0.20)).sum()
        intermediate = ((df['daf'] >= 0.20) & (df['daf'] < 0.50)).sum()
        common = ((df['daf'] >= 0.50) & (df['daf'] < 0.95)).sum()
        fixed = (df['daf'] >= 0.95).sum()
        
        total = len(df)
        
        logger.info(f"  Rare (DAF < 0.05):           {rare:5d} ({rare/total*100:5.1f}%)")
        logger.info(f"  Low (0.05 ≤ DAF < 0.20):     {low:5d} ({low/total*100:5.1f}%)")
        logger.info(f"  Intermediate (0.20 ≤ DAF < 0.50): {intermediate:5d} ({intermediate/total*100:5.1f}%)")
        logger.info(f"  Common (0.50 ≤ DAF < 0.95):  {common:5d} ({common/total*100:5.1f}%)")
        logger.info(f"  Fixed (DAF ≥ 0.95):          {fixed:5d} ({fixed/total*100:5.1f}%)")
        
        logger.info(f"\n  Mean DAF: {df['daf'].mean():.4f}")
        logger.info(f"  Median DAF: {df['daf'].median():.4f}")
        logger.info(f"  SD DAF: {df['daf'].std():.4f}")
        
        results[name] = {
            'total': total,
            'rare': rare,
            'low': low,
            'intermediate': intermediate,
            'common': common,
            'fixed': fixed,
            'mean_daf': df['daf'].mean(),
            'median_daf': df['daf'].median(),
            'sd_daf': df['daf'].std()
        }
    
    # Save statistics
    stats_df = pd.DataFrame(results).T
    stats_df.to_csv(OUTPUT_DIR / "frequency_statistics.csv")
    
    return results


def compare_with_hypotheses(results):
    """Compare observed frequency spectrum with hypothesis predictions."""
    logger.info("\n" + "=" * 70)
    logger.info("HYPOTHESIS TESTING")
    logger.info("=" * 70)
    
    logger.info("\nHypothesis A (Recent Inbreeding Load):")
    logger.info("  Prediction: High-frequency alleles dominate (>50% in intermediate-common range)")
    logger.info("  Prediction: Few fixed alleles (<5%)")
    
    logger.info("\nHypothesis B (Historical Bottleneck / Drift Load):")
    logger.info("  Prediction: Low-frequency alleles dominate OR")
    logger.info("  Prediction: Substantial fixed alleles (>10%)")
    
    logger.info("\n" + "-" * 70)
    logger.info("OBSERVED RESULTS:")
    logger.info("-" * 70)
    
    for name, stats in results.items():
        logger.info(f"\n{name}:")
        
        # Calculate combined categories for hypothesis testing
        low_freq = stats['rare'] + stats['low']
        high_freq = stats['intermediate'] + stats['common']
        fixed_freq = stats['fixed']
        total = stats['total']
        
        low_pct = low_freq / total * 100
        high_pct = high_freq / total * 100
        fixed_pct = fixed_freq / total * 100
        
        logger.info(f"  Low-frequency (DAF < 0.20):  {low_pct:.1f}%")
        logger.info(f"  High-frequency (0.20 ≤ DAF < 0.95): {high_pct:.1f}%")
        logger.info(f"  Fixed (DAF ≥ 0.95):          {fixed_pct:.1f}%")
        
        # Interpret
        logger.info(f"\n  Interpretation:")
        if high_pct > 50:
            logger.info(f"    → High-frequency alleles dominate ({high_pct:.1f}%)")
            logger.info(f"    → Consistent with Hypothesis A (recent inbreeding)")
        else:
            logger.info(f"    → Low-frequency alleles dominate ({low_pct:.1f}%)")
            if fixed_pct > 10:
                logger.info(f"    → Substantial fixed alleles ({fixed_pct:.1f}%)")
                logger.info(f"    → STRONGLY supports Hypothesis B (drift load)")
            else:
                logger.info(f"    → Few fixed alleles ({fixed_pct:.1f}%)")
                logger.info(f"    → Mixed pattern, but leans toward Hypothesis B")


def create_visualizations(lof_df, missense_df, results):
    """Create frequency spectrum visualizations."""
    logger.info("\n" + "=" * 70)
    logger.info("GENERATING VISUALIZATIONS")
    logger.info("=" * 70)
    
    # Figure 1: Frequency spectrum histograms
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    bins = np.arange(0, 1.05, 0.05)
    
    # LoF variants
    ax = axes[0, 0]
    if len(lof_df) > 0:
        ax.hist(lof_df['daf'], bins=bins, color=CB_COLORS['vermillion'], 
                alpha=0.7, edgecolor='white', linewidth=0.5)
        ax.axvline(x=0.05, color='gray', linestyle='--', linewidth=1, alpha=0.5, label='Rare/Low boundary')
        ax.axvline(x=0.20, color='gray', linestyle='--', linewidth=1, alpha=0.5, label='Low/Intermediate')
        ax.axvline(x=0.50, color='gray', linestyle='--', linewidth=1, alpha=0.5, label='Intermediate/Common')
        ax.axvline(x=0.95, color='red', linestyle='--', linewidth=1.5, alpha=0.7, label='Fixed threshold')
        ax.set_xlabel('Derived Allele Frequency (DAF)', fontsize=11, fontweight='bold')
        ax.set_ylabel('Count', fontsize=11, fontweight='bold')
        ax.set_title(f'Loss-of-Function Variants (n={len(lof_df)})', fontsize=12, fontweight='bold')
        ax.legend(fontsize=8, loc='upper right')
        ax.grid(True, alpha=0.3)
    
    # Missense variants
    ax = axes[0, 1]
    if len(missense_df) > 0:
        ax.hist(missense_df['daf'], bins=bins, color=CB_COLORS['blue'], 
                alpha=0.7, edgecolor='white', linewidth=0.5)
        ax.axvline(x=0.05, color='gray', linestyle='--', linewidth=1, alpha=0.5)
        ax.axvline(x=0.20, color='gray', linestyle='--', linewidth=1, alpha=0.5)
        ax.axvline(x=0.50, color='gray', linestyle='--', linewidth=1, alpha=0.5)
        ax.axvline(x=0.95, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
        ax.set_xlabel('Derived Allele Frequency (DAF)', fontsize=11, fontweight='bold')
        ax.set_ylabel('Count', fontsize=11, fontweight='bold')
        ax.set_title(f'Deleterious Missense Variants (n={len(missense_df)})', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
    
    # Combined frequency categories bar plot
    ax = axes[1, 0]
    categories = ['Rare\n(<0.05)', 'Low\n(0.05-0.20)', 'Intermediate\n(0.20-0.50)', 
                  'Common\n(0.50-0.95)', 'Fixed\n(≥0.95)']
    
    if 'LoF' in results:
        stats = results['LoF']
        values = [stats['rare'], stats['low'], stats['intermediate'], 
                  stats['common'], stats['fixed']]
        percentages = [v / stats['total'] * 100 for v in values]
        
        bars = ax.bar(categories, percentages, color=CB_COLORS['vermillion'], 
                     alpha=0.7, edgecolor='white', linewidth=1.5)
        
        # Add percentage labels
        for bar, pct in zip(bars, percentages):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 1,
                   f'{pct:.1f}%', ha='center', va='bottom', fontsize=9, fontweight='bold')
        
        ax.set_ylabel('Percentage of Variants (%)', fontsize=11, fontweight='bold')
        ax.set_title('LoF Variants by Frequency Category', fontsize=12, fontweight='bold')
        ax.set_ylim(0, max(percentages) * 1.2)
        ax.grid(True, alpha=0.3, axis='y')
    
    # Missense frequency categories
    ax = axes[1, 1]
    if 'Deleterious Missense' in results:
        stats = results['Deleterious Missense']
        values = [stats['rare'], stats['low'], stats['intermediate'], 
                  stats['common'], stats['fixed']]
        percentages = [v / stats['total'] * 100 for v in values]
        
        bars = ax.bar(categories, percentages, color=CB_COLORS['blue'], 
                     alpha=0.7, edgecolor='white', linewidth=1.5)
        
        for bar, pct in zip(bars, percentages):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 1,
                   f'{pct:.1f}%', ha='center', va='bottom', fontsize=9, fontweight='bold')
        
        ax.set_ylabel('Percentage of Variants (%)', fontsize=11, fontweight='bold')
        ax.set_title('Deleterious Missense by Frequency Category', fontsize=12, fontweight='bold')
        ax.set_ylim(0, max(percentages) * 1.2)
        ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "frequency_spectrum.png", dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"  Saved: frequency_spectrum.png")
    
    # Figure 2: Hypothesis testing visualization
    fig, ax = plt.subplots(figsize=(10, 6))
    
    x_pos = np.arange(3)
    width = 0.25
    
    categories_combined = ['Low-frequency\n(DAF<0.20)', 'High-frequency\n(0.20≤DAF<0.95)', 
                          'Fixed\n(DAF≥0.95)']
    
    if 'LoF' in results:
        lof_stats = results['LoF']
        lof_total = lof_stats['total']
        lof_values = [
            (lof_stats['rare'] + lof_stats['low']) / lof_total * 100,
            (lof_stats['intermediate'] + lof_stats['common']) / lof_total * 100,
            lof_stats['fixed'] / lof_total * 100
        ]
        
        x_lof = x_pos - width / 2
        ax.bar(
            x_lof,
            lof_values,
            width,
            label='LoF',
            color=CB_COLORS['vermillion'],
            alpha=0.7,
            edgecolor='white',
            linewidth=1.5,
        )
    
    if 'Deleterious Missense' in results:
        mis_stats = results['Deleterious Missense']
        mis_total = mis_stats['total']
        mis_values = [
            (mis_stats['rare'] + mis_stats['low']) / mis_total * 100,
            (mis_stats['intermediate'] + mis_stats['common']) / mis_total * 100,
            mis_stats['fixed'] / mis_total * 100
        ]
        
        x_mis = x_pos + width / 2
        ax.bar(
            x_mis,
            mis_values,
            width,
            label='Deleterious Missense',
            color=CB_COLORS['blue'],
            alpha=0.7,
            edgecolor='white',
            linewidth=1.5,
        )
    
    # Add hypothesis prediction regions
    ax.axhspan(0, 50, alpha=0.1, color='green', label='Hypothesis B region')
    ax.axhspan(50, 100, alpha=0.1, color='red', label='Hypothesis A region')
    
    ax.set_ylabel('Percentage of Variants (%)', fontsize=11, fontweight='bold')
    ax.set_title('Frequency Spectrum: Hypothesis Testing', fontsize=13, fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(categories_combined, fontsize=10)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add interpretation text
    textstr = 'Hypothesis A: High-frequency > 50%\nHypothesis B: Low-frequency dominant OR Fixed > 10%'
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=9,
           verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "hypothesis_testing.png", dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"  Saved: hypothesis_testing.png")


def main():
    """Main execution function."""
    logger.info("Starting allele frequency spectrum analysis...")
    logger.info(f"Output directory: {OUTPUT_DIR}")
    
    # Check input files
    if not VCF_FILE.exists():
        logger.error(f"VCF file not found: {VCF_FILE}")
        return 1
    
    # Step 1: Analyze frequency spectrum
    lof_df, missense_df = analyze_frequency_spectrum()
    
    if len(lof_df) == 0 and len(missense_df) == 0:
        logger.error("No deleterious variants found!")
        return 1
    
    # Step 2: Calculate statistics
    results = calculate_frequency_statistics(lof_df, missense_df)
    
    # Step 3: Compare with hypotheses
    compare_with_hypotheses(results)
    
    # Step 4: Create visualizations
    create_visualizations(lof_df, missense_df, results)
    
    logger.info("\n" + "=" * 70)
    logger.info("ANALYSIS COMPLETE")
    logger.info("=" * 70)
    logger.info(f"\nOutput files saved to: {OUTPUT_DIR}")
    logger.info("  - lof_frequency_spectrum.csv")
    logger.info("  - missense_frequency_spectrum.csv")
    logger.info("  - frequency_statistics.csv")
    logger.info("  - frequency_spectrum.png")
    logger.info("  - hypothesis_testing.png")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
