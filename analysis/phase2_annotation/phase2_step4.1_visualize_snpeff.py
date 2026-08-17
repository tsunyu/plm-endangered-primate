#!/usr/bin/env python3
"""
Phase 3: SnpEff Annotation Visualization

Creates publication-quality figures from SnpEff annotation results.

Usage:
    python phase2_step4.1_visualize_snpeff.py [options]
    
Options:
    --input-vcf PATH          Path to annotated VCF file
    --output-dir PATH         Output directory for figures
    --dpi INT                 Figure DPI (default: 300)
    --format FORMAT           Figure format: png, pdf, svg (default: png)
    --style STYLE             Plot style: default, seaborn, ggplot (default: seaborn)
    --figsize W H             Figure size in inches (default: 10 6)
    --top-effects INT         Number of top effect types to show (default: 20)
    --top-genes INT           Number of top genes to show (default: 30)
    --normalize               Normalize bar heights to percentages
    --y-scale SCALE           Y-axis scale for bars: linear or log (default: linear)
"""

import gzip
import argparse
from collections import Counter, defaultdict
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import seaborn as sns
import pandas as pd
import numpy as np

# Set default style
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.2)


def parse_snpeff_vcf(vcf_file):
    """Parse SnpEff annotated VCF and extract annotation statistics"""
    
    print(f"Parsing VCF file: {vcf_file}")
    
    impact_counts = Counter()
    effect_counts = Counter()
    lof_counts = Counter()
    gene_impacts = defaultdict(list)
    high_impact_genes = Counter()
    
    variant_count = 0
    
    with gzip.open(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            variant_count += 1
            if variant_count % 100000 == 0:
                print(f"  Processed {variant_count:,} variants...")
            
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
                            gene_impacts[gene].append(impact)
                            high_impact_genes[gene] += 1
            
            # Check for LOF
            if 'LOF' in info_dict:
                lof_counts['total'] += 1
    
    print(f"Total variants processed: {variant_count:,}")
    
    return {
        'impact_counts': impact_counts,
        'effect_counts': effect_counts,
        'lof_counts': lof_counts,
        'gene_impacts': gene_impacts,
        'high_impact_genes': high_impact_genes,
        'variant_count': variant_count
    }


def plot_impact_distribution(impact_counts, output_dir, dpi=300, fmt='png', figsize=(10, 6), normalize=False, y_scale='linear'):
    """Plot variant impact distribution"""
    
    print("Creating impact distribution plot...")
    
    # Order impacts by severity
    impact_order = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']
    impacts = [imp for imp in impact_order if imp in impact_counts]
    counts_raw = [impact_counts[imp] for imp in impacts]
    
    # Normalize to percentages if requested
    if normalize:
        total = max(1, sum(counts_raw))
        counts = [100 * c / total for c in counts_raw]
        y_label = 'Percentage of Variants (%)'
    else:
        counts = counts_raw
        y_label = 'Number of Variants'
    
    # Color scheme
    colors = {
        'HIGH': '#d62728',      # Red
        'MODERATE': '#ff7f0e',  # Orange
        'LOW': '#2ca02c',       # Green
        'MODIFIER': '#1f77b4'   # Blue
    }
    plot_colors = [colors[imp] for imp in impacts]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Bar plot
    bars = ax1.bar(impacts, counts, color=plot_colors, alpha=0.8, edgecolor='black')
    ax1.set_ylabel(y_label, fontsize=12, fontweight='bold')
    ax1.set_xlabel('Impact Category', fontsize=12, fontweight='bold')
    ax1.set_title('Variant Impact Distribution', fontsize=14, fontweight='bold')
    ax1.tick_params(axis='x', rotation=0)
    
    # Add value labels on bars
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        label = f'{count:,.1f}%' if normalize else f'{int(count):,}'
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                 label,
                 ha='center', va='bottom', fontsize=10)

    # Y scale handling for skewed distributions
    if y_scale == 'log' and not normalize:
        ax1.set_yscale('log')
    
    # Pie chart
    # For pie chart, always use proportions from raw counts
    total_pie = max(1, sum(counts_raw))
    wedges, texts, autotexts = ax2.pie(counts_raw, labels=impacts, colors=plot_colors,
                                        autopct='%1.1f%%', startangle=90,
                                        textprops={'fontsize': 11})
    ax2.set_title('Impact Distribution (%)', fontsize=14, fontweight='bold')
    
    # Make percentage text bold
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
    
    plt.tight_layout()
    
    output_file = Path(output_dir) / f'impact_distribution.{fmt}'
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_file}")
    return output_file


def plot_top_effects(effect_counts, output_dir, top_n=20, dpi=300, fmt='png', figsize=(12, 8)):
    """Plot top variant effect types"""
    
    print(f"Creating top {top_n} effect types plot...")
    
    # Get top effects
    top_effects = effect_counts.most_common(top_n)
    effects, counts = zip(*top_effects)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create horizontal bar plot
    y_pos = np.arange(len(effects))
    bars = ax.barh(y_pos, counts, color='steelblue', alpha=0.8, edgecolor='black')
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(effects, fontsize=10)
    ax.set_xlabel('Number of Variants', fontsize=12, fontweight='bold')
    ax.set_title(f'Top {top_n} Variant Effect Types', fontsize=14, fontweight='bold')
    ax.invert_yaxis()  # Highest at top
    
    # Add value labels
    for i, (bar, count) in enumerate(zip(bars, counts)):
        width = bar.get_width()
        ax.text(width, i, f' {count:,}',
                ha='left', va='center', fontsize=9)
    
    plt.tight_layout()
    
    output_file = Path(output_dir) / f'top_effect_types.{fmt}'
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_file}")
    return output_file


def plot_high_impact_genes(high_impact_genes, output_dir, top_n=30, dpi=300, fmt='png', figsize=(12, 10)):
    """Plot genes with most HIGH impact variants"""
    
    print(f"Creating top {top_n} high impact genes plot...")
    
    if not high_impact_genes:
        print("  Warning: No high impact genes found")
        return None
    
    # Get top genes
    top_genes = high_impact_genes.most_common(top_n)
    genes, counts = zip(*top_genes)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create horizontal bar plot
    y_pos = np.arange(len(genes))
    bars = ax.barh(y_pos, counts, color='crimson', alpha=0.8, edgecolor='black')
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(genes, fontsize=9)
    ax.set_xlabel('Number of HIGH Impact Variants', fontsize=12, fontweight='bold')
    ax.set_title(f'Top {top_n} Genes with HIGH Impact Variants', fontsize=14, fontweight='bold')
    ax.invert_yaxis()
    
    # Add value labels
    for i, (bar, count) in enumerate(zip(bars, counts)):
        width = bar.get_width()
        ax.text(width, i, f' {count}',
                ha='left', va='center', fontsize=9)
    
    plt.tight_layout()
    
    output_file = Path(output_dir) / f'high_impact_genes_top{top_n}.{fmt}'
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_file}")
    return output_file


def plot_summary_statistics(stats, output_dir, dpi=300, fmt='png', figsize=(14, 6)):
    """Create summary statistics figure"""
    
    print("Creating summary statistics plot...")
    
    impact_counts = stats['impact_counts']
    effect_counts = stats['effect_counts']
    lof_count = stats['lof_counts'].get('total', 0)
    high_impact_genes_count = len(stats['gene_impacts'])
    
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    
    # Panel 1: Impact summary
    impact_order = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']
    impacts = [imp for imp in impact_order if imp in impact_counts]
    counts = [impact_counts[imp] for imp in impacts]
    
    colors_map = {'HIGH': '#d62728', 'MODERATE': '#ff7f0e', 
                  'LOW': '#2ca02c', 'MODIFIER': '#1f77b4'}
    colors = [colors_map[imp] for imp in impacts]
    
    axes[0].bar(range(len(impacts)), counts, color=colors, alpha=0.8, edgecolor='black')
    axes[0].set_xticks(range(len(impacts)))
    axes[0].set_xticklabels(impacts, rotation=45, ha='right')
    axes[0].set_ylabel('Count', fontweight='bold')
    axes[0].set_title('Impact Categories', fontweight='bold')
    axes[0].set_yscale('log')
    
    # Panel 2: Key metrics
    metrics = ['Total\nVariants', 'Effect\nTypes', 'HIGH Impact\nGenes', 'LOF\nVariants']
    values = [stats['variant_count'], len(effect_counts), high_impact_genes_count, lof_count]
    
    bars = axes[1].bar(range(len(metrics)), values, color='teal', alpha=0.8, edgecolor='black')
    axes[1].set_xticks(range(len(metrics)))
    axes[1].set_xticklabels(metrics, fontsize=9)
    axes[1].set_ylabel('Count', fontweight='bold')
    axes[1].set_title('Key Metrics', fontweight='bold')
    axes[1].set_yscale('log')
    
    for bar, val in zip(bars, values):
        height = bar.get_height()
        axes[1].text(bar.get_x() + bar.get_width()/2., height,
                    f'{val:,}', ha='center', va='bottom', fontsize=8)
    
    # Panel 3: Impact proportions
    total_impacts = sum(counts)
    proportions = [100 * c / total_impacts for c in counts]
    
    axes[2].pie(proportions, labels=impacts, colors=colors, autopct='%1.1f%%',
                startangle=90, textprops={'fontsize': 10})
    axes[2].set_title('Impact Proportions', fontweight='bold')
    
    plt.tight_layout()
    
    output_file = Path(output_dir) / f'annotation_summary.{fmt}'
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_file}")
    return output_file


def main():
    parser = argparse.ArgumentParser(
        description='Visualize SnpEff annotation results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument('--input-vcf', 
                       default=str(get_base_dir() / "output/phase2_annotation/snpeff_annotation/annotated_variants.vcf.gz"),
                       help='Path to annotated VCF file')
    parser.add_argument('--output-dir',
                       default=str(get_base_dir() / "output/phase2_annotation/snpeff_annotation"),
                       help='Output directory for figures')
    parser.add_argument('--dpi', type=int, default=300,
                       help='Figure DPI (default: 300)')
    parser.add_argument('--format', default='png', choices=['png', 'pdf', 'svg'],
                       help='Figure format (default: png)')
    parser.add_argument('--style', default='seaborn', choices=['default', 'seaborn', 'ggplot'],
                       help='Plot style (default: seaborn)')
    parser.add_argument('--figsize', nargs=2, type=float, default=[10, 6],
                       metavar=('WIDTH', 'HEIGHT'),
                       help='Figure size in inches (default: 10 6)')
    parser.add_argument('--top-effects', type=int, default=20,
                       help='Number of top effect types to show (default: 20)')
    parser.add_argument('--top-genes', type=int, default=30,
                       help='Number of top HIGH impact genes to show (default: 30)')
    parser.add_argument('--normalize', action='store_true',
                       help='Normalize bar heights to percentages')
    parser.add_argument('--y-scale', default='linear', choices=['linear', 'log'],
                       help='Y-axis scale for bars (default: linear)')
    
    args = parser.parse_args()
    
    # Set plot style
    if args.style == 'seaborn':
        sns.set_style("whitegrid")
    elif args.style == 'ggplot':
        plt.style.use('ggplot')
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*70)
    print("SnpEff Annotation Visualization")
    print("="*70)
    print(f"\nInput VCF: {args.input_vcf}")
    print(f"Output directory: {output_dir}")
    print(f"Figure format: {args.format} (DPI: {args.dpi})")
    print()
    
    # Check input file exists
    if not Path(args.input_vcf).exists():
        print(f"ERROR: Input VCF file not found: {args.input_vcf}")
        sys.exit(1)
    
    # Parse VCF
    stats = parse_snpeff_vcf(args.input_vcf)
    
    print("\n" + "="*70)
    print("Generating figures...")
    print("="*70)
    print()
    
    # Generate plots
    figsize = tuple(args.figsize)
    
    plot_impact_distribution(
        stats['impact_counts'], 
        output_dir, 
        dpi=args.dpi, 
        fmt=args.format,
        figsize=figsize,
        normalize=args.normalize,
        y_scale=args.y_scale
    )
    
    plot_top_effects(
        stats['effect_counts'],
        output_dir,
        top_n=args.top_effects,
        dpi=args.dpi,
        fmt=args.format,
        figsize=(12, 8)
    )
    
    plot_high_impact_genes(
        stats['high_impact_genes'],
        output_dir,
        top_n=args.top_genes,
        dpi=args.dpi,
        fmt=args.format,
        figsize=(12, 10)
    )
    
    plot_summary_statistics(
        stats,
        output_dir,
        dpi=args.dpi,
        fmt=args.format,
        figsize=(14, 6)
    )
    
    print("\n" + "="*70)
    print("Visualization complete!")
    print("="*70)
    print(f"\nOutput files saved to: {output_dir}")
    print(f"  - impact_distribution.{args.format}")
    print(f"  - top_effect_types.{args.format}")
    print(f"  - high_impact_genes_top{args.top_genes}.{args.format}")
    print(f"  - annotation_summary.{args.format}")
    print()


if __name__ == '__main__':
    main()

