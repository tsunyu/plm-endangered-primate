#!/usr/bin/env python3
"""
Phase 3: Functional Annotation Visualization

Creates publication-quality figures from functional annotation gene lists.

Usage:
    python phase2_step4.2_visualize_functional.py [options]
    
Options:
    --input-dir PATH          Directory with gene list files
    --output-dir PATH         Output directory for figures (default: same as input)
    --gene-type TYPE          Gene type to visualize: all, genes, known (default: genes)
    --dpi INT                 Figure DPI (default: 300)
    --format FORMAT           Figure format: png, pdf, svg (default: png)
    --style STYLE             Plot style: default, seaborn, ggplot (default: seaborn)
    --show-counts             Show count labels on bars
    --normalize               Normalize bar heights to percentages
    --y-scale SCALE           Y-axis scale for bars: linear or log (default: log)
    --biotype-top N           Show top N biotypes in biotype plot (default: 15)
"""

import argparse
from pathlib import Path
from collections import Counter
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import pandas as pd
import numpy as np

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.2)


def count_genes_by_impact(input_dir, gene_type='genes'):
    """Count genes in each impact category"""
    
    print(f"Counting genes from {gene_type} files...")
    
    impact_categories = ['high', 'moderate', 'low', 'modifier']
    gene_counts = {}
    
    input_dir = Path(input_dir)
    
    for impact in impact_categories:
        filename = f'genes_{impact}_impact_{gene_type}.txt'
        filepath = input_dir / filename
        
        if filepath.exists():
            with open(filepath, 'r') as f:
                genes = [line.strip() for line in f if line.strip()]
                gene_counts[impact.upper()] = len(genes)
                print(f"  {impact.upper()}: {len(genes):,} genes")
        else:
            print(f"  WARNING: File not found: {filepath}")
            gene_counts[impact.upper()] = 0
    
    return gene_counts


def plot_gene_counts(gene_counts, output_dir, gene_type='genes', 
                     dpi=300, fmt='png', figsize=(10, 6), show_counts=True,
                     normalize=False, y_scale='log'):
    """Plot gene counts by impact category"""
    
    print(f"\nCreating gene count plot for {gene_type} genes...")
    
    # Order by severity
    impact_order = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']
    impacts = [imp for imp in impact_order if imp in gene_counts]
    raw_counts = [gene_counts[imp] for imp in impacts]
    if normalize:
        total = max(1, sum(raw_counts))
        counts = [100 * c / total for c in raw_counts]
        y_label = 'Percentage of Genes (%)'
    else:
        counts = raw_counts
        y_label = 'Number of Genes'
    
    # Color scheme
    colors_map = {
        'HIGH': '#d62728',      # Red
        'MODERATE': '#ff7f0e',  # Orange
        'LOW': '#2ca02c',       # Green
        'MODIFIER': '#1f77b4'   # Blue
    }
    colors = [colors_map[imp] for imp in impacts]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Bar plot
    bars = ax1.bar(impacts, counts, color=colors, alpha=0.8, edgecolor='black', linewidth=1.5)
    ax1.set_ylabel(y_label, fontsize=12, fontweight='bold')
    ax1.set_xlabel('Impact Category', fontsize=12, fontweight='bold')
    
    gene_type_label = {
        'all': 'All Entries',
        'genes': 'Real Genes (incl. LOC)',
        'known': 'Known Gene Symbols'
    }.get(gene_type, gene_type.title())
    
    ax1.set_title(f'Gene Counts by Impact - {gene_type_label}', 
                 fontsize=13, fontweight='bold')
    if y_scale == 'log' and not normalize:
        ax1.set_yscale('log')
    ax1.grid(axis='y', alpha=0.3)
    
    # Add value labels on bars
    if show_counts:
        for bar, count in zip(bars, counts):
            height = bar.get_height()
            label = f'{count:,.1f}%' if normalize else f'{int(count):,}'
            ax1.text(bar.get_x() + bar.get_width()/2., height * (1.05 if y_scale=='log' and not normalize else 1.02),
                     label,
                     ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Pie chart
    # Use legend to avoid label overlap: show only percentages on pie
    wedges, texts, autotexts = ax2.pie(
        raw_counts,
        labels=None,
        colors=colors,
        autopct='%1.1f%%',
        pctdistance=0.7,
        startangle=90,
        textprops={'fontsize': 9}
    )

    ax2.set_title(f'Gene Distribution - {gene_type_label}',
                  fontsize=13, fontweight='bold')

    # Percentage text styling
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(9)

    # Add legend with labels and counts to avoid overlapping wedge labels
    legend_labels = [f"{label}: {count:,}" for label, count in zip(impacts, raw_counts)]
    ax2.legend(wedges, legend_labels,
               title='Impacts',
               loc='center left',
               bbox_to_anchor=(1.02, 0.5),
               fontsize=9,
               frameon=False)
    
    plt.tight_layout()
    
    output_file = Path(output_dir) / f'gene_counts_{gene_type}.{fmt}'
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_file}")
    return output_file


def plot_comparison(input_dir, output_dir, dpi=300, fmt='png', figsize=(14, 6)):
    """Compare gene counts across all three gene types"""
    
    print("\nCreating comparison plot across gene types...")
    
    gene_types = ['all', 'genes', 'known']
    impact_order = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']
    
    # Collect data
    data = {gt: count_genes_by_impact(input_dir, gt) for gt in gene_types}
    
    # Prepare data for grouped bar plot
    gene_type_labels = {
        'all': 'All Entries',
        'genes': 'Real Genes\n(incl. LOC)',
        'known': 'Known Symbols'
    }
    
    fig, ax = plt.subplots(figsize=figsize)
    
    x = np.arange(len(impact_order))
    width = 0.25
    
    colors = ['#8c564b', '#17becf', '#bcbd22']
    
    for i, gene_type in enumerate(gene_types):
        counts = [data[gene_type].get(imp, 0) for imp in impact_order]
        offset = width * (i - 1)
        bars = ax.bar(x + offset, counts, width, 
                     label=gene_type_labels[gene_type],
                     color=colors[i], alpha=0.8, edgecolor='black')
        
        # Add value labels
        for bar, count in zip(bars, counts):
            if count > 0:
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{count:,}',
                       ha='center', va='bottom', fontsize=8, rotation=90)
    
    ax.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
    ax.set_xlabel('Impact Category', fontsize=12, fontweight='bold')
    ax.set_title('Gene Count Comparison Across Classification Levels', 
                fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(impact_order)
    ax.set_yscale('log')
    ax.legend(loc='upper left', frameon=True, fancybox=True, shadow=True)
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    
    output_file = Path(output_dir) / f'gene_counts_comparison.{fmt}'
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_file}")
    return output_file


def plot_gene_categories(input_dir, output_dir, dpi=300, fmt='png', figsize=(10, 6)):
    """Plot distribution of gene categories from detailed annotation file"""
    
    print("\nCreating gene category distribution plot...")
    
    detail_file = Path(input_dir) / 'gene_annotations_detailed.tsv'
    
    if not detail_file.exists():
        print(f"  WARNING: Detailed annotation file not found: {detail_file}")
        return None
    
    # Read detailed annotations
    df = pd.read_csv(detail_file, sep='\t')
    
    # Count categories
    category_counts = df['Category'].value_counts()
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Bar plot
    colors = {'known': '#2ca02c', 'predicted': '#ff7f0e', 
              'feature_id': '#d62728', 'other': '#7f7f7f'}
    
    plot_colors = [colors.get(cat, '#7f7f7f') for cat in category_counts.index]
    
    bars = ax1.bar(range(len(category_counts)), category_counts.values, 
                   color=plot_colors, alpha=0.8, edgecolor='black')
    ax1.set_xticks(range(len(category_counts)))
    ax1.set_xticklabels(category_counts.index, rotation=45, ha='right')
    ax1.set_ylabel('Number of Entries', fontsize=12, fontweight='bold')
    ax1.set_title('Gene Entry Categories', fontsize=13, fontweight='bold')
    ax1.set_yscale('log')
    
    # Add value labels
    for bar, count in zip(bars, category_counts.values):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height * 1.1,
                f'{count:,}',
                ha='center', va='bottom', fontsize=10)
    
    # Pie chart
    wedges, texts, autotexts = ax2.pie(
        category_counts.values,
        labels=category_counts.index,
        colors=plot_colors,
        autopct='%1.1f%%',
        startangle=90,
        textprops={'fontsize': 10}
    )
    
    ax2.set_title('Category Distribution', fontsize=13, fontweight='bold')
    
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
    
    plt.tight_layout()
    
    output_file = Path(output_dir) / f'gene_categories.{fmt}'
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_file}")
    return output_file


def plot_biotype_distribution(input_dir, output_dir, top_n=15, dpi=300, fmt='png', figsize=(10, 6)):
    """Plot distribution of genes by biotype using gene_annotations_detailed.tsv or biotype_counts.tsv"""
    print("\nCreating biotype distribution plot...")

    detail_file = Path(input_dir) / 'gene_annotations_detailed.tsv'
    counts_file = Path(input_dir) / 'biotype_counts.tsv'

    df = None
    if detail_file.exists():
        df = pd.read_csv(detail_file, sep='\t')
        # Each row may contain multiple biotypes comma-separated; explode them
        if 'Biotypes' in df.columns:
            df = df[['Gene', 'Biotypes']].copy()
            df['Biotypes'] = df['Biotypes'].fillna('')
            df = df.assign(Biotype=df['Biotypes'].str.split(',')).explode('Biotype')
            df['Biotype'] = df['Biotype'].str.strip()
            df = df[df['Biotype'] != '']
            biotype_counts = df['Biotype'].value_counts()
        else:
            df = None

    if df is None and counts_file.exists():
        cdf = pd.read_csv(counts_file, sep='\t')
        biotype_counts = pd.Series(cdf['Total'].values, index=cdf['Biotype'])

    if df is None and not counts_file.exists():
        print(f"  WARNING: No biotype source file found: {detail_file} or {counts_file}")
        return None

    # Take top N biotypes
    biotype_counts = biotype_counts.sort_values(ascending=False).head(top_n)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # Horizontal bar chart
    y = np.arange(len(biotype_counts))
    bars = ax1.barh(y, biotype_counts.values, color='#9467bd', alpha=0.85, edgecolor='black')
    ax1.set_yticks(y)
    ax1.set_yticklabels(biotype_counts.index, fontsize=9)
    ax1.invert_yaxis()
    ax1.set_xlabel('Number of Genes', fontsize=12, fontweight='bold')
    ax1.set_title(f'Top {top_n} Biotypes by Gene Count', fontsize=13, fontweight='bold')
    for i, (bar, val) in enumerate(zip(bars, biotype_counts.values)):
        ax1.text(val, i, f' {int(val):,}', va='center', fontsize=9)

    # Pie chart with legend to avoid overlap
    wedges, texts, autotexts = ax2.pie(
        biotype_counts.values,
        labels=None,
        autopct='%1.1f%%',
        startangle=90,
        textprops={'fontsize': 9}
    )
    ax2.set_title('Biotype Proportions (Top)', fontsize=13, fontweight='bold')
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(9)
    legend_labels = [f"{bt}: {int(cnt):,}" for bt, cnt in biotype_counts.items()]
    ax2.legend(wedges, legend_labels, title='Biotypes', loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=9, frameon=False)

    plt.tight_layout()
    output_file = Path(output_dir) / f'gene_biotypes_top{top_n}.{fmt}'
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file}")
    return output_file


def main():
    parser = argparse.ArgumentParser(
        description='Visualize functional annotation gene lists',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument('--input-dir',
                       default=str(get_base_dir() / "output/phase2_annotation/functional_annotation"),
                       help='Directory with gene list files')
    parser.add_argument('--output-dir',
                       help='Output directory for figures (default: same as input)')
    parser.add_argument('--gene-type', default='genes',
                       choices=['all', 'genes', 'known'],
                       help='Gene type to visualize (default: genes)')
    parser.add_argument('--dpi', type=int, default=300,
                       help='Figure DPI (default: 300)')
    parser.add_argument('--format', default='png', choices=['png', 'pdf', 'svg'],
                       help='Figure format (default: png)')
    parser.add_argument('--style', default='seaborn', choices=['default', 'seaborn', 'ggplot'],
                       help='Plot style (default: seaborn)')
    parser.add_argument('--show-counts', action='store_true',
                       help='Show count labels on bars')
    parser.add_argument('--normalize', action='store_true',
                       help='Normalize bar heights to percentages')
    parser.add_argument('--y-scale', default='log', choices=['linear', 'log'],
                       help='Y-axis scale for bars (default: log)')
    parser.add_argument('--all-types', action='store_true',
                       help='Generate plots for all gene types')
    parser.add_argument('--comparison', action='store_true',
                       help='Generate comparison plot across gene types')
    parser.add_argument('--biotype-top', type=int, default=15,
                       help='Number of top biotypes to display (default: 15)')
    
    args = parser.parse_args()
    
    # Set output directory
    output_dir = Path(args.output_dir) if args.output_dir else Path(args.input_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Set plot style
    if args.style == 'seaborn':
        sns.set_style("whitegrid")
    elif args.style == 'ggplot':
        plt.style.use('ggplot')
    
    print("="*70)
    print("Functional Annotation Visualization")
    print("="*70)
    print(f"\nInput directory: {args.input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Figure format: {args.format} (DPI: {args.dpi})")
    print()
    
    # Check input directory exists
    if not Path(args.input_dir).exists():
        print(f"ERROR: Input directory not found: {args.input_dir}")
        sys.exit(1)
    
    # Generate plots
    if args.all_types:
        print("Generating plots for all gene types...")
        for gene_type in ['all', 'genes', 'known']:
            print(f"\nProcessing {gene_type}...")
            gene_counts = count_genes_by_impact(args.input_dir, gene_type)
            plot_gene_counts(gene_counts, output_dir, gene_type=gene_type,
                           dpi=args.dpi, fmt=args.format, show_counts=args.show_counts,
                           normalize=args.normalize, y_scale=args.y_scale)
    else:
        gene_counts = count_genes_by_impact(args.input_dir, args.gene_type)
        plot_gene_counts(gene_counts, output_dir, gene_type=args.gene_type,
                        dpi=args.dpi, fmt=args.format, show_counts=args.show_counts,
                        normalize=args.normalize, y_scale=args.y_scale)
    
    # Comparison plot
    if args.comparison or args.all_types:
        plot_comparison(args.input_dir, output_dir, dpi=args.dpi, fmt=args.format)
    
    # Category distribution
    plot_gene_categories(args.input_dir, output_dir, dpi=args.dpi, fmt=args.format)

    # Biotype distribution
    plot_biotype_distribution(args.input_dir, output_dir, top_n=args.biotype_top, dpi=args.dpi, fmt=args.format)
    
    print("\n" + "="*70)
    print("Visualization complete!")
    print("="*70)
    print(f"\nOutput files saved to: {output_dir}")
    print()


if __name__ == '__main__':
    main()

