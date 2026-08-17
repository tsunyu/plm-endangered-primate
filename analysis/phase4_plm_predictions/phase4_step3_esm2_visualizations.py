#!/usr/bin/env python3
"""
Phase 4: ESM-2 Predictions Visualization

Generates comprehensive visualizations of ESM-2 mutation effect predictions.

Requires: pandas, matplotlib, seaborn, numpy
"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils import setup_logger

# ============================================================================
# SETUP
# ============================================================================

BASE_DIR = get_base_dir()
OUTPUT_DIR = f"{BASE_DIR}/output/phase4_plm_predictions/esm2"
VISUALIZATIONS_DIR = f"{OUTPUT_DIR}/visualizations"
INPUT_FILE = f"{OUTPUT_DIR}/esm2_predictions.csv"

# Create visualizations directory
Path(VISUALIZATIONS_DIR).mkdir(parents=True, exist_ok=True)

logger = setup_logger("esm2_viz", f"{VISUALIZATIONS_DIR}/visualization.log")

# Set style
plt.style.use('seaborn-v0_8-darkgrid')

# Colorblind-friendly palette (Wong palette)
# Reference: https://www.nature.com/articles/nmeth.1618
COLORBLIND_PALETTE = {
    'blue': '#0072B2',        # Dark blue
    'orange': '#E69F00',      # Orange/Gold
    'teal': '#009E73',        # Teal/Bluish green
    'yellow': '#F0E442',      # Yellow
    'sky_blue': '#56B4E9',    # Sky blue
    'vermillion': '#D55E00',  # Vermillion/Red-orange
    'pink': '#CC79A7',        # Reddish purple
    'black': '#000000',       # Black
}

# Category colors (colorblind-safe)
CATEGORY_COLORS = {
    'deleterious': COLORBLIND_PALETTE['vermillion'],       # Vermillion for harmful
    'possibly_deleterious': COLORBLIND_PALETTE['orange'],  # Orange for possibly harmful
    'benign': COLORBLIND_PALETTE['teal']                   # Teal for safe
}

sns.set_palette([COLORBLIND_PALETTE['blue'], COLORBLIND_PALETTE['orange'], 
                 COLORBLIND_PALETTE['teal'], COLORBLIND_PALETTE['vermillion']])

# ============================================================================
# FUNCTIONS
# ============================================================================

def format_category_name(category):
    """Format category names for display (remove underscores)"""
    category_map = {
        'deleterious': 'Deleterious',
        'possibly_deleterious': 'Possibly Deleterious',
        'benign': 'Benign'
    }
    return category_map.get(category, category.replace('_', ' ').title())


def load_data():
    """Load ESM-2 predictions"""
    logger.info("Loading ESM-2 predictions...")
    df = pd.read_csv(INPUT_FILE)
    logger.info(f"Loaded {len(df)} variants")
    
    # Filter out NaN scores for analysis
    valid_df = df[~df['esm2_score'].isna()].copy()
    logger.info(f"Valid predictions: {len(valid_df)} / {len(df)}")
    
    return df, valid_df


def plot_score_distribution(df, output_dir):
    """Plot distribution of ESM-2 scores"""
    logger.info("Creating score distribution plot...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. Histogram
    ax = axes[0, 0]
    ax.hist(df['esm2_score'], bins=50, color=COLORBLIND_PALETTE['blue'], alpha=0.7, edgecolor='black')
    mean_score = df['esm2_score'].mean()
    median_score = df['esm2_score'].median()
    ax.axvline(mean_score, color=COLORBLIND_PALETTE['vermillion'], linestyle='--', linewidth=2, label=f'Mean: {mean_score:.3f}')
    ax.axvline(median_score, color=COLORBLIND_PALETTE['teal'], linestyle='-.', linewidth=2, label=f'Median: {median_score:.3f}')
    ax.axvline(-2, color=COLORBLIND_PALETTE['orange'], linestyle=':', linewidth=1.5, label='Deleterious threshold (-2)')
    ax.axvline(0, color=COLORBLIND_PALETTE['black'], linestyle=':', linewidth=1, label='Neutral (0)')
    ax.set_xlabel('ESM-2 Score (log-likelihood ratio)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Count', fontsize=11, fontweight='bold')
    ax.set_title('Distribution of ESM-2 Scores', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    # 2. Density plot
    ax = axes[0, 1]
    ax.hist(df['esm2_score'], bins=50, density=True, color=COLORBLIND_PALETTE['sky_blue'], alpha=0.5, edgecolor='black')
    df['esm2_score'].plot.density(ax=ax, color=COLORBLIND_PALETTE['blue'], linewidth=2)
    ax.axvline(-2, color=COLORBLIND_PALETTE['vermillion'], linestyle='--', linewidth=1.5, label='Deleterious (-2)')
    ax.axvline(0, color=COLORBLIND_PALETTE['black'], linestyle='-.', linewidth=1.5, label='Neutral (0)')
    ax.set_xlabel('ESM-2 Score', fontsize=11, fontweight='bold')
    ax.set_ylabel('Density', fontsize=11, fontweight='bold')
    ax.set_title('Probability Density of ESM-2 Scores', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    # 3. Box plot
    ax = axes[1, 0]
    box_data = [df['esm2_score']]
    bp = ax.boxplot(box_data, vert=True, patch_artist=True, tick_labels=['ESM-2 Scores'])
    bp['boxes'][0].set_facecolor(COLORBLIND_PALETTE['sky_blue'])
    bp['boxes'][0].set_alpha(0.7)
    ax.axhline(-2, color=COLORBLIND_PALETTE['vermillion'], linestyle='--', linewidth=1.5, label='Deleterious threshold')
    ax.axhline(0, color=COLORBLIND_PALETTE['black'], linestyle='-.', linewidth=1.5, label='Neutral')
    ax.set_ylabel('ESM-2 Score', fontsize=11, fontweight='bold')
    ax.set_title('Box Plot of ESM-2 Scores', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')
    
    # 4. Cumulative distribution
    ax = axes[1, 1]
    sorted_scores = np.sort(df['esm2_score'])
    cumulative = np.arange(1, len(sorted_scores) + 1) / len(sorted_scores)
    ax.plot(sorted_scores, cumulative, linewidth=2, color=COLORBLIND_PALETTE['blue'])
    ax.axvline(-2, color=COLORBLIND_PALETTE['vermillion'], linestyle='--', linewidth=1.5, label='Deleterious (-2)')
    ax.axvline(0, color=COLORBLIND_PALETTE['black'], linestyle='-.', linewidth=1.5, label='Neutral (0)')
    ax.set_xlabel('ESM-2 Score', fontsize=11, fontweight='bold')
    ax.set_ylabel('Cumulative Probability', fontsize=11, fontweight='bold')
    ax.set_title('Cumulative Distribution of ESM-2 Scores', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = f"{output_dir}/esm2_score_distribution.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved: {output_file}")


def plot_prediction_categories(df, output_dir):
    """Plot distribution of prediction categories"""
    logger.info("Creating prediction categories plot...")
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # 1. Bar plot
    ax = axes[0]
    category_counts = df['esm2_prediction'].value_counts()
    bar_colors = [CATEGORY_COLORS.get(cat, COLORBLIND_PALETTE['black']) for cat in category_counts.index]
    
    # Format category names for display (remove underscores)
    formatted_labels = [format_category_name(cat) for cat in category_counts.index]
    
    bars = ax.bar(formatted_labels, category_counts.values, color=bar_colors, alpha=0.7, edgecolor='black')
    ax.set_xlabel('Prediction Category', fontsize=11, fontweight='bold')
    ax.set_ylabel('Number of Variants', fontsize=11, fontweight='bold')
    ax.set_title('Distribution of Prediction Categories', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    # Rotate x labels slightly for better readability
    ax.tick_params(axis='x', rotation=15)
    
    # Add count labels on bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}\n({height/len(df)*100:.1f}%)',
                ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    # 2. Pie chart
    ax = axes[1]
    wedges, texts, autotexts = ax.pie(category_counts.values, labels=formatted_labels, 
                                       autopct='%1.1f%%', colors=bar_colors, startangle=90)
    ax.set_title('Prediction Categories (Pie Chart)', fontsize=13, fontweight='bold')
    
    # Make percentage text bold
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(10)
    
    plt.tight_layout()
    output_file = f"{output_dir}/esm2_prediction_categories.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved: {output_file}")


def plot_percentile_distribution(df, output_dir):
    """Plot percentile distribution"""
    logger.info("Creating percentile distribution plot...")
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.hist(df['esm2_percentile'], bins=50, color=COLORBLIND_PALETTE['pink'], alpha=0.7, edgecolor='black')
    ax.axvline(0.1, color=COLORBLIND_PALETTE['vermillion'], linestyle='--', linewidth=2, label='Top 10% (most deleterious)')
    ax.axvline(0.5, color=COLORBLIND_PALETTE['black'], linestyle='-.', linewidth=1.5, label='Median (50%)')
    ax.axvline(0.9, color=COLORBLIND_PALETTE['blue'], linestyle=':', linewidth=2, label='Top 90%')
    ax.set_xlabel('ESM-2 Percentile', fontsize=11, fontweight='bold')
    ax.set_ylabel('Count', fontsize=11, fontweight='bold')
    ax.set_title('Distribution of ESM-2 Percentiles', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = f"{output_dir}/esm2_percentile_distribution.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved: {output_file}")


def plot_score_vs_percentile(df, output_dir):
    """Plot score vs percentile relationship"""
    logger.info("Creating score vs percentile plot...")
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.scatter(df['esm2_percentile'], df['esm2_score'], alpha=0.5, s=10, color=COLORBLIND_PALETTE['blue'])
    ax.axhline(-2, color=COLORBLIND_PALETTE['vermillion'], linestyle='--', linewidth=1.5, label='Deleterious threshold (-2)')
    ax.axhline(0, color=COLORBLIND_PALETTE['black'], linestyle='-.', linewidth=1.5, label='Neutral (0)')
    ax.axvline(0.1, color=COLORBLIND_PALETTE['orange'], linestyle=':', linewidth=1.5, label='Top 10%')
    ax.set_xlabel('ESM-2 Percentile', fontsize=11, fontweight='bold')
    ax.set_ylabel('ESM-2 Score', fontsize=11, fontweight='bold')
    ax.set_title('ESM-2 Score vs Percentile', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = f"{output_dir}/esm2_score_vs_percentile.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved: {output_file}")


def plot_summary_statistics(df, valid_df, output_dir):
    """Create summary statistics visualization"""
    logger.info("Creating summary statistics plot...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Summary statistics
    stats = {
        'Total Variants': len(df),
        'Valid Predictions': len(valid_df),
        'Deleterious': len(valid_df[valid_df['esm2_prediction'] == 'deleterious']),
        'Possibly Deleterious': len(valid_df[valid_df['esm2_prediction'] == 'possibly_deleterious']),
        'Benign': len(valid_df[valid_df['esm2_prediction'] == 'benign'])
    }
    
    # Text summary
    ax = axes[0, 0]
    ax.axis('off')
    summary_text = f"""
    ESM-2 PREDICTION SUMMARY
    
    Total Variants: {stats['Total Variants']:,}
    Valid Predictions: {stats['Valid Predictions']:,}
    
    Prediction Distribution:
    • Deleterious: {stats['Deleterious']:,} ({stats['Deleterious']/len(valid_df)*100:.1f}%)
    • Possibly Deleterious: {stats['Possibly Deleterious']:,} ({stats['Possibly Deleterious']/len(valid_df)*100:.1f}%)
    • Benign: {stats['Benign']:,} ({stats['Benign']/len(valid_df)*100:.1f}%)
    
    Score Statistics:
    • Mean Score: {valid_df['esm2_score'].mean():.3f}
    • Median Score: {valid_df['esm2_score'].median():.3f}
    • Std Dev: {valid_df['esm2_score'].std():.3f}
    • Min Score: {valid_df['esm2_score'].min():.3f}
    • Max Score: {valid_df['esm2_score'].max():.3f}
    """
    ax.text(0.1, 0.5, summary_text, fontsize=11, family='monospace',
            verticalalignment='center', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Score statistics bar chart
    ax = axes[0, 1]
    stat_names = ['Mean', 'Median', 'Std Dev']
    stat_values = [valid_df['esm2_score'].mean(), valid_df['esm2_score'].median(), valid_df['esm2_score'].std()]
    bars = ax.bar(stat_names, stat_values, color=[COLORBLIND_PALETTE['blue'], COLORBLIND_PALETTE['teal'], COLORBLIND_PALETTE['orange']], alpha=0.7, edgecolor='black')
    ax.set_ylabel('Score Value', fontsize=11, fontweight='bold')
    ax.set_title('Score Statistics', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    for bar, val in zip(bars, stat_values):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:.3f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    # Category counts
    ax = axes[1, 0]
    categories = ['Deleterious', 'Possibly Deleterious', 'Benign']
    counts = [stats['Deleterious'], stats['Possibly Deleterious'], stats['Benign']]
    colors_cat = [CATEGORY_COLORS['deleterious'], CATEGORY_COLORS['possibly_deleterious'], CATEGORY_COLORS['benign']]
    bars = ax.bar(categories, counts, color=colors_cat, alpha=0.7, edgecolor='black')
    ax.set_ylabel('Count', fontsize=11, fontweight='bold')
    ax.set_title('Prediction Category Counts', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    # Rotate x labels to fit longer category names
    ax.tick_params(axis='x', rotation=15)
    
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{count:,}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    # Percentile distribution summary
    ax = axes[1, 1]
    percentile_ranges = ['0-10%', '10-25%', '25-50%', '50-75%', '75-90%', '90-100%']
    percentile_counts = [
        len(valid_df[(valid_df['esm2_percentile'] >= 0) & (valid_df['esm2_percentile'] < 0.1)]),
        len(valid_df[(valid_df['esm2_percentile'] >= 0.1) & (valid_df['esm2_percentile'] < 0.25)]),
        len(valid_df[(valid_df['esm2_percentile'] >= 0.25) & (valid_df['esm2_percentile'] < 0.5)]),
        len(valid_df[(valid_df['esm2_percentile'] >= 0.5) & (valid_df['esm2_percentile'] < 0.75)]),
        len(valid_df[(valid_df['esm2_percentile'] >= 0.75) & (valid_df['esm2_percentile'] < 0.9)]),
        len(valid_df[(valid_df['esm2_percentile'] >= 0.9) & (valid_df['esm2_percentile'] <= 1.0)])
    ]
    bars = ax.bar(percentile_ranges, percentile_counts, color=COLORBLIND_PALETTE['teal'], alpha=0.7, edgecolor='black')
    ax.set_xlabel('Percentile Range', fontsize=11, fontweight='bold')
    ax.set_ylabel('Count', fontsize=11, fontweight='bold')
    ax.set_title('Variants by Percentile Range', fontsize=13, fontweight='bold')
    ax.tick_params(axis='x', rotation=45)
    ax.grid(True, alpha=0.3, axis='y')
    
    for bar, count in zip(bars, percentile_counts):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{count:,}', ha='center', va='bottom', fontsize=8, fontweight='bold')
    
    plt.tight_layout()
    output_file = f"{output_dir}/esm2_summary_statistics.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved: {output_file}")


def main():
    logger.info("="*70)
    logger.info("ESM-2 PREDICTIONS VISUALIZATION")
    logger.info("="*70)
    logger.info("")
    
    # Load data
    df, valid_df = load_data()
    
    if len(valid_df) == 0:
        logger.error("No valid predictions found! Cannot create visualizations.")
        return 1
    
    logger.info("")
    logger.info("Generating visualizations...")
    logger.info("")
    
    # Generate all plots
    plot_score_distribution(valid_df, VISUALIZATIONS_DIR)
    plot_prediction_categories(valid_df, VISUALIZATIONS_DIR)
    plot_percentile_distribution(valid_df, VISUALIZATIONS_DIR)
    plot_score_vs_percentile(valid_df, VISUALIZATIONS_DIR)
    plot_summary_statistics(df, valid_df, VISUALIZATIONS_DIR)
    
    logger.info("")
    logger.info("="*70)
    logger.info("VISUALIZATION COMPLETE")
    logger.info("="*70)
    logger.info(f"Visualizations saved to: {VISUALIZATIONS_DIR}")
    logger.info("")
    logger.info("Generated plots:")
    logger.info("  - esm2_score_distribution.png")
    logger.info("  - esm2_prediction_categories.png")
    logger.info("  - esm2_percentile_distribution.png")
    logger.info("  - esm2_score_vs_percentile.png")
    logger.info("  - esm2_summary_statistics.png")
    logger.info("="*70)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

