#!/usr/bin/env python3
"""
Phase 5: Population-Level Genetic Load Analysis (Python)

Replaces the R implementation with pandas + matplotlib/seaborn.
Reads individual-level load outputs and produces population summaries and
visualizations focused on POPULATION-LEVEL patterns (no overlap with individual plots).

Outputs:
- CSV summaries under output/phase5_genetic_load/population_load
- Visualizations under output/phase5_genetic_load/population_load/visualizations
"""

import os
import sys
import math
import pandas as pd
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils import setup_logger


# =============================================================================
# SETUP
# =============================================================================

BASE_DIR = get_base_dir()
OUTPUT_DIR = f"{BASE_DIR}/output/phase5_genetic_load/population_load"
INPUT_FILE = f"{BASE_DIR}/output/phase5_genetic_load/individual_load/individual_genetic_load.csv"
VIS_DIR = f"{OUTPUT_DIR}/visualizations"

Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
Path(VIS_DIR).mkdir(parents=True, exist_ok=True)

logger = setup_logger("phase5_population_load", f"{OUTPUT_DIR}/population_load.log")


# =============================================================================
# PLOTTING STYLE
# =============================================================================

# Clean white background without grid background color
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3
plt.rcParams['grid.color'] = '#cccccc'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8

# Colorblind-friendly palette (Wong palette)
# Reference: https://www.nature.com/articles/nmeth.1618
CB_COLORS = {
    'blue': '#0072B2',        # Strong blue
    'orange': '#E69F00',      # Orange
    'sky_blue': '#56B4E9',    # Sky blue
    'green': '#009E73',       # Bluish green
    'yellow': '#F0E442',      # Yellow
    'vermillion': '#D55E00',  # Vermillion (red-ish)
    'purple': '#CC79A7',      # Reddish purple
    'black': '#000000',       # Black
}

# Set colorblind-friendly palette
sns.set_palette([CB_COLORS['blue'], CB_COLORS['orange'], CB_COLORS['green'], 
                 CB_COLORS['vermillion'], CB_COLORS['purple'], CB_COLORS['sky_blue']])


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_sig_stars(p_val):
    """Return significance stars based on p-value"""
    if p_val < 0.001:
        return '***'
    elif p_val < 0.01:
        return '**'
    elif p_val < 0.05:
        return '*'
    else:
        return 'ns'


# =============================================================================
# SUMMARY HELPERS
# =============================================================================

def summarize_population(load_df: pd.DataFrame):
    logger.info("Calculating population-level summaries...")

    # Variant count summary
    variant_stats = pd.DataFrame({
        'Metric': [
            'LOF Heterozygous',
            'LOF Homozygous',
            'Deleterious Missense Heterozygous',
            'Deleterious Missense Homozygous',
            'Total Deleterious Variants',
        ],
        'Mean': [
            load_df['LOF_Het'].mean(skipna=True),
            load_df['LOF_Hom'].mean(skipna=True),
            load_df['DelMis_Het'].mean(skipna=True),
            load_df['DelMis_Hom'].mean(skipna=True),
            load_df['Total_Deleterious'].mean(skipna=True),
        ],
        'SD': [
            load_df['LOF_Het'].std(skipna=True),
            load_df['LOF_Hom'].std(skipna=True),
            load_df['DelMis_Het'].std(skipna=True),
            load_df['DelMis_Hom'].std(skipna=True),
            load_df['Total_Deleterious'].std(skipna=True),
        ],
        'Min': [
            load_df['LOF_Het'].min(skipna=True),
            load_df['LOF_Hom'].min(skipna=True),
            load_df['DelMis_Het'].min(skipna=True),
            load_df['DelMis_Hom'].min(skipna=True),
            load_df['Total_Deleterious'].min(skipna=True),
        ],
        'Max': [
            load_df['LOF_Het'].max(skipna=True),
            load_df['LOF_Hom'].max(skipna=True),
            load_df['DelMis_Het'].max(skipna=True),
            load_df['DelMis_Hom'].max(skipna=True),
            load_df['Total_Deleterious'].max(skipna=True),
        ]
    })

    # Load metric summary
    # Check for new vs old column names
    if 'Total_Genetic_Load' in load_df.columns:
        # Latest columns (v4.1)
        cols = ['Total_Genetic_Load', 'Realized_Load', 'Potential_Load', 'Hom_Realized_Load', 'Het_Realized_Load']
        metric_names = [
            'Total Genetic Load',
            'Realized Load (Expressed)',
            'Potential Load (Hidden)',
            'Homozygous Realized Load',
            'Heterozygous Realized Load',
        ]
    elif 'Genetic_Load' in load_df.columns:
        # Old simplified columns
        cols = ['Genetic_Load', 'Realized_Load', 'Masked_Load']
        metric_names = [
            'Genetic Load (Total)',
            'Realized Load (Homozygous)',
            'Masked Load (Heterozygous)',
        ]
    else:
        # Very old columns (for backwards compatibility)
        cols = ['Genetic_Load_Multiplicative', 'Genetic_Load_Additive', 'Realized_Load', 'Masked_Load', 'Fitness']
        metric_names = [
            'Genetic Load (Multiplicative)',
            'Genetic Load (Additive)',
            'Realized Load (Homozygous)',
            'Masked Load (Heterozygous)',
            'Mean Fitness',
        ]
    
    # Filter to existing columns
    cols = [c for c in cols if c in load_df.columns]
    metric_names = metric_names[:len(cols)]
    
    load_stats = pd.DataFrame({
        'Metric': metric_names,
        'Mean': [load_df[c].mean(skipna=True) for c in cols],
        'SD': [load_df[c].std(skipna=True) for c in cols],
        'Min': [load_df[c].min(skipna=True) for c in cols],
        'Max': [load_df[c].max(skipna=True) for c in cols],
    })

    # ROH-associated summary (match updated categories)
    roh_stats = None
    required = ['Load_in_ROH', 'Load_outside_ROH', 'Load_in_Short_ROH', 'Load_in_Medium_ROH', 'Load_in_Long_ROH',
                'Variants_in_Short_ROH', 'Variants_in_Medium_ROH', 'Variants_in_Long_ROH']
    if all(c in load_df.columns for c in required):
        roh_stats = pd.DataFrame({
            'Metric': [
                'Load in ROH (All)',
                'Load outside ROH',
                'Load in Short ROH (<1 Mb, Ancient)',
                'Load in Medium ROH (1-5 Mb, Intermediate)',
                'Load in Long ROH (>5 Mb, Recent)',
                'Variants in Short ROH',
                'Variants in Medium ROH',
                'Variants in Long ROH',
            ],
            'Mean': [
                load_df['Load_in_ROH'].mean(skipna=True),
                load_df['Load_outside_ROH'].mean(skipna=True),
                load_df['Load_in_Short_ROH'].mean(skipna=True),
                load_df['Load_in_Medium_ROH'].mean(skipna=True),
                load_df['Load_in_Long_ROH'].mean(skipna=True),
                load_df['Variants_in_Short_ROH'].mean(skipna=True),
                load_df['Variants_in_Medium_ROH'].mean(skipna=True),
                load_df['Variants_in_Long_ROH'].mean(skipna=True),
            ],
            'SD': [
                load_df['Load_in_ROH'].std(skipna=True),
                load_df['Load_outside_ROH'].std(skipna=True),
                load_df['Load_in_Short_ROH'].std(skipna=True),
                load_df['Load_in_Medium_ROH'].std(skipna=True),
                load_df['Load_in_Long_ROH'].std(skipna=True),
                load_df['Variants_in_Short_ROH'].std(skipna=True),
                load_df['Variants_in_Medium_ROH'].std(skipna=True),
                load_df['Variants_in_Long_ROH'].std(skipna=True),
            ],
            'Min': [
                load_df['Load_in_ROH'].min(skipna=True),
                load_df['Load_outside_ROH'].min(skipna=True),
                load_df['Load_in_Short_ROH'].min(skipna=True),
                load_df['Load_in_Medium_ROH'].min(skipna=True),
                load_df['Load_in_Long_ROH'].min(skipna=True),
                load_df['Variants_in_Short_ROH'].min(skipna=True),
                load_df['Variants_in_Medium_ROH'].min(skipna=True),
                load_df['Variants_in_Long_ROH'].min(skipna=True),
            ],
            'Max': [
                load_df['Load_in_ROH'].max(skipna=True),
                load_df['Load_outside_ROH'].max(skipna=True),
                load_df['Load_in_Short_ROH'].max(skipna=True),
                load_df['Load_in_Medium_ROH'].max(skipna=True),
                load_df['Load_in_Long_ROH'].max(skipna=True),
                load_df['Variants_in_Short_ROH'].max(skipna=True),
                load_df['Variants_in_Medium_ROH'].max(skipna=True),
                load_df['Variants_in_Long_ROH'].max(skipna=True),
            ],
        })
        logger.info("ROH-associated load metrics calculated")
    else:
        logger.warning("ROH-associated columns not found; skipping ROH summary")

    # Save CSVs
    variant_stats.to_csv(f"{OUTPUT_DIR}/variant_count_summary.csv", index=False)
    load_stats.to_csv(f"{OUTPUT_DIR}/genetic_load_summary.csv", index=False)
    if roh_stats is not None:
        roh_stats.to_csv(f"{OUTPUT_DIR}/roh_load_summary.csv", index=False)
    
    # Calculate key correlations
    logger.info("Calculating key correlations...")
    correlation_results = []
    
    try:
        from scipy.stats import spearmanr, pearsonr
        
        # Key correlation pairs to test
        test_pairs = [
            ('Total_Het', 'Total_Hom', 'Heterozygous vs Homozygous burden'),
            ('F_ROH', 'Total_Genetic_Load', 'F_ROH vs Total Load'),
            ('F_ROH', 'Realized_Load', 'F_ROH vs Realized Load'),
            ('F_ROH', 'Potential_Load', 'F_ROH vs Potential Load'),
            ('F_ROH', 'Load_in_ROH', 'F_ROH vs ROH Load'),
            ('Total_Genetic_Load', 'Total_Deleterious', 'Load vs Variant Count'),
        ]
        
        for var1, var2, description in test_pairs:
            if var1 in load_df.columns and var2 in load_df.columns:
                # Remove NaN values
                valid_data = load_df[[var1, var2]].dropna()
                
                if len(valid_data) > 2:
                    # Spearman correlation (non-parametric)
                    rho, p_spearman = spearmanr(valid_data[var1], valid_data[var2])
                    
                    # Pearson correlation
                    r, p_pearson = pearsonr(valid_data[var1], valid_data[var2])
                    
                    # Significance indicators
                    spearman_sig = get_sig_stars(p_spearman)
                    pearson_sig = get_sig_stars(p_pearson)
                    
                    correlation_results.append({
                        'Variable_1': var1,
                        'Variable_2': var2,
                        'Description': description,
                        'N': len(valid_data),
                        'Spearman_rho': rho,
                        'Spearman_p': p_spearman,
                        'Spearman_sig': spearman_sig,
                        'Pearson_r': r,
                        'Pearson_p': p_pearson,
                        'Pearson_sig': pearson_sig
                    })
                    
                    logger.info(f"  {description}: ρ={rho:.3f}, p={p_spearman:.2e} {spearman_sig}")
        
        # Save correlation results
        if correlation_results:
            corr_df = pd.DataFrame(correlation_results)
            corr_df.to_csv(f"{OUTPUT_DIR}/correlation_summary.csv", index=False)
            logger.info(f"Saved correlation results to correlation_summary.csv")
            
            # Print significant correlations
            sig_corr = corr_df[corr_df['Spearman_p'] < 0.05]
            if len(sig_corr) > 0:
                logger.info(f"\n{len(sig_corr)} significant correlations (p < 0.05):")
                for _, row in sig_corr.iterrows():
                    logger.info(f"  {row['Description']}: ρ={row['Spearman_rho']:.3f}, "
                              f"p={row['Spearman_p']:.4f} {row['Spearman_sig']}")
            else:
                logger.info("No significant correlations found (p < 0.05)")
    
    except Exception as e:
        logger.warning(f"Could not calculate correlations: {e}")
        correlation_results = None

    return variant_stats, load_stats, roh_stats


# =============================================================================
# VISUALIZATIONS (POPULATION-LEVEL ONLY, NO OVERLAP WITH INDIVIDUAL PLOTS)
# =============================================================================

def savefig(path: str):
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"  Saved: {path}")


def generate_plots(load_df: pd.DataFrame):
    logger.info("Generating population-level visualizations...")
    
    # Colorblind-friendly colors for bar plots
    cb_bar_colors = [CB_COLORS['blue'], CB_COLORS['orange'], CB_COLORS['green'], CB_COLORS['purple']]
    
    # 1) Population summary bar plot (means + error bars)
    logger.info("  Creating population summary barplot...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Variant counts
    ax = axes[0]
    categories = ['LOF\nHet', 'LOF\nHom', 'DelMis\nHet', 'DelMis\nHom']
    means = [load_df['LOF_Het'].mean(), load_df['LOF_Hom'].mean(), 
             load_df['DelMis_Het'].mean(), load_df['DelMis_Hom'].mean()]
    stds = [load_df['LOF_Het'].std(), load_df['LOF_Hom'].std(),
            load_df['DelMis_Het'].std(), load_df['DelMis_Hom'].std()]
    bars = ax.bar(categories, means, yerr=stds, capsize=5, color=cb_bar_colors, 
                  alpha=0.8, edgecolor='white', linewidth=1.5)
    for bar, mean, std in zip(bars, means, stds):
        # Place text above error bar
        y_pos = bar.get_height() + std + max(means) * 0.02
        ax.text(bar.get_x() + bar.get_width()/2., y_pos,
                f'{mean:.1f}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    ax.set_ylabel('Mean Count per Individual', fontsize=11, fontweight='bold')
    ax.set_title('Population Mean Deleterious Variant Burden', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Load metrics - handle different column versions
    ax = axes[1]
    if 'Total_Genetic_Load' in load_df.columns:
        # Latest columns (v4.1)
        categories = ['Total\nLoad', 'Realized\nLoad', 'Potential\nLoad']
        means = [load_df['Total_Genetic_Load'].mean(), load_df['Realized_Load'].mean(),
                 load_df['Potential_Load'].mean()]
        stds = [load_df['Total_Genetic_Load'].std(), load_df['Realized_Load'].std(),
                load_df['Potential_Load'].std()]
        colors_bar = [CB_COLORS['blue'], CB_COLORS['vermillion'], CB_COLORS['sky_blue']]
    elif 'Genetic_Load' in load_df.columns:
        # Old simplified columns
        categories = ['Total\nLoad', 'Realized\nLoad', 'Masked\nLoad']
        means = [load_df['Genetic_Load'].mean(), load_df['Realized_Load'].mean(),
                 load_df['Masked_Load'].mean()]
        stds = [load_df['Genetic_Load'].std(), load_df['Realized_Load'].std(),
                load_df['Masked_Load'].std()]
        colors_bar = [CB_COLORS['blue'], CB_COLORS['vermillion'], CB_COLORS['sky_blue']]
    else:
        # Very old columns
        categories = ['Multiplicative\nLoad', 'Realized\nLoad', 'Masked\nLoad', 'Mean\nFitness']
        means = [load_df['Genetic_Load_Multiplicative'].mean(), load_df['Realized_Load'].mean(),
                 load_df['Masked_Load'].mean(), load_df['Fitness'].mean()]
        stds = [load_df['Genetic_Load_Multiplicative'].std(), load_df['Realized_Load'].std(),
                load_df['Masked_Load'].std(), load_df['Fitness'].std()]
        colors_bar = cb_bar_colors
    
    bars = ax.bar(categories, means, yerr=stds, capsize=5, color=colors_bar,
                  alpha=0.8, edgecolor='white', linewidth=1.5)
    for bar, mean, std in zip(bars, means, stds):
        # Place text above error bar
        y_pos = bar.get_height() + std + max(means) * 0.02
        ax.text(bar.get_x() + bar.get_width()/2., y_pos,
                f'{mean:.3f}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    ax.set_ylabel('Mean Value', fontsize=11, fontweight='bold')
    ax.set_title('Population Load Metrics Summary', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    savefig(f"{VIS_DIR}/population_summary_barplot.png")
    
    # 2) Individual load ranking (top & bottom)
    logger.info("  Creating individual load ranking...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Determine which load column to use
    if 'Total_Genetic_Load' in load_df.columns:
        load_col = 'Total_Genetic_Load'
    elif 'Genetic_Load' in load_df.columns:
        load_col = 'Genetic_Load'
    else:
        load_col = 'Genetic_Load_Multiplicative'
    
    # Top 20
    ax = axes[0]
    top_df = load_df.sort_values(load_col, ascending=False).head(20)
    y_pos = np.arange(len(top_df))
    ax.barh(y_pos, top_df[load_col], color=CB_COLORS['vermillion'], alpha=0.8, edgecolor='white')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top_df['IID'], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel('Genetic Load', fontsize=11, fontweight='bold')
    ax.set_title('Top 20 Highest Load Individuals', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    
    # Bottom 20
    ax = axes[1]
    bottom_df = load_df.sort_values(load_col, ascending=True).head(20)
    y_pos = np.arange(len(bottom_df))
    ax.barh(y_pos, bottom_df[load_col], color=CB_COLORS['blue'], alpha=0.8, edgecolor='white')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(bottom_df['IID'], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel('Genetic Load', fontsize=11, fontweight='bold')
    ax.set_title('Bottom 20 Lowest Load Individuals', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    savefig(f"{VIS_DIR}/individual_load_ranking.png")
    
    # 3) Het vs Hom burden across population
    logger.info("  Creating het vs hom burden comparison...")
    fig, ax = plt.subplots(figsize=(10, 7))
    x = load_df['Total_Het']
    y = load_df['Total_Hom']
    # Use colorblind-friendly colormap (viridis)
    sc = ax.scatter(x, y, alpha=0.6, s=50, c=load_df[load_col],
                    cmap='viridis', edgecolors='white', linewidth=0.5)
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('Genetic Load', fontsize=11, fontweight='bold')
    try:
        from scipy import stats
        slope, intercept, r, p, _ = stats.linregress(x, y)
        xs = np.linspace(x.min(), x.max(), 100)
        ax.plot(xs, slope * xs + intercept, color=CB_COLORS['vermillion'], linestyle='--', linewidth=2)
        
        # Add significance stars
        sig_stars = get_sig_stars(p)
        ax.text(0.02, 0.95, f"r = {r:.3f}\np = {p:.2e} {sig_stars}", transform=ax.transAxes,
                va='top', ha='left', fontsize=10, fontweight='bold',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray', linewidth=1.5))
        
        # Log the correlation result
        logger.info(f"    Het vs Hom correlation: r={r:.3f}, p={p:.2e} {sig_stars}")
    except Exception as e:
        logger.warning(f"    Could not calculate Het vs Hom correlation: {e}")
        pass
    ax.set_xlabel('Total Heterozygous Deleterious', fontsize=11, fontweight='bold')
    ax.set_ylabel('Total Homozygous Deleterious', fontsize=11, fontweight='bold')
    ax.set_title('Heterozygous vs Homozygous Burden Across Population', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    savefig(f"{VIS_DIR}/het_vs_hom_burden.png")
    
    # 4) Load components stacked (realized + potential)
    logger.info("  Creating load composition stacked view...")
    fig, ax = plt.subplots(figsize=(14, 6))
    sorted_df = load_df.sort_values(load_col, ascending=False)
    x_pos = np.arange(len(sorted_df))
    
    # Use Potential_Load if available, otherwise fall back to Masked_Load
    potential_col = 'Potential_Load' if 'Potential_Load' in load_df.columns else 'Masked_Load'
    potential_label = 'Potential Load (Hidden)' if 'Potential_Load' in load_df.columns else 'Masked Load (Het)'
    
    ax.bar(x_pos, sorted_df['Realized_Load'], label='Realized Load (Expressed)', 
           color=CB_COLORS['vermillion'], alpha=0.8, edgecolor='white', linewidth=0.3)
    ax.bar(x_pos, sorted_df[potential_col], bottom=sorted_df['Realized_Load'],
           label=potential_label, color=CB_COLORS['sky_blue'], alpha=0.8, edgecolor='white', linewidth=0.3)
    ax.set_xlabel('Individuals (sorted by total load)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Genetic Load', fontsize=11, fontweight='bold')
    ax.set_title('Population Load Composition: Realized vs Potential', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.set_xticks([])
    savefig(f"{VIS_DIR}/population_load_composition.png")
    
    # Colorblind-friendly colors for ROH categories (sequential: blue → orange → vermillion)
    roh_colors = [CB_COLORS['sky_blue'], CB_COLORS['orange'], CB_COLORS['vermillion']]
    
    # 5) ROH temporal burden (violin plot + bar chart) with statistical tests
    if all(c in load_df.columns for c in ['Load_in_Short_ROH', 'Load_in_Medium_ROH', 'Load_in_Long_ROH']):
        logger.info("  Creating ROH temporal burden analysis...")
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # Violin plot
        ax = axes[0]
        roh_data = [load_df['Load_in_Short_ROH'].dropna(),
                    load_df['Load_in_Medium_ROH'].dropna(),
                    load_df['Load_in_Long_ROH'].dropna()]
        parts = ax.violinplot(roh_data, positions=[0, 1, 2], showmeans=True, showmedians=True)
        for pc, color in zip(parts['bodies'], roh_colors):
            pc.set_facecolor(color)
            pc.set_alpha(0.7)
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(['Short\n(<1 Mb)\nAncient', 'Medium\n(1-5 Mb)\nIntermediate', 'Long\n(>5 Mb)\nRecent'], fontsize=10)
        ax.set_ylabel('Genetic Load', fontsize=11, fontweight='bold')
        
        # Add Kruskal-Wallis test (non-parametric ANOVA)
        try:
            from scipy.stats import kruskal
            stat, p_kw = kruskal(roh_data[0], roh_data[1], roh_data[2])
            sig_stars = get_sig_stars(p_kw)
            ax.set_title(f'Load Distribution by ROH Age Category\n(Kruskal-Wallis p={p_kw:.2e} {sig_stars})', 
                        fontsize=12, fontweight='bold')
            logger.info(f"    Kruskal-Wallis test (ROH categories): p={p_kw:.2e} {sig_stars}")
        except Exception as e:
            ax.set_title('Load Distribution by ROH Age Category', fontsize=13, fontweight='bold')
            logger.warning(f"    Could not perform Kruskal-Wallis test: {e}")
        
        ax.grid(True, alpha=0.3, axis='y')
        
        # Perform pairwise comparisons (Mann-Whitney U tests)
        try:
            from scipy.stats import mannwhitneyu
            comparisons = [
                ('Short', 'Medium', 0, 1),
                ('Short', 'Long', 0, 2),
                ('Medium', 'Long', 1, 2)
            ]
            
            logger.info("    Pairwise comparisons (Mann-Whitney U):")
            for name1, name2, idx1, idx2 in comparisons:
                _, p_val = mannwhitneyu(roh_data[idx1], roh_data[idx2], alternative='two-sided')
                sig_stars = get_sig_stars(p_val)
                logger.info(f"      {name1} vs {name2}: p={p_val:.4f} {sig_stars}")
        except Exception as e:
            logger.warning(f"    Could not perform pairwise comparisons: {e}")
        
        # Variant counts
        ax = axes[1]
        var_means = [load_df['Variants_in_Short_ROH'].mean(), 
                     load_df['Variants_in_Medium_ROH'].mean(),
                     load_df['Variants_in_Long_ROH'].mean()]
        var_stds = [load_df['Variants_in_Short_ROH'].std(),
                    load_df['Variants_in_Medium_ROH'].std(),
                    load_df['Variants_in_Long_ROH'].std()]
        labels = ['Short\n(<1 Mb)', 'Medium\n(1-5 Mb)', 'Long\n(>5 Mb)']
        bars = ax.bar(labels, var_means, yerr=var_stds, capsize=5, color=roh_colors,
                     alpha=0.8, edgecolor='white', linewidth=1.5)
        for bar, mean in zip(bars, var_means):
            ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + bar.get_height()*0.05,
                    f'{mean:.1f}', ha='center', va='bottom', fontsize=10, fontweight='bold')
        ax.set_ylabel('Mean Variant Count', fontsize=11, fontweight='bold')
        ax.set_title('Deleterious Variants by Inbreeding Era', fontsize=13, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        savefig(f"{VIS_DIR}/roh_temporal_burden.png")
    
    # 6) ROH load contribution (pie charts)
    if all(c in load_df.columns for c in ['Load_in_ROH', 'Load_outside_ROH', 'Load_in_Short_ROH', 'Load_in_Medium_ROH', 'Load_in_Long_ROH']):
        logger.info("  Creating ROH load contribution analysis...")
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # Overall ROH vs non-ROH (colorblind-friendly: vermillion vs blue)
        ax = axes[0]
        roh_contrib = load_df['Load_in_ROH'].mean()
        non_roh_contrib = load_df['Load_outside_ROH'].mean()
        sizes = [roh_contrib, non_roh_contrib]
        labels = [f'In ROH\n{roh_contrib:.3f}', f'Outside ROH\n{non_roh_contrib:.3f}']
        colors_pie = [CB_COLORS['vermillion'], CB_COLORS['blue']]
        wedges, texts, autotexts = ax.pie(sizes, labels=labels, colors=colors_pie, autopct='%1.1f%%',
                                           startangle=90, textprops={'fontsize': 11, 'fontweight': 'bold'})
        ax.set_title('Population Mean Load: ROH vs Non-ROH', fontsize=13, fontweight='bold')
        
        # ROH age categories (colorblind-friendly sequential)
        ax = axes[1]
        short_mean = load_df['Load_in_Short_ROH'].mean()
        medium_mean = load_df['Load_in_Medium_ROH'].mean()
        long_mean = load_df['Load_in_Long_ROH'].mean()
        sizes = [short_mean, medium_mean, long_mean]
        labels = [f'Short\n{short_mean:.3f}', f'Medium\n{medium_mean:.3f}', f'Long\n{long_mean:.3f}']
        wedges, texts, autotexts = ax.pie(sizes, labels=labels, colors=roh_colors, autopct='%1.1f%%',
                                           startangle=90, textprops={'fontsize': 11, 'fontweight': 'bold'})
        ax.set_title('ROH Load by Inbreeding Era', fontsize=13, fontweight='bold')
        savefig(f"{VIS_DIR}/roh_load_contribution.png")
    
    # 7) FROH stratification by load quartiles
    if 'F_ROH' in load_df.columns:
        logger.info("  Creating FROH stratification analysis...")
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # Colorblind-friendly quartile colors (sequential: blue → sky blue → orange → vermillion)
        quartile_colors_cb = {
            'Q1 (Lowest)': CB_COLORS['blue'], 
            'Q2': CB_COLORS['sky_blue'], 
            'Q3': CB_COLORS['orange'], 
            'Q4 (Highest)': CB_COLORS['vermillion']
        }
        
        # FROH histogram by load quartile
        ax = axes[0]
        load_df['Load_Quartile'] = pd.qcut(load_df[load_col], q=4, 
                                            labels=['Q1 (Lowest)', 'Q2', 'Q3', 'Q4 (Highest)'])
        for quartile in ['Q1 (Lowest)', 'Q2', 'Q3', 'Q4 (Highest)']:
            sub = load_df[load_df['Load_Quartile'] == quartile]['F_ROH']
            ax.hist(sub, bins=15, alpha=0.5, label=quartile, edgecolor='white',
                   color=quartile_colors_cb[quartile])
        ax.set_xlabel(r'F$_{ROH}$', fontsize=11, fontweight='bold')
        ax.set_ylabel('Count', fontsize=11, fontweight='bold')
        ax.set_title(r'F$_{ROH}$ Distribution by Load Quartile', fontsize=13, fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        
        # FROH vs Load colored by quartile
        ax = axes[1]
        for quartile in ['Q1 (Lowest)', 'Q2', 'Q3', 'Q4 (Highest)']:
            sub = load_df[load_df['Load_Quartile'] == quartile]
            ax.scatter(sub['F_ROH'], sub[load_col], 
                      alpha=0.6, s=40, label=quartile, color=quartile_colors_cb[quartile], 
                      edgecolors='white', linewidth=0.5)
        try:
            from scipy import stats
            sub_all = load_df[['F_ROH', load_col]].dropna()
            slope, intercept, r, p, _ = stats.linregress(sub_all['F_ROH'], sub_all[load_col])
            xs = np.linspace(sub_all['F_ROH'].min(), sub_all['F_ROH'].max(), 100)
            
            # Add significance stars
            sig_stars = get_sig_stars(p)
            ax.plot(xs, slope * xs + intercept, 'k--', linewidth=2, 
                   label=f'r={r:.3f}, p={p:.2e} {sig_stars}')
            
            # Log the correlation result
            logger.info(f"    F_ROH vs Load correlation: r={r:.3f}, p={p:.2e} {sig_stars}")
        except Exception as e:
            logger.warning(f"    Could not calculate F_ROH vs Load correlation: {e}")
            pass
        ax.set_xlabel(r'F$_{ROH}$', fontsize=11, fontweight='bold')
        ax.set_ylabel('Genetic Load', fontsize=11, fontweight='bold')
        ax.set_title(r'Inbreeding (F$_{ROH}$) vs Load by Quartile', fontsize=13, fontweight='bold')
        ax.legend(fontsize=8, loc='best')
        ax.grid(True, alpha=0.3)
        savefig(f"{VIS_DIR}/froh_load_stratification.png")


def main():
    logger.info("========================================")
    logger.info("POPULATION-LEVEL GENETIC LOAD ANALYSIS (Python)")
    logger.info("========================================")

    if not os.path.exists(INPUT_FILE):
        logger.error(f"Input file not found: {INPUT_FILE}")
        return 1

    logger.info("Loading individual genetic load data...")
    load_df = pd.read_csv(INPUT_FILE)
    logger.info(f"Loaded data for {len(load_df)} individuals")

    # Summaries
    summarize_population(load_df)

    # Visualizations
    generate_plots(load_df)

    logger.info("\nAnalysis complete.")
    logger.info(f"Output directory: {OUTPUT_DIR}")
    logger.info(f"Visualizations: {VIS_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
