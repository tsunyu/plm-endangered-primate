#!/usr/bin/env python3
"""
Visualize Demographic History from fastsimcoal2 Results
========================================================

This script creates publication-quality visualizations of the inferred
demographic history.

Input:  fastsimcoal2 parameter estimates
Output: Demographic history plots, SFS fit plots

Author: Demographic Analysis Pipeline
Date: 2026-01-26
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FuncFormatter
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir
import logging

# Configuration
BASE_DIR = get_base_dir()
OUTPUT_DIR = BASE_DIR / "output/phase3b_fastsimcoal2"
PLOT_DIR = OUTPUT_DIR / "plots"
PLOT_DIR.mkdir(parents=True, exist_ok=True)

GENERATION_TIME = 10  # years
MUTATION_RATE = 1.36e-8

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

MODEL_LABELS = {
    'constant_ne': 'Constant Ne',
    'single_bottleneck': 'Single Bottleneck',
    'two_consecutive_bottlenecks': 'Two Consecutive Bottlenecks',
    'bottleneck_continuous_decline': 'Bottleneck + Continuous Decline',
    'bottleneck_recent_contraction': 'Bottleneck + Recent Contraction',
    'complex_multi_event': 'Complex Multi-event'
}


def configure_plot_style():
    """Configure a cleaner plotting style for publication-ready figures."""
    plt.rcParams.update({
        'figure.facecolor': 'white',
        'axes.facecolor': '#f8fafc',
        'axes.edgecolor': '#334155',
        'axes.grid': True,
        'grid.color': '#cbd5e1',
        'grid.alpha': 0.35,
        'grid.linestyle': '--',
        'font.size': 11,
        'axes.titlesize': 15,
        'axes.labelsize': 12,
        'legend.frameon': True,
        'legend.framealpha': 0.9
    })


def load_fastsimcoal2_parameters():
    """Load best-fit parameters from fastsimcoal2."""
    param_file = OUTPUT_DIR / "parameter_estimates.csv"
    
    if not param_file.exists():
        logger.error(f"Parameter file not found: {param_file}")
        logger.error("Please run phase3b_step5_analyze_results.py first")
        return None, None
    
    df = pd.read_csv(param_file)
    
    # Convert to dictionary
    params = {}
    for _, row in df.iterrows():
        params[row['Parameter']] = row['Estimate']
    
    # Determine model type from parameters
    if {'T1', 'T2', 'T3', 'T4'}.issubset(params):
        model_type = 'complex_multi_event'
    elif {'T_REC', 'T_BOT', 'TRECENT'}.issubset(params) or {'TRECOVERY_OLD', 'TBOT_OLD', 'TRECENT'}.issubset(params):
        model_type = 'bottleneck_recent_contraction'
    elif {'TRECENT_BOT', 'TOLD_BOT'}.issubset(params):
        model_type = 'two_consecutive_bottlenecks'
    elif {'TRECENT', 'TBOT', 'TANC'}.issubset(params):
        model_type = 'bottleneck_continuous_decline'
    elif {'TBOT', 'TRECOVERY'}.issubset(params):
        model_type = 'single_bottleneck'
    else:
        model_type = 'constant_ne'
    
    logger.info(f"Model type: {model_type}")
    logger.info(f"Parameters: {params}")
    
    return params, model_type


def generate_demographic_trajectory(params, model_type, max_time_gen=20000):
    """
    Generate Ne trajectory over time from fastsimcoal2 parameters.
    
    Returns:
        time_years: Array of time points (years ago)
        ne_values: Array of Ne values
    """
    time_gen = np.linspace(0, max_time_gen, 1000)
    ne_values = np.zeros_like(time_gen)
    
    if model_type == 'constant_ne':
        ne_values[:] = params['NCUR']
    
    elif model_type == 'single_bottleneck':
        ncur = params['NCUR']
        nbot = params['NBOT']
        nanc = params['NANC']
        tbot = params['TBOT']
        trecover = params['TRECOVERY']
        
        ne_values[time_gen <= trecover] = ncur
        mask_bottleneck = (time_gen > trecover) & (time_gen <= tbot)
        ne_values[mask_bottleneck] = nbot
        ne_values[time_gen > tbot] = nanc

    elif model_type == 'two_consecutive_bottlenecks':
        ncur = params['NCUR']
        nbot2 = params['NBOT2']
        ninter = params['NINTER']
        nbot1 = params['NBOT1']
        nanc = params['NANC']
        trecent_recovery = params['TRECENT_RECOVERY']
        trecent_bot = params['TRECENT_BOT']
        told_recovery = params['TOLD_RECOVERY']
        told_bot = params['TOLD_BOT']

        ne_values[time_gen <= trecent_recovery] = ncur
        ne_values[(time_gen > trecent_recovery) & (time_gen <= trecent_bot)] = nbot2
        ne_values[(time_gen > trecent_bot) & (time_gen <= told_recovery)] = ninter
        ne_values[(time_gen > told_recovery) & (time_gen <= told_bot)] = nbot1
        ne_values[time_gen > told_bot] = nanc

    elif model_type == 'bottleneck_continuous_decline':
        ncur = params['NCUR']
        nmid = params['NMID']
        nbot = params['NBOT']
        nanc = params['NANC']
        trecent = params['TRECENT']
        tbot = params['TBOT']
        tanc = params['TANC']

        ne_values[time_gen <= trecent] = ncur
        ne_values[(time_gen > trecent) & (time_gen <= tbot)] = nmid
        ne_values[(time_gen > tbot) & (time_gen <= tanc)] = nbot
        ne_values[time_gen > tanc] = nanc

    elif model_type == 'bottleneck_recent_contraction':
        ncur = params['NCUR']
        nrecover = params['NRECOVER']
        nbot = params['NBOT']
        nanc = params['NANC']
        trecent = params['TRECENT']
        trecovery_old = params.get('T_REC', params.get('TRECOVERY_OLD', 0))
        tbot_old = params.get('T_BOT', params.get('TBOT_OLD', 0))

        ne_values[time_gen <= trecent] = ncur
        ne_values[(time_gen > trecent) & (time_gen <= trecovery_old)] = nrecover
        ne_values[(time_gen > trecovery_old) & (time_gen <= tbot_old)] = nbot
        ne_values[time_gen > tbot_old] = nanc

    elif model_type == 'complex_multi_event':
        ncur = params['NCUR']
        n1 = params['N1']
        n2 = params['N2']
        n3 = params['N3']
        nanc = params['NANC']
        t1 = params['T1']
        t2 = params['T2']
        t3 = params['T3']
        t4 = params['T4']

        ne_values[time_gen <= t1] = ncur
        ne_values[(time_gen > t1) & (time_gen <= t2)] = n1
        ne_values[(time_gen > t2) & (time_gen <= t3)] = n2
        ne_values[(time_gen > t3) & (time_gen <= t4)] = n3
        ne_values[time_gen > t4] = nanc
    
    time_years = time_gen * GENERATION_TIME
    
    return time_years, ne_values




def plot_demographic_history(params, model_type):
    """Create main demographic history plot."""
    logger.info("Creating demographic history plot")
    configure_plot_style()
    
    # Generate fastsimcoal2 trajectory
    time_fsc, ne_fsc = generate_demographic_trajectory(params, model_type)
    
    # Create plot
    fig, ax = plt.subplots(figsize=(13, 7.5))
    time_kya = time_fsc / 1000.0
    
    # Use step plot to reflect epoch-wise demographic changes.
    ax.step(time_kya, ne_fsc, where='post', linewidth=2.8, color='#0f766e',
            label='Inferred Ne trajectory', zorder=3)
    ax.fill_between(time_kya, ne_fsc, np.min(ne_fsc[ne_fsc > 0]) * 0.9,
                    color='#14b8a6', alpha=0.12, step='post', zorder=2)
    
    # Mark key events
    event_lines = []
    if model_type == 'single_bottleneck':
        tbot_years = params['TBOT'] * GENERATION_TIME
        trecover_years = params['TRECOVERY'] * GENERATION_TIME
        event_lines = [
            (tbot_years / 1000, '#dc2626', f'Bottleneck onset ({tbot_years:.0f} ya)'),
            (trecover_years / 1000, '#2563eb', f'Recovery start ({trecover_years:.0f} ya)')
        ]
    
    elif model_type == 'bottleneck_recent_contraction':
        trecent_years = params['TRECENT'] * GENERATION_TIME
        event_lines = [(trecent_years / 1000, '#dc2626', f'Recent contraction ({trecent_years:.0f} ya)')]

    elif model_type == 'two_consecutive_bottlenecks':
        event_lines = [
            (params['TRECENT_RECOVERY'] * GENERATION_TIME / 1000, '#2563eb', 'Recent recovery'),
            (params['TRECENT_BOT'] * GENERATION_TIME / 1000, '#dc2626', 'Recent bottleneck'),
            (params['TOLD_RECOVERY'] * GENERATION_TIME / 1000, '#7c3aed', 'Old recovery'),
            (params['TOLD_BOT'] * GENERATION_TIME / 1000, '#b91c1c', 'Old bottleneck')
        ]
    elif model_type == 'bottleneck_continuous_decline':
        event_lines = [
            (params['TRECENT'] * GENERATION_TIME / 1000, '#dc2626', 'Recent decline'),
            (params['TBOT'] * GENERATION_TIME / 1000, '#7c3aed', 'Older shift'),
            (params['TANC'] * GENERATION_TIME / 1000, '#2563eb', 'Ancestral boundary')
        ]
    elif model_type == 'complex_multi_event':
        event_lines = [
            (params['T1'] * GENERATION_TIME / 1000, '#2563eb', 'T1'),
            (params['T2'] * GENERATION_TIME / 1000, '#7c3aed', 'T2'),
            (params['T3'] * GENERATION_TIME / 1000, '#b45309', 'T3'),
            (params['T4'] * GENERATION_TIME / 1000, '#dc2626', 'T4')
        ]

    for x_event, color, label in event_lines:
        ax.axvline(x_event, color=color, linestyle=':', linewidth=1.8, alpha=0.75, label=label, zorder=4)

    # Explicitly mark current Ne to improve readability.
    ncur = params.get('NCUR')
    if ncur is not None and ncur > 0:
        ax.axhline(ncur, color='#0f766e', linestyle='--', linewidth=1.2, alpha=0.6, zorder=1)
        x_present = max(0.05, np.min(time_kya[time_kya > 0]))
        ax.scatter([x_present], [ncur], color='#0f766e', s=28, zorder=5)
        ax.annotate(
            f'Current diploid Ne = {ncur:.0f}',
            xy=(x_present, ncur),
            xytext=(x_present * 1.8, ncur * 1.35),
            textcoords='data',
            fontsize=10,
            fontweight='bold',
            color='#134e4a',
            arrowprops=dict(arrowstyle='->', color='#134e4a', lw=1.1, alpha=0.9),
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='#99f6e4', alpha=0.95),
            zorder=6
        )
    
    # Formatting
    ax.set_xlabel('Time before present (kya)', fontweight='bold', fontsize=13)
    ax.set_ylabel('Diploid effective population size (Ne)', fontweight='bold', fontsize=13)
    ax.set_title('Demographic History of Shennongjia Golden Snub-Nosed Monkeys\n' +
                f'fastsimcoal2 ({MODEL_LABELS.get(model_type, model_type)}, generation time = {GENERATION_TIME} years)',
                fontweight='bold', fontsize=15)
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(max(0.05, np.min(time_kya[time_kya > 0])), time_kya[-1])
    ax.set_ylim(np.min(ne_fsc[ne_fsc > 0]) * 0.9, np.max(ne_fsc) * 1.2)
    ax.grid(True, which='major', alpha=0.35)
    ax.grid(True, which='minor', alpha=0.15)
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:g}'))
    ax.legend(loc='upper left', fontsize=9, ncol=2)
    
    # Add information box
    textstr = f'Model: {MODEL_LABELS.get(model_type, model_type)}\n'
    
    if 'NCUR' in params:
        textstr += f'Current diploid Ne: {params["NCUR"]:.0f}\n'
    if 'NBOT' in params:
        textstr += f'Bottleneck diploid Ne: {params["NBOT"]:.0f}\n'
    if 'TBOT' in params:
        tbot_years = params['TBOT'] * GENERATION_TIME
        textstr += f'Bottleneck: {tbot_years:.0f} ya\n'
    
    props = dict(boxstyle='round,pad=0.35', facecolor='white', edgecolor='#94a3b8', alpha=0.9)
    ax.text(0.015, 0.98, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    
    # Save
    plot_file = PLOT_DIR / "demographic_history.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    logger.info(f"Saved: {plot_file}")
    
    plot_file_pdf = PLOT_DIR / "demographic_history.pdf"
    plt.savefig(plot_file_pdf, bbox_inches='tight')
    logger.info(f"Saved: {plot_file_pdf}")
    
    plt.close()


def plot_sfs_fit():
    """Plot observed vs expected SFS under best-fit model."""
    logger.info("Creating SFS fit plot")
    configure_plot_style()
    
    # Load observed SFS
    sfs_file = OUTPUT_DIR / "sfs" / "SNJ_DAFpop0.obs"
    
    if not sfs_file.exists():
        logger.warning(f"SFS file not found: {sfs_file}")
        return
    
    # Read SFS
    with open(sfs_file, 'r') as f:
        lines = f.readlines()
        if len(lines) >= 3:
            data_line = lines[2].strip().split()
            obs_sfs = np.array([float(x) for x in data_line])
        else:
            logger.error("Invalid SFS file format")
            return
    
    # For expected SFS, we would need to run fastsimcoal2 simulation
    # For now, just plot observed SFS
    
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    
    # Plot 1: Full SFS
    ax = axes[0]
    x = np.arange(len(obs_sfs))
    ax.bar(x, obs_sfs, color='#2563eb', alpha=0.78, edgecolor='white', linewidth=0.3, label='Observed')
    ax.set_xlabel('Derived Allele Count', fontweight='bold')
    ax.set_ylabel('Number of Sites', fontweight='bold')
    ax.set_title('Observed Site Frequency Spectrum', fontweight='bold')
    ax.set_yscale('log')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    # Plot 2: Low frequency spectrum (excluding singletons)
    ax = axes[1]
    max_freq = min(30, len(obs_sfs))
    x_low = np.arange(2, max_freq)
    sfs_low = obs_sfs[2:max_freq]
    ax.bar(x_low, sfs_low, color='#fb7185', alpha=0.82, edgecolor='white', linewidth=0.3)
    ax.set_xlabel('Derived Allele Count', fontweight='bold')
    ax.set_ylabel('Number of Sites', fontweight='bold')
    ax.set_title('Low Frequency Spectrum (DAF 2-30)', fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    plot_file = PLOT_DIR / "sfs_observed.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    logger.info(f"Saved: {plot_file}")
    
    plt.close()


def main():
    """Main workflow."""
    logger.info("=" * 80)
    logger.info("DEMOGRAPHIC HISTORY VISUALIZATION")
    logger.info("=" * 80)
    logger.info("")
    
    # Load parameters
    params, model_type = load_fastsimcoal2_parameters()
    
    if params is None:
        sys.exit(1)
    
    # Create demographic history plot
    plot_demographic_history(params, model_type)
    
    # Plot SFS
    plot_sfs_fit()
    
    logger.info("\n" + "=" * 80)
    logger.info("VISUALIZATION COMPLETE")
    logger.info("=" * 80)
    logger.info(f"\nPlots saved to: {PLOT_DIR}")
    logger.info("\nGenerated plots:")
    logger.info("  - demographic_history.png/pdf")
    logger.info("  - sfs_observed.png")
    logger.info("")


if __name__ == "__main__":
    main()
