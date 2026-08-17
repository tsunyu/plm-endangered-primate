#!/usr/bin/env python3
"""
Analyze fastsimcoal2 Demographic Inference Results
===================================================

This script extracts and analyzes demographic parameters from fastsimcoal2
results, including bootstrap confidence intervals.

Input:  fastsimcoal2 best-fit parameters and bootstrap results
Output: Parameter estimates with 95% CI, demographic summary

Author: Demographic Analysis Pipeline
Date: 2026-01-26
"""

import sys
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir
import logging

# Configuration
BASE_DIR = get_base_dir()
MODEL_DIR = BASE_DIR / "output/phase3b_fastsimcoal2/models"
BOOT_DIR = BASE_DIR / "output/phase3b_fastsimcoal2/bootstrap"
OUTPUT_DIR = BASE_DIR / "output/phase3b_fastsimcoal2"

GENERATION_TIME = 10  # years per generation
MUTATION_RATE = 1.36e-8
DEFAULT_NE_SCALING = "Ne"  # fastsimcoal2 size parameters are reported as diploid Ne

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def get_best_model():
    """Determine best model from AIC comparison."""
    comparison_file = BASE_DIR / "output/phase3b_fastsimcoal2/model_comparison/model_comparison.csv"
    
    if not comparison_file.exists():
        logger.error("Model comparison file not found. Please run step 3 first.")
        return None
    
    df = pd.read_csv(comparison_file)
    best_model = df.iloc[0]['Model']
    logger.info(f"Best model: {best_model}")
    
    return best_model


def parse_parameter_line(line, model_name, header=None):
    """
    Parse parameter values from .bestlhoods file.
    
    Returns:
        Dictionary of parameter values
    """
    values = line.strip().split()

    # Known parameter names for each model (used when header is missing)
    param_names_map = {
        'constant_ne': ['NCUR'],
        'single_bottleneck': ['NCUR', 'NBOT', 'NANC', 'TRECOVERY', 'TBOT'],
        'two_consecutive_bottlenecks': [
            'NCUR', 'NBOT2', 'NINTER', 'NBOT1', 'NANC',
            'TRECENT_RECOVERY', 'TRECENT_BOT', 'TOLD_RECOVERY', 'TOLD_BOT'
        ],
        'bottleneck_continuous_decline': ['NCUR', 'NMID', 'NBOT', 'NANC', 'TRECENT', 'TBOT', 'TANC'],
        'bottleneck_recent_contraction': ['NCUR', 'NRECOVER', 'NBOT', 'NANC', 'TRECENT', 'T_REC', 'T_BOT', 'TRECOVERY_OLD', 'TBOT_OLD'],
        'complex_multi_event': ['NCUR', 'N1', 'N2', 'N3', 'NANC', 'T1', 'T2', 'T3', 'T4']
    }

    if header:
        # Header-aware parsing: use column names to extract parameters and MaxEstLhood
        name_to_idx = {name: i for i, name in enumerate(header)}
        params = {}

        # Extract known parameters when present in header
        for pname in param_names_map.get(model_name, []):
            if pname in name_to_idx and name_to_idx[pname] < len(values):
                try:
                    params[pname] = float(values[name_to_idx[pname]])
                except ValueError:
                    continue

        # Likelihood column (prefer MaxEstLhood)
        lhood_idx = None
        for key in name_to_idx:
            if key.lower().startswith("maxestlhood"):
                lhood_idx = name_to_idx[key]
                break
        if lhood_idx is None:
            lhood_idx = len(values) - 1  # fallback

        try:
            params['MaxEstLhood'] = float(values[lhood_idx])
        except ValueError:
            params['MaxEstLhood'] = None

        return params

    # Fallback for older bootstrap files that store only the raw data line.
    # fsc output is normally: parameters..., MaxEstLhood, MaxObsLhood.
    param_names = param_names_map.get(model_name)
    if param_names is None:
        # Reserve up to two trailing likelihood columns when parameter names are unknown.
        n_unknown = max(len(values) - 2, 0)
        param_names = [f'param_{i}' for i in range(n_unknown)]
    n_params = len(param_names)
    param_values = [float(v) for v in values[:n_params]]
    params = dict(zip(param_names, param_values))
    if len(values) > n_params:
        params['MaxEstLhood'] = float(values[n_params])
    else:
        params['MaxEstLhood'] = None
    return params


def get_best_parameters(model_name):
    """Extract best-fit parameters for a model."""
    logger.info(f"Extracting best-fit parameters for {model_name}")
    
    # Try both possible paths (serial and parallel script outputs differ)
    bestlhoods_file = MODEL_DIR / model_name / "best_run" / model_name / f"{model_name}.bestlhoods"
    if not bestlhoods_file.exists():
        # Parallel script puts files directly in best_run/
        bestlhoods_file = MODEL_DIR / model_name / "best_run" / f"{model_name}.bestlhoods"
    
    if not bestlhoods_file.exists():
        logger.error(f"Best run file not found: {bestlhoods_file}")
        return None
    
    with open(bestlhoods_file, 'r') as f:
        lines = [ln.strip() for ln in f.readlines() if ln.strip()]

    if len(lines) < 2:
        logger.error(f"Unexpected .bestlhoods format: {bestlhoods_file}")
        return None

    header = lines[0].split()
    last_line = lines[-1]

    params = parse_parameter_line(last_line, model_name, header=header)
    return params


def get_bootstrap_parameters(model_name):
    """Extract bootstrap parameter distributions."""
    logger.info(f"Loading bootstrap results for {model_name}")
    
    boot_file = BOOT_DIR / model_name / "bootstrap_results.csv"
    
    if not boot_file.exists():
        logger.warning(f"Bootstrap file not found: {boot_file}")
        return None
    
    # Read bootstrap results
    df = pd.read_csv(boot_file)
    
    # Remove failed replicates
    df = df.dropna()
    
    logger.info(f"Loaded {len(df)} successful bootstrap replicates")
    
    # Parse parameter values from each replicate
    boot_params = []
    
    for _, row in df.iterrows():
        param_str = row['parameters']
        if isinstance(param_str, str):
            try:
                # Bootstrap CSV currently stores raw .bestlhoods data line without header;
                # fall back to header-less parsing (same version as original analysis).
                params = parse_parameter_line(param_str, model_name, header=None)
                boot_params.append(params)
            except:
                continue
    
    if len(boot_params) == 0:
        logger.warning("No valid bootstrap parameters found")
        return None
    
    # Convert to DataFrame
    boot_df = pd.DataFrame(boot_params)
    
    return boot_df


def calculate_confidence_intervals(boot_df, param_name, percentile=95):
    """Calculate confidence intervals from bootstrap distribution."""
    if boot_df is None or param_name not in boot_df.columns:
        return None, None
    
    values = boot_df[param_name].values
    
    # Remove outliers (beyond 3 standard deviations)
    mean = np.mean(values)
    std = np.std(values)
    filtered = values[np.abs(values - mean) < 3 * std]
    
    if len(filtered) < 10:
        filtered = values  # Use all if too few remain
    
    alpha = (100 - percentile) / 2
    lower = np.percentile(filtered, alpha)
    upper = np.percentile(filtered, 100 - alpha)
    
    return lower, upper


def convert_time_to_years(time_generations):
    """Convert time in generations to years."""
    return time_generations * GENERATION_TIME


def create_parameter_summary(model_name, best_params, boot_df):
    """Create comprehensive parameter summary."""
    logger.info("Creating parameter summary")
    
    summary_file = OUTPUT_DIR / "parameter_estimates.txt"
    
    with open(summary_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("DEMOGRAPHIC PARAMETER ESTIMATES\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"Best-supported model: {model_name}\n")
        f.write(f"Generation time: {GENERATION_TIME} years\n")
        f.write(f"Mutation rate: {MUTATION_RATE:.2e} per bp per generation\n\n")
        f.write("Size parameter convention: N parameters are reported as diploid effective population size (Ne)\n")
        f.write("Likelihood convention: model selection and best-run choice use MaxEstLhood\n\n")
        
        f.write("-" * 80 + "\n")
        f.write("PARAMETER ESTIMATES\n")
        f.write("-" * 80 + "\n\n")
        
        f.write(f"{'Parameter':<20} {'Estimate':<15} {'95% CI Lower':<15} {'95% CI Upper':<15} {'Unit':<10}\n")
        f.write("-" * 80 + "\n")
        
        # Extract and display each parameter
        for param_name, value in best_params.items():
            if param_name == 'MaxEstLhood':
                continue
            
            # Get confidence intervals from bootstrap
            if boot_df is not None:
                ci_lower, ci_upper = calculate_confidence_intervals(boot_df, param_name)
            else:
                ci_lower, ci_upper = None, None
            
            # Determine unit
            if param_name.startswith('N'):
                unit = 'diploid Ne'
                ci_l_str = f"{ci_lower:.0f}" if ci_lower else "N/A"
                ci_u_str = f"{ci_upper:.0f}" if ci_upper else "N/A"
                val_str = f"{value:.0f}"
            elif param_name.startswith('T'):
                unit = 'generations'
                ci_l_str = f"{ci_lower:.0f}" if ci_lower else "N/A"
                ci_u_str = f"{ci_upper:.0f}" if ci_upper else "N/A"
                val_str = f"{value:.0f}"
            else:
                unit = 'unitless'
                ci_l_str = f"{ci_lower:.6f}" if ci_lower else "N/A"
                ci_u_str = f"{ci_upper:.6f}" if ci_upper else "N/A"
                val_str = f"{value:.6f}"
            
            f.write(f"{param_name:<20} {val_str:<15} {ci_l_str:<15} {ci_u_str:<15} {unit:<10}\n")
        
        f.write("\n")
        f.write("-" * 80 + "\n")
        f.write("TIME CONVERSION (GENERATIONS → YEARS)\n")
        f.write("-" * 80 + "\n\n")
        
        # Convert time parameters to years
        for param_name, value in best_params.items():
            if param_name.startswith('T'):
                years = convert_time_to_years(value)
                
                if boot_df is not None:
                    ci_lower, ci_upper = calculate_confidence_intervals(boot_df, param_name)
                    if ci_lower and ci_upper:
                        years_lower = convert_time_to_years(ci_lower)
                        years_upper = convert_time_to_years(ci_upper)
                        f.write(f"{param_name}: {years:.0f} years ago (95% CI: {years_lower:.0f} - {years_upper:.0f})\n")
                    else:
                        f.write(f"{param_name}: {years:.0f} years ago\n")
                else:
                    f.write(f"{param_name}: {years:.0f} years ago\n")
        
        f.write("\n")
        f.write("-" * 80 + "\n")
        f.write("BIOLOGICAL INTERPRETATION\n")
        f.write("-" * 80 + "\n\n")
        
        # Model-specific interpretation
        if model_name == 'single_bottleneck':
            ncur = best_params.get('NCUR', 0)
            nbot = best_params.get('NBOT', 0)
            nanc = best_params.get('NANC', 0)
            tbot = best_params.get('TBOT', 0)
            trecover = best_params.get('TRECOVERY', 0)
            
            f.write("Single-bottleneck model (time in generations before present):\n")
            f.write(f"  Current Ne: {ncur:.0f} individuals\n")
            f.write(f"  Bottleneck Ne: {nbot:.0f} individuals\n")
            f.write(f"  Ancestral Ne: {nanc:.0f} individuals\n")
            f.write(f"  Bottleneck onset (TBOT): {convert_time_to_years(tbot):.0f} years ago\n")
            f.write(f"  Recovery start (TRECOVERY): {convert_time_to_years(trecover):.0f} years ago\n\n")
            
            severity = nanc / nbot if nbot > 0 else 0
            recovery = ncur / nbot if nbot > 0 else 0
            
            f.write(f"  Bottleneck severity: {severity:.2f}× decline\n")
            f.write(f"  Recovery magnitude: {recovery:.2f}× increase from bottleneck\n")
            f.write(f"  Current vs ancestral: {ncur/nanc*100:.1f}% of ancestral size\n")
        
        elif model_name == 'bottleneck_recent_contraction':
            ncur = best_params.get('NCUR', 0)
            nrecover = best_params.get('NRECOVER', 0)
            nbot = best_params.get('NBOT', 0)
            nanc = best_params.get('NANC', 0)
            trecent = best_params.get('TRECENT', 0)
            trecovery_old = best_params.get('T_REC', best_params.get('TRECOVERY_OLD', 0))
            tbot_old = best_params.get('T_BOT', best_params.get('TBOT_OLD', 0))
            
            f.write("Bottleneck + recent contraction model:\n")
            f.write(f"  Current Ne: {ncur:.0f} individuals\n")
            f.write(f"  Pre-recent-contraction Ne: {nrecover:.0f} individuals\n")
            f.write(f"  Bottleneck Ne: {nbot:.0f} individuals\n")
            f.write(f"  Ancestral Ne: {nanc:.0f} individuals\n")
            f.write(f"  Recent contraction time (TRECENT): {convert_time_to_years(trecent):.0f} years ago\n")
            f.write(f"  Ancient recovery start (TRECOVERY_OLD): {convert_time_to_years(trecovery_old):.0f} years ago\n")
            f.write(f"  Ancient bottleneck onset (TBOT_OLD): {convert_time_to_years(tbot_old):.0f} years ago\n\n")

            if nrecover > 0:
                f.write(f"  Recent contraction magnitude: {ncur/nrecover*100:.1f}% of pre-recent-contraction size\n")
            if nbot > 0:
                f.write(f"  Ancient bottleneck severity: {nanc/nbot:.2f}× decline\n")
        
        f.write("\n" + "=" * 80 + "\n")
    
    logger.info(f"Parameter summary saved: {summary_file}")


def plot_bootstrap_distributions(model_name, best_params, boot_df):
    """Plot bootstrap parameter distributions."""
    if boot_df is None:
        logger.warning("No bootstrap data available for plotting")
        return
    
    logger.info("Creating bootstrap distribution plots")
    
    # Get parameter names (exclude likelihood)
    param_names = [p for p in best_params.keys() if p != 'MaxEstLhood']
    n_params = len(param_names)
    
    # Create subplot grid
    n_cols = min(3, n_params)
    n_rows = (n_params + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    if n_params == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    for idx, param_name in enumerate(param_names):
        ax = axes[idx]
        
        if param_name not in boot_df.columns:
            continue
        
        values = boot_df[param_name].values
        best_value = best_params[param_name]
        
        # Calculate CI
        ci_lower, ci_upper = calculate_confidence_intervals(boot_df, param_name)
        
        # Plot histogram
        ax.hist(values, bins=30, color='steelblue', alpha=0.7, edgecolor='black')
        
        # Mark best estimate
        ax.axvline(best_value, color='red', linewidth=2, linestyle='--', 
                  label=f'Best: {best_value:.1f}')
        
        # Mark CI
        if ci_lower and ci_upper:
            ax.axvline(ci_lower, color='orange', linewidth=1.5, linestyle=':', 
                      label=f'95% CI: [{ci_lower:.1f}, {ci_upper:.1f}]')
            ax.axvline(ci_upper, color='orange', linewidth=1.5, linestyle=':')
        
        ax.set_xlabel(f'{param_name}', fontweight='bold')
        ax.set_ylabel('Frequency', fontweight='bold')
        ax.set_title(f'Bootstrap Distribution: {param_name}', fontweight='bold')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3, axis='y')
    
    # Hide unused subplots
    for idx in range(n_params, len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    plot_file = OUTPUT_DIR / "plots" / f"{model_name}_bootstrap_distributions.png"
    plot_file.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Bootstrap plot saved: {plot_file}")


def main():
    """Main workflow."""
    logger.info("=" * 80)
    logger.info("DEMOGRAPHIC PARAMETER ANALYSIS")
    logger.info("=" * 80)
    logger.info("")
    
    # Get best model
    model_name = get_best_model()
    if model_name is None:
        sys.exit(1)
    
    # Get best-fit parameters (raw, as in fastsimcoal2 output)
    best_params = get_best_parameters(model_name)
    if best_params is None:
        logger.error("Could not extract best-fit parameters")
        sys.exit(1)
    
    logger.info(f"\nBest-fit parameters for {model_name}:")
    for param, value in best_params.items():
        logger.info(f"  {param}: {value}")
    
    # Get bootstrap results
    boot_df = get_bootstrap_parameters(model_name)

    # Optionally apply Ne scaling if NE_SCALING_NOTE.txt indicates 2Ne
    scaling_note = BASE_DIR / "output/phase3b_fastsimcoal2/NE_SCALING_NOTE.txt"
    scaling_type = DEFAULT_NE_SCALING
    if scaling_note.exists():
        with open(scaling_note, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith("SCALING="):
                    scaling_type = line.split("=", 1)[1].strip()
                    break

    if scaling_type == "2Ne":
        logger.info("Applying Ne scaling: interpreting NCUR etc. as 2*Ne → dividing N-parameters by 2")

        def _scale_ne_params(param_dict):
            scaled = param_dict.copy()
            for key, val in param_dict.items():
                if key.startswith("N") and isinstance(val, (int, float)):
                    scaled[key] = val / 2.0
            return scaled

        best_params = _scale_ne_params(best_params)

        if boot_df is not None:
            for col in boot_df.columns:
                if col.startswith("N"):
                    boot_df[col] = boot_df[col] / 2.0
    
    # Create summary
    create_parameter_summary(model_name, best_params, boot_df)
    
    # Plot bootstrap distributions
    if boot_df is not None:
        plot_bootstrap_distributions(model_name, best_params, boot_df)
    
    # Save parameter estimates to CSV
    param_csv = OUTPUT_DIR / "parameter_estimates.csv"
    
    param_data = []
    for param_name, value in best_params.items():
        if param_name == 'MaxEstLhood':
            continue
        
        row = {'Parameter': param_name, 'Estimate': value}
        
        if boot_df is not None:
            ci_lower, ci_upper = calculate_confidence_intervals(boot_df, param_name)
            row['CI_Lower'] = ci_lower
            row['CI_Upper'] = ci_upper
        else:
            row['CI_Lower'] = None
            row['CI_Upper'] = None
        
        # Add time conversion for T parameters
        if param_name.startswith('T'):
            row['Estimate_Years'] = convert_time_to_years(value)
            if row['CI_Lower']:
                row['CI_Lower_Years'] = convert_time_to_years(row['CI_Lower'])
                row['CI_Upper_Years'] = convert_time_to_years(row['CI_Upper'])
        
        param_data.append(row)
    
    pd.DataFrame(param_data).to_csv(param_csv, index=False)
    logger.info(f"\nParameter estimates saved to: {param_csv}")
    
    logger.info("\n" + "=" * 80)
    logger.info("ANALYSIS COMPLETE")
    logger.info("=" * 80)
    logger.info(f"\nOutput files:")
    logger.info(f"  - {OUTPUT_DIR / 'parameter_estimates.txt'}")
    logger.info(f"  - {OUTPUT_DIR / 'parameter_estimates.csv'}")
    if boot_df is not None:
        logger.info(f"  - {OUTPUT_DIR / 'plots' / f'{model_name}_bootstrap_distributions.png'}")
    logger.info("")


if __name__ == "__main__":
    main()
