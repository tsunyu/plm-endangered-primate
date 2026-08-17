#!/usr/bin/env python3
"""
Compare Demographic Models Using AIC
=====================================

This script compares different demographic models estimated by fastsimcoal2
using Akaike Information Criterion (AIC) to determine the best-supported model.

Input:  fastsimcoal2 results for multiple models
Output: Model comparison table, AIC weights, best model selection

Author: Demographic Analysis Pipeline
Date: 2026-01-26
"""

import sys
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
OUTPUT_DIR = BASE_DIR / "output/phase3b_fastsimcoal2/model_comparison"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(OUTPUT_DIR / "model_comparison.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Model information (name: number of parameters)
MODELS = {
    'constant_ne': 1,
    'single_bottleneck': 5,
    'two_consecutive_bottlenecks': 9,
    'bottleneck_continuous_decline': 7,
    'bottleneck_recent_contraction': 7,
    'complex_multi_event': 9
}


def extract_likelihood(model_name):
    """
    Extract maximum likelihood from fastsimcoal2 output.
    
    Returns:
        float: Maximum likelihood value, or None if not found
    """
    # Try both possible paths (serial and parallel script outputs differ)
    bestlhoods_file = MODEL_DIR / model_name / "best_run" / model_name / f"{model_name}.bestlhoods"
    if not bestlhoods_file.exists():
        # Parallel script puts files directly in best_run/
        bestlhoods_file = MODEL_DIR / model_name / "best_run" / f"{model_name}.bestlhoods"
    
    if not bestlhoods_file.exists():
        logger.warning(f"No .bestlhoods file found for {model_name}")
        return None
    
    try:
        with open(bestlhoods_file, 'r') as f:
            lines = [ln.strip() for ln in f.readlines() if ln.strip()]
            if len(lines) < 2:
                logger.error(f"Unexpected .bestlhoods format for {model_name}")
                return None
            header = lines[0].split()
            last_line = lines[-1].split()

            # Prefer column named MaxEstLhood (robust to column order)
            idx = None
            for i, name in enumerate(header):
                if name.lower().startswith("maxestlhood"):
                    idx = i
                    break

            if idx is None:
                # Fallback for older files: fsc usually writes MaxEstLhood
                # immediately before the final MaxObsLhood column.
                idx = len(last_line) - 2 if len(last_line) >= 2 else len(last_line) - 1

            max_lhood = float(last_line[idx])
            return max_lhood
    except Exception as e:
        logger.error(f"Error extracting likelihood for {model_name}: {e}")
        return None


def calculate_aic(log_likelihood, n_params):
    """
    Calculate Akaike Information Criterion.
    
    AIC = 2k - 2ln(L)
    where k = number of parameters, L = likelihood
    """
    return 2 * n_params - 2 * log_likelihood


def calculate_aic_weights(aic_values):
    """
    Calculate Akaike weights for model comparison.
    
    Weight_i = exp(-0.5 * ΔAIC_i) / Σ exp(-0.5 * ΔAIC_j)
    """
    # Calculate ΔAIC (difference from minimum AIC)
    min_aic = np.min(aic_values)
    delta_aic = aic_values - min_aic
    
    # Calculate relative likelihoods
    rel_likelihoods = np.exp(-0.5 * delta_aic)
    
    # Normalize to get weights
    weights = rel_likelihoods / np.sum(rel_likelihoods)
    
    return delta_aic, weights


def create_comparison_table(results_df):
    """Create formatted model comparison table."""
    logger.info("Creating model comparison table")
    
    table_path = OUTPUT_DIR / "model_comparison_table.txt"
    
    with open(table_path, 'w') as f:
        f.write("=" * 100 + "\n")
        f.write("DEMOGRAPHIC MODEL COMPARISON - AIC ANALYSIS\n")
        f.write("=" * 100 + "\n\n")
        
        f.write("Models compared:\n")
        for model in results_df['Model']:
            f.write(f"  - {model}\n")
        f.write(f"\nTotal models: {len(results_df)}\n\n")
        
        f.write("-" * 100 + "\n")
        f.write("MODEL COMPARISON RESULTS\n")
        f.write("-" * 100 + "\n\n")
        
        # Create formatted table
        f.write(f"{'Model':<25} {'Params':<8} {'MaxEstLhood':<15} {'AIC':<15} {'ΔAIC':<12} {'Weight':<10}\n")
        f.write("-" * 100 + "\n")
        
        for _, row in results_df.iterrows():
            f.write(f"{row['Model']:<25} "
                   f"{row['Parameters']:<8} "
                   f"{row['MaxEstLhood']:<15.2f} "
                   f"{row['AIC']:<15.2f} "
                   f"{row['DeltaAIC']:<12.2f} "
                   f"{row['Weight']:<10.4f}\n")
        
        f.write("\n")
        f.write("-" * 100 + "\n")
        f.write("INTERPRETATION\n")
        f.write("-" * 100 + "\n\n")
        
        # Best model
        best_model = results_df.iloc[0]
        f.write(f"Best-supported model: {best_model['Model']}\n")
        f.write(f"  AIC weight: {best_model['Weight']:.4f}\n")
        f.write(f"  Maximum estimated likelihood: {best_model['MaxEstLhood']:.2f}\n")
        f.write(f"  Number of parameters: {best_model['Parameters']}\n\n")
        
        # Evidence ratios
        f.write("Evidence strength (Burnham & Anderson guidelines):\n")
        f.write("  ΔAIC < 2:   Substantial support\n")
        f.write("  ΔAIC 4-7:   Considerably less support\n")
        f.write("  ΔAIC > 10:  Essentially no support\n\n")
        
        f.write("Model support summary:\n")
        for _, row in results_df.iterrows():
            delta = row['DeltaAIC']
            if delta < 2:
                support = "STRONG support"
            elif delta < 4:
                support = "Moderate support"
            elif delta < 7:
                support = "Weak support"
            elif delta < 10:
                support = "Very weak support"
            else:
                support = "No support"
            
            f.write(f"  {row['Model']:<25} ΔAIC = {delta:>6.2f}  →  {support}\n")
        
        f.write("\n")
        f.write("-" * 100 + "\n")
        f.write("BIOLOGICAL INTERPRETATION\n")
        f.write("-" * 100 + "\n\n")
        
        model_descriptions = {
            'constant_ne': 'No demographic change over time',
            'single_bottleneck': 'Single severe bottleneck followed by recovery',
            'two_consecutive_bottlenecks': 'Two independent bottlenecks from repeated deterioration events',
            'bottleneck_continuous_decline': 'Bottleneck followed by persistent decline without recovery',
            'bottleneck_recent_contraction': 'Ancient bottleneck with short recovery and recent (~1.3 ka) contraction',
            'complex_multi_event': 'Complex model with >=3 demographic size changes'
        }
        
        f.write(f"Best model ({best_model['Model']}):\n")
        f.write(f"  {model_descriptions.get(best_model['Model'], 'No description')}\n\n")
        
        if best_model['Weight'] > 0.9:
            f.write("Inference: VERY STRONG support for this demographic scenario\n")
        elif best_model['Weight'] > 0.7:
            f.write("Inference: STRONG support for this demographic scenario\n")
        elif best_model['Weight'] > 0.5:
            f.write("Inference: MODERATE support for this demographic scenario\n")
        else:
            f.write("Inference: Model uncertainty remains - consider model averaging\n")
        
        f.write("\n" + "=" * 100 + "\n")
    
    logger.info(f"Table saved: {table_path}")


def plot_model_comparison(results_df):
    """Create visualization of model comparison results."""
    logger.info("Creating model comparison plots")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Plot 1: AIC values
    ax = axes[0, 0]
    colors = ['green' if i == 0 else 'steelblue' for i in range(len(results_df))]
    ax.barh(results_df['Model'], results_df['AIC'], color=colors, alpha=0.7, edgecolor='black')
    ax.set_xlabel('AIC (lower is better)', fontweight='bold')
    ax.set_ylabel('Model', fontweight='bold')
    ax.set_title('Akaike Information Criterion', fontweight='bold')
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3, axis='x')
    
    # Plot 2: ΔAIC
    ax = axes[0, 1]
    delta_original = results_df['DeltaAIC'].values
    delta_plot = np.log10(delta_original + 1.0)  # compress extreme values
    colors_delta = ['green' if d < 2 else 'orange' if d < 7 else 'red' 
                    for d in delta_original]
    ax.barh(results_df['Model'], delta_plot, color=colors_delta, alpha=0.7, edgecolor='black')
    ax.axvline(x=np.log10(2 + 1.0), color='green', linestyle='--', linewidth=2, alpha=0.5, label='ΔAIC=2')
    ax.axvline(x=np.log10(7 + 1.0), color='orange', linestyle='--', linewidth=2, alpha=0.5, label='ΔAIC=7')
    ax.set_xlabel('log10(ΔAIC + 1)', fontweight='bold')
    ax.set_ylabel('Model', fontweight='bold')
    ax.set_title('ΔAIC (compressed scale for visibility)', fontweight='bold')
    ax.invert_yaxis()
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis='x')

    # Show original ΔAIC values at bar ends
    for i, d in enumerate(delta_original):
        x = np.log10(d + 1.0)
        label = f'{d:.1f}' if d < 1000 else f'{d:.1e}'
        ax.text(x + 0.03, i, label, va='center', fontsize=8)
    
    # Plot 3: Akaike weights
    ax = axes[1, 0]
    weights = results_df['Weight'].values
    weights_plot = np.clip(weights, 1e-12, None)
    ax.bar(range(len(results_df)), weights_plot, color='purple', alpha=0.7, edgecolor='black')
    ax.set_xticks(range(len(results_df)))
    ax.set_xticklabels(results_df['Model'], rotation=45, ha='right')
    ax.set_ylabel('Akaike Weight (log scale)', fontweight='bold')
    ax.set_title('Model Probability (small weights now visible)', fontweight='bold')
    ax.set_yscale('log')
    ax.set_ylim(1e-12, 1)
    ax.grid(True, alpha=0.3, axis='y')

    # Add readable labels for each weight
    for i, w in enumerate(weights):
        label = f'{w:.3f}' if w >= 0.01 else f'{w:.1e}'
        if w >= 0.9:
            # Keep top labels inside the plotting area to avoid title overlap.
            y = max(weights_plot[i] / 1.35, 1.5e-12)
            ax.text(i, y, label, ha='center', va='top', fontsize=8)
        else:
            y = max(weights_plot[i] * 1.25, 1.5e-12)
            ax.text(i, y, label, ha='center', va='bottom', fontsize=8)
    
    # Plot 4: Likelihood vs Parameters (overfitting check)
    ax = axes[1, 1]
    ax.scatter(results_df['Parameters'], results_df['MaxEstLhood'],
              s=results_df['Weight']*500, c=results_df['DeltaAIC'], 
              cmap='RdYlGn_r', alpha=0.7, edgecolor='black', linewidth=2)
    
    # Add labels for each point with deterministic offsets to reduce overlap
    offsets = [(-18, 8), (10, 12), (-16, -12), (12, -10), (18, 6), (-12, 14), (8, -14), (-20, 0)]
    for _, row in results_df.iterrows():
        idx = int(_)
        dx, dy = offsets[idx % len(offsets)]
        ax.annotate(row['Model'], 
                   (row['Parameters'], row['MaxEstLhood']),
                   textcoords='offset points',
                   xytext=(dx, dy),
                   fontsize=8,
                   bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.7))
    
    ax.set_xlabel('Number of Parameters', fontweight='bold')
    ax.set_ylabel('MaxEstLhood', fontweight='bold')
    ax.set_title('Likelihood vs Model Complexity\n(size = weight, color = ΔAIC)', fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    cbar = plt.colorbar(ax.collections[0], ax=ax)
    cbar.set_label('ΔAIC', fontweight='bold')
    
    plt.tight_layout()
    plot_path = OUTPUT_DIR / "model_comparison_plots.png"
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Plot saved: {plot_path}")


def main():
    """Main workflow."""
    logger.info("=" * 80)
    logger.info("DEMOGRAPHIC MODEL COMPARISON")
    logger.info("=" * 80)
    
    # Extract likelihoods for all models
    results = []
    
    for model_name, n_params in MODELS.items():
        logger.info(f"Processing model: {model_name}")
        
        max_lhood = extract_likelihood(model_name)
        
        if max_lhood is not None:
            aic = calculate_aic(max_lhood, n_params)
            results.append({
                'Model': model_name,
                'Parameters': n_params,
                'MaxEstLhood': max_lhood,
                'AIC': aic
            })
            logger.info(f"  MaxEstLhood: {max_lhood:.2f}, AIC: {aic:.2f}")
        else:
            logger.warning(f"  Skipping {model_name} - no results found")
    
    if len(results) == 0:
        logger.error("No model results found. Please run phase3b_step2_run_fastsimcoal2.sh first")
        sys.exit(1)
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Calculate ΔAIC and weights
    delta_aic, weights = calculate_aic_weights(df['AIC'].values)
    df['DeltaAIC'] = delta_aic
    df['Weight'] = weights
    
    # Sort by AIC (best model first)
    df = df.sort_values('AIC').reset_index(drop=True)
    
    # Save to CSV
    csv_path = OUTPUT_DIR / "model_comparison.csv"
    df.to_csv(csv_path, index=False, float_format='%.4f')
    logger.info(f"\nResults saved to: {csv_path}")
    
    # Create formatted table
    create_comparison_table(df)
    
    # Create plots
    plot_model_comparison(df)
    
    # Display summary
    logger.info("\n" + "=" * 80)
    logger.info("MODEL COMPARISON SUMMARY")
    logger.info("=" * 80)
    
    best_model = df.iloc[0]
    logger.info(f"\nBest model: {best_model['Model']}")
    logger.info(f"  AIC: {best_model['AIC']:.2f}")
    logger.info(f"  Weight: {best_model['Weight']:.4f}")
    logger.info(f"  Parameters: {best_model['Parameters']}")
    
    logger.info(f"\nAll models:")
    for _, row in df.iterrows():
        logger.info(f"  {row['Model']:<25} ΔAIC={row['DeltaAIC']:>6.2f}  Weight={row['Weight']:>6.4f}")
    
    logger.info("\n" + "=" * 80)
    logger.info("ANALYSIS COMPLETE")
    logger.info("=" * 80)
    logger.info(f"\nOutput directory: {OUTPUT_DIR}")
    logger.info("\nNext steps:")
    logger.info("  1. Review model_comparison_table.txt")
    logger.info("  2. Run phase3b_step4_bootstrap_ci.sh for confidence intervals")
    logger.info("  3. Run phase3b_step5_analyze_results.py for parameter interpretation")
    logger.info("")


if __name__ == "__main__":
    main()
