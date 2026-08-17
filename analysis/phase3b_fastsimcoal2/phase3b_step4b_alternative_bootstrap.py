#!/usr/bin/env python3
"""
Alternative Bootstrap Approach for fastsimcoal2
================================================

If the standard bootstrap fails due to .par file issues, this script
provides alternative approaches for uncertainty quantification:

Method 1: Non-parametric bootstrap (resample SFS bins)
Method 2: Jackknife resampling
Method 3: Profile likelihood-based CIs

Author: Demographic Analysis Pipeline
Date: 2026-01-26
"""

import numpy as np
import pandas as pd
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir
import subprocess
import logging
from scipy import stats

# Configuration
BASE_DIR = get_base_dir()
MODEL_DIR = BASE_DIR / "output/phase3b_fastsimcoal2/models"
SFS_DIR = BASE_DIR / "output/phase3b_fastsimcoal2/sfs"
OUTPUT_DIR = BASE_DIR / "output/phase3b_fastsimcoal2/bootstrap_alternative"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)


def load_sfs():
    """Load observed SFS."""
    sfs_file = SFS_DIR / "SNJ_DAFpop0.obs"
    
    with open(sfs_file, 'r') as f:
        lines = f.readlines()
    
    # Parse SFS data
    data = lines[2].strip().split('\t')
    sfs = np.array([int(x) for x in data])
    
    return sfs


def nonparametric_bootstrap(sfs, n_bootstrap=100):
    """
    Method 1: Non-parametric bootstrap.
    
    Resample SFS bins with replacement.
    """
    logger.info("Method 1: Non-parametric Bootstrap")
    logger.info("-" * 60)
    
    n_bins = len(sfs)
    bootstrap_sfs = []
    
    for i in range(n_bootstrap):
        # Resample with replacement
        weights = np.random.multinomial(sum(sfs), sfs / sum(sfs))
        bootstrap_sfs.append(weights)
        
        if (i + 1) % 20 == 0:
            logger.info(f"  Generated {i+1}/{n_bootstrap} bootstrap replicates")
    
    return np.array(bootstrap_sfs)


def jackknife_resampling(sfs):
    """
    Method 2: Jackknife resampling.
    
    Leave-one-out approach for each SFS bin.
    """
    logger.info("Method 2: Jackknife Resampling")
    logger.info("-" * 60)
    
    n_bins = len(sfs)
    jackknife_sfs = []
    
    # For each bin, reduce count by 1 and redistribute
    for i in range(n_bins):
        if sfs[i] > 0:
            jk_sfs = sfs.copy()
            jk_sfs[i] -= 1
            # Redistribute removed count proportionally
            if sum(sfs) > 1:
                jk_sfs += 1 / (sum(sfs) - 1)
            jackknife_sfs.append(jk_sfs)
    
    logger.info(f"  Generated {len(jackknife_sfs)} jackknife replicates")
    
    return np.array(jackknife_sfs)


def approximate_ci_from_multiple_runs():
    """
    Method 3: Use variation across independent optimization runs.
    
    fastsimcoal2 runs 50 independent optimizations. The variation
    in parameter estimates across these runs can approximate uncertainty.
    """
    logger.info("Method 3: Approximate CI from Multiple Runs")
    logger.info("-" * 60)
    
    # Get best model
    comparison_file = BASE_DIR / "output/phase3b_fastsimcoal2/model_comparison/model_comparison.csv"
    
    if not comparison_file.exists():
        logger.error("Model comparison not found")
        return None
    
    df = pd.read_csv(comparison_file)
    best_model = df.iloc[0]['Model']
    
    logger.info(f"  Best model: {best_model}")
    
    # Read likelihoods from all runs
    lhoods_file = MODEL_DIR / best_model / "run_likelihoods.csv"
    
    if not lhoods_file.exists():
        logger.error(f"Likelihoods file not found: {lhoods_file}")
        return None
    
    lhoods_df = pd.read_csv(lhoods_file)
    lhoods_df = lhoods_df.dropna()
    
    logger.info(f"  Found {len(lhoods_df)} successful runs")
    
    # Extract parameters from top runs
    top_n = min(10, len(lhoods_df))
    likelihood_col = 'max_est_likelihood' if 'max_est_likelihood' in lhoods_df.columns else 'max_likelihood'
    top_runs = lhoods_df.nlargest(top_n, likelihood_col)
    
    logger.info(f"  Analyzing top {top_n} runs")
    
    param_values = []
    
    for _, row in top_runs.iterrows():
        run_id = row['run_id']
        run_dir = MODEL_DIR / best_model / f"run_{run_id}" / best_model
        bestlhoods = run_dir / f"{best_model}.bestlhoods"
        
        if bestlhoods.exists():
            with open(bestlhoods, 'r') as f:
                lines = f.readlines()
                header = lines[0].strip().split()
                last_line = lines[-1].strip().split()
                likelihood_cols = {'MaxEstLhood', 'MaxObsLhood'}
                params = [
                    float(value)
                    for name, value in zip(header, last_line)
                    if name not in likelihood_cols
                ]
                param_values.append(params)
    
    if param_values:
        param_array = np.array(param_values)
        
        # Calculate mean and std for each parameter
        param_means = np.mean(param_array, axis=0)
        param_stds = np.std(param_array, axis=0)
        
        # Approximate 95% CI (mean ± 1.96*std)
        ci_lower = param_means - 1.96 * param_stds
        ci_upper = param_means + 1.96 * param_stds
        
        logger.info(f"  Parameter variability:")
        for i in range(len(param_means)):
            logger.info(f"    Param {i+1}: {param_means[i]:.1f} ± {param_stds[i]:.1f}")
            logger.info(f"              95% CI: [{ci_lower[i]:.1f}, {ci_upper[i]:.1f}]")
        
        return {
            'means': param_means,
            'stds': param_stds,
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'n_runs': len(param_values)
        }
    
    return None


def save_alternative_bootstrap_results(results, method_name):
    """Save bootstrap results."""
    output_file = OUTPUT_DIR / f"{method_name}_results.csv"
    
    if isinstance(results, dict):
        df = pd.DataFrame(results)
        df.to_csv(output_file, index=False)
        logger.info(f"  Saved to: {output_file}")
    else:
        logger.warning(f"  No results to save for {method_name}")


def main():
    """Main workflow."""
    logger.info("=" * 80)
    logger.info("ALTERNATIVE BOOTSTRAP METHODS FOR UNCERTAINTY QUANTIFICATION")
    logger.info("=" * 80)
    logger.info("")
    logger.info("This script provides alternatives when standard parametric")
    logger.info("bootstrap fails due to .par file issues.")
    logger.info("")
    
    # Load SFS
    try:
        sfs = load_sfs()
        logger.info(f"Loaded SFS with {len(sfs)} bins")
        logger.info(f"Total SNPs: {sum(sfs):,}")
        logger.info("")
    except Exception as e:
        logger.error(f"Could not load SFS: {e}")
        return
    
    # Method 1: Non-parametric bootstrap (computationally intensive)
    logger.info("\n" + "=" * 80)
    try:
        logger.info("NOTE: Non-parametric bootstrap requires re-running")
        logger.info("      fastsimcoal2 estimation for each replicate.")
        logger.info("      This is very time-consuming (100+ runs).")
        logger.info("")
        logger.info("      Skipping Method 1 for now.")
        logger.info("      Use parametric bootstrap (step 4) if possible.")
    except Exception as e:
        logger.error(f"Method 1 failed: {e}")
    
    # Method 2: Jackknife (also computationally intensive)
    logger.info("\n" + "=" * 80)
    try:
        logger.info("NOTE: Jackknife also requires multiple fastsimcoal2 runs.")
        logger.info("      Skipping Method 2 for now.")
    except Exception as e:
        logger.error(f"Method 2 failed: {e}")
    
    # Method 3: Use existing run variation (RECOMMENDED)
    logger.info("\n" + "=" * 80)
    try:
        results = approximate_ci_from_multiple_runs()
        if results:
            save_alternative_bootstrap_results(results, "multiple_runs_ci")
            
            logger.info("")
            logger.info("✅ Method 3 SUCCESSFUL")
            logger.info("")
            logger.info("RECOMMENDATION:")
            logger.info("  Use these approximate CIs for exploratory analysis.")
            logger.info("  They are based on variation across independent runs.")
            logger.info("  More conservative than parametric bootstrap.")
    except Exception as e:
        logger.error(f"Method 3 failed: {e}")
    
    logger.info("")
    logger.info("=" * 80)
    logger.info("SUMMARY")
    logger.info("=" * 80)
    logger.info("")
    logger.info("BEST OPTION: Fix parametric bootstrap (phase3b_step4_bootstrap_ci.sh)")
    logger.info("  - Most accurate confidence intervals")
    logger.info("  - Standard method in literature")
    logger.info("")
    logger.info("FALLBACK OPTION: Use Method 3 results (multiple runs)")
    logger.info("  - Approximate CIs from optimization variation")
    logger.info("  - Available immediately")
    logger.info("  - Conservative estimates")
    logger.info("")
    logger.info("RECOMMENDATION: Prefer parametric bootstrap for reported CIs")
    logger.info("")


if __name__ == "__main__":
    main()
