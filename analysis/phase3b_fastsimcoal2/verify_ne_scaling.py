#!/usr/bin/env python3
"""
Verify Ne Parameter Scaling in fastsimcoal2
============================================

This script helps verify whether fastsimcoal2 outputs Ne or 2*Ne.

Approach:
1. Compare with known values (census size)
2. Suggest additional verification methods

Author: Demographic Analysis Pipeline
Date: 2026-01-26
"""

import numpy as np
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir
import logging

# Configuration
BASE_DIR = get_base_dir()
FSC_RESULTS = BASE_DIR / "output/phase3b_fastsimcoal2/models"

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)


def get_fsc_best_ncur():
    """Extract best NCUR from fastsimcoal2 results using header-aware parsing."""
    # Find best model
    comparison_file = BASE_DIR / "output/phase3b_fastsimcoal2/model_comparison/model_comparison.csv"
    
    if not comparison_file.exists():
        logger.error("Model comparison not found. Run pipeline first.")
        return None
    
    import pandas as pd
    df = pd.read_csv(comparison_file)
    best_model = df.iloc[0]['Model']
    
    logger.info(f"Best model: {best_model}")
    
    # Read bestlhoods file
    bestlhoods = FSC_RESULTS / best_model / "best_run" / best_model / f"{best_model}.bestlhoods"
    if not bestlhoods.exists():
        bestlhoods = FSC_RESULTS / best_model / "best_run" / f"{best_model}.bestlhoods"
    
    if not bestlhoods.exists():
        logger.error(f"Bestlhoods file not found: {bestlhoods}")
        return None
    
    with open(bestlhoods, 'r') as f:
        lines = [ln.strip() for ln in f.readlines() if ln.strip()]
        if len(lines) < 2:
            logger.error(f"Unexpected .bestlhoods format: {bestlhoods}")
            return None
        header = lines[0].split()
        data = lines[-1].split()

        # Try to find NCUR column by name
        idx = None
        for i, name in enumerate(header):
            if name.upper() == "NCUR":
                idx = i
                break

        if idx is None:
            # Fallback: assume first column is NCUR (old behaviour)
            idx = 0

        try:
            ncur = float(data[idx])
        except (IndexError, ValueError):
            logger.error(f"Could not parse NCUR from {bestlhoods}")
            return None
    
    logger.info(f"fastsimcoal2 NCUR: {ncur:,.0f}")
    return ncur




def compare_with_census():
    """Compare with known census size."""
    logger.info("\n" + "="*60)
    logger.info("COMPARISON WITH CENSUS SIZE")
    logger.info("="*60)
    
    # Known census size from literature
    census_size = 1200
    logger.info(f"Known census size: {census_size}")
    
    # Expected Ne/Nc ratio for mammals: 0.1 - 0.5
    expected_ne_low = census_size * 0.1
    expected_ne_high = census_size * 0.5
    
    logger.info(f"Expected Ne range: {expected_ne_low:.0f} - {expected_ne_high:.0f}")
    logger.info(f"  (assuming Ne/Nc = 0.1 to 0.5)")
    
    return expected_ne_low, expected_ne_high


def analyze_scaling():
    """Analyze whether results are Ne or 2*Ne."""
    logger.info("=" * 80)
    logger.info("Ne PARAMETER SCALING VERIFICATION")
    logger.info("=" * 80)
    logger.info("")
    
    # Get fastsimcoal2 estimate
    fsc_ncur = get_fsc_best_ncur()
    
    if not fsc_ncur:
        logger.error("Cannot proceed without fastsimcoal2 results")
        return
    
    logger.info("")
    
    # Compare with census
    expected_low, expected_high = compare_with_census()
    
    logger.info("")
    logger.info("=" * 80)
    logger.info("INTERPRETATION")
    logger.info("=" * 80)
    logger.info("")
    
    # Scenario 1: fsc_ncur is Ne (diploid)
    logger.info(f"SCENARIO 1: fastsimcoal2 outputs Ne (diploid)")
    logger.info(f"  NCUR = {fsc_ncur:,.0f} is Ne")
    logger.info(f"  Expected range: {expected_low:.0f} - {expected_high:.0f}")
    
    if expected_low <= fsc_ncur <= expected_high:
        logger.info(f"  ✅ MATCHES EXPECTED RANGE")
        logger.info(f"  → Interpretation: NCUR is Ne (diploid)")
    elif fsc_ncur < expected_low:
        logger.info(f"  ⚠️  Lower than expected")
    elif fsc_ncur > expected_high:
        logger.info(f"  ⚠️  Higher than expected")
    
    logger.info("")
    
    # Scenario 2: fsc_ncur is 2*Ne (haploid gene copies)
    logger.info(f"SCENARIO 2: fastsimcoal2 outputs 2*Ne (gene copies)")
    logger.info(f"  NCUR = {fsc_ncur:,.0f} is 2*Ne")
    logger.info(f"  → Ne (diploid) = {fsc_ncur/2:,.0f}")
    logger.info(f"  Expected range: {expected_low:.0f} - {expected_high:.0f}")
    
    if expected_low <= fsc_ncur/2 <= expected_high:
        logger.info(f"  ✅ MATCHES EXPECTED RANGE")
        logger.info(f"  → Interpretation: NCUR is 2*Ne, divide by 2 for Ne")
    elif fsc_ncur/2 < expected_low:
        logger.info(f"  ⚠️  Lower than expected")
    elif fsc_ncur/2 > expected_high:
        logger.info(f"  ⚠️  Higher than expected")
    
    logger.info("")
    logger.info("=" * 80)
    logger.info("RECOMMENDATION")
    logger.info("=" * 80)
    logger.info("")
    
    # Additional verification methods
    logger.info("")
    logger.info("ADDITIONAL VERIFICATION METHODS:")
    logger.info("-" * 80)
    logger.info("")
    logger.info("1. Literature comparison:")
    logger.info("   Search for Ne estimates in similar species")
    logger.info("   Typical endangered primates: Ne = 100-2,000")
    logger.info("")
    logger.info("2. Heterozygosity-based estimate:")
    logger.info("   Expected Ne ≈ heterozygosity / (4 × mutation_rate × generation_time)")
    logger.info("   If you have genome-wide heterozygosity, calculate this")
    logger.info("")
    logger.info("3. LD-based Ne estimation:")
    logger.info("   Use tools like NeEstimator, GONE, or LDNe")
    logger.info("   Compare with fastsimcoal2 current Ne")
    logger.info("")
    logger.info("4. ROH analysis:")
    logger.info("   Runs of homozygosity correlate with small Ne")
    logger.info("   Extensive ROH suggests Ne < 1,000")
    
    logger.info("")
    
    # Determine most likely scenario
    scenario1_match = expected_low <= fsc_ncur <= expected_high
    scenario2_match = expected_low <= fsc_ncur/2 <= expected_high
    
    if scenario1_match and not scenario2_match:
        logger.info("✅ MOST LIKELY: NCUR represents Ne (diploid)")
        logger.info("")
        logger.info("INTERPRETATION:")
        logger.info(f"  Current Ne = {fsc_ncur:,.0f} individuals")
        logger.info(f"  No conversion needed")
        logger.info("")
        logger.info("ACTION:")
        logger.info("  Use values directly from fastsimcoal2 output")
        
        return "Ne", fsc_ncur
        
    elif scenario2_match and not scenario1_match:
        logger.info("✅ MOST LIKELY: NCUR represents 2*Ne (gene copies)")
        logger.info("")
        logger.info("INTERPRETATION:")
        logger.info(f"  NCUR = {fsc_ncur:,.0f} gene copies")
        logger.info(f"  Current Ne = {fsc_ncur/2:,.0f} individuals")
        logger.info("")
        logger.info("ACTION:")
        logger.info("  Divide all Ne values by 2 in analysis scripts")
        logger.info("  Update phase3b_step5_analyze_results.py")
        
        return "2Ne", fsc_ncur/2
        
    else:
        logger.info("⚠️  AMBIGUOUS - Both scenarios possible or neither fits well")
        logger.info("")
        logger.info("RECOMMENDATION:")
        logger.info("  1. Check fastsimcoal2 documentation")
        logger.info("  2. Run test simulation with known Ne")
        logger.info("  3. Compare with independent estimates (ROH, LD)")
        logger.info("")
        logger.info("IF STILL AMBIGUOUS:")
        logger.info(f"  Assume NCUR is Ne: {fsc_ncur:,.0f}")
        logger.info(f"  Or assume NCUR is 2*Ne: {fsc_ncur/2:,.0f}")
        logger.info("  Document uncertainty in results")
        
        return "uncertain", fsc_ncur
    
    logger.info("")


def create_scaling_note(result):
    """Create note file with scaling information and recommended interpretation."""
    output_file = BASE_DIR / "output/phase3b_fastsimcoal2/NE_SCALING_NOTE.txt"
    
    scaling_type, scaled_ne = result if result is not None else ("unknown", None)

    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("Ne PARAMETER SCALING IN FASTSIMCOAL2\n")
        f.write("=" * 80 + "\n\n")

        # Machine-readable flag for downstream scripts
        f.write(f"SCALING={scaling_type}\n")
        if scaled_ne is not None:
            f.write(f"NCUR_INTERPRETED_NE={scaled_ne:.4f}\n")
        f.write("\n")
        
        f.write("ISSUE:\n")
        f.write("------\n")
        f.write("fastsimcoal2 templates use the pipeline convention that N parameters are reported as diploid Ne.\n")
        f.write("This could mean:\n")
        f.write("  - Ne (diploid effective population size), OR\n")
        f.write("  - 2*Ne (number of gene copies = haploid)\n\n")
        
        f.write("VERIFICATION NEEDED:\n")
        f.write("-------------------\n")
        f.write("Run: python3 verify_ne_scaling.py\n\n")
        
        f.write("This script will:\n")
        f.write("  1. Compare with census size (Nc ~1,200)\n")
        f.write("  2. Suggest additional verification methods\n")
        f.write("  3. Determine most likely scaling\n\n")
        
        f.write("TYPICAL Ne/Nc RATIOS:\n")
        f.write("--------------------\n")
        f.write("  Endangered species: 0.1 - 0.3\n")
        f.write("  Most mammals: 0.2 - 0.5\n")
        f.write("  Expected Ne for Nc=1,200: 120 - 600\n\n")
        
        f.write("UNTIL VERIFIED:\n")
        f.write("--------------\n")
        f.write("Report results with caveat:\n")
        f.write("  'Ne estimates subject to ×1 or ×2 scaling verification'\n\n")
    
    logger.info(f"Scaling note saved: {output_file}")


if __name__ == "__main__":
    result = analyze_scaling()
    create_scaling_note(result)
    
    logger.info("\n" + "=" * 80)
    logger.info("NEXT STEPS:")
    logger.info("=" * 80)
    logger.info("1. Review the interpretation above")
    logger.info("2. If scaling needs adjustment, update analysis scripts")
    logger.info("3. Document the scaling in your results")
    logger.info("")
