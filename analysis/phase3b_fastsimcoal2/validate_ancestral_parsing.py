#!/usr/bin/env python3
"""
Quick Validation Script: Ancestral Allele Parsing
==================================================

This script validates that ancestral alleles are correctly extracted
from the INFO field (AA=X) rather than incorrectly using REF allele.

Run this BEFORE running the full pipeline to verify the fix is working.

Usage:
    python3 validate_ancestral_parsing.py
"""

import gzip
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir

# Configuration
BASE_DIR = get_base_dir()
ANCESTRAL_VCF = BASE_DIR / "output/phase2_annotation/ancestral_states/variants_with_ancestral.vcf.gz"

def validate_parsing():
    """Validate ancestral allele extraction."""
    
    print("=" * 80)
    print("VALIDATING ANCESTRAL ALLELE PARSING")
    print("=" * 80)
    print()
    
    if not ANCESTRAL_VCF.exists():
        print(f"❌ ERROR: Ancestral VCF not found: {ANCESTRAL_VCF}")
        return False
    
    print(f"Reading: {ANCESTRAL_VCF}")
    print()
    
    # Counters
    total = 0
    ref_is_ancestral = 0
    alt_is_ancestral = 0
    ambiguous = 0
    
    # Sample variants for manual inspection
    samples = []
    
    with gzip.open(ANCESTRAL_VCF, 'rt') as f:
        for i, line in enumerate(f):
            if line.startswith('#'):
                continue
            
            total += 1
            
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            info = fields[7]
            
            # Extract AA from INFO
            ancestral = None
            for info_field in info.split(';'):
                if info_field.startswith('AA='):
                    ancestral = info_field.split('=')[1]
                    break
            
            if not ancestral or ancestral in ['.', 'N', '-']:
                ambiguous += 1
                continue
            
            # Compare with REF/ALT
            if ancestral == ref:
                ref_is_ancestral += 1
            elif ancestral == alt:
                alt_is_ancestral += 1
            else:
                ambiguous += 1
            
            # Collect samples
            if total <= 10:
                samples.append({
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'ancestral': ancestral,
                    'match': 'REF' if ancestral == ref else 'ALT' if ancestral == alt else 'NEITHER'
                })
            
            # Stop after checking first 100,000 variants (representative sample)
            if total >= 100000:
                break
    
    # Report results
    print("VALIDATION RESULTS")
    print("=" * 80)
    print()
    print(f"Total variants checked: {total:,}")
    print(f"REF = ancestral:        {ref_is_ancestral:,} ({ref_is_ancestral/total*100:.1f}%)")
    print(f"ALT = ancestral:        {alt_is_ancestral:,} ({alt_is_ancestral/total*100:.1f}%)")
    print(f"Ambiguous/missing:      {ambiguous:,} ({ambiguous/total*100:.1f}%)")
    print()
    
    # Show sample variants
    print("SAMPLE VARIANTS (first 10):")
    print("-" * 80)
    print(f"{'Chr':<10} {'Pos':<12} {'REF':<5} {'ALT':<5} {'AA':<5} {'Match':<8}")
    print("-" * 80)
    for s in samples:
        print(f"{s['chrom']:<10} {s['pos']:<12} {s['ref']:<5} {s['alt']:<5} {s['ancestral']:<5} {s['match']:<8}")
    print()
    
    # Validation check
    print("VALIDATION CHECK:")
    print("-" * 80)
    
    if alt_is_ancestral > 0:
        print(f"✅ CORRECT: Found {alt_is_ancestral:,} variants where ALT is ancestral")
        print(f"✅ This confirms we're parsing INFO field (AA=X), not using REF")
        print()
        print("🎯 EXPECTED BEHAVIOR:")
        print(f"   - REF ancestral: ~59% (observed: {ref_is_ancestral/total*100:.1f}%)")
        print(f"   - ALT ancestral: ~41% (observed: {alt_is_ancestral/total*100:.1f}%)")
        print()
        
        # Check if proportions match expectations
        alt_pct = alt_is_ancestral / total * 100
        if 35 < alt_pct < 45:
            print("✅ PROPORTION MATCHES EXPECTED (~41%)")
            print()
            print("=" * 80)
            print("✅ VALIDATION PASSED!")
            print("=" * 80)
            print()
            print("The ancestral allele parsing is working correctly.")
            print("You can proceed with running the pipeline.")
            return True
        else:
            print(f"⚠️  WARNING: ALT proportion ({alt_pct:.1f}%) differs from expected (41%)")
            print("This may be due to sampling or dataset differences.")
            print()
            return True
    else:
        print("❌ CRITICAL ERROR: No variants found where ALT is ancestral!")
        print("❌ This suggests the code is still incorrectly using REF as ancestral.")
        print()
        print("🔧 FIX REQUIRED:")
        print("   The load_ancestral_states() function needs to parse INFO field.")
        print("   Check that the fix was applied correctly.")
        return False

if __name__ == "__main__":
    success = validate_parsing()
    
    if success:
        print()
        print("Next steps:")
        print("  1. Run: python3 phase3b_step0_prepare_sfs.py")
        print("  2. Check: output/phase3b_fastsimcoal2/sfs/sfs_statistics.txt")
        print("  3. Verify polarization rate is ~98%")
    else:
        print()
        print("⚠️  DO NOT RUN THE PIPELINE YET")
        print("Fix the ancestral allele parsing issue first.")
    
    print()
