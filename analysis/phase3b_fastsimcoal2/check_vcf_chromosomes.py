#!/usr/bin/env python3
"""
Check VCF Chromosomes
=====================

Quick utility to check which chromosomes are present in the VCF file
and determine appropriate reference genome size.

Key Features:
  - Detects number of chromosomes in VCF
  - Identifies sex chromosomes by variant count analysis
  - If 22 chromosomes detected, chromosome 22 is identified as sex chromosome
  - Provides actionable recommendations for filtering

Usage:
    python3 check_vcf_chromosomes.py                    # Check default VCF
    python3 check_vcf_chromosomes.py path/to/vcf.gz     # Check specific VCF
"""

import gzip
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir
from collections import Counter

# Configuration
BASE_DIR = get_base_dir()
MAIN_VCF = BASE_DIR / "data/hardfilted.snp.pass.autosomes.vcf.gz"

# Reference genome sizes
REF_21_AUTOSOMES = 2_948_446_826
REF_22_CHROMOSOMES = 2_964_016_442


def check_chromosomes(vcf_path):
    """Check which chromosomes are in the VCF."""
    print("=" * 80)
    print("VCF CHROMOSOME CHECK")
    print("=" * 80)
    print()
    print(f"Checking: {vcf_path}")
    print()
    
    chromosomes = set()
    chr_counts = Counter()
    
    open_func = gzip.open if str(vcf_path).endswith('.gz') else open
    
    with open_func(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            chrom = line.split('\t')[0]
            chromosomes.add(chrom)
            chr_counts[chrom] += 1
    
    # Sort chromosomes
    autosome_pattern = set()
    sex_chromosomes = set()
    
    for chrom in chromosomes:
        if chrom.upper() in ['X', 'Y', 'CHR_X', 'CHR_Y', 'CHRX', 'CHRY']:
            sex_chromosomes.add(chrom)
        elif chrom.isdigit() or chrom.startswith('chr'):
            autosome_pattern.add(chrom)
        else:
            print(f"  Unknown chromosome format: {chrom}")
    
    # Convert to numbers for sorting
    autosome_nums = []
    for chrom in autosome_pattern:
        try:
            if chrom.isdigit():
                autosome_nums.append(int(chrom))
            elif chrom.lower().startswith('chr'):
                autosome_nums.append(int(chrom[3:]))
        except:
            pass
    
    autosome_nums.sort()
    
    print("RESULTS:")
    print("-" * 80)
    print(f"Total chromosomes found: {len(chromosomes)}")
    print()
    
    if autosome_nums:
        print(f"Autosomes: {len(autosome_nums)}")
        print(f"  Range: {min(autosome_nums)} - {max(autosome_nums)}")
        print(f"  List: {', '.join(map(str, autosome_nums))}")
        print()
    
    if sex_chromosomes:
        print(f"Sex chromosomes: {len(sex_chromosomes)}")
        print(f"  Found: {', '.join(sorted(sex_chromosomes))}")
        print()
    
    # Variant counts
    print("Variant counts per chromosome:")
    for chrom in sorted(chromosomes, key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else x)):
        print(f"  {chrom:>3}: {chr_counts[chrom]:>12,}")
    
    print()
    print("=" * 80)
    print("CHROMOSOME INTERPRETATION")
    print("=" * 80)
    print()
    
    total_chroms = len(chromosomes)
    
    # Check if chr 22 has significantly fewer variants (indicating sex chromosome)
    chr_22_is_sex = False
    if total_chroms == 22 and '22' in chr_counts:
        chr_22_count = chr_counts['22']
        # Get average count of chr 1-21
        chr_1_21_counts = [chr_counts[str(i)] for i in range(1, 22) if str(i) in chr_counts]
        if chr_1_21_counts:
            avg_autosome_count = sum(chr_1_21_counts) / len(chr_1_21_counts)
            # If chr 22 has less than 50% of average, it's likely a sex chromosome
            if chr_22_count < 0.5 * avg_autosome_count:
                chr_22_is_sex = True
                print(f"📊 Variant Analysis:")
                print(f"   Chromosomes 1-21 average: {avg_autosome_count:,.0f} variants")
                print(f"   Chromosome 22: {chr_22_count:,} variants ({100*chr_22_count/avg_autosome_count:.1f}% of average)")
                print(f"   → Chromosome 22 is likely a SEX CHROMOSOME")
                print()
    
    print("=" * 80)
    print("REFERENCE GENOME RECOMMENDATION")
    print("=" * 80)
    print()
    
    if total_chroms == 21:
        print(f"✅ Your VCF has 21 chromosomes (autosomes only)")
        print(f"   Recommended reference genome: {REF_21_AUTOSOMES:,} bp")
        print(f"   Status: CORRECT for fastsimcoal2 demographic inference")
        recommended_size = REF_21_AUTOSOMES
        
    elif total_chroms == 22:
        if chr_22_is_sex:
            print(f"⚠️  Your VCF has 22 chromosomes")
            print(f"   Chromosomes 1-21: AUTOSOMES")
            print(f"   Chromosome 22: SEX CHROMOSOME (detected by low variant count)")
            print()
            print(f"   For demographic inference:")
            print(f"   ❌ Current VCF includes sex chromosome (NOT recommended)")
            print(f"   ✅ Should filter to autosomes only (chromosomes 1-21)")
            print()
            print(f"   Reference genome if using 22 chr: {REF_22_CHROMOSOMES:,} bp")
            print(f"   Reference genome for autosomes only: {REF_21_AUTOSOMES:,} bp")
            print()
            print(f"   ACTION REQUIRED: Run filter_autosomes_only.sh")
            recommended_size = REF_22_CHROMOSOMES  # Current state
        elif sex_chromosomes:
            print(f"✅ Your VCF has 22 chromosomes")
            print(f"   Includes sex chromosome(s): {', '.join(sex_chromosomes)}")
            print(f"   Recommended reference genome: {REF_22_CHROMOSOMES:,} bp")
            print()
            print(f"   ACTION REQUIRED: Filter to autosomes only for demographic inference")
            recommended_size = REF_22_CHROMOSOMES
        else:
            print(f"⚠️  Your VCF has 22 chromosomes (no clear sex chromosome detected)")
            print(f"   Recommended reference genome: {REF_22_CHROMOSOMES:,} bp")
            print()
            print(f"   Note: Standard practice is to use autosomes only")
            recommended_size = REF_22_CHROMOSOMES
            
    else:
        print(f"⚠️  Unexpected chromosome count: {total_chroms}")
        if total_chroms < 21:
            print(f"   Using 21-autosome reference: {REF_21_AUTOSOMES:,} bp")
            recommended_size = REF_21_AUTOSOMES
        else:
            print(f"   Using 22-chromosome reference: {REF_22_CHROMOSOMES:,} bp")
            recommended_size = REF_22_CHROMOSOMES
    
    print()
    
    if total_chroms == 22 and chr_22_is_sex:
        print("=" * 80)
        print("⚠️  ACTION REQUIRED: FILTER TO AUTOSOMES ONLY")
        print("=" * 80)
        print()
        print("Standard practice for fastsimcoal2: Use AUTOSOMES ONLY")
        print()
        print("Your VCF includes chromosome 22 (sex chromosome).")
        print("For accurate demographic inference, you should:")
        print()
        print("1. Filter VCF to autosomes only (chromosomes 1-21):")
        print("   bash filter_autosomes_only.sh")
        print()
        print("2. Update all scripts to use the filtered VCF:")
        print("   bash update_to_autosomes_only.sh")
        print()
        print("3. Verify the change:")
        print("   python3 check_vcf_chromosomes.py")
        print()
        print("Why exclude sex chromosomes?")
        print("  - Different effective population size (X = 3/4 autosomal Ne)")
        print("  - Different inheritance patterns (hemizygous in males)")
        print("  - fastsimcoal2 assumes autosomal diploid inheritance")
        print("  - Standard practice in published studies")
        print()
        print("See AUTOSOMES_GUIDE.txt for detailed explanation.")
        print()
    elif total_chroms == 21:
        print("=" * 80)
        print("✅ CORRECT CONFIGURATION")
        print("=" * 80)
        print()
        print("Your VCF contains autosomes only (21 chromosomes).")
        print("This is the correct configuration for demographic inference.")
        print()
        print("Reference genome size:")
        print(f"  {REF_21_AUTOSOMES:,} bp (21 autosomes)")
        print()
        print("You can proceed with the analysis:")
        print("  python3 phase3b_step0_prepare_sfs.py")
        print("  bash run_complete_pipeline.sh")
        print()
    
    return {
        'total_chromosomes': total_chroms,
        'autosomes': autosome_nums,
        'sex_chromosomes': list(sex_chromosomes),
        'chr_22_is_sex': chr_22_is_sex,
        'recommended_size': recommended_size,
        'variant_counts': dict(chr_counts)
    }


def main():
    # Allow command-line argument to specify VCF path
    if len(sys.argv) > 1:
        vcf_path = Path(sys.argv[1])
        if not vcf_path.is_absolute():
            vcf_path = Path.cwd() / vcf_path
    else:
        vcf_path = MAIN_VCF
    
    if not vcf_path.exists():
        print(f"ERROR: VCF file not found: {vcf_path}")
        print()
        print("Usage:")
        print(f"  python3 {Path(__file__).name}")
        print(f"  python3 {Path(__file__).name} <path/to/vcf.gz>")
        print()
        print(f"Default VCF: {MAIN_VCF}")
        return
    
    result = check_chromosomes(vcf_path)
    
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print()
    print(f"Total variants: {sum(result['variant_counts'].values()):,}")
    print(f"Chromosomes: {result['total_chromosomes']}")
    
    if result['chr_22_is_sex']:
        print(f"Chromosome 22: SEX CHROMOSOME (detected)")
        print()
        print(f"⚠️  STATUS: NEEDS FILTERING")
        print(f"   Current VCF includes sex chromosome")
        print(f"   Action: Run filter_autosomes_only.sh")
        print()
        print(f"After filtering to 21 autosomes:")
        # Calculate for 21 autosomes
        ref_21 = REF_21_AUTOSOMES
        callable_21 = int(ref_21 * 0.85)
        # Approximate variants for chr 1-21 (exclude chr 22)
        chr_22_variants = result['variant_counts'].get('22', 0)
        variants_21 = sum(result['variant_counts'].values()) - chr_22_variants
        monomorphic_21 = callable_21 - variants_21
        print(f"  Reference genome: {ref_21:,} bp (21 autosomes)")
        print(f"  Callable (85%): {callable_21:,} bp")
        print(f"  Variants (chr 1-21): ~{variants_21:,}")
        print(f"  Monomorphic sites: ~{monomorphic_21:,}")
    else:
        print(f"Recommended genome size: {result['recommended_size']:,} bp")
        print()
        
        if result['total_chromosomes'] == 21:
            print("✅ STATUS: CORRECT (autosomes only)")
        else:
            print("⚠️  STATUS: Review chromosome composition")
        
        print()
        
        # Calculate callable genome
        callable_prop = 0.85
        callable = int(result['recommended_size'] * callable_prop)
        total_variants = sum(result['variant_counts'].values())
        monomorphic = callable - total_variants
        
        print("ESTIMATED MONOMORPHIC SITES:")
        print(f"  Reference genome: {result['recommended_size']:,} bp")
        print(f"  Callable (85%): {callable:,} bp")
        print(f"  Total variants: {total_variants:,}")
        print(f"  Monomorphic sites: {monomorphic:,}")
    
    print()


if __name__ == "__main__":
    main()
