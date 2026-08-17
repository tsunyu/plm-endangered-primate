#!/usr/bin/env python3
"""
Estimate Monomorphic Sites Count for fastsimcoal2
==================================================

This script estimates the number of monomorphic ancestral sites that should
be included in the SFS d0_0 bin for more accurate likelihood calculations.

Three methods provided:
  Method A: Calculate from reference genome + callable regions
  Method B: Estimate from SNP density
  Method C: Use proportion-based estimate

Reference Genome Information:
  - 21 autosomes: 2,948,446,826 bp
  - 21 autosomes + 1 sex chromosome: 2,964,016,442 bp
  - Recommended: Use autosomal length for analysis

Author: Demographic Analysis Pipeline
Date: 2026-01-26
"""

import gzip
import numpy as np
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir
from collections import defaultdict
import logging

# Configuration
BASE_DIR = get_base_dir()
MAIN_VCF = BASE_DIR / "data/hardfilted.snp.pass.autosomes.vcf.gz"
REF_GENOME = BASE_DIR / "data/reference"  # Adjust if needed
OUTPUT_DIR = BASE_DIR / "output/phase3b_fastsimcoal2/sfs"

# Reference genome sizes (bp)
REF_GENOME_AUTOSOMES = 2_948_446_826  # 21 autosomes
REF_GENOME_WITH_SEX = 2_964_016_442   # 21 autosomes + 1 sex chromosome

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)


def detect_chromosomes_in_vcf(vcf_path):
    """Detect number of chromosomes in VCF."""
    chromosomes = set()
    
    open_func = gzip.open if str(vcf_path).endswith('.gz') else open
    
    try:
        with open_func(vcf_path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                chrom = line.split('\t')[0]
                chromosomes.add(chrom)
        
        return len(chromosomes)
    except Exception as e:
        logger.warning(f"Could not detect chromosomes: {e}")
        return None


def method_a_from_reference(ref_genome_path=None):
    """
    Method A: Calculate from reference genome length.
    
    Auto-detects chromosome count from VCF and uses appropriate genome size:
      - 21 chromosomes: 2,948,446,826 bp
      - 22 chromosomes: 2,964,016,442 bp
    """
    logger.info("Method A: Calculating from reference genome")
    
    # Auto-detect chromosome count
    num_chroms = detect_chromosomes_in_vcf(MAIN_VCF)
    
    # Select appropriate reference genome size
    if num_chroms == 22:
        total_length = REF_GENOME_WITH_SEX
        logger.info(f"  Detected: 22 chromosomes in VCF")
        logger.info(f"  Reference genome size (22 chromosomes): {total_length:,} bp")
    elif num_chroms == 21:
        total_length = REF_GENOME_AUTOSOMES
        logger.info(f"  Detected: 21 chromosomes in VCF")
        logger.info(f"  Reference genome size (21 autosomes): {total_length:,} bp")
    else:
        # Default to 22 chromosomes if detection failed or unexpected count
        total_length = REF_GENOME_WITH_SEX
        if num_chroms:
            logger.warning(f"  Unexpected chromosome count: {num_chroms}")
            logger.warning(f"  Using 22-chromosome reference: {total_length:,} bp")
        else:
            logger.warning(f"  Could not detect chromosome count")
            logger.warning(f"  Defaulting to 22-chromosome reference: {total_length:,} bp")
    
    # Estimate callable proportion (typically 80-95% of genome)
    # This accounts for repeats, low-quality regions, etc.
    callable_proportion = 0.85  # Conservative estimate
    callable_length = int(total_length * callable_proportion)
    
    logger.info(f"  Estimated callable: {callable_length:,} bp ({callable_proportion*100:.0f}%)")
    logger.info(f"  Non-callable (repeats, low quality): {total_length - callable_length:,} bp")
    
    return callable_length


def method_b_from_snp_density(vcf_path):
    """
    Method B: Estimate from SNP density.
    
    Logic: If we have X variants, and genome is Y bp, 
           monomorphic sites ≈ Y - X
    """
    logger.info("Method B: Estimating from SNP density")
    
    # Count variants per chromosome
    chrom_variants = defaultdict(int)
    chrom_span = {}  # Track genomic span
    
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            
            chrom_variants[chrom] += 1
            
            if chrom not in chrom_span:
                chrom_span[chrom] = [pos, pos]
            else:
                chrom_span[chrom][1] = pos
    
    # Calculate total span and variants
    total_variants = sum(chrom_variants.values())
    total_span = sum(end - start + 1 for start, end in chrom_span.values())
    
    logger.info(f"  Total variants: {total_variants:,}")
    logger.info(f"  Total span: {total_span:,} bp")
    
    # Estimate monomorphic sites
    # Assume variants cover ~90% of callable genome
    estimated_callable = total_span / 0.9
    monomorphic = int(estimated_callable - total_variants)
    
    logger.info(f"  Estimated callable genome: {estimated_callable:,.0f} bp")
    logger.info(f"  Estimated monomorphic sites: {monomorphic:,}")
    
    return monomorphic, total_variants


def method_c_proportion_based(vcf_path):
    """
    Method C: Use typical SNP/genome ratio.
    
    For mammals, typical heterozygosity is 0.001-0.005
    So monomorphic sites ≈ variants / heterozygosity
    """
    logger.info("Method C: Proportion-based estimate")
    
    # Count variants
    total_variants = 0
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            total_variants += 1
    
    logger.info(f"  Total variants: {total_variants:,}")
    
    # Use typical mammalian heterozygosity
    # For endangered species, often lower: 0.0005-0.002
    heterozygosity_low = 0.0005
    heterozygosity_high = 0.002
    
    # Estimate genome length
    genome_low = total_variants / heterozygosity_high
    genome_high = total_variants / heterozygosity_low
    
    monomorphic_low = int(genome_low - total_variants)
    monomorphic_high = int(genome_high - total_variants)
    
    logger.info(f"  Estimated genome length: {genome_low:,.0f} - {genome_high:,.0f} bp")
    logger.info(f"  Estimated heterozygosity: {heterozygosity_low} - {heterozygosity_high}")
    logger.info(f"  Estimated monomorphic: {monomorphic_low:,} - {monomorphic_high:,}")
    
    # Use midpoint
    monomorphic = int((monomorphic_low + monomorphic_high) / 2)
    
    return monomorphic, total_variants


def update_sfs_file(sfs_file, monomorphic_count):
    """Update existing SFS file with monomorphic count."""
    logger.info(f"Updating SFS file: {sfs_file}")
    
    # Read existing SFS
    with open(sfs_file, 'r') as f:
        lines = f.readlines()
    
    if len(lines) < 3:
        logger.error("Invalid SFS file format")
        return False
    
    # Parse data line (strip() removes leading tab, so data[0] is d0_0)
    data = lines[2].strip().split('\t')
    
    # Update d0_0 (first data column)
    data[0] = str(monomorphic_count)
    
    # Write updated file
    backup_file = sfs_file.with_suffix('.obs.backup')
    import shutil
    shutil.copy(sfs_file, backup_file)
    
    with open(sfs_file, 'w') as f:
        f.write(lines[0])  # observations line
        f.write(lines[1])  # header line
        f.write('\t' + '\t'.join(data) + '\n')  # updated data
    
    logger.info(f"  Backup saved: {backup_file}")
    logger.info(f"  Updated d0_0 to: {monomorphic_count:,}")
    
    return True


def main():
    """Main workflow."""
    logger.info("=" * 80)
    logger.info("ESTIMATING MONOMORPHIC SITES COUNT")
    logger.info("=" * 80)
    logger.info("")
    logger.info("REFERENCE GENOME INFORMATION:")
    logger.info(f"  21 chromosomes: {REF_GENOME_AUTOSOMES:,} bp")
    logger.info(f"  22 chromosomes: {REF_GENOME_WITH_SEX:,} bp")
    logger.info(f"  (Auto-detecting from VCF)")
    logger.info("")
    
    results = {}
    
    # Try Method A (using known reference genome size)
    callable_length = method_a_from_reference()
    if callable_length:
        results['method_a'] = callable_length
    
    logger.info("")
    
    # Method B (SNP density)
    try:
        monomorphic_b, variants_b = method_b_from_snp_density(MAIN_VCF)
        results['method_b'] = {'monomorphic': monomorphic_b, 'variants': variants_b}
    except Exception as e:
        logger.error(f"Method B failed: {e}")
    
    logger.info("")
    
    # Method C (proportion-based)
    try:
        monomorphic_c, variants_c = method_c_proportion_based(MAIN_VCF)
        results['method_c'] = {'monomorphic': monomorphic_c, 'variants': variants_c}
    except Exception as e:
        logger.error(f"Method C failed: {e}")
    
    logger.info("")
    logger.info("=" * 80)
    logger.info("SUMMARY OF ESTIMATES")
    logger.info("=" * 80)
    logger.info("")
    
    estimates = []
    
    if 'method_a' in results:
        logger.info(f"Method A (reference): {results['method_a']:,} callable bp")
        if 'method_b' in results:
            mono_a = results['method_a'] - results['method_b']['variants']
            estimates.append(mono_a)
            logger.info(f"         → Monomorphic: {mono_a:,}")
    
    if 'method_b' in results:
        mono_b = results['method_b']['monomorphic']
        estimates.append(mono_b)
        logger.info(f"Method B (SNP density): {mono_b:,} monomorphic sites")
    
    if 'method_c' in results:
        mono_c = results['method_c']['monomorphic']
        estimates.append(mono_c)
        logger.info(f"Method C (proportion): {mono_c:,} monomorphic sites")
    
    logger.info("")
    
    if estimates:
        recommended = int(np.median(estimates))
        logger.info(f"RECOMMENDED VALUE (median): {recommended:,}")
        logger.info("")
        
        # Ask to update SFS file
        sfs_file = OUTPUT_DIR / "SNJ_DAFpop0.obs"
        
        if sfs_file.exists():
            logger.info(f"SFS file found: {sfs_file}")
            logger.info("")
            logger.info(f"Automatically updating SFS with recommended monomorphic count: {recommended:,}")
            
            # Save recommendation to file
            rec_file = OUTPUT_DIR / "monomorphic_recommendation.txt"
            with open(rec_file, 'w') as f:
                f.write(f"MONOMORPHIC SITES ESTIMATE\n")
                f.write(f"=" * 60 + "\n\n")
                f.write(f"Recommended value: {recommended:,}\n\n")
                f.write(f"Individual estimates:\n")
                for i, est in enumerate(estimates, 1):
                    f.write(f"  Method {i}: {est:,}\n")
                f.write(f"\nSFS file automatically updated with this value.\n")
                f.write(f"Original SFS backed up as SNJ_DAFpop0.obs.backup\n")
            
            logger.info(f"Recommendation saved to: {rec_file}")
            
            # Automatically update SFS in-place (with backup handled by update_sfs_file)
            updated = update_sfs_file(sfs_file, recommended)
            if updated:
                logger.info("SFS update completed successfully.")
            else:
                logger.warning("SFS update reported failure; please check the SFS file format.")
        else:
            logger.warning(f"SFS file not found. Generate it first with phase3b_step0_prepare_sfs.py")
    else:
        logger.error("No valid estimates obtained")
    
    logger.info("")
    logger.info("=" * 80)
    logger.info("ANALYSIS COMPLETE")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
