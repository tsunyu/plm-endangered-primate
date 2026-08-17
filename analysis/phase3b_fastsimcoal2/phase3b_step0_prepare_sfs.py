#!/usr/bin/env python3
"""
Prepare Site Frequency Spectrum (SFS) for fastsimcoal2
=======================================================

This script generates an unfolded (derived) allele frequency spectrum from VCF data
using ancestral state information to polarize variants.

Input:
  - Main VCF: data/monkey_snp_autosomes_only.vcf.gz (21 autosomes)
  - Ancestral VCF: output/phase2_annotation/ancestral_states/variants_with_ancestral.vcf.gz
    (REF allele = ancestral allele)

Output:
  - SNJ_DAFpop0.obs: Unfolded SFS for fastsimcoal2
  - sfs_statistics.txt: Summary statistics
  - sfs_plot.png: Visualization

Author: Demographic Analysis Pipeline
Date: 2026-01-26
"""

import sys
import gzip
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir
from collections import defaultdict
import logging

# Configuration
BASE_DIR = get_base_dir()
MAIN_VCF = BASE_DIR / "data/hardfilted.snp.pass.autosomes.vcf.gz"
ANCESTRAL_VCF = BASE_DIR / "output/phase2_annotation/ancestral_states/variants_with_ancestral.vcf.gz"
OUTPUT_DIR = BASE_DIR / "output/phase3b_fastsimcoal2/sfs"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(OUTPUT_DIR / "sfs_generation.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def parse_vcf_line(line):
    """Parse a VCF data line and return key information."""
    fields = line.strip().split('\t')
    chrom = fields[0]
    pos = int(fields[1])
    ref = fields[3]
    alt = fields[4]
    genotypes = fields[9:]  # Sample genotypes
    return chrom, pos, ref, alt, genotypes


def count_alleles(genotypes, target_allele_idx):
    """
    Count occurrences of target allele and total called alleles in genotype fields.
    
    Args:
        genotypes: List of genotype strings (e.g., '0/0', '0/1', '1/1')
        target_allele_idx: 0 for REF, 1 for ALT
    
    Returns:
        (target_count, total_called): Count of target allele and total non-missing alleles
    """
    target_count = 0
    total_called = 0
    for gt in genotypes:
        # Extract genotype (GT field is first, before ':' if present)
        gt_field = gt.split(':')[0]
        
        # Handle different genotype formats
        if '/' in gt_field:
            alleles = gt_field.split('/')
        elif '|' in gt_field:
            alleles = gt_field.split('|')
        else:
            continue  # Skip malformed genotypes
        
        # Count target allele and total called alleles
        for allele in alleles:
            if allele == '.':
                continue  # Missing data
            try:
                total_called += 1
                if int(allele) == target_allele_idx:
                    target_count += 1
            except ValueError:
                continue
    
    return target_count, total_called


def load_ancestral_states(ancestral_vcf_path):
    """
    Load ancestral states from ancestral VCF.
    
    The ancestral allele is stored in the INFO field as AA=X
    
    Returns:
        Dictionary: {(chrom, pos): ancestral_allele}
    """
    logger.info(f"Loading ancestral states from {ancestral_vcf_path}")
    ancestral_dict = {}
    
    open_func = gzip.open if str(ancestral_vcf_path).endswith('.gz') else open
    
    with open_func(ancestral_vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            info = fields[7]  # INFO field
            
            # Extract ancestral allele from INFO field (AA=X)
            ancestral_allele = None
            for info_field in info.split(';'):
                if info_field.startswith('AA='):
                    ancestral_allele = info_field.split('=')[1]
                    break
            
            if ancestral_allele and ancestral_allele not in ['.', 'N', '-']:
                ancestral_dict[(chrom, pos)] = ancestral_allele
    
    logger.info(f"Loaded {len(ancestral_dict):,} ancestral states")
    return ancestral_dict


def generate_sfs(main_vcf_path, ancestral_dict, n_samples):
    """
    Generate unfolded (derived) allele frequency spectrum.
    
    Args:
        main_vcf_path: Path to main VCF file
        ancestral_dict: Dictionary of ancestral states
        n_samples: Number of diploid individuals
    
    Returns:
        sfs: Array of SFS counts (0 to 2*n_samples)
        stats: Dictionary of statistics
    """
    n_chromosomes = 2 * n_samples
    sfs = np.zeros(n_chromosomes + 1, dtype=int)
    
    stats = {
        'total_variants': 0,
        'biallelic': 0,
        'polarized': 0,
        'unpolarized': 0,
        'no_ancestral': 0,
        'ambiguous': 0,
        'monomorphic': 0,
        'polymorphic': 0
    }
    
    logger.info(f"Processing main VCF: {main_vcf_path}")
    logger.info(f"Sample size: {n_samples} diploids ({n_chromosomes} chromosomes)")
    
    open_func = gzip.open if str(main_vcf_path).endswith('.gz') else open
    
    with open_func(main_vcf_path, 'rt') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue
            
            if line_num % 10000 == 0:
                logger.info(f"Processed {line_num:,} lines, polarized {stats['polarized']:,} variants")
            
            stats['total_variants'] += 1
            
            chrom, pos, ref, alt, genotypes = parse_vcf_line(line)
            
            # Skip multi-allelic sites
            if ',' in alt:
                continue
            
            stats['biallelic'] += 1
            
            # Look up ancestral allele
            ancestral_allele = ancestral_dict.get((chrom, pos))
            
            if ancestral_allele is None:
                stats['no_ancestral'] += 1
                continue
            
            # Determine derived allele and count it
            derived_count = None
            
            if ancestral_allele == ref:
                # ALT is derived; count ALT alleles
                alt_count, total_called = count_alleles(genotypes, 1)
                derived_count = alt_count
                stats['polarized'] += 1
                
            elif ancestral_allele == alt:
                # REF is derived; derived_count = total_called - ALT_count
                alt_count, total_called = count_alleles(genotypes, 1)
                derived_count = total_called - alt_count
                stats['polarized'] += 1
                
            else:
                # Ancestral allele doesn't match REF or ALT
                stats['ambiguous'] += 1
                continue
            
            # Skip sites with too much missing data (require >=50% called)
            if total_called < n_chromosomes * 0.5:
                stats.setdefault('low_call_rate', 0)
                stats['low_call_rate'] += 1
                continue
            
            # Add to SFS (use total_called as effective sample size)
            # For fastsimcoal2, SFS must have fixed dimension (2*n_samples + 1),
            # so we project to the full sample size only when all sites are called.
            # Sites with missing data are skipped to avoid bias.
            if total_called != n_chromosomes:
                stats.setdefault('partial_missing', 0)
                stats['partial_missing'] += 1
                continue  # Only use fully called sites for unbiased SFS
            
            if 0 <= derived_count <= n_chromosomes:
                sfs[derived_count] += 1
                
                if derived_count == 0 or derived_count == n_chromosomes:
                    stats['monomorphic'] += 1
                else:
                    stats['polymorphic'] += 1
    
    logger.info(f"SFS generation complete. Polarized {stats['polarized']:,} variants")
    return sfs, stats


def write_fastsimcoal2_sfs(sfs, output_path, pop_name="SNJ"):
    """
    Write SFS in fastsimcoal2 format.
    
    Format for unfolded SFS:
    1 observations
            d0_0    d0_1    d0_2    ...    d0_n
            count0  count1  count2  ...    countn
    
    Note: For fastsimcoal2, monomorphic sites (d0_0) should ideally represent
    the total number of monomorphic ancestral sites in the genome.
    Since we only have variant sites, d0_0 will be 0 unless we estimate
    the total callable genome length.
    """
    logger.info(f"Writing fastsimcoal2 SFS to {output_path}")
    
    with open(output_path, 'w') as f:
        f.write("1 observations\n")
        
        # Header line with column names
        n_bins = len(sfs)
        header = '\t'.join([f"d0_{i}" for i in range(n_bins)])
        f.write('\t' + header + '\n')
        
        # Data line with counts
        data = '\t'.join([str(int(count)) for count in sfs])
        f.write('\t' + data + '\n')
    
    logger.info(f"SFS file written: {output_path}")
    logger.warning("Note: Monomorphic sites (d0_0) are set to 0 since only variant sites are available")
    logger.warning("For more accurate likelihood, consider estimating total callable genome length")


def write_statistics(sfs, stats, output_path):
    """Write detailed SFS statistics."""
    logger.info(f"Writing statistics to {output_path}")
    
    n_chromosomes = len(sfs) - 1
    
    with open(output_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("SITE FREQUENCY SPECTRUM (SFS) GENERATION SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("Input Files:\n")
        f.write(f"  Main VCF: {MAIN_VCF}\n")
        f.write(f"  Ancestral VCF: {ANCESTRAL_VCF}\n\n")
        
        f.write("Sample Information:\n")
        f.write(f"  Diploid individuals: {n_chromosomes // 2}\n")
        f.write(f"  Total chromosomes: {n_chromosomes}\n\n")
        
        f.write("-" * 80 + "\n")
        f.write("VARIANT PROCESSING STATISTICS\n")
        f.write("-" * 80 + "\n\n")
        
        f.write(f"Total variants in main VCF: {stats['total_variants']:,}\n")
        f.write(f"Biallelic variants: {stats['biallelic']:,} ({stats['biallelic']/max(stats['total_variants'],1)*100:.1f}%)\n")
        f.write(f"Successfully polarized: {stats['polarized']:,} ({stats['polarized']/max(stats['biallelic'],1)*100:.1f}%)\n")
        f.write(f"No ancestral state: {stats['no_ancestral']:,}\n")
        f.write(f"Ambiguous polarization: {stats['ambiguous']:,}\n\n")
        
        f.write("-" * 80 + "\n")
        f.write("DERIVED ALLELE FREQUENCY SPECTRUM\n")
        f.write("-" * 80 + "\n\n")
        
        # Frequency categories
        singletons = sfs[1]
        doubletons = sfs[2]
        low_freq = np.sum(sfs[1:10])
        intermediate = np.sum(sfs[10:n_chromosomes-10])
        high_freq = np.sum(sfs[n_chromosomes-10:n_chromosomes])
        fixed_derived = sfs[n_chromosomes]
        
        f.write(f"Monomorphic ancestral (DAF=0): {sfs[0]:,}\n")
        f.write(f"Singletons (DAF=1): {singletons:,}\n")
        f.write(f"Doubletons (DAF=2): {doubletons:,}\n")
        f.write(f"Low frequency (1≤DAF<10): {low_freq:,}\n")
        f.write(f"Intermediate (10≤DAF<{n_chromosomes-10}): {intermediate:,}\n")
        f.write(f"High frequency ({n_chromosomes-10}≤DAF<{n_chromosomes}): {high_freq:,}\n")
        f.write(f"Fixed derived (DAF={n_chromosomes}): {fixed_derived:,}\n\n")
        
        f.write(f"Total polymorphic sites: {stats['polymorphic']:,}\n")
        f.write(f"Total monomorphic sites: {stats['monomorphic']:,}\n\n")
        
        # Summary metrics
        if stats['polymorphic'] > 0:
            f.write("-" * 80 + "\n")
            f.write("SFS SUMMARY METRICS\n")
            f.write("-" * 80 + "\n\n")
            
            # Tajima's D-related metrics
            pi = np.sum([2 * i * (n_chromosomes - i) / (n_chromosomes * (n_chromosomes - 1)) * sfs[i] 
                         for i in range(1, n_chromosomes)])
            
            f.write(f"Nucleotide diversity (π): {pi:.6f}\n")
            f.write(f"Watterson's θ estimate: {stats['polymorphic'] / sum(1/i for i in range(1, n_chromosomes)):.2f}\n")
            f.write(f"Singleton/Doubleton ratio: {singletons / max(doubletons, 1):.2f}\n")
            
            # Proportion in each category
            f.write(f"\nFrequency distribution:\n")
            f.write(f"  Rare variants (DAF<5%): {np.sum(sfs[1:max(int(n_chromosomes*0.05),1)]):,} ")
            f.write(f"({np.sum(sfs[1:max(int(n_chromosomes*0.05),1)])/stats['polymorphic']*100:.1f}%)\n")
            f.write(f"  Intermediate (5-95%): {np.sum(sfs[max(int(n_chromosomes*0.05),1):int(n_chromosomes*0.95)]):,} ")
            f.write(f"({np.sum(sfs[max(int(n_chromosomes*0.05),1):int(n_chromosomes*0.95)])/stats['polymorphic']*100:.1f}%)\n")
            f.write(f"  High frequency (>95%): {np.sum(sfs[int(n_chromosomes*0.95):n_chromosomes]):,} ")
            f.write(f"({np.sum(sfs[int(n_chromosomes*0.95):n_chromosomes])/stats['polymorphic']*100:.1f}%)\n")
        
        f.write("\n" + "=" * 80 + "\n")
        f.write("OUTPUT FILES\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"fastsimcoal2 SFS file: {OUTPUT_DIR / 'SNJ_DAFpop0.obs'}\n")
        f.write(f"Statistics: {output_path}\n")
        f.write(f"Plot: {OUTPUT_DIR / 'sfs_plot.png'}\n\n")


def plot_sfs(sfs, output_path):
    """Create visualization of the SFS."""
    logger.info(f"Creating SFS plot: {output_path}")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Plot 1: Full SFS (log scale)
    ax = axes[0, 0]
    n_bins = len(sfs)
    x = np.arange(n_bins)
    ax.bar(x, sfs, color='steelblue', alpha=0.7, edgecolor='black')
    ax.set_xlabel('Derived Allele Count', fontweight='bold')
    ax.set_ylabel('Number of SNPs', fontweight='bold')
    ax.set_title('Unfolded Site Frequency Spectrum', fontweight='bold')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Plot 2: SFS excluding monomorphic and fixed
    ax = axes[0, 1]
    x_poly = np.arange(1, n_bins - 1)
    sfs_poly = sfs[1:-1]
    ax.bar(x_poly, sfs_poly, color='darkgreen', alpha=0.7, edgecolor='black')
    ax.set_xlabel('Derived Allele Count', fontweight='bold')
    ax.set_ylabel('Number of SNPs', fontweight='bold')
    ax.set_title('Polymorphic Sites Only (excluding 0 and n)', fontweight='bold')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Plot 3: Low frequency variants (DAF 1-20)
    ax = axes[1, 0]
    low_freq_range = min(20, n_bins - 1)
    x_low = np.arange(1, low_freq_range + 1)
    sfs_low = sfs[1:low_freq_range + 1]
    ax.bar(x_low, sfs_low, color='coral', alpha=0.7, edgecolor='black')
    ax.set_xlabel('Derived Allele Count', fontweight='bold')
    ax.set_ylabel('Number of SNPs', fontweight='bold')
    ax.set_title('Low Frequency Variants (DAF 1-20)', fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Plot 4: Folded vs Unfolded comparison
    ax = axes[1, 1]
    # Create folded SFS for comparison
    n_fold = (n_bins - 1) // 2 + 1
    sfs_folded = np.zeros(n_fold)
    for i in range(1, n_bins):
        if i <= n_bins // 2:
            sfs_folded[i] += sfs[i]
        else:
            sfs_folded[n_bins - i] += sfs[i]
    
    x_fold = np.arange(n_fold)
    ax.bar(x_fold, sfs_folded, color='purple', alpha=0.7, edgecolor='black', label='Folded (for reference)')
    ax.set_xlabel('Minor Allele Count', fontweight='bold')
    ax.set_ylabel('Number of SNPs', fontweight='bold')
    ax.set_title('Folded SFS (for comparison)', fontweight='bold')
    ax.set_yscale('log')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Plot saved: {output_path}")


def main():
    """Main workflow."""
    logger.info("=" * 80)
    logger.info("SITE FREQUENCY SPECTRUM GENERATION FOR FASTSIMCOAL2")
    logger.info("=" * 80)
    
    # Check input files exist
    if not MAIN_VCF.exists():
        logger.error(f"Main VCF not found: {MAIN_VCF}")
        sys.exit(1)
    
    if not ANCESTRAL_VCF.exists():
        logger.error(f"Ancestral VCF not found: {ANCESTRAL_VCF}")
        sys.exit(1)
    
    # Load ancestral states
    ancestral_dict = load_ancestral_states(ANCESTRAL_VCF)
    
    # Determine number of samples from VCF header
    n_samples = None
    
    open_func = gzip.open if str(MAIN_VCF).endswith('.gz') else open
    with open_func(MAIN_VCF, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                # Sample names start at column 10 (index 9)
                samples = line.strip().split('\t')[9:]
                n_samples = len(samples)
                logger.info(f"Auto-detected {n_samples} diploid samples from VCF header")
                break
            if not line.startswith('#'):
                break
    
    # Fallback to known value if detection fails
    if n_samples is None:
        n_samples = 68  # Known from data
        logger.warning(f"Could not auto-detect sample size, using default: {n_samples}")
    
    # Generate SFS
    sfs, stats = generate_sfs(MAIN_VCF, ancestral_dict, n_samples)
    
    # Write fastsimcoal2 format
    sfs_file = OUTPUT_DIR / "SNJ_DAFpop0.obs"
    write_fastsimcoal2_sfs(sfs, sfs_file)
    
    # Write statistics
    stats_file = OUTPUT_DIR / "sfs_statistics.txt"
    write_statistics(sfs, stats, stats_file)
    
    # Create plot
    plot_file = OUTPUT_DIR / "sfs_plot.png"
    plot_sfs(sfs, plot_file)
    
    logger.info("\n" + "=" * 80)
    logger.info("SFS GENERATION COMPLETE")
    logger.info("=" * 80)
    logger.info(f"\nOutput files:")
    logger.info(f"  SFS: {sfs_file}")
    logger.info(f"  Statistics: {stats_file}")
    logger.info(f"  Plot: {plot_file}")
    logger.info(f"\nPolarized {stats['polarized']:,} variants out of {stats['biallelic']:,} biallelic sites")
    logger.info(f"Polarization success rate: {stats['polarized']/max(stats['biallelic'],1)*100:.1f}%\n")


if __name__ == "__main__":
    main()
