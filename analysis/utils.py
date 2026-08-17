#!/usr/bin/env python3
"""
Utility functions for monkey inbreeding analysis
Common functions used across multiple phases
"""

import os
import sys
import logging
import yaml
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
import gzip
from typing import Dict, List, Tuple, Optional, Union

# ============================================================================
# CONFIGURATION LOADING
# ============================================================================

def load_config(config_path: str = "config.yaml") -> Dict:
    """
    Load configuration from YAML file
    
    Args:
        config_path: Path to config.yaml file
        
    Returns:
        Dictionary containing configuration
    """
    if not os.path.exists(config_path):
        # Try in scripts directory
        script_dir = os.path.dirname(os.path.abspath(__file__))
        config_path = os.path.join(script_dir, "config.yaml")
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    root = os.environ.get("PLM_BASE_DIR", "").strip()
    if not root:
        env_file = Path(__file__).resolve().parent / "base_dir.env"
        if env_file.is_file():
            for line in env_file.read_text(encoding="utf-8").splitlines():
                stripped = line.strip()
                if stripped.startswith("export "):
                    stripped = stripped[len("export "):].strip()
                if stripped.startswith("PLM_BASE_DIR="):
                    root = stripped.split("=", 1)[1].strip().strip("'").strip('"')
                    break
    if root:
        config = _expand_analysis_root(config, root)

    return config


def _expand_analysis_root(obj, root: str):
    """Replace /path/to/analysis_root placeholders with the configured root."""
    if isinstance(obj, str):
        return obj.replace("/path/to/analysis_root", root)
    if isinstance(obj, dict):
        return {k: _expand_analysis_root(v, root) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_expand_analysis_root(v, root) for v in obj]
    return obj

# ============================================================================
# LOGGING SETUP
# ============================================================================

def setup_logger(name: str, log_file: Optional[str] = None, 
                 level: str = "INFO") -> logging.Logger:
    """
    Setup logger with consistent formatting
    
    Args:
        name: Logger name
        log_file: Optional path to log file
        level: Logging level (DEBUG, INFO, WARNING, ERROR)
        
    Returns:
        Configured logger
    """
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, level.upper()))
    
    # Remove existing handlers to avoid duplicates
    logger.handlers.clear()
    
    # Prevent propagation to root logger
    logger.propagate = False
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(getattr(logging, level.upper()))
    
    # Formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler if specified
    if log_file:
        log_dir = os.path.dirname(log_file)
        if log_dir:  # Only create directory if there is a directory component
            os.makedirs(log_dir, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(getattr(logging, level.upper()))
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger

# ============================================================================
# FILE I/O UTILITIES
# ============================================================================

def read_vcf_info(vcf_path: str) -> Tuple[int, int, List[str]]:
    """
    Read basic VCF information without loading entire file
    
    Args:
        vcf_path: Path to VCF file (can be gzipped)
        
    Returns:
        Tuple of (num_variants, num_samples, sample_ids)
    """
    opener = gzip.open if vcf_path.endswith('.gz') else open
    
    num_variants = 0
    sample_ids = []
    
    with opener(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                continue
            elif line.startswith('#CHROM'):
                # Header line with sample IDs
                fields = line.strip().split('\t')
                sample_ids = fields[9:]  # Everything after FORMAT
                break
        
        # Count remaining lines (variants)
        num_variants = sum(1 for _ in f)
    
    return num_variants, len(sample_ids), sample_ids

def read_plink_fam(fam_path: str) -> pd.DataFrame:
    """
    Read PLINK .fam file
    
    Args:
        fam_path: Path to .fam file
        
    Returns:
        DataFrame with family information
    """
    cols = ['FID', 'IID', 'Father', 'Mother', 'Sex', 'Phenotype']
    df = pd.read_csv(fam_path, sep=r'\s+', header=None, names=cols)
    return df

def read_plink_bim(bim_path: str) -> pd.DataFrame:
    """
    Read PLINK .bim file
    
    Args:
        bim_path: Path to .bim file
        
    Returns:
        DataFrame with variant information
    """
    cols = ['CHR', 'SNP', 'cM', 'POS', 'A1', 'A2']
    df = pd.read_csv(bim_path, sep=r'\s+', header=None, names=cols)
    return df

def parse_vcf_header(vcf_path: str) -> Dict:
    """
    Parse VCF header information
    
    Args:
        vcf_path: Path to VCF file
        
    Returns:
        Dictionary with header information
    """
    opener = gzip.open if vcf_path.endswith('.gz') else open
    
    header_info = {
        'contigs': [],
        'info_fields': {},
        'format_fields': {},
        'filters': {}
    }
    
    with opener(vcf_path, 'rt') as f:
        for line in f:
            if not line.startswith('##'):
                break
            
            if line.startswith('##contig='):
                contig = line.split('ID=')[1].split(',')[0].split('>')[0]
                header_info['contigs'].append(contig)
            elif line.startswith('##INFO='):
                # Parse INFO field
                pass  # Add parsing logic if needed
            elif line.startswith('##FORMAT='):
                # Parse FORMAT field
                pass
    
    return header_info

# ============================================================================
# VCF PARSING UTILITIES
# ============================================================================

def parse_vcf_line(line: str) -> Dict:
    """
    Parse a single VCF variant line
    
    Args:
        line: VCF line string
        
    Returns:
        Dictionary with variant information
    """
    fields = line.strip().split('\t')
    
    variant = {
        'CHROM': fields[0],
        'POS': int(fields[1]),
        'ID': fields[2],
        'REF': fields[3],
        'ALT': fields[4],
        'QUAL': fields[5],
        'FILTER': fields[6],
        'INFO': fields[7]
    }
    
    return variant

def parse_info_field(info_str: str) -> Dict:
    """
    Parse VCF INFO field
    
    Args:
        info_str: INFO field string
        
    Returns:
        Dictionary of INFO tags and values
    """
    info_dict = {}
    
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[item] = True
    
    return info_dict

def extract_snpeff_annotation(info_str: str) -> List[Dict]:
    """
    Extract SnpEff annotations from INFO field
    
    Args:
        info_str: INFO field string
        
    Returns:
        List of annotation dictionaries
    """
    annotations = []
    
    info_dict = parse_info_field(info_str)
    
    if 'ANN' in info_dict:
        # SnpEff annotation format:
        # Allele|Annotation|Impact|Gene_Name|Gene_ID|Feature_Type|Feature_ID|...
        for ann in info_dict['ANN'].split(','):
            fields = ann.split('|')
            if len(fields) >= 15:
                annotations.append({
                    'Allele': fields[0],
                    'Annotation': fields[1],
                    'Impact': fields[2],
                    'Gene_Name': fields[3],
                    'Gene_ID': fields[4],
                    'Feature_Type': fields[5],
                    'Feature_ID': fields[6],
                    'Transcript_BioType': fields[7],
                    'Rank': fields[8],
                    'HGVS_c': fields[9],
                    'HGVS_p': fields[10],
                    'cDNA_pos': fields[11],
                    'CDS_pos': fields[12],
                    'Protein_pos': fields[13],
                    'Distance': fields[14]
                })
    
    return annotations

# ============================================================================
# GENOTYPE UTILITIES
# ============================================================================

def calculate_allele_frequency(genotypes: List[str]) -> Tuple[float, int, int]:
    """
    Calculate allele frequency from genotype list
    
    Args:
        genotypes: List of genotype strings (e.g., ['0/0', '0/1', '1/1'])
        
    Returns:
        Tuple of (allele_frequency, allele_count, total_alleles)
    """
    alt_count = 0
    total_alleles = 0
    
    for gt in genotypes:
        if gt in ['./.', '.|.', '.']:
            continue
        
        alleles = gt.replace('|', '/').split('/')
        
        for allele in alleles:
            if allele != '.':
                total_alleles += 1
                if allele != '0':
                    alt_count += 1
    
    if total_alleles == 0:
        return 0.0, 0, 0
    
    freq = alt_count / total_alleles
    return freq, alt_count, total_alleles

def genotype_to_012(genotype: str, missing_value: int = -9) -> int:
    """
    Convert genotype string to 0/1/2 coding
    
    Args:
        genotype: Genotype string (e.g., '0/1', '1|1')
        missing_value: Value for missing genotypes
        
    Returns:
        Integer genotype (0, 1, 2, or missing_value)
    """
    if genotype in ['./.', '.|.', '.']:
        return missing_value
    
    alleles = genotype.replace('|', '/').split('/')
    
    try:
        return sum(int(a) for a in alleles if a != '.')
    except:
        return missing_value

def is_homozygous_alt(genotype: str) -> bool:
    """Check if genotype is homozygous alternate"""
    if genotype in ['./.', '.|.', '.']:
        return False
    
    alleles = genotype.replace('|', '/').split('/')
    return len(set(alleles)) == 1 and alleles[0] not in ['0', '.']

def is_heterozygous(genotype: str) -> bool:
    """Check if genotype is heterozygous"""
    if genotype in ['./.', '.|.', '.']:
        return False
    
    alleles = genotype.replace('|', '/').split('/')
    alleles = [a for a in alleles if a != '.']
    
    if len(alleles) < 2:
        return False
    
    return len(set(alleles)) > 1

# ============================================================================
# STATISTICAL UTILITIES
# ============================================================================

def calculate_heterozygosity(genotypes: List[int], missing_value: int = -9) -> float:
    """
    Calculate observed heterozygosity
    
    Args:
        genotypes: List of genotypes in 0/1/2 format
        missing_value: Value indicating missing data
        
    Returns:
        Observed heterozygosity
    """
    valid_genotypes = [g for g in genotypes if g != missing_value]
    
    if len(valid_genotypes) == 0:
        return 0.0
    
    het_count = sum(1 for g in valid_genotypes if g == 1)
    return het_count / len(valid_genotypes)

def calculate_tajimas_d(pi: float, theta_w: float, n: int) -> float:
    """
    Calculate Tajima's D
    
    Args:
        pi: Nucleotide diversity (average pairwise differences)
        theta_w: Watterson's theta (from number of segregating sites)
        n: Sample size (number of sequences)
        
    Returns:
        Tajima's D value
    """
    if n < 2:
        return np.nan
    
    # Calculate a1 and a2
    a1 = sum(1.0 / i for i in range(1, n))
    a2 = sum(1.0 / (i * i) for i in range(1, n))
    
    # Calculate b1 and b2
    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n * n + n + 3) / (9 * n * (n - 1))
    
    # Calculate c1 and c2
    c1 = b1 - 1 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1 * a1)
    
    # Calculate e1 and e2
    e1 = c1 / a1
    e2 = c2 / (a1 * a1 + a2)
    
    # Calculate variance
    var_d = e1 * theta_w + e2 * theta_w * theta_w
    
    if var_d <= 0:
        return np.nan
    
    # Calculate D
    d = (pi - theta_w) / np.sqrt(var_d)
    
    return d

# ============================================================================
# ROH UTILITIES
# ============================================================================

def classify_roh_length(length_bp: int, 
                       short_max: int = 1000000,
                       medium_max: int = 3000000) -> str:
    """
    Classify ROH by length
    
    Args:
        length_bp: ROH length in base pairs
        short_max: Maximum length for short ROH (default 1 Mb)
        medium_max: Maximum length for medium ROH (default 3 Mb)
        
    Returns:
        Category: 'short', 'medium', or 'long'
    """
    if length_bp < short_max:
        return 'short'
    elif length_bp < medium_max:
        return 'medium'
    else:
        return 'long'

def calculate_froh(roh_total_length: float, genome_length: float) -> float:
    """
    Calculate inbreeding coefficient from ROH
    
    Args:
        roh_total_length: Total length of ROH segments (bp)
        genome_length: Autosomal genome length (bp)
        
    Returns:
        F_ROH
    """
    return roh_total_length / genome_length

# ============================================================================
# DATA FORMATTING UTILITIES
# ============================================================================

def format_pvalue(pval: float, threshold: float = 0.001) -> str:
    """
    Format p-value for display
    
    Args:
        pval: P-value
        threshold: Threshold for scientific notation
        
    Returns:
        Formatted string
    """
    if pval < threshold:
        return f"{pval:.2e}"
    else:
        return f"{pval:.4f}"

def format_number(num: float, decimals: int = 2) -> str:
    """Format number with comma separators"""
    if abs(num) >= 1000:
        return f"{num:,.{decimals}f}"
    else:
        return f"{num:.{decimals}f}"

def create_output_dir(path: str) -> None:
    """Create output directory if it doesn't exist"""
    os.makedirs(path, exist_ok=True)

# ============================================================================
# PROGRESS TRACKING
# ============================================================================

class ProgressTracker:
    """Simple progress tracker for long-running operations"""
    
    def __init__(self, total: int, desc: str = "Progress", 
                 update_interval: int = 100):
        self.total = total
        self.current = 0
        self.desc = desc
        self.update_interval = update_interval
        self.start_time = datetime.now()
    
    def update(self, n: int = 1):
        """Update progress"""
        self.current += n
        
        if self.current % self.update_interval == 0 or self.current == self.total:
            elapsed = (datetime.now() - self.start_time).total_seconds()
            rate = self.current / elapsed if elapsed > 0 else 0
            eta = (self.total - self.current) / rate if rate > 0 else 0
            
            pct = 100 * self.current / self.total
            print(f"{self.desc}: {self.current}/{self.total} ({pct:.1f}%) "
                  f"- {rate:.1f} it/s - ETA: {eta:.0f}s", flush=True)
    
    def finish(self):
        """Mark as complete"""
        elapsed = (datetime.now() - self.start_time).total_seconds()
        print(f"{self.desc}: Complete! Total time: {elapsed:.1f}s")

# ============================================================================
# VALIDATION UTILITIES
# ============================================================================

def validate_file_exists(filepath: str, description: str = "File") -> None:
    """
    Validate that a file exists
    
    Args:
        filepath: Path to file
        description: Description for error message
        
    Raises:
        FileNotFoundError if file doesn't exist
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"{description} not found: {filepath}")

def validate_output_dir(dirpath: str) -> None:
    """
    Validate and create output directory
    
    Args:
        dirpath: Directory path
    """
    os.makedirs(dirpath, exist_ok=True)
    
    # Test write permissions
    test_file = os.path.join(dirpath, '.write_test')
    try:
        with open(test_file, 'w') as f:
            f.write('test')
        os.remove(test_file)
    except:
        raise PermissionError(f"Cannot write to directory: {dirpath}")

# ============================================================================
# MAIN (for testing)
# ============================================================================

if __name__ == "__main__":
    # Test configuration loading
    try:
        config = load_config()
        print("Configuration loaded successfully")
        print(f"Project directory: {config['project']['base_dir']}")
    except Exception as e:
        print(f"Error loading configuration: {e}")


