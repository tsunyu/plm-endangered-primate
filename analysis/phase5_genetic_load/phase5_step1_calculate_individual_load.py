#!/usr/bin/env python3
"""
Phase 5: Calculate Individual Genetic Load

Calculates per-individual genetic load using sigmoid-calibrated pathogenicity probabilities:

Method:
- Uses raw ESM-2 LLR scores (no Z-normalization needed)
- Converts to pathogenicity probability via calibrated sigmoid function:
  P(Pathogenic) = 1 / (1 + exp(0.5287 * (x + 6.8920)))
- Pathogenicity probability directly represents fitness reduction

Genetic Load Metrics:
- Total Genetic Load = Realized Load + Potential Load
  (Maximum theoretical burden if all variants were fully expressed/homozygous)
  (Formula: Σ(P for het) + Σ(P for hom))
  
- Realized Load = Homozygous Realized Load + Heterozygous Realized Load
  (Currently expressed fitness reduction given dominance coefficient h)
  (Formula: Σ(h×P for het) + Σ(P for hom))
  
- Potential Load = Σ[(1-h)×P for het]
  (Hidden load masked by heterozygosity, would be exposed by inbreeding)
  
- Homozygous Realized Load = Σ(P for hom)
  (Fully expressed load from homozygous deleterious variants)
  
- Heterozygous Realized Load = Σ(h×P for het)
  (Partially expressed load from heterozygous deleterious variants)

Where h = dominance coefficient (default 0.25 for partially recessive model)

Additional metrics:
- LOF and deleterious missense counts (het/hom)
- ROH-associated load by ROH length category

Usage: python3 phase5_step1_calculate_individual_load.py
"""

import sys
import os
import gzip
import csv
import pandas as pd
import numpy as np
from collections import defaultdict
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir

# Visualization imports
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils import setup_logger, load_config

# ============================================================================
# SETUP
# ============================================================================

BASE_DIR = get_base_dir()
OUTPUT_DIR = f"{BASE_DIR}/output/phase5_genetic_load/individual_load"

# Input files
ANNOTATED_VCF = f"{BASE_DIR}/output/phase2_annotation/snpeff_annotation/annotated_variants.vcf.gz"
PLM_PREDICTIONS = f"{BASE_DIR}/output/phase4_plm_predictions/ensemble/ensemble_predictions.csv"
ROH_DATA = f"{BASE_DIR}/output/phase3a_population_genomics/roh_analysis/plink_roh.hom"
FROH_DATA = f"{BASE_DIR}/output/phase3a_population_genomics/roh_analysis/roh_summary_per_individual.csv"

os.makedirs(OUTPUT_DIR, exist_ok=True)

logger = setup_logger("phase5_load", f"{OUTPUT_DIR}/individual_load.log")

# ============================================================================
# PARAMETERS
# ============================================================================

# Sigmoid function parameters for pathogenicity probability
# P(Pathogenic) = 1 / (1 + exp(k * (x - x0)))
# where k = 0.5287, x0 = -6.8920
SIGMOID_K = 0.5287       # Steepness/scaling factor
SIGMOID_X0 = -6.8920     # Midpoint (inflection point where P = 0.5)

# Pathogenicity probability for LOF variants (assumed highly deleterious)
LOF_PATHOGENICITY = 0.95

# Dominance coefficient for heterozygous effects
# h=0 means fully recessive (het has no effect)
# h=0.5 means additive (het has half the effect of hom)
# h=1 means fully dominant (het has same effect as hom)
H_DEFAULT = 0.25  # Partially recessive (het has 25% of hom effect)

# ============================================================================
# FUNCTIONS
# ============================================================================

def load_plm_predictions():
    """
    Load PLM predictions and create variant lookup
    
    Uses raw ESM-2 LLR scores directly (no Z-normalization needed)
    
    Returns dict with variant_id -> {esm2_score, ensemble_score, is_deleterious}
    """
    logger.info("Loading PLM predictions...")
    
    try:
        plm_df = pd.read_csv(PLM_PREDICTIONS)
        
        # Log ESM-2 score statistics
        esm2_scores = plm_df['esm2_score'].dropna()
        
        if len(esm2_scores) > 0:
            logger.info(f"  ESM-2 LLR score statistics:")
            logger.info(f"    Total variants: {len(esm2_scores)}")
            logger.info(f"    Range: {esm2_scores.min():.2f} to {esm2_scores.max():.2f}")
            logger.info(f"    Mean: {esm2_scores.mean():.3f}")
            logger.info(f"    Median: {esm2_scores.median():.3f}")
            
            # Show pathogenicity probability at key thresholds
            logger.info(f"  Sigmoid pathogenicity calibration (k={SIGMOID_K}, x0={SIGMOID_X0}):")
            test_scores = [-15, -10, -6.89, -5, 0, 5]
            for ts in test_scores:
                prob = sigmoid_pathogenicity(ts)
                logger.info(f"    ESM-2 = {ts:6.2f} -> P(Pathogenic) = {prob:.4f}")
        
        # Create variant ID to score mapping (using raw scores)
        variant_scores = {}
        for _, row in plm_df.iterrows():
            var_id = row['variant_id']
            raw_score = row.get('esm2_score', None)
            
            variant_scores[var_id] = {
                'esm2_score': raw_score,                     # Raw ESM-2 LLR score
                'ensemble_score': row['ensemble_score'],     # 0-1 normalized (for compatibility)
                'is_deleterious': row['is_deleterious']
            }
        
        logger.info(f"  Loaded predictions for {len(variant_scores)} variants")
        
        return variant_scores
        
    except Exception as e:
        logger.warning(f"Could not load PLM predictions: {e}")
        return {}

def load_roh_data():
    """
    Load ROH segments per individual
    
    Converts PLINK chromosome format (1, 2, ...) to RefSeq format (NC_044549.1, ...)
    """
    logger.info("Loading ROH data...")
    
    # Chromosome mapping: PLINK format -> RefSeq format
    chrom_map = {
        '1': 'NC_044549.1', '2': 'NC_044550.1', '3': 'NC_044551.1',
        '4': 'NC_044552.1', '5': 'NC_044553.1', '6': 'NC_044554.1',
        '7': 'NC_044555.1', '8': 'NC_044556.1', '9': 'NC_044557.1',
        '10': 'NC_044558.1', '11': 'NC_044559.1', '12': 'NC_044560.1',
        '13': 'NC_044561.1', '14': 'NC_044562.1', '15': 'NC_044563.1',
        '16': 'NC_044564.1', '17': 'NC_044565.1', '18': 'NC_044566.1',
        '19': 'NC_044567.1', '20': 'NC_044568.1', '21': 'NC_044569.1',
        '22': 'NC_044570.1'
    }
    
    try:
        roh_df = pd.read_csv(ROH_DATA, sep=r'\s+')
        
        # Index ROH by individual
        roh_by_ind = defaultdict(list)
        conversion_warnings = set()
        
        for _, row in roh_df.iterrows():
            iid = row['IID']
            plink_chr = str(row['CHR'])
            
            # Convert PLINK chromosome to RefSeq format
            if plink_chr in chrom_map:
                refseq_chr = chrom_map[plink_chr]
            else:
                # Handle non-standard chromosomes
                if plink_chr not in conversion_warnings:
                    logger.warning(f"  Unknown chromosome in ROH data: {plink_chr}")
                    conversion_warnings.add(plink_chr)
                continue
            
            roh_by_ind[iid].append({
                'chr': refseq_chr,  # Now in RefSeq format
                'start': row['POS1'],
                'end': row['POS2'],
                'length': row['KB'] * 1000
            })
        
        logger.info(f"  Loaded ROH for {len(roh_by_ind)} individuals")
        total_roh = sum(len(segs) for segs in roh_by_ind.values())
        logger.info(f"  Total ROH segments: {total_roh}")
        return roh_by_ind
        
    except Exception as e:
        logger.warning(f"Could not load ROH data: {e}")
        return {}

def sigmoid_pathogenicity(esm2_score):
    """
    Calculate pathogenicity probability from raw ESM-2 LLR score using sigmoid function
    
    The calibrated probability of pathogenicity is computed as:
        P(Pathogenic) = 1 / (1 + exp(k * (x - x0)))
    
    With calibrated parameters:
        k = 0.5287 (steepness)
        x0 = -6.8920 (midpoint where P = 0.5)
    
    Simplified formula:
        P(Pathogenic) = 1 / (1 + exp(0.5287 * (x + 6.8920)))
    
    Args:
        esm2_score: Raw ESM-2 LLR score (log-likelihood ratio)
                    Negative values indicate deleterious mutations
                    More negative = more likely pathogenic
        
    Returns:
        float: Probability of pathogenicity (0 to 1)
               Higher values indicate more likely to reduce fitness
    """
    if esm2_score is None:
        return 0.0
    
    # Sigmoid function with calibrated parameters
    # P = 1 / (1 + exp(k * (x - x0))) = 1 / (1 + exp(0.5287 * (x + 6.8920)))
    exponent = SIGMOID_K * (esm2_score - SIGMOID_X0)
    prob = 1.0 / (1.0 + np.exp(exponent))
    
    return prob

def is_in_roh(chrom, pos, roh_list):
    """
    Check if a position is within any ROH segment and classify by length
    
    Returns:
        tuple: (in_roh, roh_length, roh_category)
        
    ROH categories:
    - 'short': < 1 Mb (ancient inbreeding, >50 generations)
    - 'medium': 1-5 Mb (intermediate inbreeding, 10-50 generations)  
    - 'long': > 5 Mb (recent inbreeding, <10 generations)
    """
    for roh in roh_list:
        if roh['chr'] == chrom and roh['start'] <= pos <= roh['end']:
            length_mb = roh['length'] / 1e6
            
            # Classify ROH by length
            if length_mb < 1.0:
                category = 'short'
            elif length_mb < 5.0:
                category = 'medium'
            else:
                category = 'long'
                
            return True, roh['length'], category
    
    return False, 0, None

def classify_variant_impact(info_dict, plm_scores, variant_id):
    """
    Classify variant as LOF, deleterious missense, or other
    
    Uses sigmoid function to calculate pathogenicity probability from raw ESM-2 LLR scores:
        P(Pathogenic) = 1 / (1 + exp(0.5287 * (x + 6.8920)))
    
    The pathogenicity probability directly represents the fitness reduction caused by the mutation.
    
    Args:
        info_dict: Parsed VCF INFO field
        plm_scores: Dict of variant_id -> PLM prediction data
        variant_id: Variant identifier (chrom:pos:ref:alt)
    
    Returns:
        dict: Classification information with pathogenicity probability
    """
    classification = {
        'is_lof': False,
        'is_deleterious_missense': False,
        'impact': 'MODIFIER',
        'pathogenicity_prob': 0.0,  # Probability of being pathogenic (fitness reduction)
        'dominance': H_DEFAULT
    }
    
    # Check for LOF
    if 'LOF' in info_dict:
        classification['is_lof'] = True
        classification['impact'] = 'HIGH'
        classification['pathogenicity_prob'] = LOF_PATHOGENICITY
        return classification
    
    # Check SnpEff annotation
    if 'ANN' in info_dict:
        anns = info_dict['ANN'].split(',')
        for ann in anns:
            fields = ann.split('|')
            if len(fields) >= 3:
                effect = fields[1]
                impact = fields[2]
                
                if impact == 'HIGH':
                    classification['is_lof'] = True
                    classification['impact'] = 'HIGH'
                    classification['pathogenicity_prob'] = LOF_PATHOGENICITY
                    return classification
                
                if 'missense' in effect.lower():
                    # Check PLM prediction - ONLY classify as deleterious if PLM says so
                    if variant_id in plm_scores and plm_scores[variant_id]['is_deleterious']:
                        classification['is_deleterious_missense'] = True
                        classification['impact'] = 'MODERATE'
                        
                        # Use raw ESM-2 LLR score with sigmoid function
                        esm2_score = plm_scores[variant_id].get('esm2_score', None)
                        
                        if esm2_score is not None:
                            # Calculate pathogenicity probability using sigmoid
                            classification['pathogenicity_prob'] = sigmoid_pathogenicity(esm2_score)
                        else:
                            # Fallback: use ensemble_score if ESM-2 score unavailable
                            ensemble_score = plm_scores[variant_id].get('ensemble_score', 0.5)
                            classification['pathogenicity_prob'] = ensemble_score
                    # If missense but not predicted deleterious by PLM, don't count it
    
    return classification

def calculate_individual_load(vcf_file, plm_scores, roh_by_ind):
    """
    Calculate genetic load for each individual
    
    Uses sigmoid-based pathogenicity probability model:
    - Each mutation has a pathogenicity probability P calculated via sigmoid function:
      P(Pathogenic) = 1 / (1 + exp(0.5287 * (x + 6.8920)))
    
    Genetic Load Metrics:
    - Total Genetic Load = Σ(P for het) + Σ(P for hom)
      Maximum theoretical burden if all variants were fully expressed (e.g., all homozygous)
      
    - Realized Load = Σ(h×P for het) + Σ(P for hom)
      Currently expressed fitness reduction given dominance coefficient h
      
    - Potential Load = Σ((1-h)×P for het)
      Hidden load masked by heterozygosity, would be exposed by inbreeding
      
    - Homozygous Realized Load = Σ(P for hom)
    - Heterozygous Realized Load = Σ(h×P for het)
    
    Relationships:
    - Total Genetic Load = Realized Load + Potential Load
    - Realized Load = Homozygous Realized Load + Heterozygous Realized Load
    
    Returns:
        DataFrame with per-individual load statistics
    """
    logger.info("Calculating individual genetic load...")
    logger.info(f"  Using sigmoid pathogenicity: P = 1 / (1 + exp({SIGMOID_K} * (x - ({SIGMOID_X0}))))")
    logger.info(f"  Dominance coefficient (h) = {H_DEFAULT} (partially recessive)")
    logger.info(f"  Total Genetic Load = Realized Load + Potential Load (max theoretical burden)")
    logger.info(f"  Realized Load = Hom_Realized + Het_Realized (currently expressed)")
    
    # Initialize load counters
    individual_load = defaultdict(lambda: {
        'lof_het': 0,
        'lof_hom': 0,
        'del_mis_het': 0,
        'del_mis_hom': 0,
        'synonymous': 0,
        'lof_in_roh': 0,
        'del_mis_in_roh': 0,
        # Pathogenicity-based load tracking
        # Realized Load components:
        'sum_het_realized': 0.0,   # Heterozygous Realized Load: Σ(h×P for het)
        'sum_hom_realized': 0.0,   # Homozygous Realized Load: Σ(P for hom)
        # For calculating Potential Load:
        'sum_het_full_path': 0.0,  # Full pathogenicity of het variants: Σ(P for het)
        # ROH-specific tracking
        'sum_path_in_roh': 0.0,    # Pathogenicity sum from variants in ROH
        'sum_path_outside_roh': 0.0,  # Pathogenicity sum from variants outside ROH
        'variants_in_short_roh': 0,   # < 1 Mb (ancient)
        'variants_in_medium_roh': 0,  # 1-5 Mb (intermediate)
        'variants_in_long_roh': 0,    # > 5 Mb (recent inbreeding)
        'sum_path_in_short_roh': 0.0,
        'sum_path_in_medium_roh': 0.0,
        'sum_path_in_long_roh': 0.0,
    })
    
    sample_ids = []
    
    with gzip.open(vcf_file, 'rt') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('##'):
                continue
            
            if line.startswith('#CHROM'):
                # Extract sample IDs
                fields = line.strip().split('\t')
                sample_ids = fields[9:]
                logger.info(f"  Found {len(sample_ids)} samples")
                continue
            
            # Parse variant
            fields = line.strip().split('\t')
            chrom, pos, vid, ref, alt, qual, filt, info = fields[:8]
            genotypes = fields[9:] if len(fields) > 9 else []
            
            pos = int(pos)
            
            # Parse INFO
            info_dict = {}
            for item in info.split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    info_dict[key] = value
            
            # Classify variant
            var_id = f"{chrom}:{pos}:{ref}:{alt}"
            classification = classify_variant_impact(info_dict, plm_scores, var_id)
            
            # Skip if not deleterious
            if not classification['is_lof'] and not classification['is_deleterious_missense']:
                continue
            
            # Process genotypes
            for sample_id, gt_field in zip(sample_ids, genotypes):
                gt = gt_field.split(':')[0]  # Extract genotype
                
                if gt in ['./.', '.|.', '.']:
                    continue
                
                # Parse genotype
                alleles = gt.replace('|', '/').split('/')
                
                try:
                    allele_sum = sum(int(a) for a in alleles if a != '.')
                except:
                    continue
                
                is_het = (allele_sum == 1)
                is_hom = (allele_sum == 2)
                
                if not is_het and not is_hom:
                    continue
                
                # Check if in ROH
                in_roh, roh_length, roh_category = is_in_roh(chrom, pos, roh_by_ind.get(sample_id, []))
                
                # Get pathogenicity probability and dominance
                path_prob = classification['pathogenicity_prob']
                h = classification['dominance']
                
                # Calculate realized load contribution
                # Heterozygous: h × P (partially expressed)
                # Homozygous: P (fully expressed)
                if is_het:
                    realized_effect = h * path_prob  # Heterozygous Realized Load component
                    # Track full pathogenicity for Potential Load calculation
                    individual_load[sample_id]['sum_het_full_path'] += path_prob
                    individual_load[sample_id]['sum_het_realized'] += realized_effect
                elif is_hom:
                    realized_effect = path_prob  # Homozygous Realized Load component
                    individual_load[sample_id]['sum_hom_realized'] += realized_effect
                else:
                    realized_effect = 0.0
                
                # Update variant counts
                if classification['is_lof']:
                    if is_het:
                        individual_load[sample_id]['lof_het'] += 1
                    elif is_hom:
                        individual_load[sample_id]['lof_hom'] += 1
                        if in_roh:
                            individual_load[sample_id]['lof_in_roh'] += 1
                
                elif classification['is_deleterious_missense']:
                    if is_het:
                        individual_load[sample_id]['del_mis_het'] += 1
                    elif is_hom:
                        individual_load[sample_id]['del_mis_hom'] += 1
                        if in_roh:
                            individual_load[sample_id]['del_mis_in_roh'] += 1
                
                # Track ROH-specific load (only for homozygous realized load)
                if is_hom:
                    if in_roh:
                        # Variant in ROH
                        individual_load[sample_id]['sum_path_in_roh'] += realized_effect
                        individual_load[sample_id]['sum_path_in_' + roh_category + '_roh'] += realized_effect
                        individual_load[sample_id]['variants_in_' + roh_category + '_roh'] += 1
                    else:
                        # Variant outside ROH
                        individual_load[sample_id]['sum_path_outside_roh'] += realized_effect
            
            if line_num % 10000 == 0:
                logger.info(f"  Processed {line_num} variants")
    
    # Convert to DataFrame
    load_data = []
    for sample_id in sample_ids:
        load = individual_load[sample_id]
        
        # Calculate derived load metrics
        het_realized = load['sum_het_realized']      # Σ(h×P for het)
        hom_realized = load['sum_hom_realized']      # Σ(P for hom)
        het_full_path = load['sum_het_full_path']    # Σ(P for het)
        
        realized_load = het_realized + hom_realized  # Σ(h×P for het) + Σ(P for hom)
        potential_load = het_full_path - het_realized  # Σ((1-h)×P for het)
        total_genetic_load = het_full_path + hom_realized  # Σ(P for het) + Σ(P for hom)
        
        load_data.append({
            'IID': sample_id,
            # Variant counts
            'LOF_Het': load['lof_het'],
            'LOF_Hom': load['lof_hom'],
            'LOF_Total': load['lof_het'] + load['lof_hom'],
            'DelMis_Het': load['del_mis_het'],
            'DelMis_Hom': load['del_mis_hom'],
            'DelMis_Total': load['del_mis_het'] + load['del_mis_hom'],
            'Total_Het': load['lof_het'] + load['del_mis_het'],
            'Total_Hom': load['lof_hom'] + load['del_mis_hom'],
            'Total_Deleterious': load['lof_het'] + load['lof_hom'] + load['del_mis_het'] + load['del_mis_hom'],
            'LOF_in_ROH': load['lof_in_roh'],
            'DelMis_in_ROH': load['del_mis_in_roh'],
            # Genetic load metrics (pathogenicity-based)
            # Total Genetic Load = Realized Load + Potential Load
            'Total_Genetic_Load': total_genetic_load,         # Σ(P for het) + Σ(P for hom) - max theoretical burden
            'Realized_Load': realized_load,                   # Σ(h×P for het) + Σ(P for hom) - currently expressed
            'Potential_Load': potential_load,                 # Σ((1-h)×P for het) - hidden, exposed by inbreeding
            'Hom_Realized_Load': hom_realized,                # Σ(P for hom) - fully expressed
            'Het_Realized_Load': het_realized,                # Σ(h×P for het) - partially expressed (h=0.25)
            # ROH-specific metrics (based on homozygous realized load in ROH)
            'Load_in_ROH': load['sum_path_in_roh'],           # Hom realized load from variants in any ROH
            'Load_outside_ROH': load['sum_path_outside_roh'], # Hom realized load from variants outside ROH
            'Load_in_Short_ROH': load['sum_path_in_short_roh'],    # From ancient inbreeding (<1 Mb)
            'Load_in_Medium_ROH': load['sum_path_in_medium_roh'],  # From intermediate inbreeding (1-5 Mb)
            'Load_in_Long_ROH': load['sum_path_in_long_roh'],      # From recent inbreeding (>5 Mb)
            'Variants_in_Short_ROH': load['variants_in_short_roh'],
            'Variants_in_Medium_ROH': load['variants_in_medium_roh'],
            'Variants_in_Long_ROH': load['variants_in_long_roh'],
        })
    
    load_df = pd.DataFrame(load_data)
    
    logger.info(f"Calculated load for {len(load_df)} individuals")
    
    return load_df

# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

# Setup visualization style - clean white background without grid background color
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3
plt.rcParams['grid.color'] = '#cccccc'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8

# Colorblind-friendly palette (Wong palette)
# Reference: https://www.nature.com/articles/nmeth.1618
CB_COLORS = {
    'blue': '#0072B2',        # Strong blue
    'orange': '#E69F00',      # Orange
    'sky_blue': '#56B4E9',    # Sky blue
    'green': '#009E73',       # Bluish green
    'yellow': '#F0E442',      # Yellow
    'vermillion': '#D55E00',  # Vermillion (red-ish)
    'purple': '#CC79A7',      # Reddish purple
    'black': '#000000',       # Black
}

# Set colorblind-friendly palette
sns.set_palette([CB_COLORS['blue'], CB_COLORS['orange'], CB_COLORS['green'], 
                 CB_COLORS['vermillion'], CB_COLORS['purple'], CB_COLORS['sky_blue']])

def plot_load_distributions(df, output_dir):
    """Plot distributions of deleterious variant counts"""
    logger.info("  Creating variant count distributions...")

    metrics = [
        ('LOF_Het', 'LOF Heterozygous'),
        ('LOF_Hom', 'LOF Homozygous'),
        ('DelMis_Het', 'Deleterious Missense Heterozygous'),
        ('DelMis_Hom', 'Deleterious Missense Homozygous'),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    for ax, (col, title) in zip(axes, metrics):
        vals = df[col].fillna(0)
        ax.hist(vals, bins=40, color=CB_COLORS['blue'], alpha=0.7, edgecolor='white')
        ax.set_xlabel(title, fontsize=11, fontweight='bold')
        ax.set_ylabel('Count', fontsize=11, fontweight='bold')
        ax.set_title(f'Distribution of {title}', fontsize=13, fontweight='bold')
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out = f"{output_dir}/individual_load_distributions.png"
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"    Saved: {out}")


def plot_total_deleterious(df, output_dir):
    """Plot total deleterious distribution and top individuals"""
    logger.info("  Creating total deleterious overview...")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Histogram
    ax = axes[0]
    vals = df['Total_Deleterious'].fillna(0)
    ax.hist(vals, bins=40, color=CB_COLORS['purple'], alpha=0.7, edgecolor='white')
    ax.set_xlabel('Total Deleterious Variants', fontsize=11, fontweight='bold')
    ax.set_ylabel('Count', fontsize=11, fontweight='bold')
    ax.set_title('Distribution of Total Deleterious Variants', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)

    # Top 20 individuals
    ax = axes[1]
    top_df = df.sort_values('Total_Deleterious', ascending=False).head(20)
    ax.barh(top_df['IID'], top_df['Total_Deleterious'], color=CB_COLORS['blue'], alpha=0.8, edgecolor='white')
    ax.set_xlabel('Total Deleterious', fontsize=11, fontweight='bold')
    ax.set_ylabel('Individual (IID)', fontsize=11, fontweight='bold')
    ax.set_title('Top 20 Individuals by Total Deleterious Variants', fontsize=13, fontweight='bold')
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3, axis='x')

    plt.tight_layout()
    out = f"{output_dir}/total_deleterious_overview.png"
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"    Saved: {out}")


def plot_genetic_load(df, output_dir):
    """Plot genetic load distributions (total, realized, potential)"""
    logger.info("  Creating genetic load distributions...")

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()
    
    entries = [
        ('Total_Genetic_Load', 'Total Genetic Load'),
        ('Realized_Load', 'Realized Load (Expressed)'),
        ('Potential_Load', 'Potential Load (Hidden)'),
        ('Hom_Realized_Load', 'Homozygous Realized Load'),
        ('Het_Realized_Load', 'Heterozygous Realized Load'),
        ('Total_Deleterious', 'Total Deleterious Variant Count'),
    ]

    for ax, (col, title) in zip(axes, entries):
        if col not in df.columns:
            continue
        vals = df[col].fillna(0.0)
        ax.hist(vals, bins=40, color=CB_COLORS['blue'], alpha=0.7, edgecolor='white')
        ax.set_xlabel(title, fontsize=11, fontweight='bold')
        ax.set_ylabel('Count', fontsize=11, fontweight='bold')
        ax.set_title(f'Distribution of {title}', fontsize=13, fontweight='bold')
        ax.axvline(vals.mean(), color=CB_COLORS['vermillion'], linestyle='--', linewidth=2, 
                   label=f'Mean: {vals.mean():.4f}')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out = f"{output_dir}/genetic_load_distributions.png"
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"    Saved: {out}")


def plot_roh_load_analysis(df, output_dir):
    """Plot ROH-specific genetic load patterns"""
    logger.info("  Creating ROH load analysis plots...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    
    # Colorblind-friendly colors for ROH categories (sequential: blue → orange → vermillion)
    roh_colors = [CB_COLORS['sky_blue'], CB_COLORS['orange'], CB_COLORS['vermillion']]
    
    # 1. Load by ROH length category
    ax = axes[0]
    roh_cols = ['Load_in_Short_ROH', 'Load_in_Medium_ROH', 'Load_in_Long_ROH']
    labels = ['Short\n(<1 Mb)\nAncient', 'Medium\n(1-5 Mb)\nIntermediate', 'Long\n(>5 Mb)\nRecent']
    
    if all(col in df.columns for col in roh_cols):
        data = [df[col].dropna() for col in roh_cols]
        bp = ax.boxplot(data, labels=labels, patch_artist=True)
        for patch, color in zip(bp['boxes'], roh_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        ax.set_ylabel('Genetic Load', fontsize=11, fontweight='bold')
        ax.set_title('Genetic Load by ROH Length Category', fontsize=13, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
    
    # 2. Variants in ROH by category
    ax = axes[1]
    var_cols = ['Variants_in_Short_ROH', 'Variants_in_Medium_ROH', 'Variants_in_Long_ROH']
    if all(col in df.columns for col in var_cols):
        means = [df[col].mean() for col in var_cols]
        bars = ax.bar(labels, means, color=roh_colors, alpha=0.7, edgecolor='white')
        ax.set_ylabel('Number of Variants', fontsize=11, fontweight='bold')
        ax.set_title('Deleterious Variants by ROH Category', fontsize=13, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        for bar, mean in zip(bars, means):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                    f'{mean:.1f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    # 3. Load in ROH vs outside ROH
    ax = axes[2]
    if 'Load_in_ROH' in df.columns and 'Load_outside_ROH' in df.columns:
        data = [df['Load_in_ROH'].dropna(), df['Load_outside_ROH'].dropna()]
        bp = ax.boxplot(data, labels=['In ROH', 'Outside ROH'], patch_artist=True)
        bp['boxes'][0].set_facecolor(CB_COLORS['vermillion'])
        bp['boxes'][1].set_facecolor(CB_COLORS['green'])
        for box in bp['boxes']:
            box.set_alpha(0.7)
        ax.set_ylabel('Genetic Load', fontsize=11, fontweight='bold')
        ax.set_title('Load: ROH vs Non-ROH Regions', fontsize=13, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
    
    # 4. Correlation: Realized load vs Long ROH load (recent inbreeding)
    ax = axes[3]
    if 'Realized_Load' in df.columns and 'Load_in_Long_ROH' in df.columns:
        sub = df[['Realized_Load', 'Load_in_Long_ROH']].dropna()
        ax.scatter(sub['Load_in_Long_ROH'], sub['Realized_Load'], 
                   alpha=0.6, s=30, color=CB_COLORS['vermillion'])
        ax.set_xlabel('Load in Long ROH (Recent Inbreeding)', fontsize=11, fontweight='bold')
        ax.set_ylabel('Realized Load', fontsize=11, fontweight='bold')
        ax.set_title('Recent Inbreeding vs Realized Load', fontsize=13, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # Add trend line
        try:
            from scipy import stats
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                sub['Load_in_Long_ROH'], sub['Realized_Load'])
            x_line = np.linspace(sub['Load_in_Long_ROH'].min(), 
                                sub['Load_in_Long_ROH'].max(), 100)
            y_line = slope * x_line + intercept
            ax.plot(x_line, y_line, color=CB_COLORS['blue'], linestyle='--', linewidth=2, 
                   label=f'r={r_value:.3f}, p={p_value:.2e}')
            ax.legend(fontsize=9)
        except Exception as e:
            logger.warning(f"    Could not add regression line: {e}")
    
    plt.tight_layout()
    out = f"{output_dir}/roh_load_analysis.png"
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"    Saved: {out}")


def plot_froh_relationships(df, output_dir):
    """Plot F_ROH relationships with genetic load"""
    if 'F_ROH' not in df.columns:
        logger.info("  F_ROH not found; skipping F_ROH relationship plots")
        return

    logger.info("  Creating F_ROH relationship plots...")

    fig, axes = plt.subplots(3, 2, figsize=(14, 18))
    axes = axes.flatten()

    # F_ROH vs Total_Hom
    ax = axes[0]
    sub = df[['F_ROH', 'Total_Hom']].dropna()
    ax.scatter(sub['F_ROH'], sub['Total_Hom'], alpha=0.6, s=30, color=CB_COLORS['blue'])
    ax.set_xlabel(r'F$_{ROH}$', fontsize=11, fontweight='bold')
    ax.set_ylabel('Total Homozygous Deleterious', fontsize=11, fontweight='bold')
    ax.set_title(r'F$_{ROH}$ vs Total Homozygous Deleterious', fontsize=13, fontweight='bold')
    try:
        from scipy import stats
        r1, p1 = stats.pearsonr(sub['F_ROH'], sub['Total_Hom'])
        ax.text(0.02, 0.95, f'r = {r1:.3f}\np = {p1:.2e}', transform=ax.transAxes,
                va='top', ha='left', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='gray'))
    except Exception as e:
        logger.warning(f"    Could not compute correlation (FROH vs Total_Hom): {e}")
    ax.grid(True, alpha=0.3)

    # F_ROH vs Total_Genetic_Load
    ax = axes[1]
    if 'Total_Genetic_Load' in df.columns:
        sub = df[['F_ROH', 'Total_Genetic_Load']].dropna()
        ax.scatter(sub['F_ROH'], sub['Total_Genetic_Load'], alpha=0.6, s=30, color=CB_COLORS['orange'])
        ax.set_xlabel(r'F$_{ROH}$', fontsize=11, fontweight='bold')
        ax.set_ylabel('Total Genetic Load', fontsize=11, fontweight='bold')
        ax.set_title(r'F$_{ROH}$ vs Total Genetic Load', fontsize=13, fontweight='bold')
        try:
            from scipy import stats
            r2, p2 = stats.pearsonr(sub['F_ROH'], sub['Total_Genetic_Load'])
            ax.text(0.02, 0.95, f'r = {r2:.3f}\np = {p2:.2e}', transform=ax.transAxes,
                    va='top', ha='left', fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='gray'))
        except Exception as e:
            logger.warning(f"    Could not compute correlation (FROH vs Total_Genetic_Load): {e}")
        ax.grid(True, alpha=0.3)

    # F_ROH vs Load_in_ROH
    ax = axes[2]
    if 'Load_in_ROH' in df.columns:
        sub = df[['F_ROH', 'Load_in_ROH']].dropna()
        ax.scatter(sub['F_ROH'], sub['Load_in_ROH'], alpha=0.6, s=30, color=CB_COLORS['purple'])
        ax.set_xlabel(r'F$_{ROH}$', fontsize=11, fontweight='bold')
        ax.set_ylabel('Load in ROH', fontsize=11, fontweight='bold')
        ax.set_title(r'F$_{ROH}$ vs Load from Variants in ROH', fontsize=13, fontweight='bold')
        try:
            from scipy import stats
            r5, p5 = stats.pearsonr(sub['F_ROH'], sub['Load_in_ROH'])
            ax.text(0.02, 0.95, f'r = {r5:.3f}\np = {p5:.2e}', transform=ax.transAxes,
                    va='top', ha='left', fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='gray'))
        except Exception as e:
            logger.warning(f"    Could not compute correlation (FROH vs Load_in_ROH): {e}")
        ax.grid(True, alpha=0.3)

    # F_ROH vs Realized_Load
    ax = axes[3]
    if 'Realized_Load' in df.columns:
        sub = df[['F_ROH', 'Realized_Load']].dropna()
        ax.scatter(sub['F_ROH'], sub['Realized_Load'], alpha=0.6, s=30, color=CB_COLORS['vermillion'])
        ax.set_xlabel(r'F$_{ROH}$', fontsize=11, fontweight='bold')
        ax.set_ylabel('Realized Load (Expressed)', fontsize=11, fontweight='bold')
        ax.set_title(r'F$_{ROH}$ vs Realized Load', fontsize=13, fontweight='bold')
        try:
            from scipy import stats
            r3, p3 = stats.pearsonr(sub['F_ROH'], sub['Realized_Load'])
            ax.text(0.02, 0.95, f'r = {r3:.3f}\np = {p3:.2e}', transform=ax.transAxes,
                    va='top', ha='left', fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='gray'))
        except Exception as e:
            logger.warning(f"    Could not compute correlation (FROH vs Realized_Load): {e}")
        ax.grid(True, alpha=0.3)

    # F_ROH vs Potential_Load
    ax = axes[4]
    if 'Potential_Load' in df.columns:
        sub = df[['F_ROH', 'Potential_Load']].dropna()
        ax.scatter(sub['F_ROH'], sub['Potential_Load'], alpha=0.6, s=30, color=CB_COLORS['green'])
        ax.set_xlabel(r'F$_{ROH}$', fontsize=11, fontweight='bold')
        ax.set_ylabel('Potential Load (Hidden)', fontsize=11, fontweight='bold')
        ax.set_title(r'F$_{ROH}$ vs Potential Load', fontsize=13, fontweight='bold')
        try:
            from scipy import stats
            r4, p4 = stats.pearsonr(sub['F_ROH'], sub['Potential_Load'])
            ax.text(0.02, 0.95, f'r = {r4:.3f}\np = {p4:.2e}', transform=ax.transAxes,
                    va='top', ha='left', fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='gray'))
        except Exception as e:
            logger.warning(f"    Could not compute correlation (FROH vs Potential_Load): {e}")
        ax.grid(True, alpha=0.3)

    # Hide the 6th subplot (bottom right) if not needed, or add another comparison
    ax = axes[5]
    ax.axis('off')  # Hide the last subplot

    plt.tight_layout()
    out = f"{output_dir}/froh_relationships.png"
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"    Saved: {out}")


def plot_correlation_heatmap(df, output_dir):
    """Plot correlation heatmap of load metrics"""
    logger.info("  Creating correlation heatmap...")

    cols = [
        'LOF_Het', 'LOF_Hom', 'DelMis_Het', 'DelMis_Hom',
        'Total_Deleterious', 'Total_Genetic_Load', 
        'Realized_Load', 'Potential_Load', 
        'Hom_Realized_Load', 'Het_Realized_Load', 'Load_in_ROH'
    ]
    if 'F_ROH' in df.columns:
        cols.append('F_ROH')
    
    cols = [c for c in cols if c in df.columns]
    corr = df[cols].corr()

    plt.figure(figsize=(12, 10))
    # Use colorblind-friendly diverging colormap (cividis or RdBu are good alternatives)
    sns.heatmap(corr, annot=True, fmt='.2f', cmap='cividis', square=True,
                cbar_kws={'shrink': 0.8}, linewidths=0.5, linecolor='white')
    plt.title('Correlation Between Load Metrics', fontsize=13, fontweight='bold')
    plt.tight_layout()
    out = f"{output_dir}/load_metrics_correlation.png"
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"    Saved: {out}")


def plot_load_components(df, output_dir):
    """Plot genetic load components breakdown"""
    logger.info("  Creating load components plot...")
    
    if 'Total_Genetic_Load' not in df.columns:
        return
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # 1. Realized vs Potential load scatter
    ax = axes[0]
    if 'Realized_Load' in df.columns and 'Potential_Load' in df.columns:
        sub = df[['Realized_Load', 'Potential_Load']].dropna()
        ax.scatter(sub['Potential_Load'], sub['Realized_Load'], 
                   alpha=0.6, s=30, color=CB_COLORS['green'])
        ax.set_xlabel('Potential Load (Hidden)', fontsize=11, fontweight='bold')
        ax.set_ylabel('Realized Load (Expressed)', fontsize=11, fontweight='bold')
        ax.set_title('Realized vs Potential Genetic Load', fontsize=13, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # Add correlation
        try:
            from scipy import stats
            r, p = stats.pearsonr(sub['Potential_Load'], sub['Realized_Load'])
            ax.text(0.02, 0.95, f'r = {r:.3f}\np = {p:.2e}', transform=ax.transAxes,
                    va='top', ha='left', fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='gray'))
        except Exception as e:
            logger.warning(f"    Could not compute correlation: {e}")
    
    # 2. Stacked bar: Total = Realized + Potential
    ax = axes[1]
    if all(c in df.columns for c in ['Realized_Load', 'Potential_Load']):
        categories = ['Mean Load']
        realized_mean = df['Realized_Load'].mean()
        potential_mean = df['Potential_Load'].mean()
        
        ax.bar(categories, [realized_mean], label='Realized Load', color=CB_COLORS['vermillion'], alpha=0.8)
        ax.bar(categories, [potential_mean], bottom=[realized_mean], label='Potential Load', 
               color=CB_COLORS['sky_blue'], alpha=0.8)
        ax.set_ylabel('Genetic Load', fontsize=11, fontweight='bold')
        ax.set_title('Total = Realized + Potential', fontsize=13, fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add text annotations
        ax.text(0, realized_mean/2, f'{realized_mean:.4f}', ha='center', va='center', fontsize=10, fontweight='bold')
        ax.text(0, realized_mean + potential_mean/2, f'{potential_mean:.4f}', ha='center', va='center', fontsize=10, fontweight='bold')
    
    # 3. Stacked bar: Realized = Hom + Het
    ax = axes[2]
    if all(c in df.columns for c in ['Hom_Realized_Load', 'Het_Realized_Load']):
        categories = ['Mean Realized Load']
        hom_mean = df['Hom_Realized_Load'].mean()
        het_mean = df['Het_Realized_Load'].mean()
        
        ax.bar(categories, [hom_mean], label='Homozygous Realized', color=CB_COLORS['orange'], alpha=0.8)
        ax.bar(categories, [het_mean], bottom=[hom_mean], label='Heterozygous Realized', 
               color=CB_COLORS['green'], alpha=0.8)
        ax.set_ylabel('Realized Load', fontsize=11, fontweight='bold')
        ax.set_title('Realized = Hom + Het', fontsize=13, fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add text annotations
        ax.text(0, hom_mean/2, f'{hom_mean:.4f}', ha='center', va='center', fontsize=10, fontweight='bold')
        ax.text(0, hom_mean + het_mean/2, f'{het_mean:.4f}', ha='center', va='center', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    out = f"{output_dir}/load_components.png"
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"    Saved: {out}")


def generate_visualizations(df, output_dir):
    """Generate all visualizations"""
    logger.info("")
    logger.info("="*70)
    logger.info("GENERATING VISUALIZATIONS")
    logger.info("="*70)
    logger.info("")
    
    # Create visualization directory
    vis_dir = f"{output_dir}/visualizations"
    Path(vis_dir).mkdir(parents=True, exist_ok=True)
    
    # Generate all plots
    try:
        plot_load_distributions(df, vis_dir)
        plot_total_deleterious(df, vis_dir)
        plot_genetic_load(df, vis_dir)
        plot_load_components(df, vis_dir)
        plot_roh_load_analysis(df, vis_dir)
        plot_froh_relationships(df, vis_dir)
        plot_correlation_heatmap(df, vis_dir)
        
        logger.info("")
        logger.info(f"All visualizations saved to: {vis_dir}")
        
    except Exception as e:
        logger.error(f"Error generating visualizations: {e}")
        import traceback
        logger.error(traceback.format_exc())

# ============================================================================
# MAIN
# ============================================================================

def main():
    logger.info("="*70)
    logger.info("INDIVIDUAL GENETIC LOAD CALCULATION")
    logger.info("="*70)
    logger.info("")
    
    # Load data
    plm_scores = load_plm_predictions()
    roh_by_ind = load_roh_data()
    logger.info("")
    
    # Calculate load
    load_df = calculate_individual_load(ANNOTATED_VCF, plm_scores, roh_by_ind)
    
    # Merge with F_ROH if available
    try:
        froh_df = pd.read_csv(FROH_DATA)
        load_df = load_df.merge(froh_df[['IID', 'F_ROH', 'Num_ROH']], on='IID', how='left')
        logger.info("Merged with F_ROH data")
    except:
        logger.warning("Could not merge F_ROH data")
    
    # Save results
    output_file = f"{OUTPUT_DIR}/individual_genetic_load.csv"
    load_df.to_csv(output_file, index=False)
    logger.info(f"\nSaved: {output_file}")
    
    # Summary statistics
    logger.info("")
    logger.info("="*70)
    logger.info("SUMMARY STATISTICS")
    logger.info("="*70)
    
    logger.info(f"\nSigmoid pathogenicity calibration:")
    logger.info(f"  P(Pathogenic) = 1 / (1 + exp({SIGMOID_K} * (x + {-SIGMOID_X0})))")
    logger.info(f"  Midpoint (P=0.5): ESM-2 LLR = {SIGMOID_X0:.4f}")
    logger.info(f"  LOF pathogenicity: {LOF_PATHOGENICITY}")
    logger.info(f"  Dominance (h): {H_DEFAULT}")
    
    logger.info(f"\nMean deleterious variant counts per individual:")
    logger.info(f"  LOF heterozygous: {load_df['LOF_Het'].mean():.2f} ± {load_df['LOF_Het'].std():.2f}")
    logger.info(f"  LOF homozygous: {load_df['LOF_Hom'].mean():.2f} ± {load_df['LOF_Hom'].std():.2f}")
    logger.info(f"  Deleterious missense het: {load_df['DelMis_Het'].mean():.2f} ± {load_df['DelMis_Het'].std():.2f}")
    logger.info(f"  Deleterious missense hom: {load_df['DelMis_Hom'].mean():.2f} ± {load_df['DelMis_Hom'].std():.2f}")
    logger.info(f"  Total deleterious: {load_df['Total_Deleterious'].mean():.2f} ± {load_df['Total_Deleterious'].std():.2f}")
    
    logger.info(f"\nGenetic load metrics (sum of pathogenicity probabilities):")
    logger.info(f"  Total Genetic Load = Realized Load + Potential Load")
    logger.info(f"  Realized Load = Homozygous Realized + Heterozygous Realized")
    logger.info(f"")
    logger.info(f"  Total Genetic Load: {load_df['Total_Genetic_Load'].mean():.4f} ± {load_df['Total_Genetic_Load'].std():.4f}")
    logger.info(f"    Range: {load_df['Total_Genetic_Load'].min():.4f} to {load_df['Total_Genetic_Load'].max():.4f}")
    logger.info(f"  Realized Load (expressed): {load_df['Realized_Load'].mean():.4f} ± {load_df['Realized_Load'].std():.4f}")
    logger.info(f"  Potential Load (hidden): {load_df['Potential_Load'].mean():.4f} ± {load_df['Potential_Load'].std():.4f}")
    logger.info(f"  Homozygous Realized Load: {load_df['Hom_Realized_Load'].mean():.4f} ± {load_df['Hom_Realized_Load'].std():.4f}")
    logger.info(f"  Heterozygous Realized Load: {load_df['Het_Realized_Load'].mean():.4f} ± {load_df['Het_Realized_Load'].std():.4f}")
    
    # Verify relationships
    total_check = load_df['Realized_Load'] + load_df['Potential_Load']
    realized_check = load_df['Hom_Realized_Load'] + load_df['Het_Realized_Load']
    logger.info(f"\n  Verification:")
    logger.info(f"    Realized + Potential ≈ Total: {np.allclose(total_check, load_df['Total_Genetic_Load'])}")
    logger.info(f"    Hom + Het ≈ Realized: {np.allclose(realized_check, load_df['Realized_Load'])}")
    
    logger.info(f"\nROH-associated genetic load (homozygous realized only):")
    logger.info(f"  Load in ROH: {load_df['Load_in_ROH'].mean():.4f} ± {load_df['Load_in_ROH'].std():.4f}")
    logger.info(f"  Load outside ROH: {load_df['Load_outside_ROH'].mean():.4f} ± {load_df['Load_outside_ROH'].std():.4f}")
    logger.info(f"  Load in Short ROH (<1 Mb, ancient): {load_df['Load_in_Short_ROH'].mean():.4f} ± {load_df['Load_in_Short_ROH'].std():.4f}")
    logger.info(f"  Load in Medium ROH (1-5 Mb): {load_df['Load_in_Medium_ROH'].mean():.4f} ± {load_df['Load_in_Medium_ROH'].std():.4f}")
    logger.info(f"  Load in Long ROH (>5 Mb, recent): {load_df['Load_in_Long_ROH'].mean():.4f} ± {load_df['Load_in_Long_ROH'].std():.4f}")
    
    logger.info(f"\nVariants in ROH by length category:")
    logger.info(f"  Short ROH (<1 Mb): {load_df['Variants_in_Short_ROH'].mean():.1f} ± {load_df['Variants_in_Short_ROH'].std():.1f}")
    logger.info(f"  Medium ROH (1-5 Mb): {load_df['Variants_in_Medium_ROH'].mean():.1f} ± {load_df['Variants_in_Medium_ROH'].std():.1f}")
    logger.info(f"  Long ROH (>5 Mb): {load_df['Variants_in_Long_ROH'].mean():.1f} ± {load_df['Variants_in_Long_ROH'].std():.1f}")
    
    if 'F_ROH' in load_df.columns:
        # Correlation analysis
        from scipy import stats
        valid_idx = load_df['F_ROH'].notna()
        r1, p1 = stats.pearsonr(load_df.loc[valid_idx, 'F_ROH'], load_df.loc[valid_idx, 'Total_Hom'])
        r2, p2 = stats.pearsonr(load_df.loc[valid_idx, 'F_ROH'], load_df.loc[valid_idx, 'Total_Genetic_Load'])
        r3, p3 = stats.pearsonr(load_df.loc[valid_idx, 'F_ROH'], load_df.loc[valid_idx, 'Realized_Load'])
        r4, p4 = stats.pearsonr(load_df.loc[valid_idx, 'F_ROH'], load_df.loc[valid_idx, 'Potential_Load'])
        r5, p5 = stats.pearsonr(load_df.loc[valid_idx, 'F_ROH'], load_df.loc[valid_idx, 'Hom_Realized_Load'])
        logger.info(f"\nCorrelations with F_ROH:")
        logger.info(f"  F_ROH vs Total_Hom (count): r = {r1:.3f}, p = {p1:.3e}")
        logger.info(f"  F_ROH vs Total_Genetic_Load: r = {r2:.3f}, p = {p2:.3e}")
        logger.info(f"  F_ROH vs Realized_Load: r = {r3:.3f}, p = {p3:.3e}")
        logger.info(f"  F_ROH vs Potential_Load: r = {r4:.3f}, p = {p4:.3e}")
        logger.info(f"  F_ROH vs Hom_Realized_Load: r = {r5:.3f}, p = {p5:.3e}")
    
    logger.info("")
    logger.info("="*70)
    
    # Generate visualizations
    generate_visualizations(load_df, OUTPUT_DIR)
    
    logger.info("")
    logger.info("="*70)
    logger.info("ANALYSIS COMPLETE")
    logger.info("="*70)
    logger.info(f"\nOutputs:")
    logger.info(f"  Data: {OUTPUT_DIR}/individual_genetic_load.csv")
    logger.info(f"  Visualizations: {OUTPUT_DIR}/visualizations/")
    logger.info(f"  Log: {OUTPUT_DIR}/individual_load.log")
    logger.info("")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())


