#!/usr/bin/env python3
"""
Genotype-phenotype association analysis for the golden snub-nosed monkey.

Integrates field morbidity records with genomic data (ROH, genetic load)
to test associations between inbreeding, realized load, and health scores (CHS).
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import spearmanr, pearsonr, mannwhitneyu, kruskal
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import warnings
import os
import shutil
import subprocess
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir
from datetime import datetime

warnings.filterwarnings('ignore')

# Function to format variable names for display
def format_var_name(var_name):
    """Convert variable names to readable format without underscores"""
    # Mapping for common variable names
    name_mapping = {
        'F_ROH': 'F-ROH',
        'F_ROH_SHORT': 'F-ROH Short',
        'F_ROH_MEDIUM': 'F-ROH Medium',
        'F_ROH_LONG': 'F-ROH Long',
        'Total_Genetic_Load': 'Total Genetic Load',
        'Realized_Load': 'Realized Load',
        'Potential_Load': 'Potential Load',
        'Masked_Load': 'Masked Load',
        'Load_in_ROH': 'Load in ROH',
        'Load_in_Short_ROH': 'Load in Short ROH',
        'Load_in_Medium_ROH': 'Load in Medium ROH',
        'Load_in_Long_ROH': 'Load in Long ROH',
        'Total_Deleterious': 'Total Deleterious',
        'Total_Hom': 'Total Homozygous',
        'Total_Het': 'Total Heterozygous',
        'LOF_Het': 'LOF Heterozygous',
        'LOF_Hom': 'LOF Homozygous',
        'DelMis_Het': 'Deleterious Missense Het',
        'DelMis_Hom': 'Deleterious Missense Hom',
        'Hom_Realized_Load': 'Homozygous Realized Load',
        'Het_Realized_Load': 'Heterozygous Realized Load',
        'N_ROH': 'Number of ROH',
        'KB_MEDIUM_ROH': 'KB Medium ROH',
        'KB_LONG_ROH': 'KB Long ROH',
        'Has_Disease': 'Has Disease',
        'Num_Disease_Types': 'Number of Disease Types',
    }
    
    # Return mapped name or replace underscores with spaces
    return name_mapping.get(var_name, var_name.replace('_', ' '))

# Set plot style - clean white background
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

# Paths
BASE_DIR = get_base_dir()
OUTPUT_DIR = BASE_DIR / "output" / "phenotype_genotype_analysis"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def _resolve_tool(env_name: str, *candidates: str) -> str:
    env = os.environ.get(env_name, "").strip()
    if env:
        return env
    for name in candidates:
        found = shutil.which(name)
        if found:
            return found
    return candidates[0]


# =============================================================================
# 1. PHENOTYPE SCORING SYSTEM
# =============================================================================

def load_and_process_disease_data():
    """
    Load disease records and build the Composite Health Score (CHS).
    """
    print("=" * 70)
    print("STEP 1: Loading and Processing Disease Records")
    print("=" * 70)
    
    disease_file = BASE_DIR / "data" / "disease_record_noredundancy.csv"
    disease_df = pd.read_csv(disease_file)
    
    print(f"\nTotal disease records: {len(disease_df)}")
    print(f"Unique individuals with diseases: {disease_df['ID'].nunique()}")
    print(f"\nDisease types distribution:")
    print(disease_df['Disease Type'].value_counts())
    
    # Severity weights for CHS (higher = more severe fitness impact)
    disease_weights = {
        'Eye Disease': 1.0,           # Vision critical for survival
        'Respiratory Disease': 1.0,   # Respiratory function essential
        'Skin Disease': 0.8,          # Less immediately life-threatening
        'Finger Joint Abnormality': 0.9  # Affects mobility/dexterity
    }
    
    # Define condition severity scores
    condition_severity = {
        # Eye conditions
        'Cataract': 1.5,
        'Injury/Disability': 2.0,     # Most severe - permanent damage
        'Tissue Hyperplasia': 1.0,
        
        # Respiratory conditions
        'Tracheal Infection': 2.0,
        
        # Skin conditions
        'Skin Infection': 1.5,
        'Hair Loss': 1.0,
        
        # Joint conditions
        'Flexion/Extension Abnormality': 1.5
    }
    
    # Calculate severity score for each record
    disease_df['Severity_Score'] = disease_df.apply(
        lambda x: disease_weights.get(x['Disease Type'], 0.5) * 
                  condition_severity.get(x['Suspected Condition'], 1.0),
        axis=1
    )
    
    # Aggregate by individual
    individual_disease = disease_df.groupby('ID').agg({
        'Disease Type': lambda x: list(x.unique()),
        'Suspected Condition': lambda x: list(x.unique()),
        'Severity_Score': 'sum',
        'Age': 'first'  # Age at first disease discovery
    }).reset_index()
    
    individual_disease['Num_Disease_Types'] = individual_disease['Disease Type'].apply(len)
    individual_disease['Num_Conditions'] = individual_disease['Suspected Condition'].apply(len)
    
    # Create binary disease indicators
    for disease_type in disease_df['Disease Type'].unique():
        individual_disease[f'Has_{disease_type.replace(" ", "_")}'] = \
            individual_disease['Disease Type'].apply(lambda x: 1 if disease_type in x else 0)
    
    # Calculate Composite Health Score (CHS)
    # Normalize to 0-10 scale where 10 is worst
    max_severity = individual_disease['Severity_Score'].max()
    individual_disease['CHS'] = (individual_disease['Severity_Score'] / max_severity) * 10
    
    print(f"\nIndividual disease summary:")
    print(individual_disease[['ID', 'Num_Disease_Types', 'Severity_Score', 'CHS']].to_string())
    
    return disease_df, individual_disease


def create_full_phenotype_dataset(individual_disease):
    """
    Create a complete phenotype dataset for all 68 individuals,
    with healthy individuals having CHS = 0.
    """
    # Load FAM file to get all individuals
    fam_df = pd.read_csv(
        BASE_DIR / "data" / "monkey_snp_sex_qc.fam",
        sep='\t',
        header=None,
        names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO']
    )
    
    # Sex mapping (1=male, 2=female in PLINK)
    fam_df['Sex'] = fam_df['SEX'].map({1: 'Male', 2: 'Female'})
    
    # Merge with disease data
    phenotype_df = fam_df[['IID', 'Sex']].merge(
        individual_disease[['ID', 'Num_Disease_Types', 'Num_Conditions', 
                            'Severity_Score', 'CHS', 'Age',
                            'Has_Eye_Disease', 'Has_Respiratory_Disease', 
                            'Has_Skin_Disease', 'Has_Finger_Joint_Abnormality']],
        left_on='IID', right_on='ID', how='left'
    )
    
    # Fill NaN for healthy individuals
    phenotype_df = phenotype_df.drop(columns=['ID'], errors='ignore')
    phenotype_df = phenotype_df.fillna({
        'Num_Disease_Types': 0,
        'Num_Conditions': 0,
        'Severity_Score': 0,
        'CHS': 0,
        'Has_Eye_Disease': 0,
        'Has_Respiratory_Disease': 0,
        'Has_Skin_Disease': 0,
        'Has_Finger_Joint_Abnormality': 0
    })
    
    # Create binary disease status
    phenotype_df['Has_Disease'] = (phenotype_df['CHS'] > 0).astype(int)
    
    print(f"\nFull phenotype dataset: {len(phenotype_df)} individuals")
    print(f"Individuals with disease: {phenotype_df['Has_Disease'].sum()}")
    print(f"Healthy individuals: {(phenotype_df['Has_Disease'] == 0).sum()}")
    
    return phenotype_df


# =============================================================================
# 2. MERGE WITH GENOMIC DATA
# =============================================================================

def merge_with_genomic_data(phenotype_df):
    """
    Merge phenotype data with ROH and genetic load data.
    """
    print("\n" + "=" * 70)
    print("STEP 2: Merging Phenotype with Genomic Data")
    print("=" * 70)
    
    # Load genetic load data (already contains F_ROH from ROH analysis)
    genetic_load = pd.read_csv(
        BASE_DIR / "output" / "phase5_genetic_load" / "individual_load" / "individual_genetic_load.csv"
    )
    
    # Load ROH data for additional ROH metrics
    roh_stats = pd.read_csv(
        BASE_DIR / "output" / "phase3a_population_genomics" / "roh_analysis" / "roh_individual_stats.csv",
        sep='\t'
    )
    
    # Merge phenotype with genetic load (which already has F_ROH)
    merged_df = phenotype_df.merge(genetic_load, on='IID', how='left')
    
    # Only merge ROH metrics not already in genetic_load (avoid duplicate F_ROH)
    roh_cols_to_add = ['IID', 'N_ROH', 'TOTAL_KB', 'N_SHORT_ROH', 'KB_SHORT_ROH', 
                       'N_MEDIUM_ROH', 'KB_MEDIUM_ROH', 'N_LONG_ROH', 'KB_LONG_ROH',
                       'F_ROH_SHORT', 'F_ROH_MEDIUM', 'F_ROH_LONG']
    # Remove columns that are already in merged_df (except IID for merge)
    roh_cols_to_add = [c for c in roh_cols_to_add if c not in merged_df.columns or c == 'IID']
    
    merged_df = merged_df.merge(
        roh_stats[roh_cols_to_add],
        on='IID', how='left'
    )
    
    print(f"\nMerged dataset: {len(merged_df)} individuals with {merged_df.shape[1]} variables")
    
    # Key variables for analysis
    # Updated for V4.1: Total_Genetic_Load, Realized_Load, Potential_Load, Hom/Het Realized Load
    key_vars = ['CHS', 'Has_Disease', 'F_ROH', 'Total_Genetic_Load', 'Realized_Load', 'Potential_Load', 
                'Hom_Realized_Load', 'Het_Realized_Load', 'Load_in_ROH', 
                'Total_Deleterious', 'LOF_Hom', 'DelMis_Hom']
    
    print("\nKey variable summary (all individuals):")
    print(merged_df[key_vars].describe().round(4))
    
    # Save merged dataset
    merged_df.to_csv(OUTPUT_DIR / "merged_phenotype_genotype.csv", index=False)
    
    return merged_df


# =============================================================================
# 3. GENOTYPE-PHENOTYPE ASSOCIATION ANALYSIS
# =============================================================================

def correlation_analysis(merged_df):
    """
    Perform correlation analysis between genomic metrics and health phenotypes.
    """
    print("\n" + "=" * 70)
    print("STEP 3: Correlation Analysis")
    print("=" * 70)
    
    # Genomic variables of interest
    # Note: Using sigmoid pathogenicity method V4.1 with h=0.25
    # Column names: Total_Genetic_Load, Realized_Load, Potential_Load, Hom/Het Realized Load
    genomic_vars = [
        'F_ROH', 'F_ROH_SHORT', 'F_ROH_MEDIUM', 'F_ROH_LONG',
        'Total_Genetic_Load', 'Realized_Load', 'Potential_Load', 
        'Hom_Realized_Load', 'Het_Realized_Load', 'Load_in_ROH',
        'Total_Deleterious', 'Total_Hom', 'LOF_Hom', 'DelMis_Hom'
    ]
    
    # Phenotype variables
    phenotype_vars = ['CHS', 'Num_Disease_Types', 'Severity_Score']
    
    # Calculate correlations
    results = []
    
    for gvar in genomic_vars:
        for pvar in phenotype_vars:
            if gvar in merged_df.columns and pvar in merged_df.columns:
                # Spearman correlation (non-parametric)
                rho, p_spearman = spearmanr(merged_df[gvar], merged_df[pvar])
                
                # Pearson correlation
                r, p_pearson = pearsonr(merged_df[gvar], merged_df[pvar])
                
                results.append({
                    'Genomic_Variable': gvar,
                    'Phenotype_Variable': pvar,
                    'Spearman_rho': rho,
                    'Spearman_p': p_spearman,
                    'Pearson_r': r,
                    'Pearson_p': p_pearson
                })
    
    corr_df = pd.DataFrame(results)
    
    # Multiple testing correction
    corr_df['Spearman_p_adj'] = multipletests(corr_df['Spearman_p'], method='fdr_bh')[1]
    corr_df['Pearson_p_adj'] = multipletests(corr_df['Pearson_p'], method='fdr_bh')[1]
    
    # Add significance indicators
    def get_sig_indicator(p_val):
        """Return significance indicator based on p-value"""
        if p_val < 0.001:
            return '***'
        elif p_val < 0.01:
            return '**'
        elif p_val < 0.05:
            return '*'
        else:
            return ''
    
    corr_df['Spearman_sig'] = corr_df['Spearman_p'].apply(get_sig_indicator)
    corr_df['Spearman_sig_adj'] = corr_df['Spearman_p_adj'].apply(get_sig_indicator)
    corr_df['Pearson_sig'] = corr_df['Pearson_p'].apply(get_sig_indicator)
    corr_df['Pearson_sig_adj'] = corr_df['Pearson_p_adj'].apply(get_sig_indicator)
    
    # Sort by significance
    corr_df = corr_df.sort_values('Spearman_p')
    
    print("\nCorrelation Results (sorted by p-value):")
    print("Significance levels: *** p<0.001, ** p<0.01, * p<0.05")
    print(corr_df.to_string(index=False))
    
    # Print significant correlations (uncorrected p < 0.05)
    sig_uncorrected = corr_df[corr_df['Spearman_p'] < 0.05]
    if len(sig_uncorrected) > 0:
        print(f"\n{len(sig_uncorrected)} significant correlations (uncorrected p < 0.05):")
        for _, row in sig_uncorrected.iterrows():
            print(f"  {row['Genomic_Variable']} vs {row['Phenotype_Variable']}: "
                  f"ρ={row['Spearman_rho']:.3f}, p={row['Spearman_p']:.4f}{row['Spearman_sig']}")
    
    # Print FDR-corrected significant correlations
    sig_corrected = corr_df[corr_df['Spearman_p_adj'] < 0.05]
    if len(sig_corrected) > 0:
        print(f"\n{len(sig_corrected)} significant correlations (FDR-corrected p < 0.05):")
        for _, row in sig_corrected.iterrows():
            print(f"  {row['Genomic_Variable']} vs {row['Phenotype_Variable']}: "
                  f"ρ={row['Spearman_rho']:.3f}, p_adj={row['Spearman_p_adj']:.4f}{row['Spearman_sig_adj']}")
    else:
        print("\nNo significant correlations after FDR correction (p_adj < 0.05)")
    
    # Save results
    corr_df.to_csv(OUTPUT_DIR / "correlation_results.csv", index=False)
    
    return corr_df


def case_control_analysis(merged_df):
    """
    Compare genomic metrics between diseased and healthy individuals.
    """
    print("\n" + "=" * 70)
    print("STEP 4: Case-Control Analysis (Diseased vs Healthy)")
    print("=" * 70)
    
    diseased = merged_df[merged_df['Has_Disease'] == 1]
    healthy = merged_df[merged_df['Has_Disease'] == 0]
    
    print(f"\nDiseased individuals: {len(diseased)}")
    print(f"Healthy individuals: {len(healthy)}")
    
    # Variables to compare
    # Note: Using sigmoid pathogenicity method V4.1 with h=0.25
    compare_vars = [
        'F_ROH', 'Total_Genetic_Load', 'Realized_Load', 'Potential_Load', 
        'Hom_Realized_Load', 'Het_Realized_Load', 'Load_in_ROH',
        'Total_Deleterious', 'Total_Hom', 'LOF_Hom', 'DelMis_Hom',
        'N_ROH', 'KB_MEDIUM_ROH', 'KB_LONG_ROH'
    ]
    
    results = []
    
    for var in compare_vars:
        if var in merged_df.columns:
            d_vals = diseased[var].dropna()
            h_vals = healthy[var].dropna()
            
            # Mann-Whitney U test
            stat, p = mannwhitneyu(d_vals, h_vals, alternative='two-sided')
            
            # Effect size (Cohen's d)
            pooled_std = np.sqrt((d_vals.std()**2 + h_vals.std()**2) / 2)
            cohens_d = (d_vals.mean() - h_vals.mean()) / pooled_std if pooled_std > 0 else 0
            
            results.append({
                'Variable': var,
                'Diseased_Mean': d_vals.mean(),
                'Diseased_SD': d_vals.std(),
                'Healthy_Mean': h_vals.mean(),
                'Healthy_SD': h_vals.std(),
                'Difference': d_vals.mean() - h_vals.mean(),
                'Cohens_d': cohens_d,
                'Mann_Whitney_U': stat,
                'P_value': p
            })
    
    cc_df = pd.DataFrame(results)
    cc_df['P_adj'] = multipletests(cc_df['P_value'], method='fdr_bh')[1]
    cc_df = cc_df.sort_values('P_value')
    
    print("\nCase-Control Comparison:")
    print(cc_df.to_string(index=False))
    
    # Save results
    cc_df.to_csv(OUTPUT_DIR / "case_control_analysis.csv", index=False)
    
    return cc_df


def regression_analysis(merged_df):
    """
    Perform regression analysis to test genotype-phenotype associations.
    """
    print("\n" + "=" * 70)
    print("STEP 5: Regression Analysis (GLM Models)")
    print("=" * 70)
    
    # Prepare data for regression
    reg_df = merged_df.copy()
    reg_df['Sex_numeric'] = reg_df['Sex'].map({'Male': 1, 'Female': 0})
    
    # Model 1: CHS ~ Realized_Load + Potential_Load + F_ROH + Sex
    print("\n--- Model 1: Continuous CHS ---")
    try:
        model1 = smf.ols(
            'CHS ~ Realized_Load + Potential_Load + F_ROH + Sex_numeric',
            data=reg_df
        ).fit()
        print(model1.summary())
        
        # Save model summary
        with open(OUTPUT_DIR / "model1_continuous_CHS.txt", 'w') as f:
            f.write(str(model1.summary()))
    except Exception as e:
        print(f"Model 1 error: {e}")
    
    # Model 2: Has_Disease (binary) ~ genomic variables (Logistic Regression)
    print("\n--- Model 2: Binary Disease Status (Logistic Regression) ---")
    try:
        model2 = smf.logit(
            'Has_Disease ~ Realized_Load + Load_in_ROH + F_ROH + Sex_numeric',
            data=reg_df
        ).fit(disp=0)
        print(model2.summary())
        
        # Calculate odds ratios
        print("\nOdds Ratios (exp(coef)):")
        odds_ratios = np.exp(model2.params)
        print(odds_ratios)
        
        with open(OUTPUT_DIR / "model2_logistic_disease.txt", 'w') as f:
            f.write(str(model2.summary()))
            f.write("\n\nOdds Ratios:\n")
            f.write(str(odds_ratios))
    except Exception as e:
        print(f"Model 2 error: {e}")
    
    # Model 3: CHS ~ Load_in_ROH (ROH-associated load as key predictor)
    print("\n--- Model 3: CHS ~ ROH-Associated Load (RDML) ---")
    try:
        model3 = smf.ols(
            'CHS ~ Load_in_ROH + Sex_numeric',
            data=reg_df
        ).fit()
        print(model3.summary())
        
        with open(OUTPUT_DIR / "model3_ROH_load.txt", 'w') as f:
            f.write(str(model3.summary()))
    except Exception as e:
        print(f"Model 3 error: {e}")
    
    # Model 4: CHS ~ Hom_Realized_Load + Het_Realized_Load (decompose realized load)
    print("\n--- Model 4: CHS ~ Hom + Het Realized Load ---")
    try:
        model4 = smf.ols(
            'CHS ~ Hom_Realized_Load + Het_Realized_Load + Sex_numeric',
            data=reg_df
        ).fit()
        print(model4.summary())
        
        with open(OUTPUT_DIR / "model4_hom_het_load.txt", 'w') as f:
            f.write(str(model4.summary()))
    except Exception as e:
        print(f"Model 4 error: {e}")
    
    # Model comparison using AIC
    print("\n--- Model Comparison ---")
    models_info = []
    
    # Null model
    try:
        null_model = smf.ols('CHS ~ 1', data=reg_df).fit()
        models_info.append({'Model': 'Null', 'AIC': null_model.aic, 'R2': null_model.rsquared})
    except:
        pass
    
    # Demographics only
    try:
        demo_model = smf.ols('CHS ~ Sex_numeric', data=reg_df).fit()
        models_info.append({'Model': 'Demographics', 'AIC': demo_model.aic, 'R2': demo_model.rsquared})
    except:
        pass
    
    # F_ROH only
    try:
        froh_model = smf.ols('CHS ~ F_ROH', data=reg_df).fit()
        models_info.append({'Model': 'F_ROH only', 'AIC': froh_model.aic, 'R2': froh_model.rsquared})
    except:
        pass
    
    # Genetic load only
    try:
        load_model = smf.ols('CHS ~ Realized_Load', data=reg_df).fit()
        models_info.append({'Model': 'Realized_Load only', 'AIC': load_model.aic, 'R2': load_model.rsquared})
    except:
        pass
    
    # ROH-associated load only
    try:
        roh_load_model = smf.ols('CHS ~ Load_in_ROH', data=reg_df).fit()
        models_info.append({'Model': 'ROH_Load only', 'AIC': roh_load_model.aic, 'R2': roh_load_model.rsquared})
    except:
        pass
    
    # Homozygous Realized Load only
    try:
        hom_load_model = smf.ols('CHS ~ Hom_Realized_Load', data=reg_df).fit()
        models_info.append({'Model': 'Hom_Realized_Load only', 'AIC': hom_load_model.aic, 'R2': hom_load_model.rsquared})
    except:
        pass
    
    # Heterozygous Realized Load only
    try:
        het_load_model = smf.ols('CHS ~ Het_Realized_Load', data=reg_df).fit()
        models_info.append({'Model': 'Het_Realized_Load only', 'AIC': het_load_model.aic, 'R2': het_load_model.rsquared})
    except:
        pass
    
    # Hom + Het Realized Load
    try:
        hom_het_model = smf.ols('CHS ~ Hom_Realized_Load + Het_Realized_Load', data=reg_df).fit()
        models_info.append({'Model': 'Hom+Het Realized', 'AIC': hom_het_model.aic, 'R2': hom_het_model.rsquared})
    except:
        pass
    
    # Full model
    try:
        full_model = smf.ols('CHS ~ Realized_Load + Load_in_ROH + F_ROH + Sex_numeric', data=reg_df).fit()
        models_info.append({'Model': 'Full', 'AIC': full_model.aic, 'R2': full_model.rsquared})
    except:
        pass
    
    # Full model with Hom/Het components
    try:
        full_hom_het_model = smf.ols('CHS ~ Hom_Realized_Load + Het_Realized_Load + Load_in_ROH + F_ROH + Sex_numeric', data=reg_df).fit()
        models_info.append({'Model': 'Full (Hom+Het)', 'AIC': full_hom_het_model.aic, 'R2': full_hom_het_model.rsquared})
    except:
        pass
    
    model_comp = pd.DataFrame(models_info)
    model_comp = model_comp.sort_values('AIC')
    print("\nModel Comparison (AIC - lower is better):")
    print(model_comp.to_string(index=False))
    
    model_comp.to_csv(OUTPUT_DIR / "model_comparison.csv", index=False)
    
    return model_comp


# =============================================================================
# 4. HYPOTHESIS TESTING
# =============================================================================

def test_hypotheses(merged_df, corr_df, cc_df):
    """
    Test two competing hypotheses:
    A: Inbreeding load (recent inbreeding associated with health problems)
    B: Historical bottleneck (fixed deleterious alleles from past demography)
    """
    print("\n" + "=" * 70)
    print("STEP 6: Hypothesis Testing")
    print("=" * 70)
    
    hypothesis_results = {
        'Hypothesis A (Inbreeding Load)': {},
        'Hypothesis B (Historical Bottleneck)': {}
    }
    
    # Test 1: ROH length distribution
    print("\n--- Test 1: ROH Length Distribution ---")
    
    # Proportion of long ROH (medium + long) vs short ROH - compute before subsetting
    merged_df['Long_ROH_Ratio'] = (merged_df['KB_MEDIUM_ROH'] + merged_df['KB_LONG_ROH']) / merged_df['TOTAL_KB']
    
    diseased = merged_df[merged_df['Has_Disease'] == 1]
    healthy = merged_df[merged_df['Has_Disease'] == 0]
    
    d_long_ratio = diseased['Long_ROH_Ratio'].mean() if len(diseased) > 0 else 0
    h_long_ratio = healthy['Long_ROH_Ratio'].mean() if len(healthy) > 0 else 0
    
    print(f"Mean long ROH ratio - Diseased: {d_long_ratio:.4f}, Healthy: {h_long_ratio:.4f}")
    
    # Hypothesis A predicts: diseased have more long ROH (recent inbreeding)
    # Hypothesis B predicts: no difference (all from ancient events)
    
    # Test 2: F_ROH variance
    print("\n--- Test 2: F_ROH Variance (CV) ---")
    froh_cv = merged_df['F_ROH'].std() / merged_df['F_ROH'].mean()
    print(f"F_ROH Coefficient of Variation: {froh_cv:.4f}")
    
    # Hypothesis A predicts: CV > 0.3 (high variance = differential inbreeding)
    # Hypothesis B predicts: CV < 0.15 (low variance = shared historical load)
    
    if froh_cv > 0.3:
        hypothesis_results['Hypothesis A (Inbreeding Load)']['F_ROH_CV'] = 'Supported (CV > 0.3)'
    elif froh_cv < 0.15:
        hypothesis_results['Hypothesis B (Historical Bottleneck)']['F_ROH_CV'] = 'Supported (CV < 0.15)'
    else:
        hypothesis_results['Hypothesis A (Inbreeding Load)']['F_ROH_CV'] = 'Intermediate'
        hypothesis_results['Hypothesis B (Historical Bottleneck)']['F_ROH_CV'] = 'Intermediate'
    
    # Test 3: RDML-health correlation
    print("\n--- Test 3: ROH-Associated Load vs Health ---")
    rho_rdml, p_rdml = spearmanr(merged_df['Load_in_ROH'], merged_df['CHS'])
    print(f"Spearman correlation (Load_in_ROH vs CHS): rho={rho_rdml:.4f}, p={p_rdml:.4f}")
    
    # Hypothesis A predicts: significant positive correlation (RDML drives health problems)
    # Hypothesis B predicts: weak or no correlation (health problems from fixed alleles)
    
    if p_rdml < 0.05 and rho_rdml > 0:
        hypothesis_results['Hypothesis A (Inbreeding Load)']['RDML_CHS'] = f'Supported (rho={rho_rdml:.3f}, p={p_rdml:.4f})'
    else:
        hypothesis_results['Hypothesis B (Historical Bottleneck)']['RDML_CHS'] = f'Supported (rho={rho_rdml:.3f}, p={p_rdml:.4f})'
    
    # Test 4: Individual CHS variance
    print("\n--- Test 4: CHS Variance ---")
    chs_cv = merged_df['CHS'].std() / merged_df['CHS'].mean() if merged_df['CHS'].mean() > 0 else np.nan
    print(f"CHS Coefficient of Variation: {chs_cv:.4f}")
    
    # Test 5: Realized load comparison
    print("\n--- Test 5: Realized Load Comparison ---")
    if len(cc_df) > 0:
        rl_result = cc_df[cc_df['Variable'] == 'Realized_Load']
        if len(rl_result) > 0:
            rl_diff = rl_result['Difference'].values[0]
            rl_p = rl_result['P_value'].values[0]
            print(f"Realized Load difference (Diseased - Healthy): {rl_diff:.4f}, p={rl_p:.4f}")
            
            if rl_p < 0.05 and rl_diff > 0:
                hypothesis_results['Hypothesis A (Inbreeding Load)']['Realized_Load'] = 'Supported'
    
    # Summary
    print("\n" + "=" * 70)
    print("HYPOTHESIS TESTING SUMMARY")
    print("=" * 70)
    
    for hyp, results in hypothesis_results.items():
        print(f"\n{hyp}:")
        for test, outcome in results.items():
            print(f"  - {test}: {outcome}")
    
    # Overall assessment
    print("\n--- Overall Assessment ---")
    
    # Check predominant ROH pattern
    mean_short_roh_prop = merged_df['F_ROH_SHORT'].mean() / merged_df['F_ROH'].mean()
    mean_long_roh_prop = (merged_df['F_ROH_MEDIUM'].mean() + merged_df['F_ROH_LONG'].mean()) / merged_df['F_ROH'].mean()
    
    print(f"Short ROH proportion: {mean_short_roh_prop:.2%}")
    print(f"Medium+Long ROH proportion: {mean_long_roh_prop:.2%}")
    
    if mean_short_roh_prop > 0.8:
        print("\n=> ROH pattern suggests ANCIENT inbreeding/bottleneck events dominate")
        print("   (Short ROH >80% indicates fragmentation of old autozygous segments)")
    
    if mean_long_roh_prop > 0.2:
        print("\n=> Significant medium/long ROH suggests ONGOING restricted gene flow")
        print("   contributing to inbreeding load")
    
    # Save hypothesis testing results
    with open(OUTPUT_DIR / "hypothesis_testing_results.txt", 'w') as f:
        f.write("HYPOTHESIS TESTING RESULTS\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write(f"F_ROH CV: {froh_cv:.4f}\n")
        f.write(f"RDML-CHS correlation: rho={rho_rdml:.4f}, p={p_rdml:.4f}\n")
        f.write(f"Short ROH proportion: {mean_short_roh_prop:.2%}\n")
        f.write(f"Medium+Long ROH proportion: {mean_long_roh_prop:.2%}\n\n")
        
        for hyp, results in hypothesis_results.items():
            f.write(f"\n{hyp}:\n")
            for test, outcome in results.items():
                f.write(f"  - {test}: {outcome}\n")
    
    return hypothesis_results


# =============================================================================
# 5. VISUALIZATIONS
# =============================================================================

def create_visualizations(merged_df, corr_df, cc_df):
    """
    Create publication-quality visualizations.
    """
    print("\n" + "=" * 70)
    print("STEP 7: Creating Visualizations")
    print("=" * 70)
    
    fig_dir = OUTPUT_DIR / "figures"
    fig_dir.mkdir(exist_ok=True)
    
    # Figure 1: Disease status by genomic metrics
    # Use colorblind-friendly colors: blue for healthy, vermillion for diseased
    cb_palette = [CB_COLORS['blue'], CB_COLORS['vermillion']]
    
    fig, axes = plt.subplots(3, 3, figsize=(15, 14))
    
    # Updated: Using sigmoid method V4.1 with h=0.25 - including Hom/Het Realized Load
    vars_to_plot = ['F_ROH', 'Total_Genetic_Load', 'Realized_Load', 
                    'Hom_Realized_Load', 'Het_Realized_Load', 'Potential_Load',
                    'Load_in_ROH', 'Total_Hom', 'LOF_Hom']
    titles = ['F-ROH (Inbreeding Coefficient)', 'Total Genetic Load', 'Realized Load',
              'Homozygous Realized Load', 'Heterozygous Realized Load', 'Potential Load',
              'ROH-Associated Load', 'Total Homozygous Variants', 'Homozygous LOF Variants']
    
    for ax, var, title in zip(axes.flat, vars_to_plot, titles):
        if var in merged_df.columns:
            sns.boxplot(data=merged_df, x='Has_Disease', y=var, ax=ax, palette=cb_palette)
            ax.set_xlabel('Disease Status (0=Healthy, 1=Diseased)')
            ax.set_ylabel(title)
            ax.set_title(title)
            
            # Add statistical annotation
            diseased = merged_df[merged_df['Has_Disease'] == 1][var]
            healthy = merged_df[merged_df['Has_Disease'] == 0][var]
            if len(diseased) > 0 and len(healthy) > 0:
                _, p = mannwhitneyu(diseased, healthy)
                ax.text(0.5, 0.95, f'p = {p:.4f}', transform=ax.transAxes, 
                        ha='center', fontsize=10, 
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    plt.tight_layout()
    plt.savefig(fig_dir / "01_disease_vs_genomic_metrics.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Figure 2: Correlation heatmap with significance indicators
    fig, ax = plt.subplots(figsize=(14, 11))
    
    # Updated: Using sigmoid method V4.1 with h=0.25 - including Hom/Het Realized Load
    genomic_vars = ['F_ROH', 'F_ROH_SHORT', 'F_ROH_MEDIUM', 'Total_Genetic_Load',
                    'Realized_Load', 'Potential_Load', 'Hom_Realized_Load', 'Het_Realized_Load',
                    'Load_in_ROH', 'Total_Hom']
    phenotype_vars = ['CHS', 'Has_Disease', 'Num_Disease_Types']
    
    all_vars = genomic_vars + phenotype_vars
    available_vars = [v for v in all_vars if v in merged_df.columns]
    
    # Calculate correlation matrix
    corr_matrix = merged_df[available_vars].corr(method='spearman')
    
    # Calculate p-values for each correlation
    n = len(merged_df)
    p_matrix = np.zeros_like(corr_matrix, dtype=float)
    
    for i, var1 in enumerate(available_vars):
        for j, var2 in enumerate(available_vars):
            if i != j:
                _, p = spearmanr(merged_df[var1], merged_df[var2])
                p_matrix[i, j] = p
            else:
                p_matrix[i, j] = 1.0
    
    # Create annotation matrix with significance indicators
    annot_matrix = np.empty_like(corr_matrix, dtype=object)
    for i in range(corr_matrix.shape[0]):
        for j in range(corr_matrix.shape[1]):
            r_val = corr_matrix.iloc[i, j]
            p_val = p_matrix[i, j]
            
            # Add significance stars
            if p_val < 0.001:
                sig = '***'
            elif p_val < 0.01:
                sig = '**'
            elif p_val < 0.05:
                sig = '*'
            else:
                sig = ''
            
            annot_matrix[i, j] = f'{r_val:.2f}{sig}'
    
    # Format row and column labels to remove underscores
    formatted_labels = [format_var_name(v) for v in corr_matrix.columns]
    corr_matrix.columns = formatted_labels
    corr_matrix.index = formatted_labels
    
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool), k=1)
    
    # Use colorblind-friendly colormap (cividis)
    sns.heatmap(corr_matrix, mask=mask, annot=annot_matrix, cmap='cividis', center=0,
                square=True, linewidths=0.5, fmt='', ax=ax, 
                cbar_kws={'label': 'Spearman Correlation'})
    ax.set_title('Correlation Matrix: Genomic Metrics vs Health Phenotypes\n(* p<0.05, ** p<0.01, *** p<0.001)', 
                 fontsize=14, pad=20)
    
    plt.tight_layout()
    plt.savefig(fig_dir / "02_correlation_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Figure 3: Genetic load components vs CHS scatter plots (colorblind-friendly)
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    # Use colorblind-friendly colors: blue for healthy, vermillion for diseased
    colors = [CB_COLORS['blue'] if d == 0 else CB_COLORS['vermillion'] for d in merged_df['Has_Disease']]
    
    # Helper function to get significance stars
    def get_sig_stars(p):
        if p < 0.001:
            return '***'
        elif p < 0.01:
            return '**'
        elif p < 0.05:
            return '*'
        else:
            return 'ns'
    
    # Panel A: F-ROH vs CHS
    ax = axes[0]
    ax.scatter(merged_df['F_ROH'], merged_df['CHS'], c=colors, alpha=0.7, s=60, edgecolors='white', linewidth=0.5)
    slope, intercept, r, p, se = stats.linregress(merged_df['F_ROH'], merged_df['CHS'])
    x_line = np.linspace(merged_df['F_ROH'].min(), merged_df['F_ROH'].max(), 100)
    ax.plot(x_line, slope * x_line + intercept, 'k--', alpha=0.5)
    ax.set_xlabel('F-ROH (Inbreeding Coefficient)', fontsize=11)
    ax.set_ylabel('Composite Health Score (CHS)', fontsize=11)
    sig_stars = get_sig_stars(p)
    ax.set_title(f'F-ROH vs CHS\n(r={r:.3f}, p={p:.4f} {sig_stars})', fontsize=11)
    ax.grid(True, alpha=0.3)
    
    # Panel B: Total Genetic Load vs CHS
    ax = axes[1]
    if 'Total_Genetic_Load' in merged_df.columns:
        ax.scatter(merged_df['Total_Genetic_Load'], merged_df['CHS'], c=colors, alpha=0.7, s=60, edgecolors='white', linewidth=0.5)
        slope, intercept, r, p, se = stats.linregress(merged_df['Total_Genetic_Load'], merged_df['CHS'])
        x_line = np.linspace(merged_df['Total_Genetic_Load'].min(), merged_df['Total_Genetic_Load'].max(), 100)
        ax.plot(x_line, slope * x_line + intercept, 'k--', alpha=0.5)
        ax.set_xlabel('Total Genetic Load', fontsize=11)
        ax.set_ylabel('Composite Health Score (CHS)', fontsize=11)
        sig_stars = get_sig_stars(p)
        ax.set_title(f'Total Genetic Load vs CHS\n(r={r:.3f}, p={p:.4f} {sig_stars})', fontsize=11)
        ax.grid(True, alpha=0.3)
    
    # Panel C: Realized Load vs CHS
    ax = axes[2]
    ax.scatter(merged_df['Realized_Load'], merged_df['CHS'], c=colors, alpha=0.7, s=60, edgecolors='white', linewidth=0.5)
    slope, intercept, r, p, se = stats.linregress(merged_df['Realized_Load'], merged_df['CHS'])
    x_line = np.linspace(merged_df['Realized_Load'].min(), merged_df['Realized_Load'].max(), 100)
    ax.plot(x_line, slope * x_line + intercept, 'k--', alpha=0.5)
    ax.set_xlabel('Realized Load', fontsize=11)
    ax.set_ylabel('Composite Health Score (CHS)', fontsize=11)
    sig_stars = get_sig_stars(p)
    ax.set_title(f'Realized Load vs CHS\n(r={r:.3f}, p={p:.4f} {sig_stars})', fontsize=11)
    ax.grid(True, alpha=0.3)
    
    # Panel D: Homozygous Realized Load vs CHS
    ax = axes[3]
    if 'Hom_Realized_Load' in merged_df.columns:
        ax.scatter(merged_df['Hom_Realized_Load'], merged_df['CHS'], c=colors, alpha=0.7, s=60, edgecolors='white', linewidth=0.5)
        slope, intercept, r, p, se = stats.linregress(merged_df['Hom_Realized_Load'], merged_df['CHS'])
        x_line = np.linspace(merged_df['Hom_Realized_Load'].min(), merged_df['Hom_Realized_Load'].max(), 100)
        ax.plot(x_line, slope * x_line + intercept, 'k--', alpha=0.5)
        ax.set_xlabel('Homozygous Realized Load', fontsize=11)
        ax.set_ylabel('Composite Health Score (CHS)', fontsize=11)
        sig_stars = get_sig_stars(p)
        ax.set_title(f'Hom Realized Load vs CHS\n(r={r:.3f}, p={p:.4f} {sig_stars})', fontsize=11)
        ax.grid(True, alpha=0.3)
    
    # Panel E: Heterozygous Realized Load vs CHS
    ax = axes[4]
    if 'Het_Realized_Load' in merged_df.columns:
        ax.scatter(merged_df['Het_Realized_Load'], merged_df['CHS'], c=colors, alpha=0.7, s=60, edgecolors='white', linewidth=0.5)
        slope, intercept, r, p, se = stats.linregress(merged_df['Het_Realized_Load'], merged_df['CHS'])
        x_line = np.linspace(merged_df['Het_Realized_Load'].min(), merged_df['Het_Realized_Load'].max(), 100)
        ax.plot(x_line, slope * x_line + intercept, 'k--', alpha=0.5)
        ax.set_xlabel('Heterozygous Realized Load', fontsize=11)
        ax.set_ylabel('Composite Health Score (CHS)', fontsize=11)
        sig_stars = get_sig_stars(p)
        ax.set_title(f'Het Realized Load vs CHS\n(r={r:.3f}, p={p:.4f} {sig_stars})', fontsize=11)
        ax.grid(True, alpha=0.3)
    
    # Panel F: Load in ROH vs CHS
    ax = axes[5]
    ax.scatter(merged_df['Load_in_ROH'], merged_df['CHS'], c=colors, alpha=0.7, s=60, edgecolors='white', linewidth=0.5)
    slope, intercept, r, p, se = stats.linregress(merged_df['Load_in_ROH'], merged_df['CHS'])
    x_line = np.linspace(merged_df['Load_in_ROH'].min(), merged_df['Load_in_ROH'].max(), 100)
    ax.plot(x_line, slope * x_line + intercept, 'k--', alpha=0.5)
    ax.set_xlabel('ROH-Associated Load (RDML)', fontsize=11)
    ax.set_ylabel('Composite Health Score (CHS)', fontsize=11)
    sig_stars = get_sig_stars(p)
    ax.set_title(f'ROH Load vs CHS\n(r={r:.3f}, p={p:.4f} {sig_stars})', fontsize=11)
    ax.grid(True, alpha=0.3)
    
    # Add legend (colorblind-friendly)
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=CB_COLORS['blue'], label='Healthy'),
                       Patch(facecolor=CB_COLORS['vermillion'], label='Diseased')]
    fig.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(0.98, 0.98))
    
    plt.tight_layout()
    plt.savefig(fig_dir / "03_genotype_phenotype_scatter.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Figure 4: ROH category composition by disease status (colorblind-friendly)
    fig, ax = plt.subplots(figsize=(10, 6))
    
    roh_categories = ['F_ROH_SHORT', 'F_ROH_MEDIUM', 'F_ROH_LONG']
    
    diseased = merged_df[merged_df['Has_Disease'] == 1]
    healthy = merged_df[merged_df['Has_Disease'] == 0]
    
    x = np.arange(len(roh_categories))
    width = 0.35
    
    diseased_means = [diseased[cat].mean() for cat in roh_categories]
    healthy_means = [healthy[cat].mean() for cat in roh_categories]
    diseased_sems = [diseased[cat].sem() for cat in roh_categories]
    healthy_sems = [healthy[cat].sem() for cat in roh_categories]
    
    # Use colorblind-friendly colors: blue for healthy, vermillion for diseased
    bars1 = ax.bar(x - width/2, healthy_means, width, yerr=healthy_sems, 
                   label='Healthy', color=CB_COLORS['blue'], capsize=5, edgecolor='white')
    bars2 = ax.bar(x + width/2, diseased_means, width, yerr=diseased_sems, 
                   label='Diseased', color=CB_COLORS['vermillion'], capsize=5, edgecolor='white')
    
    ax.set_ylabel('Mean F-ROH', fontsize=12)
    ax.set_xlabel('ROH Category', fontsize=12)
    ax.set_title('ROH Composition by Disease Status', fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(['Short (<1 Mb)\n(Ancient)', 'Medium (1-5 Mb)\n(Intermediate)', 'Long (>5 Mb)\n(Recent)'])
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(fig_dir / "04_roh_composition_by_disease.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Figure 5: Case-control effect sizes (colorblind-friendly)
    if len(cc_df) > 0:
        fig, ax = plt.subplots(figsize=(10, 8))
        
        cc_sorted = cc_df.sort_values('Cohens_d')
        # Use colorblind-friendly: vermillion for positive (higher in diseased), blue for negative
        colors = [CB_COLORS['vermillion'] if d > 0 else CB_COLORS['blue'] for d in cc_sorted['Cohens_d']]
        
        y_pos = np.arange(len(cc_sorted))
        ax.barh(y_pos, cc_sorted['Cohens_d'], color=colors, alpha=0.8, edgecolor='white')
        ax.set_yticks(y_pos)
        # Format variable names to remove underscores
        ax.set_yticklabels([format_var_name(v) for v in cc_sorted['Variable']])
        ax.set_xlabel("Cohen's d (Effect Size)", fontsize=12)
        ax.set_title("Effect Sizes: Diseased vs Healthy Individuals\n(Positive = higher in diseased)", fontsize=14)
        ax.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
        
        # Add significance markers
        for i, (idx, row) in enumerate(cc_sorted.iterrows()):
            if row['P_adj'] < 0.05:
                ax.text(row['Cohens_d'], i, ' *', va='center', fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(fig_dir / "05_case_control_effect_sizes.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    # Figure 6: Individual profile heatmap (diseased individuals) - colorblind-friendly
    fig, ax = plt.subplots(figsize=(14, 8))
    
    diseased_ids = merged_df[merged_df['Has_Disease'] == 1]['IID'].tolist()
    # Updated: Using sigmoid method V4.1 with h=0.25
    profile_vars = ['CHS', 'F_ROH', 'Total_Genetic_Load', 'Realized_Load', 'Potential_Load', 'Load_in_ROH', 
                    'Total_Hom', 'LOF_Hom', 'DelMis_Hom']
    # Filter to available columns
    profile_vars = [v for v in profile_vars if v in merged_df.columns]
    
    profile_data = merged_df[merged_df['Has_Disease'] == 1][['IID'] + profile_vars].set_index('IID')
    
    # Z-score normalize
    profile_normalized = (profile_data - profile_data.mean()) / profile_data.std()
    
    # Rename columns to remove underscores
    profile_normalized.columns = [format_var_name(col) for col in profile_normalized.columns]
    
    # Use colorblind-friendly colormap (cividis)
    sns.heatmap(profile_normalized, annot=True, cmap='cividis', center=0,
                fmt='.2f', linewidths=0.5, ax=ax)
    ax.set_title('Genomic Profile of Diseased Individuals (Z-scores)', fontsize=14)
    ax.set_xlabel('Variable')
    ax.set_ylabel('Individual ID')
    
    plt.tight_layout()
    plt.savefig(fig_dir / "06_diseased_individual_profiles.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\nFigures saved to: {fig_dir}")
    
    return fig_dir


def create_disease_type_analysis(merged_df, disease_df):
    """
    Analyze genomic associations for specific disease types.
    """
    print("\n" + "=" * 70)
    print("STEP 8: Disease-Specific Analysis")
    print("=" * 70)
    
    disease_types = ['Eye_Disease', 'Respiratory_Disease', 'Skin_Disease', 'Finger_Joint_Abnormality']
    genomic_vars = ['F_ROH', 'Total_Genetic_Load', 'Realized_Load', 'Hom_Realized_Load', 
                    'Het_Realized_Load', 'Load_in_ROH', 'Total_Hom']
    
    results = []
    
    for dtype in disease_types:
        col_name = f'Has_{dtype}'
        if col_name in merged_df.columns:
            affected = merged_df[merged_df[col_name] == 1]
            unaffected = merged_df[merged_df[col_name] == 0]
            
            print(f"\n{dtype.replace('_', ' ')}:")
            print(f"  Affected: {len(affected)}, Unaffected: {len(unaffected)}")
            
            for gvar in genomic_vars:
                if gvar in merged_df.columns:
                    a_vals = affected[gvar].dropna()
                    u_vals = unaffected[gvar].dropna()
                    
                    if len(a_vals) > 1 and len(u_vals) > 1:
                        stat, p = mannwhitneyu(a_vals, u_vals, alternative='two-sided')
                        
                        results.append({
                            'Disease_Type': dtype.replace('_', ' '),
                            'Genomic_Variable': gvar,
                            'Affected_Mean': a_vals.mean(),
                            'Unaffected_Mean': u_vals.mean(),
                            'Difference': a_vals.mean() - u_vals.mean(),
                            'P_value': p
                        })
    
    disease_specific = pd.DataFrame(results)
    
    if len(disease_specific) > 0:
        disease_specific['P_adj'] = multipletests(disease_specific['P_value'], method='fdr_bh')[1]
        disease_specific = disease_specific.sort_values('P_value')
        
        print("\nDisease-Specific Analysis Results:")
        print(disease_specific.to_string(index=False))
        
        disease_specific.to_csv(OUTPUT_DIR / "disease_specific_analysis.csv", index=False)
    
    return disease_specific


    # =============================================================================
    # 7. GWAS ANALYSIS (GCTA-MLMA / GEMMA LMM, GRM-based relatedness correction)
    # =============================================================================

def build_chrom_map():
    """Build NCBI accession -> chromosome number mapping from GFF region lines."""
    gff_file = BASE_DIR / "data" / "reference" / "Rhinopithecus_roxellana.gff"
    acc_to_chr = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 9 and parts[2] == 'region' and 'chromosome=' in parts[8]:
                acc_to_chr[parts[0]] = parts[8].split('chromosome=')[1].split(';')[0]
            if len(acc_to_chr) >= 25:
                break
    return acc_to_chr


def _make_grm(gwas_dir):
    """
    Compute GRM (genetic relatedness matrix) using plink2.
    Uses MAF>0.01 and LD-pruned SNPs for a stable, unbiased estimate.
    The GRM is shared by both GCTA-MLMA (CHS) and GEMMA LMM (binary traits).
    Returns the GRM prefix path (str) or None on failure.
    """
    bfile = str(BASE_DIR / "data" / "monkey_snp_sex_qc")
    grm_prefix = str(gwas_dir / "grm")
    prune_prefix = str(gwas_dir / "grm_prune")

    # Step 1: LD pruning (window=50, step=5, r2=0.2)
    r = subprocess.run(
        ['plink2', '--bfile', bfile, '--maf', '0.01',
         '--indep-pairwise', '50', '5', '0.2',
         '--out', prune_prefix, '--threads', '16'],
        capture_output=True, text=True)
    prune_in = Path(prune_prefix + '.prune.in')
    if not prune_in.exists():
        print(f"  LD pruning failed: {r.stderr[-200:]}")
        return None
    n_pruned = sum(1 for _ in open(prune_in))
    print(f"  LD-pruned SNPs for GRM: {n_pruned}")

    # Step 2: Compute GRM in binary format
    r = subprocess.run(
        ['plink2', '--bfile', bfile, '--extract', str(prune_in),
         '--make-grm-bin', '--out', grm_prefix, '--threads', '16'],
        capture_output=True, text=True)
    if not Path(grm_prefix + '.grm.bin').exists():
        print(f"  GRM computation failed: {r.stderr[-200:]}")
        return None
    print(f"  GRM written to: {grm_prefix}.grm.bin")
    return grm_prefix


def _make_sex_covar(merged_df, gwas_dir):
    """Write a minimal sex-only covariate file (GCTA format: FID IID covar, no header)."""
    covar = merged_df[['IID']].copy()
    covar.insert(0, 'FID', 0)
    covar['Sex'] = merged_df['Sex'].map({'Male': 1, 'Female': 0})
    path = str(gwas_dir / "covar_sex.txt")
    covar.to_csv(path, sep='\t', index=False, header=False)
    return path


def run_gwas_analysis(merged_df):
    """
    GWAS using linear mixed models (LMM) to jointly correct for population
    structure and cryptic relatedness via a genome-wide GRM.

    CHS (quantitative) -> GCTA-MLMA (--mlma)
      The GRM is used as the random effect; sex is the fixed covariate.
      GCTA-MLMA is the gold standard for quantitative traits in structured
      populations with known relatedness.

    Has_Eye_Disease / Has_Finger_Joint_Abnormality (binary, >=5 cases) -> GEMMA LMM
      GEMMA fits a linear mixed model with a centered relatedness matrix.
      For small binary GWAS (n<100), GEMMA LMM outperforms Firth logistic
      because the GRM random effect absorbs both structure and relatedness
      without requiring explicit PC covariates.

    Has_Respiratory_Disease (2 cases) / Has_Skin_Disease (3 cases) -> skipped.
    """
    print("\n" + "=" * 70)
    print("STEP 10: GWAS Analysis (GCTA-MLMA / GEMMA LMM)")
    print("=" * 70)

    GCTA = _resolve_tool("GCTA_BIN", "gcta64", "gcta")
    GEMMA = _resolve_tool("GEMMA_BIN", "gemma")
    bfile = str(BASE_DIR / "data" / "monkey_snp_sex_qc")

    gwas_dir = OUTPUT_DIR / "gwas"
    gwas_dir.mkdir(exist_ok=True)

    # GRM (shared by GCTA and GEMMA)
    print("\nComputing GRM (LD-pruned, MAF>0.01)...")
    grm_prefix = _make_grm(gwas_dir)
    if grm_prefix is None:
        print("  GRM failed; cannot run LMM-based GWAS")
        return {}

    # Phenotype helpers
    def _write_pheno_gcta(col, path):
        """GCTA phenotype: FID IID pheno (no header, tab-separated)."""
        df = merged_df[['IID', col]].copy()
        df.insert(0, 'FID', 0)
        df.to_csv(path, sep='\t', index=False, header=False)

    fam_orig = pd.read_csv(
        BASE_DIR / "data" / "monkey_snp_sex_qc.fam",
        sep=r'\s+', header=None,
        names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO'])

    def _write_fam_gemma(col, path):
        """
        Write a PLINK .fam with the binary phenotype for GEMMA.
        GEMMA convention: 1=control, 2=case, -9=missing.
        """
        fam = fam_orig.copy()
        pheno_map = merged_df.set_index('IID')[col]
        fam['PHENO'] = fam['IID'].map(pheno_map).apply(
            lambda x: 2 if x == 1 else (1 if x == 0 else -9))
        fam.to_csv(path, sep=' ', index=False, header=False)

    sex_covar = _make_sex_covar(merged_df, gwas_dir)
    gwas_results = {}
    MIN_CASES = 5

    # GCTA-MLMA: CHS (quantitative)
    print("\nRunning GCTA-MLMA for CHS (quantitative, n=68)...")
    pheno_chs = str(gwas_dir / "pheno_chs.txt")
    _write_pheno_gcta('CHS', pheno_chs)
    mlma_prefix = str(gwas_dir / "mlma_chs")
    r = subprocess.run(
        [GCTA, '--mlma', '--bfile', bfile,
         '--grm', grm_prefix,
         '--pheno', pheno_chs,
         '--covar', sex_covar,
         '--out', mlma_prefix,
         '--thread-num', '16'],
        capture_output=True, text=True)
    mlma_file = Path(mlma_prefix + '.mlma')
    if mlma_file.exists():
        df = pd.read_csv(mlma_file, sep='\t')
        df = df.dropna(subset=['p']).copy()
        df = df[df['p'] > 0]
        df = df.rename(columns={'Chr': '#CHROM', 'bp': 'POS', 'SNP': 'ID',
                                 'b': 'BETA', 'se': 'SE', 'p': 'P',
                                 'Freq': 'A1_FREQ'})
        df['#CHROM'] = df['#CHROM'].astype(str)
        from scipy import stats as _stats
        lam = np.median(_stats.chi2.isf(df['P'].values, 1)) / _stats.chi2.ppf(0.5, 1)
        gwas_results['CHS'] = df
        n_snp = len(df)
        p_b = 0.05 / n_snp
        n_bonf = (df['P'] < p_b).sum()
        n_sug = (df['P'] < 1e-5).sum()
        print(f"  SNPs tested: {n_snp:,}; Bonferroni p < {p_b:.2e} -> {n_bonf} hits")
        print(f"  lambda={lam:.3f}; exploratory (p<1e-5): {n_sug}")
        if lam > 1.3:
            print(f"  WARNING: lambda={lam:.3f}, residual inflation present")
        top = df.nsmallest(20, 'P')[
            [c for c in ['#CHROM', 'POS', 'ID', 'A1_FREQ', 'BETA', 'SE', 'P']
             if c in df.columns]]
        print("\n  Top 20 hits (GCTA-MLMA):")
        print(top.to_string(index=False))
        top.to_csv(gwas_dir / "top20_CHS.csv", index=False)
    else:
        print(f"  GCTA-MLMA failed: {r.stderr[-400:]}")

    # GEMMA LMM: binary disease traits
    # GEMMA needs a valid phenotype in the .fam to compute GRM.
    # Write a temporary fileset with Has_Eye_Disease as placeholder phenotype.
    print("\nComputing GEMMA centered relatedness matrix...")
    import shutil as _shutil
    gemma_base_prefix = str(gwas_dir / "tmp_gemma_base")
    fam_base_tmp = gemma_base_prefix + '.fam'
    fam_base = fam_orig.copy()
    pheno_map_eye = merged_df.set_index('IID').get('Has_Eye_Disease',
                                                    pd.Series(dtype=float))
    fam_base['PHENO'] = fam_base['IID'].map(pheno_map_eye).apply(
        lambda x: 2 if x == 1 else (1 if x == 0 else 1))  # default to 1 (ctrl)
    fam_base.to_csv(fam_base_tmp, sep=' ', index=False, header=False)
    for ext, src in [('.bim', BASE_DIR / "data" / "monkey_snp_sex_qc.bim"),
                     ('.bed', BASE_DIR / "data" / "monkey_snp_sex_qc.bed")]:
        dst = Path(gemma_base_prefix + ext)
        if not dst.exists():
            _shutil.copy(str(src), str(dst))

    gemma_grm_name = "gemma_grm"
    r_grm = subprocess.run(
        [GEMMA, '-bfile', gemma_base_prefix, '-gk', '1', '-o', gemma_grm_name],
        capture_output=True, text=True, cwd=str(gwas_dir))
    # GEMMA writes to output/ subdir relative to cwd
    gemma_k = gwas_dir / "output" / f"{gemma_grm_name}.cXX.txt"
    if not gemma_k.exists():
        print(f"  GEMMA GRM failed: {r_grm.stderr[-200:]}")
        gemma_k = None

    for col in ['Has_Eye_Disease', 'Has_Finger_Joint_Abnormality']:
        if col not in merged_df.columns:
            continue
        nc = int(merged_df[col].sum())
        if nc < MIN_CASES:
            print(f"\nSkipping {col}: only {nc} cases (< {MIN_CASES})")
            continue
        print(f"\nRunning GEMMA LMM for {col} (cases={nc}, controls={68-nc})...")
        if gemma_k is None:
            print("  GEMMA GRM unavailable, skipping")
            continue

        # Build a temporary plink fileset with this phenotype in the .fam
        import shutil
        tmp_prefix = str(gwas_dir / f"tmp_{col}")
        fam_tmp = tmp_prefix + '.fam'
        _write_fam_gemma(col, fam_tmp)
        for ext, src in [('.bim', BASE_DIR / "data" / "monkey_snp_sex_qc.bim"),
                         ('.bed', BASE_DIR / "data" / "monkey_snp_sex_qc.bed")]:
            dst = Path(tmp_prefix + ext)
            if not dst.exists():
                shutil.copy(str(src), str(dst))

        out_name = f"gemma_{col.lower()}"
        r = subprocess.run(
            [GEMMA, '-bfile', tmp_prefix,
             '-k', str(gemma_k),
             '-lmm', '4',
             '-o', out_name],
            capture_output=True, text=True, cwd=str(gwas_dir))

        assoc_file = gwas_dir / "output" / f"{out_name}.assoc.txt"
        if not assoc_file.exists():
            assoc_file = gwas_dir / f"{out_name}.assoc.txt"
        if assoc_file.exists():
            df = pd.read_csv(assoc_file, sep='\t')
            p_col = 'p_wald' if 'p_wald' in df.columns else 'p_lrt'
            df = df.dropna(subset=[p_col]).copy()
            df = df[df[p_col] > 0]
            df = df.rename(columns={'chr': '#CHROM', 'ps': 'POS', 'rs': 'ID',
                                     'af': 'A1_FREQ', 'beta': 'BETA',
                                     'se': 'SE', p_col: 'P'})
            df['#CHROM'] = df['#CHROM'].astype(str)
            from scipy import stats as _stats
            lam = np.median(_stats.chi2.isf(df['P'].values, 1)) / _stats.chi2.ppf(0.5, 1)
            gwas_results[col] = df
            n_snp = len(df)
            p_b = 0.05 / n_snp
            n_bonf = (df['P'] < p_b).sum()
            n_sug = (df['P'] < 1e-5).sum()
            print(f"  SNPs tested: {n_snp:,}; Bonferroni p < {p_b:.2e} -> {n_bonf} hits")
            print(f"  lambda={lam:.3f}; exploratory (p<1e-5): {n_sug}")
            if lam > 1.3:
                print(f"  WARNING: lambda={lam:.3f}, residual inflation present")
            top = df.nsmallest(20, 'P')[
                [c for c in ['#CHROM', 'POS', 'ID', 'A1_FREQ', 'BETA', 'SE', 'P']
                 if c in df.columns]]
            print(f"\n  Top 20 hits (GEMMA LMM):")
            print(top.to_string(index=False))
            top.to_csv(gwas_dir / f"top20_{col}.csv", index=False)
        else:
            print(f"  GEMMA failed: {r.stderr[-300:]}")

    for col in ['Has_Respiratory_Disease', 'Has_Skin_Disease']:
        if col in merged_df.columns:
            nc = int(merged_df[col].sum())
            print(f"\nSkipping {col}: only {nc} cases (< {MIN_CASES}), insufficient for GWAS")

    fig_dir = OUTPUT_DIR / "figures"
    for trait, df in gwas_results.items():
        _plot_manhattan_qq(df, trait, fig_dir)

    return gwas_results


def _plot_manhattan_qq(df, trait, fig_dir):
    """
    Publication-quality Manhattan + QQ plots (>=300 DPI raster, optional PDF).

    Manhattan horizontal threshold: Bonferroni alpha=0.05 / n_SNPs (n = rows in df
    with valid chromosome and P), not genome-wide 5e-8. Exploratory 1e-5 shown dotted.
    """
    from scipy import stats as _stats

    trait_title = trait.replace('_', ' ')
    # Publication figure: ~7 in tall at 400 dpi -> >=2800 px height; wide two-panel layout
    fig_w, fig_h = 16.0, 6.25
    dpi_raster = 400

    pub_rc = {
        'font.family': 'sans-serif',
        'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica', 'sans-serif'],
        'font.size': 10.5,
        'axes.titlesize': 12,
        'axes.labelsize': 11.5,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 9,
        'axes.linewidth': 0.9,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'lines.linewidth': 1.2,
        'figure.dpi': dpi_raster,
        'savefig.dpi': dpi_raster,
    }

    d = df.copy()
    d['CHR'] = pd.to_numeric(d['#CHROM'], errors='coerce')
    d = d.dropna(subset=['CHR']).copy()
    d['LP'] = -np.log10(np.clip(d['P'].values, 1e-323, None))

    chroms = sorted(d['CHR'].unique())
    offsets, centers, cum = {}, {}, 0
    for c in chroms:
        offsets[c] = cum
        mx = float(d.loc[d['CHR'] == c, 'POS'].max())
        centers[c] = cum + mx / 2
        cum += mx + 5e6
    d['X'] = d['POS'] + d['CHR'].map(offsets)

    # Bonferroni correction for independent tests: alpha / n_SNPs (no 5e-8 shortcut)
    n_tests = len(d)
    p_bonf = (0.05 / n_tests) if n_tests > 0 else np.nan
    y_bonf = float(-np.log10(p_bonf)) if n_tests > 0 and np.isfinite(p_bonf) else np.nan

    col_a = '#3B6E8F'   # muted blue
    col_b = '#7BA3BC'   # lighter blue-gray
    col_sig = '#C44E52' # distinct red (Wong-inspired)

    with plt.rc_context(pub_rc):
        fig, (ax1, ax2) = plt.subplots(
            1, 2, figsize=(fig_w, fig_h), constrained_layout=False)
        fig.patch.set_facecolor('white')
        fig.subplots_adjust(left=0.06, right=0.98, bottom=0.14, top=0.90, wspace=0.22)

        # --- Manhattan ---
        ax1.set_facecolor('#FAFAFA')
        for spine in ('top', 'right'):
            ax1.spines[spine].set_visible(False)
        ax1.spines['left'].set_color('#444444')
        ax1.spines['bottom'].set_color('#444444')
        ax1.grid(axis='y', linestyle='-', linewidth=0.4, alpha=0.35, color='#999999')
        ax1.set_axisbelow(True)

        for i, c in enumerate(chroms):
            sub = d[d['CHR'] == c]
            if len(sub) > 50000:
                hi = sub[sub['LP'] >= 3]
                lo = sub[sub['LP'] < 3].sample(frac=0.1, random_state=42)
                sub = pd.concat([hi, lo])
            col_c = col_a if i % 2 == 0 else col_b
            ax1.scatter(
                sub['X'], sub['LP'], c=col_c, s=3.2, alpha=0.55,
                linewidths=0, rasterized=True, zorder=1)

        sig = d[d['P'] < p_bonf] if n_tests > 0 else d.iloc[0:0]
        if len(sig) > 0:
            ax1.scatter(
                sig['X'], sig['LP'], c=col_sig, s=18, alpha=0.92,
                linewidths=0.25, edgecolors='#5C1A1C', rasterized=True, zorder=4)

        if n_tests > 0 and np.isfinite(y_bonf):
            ax1.axhline(
                y_bonf, color='#8B0000', ls='--', lw=1.0, alpha=0.75, zorder=2,
                label=f'Bonferroni alpha=0.05 (p < {p_bonf:.2e})')
        y_sug = -np.log10(1e-5)
        ax1.axhline(y_sug, color='#B8860B', ls=':', lw=1.0, alpha=0.65, zorder=2,
                    label=r'Exploratory $10^{-5}$')
        ax1.set_xlabel('Chromosome', labelpad=6)
        ax1.set_ylabel(r'$-\log_{10}(P)$', labelpad=6)
        ax1.set_title(f'Manhattan plot — {trait_title}', fontweight='semibold', pad=10)
        if n_tests > 0:
            ax1.text(
                0.012, 0.98,
                rf'$n$ = {n_tests:,} SNPs' + '\n'
                rf'Bonferroni: $p < {p_bonf:.2e}$',
                transform=ax1.transAxes, fontsize=9, verticalalignment='top',
                horizontalalignment='left',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                          edgecolor='#BBBBBB', linewidth=0.6, alpha=0.95))
        ax1.set_xticks([centers[c] for c in chroms])
        ax1.set_xticklabels([str(int(c)) for c in chroms], fontsize=8, rotation=0)
        leg1 = ax1.legend(
            loc='upper right', frameon=True, fancybox=False, framealpha=0.95,
            edgecolor='#CCCCCC', facecolor='white')
        leg1.get_frame().set_linewidth(0.6)

        ymax = float(np.nanmax(d['LP']))
        if not np.isfinite(ymax):
            ymax = 1.0
        y_top = min(ymax * 1.08 + 0.5, ymax + 3)
        if n_tests > 0 and np.isfinite(y_bonf):
            y_top = max(y_top, y_bonf + 0.8)
        ax1.set_ylim(0, y_top)

        # --- QQ (ascending expected quantiles; 95% band for order stats under null) ---
        pvals = np.clip(d['P'].values.astype(float), 1e-323, 1.0)
        m = len(pvals)
        lam = float(
            np.median(_stats.chi2.isf(d['P'].values, df=1))
            / _stats.chi2.ppf(0.5, df=1))
        # Standard GWAS QQ: x = expected -log10(P), y = observed, both increasing
        obs_asc = np.sort(-np.log10(pvals))
        j = np.arange(1, m + 1, dtype=float)
        exp_asc = -np.log10((j - 0.5) / (m + 0.5))

        ax2.set_facecolor('#FAFAFA')
        for spine in ('top', 'right'):
            ax2.spines[spine].set_visible(False)
        ax2.spines['left'].set_color('#444444')
        ax2.spines['bottom'].set_color('#444444')
        ax2.grid(linestyle='-', linewidth=0.4, alpha=0.35, color='#999999')
        ax2.set_axisbelow(True)

        lo_p = np.clip(_stats.beta.ppf(0.025, j, m - j + 1), 1e-300, 1.0)
        hi_p = np.clip(_stats.beta.ppf(0.975, j, m - j + 1), 1e-300, 1.0)
        band_lo = -np.log10(hi_p)
        band_hi = -np.log10(lo_p)
        ax2.fill_between(
            exp_asc, band_lo, band_hi, color='#D0D0D0', alpha=0.55,
            linewidth=0, zorder=0, label='95% null envelope')

        if m > 50000:
            rng = np.random.default_rng(42)
            idx = np.sort(np.concatenate([
                np.arange(min(5000, m)),
                rng.choice(range(5000, m), min(45000, m - 5000), replace=False)
            ]))
            ax2.scatter(
                exp_asc[idx], obs_asc[idx], s=4, alpha=0.45, c=col_a,
                linewidths=0, rasterized=True, zorder=2)
        else:
            ax2.scatter(
                exp_asc, obs_asc, s=4, alpha=0.45, c=col_a,
                linewidths=0, rasterized=True, zorder=2)

        mx = max(float(exp_asc.max()), float(obs_asc.max())) * 1.04
        ax2.plot([0, mx], [0, mx], color='#333333', ls='-', lw=1.15, alpha=0.9, zorder=3,
                 label='Null (y = x)')
        ax2.set_xlim(0, mx)
        ax2.set_ylim(0, mx)
        ax2.set_aspect('equal', adjustable='box')
        ax2.set_xlabel(r'Expected $-\log_{10}(P)$', labelpad=6)
        ax2.set_ylabel(r'Observed $-\log_{10}(P)$', labelpad=6)
        ax2.set_title(f'QQ plot — {trait_title}', fontweight='semibold', pad=10)
        p_bonf_qq = (0.05 / m) if m > 0 else np.nan
        lam_txt = (
            r'$\lambda_{\mathrm{GC}}$' + f' = {lam:.3f}\n'
            r'$n$' + f' = {m:,} SNPs\n'
            r'Bonferroni $p$' + f' < {p_bonf_qq:.2e}'
        )
        ax2.text(
            0.03, 0.97, lam_txt, transform=ax2.transAxes, fontsize=9.5,
            verticalalignment='top', horizontalalignment='left',
            bbox=dict(boxstyle='round,pad=0.35', facecolor='white',
                      edgecolor='#BBBBBB', linewidth=0.7, alpha=0.96),
            family='sans-serif')
        leg2 = ax2.legend(
            loc='lower right', frameon=True, fancybox=False, framealpha=0.95,
            edgecolor='#CCCCCC', facecolor='white')
        leg2.get_frame().set_linewidth(0.6)

        base = f"07_gwas_{trait}"
        out_png = fig_dir / f"{base}.png"
        out_pdf = fig_dir / f"{base}.pdf"
        plt.savefig(
            out_png, dpi=dpi_raster, bbox_inches='tight', pad_inches=0.04,
            facecolor='white', edgecolor='none', format='png')
        plt.savefig(
            out_pdf, bbox_inches='tight', pad_inches=0.04,
            facecolor='white', edgecolor='none', format='pdf')
        plt.close()

    print(f"  Saved: {out_png.name} ({dpi_raster} dpi), {out_pdf.name} (vector)")




def run_candidate_region_finemap(merged_df, gwas_results):
    """
    Candidate region fine-mapping:
    Extract SNPs within +-50kb of each optional candidate gene,
    then run association tests restricted to these regions.

    This dramatically reduces multiple testing burden (4.6M -> ~5k-20k SNPs)
    and increases power for detecting disease-relevant loci.
    Significance threshold: Bonferroni 0.05 / n_candidate_snps.
    """
    print("\n" + "=" * 70)
    print("STEP 12: Candidate Region Fine-Mapping")
    print("=" * 70)

    GCTA = _resolve_tool("GCTA_BIN", "gcta64", "gcta")
    GEMMA = _resolve_tool("GEMMA_BIN", "gemma")
    bfile = str(BASE_DIR / "data" / "monkey_snp_sex_qc")

    finemap_dir = OUTPUT_DIR / "finemap"
    finemap_dir.mkdir(exist_ok=True)

    # Optional candidate-gene table; skip if absent
    cand_file = BASE_DIR / "output" / "phase7_prioritization" / "top_100_candidates.csv"
    if not cand_file.exists():
        print("  Optional candidate-gene table not found, skipping")
        return {}

    cands = pd.read_csv(cand_file)
    print(f"  Loaded {len(cands)} candidate variants")

    # Build chrom map (accession -> chromosome number)
    acc_to_chr = build_chrom_map()

    # Parse GFF to get gene coordinates
    gff_file = BASE_DIR / "data" / "reference" / "Rhinopithecus_roxellana.gff"
    gene_coords = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.split('\t')
            if len(p) < 9 or p[2] != 'gene':
                continue
            c = acc_to_chr.get(p[0])
            if not c:
                continue
            a = p[8]
            nm = ''
            for t in ['Name=', 'gene=']:
                if t in a:
                    nm = a.split(t)[1].split(';')[0]
                    break
            if nm:
                gene_coords[nm] = (c, int(p[3]), int(p[4]))

    # Candidate gene names from the optional table
    cand_genes = set()
    if 'gene_name' in cands.columns:
        cand_genes = set(cands['gene_name'].dropna().unique())
    print(f"  Unique candidate genes: {len(cand_genes)}")

    # Build SNP extraction list: all SNPs within gene +- 50kb
    FLANK = 50000
    bim = pd.read_csv(
        BASE_DIR / "data" / "monkey_snp_sex_qc.bim", sep='\t',
        header=None, names=['CHR', 'SNP', 'CM', 'POS', 'A1', 'A2'])
    bim['CHR'] = bim['CHR'].astype(str)

    snp_list = []
    gene_windows = []
    for gene in cand_genes:
        if gene not in gene_coords:
            continue
        chrom, start, end = gene_coords[gene]
        mask = ((bim['CHR'] == chrom) &
                (bim['POS'] >= start - FLANK) &
                (bim['POS'] <= end + FLANK))
        snps = bim.loc[mask, 'SNP'].tolist()
        snp_list.extend(snps)
        gene_windows.append({'gene': gene, 'chrom': chrom,
                              'start': start - FLANK, 'end': end + FLANK,
                              'n_snps': len(snps)})

    snp_list = list(set(snp_list))
    print(f"  SNPs in candidate regions (+-50kb): {len(snp_list)}")

    if len(snp_list) == 0:
        print("  No SNPs found in candidate regions")
        return {}

    bonf_threshold = 0.05 / len(snp_list)
    print(f"  Bonferroni threshold: {bonf_threshold:.2e}")

    snp_file = str(finemap_dir / "candidate_snps.txt")
    with open(snp_file, 'w') as f:
        f.write('\n'.join(snp_list))

    # Check GRM exists (reuse from GWAS step)
    grm_prefix = str(OUTPUT_DIR / "gwas" / "grm")
    if not Path(grm_prefix + '.grm.bin').exists():
        print("  GRM not found, computing...")
        grm_prefix = _make_grm(OUTPUT_DIR / "gwas")
        if grm_prefix is None:
            return {}

    sex_covar = _make_sex_covar(merged_df, OUTPUT_DIR / "gwas")
    finemap_results = {}

    # CHS: GCTA-MLMA on candidate SNPs
    print("\n  Fine-mapping CHS with GCTA-MLMA...")
    pheno_chs = str(OUTPUT_DIR / "gwas" / "pheno_chs.txt")
    if not Path(pheno_chs).exists():
        df_tmp = merged_df[['IID', 'CHS']].copy()
        df_tmp.insert(0, 'FID', 0)
        df_tmp.to_csv(pheno_chs, sep='\t', index=False, header=False)

    mlma_prefix = str(finemap_dir / "finemap_chs")
    r = subprocess.run(
        [GCTA, '--mlma', '--bfile', bfile,
         '--grm', grm_prefix,
         '--extract', snp_file,
         '--pheno', pheno_chs,
         '--covar', sex_covar,
         '--out', mlma_prefix,
         '--thread-num', '16'],
        capture_output=True, text=True)
    mlma_file = Path(mlma_prefix + '.mlma')
    if mlma_file.exists():
        df = pd.read_csv(mlma_file, sep='\t')
        df = df.dropna(subset=['p']).copy()
        df = df[df['p'] > 0]
        df = df.rename(columns={'Chr': '#CHROM', 'bp': 'POS', 'SNP': 'ID',
                                 'b': 'BETA', 'se': 'SE', 'p': 'P',
                                 'Freq': 'A1_FREQ', 'A1': 'A1'})
        df['#CHROM'] = df['#CHROM'].astype(str)
        df['bonf_sig'] = df['P'] < bonf_threshold
        finemap_results['CHS'] = df
        n_sig = df['bonf_sig'].sum()
        print(f"  CHS: {len(df)} SNPs tested, {n_sig} Bonferroni-significant")
        top = df.nsmallest(20, 'P')
        print(top[['#CHROM', 'POS', 'ID', 'A1_FREQ', 'BETA', 'SE', 'P']].to_string(index=False))
        df.to_csv(finemap_dir / "finemap_CHS.csv", index=False)
    else:
        print(f"  GCTA failed: {r.stderr[-200:]}")

    # Eye / Finger: GEMMA LMM on candidate SNPs
    gemma_k = OUTPUT_DIR / "gwas" / "output" / "gemma_grm.cXX.txt"
    if not gemma_k.exists():
        print("  GEMMA GRM not found, skipping binary fine-mapping")
        gemma_k = None

    fam_orig = pd.read_csv(
        BASE_DIR / "data" / "monkey_snp_sex_qc.fam",
        sep=r'\s+', header=None,
        names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO'])

    for col, nc in [('Has_Eye_Disease', int(merged_df['Has_Eye_Disease'].sum())),
                    ('Has_Finger_Joint_Abnormality',
                     int(merged_df['Has_Finger_Joint_Abnormality'].sum()))]:
        if col not in merged_df.columns or nc < 5 or gemma_k is None:
            continue
        print(f"\n  Fine-mapping {col} with GEMMA LMM (cases={nc})...")
        import shutil
        tmp_prefix = str(finemap_dir / f"tmp_fm_{col}")
        fam_tmp = tmp_prefix + '.fam'
        fam = fam_orig.copy()
        pheno_map = merged_df.set_index('IID')[col]
        fam['PHENO'] = fam['IID'].map(pheno_map).apply(
            lambda x: 2 if x == 1 else (1 if x == 0 else -9))
        fam.to_csv(fam_tmp, sep=' ', index=False, header=False)
        for ext, src in [('.bim', BASE_DIR / "data" / "monkey_snp_sex_qc.bim"),
                         ('.bed', BASE_DIR / "data" / "monkey_snp_sex_qc.bed")]:
            dst = Path(tmp_prefix + ext)
            if not dst.exists():
                shutil.copy(str(src), str(dst))

        out_name = f"finemap_{col.lower()}"
        r = subprocess.run(
            [GEMMA, '-bfile', tmp_prefix,
             '-k', str(gemma_k),
             '-snps', snp_file,
             '-lmm', '4',
             '-o', out_name],
            capture_output=True, text=True, cwd=str(finemap_dir))

        assoc_file = finemap_dir / "output" / f"{out_name}.assoc.txt"
        if not assoc_file.exists():
            assoc_file = finemap_dir / f"{out_name}.assoc.txt"
        if assoc_file.exists():
            df = pd.read_csv(assoc_file, sep='\t')
            p_col = 'p_wald' if 'p_wald' in df.columns else 'p_lrt'
            df = df.dropna(subset=[p_col]).copy()
            df = df[df[p_col] > 0]
            df = df.rename(columns={'chr': '#CHROM', 'ps': 'POS', 'rs': 'ID',
                                     'af': 'A1_FREQ', 'beta': 'BETA',
                                     'se': 'SE', p_col: 'P'})
            df['#CHROM'] = df['#CHROM'].astype(str)
            df['bonf_sig'] = df['P'] < bonf_threshold
            finemap_results[col] = df
            n_sig = df['bonf_sig'].sum()
            print(f"  {col}: {len(df)} SNPs, {n_sig} Bonferroni-significant")
            top = df.nsmallest(20, 'P')[
                [c for c in ['#CHROM', 'POS', 'ID', 'A1_FREQ', 'BETA', 'SE', 'P']
                 if c in df.columns]]
            print(top.to_string(index=False))
            df.to_csv(finemap_dir / f"finemap_{col}.csv", index=False)
        else:
            print(f"  GEMMA failed: {r.stderr[-200:]}")

    # Save gene window summary
    pd.DataFrame(gene_windows).to_csv(finemap_dir / "candidate_gene_windows.csv", index=False)
    return finemap_results


def run_fisher_carrier_test(merged_df):
    """
    Fisher's exact test for candidate gene carrier analysis.

    For each candidate gene, define a 'carrier' as an individual
    with >=1 prioritized deleterious variant.
    Test: carrier status x disease status (2x2 Fisher exact test).

    This is assumption-free and well-suited for n_cases as small as 6-14.
    Run separately for each disease type (Eye, Finger) and overall Has_Disease.
    """
    print("\n" + "=" * 70)
    print("STEP 13: Fisher Exact Test - Carrier Analysis")
    print("=" * 70)

    fisher_dir = OUTPUT_DIR / "candidate_genes"
    fisher_dir.mkdir(exist_ok=True)

    # Optional table; skip if absent
    prio_file = BASE_DIR / "output" / "phase7_prioritization" / "prioritized_mutations.csv"
    if not prio_file.exists():
        print("  Optional prioritized-mutations table not found, skipping")
        return None

    prio = pd.read_csv(prio_file)
    if 'gene_name' not in prio.columns or 'chrom' not in prio.columns:
        print("  Missing gene_name or chrom column, skipping")
        return None

    # Build carrier matrix: individual x gene (1 if carries >=1 deleterious variant)
    acc_to_chr = build_chrom_map()
    prio['chr_num'] = prio['chrom'].map(acc_to_chr)
    prio = prio.dropna(subset=['chr_num', 'pos'])
    prio['pos'] = prio['pos'].astype(int)

    bim = pd.read_csv(
        BASE_DIR / "data" / "monkey_snp_sex_qc.bim", sep='\t',
        header=None, names=['CHR', 'SNP', 'CM', 'POS', 'A1', 'A2'])
    bim['CHR'] = bim['CHR'].astype(str)

    matched = prio.merge(bim[['CHR', 'SNP', 'POS']],
                         left_on=['chr_num', 'pos'],
                         right_on=['CHR', 'POS'], how='inner')
    if matched.empty:
        print("  No variants matched BIM")
        return None
    print(f"  Matched {len(matched)} prioritized variants to BIM")

    # Extract genotypes via plink2
    snp_file = str(fisher_dir / "fisher_snps.txt")
    matched['SNP'].drop_duplicates().to_csv(snp_file, index=False, header=False)

    subprocess.run(
        ['plink2', '--bfile', str(BASE_DIR / "data" / "monkey_snp_sex_qc"),
         '--extract', snp_file,
         '--export', 'A', '--out', str(fisher_dir / "fisher_geno")],
        capture_output=True, text=True)

    raw = fisher_dir / "fisher_geno.raw"
    if not raw.exists():
        print("  Genotype extraction failed")
        return None

    geno = pd.read_csv(raw, sep=r'\s+', engine='python')
    meta_cols = {'FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'}
    snp_cols = [c for c in geno.columns if c not in meta_cols]

    # Map SNP columns back to genes
    snp_to_gene = {}
    for _, row in matched.iterrows():
        snp_to_gene[row['SNP']] = row['gene_name']

    col_to_gene = {}
    for c in snp_cols:
        snp_id = c.rsplit('_', 1)[0]
        if snp_id in snp_to_gene:
            col_to_gene[c] = snp_to_gene[snp_id]

    # Build carrier matrix (1 if any alt allele dosage > 0)
    carrier = pd.DataFrame({'IID': geno['IID']})
    for gene in set(col_to_gene.values()):
        cs = [c for c, g in col_to_gene.items() if g == gene]
        if cs:
            carrier[gene] = (geno[cs].fillna(0).sum(axis=1) > 0).astype(int)

    carrier = carrier.merge(
        merged_df[['IID', 'Has_Disease', 'Has_Eye_Disease',
                   'Has_Finger_Joint_Abnormality']],
        on='IID', how='left')

    gene_cols = [c for c in carrier.columns
                 if c not in {'IID', 'Has_Disease', 'Has_Eye_Disease',
                               'Has_Finger_Joint_Abnormality'}]
    print(f"  Testing {len(gene_cols)} genes")

    from scipy.stats import fisher_exact

    results = []
    for phenotype in ['Has_Disease', 'Has_Eye_Disease', 'Has_Finger_Joint_Abnormality']:
        if phenotype not in carrier.columns:
            continue
        n_cases = int(carrier[phenotype].sum())
        if n_cases < 2:
            continue
        for gene in gene_cols:
            sub = carrier[[gene, phenotype]].dropna()
            # 2x2: rows=carrier/non-carrier, cols=case/control
            a = int(((sub[gene] == 1) & (sub[phenotype] == 1)).sum())  # carrier+case
            b = int(((sub[gene] == 1) & (sub[phenotype] == 0)).sum())  # carrier+ctrl
            c = int(((sub[gene] == 0) & (sub[phenotype] == 1)).sum())  # non-carrier+case
            d = int(((sub[gene] == 0) & (sub[phenotype] == 0)).sum())  # non-carrier+ctrl
            if a + b == 0:
                continue
            oddsratio, pval = fisher_exact([[a, b], [c, d]], alternative='two-sided')
            results.append({
                'phenotype': phenotype,
                'gene': gene,
                'n_carrier_case': a,
                'n_carrier_ctrl': b,
                'n_noncarrier_case': c,
                'n_noncarrier_ctrl': d,
                'odds_ratio': oddsratio,
                'p_fisher': pval,
            })

    if not results:
        print("  No results")
        return None

    rdf = pd.DataFrame(results)
    # FDR correction within each phenotype
    from statsmodels.stats.multitest import multipletests
    rdf_list = []
    for pheno, sub in rdf.groupby('phenotype'):
        sub = sub.copy()
        if len(sub) > 1:
            sub['p_adj'] = multipletests(sub['p_fisher'], method='fdr_bh')[1]
        else:
            sub['p_adj'] = sub['p_fisher']
        rdf_list.append(sub)
    rdf = pd.concat(rdf_list).sort_values('p_fisher')

    print("\n  Top Fisher exact test results:")
    print(rdf.head(20)[['phenotype', 'gene', 'n_carrier_case', 'n_carrier_ctrl',
                          'odds_ratio', 'p_fisher', 'p_adj']].to_string(index=False))
    rdf.to_csv(fisher_dir / "fisher_carrier_test.csv", index=False)
    return rdf


# =============================================================================
# 8. CANDIDATE GENE ASSOCIATION ANALYSIS
# =============================================================================

def run_candidate_gene_analysis(merged_df, gwas_results):
    """
    Gene-level association analysis:
    1. Aggregate GWAS p-values per gene (Simes method)
    2. Cross-reference with optional disease-associated candidate genes
    3. Gene burden test for prioritized deleterious variants (plink2)
    """
    print("\n" + "=" * 70)
    print("STEP 11: Candidate Gene Association Analysis")
    print("=" * 70)

    gene_dir = OUTPUT_DIR / "candidate_genes"
    gene_dir.mkdir(exist_ok=True)

    # Optional table; skip if absent
    da_file = BASE_DIR / "output" / "phase7_prioritization" / "disease_associations_summary.csv"
    if not da_file.exists():
        print("  Optional disease-association table not found, skipping")
        return {}

    disease_assoc = pd.read_csv(da_file)
    cand_set = set(disease_assoc['monkey_gene'].dropna().str.upper())
    print(f"Disease-associated candidate genes: {len(cand_set)}")

    acc_to_chr = build_chrom_map()

    gff_file = BASE_DIR / "data" / "reference" / "Rhinopithecus_roxellana.gff"
    genes = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.split('\t')
            if len(p) < 9 or p[2] != 'gene':
                continue
            c = acc_to_chr.get(p[0])
            if not c:
                continue
            a = p[8]
            nm = ''
            for t in ['Name=', 'gene=']:
                if t in a:
                    nm = a.split(t)[1].split(';')[0]
                    break
            genes.append((nm, c, int(p[3]), int(p[4])))
    gc = pd.DataFrame(genes, columns=['gene', 'chrom', 'start', 'end'])
    print(f"Parsed {len(gc)} genes from GFF")

    all_results = {}

    for trait, gdf in gwas_results.items():
        print(f"\nGene-level aggregation for {trait}...")
        rows = []

        for cv in gdf['#CHROM'].unique():
            cg = gdf[gdf['#CHROM'] == cv]
            gg = gc[gc['chrom'] == str(cv)]
            if gg.empty or cg.empty:
                continue

            pos_arr = cg['POS'].values
            pv_arr = cg['P'].values

            for _, g in gg.iterrows():
                mask = (pos_arr >= g['start']) & (pos_arr <= g['end'])
                if not mask.any():
                    continue
                gp = pv_arr[mask]
                n = len(gp)
                sp = np.sort(gp)
                simes = min(1.0, np.min(n * sp / np.arange(1, n + 1)))
                rows.append({
                    'gene': g['gene'], 'chrom': cv,
                    'start': g['start'], 'end': g['end'],
                    'n_snps': n, 'min_p': gp.min(), 'simes_p': simes,
                    'n_nominal': int((gp < 0.05).sum()),
                })

        if not rows:
            continue

        rdf = pd.DataFrame(rows)
        rdf['simes_p_adj'] = multipletests(rdf['simes_p'], method='fdr_bh')[1]
        rdf['is_candidate'] = rdf['gene'].str.upper().isin(cand_set)
        rdf = rdf.sort_values('simes_p')
        all_results[trait] = rdf

        print(f"\n  Top 20 genes:")
        print(rdf.head(20)[['gene', 'chrom', 'n_snps', 'min_p',
                             'simes_p', 'simes_p_adj', 'is_candidate']].to_string(index=False))

        cands = rdf[rdf['is_candidate']].sort_values('simes_p')
        if len(cands) > 0:
            print(f"\n  Disease candidate genes:")
            print(cands.head(30)[['gene', 'chrom', 'n_snps', 'min_p',
                                   'simes_p', 'simes_p_adj']].to_string(index=False))

        rdf.to_csv(gene_dir / f"gene_assoc_{trait}.csv", index=False)

    # Optional table; skip if absent
    prio_file = BASE_DIR / "output" / "phase7_prioritization" / "prioritized_mutations.csv"
    if prio_file.exists():
        prio = pd.read_csv(prio_file)
        bdf = _gene_burden_test(merged_df, prio, acc_to_chr, gene_dir)
        if bdf is not None:
            all_results['burden'] = bdf

    return all_results


def _gene_burden_test(merged_df, prio_df, acc_to_chr, gene_dir):
    """
    Per-individual gene burden test:
    extract genotypes of prioritized deleterious variants via plink2,
    aggregate per gene, test burden ~ disease association.
    """
    if 'gene_name' not in prio_df.columns or 'chrom' not in prio_df.columns:
        return None

    print("\n--- Gene Burden Test (plink2 genotype extraction) ---")
    bfile = str(BASE_DIR / "data" / "monkey_snp_sex_qc")

    prio = prio_df.copy()
    prio['chr_num'] = prio['chrom'].map(acc_to_chr)
    prio = prio.dropna(subset=['chr_num', 'pos'])
    prio['pos'] = prio['pos'].astype(int)
    prio['chr_num'] = prio['chr_num'].astype(str)

    bim = pd.read_csv(
        BASE_DIR / "data" / "monkey_snp_sex_qc.bim", sep='\t',
        header=None, names=['CHR', 'SNP', 'CM', 'POS', 'A1', 'A2'])
    bim['CHR'] = bim['CHR'].astype(str)

    matched = prio.merge(bim[['CHR', 'SNP', 'POS']],
                         left_on=['chr_num', 'pos'],
                         right_on=['CHR', 'POS'], how='inner')
    if matched.empty:
        print("  No variants matched BIM for burden test")
        return None
    print(f"  Matched {len(matched)} prioritized variants to BIM")

    snp_file = gene_dir / "burden_snps.txt"
    matched['SNP'].to_csv(snp_file, index=False, header=False)

    subprocess.run(
        ['plink2', '--bfile', bfile, '--extract', str(snp_file),
         '--export', 'A', '--out', str(gene_dir / "burden_geno")],
        capture_output=True, text=True)

    raw = gene_dir / "burden_geno.raw"
    if not raw.exists():
        print("  Genotype extraction failed")
        return None

    geno = pd.read_csv(raw, sep=r'\s+', engine='python')
    meta_cols = {'FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'}
    snp_cols = [c for c in geno.columns if c not in meta_cols]
    print(f"  Extracted {len(snp_cols)} variant dosages for {len(geno)} individuals")

    s2g = dict(zip(matched['SNP'], matched['gene_name']))
    col_gene = {}
    for c in snp_cols:
        snp_id = c.rsplit('_', 1)[0]
        if snp_id in s2g:
            col_gene[c] = s2g[snp_id]

    gene_names = sorted(set(col_gene.values()))
    burden = pd.DataFrame({'IID': geno['IID']})
    for g in gene_names:
        cs = [c for c, gn in col_gene.items() if gn == g]
        if cs:
            burden[g] = geno[cs].fillna(0).sum(axis=1)

    burden = burden.merge(merged_df[['IID', 'Has_Disease', 'CHS']], on='IID', how='left')

    results = []
    for g in gene_names:
        if g not in burden.columns:
            continue
        dv = burden[burden['Has_Disease'] == 1][g].dropna()
        hv = burden[burden['Has_Disease'] == 0][g].dropna()
        if len(dv) < 2 or len(hv) < 2 or (dv.std() == 0 and hv.std() == 0):
            continue
        try:
            _, p = mannwhitneyu(dv, hv, alternative='two-sided')
        except Exception:
            continue
        ps = np.sqrt((dv.std()**2 + hv.std()**2) / 2)
        d = (dv.mean() - hv.mean()) / ps if ps > 0 else 0
        results.append({
            'gene': g,
            'n_variants': sum(1 for c, gn in col_gene.items() if gn == g),
            'burden_diseased': dv.mean(),
            'burden_healthy': hv.mean(),
            'cohens_d': d,
            'p_value': p
        })

    if not results:
        return None

    bdf = pd.DataFrame(results).sort_values('p_value')
    if len(bdf) > 1:
        bdf['p_adj'] = multipletests(bdf['p_value'], method='fdr_bh')[1]
    else:
        bdf['p_adj'] = bdf['p_value']

    print(f"\n  Gene burden test ({len(bdf)} genes tested):")
    print(bdf.head(20).to_string(index=False))
    bdf.to_csv(gene_dir / "gene_burden_test.csv", index=False)
    return bdf


# =============================================================================
# 9. SUMMARY REPORT
# =============================================================================

def generate_summary_report(merged_df, corr_df, cc_df, hypothesis_results,
                            gwas_results=None, gene_results=None,
                            finemap_results=None, fisher_results=None):
    """
    Generate a comprehensive summary report.
    """
    if finemap_results is None:
        finemap_results = {}
    if fisher_results is None:
        fisher_results = pd.DataFrame()
    print("\n" + "=" * 70)
    print("STEP 9: Generating Summary Report")
    print("=" * 70)
    
    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append("GENOTYPE-PHENOTYPE ASSOCIATION ANALYSIS REPORT")
    report_lines.append("Golden Snub-Nosed Monkey (Rhinopithecus roxellana) - Shennongjia Population")
    report_lines.append("=" * 80)
    report_lines.append(f"\nAnalysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Sample overview
    report_lines.append("\n" + "-" * 40)
    report_lines.append("SAMPLE OVERVIEW")
    report_lines.append("-" * 40)
    report_lines.append(f"Total individuals: {len(merged_df)}")
    report_lines.append(f"Individuals with disease: {merged_df['Has_Disease'].sum()} ({merged_df['Has_Disease'].mean()*100:.1f}%)")
    report_lines.append(f"Healthy individuals: {(merged_df['Has_Disease'] == 0).sum()}")
    
    # Disease breakdown
    report_lines.append("\nDisease Types:")
    for col in merged_df.columns:
        if col.startswith('Has_') and col != 'Has_Disease':
            count = merged_df[col].sum()
            if count > 0:
                report_lines.append(f"  - {col.replace('Has_', '').replace('_', ' ')}: {int(count)}")
    
    # Statistical significance overview
    report_lines.append("\n" + "-" * 40)
    report_lines.append("STATISTICAL SIGNIFICANCE OVERVIEW")
    report_lines.append("-" * 40)
    report_lines.append("Significance levels: *** p<0.001, ** p<0.01, * p<0.05")
    report_lines.append("Multiple testing correction: Benjamini-Hochberg FDR")
    
    # Count significant results
    n_corr_uncorrected = len(corr_df[corr_df['Spearman_p'] < 0.05]) if len(corr_df) > 0 else 0
    n_corr_fdr = len(corr_df[corr_df['Spearman_p_adj'] < 0.05]) if len(corr_df) > 0 else 0
    n_cc_uncorrected = len(cc_df[cc_df['P_value'] < 0.05]) if len(cc_df) > 0 else 0
    n_cc_fdr = len(cc_df[cc_df['P_adj'] < 0.05]) if len(cc_df) > 0 else 0
    
    report_lines.append(f"\nCorrelation Analysis:")
    report_lines.append(f"  - Significant (uncorrected p<0.05): {n_corr_uncorrected}")
    report_lines.append(f"  - Significant (FDR-corrected p<0.05): {n_corr_fdr}")
    
    report_lines.append(f"\nCase-Control Analysis:")
    report_lines.append(f"  - Significant (uncorrected p<0.05): {n_cc_uncorrected}")
    report_lines.append(f"  - Significant (FDR-corrected p<0.05): {n_cc_fdr}")
    
    # Genomic summary
    report_lines.append("\n" + "-" * 40)
    report_lines.append("GENOMIC SUMMARY")
    report_lines.append("-" * 40)
    report_lines.append(f"Mean F_ROH: {merged_df['F_ROH'].mean():.4f} ± {merged_df['F_ROH'].std():.4f}")
    if 'Total_Genetic_Load' in merged_df.columns:
        report_lines.append(f"Mean Total Genetic Load: {merged_df['Total_Genetic_Load'].mean():.4f} ± {merged_df['Total_Genetic_Load'].std():.4f}")
    report_lines.append(f"Mean Realized Load: {merged_df['Realized_Load'].mean():.4f} ± {merged_df['Realized_Load'].std():.4f}")
    if 'Hom_Realized_Load' in merged_df.columns:
        report_lines.append(f"  - Homozygous Component: {merged_df['Hom_Realized_Load'].mean():.4f} ± {merged_df['Hom_Realized_Load'].std():.4f}")
    if 'Het_Realized_Load' in merged_df.columns:
        report_lines.append(f"  - Heterozygous Component: {merged_df['Het_Realized_Load'].mean():.4f} ± {merged_df['Het_Realized_Load'].std():.4f}")
    if 'Potential_Load' in merged_df.columns:
        report_lines.append(f"Mean Potential Load: {merged_df['Potential_Load'].mean():.4f} ± {merged_df['Potential_Load'].std():.4f}")
    report_lines.append(f"Mean ROH-Associated Load: {merged_df['Load_in_ROH'].mean():.4f} ± {merged_df['Load_in_ROH'].std():.4f}")
    
    # Key findings
    report_lines.append("\n" + "-" * 40)
    report_lines.append("KEY FINDINGS")
    report_lines.append("-" * 40)
    
    # Top correlations
    if len(corr_df) > 0:
        report_lines.append("\nTop Genotype-Phenotype Correlations (by p-value):")
        report_lines.append("(Significance: *** p<0.001, ** p<0.01, * p<0.05)")
        
        # Show all CHS correlations, sorted by p-value
        chs_corr = corr_df[corr_df['Phenotype_Variable'] == 'CHS'].head(10)
        for _, row in chs_corr.iterrows():
            # Determine significance level
            p_val = row['Spearman_p']
            p_adj = row['Spearman_p_adj']
            
            if p_val < 0.001:
                sig = '***'
            elif p_val < 0.01:
                sig = '**'
            elif p_val < 0.05:
                sig = '*'
            else:
                sig = ''
            
            # Add FDR indicator if still significant after correction
            fdr_sig = ' [FDR sig]' if p_adj < 0.05 else ''
            
            report_lines.append(f"  - {row['Genomic_Variable']} vs CHS: "
                              f"ρ={row['Spearman_rho']:.3f}, p={p_val:.4f}{sig}, "
                              f"p_adj={p_adj:.4f}{fdr_sig}")
        
        # Highlight FDR-significant correlations separately
        fdr_sig_corr = corr_df[corr_df['Spearman_p_adj'] < 0.05]
        if len(fdr_sig_corr) > 0:
            report_lines.append(f"\n*** {len(fdr_sig_corr)} correlations remain significant after FDR correction ***")
        else:
            report_lines.append("\nNote: No correlations remain significant after FDR correction (p_adj < 0.05)")
    
    # Case-control highlights
    if len(cc_df) > 0:
        report_lines.append("\nCase-Control Analysis (Diseased vs Healthy):")
        report_lines.append("(Significance: *** p<0.001, ** p<0.01, * p<0.05, † p<0.1)")
        
        # Show top results (p < 0.1)
        sig_results = cc_df[cc_df['P_value'] < 0.1].head(10)
        if len(sig_results) > 0:
            for _, row in sig_results.iterrows():
                direction = "higher" if row['Difference'] > 0 else "lower"
                p_val = row['P_value']
                p_adj = row['P_adj']
                
                if p_val < 0.001:
                    sig = '***'
                elif p_val < 0.01:
                    sig = '**'
                elif p_val < 0.05:
                    sig = '*'
                else:
                    sig = '†'
                
                fdr_sig = ' [FDR sig]' if p_adj < 0.05 else ''
                
                report_lines.append(f"  - {row['Variable']}: {direction} in diseased "
                                  f"(d={row['Cohens_d']:.3f}, p={p_val:.4f}{sig}, "
                                  f"p_adj={p_adj:.4f}{fdr_sig})")
        else:
            report_lines.append("  No variables show significant differences (p < 0.1)")
        
        # Count FDR-significant results
        fdr_sig_cc = cc_df[cc_df['P_adj'] < 0.05]
        if len(fdr_sig_cc) > 0:
            report_lines.append(f"\n*** {len(fdr_sig_cc)} variables remain significant after FDR correction ***")
        else:
            report_lines.append("\nNote: No variables remain significant after FDR correction (p_adj < 0.05)")
    
    # Hypothesis conclusions
    report_lines.append("\n" + "-" * 40)
    report_lines.append("HYPOTHESIS TESTING CONCLUSIONS")
    report_lines.append("-" * 40)
    
    for hyp, results in hypothesis_results.items():
        report_lines.append(f"\n{hyp}:")
        if results:
            for test, outcome in results.items():
                report_lines.append(f"  - {test}: {outcome}")
        else:
            report_lines.append("  - No strong support from current tests")
    
    # Interpretation
    report_lines.append("\n" + "-" * 40)
    report_lines.append("INTERPRETATION")
    report_lines.append("-" * 40)
    
    # Check for overall pattern
    froh_cv = merged_df['F_ROH'].std() / merged_df['F_ROH'].mean()
    short_roh_prop = merged_df['F_ROH_SHORT'].mean() / merged_df['F_ROH'].mean()
    
    if short_roh_prop > 0.8:
        report_lines.append("• ROH pattern dominated by short segments (>80%), indicating ancient")
        report_lines.append("  inbreeding events rather than recent consanguinity.")
    
    if froh_cv > 0.2 and froh_cv < 0.3:
        report_lines.append("• Moderate F_ROH variance suggests mixed contributions from both")
        report_lines.append("  historical bottleneck and ongoing restricted gene flow.")
    
    # Correlation interpretation
    if len(corr_df) > 0:
        chs_corr = corr_df[corr_df['Phenotype_Variable'] == 'CHS']
        sig_corrs = chs_corr[chs_corr['Spearman_p'] < 0.05]
        if len(sig_corrs) > 0:
            report_lines.append("• Significant correlations found between genomic metrics and health:")
            for _, row in sig_corrs.iterrows():
                if row['Spearman_rho'] > 0:
                    report_lines.append(f"  - Higher {row['Genomic_Variable']} associated with worse health")
                else:
                    report_lines.append(f"  - Higher {row['Genomic_Variable']} associated with better health")
    
    # GWAS results (LMM-based)
    if gwas_results:
        report_lines.append("\n" + "-" * 40)
        report_lines.append("GWAS RESULTS (LMM-based)")
        report_lines.append("-" * 40)
        report_lines.append("CHS: GCTA-MLMA (linear mixed model, GRM random effect, sex covariate).")
        report_lines.append("Eye/Finger: GEMMA LMM (centered relatedness matrix, lmm-4 Wald test).")
        report_lines.append("Resp/Skin: skipped (<5 cases). lambda reported per trait.")
        report_lines.append("Multiple testing: Bonferroni alpha=0.05 / n_SNPs (not 5e-8).")
        disease_binary = ['Has_Eye_Disease', 'Has_Respiratory_Disease',
                          'Has_Skin_Disease', 'Has_Finger_Joint_Abnormality']
        from scipy import stats as _stats_rep
        for trait, gdf in gwas_results.items():
            n_snp = len(gdf)
            p_b = 0.05 / n_snp if n_snp else np.nan
            n_bonf = int((gdf['P'] < p_b).sum()) if n_snp else 0
            n_sug = (gdf['P'] < 1e-5).sum()
            lam = (np.median(_stats_rep.chi2.isf(gdf['P'].values, 1))
                   / _stats_rep.chi2.ppf(0.5, 1))
            if trait in disease_binary and trait in merged_df.columns:
                nc = int(merged_df[trait].sum())
                report_lines.append(f"\n{trait} (cases={nc}, controls={68-nc}, lambda={lam:.3f}):")
            else:
                report_lines.append(f"\n{trait} (quantitative, n=68, lambda={lam:.3f}):")
            report_lines.append(f"  SNPs tested: {n_snp:,}; Bonferroni threshold p < {p_b:.2e}")
            report_lines.append(f"  Bonferroni-significant hits: {n_bonf}")
            report_lines.append(f"  Exploratory (p<1e-5): {n_sug}")
            top5 = gdf.nsmallest(5, 'P')
            for _, row in top5.iterrows():
                report_lines.append(f"  Chr{row['#CHROM']}:{int(row['POS'])} p={row['P']:.2e}")

    # Candidate gene results (Simes + burden)
    if gene_results:
        report_lines.append("\n" + "-" * 40)
        report_lines.append("CANDIDATE GENE ASSOCIATION")
        report_lines.append("-" * 40)
        for trait, gdf in gene_results.items():
            if trait == 'burden' and 'p_value' in gdf.columns:
                report_lines.append(f"\nGene Burden Test ({len(gdf)} genes):")
                for _, row in gdf.head(10).iterrows():
                    report_lines.append(
                        f"  {row['gene']}: d={row['cohens_d']:.3f}, p={row['p_value']:.4f}")
            elif 'simes_p' in gdf.columns:
                report_lines.append(f"\n{trait} - Top genes (Simes method):")
                for _, row in gdf.head(10).iterrows():
                    c = " [CANDIDATE]" if row.get('is_candidate') else ""
                    report_lines.append(f"  {row['gene']}: p={row['simes_p']:.2e}{c}")

    # Fine-mapping results
    if finemap_results:
        report_lines.append("\n" + "-" * 40)
        report_lines.append("CANDIDATE REGION FINE-MAPPING (+-50kb)")
        report_lines.append("-" * 40)
        for trait, fdf in finemap_results.items():
            n_sig = int(fdf['bonf_sig'].sum()) if 'bonf_sig' in fdf.columns else 0
            report_lines.append(f"\n{trait}: {len(fdf)} SNPs tested, "
                                 f"{n_sig} Bonferroni-significant")
            top5 = fdf.nsmallest(5, 'P')
            for _, row in top5.iterrows():
                sig_flag = " *" if row.get('bonf_sig') else ""
                report_lines.append(
                    f"  Chr{row['#CHROM']}:{int(row['POS'])} "
                    f"p={row['P']:.2e}{sig_flag}")

    # Fisher carrier test results
    if fisher_results is not None and len(fisher_results) > 0:
        report_lines.append("\n" + "-" * 40)
        report_lines.append("FISHER EXACT TEST - CARRIER ANALYSIS")
        report_lines.append("-" * 40)
        for pheno, sub in fisher_results.groupby('phenotype'):
            sig = sub[sub['p_adj'] < 0.05]
            nom = sub[sub['p_fisher'] < 0.05]
            report_lines.append(f"\n{pheno}: {len(sub)} genes tested, "
                                 f"{len(nom)} nominal (p<0.05), "
                                 f"{len(sig)} FDR-significant")
            for _, row in sub.nsmallest(5, 'p_fisher').iterrows():
                fdr = " [FDR]" if row['p_adj'] < 0.05 else ""
                report_lines.append(
                    f"  {row['gene']}: OR={row['odds_ratio']:.2f}, "
                    f"p={row['p_fisher']:.3e}{fdr}")

    report_lines.append("\n" + "=" * 80)
    report_lines.append("END OF REPORT")
    report_lines.append("=" * 80)
    
    # Save report
    report_text = "\n".join(report_lines)
    with open(OUTPUT_DIR / "analysis_summary_report.txt", 'w') as f:
        f.write(report_text)
    
    print(report_text)
    
    return report_text


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    """
    Main execution function.
    """
    print("\n" + "=" * 80)
    print("GENOTYPE-PHENOTYPE ASSOCIATION ANALYSIS")
    print("Golden Snub-Nosed Monkey (Rhinopithecus roxellana)")
    print("Shennongjia Population")
    print("=" * 80)
    
    # Step 1: Load and process disease data
    disease_df, individual_disease = load_and_process_disease_data()
    
    # Step 2: Create full phenotype dataset
    phenotype_df = create_full_phenotype_dataset(individual_disease)
    
    # Step 3: Merge with genomic data
    merged_df = merge_with_genomic_data(phenotype_df)
    
    # Step 4: Correlation analysis
    corr_df = correlation_analysis(merged_df)
    
    # Step 5: Case-control analysis
    cc_df = case_control_analysis(merged_df)
    
    # Step 6: Regression analysis
    model_comp = regression_analysis(merged_df)
    
    # Step 7: Hypothesis testing
    hypothesis_results = test_hypotheses(merged_df, corr_df, cc_df)
    
    # Step 8: Visualizations
    fig_dir = create_visualizations(merged_df, corr_df, cc_df)
    
    # Step 9: Disease-specific analysis
    disease_specific = create_disease_type_analysis(merged_df, disease_df)
    
    # Step 10: GWAS (GCTA-MLMA / GEMMA LMM)
    gwas_results = run_gwas_analysis(merged_df)

    # Step 11: Candidate gene association (Simes + burden)
    gene_results = run_candidate_gene_analysis(merged_df, gwas_results)

    # Step 12: Candidate region fine-mapping (+-50kb, LMM)
    finemap_results = run_candidate_region_finemap(merged_df, gwas_results)

    # Step 13: Fisher exact carrier test
    fisher_results = run_fisher_carrier_test(merged_df)

    # Step 14: Generate summary report
    report = generate_summary_report(merged_df, corr_df, cc_df, hypothesis_results,
                                     gwas_results, gene_results,
                                     finemap_results, fisher_results)
    
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print(f"Results saved to: {OUTPUT_DIR}")
    print("=" * 80)
    
    return merged_df, corr_df, cc_df, hypothesis_results


def redraw_gwas_figures_only():
    """
    Rebuild Manhattan + QQ figures from existing GWAS outputs (no re-running GCTA/GEMMA).
    Expects GCTA mlma_chs.mlma and GEMMA output/*.assoc.txt under OUTPUT_DIR / gwas.
    Manhattan significance line: Bonferroni p < 0.05 / n_SNPs (per trait, from rows plotted).
    """
    import sys

    gwas_dir = OUTPUT_DIR / "gwas"
    fig_dir = OUTPUT_DIR / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    jobs = []

    mlma = gwas_dir / "mlma_chs.mlma"
    if mlma.exists():
        jobs.append(("CHS", mlma, "gcta_mlma"))
    else:
        print(f"  Skip CHS: missing {mlma}", file=sys.stderr)

    gemma_specs = [
        ("Has_Eye_Disease", "gemma_has_eye_disease.assoc.txt"),
        ("Has_Finger_Joint_Abnormality", "gemma_has_finger_joint_abnormality.assoc.txt"),
    ]
    for trait, fname in gemma_specs:
        p = gwas_dir / "output" / fname
        if not p.exists():
            p = gwas_dir / fname
        if p.exists():
            jobs.append((trait, p, "gemma"))
        else:
            print(f"  Skip {trait}: missing {fname}", file=sys.stderr)

    if not jobs:
        print("No GWAS result files found; nothing to plot.", file=sys.stderr)
        return

    print("Redrawing GWAS figures from saved results...")
    print("  Significance line: Bonferroni alpha=0.05 / n_SNPs (not 5e-8).")
    for trait, path, kind in jobs:
        if kind == "gcta_mlma":
            df = pd.read_csv(path, sep='\t', low_memory=False)
            df = df.dropna(subset=['p']).copy()
            df = df[np.isfinite(df['p']) & (df['p'] > 0)]
            df = df.rename(columns={'Chr': '#CHROM', 'bp': 'POS', 'p': 'P'})
        else:
            df = pd.read_csv(path, sep='\t', low_memory=False)
            p_col = 'p_wald' if 'p_wald' in df.columns else 'p_lrt'
            df = df.dropna(subset=[p_col]).copy()
            df = df[np.isfinite(df[p_col]) & (df[p_col] > 0)]
            df = df.rename(columns={'chr': '#CHROM', 'ps': 'POS', p_col: 'P'})
        df['#CHROM'] = df['#CHROM'].astype(str)
        print(f"  {trait}: {len(df)} variants from {path.name}")
        _plot_manhattan_qq(df, trait, fig_dir)
    print(f"Done. Figures in {fig_dir}")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] in ("--redraw-gwas", "redraw-gwas"):
        redraw_gwas_figures_only()
    else:
        main()

