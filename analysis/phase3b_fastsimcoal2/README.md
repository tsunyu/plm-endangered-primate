# Phase 3b: Demographic Inference with fastsimcoal2

**Quick Links:**
- 🚀 **New user?** → [QUICKSTART.md](QUICKSTART.md)
- ⚠️ **Autosomes only?** → [AUTOSOMES_GUIDE.txt](AUTOSOMES_GUIDE.txt)

**⚠️ IMPORTANT**: Standard practice for demographic inference is to use **AUTOSOMES ONLY**. Filter your VCF first using `filter_autosomes_only.sh` if starting from an unfiltered callset. See [AUTOSOMES_GUIDE.txt](AUTOSOMES_GUIDE.txt).

**Reference Genome**: 21 autosomes = 2,948,446,826 bp (auto-detected)

---

## Overview

This directory contains scripts for inferring demographic history using **fastsimcoal2**, a coalescent-based tool that uses the Site Frequency Spectrum (SFS) to estimate population size changes over time.

### Key Features

1. **Explicit demographic models** - Test specific hypotheses (bottleneck, expansion, etc.)
2. **Model comparison** - AIC-based framework to identify best-supported scenario
3. **Unfolded SFS** - Uses ancestral state information for increased statistical power
4. **Study parameters** - Uses 10-year generation time and reports N parameters as diploid Ne
5. **Confidence intervals** - Robust uncertainty quantification via bootstrap

### Fixed conventions for these analyses

- Main VCF for SFS and monomorphic-site estimation:
  `$PLM_BASE_DIR/data/hardfilted.snp.pass.autosomes.vcf.gz`
- Generation time: 10 years per generation
- Mutation rate: 1.36e-8 per bp per generation
- Size parameters: reported as diploid effective population size (`Ne`)
- Best-run and model-comparison likelihood: `MaxEstLhood`

For a full rerun from raw inputs:

```bash
python3 phase3b_step0_prepare_sfs.py
python3 phase3b_step0b_estimate_monomorphic_sites.py
python3 phase3b_step1_create_model_templates.py
bash phase3b_step2_run_fastsimcoal2_parallel.sh
python3 phase3b_step3_model_comparison.py
bash phase3b_step4_bootstrap_ci_parallel.sh
python3 phase3b_step5_analyze_results.py
python3 phase3b_step6_visualize_demographic.py
```

---

## Directory Structure

```
analysis/phase3b_fastsimcoal2/
├── README.md                              # This file
├── phase3b_step0_prepare_sfs.py           # Generate SFS from VCF
├── phase3b_step1_create_model_templates.py # Create demographic model files
├── phase3b_step2_run_fastsimcoal2.sh      # Run parameter estimation
├── phase3b_step3_model_comparison.py       # AIC-based model selection
├── phase3b_step4_bootstrap_ci.sh           # Bootstrap confidence intervals
├── phase3b_step5_analyze_results.py        # Extract and interpret results
├── phase3b_step6_visualize_demographic.py  # Create publication plots
└── models/                                # Model template files (.tpl, .est)

output/phase3b_fastsimcoal2/
├── sfs/                                   # Site frequency spectrum
│   ├── SNJ_DAFpop0.obs                   # Unfolded SFS (fastsimcoal2 format)
│   ├── sfs_statistics.txt                # SFS summary
│   └── sfs_plot.png                      # SFS visualization
├── models/                                # Model fitting results
│   ├── constant_ne/
│   ├── single_bottleneck/
│   ├── two_consecutive_bottlenecks/
│   ├── bottleneck_continuous_decline/
│   ├── bottleneck_recent_contraction/
│   └── complex_multi_event/
├── model_comparison/                      # Model selection results
│   ├── model_comparison.csv
│   ├── model_comparison_table.txt
│   └── model_comparison_plots.png
├── bootstrap/                             # Bootstrap CI results
├── plots/                                 # Final visualizations
├── parameter_estimates.txt                # Final parameter summary
└── parameter_estimates.csv                # Parameter table
```

---

## Pipeline Steps

### Step 0: Prepare Site Frequency Spectrum

**Script**: `phase3b_step0_prepare_sfs.py`

**Purpose**: Generate unfolded (derived) allele frequency spectrum using ancestral state information.

**Input**:
- Main VCF: `$PLM_BASE_DIR/data/hardfilted.snp.pass.autosomes.vcf.gz` (autosomes only, chr 1-21)
- Ancestral VCF: `$PLM_BASE_DIR/output/phase2_annotation/ancestral_states/variants_with_ancestral.vcf.gz`

**Output**:
- `SNJ_DAFpop0.obs` - Unfolded SFS in fastsimcoal2 format
- `sfs_statistics.txt` - Detailed SFS statistics
- `sfs_plot.png` - SFS visualization

**Run**:
```bash
python3 phase3b_step0_prepare_sfs.py
```

**Key features**:
- Uses ancestral allele information to polarize variants
- Generates unfolded (not folded) SFS for maximum statistical power
- Reports polarization success rate
- Creates diagnostic plots

---

### Optional Step 0b: Estimate Monomorphic Sites (Recommended)

**Script**: `phase3b_step0b_estimate_monomorphic_sites.py`

**Purpose**: Estimate the number of monomorphic ancestral sites and update the SFS `d0_0` bin for more accurate likelihoods and Ne scaling in fastsimcoal2.

**Input**:
- Reference genome size (auto-detected from autosome count)
- Main VCF: `$PLM_BASE_DIR/data/hardfilted.snp.pass.autosomes.vcf.gz`

**Output**:
- `monomorphic_recommendation.txt` - Recommended monomorphic site count and methods summary
- Updated `SNJ_DAFpop0.obs` (if you apply the recommendation)

**Run**:
```bash
python3 phase3b_step0b_estimate_monomorphic_sites.py
```

Then optionally apply the recommended value to the SFS.

**Recommendation**:
- Exploratory runs can use Step 0 alone.
- For manuscript-quality inference, run Step 0b after Step 0 and update `SNJ_DAFpop0.obs`.

---

### Step 1: Create Demographic Model Templates

**Script**: `phase3b_step1_create_model_templates.py`

**Purpose**: Generate `.tpl` (template) and `.est` (estimation) files for different demographic scenarios.

**Output**: Model files in `models/` directory

**Models created**:

1. **constant_ne** (1 parameter)
   - Null model with constant population size
   - Baseline for model comparison

2. **single_bottleneck** (5 parameters)
   - Single severe bottleneck followed by recovery
   - Parameters: NCUR, NBOT, NANC, TRECOVERY, TBOT

3. **two_consecutive_bottlenecks** (9 parameters)
   - Two independent bottlenecks caused by repeated environmental deterioration
   - Parameters: NCUR, NBOT2, NINTER, NBOT1, NANC, TRECENT_RECOVERY, TRECENT_BOT, TOLD_RECOVERY, TOLD_BOT

4. **bottleneck_continuous_decline** (7 parameters)
   - Bottleneck followed by sustained decline without full recovery
   - Parameters: NCUR, NMID, NBOT, NANC, TRECENT, TBOT, TANC

5. **bottleneck_recent_contraction** (7 parameters)
   - Ancient bottleneck + short recovery + recent severe contraction (~1.3 ka)
   - Parameters: NCUR, NRECOVER, NBOT, NANC, TRECENT, TRECOVERY_OLD, TBOT_OLD

6. **complex_multi_event** (9 parameters)
   - Flexible model with 3+ demographic size changes
   - Parameters: NCUR, N1, N2, N3, NANC, T1, T2, T3, T4

**Run**:
```bash
python3 phase3b_step1_create_model_templates.py
```

---

### Step 2: Run fastsimcoal2 Parameter Estimation

**Script**: `phase3b_step2_run_fastsimcoal2.sh`

**Purpose**: Run fastsimcoal2 to estimate parameters for all demographic models.

**Configuration**:
- Number of runs per model: 50 (independent optimizations)
- Cores: 4
- Simulations per likelihood: 100,000
- Optimization cycles: 40

**Output**: Best-fit parameters and `MaxEstLhood` values for each model

**Run**:
```bash
bash phase3b_step2_run_fastsimcoal2.sh
```

**Duration**: ~2-6 hours (depending on model complexity and number of runs)

**Note**: Requires fastsimcoal2 (fsc27 or fsc26) installed and in PATH.
Best runs are selected by `MaxEstLhood`; `MaxObsLhood` is the likelihood of the
fixed observed SFS and should not be used to rank optimization runs.

---

### Step 3: Model Comparison

**Script**: `phase3b_step3_model_comparison.py`

**Purpose**: Compare models using Akaike Information Criterion (AIC).

**Analysis**:
- Calculate AIC for each model: `AIC = 2k - 2ln(L)`
  - Note: fastsimcoal2 uses a multinomial likelihood for a single-population 1D SFS, so AIC is appropriate. For multi-population joint SFS it reports a composite likelihood; standard AIC may not apply (use clAIC instead).
- Calculate ΔAIC and Akaike weights
- Identify best-supported model
- Assess strength of support

**Output**:
- `model_comparison.csv` - Comparison table
- `model_comparison_table.txt` - Formatted summary
- `model_comparison_plots.png` - Visualizations

**Run**:
```bash
python3 phase3b_step3_model_comparison.py
```

**Interpretation guidelines**:
- ΔAIC < 2: Substantial support
- ΔAIC 4-7: Considerably less support
- ΔAIC > 10: Essentially no support

---

### Step 4: Bootstrap Confidence Intervals

**Script**: `phase3b_step4_bootstrap_ci.sh`

**Purpose**: Generate confidence intervals using parametric bootstrap.

**Method**:
1. Simulate SFS from best-fit model
2. Re-estimate parameters for each simulated SFS
3. Calculate 95% CI from bootstrap distribution

**Configuration**:
- Bootstrap replicates: 100
- Cores: 4

**Output**: Bootstrap parameter distributions

**Run**:
```bash
bash phase3b_step4_bootstrap_ci.sh
```

**Duration**: ~2-4 hours

**Optional**: Can specify model manually:
```bash
bash phase3b_step4_bootstrap_ci.sh single_bottleneck
```

---

### Step 5: Analyze Results

**Script**: `phase3b_step5_analyze_results.py`

**Purpose**: Extract parameter estimates and calculate confidence intervals.

**Output**:
- `parameter_estimates.txt` - Comprehensive summary
- `parameter_estimates.csv` - Parameter table
- Bootstrap distribution plots

**Run**:
```bash
python3 phase3b_step5_analyze_results.py
```

**Key outputs**:
- Parameter estimates with 95% CI
- Time conversion (generations → years)
- Biological interpretation

---

### Step 6: Visualize Demographic History

**Script**: `phase3b_step6_visualize_demographic.py`

**Purpose**: Create publication-quality figures.

**Plots created**:
1. **demographic_history.png/pdf** - Ne over time
2. **sfs_observed.png** - Observed SFS

**Run**:
```bash
python3 phase3b_step6_visualize_demographic.py
```

---

## Complete Pipeline Execution

### Option 1: Automated Pipeline (Recommended)

The master script automatically uses parallel versions for 8-10x speedup:

```bash
# Run complete pipeline with parallel optimization (default)
bash run_complete_pipeline.sh

# Skip bootstrap (faster, no confidence intervals)
bash run_complete_pipeline.sh --skip-bootstrap

# Use serial version (slower, for systems with <16 cores)
bash run_complete_pipeline.sh --no-parallel

# Test mode (reduced runs for quick testing)
bash run_complete_pipeline.sh --test
```

**Performance**: With 32-core CPU, complete pipeline takes ~6-8 hours (vs ~50-60 hours serial).

### Option 2: Manual Step-by-Step Execution

For manual control over each step:

```bash
# Step 0: Prepare SFS
python3 phase3b_step0_prepare_sfs.py

# Step 0b: Estimate monomorphic sites (recommended for accurate likelihood / Ne scale)
python3 phase3b_step0b_estimate_monomorphic_sites.py

# Step 1: Create model templates
python3 phase3b_step1_create_model_templates.py

# Step 2: Run parameter estimation (longest step)
bash phase3b_step2_run_fastsimcoal2_parallel.sh      # Parallel (recommended)
# OR: bash phase3b_step2_run_fastsimcoal2.sh         # Serial (slower)

# Step 3: Model comparison
python3 phase3b_step3_model_comparison.py

# Step 4: Bootstrap confidence intervals (optional but recommended)
bash phase3b_step4_bootstrap_ci_parallel.sh          # Parallel (recommended)
# OR: bash phase3b_step4_bootstrap_ci.sh             # Serial (slower)

# Step 5: Analyze results
python3 phase3b_step5_analyze_results.py

# Step 6: Create visualizations
python3 phase3b_step6_visualize_demographic.py
```

**Note**: For systems with 16+ CPU cores, parallel versions provide 8-10x speedup. See `PARALLEL_OPTIMIZATION.md` for details.

**Total estimated time**: 4-10 hours (depending on bootstrap)

---

## Key Parameters

| Parameter | Value | Notes |
|-----------|-------|-------|
| Generation time | 10 years | Biologically appropriate for this species |
| Mutation rate | 1.36×10⁻⁸ per bp/gen | Published estimate used in this project |
| Sample size | 68 diploids (136 chromosomes) | From VCF |
| Ancestral states | Available | Enables unfolded SFS |

---

## Expected Results

### Parameter Estimates

For the **single_bottleneck** model (example interpretation):

- **NCUR**: Current Ne (~100-1,000)
- **NBOT**: Bottleneck Ne (~100-1,000)
- **NANC**: Ancestral Ne (~5,000-50,000)
- **TBOT**: Bottleneck onset, more ancient (~500-5,000 generations = 5,000-50,000 years ago)
- **TRECOVERY**: Recovery start, more recent (~50-500 generations = 500-5,000 years ago)
- Timing: TRECOVERY < TBOT (TRECOVERY is more recent; TBOT is older)

### Biological Interpretation

- **Bottleneck severity**: NANC/NBOT ratio
- **Recovery magnitude**: NCUR/NBOT ratio
- **Current status**: NCUR/NANC ratio

---

---

## Requirements

### Software

- **fastsimcoal2** (fsc27 or fsc26)
  - Download: http://cmpg.unibe.ch/software/fastsimcoal27/
  - Must be in PATH or specify full path in scripts

### Python packages

```bash
pip install numpy pandas matplotlib scipy
```

### Input data

- ✅ Main VCF: `data/monkey_snp_autosomes_only.vcf.gz` (21 autosomes)
- ✅ Ancestral VCF: `output/phase2_annotation/ancestral_states/variants_with_ancestral.vcf.gz`
- ✅ Reference genome: 2,948,446,826 bp (21 autosomes)

---

## Fixing Tools

Several tools are provided to verify fixes and improve results:

### 1. Validate Ancestral Allele Parsing (Required)
```bash
python3 validate_ancestral_parsing.py
```
Verifies the critical ancestral allele fix is working correctly. Run BEFORE the pipeline.

### 2. Update Parameter Bounds (Recommended)
```bash
python3 update_parameter_bounds.py --type expanded --apply
```
Expands parameter ranges in .est files to avoid hitting boundaries. Run BEFORE the pipeline.

### 3. Estimate Monomorphic Sites (Recommended for accurate likelihood / Ne scale)
```bash
python3 phase3b_step0b_estimate_monomorphic_sites.py
```
Estimates and adds monomorphic site counts for more accurate likelihoods and Ne scale. Run BEFORE the pipeline (highly recommended, especially for final / publication analysis).

### 4. Verify Ne Scaling (Required)
```bash
python3 verify_ne_scaling.py
```
Determines if NCUR represents Ne or 2*Ne by comparing with census size and other methods. Run AFTER the pipeline.  
The result is written to `NE_SCALING_NOTE.txt` (e.g. `SCALING=Ne` or `SCALING=2Ne`), and `phase3b_step5_analyze_results.py` automatically applies the corresponding scaling to N* parameters before reporting Ne.

### 5. Alternative Bootstrap (If Needed)
```bash
python3 phase3b_step4b_alternative_bootstrap.py
```
Provides approximate confidence intervals if standard bootstrap fails. Run AFTER the pipeline.

---

## Troubleshooting

### Issue: No SFS generated

**Solution**: Check that both VCF files exist and are readable. Verify ancestral VCF has REF = ancestral allele.

### Issue: fastsimcoal2 not found

**Solution**: Install fastsimcoal2 and ensure it's in PATH:
```bash
which fsc27  # Should show path
```

### Issue: Low polarization rate

**Solution**: Check ancestral VCF format. REF allele should be ancestral allele.

### Issue: Model estimation fails

**Solution**: 
- Reduce number of simulations (-n parameter)
- Check parameter bounds in .est files
- Verify SFS file format

### Issue: Bootstrap fails

**Solution**:
- Ensure best-fit .par file exists
- Reduce number of bootstrap replicates
- Check disk space

---

## Literature References

- **fastsimcoal2**: Excoffier et al. (2013) Robust Demographic Inference from Genomic and SNP Data. PLoS Genetics.
- **SFS-based inference**: Gutenkunst et al. (2009) Inferring the Joint Demographic History. PLoS Genetics.
- **Model selection**: Burnham & Anderson (2002) Model Selection and Multimodel Inference.

---

## Citation

If you use this pipeline, please cite:
- fastsimcoal2: Excoffier et al. (2013)
- This analysis pipeline: Protein language modelling reveals a latent drift load underlying health risk in an endangered primate

---

## Contact

For questions about this pipeline:
- Check logs in `output/phase3b_fastsimcoal2/`
- Review model descriptions in `models/MODEL_DESCRIPTIONS.txt`
- See parameter summaries in `parameter_estimates.txt`

---

## Notes

- **Unfolded SFS**: Using ancestral states provides ~30-50% more statistical power than folded SFS
- **Generation time**: 10 years is biologically accurate for golden snub-nosed monkeys
- **Model selection**: Best model is determined objectively by AIC, not a priori assumptions
- **Uncertainty**: Bootstrap provides realistic confidence intervals for all parameters

---

Last updated: 2026-01-26
