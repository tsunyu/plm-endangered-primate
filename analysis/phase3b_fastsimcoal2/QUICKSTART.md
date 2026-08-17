# fastsimcoal2 Demographic Inference - Quick Start Guide

Get from data to results in 3 steps.

## Prerequisites

1. **fastsimcoal2 installed**:
   ```bash
   # Check installation
   which fsc28  # or fsc27, fsc26
   
   # If not installed, download from:
   # http://cmpg.unibe.ch/software/fastsimcoal27/
   ```

2. **Python packages**:
   ```bash
   pip install numpy pandas matplotlib
   ```

3. **Input data available**:
   - Main VCF: `data/hardfilted.snp.pass.autosomes.vcf.gz` (autosomes only, chr 1-21) ✓
   - Ancestral VCF: `output/phase2_annotation/ancestral_states/variants_with_ancestral.vcf.gz` ✓

4. **⚠️ IMPORTANT - Filter to autosomes only** (first time only):
   ```bash
   # Standard practice for fastsimcoal2: use autosomes only (chr 1-21)
   bash filter_autosomes_only.sh          # ~5 min
   bash update_to_autosomes_only.sh       # instant
   python3 check_vcf_chromosomes.py       # verify
   ```
   See [AUTOSOMES_GUIDE.txt](AUTOSOMES_GUIDE.txt) for why this is essential.

---

## Available Fixing Tools

Before/after running the pipeline, these tools are available:

| Tool | When | Purpose |
|------|------|---------|
| `validate_ancestral_parsing.py` | **Before** | ✅ Verify critical fix (required) |
| `update_parameter_bounds.py` | Before | Expand parameter ranges (recommended) |
| `phase3b_step0b_estimate_monomorphic_sites.py` | Before | Add monomorphic counts (recommended for accurate likelihood / Ne scale) |
| `verify_ne_scaling.py` | **After** | ✅ Determine Ne interpretation (required) |
| `phase3b_step4b_alternative_bootstrap.py` | After | Alternative CIs (if needed) |

**Quick examples:**
```bash
# Validate fix (before running)
python3 validate_ancestral_parsing.py

# Update bounds (before running, recommended)
python3 update_parameter_bounds.py --type expanded --apply

# Verify Ne scaling (after running)
python3 verify_ne_scaling.py
```

---

## Quick Start (Automated - Recommended)

Use the master script for automatic parallel optimization:

```bash
cd analysis/phase3b_fastsimcoal2

# Complete pipeline with parallel optimization (6-8 hours on 32-core CPU)
bash run_complete_pipeline.sh

# Or skip bootstrap for faster results (4-5 hours)
bash run_complete_pipeline.sh --skip-bootstrap
```

## Quick Start (Manual - Step by Step)

For manual control over each step:

```bash
cd analysis/phase3b_fastsimcoal2

# 1. Generate SFS (~5 min)
python3 phase3b_step0_prepare_sfs.py

# 2. Create model templates (~1 min)
python3 phase3b_step1_create_model_templates.py

# 3. Run parameter estimation (~4-5 hours parallel, ~40 hours serial)
bash phase3b_step2_run_fastsimcoal2_parallel.sh      # Recommended
# OR: bash phase3b_step2_run_fastsimcoal2.sh         # Serial version

# 4. Compare models (~5 min)
python3 phase3b_step3_model_comparison.py

# 5. Analyze results (~5 min)
python3 phase3b_step5_analyze_results.py

# 6. Create plots (~5 min)
python3 phase3b_step6_visualize_demographic.py
```

**Total time**: ~3-7 hours (mostly step 3)

**Results**: 
- `output/phase3b_fastsimcoal2/parameter_estimates.txt`
- `output/phase3b_fastsimcoal2/plots/demographic_history.png`

---

## Full Analysis (With Bootstrap)

For complete analysis with confidence intervals:

```bash
cd analysis/phase3b_fastsimcoal2

# Steps 1-4 as above, then:

# 5. Bootstrap confidence intervals (~2-4 hours)
bash phase3b_step4_bootstrap_ci.sh

# 6. Analyze results with CI (~5 min)
python3 phase3b_step5_analyze_results.py

# 7. Create plots (~5 min)
python3 phase3b_step6_visualize_demographic.py
```

**Total time**: ~5-11 hours

---

## What to Expect

### Step 0: SFS Generation

**Output**:
```
Loading ancestral states... ✓
Loaded 1,234,567 ancestral states
Processing main VCF...
Polarized 987,654 variants
Polarization success rate: 98.5%
```

**Check**: `output/phase3b_fastsimcoal2/sfs/sfs_plot.png` should show typical SFS (high singletons)

---

### Step 2: Parameter Estimation

**Output**:
```
Running 50 independent optimizations...
Run 1/50: single_bottleneck
  MaxEstLhood: -123456.78
Run 2/50: single_bottleneck
  MaxEstLhood: -123455.12
...
Best run: 23
```

**Check**: Each model should have `best_run/` directory with results

---

### Step 3: Model Comparison

**Output**:
```
Best model: single_bottleneck
  AIC: 246913.56
  Weight: 0.8523
  Parameters: 5
```

**Check**: `output/phase3b_fastsimcoal2/model_comparison/model_comparison_table.txt`

---

### Step 5: Results

**Output**:
```
Parameter    Estimate    95% CI Lower    95% CI Upper
NCUR         2341        1850            3120
NBOT         412         245             628
TBOT         127         89              178
```

**Check**: `output/phase3b_fastsimcoal2/parameter_estimates.txt`

---

## Key Results to Look For

1. **Best-supported model**: Which demographic scenario fits best?
2. **Current Ne**: Effective population size today
3. **Bottleneck Ne**: How severe was the bottleneck?
4. **Bottleneck timing**: When did it occur (in years)?
5. **Comparison with literature**: Does it match ~1,300 years ago?

---

## Common Issues

### Issue: "fastsimcoal2 not found"

```bash
# Download from http://cmpg.unibe.ch/software/fastsimcoal27/
# Extract and add to PATH:
export PATH=/path/to/fsc27:$PATH
```

### Issue: "No ancestral states found"

Check that ancestral VCF exists:
```bash
ls -lh $PLM_BASE_DIR/output/phase2_annotation/ancestral_states/variants_with_ancestral.vcf.gz
```

### Issue: Step 2 taking too long

Reduce number of runs:
```bash
# Edit phase3b_step2_run_fastsimcoal2.sh
# Change: N_RUNS=50
# To:     N_RUNS=20
```

### Issue: Out of memory

Reduce simulations:
```bash
# In phase3b_step2_run_fastsimcoal2.sh
# Change: -n 100000
# To:     -n 50000
```

---

## Interpreting Results

### Generation Time Conversion

fastsimcoal2 outputs are in **generations**. Convert to **years**:

```
Years = Generations × 10
```

Example:
- TBOT = 130 generations → **1,300 years ago** ✓ (matches literature!)

### Ne/Nc Ratio

Effective size (Ne) is usually smaller than census size (Nc):

- Ne/Nc = 0.1-0.5 is typical for mammals
- If Ne/Nc > 1: May indicate underestimated census size

### Bottleneck Severity

```
Severity = NANC / NBOT
```

- Severity > 10: **Severe bottleneck**
- Severity 2-10: **Moderate bottleneck**
- Severity < 2: **Mild bottleneck**

---

## Next Steps After Analysis

1. **Verify Ne scaling** (REQUIRED):
   ```bash
   python3 verify_ne_scaling.py
   ```
   - Determines if NCUR is Ne or 2*Ne
   - Uses census size comparison
   - Suggests additional verification methods

2. **Validate with independent data**:
   - ROH analysis
   - Heterozygosity estimates
   - LD-based Ne estimates

3. **Reporting**:
   - Use `demographic_history.pdf` for publication figures
   - Report parameter estimates with 95% CI
   - Cite best-supported model with AIC weights

---

## File Locations

### Input
- VCF: `data/monkey_snp_autosomes_only.vcf.gz` (21 autosomes)
- Ancestral: `output/phase2_annotation/ancestral_states/variants_with_ancestral.vcf.gz`

### Output
- SFS: `output/phase3b_fastsimcoal2/sfs/SNJ_DAFpop0.obs`
- Results: `output/phase3b_fastsimcoal2/parameter_estimates.txt`
- Plots: `output/phase3b_fastsimcoal2/plots/`

### Logs
- SFS: `output/phase3b_fastsimcoal2/sfs/sfs_generation.log`
- Estimation: `output/phase3b_fastsimcoal2/models/fastsimcoal2_runs.log`
- Bootstrap: `output/phase3b_fastsimcoal2/bootstrap/bootstrap.log`

---

## Quick Diagnostics

Check if analysis completed successfully:

```bash
cd $PLM_BASE_DIR/output/phase3b_fastsimcoal2

# SFS generated?
ls -lh sfs/SNJ_DAFpop0.obs

# Models estimated?
ls -d models/*/best_run

# Best model identified?
cat model_comparison/model_comparison.csv | head -2

# Parameters extracted?
cat parameter_estimates.txt | head -30

# Plots created?
ls plots/*.png
```

---

## Time Estimates

| Step | Time | Can skip? |
|------|------|-----------|
| 0. SFS | 5 min | No |
| 1. Models | 1 min | No |
| 2. Estimation | 2-6 hours | No |
| 3. Comparison | 5 min | No |
| 4. Bootstrap | 2-4 hours | Yes* |
| 5. Analysis | 5 min | No |
| 6. Plots | 5 min | No |

*Bootstrap is recommended for publication but not required for initial exploration.

---

## Getting Help

1. **Check logs**: All scripts write detailed logs
2. **Read README.md**: Comprehensive documentation
3. **Verify inputs**: Make sure VCF files are correct
4. **Test with subset**: Reduce N_RUNS for testing

---

Last updated: 2026-01-26
