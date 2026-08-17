#!/usr/bin/env python3
"""
Create fastsimcoal2 Demographic Model Templates
================================================

This script generates .tpl (template) and .est (estimation) files for
different demographic scenarios to be tested with fastsimcoal2.

Models:
  1. constant_ne
  2. single_bottleneck
  3. two_consecutive_bottlenecks
  4. bottleneck_continuous_decline
  5. bottleneck_recent_contraction
  6. complex_multi_event

Output: Model files in analysis/phase3b_fastsimcoal2/models/

Author: Demographic Analysis Pipeline
Date: 2026-01-26
"""

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir
import logging

# Configuration
BASE_DIR = get_base_dir()
MODEL_DIR = Path(__file__).resolve().parent / "models"
MODEL_DIR.mkdir(parents=True, exist_ok=True)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Population parameters
N_SAMPLES = 68  # diploid individuals
N_CHROMOSOMES = 136
MUTATION_RATE = 1.36e-8


def write_model(model_name, tpl_content, est_content):
    """Write .tpl and .est files for a model."""
    tpl_path = MODEL_DIR / f"{model_name}.tpl"
    est_path = MODEL_DIR / f"{model_name}.est"

    with open(tpl_path, 'w') as f:
        f.write(tpl_content)

    with open(est_path, 'w') as f:
        f.write(est_content)

    logger.info(f"  Created: {tpl_path}")
    logger.info(f"  Created: {est_path}")


def create_constant_ne_model():
    """Model 1: Constant population size (null model)."""
    model_name = "constant_ne"
    logger.info(f"Creating {model_name} model")
    
    tpl_content = f"""//Number of population samples (demes)
1
//Population effective sizes (diploid Ne reported by this pipeline)
NCUR
//Sample sizes
{N_CHROMOSOMES}
//Growth rates: negative growth implies population expansion
0
//Number of migration matrices: 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix 
0  historical events
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
FREQ 1 0 {MUTATION_RATE} OUTEXP
"""
    
    # Estimation file (.est)
    est_content = """// Priors file for constant Ne model
// Issue: this model has 1 parameter to estimate
[PARAMETERS]
//#isInt? #name   #dist.#min  #max
1        NCUR    unif  500    100000   output

[COMPLEX PARAMETERS]
"""
    
    write_model(model_name, tpl_content, est_content)
    return model_name, 1  # 1 parameter


def create_single_bottleneck_model():
    """Model 2: single bottleneck with recovery."""
    model_name = "single_bottleneck"
    logger.info(f"Creating {model_name} model")

    tpl_content = f"""//Number of population samples (demes)
1
//Population effective sizes (diploid Ne reported by this pipeline)
NCUR
//Sample sizes
{N_CHROMOSOMES}
//Growth rates
0
//Number of migration matrices
0
//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix 
2 historical events
TRECOVERY 0 0 0 RESIZEBOT 0 0
TBOT 0 0 0 RESIZEANC 0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
FREQ 1 0 {MUTATION_RATE} OUTEXP
"""

    est_content = """// Priors file for single bottleneck model
[PARAMETERS]
//#isInt? #name       #dist.#min  #max
1        NCUR         unif  500    100000   output
1        NBOT         unif  50     50000    output
1        NANC         unif  500    200000   output
1        TRECOVERY    unif  10     3000     output
1        DTBOT        unif  1      17000    output

[COMPLEX PARAMETERS]
1 TBOT = TRECOVERY + DTBOT output
0 RESIZEBOT = NBOT/NCUR output
0 RESIZEANC = NANC/NBOT output
"""

    write_model(model_name, tpl_content, est_content)
    return model_name, 5


def create_two_consecutive_bottlenecks_model():
    """Model 3: two independent bottlenecks."""
    model_name = "two_consecutive_bottlenecks"
    logger.info(f"Creating {model_name} model")

    tpl_content = f"""//Number of population samples (demes)
1
//Population effective sizes (diploid Ne reported by this pipeline)
NCUR
//Sample sizes
{N_CHROMOSOMES}
//Growth rates
0
//Number of migration matrices
0
//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix 
4 historical events
TRECENT_RECOVERY 0 0 0 RESIZE_BOT2 0 0
TRECENT_BOT 0 0 0 RESIZE_INTER 0 0
TOLD_RECOVERY 0 0 0 RESIZE_BOT1 0 0
TOLD_BOT 0 0 0 RESIZE_ANC 0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
FREQ 1 0 {MUTATION_RATE} OUTEXP
"""

    est_content = """// Priors file for two consecutive bottlenecks model
[PARAMETERS]
//#isInt? #name             #dist.#min  #max
1        NCUR               unif  500    100000   output
1        NBOT2              unif  50     50000    output
1        NINTER             unif  100    100000   output
1        NBOT1              unif  50     50000    output
1        NANC               unif  500    200000   output
1        TRECENT_RECOVERY   unif  10     1000     output
1        DTRECENT_BOT       unif  1      3000     output
1        DTOLD_RECOVERY     unif  1      8000     output
1        DTOLD_BOT          unif  1      20000    output

[COMPLEX PARAMETERS]
1 TRECENT_BOT = TRECENT_RECOVERY + DTRECENT_BOT output
1 TOLD_RECOVERY = TRECENT_BOT + DTOLD_RECOVERY output
1 TOLD_BOT = TOLD_RECOVERY + DTOLD_BOT output
0 RESIZE_BOT2 = NBOT2/NCUR output
0 RESIZE_INTER = NINTER/NBOT2 output
0 RESIZE_BOT1 = NBOT1/NINTER output
0 RESIZE_ANC = NANC/NBOT1 output
"""

    write_model(model_name, tpl_content, est_content)
    return model_name, 9


def create_bottleneck_continuous_decline_model():
    """Model 4: bottleneck followed by sustained decline."""
    model_name = "bottleneck_continuous_decline"
    logger.info(f"Creating {model_name} model")

    tpl_content = f"""//Number of population samples (demes)
1
//Population effective sizes (diploid Ne reported by this pipeline)
NCUR
//Sample sizes
{N_CHROMOSOMES}
//Growth rates
0
//Number of migration matrices
0
//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix 
3 historical events
TRECENT 0 0 0 RESIZE_MID 0 0
TBOT 0 0 0 RESIZE_BOT 0 0
TANC 0 0 0 RESIZE_ANC 0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
FREQ 1 0 {MUTATION_RATE} OUTEXP
"""

    est_content = """// Priors file for bottleneck + sustained decline model
[PARAMETERS]
//#isInt? #name      #dist.#min  #max
1        NCUR        unif  200    50000    output
0        RMID        unif  1.01   5.0      output
0        RBOT        unif  1.01   5.0      output
0        RANC        unif  1.01   5.0      output
1        TRECENT     unif  20     2000     output
1        DTBOT       unif  1      8000     output
1        DTANC       unif  1      20000    output

[COMPLEX PARAMETERS]
1 NMID = NCUR*RMID output
1 NBOT = NMID*RBOT output
1 NANC = NBOT*RANC output
1 TBOT = TRECENT + DTBOT output
1 TANC = TBOT + DTANC output
0 RESIZE_MID = NMID/NCUR output
0 RESIZE_BOT = NBOT/NMID output
0 RESIZE_ANC = NANC/NBOT output
"""

    write_model(model_name, tpl_content, est_content)
    return model_name, 7


def create_bottleneck_recent_contraction_model():
    """Model 5: ancient bottleneck + short recovery + recent contraction."""
    model_name = "bottleneck_recent_contraction"
    logger.info(f"Creating {model_name} model")

    tpl_content = f"""//Number of population samples (demes)
1
//Population effective sizes (diploid Ne reported by this pipeline)
NCUR
//Sample sizes
{N_CHROMOSOMES}
//Growth rates
0
//Number of migration matrices
0
//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix 
3 historical events
TRECENT 0 0 0 RESIZE_RECOVER 0 0
TRECOVERY_OLD 0 0 0 RESIZE_BOT 0 0
TBOT_OLD 0 0 0 RESIZE_ANC 0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
FREQ 1 0 {MUTATION_RATE} OUTEXP
"""

    est_content = """// Priors file for bottleneck + recent contraction model
[PARAMETERS]
//#isInt? #name           #dist.#min  #max
1        NCUR             unif  100    50000    output
1        NRECOVER         unif  300    100000   output
1        NBOT             unif  50     40000    output
1        NANC             unif  500    200000   output
1        TRECENT          unif  80     180      output
1        DTRECOVERY_OLD   unif  1      3800     output
1        DTBOT_OLD        unif  1      16000    output

[COMPLEX PARAMETERS]
1 TRECOVERY_OLD = TRECENT + DTRECOVERY_OLD output
1 TBOT_OLD = TRECOVERY_OLD + DTBOT_OLD output
0 RESIZE_RECOVER = NRECOVER/NCUR output
0 RESIZE_BOT = NBOT/NRECOVER output
0 RESIZE_ANC = NANC/NBOT output
"""

    write_model(model_name, tpl_content, est_content)
    return model_name, 7


def create_complex_multi_event_model():
    """Model 6: complex model with >=3 historical size changes."""
    model_name = "complex_multi_event"
    logger.info(f"Creating {model_name} model")

    tpl_content = f"""//Number of population samples (demes)
1
//Population effective sizes (diploid Ne reported by this pipeline)
NCUR
//Sample sizes
{N_CHROMOSOMES}
//Growth rates
0
//Number of migration matrices
0
//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix 
4 historical events
T1 0 0 0 RESIZE1 0 0
T2 0 0 0 RESIZE2 0 0
T3 0 0 0 RESIZE3 0 0
T4 0 0 0 RESIZE4 0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
FREQ 1 0 {MUTATION_RATE} OUTEXP
"""

    est_content = """// Priors file for complex multi-event model
[PARAMETERS]
//#isInt? #name   #dist.#min  #max
1        NCUR     unif  100    50000    output
1        N1       unif  100    100000   output
1        N2       unif  100    120000   output
1        N3       unif  100    150000   output
1        NANC     unif  200    250000   output
1        T1       unif  10     500      output
1        DT2      unif  1      2500     output
1        DT3      unif  1      7000     output
1        DT4      unif  1      20000    output

[COMPLEX PARAMETERS]
1 T2 = T1 + DT2 output
1 T3 = T2 + DT3 output
1 T4 = T3 + DT4 output
0 RESIZE1 = N1/NCUR output
0 RESIZE2 = N2/N1 output
0 RESIZE3 = N3/N2 output
0 RESIZE4 = NANC/N3 output
"""

    write_model(model_name, tpl_content, est_content)
    return model_name, 9


def create_summary_file(models_info):
    """Create a summary file describing all models."""
    summary_path = MODEL_DIR / "MODEL_DESCRIPTIONS.txt"
    
    with open(summary_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("FASTSIMCOAL2 DEMOGRAPHIC MODELS SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("Generated: 2026-01-26\n")
        f.write("Species: Rhinopithecus roxellana (Shennongjia population)\n")
        f.write(f"Sample size: {N_SAMPLES} diploid individuals ({N_CHROMOSOMES} chromosomes)\n")
        f.write("Generation time: 10 years\n")
        f.write("Mutation rate: 1.36×10⁻⁸ per bp per generation\n\n")
        
        f.write("-" * 80 + "\n")
        f.write("DEMOGRAPHIC MODELS\n")
        f.write("-" * 80 + "\n\n")
        
        descriptions = {
            'constant_ne': (
                "Null model with constant population size over time.\n"
                "Parameters: NCUR (current diploid Ne)\n"
                "Use: Baseline for model comparison"
            ),
            'single_bottleneck': (
                "Single bottleneck model: ancestral -> bottleneck -> recovery/current.\n"
                "Parameters: NCUR, NBOT, NANC, TRECOVERY (younger), TBOT (older)\n"
                "Use: Test one severe bottleneck and subsequent recovery"
            ),
            'two_consecutive_bottlenecks': (
                "Two consecutive bottlenecks caused by independent events.\n"
                "Parameters: NCUR, NBOT2, NINTER, NBOT1, NANC, plus 4 ordered times\n"
                "Use: Test repeated shrinkage under multiple climate deteriorations"
            ),
            'bottleneck_continuous_decline': (
                "Bottleneck followed by sustained decline without recovery.\n"
                "Parameters: NCUR, NMID, NBOT, NANC, TRECENT, TBOT, TANC\n"
                "Use: Test post-bottleneck persistent deterioration"
            ),
            'bottleneck_recent_contraction': (
                "Ancient bottleneck + short recovery + recent severe contraction.\n"
                "Parameters: NCUR, NRECOVER, NBOT, NANC, TRECENT, TRECOVERY_OLD, TBOT_OLD\n"
                "Use: Test bottleneck history with recent contraction (~1.3 ka)"
            ),
            'complex_multi_event': (
                "Complex multi-event model with at least 4 time breaks.\n"
                "Parameters: NCUR, N1, N2, N3, NANC, T1, T2, T3, T4\n"
                "Use: Flexible model for >=3 demographic size changes"
            )
        }
        
        for model_name, n_params in models_info:
            f.write(f"Model: {model_name}\n")
            f.write(f"Parameters: {n_params}\n")
            f.write(f"Files: {model_name}.tpl, {model_name}.est\n")
            f.write(f"Description:\n{descriptions.get(model_name, 'No description')}\n")
            f.write("\n" + "-" * 80 + "\n\n")
        
        f.write("TIME CONVERSION\n")
        f.write("-" * 80 + "\n\n")
        f.write("fastsimcoal2 uses units of GENERATIONS.\n")
        f.write("To convert to YEARS: years = generations × 10\n\n")
        f.write("Examples:\n")
        f.write("  TBOT = 130 generations = 1,300 years ago\n")
        f.write("  TBOT = 500 generations = 5,000 years ago\n\n")
        
        f.write("PARAMETER INTERPRETATION\n")
        f.write("-" * 80 + "\n\n")
        f.write("Ne values:\n")
        f.write("  - Reported as diploid effective population size (not census size)\n")
        f.write("  - Typically Ne < Nc (census size)\n")
        f.write("  - Ne/Nc ratio often 0.1-0.5 for mammals\n\n")
        f.write("Likelihood values:\n")
        f.write("  - Use MaxEstLhood for best-run selection and AIC comparison\n")
        f.write("  - Do not use MaxObsLhood to rank independent optimization runs\n\n")
        
        f.write("Growth rates:\n")
        f.write("  - Negative values = population expansion\n")
        f.write("  - Positive values = population decline\n")
        f.write("  - Rate is per generation\n\n")
        
        f.write("=" * 80 + "\n")
    
    logger.info(f"Created summary file: {summary_path}")


def main():
    """Main workflow."""
    logger.info("=" * 80)
    logger.info("CREATING FASTSIMCOAL2 DEMOGRAPHIC MODEL TEMPLATES")
    logger.info("=" * 80)
    logger.info("")
    
    models_info = []
    
    # Create all models
    models_info.append(create_constant_ne_model())
    models_info.append(create_single_bottleneck_model())
    models_info.append(create_two_consecutive_bottlenecks_model())
    models_info.append(create_bottleneck_continuous_decline_model())
    models_info.append(create_bottleneck_recent_contraction_model())
    models_info.append(create_complex_multi_event_model())
    
    # Create summary
    create_summary_file(models_info)
    
    logger.info("")
    logger.info("=" * 80)
    logger.info("MODEL TEMPLATE CREATION COMPLETE")
    logger.info("=" * 80)
    logger.info(f"\nCreated {len(models_info)} demographic models in {MODEL_DIR}")
    logger.info("\nModels:")
    for model_name, n_params in models_info:
        logger.info(f"  - {model_name} ({n_params} parameters)")
    logger.info("")


if __name__ == "__main__":
    main()
