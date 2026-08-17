#!/usr/bin/env python3
"""
Phase 4: ESM-2 Structure-Aware Predictions

Runs mutation effect predictions using ESM-2 via masked language modeling
log-likelihood ratios.

Requires: transformers, torch, pandas
Install: pip install transformers torch pandas
"""

import sys
import os
import pandas as pd
import numpy as np
import torch
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from project_root import get_base_dir

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from utils import setup_logger, ProgressTracker

BASE_DIR = get_base_dir()
OUTPUT_DIR = f"{BASE_DIR}/output/phase4_plm_predictions/esm2"
# Use deduplicated variants (recommended for PLM predictions)
INPUT_VARIANTS = f"{BASE_DIR}/output/phase4_plm_predictions/missense_variants_for_prediction_unique.csv"
PROTEIN_SEQUENCES = f"{BASE_DIR}/output/phase4_plm_predictions/protein_sequences.fasta"

# Model parameters
# Set LOCAL_MODEL_PATH to use a locally downloaded model instead of downloading from HuggingFace
# Option 1: Set directly in code (uncomment and set your path):
# LOCAL_MODEL_PATH = "/path/to/esm2_t33_650M_UR50D"
# Option 2: Set via environment variable (e.g., in shell: export ESM2_LOCAL_PATH="/path/to/model")
# Example: export ESM2_LOCAL_PATH="/path/to/esm2_t33_650M_UR50D"
# If LOCAL_MODEL_PATH is set and exists, it will be used; otherwise falls back to MODEL_NAME
#LOCAL_MODEL_PATH = os.environ.get("ESM2_LOCAL_PATH", None)  # Can be set via environment variable or uncomment line above
LOCAL_MODEL_PATH = os.environ.get("ESM2_LOCAL_PATH")
MODEL_NAME = "facebook/esm2_t33_650M_UR50D"  # HuggingFace model identifier (fallback)
BATCH_SIZE = 32
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
logger = setup_logger("esm2", f"{OUTPUT_DIR}/esm2_predictions.log")


def load_fasta_sequences(fasta_file):
    """
    Load protein sequences from FASTA file
    
    Args:
        fasta_file: Path to FASTA file
        
    Returns:
        dict: Mapping of protein_id to sequence
    """
    sequences = {}
    current_id = None
    current_seq = []
    
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous sequence
                    if current_id is not None:
                        sequences[current_id] = ''.join(current_seq)
                    # Start new sequence
                    current_id = line[1:].split()[0]  # Get ID (first part after >)
                    current_seq = []
                elif line:
                    current_seq.append(line)
            
            # Save last sequence
            if current_id is not None:
                sequences[current_id] = ''.join(current_seq)
        
        logger.info(f"Loaded {len(sequences)} protein sequences from FASTA")
        return sequences
    
    except Exception as e:
        logger.error(f"Error loading FASTA file: {e}")
        return {}


def check_gpu():
    """Check GPU availability"""
    try:
        import torch
        if torch.cuda.is_available():
            logger.info(f"GPU available: {torch.cuda.get_device_name(0)}")
            logger.info(f"GPU memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
            return True
        else:
            logger.warning("GPU not available, using CPU (will be slower)")
            return False
    except ImportError:
        logger.error("PyTorch not installed!")
        return False


def load_esm2_model():
    """
    Load ESM-2 model and tokenizer
    Uses local model path if provided, otherwise downloads from HuggingFace

    Returns:
        tuple: (model, tokenizer, device)
    """
    try:
        from transformers import AutoTokenizer, AutoModelForMaskedLM
        import torch

        # Determine which model path to use
        if LOCAL_MODEL_PATH and os.path.exists(LOCAL_MODEL_PATH):
            model_path = LOCAL_MODEL_PATH
            logger.info(f"Using local ESM-2 model from: {model_path}")
        else:
            model_path = MODEL_NAME
            if LOCAL_MODEL_PATH:
                logger.warning(f"Local model path specified but not found: {LOCAL_MODEL_PATH}")
                logger.info(f"Falling back to HuggingFace model: {MODEL_NAME}")
            else:
                logger.info(f"Loading ESM-2 model from HuggingFace: {MODEL_NAME}")
            logger.info("This may take several minutes for first-time download...")

        tokenizer = AutoTokenizer.from_pretrained(model_path)
        model = AutoModelForMaskedLM.from_pretrained(model_path)

        device = torch.device(DEVICE)
        model = model.to(device)
        model.eval()

        logger.info(f"Model loaded successfully on {device}")
        return model, tokenizer, device

    except Exception as e:
        logger.error(f"Error loading ESM-2 model: {e}")
        logger.error("Please install: pip install transformers torch")
        return None, None, None


def calculate_log_likelihood_ratio(model, tokenizer, device, sequence, position, wt_aa, mut_aa):
    """
    Calculate log-likelihood ratio for a single mutation using ESM-2

    Args:
        model: ESM-2 model
        tokenizer: Tokenizer
        device: Device (cuda/cpu)
        sequence: Wild-type protein sequence
        position: Position of mutation (1-indexed)
        wt_aa: Wild-type amino acid
        mut_aa: Mutant amino acid

    Returns:
        float: Log-likelihood ratio (negative values = deleterious)
    """
    import torch

    try:
        # Tokenize sequence
        inputs = tokenizer(sequence, return_tensors="pt")
        inputs = {k: v.to(device) for k, v in inputs.items()}

        # Forward pass
        with torch.no_grad():
            outputs = model(**inputs)
            logits = outputs.logits

        # ESM tokenizers include a BOS token; residue index maps to logits[:, position, :]
        position_idx = position
        position_logits = logits[0, position_idx, :]

        probs = torch.nn.functional.softmax(position_logits, dim=0)

        wt_token = tokenizer.convert_tokens_to_ids(wt_aa)
        mut_token = tokenizer.convert_tokens_to_ids(mut_aa)

        wt_prob = probs[wt_token].item()
        mut_prob = probs[mut_token].item()

        llr = np.log(mut_prob) - np.log(wt_prob)
        return llr

    except Exception as e:
        logger.warning(f"Error calculating LLR for {wt_aa}{position}{mut_aa}: {e}")
        return np.nan


def predict_batch(model, tokenizer, device, variants_batch):
    """
    Predict scores for a batch of variants
    """
    scores = []
    for var in variants_batch:
        # Get sequence from variant record
        sequence = var.get('wt_sequence', '')
        if not sequence:
            logger.warning(f"No sequence found for variant {var.get('variant_id', 'unknown')}")
            scores.append(np.nan)
            continue
            
        score = calculate_log_likelihood_ratio(
            model, tokenizer, device,
            sequence,
            var['aa_pos'],
            var['wt_aa'],
            var['mut_aa']
        )
        scores.append(score)
    return scores

def main():
    logger.info("="*70)
    logger.info("ESM-2 MUTATION EFFECT PREDICTIONS")
    logger.info("="*70)
    logger.info("")

    # Check GPU
    _ = check_gpu()
    logger.info("")

    # Load model
    logger.info("Loading ESM-2 model...")
    model, tokenizer, device = load_esm2_model()

    if model is None:
        logger.error(
            "Failed to load ESM-2 model. Install torch/transformers and "
            "optionally set ESM2_LOCAL_PATH; refusing to write placeholder scores."
        )
        return 1

    logger.info("")

    # Load protein sequences
    logger.info("Loading protein sequences from FASTA...")
    protein_sequences = load_fasta_sequences(PROTEIN_SEQUENCES)
    if not protein_sequences:
        logger.error("Failed to load protein sequences!")
        return 1
    logger.info("")

    # Load variants
    logger.info("Loading variants for ESM-2 prediction...")
    variants_df = pd.read_csv(INPUT_VARIANTS)
    logger.info(f"Loaded {len(variants_df)} missense variants")
    logger.info("")

    # Map protein sequences to variants
    logger.info("Mapping sequences to variants...")
    if 'matched_protein_id' in variants_df.columns:
        variants_df['wt_sequence'] = variants_df['matched_protein_id'].map(protein_sequences)
        missing_sequences = variants_df['wt_sequence'].isna().sum()
        if missing_sequences > 0:
            logger.warning(f"{missing_sequences} variants missing protein sequences")
        logger.info(f"Successfully mapped sequences for {len(variants_df) - missing_sequences} variants")
    else:
        logger.error("Column 'matched_protein_id' not found in variants file!")
        return 1
    logger.info("")

    logger.info("Running ESM-2 predictions...")
    logger.info(f"Batch size: {BATCH_SIZE}")
    logger.info(f"Device: {device}")
    logger.info("")

    variants = variants_df.to_dict('records')
    esm2_scores = []

    progress = ProgressTracker(len(variants), "ESM-2 predictions", update_interval=10)

    for i in range(0, len(variants), BATCH_SIZE):
        batch = variants[i:i+BATCH_SIZE]
        batch_scores = predict_batch(model, tokenizer, device, batch)
        esm2_scores.extend(batch_scores)
        progress.update(len(batch))

    progress.finish()

    variants_df['esm2_score'] = esm2_scores
    variants_df['esm2_percentile'] = variants_df['esm2_score'].rank(pct=True)
    variants_df['esm2_prediction'] = variants_df['esm2_score'].apply(
        lambda x: 'deleterious' if x < -2 else 'possibly_deleterious' if x < 0 else 'benign'
    )

    logger.info("")
    logger.info("Saving ESM-2 results...")
    output_file = f"{OUTPUT_DIR}/esm2_predictions.csv"
    variants_df.to_csv(output_file, index=False)
    logger.info(f"Predictions saved: {output_file}")

    # Brief summary
    valid_scores = variants_df[~variants_df['esm2_score'].isna()]
    if len(valid_scores) > 0:
        logger.info(f"Valid predictions: {len(valid_scores)} / {len(variants_df)}")
        logger.info(f"Mean: {valid_scores['esm2_score'].mean():.4f} | Std: {valid_scores['esm2_score'].std():.4f}")

    logger.info("")
    logger.info("="*70)
    logger.info("ESM-2 PREDICTIONS COMPLETE")
    logger.info("="*70)

    return 0

if __name__ == "__main__":
    sys.exit(main())


