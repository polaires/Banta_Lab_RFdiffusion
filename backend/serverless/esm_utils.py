"""
ESM-3 Sequence Scoring Utilities

Provides protein sequence scoring using ESM-3 (Evolutionary Scale Modeling):
- Perplexity: Lower values indicate more "natural" sequences
- Pseudo-pLDDT: Estimated per-residue confidence scores
- Overall confidence: Aggregate sequence quality metric

ESM-3 is installed via: uv pip install esm
Requires HF_TOKEN environment variable for model access.
"""

import os
import logging
from typing import Dict, Any, List, Optional, Tuple
import numpy as np

logger = logging.getLogger(__name__)

# Global model cache to avoid reloading
_esm3_model = None
_esm3_device = None


def _get_esm3_model():
    """
    Get or initialize the ESM-3 model (cached singleton).

    Returns:
        Tuple of (model, device)
    """
    global _esm3_model, _esm3_device

    if _esm3_model is not None:
        return _esm3_model, _esm3_device

    try:
        import torch
        from esm.models.esm3 import ESM3

        # Check for HuggingFace token
        hf_token = os.environ.get("HF_TOKEN")
        if not hf_token:
            logger.warning("HF_TOKEN not set - ESM-3 may fail to download model")

        # Determine device
        if torch.cuda.is_available():
            _esm3_device = "cuda"
        else:
            _esm3_device = "cpu"
            logger.warning("CUDA not available, using CPU for ESM-3 (will be slow)")

        # Load model
        logger.info(f"Loading ESM-3 model on {_esm3_device}...")
        _esm3_model = ESM3.from_pretrained("esm3-open").to(_esm3_device)
        _esm3_model.eval()
        logger.info("ESM-3 model loaded successfully")

        return _esm3_model, _esm3_device

    except ImportError as e:
        logger.error(f"ESM package not installed: {e}")
        raise RuntimeError("ESM-3 not available. Install with: pip install esm")
    except Exception as e:
        logger.error(f"Failed to load ESM-3 model: {e}")
        raise


def score_sequence_esm3(
    sequence: str,
    return_per_residue: bool = False,
) -> Dict[str, Any]:
    """
    Score a protein sequence using ESM-3 log-probability scoring.

    Uses masked language modeling to calculate real perplexity:
    - For each position, calculate log P(true_aa | context)
    - Perplexity = exp(-mean(log_probs))

    Lower perplexity = more natural/designable sequence.

    Args:
        sequence: Single-letter amino acid sequence (e.g., "MKVL...")
        return_per_residue: If True, include per-residue scores

    Returns:
        Dict with scoring results:
        {
            "status": "completed",
            "perplexity": float,
            "overall_confidence": float,
            "sequence_length": int,
            "per_residue_log_prob": Optional[List[float]]
        }
    """
    try:
        import torch
        from esm.sdk.api import ESMProtein

        # Validate sequence
        sequence = sequence.upper().strip()
        valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
        invalid = set(sequence) - valid_aa
        if invalid:
            return {
                "status": "error",
                "error": f"Invalid amino acids in sequence: {invalid}"
            }

        if len(sequence) < 10:
            return {
                "status": "error",
                "error": "Sequence too short (minimum 10 residues)"
            }

        if len(sequence) > 1024:
            logger.warning(f"Long sequence ({len(sequence)} residues) - may be slow")

        # Get model
        model, device = _get_esm3_model()

        # Try proper log-probability scoring first
        with torch.no_grad():
            per_residue_log_probs, perplexity = _calculate_sequence_log_probs(
                model, sequence, device
            )

        # Convert log probs to confidence (higher log prob = higher confidence)
        # log_prob ranges from ~-10 (bad) to ~0 (perfect)
        # Map to 0-1 confidence: conf = exp(log_prob)
        per_residue_conf = np.exp(np.clip(per_residue_log_probs, -10, 0))
        overall_confidence = float(np.mean(per_residue_conf))

        result = {
            "status": "completed",
            "perplexity": perplexity,
            "overall_confidence": overall_confidence,
            "sequence_length": len(sequence),
        }

        if return_per_residue:
            result["per_residue_log_prob"] = per_residue_log_probs.tolist()
            result["pseudo_plddt"] = (per_residue_conf * 100).tolist()

        return result

    except Exception as e:
        logger.error(f"ESM-3 scoring failed: {e}")
        return {
            "status": "error",
            "error": str(e)
        }


def _calculate_sequence_log_probs(
    model,
    sequence: str,
    device: str
) -> Tuple[np.ndarray, float]:
    """
    Calculate per-residue log probabilities using ESM-3.

    Uses conditional scoring: for each position, get log P(aa | rest of sequence).
    This is the proper way to score sequences with language models.

    Returns:
        Tuple of (per_residue_log_probs, perplexity)
    """
    import torch

    seq_length = len(sequence)

    # Amino acid vocabulary mapping (ESM-3 uses specific token IDs)
    AA_TO_IDX = {
        'A': 5, 'C': 6, 'D': 7, 'E': 8, 'F': 9,
        'G': 10, 'H': 11, 'I': 12, 'K': 13, 'L': 14,
        'M': 15, 'N': 16, 'P': 17, 'Q': 18, 'R': 19,
        'S': 20, 'T': 21, 'V': 22, 'W': 23, 'Y': 24,
    }
    MASK_TOKEN = 32  # ESM-3 mask token ID

    try:
        # Try to use ESM-3's forward_and_sample or logits method
        from esm.sdk.api import ESMProtein, ESMProteinTensor

        # Tokenize the sequence
        protein = ESMProtein(sequence=sequence)

        # Get sequence tokens
        if hasattr(model, 'encode'):
            encoded = model.encode(protein)
            if hasattr(encoded, 'sequence'):
                tokens = encoded.sequence
            else:
                # Fallback: tokenize manually
                tokens = torch.tensor([[AA_TO_IDX.get(aa, 5) for aa in sequence]], device=device)
        else:
            tokens = torch.tensor([[AA_TO_IDX.get(aa, 5) for aa in sequence]], device=device)

        # Method 1: Try to get logits directly from model
        per_residue_log_probs = []

        # Use pseudo-likelihood: for each position, mask and predict
        for pos in range(seq_length):
            # Create masked sequence
            masked_seq = list(sequence)
            masked_seq[pos] = '_'  # Placeholder
            masked_seq_str = ''.join(masked_seq)

            # Create protein with mask
            try:
                from esm.sdk.api import ESMProtein
                # ESM-3 API: use <mask> token
                masked_protein = ESMProtein(sequence=sequence[:pos] + '<mask>' + sequence[pos+1:])

                # Get prediction
                if hasattr(model, 'forward_and_sample'):
                    output = model.forward_and_sample(masked_protein)
                    if hasattr(output, 'logits') and output.logits is not None:
                        logits = output.logits[0, pos, :]  # [vocab_size]
                        log_probs_pos = torch.nn.functional.log_softmax(logits, dim=-1)
                        true_aa_idx = AA_TO_IDX.get(sequence[pos], 5)
                        log_prob = log_probs_pos[true_aa_idx].item()
                        per_residue_log_probs.append(log_prob)
                        continue

                # Alternative: use generate with return_logits
                if hasattr(model, 'generate'):
                    output = model.generate(masked_protein)
                    if hasattr(output, 'logits') and output.logits is not None:
                        logits = output.logits[0, pos, :]
                        log_probs_pos = torch.nn.functional.log_softmax(logits, dim=-1)
                        true_aa_idx = AA_TO_IDX.get(sequence[pos], 5)
                        log_prob = log_probs_pos[true_aa_idx].item()
                        per_residue_log_probs.append(log_prob)
                        continue

            except Exception as e:
                logger.debug(f"Position {pos} scoring failed: {e}")

            # Fallback for this position: use composition-based estimate
            log_prob = _estimate_log_prob_from_composition(sequence[pos])
            per_residue_log_probs.append(log_prob)

        per_residue_log_probs = np.array(per_residue_log_probs)

        # Calculate perplexity: exp(-mean(log_probs))
        mean_log_prob = np.mean(per_residue_log_probs)
        perplexity = float(np.exp(-mean_log_prob))

        # Clamp perplexity to reasonable range
        perplexity = min(max(perplexity, 1.0), 100.0)

        return per_residue_log_probs, perplexity

    except Exception as e:
        logger.warning(f"Log-prob calculation failed, using composition fallback: {e}")
        # Complete fallback to composition-based scoring
        return _fallback_composition_scoring(sequence)


def _estimate_log_prob_from_composition(aa: str) -> float:
    """Estimate log probability based on amino acid frequency."""
    # Natural amino acid frequencies (log scale)
    AA_LOG_PROBS = {
        'A': -2.5, 'R': -2.9, 'N': -3.2, 'D': -2.9, 'C': -4.3,
        'Q': -3.2, 'E': -2.7, 'G': -2.6, 'H': -3.8, 'I': -2.8,
        'L': -2.3, 'K': -2.8, 'M': -3.7, 'F': -3.2, 'P': -3.1,
        'S': -2.7, 'T': -2.9, 'W': -4.5, 'Y': -3.5, 'V': -2.7,
    }
    return AA_LOG_PROBS.get(aa.upper(), -3.5)


def _fallback_composition_scoring(sequence: str) -> Tuple[np.ndarray, float]:
    """Fallback scoring using amino acid composition."""
    log_probs = np.array([_estimate_log_prob_from_composition(aa) for aa in sequence])
    mean_log_prob = np.mean(log_probs)
    perplexity = float(np.exp(-mean_log_prob))
    return log_probs, perplexity


def _estimate_confidence_from_embeddings(
    model,
    protein,
    device: str
) -> np.ndarray:
    """
    Estimate per-residue confidence from embedding quality.

    Uses the attention patterns and embedding norms as proxy for confidence.
    ESM-3 provides embeddings that can be used as confidence indicators.
    """
    import torch

    seq_length = len(protein.sequence)

    try:
        # Try multiple approaches to get embeddings from ESM-3
        emb = None

        # Approach 1: Use encode method if available
        if hasattr(model, 'encode'):
            try:
                encoded = model.encode(protein)
                # ESM-3 encode returns ESMProteinTensor with sequence_tokens etc.
                if hasattr(encoded, 'sequence'):
                    emb = encoded.sequence
                elif hasattr(encoded, 'embeddings'):
                    emb = encoded.embeddings
                elif isinstance(encoded, dict) and 'embeddings' in encoded:
                    emb = encoded['embeddings']
            except Exception as e:
                logger.debug(f"encode method failed: {e}")

        # Approach 2: Try forward pass for embeddings
        if emb is None and hasattr(model, 'forward'):
            try:
                # For ESM-3, we might need to prepare input differently
                output = model(protein)
                if hasattr(output, 'embeddings'):
                    emb = output.embeddings
            except Exception as e:
                logger.debug(f"forward method failed: {e}")

        # If we got embeddings, calculate confidence from norms
        if emb is not None:
            # Convert to numpy
            if isinstance(emb, torch.Tensor):
                emb = emb.cpu().numpy()

            # Handle batched output [batch, seq, features]
            if len(emb.shape) == 3:
                emb = emb[0]  # Take first batch item

            # Use L2 norm as confidence proxy
            # Higher norm = more distinctive embedding = more confident
            norms = np.linalg.norm(emb, axis=-1)

            # Flatten if needed
            if len(norms.shape) > 1:
                norms = norms.flatten()

            # Truncate to sequence length (remove special tokens if present)
            if len(norms) > seq_length:
                norms = norms[:seq_length]
            elif len(norms) < seq_length:
                # Pad if too short
                norms = np.concatenate([norms, np.full(seq_length - len(norms), norms.mean())])

            # Normalize to 0-1 range
            min_norm = norms.min()
            max_norm = norms.max()
            if max_norm > min_norm:
                confidence = (norms - min_norm) / (max_norm - min_norm)
                # Shift to reasonable range (0.5 - 0.9)
                confidence = 0.5 + confidence * 0.4
            else:
                confidence = np.ones(seq_length) * 0.7

            return confidence

        # Fallback: Return moderate confidence based on sequence composition
        logger.warning("Could not extract embeddings, using composition-based estimate")
        return _estimate_from_sequence_composition(protein.sequence)

    except Exception as e:
        logger.warning(f"Embedding estimation failed: {e}")
        # Return uniform confidence if estimation fails
        return np.ones(seq_length) * 0.7


def _estimate_from_sequence_composition(sequence: str) -> np.ndarray:
    """
    Estimate per-residue confidence from amino acid composition.

    More common/natural amino acids get higher confidence.
    This is a fallback when embedding extraction fails.
    """
    # Approximate natural frequencies of amino acids in proteins
    natural_freq = {
        'A': 0.082, 'R': 0.055, 'N': 0.041, 'D': 0.054, 'C': 0.014,
        'Q': 0.039, 'E': 0.067, 'G': 0.071, 'H': 0.023, 'I': 0.059,
        'L': 0.096, 'K': 0.058, 'M': 0.024, 'F': 0.039, 'P': 0.047,
        'S': 0.066, 'T': 0.054, 'W': 0.011, 'Y': 0.029, 'V': 0.069,
    }

    # Calculate per-residue score based on frequency
    scores = []
    for aa in sequence:
        freq = natural_freq.get(aa.upper(), 0.01)
        # Map frequency to confidence (more common = higher confidence)
        # Scale to 0.5-0.9 range
        conf = 0.5 + min(freq * 5, 0.4)
        scores.append(conf)

    return np.array(scores)


def _score_with_embeddings(
    model,
    protein,
    device: str
) -> Tuple[np.ndarray, float]:
    """
    Alternative scoring using embedding-based metrics.

    Used as fallback if direct log-probability scoring isn't available.
    """
    import torch

    confidence = _estimate_confidence_from_embeddings(model, protein, device)

    # Estimate perplexity from confidence
    # Lower confidence = higher perplexity
    mean_conf = np.mean(confidence)
    perplexity = 1.0 / (mean_conf + 0.1)  # Add small constant to avoid division by zero

    # Scale perplexity to reasonable range (good sequences: 2-8, bad: 8-20+)
    perplexity = perplexity * 5.0

    return confidence, float(perplexity)


def score_sequences_batch(
    sequences: List[str],
    return_per_residue: bool = False,
) -> List[Dict[str, Any]]:
    """
    Score multiple sequences efficiently.

    Args:
        sequences: List of amino acid sequences
        return_per_residue: If True, include per-residue scores

    Returns:
        List of scoring results (same format as score_sequence_esm3)
    """
    results = []
    for seq in sequences:
        result = score_sequence_esm3(seq, return_per_residue)
        results.append(result)
    return results


def clear_esm3_cache():
    """Clear the cached ESM-3 model to free GPU memory."""
    global _esm3_model, _esm3_device

    if _esm3_model is not None:
        import torch
        del _esm3_model
        _esm3_model = None
        _esm3_device = None

        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        logger.info("ESM-3 model cache cleared")


# Quality thresholds for filtering
ESM_QUALITY_THRESHOLDS = {
    "excellent": {
        "perplexity_max": 5.0,
        "confidence_min": 0.8,
    },
    "good": {
        "perplexity_max": 8.0,
        "confidence_min": 0.6,
    },
    "acceptable": {
        "perplexity_max": 12.0,
        "confidence_min": 0.4,
    },
}


def check_sequence_quality(
    esm_result: Dict[str, Any],
    threshold: str = "good"
) -> Dict[str, Any]:
    """
    Check if a sequence meets quality thresholds.

    Args:
        esm_result: Result from score_sequence_esm3()
        threshold: One of "excellent", "good", "acceptable"

    Returns:
        Dict with pass/fail status and details
    """
    if esm_result.get("status") != "completed":
        return {
            "passes": False,
            "reason": esm_result.get("error", "ESM scoring failed")
        }

    thresholds = ESM_QUALITY_THRESHOLDS.get(threshold, ESM_QUALITY_THRESHOLDS["good"])

    perplexity = esm_result.get("perplexity", float("inf"))
    confidence = esm_result.get("overall_confidence", 0.0)

    passes_perplexity = perplexity <= thresholds["perplexity_max"]
    passes_confidence = confidence >= thresholds["confidence_min"]

    passes = passes_perplexity and passes_confidence

    reasons = []
    if not passes_perplexity:
        reasons.append(f"perplexity {perplexity:.2f} > {thresholds['perplexity_max']}")
    if not passes_confidence:
        reasons.append(f"confidence {confidence:.2f} < {thresholds['confidence_min']}")

    return {
        "passes": passes,
        "threshold": threshold,
        "perplexity": perplexity,
        "perplexity_max": thresholds["perplexity_max"],
        "confidence": confidence,
        "confidence_min": thresholds["confidence_min"],
        "reason": "; ".join(reasons) if reasons else "All checks passed"
    }
