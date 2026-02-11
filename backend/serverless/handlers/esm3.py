"""ESM3 protein language model handlers.

Tasks: esm3_score, esm3_generate, esm3_embed
"""

from typing import Dict, Any

from inference_utils import (
    esm3_score_sequences,
    esm3_generate_sequence,
    esm3_get_embeddings,
)


def handle_esm3_score(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Score protein sequences using ESM3 perplexity.

    Input:
        sequences: list[str] - List of amino acid sequences to score

    Returns:
        scores: list of {sequence, perplexity, score}
        Lower perplexity = better sequence quality
    """
    sequences = job_input.get("sequences")
    if not sequences:
        return {"status": "failed", "error": "Missing 'sequences' parameter"}

    if not isinstance(sequences, list):
        return {"status": "failed", "error": "'sequences' must be a list of strings"}

    return esm3_score_sequences(job_input)


def handle_esm3_generate(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Generate protein sequences using ESM3 with optional function conditioning.

    Input:
        prompt: str (optional) - Partial sequence to complete
        functions: list[str] (optional) - Function keywords like ["zinc-binding", "hydrolase"]
        num_sequences: int (default: 4) - Number of sequences to generate
        temperature: float (default: 0.7) - Sampling temperature
        max_length: int (default: 200) - Maximum sequence length

    Returns:
        sequences: list[str] - Generated protein sequences
    """
    return esm3_generate_sequence(job_input)


def handle_esm3_embed(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Get ESM3 embeddings for protein sequences.

    Input:
        sequences: list[str] - List of amino acid sequences

    Returns:
        embeddings: list of {sequence, per_residue, global}
        per_residue: [seq_len, hidden_dim] per-position embeddings
        global: [hidden_dim] mean-pooled sequence embedding
    """
    sequences = job_input.get("sequences")
    if not sequences:
        return {"status": "failed", "error": "Missing 'sequences' parameter"}

    if not isinstance(sequences, list):
        return {"status": "failed", "error": "'sequences' must be a list of strings"}

    return esm3_get_embeddings(job_input)
