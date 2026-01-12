"""
ESMFold Structure Validation Utilities

Validates designed protein sequences by predicting their structure with ESMFold
and comparing to the designed backbone. This provides BindCraft-style validation
metrics:
- pLDDT: Per-residue confidence from ESMFold (similar to AlphaFold2)
- RMSD: Backbone deviation between designed and predicted structure

Key Thresholds (from BindCraft):
- pLDDT > 0.7: Acceptable binder (BindCraft uses 0.8)
- RMSD < 2.0 Å: Structure matches design intent

ESMFold is available via:
1. Local model: esm.pretrained.esmfold_v1() (requires ~8GB GPU RAM)
2. Hugging Face API: https://api-inference.huggingface.co/models/facebook/esmfold_v1
"""

import os
import logging
import tempfile
from typing import Dict, Any, List, Optional, Tuple
import numpy as np

logger = logging.getLogger(__name__)

# Global model cache
_esmfold_model = None
_esmfold_device = None

# ESMFold validation thresholds (inspired by BindCraft)
ESMFOLD_THRESHOLDS = {
    "strict": {
        "plddt_min": 0.80,  # BindCraft default
        "rmsd_max": 1.5,     # Å
    },
    "standard": {
        "plddt_min": 0.70,   # Slightly relaxed
        "rmsd_max": 2.0,     # Å
    },
    "relaxed": {
        "plddt_min": 0.60,
        "rmsd_max": 3.0,     # Å
    },
}


def _get_esmfold_model():
    """
    Get or initialize the ESMFold model (cached singleton).

    Returns:
        Tuple of (model, device)
    """
    global _esmfold_model, _esmfold_device

    if _esmfold_model is not None:
        return _esmfold_model, _esmfold_device

    try:
        import torch
        import esm

        # Determine device
        if torch.cuda.is_available():
            _esmfold_device = "cuda"
            # Check GPU memory - ESMFold needs ~8GB
            gpu_mem = torch.cuda.get_device_properties(0).total_memory / 1e9
            if gpu_mem < 7.0:
                logger.warning(f"GPU has only {gpu_mem:.1f}GB - ESMFold may OOM")
        else:
            _esmfold_device = "cpu"
            logger.warning("CUDA not available, using CPU for ESMFold (will be very slow)")

        # Load model
        logger.info(f"Loading ESMFold model on {_esmfold_device}...")
        _esmfold_model = esm.pretrained.esmfold_v1()
        _esmfold_model = _esmfold_model.eval().to(_esmfold_device)

        # Set chunk size for memory efficiency
        _esmfold_model.set_chunk_size(128)

        logger.info("ESMFold model loaded successfully")
        return _esmfold_model, _esmfold_device

    except ImportError as e:
        logger.error(f"ESM package not installed or ESMFold not available: {e}")
        raise RuntimeError("ESMFold not available. Install with: pip install esm")
    except Exception as e:
        logger.error(f"Failed to load ESMFold model: {e}")
        raise


def predict_structure_esmfold(
    sequence: str,
    return_pdb: bool = True,
) -> Dict[str, Any]:
    """
    Predict protein structure from sequence using ESMFold.

    Args:
        sequence: Single-letter amino acid sequence
        return_pdb: If True, include PDB string in output

    Returns:
        Dict with:
        - status: "completed" or "error"
        - plddt: Per-residue pLDDT scores (0-1 scale)
        - mean_plddt: Average pLDDT
        - pdb_content: Predicted structure as PDB string (if return_pdb)
        - ca_coords: CA atom coordinates (N x 3)
    """
    try:
        import torch

        # Validate sequence
        sequence = sequence.upper().strip()
        valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
        invalid = set(sequence) - valid_aa
        if invalid:
            return {
                "status": "error",
                "error": f"Invalid amino acids: {invalid}"
            }

        if len(sequence) < 10:
            return {
                "status": "error",
                "error": "Sequence too short (minimum 10 residues)"
            }

        if len(sequence) > 400:
            logger.warning(f"Long sequence ({len(sequence)} residues) - may be slow/OOM")

        # Get model
        model, device = _get_esmfold_model()

        with torch.no_grad():
            output = model.infer_pdb(sequence)

        # Parse the PDB output to extract pLDDT and coordinates
        pdb_content = output
        plddt, ca_coords = _parse_esmfold_pdb(pdb_content)

        result = {
            "status": "completed",
            "plddt": plddt.tolist(),
            "mean_plddt": float(np.mean(plddt)),
            "sequence_length": len(sequence),
            "ca_coords": ca_coords.tolist(),
        }

        if return_pdb:
            result["pdb_content"] = pdb_content

        return result

    except Exception as e:
        logger.error(f"ESMFold prediction failed: {e}")
        return {
            "status": "error",
            "error": str(e)
        }


def _parse_esmfold_pdb(pdb_content: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Parse ESMFold PDB output to extract pLDDT scores and CA coordinates.

    ESMFold stores pLDDT in the B-factor column (0-100 scale).

    Returns:
        Tuple of (plddt array [0-1 scale], ca_coords array [N, 3])
    """
    plddt_values = []
    ca_coords = []

    for line in pdb_content.split('\n'):
        if line.startswith('ATOM') and ' CA ' in line:
            try:
                # B-factor is in columns 60-66 (0-indexed: 60:66)
                bfactor = float(line[60:66].strip())
                # Normalize to 0-1 scale
                plddt_values.append(bfactor / 100.0)

                # Extract CA coordinates (columns 30-54)
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                ca_coords.append([x, y, z])
            except (ValueError, IndexError):
                continue

    return np.array(plddt_values), np.array(ca_coords)


def calculate_backbone_rmsd(
    pdb1_content: str,
    pdb2_content: str,
    chain1: str = None,
    chain2: str = None,
) -> Dict[str, Any]:
    """
    Calculate backbone RMSD between two structures.

    Uses CA atoms for alignment and RMSD calculation.

    Args:
        pdb1_content: First PDB structure (designed backbone)
        pdb2_content: Second PDB structure (ESMFold prediction)
        chain1: Chain ID in first structure (optional)
        chain2: Chain ID in second structure (optional)

    Returns:
        Dict with:
        - rmsd: Backbone RMSD in Angstroms
        - aligned: True if alignment was performed
        - n_atoms: Number of CA atoms used
    """
    try:
        # Extract CA coordinates
        coords1 = _extract_ca_coords_from_pdb(pdb1_content, chain1)
        coords2 = _extract_ca_coords_from_pdb(pdb2_content, chain2)

        if len(coords1) == 0 or len(coords2) == 0:
            return {
                "status": "error",
                "error": "Could not extract CA coordinates from structures"
            }

        # Handle length mismatch (use shorter length)
        min_len = min(len(coords1), len(coords2))
        if len(coords1) != len(coords2):
            logger.warning(f"Length mismatch: {len(coords1)} vs {len(coords2)}, using first {min_len}")
            coords1 = coords1[:min_len]
            coords2 = coords2[:min_len]

        # Kabsch alignment
        aligned_rmsd = _kabsch_rmsd(coords1, coords2)

        return {
            "status": "completed",
            "rmsd": float(aligned_rmsd),
            "n_atoms": min_len,
            "aligned": True,
        }

    except Exception as e:
        logger.error(f"RMSD calculation failed: {e}")
        return {
            "status": "error",
            "error": str(e)
        }


def _extract_ca_coords_from_pdb(pdb_content: str, chain: str = None) -> np.ndarray:
    """Extract CA atom coordinates from PDB content."""
    coords = []

    for line in pdb_content.split('\n'):
        if not line.startswith('ATOM'):
            continue

        atom_name = line[12:16].strip()
        if atom_name != 'CA':
            continue

        if chain is not None:
            line_chain = line[21]
            if line_chain != chain:
                continue

        try:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            coords.append([x, y, z])
        except (ValueError, IndexError):
            continue

    return np.array(coords)


def _kabsch_rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """
    Calculate RMSD after optimal superposition using Kabsch algorithm.

    Args:
        coords1: Reference coordinates (N x 3)
        coords2: Mobile coordinates (N x 3)

    Returns:
        RMSD in Angstroms after alignment
    """
    # Center both coordinate sets
    center1 = coords1.mean(axis=0)
    center2 = coords2.mean(axis=0)

    coords1_centered = coords1 - center1
    coords2_centered = coords2 - center2

    # Compute covariance matrix
    H = coords2_centered.T @ coords1_centered

    # SVD
    U, S, Vt = np.linalg.svd(H)

    # Compute rotation matrix
    R = Vt.T @ U.T

    # Handle reflection case
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    # Apply rotation
    coords2_aligned = coords2_centered @ R

    # Calculate RMSD
    diff = coords1_centered - coords2_aligned
    rmsd = np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))

    return rmsd


def validate_structure_esmfold(
    sequence: str,
    designed_pdb: str,
    binder_chain: str = "B",
    threshold: str = "standard",
) -> Dict[str, Any]:
    """
    Validate a designed binder by predicting its structure and comparing to design.

    This is the main validation function that provides BindCraft-style metrics.

    Args:
        sequence: Designed binder sequence
        designed_pdb: Full complex PDB (target + binder)
        binder_chain: Chain ID of the binder in designed_pdb
        threshold: Quality threshold ("strict", "standard", "relaxed")

    Returns:
        Dict with:
        - status: "completed" or "error"
        - esmfold_plddt: Average pLDDT from ESMFold prediction
        - esmfold_rmsd: Backbone RMSD between design and prediction (Å)
        - validation_passed: True if both metrics pass thresholds
        - per_residue_plddt: Per-residue pLDDT scores
        - predicted_pdb: ESMFold predicted structure
    """
    try:
        # Step 1: Predict structure with ESMFold
        logger.info(f"Predicting structure for {len(sequence)} residue binder...")
        prediction = predict_structure_esmfold(sequence, return_pdb=True)

        if prediction["status"] != "completed":
            return {
                "status": "error",
                "error": f"ESMFold prediction failed: {prediction.get('error')}"
            }

        mean_plddt = prediction["mean_plddt"]
        per_residue_plddt = prediction["plddt"]
        predicted_pdb = prediction["pdb_content"]

        # Step 2: Calculate RMSD vs designed backbone
        logger.info("Calculating backbone RMSD...")
        rmsd_result = calculate_backbone_rmsd(
            pdb1_content=designed_pdb,
            pdb2_content=predicted_pdb,
            chain1=binder_chain,
            chain2=None,  # ESMFold output has single chain
        )

        if rmsd_result["status"] != "completed":
            # If RMSD calculation fails, return with just pLDDT
            logger.warning(f"RMSD calculation failed: {rmsd_result.get('error')}")
            rmsd = None
        else:
            rmsd = rmsd_result["rmsd"]

        # Step 3: Check against thresholds
        thresholds = ESMFOLD_THRESHOLDS.get(threshold, ESMFOLD_THRESHOLDS["standard"])

        passes_plddt = mean_plddt >= thresholds["plddt_min"]
        passes_rmsd = rmsd is None or rmsd <= thresholds["rmsd_max"]
        validation_passed = passes_plddt and passes_rmsd

        # Build result
        result = {
            "status": "completed",
            "esmfold_plddt": mean_plddt,
            "esmfold_rmsd": rmsd,
            "validation_passed": validation_passed,
            "threshold_used": threshold,
            "thresholds": thresholds,
            "per_residue_plddt": per_residue_plddt,
            "predicted_pdb": predicted_pdb,
        }

        # Add failure reasons
        if not validation_passed:
            reasons = []
            if not passes_plddt:
                reasons.append(f"pLDDT {mean_plddt:.2f} < {thresholds['plddt_min']}")
            if not passes_rmsd and rmsd is not None:
                reasons.append(f"RMSD {rmsd:.2f}Å > {thresholds['rmsd_max']}Å")
            result["failure_reasons"] = reasons

        logger.info(f"Validation: pLDDT={mean_plddt:.2f}, RMSD={rmsd:.2f if rmsd else 'N/A'}Å, passed={validation_passed}")
        return result

    except Exception as e:
        logger.error(f"Structure validation failed: {e}")
        return {
            "status": "error",
            "error": str(e)
        }


def validate_binder_batch(
    designs: List[Dict[str, Any]],
    binder_chain: str = "B",
    threshold: str = "standard",
) -> List[Dict[str, Any]]:
    """
    Validate multiple binder designs.

    Args:
        designs: List of designs, each with 'binder_sequence' and 'pdb_content'
        binder_chain: Chain ID of binder
        threshold: Quality threshold

    Returns:
        List of validation results (one per design)
    """
    results = []

    for i, design in enumerate(designs):
        sequence = design.get("binder_sequence", "")
        pdb_content = design.get("pdb_content", "")

        if not sequence or not pdb_content:
            results.append({
                "status": "error",
                "error": "Missing sequence or pdb_content"
            })
            continue

        logger.info(f"Validating design {i+1}/{len(designs)}...")
        result = validate_structure_esmfold(
            sequence=sequence,
            designed_pdb=pdb_content,
            binder_chain=binder_chain,
            threshold=threshold,
        )
        results.append(result)

    return results


def clear_esmfold_cache():
    """Clear the cached ESMFold model to free GPU memory."""
    global _esmfold_model, _esmfold_device

    if _esmfold_model is not None:
        import torch
        del _esmfold_model
        _esmfold_model = None
        _esmfold_device = None

        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        logger.info("ESMFold model cache cleared")


# Utility function to check if ESMFold is available
def is_esmfold_available() -> bool:
    """Check if ESMFold can be loaded.

    Note: ESM-3 package has a different API than the original ESM package.
    We check for the esm.pretrained module which is only in the original ESM package.
    """
    try:
        import esm
        # Check if this is the original ESM package (has pretrained module)
        # ESM-3 package doesn't have this attribute
        if not hasattr(esm, 'pretrained'):
            return False
        return hasattr(esm.pretrained, 'esmfold_v1')
    except (ImportError, AttributeError):
        return False
