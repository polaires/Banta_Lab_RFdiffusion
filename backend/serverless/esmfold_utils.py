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

ESMFold is available via HuggingFace Transformers:
- transformers.EsmForProteinFolding (requires ~8GB GPU RAM)
- Model: facebook/esmfold_v1
"""

import os
import logging
import tempfile
from typing import Dict, Any, List, Optional, Tuple
import numpy as np

logger = logging.getLogger(__name__)

# Global model cache
_esmfold_model = None
_esmfold_tokenizer = None
_esmfold_device = None

# Check if ESMFold is available via transformers
def _check_esmfold_available() -> bool:
    """Check if ESMFold can be loaded via HuggingFace Transformers."""
    try:
        from transformers import EsmForProteinFolding, EsmTokenizer
        return True
    except ImportError:
        return False


ESMFOLD_AVAILABLE = _check_esmfold_available()
if ESMFOLD_AVAILABLE:
    logger.info("ESMFold available via HuggingFace Transformers")


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
    # Exploratory threshold for initial testing / debugging
    # High RMSD tolerance for designs that may not fold exactly as designed
    # but still have good predicted confidence
    "exploratory": {
        "plddt_min": 0.70,   # Still require decent confidence
        "rmsd_max": 25.0,    # Very lenient RMSD for early-stage designs
    },
}


def _get_esmfold_model():
    """
    Get or initialize the ESMFold model (cached singleton).

    Returns:
        Tuple of (model, tokenizer, device)
    """
    global _esmfold_model, _esmfold_tokenizer, _esmfold_device

    if _esmfold_model is not None:
        return _esmfold_model, _esmfold_tokenizer, _esmfold_device

    try:
        import torch
        from transformers import EsmForProteinFolding, EsmTokenizer

        # Determine device
        if torch.cuda.is_available():
            _esmfold_device = torch.device("cuda")
            # Check GPU memory - ESMFold needs ~8GB
            gpu_mem = torch.cuda.get_device_properties(0).total_memory / 1e9
            if gpu_mem < 7.0:
                logger.warning(f"GPU has only {gpu_mem:.1f}GB - ESMFold may OOM")
        else:
            _esmfold_device = torch.device("cpu")
            logger.warning("CUDA not available, using CPU for ESMFold (will be very slow)")

        # Load model and tokenizer from HuggingFace
        logger.info(f"Loading ESMFold model on {_esmfold_device}...")
        _esmfold_tokenizer = EsmTokenizer.from_pretrained("facebook/esmfold_v1")
        _esmfold_model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1")
        _esmfold_model = _esmfold_model.eval().to(_esmfold_device)

        # Set chunk size for memory efficiency (if supported)
        if hasattr(_esmfold_model, 'trunk') and hasattr(_esmfold_model.trunk, 'set_chunk_size'):
            _esmfold_model.trunk.set_chunk_size(64)

        logger.info("ESMFold model loaded successfully")
        return _esmfold_model, _esmfold_tokenizer, _esmfold_device

    except ImportError as e:
        logger.error(f"transformers package not installed or ESMFold not available: {e}")
        raise RuntimeError("ESMFold not available. Install with: pip install transformers accelerate")
    except Exception as e:
        logger.error(f"Failed to load ESMFold model: {e}")
        raise


def _convert_outputs_to_pdb(outputs, sequence: str) -> str:
    """Convert ESMFold outputs to PDB format string."""
    import torch

    # Get atom positions from model output
    # positions shape: [num_recycles, batch, seq_len, atom14, 3]
    positions = outputs["positions"][-1][0].cpu()  # Final recycling, first batch: [seq_len, 14, 3]

    # Get pLDDT - shape: [batch, seq_len, 37]
    plddt_raw = outputs["plddt"][0].cpu()  # [seq_len, 37]
    atom_mask = outputs["atom37_atom_exists"][0].cpu()  # [seq_len, 37]
    # Compute per-residue pLDDT as mean over existing atoms
    masked_plddt = plddt_raw * atom_mask
    plddt = masked_plddt.sum(dim=-1) / atom_mask.sum(dim=-1).clamp(min=1)  # [seq_len]

    # Build PDB lines
    pdb_lines = []
    atom_idx = 1

    # Atom14 representation: N, CA, C, O, CB, and sidechains
    # We'll output the first 4 backbone atoms + CB if present
    atom_names = ["N", "CA", "C", "O", "CB"]
    atom_elements = ["N", "C", "C", "O", "C"]

    # 3-letter amino acid codes
    aa_3letter = {
        'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
        'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
        'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
        'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR',
    }

    for res_idx, aa in enumerate(sequence):
        res_num = res_idx + 1
        res_plddt = float(plddt[res_idx]) * 100  # Convert to 0-100 scale
        res_name = aa_3letter.get(aa, 'UNK')

        # Get coordinates for this residue: [14, 3]
        res_coords = positions[res_idx]

        # Output backbone atoms (indices 0-3) and CB (index 4) if not Glycine
        num_atoms = 4 if aa == 'G' else 5
        for atom_local_idx in range(min(num_atoms, len(atom_names))):
            atom_name = atom_names[atom_local_idx]
            element = atom_elements[atom_local_idx]
            x, y, z = res_coords[atom_local_idx].tolist()

            # Skip if coordinates are zero (missing atom)
            if x == 0.0 and y == 0.0 and z == 0.0:
                continue

            pdb_lines.append(
                f"ATOM  {atom_idx:5d}  {atom_name:<3s} {res_name:>3s} A{res_num:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00{res_plddt:6.2f}           {element:>2s}"
            )
            atom_idx += 1

    pdb_lines.append("END")
    return "\n".join(pdb_lines)


def predict_structure_esmfold(
    sequence: str,
    return_pdb: bool = True,
) -> Dict[str, Any]:
    """
    Predict protein structure from sequence using ESMFold via HuggingFace Transformers.

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
    # Validate sequence first
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

    # Model inference using HuggingFace Transformers
    try:
        import torch

        # Get model, tokenizer, and device
        model, tokenizer, device = _get_esmfold_model()

        logger.info(f"Predicting structure for {len(sequence)} residues...")

        # Tokenize sequence
        inputs = tokenizer([sequence], return_tensors="pt", add_special_tokens=False)
        inputs = {k: v.to(device) for k, v in inputs.items()}

        # Predict structure
        with torch.no_grad():
            outputs = model(**inputs)

        # Extract pLDDT (per-residue confidence)
        # plddt shape: [batch, seq_len, 37] (37 atom types)
        # Take mean across atom types or just CA atom (index 1)
        plddt_raw = outputs["plddt"][0].cpu()  # Shape: [seq_len, 37]
        # Use atom14_atom_exists mask to only average over existing atoms
        atom_mask = outputs["atom37_atom_exists"][0].cpu()  # Shape: [seq_len, 37]
        # Compute masked mean per residue
        masked_plddt = plddt_raw * atom_mask
        plddt = (masked_plddt.sum(dim=-1) / atom_mask.sum(dim=-1).clamp(min=1)).numpy()

        # Extract CA coordinates from atom positions
        # positions shape: [num_recycles+1, batch, seq_len, atom14, 3]
        # In atom14 representation: 0=N, 1=CA, 2=C, 3=O
        positions = outputs["positions"][-1][0].cpu().numpy()  # Final recycling step, shape: [seq_len, 14, 3]
        ca_coords = positions[:, 1, :]  # CA atom coordinates (index 1)

        # Convert to PDB format if requested
        pdb_content = None
        if return_pdb:
            pdb_content = _convert_outputs_to_pdb(outputs, sequence)

        result = {
            "status": "completed",
            "plddt": plddt.tolist(),
            "mean_plddt": float(np.mean(plddt)),
            "sequence_length": len(sequence),
            "ca_coords": ca_coords.tolist(),
        }

        if return_pdb:
            result["pdb_content"] = pdb_content

        logger.info(f"Prediction complete: mean pLDDT = {result['mean_plddt']:.2f}")
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

        rmsd_str = f"{rmsd:.2f}" if rmsd is not None else "N/A"
        logger.info(f"Validation: pLDDT={mean_plddt:.2f}, RMSD={rmsd_str}Å, passed={validation_passed}")
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
    global _esmfold_model, _esmfold_tokenizer, _esmfold_device

    if _esmfold_model is not None:
        import torch
        del _esmfold_model
        del _esmfold_tokenizer
        _esmfold_model = None
        _esmfold_tokenizer = None
        _esmfold_device = None

        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        logger.info("ESMFold model cache cleared")


# Utility function to check if ESMFold is available
def is_esmfold_available() -> bool:
    """Check if ESMFold can be loaded via HuggingFace Transformers.

    Returns the pre-computed availability flag.
    """
    return ESMFOLD_AVAILABLE


# =============================================================================
# RF3 (RoseTTAFold3) Validation Support
# =============================================================================

def _extract_ca_coords_from_cif(cif_text: str, chain: str = None) -> np.ndarray:
    """Extract CA atom coordinates from mmCIF text.

    Parses the _atom_site loop in mmCIF format to find CA atoms.

    Args:
        cif_text: mmCIF file content as string
        chain: Optional chain ID to filter (matches auth_asym_id or label_asym_id)

    Returns:
        np.ndarray of shape (N, 3) with CA coordinates
    """
    coords = []
    in_atom_site = False
    col_map = {}
    col_idx = 0

    for line in cif_text.splitlines():
        if line.startswith("_atom_site."):
            in_atom_site = True
            field_name = line.strip().split(".")[1].strip()
            col_map[field_name] = col_idx
            col_idx += 1
            continue

        if in_atom_site and (line.startswith("_") or line.startswith("#")):
            in_atom_site = False
            col_map = {}
            col_idx = 0
            continue

        if in_atom_site and line.strip() and not line.startswith("_") and not line.startswith("#"):
            if not col_map:
                continue
            parts = line.split()
            if len(parts) < max(col_map.values()) + 1:
                continue

            # Get atom name
            label_atom_id = col_map.get("label_atom_id") or col_map.get("auth_atom_id")
            if label_atom_id is None:
                continue
            atom_name = parts[label_atom_id].strip('"')
            if atom_name != "CA":
                continue

            # Filter to ATOM records only (not HETATM)
            group_idx = col_map.get("group_PDB")
            if group_idx is not None and parts[group_idx] != "ATOM":
                continue

            # Filter by chain if specified
            if chain is not None:
                chain_idx = col_map.get("auth_asym_id") or col_map.get("label_asym_id")
                if chain_idx is not None and parts[chain_idx] != chain:
                    continue

            # Extract coordinates
            x_idx = col_map.get("Cartn_x")
            y_idx = col_map.get("Cartn_y")
            z_idx = col_map.get("Cartn_z")
            if x_idx is None or y_idx is None or z_idx is None:
                continue
            try:
                coords.append([float(parts[x_idx]), float(parts[y_idx]), float(parts[z_idx])])
            except ValueError:
                continue

    return np.array(coords) if coords else np.array([]).reshape(0, 3)


def calculate_backbone_rmsd_from_cif(
    pdb_content: str,
    cif_content: str,
    pdb_chain: str = None,
    cif_chain: str = None,
) -> Dict[str, Any]:
    """Calculate backbone RMSD between a PDB structure and a CIF structure.

    Uses CA atoms for Kabsch alignment and RMSD calculation.

    Args:
        pdb_content: Reference PDB structure (designed backbone)
        cif_content: Predicted CIF structure (e.g., from RF3)
        pdb_chain: Chain ID in PDB structure (optional)
        cif_chain: Chain ID in CIF structure (optional)

    Returns:
        Dict with rmsd, n_atoms, aligned status, or error
    """
    try:
        coords_pdb = _extract_ca_coords_from_pdb(pdb_content, pdb_chain)
        coords_cif = _extract_ca_coords_from_cif(cif_content, cif_chain)

        if len(coords_pdb) == 0 or len(coords_cif) == 0:
            return {
                "status": "error",
                "error": f"Could not extract CA coordinates (PDB: {len(coords_pdb)}, CIF: {len(coords_cif)})"
            }

        min_len = min(len(coords_pdb), len(coords_cif))
        if len(coords_pdb) != len(coords_cif):
            logger.warning(f"Length mismatch: PDB={len(coords_pdb)} vs CIF={len(coords_cif)}, using first {min_len}")
            coords_pdb = coords_pdb[:min_len]
            coords_cif = coords_cif[:min_len]

        aligned_rmsd = _kabsch_rmsd(coords_pdb, coords_cif)

        return {
            "status": "completed",
            "rmsd": float(aligned_rmsd),
            "n_atoms": min_len,
            "aligned": True,
        }

    except Exception as e:
        logger.error(f"CIF RMSD calculation failed: {e}")
        return {
            "status": "error",
            "error": str(e),
        }


def is_rf3_available() -> bool:
    """Check if RF3 (RoseTTAFold3) inference is available in-container.

    Returns:
        True if inference_utils can be imported and RF3 engine exists
    """
    try:
        import inference_utils
        return hasattr(inference_utils, 'run_rf3_inference')
    except ImportError:
        return False


def validate_structure_rf3(
    sequence: str,
    designed_pdb: str,
    binder_chain: str = None,
    ligand_smiles: str = None,
    threshold: str = "standard",
) -> Dict[str, Any]:
    """Validate a designed protein by predicting its structure with RF3 and comparing to design.

    Mirrors validate_structure_esmfold() signature and return format, but uses
    RoseTTAFold3 for structure prediction. RF3 outputs mmCIF format.

    Args:
        sequence: Designed protein sequence
        designed_pdb: Reference PDB structure (designed backbone)
        binder_chain: Chain ID of the binder in designed_pdb (optional)
        ligand_smiles: SMILES string for ligand-aware prediction (optional)
        threshold: Quality threshold ("strict", "standard", "relaxed", "exploratory")

    Returns:
        Dict with:
        - status: "completed" or "error"
        - rf3_plddt: Average pLDDT from RF3 prediction
        - rf3_ptm: Predicted TM score
        - rf3_iptm: Interface pTM (if ligand provided)
        - rf3_rmsd: Backbone RMSD between design and prediction (Å)
        - rf3_cif_content: RF3 predicted structure in CIF format
        - validation_passed: True if metrics pass thresholds
    """
    try:
        import inference_utils

        logger.info(f"RF3 prediction for {len(sequence)} residues (ligand={'yes' if ligand_smiles else 'no'})...")

        # Run RF3 inference
        rf3_result = inference_utils.run_rf3_inference(
            sequence=sequence,
            name="rf3_validation",
            ligand_smiles=ligand_smiles,
        )

        if rf3_result.get("status") != "completed":
            error = rf3_result.get("error", "Unknown error")
            logger.warning(f"RF3 prediction failed: {error}")
            return {
                "status": "error",
                "error": f"RF3 prediction failed: {error}",
            }

        # Extract results
        result_data = rf3_result.get("result", {})
        confidences = result_data.get("confidences") or {}
        predictions = result_data.get("predictions", [])

        rf3_plddt = confidences.get("mean_plddt")
        rf3_ptm = confidences.get("ptm")
        rf3_iptm = confidences.get("iptm")

        # Get CIF content from first prediction
        cif_content = None
        for pred in predictions:
            cif = pred.get("cif_content")
            if cif:
                cif_content = cif
                break

        # Compute RMSD vs designed backbone
        rf3_rmsd = None
        if cif_content and designed_pdb:
            rmsd_result = calculate_backbone_rmsd_from_cif(
                pdb_content=designed_pdb,
                cif_content=cif_content,
                pdb_chain=binder_chain,
            )
            if rmsd_result.get("status") == "completed":
                rf3_rmsd = rmsd_result["rmsd"]
            else:
                logger.warning(f"RF3 RMSD calculation failed: {rmsd_result.get('error')}")

        # Check against thresholds
        thresholds = ESMFOLD_THRESHOLDS.get(threshold, ESMFOLD_THRESHOLDS["standard"])
        passes_plddt = rf3_plddt is not None and rf3_plddt >= thresholds["plddt_min"]
        passes_rmsd = rf3_rmsd is None or rf3_rmsd <= thresholds["rmsd_max"]
        validation_passed = passes_plddt and passes_rmsd

        result = {
            "status": "completed",
            "rf3_plddt": rf3_plddt,
            "rf3_ptm": rf3_ptm,
            "rf3_iptm": rf3_iptm,
            "rf3_rmsd": rf3_rmsd,
            "rf3_cif_content": cif_content,
            "validation_passed": validation_passed,
            "threshold_used": threshold,
        }

        rmsd_str = f"{rf3_rmsd:.2f}" if rf3_rmsd is not None else "N/A"
        plddt_str = f"{rf3_plddt:.2f}" if rf3_plddt is not None else "N/A"
        logger.info(f"RF3 validation: pLDDT={plddt_str}, pTM={rf3_ptm}, RMSD={rmsd_str}Å, passed={validation_passed}")
        return result

    except ImportError:
        logger.warning("RF3 not available (inference_utils not importable)")
        return {
            "status": "error",
            "error": "RF3 not available",
        }
    except Exception as e:
        logger.error(f"RF3 validation failed: {e}")
        return {
            "status": "error",
            "error": str(e),
        }
