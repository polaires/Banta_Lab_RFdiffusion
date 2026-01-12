"""
SASA (Solvent Accessible Surface Area) Calculator

Calculate solvent accessible surface area for protein structures.
Uses FreeSASA library if available, with biotite fallback.
"""

import tempfile
import os
from typing import Dict, List, Any, Optional, Tuple
import numpy as np


# Try to import SASA calculation libraries
FREESASA_AVAILABLE = False
BIOTITE_SASA_AVAILABLE = False

try:
    import freesasa
    FREESASA_AVAILABLE = True
except ImportError:
    pass

try:
    from biotite.structure import sasa as biotite_sasa
    from biotite.structure.io.pdb import PDBFile
    import io
    BIOTITE_SASA_AVAILABLE = True
except ImportError:
    pass


# Van der Waals radii for common atoms (in Angstroms)
VDW_RADII = {
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "H": 1.20,
    "P": 1.80,
    "FE": 1.72,
    "ZN": 1.39,
    "CA": 1.94,
    "MG": 1.60,
    "MN": 1.61,
    "CU": 1.40,
    "CO": 1.35,
    "NI": 1.35,
}

# Burial thresholds (relative SASA)
BURIAL_THRESHOLDS = {
    "buried": 0.05,      # <5% exposed
    "partial": 0.25,     # 5-25% exposed
    "exposed": 1.0,      # >25% exposed
}

# Reference SASA values for fully exposed residues (Gly-X-Gly peptide)
REFERENCE_SASA = {
    "ALA": 115.0,
    "ARG": 249.0,
    "ASN": 158.0,
    "ASP": 150.0,
    "CYS": 141.0,
    "GLN": 189.0,
    "GLU": 183.0,
    "GLY": 85.0,
    "HIS": 194.0,
    "ILE": 182.0,
    "LEU": 180.0,
    "LYS": 211.0,
    "MET": 204.0,
    "PHE": 218.0,
    "PRO": 143.0,
    "SER": 122.0,
    "THR": 146.0,
    "TRP": 259.0,
    "TYR": 229.0,
    "VAL": 160.0,
}


def _calculate_sasa_freesasa(pdb_content: str) -> Dict[str, Any]:
    """
    Calculate SASA using FreeSASA library.

    Args:
        pdb_content: PDB file content

    Returns:
        Dict with SASA results
    """
    # Write to temp file
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode='w') as f:
        f.write(pdb_content)
        temp_path = f.name

    try:
        # Calculate SASA
        structure = freesasa.Structure(temp_path)
        result = freesasa.calc(structure)

        # Get per-residue SASA
        residue_sasa = {}
        for i in range(structure.nAtoms()):
            chain = structure.chainLabel(i)
            res_name = structure.residueName(i)
            res_num = structure.residueNumber(i)
            key = f"{chain}:{res_name}{res_num}"

            if key not in residue_sasa:
                residue_sasa[key] = {
                    "chain": chain,
                    "residue_name": res_name,
                    "residue_number": res_num,
                    "total_sasa": 0.0,
                    "atoms": {},
                }

            atom_name = structure.atomName(i)
            atom_sasa = result.atomArea(i)
            residue_sasa[key]["total_sasa"] += atom_sasa
            residue_sasa[key]["atoms"][atom_name] = atom_sasa

        # Calculate relative SASA and burial classification
        for key, data in residue_sasa.items():
            ref_sasa = REFERENCE_SASA.get(data["residue_name"], 200.0)
            data["relative_sasa"] = data["total_sasa"] / ref_sasa if ref_sasa > 0 else 0
            data["burial_classification"] = classify_burial(data["relative_sasa"])

        return {
            "success": True,
            "method": "freesasa",
            "total_sasa": result.totalArea(),
            "polar_sasa": result.polarArea() if hasattr(result, 'polarArea') else None,
            "apolar_sasa": result.apolarArea() if hasattr(result, 'apolarArea') else None,
            "per_residue": residue_sasa,
        }

    finally:
        os.unlink(temp_path)


def _calculate_sasa_biotite(pdb_content: str) -> Dict[str, Any]:
    """
    Calculate SASA using Biotite library.

    Args:
        pdb_content: PDB file content

    Returns:
        Dict with SASA results
    """
    # Parse PDB
    pdb_file = PDBFile.read(io.StringIO(pdb_content))
    structure = pdb_file.get_structure(model=1)

    # Calculate SASA for each atom
    atom_sasa = biotite_sasa(structure, vdw_radii="ProtOr")

    # Aggregate per residue
    residue_sasa = {}
    for i, atom in enumerate(structure):
        chain = atom.chain_id
        res_name = atom.res_name
        res_num = atom.res_id
        key = f"{chain}:{res_name}{res_num}"

        if key not in residue_sasa:
            residue_sasa[key] = {
                "chain": chain,
                "residue_name": res_name,
                "residue_number": res_num,
                "total_sasa": 0.0,
                "atoms": {},
            }

        atom_name = atom.atom_name
        sasa_value = atom_sasa[i]
        residue_sasa[key]["total_sasa"] += sasa_value
        residue_sasa[key]["atoms"][atom_name] = float(sasa_value)

    # Calculate relative SASA and burial classification
    total_sasa = sum(data["total_sasa"] for data in residue_sasa.values())

    for key, data in residue_sasa.items():
        ref_sasa = REFERENCE_SASA.get(data["residue_name"], 200.0)
        data["relative_sasa"] = data["total_sasa"] / ref_sasa if ref_sasa > 0 else 0
        data["burial_classification"] = classify_burial(data["relative_sasa"])

    return {
        "success": True,
        "method": "biotite",
        "total_sasa": float(total_sasa),
        "polar_sasa": None,  # Biotite doesn't separate by polarity
        "apolar_sasa": None,
        "per_residue": residue_sasa,
    }


def _calculate_sasa_simple(pdb_content: str) -> Dict[str, Any]:
    """
    Simple SASA estimation based on exposed atom counting.
    Fallback when no SASA library is available.

    Args:
        pdb_content: PDB file content

    Returns:
        Dict with estimated SASA results
    """
    # Parse atoms
    atoms = []
    for line in pdb_content.split('\n'):
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            continue

        try:
            chain = line[21]
            res_name = line[17:20].strip()
            res_num = int(line[22:26].strip())
            atom_name = line[12:16].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            element = atom_name[0].upper()
            radius = VDW_RADII.get(element, 1.7)

            atoms.append({
                "chain": chain,
                "res_name": res_name,
                "res_num": res_num,
                "atom_name": atom_name,
                "position": np.array([x, y, z]),
                "radius": radius,
            })
        except (ValueError, IndexError):
            continue

    if not atoms:
        return {
            "success": False,
            "error": "No atoms found in PDB",
        }

    # Simple neighbor counting approach
    # Atoms with fewer neighbors are more exposed
    positions = np.array([a["position"] for a in atoms])

    residue_sasa = {}
    for i, atom in enumerate(atoms):
        key = f"{atom['chain']}:{atom['res_name']}{atom['res_num']}"

        if key not in residue_sasa:
            residue_sasa[key] = {
                "chain": atom["chain"],
                "residue_name": atom["res_name"],
                "residue_number": atom["res_num"],
                "total_sasa": 0.0,
                "atoms": {},
            }

        # Count neighbors within 8 Angstroms
        distances = np.linalg.norm(positions - atom["position"], axis=1)
        n_neighbors = np.sum((distances > 0) & (distances < 8.0))

        # Estimate exposure (more neighbors = less exposed)
        max_neighbors = 30
        exposure = max(0, 1.0 - n_neighbors / max_neighbors)

        # Approximate SASA for this atom
        atom_sasa = 4 * np.pi * (atom["radius"] + 1.4) ** 2 * exposure

        residue_sasa[key]["total_sasa"] += atom_sasa
        residue_sasa[key]["atoms"][atom["atom_name"]] = atom_sasa

    # Calculate relative SASA and burial classification
    total_sasa = sum(data["total_sasa"] for data in residue_sasa.values())

    for key, data in residue_sasa.items():
        ref_sasa = REFERENCE_SASA.get(data["residue_name"], 200.0)
        data["relative_sasa"] = data["total_sasa"] / ref_sasa if ref_sasa > 0 else 0
        data["burial_classification"] = classify_burial(data["relative_sasa"])

    return {
        "success": True,
        "method": "simple_estimation",
        "total_sasa": total_sasa,
        "polar_sasa": None,
        "apolar_sasa": None,
        "per_residue": residue_sasa,
        "warning": "Using simple estimation - install freesasa or biotite for accurate results",
    }


def calculate_sasa(pdb_content: str) -> Dict[str, Any]:
    """
    Calculate Solvent Accessible Surface Area.

    Uses FreeSASA if available, falls back to Biotite, then simple estimation.

    Args:
        pdb_content: PDB file content

    Returns:
        Dict with SASA results including:
        - total_sasa: Total protein SASA
        - per_residue: Per-residue SASA values
        - burial_classification: Classification for each residue
    """
    if FREESASA_AVAILABLE:
        try:
            return _calculate_sasa_freesasa(pdb_content)
        except Exception as e:
            print(f"FreeSASA failed: {e}, trying biotite...")

    if BIOTITE_SASA_AVAILABLE:
        try:
            return _calculate_sasa_biotite(pdb_content)
        except Exception as e:
            print(f"Biotite SASA failed: {e}, using simple estimation...")

    return _calculate_sasa_simple(pdb_content)


def classify_burial(relative_sasa: float) -> str:
    """
    Classify residue burial based on relative SASA.

    Args:
        relative_sasa: Relative SASA (0-1, can be >1 for flexible residues)

    Returns:
        Classification: "buried", "partial", or "exposed"
    """
    if relative_sasa < BURIAL_THRESHOLDS["buried"]:
        return "buried"
    elif relative_sasa < BURIAL_THRESHOLDS["partial"]:
        return "partial"
    else:
        return "exposed"


def get_binding_pocket_sasa(
    pdb_content: str,
    pocket_residues: List[str]
) -> Dict[str, Any]:
    """
    Calculate SASA for a specific binding pocket.

    Args:
        pdb_content: PDB file content
        pocket_residues: List of residue identifiers (e.g., ["A:ASP15", "A:GLU20"])

    Returns:
        Dict with pocket SASA analysis
    """
    full_sasa = calculate_sasa(pdb_content)

    if not full_sasa.get("success"):
        return full_sasa

    per_residue = full_sasa["per_residue"]

    pocket_analysis = {
        "residues": [],
        "total_pocket_sasa": 0.0,
        "avg_relative_sasa": 0.0,
        "burial_summary": {"buried": 0, "partial": 0, "exposed": 0},
    }

    for res_id in pocket_residues:
        # Try different key formats
        for key in per_residue:
            if res_id in key or key.endswith(res_id.split(":")[-1]):
                data = per_residue[key]
                pocket_analysis["residues"].append({
                    "residue_id": key,
                    "sasa": data["total_sasa"],
                    "relative_sasa": data["relative_sasa"],
                    "classification": data["burial_classification"],
                })
                pocket_analysis["total_pocket_sasa"] += data["total_sasa"]
                pocket_analysis["burial_summary"][data["burial_classification"]] += 1
                break

    if pocket_analysis["residues"]:
        pocket_analysis["avg_relative_sasa"] = np.mean(
            [r["relative_sasa"] for r in pocket_analysis["residues"]]
        )

    return {
        "success": True,
        "pocket_analysis": pocket_analysis,
        "full_sasa": full_sasa,
    }


def suggest_rfd3_burial_params(
    sasa_analysis: Dict[str, Any],
    target_burial: str = "buried"
) -> Dict[str, str]:
    """
    Generate RFD3 select_buried/select_exposed parameters based on SASA analysis.

    Args:
        sasa_analysis: Output from calculate_sasa()
        target_burial: Target burial state for binding pocket

    Returns:
        Dict mapping residue IDs to atom selections for RFD3
    """
    if not sasa_analysis.get("success"):
        return {}

    per_residue = sasa_analysis["per_residue"]
    burial_params = {}

    for key, data in per_residue.items():
        # Extract chain and residue number
        parts = key.split(":")
        if len(parts) != 2:
            continue

        chain = parts[0]
        res_info = parts[1]
        res_name = res_info[:3]
        res_num = res_info[3:]

        rfd3_key = f"{chain}{res_num}"

        current_burial = data["burial_classification"]

        # If residue needs to change burial state
        if target_burial == "buried" and current_burial == "exposed":
            burial_params[rfd3_key] = "TIP"  # Bury the tip atom
        elif target_burial == "exposed" and current_burial == "buried":
            burial_params[rfd3_key] = "TIP"  # Expose the tip atom

    return burial_params
