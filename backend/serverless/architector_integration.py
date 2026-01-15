"""
Architector Integration for Lanthanide Coordination Geometry

Architector (Los Alamos National Lab) provides validated 3D structure generation
for f-block metal complexes with ~0.5 Å RMSD accuracy against experimental structures.

Reference: Nandy et al., Nature Communications 2023

Note: Architector is installed in a separate conda environment due to dependency
conflicts. We use subprocess to call it via the conda environment's Python.
"""

import json
import logging
import math
import os
import subprocess
import tempfile
import numpy as np
from typing import Any, Dict, List, Optional, Tuple, Union

logger = logging.getLogger(__name__)

# Path to Architector's Python in conda environment
# Set via ARCHITECTOR_PYTHON env var or default to conda path
ARCHITECTOR_PYTHON = os.environ.get(
    "ARCHITECTOR_PYTHON",
    "/opt/conda/envs/architector/bin/python"
)

# Check if Architector is available via subprocess
def _check_architector_available() -> bool:
    """Check if Architector is available in the conda environment."""
    if not os.path.exists(ARCHITECTOR_PYTHON):
        return False
    try:
        result = subprocess.run(
            [ARCHITECTOR_PYTHON, "-c", "from architector import build_complex; print('ok')"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        return result.returncode == 0 and "ok" in result.stdout
    except Exception:
        return False

# Try direct import first (for local development), fallback to subprocess
ARCHITECTOR_AVAILABLE = False
ARCHITECTOR_USE_SUBPROCESS = False
try:
    from architector import build_complex
    from architector.io_obabel import get_obmol_smiles
    ARCHITECTOR_AVAILABLE = True
except ImportError:
    # Check if available via subprocess (conda environment)
    if _check_architector_available():
        ARCHITECTOR_AVAILABLE = True
        ARCHITECTOR_USE_SUBPROCESS = True
        logger.info("Architector available via conda environment subprocess")


# Lanthanide metal parameters for Architector
ARCHITECTOR_METAL_PARAMS = {
    "TB": {"ox_state": 3, "spin": 6, "element": "Tb"},  # f8
    "GD": {"ox_state": 3, "spin": 7, "element": "Gd"},  # f7
    "EU": {"ox_state": 3, "spin": 6, "element": "Eu"},  # f6
    "LA": {"ox_state": 3, "spin": 0, "element": "La"},  # f0
    "CE": {"ox_state": 3, "spin": 1, "element": "Ce"},  # f1
    "YB": {"ox_state": 3, "spin": 1, "element": "Yb"},  # f13
    "PR": {"ox_state": 3, "spin": 2, "element": "Pr"},  # f2
    "ND": {"ox_state": 3, "spin": 3, "element": "Nd"},  # f3
    "SM": {"ox_state": 3, "spin": 5, "element": "Sm"},  # f5
    "DY": {"ox_state": 3, "spin": 5, "element": "Dy"},  # f9
    "HO": {"ox_state": 3, "spin": 4, "element": "Ho"},  # f10
    "ER": {"ox_state": 3, "spin": 3, "element": "Er"},  # f11
    "TM": {"ox_state": 3, "spin": 2, "element": "Tm"},  # f12
    "LU": {"ox_state": 3, "spin": 0, "element": "Lu"},  # f14
}

# Carboxylate SMILES for different coordination modes
LIGAND_SMILES = {
    "bidentate_carboxylate": "[O-]C(=O)C",  # Acetate (simplest bidentate)
    "monodentate_carboxylate": "[O-]C(=O)C",  # Same SMILES, mode determined by geometry
    "water": "O",  # Water molecule
}

# Geometry mapping for Architector
GEOMETRY_MAP = {
    "square_antiprism": "sap",
    "tricapped_trigonal_prism": "ttp",
    "octahedral": "oct",
    "pentagonal_bipyramidal": "pbp",
    "capped_square_antiprism": "csap",
    "bicapped_trigonal_prism": "btp",
}


def is_architector_available() -> bool:
    """Check if Architector is available for use."""
    return ARCHITECTOR_AVAILABLE


def _run_architector_subprocess(
    metal_element: str,
    ligands: List[str],
    coordination_number: int,
    optimize: bool,
) -> Optional[Dict[str, Any]]:
    """Run Architector via subprocess in conda environment."""
    # Python script to run in conda environment
    # Build inputDict for Architector API
    script = '''
import json
import sys
import numpy as np
from architector import build_complex

params = json.loads(sys.argv[1])
try:
    # Build inputDict for Architector
    ligand_defs = []
    for smiles in params["ligands"]:
        # Carboxylate: coordList=[0,2] for O atoms (bi-dentate)
        ligand_defs.append({
            "smiles": smiles,
            "coordList": [0, 2],  # Both carboxylate oxygens
            "ligType": "bi_cis"
        })

    inputDict = {
        "core": {
            "metal": params["metal_element"],
            "coreCN": params["coordination_number"]
        },
        "ligands": ligand_defs,
        "parameters": {
            "assemble_method": "GFN-FF" if params["optimize"] else "UFF",
            "full_method": "GFN2-xTB" if params["optimize"] else "UFF",
            "n_conformers": 1,
            "relax": params["optimize"]
        }
    }

    results = build_complex(inputDict)
    first_key = list(results.keys())[0]
    atoms = results[first_key]["ase_atoms"]
    coords = atoms.get_positions().tolist()
    symbols = atoms.get_chemical_symbols()
    print(json.dumps({"success": True, "coords": coords, "symbols": symbols}))
except Exception as e:
    import traceback
    print(json.dumps({"success": False, "error": str(e), "traceback": traceback.format_exc()}))
'''

    params = {
        "metal_element": metal_element,
        "ligands": ligands,
        "coordination_number": coordination_number,
        "optimize": optimize,
    }

    try:
        result = subprocess.run(
            [ARCHITECTOR_PYTHON, "-c", script, json.dumps(params)],
            capture_output=True,
            text=True,
            timeout=600,  # 10 min timeout for optimization
        )

        if result.returncode != 0:
            logger.error(f"Architector subprocess failed: {result.stderr}")
            return None

        # Parse JSON output (last line should be JSON)
        output_lines = result.stdout.strip().split('\n')
        json_output = output_lines[-1]  # Get last line which should be JSON
        output = json.loads(json_output)
        if not output.get("success"):
            logger.error(f"Architector build failed: {output.get('error')}")
            if output.get("traceback"):
                logger.error(f"Traceback: {output.get('traceback')}")
            return None

        return {
            "coords": np.array(output["coords"]),
            "symbols": output["symbols"],
        }
    except Exception as e:
        logger.error(f"Architector subprocess error: {e}")
        return None


def generate_architector_complex(
    metal: str,
    coordination_number: int = 8,
    num_bidentate: int = 4,
    num_monodentate: int = 0,
    num_waters: int = 0,
    geometry: str = "square_antiprism",
    optimize: bool = True,
) -> Optional[Dict[str, Any]]:
    """
    Generate a lanthanide coordination complex using Architector.

    Args:
        metal: Lanthanide metal code (TB, GD, EU, etc.)
        coordination_number: Target coordination number (6-10)
        num_bidentate: Number of bidentate carboxylates
        num_monodentate: Number of monodentate carboxylates
        num_waters: Number of coordinating water molecules
        geometry: Target coordination geometry
        optimize: Whether to optimize with GFN2-xTB

    Returns:
        Dictionary with optimized coordinates and analysis, or None if unavailable
    """
    if not ARCHITECTOR_AVAILABLE:
        logger.warning("Architector not available, falling back to manual geometry")
        return None

    metal_upper = metal.upper()
    if metal_upper not in ARCHITECTOR_METAL_PARAMS:
        logger.error(f"Unknown metal: {metal}")
        return None

    metal_params = ARCHITECTOR_METAL_PARAMS[metal_upper]

    # Build ligand list
    ligands = []
    ligand_coords_hints = []

    # Add bidentate carboxylates
    for _ in range(num_bidentate):
        ligands.append(LIGAND_SMILES["bidentate_carboxylate"])
        ligand_coords_hints.append("bidentate")

    # Add monodentate carboxylates
    for _ in range(num_monodentate):
        ligands.append(LIGAND_SMILES["monodentate_carboxylate"])
        ligand_coords_hints.append("monodentate")

    # Add water molecules
    for _ in range(num_waters):
        ligands.append(LIGAND_SMILES["water"])
        ligand_coords_hints.append(None)

    # Map geometry to Architector format
    arch_geometry = GEOMETRY_MAP.get(geometry, "sap")

    try:
        # Use subprocess if Architector is in conda environment
        if ARCHITECTOR_USE_SUBPROCESS:
            subprocess_result = _run_architector_subprocess(
                metal_element=metal_params["element"],
                ligands=ligands,
                coordination_number=coordination_number,
                optimize=optimize,
            )
            if subprocess_result is None:
                return None
            coords = subprocess_result["coords"]
            atom_symbols = subprocess_result["symbols"]
        else:
            # Direct import available - build inputDict
            ligand_defs = []
            for smiles in ligands:
                ligand_defs.append({
                    "smiles": smiles,
                    "coordList": [0, 2],
                    "ligType": "bi_cis"
                })
            inputDict = {
                "core": {
                    "metal": metal_params["element"],
                    "coreCN": coordination_number
                },
                "ligands": ligand_defs,
                "parameters": {
                    "assemble_method": "GFN-FF" if optimize else "UFF",
                    "full_method": "GFN2-xTB" if optimize else "UFF",
                    "n_conformers": 1,
                    "relax": optimize
                }
            }
            results = build_complex(inputDict)
            first_key = list(results.keys())[0]
            atoms = results[first_key]["ase_atoms"]
            coords = atoms.get_positions()
            atom_symbols = atoms.get_chemical_symbols()

        # Find metal index and coordinating oxygen indices
        metal_idx = None
        oxygen_indices = []

        for i, symbol in enumerate(atom_symbols):
            if symbol == metal_params["element"]:
                metal_idx = i
            elif symbol == "O":
                oxygen_indices.append(i)

        if metal_idx is None:
            logger.error("Metal atom not found in Architector output")
            return None

        metal_pos = coords[metal_idx]

        # Calculate coordination analysis
        oxygen_positions = [coords[i] for i in oxygen_indices]
        o_distances = [np.linalg.norm(o_pos - metal_pos) for o_pos in oxygen_positions]

        # Calculate O-O distances
        o_o_distances = []
        for i in range(len(oxygen_positions)):
            for j in range(i + 1, len(oxygen_positions)):
                dist = np.linalg.norm(oxygen_positions[i] - oxygen_positions[j])
                o_o_distances.append(dist)

        return {
            "success": True,
            "metal": metal_upper,
            "metal_position": metal_pos.tolist(),
            "oxygen_positions": [pos.tolist() for pos in oxygen_positions],
            "coordination_number": len(oxygen_indices),
            "metal_oxygen_distances": o_distances,
            "oxygen_oxygen_distances": o_o_distances,
            "min_o_o_distance": min(o_o_distances) if o_o_distances else None,
            "avg_metal_o_distance": np.mean(o_distances) if o_distances else None,
            "all_coordinates": coords.tolist(),
            "atom_symbols": atom_symbols,
            "geometry": geometry,
        }

    except Exception as e:
        logger.error(f"Architector build failed: {e}")
        return None


def convert_architector_to_pdb(
    architector_result: Dict[str, Any],
    metal: str,
    residue_assignments: List[Dict[str, Any]],
    ca_distance: float = 6.0,
) -> str:
    """
    Convert Architector output to PDB format compatible with RFdiffusion.

    This function takes the optimized coordination geometry from Architector
    and constructs a full PDB with proper amino acid residues (ASP/GLU)
    positioned to coordinate the metal ion.

    Args:
        architector_result: Output from generate_architector_complex
        metal: Lanthanide metal code
        residue_assignments: List of dicts with chain, resnum, type for each residue
        ca_distance: Distance from metal to CA atoms

    Returns:
        PDB content string
    """
    if not architector_result or not architector_result.get("success"):
        return ""

    metal_pos = np.array(architector_result["metal_position"])
    oxygen_positions = [np.array(pos) for pos in architector_result["oxygen_positions"]]

    lines = []
    atom_serial = 1

    # We need to pair oxygen positions with residues
    # For bidentate residues, we need pairs of oxygens
    # This is a simplified mapping - real implementation would be more sophisticated

    current_oxygen_idx = 0

    for residue in residue_assignments:
        chain = residue["chain"]
        resnum = residue["resnum"]
        res_type = residue["type"].upper()
        mode = residue.get("mode", "bidentate")

        if current_oxygen_idx >= len(oxygen_positions):
            break

        if mode == "bidentate":
            # Need two oxygens for bidentate
            if current_oxygen_idx + 1 >= len(oxygen_positions):
                break

            o1_pos = oxygen_positions[current_oxygen_idx]
            o2_pos = oxygen_positions[current_oxygen_idx + 1]
            current_oxygen_idx += 2

            # Calculate centroid of the two oxygens
            centroid = (o1_pos + o2_pos) / 2

        else:
            # Monodentate - single oxygen
            o1_pos = oxygen_positions[current_oxygen_idx]
            o2_pos = None
            centroid = o1_pos
            current_oxygen_idx += 1

        # Calculate CA position - away from metal, opposite to carboxylate
        direction = centroid - metal_pos
        direction = direction / np.linalg.norm(direction)
        ca_pos = metal_pos + direction * ca_distance

        # Generate residue atoms
        pdb_lines, atom_serial = _generate_residue_from_architector(
            chain=chain,
            res_num=resnum,
            res_type=res_type,
            mode=mode,
            o1_pos=o1_pos,
            o2_pos=o2_pos,
            ca_pos=ca_pos,
            metal_pos=metal_pos,
            atom_serial_start=atom_serial,
        )
        lines.extend(pdb_lines)

    # Add chain terminator
    if lines:
        lines.append("TER")

    # Add metal atom
    lines.append(
        f"HETATM{atom_serial:>5} {metal.upper():>2}   {metal.upper():>3} L   1    "
        f"{metal_pos[0]:>8.3f}{metal_pos[1]:>8.3f}{metal_pos[2]:>8.3f}"
        f"  1.00  0.00          {metal.upper():>2}"
    )

    lines.append("END")

    return "\n".join(lines)


def _generate_residue_from_architector(
    chain: str,
    res_num: int,
    res_type: str,
    mode: str,
    o1_pos: np.ndarray,
    o2_pos: Optional[np.ndarray],
    ca_pos: np.ndarray,
    metal_pos: np.ndarray,
    atom_serial_start: int,
) -> Tuple[List[str], int]:
    """
    Generate PDB ATOM lines for a single residue from Architector coordinates.

    Args:
        chain: Chain identifier
        res_num: Residue number
        res_type: ASP or GLU
        mode: bidentate or monodentate
        o1_pos: First oxygen position (OD1/OE1)
        o2_pos: Second oxygen position (OD2/OE2), None for monodentate
        ca_pos: CA position
        metal_pos: Metal position for reference
        atom_serial_start: Starting atom serial number

    Returns:
        Tuple of (PDB lines list, next atom serial number)
    """
    lines = []
    serial = atom_serial_start

    # Determine atom naming based on residue type
    if res_type == "GLU":
        o1_name, o2_name = "OE1", "OE2"
        cg_name, cd_name = "CG", "CD"
    else:  # ASP
        o1_name, o2_name = "OD1", "OD2"
        cg_name, cd_name = "CB", "CG"

    # For bidentate, use both oxygens
    # For monodentate, we need to position the second oxygen pointing away
    if o2_pos is None and mode == "monodentate":
        # Position second oxygen ~120° away from metal
        direction_to_metal = metal_pos - o1_pos
        direction_to_metal = direction_to_metal / np.linalg.norm(direction_to_metal)

        # Create perpendicular vector
        if abs(direction_to_metal[0]) < 0.9:
            perp = np.cross(direction_to_metal, np.array([1, 0, 0]))
        else:
            perp = np.cross(direction_to_metal, np.array([0, 1, 0]))
        perp = perp / np.linalg.norm(perp)

        # Position O2 ~2.2 Å from O1, pointing away from metal
        o2_pos = o1_pos - direction_to_metal * 1.1 + perp * 1.9

    # Calculate CG/CD position (carboxyl carbon) - between oxygens
    cg_pos = (o1_pos + o2_pos) / 2

    # Calculate CB position - between CA and CG
    cb_pos = ca_pos + (cg_pos - ca_pos) * 0.4

    # For GLU, add extra CG between CB and CD
    if res_type == "GLU":
        cd_pos = cg_pos
        cg_pos = ca_pos + (cd_pos - ca_pos) * 0.6

    # Calculate backbone N, C, O
    # Simple placement along CA direction
    ca_to_metal = metal_pos - ca_pos
    ca_to_metal = ca_to_metal / np.linalg.norm(ca_to_metal)

    # N is behind CA
    n_pos = ca_pos - ca_to_metal * 1.46

    # C is to the side
    if abs(ca_to_metal[0]) < 0.9:
        side = np.cross(ca_to_metal, np.array([1, 0, 0]))
    else:
        side = np.cross(ca_to_metal, np.array([0, 1, 0]))
    side = side / np.linalg.norm(side)
    c_pos = ca_pos + side * 1.52

    # Backbone O is near C
    o_pos = c_pos + side * 1.23

    def pdb_line(atom_name: str, pos: np.ndarray, element: str) -> str:
        return (
            f"ATOM  {serial:>5} {atom_name:<4} {res_type:>3} {chain}{res_num:>4}    "
            f"{pos[0]:>8.3f}{pos[1]:>8.3f}{pos[2]:>8.3f}"
            f"  1.00  0.00           {element:>2}"
        )

    # Generate ATOM lines in standard order
    lines.append(pdb_line("N", n_pos, "N"))
    serial += 1

    lines.append(pdb_line("CA", ca_pos, "C"))
    serial += 1

    lines.append(pdb_line("C", c_pos, "C"))
    serial += 1

    lines.append(pdb_line("O", o_pos, "O"))
    serial += 1

    lines.append(pdb_line("CB", cb_pos, "C"))
    serial += 1

    if res_type == "GLU":
        lines.append(pdb_line("CG", cg_pos, "C"))
        serial += 1
        lines.append(pdb_line("CD", cd_pos, "C"))
        serial += 1
    else:
        lines.append(pdb_line("CG", cg_pos, "C"))
        serial += 1

    lines.append(pdb_line(o1_name, o1_pos, "O"))
    serial += 1

    lines.append(pdb_line(o2_name, o2_pos, "O"))
    serial += 1

    return lines, serial


def generate_template_with_architector(
    template_def: Dict[str, Any],
    metal: str = "TB",
) -> Optional[str]:
    """
    Generate a template PDB using Architector for geometry optimization.

    This is the main entry point for integrating Architector with the
    template library system.

    Args:
        template_def: Template definition from TEMPLATE_LIBRARY
        metal: Lanthanide metal code

    Returns:
        PDB content string, or None if Architector is unavailable
    """
    if not ARCHITECTOR_AVAILABLE:
        return None

    # Count residue types
    residues = template_def.get("residues", [])
    num_bidentate = sum(1 for r in residues if r.get("mode") == "bidentate")
    num_monodentate = sum(1 for r in residues if r.get("mode") == "monodentate")

    # Count waters
    waters = template_def.get("waters", [])
    num_waters = len(waters)

    # Get coordination number and geometry
    cn = template_def.get("coordination_number", 8)
    geometry = template_def.get("geometry", "square_antiprism")

    # Generate complex (using UFF for fast geometry, skip GFN2-xTB optimization)
    # Note: GFN2-xTB optimization (optimize=True) can take 2+ hours - disabled for now
    result = generate_architector_complex(
        metal=metal,
        coordination_number=cn,
        num_bidentate=num_bidentate,
        num_monodentate=num_monodentate,
        num_waters=num_waters,
        geometry=geometry,
        optimize=False,  # Use UFF (fast) instead of GFN2-xTB (slow)
    )

    if not result or not result.get("success"):
        return None

    # Convert to PDB format
    pdb_content = convert_architector_to_pdb(
        architector_result=result,
        metal=metal,
        residue_assignments=residues,
    )

    return pdb_content


# Fallback: Generate approximate geometry using analytical formulas
# This is used when Architector is not available

def generate_sap_positions(
    metal_pos: np.ndarray,
    distance: float,
    cn: int = 8,
) -> List[np.ndarray]:
    """
    Generate ideal Square Antiprism (SAP) coordination positions.

    SAP geometry has two square faces rotated 45° relative to each other.
    This provides a reasonable approximation for CN=8 lanthanide complexes.

    Args:
        metal_pos: Position of the metal center
        distance: Metal-ligand bond distance
        cn: Coordination number (default 8)

    Returns:
        List of coordination positions
    """
    # SAP geometry parameters
    # Two squares at ±z, rotated 45° relative to each other
    theta = math.acos(1/3)  # ~54.7° - angle from z-axis for ideal SAP

    positions = []

    # Top square (z > 0): angles 0°, 90°, 180°, 270°
    for angle in [0, 90, 180, 270]:
        rad = math.radians(angle)
        pos = metal_pos + distance * np.array([
            math.sin(theta) * math.cos(rad),
            math.sin(theta) * math.sin(rad),
            math.cos(theta),
        ])
        positions.append(pos)

    # Bottom square (z < 0): angles 45°, 135°, 225°, 315° (rotated 45°)
    for angle in [45, 135, 225, 315]:
        rad = math.radians(angle)
        pos = metal_pos + distance * np.array([
            math.sin(theta) * math.cos(rad),
            math.sin(theta) * math.sin(rad),
            -math.cos(theta),
        ])
        positions.append(pos)

    return positions[:cn]


def generate_ttp_positions(
    metal_pos: np.ndarray,
    distance: float,
    cn: int = 9,
) -> List[np.ndarray]:
    """
    Generate ideal Tricapped Trigonal Prism (TTP) coordination positions.

    TTP geometry has two triangular faces and three capping positions.
    Common for CN=9 lanthanide complexes.

    Args:
        metal_pos: Position of the metal center
        distance: Metal-ligand bond distance
        cn: Coordination number (default 9)

    Returns:
        List of coordination positions
    """
    positions = []

    # TTP parameters
    theta_prism = math.radians(48.2)  # Angle from z-axis for prism vertices
    theta_cap = math.radians(90)  # Caps are in equatorial plane

    # Top triangle (z > 0): angles 0°, 120°, 240°
    for angle in [0, 120, 240]:
        rad = math.radians(angle)
        pos = metal_pos + distance * np.array([
            math.sin(theta_prism) * math.cos(rad),
            math.sin(theta_prism) * math.sin(rad),
            math.cos(theta_prism),
        ])
        positions.append(pos)

    # Bottom triangle (z < 0): angles 60°, 180°, 300° (rotated 60°)
    for angle in [60, 180, 300]:
        rad = math.radians(angle)
        pos = metal_pos + distance * np.array([
            math.sin(theta_prism) * math.cos(rad),
            math.sin(theta_prism) * math.sin(rad),
            -math.cos(theta_prism),
        ])
        positions.append(pos)

    # Equatorial caps: angles 30°, 150°, 270°
    for angle in [30, 150, 270]:
        rad = math.radians(angle)
        pos = metal_pos + distance * np.array([
            math.cos(rad),
            math.sin(rad),
            0,
        ])
        positions.append(pos)

    return positions[:cn]
