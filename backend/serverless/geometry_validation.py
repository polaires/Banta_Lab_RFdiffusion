"""
Geometry Validation Module for Lanthanide Coordination Templates

This module provides quality checks for generated lanthanide binding site
geometries. It validates:
1. Metal-Oxygen distances (should be 2.2-2.6 Å)
2. Oxygen-Oxygen distances (should be >2.4 Å to avoid clashes)
3. Coordination number (should match expected value)
4. Coordination geometry (RMSD to ideal polyhedra)

Publication-quality templates should pass all validation checks.
"""

import logging
import math
import numpy as np
from itertools import combinations
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


# Validation thresholds based on literature and CSD analysis
VALIDATION_THRESHOLDS = {
    # Metal-Oxygen bond distances (Å)
    "ln_o_min": 2.2,   # Minimum acceptable Ln-O distance
    "ln_o_max": 2.7,   # Maximum acceptable Ln-O distance
    "ln_o_ideal_min": 2.30,  # Ideal range lower bound
    "ln_o_ideal_max": 2.45,  # Ideal range upper bound

    # Oxygen-Oxygen distances (Å)
    "o_o_clash": 2.4,  # Below this is a steric clash
    "o_o_warn": 2.8,   # Below this triggers a warning
    "o_o_max": 5.0,    # Above this oxygens are too far apart

    # Coordination number
    "cn_min": 6,       # Minimum coordination number
    "cn_max": 10,      # Maximum coordination number

    # Geometry RMSD (Å)
    "rmsd_excellent": 0.3,  # Excellent match to ideal geometry
    "rmsd_good": 0.5,       # Good match (Architector target)
    "rmsd_acceptable": 1.0, # Still acceptable
}


def _generate_ideal_sap() -> List[List[float]]:
    """Generate ideal square antiprism positions."""
    theta = math.acos(1/3)  # ~54.7°
    positions = []

    # Top square
    for angle in [0, 90, 180, 270]:
        rad = math.radians(angle)
        positions.append([
            math.sin(theta) * math.cos(rad),
            math.sin(theta) * math.sin(rad),
            math.cos(theta),
        ])

    # Bottom square (rotated 45°)
    for angle in [45, 135, 225, 315]:
        rad = math.radians(angle)
        positions.append([
            math.sin(theta) * math.cos(rad),
            math.sin(theta) * math.sin(rad),
            -math.cos(theta),
        ])

    return positions


def _generate_ideal_ttp() -> List[List[float]]:
    """Generate ideal tricapped trigonal prism positions."""
    theta_prism = math.radians(48.2)
    positions = []

    # Top triangle
    for angle in [0, 120, 240]:
        rad = math.radians(angle)
        positions.append([
            math.sin(theta_prism) * math.cos(rad),
            math.sin(theta_prism) * math.sin(rad),
            math.cos(theta_prism),
        ])

    # Bottom triangle
    for angle in [60, 180, 300]:
        rad = math.radians(angle)
        positions.append([
            math.sin(theta_prism) * math.cos(rad),
            math.sin(theta_prism) * math.sin(rad),
            -math.cos(theta_prism),
        ])

    # Equatorial caps
    for angle in [30, 150, 270]:
        rad = math.radians(angle)
        positions.append([
            math.cos(rad),
            math.sin(rad),
            0,
        ])

    return positions


# Ideal coordination geometries (vertex positions on unit sphere)
# These are used for RMSD comparison
# Note: Defined after helper functions to avoid forward reference errors
IDEAL_GEOMETRIES = {
    "octahedral": {
        "cn": 6,
        "positions": [
            [1, 0, 0], [-1, 0, 0],
            [0, 1, 0], [0, -1, 0],
            [0, 0, 1], [0, 0, -1],
        ],
    },
    "square_antiprism": {
        "cn": 8,
        "positions": _generate_ideal_sap(),
    },
    "tricapped_trigonal_prism": {
        "cn": 9,
        "positions": _generate_ideal_ttp(),
    },
    "distorted_square_antiprism": {
        "cn": 8,
        "positions": _generate_ideal_sap(),  # Same as SAP for comparison
    },
}


def parse_pdb_for_validation(pdb_content: str) -> Dict[str, Any]:
    """
    Parse PDB content to extract metal and coordinating atom positions.

    Args:
        pdb_content: PDB file content as string

    Returns:
        Dictionary with metal_pos, oxygen_positions, atom_details
    """
    metal_pos = None
    metal_type = None
    oxygen_positions = []
    atom_details = []

    # Lanthanide metal codes
    LANTHANIDES = {"TB", "GD", "EU", "LA", "CE", "YB", "PR", "ND", "SM", "DY", "HO", "ER", "TM", "LU"}

    for line in pdb_content.split("\n"):
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue

        try:
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain = line[21] if len(line) > 21 else ""
            res_num = int(line[22:26].strip()) if line[22:26].strip() else 0
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            pos = np.array([x, y, z])

            # Check if this is a metal
            if res_name.upper() in LANTHANIDES:
                metal_pos = pos
                metal_type = res_name.upper()
                continue

            # Check if this is a coordinating oxygen
            if atom_name in ["OD1", "OD2", "OE1", "OE2"]:
                oxygen_positions.append(pos)
                atom_details.append({
                    "name": atom_name,
                    "resname": res_name,
                    "chain": chain,
                    "resnum": res_num,
                    "position": pos,
                })

        except (ValueError, IndexError):
            continue

    return {
        "metal_pos": metal_pos,
        "metal_type": metal_type,
        "oxygen_positions": oxygen_positions,
        "atom_details": atom_details,
    }


def calculate_distances(
    metal_pos: np.ndarray,
    oxygen_positions: List[np.ndarray],
    atom_details: Optional[List[Dict]] = None,
) -> Dict[str, Any]:
    """
    Calculate all relevant distances for validation.

    Distinguishes between:
    - Intra-residue O-O distances (OE1-OE2 or OD1-OD2 within same carboxylate, ~2.21 Å expected)
    - Inter-residue O-O distances (between different residues, should be >2.4 Å)

    Args:
        metal_pos: Metal atom position
        oxygen_positions: List of coordinating oxygen positions
        atom_details: Optional list of atom info dicts with chain, resnum for residue grouping

    Returns:
        Dictionary with distance metrics including separated intra/inter O-O distances
    """
    if metal_pos is None or not oxygen_positions:
        return {
            "ln_o_distances": [],
            "o_o_distances": [],
            "o_o_intra_residue": [],
            "o_o_inter_residue": [],
            "ln_o_mean": None,
            "ln_o_min": None,
            "ln_o_max": None,
            "o_o_min": None,
            "o_o_max": None,
            "o_o_mean": None,
            "o_o_inter_min": None,
            "o_o_inter_max": None,
        }

    # Metal-Oxygen distances
    ln_o_distances = [
        float(np.linalg.norm(o_pos - metal_pos))
        for o_pos in oxygen_positions
    ]

    # Oxygen-Oxygen distances - separate intra vs inter residue
    o_o_distances = []
    o_o_intra_residue = []  # Within same carboxylate (expected ~2.21 Å)
    o_o_inter_residue = []  # Between different residues (should be >2.4 Å)

    for i, pos1 in enumerate(oxygen_positions):
        for j, pos2 in enumerate(oxygen_positions):
            if i < j:
                dist = float(np.linalg.norm(pos1 - pos2))
                o_o_distances.append(dist)

                # Check if same residue (intra-carboxylate)
                if atom_details and i < len(atom_details) and j < len(atom_details):
                    detail1 = atom_details[i]
                    detail2 = atom_details[j]
                    same_residue = (
                        detail1.get("chain") == detail2.get("chain") and
                        detail1.get("resnum") == detail2.get("resnum")
                    )
                    if same_residue:
                        o_o_intra_residue.append({
                            "distance": dist,
                            "residue": f"{detail1.get('resname', '?')}{detail1.get('chain', '?')}{detail1.get('resnum', '?')}",
                            "atoms": f"{detail1.get('name', '?')}-{detail2.get('name', '?')}",
                        })
                    else:
                        o_o_inter_residue.append({
                            "distance": dist,
                            "residue1": f"{detail1.get('resname', '?')}{detail1.get('chain', '?')}{detail1.get('resnum', '?')}.{detail1.get('name', '?')}",
                            "residue2": f"{detail2.get('resname', '?')}{detail2.get('chain', '?')}{detail2.get('resnum', '?')}.{detail2.get('name', '?')}",
                        })
                else:
                    # No details available, treat all as inter-residue for safety
                    o_o_inter_residue.append({"distance": dist})

    # Extract just distances for inter-residue
    inter_distances = [d["distance"] for d in o_o_inter_residue]

    return {
        "ln_o_distances": ln_o_distances,
        "o_o_distances": o_o_distances,
        "o_o_intra_residue": o_o_intra_residue,
        "o_o_inter_residue": o_o_inter_residue,
        "ln_o_mean": float(np.mean(ln_o_distances)) if ln_o_distances else None,
        "ln_o_min": float(min(ln_o_distances)) if ln_o_distances else None,
        "ln_o_max": float(max(ln_o_distances)) if ln_o_distances else None,
        "o_o_min": float(min(o_o_distances)) if o_o_distances else None,
        "o_o_max": float(max(o_o_distances)) if o_o_distances else None,
        "o_o_mean": float(np.mean(o_o_distances)) if o_o_distances else None,
        "o_o_inter_min": float(min(inter_distances)) if inter_distances else None,
        "o_o_inter_max": float(max(inter_distances)) if inter_distances else None,
    }


def calculate_geometry_rmsd(
    oxygen_positions: List[np.ndarray],
    metal_pos: np.ndarray,
    target_geometry: str = "square_antiprism",
) -> Optional[float]:
    """
    Calculate RMSD between actual coordination and ideal geometry.

    Uses Kabsch alignment to find optimal superposition before RMSD.

    Args:
        oxygen_positions: Actual coordinating oxygen positions
        metal_pos: Metal center position
        target_geometry: Name of ideal geometry to compare against

    Returns:
        RMSD in Angstroms, or None if calculation fails
    """
    if target_geometry not in IDEAL_GEOMETRIES:
        logger.warning(f"Unknown geometry: {target_geometry}")
        return None

    ideal = IDEAL_GEOMETRIES[target_geometry]
    ideal_cn = ideal["cn"]
    ideal_positions = np.array(ideal["positions"])

    actual_cn = len(oxygen_positions)
    if actual_cn != ideal_cn:
        # If CN doesn't match, we can't do direct comparison
        logger.info(f"CN mismatch: actual={actual_cn}, ideal={ideal_cn}")
        return None

    # Convert to metal-centered, unit-distance coordinates
    actual_normalized = []
    for o_pos in oxygen_positions:
        vec = o_pos - metal_pos
        dist = np.linalg.norm(vec)
        if dist > 0:
            actual_normalized.append(vec / dist)
        else:
            return None

    actual_array = np.array(actual_normalized)

    # Kabsch alignment
    try:
        # Center both sets
        actual_centered = actual_array - np.mean(actual_array, axis=0)
        ideal_centered = ideal_positions - np.mean(ideal_positions, axis=0)

        # Compute optimal rotation matrix
        H = actual_centered.T @ ideal_centered
        U, S, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T

        # Apply rotation
        actual_rotated = actual_centered @ R

        # Calculate RMSD
        diff = actual_rotated - ideal_centered
        rmsd = float(np.sqrt(np.mean(np.sum(diff**2, axis=1))))

        return rmsd

    except Exception as e:
        logger.error(f"RMSD calculation failed: {e}")
        return None


def validate_coordination_geometry(
    pdb_content: str,
    expected_cn: Optional[int] = None,
    expected_geometry: str = "square_antiprism",
    cutoff_distance: float = 3.5,
) -> Dict[str, Any]:
    """
    Comprehensive validation of lanthanide coordination geometry.

    This is the main validation entry point. It checks all quality metrics
    and returns a detailed report.

    Args:
        pdb_content: PDB file content
        expected_cn: Expected coordination number (optional)
        expected_geometry: Expected coordination geometry
        cutoff_distance: Maximum distance to consider for coordination

    Returns:
        Validation report dictionary with all metrics and pass/fail status
    """
    # Parse PDB
    parsed = parse_pdb_for_validation(pdb_content)

    if parsed["metal_pos"] is None:
        return {
            "valid": False,
            "error": "No lanthanide metal found in PDB",
            "details": {},
        }

    metal_pos = parsed["metal_pos"]
    oxygen_positions = parsed["oxygen_positions"]

    # Filter to coordinating oxygens only
    coord_oxygens = []
    coord_details = []
    for i, o_pos in enumerate(oxygen_positions):
        dist = np.linalg.norm(o_pos - metal_pos)
        if dist <= cutoff_distance:
            coord_oxygens.append(o_pos)
            if i < len(parsed["atom_details"]):
                coord_details.append(parsed["atom_details"][i])

    # Calculate distances with atom details for residue-aware clash detection
    distances = calculate_distances(metal_pos, coord_oxygens, coord_details)

    # Calculate geometry RMSD
    rmsd = calculate_geometry_rmsd(coord_oxygens, metal_pos, expected_geometry)

    # Perform validation checks
    issues = []
    warnings = []

    coordination_number = len(coord_oxygens)

    # Check 1: Coordination number
    if expected_cn and coordination_number != expected_cn:
        issues.append(f"Coordination number {coordination_number} != expected {expected_cn}")
    if not (VALIDATION_THRESHOLDS["cn_min"] <= coordination_number <= VALIDATION_THRESHOLDS["cn_max"]):
        issues.append(f"Coordination number {coordination_number} outside valid range")

    # Check 2: Metal-Oxygen distances
    if distances["ln_o_distances"]:
        for i, dist in enumerate(distances["ln_o_distances"]):
            if dist < VALIDATION_THRESHOLDS["ln_o_min"]:
                issues.append(f"Ln-O distance {dist:.2f} Å too short (< {VALIDATION_THRESHOLDS['ln_o_min']} Å)")
            elif dist > VALIDATION_THRESHOLDS["ln_o_max"]:
                issues.append(f"Ln-O distance {dist:.2f} Å too long (> {VALIDATION_THRESHOLDS['ln_o_max']} Å)")
            elif dist < VALIDATION_THRESHOLDS["ln_o_ideal_min"]:
                warnings.append(f"Ln-O distance {dist:.2f} Å below ideal range")
            elif dist > VALIDATION_THRESHOLDS["ln_o_ideal_max"]:
                warnings.append(f"Ln-O distance {dist:.2f} Å above ideal range")

    # Check 3a: Intra-carboxylate O-O distances (expected ~2.21 Å for OE1-OE2 or OD1-OD2)
    # These are NOT clashes - they're the natural geometry of bidentate carboxylates
    INTRA_CARBOXYLATE_EXPECTED = 2.21  # Å
    INTRA_CARBOXYLATE_TOLERANCE = 0.3  # Å
    for intra in distances.get("o_o_intra_residue", []):
        dist = intra["distance"]
        if abs(dist - INTRA_CARBOXYLATE_EXPECTED) > INTRA_CARBOXYLATE_TOLERANCE:
            warnings.append(
                f"Intra-carboxylate {intra.get('residue', '?')} {intra.get('atoms', '?')} "
                f"distance {dist:.2f} Å deviates from expected ~{INTRA_CARBOXYLATE_EXPECTED} Å"
            )

    # Check 3b: Inter-residue O-O distances (steric clashes between different residues)
    # This is what we actually care about for clash detection
    inter_clash_count = 0
    if distances.get("o_o_inter_residue"):
        for inter in distances["o_o_inter_residue"]:
            dist = inter["distance"]
            if dist < VALIDATION_THRESHOLDS["o_o_clash"]:
                inter_clash_count += 1
                res1 = inter.get("residue1", "?")
                res2 = inter.get("residue2", "?")
                issues.append(f"Inter-residue O-O clash: {res1} - {res2} = {dist:.2f} Å < {VALIDATION_THRESHOLDS['o_o_clash']} Å")
            elif dist < VALIDATION_THRESHOLDS["o_o_warn"]:
                res1 = inter.get("residue1", "?")
                res2 = inter.get("residue2", "?")
                warnings.append(f"Inter-residue O-O distance {res1} - {res2} = {dist:.2f} Å is close to clash threshold")

    # Use inter-residue clash count (not total O-O clashes which includes intra-carboxylate)
    clash_count = inter_clash_count

    # Check 4: Geometry RMSD
    geometry_quality = "unknown"
    if rmsd is not None:
        if rmsd <= VALIDATION_THRESHOLDS["rmsd_excellent"]:
            geometry_quality = "excellent"
        elif rmsd <= VALIDATION_THRESHOLDS["rmsd_good"]:
            geometry_quality = "good"
        elif rmsd <= VALIDATION_THRESHOLDS["rmsd_acceptable"]:
            geometry_quality = "acceptable"
            warnings.append(f"Geometry RMSD {rmsd:.2f} Å is marginally acceptable")
        else:
            geometry_quality = "poor"
            issues.append(f"Geometry RMSD {rmsd:.2f} Å indicates poor coordination geometry")

    # Determine overall validation status
    passes_validation = len(issues) == 0

    return {
        "valid": passes_validation,
        "coordination_number": coordination_number,
        "metal_type": parsed["metal_type"],

        # Distance metrics
        "ln_o_mean": distances["ln_o_mean"],
        "ln_o_min": distances["ln_o_min"],
        "ln_o_max": distances["ln_o_max"],
        "o_o_min": distances["o_o_min"],
        "o_o_max": distances["o_o_max"],
        "o_o_mean": distances["o_o_mean"],
        "o_o_inter_min": distances.get("o_o_inter_min"),
        "o_o_inter_max": distances.get("o_o_inter_max"),
        "o_o_clash_count": clash_count,

        # Geometry metrics
        "geometry_rmsd": rmsd,
        "geometry_quality": geometry_quality,
        "target_geometry": expected_geometry,

        # Detailed lists
        "ln_o_distances": distances["ln_o_distances"],
        "o_o_distances": distances["o_o_distances"],
        "o_o_intra_residue": distances.get("o_o_intra_residue", []),
        "o_o_inter_residue": distances.get("o_o_inter_residue", []),
        "coordinating_atoms": coord_details,

        # Issues and warnings
        "issues": issues,
        "warnings": warnings,

        # Summary
        "summary": _generate_summary(
            passes_validation, coordination_number, distances, rmsd, issues, warnings
        ),
    }


def _generate_summary(
    passes: bool,
    cn: int,
    distances: Dict,
    rmsd: Optional[float],
    issues: List[str],
    warnings: List[str],
) -> str:
    """Generate a human-readable summary of validation results."""
    status = "PASS" if passes else "FAIL"

    lines = [
        f"Validation Status: {status}",
        f"Coordination Number: {cn}",
    ]

    if distances["ln_o_mean"]:
        lines.append(f"Ln-O Average: {distances['ln_o_mean']:.2f} Å")

    if distances["o_o_min"]:
        lines.append(f"O-O Minimum: {distances['o_o_min']:.2f} Å")

    if rmsd:
        lines.append(f"Geometry RMSD: {rmsd:.2f} Å")

    if issues:
        lines.append(f"Issues ({len(issues)}): " + "; ".join(issues[:3]))

    if warnings:
        lines.append(f"Warnings ({len(warnings)})")

    return " | ".join(lines)


def quick_validate(pdb_content: str) -> bool:
    """
    Quick validation check - returns True if no critical issues.

    Use this for fast filtering of templates before full validation.

    Args:
        pdb_content: PDB file content

    Returns:
        True if template passes basic validation
    """
    result = validate_coordination_geometry(pdb_content)
    return result.get("valid", False)


def get_validation_metrics(pdb_content: str) -> Dict[str, float]:
    """
    Get key metrics for comparison/logging.

    Args:
        pdb_content: PDB file content

    Returns:
        Dictionary of key numeric metrics
    """
    result = validate_coordination_geometry(pdb_content)
    return {
        "coordination_number": result.get("coordination_number", 0),
        "ln_o_mean": result.get("ln_o_mean", 0),
        "o_o_min": result.get("o_o_min", 0),
        "o_o_clash_count": result.get("o_o_clash_count", 0),
        "geometry_rmsd": result.get("geometry_rmsd", 999),
        "passes": 1.0 if result.get("valid") else 0.0,
    }
