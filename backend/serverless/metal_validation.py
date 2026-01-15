"""
Metal Coordination Validation Pipeline

Comprehensive validation of metal binding sites using coordination.py analysis.
Provides lanthanide-specific validation, TEBL readiness checks, and improvement suggestions.
"""

import sys
import os
import math
import numpy as np
from typing import Dict, List, Any, Optional, Tuple

# Add utils to path for coordination module
# Try multiple possible locations for the utils module
_utils_paths = [
    os.path.join(os.path.dirname(__file__), '..', 'utils'),  # Local development
    os.path.join(os.path.dirname(__file__), 'utils'),  # Docker container
    '/app/utils',  # Docker absolute path
]
for _path in _utils_paths:
    if os.path.exists(_path) and _path not in sys.path:
        sys.path.insert(0, _path)
        break

try:
    from coordination import (
        analyze_coordination_geometry,
        suggest_lanthanide_conversion,
        get_coordinating_atoms,
        METAL_PROPERTIES,
    )
    COORDINATION_AVAILABLE = True
except ImportError:
    COORDINATION_AVAILABLE = False


# Lanthanide success criteria based on Caldwell et al. 2020
LANTHANIDE_CRITERIA = {
    "coordination_number": {
        "min": 6,
        "good": 8,
        "excellent": 9,
    },
    "carboxylate_donors": {
        "min": 4,
        "good": 6,
        "excellent": 8,
    },
    "bond_distance": {
        "min": 2.1,
        "max": 2.8,
        "optimal": (2.3, 2.5),
    },
    "geometry_rmsd": {
        "max": 15.0,       # degrees - tightened to match crystallographic standards
        "good": 10.0,      # (was 15.0)
        "excellent": 5.0,  # (was 10.0) - crystallographic metal sites have <10 deg deviation
    },
    "chain_contribution": {
        "min_per_chain": 2,  # Each chain should contribute at least 2 donors
    },
    "trp_antenna_distance": {
        "min": 3.5,
        "max": 6.0,
        "optimal": (4.0, 5.5),
    },
}


def validate_lanthanide_site(
    pdb_content: str,
    metal: str,
    metal_chain: str = "L",
    metal_resnum: int = 1,
    target_coordination: int = 8,
    check_tebl: bool = False,
    distance_cutoff: float = 3.5,
) -> Dict[str, Any]:
    """
    Comprehensive validation of lanthanide binding site using coordination.py.

    Args:
        pdb_content: PDB content string
        metal: Metal code (TB, GD, EU, etc.)
        metal_chain: Chain ID of metal
        metal_resnum: Residue number of metal
        target_coordination: Target coordination number
        check_tebl: Whether to check for TEBL-ready Trp antenna
        distance_cutoff: Maximum distance for coordination (Å)

    Returns:
        Validation results dictionary with:
        - success: bool
        - coordination_number: int
        - geometry_type: str
        - geometry_rmsd: float
        - donor_types: dict
        - bond_distances: dict
        - chain_contribution: dict
        - quality_score: float (0-100)
        - suggestions: list
        - tebl_ready: bool (if check_tebl=True)
    """
    if not COORDINATION_AVAILABLE:
        return {
            "success": False,
            "error": "coordination.py module not available",
        }

    # Run coordination geometry analysis
    analysis = analyze_coordination_geometry(
        pdb_content=pdb_content,
        metal_chain=metal_chain,
        metal_residue=metal,
        metal_resnum=metal_resnum,
        distance_cutoff=distance_cutoff,
    )

    if not analysis.get("success"):
        return analysis

    # Extract key metrics
    coord_num = analysis["coordination"]["number"]
    geometry = analysis["coordination"]["geometry"]
    geometry_rmsd = analysis["coordination"]["geometry_rmsd"]
    donor_types = analysis["donor_analysis"]["types"]
    avg_distance = analysis["bond_analysis"]["average_distance"]

    # Count donors by chain
    chain_contribution = {"A": 0, "B": 0, "other": 0}
    for atom in analysis["coordination"]["coordinating_atoms"]:
        chain = atom["chain"]
        if chain in chain_contribution:
            chain_contribution[chain] += 1
        else:
            chain_contribution["other"] += 1

    # Count carboxylate donors
    carboxylate_count = sum(
        count for dtype, count in donor_types.items()
        if "carboxylate" in dtype
    )

    # Calculate quality score (0-100)
    quality_score = _calculate_quality_score(
        coord_num=coord_num,
        target_coordination=target_coordination,
        carboxylate_count=carboxylate_count,
        geometry_rmsd=geometry_rmsd,
        avg_distance=avg_distance,
        chain_contribution=chain_contribution,
        metal=metal,
    )

    # Generate suggestions
    suggestions = analysis.get("suggestions", [])

    # Add heterodimer-specific suggestions
    if chain_contribution["A"] < 2 or chain_contribution["B"] < 2:
        suggestions.append(
            f"Both chains should contribute ≥2 donors. "
            f"Current: A={chain_contribution['A']}, B={chain_contribution['B']}. "
            "Consider adjusting design constraints."
        )

    if coord_num < target_coordination:
        suggestions.append(
            f"Coordination ({coord_num}) below target ({target_coordination}). "
            "Use iterative refinement with suggest_lanthanide_conversion()."
        )

    # Check TEBL readiness if requested
    tebl_result = None
    if check_tebl:
        tebl_result = calculate_tebl_readiness(
            pdb_content=pdb_content,
            metal_chain=metal_chain,
            metal_resnum=metal_resnum,
        )

    result = {
        "success": True,
        "metal": metal,
        "coordination_number": coord_num,
        "target_coordination": target_coordination,
        "geometry_type": geometry,
        "geometry_rmsd": geometry_rmsd,
        "donor_types": donor_types,
        "carboxylate_count": carboxylate_count,
        "bond_distances": {
            "average": avg_distance,
            "min": analysis["bond_analysis"]["min_distance"],
            "max": analysis["bond_analysis"]["max_distance"],
        },
        "chain_contribution": chain_contribution,
        "quality_score": quality_score,
        "quality_rating": _get_quality_rating(quality_score),
        "suggestions": suggestions,
        "coordinating_atoms": analysis["coordination"]["coordinating_atoms"],
    }

    if tebl_result:
        result["tebl_ready"] = tebl_result["has_antenna"]
        result["tebl_details"] = tebl_result

    return result


def _calculate_quality_score(
    coord_num: int,
    target_coordination: int,
    carboxylate_count: int,
    geometry_rmsd: float,
    avg_distance: float,
    chain_contribution: Dict[str, int],
    metal: str,
) -> float:
    """
    Calculate quality score (0-100) based on multiple metrics.
    """
    score = 0.0

    # Coordination number (40 points max)
    coord_ratio = min(coord_num / target_coordination, 1.0)
    score += coord_ratio * 40

    # Carboxylate donors (20 points max)
    carbox_target = 6  # Good threshold
    carbox_ratio = min(carboxylate_count / carbox_target, 1.0)
    score += carbox_ratio * 20

    # Geometry RMSD (15 points max)
    criteria = LANTHANIDE_CRITERIA["geometry_rmsd"]
    if geometry_rmsd <= criteria["excellent"]:
        score += 15
    elif geometry_rmsd <= criteria["good"]:
        score += 10
    elif geometry_rmsd <= criteria["max"]:
        score += 5

    # Bond distance (15 points max)
    optimal = LANTHANIDE_CRITERIA["bond_distance"]["optimal"]
    if optimal[0] <= avg_distance <= optimal[1]:
        score += 15
    elif LANTHANIDE_CRITERIA["bond_distance"]["min"] <= avg_distance <= LANTHANIDE_CRITERIA["bond_distance"]["max"]:
        score += 8

    # Chain contribution balance (10 points max)
    min_per_chain = LANTHANIDE_CRITERIA["chain_contribution"]["min_per_chain"]
    if chain_contribution["A"] >= min_per_chain and chain_contribution["B"] >= min_per_chain:
        # Bonus for balanced contribution
        balance = min(chain_contribution["A"], chain_contribution["B"]) / max(
            chain_contribution["A"], chain_contribution["B"], 1
        )
        score += 10 * balance

    return round(min(score, 100), 1)


def _get_quality_rating(score: float) -> str:
    """Convert score to rating."""
    if score >= 90:
        return "excellent"
    elif score >= 75:
        return "good"
    elif score >= 60:
        return "acceptable"
    elif score >= 40:
        return "needs_improvement"
    else:
        return "poor"


def calculate_tebl_readiness(
    pdb_content: str,
    metal_chain: str = "L",
    metal_resnum: int = 1,
) -> Dict[str, Any]:
    """
    Check if design has proper Trp antenna for TEBL assay.

    Tryptophan-Enhanced Terbium Luminescence (TEBL) requires:
    - Trp residue within 4-6 Å of lanthanide
    - Proper orientation for energy transfer

    Args:
        pdb_content: PDB content string
        metal_chain: Chain ID of metal
        metal_resnum: Residue number of metal

    Returns:
        TEBL readiness assessment
    """
    # Find metal position
    metal_pos = None
    trp_residues = []

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
            pos = np.array([x, y, z])

            # Check for metal
            if line.startswith('HETATM'):
                if chain == metal_chain and res_num == metal_resnum:
                    metal_pos = pos

            # Check for Trp NE1 (indole nitrogen - key for energy transfer)
            if res_name == "TRP" and atom_name == "NE1":
                trp_residues.append({
                    "chain": chain,
                    "res_num": res_num,
                    "ne1_pos": pos,
                })
        except (ValueError, IndexError):
            continue

    if metal_pos is None:
        return {
            "has_antenna": False,
            "error": "Metal not found",
            "trp_count": len(trp_residues),
        }

    if not trp_residues:
        return {
            "has_antenna": False,
            "error": "No Trp residues found",
            "suggestion": "Add Trp residue within 4-6 Å of metal for TEBL",
        }

    # Find closest Trp to metal
    best_trp = None
    best_distance = float('inf')

    for trp in trp_residues:
        distance = np.linalg.norm(trp["ne1_pos"] - metal_pos)
        if distance < best_distance:
            best_distance = distance
            best_trp = trp

    # Evaluate TEBL readiness
    optimal = LANTHANIDE_CRITERIA["trp_antenna_distance"]["optimal"]
    is_optimal = optimal[0] <= best_distance <= optimal[1]
    is_acceptable = (
        LANTHANIDE_CRITERIA["trp_antenna_distance"]["min"] <= best_distance <=
        LANTHANIDE_CRITERIA["trp_antenna_distance"]["max"]
    )

    # Estimate energy transfer efficiency (simplified model)
    # FRET-like distance dependence: E = 1 / (1 + (r/R0)^6)
    # R0 ≈ 4.5 Å for Trp-Tb
    r0 = 4.5
    efficiency = 1.0 / (1.0 + (best_distance / r0) ** 6)

    result = {
        "has_antenna": is_acceptable,
        "is_optimal": is_optimal,
        "trp_chain": best_trp["chain"],
        "trp_resnum": best_trp["res_num"],
        "trp_metal_distance": round(best_distance, 2),
        "energy_transfer_efficiency": round(efficiency, 3),
        "trp_count": len(trp_residues),
    }

    if not is_acceptable:
        result["suggestion"] = (
            f"Trp at {best_distance:.1f}Å from metal. "
            f"Optimal: {optimal[0]}-{optimal[1]}Å. "
            "Reposition Trp for better energy transfer."
        )
    elif not is_optimal:
        result["note"] = (
            f"Trp at {best_distance:.1f}Å. "
            f"Could be optimized to {optimal[0]}-{optimal[1]}Å."
        )

    return result


def get_refinement_parameters(
    validation_result: Dict[str, Any],
    metal: str = "TB",
) -> Dict[str, Any]:
    """
    Get RFD3 parameters for refinement based on validation results.

    Uses suggest_lanthanide_conversion() to generate improvement parameters.

    Args:
        validation_result: Output from validate_lanthanide_site()
        metal: Target metal

    Returns:
        RFD3 refinement parameters
    """
    if not validation_result.get("success"):
        return {"error": "Cannot generate parameters from failed validation"}

    if not COORDINATION_AVAILABLE:
        return {"error": "coordination.py module not available"}

    # Build analysis dict for suggest_lanthanide_conversion
    analysis = {
        "success": True,
        "metal": {"element": metal},
        "coordination": {
            "number": validation_result["coordination_number"],
            "coordinating_atoms": validation_result.get("coordinating_atoms", []),
        },
        "donor_analysis": {
            "types": validation_result["donor_types"],
            "donor_list": [],
        },
    }

    # Get lanthanide conversion suggestions
    suggestions = suggest_lanthanide_conversion(analysis, metal)

    return suggestions.get("rfd3_parameters", {})


def batch_validate(
    designs: List[Dict[str, str]],
    metal: str,
    target_coordination: int = 8,
    check_tebl: bool = False,
) -> Dict[str, Any]:
    """
    Validate multiple designs and rank by quality.

    Args:
        designs: List of dicts with 'pdb_content' key
        metal: Metal code
        target_coordination: Target coordination
        check_tebl: Whether to check TEBL readiness

    Returns:
        Batch validation results with ranking
    """
    results = []

    for i, design in enumerate(designs):
        pdb_content = design.get("pdb_content", "")
        if not pdb_content:
            results.append({
                "design_index": i,
                "success": False,
                "error": "No PDB content",
            })
            continue

        validation = validate_lanthanide_site(
            pdb_content=pdb_content,
            metal=metal,
            target_coordination=target_coordination,
            check_tebl=check_tebl,
        )
        validation["design_index"] = i
        results.append(validation)

    # Sort by quality score
    valid_results = [r for r in results if r.get("success")]
    valid_results.sort(key=lambda x: x.get("quality_score", 0), reverse=True)

    # Calculate statistics
    if valid_results:
        scores = [r["quality_score"] for r in valid_results]
        coord_nums = [r["coordination_number"] for r in valid_results]
        stats = {
            "average_score": round(sum(scores) / len(scores), 1),
            "max_score": max(scores),
            "min_score": min(scores),
            "average_coordination": round(sum(coord_nums) / len(coord_nums), 1),
            "designs_meeting_target": sum(
                1 for r in valid_results
                if r["coordination_number"] >= target_coordination
            ),
        }
    else:
        stats = {}

    return {
        "total_designs": len(designs),
        "valid_designs": len(valid_results),
        "ranked_results": valid_results,
        "failed_results": [r for r in results if not r.get("success")],
        "statistics": stats,
        "best_design_index": valid_results[0]["design_index"] if valid_results else None,
    }
