"""
TEBL (Tryptophan-Enhanced Terbium Luminescence) Analysis Module

Optional add-on for lanthanide binding site analysis.
Based on Caldwell et al. (2020) PNAS methodology.

TEBL provides:
1. Enhanced sensitivity for detecting lanthanide binding
2. Structural validation through distance-dependent energy transfer
3. Binding affinity measurements via luminescence kinetics
"""

import math
import numpy as np
from typing import Dict, List, Any, Optional, Tuple


# Physical constants for TEBL calculations
TEBL_CONSTANTS = {
    # Förster radius (R0) for Trp → Lanthanide energy transfer
    # Based on spectral overlap and transition moments
    "R0_TB": 4.5,   # Å for Terbium
    "R0_EU": 5.0,   # Å for Europium
    "R0_GD": 4.0,   # Å for Gadolinium (paramagnetic, lower efficiency)
    "R0_YB": 3.8,   # Å for Ytterbium

    # Optimal Trp-Ln distance range
    "optimal_distance": (4.0, 5.5),  # Å

    # Trp excitation wavelength
    "trp_excitation": 280,  # nm

    # Lanthanide emission wavelengths
    "emission": {
        "TB": {"peak": 545, "range": (480, 630)},  # Green
        "EU": {"peak": 615, "range": (580, 700)},  # Red
        "GD": None,  # Paramagnetic, not luminescent
        "YB": {"peak": 980, "range": (950, 1050)},  # NIR
    },
}


def predict_tebl_signal(
    pdb_content: str,
    metal: str = "TB",
    metal_chain: str = "L",
    metal_resnum: int = 1,
) -> Dict[str, Any]:
    """
    Predict TEBL signal strength based on structure.

    Uses simplified Förster resonance energy transfer (FRET) model
    to estimate energy transfer efficiency from Trp to lanthanide.

    Args:
        pdb_content: PDB content string
        metal: Lanthanide code (TB, EU, YB)
        metal_chain: Chain ID of metal
        metal_resnum: Residue number of metal

    Returns:
        TEBL signal prediction with:
        - has_antenna: bool
        - predicted_efficiency: float (0-1)
        - signal_strength: str (strong/medium/weak/none)
        - recommendations: list
    """
    # Check if metal is luminescent
    emission = TEBL_CONSTANTS["emission"].get(metal.upper())
    if emission is None:
        return {
            "has_antenna": False,
            "error": f"{metal} is not luminescent. Use TB or EU for TEBL.",
            "suggestion": "Consider Tb³⁺ (green) or Eu³⁺ (red) for luminescence assays.",
        }

    # Find metal position
    metal_pos = _find_metal_position(pdb_content, metal_chain, metal_resnum)
    if metal_pos is None:
        return {
            "has_antenna": False,
            "error": "Metal not found in structure",
        }

    # Find all Trp residues and their positions
    trp_atoms = _find_trp_residues(pdb_content)
    if not trp_atoms:
        return {
            "has_antenna": False,
            "predicted_efficiency": 0.0,
            "signal_strength": "none",
            "recommendations": [
                "No Trp residues found. Add Trp antenna within 4-6 Å of metal.",
                "Use position_trp_antenna() from lanthanide_templates module.",
            ],
        }

    # Calculate energy transfer for each Trp
    r0 = TEBL_CONSTANTS.get(f"R0_{metal.upper()}", 4.5)
    trp_analysis = []

    for trp in trp_atoms:
        ne1_pos = trp.get("ne1_pos")
        if ne1_pos is None:
            continue

        distance = np.linalg.norm(ne1_pos - metal_pos)

        # FRET efficiency: E = 1 / (1 + (r/R0)^6)
        efficiency = 1.0 / (1.0 + (distance / r0) ** 6)

        # Orientation factor (simplified - assume random orientation)
        # In reality, indole ring orientation affects transfer
        kappa2 = 2.0 / 3.0  # Random orientation

        # Adjusted efficiency with orientation
        adjusted_efficiency = efficiency * kappa2

        trp_analysis.append({
            "chain": trp["chain"],
            "residue_number": trp["res_num"],
            "distance": round(distance, 2),
            "raw_efficiency": round(efficiency, 3),
            "adjusted_efficiency": round(adjusted_efficiency, 3),
            "is_optimal": TEBL_CONSTANTS["optimal_distance"][0] <= distance <= TEBL_CONSTANTS["optimal_distance"][1],
        })

    # Sort by efficiency
    trp_analysis.sort(key=lambda x: x["adjusted_efficiency"], reverse=True)
    best_trp = trp_analysis[0]

    # Overall signal strength
    total_efficiency = best_trp["adjusted_efficiency"]
    if total_efficiency >= 0.5:
        signal_strength = "strong"
    elif total_efficiency >= 0.2:
        signal_strength = "medium"
    elif total_efficiency >= 0.05:
        signal_strength = "weak"
    else:
        signal_strength = "very_weak"

    # Generate recommendations
    recommendations = []
    if not best_trp["is_optimal"]:
        opt_range = TEBL_CONSTANTS["optimal_distance"]
        recommendations.append(
            f"Trp at {best_trp['distance']:.1f}Å. "
            f"Move to {opt_range[0]}-{opt_range[1]}Å for optimal signal."
        )

    if total_efficiency < 0.3:
        recommendations.append(
            "Consider adding a second Trp residue for enhanced signal."
        )

    if len(trp_analysis) > 1:
        # Check for potential quenching from distant Trp
        for trp in trp_analysis[1:]:
            if trp["distance"] < 10.0 and trp["adjusted_efficiency"] < 0.1:
                recommendations.append(
                    f"Trp at {trp['chain']}{trp['residue_number']} "
                    f"({trp['distance']:.1f}Å) may cause background. "
                    "Consider mutating to Phe."
                )

    return {
        "has_antenna": True,
        "metal": metal.upper(),
        "emission_peak": emission["peak"],
        "emission_range": emission["range"],
        "best_antenna": best_trp,
        "all_antennas": trp_analysis,
        "predicted_efficiency": round(total_efficiency, 3),
        "signal_strength": signal_strength,
        "recommendations": recommendations,
        "assay_protocol": _get_tebl_protocol(metal),
    }


def _find_metal_position(
    pdb_content: str,
    metal_chain: str,
    metal_resnum: int,
) -> Optional[np.ndarray]:
    """Find metal position in PDB."""
    for line in pdb_content.split('\n'):
        if not line.startswith('HETATM'):
            continue

        try:
            chain = line[21]
            res_num = int(line[22:26].strip())

            if chain == metal_chain and res_num == metal_resnum:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                return np.array([x, y, z])
        except (ValueError, IndexError):
            continue

    return None


def _find_trp_residues(pdb_content: str) -> List[Dict[str, Any]]:
    """Find all Trp residues and their key atoms."""
    trp_dict = {}  # Key: (chain, resnum)

    for line in pdb_content.split('\n'):
        if not line.startswith('ATOM'):
            continue

        try:
            res_name = line[17:20].strip()
            if res_name != "TRP":
                continue

            chain = line[21]
            res_num = int(line[22:26].strip())
            atom_name = line[12:16].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            key = (chain, res_num)
            if key not in trp_dict:
                trp_dict[key] = {
                    "chain": chain,
                    "res_num": res_num,
                }

            # Store key atoms
            if atom_name == "NE1":
                trp_dict[key]["ne1_pos"] = np.array([x, y, z])
            elif atom_name == "CZ2":
                trp_dict[key]["cz2_pos"] = np.array([x, y, z])
            elif atom_name == "CD1":
                trp_dict[key]["cd1_pos"] = np.array([x, y, z])

        except (ValueError, IndexError):
            continue

    return list(trp_dict.values())


def _get_tebl_protocol(metal: str) -> Dict[str, Any]:
    """Get recommended TEBL assay protocol."""
    emission = TEBL_CONSTANTS["emission"].get(metal.upper())

    return {
        "excitation_wavelength": 280,
        "emission_wavelength": emission["peak"] if emission else None,
        "emission_filter": f"{emission['range'][0]}-{emission['range'][1]} nm" if emission else None,
        "delay_time": "50-100 μs (to eliminate Trp fluorescence)",
        "measurement_mode": "Time-resolved fluorescence (TRF)",
        "buffer_recommendations": [
            "Use HEPES or Tris buffer (avoid phosphate - competes for Ln binding)",
            "Include 100 mM NaCl to reduce non-specific interactions",
            "pH 7.0-7.5 optimal for carboxylate coordination",
        ],
        "controls": [
            "Apo-protein (no metal) for baseline",
            "Trp-to-Phe mutant to confirm energy transfer",
            "EDTA displacement for affinity estimation",
        ],
    }


def add_trp_antenna_to_design(
    pdb_content: str,
    metal_chain: str = "L",
    metal_resnum: int = 1,
    target_chain: str = "A",
    target_distance: float = 4.5,
) -> Dict[str, Any]:
    """
    Suggest mutation sites to add Trp antenna near metal site.

    Analyzes the structure to find suitable residues for Trp mutation
    that would place the indole ring at optimal distance from metal.

    Args:
        pdb_content: PDB content string
        metal_chain: Chain ID of metal
        metal_resnum: Residue number of metal
        target_chain: Which chain to mutate
        target_distance: Target Trp-metal distance (Å)

    Returns:
        Mutation suggestions with predicted efficiency
    """
    metal_pos = _find_metal_position(pdb_content, metal_chain, metal_resnum)
    if metal_pos is None:
        return {
            "success": False,
            "error": "Metal not found",
        }

    # Find candidate residues on target chain
    candidates = []

    for line in pdb_content.split('\n'):
        if not line.startswith('ATOM'):
            continue

        try:
            chain = line[21]
            if chain != target_chain:
                continue

            res_name = line[17:20].strip()
            res_num = int(line[22:26].strip())
            atom_name = line[12:16].strip()

            # Use CA position to estimate side chain position
            if atom_name != "CA":
                continue

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            ca_pos = np.array([x, y, z])

            # Estimate where Trp NE1 would be (~6 Å from CA toward metal)
            direction = metal_pos - ca_pos
            direction = direction / np.linalg.norm(direction)
            estimated_ne1_pos = ca_pos + direction * 6.0
            estimated_distance = np.linalg.norm(estimated_ne1_pos - metal_pos)

            # Check if this would put Trp in good range
            if 2.0 <= estimated_distance <= 8.0:
                candidates.append({
                    "chain": chain,
                    "residue_number": res_num,
                    "current_residue": res_name,
                    "ca_distance_to_metal": round(np.linalg.norm(ca_pos - metal_pos), 2),
                    "estimated_ne1_distance": round(estimated_distance, 2),
                    "is_optimal": 4.0 <= estimated_distance <= 5.5,
                })

        except (ValueError, IndexError):
            continue

    # Sort by how close to target distance
    candidates.sort(
        key=lambda x: abs(x["estimated_ne1_distance"] - target_distance)
    )

    # Filter out problematic residues
    good_candidates = []
    problematic = {"CYS", "GLY", "PRO"}  # Hard to mutate or structural
    for c in candidates:
        if c["current_residue"] not in problematic:
            good_candidates.append(c)

    if not good_candidates:
        return {
            "success": False,
            "error": "No suitable mutation sites found",
            "all_candidates": candidates[:10],
        }

    return {
        "success": True,
        "best_candidate": good_candidates[0],
        "top_candidates": good_candidates[:5],
        "mutation_command": (
            f"Mutate {good_candidates[0]['current_residue']}"
            f"{good_candidates[0]['residue_number']} to TRP"
        ),
        "expected_efficiency": _estimate_efficiency(
            good_candidates[0]["estimated_ne1_distance"]
        ),
    }


def _estimate_efficiency(distance: float, r0: float = 4.5) -> float:
    """Estimate FRET efficiency for given distance."""
    return round(1.0 / (1.0 + (distance / r0) ** 6), 3)


def estimate_binding_affinity(
    tebl_result: Dict[str, Any],
    reference_signal: float = 1.0,
    measured_signal: float = 1.0,
) -> Dict[str, Any]:
    """
    Estimate binding affinity from TEBL signal (placeholder for experimental use).

    In practice, this would use titration data to calculate Kd.
    Here we provide the framework for the calculation.

    Args:
        tebl_result: Output from predict_tebl_signal()
        reference_signal: Reference signal at saturation
        measured_signal: Measured signal

    Returns:
        Affinity estimation (placeholder)
    """
    if not tebl_result.get("has_antenna"):
        return {"error": "No antenna - cannot estimate affinity"}

    # Placeholder - real calculation requires experimental titration
    fractional_saturation = measured_signal / reference_signal

    return {
        "fractional_saturation": fractional_saturation,
        "note": "Actual Kd requires titration experiment",
        "recommended_method": "Competitive displacement with EDTA",
        "caldwell_reference": (
            "Caldwell et al. achieved subfemtomolar Kd "
            "(< 10^-18 M) with 4 Glu coordination"
        ),
    }
