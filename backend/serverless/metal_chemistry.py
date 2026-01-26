# metal_chemistry.py
"""
Metal Chemistry Database Module

Provides HSAB (Hard-Soft Acid-Base) theory-compliant metal coordination chemistry
for protein design with metal-binding sites.

Based on established coordination chemistry principles:
- Hard acids prefer hard bases (O donors like Glu, Asp, Asn)
- Soft acids prefer soft bases (S donors like Cys, Met)
- Borderline acids can coordinate both

References:
- Pearson RG. Hard and Soft Acids and Bases. JACS 1963.
- Holm RH et al. Structural and Functional Aspects of Metal Sites in Biology. Chem Rev 1996.
"""

import logging
from typing import Dict, List, Tuple, Any, Optional

logger = logging.getLogger(__name__)


# =============================================================================
# METAL DATABASE
# =============================================================================

METAL_DATABASE: Dict[str, Dict[str, Any]] = {
    # -------------------------------------------------------------------------
    # TRANSITION METALS - BORDERLINE ACIDS
    # -------------------------------------------------------------------------
    "ZN": {
        "name": "Zinc",
        "hsab_class": {2: "borderline"},
        "default_oxidation": 2,
        "preferred_donors": {
            2: {
                "catalytic": {"N": 2.0, "S": 2.0, "O": 1.0},
                "structural": {"S": 2.5, "N": 1.5, "O": 0.5},
            }
        },
        "bond_distances": {  # In Angstroms
            2: {"N": (1.95, 2.20), "S": (2.25, 2.45), "O": (1.95, 2.20)}
        },
        "coordination_numbers": {2: (4, 6)},  # (min, max)
        "common_residues": ["HIS", "CYS", "GLU", "ASP"],
        "description": "Borderline Lewis acid, common in zinc fingers and catalytic sites",
    },

    "FE": {
        "name": "Iron",
        "hsab_class": {2: "borderline", 3: "hard"},
        "default_oxidation": 2,
        "preferred_donors": {
            2: {
                "catalytic": {"N": 2.0, "S": 1.5, "O": 1.5},
                "structural": {"N": 2.0, "O": 1.5, "S": 1.0},
            },
            3: {
                "catalytic": {"O": 2.5, "N": 1.5, "S": 0.5},
                "structural": {"O": 2.5, "N": 1.0, "S": 0.0},
            }
        },
        "bond_distances": {
            2: {"N": (1.95, 2.20), "S": (2.20, 2.40), "O": (1.95, 2.20)},
            3: {"N": (1.90, 2.15), "S": (2.15, 2.35), "O": (1.90, 2.10)}
        },
        "coordination_numbers": {2: (4, 6), 3: (4, 6)},
        "common_residues": ["HIS", "CYS", "GLU", "ASP", "MET"],
        "description": "Fe2+ is borderline, Fe3+ is hard acid",
    },

    "CU": {
        "name": "Copper",
        "hsab_class": {1: "soft", 2: "borderline"},
        "default_oxidation": 2,
        "preferred_donors": {
            1: {
                "catalytic": {"S": 3.0, "N": 1.5, "O": 0.5},
                "structural": {"S": 3.0, "N": 1.0, "O": 0.0},
            },
            2: {
                "catalytic": {"N": 2.0, "S": 1.5, "O": 1.5},
                "structural": {"N": 2.0, "O": 1.5, "S": 1.0},
            }
        },
        "bond_distances": {
            1: {"S": (2.10, 2.35), "N": (1.90, 2.15), "O": (1.90, 2.10)},
            2: {"N": (1.90, 2.10), "S": (2.15, 2.40), "O": (1.90, 2.10)}
        },
        "coordination_numbers": {1: (2, 4), 2: (4, 6)},
        "common_residues": ["HIS", "CYS", "MET", "GLU", "ASP"],
        "description": "Cu1+ is soft, Cu2+ is borderline",
    },

    "MN": {
        "name": "Manganese",
        "hsab_class": {2: "borderline"},
        "default_oxidation": 2,
        "preferred_donors": {
            2: {
                "catalytic": {"O": 2.0, "N": 1.5, "S": 0.5},
                "structural": {"O": 2.0, "N": 1.0, "S": 0.0},
            }
        },
        "bond_distances": {
            2: {"N": (2.00, 2.25), "O": (2.00, 2.25), "S": (2.30, 2.50)}
        },
        "coordination_numbers": {2: (4, 6)},
        "common_residues": ["HIS", "GLU", "ASP"],
        "description": "Borderline acid, common in oxidoreductases",
    },

    "CO": {
        "name": "Cobalt",
        "hsab_class": {2: "borderline"},
        "default_oxidation": 2,
        "preferred_donors": {
            2: {
                "catalytic": {"N": 2.0, "S": 1.5, "O": 1.5},
                "structural": {"N": 2.0, "O": 1.5, "S": 1.0},
            }
        },
        "bond_distances": {
            2: {"N": (1.90, 2.15), "S": (2.20, 2.40), "O": (1.90, 2.15)}
        },
        "coordination_numbers": {2: (4, 6)},
        "common_residues": ["HIS", "CYS", "GLU", "ASP"],
        "description": "Borderline acid, common in vitamin B12",
    },

    "NI": {
        "name": "Nickel",
        "hsab_class": {2: "borderline"},
        "default_oxidation": 2,
        "preferred_donors": {
            2: {
                "catalytic": {"S": 2.0, "N": 2.0, "O": 1.0},
                "structural": {"S": 2.0, "N": 1.5, "O": 0.5},
            }
        },
        "bond_distances": {
            2: {"N": (1.85, 2.10), "S": (2.15, 2.35), "O": (1.90, 2.10)}
        },
        "coordination_numbers": {2: (4, 6)},
        "common_residues": ["HIS", "CYS", "GLU", "ASP"],
        "description": "Borderline acid, found in ureases and hydrogenases",
    },

    # -------------------------------------------------------------------------
    # ALKALINE EARTH METALS - HARD ACIDS
    # -------------------------------------------------------------------------
    "CA": {
        "name": "Calcium",
        "hsab_class": {2: "hard"},
        "default_oxidation": 2,
        "preferred_donors": {
            2: {
                "catalytic": {"O": 3.0, "N": 0.5, "S": -5.0},
                "structural": {"O": 3.0, "N": 0.0, "S": -5.0},
            }
        },
        "bond_distances": {
            2: {"O": (2.30, 2.60), "N": (2.40, 2.70)}
        },
        "coordination_numbers": {2: (6, 8)},
        "common_residues": ["GLU", "ASP", "ASN", "GLN"],
        "description": "Hard acid, prefers O donors exclusively",
    },

    "MG": {
        "name": "Magnesium",
        "hsab_class": {2: "hard"},
        "default_oxidation": 2,
        "preferred_donors": {
            2: {
                "catalytic": {"O": 3.0, "N": 0.5, "S": -5.0},
                "structural": {"O": 3.0, "N": 0.0, "S": -5.0},
            }
        },
        "bond_distances": {
            2: {"O": (2.00, 2.30), "N": (2.10, 2.40)}
        },
        "coordination_numbers": {2: (4, 6)},
        "common_residues": ["GLU", "ASP", "ASN"],
        "description": "Hard acid, essential for ATP binding",
    },

    # -------------------------------------------------------------------------
    # LANTHANIDES - HARD ACIDS (O preferred, N allowed with lower priority, exclude Cys)
    # -------------------------------------------------------------------------
    "TB": {
        "name": "Terbium",
        "hsab_class": {3: "hard"},
        "default_oxidation": 3,
        "preferred_donors": {
            3: {
                # O strongly preferred (Glu/Asp), N allowed lower priority (His), S excluded (Cys)
                "catalytic": {"O": 3.0, "N": 0.5, "S": -5.0},
                "structural": {"O": 3.0, "N": 0.5, "S": -5.0},
            }
        },
        "bond_distances": {
            3: {"O": (2.30, 2.60), "N": (2.40, 2.70)}
        },
        "coordination_numbers": {3: (8, 9)},
        "common_residues": ["GLU", "ASP", "ASN", "GLN", "HIS"],
        "description": "Hard acid lanthanide, O donors preferred, His allowed, exclude Cys",
    },

    "EU": {
        "name": "Europium",
        "hsab_class": {3: "hard"},
        "default_oxidation": 3,
        "preferred_donors": {
            3: {
                "catalytic": {"O": 3.0, "N": 0.5, "S": -5.0},
                "structural": {"O": 3.0, "N": 0.5, "S": -5.0},
            }
        },
        "bond_distances": {
            3: {"O": (2.35, 2.65), "N": (2.45, 2.75)}
        },
        "coordination_numbers": {3: (8, 9)},
        "common_residues": ["GLU", "ASP", "ASN", "GLN", "HIS"],
        "description": "Hard acid lanthanide, luminescent, O preferred, His allowed",
    },

    "GD": {
        "name": "Gadolinium",
        "hsab_class": {3: "hard"},
        "default_oxidation": 3,
        "preferred_donors": {
            3: {
                "catalytic": {"O": 3.0, "N": 0.5, "S": -5.0},
                "structural": {"O": 3.0, "N": 0.5, "S": -5.0},
            }
        },
        "bond_distances": {
            3: {"O": (2.30, 2.60), "N": (2.40, 2.70)}
        },
        "coordination_numbers": {3: (8, 9)},
        "common_residues": ["GLU", "ASP", "ASN", "GLN", "HIS"],
        "description": "Hard acid lanthanide, MRI contrast agent, O preferred, His allowed",
    },

    "LA": {
        "name": "Lanthanum",
        "hsab_class": {3: "hard"},
        "default_oxidation": 3,
        "preferred_donors": {
            3: {
                "catalytic": {"O": 3.0, "N": 0.5, "S": -5.0},
                "structural": {"O": 3.0, "N": 0.5, "S": -5.0},
            }
        },
        "bond_distances": {
            3: {"O": (2.40, 2.70), "N": (2.50, 2.80)}
        },
        "coordination_numbers": {3: (8, 10)},
        "common_residues": ["GLU", "ASP", "ASN", "GLN", "HIS"],
        "description": "Hard acid lanthanide, largest ionic radius, O preferred, His allowed",
    },

    "CE": {
        "name": "Cerium",
        "hsab_class": {3: "hard"},
        "default_oxidation": 3,
        "preferred_donors": {
            3: {
                "catalytic": {"O": 3.0, "N": 0.5, "S": -5.0},
                "structural": {"O": 3.0, "N": 0.5, "S": -5.0},
            }
        },
        "bond_distances": {
            3: {"O": (2.35, 2.65), "N": (2.45, 2.75)}
        },
        "coordination_numbers": {3: (8, 10)},
        "common_residues": ["GLU", "ASP", "ASN", "GLN", "HIS"],
        "description": "Hard acid lanthanide, O preferred, His allowed",
    },

    "SM": {
        "name": "Samarium",
        "hsab_class": {3: "hard"},
        "default_oxidation": 3,
        "preferred_donors": {
            3: {
                "catalytic": {"O": 3.0, "N": 0.5, "S": -5.0},
                "structural": {"O": 3.0, "N": 0.5, "S": -5.0},
            }
        },
        "bond_distances": {
            3: {"O": (2.30, 2.60), "N": (2.40, 2.70)}
        },
        "coordination_numbers": {3: (8, 9)},
        "common_residues": ["GLU", "ASP", "ASN", "GLN", "HIS"],
        "description": "Hard acid lanthanide, O preferred, His allowed",
    },

    "YB": {
        "name": "Ytterbium",
        "hsab_class": {3: "hard"},
        "default_oxidation": 3,
        "preferred_donors": {
            3: {
                "catalytic": {"O": 3.0, "N": 0.5, "S": -5.0},
                "structural": {"O": 3.0, "N": 0.5, "S": -5.0},
            }
        },
        "bond_distances": {
            3: {"O": (2.25, 2.55), "N": (2.35, 2.65)}
        },
        "coordination_numbers": {3: (8, 9)},
        "common_residues": ["GLU", "ASP", "ASN", "GLN", "HIS"],
        "description": "Hard acid lanthanide, smallest lanthanide, O preferred, His allowed",
    },

    "DY": {
        "name": "Dysprosium",
        "hsab_class": {3: "hard"},
        "default_oxidation": 3,
        "preferred_donors": {
            3: {
                "catalytic": {"O": 3.0, "N": 0.5, "S": -5.0},
                "structural": {"O": 3.0, "N": 0.5, "S": -5.0},
            }
        },
        "bond_distances": {
            3: {"O": (2.28, 2.58), "N": (2.38, 2.68)}
        },
        "coordination_numbers": {3: (8, 9)},
        "common_residues": ["GLU", "ASP", "ASN", "GLN", "HIS"],
        "description": "Hard acid lanthanide, SMM applications, O preferred, His allowed",
    },
}


# =============================================================================
# RESIDUE TO DONOR MAPPING
# =============================================================================

RESIDUE_DONORS: Dict[str, str] = {
    # Oxygen donors (carboxylate)
    "GLU": "O",
    "ASP": "O",
    "E": "O",
    "D": "O",
    # Oxygen donors (amide)
    "ASN": "O",
    "GLN": "O",
    "N": "O",
    "Q": "O",
    # Nitrogen donors (imidazole)
    "HIS": "N",
    "H": "N",
    # Sulfur donors (thiolate)
    "CYS": "S",
    "C": "S",
    # Sulfur donors (thioether)
    "MET": "S",
    "M": "S",
    # Hydroxyl (weak coordination)
    "SER": "O",
    "THR": "O",
    "TYR": "O",
    "S": "O",
    "T": "O",
    "Y": "O",
}

# One-letter to three-letter code mapping
AA_1TO3: Dict[str, str] = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
    "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
}

AA_3TO1: Dict[str, str] = {v: k for k, v in AA_1TO3.items()}


# =============================================================================
# PUBLIC FUNCTIONS
# =============================================================================

def get_hsab_class(metal: str, oxidation_state: int) -> str:
    """
    Get HSAB (Hard-Soft Acid-Base) classification for a metal ion.

    Args:
        metal: Metal element symbol (e.g., "ZN", "FE", "TB")
        oxidation_state: Oxidation state (e.g., 2 for Zn2+, 3 for Fe3+)

    Returns:
        HSAB class: "hard", "borderline", or "soft"

    Raises:
        ValueError: If metal is unknown or oxidation state is invalid
    """
    metal = metal.upper()

    if metal not in METAL_DATABASE:
        raise ValueError(f"Unknown metal: {metal}. Valid metals: {list(METAL_DATABASE.keys())}")

    hsab_data = METAL_DATABASE[metal]["hsab_class"]

    if oxidation_state not in hsab_data:
        valid_states = list(hsab_data.keys())
        raise ValueError(
            f"Invalid oxidation state {oxidation_state} for {metal}. "
            f"Valid states: {valid_states}"
        )

    return hsab_data[oxidation_state]


def get_preferred_donors(
    metal: str,
    oxidation_state: int,
    site_type: str = "catalytic"
) -> Dict[str, float]:
    """
    Get preferred donor atoms with weights for a metal ion.

    Args:
        metal: Metal element symbol
        oxidation_state: Oxidation state
        site_type: "catalytic" or "structural"

    Returns:
        Dict mapping donor element (N, O, S) to preference weight.
        Positive = preferred, negative = penalized.
    """
    metal = metal.upper()

    if metal not in METAL_DATABASE:
        raise ValueError(f"Unknown metal: {metal}")

    data = METAL_DATABASE[metal]
    ox_state = oxidation_state if oxidation_state in data["preferred_donors"] else data["default_oxidation"]

    if ox_state not in data["preferred_donors"]:
        raise ValueError(f"Invalid oxidation state {oxidation_state} for {metal}")

    donors_by_site = data["preferred_donors"][ox_state]
    site = site_type if site_type in donors_by_site else "catalytic"

    return donors_by_site[site].copy()


def get_bond_distance_range(
    metal: str,
    donor_element: str,
    oxidation_state: int
) -> Tuple[float, float]:
    """
    Get expected bond distance range for a metal-donor pair.

    Args:
        metal: Metal element symbol
        donor_element: Donor atom type (N, O, or S)
        oxidation_state: Oxidation state

    Returns:
        Tuple of (min_distance, max_distance) in Angstroms
    """
    metal = metal.upper()
    donor = donor_element.upper()

    if metal not in METAL_DATABASE:
        raise ValueError(f"Unknown metal: {metal}")

    data = METAL_DATABASE[metal]
    ox_state = oxidation_state if oxidation_state in data["bond_distances"] else data["default_oxidation"]

    if ox_state not in data["bond_distances"]:
        raise ValueError(f"Invalid oxidation state {oxidation_state} for {metal}")

    distances = data["bond_distances"][ox_state]

    if donor not in distances:
        # Return a default range if donor not found
        return (2.0, 3.0)

    return distances[donor]


def get_amino_acid_bias(
    metal: str,
    oxidation_state: int,
    site_type: str = "catalytic"
) -> str:
    """
    Get amino acid bias string for LigandMPNN.

    Args:
        metal: Metal element symbol
        oxidation_state: Oxidation state
        site_type: "catalytic" or "structural"

    Returns:
        Bias string in format "E:3.0,D:3.0,C:-5.0"
    """
    metal = metal.upper()
    donors = get_preferred_donors(metal, oxidation_state, site_type)
    hsab_class = get_hsab_class(metal, oxidation_state)

    # Map donor preferences to amino acid biases
    bias_parts = []

    # Oxygen donors: Glu, Asp (carboxylate), Asn, Gln (amide carbonyl)
    if "O" in donors:
        o_weight = donors["O"]
        if o_weight > 0:
            bias_parts.append(f"E:{o_weight:.1f}")
            bias_parts.append(f"D:{o_weight:.1f}")
            # Asn/Gln slightly lower for amide vs carboxylate
            bias_parts.append(f"N:{o_weight * 0.8:.1f}")
            bias_parts.append(f"Q:{o_weight * 0.8:.1f}")

    # Nitrogen donors: His
    if "N" in donors:
        n_weight = donors["N"]
        if n_weight > 0:
            bias_parts.append(f"H:{n_weight:.1f}")

    # Sulfur donors: Cys, Met
    if "S" in donors:
        s_weight = donors["S"]
        # Always include Cys with its weight (even if negative for hard acids)
        bias_parts.append(f"C:{s_weight:.1f}")
        # Met only if soft/borderline
        if s_weight >= 0:
            bias_parts.append(f"M:{s_weight * 0.5:.1f}")

    # For hard acids, ensure Cys is penalized if not already
    if hsab_class == "hard" and "S" not in donors:
        bias_parts.append("C:-5.0")

    # Add alanine penalty for all metal sites (favor coordinating residues)
    bias_parts.append("A:-2.0")

    return ",".join(bias_parts)


def get_coordination_number_range(
    metal: str,
    oxidation_state: int
) -> Tuple[int, int]:
    """
    Get expected coordination number range for a metal ion.

    Args:
        metal: Metal element symbol
        oxidation_state: Oxidation state

    Returns:
        Tuple of (min_CN, max_CN)
    """
    metal = metal.upper()

    if metal not in METAL_DATABASE:
        raise ValueError(f"Unknown metal: {metal}")

    data = METAL_DATABASE[metal]
    ox_state = oxidation_state if oxidation_state in data["coordination_numbers"] else data["default_oxidation"]

    if ox_state not in data["coordination_numbers"]:
        raise ValueError(f"Invalid oxidation state {oxidation_state} for {metal}")

    return data["coordination_numbers"][ox_state]


def validate_coordination_chemistry(
    metal: str,
    coordinating_residues: List[str],
    oxidation_state: int
) -> Dict[str, Any]:
    """
    Validate HSAB compatibility of coordinating residues with a metal ion.

    Args:
        metal: Metal element symbol
        coordinating_residues: List of residue names (3-letter or 1-letter codes)
        oxidation_state: Oxidation state

    Returns:
        Dict with keys:
            - valid: bool - overall validity
            - hsab_compatible: bool - HSAB theory compliance
            - warnings: List[str] - any warnings
            - coordination_number: int - actual CN
            - expected_cn_range: Tuple[int, int] - expected CN range
    """
    metal = metal.upper()
    warnings = []

    # Empty residue list
    if not coordinating_residues:
        return {
            "valid": False,
            "hsab_compatible": False,
            "warnings": ["No coordinating residues provided"],
            "coordination_number": 0,
            "expected_cn_range": (0, 0),
        }

    # Get metal properties
    try:
        hsab_class = get_hsab_class(metal, oxidation_state)
        cn_range = get_coordination_number_range(metal, oxidation_state)
        donors = get_preferred_donors(metal, oxidation_state)
    except ValueError as e:
        return {
            "valid": False,
            "hsab_compatible": False,
            "warnings": [str(e)],
            "coordination_number": len(coordinating_residues),
            "expected_cn_range": (0, 0),
        }

    # Normalize residue names to uppercase 3-letter codes
    normalized_residues = []
    for res in coordinating_residues:
        res = res.upper()
        if len(res) == 1:
            res = AA_1TO3.get(res, res)
        normalized_residues.append(res)

    # Check coordination number
    actual_cn = len(normalized_residues)
    min_cn, max_cn = cn_range

    if actual_cn < min_cn:
        warnings.append(f"Coordination number {actual_cn} is below minimum {min_cn}")
    elif actual_cn > max_cn:
        warnings.append(f"Coordination number {actual_cn} exceeds maximum {max_cn}")

    # Check HSAB compatibility
    hsab_compatible = True
    incompatible_residues = []

    for res in normalized_residues:
        donor = RESIDUE_DONORS.get(res)

        if donor is None:
            # Non-coordinating residue
            warnings.append(f"Residue {res} is not a typical coordinating residue")
            hsab_compatible = False
            incompatible_residues.append(res)
            continue

        # Check if this donor is compatible with metal's HSAB class
        donor_weight = donors.get(donor, 0)

        # Hard acid with soft donor (S from Cys) = incompatible
        if hsab_class == "hard" and donor == "S":
            hsab_compatible = False
            incompatible_residues.append(res)
            warnings.append(
                f"HSAB incompatible: {metal} ({hsab_class} acid) with {res} "
                f"(soft S donor). Cysteine coordination is not favorable."
            )

        # Hard acid with borderline donor (N from His)
        # For lanthanides, N is allowed with lower priority (not warned)
        # For other hard acids (Ca, Mg), N is less common but still acceptable
        # No warning needed - N donors are valid for hard acids just with lower preference

        # Soft acid with hard donor (O) = suboptimal
        if hsab_class == "soft" and donor == "O":
            warnings.append(
                f"Suboptimal: {metal} ({hsab_class} acid) with {res} "
                f"(hard O donor). Consider soft donors like Cys/Met."
            )

    # Determine overall validity
    valid = (
        hsab_compatible
        and min_cn <= actual_cn <= max_cn
        and len(incompatible_residues) == 0
    )

    # Allow valid=True with warnings for borderline cases
    if not hsab_compatible and len([w for w in warnings if "Suboptimal" in w]) == len(warnings):
        # Only suboptimal warnings, not hard incompatibilities
        valid = True
        hsab_compatible = True  # Technically compatible, just not ideal

    return {
        "valid": valid,
        "hsab_compatible": hsab_compatible,
        "warnings": warnings,
        "coordination_number": actual_cn,
        "expected_cn_range": cn_range,
        "incompatible_residues": incompatible_residues,
    }


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def is_lanthanide(metal: str) -> bool:
    """Check if a metal is a lanthanide."""
    return metal.upper() in ["TB", "EU", "GD", "LA", "CE", "SM", "YB", "DY"]


def get_hsab_class_simple(metal: str) -> str:
    """
    Get HSAB classification for a metal using default oxidation state.

    Simplified version for the design decision engine.

    Args:
        metal: Metal element symbol (e.g., "ZN", "FE", "TB")

    Returns:
        HSAB class: "hard", "borderline", or "soft"
    """
    metal = metal.upper()

    if metal not in METAL_DATABASE:
        # Fallback logic for unknown metals
        if metal in ["TB", "EU", "GD", "DY", "LA", "CE", "SM", "YB", "CA", "MG"]:
            return "hard"
        elif metal in ["CU", "AG", "AU", "HG"]:
            return "soft"
        return "borderline"

    entry = METAL_DATABASE[metal]
    default_ox = entry.get("default_oxidation", 2)
    hsab_data = entry.get("hsab_class", {})

    return hsab_data.get(default_ox, "borderline")


def get_common_residues(metal: str) -> List[str]:
    """Get list of commonly coordinating residues for a metal."""
    metal = metal.upper()
    if metal not in METAL_DATABASE:
        raise ValueError(f"Unknown metal: {metal}")
    return METAL_DATABASE[metal]["common_residues"].copy()


def get_all_metals() -> List[str]:
    """Get list of all supported metals."""
    return list(METAL_DATABASE.keys())


def get_metal_info(metal: str) -> Dict[str, Any]:
    """Get full information dictionary for a metal."""
    metal = metal.upper()
    if metal not in METAL_DATABASE:
        raise ValueError(f"Unknown metal: {metal}")
    return METAL_DATABASE[metal].copy()


# =============================================================================
# COORDINATION MOTIFS DATABASE
# =============================================================================

import math

COORDINATION_MOTIFS: Dict[str, Dict[str, Any]] = {
    "ef_hand": {
        "description": "EF-hand calcium binding motif",
        "metals": ["CA"],
        "coordination_number": 7,
        "geometry": "pentagonal_bipyramidal",
        "example_pdbs": ["1CLL", "3CLN"],
        "typical_residues": ["GLU", "ASP", "ASN", "GLN"],
        "notes": "Helix-loop-helix motif, 12-residue loop with conserved positions",
    },
    "zinc_finger": {
        "description": "Zinc finger structural motif",
        "metals": ["ZN"],
        "coordination_number": 4,
        "geometry": "tetrahedral",
        "example_pdbs": ["1ZNF", "1AAY"],
        "typical_residues": ["CYS", "HIS"],
        "notes": "C2H2, C4, or C3H variants common in DNA-binding proteins",
    },
    "his_tag": {
        "description": "Histidine tag binding site (IMAC)",
        "metals": ["NI", "CO", "ZN", "CU"],
        "coordination_number": 6,
        "geometry": "octahedral",
        "example_pdbs": [],
        "typical_residues": ["HIS"],
        "notes": "Polyhistidine sequence coordination for affinity purification",
    },
    "lanthanide_site": {
        "description": "Lanthanide binding site (engineered or natural)",
        "metals": ["TB", "EU", "GD", "LA"],
        "coordination_number": [8, 9],
        "geometry": ["square_antiprism", "tricapped_trigonal_prism"],
        "example_pdbs": [],
        "typical_residues": ["GLU", "ASP", "ASN", "GLN"],
        "notes": "High CN due to large ionic radius, O donors strongly preferred",
    },
    "rubredoxin": {
        "description": "Rubredoxin iron-sulfur center",
        "metals": ["FE"],
        "coordination_number": 4,
        "geometry": "tetrahedral",
        "example_pdbs": ["1IRO", "5RXN"],
        "typical_residues": ["CYS"],
        "notes": "FeS4 center with 4 Cys ligands, electron transfer protein",
    },
    "carbonic_anhydrase": {
        "description": "Carbonic anhydrase zinc active site",
        "metals": ["ZN"],
        "coordination_number": 4,
        "geometry": "tetrahedral",
        "example_pdbs": ["1CA2"],
        "typical_residues": ["HIS"],
        "notes": "Three His + water/hydroxide, catalyzes CO2 hydration",
    },
}


# =============================================================================
# IDEAL GEOMETRIES DATABASE
# =============================================================================

def _generate_sap_positions() -> List[Tuple[float, float, float]]:
    """
    Generate square antiprism (SAP) positions on unit sphere.

    CN=8: Two square faces rotated 45 degrees relative to each other.
    """
    positions = []
    # Top square (z positive)
    z_top = math.sin(math.radians(45))
    r_top = math.cos(math.radians(45))
    for i in range(4):
        angle = math.radians(i * 90)
        x = r_top * math.cos(angle)
        y = r_top * math.sin(angle)
        positions.append((x, y, z_top))

    # Bottom square (z negative), rotated 45 degrees
    z_bot = -z_top
    for i in range(4):
        angle = math.radians(i * 90 + 45)  # 45 degree rotation
        x = r_top * math.cos(angle)
        y = r_top * math.sin(angle)
        positions.append((x, y, z_bot))

    # Normalize to unit sphere
    normalized = []
    for x, y, z in positions:
        length = math.sqrt(x*x + y*y + z*z)
        normalized.append((x/length, y/length, z/length))

    return normalized


def _generate_ttp_positions() -> List[Tuple[float, float, float]]:
    """
    Generate tricapped trigonal prism (TTP) positions on unit sphere.

    CN=9: Trigonal prism (6 vertices) with 3 caps on the rectangular faces.
    """
    positions = []

    # Trigonal prism: two triangular faces
    z_top = 0.5
    z_bot = -0.5
    r_prism = math.sqrt(1 - z_top*z_top)  # Unit sphere radius at z_top

    # Top triangle
    for i in range(3):
        angle = math.radians(i * 120)
        x = r_prism * math.cos(angle)
        y = r_prism * math.sin(angle)
        positions.append((x, y, z_top))

    # Bottom triangle (aligned with top)
    for i in range(3):
        angle = math.radians(i * 120)
        x = r_prism * math.cos(angle)
        y = r_prism * math.sin(angle)
        positions.append((x, y, z_bot))

    # Three caps on rectangular faces (at z=0, midway between prism vertices)
    r_cap = 1.0
    for i in range(3):
        angle = math.radians(i * 120 + 60)  # Offset 60 degrees from prism vertices
        x = r_cap * math.cos(angle)
        y = r_cap * math.sin(angle)
        positions.append((x, y, 0.0))

    # Normalize to unit sphere
    normalized = []
    for x, y, z in positions:
        length = math.sqrt(x*x + y*y + z*z)
        if length > 0:
            normalized.append((x/length, y/length, z/length))
        else:
            normalized.append((0.0, 0.0, 1.0))  # Fallback for origin

    return normalized


IDEAL_GEOMETRIES: Dict[str, Dict[str, Any]] = {
    "linear": {
        "coordination_number": 2,
        "positions": [
            (0.0, 0.0, 1.0),
            (0.0, 0.0, -1.0),
        ],
        "description": "Linear geometry, 180 degree angle",
    },
    "trigonal_planar": {
        "coordination_number": 3,
        "positions": [
            (1.0, 0.0, 0.0),
            (-0.5, math.sqrt(3)/2, 0.0),
            (-0.5, -math.sqrt(3)/2, 0.0),
        ],
        "description": "Trigonal planar, 120 degree angles in xy plane",
    },
    "tetrahedral": {
        "coordination_number": 4,
        "positions": [
            (1.0, 1.0, 1.0),
            (1.0, -1.0, -1.0),
            (-1.0, 1.0, -1.0),
            (-1.0, -1.0, 1.0),
        ],
        "description": "Tetrahedral geometry, 109.5 degree angles",
    },
    "square_planar": {
        "coordination_number": 4,
        "positions": [
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (-1.0, 0.0, 0.0),
            (0.0, -1.0, 0.0),
        ],
        "description": "Square planar, 90 degree angles in xy plane",
    },
    "trigonal_bipyramidal": {
        "coordination_number": 5,
        "positions": [
            (0.0, 0.0, 1.0),   # Axial top
            (0.0, 0.0, -1.0),  # Axial bottom
            (1.0, 0.0, 0.0),   # Equatorial
            (-0.5, math.sqrt(3)/2, 0.0),
            (-0.5, -math.sqrt(3)/2, 0.0),
        ],
        "description": "Trigonal bipyramidal, axial 180, equatorial 120",
    },
    "square_pyramidal": {
        "coordination_number": 5,
        "positions": [
            (0.0, 0.0, 1.0),   # Apex
            (1.0, 0.0, 0.0),   # Base
            (0.0, 1.0, 0.0),
            (-1.0, 0.0, 0.0),
            (0.0, -1.0, 0.0),
        ],
        "description": "Square pyramidal, apex + square base",
    },
    "octahedral": {
        "coordination_number": 6,
        "positions": [
            (1.0, 0.0, 0.0),
            (-1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, -1.0, 0.0),
            (0.0, 0.0, 1.0),
            (0.0, 0.0, -1.0),
        ],
        "description": "Octahedral, 90 degree angles",
    },
    "pentagonal_bipyramidal": {
        "coordination_number": 7,
        "positions": [
            (0.0, 0.0, 1.0),   # Axial top
            (0.0, 0.0, -1.0),  # Axial bottom
            (1.0, 0.0, 0.0),   # Equatorial pentagon
            (math.cos(math.radians(72)), math.sin(math.radians(72)), 0.0),
            (math.cos(math.radians(144)), math.sin(math.radians(144)), 0.0),
            (math.cos(math.radians(216)), math.sin(math.radians(216)), 0.0),
            (math.cos(math.radians(288)), math.sin(math.radians(288)), 0.0),
        ],
        "description": "Pentagonal bipyramidal, common for Ca2+ in EF-hands",
    },
    "square_antiprism": {
        "coordination_number": 8,
        "positions": _generate_sap_positions(),
        "description": "Square antiprism, two squares rotated 45 degrees",
    },
    "tricapped_trigonal_prism": {
        "coordination_number": 9,
        "positions": _generate_ttp_positions(),
        "description": "Tricapped trigonal prism, common for lanthanides",
    },
}

# Normalize tetrahedral positions to unit sphere
_tet_norm = math.sqrt(3)
IDEAL_GEOMETRIES["tetrahedral"]["positions"] = [
    (x/_tet_norm, y/_tet_norm, z/_tet_norm)
    for x, y, z in IDEAL_GEOMETRIES["tetrahedral"]["positions"]
]


# =============================================================================
# COORDINATION TEMPLATE FUNCTIONS
# =============================================================================

def get_coordination_template(
    metal: str,
    coordination_number: Optional[int] = None
) -> Dict[str, Any]:
    """
    Get coordination template for a metal with optional coordination number.

    Returns information about expected geometry, typical residues, and
    example structures for the given metal and coordination number.

    Args:
        metal: Metal element symbol (e.g., "ZN", "CA", "TB")
        coordination_number: Optional specific coordination number. If not provided,
                           returns template for the metal's default coordination.

    Returns:
        Dict with keys:
            - metal: The metal symbol
            - coordination_number: The CN used
            - geometries: List of possible geometry names
            - typical_residues: List of residue codes
            - example_pdbs: List of example PDB IDs
            - motifs: List of matching coordination motifs
    """
    metal = metal.upper()

    if metal not in METAL_DATABASE:
        raise ValueError(f"Unknown metal: {metal}. Valid metals: {list(METAL_DATABASE.keys())}")

    metal_data = METAL_DATABASE[metal]
    default_ox = metal_data["default_oxidation"]
    cn_range = metal_data["coordination_numbers"].get(default_ox, (4, 6))

    # Determine the coordination number to use
    if coordination_number is None:
        # Use middle of range
        cn = (cn_range[0] + cn_range[1]) // 2
    else:
        cn = coordination_number

    # Find matching motifs for this metal and CN
    matching_motifs = []
    geometries = []
    example_pdbs = []

    for motif_name, motif_data in COORDINATION_MOTIFS.items():
        if metal in motif_data["metals"]:
            motif_cn = motif_data["coordination_number"]
            # Handle list of CNs (e.g., lanthanide_site)
            if isinstance(motif_cn, list):
                if cn in motif_cn:
                    matching_motifs.append(motif_name)
                    geom = motif_data["geometry"]
                    if isinstance(geom, list):
                        geometries.extend(geom)
                    else:
                        geometries.append(geom)
                    example_pdbs.extend(motif_data.get("example_pdbs", []))
            elif motif_cn == cn:
                matching_motifs.append(motif_name)
                geom = motif_data["geometry"]
                if isinstance(geom, list):
                    geometries.extend(geom)
                else:
                    geometries.append(geom)
                example_pdbs.extend(motif_data.get("example_pdbs", []))

    # If no motifs found, infer geometry from CN
    if not geometries:
        cn_to_geom = {
            2: ["linear"],
            3: ["trigonal_planar"],
            4: ["tetrahedral", "square_planar"],
            5: ["trigonal_bipyramidal", "square_pyramidal"],
            6: ["octahedral"],
            7: ["pentagonal_bipyramidal"],
            8: ["square_antiprism"],
            9: ["tricapped_trigonal_prism"],
        }
        geometries = cn_to_geom.get(cn, ["unknown"])

    # Remove duplicates while preserving order
    seen = set()
    unique_geometries = []
    for g in geometries:
        if g not in seen:
            seen.add(g)
            unique_geometries.append(g)

    return {
        "metal": metal,
        "coordination_number": cn,
        "geometries": unique_geometries,
        "typical_residues": metal_data.get("common_residues", []),
        "example_pdbs": list(set(example_pdbs)),
        "motifs": matching_motifs,
    }


def get_ideal_geometry(geometry_name: str) -> List[Tuple[float, float, float]]:
    """
    Get ideal ligand positions for a coordination geometry on unit sphere.

    Args:
        geometry_name: Name of geometry (e.g., "octahedral", "tetrahedral")

    Returns:
        List of (x, y, z) tuples representing ligand positions on unit sphere.
        Metal is assumed to be at origin (0, 0, 0).

    Raises:
        ValueError: If geometry_name is unknown.
    """
    geometry_name = geometry_name.lower()

    if geometry_name not in IDEAL_GEOMETRIES:
        valid = list(IDEAL_GEOMETRIES.keys())
        raise ValueError(f"Unknown geometry: {geometry_name}. Valid geometries: {valid}")

    return IDEAL_GEOMETRIES[geometry_name]["positions"].copy()


def validate_coordination(
    positions: List[Tuple[float, float, float]],
    metal_pos: Tuple[float, float, float],
    metal: str,
    expected_cn: Optional[int] = None
) -> Dict[str, Any]:
    """
    Validate coordination geometry around a metal center.

    Checks that ligand positions are within expected distance ranges for the
    metal and flags any geometric issues.

    Args:
        positions: List of (x, y, z) ligand atom positions in Angstroms
        metal_pos: (x, y, z) position of metal center in Angstroms
        metal: Metal element symbol
        expected_cn: Expected coordination number (optional)

    Returns:
        Dict with keys:
            - valid: bool - whether geometry is valid
            - coordination_number: int - actual number of coordinating atoms
            - issues: List[str] - list of any detected issues
            - distances: List[float] - distances to each ligand
            - geometry_match: str or None - best matching ideal geometry
    """
    metal = metal.upper()
    issues = []

    if not positions:
        return {
            "valid": False,
            "coordination_number": 0,
            "issues": ["No coordinating positions provided"],
            "distances": [],
            "geometry_match": None,
        }

    if metal not in METAL_DATABASE:
        return {
            "valid": False,
            "coordination_number": len(positions),
            "issues": [f"Unknown metal: {metal}"],
            "distances": [],
            "geometry_match": None,
        }

    metal_data = METAL_DATABASE[metal]
    default_ox = metal_data["default_oxidation"]

    # Calculate distances from metal to each ligand
    distances = []
    mx, my, mz = metal_pos
    for lx, ly, lz in positions:
        dist = math.sqrt((lx - mx)**2 + (ly - my)**2 + (lz - mz)**2)
        distances.append(dist)

    actual_cn = len(positions)

    # Check coordination number
    cn_range = metal_data["coordination_numbers"].get(default_ox, (4, 6))
    if expected_cn is not None:
        if actual_cn != expected_cn:
            issues.append(f"coordination_number: expected {expected_cn}, found {actual_cn}")

    if actual_cn < cn_range[0]:
        issues.append(f"coordination_number: {actual_cn} below minimum {cn_range[0]}")
    elif actual_cn > cn_range[1]:
        issues.append(f"coordination_number: {actual_cn} above maximum {cn_range[1]}")

    if actual_cn > 9:
        issues.append(f"Coordination number {actual_cn} exceeds common ranges (typical: 2-9)")

    # Check distances against expected ranges
    # Use O donor range as default (most common)
    bond_distances = metal_data.get("bond_distances", {}).get(default_ox, {})

    # Get the broadest acceptable range from all donor types
    min_dist = float('inf')
    max_dist = 0.0
    for donor, (d_min, d_max) in bond_distances.items():
        min_dist = min(min_dist, d_min)
        max_dist = max(max_dist, d_max)

    # Default range if no specific data
    if min_dist == float('inf'):
        logger.debug(f"Using default distance range for metal {metal}")
        min_dist = 1.8
        max_dist = 3.0

    # Add some tolerance
    min_dist_tol = min_dist - 0.2
    max_dist_tol = max_dist + 0.3

    for i, dist in enumerate(distances):
        if dist < min_dist_tol:
            issues.append(f"distance: ligand {i} at {dist:.2f}A is too close (min {min_dist_tol:.2f}A)")
        elif dist > max_dist_tol:
            issues.append(f"distance: ligand {i} at {dist:.2f}A is too far (max {max_dist_tol:.2f}A)")

    # Try to match to ideal geometry
    geometry_match = None
    if actual_cn in [2, 3, 4, 5, 6, 7, 8, 9]:
        cn_to_geom = {
            2: "linear",
            3: "trigonal_planar",
            4: "tetrahedral",  # Default; could also be square_planar
            5: "trigonal_bipyramidal",
            6: "octahedral",
            7: "pentagonal_bipyramidal",
            8: "square_antiprism",
            9: "tricapped_trigonal_prism",
        }
        geometry_match = cn_to_geom.get(actual_cn)

    valid = len(issues) == 0

    return {
        "valid": valid,
        "coordination_number": actual_cn,
        "issues": issues,
        "distances": distances,
        "geometry_match": geometry_match,
    }
