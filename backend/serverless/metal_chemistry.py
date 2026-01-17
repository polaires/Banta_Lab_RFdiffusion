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

from typing import Dict, List, Tuple, Any, Optional


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
    # LANTHANIDES - HARD ACIDS (ONLY O DONORS)
    # -------------------------------------------------------------------------
    "TB": {
        "name": "Terbium",
        "hsab_class": {3: "hard"},
        "default_oxidation": 3,
        "preferred_donors": {
            3: {
                "catalytic": {"O": 3.0, "N": -2.0, "S": -5.0},
                "structural": {"O": 3.0, "N": -2.0, "S": -5.0},
            }
        },
        "bond_distances": {
            3: {"O": (2.30, 2.60)}
        },
        "coordination_numbers": {3: (8, 9)},
        "common_residues": ["GLU", "ASP", "ASN", "GLN"],
        "description": "Hard acid lanthanide, ONLY O donors, exclude Cys",
    },

    "EU": {
        "name": "Europium",
        "hsab_class": {3: "hard"},
        "default_oxidation": 3,
        "preferred_donors": {
            3: {
                "catalytic": {"O": 3.0, "N": -2.0, "S": -5.0},
                "structural": {"O": 3.0, "N": -2.0, "S": -5.0},
            }
        },
        "bond_distances": {
            3: {"O": (2.35, 2.65)}
        },
        "coordination_numbers": {3: (8, 9)},
        "common_residues": ["GLU", "ASP", "ASN", "GLN"],
        "description": "Hard acid lanthanide, luminescent, ONLY O donors",
    },

    "GD": {
        "name": "Gadolinium",
        "hsab_class": {3: "hard"},
        "default_oxidation": 3,
        "preferred_donors": {
            3: {
                "catalytic": {"O": 3.0, "N": -2.0, "S": -5.0},
                "structural": {"O": 3.0, "N": -2.0, "S": -5.0},
            }
        },
        "bond_distances": {
            3: {"O": (2.30, 2.60)}
        },
        "coordination_numbers": {3: (8, 9)},
        "common_residues": ["GLU", "ASP", "ASN", "GLN"],
        "description": "Hard acid lanthanide, MRI contrast agent, ONLY O donors",
    },

    "LA": {
        "name": "Lanthanum",
        "hsab_class": {3: "hard"},
        "default_oxidation": 3,
        "preferred_donors": {
            3: {
                "catalytic": {"O": 3.0, "N": -2.0, "S": -5.0},
                "structural": {"O": 3.0, "N": -2.0, "S": -5.0},
            }
        },
        "bond_distances": {
            3: {"O": (2.40, 2.70)}
        },
        "coordination_numbers": {3: (8, 10)},
        "common_residues": ["GLU", "ASP", "ASN", "GLN"],
        "description": "Hard acid lanthanide, largest ionic radius, ONLY O donors",
    },

    "CE": {
        "name": "Cerium",
        "hsab_class": {3: "hard"},
        "default_oxidation": 3,
        "preferred_donors": {
            3: {
                "catalytic": {"O": 3.0, "N": -2.0, "S": -5.0},
                "structural": {"O": 3.0, "N": -2.0, "S": -5.0},
            }
        },
        "bond_distances": {
            3: {"O": (2.35, 2.65)}
        },
        "coordination_numbers": {3: (8, 10)},
        "common_residues": ["GLU", "ASP", "ASN", "GLN"],
        "description": "Hard acid lanthanide, ONLY O donors",
    },

    "SM": {
        "name": "Samarium",
        "hsab_class": {3: "hard"},
        "default_oxidation": 3,
        "preferred_donors": {
            3: {
                "catalytic": {"O": 3.0, "N": -2.0, "S": -5.0},
                "structural": {"O": 3.0, "N": -2.0, "S": -5.0},
            }
        },
        "bond_distances": {
            3: {"O": (2.30, 2.60)}
        },
        "coordination_numbers": {3: (8, 9)},
        "common_residues": ["GLU", "ASP", "ASN", "GLN"],
        "description": "Hard acid lanthanide, ONLY O donors",
    },

    "YB": {
        "name": "Ytterbium",
        "hsab_class": {3: "hard"},
        "default_oxidation": 3,
        "preferred_donors": {
            3: {
                "catalytic": {"O": 3.0, "N": -2.0, "S": -5.0},
                "structural": {"O": 3.0, "N": -2.0, "S": -5.0},
            }
        },
        "bond_distances": {
            3: {"O": (2.25, 2.55)}
        },
        "coordination_numbers": {3: (8, 9)},
        "common_residues": ["GLU", "ASP", "ASN", "GLN"],
        "description": "Hard acid lanthanide, smallest lanthanide, ONLY O donors",
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

        # Hard acid with borderline donor (N from His) = suboptimal but tolerated
        if hsab_class == "hard" and donor == "N":
            # For lanthanides specifically, N is penalized
            if metal in ["TB", "EU", "GD", "LA", "CE", "SM", "YB"]:
                warnings.append(
                    f"Suboptimal: {metal} (lanthanide) prefers O donors over N donors ({res}). "
                    f"Consider using Glu/Asp instead of His."
                )

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
    return metal.upper() in ["TB", "EU", "GD", "LA", "CE", "SM", "YB"]


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
