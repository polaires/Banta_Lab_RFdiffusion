# ligand_donors.py
"""
Ligand Donor Atom Identification Module

Parses ligand SMILES to identify coordinating atoms and score metal-ligand
compatibility using HSAB (Hard-Soft Acid-Base) theory.

Key features:
- Identify coordinating atoms (O, N, S, P) from SMILES
- Classify donor types (carboxylate, amine, thiol, etc.)
- HSAB compatibility scoring with metals
- Denticity analysis for chelation
- Coordination mode recommendations
"""

import logging
from typing import Dict, List, Any, Optional
from dataclasses import dataclass

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    logger.warning("RDKit not available - ligand analysis limited")

from metal_chemistry import get_hsab_class


# =============================================================================
# DONOR TYPE DEFINITIONS
# =============================================================================

DONOR_TYPES: Dict[str, Dict[str, Any]] = {
    "carboxylate": {"element": "O", "hsab": "hard", "denticity": 2},
    "carbonyl": {"element": "O", "hsab": "hard", "denticity": 1},
    "hydroxyl": {"element": "O", "hsab": "hard", "denticity": 1},
    "ether": {"element": "O", "hsab": "borderline", "denticity": 1},
    "amine": {"element": "N", "hsab": "borderline", "denticity": 1},
    "imine": {"element": "N", "hsab": "borderline", "denticity": 1},
    "pyridine": {"element": "N", "hsab": "borderline", "denticity": 1},
    "amide": {"element": "N", "hsab": "hard", "denticity": 1},
    "imidazole": {"element": "N", "hsab": "borderline", "denticity": 1},
    "thiol": {"element": "S", "hsab": "soft", "denticity": 1},
    "thioether": {"element": "S", "hsab": "soft", "denticity": 1},
    "thiolate": {"element": "S", "hsab": "soft", "denticity": 1},
    "phosphate": {"element": "P", "hsab": "hard", "denticity": 2},
}

# HSAB compatibility matrix
# Scores: matching = 1.0, borderline involved = 0.7, mismatch = 0.1
HSAB_COMPATIBILITY: Dict[str, Dict[str, float]] = {
    "hard": {"hard": 1.0, "borderline": 0.7, "soft": 0.1},
    "borderline": {"hard": 0.7, "borderline": 1.0, "soft": 0.7},
    "soft": {"hard": 0.1, "borderline": 0.7, "soft": 1.0},
}

# SMARTS patterns for donor classification
DONOR_SMARTS: Dict[str, str] = {
    # Carboxylate: -C(=O)[O-] or -C(=O)O (acid form)
    "carboxylate_anion": "[OX1-][CX3]=[OX1]",
    "carboxylic_acid": "[OX2H1][CX3]=[OX1]",
    "carboxylate_generic": "[OX1,OX2H1][CX3]=[OX1]",
    # Carbonyl (not carboxylate)
    "carbonyl": "[CX3]=[OX1]",
    # Hydroxyl
    "hydroxyl": "[OX2H1][#6]",
    # Ether
    "ether": "[OX2]([#6])[#6]",
    # Amine (primary, secondary, tertiary)
    "amine_primary": "[NX3H2][#6]",
    "amine_secondary": "[NX3H1]([#6])[#6]",
    "amine_tertiary": "[NX3]([#6])([#6])[#6]",
    # Imine
    "imine": "[NX2]=[CX3]",
    # Amide N (N of amide)
    "amide_n": "[NX3][CX3]=[OX1]",
    # Pyridine
    "pyridine": "[nX2]",
    # Imidazole
    "imidazole": "[nR1]1[cR1][nR1][cR1][cR1]1",
    # Thiol
    "thiol": "[SX2H1]",
    # Thioether
    "thioether": "[SX2]([#6])[#6]",
    # Phosphate
    "phosphate": "[PX4](=[OX1])([OX2,OX1-])([OX2,OX1-])[OX2,OX1-]",
}


# =============================================================================
# PUBLIC FUNCTIONS
# =============================================================================

def identify_donors_from_smiles(smiles: str) -> List[Dict[str, Any]]:
    """
    Parse SMILES and identify potential coordinating atoms.

    Args:
        smiles: SMILES string of the ligand

    Returns:
        List of donor dictionaries with keys:
            - atom_idx: int - Index of atom in molecule
            - element: str - Element symbol (O, N, S, P)
            - type: str - Donor type classification
            - hsab: str - HSAB class of the donor
            - coords: Optional tuple - 3D coordinates if available
    """
    if not HAS_RDKIT:
        logger.warning("RDKit not available - returning empty donor list")
        return []

    if not smiles or not isinstance(smiles, str):
        return []

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.warning(f"Failed to parse SMILES: {smiles}")
        return []

    # Add hydrogens for better classification
    mol = Chem.AddHs(mol)

    donors = []
    processed_atoms = set()

    # Iterate through potential donor atoms (O, N, S, P)
    for atom in mol.GetAtoms():
        atom_idx = atom.GetIdx()
        element = atom.GetSymbol()

        if element not in ["O", "N", "S", "P"]:
            continue

        if atom_idx in processed_atoms:
            continue

        # Classify the donor type
        donor_type, hsab = _classify_donor_type(atom, mol)

        if donor_type is not None:
            donors.append({
                "atom_idx": atom_idx,
                "element": element,
                "type": donor_type,
                "hsab": hsab,
                "coords": None,  # Would need 3D conformer generation
            })
            processed_atoms.add(atom_idx)

    return donors


def _classify_donor_type(atom, mol) -> tuple:
    """
    Classify the donor type of an atom based on its chemical environment.

    Args:
        atom: RDKit atom object
        mol: RDKit molecule object

    Returns:
        Tuple of (donor_type, hsab_class) or (None, None) if not a donor
    """
    element = atom.GetSymbol()
    atom_idx = atom.GetIdx()

    # Get neighbors
    neighbors = atom.GetNeighbors()
    neighbor_symbols = [n.GetSymbol() for n in neighbors]

    # =============================================================================
    # OXYGEN CLASSIFICATION
    # =============================================================================
    if element == "O":
        # Check for carboxylate pattern
        # O connected to C which is double bonded to another O
        for neighbor in neighbors:
            if neighbor.GetSymbol() == "C":
                c_neighbors = neighbor.GetNeighbors()
                # Look for C=O pattern
                for c_neighbor in c_neighbors:
                    if c_neighbor.GetSymbol() == "O" and c_neighbor.GetIdx() != atom_idx:
                        bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), c_neighbor.GetIdx())
                        if bond and bond.GetBondTypeAsDouble() == 2.0:
                            # This is part of a carboxylate group
                            return ("carboxylate", "hard")
                        # Also check if we ARE the double-bonded oxygen
                bond_to_c = mol.GetBondBetweenAtoms(atom_idx, neighbor.GetIdx())
                if bond_to_c and bond_to_c.GetBondTypeAsDouble() == 2.0:
                    # This is C=O, check if connected to another O
                    for c_neighbor in c_neighbors:
                        if c_neighbor.GetSymbol() == "O" and c_neighbor.GetIdx() != atom_idx:
                            return ("carboxylate", "hard")
                    # Plain carbonyl
                    return ("carbonyl", "hard")

        # Check for hydroxyl (O-H)
        total_h = atom.GetTotalNumHs()
        if total_h > 0 and len(neighbors) == 1:
            # Check if attached to carbon (alcohol) vs other
            if neighbor_symbols[0] == "C":
                return ("hydroxyl", "hard")

        # Ether: O bonded to two carbons
        if len(neighbors) == 2 and all(n.GetSymbol() == "C" for n in neighbors):
            return ("ether", "borderline")

        # Default oxygen as hard donor
        return ("hydroxyl", "hard")

    # =============================================================================
    # NITROGEN CLASSIFICATION
    # =============================================================================
    if element == "N":
        # Check for amide N first (N bonded to C=O)
        for neighbor in neighbors:
            if neighbor.GetSymbol() == "C":
                c_neighbors = neighbor.GetNeighbors()
                for c_neighbor in c_neighbors:
                    if c_neighbor.GetSymbol() == "O":
                        bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), c_neighbor.GetIdx())
                        if bond and bond.GetBondTypeAsDouble() == 2.0:
                            return ("amide", "hard")

        # Check for imine (N=C)
        for neighbor in neighbors:
            if neighbor.GetSymbol() == "C":
                bond = mol.GetBondBetweenAtoms(atom_idx, neighbor.GetIdx())
                if bond and bond.GetBondTypeAsDouble() == 2.0:
                    return ("imine", "borderline")

        # Check if in aromatic ring (pyridine-like)
        if atom.GetIsAromatic():
            return ("pyridine", "borderline")

        # Regular amine
        total_h = atom.GetTotalNumHs()
        if total_h == 2:
            return ("amine", "borderline")  # Primary
        elif total_h == 1:
            return ("amine", "borderline")  # Secondary
        else:
            return ("amine", "borderline")  # Tertiary

    # =============================================================================
    # SULFUR CLASSIFICATION
    # =============================================================================
    if element == "S":
        # Check for thiol (S-H)
        total_h = atom.GetTotalNumHs()
        if total_h > 0:
            return ("thiol", "soft")

        # Thioether: S bonded to two carbons
        if len(neighbors) == 2 and all(n.GetSymbol() == "C" for n in neighbors):
            return ("thioether", "soft")

        # Default sulfur as soft donor
        return ("thiol", "soft")

    # =============================================================================
    # PHOSPHORUS CLASSIFICATION
    # =============================================================================
    if element == "P":
        # Check for phosphate (P bonded to multiple O)
        o_count = sum(1 for n in neighbors if n.GetSymbol() == "O")
        if o_count >= 3:
            return ("phosphate", "hard")

        # Default phosphorus
        return ("phosphate", "hard")

    return (None, None)


def score_ligand_metal_compatibility(
    smiles: str,
    metal: str,
    oxidation_state: Optional[int] = None
) -> float:
    """
    Score ligand-metal compatibility using HSAB theory.

    Args:
        smiles: SMILES string of the ligand
        metal: Metal element symbol (e.g., "ZN", "TB", "FE")
        oxidation_state: Optional oxidation state. If None, uses default.

    Returns:
        Compatibility score from 0.0 (incompatible) to 1.0 (highly compatible)
    """
    if not HAS_RDKIT:
        logger.warning("RDKit not available - returning 0.0 score")
        return 0.0

    # Get donors
    donors = identify_donors_from_smiles(smiles)

    if not donors:
        return 0.0

    # Get metal HSAB class
    try:
        # Import metal database to get default oxidation state
        from metal_chemistry import METAL_DATABASE
        metal_upper = metal.upper()

        if oxidation_state is None:
            if metal_upper in METAL_DATABASE:
                oxidation_state = METAL_DATABASE[metal_upper]["default_oxidation"]
            else:
                oxidation_state = 2  # Default guess

        metal_hsab = get_hsab_class(metal, oxidation_state)
    except ValueError:
        logger.warning(f"Unknown metal {metal}, assuming borderline")
        metal_hsab = "borderline"

    # Calculate weighted average compatibility
    total_score = 0.0
    total_weight = 0.0

    for donor in donors:
        donor_hsab = donor["hsab"]
        compatibility = HSAB_COMPATIBILITY.get(metal_hsab, {}).get(donor_hsab, 0.5)

        # Weight by donor type (carboxylates are generally better ligands)
        weight = 1.0
        if donor["type"] == "carboxylate":
            weight = 1.5
        elif donor["type"] in ["thiol", "thiolate"]:
            weight = 1.5

        total_score += compatibility * weight
        total_weight += weight

    if total_weight == 0:
        return 0.0

    return total_score / total_weight


def analyze_denticity(smiles: str) -> List[Dict[str, Any]]:
    """
    Analyze potential denticity of coordinating groups in a ligand.

    Args:
        smiles: SMILES string of the ligand

    Returns:
        List of dictionaries describing chelating groups:
            - group_type: str - Type of chelating group
            - atom_indices: List[int] - Indices of donor atoms in the group
            - max_denticity: int - Maximum coordination points
            - elements: List[str] - Elements of donor atoms
    """
    if not HAS_RDKIT:
        logger.warning("RDKit not available - returning empty denticity analysis")
        return []

    if not smiles or not isinstance(smiles, str):
        return []

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.warning(f"Failed to parse SMILES: {smiles}")
        return []

    mol = Chem.AddHs(mol)

    denticity_groups = []

    # Find carboxylate groups (bidentate potential)
    carboxylate_pattern = Chem.MolFromSmarts("[OX1,OX2H1][CX3](=[OX1])")
    if carboxylate_pattern:
        matches = mol.GetSubstructMatches(carboxylate_pattern)

        # Group by the central carbon
        seen_carbons = set()
        for match in matches:
            # match[0] is the O, match[1] is the C, match[2] is the =O
            carbon_idx = match[1]
            if carbon_idx in seen_carbons:
                continue
            seen_carbons.add(carbon_idx)

            # Find both oxygens of this carboxylate
            carbon = mol.GetAtomWithIdx(carbon_idx)
            oxygen_indices = []
            for neighbor in carbon.GetNeighbors():
                if neighbor.GetSymbol() == "O":
                    oxygen_indices.append(neighbor.GetIdx())

            if len(oxygen_indices) >= 2:
                denticity_groups.append({
                    "group_type": "carboxylate",
                    "atom_indices": oxygen_indices,
                    "max_denticity": 2,
                    "elements": ["O", "O"],
                })

    # Find phosphate groups (multi-dentate potential)
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=[OX1])([OX2,OX1-])([OX2,OX1-])[OX2,OX1-]")
    if phosphate_pattern:
        matches = mol.GetSubstructMatches(phosphate_pattern)
        for match in matches:
            # Phosphate can be tridentate or more
            oxygen_indices = [idx for idx in match if mol.GetAtomWithIdx(idx).GetSymbol() == "O"]
            denticity_groups.append({
                "group_type": "phosphate",
                "atom_indices": oxygen_indices,
                "max_denticity": len(oxygen_indices),
                "elements": ["O"] * len(oxygen_indices),
            })

    # Find hydroxyl groups (single O but can pair with carboxylate for chelation)
    # Alpha-hydroxy acids are common chelators

    return denticity_groups


def get_recommended_coordination_mode(
    smiles: str,
    metal: str,
    target_cn: int = 6
) -> Dict[str, Any]:
    """
    Get recommended coordination mode for a ligand with a specific metal.

    Args:
        smiles: SMILES string of the ligand
        metal: Metal element symbol
        target_cn: Target coordination number

    Returns:
        Dictionary with recommendations:
            - recommended_donors: List of donors to use
            - coordination_mode: str - "monodentate", "bidentate", "polydentate"
            - estimated_cn_contribution: int - How many coordination sites this ligand fills
            - compatibility_score: float - HSAB-based score
    """
    if not HAS_RDKIT:
        return {
            "recommended_donors": [],
            "coordination_mode": "unknown",
            "estimated_cn_contribution": 0,
            "compatibility_score": 0.0,
            "error": "RDKit not available",
        }

    donors = identify_donors_from_smiles(smiles)
    denticity = analyze_denticity(smiles)
    compatibility = score_ligand_metal_compatibility(smiles, metal)

    # Determine coordination mode based on denticity groups
    if not denticity:
        coordination_mode = "monodentate"
        cn_contribution = min(len(donors), 1) if donors else 0
    else:
        max_dent = max(d["max_denticity"] for d in denticity)
        if max_dent >= 3:
            coordination_mode = "polydentate"
        elif max_dent == 2:
            coordination_mode = "bidentate"
        else:
            coordination_mode = "monodentate"

        # Estimate CN contribution
        cn_contribution = sum(d["max_denticity"] for d in denticity)

    # Filter donors by compatibility with metal
    try:
        from metal_chemistry import METAL_DATABASE
        metal_upper = metal.upper()
        if metal_upper in METAL_DATABASE:
            oxidation_state = METAL_DATABASE[metal_upper]["default_oxidation"]
            metal_hsab = get_hsab_class(metal, oxidation_state)
        else:
            metal_hsab = "borderline"
    except (ValueError, ImportError):
        metal_hsab = "borderline"

    # Rank donors by compatibility
    recommended_donors = []
    for donor in donors:
        donor_hsab = donor["hsab"]
        compat = HSAB_COMPATIBILITY.get(metal_hsab, {}).get(donor_hsab, 0.5)
        if compat >= 0.5:  # Only recommend compatible donors
            recommended_donors.append({
                **donor,
                "compatibility": compat,
            })

    # Sort by compatibility
    recommended_donors.sort(key=lambda x: x.get("compatibility", 0), reverse=True)

    return {
        "recommended_donors": recommended_donors,
        "coordination_mode": coordination_mode,
        "estimated_cn_contribution": cn_contribution,
        "compatibility_score": compatibility,
        "total_donors": len(donors),
        "denticity_groups": denticity,
    }


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def get_ligand_summary(smiles: str) -> Dict[str, Any]:
    """
    Get a summary of ligand donor properties.

    Args:
        smiles: SMILES string of the ligand

    Returns:
        Summary dictionary with donor counts and types
    """
    donors = identify_donors_from_smiles(smiles)
    denticity = analyze_denticity(smiles)

    # Count by element
    element_counts = {}
    for donor in donors:
        elem = donor["element"]
        element_counts[elem] = element_counts.get(elem, 0) + 1

    # Count by type
    type_counts = {}
    for donor in donors:
        dtype = donor["type"]
        type_counts[dtype] = type_counts.get(dtype, 0) + 1

    # Determine overall HSAB character
    hsab_counts = {"hard": 0, "borderline": 0, "soft": 0}
    for donor in donors:
        hsab = donor["hsab"]
        hsab_counts[hsab] = hsab_counts.get(hsab, 0) + 1

    if not donors:
        overall_hsab = "unknown"
    elif hsab_counts["hard"] >= hsab_counts["soft"] and hsab_counts["hard"] >= hsab_counts["borderline"]:
        overall_hsab = "hard"
    elif hsab_counts["soft"] >= hsab_counts["hard"] and hsab_counts["soft"] >= hsab_counts["borderline"]:
        overall_hsab = "soft"
    else:
        overall_hsab = "borderline"

    return {
        "smiles": smiles,
        "total_donors": len(donors),
        "element_counts": element_counts,
        "type_counts": type_counts,
        "hsab_character": overall_hsab,
        "hsab_counts": hsab_counts,
        "denticity_groups": len(denticity),
        "max_denticity": max((d["max_denticity"] for d in denticity), default=1),
    }
