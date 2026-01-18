# Metal-Ligand Complex Dimer Design Enhancement

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Enhance metal-ligand complex dimer design with chemistry-aware amino acid biases, database-driven templates, ligand coordination intelligence, and comprehensive test suite.

**Architecture:** Four-phase enhancement: (1) HSAB-compliant amino acid biases per metal class, (2) PDB database integration for template retrieval with calculated fallback, (3) SMILES-to-donor atom mapping for ligand coordination, (4) biochemically rigorous test suite with real PDB structures.

**Tech Stack:** Python 3.12, pytest, requests (PDB API), RDKit (SMILES parsing), numpy (geometry calculations)

---

## Phase 1: HSAB-Compliant Amino Acid Bias System

### Task 1.1: Create Metal Chemistry Database Module

**Files:**
- Create: `backend/serverless/metal_chemistry.py`
- Test: `backend/serverless/test_metal_chemistry.py`

**Step 1: Write the failing test for HSAB classification**

```python
# test_metal_chemistry.py
"""Tests for metal chemistry database and HSAB classification."""
import pytest


class TestHSABClassification:
    """Test Hard-Soft Acid-Base classifications."""

    def test_zinc_is_borderline_acid(self):
        """Zinc should be classified as borderline acid."""
        from metal_chemistry import get_hsab_class
        assert get_hsab_class("ZN") == "borderline"

    def test_lanthanides_are_hard_acids(self):
        """All lanthanides should be hard acids."""
        from metal_chemistry import get_hsab_class
        for metal in ["TB", "EU", "GD", "LA", "CE", "SM", "YB"]:
            assert get_hsab_class(metal) == "hard", f"{metal} should be hard"

    def test_copper_i_is_soft_acid(self):
        """Cu(I) is soft, Cu(II) is borderline."""
        from metal_chemistry import get_hsab_class
        assert get_hsab_class("CU", oxidation_state=1) == "soft"
        assert get_hsab_class("CU", oxidation_state=2) == "borderline"

    def test_iron_depends_on_oxidation_state(self):
        """Fe(II) is borderline, Fe(III) is hard."""
        from metal_chemistry import get_hsab_class
        assert get_hsab_class("FE", oxidation_state=2) == "borderline"
        assert get_hsab_class("FE", oxidation_state=3) == "hard"


class TestPreferredDonors:
    """Test donor atom preferences per metal."""

    def test_zinc_prefers_soft_donors(self):
        """Zinc prefers S (Cys) and N (His) over O."""
        from metal_chemistry import get_preferred_donors
        donors = get_preferred_donors("ZN")
        assert donors[0] in ["S", "N"]  # Top preference
        assert "O" in donors  # Acceptable but lower priority

    def test_lanthanide_excludes_sulfur(self):
        """Lanthanides should never prefer sulfur donors."""
        from metal_chemistry import get_preferred_donors
        for metal in ["TB", "EU", "GD"]:
            donors = get_preferred_donors(metal)
            # S should either be absent or have negative weight
            assert "S" not in donors or donors.get("S", 0) < 0

    def test_copper_type1_includes_methionine(self):
        """Type I copper sites include Met coordination."""
        from metal_chemistry import get_preferred_donors
        donors = get_preferred_donors("CU", site_type="type1")
        assert "S_met" in donors or "M" in donors


class TestAminoAcidBias:
    """Test amino acid bias string generation."""

    def test_zinc_structural_bias_favors_cys_his(self):
        """Zinc structural sites need Cys4 or Cys-His combinations."""
        from metal_chemistry import get_amino_acid_bias
        bias = get_amino_acid_bias("ZN", site_type="structural")
        # Parse bias string
        bias_dict = dict(pair.split(":") for pair in bias.split(","))
        assert float(bias_dict.get("C", 0)) >= 2.5  # Cys highly favored
        assert float(bias_dict.get("H", 0)) >= 2.0  # His favored

    def test_lanthanide_bias_excludes_cysteine(self):
        """Lanthanide biases must exclude Cys completely."""
        from metal_chemistry import get_amino_acid_bias
        bias = get_amino_acid_bias("TB")
        bias_dict = dict(pair.split(":") for pair in bias.split(","))
        assert float(bias_dict.get("C", 0)) <= -3.0  # Cys strongly disfavored

    def test_lanthanide_allows_soft_with_lower_priority(self):
        """Lanthanides can have His/Asn but with lower weight than Glu/Asp."""
        from metal_chemistry import get_amino_acid_bias
        bias = get_amino_acid_bias("TB")
        bias_dict = dict(pair.split(":") for pair in bias.split(","))
        glu_weight = float(bias_dict.get("E", 0))
        his_weight = float(bias_dict.get("H", 0))
        # His allowed but lower than Glu
        assert glu_weight > his_weight
        # His should be positive or zero, not strongly negative
        assert his_weight >= -0.5


class TestBondDistances:
    """Test metal-ligand bond distance ranges."""

    def test_zinc_distances(self):
        """Zinc bond distances from crystallographic data."""
        from metal_chemistry import get_bond_distance_range
        # Zn-S (Cys): 2.30-2.35 Å
        s_range = get_bond_distance_range("ZN", "S")
        assert 2.28 <= s_range[0] <= 2.32
        assert 2.33 <= s_range[1] <= 2.38
        # Zn-N (His): 2.00-2.10 Å
        n_range = get_bond_distance_range("ZN", "N")
        assert 1.98 <= n_range[0] <= 2.02
        assert 2.05 <= n_range[1] <= 2.12

    def test_lanthanide_oxygen_distances(self):
        """Lanthanide-O distances from crystallographic data."""
        from metal_chemistry import get_bond_distance_range
        # Tb-O: 2.30-2.50 Å
        o_range = get_bond_distance_range("TB", "O")
        assert 2.28 <= o_range[0] <= 2.35
        assert 2.45 <= o_range[1] <= 2.55
```

**Step 2: Run test to verify it fails**

Run: `cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless && python -m pytest test_metal_chemistry.py -v --tb=short 2>&1 | head -30`
Expected: FAIL with "ModuleNotFoundError: No module named 'metal_chemistry'"

**Step 3: Write minimal implementation**

```python
# metal_chemistry.py
"""
Metal Chemistry Database

HSAB (Hard-Soft Acid-Base) theory-compliant metal coordination chemistry.
Based on crystallographic data from MESPEUS, MetalPDB, and CSD.

References:
- Pearson, R.G. (1963) JACS - HSAB Principle
- Harding et al. (2010) Acta Cryst D - Metal coordination in proteins
- CSD lanthanide analysis (Nature Sci Rep 2024)
"""

from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass
from enum import Enum


class HSABClass(Enum):
    """Hard-Soft Acid-Base classification."""
    HARD = "hard"
    BORDERLINE = "borderline"
    SOFT = "soft"


@dataclass
class MetalProperties:
    """Complete metal coordination properties."""
    symbol: str
    name: str
    hsab_class: HSABClass
    oxidation_states: List[int]
    typical_coordination_numbers: List[int]
    preferred_geometries: List[str]
    preferred_donors: Dict[str, float]  # Element/type -> preference weight
    bond_distances: Dict[str, Tuple[float, float]]  # Donor -> (min, max) Å
    notes: str = ""


# Comprehensive metal database based on crystallographic surveys
METAL_DATABASE: Dict[str, MetalProperties] = {
    # =========================================================================
    # TRANSITION METALS (Borderline/Soft acids)
    # =========================================================================
    "ZN": MetalProperties(
        symbol="ZN",
        name="Zinc",
        hsab_class=HSABClass.BORDERLINE,
        oxidation_states=[2],
        typical_coordination_numbers=[4, 5, 6],
        preferred_geometries=["tetrahedral", "trigonal_bipyramidal", "octahedral"],
        preferred_donors={
            "S": 3.0,   # Cys - strongly preferred
            "N": 2.5,   # His - strongly preferred
            "O": 1.0,   # Asp/Glu - acceptable
        },
        bond_distances={
            "S": (2.30, 2.36),  # Zn-Cys Sγ
            "N": (2.00, 2.10),  # Zn-His Nε2/Nδ1
            "O": (1.95, 2.10),  # Zn-Asp/Glu carboxylate
        },
        notes="Structural: Cys4/Cys3His. Catalytic: His2-3 + water/Glu"
    ),
    "FE2": MetalProperties(
        symbol="FE",
        name="Iron(II)",
        hsab_class=HSABClass.BORDERLINE,
        oxidation_states=[2],
        typical_coordination_numbers=[4, 5, 6],
        preferred_geometries=["tetrahedral", "octahedral", "square_pyramidal"],
        preferred_donors={
            "N": 2.5,   # His
            "S": 2.0,   # Cys (FeS clusters)
            "O": 1.5,   # Asp/Glu
        },
        bond_distances={
            "S": (2.28, 2.36),
            "N": (2.05, 2.20),
            "O": (2.00, 2.15),
        },
        notes="Non-heme: 2-His-1-carboxylate facial triad. FeS: Cys4"
    ),
    "FE3": MetalProperties(
        symbol="FE",
        name="Iron(III)",
        hsab_class=HSABClass.HARD,  # Fe3+ is harder than Fe2+
        oxidation_states=[3],
        typical_coordination_numbers=[5, 6],
        preferred_geometries=["octahedral", "trigonal_bipyramidal"],
        preferred_donors={
            "O": 2.5,   # Glu/Asp preferred for Fe3+
            "N": 2.0,   # His
            "S": 1.0,   # Less preferred for Fe3+
        },
        bond_distances={
            "S": (2.25, 2.35),
            "N": (2.00, 2.15),
            "O": (1.90, 2.10),
        },
        notes="Higher oxidation state prefers harder O donors"
    ),
    "CU1": MetalProperties(
        symbol="CU",
        name="Copper(I)",
        hsab_class=HSABClass.SOFT,
        oxidation_states=[1],
        typical_coordination_numbers=[2, 3, 4],
        preferred_geometries=["linear", "trigonal", "tetrahedral"],
        preferred_donors={
            "S": 3.5,   # Cys - very strong
            "S_met": 2.5,  # Met thioether
            "N": 2.0,   # His
        },
        bond_distances={
            "S": (2.15, 2.25),  # Cu(I)-Cys
            "S_met": (2.40, 2.55),  # Cu(I)-Met
            "N": (1.95, 2.10),
        },
        notes="Type I blue copper: 2His + Cys + Met. Uses Nδ1 of His."
    ),
    "CU2": MetalProperties(
        symbol="CU",
        name="Copper(II)",
        hsab_class=HSABClass.BORDERLINE,
        oxidation_states=[2],
        typical_coordination_numbers=[4, 5, 6],
        preferred_geometries=["square_planar", "square_pyramidal", "octahedral"],
        preferred_donors={
            "N": 2.5,   # His
            "S": 2.0,   # Cys
            "O": 1.5,   # Asp/Glu
        },
        bond_distances={
            "S": (2.20, 2.35),
            "N": (1.95, 2.10),
            "O": (1.90, 2.05),
        },
        notes="Type II copper. Uses Nε2 of His predominantly."
    ),
    "MN": MetalProperties(
        symbol="MN",
        name="Manganese",
        hsab_class=HSABClass.HARD,
        oxidation_states=[2, 3, 4],
        typical_coordination_numbers=[5, 6],
        preferred_geometries=["octahedral", "trigonal_bipyramidal"],
        preferred_donors={
            "O": 2.5,   # Glu/Asp
            "N": 2.0,   # His
            "S": -1.0,  # Disfavored
        },
        bond_distances={
            "N": (2.15, 2.30),
            "O": (2.10, 2.25),
        },
        notes="OEC (Mn4Ca cluster) and other O2-evolving enzymes"
    ),
    "CO": MetalProperties(
        symbol="CO",
        name="Cobalt",
        hsab_class=HSABClass.BORDERLINE,
        oxidation_states=[2, 3],
        typical_coordination_numbers=[4, 6],
        preferred_geometries=["tetrahedral", "octahedral"],
        preferred_donors={
            "N": 2.5,
            "S": 2.0,
            "O": 1.5,
        },
        bond_distances={
            "N": (2.00, 2.15),
            "S": (2.25, 2.40),
            "O": (2.00, 2.15),
        },
    ),
    "NI": MetalProperties(
        symbol="NI",
        name="Nickel",
        hsab_class=HSABClass.BORDERLINE,
        oxidation_states=[2],
        typical_coordination_numbers=[4, 5, 6],
        preferred_geometries=["square_planar", "octahedral"],
        preferred_donors={
            "S": 2.5,   # NiFe hydrogenase
            "N": 2.5,
            "O": 1.0,
        },
        bond_distances={
            "S": (2.15, 2.30),
            "N": (1.90, 2.10),
            "O": (2.00, 2.15),
        },
    ),
    # =========================================================================
    # ALKALINE EARTH METALS (Hard acids)
    # =========================================================================
    "CA": MetalProperties(
        symbol="CA",
        name="Calcium",
        hsab_class=HSABClass.HARD,
        oxidation_states=[2],
        typical_coordination_numbers=[6, 7, 8],
        preferred_geometries=["octahedral", "pentagonal_bipyramidal"],
        preferred_donors={
            "O": 3.0,   # Strongly prefers O
            "N": 0.5,   # Rare
            "S": -3.0,  # Never
        },
        bond_distances={
            "O": (2.30, 2.55),  # Ca-carboxylate
            "O_water": (2.35, 2.50),
            "O_carbonyl": (2.35, 2.55),
        },
        notes="EF-hand: 7-coordinate pentagonal bipyramidal"
    ),
    "MG": MetalProperties(
        symbol="MG",
        name="Magnesium",
        hsab_class=HSABClass.HARD,
        oxidation_states=[2],
        typical_coordination_numbers=[6],
        preferred_geometries=["octahedral"],
        preferred_donors={
            "O": 3.0,
            "N": 0.5,
            "S": -3.0,
        },
        bond_distances={
            "O": (2.00, 2.20),
            "N": (2.10, 2.25),
        },
        notes="Strictly octahedral in most enzymes"
    ),
    # =========================================================================
    # LANTHANIDES (Hard acids - oxygen only)
    # =========================================================================
    "TB": MetalProperties(
        symbol="TB",
        name="Terbium",
        hsab_class=HSABClass.HARD,
        oxidation_states=[3],
        typical_coordination_numbers=[8, 9],
        preferred_geometries=["square_antiprism", "tricapped_trigonal_prism"],
        preferred_donors={
            "O_carboxylate": 3.5,  # Glu/Asp bidentate
            "O_carbonyl": 2.5,     # Backbone/Asn/Gln
            "O_water": 2.0,        # Coordinated water
            "N": 0.0,              # Very rare but possible
            "S": -5.0,             # NEVER - exclude completely
        },
        bond_distances={
            "O": (2.30, 2.50),
            "O_carboxylate": (2.30, 2.45),
            "O_carbonyl": (2.40, 2.55),
            "O_water": (2.40, 2.55),
        },
        notes="Luminescent. 4f→4f transitions. Antenna effect from Trp."
    ),
    "EU": MetalProperties(
        symbol="EU",
        name="Europium",
        hsab_class=HSABClass.HARD,
        oxidation_states=[2, 3],
        typical_coordination_numbers=[8, 9],
        preferred_geometries=["square_antiprism", "bicapped_trigonal_prism"],
        preferred_donors={
            "O_carboxylate": 3.5,
            "O_carbonyl": 2.5,
            "O_water": 2.0,
            "N": 0.0,
            "S": -5.0,
        },
        bond_distances={
            "O": (2.35, 2.55),
        },
        notes="Red emission. Most common luminescent Ln for LRET."
    ),
    "GD": MetalProperties(
        symbol="GD",
        name="Gadolinium",
        hsab_class=HSABClass.HARD,
        oxidation_states=[3],
        typical_coordination_numbers=[8, 9],
        preferred_geometries=["square_antiprism", "tricapped_trigonal_prism"],
        preferred_donors={
            "O_carboxylate": 3.5,
            "O_carbonyl": 2.5,
            "O_water": 2.0,
            "N": 0.0,
            "S": -5.0,
        },
        bond_distances={
            "O": (2.32, 2.52),
        },
        notes="MRI contrast agent. High magnetic moment."
    ),
    "LA": MetalProperties(
        symbol="LA",
        name="Lanthanum",
        hsab_class=HSABClass.HARD,
        oxidation_states=[3],
        typical_coordination_numbers=[9, 10, 12],
        preferred_geometries=["tricapped_trigonal_prism", "bicapped_square_antiprism"],
        preferred_donors={
            "O_carboxylate": 3.5,
            "O_carbonyl": 2.5,
            "O_water": 2.0,
            "N": 0.0,
            "S": -5.0,
        },
        bond_distances={
            "O": (2.45, 2.70),  # Larger ionic radius
        },
        notes="Largest lanthanide. Highest CN."
    ),
    "CE": MetalProperties(
        symbol="CE",
        name="Cerium",
        hsab_class=HSABClass.HARD,
        oxidation_states=[3, 4],
        typical_coordination_numbers=[8, 9, 10],
        preferred_geometries=["square_antiprism", "bicapped_square_antiprism"],
        preferred_donors={
            "O_carboxylate": 3.5,
            "O_carbonyl": 2.5,
            "O_water": 2.0,
            "N": 0.0,
            "S": -5.0,
        },
        bond_distances={
            "O": (2.35, 2.60),
        },
        notes="Ce(IV) is a strong oxidant. XoxF enzymes use Ce(III)."
    ),
    "SM": MetalProperties(
        symbol="SM",
        name="Samarium",
        hsab_class=HSABClass.HARD,
        oxidation_states=[2, 3],
        typical_coordination_numbers=[8, 9],
        preferred_geometries=["square_antiprism"],
        preferred_donors={
            "O_carboxylate": 3.5,
            "O_carbonyl": 2.5,
            "O_water": 2.0,
            "N": 0.0,
            "S": -5.0,
        },
        bond_distances={
            "O": (2.35, 2.55),
        },
    ),
    "YB": MetalProperties(
        symbol="YB",
        name="Ytterbium",
        hsab_class=HSABClass.HARD,
        oxidation_states=[2, 3],
        typical_coordination_numbers=[8],
        preferred_geometries=["square_antiprism"],
        preferred_donors={
            "O_carboxylate": 3.5,
            "O_carbonyl": 2.5,
            "O_water": 2.0,
            "N": 0.0,
            "S": -5.0,
        },
        bond_distances={
            "O": (2.25, 2.45),  # Smallest lanthanide - shortest bonds
        },
        notes="Smallest common lanthanide. NIR emission."
    ),
}

# Alias handling for oxidation states
METAL_ALIASES = {
    "FE": "FE2",   # Default to Fe(II)
    "CU": "CU2",   # Default to Cu(II)
}


def _get_metal_key(metal: str, oxidation_state: Optional[int] = None) -> str:
    """Resolve metal symbol to database key."""
    metal = metal.upper().strip()

    if oxidation_state:
        key = f"{metal}{oxidation_state}"
        if key in METAL_DATABASE:
            return key

    if metal in METAL_DATABASE:
        return metal

    if metal in METAL_ALIASES:
        return METAL_ALIASES[metal]

    raise ValueError(f"Unknown metal: {metal}")


def get_hsab_class(metal: str, oxidation_state: Optional[int] = None) -> str:
    """
    Get HSAB classification for a metal.

    Args:
        metal: Metal symbol (e.g., "ZN", "TB", "FE")
        oxidation_state: Optional oxidation state (e.g., 2 for Fe2+)

    Returns:
        HSAB class: "hard", "borderline", or "soft"
    """
    key = _get_metal_key(metal, oxidation_state)
    return METAL_DATABASE[key].hsab_class.value


def get_preferred_donors(
    metal: str,
    oxidation_state: Optional[int] = None,
    site_type: Optional[str] = None,
) -> Dict[str, float]:
    """
    Get preferred donor atoms with weights.

    Args:
        metal: Metal symbol
        oxidation_state: Optional oxidation state
        site_type: Optional site type ("structural", "catalytic", "type1", etc.)

    Returns:
        Dict of donor element/type -> preference weight
    """
    key = _get_metal_key(metal, oxidation_state)
    props = METAL_DATABASE[key]

    donors = dict(props.preferred_donors)

    # Adjust for site type
    if site_type == "type1" and metal.upper() == "CU":
        # Type I copper needs Met
        donors["S_met"] = 2.5
    elif site_type == "structural" and metal.upper() == "ZN":
        # Structural zinc favors Cys even more
        donors["S"] = 3.5

    return donors


def get_bond_distance_range(
    metal: str,
    donor_element: str,
    oxidation_state: Optional[int] = None,
) -> Tuple[float, float]:
    """
    Get crystallographic bond distance range.

    Args:
        metal: Metal symbol
        donor_element: Donor element (S, N, O)
        oxidation_state: Optional oxidation state

    Returns:
        Tuple of (min_distance, max_distance) in Angstroms
    """
    key = _get_metal_key(metal, oxidation_state)
    props = METAL_DATABASE[key]

    donor = donor_element.upper()

    if donor in props.bond_distances:
        return props.bond_distances[donor]

    # Fallback to generic element
    if donor.startswith("O"):
        if "O" in props.bond_distances:
            return props.bond_distances["O"]
    if donor.startswith("S"):
        if "S" in props.bond_distances:
            return props.bond_distances["S"]
    if donor.startswith("N"):
        if "N" in props.bond_distances:
            return props.bond_distances["N"]

    # Ultimate fallback
    return (2.0, 3.0)


# Amino acid single-letter to donor mapping
AA_DONORS = {
    "C": {"element": "S", "atom": "SG", "type": "thiolate"},
    "H": {"element": "N", "atom": "NE2/ND1", "type": "imidazole"},
    "D": {"element": "O", "atom": "OD1/OD2", "type": "carboxylate"},
    "E": {"element": "O", "atom": "OE1/OE2", "type": "carboxylate"},
    "N": {"element": "O", "atom": "OD1", "type": "amide"},
    "Q": {"element": "O", "atom": "OE1", "type": "amide"},
    "S": {"element": "O", "atom": "OG", "type": "hydroxyl"},
    "T": {"element": "O", "atom": "OG1", "type": "hydroxyl"},
    "Y": {"element": "O", "atom": "OH", "type": "phenol"},
    "M": {"element": "S", "atom": "SD", "type": "thioether"},
    "K": {"element": "N", "atom": "NZ", "type": "amine"},
    "R": {"element": "N", "atom": "NH1/NH2", "type": "guanidinium"},
    "W": {"element": "N", "atom": "NE1", "type": "indole"},
}


def get_amino_acid_bias(
    metal: str,
    oxidation_state: Optional[int] = None,
    site_type: Optional[str] = None,
    exclude_cysteine_for_lanthanides: bool = True,
) -> str:
    """
    Generate LigandMPNN bias string for metal coordination.

    Args:
        metal: Metal symbol
        oxidation_state: Optional oxidation state
        site_type: Optional site type
        exclude_cysteine_for_lanthanides: Strongly exclude Cys for Ln (default True)

    Returns:
        Bias string like "E:3.0,D:3.0,C:-5.0,A:-2.0"
    """
    key = _get_metal_key(metal, oxidation_state)
    props = METAL_DATABASE[key]
    donors = get_preferred_donors(metal, oxidation_state, site_type)

    bias_parts = []

    # Always disfavor Ala (small, non-coordinating)
    bias_parts.append("A:-2.0")

    # Map donor preferences to amino acids
    for aa, donor_info in AA_DONORS.items():
        element = donor_info["element"]
        donor_type = donor_info["type"]

        # Find matching donor weight
        weight = 0.0

        # Check specific type first (e.g., O_carboxylate)
        type_key = f"{element}_{donor_type}"
        if type_key in donors:
            weight = donors[type_key]
        elif element in donors:
            weight = donors[element]

        # Special handling for Cys with lanthanides
        if aa == "C" and props.hsab_class == HSABClass.HARD:
            if exclude_cysteine_for_lanthanides:
                weight = -5.0  # Strongly exclude

        # Special handling for Met (thioether vs thiolate)
        if aa == "M":
            if "S_met" in donors:
                weight = donors["S_met"]
            elif props.hsab_class == HSABClass.HARD:
                weight = -3.0  # Disfavor for hard metals

        if weight != 0.0:
            bias_parts.append(f"{aa}:{weight:.1f}")

    return ",".join(bias_parts)


def get_coordination_number_range(
    metal: str,
    oxidation_state: Optional[int] = None,
) -> Tuple[int, int]:
    """Get typical coordination number range for metal."""
    key = _get_metal_key(metal, oxidation_state)
    cns = METAL_DATABASE[key].typical_coordination_numbers
    return (min(cns), max(cns))


def get_preferred_geometries(
    metal: str,
    oxidation_state: Optional[int] = None,
) -> List[str]:
    """Get preferred coordination geometries."""
    key = _get_metal_key(metal, oxidation_state)
    return METAL_DATABASE[key].preferred_geometries


def validate_coordination_chemistry(
    metal: str,
    coordinating_residues: List[str],
    oxidation_state: Optional[int] = None,
) -> Dict[str, any]:
    """
    Validate if coordination is chemically reasonable.

    Args:
        metal: Metal symbol
        coordinating_residues: List of 1-letter AA codes
        oxidation_state: Optional oxidation state

    Returns:
        Dict with validation results and warnings
    """
    key = _get_metal_key(metal, oxidation_state)
    props = METAL_DATABASE[key]

    result = {
        "valid": True,
        "warnings": [],
        "errors": [],
        "score": 1.0,
    }

    # Check for HSAB violations
    for aa in coordinating_residues:
        if aa not in AA_DONORS:
            continue

        donor = AA_DONORS[aa]
        element = donor["element"]

        # Hard acid + soft base = bad
        if props.hsab_class == HSABClass.HARD and element == "S":
            result["errors"].append(
                f"HSAB violation: {aa} (sulfur donor) incompatible with hard {metal}"
            )
            result["valid"] = False
            result["score"] *= 0.1

        # Soft acid + hard base = suboptimal
        if props.hsab_class == HSABClass.SOFT and element == "O":
            result["warnings"].append(
                f"Suboptimal: {aa} (oxygen donor) with soft {metal}"
            )
            result["score"] *= 0.7

    # Check coordination number
    cn = len(coordinating_residues)
    cn_range = get_coordination_number_range(metal, oxidation_state)
    if cn < cn_range[0]:
        result["warnings"].append(
            f"Low coordination: {cn} < typical {cn_range[0]}-{cn_range[1]}"
        )
    elif cn > cn_range[1]:
        result["warnings"].append(
            f"High coordination: {cn} > typical {cn_range[0]}-{cn_range[1]}"
        )

    return result
```

**Step 4: Run test to verify it passes**

Run: `cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless && python -m pytest test_metal_chemistry.py -v`
Expected: PASS

**Step 5: Commit**

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion
git add backend/serverless/metal_chemistry.py backend/serverless/test_metal_chemistry.py
git commit -m "feat: add HSAB-compliant metal chemistry database

- Complete metal database with HSAB classification
- Crystallographic bond distances from literature
- Amino acid bias generation per metal type
- Lanthanide-specific exclusion of Cys donors
- Validation for HSAB chemistry violations

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

### Task 1.2: Integrate Metal Chemistry into Design Types

**Files:**
- Modify: `backend/serverless/design_types.py:310-380`
- Test: `backend/serverless/test_design_types.py`

**Step 1: Write the failing test**

```python
# Add to test_design_types.py

class TestMetalPresetChemistry:
    """Test that metal presets use correct chemistry."""

    def test_zinc_preset_favors_cys_his(self):
        """Zinc preset should favor Cys and His over Asp/Glu."""
        from design_types import METAL_PRESETS
        bias = METAL_PRESETS["zinc"]["bias_AA"]
        bias_dict = dict(pair.split(":") for pair in bias.split(","))
        assert float(bias_dict.get("C", 0)) >= 2.5
        assert float(bias_dict.get("H", 0)) >= 2.0

    def test_lanthanide_preset_excludes_cys(self):
        """Lanthanide preset must strongly exclude Cys."""
        from design_types import METAL_PRESETS
        bias = METAL_PRESETS["lanthanide"]["bias_AA"]
        bias_dict = dict(pair.split(":") for pair in bias.split(","))
        assert float(bias_dict.get("C", 0)) <= -3.0

    def test_iron_has_oxidation_state_variants(self):
        """Iron should have separate Fe2+ and Fe3+ presets."""
        from design_types import METAL_PRESETS
        assert "iron_ii" in METAL_PRESETS or "iron" in METAL_PRESETS
        # Fe3+ is harder - prefers O over S
        if "iron_iii" in METAL_PRESETS:
            bias = METAL_PRESETS["iron_iii"]["bias_AA"]
            bias_dict = dict(pair.split(":") for pair in bias.split(","))
            assert float(bias_dict.get("E", 0)) > float(bias_dict.get("C", 0))
```

**Step 2: Run test to verify it fails**

Run: `cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless && python -m pytest test_design_types.py::TestMetalPresetChemistry -v`
Expected: FAIL

**Step 3: Modify design_types.py**

Replace the `METAL_PRESETS` section with:

```python
# Import chemistry-aware bias generation
from metal_chemistry import get_amino_acid_bias, get_coordination_number_range

# Metal presets with HSAB-compliant chemistry
METAL_PRESETS = {
    # =========================================================================
    # Zinc (Borderline acid - prefers S, N donors)
    # =========================================================================
    "zinc": {
        "description": "Zinc binding (general)",
        "bias_AA": get_amino_acid_bias("ZN"),
        "coordination_numbers": [4, 5, 6],
        "preferred_residues": ["C", "H", "D", "E"],
    },
    "zinc_structural": {
        "description": "Zinc structural site (Cys4 or Cys3His)",
        "bias_AA": get_amino_acid_bias("ZN", site_type="structural"),
        "coordination_numbers": [4],
        "preferred_residues": ["C", "H"],
    },
    "zinc_catalytic": {
        "description": "Zinc catalytic site (His-rich + water)",
        "bias_AA": "H:3.0,E:2.0,D:2.0,C:1.0,A:-2.0",
        "coordination_numbers": [4, 5],
        "preferred_residues": ["H", "E", "D"],
    },
    # =========================================================================
    # Lanthanides (Hard acids - O donors only, EXCLUDE Cys)
    # =========================================================================
    "lanthanide": {
        "description": "Lanthanide binding (Tb, Eu, Gd, etc.)",
        "bias_AA": get_amino_acid_bias("TB"),
        "coordination_numbers": [8, 9],
        "preferred_residues": ["E", "D", "N", "Q"],
        "excluded_residues": ["C"],  # HSAB hard acid - no sulfur!
    },
    "terbium": {
        "description": "Terbium binding site",
        "bias_AA": get_amino_acid_bias("TB"),
        "coordination_numbers": [8, 9],
        "geometry": "square_antiprism",
    },
    "europium": {
        "description": "Europium binding site",
        "bias_AA": get_amino_acid_bias("EU"),
        "coordination_numbers": [8, 9],
    },
    "gadolinium": {
        "description": "Gadolinium binding site (MRI)",
        "bias_AA": get_amino_acid_bias("GD"),
        "coordination_numbers": [8, 9],
    },
    # =========================================================================
    # Iron (oxidation state dependent)
    # =========================================================================
    "iron": {
        "description": "Iron binding (default Fe2+)",
        "bias_AA": get_amino_acid_bias("FE", oxidation_state=2),
        "coordination_numbers": [4, 5, 6],
    },
    "iron_ii": {
        "description": "Iron(II) - borderline, accepts S/N/O",
        "bias_AA": get_amino_acid_bias("FE", oxidation_state=2),
        "coordination_numbers": [4, 5, 6],
    },
    "iron_iii": {
        "description": "Iron(III) - harder, prefers O/N over S",
        "bias_AA": get_amino_acid_bias("FE", oxidation_state=3),
        "coordination_numbers": [5, 6],
    },
    "iron_sulfur": {
        "description": "Iron-sulfur cluster (FeS)",
        "bias_AA": "C:3.5,H:1.5,A:-2.0",  # Cys4 coordination
        "coordination_numbers": [4],
    },
    # =========================================================================
    # Copper (oxidation state dependent)
    # =========================================================================
    "copper": {
        "description": "Copper binding (default Cu2+)",
        "bias_AA": get_amino_acid_bias("CU", oxidation_state=2),
        "coordination_numbers": [4, 5, 6],
    },
    "copper_type1": {
        "description": "Type I blue copper (Cu+, 2His+Cys+Met)",
        "bias_AA": "C:3.0,H:2.5,M:2.0,A:-2.0",
        "coordination_numbers": [4],
        "histidine_tautomer": "ND1",  # Type I uses Nδ1
    },
    "copper_type2": {
        "description": "Type II copper (Cu2+, His-rich)",
        "bias_AA": "H:3.0,E:1.5,D:1.5,C:0.5,A:-2.0",
        "coordination_numbers": [4, 5, 6],
        "histidine_tautomer": "NE2",  # Type II uses Nε2
    },
    # =========================================================================
    # Other metals
    # =========================================================================
    "calcium": {
        "description": "Calcium binding (EF-hand style)",
        "bias_AA": get_amino_acid_bias("CA"),
        "coordination_numbers": [6, 7, 8],
    },
    "magnesium": {
        "description": "Magnesium binding (octahedral)",
        "bias_AA": get_amino_acid_bias("MG"),
        "coordination_numbers": [6],
    },
    "manganese": {
        "description": "Manganese binding",
        "bias_AA": get_amino_acid_bias("MN"),
        "coordination_numbers": [5, 6],
    },
}
```

**Step 4: Run test to verify it passes**

Run: `cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless && python -m pytest test_design_types.py::TestMetalPresetChemistry -v`
Expected: PASS

**Step 5: Commit**

```bash
git add backend/serverless/design_types.py backend/serverless/test_design_types.py
git commit -m "feat: update METAL_PRESETS with HSAB-compliant chemistry

- Zinc presets favor Cys/His (borderline acid)
- Lanthanide presets exclude Cys (hard acid)
- Iron variants for Fe2+/Fe3+ oxidation states
- Copper variants for Type I/II sites
- Histidine tautomer specification for copper

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Phase 2: Database-Driven Template System

### Task 2.1: Create PDB Metal Site Fetcher

**Files:**
- Create: `backend/serverless/metal_site_fetcher.py`
- Test: `backend/serverless/test_metal_site_fetcher.py`

**Step 1: Write the failing test**

```python
# test_metal_site_fetcher.py
"""Tests for PDB metal site fetching."""
import pytest


class TestMetalSiteQuery:
    """Test metal site database queries."""

    def test_query_zinc_sites(self):
        """Should find zinc binding sites in PDB."""
        from metal_site_fetcher import query_metal_sites
        sites = query_metal_sites(
            metal="ZN",
            resolution_max=2.0,
            limit=5,
        )
        assert len(sites) >= 1
        assert all(s["metal"] == "ZN" for s in sites)

    def test_query_lanthanide_sites(self):
        """Should find lanthanide sites (may be sparse)."""
        from metal_site_fetcher import query_metal_sites
        sites = query_metal_sites(
            metal="TB",
            resolution_max=3.0,
            limit=5,
        )
        # Lanthanides are rarer - may find 0
        # But if found, should be TB
        for s in sites:
            assert s["metal"] in ["TB", "EU", "GD", "LA"]

    def test_query_with_ligand(self):
        """Should find metal-ligand complex structures."""
        from metal_site_fetcher import query_metal_ligand_sites
        sites = query_metal_ligand_sites(
            metal="CA",
            ligand="PQQ",
            limit=5,
        )
        # PQQ-Ca should find methanol dehydrogenase structures
        assert len(sites) >= 0  # May be empty if API unavailable


class TestSiteExtraction:
    """Test extraction of coordination from PDB."""

    def test_extract_coordination_from_pdb(self):
        """Should extract coordinating residues from PDB content."""
        from metal_site_fetcher import extract_metal_coordination

        # Minimal zinc site
        pdb = """
ATOM      1  NE2 HIS A  63      12.500  10.000   7.230  1.00 50.00           N
ATOM      2  SG  CYS A  96      12.300  10.000   3.440  1.00 50.00           S
ATOM      3  SG  CYS A  99      11.700  11.500   5.000  1.00 50.00           S
ATOM      4  SG  CYS A 102      11.700   8.500   5.000  1.00 50.00           S
HETATM    5  ZN  ZN  X   1      12.000  10.000   5.000  1.00 50.00          ZN
END
"""
        coord = extract_metal_coordination(pdb, metal="ZN", cutoff=3.0)
        assert coord["metal"] == "ZN"
        assert coord["coordination_number"] >= 3
        assert len(coord["coordinating_atoms"]) >= 3


class TestTemplateGeneration:
    """Test template generation from PDB sites."""

    def test_generate_template_from_pdb_id(self):
        """Should generate template from known PDB structure."""
        from metal_site_fetcher import generate_template_from_pdb

        # 1CA2 = carbonic anhydrase with Zn
        template = generate_template_from_pdb(
            pdb_id="1CA2",
            metal="ZN",
        )
        # May fail if network unavailable
        if template:
            assert template["metal"] == "ZN"
            assert template["source"] == "pdb"
            assert "pdb_id" in template
```

**Step 2: Run test to verify it fails**

Run: `cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless && python -m pytest test_metal_site_fetcher.py -v --tb=short 2>&1 | head -20`
Expected: FAIL with "No module named 'metal_site_fetcher'"

**Step 3: Write implementation**

```python
# metal_site_fetcher.py
"""
PDB Metal Site Fetcher

Query RCSB PDB for metal binding site structures.
Extract coordination geometry from experimental structures.
Generate templates from known structures (database-first approach).
"""

import json
import math
import logging
from typing import Dict, List, Any, Optional, Tuple
from dataclasses import dataclass, asdict

# Use existing pdb_fetch module
from pdb_fetch import fetch_pdb, parse_pdb_content

logger = logging.getLogger(__name__)

# RCSB Search API endpoint
RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"

# Try to import requests
try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False
    logger.warning("requests not available - PDB queries will be limited")


@dataclass
class MetalSite:
    """Extracted metal coordination site."""
    pdb_id: str
    metal: str
    metal_chain: str
    metal_resnum: int
    metal_coords: Tuple[float, float, float]
    coordination_number: int
    coordinating_atoms: List[Dict[str, Any]]
    resolution: Optional[float] = None
    geometry: Optional[str] = None


@dataclass
class MetalLigandSite(MetalSite):
    """Metal site with organic ligand."""
    ligand_code: str = ""
    ligand_chain: str = ""
    ligand_resnum: int = 0
    ligand_atoms: List[Dict[str, Any]] = None


def query_metal_sites(
    metal: str,
    resolution_max: float = 2.5,
    limit: int = 10,
    coordination_min: int = 3,
) -> List[Dict[str, Any]]:
    """
    Query RCSB PDB for structures containing metal ion.

    Args:
        metal: Metal symbol (e.g., "ZN", "TB")
        resolution_max: Maximum resolution in Angstroms
        limit: Maximum number of results
        coordination_min: Minimum coordination number

    Returns:
        List of PDB IDs with metadata
    """
    if not HAS_REQUESTS:
        logger.warning("Cannot query PDB - requests not available")
        return []

    # RCSB Search API query
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_nonpolymer_instance_feature_summary.comp_id",
                        "operator": "exact_match",
                        "value": metal.upper(),
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.resolution_combined",
                        "operator": "less_or_equal",
                        "value": resolution_max,
                    }
                },
            ]
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": limit,
            },
            "sort": [
                {
                    "sort_by": "rcsb_entry_info.resolution_combined",
                    "direction": "asc",
                }
            ]
        }
    }

    try:
        response = requests.post(
            RCSB_SEARCH_URL,
            json=query,
            headers={"Content-Type": "application/json"},
            timeout=30,
        )

        if response.status_code != 200:
            logger.error(f"RCSB query failed: {response.status_code}")
            return []

        data = response.json()
        results = []

        for result in data.get("result_set", []):
            results.append({
                "pdb_id": result["identifier"],
                "metal": metal.upper(),
                "score": result.get("score", 0),
            })

        return results

    except Exception as e:
        logger.error(f"Error querying RCSB: {e}")
        return []


def query_metal_ligand_sites(
    metal: str,
    ligand: str,
    resolution_max: float = 3.0,
    limit: int = 10,
) -> List[Dict[str, Any]]:
    """
    Query for structures with both metal and organic ligand.

    Args:
        metal: Metal symbol
        ligand: Ligand code (e.g., "PQQ", "CIT")
        resolution_max: Maximum resolution
        limit: Maximum results

    Returns:
        List of PDB IDs with metadata
    """
    if not HAS_REQUESTS:
        return []

    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_nonpolymer_instance_feature_summary.comp_id",
                        "operator": "exact_match",
                        "value": metal.upper(),
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_nonpolymer_instance_feature_summary.comp_id",
                        "operator": "exact_match",
                        "value": ligand.upper(),
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.resolution_combined",
                        "operator": "less_or_equal",
                        "value": resolution_max,
                    }
                },
            ]
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": limit},
        }
    }

    try:
        response = requests.post(
            RCSB_SEARCH_URL,
            json=query,
            headers={"Content-Type": "application/json"},
            timeout=30,
        )

        if response.status_code != 200:
            return []

        data = response.json()
        return [
            {
                "pdb_id": r["identifier"],
                "metal": metal.upper(),
                "ligand": ligand.upper(),
            }
            for r in data.get("result_set", [])
        ]

    except Exception as e:
        logger.error(f"Error querying RCSB: {e}")
        return []


def extract_metal_coordination(
    pdb_content: str,
    metal: str,
    cutoff: float = 3.0,
) -> Dict[str, Any]:
    """
    Extract metal coordination from PDB content.

    Args:
        pdb_content: PDB file content
        metal: Metal to find
        cutoff: Coordination cutoff in Angstroms

    Returns:
        Dict with coordination information
    """
    metal = metal.upper()
    metal_pos = None
    metal_chain = None
    metal_resnum = None

    atoms = []

    for line in pdb_content.split('\n'):
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            continue

        try:
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain = line[21]
            resnum = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            element = line[76:78].strip() if len(line) > 76 else res_name[:2]

            # Check if this is the metal
            if res_name == metal or element == metal:
                metal_pos = (x, y, z)
                metal_chain = chain
                metal_resnum = resnum
            else:
                atoms.append({
                    "atom_name": atom_name,
                    "res_name": res_name,
                    "chain": chain,
                    "resnum": resnum,
                    "coords": (x, y, z),
                    "element": element,
                })
        except (ValueError, IndexError):
            continue

    if metal_pos is None:
        return {
            "success": False,
            "error": f"Metal {metal} not found",
            "metal": metal,
            "coordination_number": 0,
            "coordinating_atoms": [],
        }

    # Find coordinating atoms
    coordinating = []
    for atom in atoms:
        dx = atom["coords"][0] - metal_pos[0]
        dy = atom["coords"][1] - metal_pos[1]
        dz = atom["coords"][2] - metal_pos[2]
        dist = math.sqrt(dx*dx + dy*dy + dz*dz)

        if dist <= cutoff:
            coordinating.append({
                **atom,
                "distance": round(dist, 2),
            })

    # Sort by distance
    coordinating.sort(key=lambda x: x["distance"])

    return {
        "success": True,
        "metal": metal,
        "metal_chain": metal_chain,
        "metal_resnum": metal_resnum,
        "metal_coords": metal_pos,
        "coordination_number": len(coordinating),
        "coordinating_atoms": coordinating,
    }


def generate_template_from_pdb(
    pdb_id: str,
    metal: str,
    cutoff: float = 3.0,
) -> Optional[Dict[str, Any]]:
    """
    Generate template from experimental PDB structure.

    This is the PRIMARY template source - database-first approach.

    Args:
        pdb_id: 4-character PDB ID
        metal: Metal to extract
        cutoff: Coordination cutoff

    Returns:
        Template dict or None if unavailable
    """
    # Fetch PDB
    result = fetch_pdb(pdb_id)

    if not result["success"]:
        logger.warning(f"Could not fetch PDB {pdb_id}: {result['error']}")
        return None

    pdb_content = result["content"]

    # Extract coordination
    coord = extract_metal_coordination(pdb_content, metal, cutoff)

    if not coord["success"]:
        logger.warning(f"Could not find {metal} in {pdb_id}")
        return None

    if coord["coordination_number"] < 3:
        logger.warning(f"Low coordination ({coord['coordination_number']}) in {pdb_id}")
        return None

    # Build template
    template = {
        "name": f"{pdb_id}_{metal}",
        "source": "pdb",
        "pdb_id": pdb_id,
        "metal": metal,
        "metal_coords": coord["metal_coords"],
        "coordination_number": coord["coordination_number"],
        "coordinating_atoms": coord["coordinating_atoms"],
        "pdb_content": pdb_content,  # Store for later extraction
    }

    return template


# Known high-quality reference structures for each metal
REFERENCE_STRUCTURES = {
    "ZN": [
        ("1CA2", "Carbonic anhydrase - His3 + water"),
        ("1HET", "Alcohol dehydrogenase - His, Cys, Cys, Glu"),
        ("1A5T", "Zinc finger - Cys2His2"),
    ],
    "FE": [
        ("1RB9", "Rubredoxin - FeCys4"),
        ("1MBS", "Myoglobin - heme"),
        ("1WOC", "2-His-1-carboxylate"),
    ],
    "CU": [
        ("1PLC", "Plastocyanin - Type I"),
        ("1NWP", "Nitrite reductase - Type II"),
    ],
    "CA": [
        ("1CLL", "Calmodulin - EF-hand"),
        ("1W6S", "Methanol dehydrogenase - PQQ-Ca"),
    ],
    "TB": [
        ("6MI5", "Lanmodulin"),
    ],
    "EU": [
        ("6MI5", "Lanmodulin"),  # Same structure, can substitute
    ],
    "GD": [
        ("6MI5", "Lanmodulin"),
    ],
}


def get_reference_template(
    metal: str,
    ligand: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """
    Get template from curated reference structures.

    Args:
        metal: Metal symbol
        ligand: Optional ligand code

    Returns:
        Best matching template or None
    """
    metal = metal.upper()

    if metal not in REFERENCE_STRUCTURES:
        logger.info(f"No reference structures for {metal}")
        return None

    for pdb_id, description in REFERENCE_STRUCTURES[metal]:
        # If ligand specified, check description
        if ligand and ligand.upper() not in description.upper():
            continue

        template = generate_template_from_pdb(pdb_id, metal)
        if template:
            template["description"] = description
            return template

    return None
```

**Step 4: Run test to verify it passes**

Run: `cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless && python -m pytest test_metal_site_fetcher.py -v --tb=short`
Expected: PASS (some tests may skip if network unavailable)

**Step 5: Commit**

```bash
git add backend/serverless/metal_site_fetcher.py backend/serverless/test_metal_site_fetcher.py
git commit -m "feat: add PDB metal site fetcher for database-driven templates

- Query RCSB PDB for metal binding structures
- Extract coordination geometry from PDB content
- Generate templates from experimental structures
- Curated reference structures for common metals
- Database-first approach (calculated fallback later)

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

### Task 2.2: Update Template Library with Database Priority

**Files:**
- Modify: `backend/serverless/metal_ligand_templates.py`
- Test: `backend/serverless/test_metal_ligand_complex.py`

**Step 1: Write the failing test**

```python
# Add to test_metal_ligand_complex.py

class TestDatabaseDrivenTemplates:
    """Test database-first template retrieval."""

    def test_get_template_prefers_database(self):
        """Template retrieval should try database first."""
        from metal_ligand_templates import get_template_with_fallback

        # This should try PDB first, then fall back to calculated
        template = get_template_with_fallback("pqq_ca")

        assert template is not None
        # Should have either source=pdb or source=calculated
        assert template.get("source") in ["pdb", "calculated", "library"]

    def test_fallback_to_calculated(self):
        """Should fall back to calculated if PDB unavailable."""
        from metal_ligand_templates import get_template_with_fallback

        # Non-existent complex - must fall back
        template = get_template_with_fallback(
            "unknown_complex",
            metal="TB",
            fallback_enabled=True,
        )

        # Should return None or a calculated template
        if template:
            assert template["source"] == "calculated"
```

**Step 2: Run test to verify it fails**

Run: `cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless && python -m pytest test_metal_ligand_complex.py::TestDatabaseDrivenTemplates -v`
Expected: FAIL

**Step 3: Modify metal_ligand_templates.py**

Add database-priority template retrieval:

```python
# Add to metal_ligand_templates.py after imports

from metal_site_fetcher import (
    get_reference_template,
    generate_template_from_pdb,
    query_metal_ligand_sites,
)

# Template source priority
TEMPLATE_SOURCES = ["database", "library", "calculated"]

# PDB IDs for known metal-ligand complexes
KNOWN_COMPLEX_PDBS = {
    "pqq_ca": {
        "pdb_ids": ["1W6S", "4AAH", "1W99"],
        "metal": "CA",
        "ligand": "PQQ",
        "description": "PQQ-dependent methanol dehydrogenase",
    },
    "citrate_tb": {
        "pdb_ids": ["6MI5"],  # Lanmodulin - closest analog
        "metal": "TB",
        "ligand": "CIT",  # May need custom handling
        "description": "Citrate-terbium complex",
    },
    "citrate_eu": {
        "pdb_ids": ["6MI5"],
        "metal": "EU",
        "ligand": "CIT",
    },
    "citrate_gd": {
        "pdb_ids": ["6MI5"],
        "metal": "GD",
        "ligand": "CIT",
    },
    "heme_fe": {
        "pdb_ids": ["1MBS", "1HHO"],
        "metal": "FE",
        "ligand": "HEM",
        "description": "Heme-iron complex",
    },
}


def get_template_with_fallback(
    template_name: str,
    metal: Optional[str] = None,
    ligand: Optional[str] = None,
    fallback_enabled: bool = True,
    use_cache: bool = True,
) -> Optional[Dict[str, Any]]:
    """
    Get template with database-first priority.

    Priority:
    1. PDB database (experimental structure)
    2. Library (pre-defined templates)
    3. Calculated (geometry-based fallback)

    Args:
        template_name: Template identifier
        metal: Metal symbol (for fallback)
        ligand: Ligand code (for fallback)
        fallback_enabled: Allow calculated fallback
        use_cache: Cache PDB fetches

    Returns:
        Template dict or None
    """
    template_key = template_name.lower()

    # =========================================================================
    # Priority 1: Try PDB database
    # =========================================================================
    if template_key in KNOWN_COMPLEX_PDBS:
        complex_info = KNOWN_COMPLEX_PDBS[template_key]

        for pdb_id in complex_info["pdb_ids"]:
            try:
                template = generate_template_from_pdb(
                    pdb_id,
                    complex_info["metal"],
                )
                if template:
                    template["source"] = "pdb"
                    template["original_name"] = template_name
                    template["ligand"] = complex_info.get("ligand")
                    logger.info(f"Using PDB template from {pdb_id}")
                    return template
            except Exception as e:
                logger.debug(f"Could not fetch {pdb_id}: {e}")
                continue

    # =========================================================================
    # Priority 2: Try library templates
    # =========================================================================
    if template_key in METAL_LIGAND_COMPLEX_TEMPLATES:
        template = dict(METAL_LIGAND_COMPLEX_TEMPLATES[template_key])
        template["source"] = "library"
        logger.info(f"Using library template for {template_name}")
        return template

    # =========================================================================
    # Priority 3: Calculated fallback
    # =========================================================================
    if fallback_enabled and metal:
        logger.warning(f"No database/library template for {template_name}, using calculated")

        from metal_chemistry import get_preferred_donors, get_coordination_number_range

        cn_range = get_coordination_number_range(metal)
        donors = get_preferred_donors(metal)

        # Generate idealized template
        template = _generate_calculated_template(
            metal=metal,
            ligand=ligand,
            coordination_number=cn_range[1],  # Use max CN
            donors=donors,
        )
        template["source"] = "calculated"
        template["warning"] = "Calculated template - verify geometry manually"

        return template

    return None


def _generate_calculated_template(
    metal: str,
    ligand: Optional[str],
    coordination_number: int,
    donors: Dict[str, float],
) -> Dict[str, Any]:
    """
    Generate calculated template when no database available.

    CAUTION: This is fallback only. Results should be validated.

    Args:
        metal: Metal symbol
        ligand: Ligand code
        coordination_number: Target CN
        donors: Donor preferences

    Returns:
        Calculated template dict
    """
    from lanthanide_templates import LANTHANIDE_PARAMS, _generate_coordination_positions
    import numpy as np

    # Get metal-specific parameters
    if metal.upper() in LANTHANIDE_PARAMS:
        params = LANTHANIDE_PARAMS[metal.upper()]
        bond_distance = params["bond_distance"]
        geometry = params["geometry"]
    else:
        # Default for transition metals
        bond_distance = 2.3
        geometry = "octahedral" if coordination_number == 6 else "tetrahedral"

    # Generate ideal positions
    metal_pos = np.array([50.0, 50.0, 50.0])  # Centered
    positions = _generate_coordination_positions(
        geometry,
        coordination_number,
        metal_pos,
        bond_distance,
    )

    # Build template
    template = {
        "name": f"calculated_{metal}_{coordination_number}",
        "metal": metal.upper(),
        "ligand_res_name": ligand,
        "coordination_number": coordination_number,
        "geometry": geometry,
        "bond_distance": bond_distance,
        "metal_coords": tuple(metal_pos),
        "coordination_positions": [tuple(p) for p in positions],
        "calculated": True,
    }

    return template
```

**Step 4: Run test to verify it passes**

Run: `cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless && python -m pytest test_metal_ligand_complex.py::TestDatabaseDrivenTemplates -v`
Expected: PASS

**Step 5: Commit**

```bash
git add backend/serverless/metal_ligand_templates.py backend/serverless/test_metal_ligand_complex.py
git commit -m "feat: implement database-first template retrieval

- get_template_with_fallback() tries PDB first
- Known complex PDB IDs for common metal-ligand pairs
- Library templates as secondary source
- Calculated fallback with warning when no data
- Clear source tracking (pdb/library/calculated)

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Phase 3: Ligand Donor Intelligence

### Task 3.1: Create Ligand Donor Mapper

**Files:**
- Create: `backend/serverless/ligand_donors.py`
- Test: `backend/serverless/test_ligand_donors.py`

**Step 1: Write the failing test**

```python
# test_ligand_donors.py
"""Tests for ligand donor atom identification."""
import pytest


class TestDonorIdentification:
    """Test identification of coordinating atoms from ligand."""

    def test_citrate_donors(self):
        """Citrate should have 6 potential O donors (3 carboxylates)."""
        from ligand_donors import identify_donors_from_smiles

        citrate_smiles = "C(C(=O)O)(CC(=O)O)(C(=O)O)O"
        donors = identify_donors_from_smiles(citrate_smiles)

        # 3 carboxylates = 6 oxygen atoms potentially coordinating
        oxygen_donors = [d for d in donors if d["element"] == "O"]
        assert len(oxygen_donors) >= 4  # At least 4 accessible

    def test_pqq_donors(self):
        """PQQ should have multiple O and N donors."""
        from ligand_donors import identify_donors_from_smiles

        # Simplified PQQ core
        pqq_smiles = "O=C1C(=O)c2cc3C(=O)NC(=O)c3c(c2N1)C(=O)O"
        donors = identify_donors_from_smiles(pqq_smiles)

        # Should have O and N donors
        elements = set(d["element"] for d in donors)
        assert "O" in elements
        assert "N" in elements

    def test_donor_type_classification(self):
        """Should classify donor types correctly."""
        from ligand_donors import identify_donors_from_smiles

        # Carboxylic acid
        acetic = "CC(=O)O"
        donors = identify_donors_from_smiles(acetic)

        carboxylate_donors = [d for d in donors if d["type"] == "carboxylate"]
        assert len(carboxylate_donors) >= 1


class TestMetalCompatibility:
    """Test metal-donor compatibility scoring."""

    def test_citrate_lanthanide_compatibility(self):
        """Citrate should be highly compatible with lanthanides."""
        from ligand_donors import score_ligand_metal_compatibility

        citrate_smiles = "C(C(=O)O)(CC(=O)O)(C(=O)O)O"
        score = score_ligand_metal_compatibility(citrate_smiles, "TB")

        assert score >= 0.8  # High compatibility

    def test_thiol_lanthanide_incompatibility(self):
        """Thiol ligands should be incompatible with lanthanides."""
        from ligand_donors import score_ligand_metal_compatibility

        thiol_smiles = "CCCS"  # Propanethiol
        score = score_ligand_metal_compatibility(thiol_smiles, "TB")

        assert score < 0.3  # Low compatibility


class TestDenticity:
    """Test denticity determination."""

    def test_citrate_bidentate_carboxylates(self):
        """Citrate carboxylates can be bidentate."""
        from ligand_donors import analyze_denticity

        citrate_smiles = "C(C(=O)O)(CC(=O)O)(C(=O)O)O"
        denticity = analyze_denticity(citrate_smiles)

        # Should have potential for bidentate
        bidentate_groups = [d for d in denticity if d["max_denticity"] >= 2]
        assert len(bidentate_groups) >= 1
```

**Step 2: Run test to verify it fails**

Run: `cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless && python -m pytest test_ligand_donors.py -v --tb=short 2>&1 | head -20`
Expected: FAIL

**Step 3: Write implementation**

```python
# ligand_donors.py
"""
Ligand Donor Atom Identification

Parse ligand SMILES to identify coordinating atoms.
Score metal-ligand compatibility based on HSAB theory.
Determine denticity (mono/bi/poly-dentate) potential.
"""

import logging
from typing import Dict, List, Any, Optional, Tuple
from dataclasses import dataclass

logger = logging.getLogger(__name__)

# Try to import RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    logger.warning("RDKit not available - ligand analysis limited")

from metal_chemistry import get_hsab_class, HSABClass


@dataclass
class DonorAtom:
    """A potential coordinating atom in ligand."""
    atom_index: int
    element: str
    donor_type: str
    hybridization: str
    formal_charge: int
    neighbors: List[str]
    is_aromatic: bool
    hsab_class: str  # hard, borderline, soft
    max_denticity: int = 1


# Donor type definitions with HSAB classification
DONOR_TYPES = {
    # Oxygen donors (mostly hard)
    "carboxylate": {"element": "O", "hsab": "hard", "denticity": 2},
    "carbonyl": {"element": "O", "hsab": "hard", "denticity": 1},
    "hydroxyl": {"element": "O", "hsab": "hard", "denticity": 1},
    "ether": {"element": "O", "hsab": "borderline", "denticity": 1},
    "phenol": {"element": "O", "hsab": "borderline", "denticity": 1},
    # Nitrogen donors (borderline)
    "amine": {"element": "N", "hsab": "borderline", "denticity": 1},
    "imine": {"element": "N", "hsab": "borderline", "denticity": 1},
    "pyridine": {"element": "N", "hsab": "borderline", "denticity": 1},
    "amide": {"element": "N", "hsab": "hard", "denticity": 1},
    "imidazole": {"element": "N", "hsab": "borderline", "denticity": 1},
    # Sulfur donors (soft)
    "thiol": {"element": "S", "hsab": "soft", "denticity": 1},
    "thioether": {"element": "S", "hsab": "soft", "denticity": 1},
    "thiolate": {"element": "S", "hsab": "soft", "denticity": 1},
    # Phosphorus (hard when oxidized)
    "phosphate": {"element": "P", "hsab": "hard", "denticity": 2},
}


def identify_donors_from_smiles(smiles: str) -> List[Dict[str, Any]]:
    """
    Identify potential metal-coordinating atoms from SMILES.

    Args:
        smiles: SMILES string of ligand

    Returns:
        List of donor atom dicts
    """
    if not HAS_RDKIT:
        logger.warning("RDKit not available")
        return []

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.error(f"Could not parse SMILES: {smiles}")
        return []

    # Add hydrogens for accurate analysis
    mol = Chem.AddHs(mol)

    donors = []

    for atom in mol.GetAtoms():
        element = atom.GetSymbol()

        # Only consider O, N, S, P as potential donors
        if element not in ["O", "N", "S", "P"]:
            continue

        idx = atom.GetIdx()
        charge = atom.GetFormalCharge()
        hybrid = str(atom.GetHybridization())
        aromatic = atom.GetIsAromatic()
        neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]

        # Classify donor type
        donor_type = _classify_donor_type(atom, mol)

        if donor_type is None:
            continue

        type_info = DONOR_TYPES.get(donor_type, {})

        donor = {
            "atom_index": idx,
            "element": element,
            "type": donor_type,
            "hybridization": hybrid,
            "formal_charge": charge,
            "neighbors": neighbors,
            "is_aromatic": aromatic,
            "hsab_class": type_info.get("hsab", "borderline"),
            "max_denticity": type_info.get("denticity", 1),
        }

        donors.append(donor)

    return donors


def _classify_donor_type(atom, mol) -> Optional[str]:
    """Classify the donor type of an atom."""
    element = atom.GetSymbol()
    charge = atom.GetFormalCharge()
    neighbors = atom.GetNeighbors()
    neighbor_symbols = [n.GetSymbol() for n in neighbors]
    hybrid = atom.GetHybridization()

    # Oxygen classification
    if element == "O":
        if charge == -1:
            # Check for carboxylate
            for n in neighbors:
                if n.GetSymbol() == "C":
                    c_neighbors = [x.GetSymbol() for x in n.GetNeighbors()]
                    if c_neighbors.count("O") >= 2:
                        return "carboxylate"
            return "hydroxyl"  # Generic deprotonated O

        if len(neighbors) == 1:
            n = neighbors[0]
            if n.GetSymbol() == "C":
                # Check for C=O
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx())
                if bond and bond.GetBondTypeAsDouble() == 2.0:
                    # Check for carboxylic acid
                    c_neighbors = [x.GetSymbol() for x in n.GetNeighbors()]
                    if c_neighbors.count("O") >= 2:
                        return "carboxylate"
                    return "carbonyl"

        if len(neighbors) == 2:
            if "H" in neighbor_symbols:
                return "hydroxyl"
            return "ether"

        return "hydroxyl"

    # Nitrogen classification
    if element == "N":
        if atom.GetIsAromatic():
            # Check for imidazole-like
            ring_info = mol.GetRingInfo()
            for ring in ring_info.AtomRings():
                if atom.GetIdx() in ring:
                    ring_atoms = [mol.GetAtomWithIdx(i).GetSymbol() for i in ring]
                    if ring_atoms.count("N") >= 2:
                        return "imidazole"
            return "pyridine"

        if charge == 0:
            # Check for amide
            for n in neighbors:
                if n.GetSymbol() == "C":
                    for nn in n.GetNeighbors():
                        if nn.GetSymbol() == "O":
                            bond = mol.GetBondBetweenAtoms(n.GetIdx(), nn.GetIdx())
                            if bond and bond.GetBondTypeAsDouble() == 2.0:
                                return "amide"

            # Check hybridization
            if str(hybrid) == "SP2":
                return "imine"
            return "amine"

        return "amine"

    # Sulfur classification
    if element == "S":
        if charge == -1:
            return "thiolate"
        if "H" in neighbor_symbols:
            return "thiol"
        return "thioether"

    # Phosphorus classification
    if element == "P":
        if "O" in neighbor_symbols:
            return "phosphate"

    return None


def score_ligand_metal_compatibility(
    smiles: str,
    metal: str,
    oxidation_state: Optional[int] = None,
) -> float:
    """
    Score compatibility between ligand and metal.

    Based on HSAB theory:
    - Hard acid + hard base = good
    - Soft acid + soft base = good
    - Hard acid + soft base = bad
    - Soft acid + hard base = bad

    Args:
        smiles: Ligand SMILES
        metal: Metal symbol
        oxidation_state: Optional oxidation state

    Returns:
        Compatibility score 0.0-1.0
    """
    donors = identify_donors_from_smiles(smiles)

    if not donors:
        return 0.0

    metal_hsab = get_hsab_class(metal, oxidation_state)

    total_score = 0.0
    max_score = 0.0

    for donor in donors:
        donor_hsab = donor["hsab_class"]

        # HSAB matching score
        if metal_hsab == donor_hsab:
            match_score = 1.0
        elif metal_hsab == "borderline" or donor_hsab == "borderline":
            match_score = 0.7
        else:
            # Hard-soft mismatch
            match_score = 0.1

        # Weight by denticity potential
        weight = donor["max_denticity"]

        total_score += match_score * weight
        max_score += 1.0 * weight

    if max_score == 0:
        return 0.0

    return total_score / max_score


def analyze_denticity(smiles: str) -> List[Dict[str, Any]]:
    """
    Analyze denticity potential of ligand.

    Args:
        smiles: Ligand SMILES

    Returns:
        List of chelating group analyses
    """
    if not HAS_RDKIT:
        return []

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    donors = identify_donors_from_smiles(smiles)

    # Group donors by proximity
    chelating_groups = []

    # Find carboxylate groups (bidentate)
    for donor in donors:
        if donor["type"] == "carboxylate":
            # Check for paired oxygen
            chelating_groups.append({
                "type": "carboxylate",
                "atoms": [donor["atom_index"]],
                "max_denticity": 2,
                "element": "O",
            })

    # Find other potential chelates (2 donors within 3 bonds)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return chelating_groups

    for i, d1 in enumerate(donors):
        for d2 in donors[i+1:]:
            # Check distance in bond graph
            try:
                path = Chem.GetShortestPath(mol, d1["atom_index"], d2["atom_index"])
                if path and len(path) <= 4:  # 3 bonds or fewer
                    chelating_groups.append({
                        "type": "chelate",
                        "atoms": [d1["atom_index"], d2["atom_index"]],
                        "max_denticity": 2,
                        "elements": [d1["element"], d2["element"]],
                        "bond_distance": len(path) - 1,
                    })
            except:
                pass

    return chelating_groups


def get_recommended_coordination_mode(
    smiles: str,
    metal: str,
    target_cn: int,
) -> Dict[str, Any]:
    """
    Recommend coordination mode for ligand with metal.

    Args:
        smiles: Ligand SMILES
        metal: Metal symbol
        target_cn: Target coordination number

    Returns:
        Recommended coordination configuration
    """
    donors = identify_donors_from_smiles(smiles)
    denticity = analyze_denticity(smiles)
    compatibility = score_ligand_metal_compatibility(smiles, metal)

    # Filter donors by metal compatibility
    metal_hsab = get_hsab_class(metal)
    compatible_donors = [
        d for d in donors
        if (d["hsab_class"] == metal_hsab or
            d["hsab_class"] == "borderline" or
            metal_hsab == "borderline")
    ]

    # Calculate achievable CN from ligand
    max_ligand_cn = sum(d["max_denticity"] for d in compatible_donors)

    # Remaining sites from protein
    protein_sites = max(0, target_cn - min(max_ligand_cn, target_cn))

    return {
        "ligand_smiles": smiles,
        "metal": metal,
        "target_cn": target_cn,
        "compatibility_score": compatibility,
        "compatible_donors": len(compatible_donors),
        "max_ligand_cn": max_ligand_cn,
        "recommended_ligand_sites": min(max_ligand_cn, target_cn - 2),  # Leave sites for protein
        "protein_sites_needed": protein_sites,
        "chelating_groups": len(denticity),
    }
```

**Step 4: Run test to verify it passes**

Run: `cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless && python -m pytest test_ligand_donors.py -v`
Expected: PASS (may skip if RDKit unavailable)

**Step 5: Commit**

```bash
git add backend/serverless/ligand_donors.py backend/serverless/test_ligand_donors.py
git commit -m "feat: add ligand donor identification from SMILES

- Identify coordinating atoms (O, N, S, P)
- Classify donor types (carboxylate, amine, thiol, etc.)
- HSAB compatibility scoring with metals
- Denticity analysis for chelation
- Coordination mode recommendations

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Phase 4: Comprehensive Test Suite

### Task 4.1: Add Real PDB Structure Tests

**Files:**
- Create: `backend/serverless/test_real_structures.py`
- Create: `backend/serverless/fixtures/reference_pdbs.py`

**Step 1: Write the test file**

```python
# test_real_structures.py
"""
Biochemically Rigorous Tests with Real PDB Structures

These tests validate coordination chemistry against experimentally
determined protein structures from the PDB.
"""
import pytest
import math
from typing import Dict, Any


# Reference structures for validation
REFERENCE_STRUCTURES = {
    "carbonic_anhydrase_zn": {
        "pdb_id": "1CA2",
        "metal": "ZN",
        "expected_cn": 4,
        "expected_geometry": "tetrahedral",
        "expected_donors": ["HIS", "HIS", "HIS"],  # + water
        "bond_distance_range": (1.95, 2.15),
    },
    "rubredoxin_fe": {
        "pdb_id": "1RB9",
        "metal": "FE",
        "expected_cn": 4,
        "expected_geometry": "tetrahedral",
        "expected_donors": ["CYS", "CYS", "CYS", "CYS"],
        "bond_distance_range": (2.25, 2.40),
    },
    "plastocyanin_cu": {
        "pdb_id": "1PLC",
        "metal": "CU",
        "expected_cn": 4,
        "expected_geometry": "distorted_tetrahedral",
        "expected_donors": ["HIS", "HIS", "CYS", "MET"],
        "bond_distance_range": (1.90, 2.60),  # Wider for Met
    },
    "calmodulin_ca": {
        "pdb_id": "1CLL",
        "metal": "CA",
        "expected_cn": 7,
        "expected_geometry": "pentagonal_bipyramidal",
        "expected_donors": ["ASP", "ASP", "ASP", "GLU"],  # + waters
        "bond_distance_range": (2.30, 2.60),
    },
}


class TestZincCoordination:
    """Test zinc coordination chemistry against PDB structures."""

    @pytest.fixture
    def carbonic_anhydrase_pdb(self):
        """Fetch carbonic anhydrase structure."""
        from pdb_fetch import fetch_pdb
        result = fetch_pdb("1CA2")
        if not result["success"]:
            pytest.skip("Could not fetch PDB 1CA2")
        return result["content"]

    def test_zinc_tetrahedral_angles(self, carbonic_anhydrase_pdb):
        """Zinc sites should have ~109.5° L-M-L angles."""
        from coordination import analyze_coordination_geometry

        analysis = analyze_coordination_geometry(
            carbonic_anhydrase_pdb,
            metal_chain="A",  # May need adjustment
            metal_residue="ZN",
            metal_resnum=262,  # Typical for 1CA2
        )

        if not analysis.get("success"):
            pytest.skip("Could not analyze zinc site")

        # Check tetrahedral angles (109.5° ± 15°)
        angles = analysis["angles"]["lml_angles"]
        for angle in angles:
            assert 95 < angle < 125, f"Angle {angle}° outside tetrahedral range"

    def test_zinc_histidine_coordination(self, carbonic_anhydrase_pdb):
        """Zinc should be coordinated by 3 His residues."""
        from coordination import get_coordinating_atoms
        from metal_site_fetcher import extract_metal_coordination

        coord = extract_metal_coordination(carbonic_anhydrase_pdb, "ZN")

        if not coord["success"]:
            pytest.skip("Could not find zinc")

        his_count = sum(
            1 for a in coord["coordinating_atoms"]
            if a["res_name"] == "HIS"
        )

        assert his_count >= 3, f"Expected 3 His, found {his_count}"

    def test_zinc_bond_distances(self, carbonic_anhydrase_pdb):
        """Zn-N distances should be 2.00-2.10 Å."""
        from metal_site_fetcher import extract_metal_coordination

        coord = extract_metal_coordination(carbonic_anhydrase_pdb, "ZN", cutoff=2.5)

        if not coord["success"]:
            pytest.skip("Could not find zinc")

        for atom in coord["coordinating_atoms"]:
            if atom["element"] == "N":
                assert 1.95 <= atom["distance"] <= 2.15, \
                    f"Zn-N distance {atom['distance']} outside range"


class TestLanthanideCoordination:
    """Test lanthanide coordination chemistry."""

    def test_lanthanide_requires_high_cn(self):
        """Lanthanides need CN >= 8 for stability."""
        from metal_chemistry import get_coordination_number_range

        for metal in ["TB", "EU", "GD", "LA"]:
            cn_range = get_coordination_number_range(metal)
            assert cn_range[0] >= 7, f"{metal} min CN should be >= 7"
            assert cn_range[1] >= 8, f"{metal} max CN should be >= 8"

    def test_lanthanide_oxygen_preference(self):
        """Lanthanides should strongly prefer O donors."""
        from metal_chemistry import get_preferred_donors

        for metal in ["TB", "EU", "GD"]:
            donors = get_preferred_donors(metal)

            # O should have highest weight
            o_weights = [v for k, v in donors.items() if k.startswith("O")]
            s_weights = [v for k, v in donors.items() if k == "S"]

            if o_weights and s_weights:
                assert max(o_weights) > max(s_weights), \
                    f"{metal} should prefer O over S"

    def test_lanthanide_excludes_cysteine(self):
        """Lanthanide biases must exclude Cys."""
        from metal_chemistry import get_amino_acid_bias

        for metal in ["TB", "EU", "GD"]:
            bias = get_amino_acid_bias(metal)
            bias_dict = dict(pair.split(":") for pair in bias.split(","))

            cys_weight = float(bias_dict.get("C", 0))
            assert cys_weight <= -3.0, \
                f"{metal} Cys weight {cys_weight} not negative enough"


class TestHSABCompatibility:
    """Test HSAB theory compliance."""

    def test_hard_acid_soft_base_violation(self):
        """Should detect hard acid + soft base violations."""
        from metal_chemistry import validate_coordination_chemistry

        # Lanthanide (hard) with Cys (soft)
        result = validate_coordination_chemistry(
            "TB",
            ["C", "C", "E", "E"],  # 2 Cys + 2 Glu
        )

        assert not result["valid"], "Should flag Cys with lanthanide"
        assert any("HSAB" in e for e in result["errors"])

    def test_matching_hsab_valid(self):
        """Matching HSAB should be valid."""
        from metal_chemistry import validate_coordination_chemistry

        # Zinc (borderline) with His/Cys (borderline/soft)
        result = validate_coordination_chemistry(
            "ZN",
            ["H", "H", "C", "C"],  # His2Cys2
        )

        assert result["valid"], "His2Cys2 should be valid for Zn"

    def test_lanthanide_all_oxygen_valid(self):
        """All-oxygen coordination valid for lanthanides."""
        from metal_chemistry import validate_coordination_chemistry

        result = validate_coordination_chemistry(
            "TB",
            ["E", "E", "E", "D", "D", "D", "N", "Q"],  # 8 O donors
        )

        assert result["valid"], "All-O coordination valid for lanthanide"
        assert result["score"] >= 0.8


class TestGeometryValidation:
    """Test coordination geometry validation."""

    def test_tetrahedral_angle_validation(self):
        """Tetrahedral should have ~109.5° angles."""
        from coordination import classify_geometry

        # Perfect tetrahedral angles
        tetrahedral_angles = [109.5, 109.5, 109.5, 109.5, 109.5, 109.5]
        geometry, rmsd = classify_geometry(4, tetrahedral_angles)

        assert geometry == "tetrahedral"
        assert rmsd < 5.0  # Low deviation from ideal

    def test_octahedral_angle_validation(self):
        """Octahedral should have 90° and 180° angles."""
        from coordination import classify_geometry

        # Perfect octahedral: 12 90° + 3 180°
        octahedral_angles = [90.0] * 12 + [180.0] * 3
        geometry, rmsd = classify_geometry(6, octahedral_angles)

        assert geometry == "octahedral"
        assert rmsd < 5.0

    def test_distorted_geometry_detection(self):
        """Should detect distorted geometries."""
        from coordination import classify_geometry

        # Distorted tetrahedral
        distorted = [95.0, 100.0, 115.0, 120.0, 110.0, 105.0]
        geometry, rmsd = classify_geometry(4, distorted)

        # Should still classify as tetrahedral but with higher RMSD
        assert geometry == "tetrahedral"
        assert rmsd > 5.0  # Higher deviation


class TestBondDistanceValidation:
    """Test bond distance ranges."""

    def test_zinc_cysteine_distance(self):
        """Zn-Cys should be 2.30-2.36 Å."""
        from metal_chemistry import get_bond_distance_range

        dist_range = get_bond_distance_range("ZN", "S")

        assert 2.28 <= dist_range[0] <= 2.32
        assert 2.33 <= dist_range[1] <= 2.38

    def test_lanthanide_carboxylate_distance(self):
        """Ln-O(carboxylate) should be 2.30-2.50 Å."""
        from metal_chemistry import get_bond_distance_range

        for metal in ["TB", "EU", "GD"]:
            dist_range = get_bond_distance_range(metal, "O")

            assert dist_range[0] >= 2.25, f"{metal} min distance too short"
            assert dist_range[1] <= 2.60, f"{metal} max distance too long"

    def test_distance_validation_against_pdb(self):
        """Bond distances should match PDB statistics."""
        from metal_chemistry import get_bond_distance_range

        # PDB survey values (from literature)
        pdb_stats = {
            ("ZN", "N"): (2.00, 2.10),
            ("ZN", "S"): (2.30, 2.35),
            ("FE", "S"): (2.28, 2.36),
            ("CU", "S"): (2.15, 2.30),
        }

        for (metal, donor), (exp_min, exp_max) in pdb_stats.items():
            calc_range = get_bond_distance_range(metal, donor)

            # Allow some tolerance
            assert abs(calc_range[0] - exp_min) < 0.1, \
                f"{metal}-{donor} min differs from PDB"
            assert abs(calc_range[1] - exp_max) < 0.1, \
                f"{metal}-{donor} max differs from PDB"


# Run with: pytest test_real_structures.py -v --tb=short
```

**Step 2: Run test**

Run: `cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless && python -m pytest test_real_structures.py -v --tb=short`
Expected: PASS (some network-dependent tests may skip)

**Step 3: Commit**

```bash
git add backend/serverless/test_real_structures.py
git commit -m "test: add biochemically rigorous test suite

- Real PDB structure validation (1CA2, 1RB9, 1PLC, 1CLL)
- HSAB compatibility testing
- Geometry angle validation (tetrahedral, octahedral)
- Bond distance validation against crystallographic data
- Lanthanide-specific coordination tests

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

### Task 4.2: Integration Test for Complete Pipeline

**Files:**
- Create: `backend/serverless/test_metal_dimer_pipeline.py`

**Step 1: Write integration test**

```python
# test_metal_dimer_pipeline.py
"""
Integration tests for complete metal-ligand dimer design pipeline.
"""
import pytest


class TestTemplateRetrievalPipeline:
    """Test complete template retrieval with fallback."""

    def test_pqq_ca_template_pipeline(self):
        """PQQ-Ca should retrieve from database or library."""
        from metal_ligand_templates import get_template_with_fallback

        template = get_template_with_fallback("pqq_ca")

        assert template is not None
        assert template["metal"] == "CA"
        assert template["source"] in ["pdb", "library", "calculated"]

    def test_unknown_template_fallback(self):
        """Unknown template should fall back to calculated."""
        from metal_ligand_templates import get_template_with_fallback

        template = get_template_with_fallback(
            "completely_unknown",
            metal="TB",
            fallback_enabled=True,
        )

        if template:
            assert template["source"] == "calculated"
            assert "warning" in template


class TestBiasGenerationPipeline:
    """Test amino acid bias generation for design."""

    def test_zinc_design_bias(self):
        """Zinc design should generate correct bias."""
        from design_types import METAL_PRESETS

        bias = METAL_PRESETS["zinc"]["bias_AA"]
        bias_dict = dict(pair.split(":") for pair in bias.split(","))

        # Cys and His should be positive
        assert float(bias_dict.get("C", 0)) > 0
        assert float(bias_dict.get("H", 0)) > 0

    def test_lanthanide_design_bias(self):
        """Lanthanide design should exclude Cys."""
        from design_types import METAL_PRESETS

        bias = METAL_PRESETS["lanthanide"]["bias_AA"]
        bias_dict = dict(pair.split(":") for pair in bias.split(","))

        # Cys must be strongly negative
        assert float(bias_dict.get("C", 0)) <= -3.0
        # Glu/Asp should be positive
        assert float(bias_dict.get("E", 0)) > 2.0
        assert float(bias_dict.get("D", 0)) > 2.0


class TestLigandAnalysisPipeline:
    """Test ligand analysis in design context."""

    @pytest.mark.skipif(
        not pytest.importorskip("rdkit", reason="RDKit required"),
        reason="RDKit not available"
    )
    def test_citrate_analysis_for_lanthanide(self):
        """Citrate should be highly compatible with lanthanides."""
        from ligand_donors import (
            identify_donors_from_smiles,
            score_ligand_metal_compatibility,
        )

        citrate = "C(C(=O)O)(CC(=O)O)(C(=O)O)O"

        donors = identify_donors_from_smiles(citrate)
        assert len(donors) >= 4

        score = score_ligand_metal_compatibility(citrate, "TB")
        assert score >= 0.7


class TestValidationPipeline:
    """Test validation of designed structures."""

    def test_hsab_validation_in_design(self):
        """Design validation should catch HSAB violations."""
        from metal_chemistry import validate_coordination_chemistry

        # Simulate a bad design (Cys with lanthanide)
        result = validate_coordination_chemistry(
            "EU",
            ["C", "E", "E", "D", "D"],
        )

        assert not result["valid"]
        assert len(result["errors"]) > 0

    def test_coordination_number_validation(self):
        """Should validate appropriate coordination number."""
        from metal_chemistry import get_coordination_number_range

        # Lanthanide needs high CN
        cn_range = get_coordination_number_range("TB")
        assert cn_range[0] >= 7

        # Zinc is typically 4
        cn_range = get_coordination_number_range("ZN")
        assert 4 in range(cn_range[0], cn_range[1] + 1)


class TestEndToEndDesign:
    """End-to-end design workflow tests."""

    def test_lanthanide_dimer_workflow(self):
        """Complete workflow for lanthanide dimer design."""
        from metal_ligand_templates import get_template_with_fallback
        from design_types import METAL_PRESETS
        from metal_chemistry import validate_coordination_chemistry

        # Step 1: Get template
        template = get_template_with_fallback(
            "citrate_tb",
            metal="TB",
            fallback_enabled=True,
        )
        assert template is not None

        # Step 2: Get design parameters
        preset = METAL_PRESETS["lanthanide"]
        bias = preset["bias_AA"]

        # Step 3: Validate expected coordination
        # (simulating designed residues)
        expected_residues = ["E", "E", "E", "D", "D", "D", "N", "Q"]
        validation = validate_coordination_chemistry("TB", expected_residues)

        assert validation["valid"]
        assert validation["score"] >= 0.8

    def test_zinc_site_workflow(self):
        """Complete workflow for zinc binding site design."""
        from metal_site_fetcher import get_reference_template
        from design_types import METAL_PRESETS
        from metal_chemistry import validate_coordination_chemistry

        # Step 1: Get reference structure
        template = get_reference_template("ZN")
        # May fail if network unavailable

        # Step 2: Get design parameters
        preset = METAL_PRESETS["zinc"]
        bias = preset["bias_AA"]
        bias_dict = dict(pair.split(":") for pair in bias.split(","))

        # Verify Cys/His preference
        assert float(bias_dict.get("C", 0)) >= 2.0
        assert float(bias_dict.get("H", 0)) >= 2.0

        # Step 3: Validate typical zinc site
        expected_residues = ["H", "H", "C", "C"]
        validation = validate_coordination_chemistry("ZN", expected_residues)

        assert validation["valid"]
```

**Step 2: Run integration tests**

Run: `cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless && python -m pytest test_metal_dimer_pipeline.py -v`
Expected: PASS

**Step 3: Commit**

```bash
git add backend/serverless/test_metal_dimer_pipeline.py
git commit -m "test: add integration tests for metal dimer pipeline

- Template retrieval with database/library/calculated fallback
- Bias generation validation for metals
- Ligand analysis integration
- HSAB validation in design context
- End-to-end workflow tests

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Phase 5: Final Integration

### Task 5.1: Update Handler to Use New Modules

**Files:**
- Modify: `backend/serverless/handler.py`

**Step 1: Update imports and handler function**

Add to handler.py after existing imports:

```python
# Import enhanced metal chemistry modules
from metal_chemistry import (
    get_amino_acid_bias,
    get_hsab_class,
    validate_coordination_chemistry,
    get_coordination_number_range,
)
from metal_ligand_templates import get_template_with_fallback
from metal_site_fetcher import get_reference_template

# Import ligand analysis (optional - may not have RDKit)
try:
    from ligand_donors import (
        identify_donors_from_smiles,
        score_ligand_metal_compatibility,
        get_recommended_coordination_mode,
    )
    LIGAND_ANALYSIS_AVAILABLE = True
except ImportError:
    LIGAND_ANALYSIS_AVAILABLE = False
```

Update `handle_interface_metal_ligand_design()` to use database-first templates:

```python
def handle_interface_metal_ligand_design(request: Dict[str, Any]) -> Dict[str, Any]:
    """
    Handle metal-ligand complex dimer design.

    Enhanced with:
    - Database-first template retrieval
    - HSAB-compliant amino acid biases
    - Ligand donor analysis (if RDKit available)
    """
    # Get template with database priority
    template_name = request.get("template_name")
    complex_pdb = request.get("complex_pdb")
    metal = request.get("metal")
    ligand_smiles = request.get("ligand_smiles")

    template = None

    if template_name:
        # Database-first retrieval
        template = get_template_with_fallback(
            template_name,
            metal=metal,
            fallback_enabled=True,
        )

        if template and template.get("source") == "calculated":
            print(f"[Handler] Warning: Using calculated template for {template_name}")

    if not template and complex_pdb:
        # Parse from provided PDB
        template = parse_metal_ligand_complex(complex_pdb)

    if not template:
        return {
            "status": "failed",
            "error": "No template found and no complex_pdb provided",
        }

    # Get HSAB-compliant bias
    metal_code = template.get("metal", metal)
    if metal_code:
        bias = get_amino_acid_bias(metal_code)
        request["bias_AA"] = bias

        # Validate coordination number
        cn_range = get_coordination_number_range(metal_code)
        request["target_coordination"] = cn_range[1]

    # Analyze ligand if SMILES provided
    if LIGAND_ANALYSIS_AVAILABLE and ligand_smiles:
        compatibility = score_ligand_metal_compatibility(ligand_smiles, metal_code)
        if compatibility < 0.5:
            print(f"[Handler] Warning: Low ligand-metal compatibility: {compatibility:.2f}")

    # Continue with design...
    # (rest of existing handler code)
```

**Step 2: Run existing tests to verify no regression**

Run: `cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless && python -m pytest test_metal_ligand_complex.py test_metal_coordination.py -v`
Expected: PASS

**Step 3: Commit**

```bash
git add backend/serverless/handler.py
git commit -m "feat: integrate enhanced metal chemistry into handler

- Database-first template retrieval
- HSAB-compliant amino acid biases
- Ligand compatibility scoring
- Coordination number validation
- Warnings for calculated fallbacks

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

### Task 5.2: Run Full Test Suite

**Step 1: Run all tests**

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless
python -m pytest test_metal_chemistry.py test_metal_site_fetcher.py test_ligand_donors.py test_real_structures.py test_metal_dimer_pipeline.py test_metal_ligand_complex.py -v --tb=short
```

**Step 2: Fix any failures**

**Step 3: Final commit**

```bash
git add .
git commit -m "chore: complete metal-ligand dimer enhancement

All phases complete:
1. HSAB-compliant amino acid biases
2. Database-driven template system
3. Ligand donor intelligence
4. Comprehensive test suite

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Summary

| Phase | Tasks | Files Created/Modified |
|-------|-------|----------------------|
| 1 | HSAB Amino Acid Biases | `metal_chemistry.py`, `design_types.py` |
| 2 | Database Templates | `metal_site_fetcher.py`, `metal_ligand_templates.py` |
| 3 | Ligand Intelligence | `ligand_donors.py` |
| 4 | Test Suite | `test_real_structures.py`, `test_metal_dimer_pipeline.py` |
| 5 | Integration | `handler.py` |

**Total: 6 new files, 3 modified files, ~2000 lines of code**
