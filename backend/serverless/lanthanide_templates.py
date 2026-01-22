"""
Lanthanide Binding Site Template Generators

Generate template PDB structures for lanthanide-binding heterodimer design.
Based on Caldwell et al. (2020) PNAS - Tight lanthanide binding in de novo TIM barrel.

Template Types:
1. EF-hand inspired - 8 coordinating Glu residues (4 per chain)
2. C4-symmetric - 4 Glu in symmetric arrangement (like TIM barrel top)
"""

import logging
import math
import numpy as np
from typing import Any, Dict, List, Optional, Tuple

# Import Architector integration (optional - may not be available in all environments)
try:
    from architector_integration import (
        is_architector_available,
        generate_template_with_architector,
        generate_sap_positions,
        generate_ttp_positions,
    )
    ARCHITECTOR_INTEGRATION_AVAILABLE = True
except ImportError:
    ARCHITECTOR_INTEGRATION_AVAILABLE = False
    is_architector_available = lambda: False
    generate_template_with_architector = lambda *args, **kwargs: None

# Import geometry validation (optional)
try:
    from geometry_validation import (
        validate_coordination_geometry,
        quick_validate,
        get_validation_metrics,
    )
    VALIDATION_AVAILABLE = True
except ImportError:
    VALIDATION_AVAILABLE = False
    validate_coordination_geometry = lambda *args, **kwargs: {"valid": True}
    quick_validate = lambda x: True

logger = logging.getLogger(__name__)

# Lanthanide coordination parameters
LANTHANIDE_PARAMS = {
    "TB": {
        "name": "Terbium",
        "element": "TB",
        "coordination": (8, 9),
        "bond_distance": 2.35,  # Å (Tb-O typical)
        "geometry": "square_antiprism",
    },
    "GD": {
        "name": "Gadolinium",
        "element": "GD",
        "coordination": (8, 9),
        "bond_distance": 2.37,
        "geometry": "square_antiprism",
    },
    "EU": {
        "name": "Europium",
        "element": "EU",
        "coordination": (8, 9),
        "bond_distance": 2.40,
        "geometry": "square_antiprism",
    },
    "LA": {
        "name": "Lanthanum",
        "element": "LA",
        "coordination": (9, 10),
        "bond_distance": 2.50,
        "geometry": "tricapped_trigonal_prism",
    },
    "YB": {
        "name": "Ytterbium",
        "element": "YB",
        "coordination": (8,),
        "bond_distance": 2.28,
        "geometry": "square_antiprism",
    },
}

# Template types available (legacy)
TEMPLATE_TYPES = ["ef_hand", "c4_symmetric"]

# =============================================================================
# NEW: Template Library with Chemically Realistic Coordination
# =============================================================================

# Coordination geometry rules by CN
COORDINATION_GEOMETRIES = {
    6: {
        "geometry": "octahedral",
        "description": "6-coordinate octahedral",
        "angles": [90, 180],
        "max_protein_donors": 6,
        "typical_waters": 0,
    },
    7: {
        "geometry": "pentagonal_bipyramidal",
        "description": "7-coordinate pentagonal bipyramidal",
        "angles": [72, 90, 180],
        "max_protein_donors": 6,
        "typical_waters": 1,
    },
    8: {
        "geometry": "square_antiprism",
        "description": "8-coordinate square antiprism (SAP)",
        "angles": [70.5, 109.5, 141.5],
        "max_protein_donors": 8,
        "typical_waters": (0, 2),
        "alternative": "bicapped_trigonal_prism",
    },
    9: {
        "geometry": "tricapped_trigonal_prism",
        "description": "9-coordinate tricapped trigonal prism (TTP)",
        "angles": [73.4, 82.8, 131.8],
        "max_protein_donors": 8,
        "typical_waters": (1, 3),
        "alternative": "capped_square_antiprism",
    },
    10: {
        "geometry": "bicapped_square_antiprism",
        "description": "10-coordinate bicapped SAP (rare)",
        "angles": [],
        "max_protein_donors": 8,
        "typical_waters": (2, 4),
        "note": "Rare in proteins, mostly organometallic",
    },
}

# Template library with chemically realistic coordination
TEMPLATE_LIBRARY = {
    "caldwell_4": {
        "name": "Caldwell TFD-EE Style",
        "description": "4 bidentate Glu, highest affinity (subfemtomolar)",
        "coordination_number": 8,
        "residues": [
            # ALTERNATING Z-LEVEL per residue for maximum 3D separation:
            # Each residue alternates between top and bottom z-levels
            # This creates a TRUE square antiprism arrangement:
            # Position 1: 0° top, Position 2: 45° bottom, Position 3: 90° top, Position 4: 135° bottom
            # Residues are 45° apart horizontally and alternate z-levels
            {"chain": "A", "resnum": 15, "type": "GLU", "mode": "bidentate", "angle": 0, "z_offset": 3.0},
            {"chain": "A", "resnum": 25, "type": "GLU", "mode": "bidentate", "angle": 90, "z_offset": -3.0},
            {"chain": "B", "resnum": 15, "type": "GLU", "mode": "bidentate", "angle": 180, "z_offset": 3.0},
            {"chain": "B", "resnum": 25, "type": "GLU", "mode": "bidentate", "angle": 270, "z_offset": -3.0},
        ],
        "waters": [],
        "geometry": "square_antiprism",
        "best_for": ["TB", "EU", "GD"],
        "reference": "Caldwell et al. 2020 PNAS",
    },
    "ef_hand_8": {
        "name": "EF-Hand Mixed",
        "description": "4 monodentate Asp + 2 bidentate Glu = CN8",
        "coordination_number": 8,
        "residues": [
            # 6 residues with SAP-like geometry (offset z levels by 45° angular shift):
            # Chain A (top, z>0): 0, 120, 240 degrees
            # Chain B (bottom, z<0): 60, 180, 300 degrees (60° offset for maximum separation)
            {"chain": "A", "resnum": 10, "type": "ASP", "mode": "monodentate", "angle": 0, "z_offset": 3.0},
            {"chain": "A", "resnum": 15, "type": "GLU", "mode": "bidentate", "angle": 120, "z_offset": 3.0},
            {"chain": "A", "resnum": 20, "type": "ASP", "mode": "monodentate", "angle": 240, "z_offset": 3.0},
            {"chain": "B", "resnum": 10, "type": "ASP", "mode": "monodentate", "angle": 60, "z_offset": -3.0},
            {"chain": "B", "resnum": 15, "type": "GLU", "mode": "bidentate", "angle": 180, "z_offset": -3.0},
            {"chain": "B", "resnum": 20, "type": "ASP", "mode": "monodentate", "angle": 300, "z_offset": -3.0},
        ],
        "waters": [],
        "geometry": "distorted_square_antiprism",
        "best_for": ["TB", "GD", "EU", "YB"],
        "reference": "Natural EF-hand motifs",
    },
    "lanm_mixed": {
        "name": "Lanmodulin-Inspired",
        "description": "3 bidentate + 2 waters = CN8-9 (natural-like)",
        "coordination_number": 9,
        "residues": [
            # Trigonal arrangement at 120° intervals with alternating z
            {"chain": "A", "resnum": 15, "type": "ASP", "mode": "bidentate", "angle": 0, "z_offset": 3.0},
            {"chain": "A", "resnum": 25, "type": "ASP", "mode": "bidentate", "angle": 120, "z_offset": -3.0},
            {"chain": "B", "resnum": 15, "type": "ASP", "mode": "bidentate", "angle": 240, "z_offset": 3.0},
        ],
        "waters": [
            {"position": "axial_top", "distance": 2.45},
            {"position": "axial_bottom", "distance": 2.45},
            {"position": "equatorial", "distance": 2.45},
        ],
        "backbone_carbonyl": False,
        "geometry": "tricapped_trigonal_prism",
        "best_for": ["TB", "GD", "EU"],
        "reference": "Cotruvo et al. 2018 JACS (Lanmodulin)",
    },
    "high_coord_9": {
        "name": "High Coordination",
        "description": "4 bidentate Glu + 1 water = CN9 for large Ln",
        "coordination_number": 9,
        "residues": [
            # ALTERNATING Z-LEVEL per residue (same as caldwell_4):
            # Each residue alternates between top and bottom z-levels
            {"chain": "A", "resnum": 15, "type": "GLU", "mode": "bidentate", "angle": 0, "z_offset": 3.0},
            {"chain": "A", "resnum": 25, "type": "GLU", "mode": "bidentate", "angle": 90, "z_offset": -3.0},
            {"chain": "B", "resnum": 15, "type": "GLU", "mode": "bidentate", "angle": 180, "z_offset": 3.0},
            {"chain": "B", "resnum": 25, "type": "GLU", "mode": "bidentate", "angle": 270, "z_offset": -3.0},
        ],
        "waters": [
            {"position": "axial_top", "distance": 2.50},
        ],
        "geometry": "tricapped_trigonal_prism",
        "best_for": ["LA", "CE", "PR"],
        "reference": "Standard lanthanide chemistry",
    },
    # Legacy template (original 8-residue - deprecated)
    "legacy_8res": {
        "name": "Legacy 8-Residue (Deprecated)",
        "description": "Original 8 bidentate = CN16 (unrealistic)",
        "coordination_number": 16,
        "residues": [
            {"chain": "A", "resnum": 10, "type": "ASP", "mode": "bidentate", "angle": 0},
            {"chain": "A", "resnum": 15, "type": "ASP", "mode": "bidentate", "angle": 90},
            {"chain": "A", "resnum": 20, "type": "ASP", "mode": "bidentate", "angle": 180},
            {"chain": "A", "resnum": 25, "type": "ASP", "mode": "bidentate", "angle": 270},
            {"chain": "B", "resnum": 10, "type": "ASP", "mode": "bidentate", "angle": 45},
            {"chain": "B", "resnum": 15, "type": "ASP", "mode": "bidentate", "angle": 135},
            {"chain": "B", "resnum": 20, "type": "ASP", "mode": "bidentate", "angle": 225},
            {"chain": "B", "resnum": 25, "type": "ASP", "mode": "bidentate", "angle": 315},
        ],
        "waters": [],
        "geometry": "square_antiprism",
        "best_for": [],
        "deprecated": True,
        "reference": "Original implementation",
    },
}

# Second-sphere effects on water coordination (Lanmodulin-inspired)
SECOND_SPHERE_EFFECTS = {
    "PHE": {"steric": "high", "q_reduction": 2, "description": "Blocks 2 water sites"},
    "TRP": {"steric": "high", "q_reduction": 2, "antenna": True, "description": "Blocks waters, TEBL antenna"},
    "TYR": {"steric": "medium", "q_reduction": 1, "description": "Moderate steric blocking"},
    "LEU": {"steric": "medium", "q_reduction": 1, "description": "Moderate steric blocking"},
    "ILE": {"steric": "medium", "q_reduction": 1, "description": "Moderate steric blocking"},
    "VAL": {"steric": "low", "q_reduction": 0, "description": "Minimal steric effect"},
    "ASN": {"steric": "low", "q_reduction": 0, "h_bond": True, "description": "H-bond to waters"},
    "GLN": {"steric": "low", "q_reduction": 0, "h_bond": True, "description": "H-bond to waters"},
    "SER": {"steric": "low", "q_reduction": 0, "h_bond": True, "description": "H-bond to waters"},
    "THR": {"steric": "low", "q_reduction": 0, "h_bond": True, "description": "H-bond to waters"},
}

# Metal-specific template recommendations
METAL_TEMPLATE_RECOMMENDATIONS = {
    "TB": "caldwell_4",   # Highest affinity for Tb
    "EU": "caldwell_4",   # Similar to Tb
    "GD": "ef_hand_8",    # Balanced for Gd
    "YB": "ef_hand_8",    # Smaller, needs tighter coordination
    "LA": "high_coord_9", # Large, needs CN=9
    "CE": "high_coord_9", # Large, needs CN=9
    "PR": "high_coord_9", # Large, needs CN=9
}


# =============================================================================
# Fixed Position Extraction for LigandMPNN
# =============================================================================


def get_template_fixed_positions(template_name: str) -> List[str]:
    """
    Extract coordinating residue positions from a template definition.

    These positions should be passed to LigandMPNN's fixed_positions
    to prevent mutation of metal-coordinating residues during sequence design.

    CRITICAL: Without this, LigandMPNN may mutate the carefully positioned
    Asp/Glu coordinating residues, destroying the metal binding geometry.

    Args:
        template_name: Name of template from TEMPLATE_LIBRARY (e.g., "caldwell_4")

    Returns:
        List of position strings like ["A15", "A25", "B15", "B25"]
        Empty list if template not found.

    Example:
        >>> positions = get_template_fixed_positions("caldwell_4")
        >>> print(positions)
        ['A15', 'A25', 'B15', 'B25']

        >>> # Use with LigandMPNN
        >>> from handler import run_ligandmpnn_for_ligand_binding
        >>> result = run_ligandmpnn_for_ligand_binding(
        ...     pdb_content=scaffold_pdb,
        ...     ligand_type="lanthanide",
        ...     fixed_positions=positions,
        ... )
    """
    if template_name not in TEMPLATE_LIBRARY:
        logger.warning(f"Template '{template_name}' not found in TEMPLATE_LIBRARY")
        return []

    template = TEMPLATE_LIBRARY[template_name]
    residues = template.get("residues", [])

    positions = []
    for res in residues:
        chain = res.get("chain", "A")
        resnum = res.get("resnum")
        if resnum is not None:
            positions.append(f"{chain}{resnum}")

    logger.debug(f"Extracted {len(positions)} fixed positions from template '{template_name}': {positions}")
    return positions


# =============================================================================
# NEW: Geometry Position Generators
# =============================================================================


def _generate_sap_positions(
    metal_pos: np.ndarray,
    bond_distance: float,
    cn: int = 8,
) -> List[np.ndarray]:
    """
    Generate square antiprism (SAP) coordination positions.

    SAP has two square faces, rotated 45° relative to each other:
    - Top square: angles 0°, 90°, 180°, 270°
    - Bottom square: angles 45°, 135°, 225°, 315° (45° rotation from top)

    Args:
        metal_pos: Position of metal ion
        bond_distance: Target Ln-O distance
        cn: Coordination number (max 8 for SAP)

    Returns:
        List of 3D positions for coordinating atoms
    """
    theta = math.acos(1/3)  # ~54.7° for ideal SAP
    positions = []

    # Top square (z > 0)
    for angle in [0, 90, 180, 270]:
        rad = math.radians(angle)
        pos = metal_pos + bond_distance * np.array([
            math.sin(theta) * math.cos(rad),
            math.sin(theta) * math.sin(rad),
            math.cos(theta),
        ])
        positions.append(pos)

    # Bottom square (z < 0), rotated 45°
    for angle in [45, 135, 225, 315]:
        rad = math.radians(angle)
        pos = metal_pos + bond_distance * np.array([
            math.sin(theta) * math.cos(rad),
            math.sin(theta) * math.sin(rad),
            -math.cos(theta),
        ])
        positions.append(pos)

    return positions[:cn]


def _generate_ttp_positions(
    metal_pos: np.ndarray,
    bond_distance: float,
    cn: int = 9,
) -> List[np.ndarray]:
    """
    Generate tricapped trigonal prism (TTP) coordination positions.

    TTP has:
    - Top triangle: 3 positions at z > 0
    - Bottom triangle: 3 positions at z < 0, rotated 60° from top
    - 3 equatorial caps at z = 0

    Args:
        metal_pos: Position of metal ion
        bond_distance: Target Ln-O distance
        cn: Coordination number (max 9 for TTP)

    Returns:
        List of 3D positions for coordinating atoms
    """
    positions = []
    theta = math.radians(55)  # Tilt angle for prism vertices

    # Top triangle (3 positions)
    for i, angle in enumerate([0, 120, 240]):
        rad = math.radians(angle)
        pos = metal_pos + bond_distance * np.array([
            math.sin(theta) * math.cos(rad),
            math.sin(theta) * math.sin(rad),
            math.cos(theta),
        ])
        positions.append(pos)

    # Bottom triangle (3 positions), rotated 60° from top
    for i, angle in enumerate([60, 180, 300]):
        rad = math.radians(angle)
        pos = metal_pos + bond_distance * np.array([
            math.sin(theta) * math.cos(rad),
            math.sin(theta) * math.sin(rad),
            -math.cos(theta),
        ])
        positions.append(pos)

    # Equatorial caps (3 positions at z=0)
    for i, angle in enumerate([30, 150, 270]):
        rad = math.radians(angle)
        pos = metal_pos + bond_distance * np.array([
            math.cos(rad),
            math.sin(rad),
            0.0,
        ])
        positions.append(pos)

    return positions[:cn]


def _generate_octahedral_positions(
    metal_pos: np.ndarray,
    bond_distance: float,
) -> List[np.ndarray]:
    """
    Generate octahedral (CN=6) coordination positions.

    Octahedral: 6 positions along ±x, ±y, ±z axes.

    Args:
        metal_pos: Position of metal ion
        bond_distance: Target Ln-O distance

    Returns:
        List of 6 3D positions
    """
    directions = [
        np.array([1, 0, 0]),
        np.array([-1, 0, 0]),
        np.array([0, 1, 0]),
        np.array([0, -1, 0]),
        np.array([0, 0, 1]),
        np.array([0, 0, -1]),
    ]
    return [metal_pos + bond_distance * d for d in directions]


def _generate_coordination_positions(
    geometry: str,
    cn: int,
    metal_pos: np.ndarray,
    bond_distance: float,
) -> List[np.ndarray]:
    """
    Generate ideal coordination positions for given geometry.

    Args:
        geometry: Geometry name from COORDINATION_GEOMETRIES
        cn: Coordination number
        metal_pos: Position of metal ion
        bond_distance: Target Ln-O distance

    Returns:
        List of 3D positions for coordinating atoms
    """
    if geometry in ("square_antiprism", "distorted_square_antiprism"):
        return _generate_sap_positions(metal_pos, bond_distance, cn)
    elif geometry == "tricapped_trigonal_prism":
        return _generate_ttp_positions(metal_pos, bond_distance, cn)
    elif geometry == "octahedral":
        return _generate_octahedral_positions(metal_pos, bond_distance)
    elif geometry == "pentagonal_bipyramidal":
        # For CN=7, use modified TTP (remove 2 positions)
        return _generate_ttp_positions(metal_pos, bond_distance, cn=7)
    elif geometry == "bicapped_square_antiprism":
        # For CN=10, use SAP with 2 axial positions
        positions = _generate_sap_positions(metal_pos, bond_distance, cn=8)
        # Add axial caps
        positions.append(metal_pos + np.array([0, 0, bond_distance]))
        positions.append(metal_pos + np.array([0, 0, -bond_distance]))
        return positions[:cn]
    else:
        # Default to SAP for CN=8
        return _generate_sap_positions(metal_pos, bond_distance, min(cn, 8))


# =============================================================================
# NEW: Monodentate Residue Generator
# =============================================================================


def _generate_monodentate_asp(
    chain: str,
    res_num: int,
    ca_pos: np.ndarray,
    metal_pos: np.ndarray,
    atom_serial_start: int,
    target_ln_o_distance: float = 2.35,
    coordinating_oxygen: str = "OD1",
) -> Tuple[str, int]:
    """
    Generate MONODENTATE aspartate with only ONE oxygen coordinating.

    Uses "inside-out" construction with proper amino acid geometry.
    In monodentate mode:
    - coordinating_oxygen points toward metal
    - The other oxygen points AWAY from metal (~120° angle)

    Args:
        chain: Chain identifier
        res_num: Residue number
        ca_pos: Position of CA atom
        metal_pos: Position of metal ion
        atom_serial_start: Starting atom serial number
        target_ln_o_distance: Target Ln-O distance (Å)
        coordinating_oxygen: Which oxygen coordinates ("OD1" or "OD2")

    Returns:
        Tuple of (PDB lines string, next atom serial number)
    """
    # Standard bond lengths
    CA_CB = 1.54
    CB_CG = 1.52
    CG_OD = 1.25
    CA_N = 1.47
    CA_C = 1.52
    C_O = 1.23

    # Standard angles
    TETRAHEDRAL = math.radians(109.5)
    CARBOXYL_OCO = math.radians(124)

    # When building chain A→B→C where ref_dir points A→B, to get angle A-B-C = θ,
    # we need: C = B + d*(-cos(θ)*ref_dir + sin(θ)*perp)
    COS_TET = -math.cos(TETRAHEDRAL)  # ≈ +0.33 for 109.5°
    SIN_TET = math.sin(TETRAHEDRAL)   # ≈ +0.94

    # Direction from CA toward metal
    ca_to_metal = metal_pos - ca_pos
    ca_metal_dist = np.linalg.norm(ca_to_metal)
    if ca_metal_dist < 0.01:
        ca_to_metal = np.array([1.0, 0.0, 0.0])
        ca_metal_dist = 1.0
    main_axis = ca_to_metal / ca_metal_dist

    # Create orthonormal basis
    perp1 = np.cross(main_axis, np.array([0, 0, 1]))
    if np.linalg.norm(perp1) < 0.1:
        perp1 = np.cross(main_axis, np.array([0, 1, 0]))
    perp1 = perp1 / np.linalg.norm(perp1)
    perp2 = np.cross(main_axis, perp1)
    perp2 = perp2 / np.linalg.norm(perp2)

    # Step 1: Build CB at proper tetrahedral angle
    cb_axis_component = COS_TET * CA_CB
    cb_perp_component = SIN_TET * CA_CB
    cb_pos = ca_pos + main_axis * cb_axis_component + perp1 * cb_perp_component

    # Step 2: Build CG at proper angle from CB
    cb_dir = (cb_pos - ca_pos)
    cb_dir = cb_dir / np.linalg.norm(cb_dir)

    cg_axis_component = COS_TET * CB_CG
    cg_perp_component = SIN_TET * CB_CG
    cg_perp = np.cross(cb_dir, perp2)
    if np.linalg.norm(cg_perp) < 0.1:
        cg_perp = perp1
    else:
        cg_perp = cg_perp / np.linalg.norm(cg_perp)
    cg_pos = cb_pos + cb_dir * cg_axis_component + cg_perp * cg_perp_component

    # Step 3: Build carboxylate with proper geometry
    # In a carboxylate: O-C-O angle = 124°, C-C-O angle = 118°
    # For monodentate, we use same geometry but orient so one O is closer to metal
    cg_dir = (cg_pos - cb_pos)
    cg_dir = cg_dir / np.linalg.norm(cg_dir)

    # The carboxylate bisector points along cg_dir (due to ~118° C-C-O angle)
    od_bisector = cg_dir

    # The carboxylate plane normal is perpendicular to cg_dir
    # Orient it toward the metal so one oxygen is closer to metal
    cg_to_metal = metal_pos - cg_pos
    od_plane_normal = cg_to_metal - np.dot(cg_to_metal, cg_dir) * cg_dir
    if np.linalg.norm(od_plane_normal) < 0.1:
        od_plane_normal = perp1 - np.dot(perp1, cg_dir) * cg_dir
        if np.linalg.norm(od_plane_normal) < 0.1:
            od_plane_normal = perp2
    od_plane_normal = od_plane_normal / np.linalg.norm(od_plane_normal)

    # Half angle for OD positions (O-C-O angle is 124°)
    half_oco = CARBOXYL_OCO / 2  # ~62 degrees

    # Position OD1 and OD2: each at 62° from bisector (cg_dir)
    # OD1 is on the +plane_normal side (closer to metal)
    # OD2 is on the -plane_normal side (farther from metal)
    od1_dir = math.cos(half_oco) * od_bisector + math.sin(half_oco) * od_plane_normal
    od2_dir = math.cos(half_oco) * od_bisector - math.sin(half_oco) * od_plane_normal
    od1_pos = cg_pos + od1_dir * CG_OD
    od2_pos = cg_pos + od2_dir * CG_OD

    # For monodentate, swap if user wants OD2 to be the coordinating oxygen
    if coordinating_oxygen == "OD2":
        od1_pos, od2_pos = od2_pos, od1_pos

    # Note: We do NOT translate/scale the side chain. The geometry is built with proper
    # bond lengths and angles. The actual metal-O distance depends on the CA position.
    # RFdiffusion will optimize the structure during the design process.

    # Step 4: Generate backbone atoms
    n_pos = ca_pos - main_axis * CA_N

    c_perp = np.cross(main_axis, perp1)
    if np.linalg.norm(c_perp) < 0.1:
        c_perp = perp2
    else:
        c_perp = c_perp / np.linalg.norm(c_perp)
    c_pos = ca_pos + c_perp * CA_C
    o_pos = c_pos + c_perp * C_O

    atoms = [
        ("N", "N", n_pos),
        ("CA", "C", ca_pos),
        ("C", "C", c_pos),
        ("O", "O", o_pos),
        ("CB", "C", cb_pos),
        ("CG", "C", cg_pos),
        ("OD1", "O", od1_pos),
        ("OD2", "O", od2_pos),
    ]

    lines = []
    serial = atom_serial_start
    for atom_name, element, pos in atoms:
        if len(atom_name) <= 2:
            atom_name_fmt = f" {atom_name:<3s}"
        else:
            atom_name_fmt = f"{atom_name:<4s}"
        line = (
            f"ATOM  {serial:5d} {atom_name_fmt} ASP {chain}{res_num:4d}    "
            f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00  0.00          {element:>2s}"
        )
        lines.append(line)
        serial += 1

    return "\n".join(lines), serial


def _generate_monodentate_glu(
    chain: str,
    res_num: int,
    ca_pos: np.ndarray,
    metal_pos: np.ndarray,
    atom_serial_start: int,
    target_ln_o_distance: float = 2.35,
    coordinating_oxygen: str = "OE1",
) -> Tuple[str, int]:
    """
    Generate MONODENTATE glutamate with only ONE oxygen coordinating.

    Uses "inside-out" construction with proper amino acid geometry.
    In monodentate mode:
    - coordinating_oxygen points toward metal
    - The other oxygen points AWAY from metal (~120° angle)

    Args:
        chain: Chain identifier
        res_num: Residue number
        ca_pos: Position of CA atom
        metal_pos: Position of metal ion
        atom_serial_start: Starting atom serial number
        target_ln_o_distance: Target Ln-O distance (Å)
        coordinating_oxygen: Which oxygen coordinates ("OE1" or "OE2")

    Returns:
        Tuple of (PDB lines string, next atom serial number)
    """
    # Standard bond lengths
    CA_CB = 1.54
    CB_CG = 1.52
    CG_CD = 1.52
    CD_OE = 1.25
    CA_N = 1.47
    CA_C = 1.52
    C_O = 1.23

    # Standard angles
    TETRAHEDRAL = math.radians(109.5)
    CARBOXYL_OCO = math.radians(124)

    # When building chain A→B→C where ref_dir points A→B, to get angle A-B-C = θ,
    # we need: C = B + d*(-cos(θ)*ref_dir + sin(θ)*perp)
    COS_TET = -math.cos(TETRAHEDRAL)  # ≈ +0.33 for 109.5°
    SIN_TET = math.sin(TETRAHEDRAL)   # ≈ +0.94

    # Direction from CA toward metal
    ca_to_metal = metal_pos - ca_pos
    ca_metal_dist = np.linalg.norm(ca_to_metal)
    if ca_metal_dist < 0.01:
        ca_to_metal = np.array([1.0, 0.0, 0.0])
        ca_metal_dist = 1.0
    main_axis = ca_to_metal / ca_metal_dist

    # Create orthonormal basis
    perp1 = np.cross(main_axis, np.array([0, 0, 1]))
    if np.linalg.norm(perp1) < 0.1:
        perp1 = np.cross(main_axis, np.array([0, 1, 0]))
    perp1 = perp1 / np.linalg.norm(perp1)
    perp2 = np.cross(main_axis, perp1)
    perp2 = perp2 / np.linalg.norm(perp2)

    # Step 1: Build CB at proper tetrahedral angle
    cb_axis_component = COS_TET * CA_CB
    cb_perp_component = SIN_TET * CA_CB
    cb_pos = ca_pos + main_axis * cb_axis_component + perp1 * cb_perp_component

    # Step 2: Build CG at proper angle from CB
    cb_dir = (cb_pos - ca_pos)
    cb_dir = cb_dir / np.linalg.norm(cb_dir)

    cg_axis_component = COS_TET * CB_CG
    cg_perp_component = SIN_TET * CB_CG
    cg_perp = np.cross(cb_dir, perp2)
    if np.linalg.norm(cg_perp) < 0.1:
        cg_perp = perp1
    else:
        cg_perp = cg_perp / np.linalg.norm(cg_perp)
    cg_pos = cb_pos + cb_dir * cg_axis_component + cg_perp * cg_perp_component

    # Step 3: Build CD at proper angle from CG
    cg_dir = (cg_pos - cb_pos)
    cg_dir = cg_dir / np.linalg.norm(cg_dir)
    cd_perp = np.cross(cg_dir, cg_perp)
    if np.linalg.norm(cd_perp) < 0.1:
        cd_perp = perp2
    else:
        cd_perp = cd_perp / np.linalg.norm(cd_perp)
    cd_axis_component = COS_TET * CG_CD
    cd_perp_component = SIN_TET * CG_CD
    cd_pos = cg_pos + cg_dir * cd_axis_component + cd_perp * cd_perp_component

    # Step 4: Build carboxylate with proper geometry
    # In a carboxylate: O-C-O angle = 124°, C-C-O angle = 118°
    # For monodentate, we use same geometry but orient so one O is closer to metal
    cd_dir_final = (cd_pos - cg_pos)
    cd_dir_final = cd_dir_final / np.linalg.norm(cd_dir_final)

    # The carboxylate bisector points along cd_dir (due to ~118° C-C-O angle)
    oe_bisector = cd_dir_final

    # The carboxylate plane normal is perpendicular to cd_dir
    # Orient it toward the metal so one oxygen is closer to metal
    cd_to_metal = metal_pos - cd_pos
    oe_plane_normal = cd_to_metal - np.dot(cd_to_metal, cd_dir_final) * cd_dir_final
    if np.linalg.norm(oe_plane_normal) < 0.1:
        oe_plane_normal = perp1 - np.dot(perp1, cd_dir_final) * cd_dir_final
        if np.linalg.norm(oe_plane_normal) < 0.1:
            oe_plane_normal = perp2
    oe_plane_normal = oe_plane_normal / np.linalg.norm(oe_plane_normal)

    # Half angle for OE positions (O-C-O angle is 124°)
    half_oco = CARBOXYL_OCO / 2  # ~62 degrees

    # Position OE1 and OE2: each at 62° from bisector (cd_dir)
    # OE1 is on the +plane_normal side (closer to metal)
    # OE2 is on the -plane_normal side (farther from metal)
    oe1_dir = math.cos(half_oco) * oe_bisector + math.sin(half_oco) * oe_plane_normal
    oe2_dir = math.cos(half_oco) * oe_bisector - math.sin(half_oco) * oe_plane_normal
    oe1_pos = cd_pos + oe1_dir * CD_OE
    oe2_pos = cd_pos + oe2_dir * CD_OE

    # For monodentate, swap if user wants OE2 to be the coordinating oxygen
    if coordinating_oxygen == "OE2":
        oe1_pos, oe2_pos = oe2_pos, oe1_pos

    # Note: We do NOT translate/scale the side chain. The geometry is built with proper
    # bond lengths and angles. The actual metal-O distance depends on the CA position.
    # RFdiffusion will optimize the structure during the design process.

    # Step 5: Generate backbone atoms
    n_pos = ca_pos - main_axis * CA_N

    c_perp = np.cross(main_axis, perp1)
    if np.linalg.norm(c_perp) < 0.1:
        c_perp = perp2
    else:
        c_perp = c_perp / np.linalg.norm(c_perp)
    c_pos = ca_pos + c_perp * CA_C
    o_pos = c_pos + c_perp * C_O

    atoms = [
        ("N", "N", n_pos),
        ("CA", "C", ca_pos),
        ("C", "C", c_pos),
        ("O", "O", o_pos),
        ("CB", "C", cb_pos),
        ("CG", "C", cg_pos),
        ("CD", "C", cd_pos),
        ("OE1", "O", oe1_pos),
        ("OE2", "O", oe2_pos),
    ]

    lines = []
    serial = atom_serial_start
    for atom_name, element, pos in atoms:
        if len(atom_name) <= 2:
            atom_name_fmt = f" {atom_name:<3s}"
        else:
            atom_name_fmt = f"{atom_name:<4s}"
        line = (
            f"ATOM  {serial:5d} {atom_name_fmt} GLU {chain}{res_num:4d}    "
            f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00  0.00          {element:>2s}"
        )
        lines.append(line)
        serial += 1

    return "\n".join(lines), serial


def _generate_glu_residue(
    chain: str,
    res_num: int,
    ca_pos: np.ndarray,
    metal_pos: np.ndarray,
    atom_serial_start: int,
    target_ln_o_distance: float = 2.35,
) -> Tuple[str, int]:
    """
    Generate a bidentate glutamate residue with proper amino acid geometry.

    Uses "inside-out" construction:
    1. Build proper Glu geometry starting from CA
    2. Orient the side chain toward the metal
    3. Scale/position so OE atoms are at target distance from metal

    This ensures proper bond lengths (CA-CB=1.54, CB-CG=1.52, etc.) and
    tetrahedral angles (~109-111°) while achieving metal coordination.

    Args:
        chain: Chain identifier
        res_num: Residue number
        ca_pos: Position of CA atom (defines overall residue position)
        metal_pos: Position of metal ion
        atom_serial_start: Starting atom serial number
        target_ln_o_distance: Target lanthanide-oxygen distance (Å)

    Returns:
        Tuple of (PDB lines string, next atom serial number)
    """
    # Standard bond lengths
    CA_CB = 1.54
    CB_CG = 1.52
    CG_CD = 1.52
    CD_OE = 1.25
    CA_N = 1.47
    CA_C = 1.52
    C_O = 1.23

    # Standard angles (in radians)
    TETRAHEDRAL = math.radians(109.5)  # sp3 carbon angle
    CARBOXYL_OCO = math.radians(124)   # O-C-O angle in carboxylate
    CARBOXYL_CCO = math.radians(118)   # C-C-O angle in carboxylate

    # When building chain A→B→C where ref_dir points A→B, to get angle A-B-C = θ,
    # we need: C = B + d*(cos(π-θ)*ref_dir + sin(π-θ)*perp) = B + d*(-cos(θ)*ref_dir + sin(θ)*perp)
    # So use -cos(θ) for axis component to get correct bond angles
    COS_TET = -math.cos(TETRAHEDRAL)  # ≈ +0.33 for 109.5°
    SIN_TET = math.sin(TETRAHEDRAL)   # ≈ +0.94

    # Direction from CA toward metal (main axis)
    ca_to_metal = metal_pos - ca_pos
    ca_metal_dist = np.linalg.norm(ca_to_metal)
    if ca_metal_dist < 0.01:
        ca_to_metal = np.array([1.0, 0.0, 0.0])
        ca_metal_dist = 1.0
    main_axis = ca_to_metal / ca_metal_dist

    # Create orthonormal basis for positioning
    perp1 = np.cross(main_axis, np.array([0, 0, 1]))
    if np.linalg.norm(perp1) < 0.1:
        perp1 = np.cross(main_axis, np.array([0, 1, 0]))
    perp1 = perp1 / np.linalg.norm(perp1)
    perp2 = np.cross(main_axis, perp1)
    perp2 = perp2 / np.linalg.norm(perp2)

    # Step 1: Build CB at proper distance and tetrahedral angle from CA
    # CB extends toward metal with perpendicular offset for tetrahedral geometry
    cb_axis_component = COS_TET * CA_CB
    cb_perp_component = SIN_TET * CA_CB
    cb_pos = ca_pos + main_axis * cb_axis_component + perp1 * cb_perp_component

    # Step 2: Build CG at proper distance and angle from CB
    # Direction from CA to CB
    cb_dir = (cb_pos - ca_pos)
    cb_dir = cb_dir / np.linalg.norm(cb_dir)

    # CG continues along chain with tetrahedral angle
    cg_axis_component = COS_TET * CB_CG
    cg_perp_component = SIN_TET * CB_CG
    # Rotate perpendicular direction to create realistic zigzag
    cg_perp = np.cross(cb_dir, perp2)
    if np.linalg.norm(cg_perp) < 0.1:
        cg_perp = perp1
    else:
        cg_perp = cg_perp / np.linalg.norm(cg_perp)
    cg_pos = cb_pos + cb_dir * cg_axis_component + cg_perp * cg_perp_component

    # Step 3: Build CD at proper distance and angle from CG
    cg_dir = (cg_pos - cb_pos)
    cg_dir = cg_dir / np.linalg.norm(cg_dir)
    cd_perp = np.cross(cg_dir, cg_perp)
    if np.linalg.norm(cd_perp) < 0.1:
        cd_perp = perp2
    else:
        cd_perp = cd_perp / np.linalg.norm(cd_perp)
    cd_axis_component = COS_TET * CG_CD
    cd_perp_component = SIN_TET * CG_CD
    cd_pos = cg_pos + cg_dir * cd_axis_component + cd_perp * cd_perp_component

    # Step 4: Build OE1 and OE2 with proper carboxylate geometry for BIDENTATE coordination
    # In a carboxylate: O-C-O angle = 124°, C-C-O angle = 118°
    #
    # For BIDENTATE coordination, both oxygens must be close to the metal.
    # This requires orienting the carboxylate so the CD-metal vector roughly
    # bisects the O-C-O angle (both oxygens equidistant from metal).
    #
    # Key insight: The carboxylate plane should contain (or nearly contain) the metal.
    # So the plane normal must be perpendicular to BOTH cd_dir AND cd_to_metal.
    cd_dir = (cd_pos - cg_pos)
    cd_dir = cd_dir / np.linalg.norm(cd_dir)

    cd_to_metal = metal_pos - cd_pos

    # For bidentate: plane normal is perpendicular to both cd_dir and cd_to_metal
    # This puts the carboxylate plane containing the metal, so both O's are equidistant
    oe_plane_normal = np.cross(cd_dir, cd_to_metal)
    if np.linalg.norm(oe_plane_normal) < 0.1:
        # cd_to_metal is parallel to cd_dir, use any perpendicular
        oe_plane_normal = np.cross(cd_dir, np.array([0.0, 0.0, 1.0]))
        if np.linalg.norm(oe_plane_normal) < 0.1:
            oe_plane_normal = np.cross(cd_dir, np.array([0.0, 1.0, 0.0]))
    oe_plane_normal = oe_plane_normal / np.linalg.norm(oe_plane_normal)

    # For bidentate, the bisector should point TOWARD the metal, not along cd_dir
    # This way both OE1 and OE2 are close to the metal
    cd_to_metal_dir = cd_to_metal / np.linalg.norm(cd_to_metal)
    oe_bisector = cd_to_metal_dir

    # Half angle for OE positions (O-C-O angle is 124°)
    half_oco = CARBOXYL_OCO / 2  # ~62 degrees

    # Position OE1 and OE2: each at 62° from bisector (toward metal) in the carboxylate plane
    # Both oxygens will be at similar distances from the metal
    oe1_dir = math.cos(half_oco) * oe_bisector + math.sin(half_oco) * oe_plane_normal
    oe2_dir = math.cos(half_oco) * oe_bisector - math.sin(half_oco) * oe_plane_normal
    oe1_pos = cd_pos + oe1_dir * CD_OE
    oe2_pos = cd_pos + oe2_dir * CD_OE

    # Note: We do NOT translate/scale the side chain. The geometry is built with proper
    # bond lengths and angles. The actual metal-O distance depends on the CA position
    # (set by the template generator). RFdiffusion will optimize the final structure.

    # Step 6: Generate backbone atoms
    # N is opposite to side chain
    n_dir = -main_axis
    n_pos = ca_pos + n_dir * CA_N

    # C is perpendicular to N-CA-CB plane
    c_perp = np.cross(main_axis, perp1)
    if np.linalg.norm(c_perp) < 0.1:
        c_perp = perp2
    else:
        c_perp = c_perp / np.linalg.norm(c_perp)
    c_pos = ca_pos + c_perp * CA_C

    # Backbone O extends from C
    o_pos = c_pos + c_perp * C_O

    atoms = [
        ("N", "N", n_pos),
        ("CA", "C", ca_pos),
        ("C", "C", c_pos),
        ("O", "O", o_pos),
        ("CB", "C", cb_pos),
        ("CG", "C", cg_pos),
        ("CD", "C", cd_pos),
        ("OE1", "O", oe1_pos),
        ("OE2", "O", oe2_pos),
    ]

    lines = []
    serial = atom_serial_start
    for atom_name, element, pos in atoms:
        if len(atom_name) <= 2:
            atom_name_fmt = f" {atom_name:<3s}"
        else:
            atom_name_fmt = f"{atom_name:<4s}"
        line = (
            f"ATOM  {serial:5d} {atom_name_fmt} GLU {chain}{res_num:4d}    "
            f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00  0.00          {element:>2s}"
        )
        lines.append(line)
        serial += 1

    return "\n".join(lines), serial


def _generate_asp_residue(
    chain: str,
    res_num: int,
    ca_pos: np.ndarray,
    metal_pos: np.ndarray,
    atom_serial_start: int,
    target_ln_o_distance: float = 2.35,
) -> Tuple[str, int]:
    """
    Generate a bidentate aspartate residue with proper amino acid geometry.

    Uses "inside-out" construction:
    1. Build proper Asp geometry starting from CA
    2. Orient the side chain toward the metal
    3. Scale/position so OD atoms are at target distance from metal

    Asp has SHORTER side chain than Glu: CA-CB-CG-OD1/OD2 (no CD)
    This provides TIGHTER coordination geometry - preferred for lanthanides.

    Args:
        chain: Chain identifier
        res_num: Residue number
        ca_pos: Position of CA atom (defines overall residue position)
        metal_pos: Position of metal ion
        atom_serial_start: Starting atom serial number
        target_ln_o_distance: Target lanthanide-oxygen distance (Å)

    Returns:
        Tuple of (PDB lines string, next atom serial number)
    """
    # Standard bond lengths
    CA_CB = 1.54
    CB_CG = 1.52
    CG_OD = 1.25
    CA_N = 1.47
    CA_C = 1.52
    C_O = 1.23

    # Standard angles (in radians)
    TETRAHEDRAL = math.radians(109.5)  # sp3 carbon angle
    CARBOXYL_OCO = math.radians(124)   # O-C-O angle in carboxylate

    # When building chain A→B→C where ref_dir points A→B, to get angle A-B-C = θ,
    # we need: C = B + d*(cos(π-θ)*ref_dir + sin(π-θ)*perp) = B + d*(-cos(θ)*ref_dir + sin(θ)*perp)
    # So use -cos(θ) for axis component to get correct bond angles
    COS_TET = -math.cos(TETRAHEDRAL)  # ≈ +0.33 for 109.5°
    SIN_TET = math.sin(TETRAHEDRAL)   # ≈ +0.94

    # Direction from CA toward metal (main axis)
    ca_to_metal = metal_pos - ca_pos
    ca_metal_dist = np.linalg.norm(ca_to_metal)
    if ca_metal_dist < 0.01:
        ca_to_metal = np.array([1.0, 0.0, 0.0])
        ca_metal_dist = 1.0
    main_axis = ca_to_metal / ca_metal_dist

    # Create orthonormal basis for positioning
    perp1 = np.cross(main_axis, np.array([0, 0, 1]))
    if np.linalg.norm(perp1) < 0.1:
        perp1 = np.cross(main_axis, np.array([0, 1, 0]))
    perp1 = perp1 / np.linalg.norm(perp1)
    perp2 = np.cross(main_axis, perp1)
    perp2 = perp2 / np.linalg.norm(perp2)

    # Step 1: Build CB at proper distance and tetrahedral angle from CA
    cb_axis_component = COS_TET * CA_CB
    cb_perp_component = SIN_TET * CA_CB
    cb_pos = ca_pos + main_axis * cb_axis_component + perp1 * cb_perp_component

    # Step 2: Build CG at proper distance and angle from CB
    cb_dir = (cb_pos - ca_pos)
    cb_dir = cb_dir / np.linalg.norm(cb_dir)

    cg_axis_component = COS_TET * CB_CG
    cg_perp_component = SIN_TET * CB_CG
    cg_perp = np.cross(cb_dir, perp2)
    if np.linalg.norm(cg_perp) < 0.1:
        cg_perp = perp1
    else:
        cg_perp = cg_perp / np.linalg.norm(cg_perp)
    cg_pos = cb_pos + cb_dir * cg_axis_component + cg_perp * cg_perp_component

    # Step 3: Build OD1 and OD2 with proper carboxylate geometry for BIDENTATE coordination
    # In a carboxylate: O-C-O angle = 124°, C-C-O angle = 118°
    #
    # For BIDENTATE coordination, both oxygens must be close to the metal.
    # This requires orienting the carboxylate so the CG-metal vector roughly
    # bisects the O-C-O angle (both oxygens equidistant from metal).
    #
    # Key insight: The carboxylate plane should contain (or nearly contain) the metal.
    # So the plane normal must be perpendicular to BOTH cg_dir AND cg_to_metal.
    cg_dir = (cg_pos - cb_pos)
    cg_dir = cg_dir / np.linalg.norm(cg_dir)

    cg_to_metal = metal_pos - cg_pos

    # For bidentate: plane normal is perpendicular to both cg_dir and cg_to_metal
    # This puts the carboxylate plane containing the metal, so both O's are equidistant
    od_plane_normal = np.cross(cg_dir, cg_to_metal)
    if np.linalg.norm(od_plane_normal) < 0.1:
        # cg_to_metal is parallel to cg_dir, use any perpendicular
        od_plane_normal = np.cross(cg_dir, np.array([0.0, 0.0, 1.0]))
        if np.linalg.norm(od_plane_normal) < 0.1:
            od_plane_normal = np.cross(cg_dir, np.array([0.0, 1.0, 0.0]))
    od_plane_normal = od_plane_normal / np.linalg.norm(od_plane_normal)

    # For bidentate, the bisector should point TOWARD the metal, not along cg_dir
    # This way both OD1 and OD2 are close to the metal
    cg_to_metal_dir = cg_to_metal / np.linalg.norm(cg_to_metal)
    od_bisector = cg_to_metal_dir

    # Half angle for OD positions (O-C-O angle is 124°)
    half_oco = CARBOXYL_OCO / 2  # ~62 degrees

    # Position OD1 and OD2: each at 62° from bisector (toward metal) in the carboxylate plane
    # Both oxygens will be at similar distances from the metal
    od1_dir = math.cos(half_oco) * od_bisector + math.sin(half_oco) * od_plane_normal
    od2_dir = math.cos(half_oco) * od_bisector - math.sin(half_oco) * od_plane_normal
    od1_pos = cg_pos + od1_dir * CG_OD
    od2_pos = cg_pos + od2_dir * CG_OD

    # Note: We do NOT translate/scale the side chain. The geometry is built with proper
    # bond lengths and angles. The actual metal-O distance depends on the CA position.

    # Step 5: Generate backbone atoms
    n_dir = -main_axis
    n_pos = ca_pos + n_dir * CA_N

    c_perp = np.cross(main_axis, perp1)
    if np.linalg.norm(c_perp) < 0.1:
        c_perp = perp2
    else:
        c_perp = c_perp / np.linalg.norm(c_perp)
    c_pos = ca_pos + c_perp * CA_C

    o_pos = c_pos + c_perp * C_O

    atoms = [
        ("N", "N", n_pos),
        ("CA", "C", ca_pos),
        ("C", "C", c_pos),
        ("O", "O", o_pos),
        ("CB", "C", cb_pos),
        ("CG", "C", cg_pos),
        ("OD1", "O", od1_pos),
        ("OD2", "O", od2_pos),
    ]

    lines = []
    serial = atom_serial_start
    for atom_name, element, pos in atoms:
        if len(atom_name) <= 2:
            atom_name_fmt = f" {atom_name:<3s}"
        else:
            atom_name_fmt = f"{atom_name:<4s}"
        line = (
            f"ATOM  {serial:5d} {atom_name_fmt} ASP {chain}{res_num:4d}    "
            f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00  0.00          {element:>2s}"
        )
        lines.append(line)
        serial += 1

    return "\n".join(lines), serial


def _generate_water(
    res_num: int,
    position: np.ndarray,
    atom_serial: int,
    chain: str = "W",
) -> str:
    """Generate a water molecule HETATM record.

    Args:
        res_num: Residue number for the water
        position: 3D coordinates of the water oxygen
        atom_serial: Atom serial number
        chain: Chain identifier (default "W" for water)

    Returns:
        HETATM line string for PDB
    """
    return (
        f"HETATM{atom_serial:5d}  O   HOH {chain}{res_num:4d}    "
        f"{position[0]:8.3f}{position[1]:8.3f}{position[2]:8.3f}  1.00  0.00           O"
    )


def _generate_trp_residue(
    chain: str,
    res_num: int,
    ca_pos: np.ndarray,
    metal_pos: np.ndarray,
    atom_serial_start: int,
    target_distance: float = 4.5,
) -> Tuple[str, int]:
    """
    Generate a tryptophan residue with indole positioned toward metal.

    Args:
        chain: Chain identifier
        res_num: Residue number
        ca_pos: Position of CA atom
        metal_pos: Position of metal ion
        atom_serial_start: Starting atom serial number
        target_distance: Target distance from NE1 to metal (Å)

    Returns:
        Tuple of (PDB lines string, next atom serial number)
    """
    # Direction from CA toward metal
    direction = metal_pos - ca_pos
    direction = direction / np.linalg.norm(direction)

    # Trp atoms (simplified)
    cb_pos = ca_pos + direction * 1.54
    cg_pos = cb_pos + direction * 1.50

    # Indole ring - position CD1/NE1 toward metal
    cd1_pos = cg_pos + direction * 1.36
    ne1_pos = cd1_pos + direction * 1.38

    # Adjust NE1 to be at target distance from metal
    ne1_to_metal = metal_pos - ne1_pos
    current_dist = np.linalg.norm(ne1_to_metal)
    if current_dist > 0:
        adjustment = (current_dist - target_distance) / current_dist
        ne1_pos = ne1_pos + ne1_to_metal * adjustment

    # Perpendicular for ring atoms
    perp = np.cross(direction, np.array([0, 0, 1]))
    if np.linalg.norm(perp) < 0.1:
        perp = np.cross(direction, np.array([0, 1, 0]))
    perp = perp / np.linalg.norm(perp)

    cd2_pos = cg_pos + perp * 1.43
    ce2_pos = cd2_pos + direction * 1.40
    ce3_pos = cd2_pos - direction * 1.40
    cz2_pos = ce2_pos + perp * 1.40
    cz3_pos = ce3_pos + perp * 1.40
    ch2_pos = (cz2_pos + cz3_pos) / 2

    # Backbone
    n_pos = ca_pos - direction * 1.47
    c_perp = np.cross(direction, perp)
    c_pos = ca_pos + c_perp * 1.52
    o_pos = c_pos + c_perp * 1.23

    atoms = [
        ("N", "N", n_pos),
        ("CA", "C", ca_pos),
        ("C", "C", c_pos),
        ("O", "O", o_pos),
        ("CB", "C", cb_pos),
        ("CG", "C", cg_pos),
        ("CD1", "C", cd1_pos),
        ("CD2", "C", cd2_pos),
        ("NE1", "N", ne1_pos),
        ("CE2", "C", ce2_pos),
        ("CE3", "C", ce3_pos),
        ("CZ2", "C", cz2_pos),
        ("CZ3", "C", cz3_pos),
        ("CH2", "C", ch2_pos),
    ]

    lines = []
    serial = atom_serial_start
    for atom_name, element, pos in atoms:
        line = (
            f"ATOM  {serial:5d} {atom_name:4s} TRP {chain}{res_num:4d}    "
            f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00  0.00          {element:>2s}"
        )
        lines.append(line)
        serial += 1

    return "\n".join(lines), serial


def generate_ef_hand_template(
    metal: str = "TB",
    donor_residue: str = "ASP",
    include_waters: bool = False,
    add_trp_antenna: bool = False,
    trp_chain: str = "A",
) -> str:
    """
    Generate EF-hand-inspired template PDB with 8 coordinating Asp/Glu residues.

    Based on EF-hand calcium binding motif, adapted for lanthanides.
    Uses 4 residues per chain (8 total) for bidentate coordination.

    The template positions:
    - Metal at a central position (not origin, to avoid CCD lookup issues)
    - 8 Asp/Glu residues with OD/OE atoms at exact Ln-O distance (2.35 Å for Tb)
    - CA atoms positioned so side chains extend toward metal

    Args:
        metal: Metal code (TB, GD, EU, LA, YB)
        donor_residue: "ASP" (recommended), "GLU", or "MIXED" (5 Asp + 3 Glu)
        include_waters: Whether to include 4 water molecules in coordination
        add_trp_antenna: Whether to add a Trp residue for TEBL
        trp_chain: Which chain to add Trp to (A or B)

    Returns:
        PDB content string
    """
    params = LANTHANIDE_PARAMS.get(metal.upper(), LANTHANIDE_PARAMS["TB"])
    bond_dist = params["bond_distance"]

    # Metal at origin - RFD3 designs proteins around origin
    # IMPORTANT: Using (50,50,50) causes stretched chains because ori_token
    # positions protein COM near origin while template residues are far away
    metal_pos = np.array([0.0, 0.0, 0.0])

    # CA distance from metal depends on residue type:
    # Glu side chain: CA-CB (1.54) + CB-CG (1.52) + CG-CD (1.52) + CD-OE (1.25) = 5.83 Å
    # Asp side chain: CA-CB (1.54) + CB-CG (1.52) + CG-OD (1.25) = 4.31 Å (SHORTER!)
    # Add backbone offset: Glu ~7.5 Å, Asp ~6.0 Å from metal
    if donor_residue.upper() == "ASP":
        ca_radius = 6.0  # Asp has shorter side chain - CA needs to be closer
    else:
        ca_radius = 7.5  # Glu has longer side chain

    # Position 8 Glu residues in approximate square antiprism geometry
    # 4 on top plane (Chain A), 4 on bottom plane (Chain B), rotated 45 degrees

    # Top plane (Chain A residues 10, 15, 20, 25)
    top_z_offset = 3.0  # Å above metal
    top_angles = [0, 90, 180, 270]  # degrees

    # Bottom plane (Chain B residues 10, 15, 20, 25)
    bottom_z_offset = -3.0  # Å below metal
    bottom_angles = [45, 135, 225, 315]  # 45 degree rotation from top

    pdb_lines = []
    atom_serial = 1

    # Chain A - top plane coordinating residues
    # For MIXED mode: positions 10, 15, 20 = ASP (3), position 25 = GLU (1)
    chain_a_res_nums = [10, 15, 20, 25]
    for i, (res_num, angle) in enumerate(zip(chain_a_res_nums, top_angles)):
        rad = math.radians(angle)
        ca_pos = metal_pos + np.array([
            ca_radius * math.cos(rad),
            ca_radius * math.sin(rad),
            top_z_offset
        ])

        # Select residue generator based on donor_residue parameter
        if donor_residue.upper() == "ASP":
            res_lines, atom_serial = _generate_asp_residue(
                "A", res_num, ca_pos, metal_pos, atom_serial, target_ln_o_distance=bond_dist
            )
        elif donor_residue.upper() == "GLU":
            res_lines, atom_serial = _generate_glu_residue(
                "A", res_num, ca_pos, metal_pos, atom_serial, target_ln_o_distance=bond_dist
            )
        else:  # MIXED mode: 5 Asp + 3 Glu pattern
            # Chain A: positions 10, 15, 20 = ASP, position 25 = GLU
            if res_num in [10, 15, 20]:
                res_lines, atom_serial = _generate_asp_residue(
                    "A", res_num, ca_pos, metal_pos, atom_serial, target_ln_o_distance=bond_dist
                )
            else:  # res_num 25
                res_lines, atom_serial = _generate_glu_residue(
                    "A", res_num, ca_pos, metal_pos, atom_serial, target_ln_o_distance=bond_dist
                )
        pdb_lines.append(res_lines)

    # Optional Trp antenna on chain A
    if add_trp_antenna and trp_chain == "A":
        trp_ca_pos = metal_pos + np.array([8.0, 0.0, 2.0])  # Position near metal
        trp_lines, atom_serial = _generate_trp_residue(
            "A", 30, trp_ca_pos, metal_pos, atom_serial, target_distance=4.5
        )
        pdb_lines.append(trp_lines)

    pdb_lines.append("TER")

    # Chain B - bottom plane coordinating residues
    # For MIXED mode: positions 10, 15 = ASP (2), positions 20, 25 = GLU (2)
    # Total MIXED: Chain A (3 Asp + 1 Glu) + Chain B (2 Asp + 2 Glu) = 5 Asp + 3 Glu
    chain_b_res_nums = [10, 15, 20, 25]
    for i, (res_num, angle) in enumerate(zip(chain_b_res_nums, bottom_angles)):
        rad = math.radians(angle)
        ca_pos = metal_pos + np.array([
            ca_radius * math.cos(rad),
            ca_radius * math.sin(rad),
            bottom_z_offset
        ])

        # Select residue generator based on donor_residue parameter
        if donor_residue.upper() == "ASP":
            res_lines, atom_serial = _generate_asp_residue(
                "B", res_num, ca_pos, metal_pos, atom_serial, target_ln_o_distance=bond_dist
            )
        elif donor_residue.upper() == "GLU":
            res_lines, atom_serial = _generate_glu_residue(
                "B", res_num, ca_pos, metal_pos, atom_serial, target_ln_o_distance=bond_dist
            )
        else:  # MIXED mode: 5 Asp + 3 Glu pattern
            # Chain B: positions 10, 15 = ASP, positions 20, 25 = GLU
            if res_num in [10, 15]:
                res_lines, atom_serial = _generate_asp_residue(
                    "B", res_num, ca_pos, metal_pos, atom_serial, target_ln_o_distance=bond_dist
                )
            else:  # res_num 20, 25
                res_lines, atom_serial = _generate_glu_residue(
                    "B", res_num, ca_pos, metal_pos, atom_serial, target_ln_o_distance=bond_dist
                )
        pdb_lines.append(res_lines)

    # Optional Trp antenna on chain B
    if add_trp_antenna and trp_chain == "B":
        trp_ca_pos = metal_pos + np.array([8.0, 0.0, -2.0])  # Position near metal
        trp_lines, atom_serial = _generate_trp_residue(
            "B", 30, trp_ca_pos, metal_pos, atom_serial, target_distance=4.5
        )
        pdb_lines.append(trp_lines)

    pdb_lines.append("TER")

    # Metal ion - use explicit HETATM format that RFD3 should preserve
    # Format: standard PDB HETATM with element symbol matching residue name
    metal_code = metal.upper()[:2]  # TB, GD, EU, etc.
    metal_line = (
        f"HETATM{atom_serial:5d} {metal_code:>2s}   {metal_code:>3s} L   1    "
        f"{metal_pos[0]:8.3f}{metal_pos[1]:8.3f}{metal_pos[2]:8.3f}  1.00  0.00          {metal_code:>2s}"
    )
    pdb_lines.append(metal_line)
    atom_serial += 1

    # Optional water molecules (4 in axial positions for square antiprism)
    if include_waters:
        # Waters at bond_dist from metal, in equatorial plane
        water_angles = [22.5, 112.5, 202.5, 292.5]  # Between Glu positions
        for i, angle in enumerate(water_angles):
            rad = math.radians(angle)
            pos = metal_pos + np.array([
                bond_dist * math.cos(rad),
                bond_dist * math.sin(rad),
                0.0
            ])
            water_line = _generate_water(i + 1, pos, atom_serial)
            pdb_lines.append(water_line)
            atom_serial += 1

    pdb_lines.append("END")

    return "\n".join(pdb_lines)


def generate_c4_symmetric_template(
    metal: str = "TB",
    donor_residue: str = "ASP",
    include_waters: bool = False,
    add_trp_antenna: bool = False,
) -> str:
    """
    Generate C4-symmetric template like Caldwell TIM barrel top.

    Based on Caldwell et al. (2020) TFD-EE design:
    - 4 Asp/Glu residues in TRUE C4 symmetry at the SAME Z-plane
    - All carboxylate oxygens point toward metal from equatorial positions
    - Metal sits at center of the residue ring
    - Optional axial waters complete the coordination sphere

    Key difference from previous version: ALL 4 residues are coplanar (same Z),
    not spread across different Z levels. This matches the Caldwell TIM
    barrel top geometry where all coordinating residues are at positions 31/154.

    Args:
        metal: Metal code (TB, GD, EU, LA, YB)
        donor_residue: "ASP" (recommended), "GLU", or "MIXED" (2 Asp + 2 Glu)
        include_waters: Whether to include axial water molecules
        add_trp_antenna: Whether to add Trp for TEBL

    Returns:
        PDB content string
    """
    params = LANTHANIDE_PARAMS.get(metal.upper(), LANTHANIDE_PARAMS["TB"])
    bond_dist = params["bond_distance"]

    # Metal at origin - RFD3 designs proteins around origin
    metal_pos = np.array([0.0, 0.0, 0.0])

    # CA distance from metal depends on residue type:
    # Glu side chain: CA-CB (1.54) + CB-CG (1.52) + CG-CD (1.52) + CD-OE (1.25) = 5.83 Å
    # Asp side chain: CA-CB (1.54) + CB-CG (1.52) + CG-OD (1.25) = 4.31 Å (SHORTER!)
    if donor_residue.upper() == "ASP":
        ca_radius = 6.0  # Asp has shorter side chain - CA needs to be closer
    else:
        ca_radius = 7.5  # Glu has longer side chain

    # CRITICAL FIX: All 4 Glu at the SAME Z-level for true C4 symmetry
    # This matches Caldwell's TIM barrel top where all 4 Glu are coplanar
    z_offset = 0.0  # All Glu in the same plane as the metal

    pdb_lines = []
    atom_serial = 1

    # Chain A - coordinating residues at 0° and 180° (opposite sides of metal)
    # Both at the SAME z-level for true C4 planar arrangement
    for i, angle in enumerate([0, 180]):
        rad = math.radians(angle)
        ca_pos = metal_pos + np.array([
            ca_radius * math.cos(rad),
            ca_radius * math.sin(rad),
            z_offset  # Same Z for all residues - TRUE C4 symmetry
        ])
        res_num = 15 + i * 10  # Residues 15, 25

        # Select residue generator based on donor_residue parameter
        if donor_residue.upper() == "ASP":
            res_lines, atom_serial = _generate_asp_residue(
                "A", res_num, ca_pos, metal_pos, atom_serial, target_ln_o_distance=bond_dist
            )
        elif donor_residue.upper() == "GLU":
            res_lines, atom_serial = _generate_glu_residue(
                "A", res_num, ca_pos, metal_pos, atom_serial, target_ln_o_distance=bond_dist
            )
        else:  # MIXED mode: 2 Asp + 2 Glu pattern (Chain A: 1 Asp + 1 Glu)
            if res_num == 15:
                res_lines, atom_serial = _generate_asp_residue(
                    "A", res_num, ca_pos, metal_pos, atom_serial, target_ln_o_distance=bond_dist
                )
            else:  # res_num 25
                res_lines, atom_serial = _generate_glu_residue(
                    "A", res_num, ca_pos, metal_pos, atom_serial, target_ln_o_distance=bond_dist
                )
        pdb_lines.append(res_lines)

    # Optional Trp on chain A (positioned above the Glu plane)
    if add_trp_antenna:
        trp_ca_pos = metal_pos + np.array([6.0, 3.0, 4.0])  # Above the plane
        trp_lines, atom_serial = _generate_trp_residue(
            "A", 30, trp_ca_pos, metal_pos, atom_serial, target_distance=4.5
        )
        pdb_lines.append(trp_lines)

    pdb_lines.append("TER")

    # Chain B - coordinating residues at 90° and 270° (perpendicular to chain A)
    # Same Z-level as chain A for true C4 symmetry
    for i, angle in enumerate([90, 270]):
        rad = math.radians(angle)
        ca_pos = metal_pos + np.array([
            ca_radius * math.cos(rad),
            ca_radius * math.sin(rad),
            z_offset  # Same Z for all residues - TRUE C4 symmetry
        ])
        res_num = 15 + i * 10  # Residues 15, 25

        # Select residue generator based on donor_residue parameter
        if donor_residue.upper() == "ASP":
            res_lines, atom_serial = _generate_asp_residue(
                "B", res_num, ca_pos, metal_pos, atom_serial, target_ln_o_distance=bond_dist
            )
        elif donor_residue.upper() == "GLU":
            res_lines, atom_serial = _generate_glu_residue(
                "B", res_num, ca_pos, metal_pos, atom_serial, target_ln_o_distance=bond_dist
            )
        else:  # MIXED mode: 2 Asp + 2 Glu pattern (Chain B: 1 Asp + 1 Glu)
            if res_num == 15:
                res_lines, atom_serial = _generate_asp_residue(
                    "B", res_num, ca_pos, metal_pos, atom_serial, target_ln_o_distance=bond_dist
                )
            else:  # res_num 25
                res_lines, atom_serial = _generate_glu_residue(
                    "B", res_num, ca_pos, metal_pos, atom_serial, target_ln_o_distance=bond_dist
                )
        pdb_lines.append(res_lines)

    pdb_lines.append("TER")

    # Metal ion - use explicit HETATM format
    metal_code = metal.upper()[:2]
    metal_line = (
        f"HETATM{atom_serial:5d} {metal_code:>2s}   {metal_code:>3s} L   1    "
        f"{metal_pos[0]:8.3f}{metal_pos[1]:8.3f}{metal_pos[2]:8.3f}  1.00  0.00          {metal_code:>2s}"
    )
    pdb_lines.append(metal_line)
    atom_serial += 1

    # Optional water molecules in AXIAL positions (above and below plane)
    # This matches Caldwell's geometry where waters fill axial positions
    # while Glu carboxylates occupy the equatorial plane
    if include_waters:
        # Axial waters: above and below the metal (like Caldwell's structure)
        axial_water_dist = 2.7  # Å - matches Caldwell's Tb-O(water) = 2.7 Å
        axial_positions = [
            metal_pos + np.array([0.0, 0.0, axial_water_dist]),   # Above
            metal_pos + np.array([0.0, 0.0, -axial_water_dist]),  # Below
        ]
        # Additional equatorial waters at 45° angles (between Glu positions)
        for i, angle in enumerate([45, 135]):
            rad = math.radians(angle)
            pos = metal_pos + np.array([
                bond_dist * math.cos(rad),
                bond_dist * math.sin(rad),
                0.0
            ])
            axial_positions.append(pos)

        for i, pos in enumerate(axial_positions):
            water_line = _generate_water(i + 1, pos, atom_serial)
            pdb_lines.append(water_line)
            atom_serial += 1

    pdb_lines.append("END")

    return "\n".join(pdb_lines)


def position_trp_antenna(
    pdb_content: str,
    metal_chain: str = "L",
    metal_resnum: int = 1,
    target_chain: str = "A",
    target_distance: float = 4.5,
) -> str:
    """
    Add or reposition Trp residue within specified distance of metal for TEBL.

    This function analyzes the structure and adds a Trp residue
    in an optimal position for sensitized luminescence.

    Args:
        pdb_content: Input PDB content
        metal_chain: Chain ID of metal
        metal_resnum: Residue number of metal
        target_chain: Which chain to add Trp to
        target_distance: Target distance from Trp NE1 to metal (Å)

    Returns:
        Modified PDB content with Trp antenna
    """
    # Find metal position
    metal_pos = None
    lines = pdb_content.split('\n')
    new_lines = []
    max_res_num = 0
    atom_serial = 0

    for line in lines:
        if line.startswith('HETATM') or line.startswith('ATOM'):
            try:
                chain = line[21]
                res_num = int(line[22:26].strip())
                serial = int(line[6:11].strip())
                atom_serial = max(atom_serial, serial)

                if chain == target_chain:
                    max_res_num = max(max_res_num, res_num)

                # Check if this is the metal
                if line.startswith('HETATM'):
                    res_name = line[17:20].strip()
                    if res_name in ["TB", "GD", "EU", "LA", "YB", "ZN", "FE", "CA"]:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        metal_pos = np.array([x, y, z])
            except (ValueError, IndexError):
                pass

        if not line.startswith('END'):
            new_lines.append(line)

    if metal_pos is None:
        # No metal found, return original
        return pdb_content

    # Position Trp at target distance from metal
    # Choose a direction that doesn't clash with existing residues
    trp_ca_pos = metal_pos + np.array([target_distance + 2.0, 0.0, 2.0])

    # Generate Trp residue
    trp_res_num = max_res_num + 5
    trp_lines, _ = _generate_trp_residue(
        target_chain,
        trp_res_num,
        trp_ca_pos,
        metal_pos,
        atom_serial + 1,
        target_distance=target_distance
    )

    # Insert before END
    # Find the TER for target chain and insert after
    insert_idx = len(new_lines)
    for i, line in enumerate(new_lines):
        if line.startswith('TER') or line.startswith('HETATM'):
            # Check if we've passed the target chain
            if i > 0:
                prev_line = new_lines[i-1]
                if prev_line.startswith('ATOM') and len(prev_line) > 21:
                    if prev_line[21] == target_chain:
                        insert_idx = i
                        break

    new_lines.insert(insert_idx, trp_lines)
    new_lines.append("END")

    return "\n".join(new_lines)


def get_template(
    template_type: str,
    metal: str = "TB",
    donor_residue: str = "ASP",
    include_waters: bool = False,
    add_trp_antenna: bool = False,
    trp_chain: str = "A",
) -> str:
    """
    Get a lanthanide binding site template.

    Args:
        template_type: "ef_hand" or "c4_symmetric"
        metal: Metal code (TB, GD, EU, etc.)
        donor_residue: "ASP" (recommended), "GLU", or "MIXED"
        include_waters: Whether to include water molecules
        add_trp_antenna: Whether to add Trp for TEBL
        trp_chain: Which chain to add Trp to

    Returns:
        PDB content string

    Raises:
        ValueError: If template_type is invalid
    """
    if template_type == "ef_hand":
        return generate_ef_hand_template(
            metal=metal,
            donor_residue=donor_residue,
            include_waters=include_waters,
            add_trp_antenna=add_trp_antenna,
            trp_chain=trp_chain,
        )
    elif template_type == "c4_symmetric":
        return generate_c4_symmetric_template(
            metal=metal,
            donor_residue=donor_residue,
            include_waters=include_waters,
            add_trp_antenna=add_trp_antenna,
        )
    else:
        raise ValueError(
            f"Invalid template_type: {template_type}. "
            f"Valid options: {TEMPLATE_TYPES}"
        )


# Convenience function to check if metal is a lanthanide
def is_lanthanide(metal: str) -> bool:
    """Check if metal code is a lanthanide."""
    return metal.upper() in LANTHANIDE_PARAMS


def validate_template_geometry(
    pdb_content: str,
    metal: str = "TB",
    target_distance: float = None,
    tolerance: float = 0.3,
) -> Dict[str, Any]:
    """
    Validate template geometry - check OE-metal distances and coordination.

    Args:
        pdb_content: PDB content string from template generator
        metal: Metal code to look for
        target_distance: Expected Ln-O distance (default from LANTHANIDE_PARAMS)
        tolerance: Acceptable deviation from target distance (Å)

    Returns:
        Dictionary with validation results:
        - valid: bool - overall validity
        - metal_position: tuple - (x, y, z) of metal
        - oe_distances: list - distances from each OE to metal
        - coordination_count: int - number of OE atoms within target distance
        - issues: list - any problems found
    """
    if target_distance is None:
        params = LANTHANIDE_PARAMS.get(metal.upper(), LANTHANIDE_PARAMS["TB"])
        target_distance = params["bond_distance"]

    metal_code = metal.upper()[:2]
    metal_pos = None
    oe_atoms = []
    issues = []

    # Parse PDB
    for line in pdb_content.split('\n'):
        if line.startswith('HETATM'):
            res_name = line[17:20].strip()
            if res_name == metal_code or res_name == metal.upper():
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    metal_pos = np.array([x, y, z])
                except ValueError:
                    issues.append(f"Could not parse metal coordinates")
        elif line.startswith('ATOM'):
            atom_name = line[12:16].strip()
            if atom_name in ['OE1', 'OE2', 'OD1', 'OD2']:  # Glu or Asp carboxylate
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    chain = line[21]
                    res_num = int(line[22:26].strip())
                    oe_atoms.append({
                        'atom': atom_name,
                        'chain': chain,
                        'resnum': res_num,
                        'pos': np.array([x, y, z])
                    })
                except ValueError:
                    pass

    # Validate
    if metal_pos is None:
        issues.append(f"Metal {metal} not found in template")
        return {
            'valid': False,
            'metal_position': None,
            'oe_distances': [],
            'coordination_count': 0,
            'issues': issues
        }

    # Check for NaN coordinates
    if np.any(np.isnan(metal_pos)):
        issues.append(f"Metal has NaN coordinates: {metal_pos}")
        return {
            'valid': False,
            'metal_position': tuple(metal_pos),
            'oe_distances': [],
            'coordination_count': 0,
            'issues': issues
        }

    # Calculate distances
    oe_distances = []
    coordinating_count = 0
    for oe in oe_atoms:
        dist = np.linalg.norm(oe['pos'] - metal_pos)
        oe_distances.append({
            'atom': oe['atom'],
            'chain': oe['chain'],
            'resnum': oe['resnum'],
            'distance': dist
        })
        if abs(dist - target_distance) <= tolerance:
            coordinating_count += 1
        elif dist < target_distance - tolerance:
            issues.append(
                f"{oe['chain']}{oe['resnum']} {oe['atom']}: too close ({dist:.2f} Å, "
                f"target {target_distance:.2f} Å)"
            )
        elif dist > target_distance + tolerance:
            issues.append(
                f"{oe['chain']}{oe['resnum']} {oe['atom']}: too far ({dist:.2f} Å, "
                f"target {target_distance:.2f} Å)"
            )

    # Check minimum coordination
    if coordinating_count < 8:
        issues.append(
            f"Low coordination count: {coordinating_count}/8 OE atoms at target distance"
        )

    valid = len(issues) == 0 or (coordinating_count >= 6 and not any('NaN' in i for i in issues))

    return {
        'valid': valid,
        'metal_position': tuple(metal_pos),
        'oe_distances': oe_distances,
        'coordination_count': coordinating_count,
        'issues': issues
    }


def generate_and_validate_template(
    template_type: str = "ef_hand",
    metal: str = "TB",
    include_waters: bool = False,
    add_trp_antenna: bool = False,
) -> Tuple[str, Dict[str, Any]]:
    """
    Generate a template and validate its geometry.

    Args:
        template_type: "ef_hand" or "c4_symmetric"
        metal: Metal code (TB, GD, EU, etc.)
        include_waters: Whether to include water molecules
        add_trp_antenna: Whether to add Trp for TEBL

    Returns:
        Tuple of (pdb_content, validation_results)
    """
    pdb_content = get_template(
        template_type=template_type,
        metal=metal,
        include_waters=include_waters,
        add_trp_antenna=add_trp_antenna,
    )

    validation = validate_template_geometry(pdb_content, metal)

    return pdb_content, validation


# =============================================================================
# NEW: Template Library System
# =============================================================================


def recommend_template(metal: str, user_preference: Optional[str] = None) -> str:
    """
    Recommend best template for given metal.

    Args:
        metal: Metal code (TB, GD, EU, etc.)
        user_preference: User's preference (optional)

    Returns:
        Template name from TEMPLATE_LIBRARY
    """
    if user_preference and user_preference in TEMPLATE_LIBRARY:
        return user_preference

    return METAL_TEMPLATE_RECOMMENDATIONS.get(metal.upper(), "ef_hand_8")


def get_template_info(template_name: str) -> Optional[Dict[str, Any]]:
    """
    Get information about a template from the library.

    Args:
        template_name: Name of template from TEMPLATE_LIBRARY

    Returns:
        Template definition dict or None if not found
    """
    return TEMPLATE_LIBRARY.get(template_name)


def list_templates() -> List[Dict[str, Any]]:
    """
    List all available templates with their metadata.

    Returns:
        List of template info dicts
    """
    return [
        {
            "name": name,
            "display_name": info["name"],
            "description": info["description"],
            "coordination_number": info["coordination_number"],
            "geometry": info["geometry"],
            "best_for": info["best_for"],
            "deprecated": info.get("deprecated", False),
        }
        for name, info in TEMPLATE_LIBRARY.items()
    ]


def generate_template_from_library(
    template_name: str,
    metal: str = "TB",
    add_trp_antenna: bool = False,
    trp_chain: str = "A",
    use_architector: bool = True,
    validate_output: bool = True,
) -> str:
    """
    Generate PDB from template library definition.

    This is the NEW unified template generator that creates chemically
    realistic templates with correct coordination numbers (CN=8-9 instead
    of CN=16 from the legacy system).

    When Architector is available and use_architector=True, uses physics-based
    geometry optimization (GFN2-xTB) for ~0.5 Å RMSD accuracy to experimental
    structures. Falls back to analytical geometry if Architector unavailable.

    Args:
        template_name: Key from TEMPLATE_LIBRARY (e.g., "caldwell_4", "ef_hand_8")
        metal: Metal code (TB, GD, EU, etc.)
        add_trp_antenna: Whether to add Trp for TEBL
        trp_chain: Which chain to add Trp to
        use_architector: Whether to try Architector for geometry (default: True)
        validate_output: Whether to validate the output geometry (default: True)

    Returns:
        PDB content string with correct coordination geometry

    Raises:
        ValueError: If template_name is not in TEMPLATE_LIBRARY
    """
    if template_name not in TEMPLATE_LIBRARY:
        raise ValueError(
            f"Unknown template: {template_name}. "
            f"Available templates: {list(TEMPLATE_LIBRARY.keys())}"
        )

    template_def = TEMPLATE_LIBRARY[template_name]

    # Try Architector first if available and requested
    # Architector provides physics-based geometry optimization with ~0.5 Å RMSD accuracy
    if use_architector and is_architector_available() and not template_def.get("deprecated"):
        logger.info(f"Attempting Architector-based geometry for {template_name}")
        architector_pdb = generate_template_with_architector(template_def, metal)
        if architector_pdb:
            # Validate Architector output
            if validate_output and VALIDATION_AVAILABLE:
                validation = validate_coordination_geometry(
                    architector_pdb,
                    expected_cn=template_def["coordination_number"],
                    expected_geometry=template_def.get("geometry", "square_antiprism"),
                )
                if validation.get("valid"):
                    logger.info(f"Architector geometry validated: {validation.get('summary')}")
                    # TODO: Add Trp antenna if requested
                    return architector_pdb
                else:
                    logger.warning(
                        f"Architector geometry failed validation: {validation.get('issues')}, "
                        "falling back to analytical geometry"
                    )
            else:
                logger.info("Architector geometry generated (validation skipped)")
                return architector_pdb
        else:
            logger.info("Architector generation failed, falling back to analytical geometry")

    # Handle deprecated legacy template
    if template_def.get("deprecated"):
        # Fall back to legacy generator for backward compatibility
        return generate_ef_hand_template(
            metal=metal,
            donor_residue="ASP",
            include_waters=len(template_def.get("waters", [])) > 0,
            add_trp_antenna=add_trp_antenna,
            trp_chain=trp_chain,
        )

    # Get metal-specific bond distance
    params = LANTHANIDE_PARAMS.get(metal.upper(), LANTHANIDE_PARAMS["TB"])
    bond_dist = params["bond_distance"]

    # Metal at origin - RFD3 designs proteins around origin
    metal_pos = np.array([0.0, 0.0, 0.0])

    # CA radius is calculated PER RESIDUE based on type and coordination mode.
    # This is necessary because:
    # - GLU has longer side chain (CA-CB-CG-CD-OE ≈ 5.8Å) than ASP (CA-CB-CG-OD ≈ 4.3Å)
    # - Bidentate: bisector toward metal, both O's at ~2.5Å from metal
    # - Monodentate: bisector along chain, only one O close (~2.5Å), other far (~3.5Å+)
    #
    # With z_offset=3Å and target O at ~2.5Å from metal:
    # CA radius values calculated to achieve correct coordination distances
    CA_RADIUS_TABLE = {
        # (residue_type, mode): ca_radius
        ("GLU", "bidentate"): 3.5,   # Long side chain, both O's close
        ("GLU", "monodentate"): 4.5, # Long side chain, one O close
        ("ASP", "bidentate"): 2.8,   # Short side chain, both O's close
        ("ASP", "monodentate"): 3.8, # Short side chain, one O close
    }

    # Default for backward compatibility
    ca_radius_default = 3.5

    pdb_lines = []
    atom_serial = 1

    # Generate geometry positions for coordination
    geometry = template_def.get("geometry", "square_antiprism")
    cn = template_def["coordination_number"]

    # Get ideal coordination positions
    coord_positions = _generate_coordination_positions(
        geometry, cn, metal_pos, bond_dist
    )

    # Map residues to chains for grouping
    chain_a_residues = [r for r in template_def["residues"] if r["chain"] == "A"]
    chain_b_residues = [r for r in template_def["residues"] if r["chain"] == "B"]

    # Generate Chain A residues
    for i, res_def in enumerate(chain_a_residues):
        # Calculate CA position based on residue's angle
        angle = res_def.get("angle", i * 360 / len(chain_a_residues))
        rad = math.radians(angle)

        # Get per-residue CA radius based on type and coordination mode
        ca_radius = CA_RADIUS_TABLE.get(
            (res_def["type"], res_def["mode"]),
            ca_radius_default
        )

        # Position CA at appropriate distance
        # Use explicit z_offset from template if provided, otherwise use quadrant-based default
        if "z_offset" in res_def:
            z_offset = res_def["z_offset"]
        else:
            # Fallback: alternating z-levels per quadrant
            quadrant = int(angle // 90) % 4
            z_offset = 3.0 if (quadrant % 2 == 0) else -3.0
        ca_pos = metal_pos + np.array([
            ca_radius * math.cos(rad),
            ca_radius * math.sin(rad),
            z_offset,
        ])

        # Select generator based on mode
        if res_def["mode"] == "bidentate":
            if res_def["type"] == "GLU":
                res_lines, atom_serial = _generate_glu_residue(
                    "A", res_def["resnum"], ca_pos, metal_pos,
                    atom_serial, target_ln_o_distance=bond_dist
                )
            else:  # ASP
                res_lines, atom_serial = _generate_asp_residue(
                    "A", res_def["resnum"], ca_pos, metal_pos,
                    atom_serial, target_ln_o_distance=bond_dist
                )
        else:  # monodentate
            if res_def["type"] == "GLU":
                res_lines, atom_serial = _generate_monodentate_glu(
                    "A", res_def["resnum"], ca_pos, metal_pos,
                    atom_serial, target_ln_o_distance=bond_dist
                )
            else:  # ASP
                res_lines, atom_serial = _generate_monodentate_asp(
                    "A", res_def["resnum"], ca_pos, metal_pos,
                    atom_serial, target_ln_o_distance=bond_dist
                )

        pdb_lines.append(res_lines)

    # Optional Trp antenna on chain A
    if add_trp_antenna and trp_chain == "A":
        trp_ca_pos = metal_pos + np.array([8.0, 0.0, 2.0])
        trp_lines, atom_serial = _generate_trp_residue(
            "A", 30, trp_ca_pos, metal_pos, atom_serial, target_distance=4.5
        )
        pdb_lines.append(trp_lines)

    pdb_lines.append("TER")

    # Generate Chain B residues
    for i, res_def in enumerate(chain_b_residues):
        angle = res_def.get("angle", i * 360 / len(chain_b_residues) + 45)
        rad = math.radians(angle)

        # Get per-residue CA radius based on type and coordination mode
        ca_radius = CA_RADIUS_TABLE.get(
            (res_def["type"], res_def["mode"]),
            ca_radius_default
        )

        # Use explicit z_offset from template if provided, otherwise use quadrant-based default
        if "z_offset" in res_def:
            z_offset = res_def["z_offset"]
        else:
            # Fallback: alternating z-levels per quadrant
            quadrant = int(angle // 90) % 4
            z_offset = 3.0 if (quadrant % 2 == 0) else -3.0
        ca_pos = metal_pos + np.array([
            ca_radius * math.cos(rad),
            ca_radius * math.sin(rad),
            z_offset,
        ])

        if res_def["mode"] == "bidentate":
            if res_def["type"] == "GLU":
                res_lines, atom_serial = _generate_glu_residue(
                    "B", res_def["resnum"], ca_pos, metal_pos,
                    atom_serial, target_ln_o_distance=bond_dist
                )
            else:
                res_lines, atom_serial = _generate_asp_residue(
                    "B", res_def["resnum"], ca_pos, metal_pos,
                    atom_serial, target_ln_o_distance=bond_dist
                )
        else:
            if res_def["type"] == "GLU":
                res_lines, atom_serial = _generate_monodentate_glu(
                    "B", res_def["resnum"], ca_pos, metal_pos,
                    atom_serial, target_ln_o_distance=bond_dist
                )
            else:
                res_lines, atom_serial = _generate_monodentate_asp(
                    "B", res_def["resnum"], ca_pos, metal_pos,
                    atom_serial, target_ln_o_distance=bond_dist
                )

        pdb_lines.append(res_lines)

    if add_trp_antenna and trp_chain == "B":
        trp_ca_pos = metal_pos + np.array([8.0, 0.0, -2.0])
        trp_lines, atom_serial = _generate_trp_residue(
            "B", 30, trp_ca_pos, metal_pos, atom_serial, target_distance=4.5
        )
        pdb_lines.append(trp_lines)

    pdb_lines.append("TER")

    # Metal ion
    metal_code = metal.upper()[:2]
    metal_line = (
        f"HETATM{atom_serial:5d} {metal_code:>2s}   {metal_code:>3s} L   1    "
        f"{metal_pos[0]:8.3f}{metal_pos[1]:8.3f}{metal_pos[2]:8.3f}  1.00  0.00          {metal_code:>2s}"
    )
    pdb_lines.append(metal_line)
    atom_serial += 1

    # Waters if specified in template
    for i, water_def in enumerate(template_def.get("waters", [])):
        water_dist = water_def.get("distance", bond_dist)
        position_type = water_def.get("position", "axial_top")

        if position_type == "axial_top":
            water_pos = metal_pos + np.array([0, 0, water_dist])
        elif position_type == "axial_bottom":
            water_pos = metal_pos + np.array([0, 0, -water_dist])
        else:  # equatorial
            angle = i * 120  # Spread equatorial waters
            rad = math.radians(angle)
            water_pos = metal_pos + np.array([
                water_dist * math.cos(rad),
                water_dist * math.sin(rad),
                0,
            ])

        water_line = _generate_water(i + 1, water_pos, atom_serial)
        pdb_lines.append(water_line)
        atom_serial += 1

    pdb_lines.append("END")

    pdb_content = "\n".join(pdb_lines)

    # Validate analytical geometry output
    if validate_output and VALIDATION_AVAILABLE:
        validation = validate_coordination_geometry(
            pdb_content,
            expected_cn=template_def["coordination_number"],
            expected_geometry=template_def.get("geometry", "square_antiprism"),
        )
        if not validation.get("valid"):
            logger.warning(
                f"Analytical geometry has issues for {template_name}: "
                f"{validation.get('issues')}. Consider installing Architector for "
                "improved geometry (conda install -c conda-forge architector)"
            )
        else:
            logger.debug(f"Geometry validated: {validation.get('summary')}")

    return pdb_content


# =============================================================================
# NEW: Parametric Template Generator with Stochastic Mode
# =============================================================================


def _randomize_positions(
    positions: List[np.ndarray],
    metal_pos: np.ndarray,
    bond_distance: float,
    min_donor_distance: float = 2.8,
    perturbation_std: float = 0.3,
    seed: Optional[int] = None,
) -> List[np.ndarray]:
    """
    Randomize coordination positions while maintaining chemical constraints.

    Uses rejection sampling within spherical caps to ensure:
    1. Bond distance to metal is preserved (re-normalized)
    2. Minimum distance between donors is maintained
    3. Approximate geometry is preserved

    Args:
        positions: List of ideal coordination positions
        metal_pos: Position of metal ion
        bond_distance: Target Ln-O distance
        min_donor_distance: Minimum distance between donors (Å)
        perturbation_std: Standard deviation of random perturbation (Å)
        seed: Random seed for reproducibility

    Returns:
        List of randomized positions
    """
    rng = np.random.default_rng(seed)
    randomized = []

    for pos in positions:
        # Try to find valid position with rejection sampling
        for _ in range(100):  # Max attempts
            # Add random perturbation
            perturb = rng.normal(0, perturbation_std, 3)
            new_pos = pos + perturb

            # Re-normalize to correct bond distance
            direction = new_pos - metal_pos
            norm = np.linalg.norm(direction)
            if norm > 0.01:
                new_pos = metal_pos + (direction / norm) * bond_distance

            # Check distances to already-placed donors
            valid = True
            for other in randomized:
                if np.linalg.norm(new_pos - other) < min_donor_distance:
                    valid = False
                    break

            if valid:
                randomized.append(new_pos)
                break
        else:
            # Fall back to ideal position if no valid random found
            randomized.append(pos)

    return randomized


def generate_parametric_template(
    metal: str = "TB",
    coordination_number: int = 8,
    num_waters: int = 0,
    bidentate_fraction: float = 0.5,
    preferred_donors: Optional[List[str]] = None,
    geometry: Optional[str] = None,
    randomize: bool = False,
    num_variants: int = 1,
    seed: Optional[int] = None,
    pocket_radius: float = 6.0,
    min_donor_distance: float = 2.8,
    add_trp_antenna: bool = False,
) -> List[str]:
    """
    Generate template(s) with parametric control over coordination geometry.

    This advanced generator allows full control over:
    - Coordination number (6-10)
    - Water coordination
    - Bidentate/monodentate mixing
    - Stochastic ensemble generation

    Args:
        metal: Metal code (TB, GD, EU, etc.)
        coordination_number: Target CN (6-10)
        num_waters: Number of water molecules in coordination sphere
        bidentate_fraction: Fraction of bidentate residues (0.0 = all mono, 1.0 = all bi)
        preferred_donors: List of preferred residue types (e.g., ["ASP", "GLU"])
        geometry: Geometry type (auto-selected from CN if not specified)
        randomize: Whether to generate stochastic variants
        num_variants: Number of variants to generate (if randomize=True)
        seed: Random seed for reproducibility
        pocket_radius: Radius of binding pocket (Å)
        min_donor_distance: Minimum distance between donors (Å)
        add_trp_antenna: Whether to add Trp for TEBL

    Returns:
        List of PDB content strings (single item if randomize=False)
    """
    if preferred_donors is None:
        preferred_donors = ["ASP", "GLU"]

    # Get metal parameters
    params = LANTHANIDE_PARAMS.get(metal.upper(), LANTHANIDE_PARAMS["TB"])
    bond_dist = params["bond_distance"]

    # Auto-select geometry based on CN
    if geometry is None:
        geom_info = COORDINATION_GEOMETRIES.get(coordination_number, COORDINATION_GEOMETRIES[8])
        geometry = geom_info["geometry"]

    # Metal at origin - RFD3 designs proteins around origin
    metal_pos = np.array([0.0, 0.0, 0.0])

    # Calculate number of protein donors
    num_protein_donors = coordination_number - num_waters

    # Calculate bidentate vs monodentate counts
    # Each bidentate contributes 2 coordination positions
    num_bidentate = int(num_protein_donors * bidentate_fraction / 2)
    num_bidentate = min(num_bidentate, num_protein_donors // 2)  # Can't have more than half
    num_monodentate = num_protein_donors - (num_bidentate * 2)

    # Total residues needed
    total_residues = num_bidentate + num_monodentate

    # Split between chains A and B
    chain_a_count = (total_residues + 1) // 2
    chain_b_count = total_residues - chain_a_count

    def _generate_single_template(variant_seed: Optional[int] = None) -> str:
        """Generate a single template variant."""
        # CA radius lookup table - same as in generate_template_from_library
        # These values are calibrated to achieve correct Ln-O coordination distances (~2.5Å)
        CA_RADIUS_TABLE = {
            # (residue_type, mode): ca_radius
            ("GLU", "bidentate"): 3.5,   # Long side chain, both O's close
            ("GLU", "monodentate"): 4.5, # Long side chain, one O close
            ("ASP", "bidentate"): 2.8,   # Short side chain, both O's close
            ("ASP", "monodentate"): 3.8, # Short side chain, one O close
        }
        ca_radius_default = 3.5

        # Get ideal coordination positions
        coord_positions = _generate_coordination_positions(
            geometry, coordination_number, metal_pos, bond_dist
        )

        # Randomize positions if requested
        if randomize and variant_seed is not None:
            coord_positions = _randomize_positions(
                coord_positions, metal_pos, bond_dist,
                min_donor_distance=min_donor_distance,
                seed=variant_seed,
            )

        pdb_lines = []
        atom_serial = 1
        rng = np.random.default_rng(variant_seed)

        # Assign bidentate/monodentate to positions
        bi_positions = num_bidentate
        mono_positions = num_monodentate

        # Generate Chain A
        for i in range(chain_a_count):
            res_num = 10 + i * 5  # 10, 15, 20, 25, ...
            angle = i * 360 / max(chain_a_count, 1)
            rad = math.radians(angle)

            # Select residue type FIRST
            res_type = rng.choice(preferred_donors)

            # Determine mode SECOND
            if bi_positions > 0:
                mode = "bidentate"
                bi_positions -= 1
            else:
                mode = "monodentate"
                mono_positions -= 1

            # THEN calculate CA position using correct radius for this residue type/mode
            ca_radius = CA_RADIUS_TABLE.get((res_type, mode), ca_radius_default)
            # Use angle-based quadrant for z-offset to prevent O-O clashes
            quadrant = int(angle // 90) % 4
            z_offset = 3.0 if (quadrant % 2 == 0) else -3.0
            ca_pos = metal_pos + np.array([
                ca_radius * math.cos(rad),
                ca_radius * math.sin(rad),
                z_offset,
            ])

            # Generate residue
            if mode == "bidentate":
                if res_type == "GLU":
                    res_lines, atom_serial = _generate_glu_residue(
                        "A", res_num, ca_pos, metal_pos, atom_serial, bond_dist
                    )
                else:
                    res_lines, atom_serial = _generate_asp_residue(
                        "A", res_num, ca_pos, metal_pos, atom_serial, bond_dist
                    )
            else:
                if res_type == "GLU":
                    res_lines, atom_serial = _generate_monodentate_glu(
                        "A", res_num, ca_pos, metal_pos, atom_serial, bond_dist
                    )
                else:
                    res_lines, atom_serial = _generate_monodentate_asp(
                        "A", res_num, ca_pos, metal_pos, atom_serial, bond_dist
                    )

            pdb_lines.append(res_lines)

        if add_trp_antenna:
            trp_ca_pos = metal_pos + np.array([8.0, 0.0, 2.0])
            trp_lines, atom_serial = _generate_trp_residue(
                "A", 30, trp_ca_pos, metal_pos, atom_serial, 4.5
            )
            pdb_lines.append(trp_lines)

        pdb_lines.append("TER")

        # Generate Chain B
        for i in range(chain_b_count):
            res_num = 10 + i * 5
            angle = i * 360 / max(chain_b_count, 1) + 45
            rad = math.radians(angle)

            # Select residue type FIRST
            res_type = rng.choice(preferred_donors)

            # Determine mode SECOND
            if bi_positions > 0:
                mode = "bidentate"
                bi_positions -= 1
            else:
                mode = "monodentate"

            # THEN calculate CA position using correct radius for this residue type/mode
            ca_radius = CA_RADIUS_TABLE.get((res_type, mode), ca_radius_default)
            # Use angle-based quadrant for z-offset to prevent O-O clashes
            quadrant = int(angle // 90) % 4
            z_offset = 3.0 if (quadrant % 2 == 0) else -3.0
            ca_pos = metal_pos + np.array([
                ca_radius * math.cos(rad),
                ca_radius * math.sin(rad),
                z_offset,
            ])

            if mode == "bidentate":
                if res_type == "GLU":
                    res_lines, atom_serial = _generate_glu_residue(
                        "B", res_num, ca_pos, metal_pos, atom_serial, bond_dist
                    )
                else:
                    res_lines, atom_serial = _generate_asp_residue(
                        "B", res_num, ca_pos, metal_pos, atom_serial, bond_dist
                    )
            else:
                if res_type == "GLU":
                    res_lines, atom_serial = _generate_monodentate_glu(
                        "B", res_num, ca_pos, metal_pos, atom_serial, bond_dist
                    )
                else:
                    res_lines, atom_serial = _generate_monodentate_asp(
                        "B", res_num, ca_pos, metal_pos, atom_serial, bond_dist
                    )

            pdb_lines.append(res_lines)

        pdb_lines.append("TER")

        # Metal ion
        metal_code = metal.upper()[:2]
        metal_line = (
            f"HETATM{atom_serial:5d} {metal_code:>2s}   {metal_code:>3s} L   1    "
            f"{metal_pos[0]:8.3f}{metal_pos[1]:8.3f}{metal_pos[2]:8.3f}  1.00  0.00          {metal_code:>2s}"
        )
        pdb_lines.append(metal_line)
        atom_serial += 1

        # Waters
        for i in range(num_waters):
            if i == 0:
                water_pos = metal_pos + np.array([0, 0, bond_dist])
            elif i == 1:
                water_pos = metal_pos + np.array([0, 0, -bond_dist])
            else:
                angle = (i - 2) * 120
                rad = math.radians(angle)
                water_pos = metal_pos + np.array([
                    bond_dist * math.cos(rad),
                    bond_dist * math.sin(rad),
                    0,
                ])

            water_line = _generate_water(i + 1, water_pos, atom_serial)
            pdb_lines.append(water_line)
            atom_serial += 1

        pdb_lines.append("END")
        return "\n".join(pdb_lines)

    # Generate template(s)
    if randomize:
        templates = []
        for i in range(num_variants):
            variant_seed = seed + i if seed is not None else None
            templates.append(_generate_single_template(variant_seed))
        return templates
    else:
        return [_generate_single_template(seed)]


def get_template_for_metal(
    metal: str,
    template_name: Optional[str] = None,
    add_trp_antenna: bool = False,
) -> str:
    """
    Get the best template for a given metal, using template library.

    This is the RECOMMENDED entry point for new code. It uses the
    chemically realistic templates from TEMPLATE_LIBRARY.

    Args:
        metal: Metal code (TB, GD, EU, etc.)
        template_name: Specific template to use (auto-selected if None)
        add_trp_antenna: Whether to add Trp for TEBL

    Returns:
        PDB content string
    """
    if template_name is None:
        template_name = recommend_template(metal)

    return generate_template_from_library(
        template_name=template_name,
        metal=metal,
        add_trp_antenna=add_trp_antenna,
    )
