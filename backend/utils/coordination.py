"""
Metal Coordination Geometry Analyzer

Analyze metal binding sites in protein structures:
- Coordination number calculation
- Geometry classification (tetrahedral, octahedral, square antiprism, etc.)
- Donor atom type classification
- Bond distance and angle analysis
- Lanthanide optimization suggestions
"""

import math
import numpy as np
from typing import Dict, List, Any, Optional, Tuple, Set
from dataclasses import dataclass


# Metal ion properties database
METAL_PROPERTIES = {
    # Transition metals
    "FE": {
        "name": "Iron",
        "typical_coordination": [4, 6],
        "typical_geometry": ["tetrahedral", "octahedral"],
        "preferred_donors": ["S", "N", "O"],
        "ionic_radius": {"2+": 0.78, "3+": 0.64},
        "bond_distance_range": (1.9, 2.5),
    },
    "ZN": {
        "name": "Zinc",
        "typical_coordination": [4, 5, 6],
        "typical_geometry": ["tetrahedral", "trigonal_bipyramidal", "octahedral"],
        "preferred_donors": ["S", "N", "O"],
        "ionic_radius": {"2+": 0.74},
        "bond_distance_range": (1.9, 2.4),
    },
    "CA": {
        "name": "Calcium",
        "typical_coordination": [6, 7, 8],
        "typical_geometry": ["octahedral", "pentagonal_bipyramidal", "square_antiprism"],
        "preferred_donors": ["O"],
        "ionic_radius": {"2+": 1.00},
        "bond_distance_range": (2.3, 2.6),
    },
    "MG": {
        "name": "Magnesium",
        "typical_coordination": [6],
        "typical_geometry": ["octahedral"],
        "preferred_donors": ["O", "N"],
        "ionic_radius": {"2+": 0.72},
        "bond_distance_range": (2.0, 2.3),
    },
    "MN": {
        "name": "Manganese",
        "typical_coordination": [4, 6],
        "typical_geometry": ["tetrahedral", "octahedral"],
        "preferred_donors": ["O", "N"],
        "ionic_radius": {"2+": 0.83, "3+": 0.65},
        "bond_distance_range": (2.0, 2.4),
    },
    "CU": {
        "name": "Copper",
        "typical_coordination": [4, 5, 6],
        "typical_geometry": ["square_planar", "square_pyramidal", "octahedral"],
        "preferred_donors": ["N", "S", "O"],
        "ionic_radius": {"1+": 0.77, "2+": 0.73},
        "bond_distance_range": (1.9, 2.3),
    },
    "CO": {
        "name": "Cobalt",
        "typical_coordination": [4, 6],
        "typical_geometry": ["tetrahedral", "octahedral"],
        "preferred_donors": ["N", "O", "S"],
        "ionic_radius": {"2+": 0.75, "3+": 0.61},
        "bond_distance_range": (1.9, 2.3),
    },
    "NI": {
        "name": "Nickel",
        "typical_coordination": [4, 6],
        "typical_geometry": ["square_planar", "octahedral"],
        "preferred_donors": ["N", "S", "O"],
        "ionic_radius": {"2+": 0.69},
        "bond_distance_range": (1.9, 2.2),
    },
    # Lanthanides (high coordination, oxygen preference)
    "LA": {
        "name": "Lanthanum",
        "typical_coordination": [9, 10, 12],
        "typical_geometry": ["tricapped_trigonal_prism", "capped_square_antiprism"],
        "preferred_donors": ["O"],
        "ionic_radius": {"3+": 1.16},
        "bond_distance_range": (2.4, 2.8),
    },
    "TB": {
        "name": "Terbium",
        "typical_coordination": [8, 9],
        "typical_geometry": ["square_antiprism", "tricapped_trigonal_prism"],
        "preferred_donors": ["O"],
        "ionic_radius": {"3+": 1.04},
        "bond_distance_range": (2.3, 2.6),
    },
    "GD": {
        "name": "Gadolinium",
        "typical_coordination": [8, 9],
        "typical_geometry": ["square_antiprism", "tricapped_trigonal_prism"],
        "preferred_donors": ["O"],
        "ionic_radius": {"3+": 1.05},
        "bond_distance_range": (2.3, 2.6),
    },
    "EU": {
        "name": "Europium",
        "typical_coordination": [8, 9],
        "typical_geometry": ["square_antiprism", "bicapped_trigonal_prism"],
        "preferred_donors": ["O"],
        "ionic_radius": {"3+": 1.07},
        "bond_distance_range": (2.3, 2.6),
    },
    "CE": {
        "name": "Cerium",
        "typical_coordination": [8, 9, 10],
        "typical_geometry": ["square_antiprism", "bicapped_square_antiprism"],
        "preferred_donors": ["O"],
        "ionic_radius": {"3+": 1.14, "4+": 0.97},
        "bond_distance_range": (2.3, 2.7),
    },
    "SM": {
        "name": "Samarium",
        "typical_coordination": [8, 9],
        "typical_geometry": ["square_antiprism", "tricapped_trigonal_prism"],
        "preferred_donors": ["O"],
        "ionic_radius": {"3+": 1.08},
        "bond_distance_range": (2.3, 2.6),
    },
    "YB": {
        "name": "Ytterbium",
        "typical_coordination": [8],
        "typical_geometry": ["square_antiprism"],
        "preferred_donors": ["O"],
        "ionic_radius": {"3+": 0.99},
        "bond_distance_range": (2.2, 2.5),
    },
}

# Ideal geometry angles for classification
IDEAL_GEOMETRIES = {
    "tetrahedral": {
        "coordination": 4,
        "ideal_angles": [109.5] * 6,  # 6 L-M-L angles
        "angle_tolerance": 15.0,
    },
    "square_planar": {
        "coordination": 4,
        "ideal_angles": [90.0, 90.0, 90.0, 90.0, 180.0, 180.0],
        "angle_tolerance": 15.0,
    },
    "trigonal_bipyramidal": {
        "coordination": 5,
        "ideal_angles": [90.0, 90.0, 90.0, 120.0, 120.0, 120.0, 180.0],
        "angle_tolerance": 15.0,
    },
    "square_pyramidal": {
        "coordination": 5,
        "ideal_angles": [90.0, 90.0, 90.0, 90.0, 180.0],
        "angle_tolerance": 15.0,
    },
    "octahedral": {
        "coordination": 6,
        "ideal_angles": [90.0] * 12 + [180.0] * 3,
        "angle_tolerance": 15.0,
    },
    "pentagonal_bipyramidal": {
        "coordination": 7,
        "ideal_angles": [72.0] * 5 + [90.0] * 10 + [180.0],
        "angle_tolerance": 15.0,
    },
    "square_antiprism": {
        "coordination": 8,
        "ideal_angles": [52.4, 52.4, 52.4, 52.4, 52.4, 52.4, 52.4, 52.4],  # Simplified
        "angle_tolerance": 15.0,
    },
    "tricapped_trigonal_prism": {
        "coordination": 9,
        "ideal_angles": [70.5] * 9,  # Simplified
        "angle_tolerance": 18.0,
    },
}

# Donor atom classification
DONOR_TYPES = {
    # Oxygen donors
    ("ASP", "OD1"): "O_carboxylate",
    ("ASP", "OD2"): "O_carboxylate",
    ("GLU", "OE1"): "O_carboxylate",
    ("GLU", "OE2"): "O_carboxylate",
    ("ASN", "OD1"): "O_amide",
    ("GLN", "OE1"): "O_amide",
    ("SER", "OG"): "O_hydroxyl",
    ("THR", "OG1"): "O_hydroxyl",
    ("TYR", "OH"): "O_phenol",
    ("HOH", "O"): "O_water",
    ("WAT", "O"): "O_water",
    # Backbone oxygen
    ("*", "O"): "O_carbonyl",
    # Nitrogen donors
    ("HIS", "ND1"): "N_imidazole",
    ("HIS", "NE2"): "N_imidazole",
    ("LYS", "NZ"): "N_amine",
    ("ARG", "NH1"): "N_guanidinium",
    ("ARG", "NH2"): "N_guanidinium",
    ("TRP", "NE1"): "N_indole",
    # Backbone nitrogen
    ("*", "N"): "N_backbone",
    # Sulfur donors
    ("CYS", "SG"): "S_thiolate",
    ("MET", "SD"): "S_thioether",
}


@dataclass
class CoordinatingAtom:
    """Information about a coordinating atom."""
    chain: str
    residue_name: str
    residue_number: int
    atom_name: str
    element: str
    distance: float
    position: np.ndarray
    donor_type: str


def classify_donor_type(residue_name: str, atom_name: str) -> str:
    """
    Classify the donor atom type.

    Args:
        residue_name: 3-letter residue code
        atom_name: PDB atom name

    Returns:
        Donor type classification string
    """
    # Check specific residue-atom combination
    key = (residue_name.upper(), atom_name.strip())
    if key in DONOR_TYPES:
        return DONOR_TYPES[key]

    # Check wildcard (backbone atoms)
    wildcard_key = ("*", atom_name.strip())
    if wildcard_key in DONOR_TYPES:
        return DONOR_TYPES[wildcard_key]

    # Classify by element
    element = atom_name.strip()[0].upper()
    if element == "O":
        return "O_other"
    elif element == "N":
        return "N_other"
    elif element == "S":
        return "S_other"
    else:
        return f"{element}_other"


def get_coordinating_atoms(
    pdb_content: str,
    metal_chain: str,
    metal_residue: str,
    metal_resnum: int,
    distance_cutoff: float = 3.0
) -> List[CoordinatingAtom]:
    """
    Find all atoms coordinating to a metal ion.

    Args:
        pdb_content: PDB file content
        metal_chain: Chain ID of metal
        metal_residue: Residue name of metal (e.g., "FE", "ZN")
        metal_resnum: Residue number of metal
        distance_cutoff: Maximum distance for coordination (Angstroms)

    Returns:
        List of CoordinatingAtom objects
    """
    # Parse metal position
    metal_pos = None
    protein_atoms = []

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
            element = line[76:78].strip() if len(line) > 76 else atom_name[0]

            pos = np.array([x, y, z])

            # Check if this is the metal
            if (chain == metal_chain and
                res_name.upper() == metal_residue.upper() and
                res_num == metal_resnum):
                metal_pos = pos
            else:
                protein_atoms.append({
                    "chain": chain,
                    "res_name": res_name,
                    "res_num": res_num,
                    "atom_name": atom_name,
                    "element": element,
                    "position": pos,
                })
        except (ValueError, IndexError):
            continue

    if metal_pos is None:
        return []

    # Find coordinating atoms
    coordinating = []
    for atom in protein_atoms:
        distance = np.linalg.norm(atom["position"] - metal_pos)
        if distance <= distance_cutoff:
            donor_type = classify_donor_type(atom["res_name"], atom["atom_name"])
            coordinating.append(CoordinatingAtom(
                chain=atom["chain"],
                residue_name=atom["res_name"],
                residue_number=atom["res_num"],
                atom_name=atom["atom_name"],
                element=atom["element"],
                distance=distance,
                position=atom["position"],
                donor_type=donor_type,
            ))

    # Sort by distance
    coordinating.sort(key=lambda x: x.distance)
    return coordinating


def calculate_coordination_number(
    pdb_content: str,
    metal_chain: str,
    metal_residue: str,
    metal_resnum: int,
    distance_cutoff: float = 3.0
) -> int:
    """
    Calculate the coordination number of a metal ion.

    Args:
        pdb_content: PDB file content
        metal_chain: Chain ID of metal
        metal_residue: Residue name of metal
        metal_resnum: Residue number of metal
        distance_cutoff: Maximum distance for coordination

    Returns:
        Coordination number
    """
    atoms = get_coordinating_atoms(
        pdb_content, metal_chain, metal_residue, metal_resnum, distance_cutoff
    )
    return len(atoms)


def calculate_lml_angles(metal_pos: np.ndarray, ligand_positions: List[np.ndarray]) -> List[float]:
    """
    Calculate all ligand-metal-ligand angles.

    Args:
        metal_pos: Position of metal center
        ligand_positions: List of ligand positions

    Returns:
        List of L-M-L angles in degrees
    """
    n = len(ligand_positions)
    if n < 2:
        return []

    angles = []
    for i in range(n):
        for j in range(i + 1, n):
            # Vectors from metal to ligands
            v1 = ligand_positions[i] - metal_pos
            v2 = ligand_positions[j] - metal_pos

            # Calculate angle
            cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
            cos_angle = np.clip(cos_angle, -1.0, 1.0)
            angle = np.degrees(np.arccos(cos_angle))
            angles.append(angle)

    return angles


def classify_geometry(
    coordination_number: int,
    angles: List[float]
) -> Tuple[str, float]:
    """
    Classify coordination geometry based on angles.

    Args:
        coordination_number: Number of coordinating atoms
        angles: List of L-M-L angles

    Returns:
        Tuple of (geometry_name, rmsd_from_ideal)
    """
    if coordination_number < 4:
        return "underdetermined", 0.0

    best_geometry = "unknown"
    best_rmsd = float('inf')

    for geom_name, geom_info in IDEAL_GEOMETRIES.items():
        if geom_info["coordination"] != coordination_number:
            continue

        # Calculate RMSD from ideal angles
        ideal = geom_info["ideal_angles"]
        if len(angles) != len(ideal):
            continue

        # Sort both for comparison
        sorted_angles = sorted(angles)
        sorted_ideal = sorted(ideal)

        rmsd = np.sqrt(np.mean([(a - i) ** 2 for a, i in zip(sorted_angles, sorted_ideal)]))

        if rmsd < best_rmsd:
            best_rmsd = rmsd
            best_geometry = geom_name

    return best_geometry, best_rmsd


def analyze_coordination_geometry(
    pdb_content: str,
    metal_chain: str,
    metal_residue: str,
    metal_resnum: int,
    distance_cutoff: float = 3.0
) -> Dict[str, Any]:
    """
    Comprehensive analysis of metal coordination geometry.

    Args:
        pdb_content: PDB file content
        metal_chain: Chain ID of metal
        metal_residue: Residue name of metal (e.g., "FE", "TB")
        metal_resnum: Residue number of metal
        distance_cutoff: Maximum distance for coordination

    Returns:
        Dict with comprehensive coordination analysis
    """
    # Get coordinating atoms
    coordinating_atoms = get_coordinating_atoms(
        pdb_content, metal_chain, metal_residue, metal_resnum, distance_cutoff
    )

    if not coordinating_atoms:
        return {
            "success": False,
            "error": f"No metal found at {metal_chain}:{metal_residue}{metal_resnum}",
        }

    # Get metal position (from first coordinating atom's distance)
    # We need to recalculate metal position
    metal_pos = None
    for line in pdb_content.split('\n'):
        if line.startswith('HETATM') or line.startswith('ATOM'):
            try:
                chain = line[21]
                res_name = line[17:20].strip()
                res_num = int(line[22:26].strip())
                if (chain == metal_chain and
                    res_name.upper() == metal_residue.upper() and
                    res_num == metal_resnum):
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    metal_pos = np.array([x, y, z])
                    break
            except (ValueError, IndexError):
                continue

    # Coordination number
    coordination_number = len(coordinating_atoms)

    # Calculate L-M-L angles
    ligand_positions = [atom.position for atom in coordinating_atoms]
    angles = calculate_lml_angles(metal_pos, ligand_positions)

    # Classify geometry
    geometry_type, geometry_rmsd = classify_geometry(coordination_number, angles)

    # Analyze donor types
    donor_types = [atom.donor_type for atom in coordinating_atoms]
    donor_summary = {}
    for dt in donor_types:
        donor_summary[dt] = donor_summary.get(dt, 0) + 1

    # Bond distances statistics
    distances = [atom.distance for atom in coordinating_atoms]
    avg_distance = np.mean(distances)
    min_distance = np.min(distances)
    max_distance = np.max(distances)

    # Get metal properties
    metal_elem = metal_residue.upper()
    metal_props = METAL_PROPERTIES.get(metal_elem, {})

    # Generate suggestions for improvement
    suggestions = []

    # Check if coordination matches expected for this metal
    typical_coord = metal_props.get("typical_coordination", [])
    if typical_coord and coordination_number not in typical_coord:
        if coordination_number < min(typical_coord):
            suggestions.append(
                f"Low coordination ({coordination_number}). "
                f"Typical for {metal_elem}: {typical_coord}. "
                "Consider adding more coordinating residues."
            )
        elif coordination_number > max(typical_coord):
            suggestions.append(
                f"High coordination ({coordination_number}). "
                f"Typical for {metal_elem}: {typical_coord}."
            )

    # Check donor type preferences
    preferred_donors = metal_props.get("preferred_donors", [])
    if preferred_donors:
        for dt, count in donor_summary.items():
            donor_element = dt.split("_")[0]
            if donor_element not in preferred_donors:
                suggestions.append(
                    f"Donor type '{dt}' may not be optimal for {metal_elem}. "
                    f"Preferred donors: {preferred_donors}"
                )

    # Check bond distances
    expected_range = metal_props.get("bond_distance_range", (2.0, 3.0))
    if avg_distance < expected_range[0] or avg_distance > expected_range[1]:
        suggestions.append(
            f"Average bond distance ({avg_distance:.2f}Å) outside typical range "
            f"({expected_range[0]}-{expected_range[1]}Å) for {metal_elem}."
        )

    # Lanthanide-specific suggestions
    if metal_elem in ["LA", "TB", "GD", "EU", "CE", "SM", "YB"]:
        # Check for carboxylate donors
        carboxylate_count = sum(1 for dt in donor_types if "carboxylate" in dt)
        if carboxylate_count < 4:
            suggestions.append(
                f"Lanthanides prefer carboxylate donors. "
                f"Current carboxylate count: {carboxylate_count}. "
                "Consider adding more Asp/Glu residues."
            )

        # Check for sulfur donors (not preferred for lanthanides)
        sulfur_count = sum(1 for dt in donor_types if dt.startswith("S_"))
        if sulfur_count > 0:
            suggestions.append(
                f"Lanthanides prefer O-donors over S-donors. "
                f"Found {sulfur_count} sulfur donors. "
                "Consider replacing Cys/Met with Asp/Glu/Asn."
            )

    # Format coordinating residues for output
    coord_residues = []
    for atom in coordinating_atoms:
        coord_residues.append({
            "chain": atom.chain,
            "residue": atom.residue_name,
            "residue_number": atom.residue_number,
            "atom": atom.atom_name,
            "distance": round(atom.distance, 2),
            "donor_type": atom.donor_type,
        })

    return {
        "success": True,
        "metal": {
            "chain": metal_chain,
            "residue": metal_residue,
            "residue_number": metal_resnum,
            "element": metal_elem,
            "properties": metal_props,
        },
        "coordination": {
            "number": coordination_number,
            "geometry": geometry_type,
            "geometry_rmsd": round(geometry_rmsd, 2),
            "coordinating_atoms": coord_residues,
        },
        "bond_analysis": {
            "average_distance": round(avg_distance, 2),
            "min_distance": round(min_distance, 2),
            "max_distance": round(max_distance, 2),
            "distances": [round(d, 2) for d in distances],
        },
        "donor_analysis": {
            "types": donor_summary,
            "donor_list": donor_types,
        },
        "angles": {
            "lml_angles": [round(a, 1) for a in angles],
        },
        "suggestions": suggestions,
    }


def suggest_lanthanide_conversion(
    current_analysis: Dict[str, Any],
    target_lanthanide: str = "TB"
) -> Dict[str, Any]:
    """
    Suggest modifications for converting a metal site to bind lanthanides.

    Args:
        current_analysis: Output from analyze_coordination_geometry
        target_lanthanide: Target lanthanide element (e.g., "TB", "GD", "EU")

    Returns:
        Dict with conversion suggestions and RFD3 parameters
    """
    if not current_analysis.get("success"):
        return current_analysis

    target_props = METAL_PROPERTIES.get(target_lanthanide.upper(), {})
    target_coordination = target_props.get("typical_coordination", [8, 9])

    current_coord = current_analysis["coordination"]["number"]
    coord_delta = min(target_coordination) - current_coord

    suggestions = []
    rfd3_params = {}

    # Analyze current donor types
    donor_analysis = current_analysis["donor_analysis"]
    donor_types = donor_analysis.get("donor_list", [])

    # Count sulfur donors (need to be replaced)
    sulfur_donors = [d for d in donor_types if d.startswith("S_")]
    oxygen_donors = [d for d in donor_types if d.startswith("O_")]

    if sulfur_donors:
        suggestions.append({
            "type": "replace_donors",
            "description": f"Replace {len(sulfur_donors)} sulfur donors with oxygen donors (Asp/Glu)",
            "residues_to_modify": [
                atom for atom in current_analysis["coordination"]["coordinating_atoms"]
                if atom["donor_type"].startswith("S_")
            ],
        })

    if coord_delta > 0:
        suggestions.append({
            "type": "increase_coordination",
            "description": f"Increase coordination from {current_coord} to {min(target_coordination)}-{max(target_coordination)}",
            "needed_donors": coord_delta,
            "recommendation": "Add Asp/Glu residues or backbone carbonyl interactions",
        })

    # Generate RFD3 parameters
    coord_residues = current_analysis["coordination"]["coordinating_atoms"]
    residue_ids = [f"A{atom['residue_number']}" for atom in coord_residues]

    rfd3_params = {
        "ligand": target_lanthanide.upper(),
        "partial_t": 12.0 if coord_delta <= 2 else 15.0,
        "unindex": ",".join(residue_ids),
        "select_fixed_atoms": {
            res_id: "BKBN" for res_id in residue_ids  # Fix backbone, redesign side chains
        },
        "num_timesteps": 200,
        "step_scale": 1.5,
        "gamma_0": 0.6,
        "num_designs": 10,
    }

    # Suggest evaluation criteria
    evaluation_criteria = {
        "target_coordination": target_coordination,
        "target_distance_range": target_props.get("bond_distance_range", (2.3, 2.6)),
        "preferred_geometry": target_props.get("typical_geometry", ["square_antiprism"]),
        "preferred_donors": ["O_carboxylate", "O_carbonyl", "O_amide"],
        "min_carboxylate_donors": 4,
    }

    return {
        "success": True,
        "current_metal": current_analysis["metal"]["element"],
        "target_metal": target_lanthanide.upper(),
        "target_properties": target_props,
        "conversion_suggestions": suggestions,
        "rfd3_parameters": rfd3_params,
        "evaluation_criteria": evaluation_criteria,
    }
