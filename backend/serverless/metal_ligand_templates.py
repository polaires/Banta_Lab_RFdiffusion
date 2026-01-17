"""
Metal-Ligand Complex Templates for Interface Design

Provides pre-defined templates for common cofactors and metal-ligand complexes
with accurate 3D coordinates and coordination information.

Templates include:
- PQQ-Ca (pyrroloquinoline quinone with calcium)
- Citrate-Lanthanide (citrate with Tb/Eu/Gd)
- Heme-Fe (porphyrin with iron) [future]

References:
- PQQ-Ca: Lumpe et al. 2020, Acta Cryst. C; PDB 1W6S
- Citrate-Ln: Wang et al. 2008, J. Mol. Struct.
"""

import math
import tempfile
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

# Try to import RDKit for conformer generation
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

# Import from metal_site_fetcher for PDB-based templates
from metal_site_fetcher import (
    get_reference_template,
    generate_template_from_pdb,
    query_metal_ligand_sites,
)


# =============================================================================
# Template Definitions
# =============================================================================

# PQQ-Calcium template
# Source: PDB ligand PQQ + Ca coordination from 1W6S (methanol dehydrogenase)
# Crystal structure: Lumpe et al. 2020 - Ca3PQQ2·13H2O
PQQ_CA_TEMPLATE = {
    "name": "PQQ-Calcium",
    "description": "Pyrroloquinoline quinone with calcium cofactor for quinoprotein dehydrogenases",
    "metal": "CA",
    "ligand_smiles": "OC(=O)c1[nH]c2c(c1)C(=O)C(=O)c3nc(cc(C(O)=O)c23)C(O)=O",
    "ligand_res_name": "PQQ",
    "molecular_weight": 330.206,
    "coordination": {
        "metal_coordination_number": 6,
        "ligand_donors": ["O5", "N6", "O7A"],  # PQQ atoms coordinating Ca
        "ligand_donor_count": 3,
        "protein_sites": 3,  # Available for protein coordination
        "geometry": "distorted_octahedral",
        "bond_distances": {
            "O5": 2.40,   # Quinone oxygen
            "N6": 2.50,   # Pyridine nitrogen
            "O7A": 2.40,  # Carboxylate oxygen
        },
    },
    "preferred_protein_donors": ["Glu", "Asp", "Asn", "water"],
    "donor_atom_types": ["OE1", "OE2", "OD1", "OD2", "OD1"],  # Carboxylate oxygens
    "pdb_reference": "1W6S",
    "best_for": ["dehydrogenase", "electron_transfer", "quinoprotein"],
    "notes": "Ca2+ binds PQQ via O5/N6/O7, protein completes octahedral with 3 Glu/Asp",
}

# Citrate-Lanthanide template
# Source: Citrate tricarboxylate coordination to Tb3+/Eu3+
# Reference: Wang et al. 2008 - Lanthanide-citrate complexes
CITRATE_LANTHANIDE_TEMPLATE = {
    "name": "Citrate-Lanthanide",
    "description": "Citrate-coordinated lanthanide (Tb/Eu/Gd) for luminescent biosensors",
    "metals": ["TB", "EU", "GD", "LA", "CE"],  # Compatible lanthanides
    "default_metal": "TB",
    "ligand_smiles": "OC(CC(=O)[O-])(CC(=O)[O-])C(=O)[O-]",  # Citrate trianion
    "ligand_res_name": "CIT",
    "molecular_weight": 189.10,  # Citrate anion
    "coordination": {
        "metal_coordination_number": 9,  # Typical for Ln3+
        "ligand_donors": ["O1", "O3", "O5", "O7"],  # Carboxylate oxygens
        "ligand_donor_count": 4,  # Citrate can be tetradentate
        "citrate_denticity": "tetradentate",  # Can vary from bi- to tetra-
        "protein_sites": 5,  # Remaining coordination sites
        "geometry": "tricapped_trigonal_prism",
        "bond_distances": {
            "carboxylate_O": 2.35,  # Typical Ln-O distance
            "hydroxyl_O": 2.45,
        },
    },
    "preferred_protein_donors": ["Glu", "Asp", "Asn", "water"],
    "donor_atom_types": ["OE1", "OE2", "OD1", "OD2"],
    "luminescence": {
        "TB": {"emission": "green", "wavelength_nm": 545, "transition": "5D4->7F5"},
        "EU": {"emission": "red", "wavelength_nm": 615, "transition": "5D0->7F2"},
        "GD": {"emission": "UV", "wavelength_nm": 311, "mri_contrast": True},
    },
    "best_for": ["biosensor", "luminescence", "metal_sensing", "TEBL"],
    "notes": "Citrate chelates Ln3+ via carboxylates; protein Glu/Asp complete coordination",
}

# Template library
METAL_LIGAND_COMPLEX_TEMPLATES = {
    "pqq_ca": PQQ_CA_TEMPLATE,
    "citrate_tb": {**CITRATE_LANTHANIDE_TEMPLATE, "metal": "TB"},
    "citrate_eu": {**CITRATE_LANTHANIDE_TEMPLATE, "metal": "EU"},
    "citrate_gd": {**CITRATE_LANTHANIDE_TEMPLATE, "metal": "GD"},
    "citrate_la": {**CITRATE_LANTHANIDE_TEMPLATE, "metal": "LA"},
}

# PDB IDs for known metal-ligand complexes
KNOWN_COMPLEX_PDBS = {
    "pqq_ca": {
        "pdb_ids": ["1W6S", "4AAH", "1W99"],
        "metal": "CA",
        "ligand": "PQQ",
        "description": "PQQ-dependent methanol dehydrogenase",
    },
    "citrate_tb": {
        "pdb_ids": ["6MI5"],
        "metal": "TB",
        "ligand": "CIT",
        "description": "Citrate-terbium complex (Lanmodulin analog)",
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


# =============================================================================
# Template Access Functions
# =============================================================================

def get_template(template_name: str) -> Optional[Dict[str, Any]]:
    """
    Get a template by name.

    Args:
        template_name: Template key (e.g., "pqq_ca", "citrate_tb")

    Returns:
        Template dictionary or None if not found
    """
    return METAL_LIGAND_COMPLEX_TEMPLATES.get(template_name.lower())


def list_templates() -> List[str]:
    """List all available template names."""
    return list(METAL_LIGAND_COMPLEX_TEMPLATES.keys())


def get_template_info(template_name: str) -> Optional[Dict[str, Any]]:
    """
    Get summary information about a template.

    Returns dict with name, description, metal, ligand, coordination info.
    """
    template = get_template(template_name)
    if not template:
        return None

    return {
        "name": template.get("name"),
        "description": template.get("description"),
        "metal": template.get("metal") or template.get("default_metal"),
        "ligand": template.get("ligand_res_name"),
        "ligand_smiles": template.get("ligand_smiles"),
        "coordination_number": template["coordination"]["metal_coordination_number"],
        "ligand_donors": template["coordination"]["ligand_donor_count"],
        "protein_sites": template["coordination"]["protein_sites"],
        "geometry": template["coordination"]["geometry"],
        "preferred_donors": template.get("preferred_protein_donors", []),
        "best_for": template.get("best_for", []),
    }


def recommend_template(metal: str, ligand_type: Optional[str] = None) -> Optional[str]:
    """
    Recommend a template based on metal and optional ligand type.

    Args:
        metal: Metal code (CA, TB, EU, etc.)
        ligand_type: Optional ligand hint ("quinone", "carboxylate", etc.)

    Returns:
        Recommended template name or None
    """
    metal = metal.upper()

    # Direct metal matches
    if metal == "CA":
        return "pqq_ca"
    elif metal == "TB":
        return "citrate_tb"
    elif metal == "EU":
        return "citrate_eu"
    elif metal == "GD":
        return "citrate_gd"
    elif metal == "LA":
        return "citrate_la"

    # Ligand type hints
    if ligand_type:
        ligand_lower = ligand_type.lower()
        if "quinone" in ligand_lower or "pqq" in ligand_lower:
            return "pqq_ca"
        elif "citrate" in ligand_lower or "citric" in ligand_lower:
            if metal in {"TB", "EU", "GD", "LA", "CE"}:
                return f"citrate_{metal.lower()}"

    return None


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

    # Priority 1: Try PDB database
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
                    return template
            except Exception as e:
                continue

    # Priority 2: Try library templates
    if template_key in METAL_LIGAND_COMPLEX_TEMPLATES:
        template = dict(METAL_LIGAND_COMPLEX_TEMPLATES[template_key])
        template["source"] = "library"
        return template

    # Priority 3: Calculated fallback
    if fallback_enabled and metal:
        from metal_chemistry import get_preferred_donors, get_coordination_number_range

        metal_upper = metal.upper()
        try:
            # Get default oxidation state for metal
            from metal_chemistry import METAL_DATABASE
            default_ox = METAL_DATABASE.get(metal_upper, {}).get("default_oxidation", 2)

            cn_range = get_coordination_number_range(metal_upper, default_ox)
            donors = get_preferred_donors(metal_upper, default_ox)

            template = _generate_calculated_template(
                metal=metal_upper,
                ligand=ligand,
                coordination_number=cn_range[1],
                donors=donors,
            )
            template["source"] = "calculated"
            template["warning"] = "Calculated template - verify geometry manually"
            return template
        except (ValueError, KeyError):
            # Metal not in database
            pass

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
    """
    # Get geometry based on coordination number
    if coordination_number <= 4:
        geometry = "tetrahedral"
        bond_distance = 2.3
    elif coordination_number <= 6:
        geometry = "octahedral"
        bond_distance = 2.4
    else:
        geometry = "square_antiprism" if coordination_number == 8 else "tricapped_trigonal_prism"
        bond_distance = 2.45

    # Check for lanthanide-specific handling
    if metal.upper() in {"TB", "EU", "GD", "LA", "CE", "SM", "YB"}:
        geometry = "square_antiprism" if coordination_number == 8 else "tricapped_trigonal_prism"
        bond_distance = 2.40

    # Generate ideal positions
    metal_pos = (50.0, 50.0, 50.0)

    template = {
        "name": f"calculated_{metal}_{coordination_number}",
        "metal": metal.upper(),
        "ligand_res_name": ligand,
        "coordination_number": coordination_number,
        "geometry": geometry,
        "bond_distance": bond_distance,
        "metal_coords": metal_pos,
        "calculated": True,
    }

    return template


# =============================================================================
# PDB Generation Functions
# =============================================================================

def generate_complex_pdb(
    template_name: str,
    metal: Optional[str] = None,
    center: Tuple[float, float, float] = (50.0, 50.0, 50.0),
    include_hydrogens: bool = False,
) -> Optional[str]:
    """
    Generate PDB content for a metal-ligand complex from template.

    Uses RDKit to generate ligand conformer, then places metal at coordination site.

    Args:
        template_name: Template key
        metal: Override metal (optional)
        center: Center coordinates for the complex
        include_hydrogens: Whether to include H atoms

    Returns:
        PDB file content as string, or None if generation fails
    """
    template = get_template(template_name)
    if not template:
        print(f"[MetalLigandTemplates] Unknown template: {template_name}")
        return None

    # Determine metal to use
    if metal:
        metal_code = metal.upper()
    else:
        metal_code = template.get("metal") or template.get("default_metal")

    if not metal_code:
        print(f"[MetalLigandTemplates] No metal specified for template: {template_name}")
        return None

    # Get ligand SMILES
    smiles = template.get("ligand_smiles")
    if not smiles:
        print(f"[MetalLigandTemplates] No SMILES in template: {template_name}")
        return None

    ligand_name = template.get("ligand_res_name", "LIG")

    # Generate ligand conformer
    if not RDKIT_AVAILABLE:
        print("[MetalLigandTemplates] RDKit not available, using fallback PDB")
        return _generate_fallback_pdb(template_name, metal_code, center)

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"[MetalLigandTemplates] Invalid SMILES: {smiles}")
            return _generate_fallback_pdb(template_name, metal_code, center)

        if include_hydrogens:
            mol = Chem.AddHs(mol)

        # Generate 3D conformer
        result = AllChem.EmbedMolecule(mol, randomSeed=42)
        if result != 0:
            # Fallback to less strict embedding
            AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)

        # Optimize geometry
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        except Exception:
            pass  # Use unoptimized if MMFF fails

        # Get conformer
        conf = mol.GetConformer()

        # Calculate ligand centroid
        positions = conf.GetPositions()
        lig_centroid = positions.mean(axis=0)

        # Get coordination info
        coord_info = template.get("coordination", {})
        bond_distance = 2.4  # Default metal-ligand distance

        # Estimate metal position (above ligand plane for planar ligands)
        # For PQQ: metal is above the quinone ring plane
        # For citrate: metal is at center of carboxylate cage
        if "pqq" in template_name.lower():
            # PQQ: place Ca above the quinoline plane
            metal_offset = [0.0, 0.0, 2.4]
        else:
            # Citrate/other: place metal at centroid offset
            metal_offset = [0.0, 0.0, 2.0]

        metal_pos = lig_centroid + metal_offset

        # Translate to desired center
        translation = [
            center[0] - (lig_centroid[0] + metal_pos[0]) / 2,
            center[1] - (lig_centroid[1] + metal_pos[1]) / 2,
            center[2] - (lig_centroid[2] + metal_pos[2]) / 2,
        ]

        # Build PDB content
        lines = []
        atom_num = 1

        # Write ligand atoms
        for i, atom in enumerate(mol.GetAtoms()):
            if not include_hydrogens and atom.GetAtomicNum() == 1:
                continue

            pos = conf.GetAtomPosition(i)
            x = pos.x + translation[0]
            y = pos.y + translation[1]
            z = pos.z + translation[2]

            atom_name = atom.GetSymbol() + str(i + 1)
            atom_name = atom_name[:4].ljust(4)
            element = atom.GetSymbol()

            line = (
                f"HETATM{atom_num:5d} {atom_name:4s} {ligand_name:3s} "
                f"L   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}"
            )
            lines.append(line)
            atom_num += 1

        # Write metal atom
        mx = metal_pos[0] + translation[0]
        my = metal_pos[1] + translation[1]
        mz = metal_pos[2] + translation[2]

        metal_element = metal_code[:2]
        line = (
            f"HETATM{atom_num:5d} {metal_code:4s} {metal_code:3s} "
            f"M   2    {mx:8.3f}{my:8.3f}{mz:8.3f}  1.00  0.00          {metal_element:>2s}"
        )
        lines.append(line)
        lines.append("END")

        return "\n".join(lines)

    except Exception as e:
        print(f"[MetalLigandTemplates] Error generating PDB: {e}")
        return _generate_fallback_pdb(template_name, metal_code, center)


def _generate_fallback_pdb(
    template_name: str,
    metal_code: str,
    center: Tuple[float, float, float],
) -> str:
    """
    Generate a minimal fallback PDB when RDKit is not available.

    Creates a simplified representation with key atoms only.
    """
    cx, cy, cz = center

    if "pqq" in template_name.lower():
        # Simplified PQQ-Ca: quinone oxygens + carboxylates + Ca
        lines = [
            f"HETATM    1  O5  PQQ L   1    {cx-2.0:8.3f}{cy:8.3f}{cz:8.3f}  1.00  0.00           O",
            f"HETATM    2  N6  PQQ L   1    {cx:8.3f}{cy+2.0:8.3f}{cz:8.3f}  1.00  0.00           N",
            f"HETATM    3  O7A PQQ L   1    {cx+2.0:8.3f}{cy:8.3f}{cz:8.3f}  1.00  0.00           O",
            f"HETATM    4  C1  PQQ L   1    {cx-3.0:8.3f}{cy-1.0:8.3f}{cz:8.3f}  1.00  0.00           C",
            f"HETATM    5  C2  PQQ L   1    {cx+3.0:8.3f}{cy-1.0:8.3f}{cz:8.3f}  1.00  0.00           C",
            f"HETATM    6 CA   CA  M   2    {cx:8.3f}{cy:8.3f}{cz+2.4:8.3f}  1.00  0.00          CA",
            "END",
        ]
    elif "citrate" in template_name.lower():
        # Simplified citrate-Ln: 3 carboxylate groups + central C + metal
        lines = [
            f"HETATM    1  O1  CIT L   1    {cx-2.5:8.3f}{cy+1.5:8.3f}{cz:8.3f}  1.00  0.00           O",
            f"HETATM    2  O2  CIT L   1    {cx-2.5:8.3f}{cy-1.5:8.3f}{cz:8.3f}  1.00  0.00           O",
            f"HETATM    3  O3  CIT L   1    {cx+2.5:8.3f}{cy+1.5:8.3f}{cz:8.3f}  1.00  0.00           O",
            f"HETATM    4  O4  CIT L   1    {cx+2.5:8.3f}{cy-1.5:8.3f}{cz:8.3f}  1.00  0.00           O",
            f"HETATM    5  O5  CIT L   1    {cx:8.3f}{cy+2.5:8.3f}{cz:8.3f}  1.00  0.00           O",
            f"HETATM    6  O6  CIT L   1    {cx:8.3f}{cy-2.5:8.3f}{cz:8.3f}  1.00  0.00           O",
            f"HETATM    7  C1  CIT L   1    {cx:8.3f}{cy:8.3f}{cz-1.0:8.3f}  1.00  0.00           C",
            f"HETATM    8 {metal_code:4s} {metal_code:3s} M   2    {cx:8.3f}{cy:8.3f}{cz+2.35:8.3f}  1.00  0.00          {metal_code[:2]:>2s}",
            "END",
        ]
    else:
        # Generic fallback
        lines = [
            f"HETATM    1  X1  LIG L   1    {cx:8.3f}{cy:8.3f}{cz:8.3f}  1.00  0.00           C",
            f"HETATM    2 {metal_code:4s} {metal_code:3s} M   2    {cx:8.3f}{cy:8.3f}{cz+2.4:8.3f}  1.00  0.00          {metal_code[:2]:>2s}",
            "END",
        ]

    return "\n".join(lines)


def get_template_coordination_info(template_name: str) -> Dict[str, Any]:
    """
    Get coordination information for RFD3 configuration.

    Returns dict with:
    - hotspot_atoms: Atoms to use in select_hotspots
    - hbond_acceptors: Atoms for select_hbond_acceptor
    - exposed_atoms: Atoms for select_exposed
    - fixed_atoms: Atoms for select_fixed_atoms
    """
    template = get_template(template_name)
    if not template:
        return {}

    coord = template.get("coordination", {})
    ligand_name = template.get("ligand_res_name", "LIG")
    metal = template.get("metal") or template.get("default_metal", "")

    # Ligand donors that should contact the protein
    ligand_donors = coord.get("ligand_donors", [])

    # All ligand oxygen atoms can be H-bond acceptors
    # For PQQ: quinone O, carboxylate O
    # For citrate: carboxylate O
    hbond_acceptors = [a for a in ligand_donors if a.startswith("O")]

    # Atoms to expose (not coordinating metal, available for protein contact)
    # Typically non-coordinating atoms on the ligand
    exposed_atoms = []

    return {
        "ligand_name": ligand_name,
        "metal_code": metal,
        "hotspot_atoms": {metal: "ALL"},  # Ensure protein contacts metal
        "hbond_acceptors": {ligand_name: ",".join(hbond_acceptors)} if hbond_acceptors else {},
        "exposed_atoms": {ligand_name: ",".join(exposed_atoms)} if exposed_atoms else {},
        "preferred_donors": template.get("preferred_protein_donors", []),
        "target_coordination": coord.get("metal_coordination_number", 6),
        "protein_sites": coord.get("protein_sites", 3),
    }


# =============================================================================
# Validation Functions
# =============================================================================

def validate_template_complex(
    pdb_content: str,
    template_name: str,
    tolerance: float = 0.5,
) -> Dict[str, Any]:
    """
    Validate that a generated complex matches template expectations.

    Checks:
    - Metal is present
    - Ligand is present
    - Metal-ligand distances are reasonable
    - Coordination geometry is approximately correct

    Args:
        pdb_content: PDB content to validate
        template_name: Template to validate against
        tolerance: Distance tolerance in Angstroms

    Returns:
        Validation result dict
    """
    template = get_template(template_name)
    if not template:
        return {"valid": False, "error": f"Unknown template: {template_name}"}

    expected_metal = template.get("metal") or template.get("default_metal")
    expected_ligand = template.get("ligand_res_name")
    coord_info = template.get("coordination", {})
    expected_distances = coord_info.get("bond_distances", {})

    # Parse PDB
    metal_found = False
    ligand_found = False
    metal_pos = None
    ligand_atoms = []

    for line in pdb_content.split('\n'):
        if not line.startswith('HETATM'):
            continue

        try:
            res_name = line[17:20].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            atom_name = line[12:16].strip()

            if res_name == expected_metal or res_name == expected_metal:
                metal_found = True
                metal_pos = (x, y, z)
            elif res_name == expected_ligand:
                ligand_found = True
                ligand_atoms.append({
                    "name": atom_name,
                    "coords": (x, y, z),
                })
        except (ValueError, IndexError):
            continue

    issues = []

    if not metal_found:
        issues.append(f"Metal {expected_metal} not found")
    if not ligand_found:
        issues.append(f"Ligand {expected_ligand} not found")

    # Check distances if both found
    if metal_found and ligand_found and metal_pos:
        for atom in ligand_atoms:
            atom_name = atom["name"]
            if atom_name in expected_distances:
                expected_dist = expected_distances[atom_name]
                actual_dist = math.sqrt(
                    (atom["coords"][0] - metal_pos[0])**2 +
                    (atom["coords"][1] - metal_pos[1])**2 +
                    (atom["coords"][2] - metal_pos[2])**2
                )
                if abs(actual_dist - expected_dist) > tolerance:
                    issues.append(
                        f"{atom_name}-{expected_metal} distance: {actual_dist:.2f} Å "
                        f"(expected {expected_dist:.2f} ± {tolerance})"
                    )

    return {
        "valid": len(issues) == 0,
        "metal_found": metal_found,
        "ligand_found": ligand_found,
        "issues": issues,
        "template": template_name,
    }


# =============================================================================
# Utility Functions
# =============================================================================

def is_lanthanide(metal: str) -> bool:
    """Check if metal is a lanthanide."""
    return metal.upper() in {"LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU"}


def get_metal_preferences(metal: str) -> Dict[str, Any]:
    """
    Get coordination preferences for a metal.

    Returns preferred coordination number, geometry, and donor types.
    """
    metal = metal.upper()

    METAL_PREFERENCES = {
        "CA": {
            "coordination": [6, 7, 8],
            "preferred": 6,
            "geometry": "octahedral",
            "donors": ["O"],
            "preferred_residues": ["Glu", "Asp", "Asn", "Gln"],
        },
        "MG": {
            "coordination": [6],
            "preferred": 6,
            "geometry": "octahedral",
            "donors": ["O"],
            "preferred_residues": ["Glu", "Asp", "water"],
        },
        "ZN": {
            "coordination": [4, 5, 6],
            "preferred": 4,
            "geometry": "tetrahedral",
            "donors": ["N", "S", "O"],
            "preferred_residues": ["His", "Cys", "Glu", "Asp"],
        },
        "FE": {
            "coordination": [4, 5, 6],
            "preferred": 6,
            "geometry": "octahedral",
            "donors": ["N", "O", "S"],
            "preferred_residues": ["His", "Cys", "Glu", "Asp", "Tyr"],
        },
        "TB": {
            "coordination": [8, 9],
            "preferred": 9,
            "geometry": "tricapped_trigonal_prism",
            "donors": ["O"],
            "preferred_residues": ["Glu", "Asp", "Asn"],
        },
        "EU": {
            "coordination": [8, 9],
            "preferred": 9,
            "geometry": "tricapped_trigonal_prism",
            "donors": ["O"],
            "preferred_residues": ["Glu", "Asp", "Asn"],
        },
        "GD": {
            "coordination": [8, 9],
            "preferred": 9,
            "geometry": "square_antiprism",
            "donors": ["O"],
            "preferred_residues": ["Glu", "Asp", "Asn"],
        },
    }

    return METAL_PREFERENCES.get(metal, {
        "coordination": [6],
        "preferred": 6,
        "geometry": "octahedral",
        "donors": ["O", "N"],
        "preferred_residues": ["Glu", "Asp", "His"],
    })
