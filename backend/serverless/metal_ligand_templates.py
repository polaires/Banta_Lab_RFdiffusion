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
#
# Coordination accounting:
#   Metal (Tb³⁺) total capacity: 9
#   Ligand (citrate) provides: 3 (O2, O5, O7 - tridentate based on PDB 3C9H)
#   Protein must provide: 6 (9 - 3 = 6 remaining sites)
#
# Ligand atom roles (based on PDB 3C9H citrate-Ca geometry):
#   Metal-coordinating: O2, O5, O7 (tridentate, 2.36-2.61Å from metal)
#   Protein H-bond acceptors: O1, O3, O4, O6 (non-coordinating carboxylate O)
#
CITRATE_LANTHANIDE_TEMPLATE = {
    "name": "Citrate-Lanthanide",
    "description": "Citrate-coordinated lanthanide (Tb/Eu/Gd) for luminescent biosensors",
    "metals": ["TB", "EU", "GD", "LA", "CE"],  # Compatible lanthanides
    "default_metal": "TB",
    "ligand_smiles": "OC(CC(=O)[O-])(CC(=O)[O-])C(=O)[O-]",  # Citrate trianion
    "ligand_res_name": "CIT",
    "molecular_weight": 189.10,  # Citrate anion
    "coordination": {
        # Metal coordination budget
        "metal_coordination_number": 9,  # Total capacity for Ln³⁺
        "geometry": "tricapped_trigonal_prism",

        # Ligand contribution (FIXED - pre-formed complex)
        # Based on PDB 3C9H citrate-Ca geometry (representative of citrate-metal binding)
        # O2 (carboxylate), O5 (carboxylate), O7 (hydroxyl) coordinate metal
        "ligand_metal_donors": ["O2", "O5", "O7"],  # Atoms that bind metal
        "ligand_donor_count": 3,  # Citrate is TRIDENTATE based on 3C9H
        "ligand_denticity": "tridentate",
        "ligand_metal_bond_distance": 2.48,  # Based on 3C9H: O2=2.48, O5=2.61, O7=2.36Å

        # Protein contribution (NEEDED - what design must achieve)
        "protein_sites_needed": 6,  # Remaining sites: 9 - 3 = 6
        "protein_donor_distance": 2.4,  # Target Ln-O distance for protein (Å)

        # Ligand atoms available for protein interaction (H-bonds, not metal coord)
        # O1, O3, O4, O6 are non-coordinating carboxylate oxygens
        "ligand_hbond_acceptors": ["O1", "O3", "O4", "O6"],  # Non-coordinating O atoms
        "ligand_hbond_distance": 2.8,  # Typical H-bond distance (Å)
    },
    # Suggested distribution for dimer interface
    "dimer_coordination_split": {
        "chain_a_target": 3,  # 3 donors from chain A
        "chain_b_target": 3,  # 3 donors from chain B
        "symmetric_ok": True,  # 3+3 symmetric
    },
    "preferred_protein_donors": ["Glu", "Asp", "Asn", "His", "water"],
    "donor_atom_types": ["OE1", "OE2", "OD1", "OD2", "NE2"],
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

    coord = template.get("coordination", {})
    metal_cn = coord.get("metal_coordination_number", 6)
    ligand_donors = coord.get("ligand_donor_count", 0)

    # Handle both field names: "protein_sites" (PQQ) and "protein_sites_needed" (citrate)
    protein_sites = coord.get("protein_sites", coord.get("protein_sites_needed", metal_cn - ligand_donors))

    return {
        "name": template.get("name"),
        "description": template.get("description"),
        "metal": template.get("metal") or template.get("default_metal"),
        "ligand": template.get("ligand_res_name"),
        "ligand_smiles": template.get("ligand_smiles"),
        "coordination_number": metal_cn,
        "ligand_donors": ligand_donors,
        "protein_sites": protein_sites,
        "geometry": coord.get("geometry", "octahedral"),
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

    # Priority 1: Library templates (curated with matching atom names)
    # IMPORTANT: Library templates have coordination info (ligand_donors, etc.) that
    # match our fallback PDB atom names. PDB-derived templates use CCD atom names
    # which may differ, causing select_buried/hotspots to fail.
    if template_key in METAL_LIGAND_COMPLEX_TEMPLATES:
        template = dict(METAL_LIGAND_COMPLEX_TEMPLATES[template_key])
        template["source"] = "library"
        return template

    # Priority 2: Try PDB database (for templates not in library)
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

    # For citrate-lanthanide complexes, use fallback with proper coordination geometry
    # RDKit doesn't understand metal coordination, so the fallback has correct Ln-O distances
    if "citrate" in template_name.lower():
        print("[MetalLigandTemplates] Using coordination-aware geometry for citrate-Ln")
        return _generate_fallback_pdb(template_name, metal_code, center)

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
        # PQQ-Ca using REAL geometry from PDB 1W6S (methanol dehydrogenase)
        # Resolution: 1.2 Angstrom - high quality experimental structure
        #
        # Reference: Williams et al. (2005) Acta Cryst D61:75-79
        # "Methanol dehydrogenase, a PQQ-containing quinoprotein dehydrogenase"
        #
        # Original metal position in 1W6S chain A: (6.082, 22.564, -14.869)
        #
        # Coordination geometry (tridentate):
        # - O5 (quinone C=O): 2.25 Angstrom from Ca
        # - N6 (pyridine N): 2.31 Angstrom from Ca
        # - O7A (carboxylate O): 2.44 Angstrom from Ca
        #
        # All carbon atoms are >3 Angstrom from metal (RFD3 compatible)
        # Full tricyclic quinoline ring preserved for accurate geometry

        # Translation offset from original metal pos to target center
        # 1W6S metal: (6.082, 22.564, -14.869) -> target: (cx, cy, cz)
        tx = cx - 6.082
        ty = cy - 22.564
        tz = cz - (-14.869)

        # Real PQQ atom positions from PDB 1W6S chain A residue 1596
        # These maintain correct bond lengths and ring geometry
        atoms = {
            # Pyrrole ring (ring A)
            'N1':  (3.250 + tx, 16.518 + ty, -12.689 + tz),   # 7.02 Angstrom from metal
            'C2':  (1.971 + tx, 16.709 + ty, -11.910 + tz),   # 7.74 Angstrom from metal
            'C3':  (1.631 + tx, 18.117 + ty, -11.826 + tz),   # 6.99 Angstrom from metal
            'C3A': (2.703 + tx, 18.759 + ty, -12.535 + tz),   # 5.60 Angstrom from metal
            'C1A': (3.727 + tx, 17.794 + ty, -13.042 + tz),   # 5.62 Angstrom from metal
            # C2 carboxylate (pendant)
            'C2X': (1.280 + tx, 15.500 + ty, -11.390 + tz),   # 9.22 Angstrom from metal
            'O2A': (1.645 + tx, 14.344 + ty, -11.761 + tz),   # 9.84 Angstrom from metal
            'O2B': (0.280 + tx, 15.701 + ty, -10.558 + tz),   # 9.97 Angstrom from metal
            # Quinone ring (ring B) - contains coordinating O5
            'C4':  (2.818 + tx, 20.191 + ty, -12.853 + tz),   # 4.51 Angstrom from metal
            'O4':  (1.839 + tx, 21.051 + ty, -12.585 + tz),   # 5.05 Angstrom from metal (ketone)
            'C5':  (3.900 + tx, 20.662 + ty, -13.461 + tz),   # 3.22 Angstrom from metal
            'O5':  (4.091 + tx, 21.938 + ty, -14.022 + tz),   # 2.25 Angstrom COORDINATING
            'C6A': (4.948 + tx, 19.731 + ty, -13.992 + tz),   # 3.18 Angstrom from metal
            # Pyridine ring (ring C) - contains coordinating N6
            'N6':  (5.880 + tx, 20.261 + ty, -14.866 + tz),   # 2.31 Angstrom COORDINATING
            'C7':  (6.890 + tx, 19.466 + ty, -15.499 + tz),   # 3.26 Angstrom from metal
            'C8':  (6.875 + tx, 18.045 + ty, -15.256 + tz),   # 4.60 Angstrom from metal
            'C9':  (5.881 + tx, 17.331 + ty, -14.460 + tz),   # 5.25 Angstrom from metal
            'C9A': (4.853 + tx, 18.236 + ty, -13.817 + tz),   # 4.62 Angstrom from metal
            # C7 carboxylate - contains coordinating O7A
            'C7X': (7.810 + tx, 20.257 + ty, -16.402 + tz),   # 3.26 Angstrom from metal
            'O7A': (7.678 + tx, 21.485 + ty, -16.371 + tz),   # 2.44 Angstrom COORDINATING
            'O7B': (8.661 + tx, 19.640 + ty, -17.042 + tz),   # 4.46 Angstrom from metal
            # C9 carboxylate (pendant)
            'C9X': (5.926 + tx, 15.807 + ty, -14.467 + tz),   # 6.77 Angstrom from metal
            'O9A': (6.904 + tx, 15.311 + ty, -15.202 + tz),   # 7.31 Angstrom from metal
            'O9B': (5.157 + tx, 15.120 + ty, -13.748 + tz),   # 7.58 Angstrom from metal
        }

        lines = [
            # Ring A - Pyrrole
            f"HETATM    1  N1  PQQ L   1    {atoms['N1'][0]:8.3f}{atoms['N1'][1]:8.3f}{atoms['N1'][2]:8.3f}  1.00  0.00           N",
            f"HETATM    2  C2  PQQ L   1    {atoms['C2'][0]:8.3f}{atoms['C2'][1]:8.3f}{atoms['C2'][2]:8.3f}  1.00  0.00           C",
            f"HETATM    3  C3  PQQ L   1    {atoms['C3'][0]:8.3f}{atoms['C3'][1]:8.3f}{atoms['C3'][2]:8.3f}  1.00  0.00           C",
            f"HETATM    4  C3A PQQ L   1    {atoms['C3A'][0]:8.3f}{atoms['C3A'][1]:8.3f}{atoms['C3A'][2]:8.3f}  1.00  0.00           C",
            f"HETATM    5  C1A PQQ L   1    {atoms['C1A'][0]:8.3f}{atoms['C1A'][1]:8.3f}{atoms['C1A'][2]:8.3f}  1.00  0.00           C",
            # C2 carboxylate
            f"HETATM    6  C2X PQQ L   1    {atoms['C2X'][0]:8.3f}{atoms['C2X'][1]:8.3f}{atoms['C2X'][2]:8.3f}  1.00  0.00           C",
            f"HETATM    7  O2A PQQ L   1    {atoms['O2A'][0]:8.3f}{atoms['O2A'][1]:8.3f}{atoms['O2A'][2]:8.3f}  1.00  0.00           O",
            f"HETATM    8  O2B PQQ L   1    {atoms['O2B'][0]:8.3f}{atoms['O2B'][1]:8.3f}{atoms['O2B'][2]:8.3f}  1.00  0.00           O",
            # Ring B - Quinone
            f"HETATM    9  C4  PQQ L   1    {atoms['C4'][0]:8.3f}{atoms['C4'][1]:8.3f}{atoms['C4'][2]:8.3f}  1.00  0.00           C",
            f"HETATM   10  O4  PQQ L   1    {atoms['O4'][0]:8.3f}{atoms['O4'][1]:8.3f}{atoms['O4'][2]:8.3f}  1.00  0.00           O",
            f"HETATM   11  C5  PQQ L   1    {atoms['C5'][0]:8.3f}{atoms['C5'][1]:8.3f}{atoms['C5'][2]:8.3f}  1.00  0.00           C",
            f"HETATM   12  O5  PQQ L   1    {atoms['O5'][0]:8.3f}{atoms['O5'][1]:8.3f}{atoms['O5'][2]:8.3f}  1.00  0.00           O",
            f"HETATM   13  C6A PQQ L   1    {atoms['C6A'][0]:8.3f}{atoms['C6A'][1]:8.3f}{atoms['C6A'][2]:8.3f}  1.00  0.00           C",
            # Ring C - Pyridine
            f"HETATM   14  N6  PQQ L   1    {atoms['N6'][0]:8.3f}{atoms['N6'][1]:8.3f}{atoms['N6'][2]:8.3f}  1.00  0.00           N",
            f"HETATM   15  C7  PQQ L   1    {atoms['C7'][0]:8.3f}{atoms['C7'][1]:8.3f}{atoms['C7'][2]:8.3f}  1.00  0.00           C",
            f"HETATM   16  C8  PQQ L   1    {atoms['C8'][0]:8.3f}{atoms['C8'][1]:8.3f}{atoms['C8'][2]:8.3f}  1.00  0.00           C",
            f"HETATM   17  C9  PQQ L   1    {atoms['C9'][0]:8.3f}{atoms['C9'][1]:8.3f}{atoms['C9'][2]:8.3f}  1.00  0.00           C",
            f"HETATM   18  C9A PQQ L   1    {atoms['C9A'][0]:8.3f}{atoms['C9A'][1]:8.3f}{atoms['C9A'][2]:8.3f}  1.00  0.00           C",
            # C7 carboxylate
            f"HETATM   19  C7X PQQ L   1    {atoms['C7X'][0]:8.3f}{atoms['C7X'][1]:8.3f}{atoms['C7X'][2]:8.3f}  1.00  0.00           C",
            f"HETATM   20  O7A PQQ L   1    {atoms['O7A'][0]:8.3f}{atoms['O7A'][1]:8.3f}{atoms['O7A'][2]:8.3f}  1.00  0.00           O",
            f"HETATM   21  O7B PQQ L   1    {atoms['O7B'][0]:8.3f}{atoms['O7B'][1]:8.3f}{atoms['O7B'][2]:8.3f}  1.00  0.00           O",
            # C9 carboxylate
            f"HETATM   22  C9X PQQ L   1    {atoms['C9X'][0]:8.3f}{atoms['C9X'][1]:8.3f}{atoms['C9X'][2]:8.3f}  1.00  0.00           C",
            f"HETATM   23  O9A PQQ L   1    {atoms['O9A'][0]:8.3f}{atoms['O9A'][1]:8.3f}{atoms['O9A'][2]:8.3f}  1.00  0.00           O",
            f"HETATM   24  O9B PQQ L   1    {atoms['O9B'][0]:8.3f}{atoms['O9B'][1]:8.3f}{atoms['O9B'][2]:8.3f}  1.00  0.00           O",
            # Metal at coordination center
            f"HETATM   25 CA   CA  M   2    {cx:8.3f}{cy:8.3f}{cz:8.3f}  1.00  0.00          CA",
            "END",
        ]
    elif "citrate" in template_name.lower():
        # Citrate-Ln using REAL geometry from PDB 3C9H (citrate-Ca structure)
        # This ensures proper Y-shaped citrate with ALL CARBONS >3Å from metal
        #
        # Real citrate structure from 3C9H:
        # - C3 is central quaternary carbon with hydroxyl O7
        # - C1-O1/O2 is terminal carboxylate (via C2 methylene)
        # - C5-O3/O4 is another carboxylate arm (via C4 methylene)
        # - C6-O5/O6 is third carboxylate arm
        #
        # Coordination: Metal binds O2, O5, O7 (tridentate)
        # Key: All carbons must be >3Å from metal!
        #
        # Reference coordinates from 3C9H, translated so metal is at (cx, cy, cz)
        # Original metal position in 3C9H: (48.044, 3.540, 11.340)

        # Translation offset from original metal pos to target center
        # 3C9H metal: (48.044, 3.540, 11.340) -> target: (cx, cy, cz)
        tx = cx - 48.044
        ty = cy - 3.540
        tz = cz - 11.340

        # Real citrate atom positions from PDB 3C9H chain A, translated
        # These maintain correct bond lengths (~1.5Å C-C, ~1.25Å C=O)
        # and keep all carbons >3Å from metal
        atoms = {
            # Terminal carboxylate 1
            'C1': (50.774 + tx, 5.151 + ty, 11.017 + tz),   # 3.19Å from metal
            'O1': (51.750 + tx, 4.983 + ty, 11.801 + tz),   # 4.00Å (non-coord)
            'O2': (50.263 + tx, 4.175 + ty, 10.429 + tz),   # 2.48Å COORDINATING
            # Methylene connecting to central C
            'C2': (50.238 + tx, 6.531 + ty, 10.743 + tz),   # 3.76Å from metal
            # Central quaternary carbon
            'C3': (48.708 + tx, 6.584 + ty, 10.576 + tz),   # 3.21Å from metal
            'O7': (48.063 + tx, 5.873 + ty, 11.676 + tz),   # 2.36Å COORDINATING (hydroxyl)
            # Carboxylate arm 2 (via C4)
            'C4': (48.209 + tx, 8.025 + ty, 10.587 + tz),   # 4.55Å from metal
            'C5': (48.428 + tx, 8.713 + ty, 11.897 + tz),   # 5.22Å from metal
            'O3': (48.912 + tx, 9.889 + ty, 11.949 + tz),   # 6.44Å (non-coord)
            'O4': (48.055 + tx, 8.096 + ty, 12.906 + tz),   # 4.82Å (non-coord)
            # Carboxylate arm 3
            'C6': (48.248 + tx, 5.943 + ty, 9.256 + tz),    # 3.19Å from metal
            'O5': (47.458 + tx, 4.961 + ty, 9.225 + tz),    # 2.61Å COORDINATING
            'O6': (48.661 + tx, 6.370 + ty, 8.171 + tz),    # 4.29Å (non-coord)
        }

        lines = [
            f"HETATM    1  C1  CIT L   1    {atoms['C1'][0]:8.3f}{atoms['C1'][1]:8.3f}{atoms['C1'][2]:8.3f}  1.00  0.00           C",
            f"HETATM    2  O1  CIT L   1    {atoms['O1'][0]:8.3f}{atoms['O1'][1]:8.3f}{atoms['O1'][2]:8.3f}  1.00  0.00           O",
            f"HETATM    3  O2  CIT L   1    {atoms['O2'][0]:8.3f}{atoms['O2'][1]:8.3f}{atoms['O2'][2]:8.3f}  1.00  0.00           O",
            f"HETATM    4  C2  CIT L   1    {atoms['C2'][0]:8.3f}{atoms['C2'][1]:8.3f}{atoms['C2'][2]:8.3f}  1.00  0.00           C",
            f"HETATM    5  C3  CIT L   1    {atoms['C3'][0]:8.3f}{atoms['C3'][1]:8.3f}{atoms['C3'][2]:8.3f}  1.00  0.00           C",
            f"HETATM    6  O7  CIT L   1    {atoms['O7'][0]:8.3f}{atoms['O7'][1]:8.3f}{atoms['O7'][2]:8.3f}  1.00  0.00           O",
            f"HETATM    7  C4  CIT L   1    {atoms['C4'][0]:8.3f}{atoms['C4'][1]:8.3f}{atoms['C4'][2]:8.3f}  1.00  0.00           C",
            f"HETATM    8  C5  CIT L   1    {atoms['C5'][0]:8.3f}{atoms['C5'][1]:8.3f}{atoms['C5'][2]:8.3f}  1.00  0.00           C",
            f"HETATM    9  O3  CIT L   1    {atoms['O3'][0]:8.3f}{atoms['O3'][1]:8.3f}{atoms['O3'][2]:8.3f}  1.00  0.00           O",
            f"HETATM   10  O4  CIT L   1    {atoms['O4'][0]:8.3f}{atoms['O4'][1]:8.3f}{atoms['O4'][2]:8.3f}  1.00  0.00           O",
            f"HETATM   11  C6  CIT L   1    {atoms['C6'][0]:8.3f}{atoms['C6'][1]:8.3f}{atoms['C6'][2]:8.3f}  1.00  0.00           C",
            f"HETATM   12  O5  CIT L   1    {atoms['O5'][0]:8.3f}{atoms['O5'][1]:8.3f}{atoms['O5'][2]:8.3f}  1.00  0.00           O",
            f"HETATM   13  O6  CIT L   1    {atoms['O6'][0]:8.3f}{atoms['O6'][1]:8.3f}{atoms['O6'][2]:8.3f}  1.00  0.00           O",
            # Metal at coordination center
            f"HETATM   14 {metal_code:4s} {metal_code:3s} M   2    {cx:8.3f}{cy:8.3f}{cz:8.3f}  1.00  0.00          {metal_code[:2]:>2s}",
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
    Get coordination information for RFD3 configuration and validation.

    Returns dict with:
    - Coordination accounting (metal capacity, ligand donors, protein needed)
    - RFD3 config (hotspot_atoms, hbond_acceptors, exposed_atoms)
    - Validation parameters (ligand H-bond acceptors)
    """
    template = get_template(template_name)
    if not template:
        return {}

    coord = template.get("coordination", {})
    ligand_name = template.get("ligand_res_name", "LIG")
    metal = template.get("metal") or template.get("default_metal", "")

    # Ligand atoms that coordinate metal (fixed in pre-formed complex)
    ligand_metal_donors = coord.get("ligand_metal_donors", coord.get("ligand_donors", []))

    # Ligand atoms available for protein H-bonds (not coordinating metal)
    ligand_hbond_acceptors = coord.get("ligand_hbond_acceptors", [])

    # Coordination accounting
    metal_coord_number = coord.get("metal_coordination_number", 6)
    ligand_donor_count = coord.get("ligand_donor_count", len(ligand_metal_donors))
    protein_sites_needed = coord.get("protein_sites_needed", metal_coord_number - ligand_donor_count)

    # Dimer coordination split
    dimer_split = template.get("dimer_coordination_split", {})

    # Atoms to expose (not coordinating metal, available for protein contact)
    exposed_atoms = []

    return {
        # Basic info
        "ligand_name": ligand_name,
        "metal_code": metal,

        # Coordination accounting (for validation)
        "coordination": {
            "metal_coordination_number": metal_coord_number,
            "ligand_metal_donors": ligand_metal_donors,
            "ligand_donor_count": ligand_donor_count,
            "protein_sites_needed": protein_sites_needed,
            "ligand_hbond_acceptors": ligand_hbond_acceptors,
            "ligand_metal_bond_distance": coord.get("ligand_metal_bond_distance", 2.35),
            "geometry": coord.get("geometry", "octahedral"),
        },

        # Dimer interface targets
        "dimer_split": {
            "chain_a_target": dimer_split.get("chain_a_target", protein_sites_needed // 2),
            "chain_b_target": dimer_split.get("chain_b_target", protein_sites_needed - protein_sites_needed // 2),
        },

        # RFD3 configuration
        "hotspot_atoms": {metal: "ALL"},  # Ensure protein contacts metal
        "hbond_acceptors": {ligand_name: ",".join(ligand_hbond_acceptors)} if ligand_hbond_acceptors else {},
        "exposed_atoms": {ligand_name: ",".join(exposed_atoms)} if exposed_atoms else {},
        "preferred_donors": template.get("preferred_protein_donors", []),

        # Legacy compatibility
        "target_coordination": metal_coord_number,
        "protein_sites": protein_sites_needed,
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
