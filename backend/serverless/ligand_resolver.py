"""
Ligand Resolver Module

Resolves ligand structures from multiple sources and analyzes chemistry
for RFD3 configuration.

Resolution priority:
1. Pre-built templates (metal_ligand_templates.py)
2. PDB experimental structures (structure_discovery.py)
3. PubChem SMILES (database_adapters/pubchem_adapter.py)

Also analyzes ligand chemistry to recommend:
- Hotspot atoms for RFD3 select_hotspots
- H-bond donors/acceptors for conditioning
- Burial recommendations
- CFG scale suggestions
"""

import logging
from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional, Tuple

# Import template library
from metal_ligand_templates import (
    get_template,
    get_template_info,
    get_template_with_fallback,
    get_template_coordination_info,
    generate_complex_pdb,
    METAL_LIGAND_COMPLEX_TEMPLATES,
    KNOWN_COMPLEX_PDBS,
)

# Import metal chemistry
from metal_chemistry import (
    METAL_DATABASE,
    is_lanthanide,
    get_coordination_number_range,
)

logger = logging.getLogger(__name__)


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class LigandChemistry:
    """
    Chemical analysis of a ligand for RFD3 configuration.
    """
    # Atom classification
    donor_atoms: List[str] = field(default_factory=list)       # Metal-coordinating atoms
    hbond_donors: List[str] = field(default_factory=list)      # H-bond donor atoms
    hbond_acceptors: List[str] = field(default_factory=list)   # H-bond acceptor atoms
    hydrophobic_atoms: List[str] = field(default_factory=list) # Hydrophobic atoms

    # NEW: Additional fields for decision engine
    coordinating_atoms: List[str] = field(default_factory=list)  # O, N, S atoms that can coordinate metals
    heavy_atom_count: int = 0                                     # Total heavy atoms (non-H)

    # Counts
    total_atoms: int = 0
    polar_atoms: int = 0
    nonpolar_atoms: int = 0

    # Character classification
    ligand_character: str = "mixed"  # polar | hydrophobic | mixed

    # RFD3 recommendations
    recommended_burial: str = "partial"  # buried | partial | exposed
    hotspot_atoms: List[str] = field(default_factory=list)
    cfg_scale_suggestion: float = 2.5

    # Warnings
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return asdict(self)


@dataclass
class ResolvedLigand:
    """
    Fully resolved ligand with structure and chemistry.
    """
    # Identity
    name: str = ""                    # Normalized ligand name
    residue_code: str = "LIG"         # 3-letter PDB residue code
    smiles: Optional[str] = None      # SMILES string if available

    # Structure
    pdb_content: Optional[str] = None # PDB content for the ligand+metal complex
    source: str = "unknown"           # template | pdb | pubchem | calculated

    # Metal association
    metal_type: Optional[str] = None  # Associated metal
    metal_coords: Optional[Tuple[float, float, float]] = None

    # Chemistry analysis
    chemistry: Optional[LigandChemistry] = None

    # Coordination info (from template)
    coordination: Dict[str, Any] = field(default_factory=dict)

    # Template info
    template_name: Optional[str] = None
    template_info: Dict[str, Any] = field(default_factory=dict)

    # Status
    resolved: bool = False
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        result = asdict(self)
        if self.chemistry:
            result["chemistry"] = self.chemistry.to_dict()
        return result


# =============================================================================
# Ligand Name to Template Mapping
# =============================================================================

LIGAND_TO_TEMPLATE: Dict[str, Dict[str, str]] = {
    # Citrate with various metals
    "citrate": {
        "TB": "citrate_tb",
        "EU": "citrate_eu",
        "GD": "citrate_gd",
        "LA": "citrate_la",
        "default": "citrate_tb",  # Default lanthanide
    },
    # PQQ
    "pqq": {
        "CA": "pqq_ca",
        "default": "pqq_ca",
    },
}


# Isomeric SMILES for ligands with known cis/trans variants
ISOMERIC_SMILES: Dict[str, Dict[str, str]] = {
    "azobenzene": {
        "trans": "c1ccc(/N=N/c2ccccc2)cc1",
        "cis": r"c1ccc(/N=N\c2ccccc2)cc1",
        "default": "c1ccc(/N=N/c2ccccc2)cc1",
    },
    "stilbene": {
        "trans": "C(/c1ccccc1)=C/c1ccccc1",
        "cis": r"C(/c1ccccc1)=C\c1ccccc1",
        "default": "C(/c1ccccc1)=C/c1ccccc1",
    },
}


def get_template_name_for_ligand(ligand: str, metal: Optional[str] = None) -> Optional[str]:
    """
    Get template name for a ligand-metal combination.

    Args:
        ligand: Ligand name (normalized)
        metal: Metal symbol (optional)

    Returns:
        Template name or None if not found
    """
    ligand_lower = ligand.lower()

    if ligand_lower not in LIGAND_TO_TEMPLATE:
        return None

    metal_map = LIGAND_TO_TEMPLATE[ligand_lower]

    if metal and metal.upper() in metal_map:
        return metal_map[metal.upper()]

    return metal_map.get("default")


# =============================================================================
# Chemistry Analyzer
# =============================================================================

# Atom type classification
POLAR_ATOMS = {"O", "N", "S"}
HBOND_DONOR_ATOMS = {"N", "O"}  # When bonded to H
HBOND_ACCEPTOR_ATOMS = {"O", "N", "S"}
HYDROPHOBIC_ATOMS = {"C"}


def analyze_ligand_chemistry(
    pdb_content: str,
    ligand_res_name: str = "LIG",
    template_info: Optional[Dict[str, Any]] = None,
) -> LigandChemistry:
    """
    Analyze ligand chemistry from PDB content.

    Args:
        pdb_content: PDB content containing the ligand
        ligand_res_name: Residue name of the ligand
        template_info: Optional template info with coordination data

    Returns:
        LigandChemistry with analysis results
    """
    chem = LigandChemistry()

    # Parse ligand atoms from PDB
    ligand_atoms = []
    metal_atoms = []

    for line in pdb_content.strip().split('\n'):
        if not line.startswith('HETATM'):
            continue

        try:
            res_name = line[17:20].strip()
            atom_name = line[12:16].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            element = line[76:78].strip() if len(line) >= 78 else atom_name[0]

            if res_name == ligand_res_name:
                ligand_atoms.append({
                    "name": atom_name,
                    "element": element,
                    "coords": (x, y, z),
                })
            elif res_name in METAL_DATABASE or element in METAL_DATABASE:
                metal_atoms.append({
                    "name": atom_name,
                    "element": element,
                    "coords": (x, y, z),
                })
        except (ValueError, IndexError):
            continue

    chem.total_atoms = len(ligand_atoms)

    if not ligand_atoms:
        chem.warnings.append(f"No atoms found for ligand {ligand_res_name}")
        return chem

    # Classify atoms by element
    for atom in ligand_atoms:
        element = atom["element"].upper()
        atom_name = atom["name"]

        if element in POLAR_ATOMS:
            chem.polar_atoms += 1
            if element in HBOND_ACCEPTOR_ATOMS:
                chem.hbond_acceptors.append(atom_name)
            # NEW: Track coordinating atoms (O, N, S can coordinate metals)
            chem.coordinating_atoms.append(atom_name)
        else:
            chem.nonpolar_atoms += 1
            if element == "C":
                chem.hydrophobic_atoms.append(atom_name)

    # NEW: Calculate heavy atom count for decision engine
    chem.heavy_atom_count = chem.polar_atoms + chem.nonpolar_atoms

    # Use template coordination info if available
    if template_info:
        coord = template_info.get("coordination", {})

        # Metal-coordinating atoms (from template)
        chem.donor_atoms = coord.get("ligand_metal_donors", [])

        # H-bond acceptors for protein contacts (from template)
        template_hbond = coord.get("ligand_hbond_acceptors", [])
        if template_hbond:
            chem.hbond_acceptors = template_hbond

        # H-bond donors (atoms that can donate H to protein)
        chem.hbond_donors = coord.get("ligand_hbond_donors", [])

    else:
        # Infer from structure if no template
        # Calculate distances to metal
        if metal_atoms:
            metal_coord = metal_atoms[0]["coords"]
            for atom in ligand_atoms:
                if atom["element"].upper() in {"O", "N", "S"}:
                    dist = _distance(atom["coords"], metal_coord)
                    if dist < 3.0:  # Coordination distance
                        chem.donor_atoms.append(atom["name"])

    # Classify ligand character
    if chem.total_atoms > 0:
        polar_ratio = chem.polar_atoms / chem.total_atoms
        if polar_ratio > 0.5:
            chem.ligand_character = "polar"
        elif polar_ratio < 0.2:
            chem.ligand_character = "hydrophobic"
        else:
            chem.ligand_character = "mixed"

    # Recommend burial based on character
    if chem.ligand_character == "polar":
        chem.recommended_burial = "partial"  # Some exposure for solvation
        chem.cfg_scale_suggestion = 2.5
    elif chem.ligand_character == "hydrophobic":
        chem.recommended_burial = "buried"  # Full burial
        chem.cfg_scale_suggestion = 2.0
    else:
        chem.recommended_burial = "partial"
        chem.cfg_scale_suggestion = 2.5

    # Recommend hotspot atoms
    # For polar ligands: use H-bond acceptors
    # For hydrophobic: use carbon atoms
    if chem.ligand_character == "hydrophobic":
        chem.hotspot_atoms = chem.hydrophobic_atoms[:10]  # Limit to 10
    else:
        # Prefer atoms that are NOT metal donors (those are fixed)
        chem.hotspot_atoms = [
            a for a in chem.hbond_acceptors
            if a not in chem.donor_atoms
        ][:10]

    return chem


def _distance(c1: Tuple[float, float, float], c2: Tuple[float, float, float]) -> float:
    """Calculate Euclidean distance between two 3D points."""
    return ((c1[0]-c2[0])**2 + (c1[1]-c2[1])**2 + (c1[2]-c2[2])**2) ** 0.5


# =============================================================================
# Main Resolver Class
# =============================================================================

class LigandResolver:
    """
    Resolves ligand structures and analyzes chemistry.

    Priority:
    1. Pre-built templates (fast, curated)
    2. PDB experimental structures (accurate)
    3. PubChem SMILES (comprehensive)
    4. Calculated fallback (last resort)
    """

    def __init__(
        self,
        use_pdb: bool = True,
        use_pubchem: bool = True,
        default_center: Tuple[float, float, float] = (50.0, 50.0, 50.0),
    ):
        """
        Initialize resolver.

        Args:
            use_pdb: Enable PDB structure lookup
            use_pubchem: Enable PubChem SMILES lookup
            default_center: Default coordinates for generated structures
        """
        self.use_pdb = use_pdb
        self.use_pubchem = use_pubchem
        self.default_center = default_center

    def resolve(
        self,
        ligand_name: str,
        metal_type: Optional[str] = None,
        center: Optional[Tuple[float, float, float]] = None,
        isomer_spec: Optional[str] = None,
    ) -> ResolvedLigand:
        """
        Resolve a ligand by name, optionally with associated metal.

        Args:
            ligand_name: Ligand name (e.g., "citrate", "pqq")
            metal_type: Metal symbol (e.g., "TB", "CA")
            center: Coordinates for structure placement
            isomer_spec: Isomer specification ("cis", "trans") for ligands with known variants

        Returns:
            ResolvedLigand with structure and chemistry
        """
        if center is None:
            center = self.default_center

        result = ResolvedLigand(
            name=ligand_name.lower(),
            metal_type=metal_type.upper() if metal_type else None,
        )

        # Store isomer spec for use in PubChem resolution
        self._current_isomer_spec = isomer_spec

        # Try resolution sources in priority order

        # 1. Pre-built templates (PREFERRED - curated with matching atom names)
        if self._try_template_resolution(result, center):
            result.resolved = True
            result.source = "template"
            return result

        # 2. PDB experimental structures
        if self.use_pdb and self._try_pdb_resolution(result, center):
            result.resolved = True
            result.source = "pdb"
            return result

        # 3. PubChem SMILES
        if self.use_pubchem and self._try_pubchem_resolution(result, center):
            result.resolved = True
            result.source = "pubchem"
            return result

        # 4. Calculated fallback
        if self._try_calculated_fallback(result, center):
            result.resolved = True
            result.source = "calculated"
            result.warnings.append("Using calculated fallback - verify geometry")
            return result

        # Failed to resolve
        result.warnings.append(f"Could not resolve ligand: {ligand_name}")
        return result

    def _try_template_resolution(
        self,
        result: ResolvedLigand,
        center: Tuple[float, float, float],
    ) -> bool:
        """Try to resolve from pre-built templates."""
        template_name = get_template_name_for_ligand(
            result.name,
            result.metal_type
        )

        if not template_name:
            return False

        # Get template info
        template_info = get_template_info(template_name)
        if not template_info:
            return False

        # Generate PDB content
        pdb_content = generate_complex_pdb(
            template_name,
            metal=result.metal_type,
            center=center,
        )

        if not pdb_content:
            return False

        # Get coordination info
        coord_info = get_template_coordination_info(template_name)

        # Store results
        result.template_name = template_name
        result.template_info = template_info
        result.coordination = coord_info.get("coordination", {})
        result.pdb_content = pdb_content
        result.residue_code = template_info.get("ligand", "LIG")
        result.smiles = template_info.get("ligand_smiles")
        result.metal_coords = center  # Metal placed at center

        # Analyze chemistry
        result.chemistry = analyze_ligand_chemistry(
            pdb_content,
            result.residue_code,
            coord_info,
        )

        logger.info(f"Resolved {result.name} from template: {template_name}")
        return True

    def _try_pdb_resolution(
        self,
        result: ResolvedLigand,
        center: Tuple[float, float, float],
    ) -> bool:
        """Try to resolve from PDB structures."""
        # Check if we have known PDB entries
        ligand_lower = result.name.lower()
        metal = result.metal_type

        # Build lookup key
        lookup_key = None
        for key, info in KNOWN_COMPLEX_PDBS.items():
            if ligand_lower in key:
                if metal and info.get("metal") == metal:
                    lookup_key = key
                    break
                elif not metal:
                    lookup_key = key
                    break

        if not lookup_key:
            return False

        try:
            template = get_template_with_fallback(
                lookup_key,
                metal=metal,
                ligand=result.name,
                fallback_enabled=False,
            )

            if template and template.get("source") == "pdb":
                # Use PDB-derived template
                result.template_info = template
                # Would need to generate PDB content from PDB...
                # For now, fall through to template
                return False

        except Exception as e:
            logger.debug(f"PDB resolution failed: {e}")

        return False

    def _try_pubchem_resolution(
        self,
        result: ResolvedLigand,
        center: Tuple[float, float, float],
    ) -> bool:
        """Try to resolve from PubChem SMILES (with isomeric SMILES support)."""
        ligand_lower = result.name.lower()
        isomer_spec = getattr(self, '_current_isomer_spec', None)

        # 1. Check ISOMERIC_SMILES table first (for known isomers)
        if ligand_lower in ISOMERIC_SMILES:
            isomer_table = ISOMERIC_SMILES[ligand_lower]
            if isomer_spec and isomer_spec in isomer_table:
                smiles = isomer_table[isomer_spec]
                logger.info(f"Using {isomer_spec}-isomer SMILES for {ligand_lower}")
            else:
                smiles = isomer_table.get("default", "")
                if isomer_spec:
                    result.warnings.append(
                        f"Isomer '{isomer_spec}' not found for {ligand_lower}, using default"
                    )
        else:
            # 2. Try PubChem API (already returns IsomericSMILES)
            smiles = None
            try:
                from database_adapters.pubchem_adapter import PubChemAdapter
                adapter = PubChemAdapter()
                smiles = adapter.get_smiles(result.name)
            except Exception as e:
                logger.debug(f"PubChem resolution failed for {result.name}: {e}")
                return False

        if not smiles:
            return False

        result.smiles = smiles

        # 3. Generate 3D conformer from SMILES
        try:
            from conformer_utils import generate_conformer
            pdb_content = generate_conformer(
                smiles,
                name=result.residue_code or "LIG",
                center=center,
            )
            if pdb_content:
                result.pdb_content = pdb_content
        except Exception as e:
            logger.warning(f"Conformer generation failed for {result.name}: {e}")
            # Continue without 3D structure — SMILES alone is still useful

        # 4. Analyze chemistry if we have PDB content
        if result.pdb_content:
            result.chemistry = analyze_ligand_chemistry(
                result.pdb_content,
                result.residue_code or "LIG",
            )

        # 5. Set recommended fixing strategy
        # PubChem ligands have no metal coordination info, so default to diffuse_all
        # (ISOMERIC_SMILES ligands are typically organic — same logic)
        result.coordination["recommended_fixing_strategy"] = "diffuse_all"

        logger.info(f"Resolved {result.name} via PubChem/isomeric table (SMILES: {smiles[:40]}...)")
        return True

    def _try_calculated_fallback(
        self,
        result: ResolvedLigand,
        center: Tuple[float, float, float],
    ) -> bool:
        """Generate a calculated fallback when other methods fail."""
        if not result.metal_type:
            return False

        try:
            template = get_template_with_fallback(
                f"calculated_{result.metal_type}",
                metal=result.metal_type,
                ligand=result.name,
                fallback_enabled=True,
            )

            if template and template.get("source") == "calculated":
                result.template_info = template
                result.coordination = {
                    "metal_coordination_number": template.get("coordination_number", 6),
                    "geometry": template.get("geometry", "octahedral"),
                }

                # Generate minimal PDB
                cx, cy, cz = center
                metal = result.metal_type
                result.pdb_content = (
                    f"HETATM    1  C1  LIG L   1    {cx:8.3f}{cy:8.3f}{cz:8.3f}  1.00  0.00           C\n"
                    f"HETATM    2 {metal:4s} {metal:3s} X   1    {cx:8.3f}{cy:8.3f}{cz+2.4:8.3f}  1.00  0.00          {metal[:2]:>2s}\n"
                    "END"
                )
                result.metal_coords = (cx, cy, cz + 2.4)

                # Minimal chemistry
                result.chemistry = LigandChemistry(
                    total_atoms=1,
                    nonpolar_atoms=1,
                    ligand_character="unknown",
                    cfg_scale_suggestion=2.5,
                    warnings=["Calculated fallback - no real ligand structure"],
                )

                return True

        except Exception as e:
            logger.debug(f"Calculated fallback failed: {e}")

        return False

    def analyze_chemistry(
        self,
        pdb_content: str,
        ligand_res_name: str = "LIG",
    ) -> LigandChemistry:
        """
        Analyze ligand chemistry from PDB content.

        Args:
            pdb_content: PDB content containing the ligand
            ligand_res_name: Residue name of the ligand

        Returns:
            LigandChemistry analysis
        """
        return analyze_ligand_chemistry(pdb_content, ligand_res_name)


# =============================================================================
# Convenience Functions
# =============================================================================

def resolve_ligand(
    ligand_name: str,
    metal_type: Optional[str] = None,
    isomer_spec: Optional[str] = None,
) -> ResolvedLigand:
    """
    Convenience function to resolve a ligand.

    Args:
        ligand_name: Ligand name
        metal_type: Metal symbol
        isomer_spec: Isomer specification ("cis", "trans")

    Returns:
        ResolvedLigand
    """
    resolver = LigandResolver()
    return resolver.resolve(ligand_name, metal_type, isomer_spec=isomer_spec)


def get_recommended_rfd3_params(
    resolved: ResolvedLigand,
) -> Dict[str, Any]:
    """
    Get recommended RFD3 parameters based on resolved ligand.

    Args:
        resolved: Resolved ligand with chemistry analysis

    Returns:
        Dict with RFD3 parameter recommendations
    """
    params = {
        "select_fixed_atoms": {},
        "select_hotspots": {},
        "select_buried": {},
        "select_hbond_acceptor": {},
        "cfg_scale": 2.5,
    }

    if not resolved.resolved:
        return params

    # Metal chain (typically X1)
    if resolved.metal_type:
        params["select_fixed_atoms"]["X1"] = "all"

    # Ligand chain (typically L1)
    params["select_fixed_atoms"]["L1"] = "all"

    # Hotspots based on chemistry
    if resolved.chemistry:
        chem = resolved.chemistry

        # Use hotspot atoms from analysis
        if chem.hotspot_atoms:
            params["select_hotspots"]["L1"] = ",".join(chem.hotspot_atoms)

        # Burial recommendation
        if chem.recommended_burial == "buried":
            params["select_buried"]["L1"] = "all"
        elif chem.recommended_burial == "partial":
            # Bury hydrophobic atoms
            if chem.hydrophobic_atoms:
                params["select_buried"]["L1"] = ",".join(chem.hydrophobic_atoms[:5])

        # H-bond acceptors for protein contact
        if chem.hbond_acceptors:
            # Exclude metal-coordinating atoms
            available = [a for a in chem.hbond_acceptors if a not in chem.donor_atoms]
            if available:
                params["select_hbond_acceptor"]["L1"] = ",".join(available[:4])

        # CFG scale from analysis
        params["cfg_scale"] = chem.cfg_scale_suggestion

    # Use template coordination info if available
    if resolved.coordination:
        coord = resolved.coordination

        # Override with template-specific recommendations
        if "ligand_hbond_acceptors" in coord:
            acceptors = coord["ligand_hbond_acceptors"]
            if acceptors:
                params["select_hbond_acceptor"]["L1"] = ",".join(acceptors)

    return params


# =============================================================================
# Test Function
# =============================================================================

def test_resolver():
    """Test the ligand resolver."""
    test_cases = [
        ("citrate", "TB"),
        ("citrate", "EU"),
        ("pqq", "CA"),
        ("unknown_ligand", "ZN"),
    ]

    resolver = LigandResolver()

    for ligand, metal in test_cases:
        print(f"\nResolving: {ligand} with {metal}")
        result = resolver.resolve(ligand, metal)

        print(f"  Resolved: {result.resolved}")
        print(f"  Source: {result.source}")
        print(f"  Template: {result.template_name}")

        if result.chemistry:
            print(f"  Character: {result.chemistry.ligand_character}")
            print(f"  Donor atoms: {result.chemistry.donor_atoms}")
            print(f"  Hotspots: {result.chemistry.hotspot_atoms}")
            print(f"  CFG scale: {result.chemistry.cfg_scale_suggestion}")

        if result.warnings:
            print(f"  Warnings: {result.warnings}")

        # Get RFD3 params
        params = get_recommended_rfd3_params(result)
        print(f"  RFD3 params: {params}")


if __name__ == "__main__":
    test_resolver()
