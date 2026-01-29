"""
Scaffolding Workflow for AI Design Pipeline.

Wraps existing handler.py scaffolding functions for use with AI-driven design.
Enables queries like "scaffold the PQQ-Ca pocket of 4CVB".

Supports generalizable theozyme scaffolding:
- Ligand code normalization (human names → PDB 3-letter codes)
- Metal substitution (extract MG site, rewrite as TB for lanthanide design)
- HSAB-based compatibility checking for metal substitution

This module provides:
- ScaffoldingWorkflow: Main workflow class for motif scaffolding
- ScaffoldResult: Result dataclass with all scaffolding info
- Integration with metal_site_fetcher for PDB fetching and active site extraction
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any, Tuple, Set
import logging
import asyncio
import re

# Import existing infrastructure
from metal_site_fetcher import (
    find_metal_ligand_active_site,
    _fetch_pdb_content,
    extract_metal_coordination,
    REFERENCE_STRUCTURES,
)

# Import metal chemistry for substitution compatibility
from metal_chemistry import (
    METAL_DATABASE,
    get_hsab_class,
    get_coordination_number_range,
)

logger = logging.getLogger(__name__)


# =============================================================================
# LIGAND NAME → PDB 3-LETTER CODE MAPPING
# =============================================================================

LIGAND_NAME_TO_PDB_CODE: Dict[str, str] = {
    # Citrate
    "citrate": "CIT",
    "citric acid": "CIT",
    "citric": "CIT",
    # PQQ
    "pqq": "PQQ",
    "pyrroloquinoline quinone": "PQQ",
    "pyrroloquinoline": "PQQ",
    # Heme / Porphyrin
    "heme": "HEM",
    "hem": "HEM",
    "protoporphyrin": "HEM",
    # Nucleotides
    "atp": "ATP",
    "adp": "ADP",
    "amp": "AMP",
    "gtp": "GTP",
    "gdp": "GDP",
    # Cofactors
    "nad": "NAD",
    "nad+": "NAD",
    "nadh": "NAI",
    "nadp": "NAP",
    "nadp+": "NAP",
    "nadph": "NDP",
    "fad": "FAD",
    "fmn": "FMN",
    "plp": "PLP",
    "pyridoxal phosphate": "PLP",
    # Common small molecules
    "acetate": "ACT",
    "phosphate": "PO4",
    "sulfate": "SO4",
    "carbonate": "CO3",
    "oxalate": "OXL",
    "malonate": "MLN",
    "succinate": "SIN",
    "tartrate": "TAR",
    "edta": "EDO",
    "glycerol": "GOL",
}


def normalize_ligand_code(ligand_code: str) -> str:
    """Normalize a human-readable ligand name to its PDB 3-letter code.

    If the input is already a valid 3-letter code (uppercase, 1-3 chars),
    it is returned as-is. Otherwise, looks up the mapping table.

    Args:
        ligand_code: Human name ("citrate") or PDB code ("CIT")

    Returns:
        PDB 3-letter code (e.g., "CIT")
    """
    if not ligand_code:
        return ligand_code

    code_upper = ligand_code.strip().upper()

    # If it's already a short uppercase code (1-3 chars), return as-is
    # This handles cases where the caller already passes "CIT", "PQQ", etc.
    if len(code_upper) <= 3 and code_upper.isalnum():
        return code_upper

    # Look up the mapping table (case-insensitive)
    code_lower = ligand_code.strip().lower()
    if code_lower in LIGAND_NAME_TO_PDB_CODE:
        return LIGAND_NAME_TO_PDB_CODE[code_lower]

    # Fallback: return uppercase (might be a 3-letter code we don't know)
    return code_upper


# =============================================================================
# METAL SUBSTITUTION COMPATIBILITY
# =============================================================================

# Known metals that appear in PDB HETATM records
KNOWN_PDB_METALS = {
    "ZN", "FE", "MG", "CA", "MN", "CU", "CO", "NI",
    "TB", "EU", "GD", "LA", "CE", "SM", "YB", "DY",
    "CD", "HG", "PB", "MO", "W", "V", "CR", "PT", "RU",
}


def scan_pdb_hetatm(pdb_content: str) -> Tuple[Set[str], Set[str]]:
    """Scan PDB content and return sets of metals and ligands found.

    Args:
        pdb_content: Raw PDB file content

    Returns:
        (metals_found, ligands_found) — sets of 3-letter residue names
    """
    metals = set()
    ligands = set()

    for line in pdb_content.split('\n'):
        if line.startswith('HETATM'):
            res_name = line[17:20].strip()
            if res_name in KNOWN_PDB_METALS:
                metals.add(res_name)
            elif len(res_name) >= 1 and res_name not in {"HOH", "WAT", "DOD"}:
                ligands.add(res_name)

    return metals, ligands


def is_metal_substitution_compatible(source_metal: str, target_metal: str) -> bool:
    """Check if a source metal site can serve as geometric template for target metal.

    Compatibility requires:
    1. Both metals in METAL_DATABASE
    2. Same HSAB class (hard↔hard, borderline↔borderline, soft↔soft)
    3. Same dominant donor preference (O-preferring, N-preferring, or S-preferring)

    Examples:
        MG → TB: both hard, both O-preferring → True
        CA → EU: both hard, both O-preferring → True
        ZN → TB: borderline vs hard → False
        ZN → FE: both borderline → True

    Args:
        source_metal: Metal found in the PDB (e.g., "MG")
        target_metal: Metal the user wants (e.g., "TB")

    Returns:
        True if the source site geometry is suitable as a template for the target
    """
    source_metal = source_metal.upper()
    target_metal = target_metal.upper()

    if source_metal == target_metal:
        return True

    if source_metal not in METAL_DATABASE or target_metal not in METAL_DATABASE:
        return False

    # Get HSAB classes
    src_info = METAL_DATABASE[source_metal]
    tgt_info = METAL_DATABASE[target_metal]

    src_ox = src_info.get("default_oxidation", 2)
    tgt_ox = tgt_info.get("default_oxidation", 2)

    try:
        src_hsab = get_hsab_class(source_metal, src_ox)
        tgt_hsab = get_hsab_class(target_metal, tgt_ox)
    except (ValueError, KeyError):
        return False

    # Rule 1: Same HSAB class
    if src_hsab != tgt_hsab:
        return False

    # Rule 2: Same dominant donor preference
    src_donors = src_info.get("preferred_donors", {}).get(src_ox, {}).get("catalytic", {})
    tgt_donors = tgt_info.get("preferred_donors", {}).get(tgt_ox, {}).get("catalytic", {})

    if src_donors and tgt_donors:
        src_best = max(src_donors, key=src_donors.get) if src_donors else None
        tgt_best = max(tgt_donors, key=tgt_donors.get) if tgt_donors else None
        if src_best != tgt_best:
            return False

    return True


# =============================================================================
# BACKBONE CONTINUITY CHECK
# =============================================================================

def check_backbone_continuity(pdb_content: str, max_ca_dist: float = 4.2) -> dict:
    """Check if a backbone PDB has continuous chain (no breaks).

    Validates that all adjacent CA-CA distances within the same chain
    are within the normal peptide bond range. Normal CA-CA distance is ~3.8A.

    Args:
        pdb_content: PDB file content as string
        max_ca_dist: Maximum allowed CA-CA distance (default 4.2A).
            Normal peptide bond CA-CA is ~3.8A; 4.2A allows for slight
            distortions while catching real chain breaks (>5A).

    Returns:
        dict with:
          - continuous (bool): True if no breaks found
          - num_breaks (int): Number of chain breaks
          - breaks (list): List of (resnum1, resnum2, distance) tuples
          - num_ca (int): Total number of CA atoms
          - mean_ca_dist (float): Mean adjacent CA-CA distance
          - max_observed (float): Maximum adjacent CA-CA distance
    """
    import math

    # Extract CA atoms per chain
    ca_by_chain = {}
    for line in pdb_content.split('\n'):
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            chain = line[21]
            try:
                resnum = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                ca_by_chain.setdefault(chain, []).append((resnum, x, y, z))
            except (ValueError, IndexError):
                continue

    breaks = []
    all_dists = []
    total_ca = 0

    for chain, cas in ca_by_chain.items():
        # Sort by residue number
        cas.sort(key=lambda c: c[0])
        total_ca += len(cas)

        for i in range(1, len(cas)):
            r1, x1, y1, z1 = cas[i - 1]
            r2, x2, y2, z2 = cas[i]
            d = math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
            all_dists.append(d)
            if d > max_ca_dist:
                breaks.append((r1, r2, round(d, 2)))

    return {
        "continuous": len(breaks) == 0,
        "num_breaks": len(breaks),
        "breaks": breaks,
        "num_ca": total_ca,
        "mean_ca_dist": round(sum(all_dists) / len(all_dists), 2) if all_dists else 0.0,
        "max_observed": round(max(all_dists), 2) if all_dists else 0.0,
    }


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class ScaffoldResult:
    """Result from scaffolding workflow."""
    # Core outputs
    motif_pdb: str                         # Extracted motif ready for RFD3
    motif_residues: List[str] = field(default_factory=list)  # Residues to preserve (e.g., ["A10-15"])
    contig: str = ""                       # Full contig string for RFD3 (DEPRECATED - use length)

    # NEW: Enzyme Scaffold Design approach (recommended)
    length: str = ""                       # Total scaffold length (e.g., "150" or "120-180")
    unindex: str = ""                      # Catalytic residues marked as flexible (e.g., "A1,A2,A3")
    ligand_codes: str = ""                 # Comma-separated ligand 3-letter codes (e.g., "PQQ,CA")

    # Conditioning parameters for RFD3
    fixed_atoms: Dict[str, str] = field(default_factory=dict)      # {"X1": "all", "L1": "all"}
    rasa_targets: Dict[str, str] = field(default_factory=dict)     # Burial targets
    hbond_acceptors: Dict[str, str] = field(default_factory=dict)  # H-bond acceptor atoms
    hbond_donors: Dict[str, str] = field(default_factory=dict)     # H-bond donor atoms

    # Conservation data (from ConSurf-style analysis)
    conservation_fixed: List[str] = field(default_factory=list)    # Highly conserved residues to fix
    conservation_grades: Dict[int, int] = field(default_factory=dict)  # Position -> grade (1-9)

    # Metadata
    source_info: Dict[str, Any] = field(default_factory=dict)      # Source structure info
    coordinating_residues: List[Dict] = field(default_factory=list)  # Coordination shell

    # Minimal motif selection result (if used)
    motif_result: Optional[Any] = None  # MinimalMotifResult when use_minimal_motif=True

    # Status
    success: bool = True
    error_message: str = ""

    def apply_conservation_to_fixed_atoms(
        self,
        threshold: int = 3,
        fix_type: str = "CA",
    ) -> None:
        """
        Add highly conserved residues to fixed_atoms based on conservation grades.

        Args:
            threshold: Maximum grade to consider fixed (1-3 = highly conserved)
            fix_type: What to fix ("CA" for backbone, "all" for everything)
        """
        for pos, grade in self.conservation_grades.items():
            if grade <= threshold:
                res_id = f"A{pos}"
                if res_id not in self.fixed_atoms:
                    self.fixed_atoms[res_id] = fix_type
                    self.conservation_fixed.append(res_id)


# =============================================================================
# SCAFFOLDING WORKFLOW
# =============================================================================

class ScaffoldingWorkflow:
    """
    Workflow for scaffolding from existing PDB structures.

    This class wraps the existing metal_site_fetcher infrastructure to provide
    a clean interface for the AI design pipeline.

    Usage:
        workflow = ScaffoldingWorkflow()
        result = await workflow.run(
            pdb_id="4CVB",
            metal="CA",
            ligand_code="PQQ",
            chain_length="100-130"
        )

        if result.success:
            # Use result.motif_pdb, result.contig, etc. for RFD3
            pass
    """

    def __init__(
        self,
        default_coordination_cutoff: float = 3.0,
        default_ligand_contact_cutoff: float = 4.5,
    ):
        """
        Initialize the scaffolding workflow.

        Args:
            default_coordination_cutoff: Default distance cutoff for metal coordination (3.0 Å)
            default_ligand_contact_cutoff: Default distance cutoff for all ligand contacts (4.5 Å)
        """
        self.default_cutoff = default_coordination_cutoff
        self.default_ligand_cutoff = default_ligand_contact_cutoff

    async def run(
        self,
        pdb_id: str,
        metal: Optional[str] = None,
        ligand_code: Optional[str] = None,
        chain_length: Optional[str] = None,  # None = auto-calculate from active site geometry
        coordination_cutoff: Optional[float] = None,
        include_second_shell: bool = False,
        include_all_ligand_contacts: bool = False,
        ligand_contact_cutoff: Optional[float] = None,
        fixed_atom_type: str = "BKBN",  # What to fix on catalytic residues: BKBN, ALL, TIP, or empty
        use_minimal_motif: bool = False,  # NEW: use evidence-based minimal motif selection
        enzyme_class: Optional[str] = None,  # NEW: enzyme class for chemistry evidence
        use_conservation: bool = False,  # NEW: include ConSurf conservation analysis
    ) -> ScaffoldResult:
        """
        Run scaffolding workflow.

        Steps:
        1. Normalize ligand code (human name → PDB 3-letter code)
        2. Fetch PDB structure from RCSB
        3. Find active site (with metal substitution fallback)
        4. Extract motif with fixed atoms
        5. Build scaffolding contig
        6. Determine conditioning parameters

        Supports theozyme scaffolding: if the requested metal is not in the
        PDB but a compatible metal is found (same HSAB class + donor preference),
        the compatible metal's site is extracted and the metal is rewritten
        in the output theozyme PDB.

        Args:
            pdb_id: 4-character PDB ID (e.g., "4CVB")
            metal: Metal element symbol (e.g., "CA", "ZN", "TB")
            ligand_code: Ligand name or 3-letter code (e.g., "citrate", "PQQ", "CIT")
            chain_length: Chain length range for designed linkers (e.g., "80-120")
            coordination_cutoff: Distance cutoff for metal coordination (default: 3.0 Å)
            include_second_shell: Include second-shell residues in motif
            include_all_ligand_contacts: If True, include ALL residues contacting ligand (default: False)
            ligand_contact_cutoff: Distance cutoff for ligand contacts (default: 4.5 Å)

        Returns:
            ScaffoldResult with motif PDB, contig, and conditioning parameters
        """
        cutoff = coordination_cutoff or self.default_cutoff
        ligand_cutoff = ligand_contact_cutoff or self.default_ligand_cutoff

        # Step 0: Normalize ligand code (human name → PDB 3-letter code)
        original_ligand_name = ligand_code
        if ligand_code:
            ligand_code = normalize_ligand_code(ligand_code)
            if ligand_code != original_ligand_name:
                logger.info(f"Ligand code normalized: '{original_ligand_name}' → '{ligand_code}'")

        logger.info(f"ScaffoldingWorkflow: PDB={pdb_id}, Metal={metal}, Ligand={ligand_code}, all_contacts={include_all_ligand_contacts}")

        # Track metal substitution for theozyme rewriting
        substituted_from = None  # Source metal found in PDB (e.g., "MG")
        target_metal = metal     # Metal the user requested (e.g., "TB")

        try:
            # Step 1: Fetch structure from RCSB
            pdb_content = await asyncio.to_thread(_fetch_pdb_content, pdb_id)
            if not pdb_content:
                return ScaffoldResult(
                    motif_pdb="",
                    success=False,
                    error_message=f"Could not fetch PDB: {pdb_id}",
                )

            # Step 2: Find active site (with metal substitution fallback)
            active_site = None

            if metal and ligand_code:
                # Try direct match first
                active_site = await asyncio.to_thread(
                    find_metal_ligand_active_site,
                    pdb_id, metal, ligand_code, cutoff
                )

                # Fallback: metal substitution if direct match fails
                if not active_site:
                    sub_result = self._find_site_with_metal_substitution(
                        pdb_content, pdb_id, metal, ligand_code, cutoff
                    )
                    if sub_result:
                        active_site = sub_result["active_site"]
                        substituted_from = sub_result["source_metal"]
                        # Use source metal for PDB parsing, target metal for output
                        metal = sub_result["source_metal"]
                        logger.info(
                            f"Metal substitution: {substituted_from} → {target_metal} "
                            f"(HSAB-compatible, site extracted from {pdb_id})"
                        )

            elif metal:
                # Extract metal coordination only
                coord_info = await asyncio.to_thread(
                    extract_metal_coordination,
                    pdb_content, metal, cutoff
                )
                if coord_info.get("coordination_number", 0) > 0:
                    active_site = {
                        "pdb_id": pdb_id.upper(),
                        "metal": metal,
                        "metal_chain": coord_info.get("metal_chain", "A"),
                        "metal_resnum": coord_info.get("metal_resnum", 0),
                        "metal_coords": coord_info.get("metal_coords"),
                        "coordinating_atoms": coord_info.get("coordinating_atoms", []),
                        "coordination_number": coord_info.get("coordination_number", 0),
                    }

                # Fallback: try substitution for metal-only queries too
                if not active_site:
                    sub_result = self._find_site_with_metal_substitution(
                        pdb_content, pdb_id, target_metal, ligand_code, cutoff
                    )
                    if sub_result:
                        active_site = sub_result["active_site"]
                        substituted_from = sub_result["source_metal"]
                        metal = sub_result["source_metal"]
            else:
                # Try to auto-detect metal from PDB
                active_site = self._auto_detect_active_site(pdb_content, cutoff)

            if not active_site:
                # Provide helpful error listing what's actually in the PDB
                metals_found, ligands_found = scan_pdb_hetatm(pdb_content)
                error_parts = [f"No active site found in {pdb_id}."]
                if metals_found:
                    error_parts.append(f"Metals in PDB: {', '.join(sorted(metals_found))}")
                if ligands_found:
                    error_parts.append(f"Ligands in PDB: {', '.join(sorted(ligands_found))}")
                if target_metal and target_metal not in metals_found:
                    error_parts.append(
                        f"Requested metal {target_metal} not found and no "
                        f"HSAB-compatible substitute available."
                    )
                return ScaffoldResult(
                    motif_pdb="",
                    success=False,
                    error_message=" ".join(error_parts),
                )

            # Step 3: Extract coordinating residues
            coord_residues = active_site.get("coordinating_atoms", [])

            # Group into unique residues (metal coordination)
            unique_residues = self._get_unique_residues(coord_residues)

            # Step 3b: If include_all_ligand_contacts, find ALL residues contacting ligand
            if include_all_ligand_contacts and ligand_code:
                ligand_contacts = self._find_all_ligand_contacts(
                    pdb_content, ligand_code, ligand_cutoff
                )
                # Merge with coordination residues (avoid duplicates)
                existing_keys = {(r["chain"], r["resnum"]) for r in unique_residues}
                for contact in ligand_contacts:
                    key = (contact["chain"], contact["resnum"])
                    if key not in existing_keys:
                        unique_residues.append(contact)
                        existing_keys.add(key)
                logger.info(f"Extended pocket: {len(unique_residues)} residues (was {len(coord_residues)} from metal coordination)")

            # Step 3c: Minimal motif selection (evidence-based)
            # When enabled, filters unique_residues down to only evidence-supported
            # residues (Tier 1-3), with per-residue atom fix types (TIP/CA).
            minimal_motif_result = None
            if use_minimal_motif:
                try:
                    from minimal_motif_selector import MinimalMotifSelector
                    selector = MinimalMotifSelector(use_conservation=use_conservation)
                    minimal_motif_result = await selector.select(
                        pdb_id=pdb_id,
                        pdb_content=pdb_content,
                        ligand_code=ligand_code,
                        metal_type=metal,
                        enzyme_class=enzyme_class,
                        chain=active_site.get("metal_chain", "A"),
                        pocket_residues=unique_residues,
                    )
                    # Replace unique_residues with only motif residues (Tier 1-3)
                    motif_residue_keys = minimal_motif_result.get_motif_residue_keys()
                    motif_unique_residues = [
                        r for r in unique_residues
                        if (r.get("chain", "A"), r.get("resnum", 0)) in motif_residue_keys
                    ]
                    logger.info(
                        f"Minimal motif: {len(motif_unique_residues)} residues "
                        f"(was {len(unique_residues)} full pocket), "
                        f"{minimal_motif_result.total_fixed_atoms} fixed atoms"
                    )
                    unique_residues = motif_unique_residues
                except Exception as e:
                    logger.warning(f"Minimal motif selection failed, using full pocket: {e}")

            # Convert to residue ranges for contig
            motif_residues, all_residue_keys = self._residues_to_ranges(unique_residues)

            # Step 4: Build theozyme PDB (catalytic residues + ligand/metal)
            # This is the RECOMMENDED approach: extract theozyme, let RFD3 design entire scaffold
            # NO blocker residues needed - RFD3 respects ligand when using proper parameters
            # If metal was substituted, rewrite HETATM to use the target metal
            motif_pdb, scaffold_residue_ids = self._extract_theozyme_pdb(
                pdb_content,
                active_site,
                unique_residues,
                metal,
                ligand_code,
                all_residue_keys=all_residue_keys,
                target_metal=target_metal if substituted_from else None,
            )

            # Step 5: Build parameters for Enzyme Scaffold Design approach
            # Use `length` (NOT contig) - let RFD3 design entire scaffold from scratch
            # Use `unindex` to mark catalytic residues as flexible during diffusion
            if chain_length is None:
                chain_length = self._calculate_scaffold_length(
                    motif_pdb, len(scaffold_residue_ids),
                )
            length = chain_length  # Total scaffold length

            # Build unindex string: comma-separated residue IDs for flexible catalytic residues
            # Format: "A1,A2,A3" (after renumbering)
            unindex = ",".join(scaffold_residue_ids)

            # Build ligand_codes string: comma-separated 3-letter codes
            # Use target_metal (user's requested metal) for output, not source metal
            output_metal = target_metal if substituted_from else metal
            ligand_codes_list = []
            if ligand_code:
                ligand_codes_list.append(ligand_code.upper())
            if output_metal:
                ligand_codes_list.append(output_metal.upper())
            ligand_codes_str = ",".join(ligand_codes_list)

            # Step 6: Determine conditioning parameters
            # For Enzyme Scaffold Design, we fix LIGAND atoms and optionally catalytic residue atoms
            fixed_atoms = {}

            # Fix ligand position (always recommended)
            if ligand_code:
                fixed_atoms[ligand_code.upper()] = "ALL"
            if output_metal:
                fixed_atoms[output_metal.upper()] = "ALL"

            # Fix catalytic residue atoms — use per-residue atoms from motif selector
            # when available, otherwise fall back to uniform fixed_atom_type
            if minimal_motif_result and use_minimal_motif:
                # Build mapping from renumbered IDs to original (chain, resnum)
                # _extract_theozyme_pdb sorts all_residue_keys and numbers A1, A2, ...
                sorted_keys = sorted(all_residue_keys)
                # Build lookup: (chain, resnum) -> MotifResidue
                motif_by_key = {
                    (r.chain, r.resnum): r
                    for r in minimal_motif_result.motif_residues
                }
                motif_fixed_count = 0
                for idx, (orig_chain, orig_resnum) in enumerate(sorted_keys, start=1):
                    res_id = f"A{idx}"
                    mr = motif_by_key.get((orig_chain, orig_resnum))
                    if mr and mr.fix_type == "TIP" and mr.functional_atoms:
                        fixed_atoms[res_id] = ",".join(mr.functional_atoms)
                        motif_fixed_count += 1
                    elif mr and mr.fix_type == "CA":
                        fixed_atoms[res_id] = "CA"
                        motif_fixed_count += 1
                    elif mr and mr.fix_type == "BKBN":
                        fixed_atoms[res_id] = "N,CA,C,O"
                        motif_fixed_count += 1
                    # fix_type == "NONE" or not in motif → not fixed
                logger.info(f"Per-residue fixed atoms from motif selector: {motif_fixed_count} entries")
            elif fixed_atom_type and scaffold_residue_ids:
                # Legacy: uniform fix type for all catalytic residues
                for res_id in scaffold_residue_ids:
                    fixed_atoms[res_id] = fixed_atom_type

            # RASA targets: bury metal and ligand to create proper pocket
            rasa_targets = self._determine_rasa_targets_renumbered(output_metal, ligand_code)

            # H-bond conditioning from ligand chemistry
            hbond_acceptors, hbond_donors = self._determine_hbond_conditioning(
                pdb_content, ligand_code
            )

            # LEGACY: Also build contig for backward compatibility
            # (but length + unindex is the recommended approach)
            contig = self._build_renumbered_contig(
                scaffold_fixed_residues=scaffold_residue_ids,
                chain_length=chain_length,
                design_mode="dimer",
            )

            return ScaffoldResult(
                motif_pdb=motif_pdb,
                motif_residues=motif_residues,
                contig=contig,  # LEGACY - use length instead
                # NEW: Enzyme Scaffold Design parameters
                length=length,
                unindex=unindex,
                ligand_codes=ligand_codes_str,
                # Conditioning
                fixed_atoms=fixed_atoms,
                rasa_targets=rasa_targets,
                hbond_acceptors=hbond_acceptors,
                hbond_donors=hbond_donors,
                coordinating_residues=unique_residues,
                # Minimal motif result (if used)
                motif_result=minimal_motif_result,
                source_info={
                    "pdb_id": pdb_id.upper(),
                    "metal": output_metal,
                    "ligand": ligand_code,
                    "original_ligand_name": original_ligand_name,
                    "metal_chain": active_site.get("metal_chain"),
                    "metal_resnum": active_site.get("metal_resnum"),
                    "coordination_number": active_site.get("coordination_number", 0),
                    "protein_donors": active_site.get("protein_donors", []),
                    "ligand_donors": active_site.get("ligand_donors", []),
                    "include_all_contacts": include_all_ligand_contacts,
                    "total_pocket_residues": len(unique_residues),
                    "catalytic_residue_ids": scaffold_residue_ids,
                    "design_approach": "enzyme_scaffold",
                    "use_minimal_motif": use_minimal_motif,
                    # Metal substitution metadata
                    "metal_substituted": substituted_from is not None,
                    "substituted_from": substituted_from,
                    "substituted_to": target_metal if substituted_from else None,
                },
                success=True,
            )

        except Exception as e:
            logger.error(f"ScaffoldingWorkflow error: {e}")
            return ScaffoldResult(
                motif_pdb="",
                success=False,
                error_message=str(e),
            )

    def _auto_detect_active_site(
        self,
        pdb_content: str,
        cutoff: float,
    ) -> Optional[Dict[str, Any]]:
        """
        Try to auto-detect active site from PDB content.

        Searches for common metals in order of likelihood.
        """
        # Common active site metals
        common_metals = ["ZN", "FE", "MG", "CA", "MN", "CU", "CO", "NI"]

        for metal in common_metals:
            coord_info = extract_metal_coordination(pdb_content, metal, cutoff)
            if coord_info.get("coordination_number", 0) >= 3:
                return {
                    "metal": metal,
                    "metal_chain": coord_info.get("metal_chain", "A"),
                    "metal_resnum": coord_info.get("metal_resnum", 0),
                    "metal_coords": coord_info.get("metal_coords"),
                    "coordinating_atoms": coord_info.get("coordinating_atoms", []),
                    "coordination_number": coord_info.get("coordination_number", 0),
                }

        return None

    def _find_site_with_metal_substitution(
        self,
        pdb_content: str,
        pdb_id: str,
        target_metal: str,
        ligand_code: Optional[str],
        cutoff: float,
    ) -> Optional[Dict[str, Any]]:
        """Find an active site using a compatible substitute metal.

        When the user requests a metal (e.g., TB) that isn't in the PDB,
        this method scans the PDB for metals with compatible HSAB chemistry
        and extracts the site using the substitute. The caller is responsible
        for rewriting the metal in the output theozyme.

        Args:
            pdb_content: Raw PDB file content (already fetched)
            pdb_id: PDB ID for logging and API calls
            target_metal: Metal the user wants (e.g., "TB")
            ligand_code: Normalized ligand 3-letter code (e.g., "CIT") or None
            cutoff: Coordination distance cutoff in Angstroms

        Returns:
            Dict with keys:
                "active_site": active site dict from find_metal_ligand_active_site
                "source_metal": the substitute metal found in PDB (e.g., "MG")
                "target_metal": the user's requested metal (e.g., "TB")
                "cn_warning": optional warning about coordination number difference
            Or None if no compatible substitute found.
        """
        metals_found, ligands_found = scan_pdb_hetatm(pdb_content)

        if not metals_found:
            logger.info(f"No metals found in {pdb_id} for substitution")
            return None

        # Check each found metal for HSAB compatibility
        compatible = []
        for source_metal in sorted(metals_found):
            if source_metal == target_metal:
                continue  # Already tried direct match
            if is_metal_substitution_compatible(source_metal, target_metal):
                compatible.append(source_metal)

        if not compatible:
            logger.info(
                f"No HSAB-compatible metals for {target_metal} in {pdb_id}. "
                f"Found: {sorted(metals_found)}"
            )
            return None

        logger.info(
            f"Found {len(compatible)} compatible substitute(s) for {target_metal}: "
            f"{compatible}"
        )

        # Try each compatible metal (prefer the first that works with the ligand)
        for source_metal in compatible:
            active_site = None

            if ligand_code:
                active_site = find_metal_ligand_active_site(
                    pdb_id, source_metal, ligand_code, cutoff
                )
            else:
                coord_info = extract_metal_coordination(
                    pdb_content, source_metal, cutoff
                )
                if coord_info.get("coordination_number", 0) > 0:
                    active_site = {
                        "pdb_id": pdb_id.upper(),
                        "metal": source_metal,
                        "metal_chain": coord_info.get("metal_chain", "A"),
                        "metal_resnum": coord_info.get("metal_resnum", 0),
                        "metal_coords": coord_info.get("metal_coords"),
                        "coordinating_atoms": coord_info.get("coordinating_atoms", []),
                        "coordination_number": coord_info.get("coordination_number", 0),
                    }

            if active_site:
                # Build CN warning if ranges differ significantly
                cn_warning = None
                try:
                    src_cn = get_coordination_number_range(
                        source_metal,
                        METAL_DATABASE[source_metal].get("default_oxidation", 2),
                    )
                    tgt_cn = get_coordination_number_range(
                        target_metal,
                        METAL_DATABASE[target_metal].get("default_oxidation", 2),
                    )
                    if src_cn and tgt_cn and tgt_cn[0] > src_cn[1]:
                        cn_warning = (
                            f"CN mismatch: {source_metal} typically {src_cn[0]}-{src_cn[1]}, "
                            f"{target_metal} needs {tgt_cn[0]}-{tgt_cn[1]}. "
                            f"RFD3/LigandMPNN will design additional coordinating residues."
                        )
                        logger.warning(cn_warning)
                except (KeyError, ValueError):
                    pass

                logger.info(
                    f"Metal substitution successful: {source_metal} → {target_metal} "
                    f"in {pdb_id} (CN={active_site.get('coordination_number', '?')})"
                )

                return {
                    "active_site": active_site,
                    "source_metal": source_metal,
                    "target_metal": target_metal,
                    "cn_warning": cn_warning,
                }

        logger.info(f"Compatible metals found but none coordinate {ligand_code} in {pdb_id}")
        return None

    def _get_unique_residues(self, coord_atoms: List[Dict]) -> List[Dict]:
        """
        Get unique residues from coordinating atoms.

        Returns list of unique residue dicts with chain, resnum, res_name.
        """
        seen = set()
        unique = []

        for atom in coord_atoms:
            key = (atom.get("chain_id", "A"), atom.get("res_seq", 0))
            if key not in seen:
                seen.add(key)
                unique.append({
                    "chain": atom.get("chain_id", "A"),
                    "resnum": atom.get("res_seq", 0),
                    "res_name": atom.get("res_name", "UNK"),
                })

        return unique

    def _find_all_ligand_contacts(
        self,
        pdb_content: str,
        ligand_code: str,
        cutoff: float = 4.5,
    ) -> List[Dict]:
        """
        Find ALL protein residues contacting the ligand within cutoff distance.

        This captures the full binding pocket, including:
        - Hydrogen bond donors/acceptors
        - Hydrophobic contacts
        - Pi-stacking interactions
        - Salt bridges

        Args:
            pdb_content: Full PDB file content
            ligand_code: 3-letter ligand code (e.g., "PQQ")
            cutoff: Distance cutoff in Angstroms (default 4.5Å)

        Returns:
            List of residue dicts with chain, resnum, res_name, min_distance
        """
        from math import sqrt

        ligand_upper = ligand_code.upper()

        # Parse ligand atoms
        ligand_atoms = []
        for line in pdb_content.split('\n'):
            if line.startswith('HETATM') and line[17:20].strip() == ligand_upper:
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    atom_name = line[12:16].strip()
                    ligand_atoms.append({'name': atom_name, 'coords': (x, y, z)})
                except (ValueError, IndexError):
                    continue

        if not ligand_atoms:
            logger.warning(f"No atoms found for ligand {ligand_code}")
            return []

        # Find protein residues within cutoff of any ligand atom
        contacts: Dict[Tuple[str, int], Dict] = {}

        for line in pdb_content.split('\n'):
            if line.startswith('ATOM'):
                try:
                    chain = line[21]
                    res_name = line[17:20].strip()
                    res_seq = int(line[22:26].strip())
                    atom_name = line[12:16].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except (ValueError, IndexError):
                    continue

                # Check distance to all ligand atoms
                for lig_atom in ligand_atoms:
                    lx, ly, lz = lig_atom['coords']
                    dist = sqrt((x - lx)**2 + (y - ly)**2 + (z - lz)**2)

                    if dist <= cutoff:
                        key = (chain, res_seq)
                        if key not in contacts:
                            contacts[key] = {
                                "chain": chain,
                                "resnum": res_seq,
                                "res_name": res_name,
                                "min_distance": dist,
                                "contact_atoms": [],
                            }
                        if dist < contacts[key]["min_distance"]:
                            contacts[key]["min_distance"] = dist
                        contacts[key]["contact_atoms"].append(
                            f"{atom_name}-{lig_atom['name']}@{dist:.2f}A"
                        )

        # Sort by residue number and return
        result = sorted(contacts.values(), key=lambda x: (x["chain"], x["resnum"]))
        logger.info(f"Found {len(result)} residues within {cutoff}A of {ligand_code}")
        return result

    def _residues_to_ranges(
        self,
        residues: List[Dict],
        gap_threshold: int = 3,
    ) -> Tuple[List[str], Set[Tuple[str, int]]]:
        """
        Convert residue list to range strings like 'A10-15'.

        IMPORTANT: Uses gap_threshold=3 to only merge consecutive or nearly
        consecutive residues. This ensures the motif PDB contains all residues
        referenced in the contig ranges.

        NOTE: Excludes ligand/metal residues (resnum >= 1000) - these should be
        handled via pdb_content and fixed_atoms, NOT via contig.

        Args:
            residues: List of residue dicts
            gap_threshold: Max gap to merge into single range (default 3)

        Returns:
            Tuple of (List of range strings, Set of all residue keys to include)
        """
        if not residues:
            return [], set()

        # Group by chain, EXCLUDING ligand/metal (resnum >= 1000)
        chains: Dict[str, List[int]] = {}
        for res in residues:
            chain = res.get("chain", "A")
            resnum = res.get("resnum", 0)
            # Skip ligand/metal residues - they go via pdb_content, not contig
            if resnum >= 1000:
                continue
            if chain not in chains:
                chains[chain] = []
            chains[chain].append(resnum)

        ranges = []
        all_residue_keys: Set[Tuple[str, int]] = set()

        for chain, resnums in chains.items():
            resnums = sorted(set(resnums))
            if not resnums:
                continue

            # Group residues within gap_threshold
            start = resnums[0]
            end = resnums[0]

            for i in range(1, len(resnums)):
                gap = resnums[i] - end
                if gap <= gap_threshold:
                    # Extend the range (includes gap residues)
                    end = resnums[i]
                else:
                    # Start new range - add all residues in range
                    ranges.append(f"{chain}{start}-{end}")
                    for r in range(start, end + 1):
                        all_residue_keys.add((chain, r))
                    start = end = resnums[i]

            # Add final range
            ranges.append(f"{chain}{start}-{end}")
            for r in range(start, end + 1):
                all_residue_keys.add((chain, r))

        return ranges, all_residue_keys

    def _build_scaffolding_contig(
        self,
        motif_residues: List[str],
        chain_length: str,
    ) -> str:
        """
        Build motif scaffolding contig with smart linker sizing.

        Strategy:
        - Parse the original residue ranges to understand gaps
        - Use gap-proportional linker sizes (smaller gaps = shorter linkers)
        - Limit total designed length to reasonable bounds

        Example output: "0 15-25/A104-220/0 10-20/A270-335/0 15-25"

        Args:
            motif_residues: List of residue range strings (e.g., ["A104-220", "A270-335"])
            chain_length: Chain length range for designed regions (e.g., "80-120")

        Returns:
            Complete contig string for RFD3
        """
        if not motif_residues:
            return chain_length

        # Filter to protein residues only (exclude ligand/metal with high resnums)
        protein_ranges = []
        for r in motif_residues:
            # Parse range like "A104-220"
            match = re.match(r'([A-Z])(\d+)-(\d+)', r)
            if match:
                chain, start, end = match.groups()
                start, end = int(start), int(end)
                # Skip if this looks like ligand (resnum > 1000)
                if start < 1000:
                    protein_ranges.append((chain, start, end, r))

        if not protein_ranges:
            return chain_length

        # Sort by start position
        protein_ranges.sort(key=lambda x: x[1])

        parts = []

        # N-terminal extension (short) - use comma notation for RFD3
        parts.append("15-30")

        # Add ranges with gap-proportional linkers
        for i, (chain, start, end, range_str) in enumerate(protein_ranges):
            parts.append(range_str)

            if i < len(protein_ranges) - 1:
                # Calculate gap to next range
                next_start = protein_ranges[i + 1][1]
                gap = next_start - end

                # Proportional linker: roughly 1/4 of the gap, min 5, max 30
                linker_min = max(5, min(15, gap // 6))
                linker_max = max(10, min(40, gap // 4))
                parts.append(f"{linker_min}-{linker_max}")

        # C-terminal extension (short)
        parts.append("15-30")

        # NOTE: Ligand/metal residues are NOT added to contig - they are
        # handled via pdb_content and fixed_atoms parameters

        # Build final contig - use comma notation for RFD3
        contig = ",".join(parts)

        return contig

    def _calculate_scaffold_length(
        self,
        theozyme_pdb: str,
        num_catalytic_residues: int,
    ) -> str:
        """Calculate optimal scaffold length from active site geometry.

        Uses three rules to determine minimum viable scaffold size:
        1. Motif fraction: catalytic residues should be <= 20% of total protein
        2. Spatial coverage: need enough protein to enclose the active site span
        3. Large-motif bonus: active sites with >15 residues need extra scaffold

        Args:
            theozyme_pdb: PDB content of the theozyme (catalytic residues + ligand/metal)
            num_catalytic_residues: Number of catalytic protein residues

        Returns:
            Length range string like "130-160" for the RFD3 length parameter
        """
        import math

        # Extract CA coordinates from theozyme
        ca_coords = []
        for line in theozyme_pdb.split('\n'):
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    ca_coords.append((x, y, z))
                except ValueError:
                    continue

        # Compute max pairwise CA-CA distance (active site span)
        max_ca_dist = 0.0
        for i in range(len(ca_coords)):
            for j in range(i + 1, len(ca_coords)):
                dx = ca_coords[i][0] - ca_coords[j][0]
                dy = ca_coords[i][1] - ca_coords[j][1]
                dz = ca_coords[i][2] - ca_coords[j][2]
                d = math.sqrt(dx*dx + dy*dy + dz*dz)
                max_ca_dist = max(max_ca_dist, d)

        n_cat = num_catalytic_residues

        # Rule 1: Motif fraction - catalytic residues should be <= 20% of total
        n_fraction = n_cat / 0.20  # e.g., 21/0.20 = 105

        # Rule 2: Spatial enclosure - need enough protein to wrap around active site
        # A globular protein of N residues has diameter ~ 5 * N^0.4 Angstroms
        # Invert: N ~ (diameter / 5)^2.5
        # Add shell of 8A beyond active site for proper folding
        if max_ca_dist > 0:
            required_diameter = max_ca_dist + 16  # 8A shell on each side
            n_spatial = (required_diameter / 5.0) ** 2.5
        else:
            n_spatial = n_cat * 4

        # Rule 3: Large motif bonus - active sites with >15 residues need extra scaffold
        # to form adequate secondary structure around the dispersed catalytic residues
        if n_cat > 15:
            n_large_motif = n_cat * 6
        elif n_cat > 10:
            n_large_motif = n_cat * 5
        else:
            n_large_motif = n_cat * 4

        # Rule 4: Absolute minimum for any foldable protein
        n_min = 80

        # Combine: take the maximum of all rules
        target = max(n_fraction, n_spatial, n_large_motif, n_min)

        # Cap at 250 (RFD3 performance degrades above this for single-chain designs)
        target = min(target, 250)

        # Build range string: round to nearest 10, allow 20% flexibility above
        low = int(round(target / 10) * 10)
        high = int(round(target * 1.2 / 10) * 10)
        high = min(high, 250)

        # Ensure low < high
        if low >= high:
            high = low + 20

        logger.info(
            f"Scaffold length auto-calculated: {low}-{high} "
            f"(n_cat={n_cat}, max_span={max_ca_dist:.1f}A, "
            f"rule1={n_fraction:.0f}, rule2={n_spatial:.0f}, rule3={n_large_motif})"
        )

        return f"{low}-{high}"

    def _extract_theozyme_pdb(
        self,
        pdb_content: str,
        active_site: Dict,
        coord_residues: List[Dict],
        metal: Optional[str],
        ligand_code: Optional[str],
        all_residue_keys: Optional[Set[Tuple[str, int]]] = None,
        target_metal: Optional[str] = None,
    ) -> Tuple[str, List[str]]:
        """
        Extract theozyme PDB for Enzyme Scaffold Design approach.

        Creates a clean theozyme with:
        - Catalytic residues (renumbered to chain A, sequential from 1)
        - Metal (on chain M, residue 1, or using 3-letter code)
        - Ligand (on chain L, residue 1, or using 3-letter code)

        NO blocker residues - the Enzyme Scaffold Design approach doesn't need them.
        RFD3 designs the entire scaffold from scratch around the theozyme.

        Args:
            pdb_content: Full PDB content
            active_site: Active site dict from find_metal_ligand_active_site
            coord_residues: List of coordinating residue dicts
            metal: Metal code as found in PDB (e.g., "MG" — the source metal)
            ligand_code: Ligand code (e.g., "PQQ")
            all_residue_keys: Set of (chain, resnum) tuples to include
            target_metal: If set, rewrite metal HETATM records to this symbol
                         (e.g., "TB" when substituting MG → TB)

        Returns:
            Tuple of (theozyme_pdb_content, list_of_residue_ids)
            where residue_ids are like ["A1", "A2", "A3"]
        """
        # Use all_residue_keys if provided, otherwise fall back to coord_residues
        if all_residue_keys is not None:
            residue_keys_to_keep = all_residue_keys
        else:
            residue_keys_to_keep = set()
            for res in coord_residues:
                key = (res.get("chain", "A"), res.get("resnum", 0))
                residue_keys_to_keep.add(key)

        metal_upper = metal.upper() if metal else None
        ligand_upper = ligand_code.upper() if ligand_code else None

        # Collect ATOM and HETATM lines
        atom_lines = []
        metal_lines = []
        ligand_lines = []

        for line in pdb_content.split('\n'):
            if line.startswith('HETATM'):
                res_name = line[17:20].strip()
                if metal_upper and res_name == metal_upper:
                    metal_lines.append(line)
                elif ligand_upper and res_name == ligand_upper:
                    ligand_lines.append(line)
            elif line.startswith('ATOM'):
                chain = line[21]
                try:
                    resnum = int(line[22:26].strip())
                except ValueError:
                    continue
                if (chain, resnum) in residue_keys_to_keep:
                    atom_lines.append((line, chain, resnum))

        # Build old->new residue mapping (sequential numbering from 1)
        old_to_new_resnum = {}
        sorted_keys = sorted(residue_keys_to_keep)
        for idx, (old_chain, old_resnum) in enumerate(sorted_keys, start=1):
            old_to_new_resnum[(old_chain, old_resnum)] = idx

        # Build residue IDs list
        residue_ids = [f"A{idx}" for idx in range(1, len(sorted_keys) + 1)]

        # Center scaffold on metal position (metal at origin)
        metal_center = None
        for line in metal_lines:
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                metal_center = (x, y, z)
                break
            except ValueError:
                pass

        def center_line(line: str) -> str:
            if not metal_center:
                return line
            try:
                x = float(line[30:38]) - metal_center[0]
                y = float(line[38:46]) - metal_center[1]
                z = float(line[46:54]) - metal_center[2]
                return line[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + line[54:]
            except ValueError:
                return line

        # Build output PDB
        output_lines = []
        atom_serial = 1

        def renumber_atom_serial(line: str, serial: int) -> str:
            return f"{line[:6]}{serial:5d}{line[11:]}"

        # Protein atoms -> chain A, renumbered
        for line, chain, resnum in atom_lines:
            new_resnum = old_to_new_resnum.get((chain, resnum), resnum)
            new_line = line[:21] + "A" + f"{new_resnum:4d}" + line[26:]
            new_line = renumber_atom_serial(new_line, atom_serial)
            atom_serial += 1
            new_line = center_line(new_line)
            output_lines.append(new_line)

        # TER after protein
        if output_lines:
            output_lines.append("TER")

        # Metal atoms — rewrite residue name if target_metal substitution is active
        for line in metal_lines:
            new_line = line
            if target_metal and target_metal.upper() != metal_upper:
                # Rewrite the residue name (cols 17-20) and atom name (cols 12-16)
                # and element (cols 76-78) to the target metal
                tgt = target_metal.upper()
                # Residue name: right-justified in cols 17:20
                new_line = new_line[:17] + f"{tgt:>3s}" + new_line[20:]
                # Atom name: cols 12:16 — metal atoms use element as name
                atom_name_field = f" {tgt:<3s}" if len(tgt) <= 2 else f"{tgt:<4s}"
                new_line = new_line[:12] + atom_name_field + new_line[16:]
                # Element symbol: cols 76:78 (right-justified)
                if len(new_line) >= 78:
                    new_line = new_line[:76] + f"{tgt:>2s}" + new_line[78:]
                logger.info(f"Theozyme metal rewritten: {metal_upper} → {tgt}")
            new_line = renumber_atom_serial(new_line, atom_serial)
            atom_serial += 1
            new_line = center_line(new_line)
            output_lines.append(new_line)

        # Ligand atoms (keep original residue name, e.g., "PQQ")
        for line in ligand_lines:
            # Keep the ligand residue name - RFD3 recognizes it
            new_line = renumber_atom_serial(line, atom_serial)
            atom_serial += 1
            new_line = center_line(new_line)
            output_lines.append(new_line)

        if metal_center:
            logger.info(f"Theozyme centered on metal at ({metal_center[0]:.1f}, {metal_center[1]:.1f}, {metal_center[2]:.1f})")

        output_lines.append("END")
        return '\n'.join(output_lines), residue_ids

    def _extract_motif_pdb(
        self,
        pdb_content: str,
        active_site: Dict,
        coord_residues: List[Dict],
        metal: Optional[str],
        ligand_code: Optional[str],
        all_residue_keys: Optional[Set[Tuple[str, int]]] = None,
        create_blockers: bool = True,
    ) -> Tuple[str, List[str], List[str]]:
        """
        Extract motif atoms and renumber for RFD3 compatibility.

        IMPORTANT: Renumbers to avoid RFD3 chain conflicts:
        - Chain A for protein residues (renumbered sequentially from 1)
        - Chain M for metal (residue 1)
        - Chain L for ligand (residue 1)

        CRITICAL: Also creates "blocker" ALA residues at ligand positions.
        RFD3 avoids protein ATOM records but ignores HETATM - blockers prevent
        backbone from passing through the ligand.

        Args:
            all_residue_keys: Set of (chain, resnum) tuples for ALL residues
                             that must be included (from _residues_to_ranges)
            create_blockers: If True, create ALA blockers at ligand positions

        Returns:
            Tuple of (motif_pdb_content, fixed_residues_list, blocker_residues_list)
        """
        # Use all_residue_keys if provided, otherwise fall back to coord_residues
        if all_residue_keys is not None:
            residue_keys_to_keep = all_residue_keys
        else:
            # Fallback: build set from coord_residues
            residue_keys_to_keep = set()
            for res in coord_residues:
                key = (res.get("chain", "A"), res.get("resnum", 0))
                residue_keys_to_keep.add(key)

        metal_upper = metal.upper() if metal else None
        ligand_upper = ligand_code.upper() if ligand_code else None

        # Collect ATOM and HETATM lines
        atom_lines = []
        metal_lines = []
        ligand_lines = []

        for line in pdb_content.split('\n'):
            if line.startswith('HETATM'):
                res_name = line[17:20].strip()

                # Keep metal
                if metal_upper and res_name == metal_upper:
                    metal_lines.append(line)
                    continue

                # Keep ligand
                if ligand_upper and res_name == ligand_upper:
                    ligand_lines.append(line)
                    continue

            elif line.startswith('ATOM'):
                chain = line[21]
                try:
                    resnum = int(line[22:26].strip())
                except ValueError:
                    continue

                # Keep residues in ranges
                if (chain, resnum) in residue_keys_to_keep:
                    atom_lines.append((line, chain, resnum))

        # Build old->new residue mapping for protein atoms
        # Renumber sequentially starting from 1
        old_to_new_resnum = {}
        sorted_keys = sorted(residue_keys_to_keep)
        for idx, (old_chain, old_resnum) in enumerate(sorted_keys, start=1):
            old_to_new_resnum[(old_chain, old_resnum)] = idx

        # Build fixed_residues list (new A-chain residue identifiers)
        fixed_residues = [f"A{idx}" for idx in range(1, len(sorted_keys) + 1)]

        # =====================================================================
        # CENTER SCAFFOLD ON METAL POSITION
        # =====================================================================
        # Critical: We must center the scaffold so the metal is at origin (0,0,0).
        # This way RFD3 designs the protein AROUND the metal-ligand complex.
        # Without centering, RFD3 generates around origin but HETATM is at
        # crystal coordinates -> designs won't coordinate the metal.

        metal_center = None
        for line in metal_lines:
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                metal_center = (x, y, z)
                break
            except ValueError:
                pass

        # Helper to center coordinates
        def center_line(line: str) -> str:
            if not metal_center:
                return line
            try:
                x = float(line[30:38]) - metal_center[0]
                y = float(line[38:46]) - metal_center[1]
                z = float(line[46:54]) - metal_center[2]
                return line[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + line[54:]
            except ValueError:
                return line

        # Build renumbered and centered PDB
        output_lines = []
        atom_serial = 1  # Sequential atom serial numbering

        # Helper to renumber atom serial
        def renumber_atom_serial(line: str, serial: int) -> str:
            return f"{line[:6]}{serial:5d}{line[11:]}"

        # Renumber protein atoms to chain A with sequential residue AND atom numbering
        for line, chain, resnum in atom_lines:
            new_resnum = old_to_new_resnum.get((chain, resnum), resnum)
            # Rebuild line with new chain (A), resnum, and atom serial
            new_line = line[:21] + "A" + f"{new_resnum:4d}" + line[26:]
            new_line = renumber_atom_serial(new_line, atom_serial)
            atom_serial += 1
            # Center on metal
            new_line = center_line(new_line)
            output_lines.append(new_line)

        # Renumber metal to chain M, resnum 1 (at origin after centering)
        for line in metal_lines:
            new_line = line[:21] + "M" + "   1" + line[26:]
            new_line = renumber_atom_serial(new_line, atom_serial)
            atom_serial += 1
            new_line = center_line(new_line)
            output_lines.append(new_line)

        # Renumber ligand to chain L, resnum 1 (centered relative to metal)
        # NOTE: Don't add ligand to output_lines yet - we need to insert blockers first
        centered_ligand_lines = []
        for line in ligand_lines:
            new_line = line[:21] + "L" + "   1" + line[26:]
            new_line = center_line(new_line)
            centered_ligand_lines.append(new_line)

        # =====================================================================
        # CREATE BLOCKER RESIDUES AT LIGAND POSITIONS
        # =====================================================================
        # Critical: RFD3 avoids ATOM records but ignores HETATM.
        # Create ALA "blocker" residues at ligand atom positions to prevent
        # the designed backbone from passing through the ligand.
        # IMPORTANT: Blockers must come BEFORE HETATM and have proper atom serials
        blocker_residues = []
        if create_blockers and ligand_code and centered_ligand_lines:
            # Blockers continue chain A numbering after scaffold residues
            blocker_start = len(fixed_residues) + 1
            blocker_pdb, blocker_residues = self._create_blocker_residues(
                ligand_lines=centered_ligand_lines,
                ligand_code=ligand_code,
                start_resnum=blocker_start,
                min_spacing=2.5,
                start_atom_serial=atom_serial,
            )
            if blocker_pdb:
                # Add blockers AFTER protein ATOM records (before HETATM)
                output_lines.append(blocker_pdb)
                # Update atom_serial counter for ligand HETATM
                atom_serial += len(blocker_residues) * 5  # 5 atoms per ALA blocker
                logger.info(f"Added {len(blocker_residues)} blocker residues: {blocker_residues[:3]}...")

        # Now add ligand HETATM (after blockers) with proper atom serials
        for line in centered_ligand_lines:
            new_line = renumber_atom_serial(line, atom_serial)
            atom_serial += 1
            output_lines.append(new_line)

        if metal_center:
            logger.info(f"Scaffold centered on metal at ({metal_center[0]:.1f}, {metal_center[1]:.1f}, {metal_center[2]:.1f})")

        output_lines.append("END")
        return '\n'.join(output_lines), fixed_residues, blocker_residues

    def _create_blocker_residues(
        self,
        ligand_lines: List[str],
        ligand_code: str,
        start_resnum: int = 1,
        min_spacing: float = 2.5,
        start_atom_serial: int = 1,
    ) -> Tuple[str, List[str]]:
        """
        Create "blocker" protein residues at ligand atom positions.

        CRITICAL: RFD3 avoids protein ATOM records but IGNORES ligand HETATM records.
        By placing ALA residues at ligand positions, we force RFD3 to design around them.

        Args:
            ligand_lines: HETATM lines for the ligand (already centered)
            ligand_code: 3-letter ligand code (e.g., "PQQ")
            start_resnum: Starting residue number for blockers
            min_spacing: Minimum spacing between blockers (Å)
            start_atom_serial: Starting atom serial number (for PDB ordering)

        Returns:
            Tuple of (blocker PDB lines, list of blocker residue IDs like ["A16", "A17", ...])
        """
        import math

        ligand_upper = ligand_code.upper()

        # Extract ligand atom positions
        ligand_atoms = []
        for line in ligand_lines:
            if line.startswith('HETATM') and line[17:20].strip() == ligand_upper:
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    atom_name = line[12:16].strip()
                    ligand_atoms.append({'atom': atom_name, 'pos': (x, y, z)})
                except (ValueError, IndexError):
                    pass

        if not ligand_atoms:
            logger.warning(f"No atoms found for ligand {ligand_code}")
            return "", []

        logger.info(f"Creating blockers from {len(ligand_atoms)} ligand atoms")

        # Cluster atoms that are too close together
        clusters = []
        used = set()

        for i, atom in enumerate(ligand_atoms):
            if i in used:
                continue

            cluster = [atom['pos']]
            used.add(i)

            for j, other in enumerate(ligand_atoms):
                if j in used:
                    continue
                dist = math.sqrt(
                    (atom['pos'][0] - other['pos'][0])**2 +
                    (atom['pos'][1] - other['pos'][1])**2 +
                    (atom['pos'][2] - other['pos'][2])**2
                )
                if dist < min_spacing:
                    cluster.append(other['pos'])
                    used.add(j)

            # Compute cluster centroid
            centroid = (
                sum(p[0] for p in cluster) / len(cluster),
                sum(p[1] for p in cluster) / len(cluster),
                sum(p[2] for p in cluster) / len(cluster),
            )
            clusters.append(centroid)

        logger.info(f"Created {len(clusters)} blocker positions")

        # Create ALA residues at each cluster centroid
        blocker_lines = []
        blocker_residues = []
        atom_serial = start_atom_serial  # Continue from protein atoms

        for res_num, centroid in enumerate(clusters, start=start_resnum):
            x, y, z = centroid

            # ALA backbone + CB
            atoms = [
                ('N',   x - 1.5, y, z),
                ('CA',  x, y, z),
                ('C',   x + 1.5, y, z),
                ('O',   x + 1.5, y + 1.2, z),
                ('CB',  x, y - 1.5, z),
            ]

            for atom_name, ax, ay, az in atoms:
                if len(atom_name) <= 2:
                    atom_field = f" {atom_name:<3s}"
                else:
                    atom_field = f"{atom_name:<4s}"

                line = f"ATOM  {atom_serial:5d} {atom_field:4s} ALA A{res_num:4d}    {ax:8.3f}{ay:8.3f}{az:8.3f}  1.00  0.00          {atom_name[0]:>2s}"
                blocker_lines.append(line)
                atom_serial += 1

            blocker_residues.append(f"A{res_num}")

        logger.info(f"Generated {len(blocker_residues)} blocker ALA residues: {blocker_residues[:5]}...")

        return '\n'.join(blocker_lines), blocker_residues

    def _build_renumbered_contig(
        self,
        scaffold_fixed_residues: List[str],
        chain_length: str,
        design_mode: str = "dimer",
    ) -> str:
        """
        Build contig for renumbered scaffold (A1-An, M1, L1).

        The scaffold now has sequential residue numbering on chain A.

        Args:
            scaffold_fixed_residues: List of residue IDs on chain A (including blockers)
            chain_length: Desired length for new chains (e.g., "60-80")
            design_mode: "dimer" (separate chain B) or "monomer" (extensions on chain A)

        design_mode="dimer" (recommended for ligand scaffolding):
            Creates chain B SEPARATE from scaffold chain A
            Contig: "A1-{n},/0,{chain_length}"
            The /0 creates a chain break, chain B is positioned via ori_token

        design_mode="monomer":
            Extensions on same chain as scaffold
            Contig: "15-30,A1-{n},15-30"
            Risk: new residues adjacent to scaffold may clash with ligand
        """
        if not scaffold_fixed_residues:
            return chain_length

        # All scaffold residues are now on chain A, sequentially numbered
        n = len(scaffold_fixed_residues)

        if design_mode == "dimer":
            # DIMER DESIGN: Scaffold (A) + new chain (B)
            # Chain B is positioned away from ligand via ori_token
            # This avoids linker regions passing through ligand space
            contig = f"A1-{n},/0,{chain_length}"
            logger.info(f"Using DIMER contig: {contig}")
        else:
            # MONOMER DESIGN: Extensions on chain A
            # Risk of clashes if linkers pass through ligand
            contig = f"15-30,A1-{n},15-30"
            logger.info(f"Using MONOMER contig: {contig}")

        return contig

    def _determine_fixed_atoms_renumbered(
        self,
        metal: Optional[str],
        ligand_code: Optional[str],
        scaffold_fixed_residues: List[str],
    ) -> Dict[str, str]:
        """
        Determine which atoms should be fixed during RFD3.

        After renumbering:
        - Metal is on chain M, residue 1 -> "M1"
        - Ligand is on chain L, residue 1 -> "L1"
        - Protein scaffold residues are A1, A2, ... -> each one fixed
        """
        fixed = {}

        # Fix metal
        if metal:
            fixed["M1"] = "ALL"

        # Fix ligand
        if ligand_code:
            fixed["L1"] = "ALL"

        # Fix scaffold protein residues
        for res_id in scaffold_fixed_residues:
            fixed[res_id] = "ALL"

        return fixed

    def _determine_rasa_targets_renumbered(
        self,
        metal: Optional[str],
        ligand_code: Optional[str] = None,
    ) -> Dict[str, str]:
        """
        Determine burial (RASA) targets for RFD3.

        Buries both metal and ligand to create a proper binding pocket.
        Uses residue NAME (e.g., "CA", "PQQ") not chain+resnum format.
        """
        targets = {}

        if metal:
            targets[metal.upper()] = "ALL"  # Bury metal

        if ligand_code:
            targets[ligand_code.upper()] = "ALL"  # Bury ligand

        return targets

    def _determine_hbond_conditioning(
        self,
        pdb_content: str,
        ligand_code: Optional[str],
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Determine H-bond conditioning from ligand chemistry.

        Uses ligand NAME (e.g., "PQQ") as key, not chain format "L1".

        Returns:
            (hbond_acceptors, hbond_donors) dicts
        """
        acceptors = {}
        donors = {}

        if not ligand_code:
            return acceptors, donors

        ligand_upper = ligand_code.upper()

        # Find oxygen atoms (acceptors) and nitrogen atoms (potential donors)
        o_atoms = []
        n_atoms = []

        for line in pdb_content.split('\n'):
            if line.startswith('HETATM') and line[17:20].strip() == ligand_upper:
                atom_name = line[12:16].strip()
                element = line[76:78].strip() if len(line) > 76 else atom_name[0]

                if element == 'O':
                    o_atoms.append(atom_name)
                elif element == 'N':
                    n_atoms.append(atom_name)

        # Use ligand NAME as key (not "L1" chain format)
        if o_atoms:
            acceptors[ligand_upper] = ",".join(o_atoms[:6])  # Limit to 6 atoms

        if n_atoms:
            donors[ligand_upper] = ",".join(n_atoms[:4])  # Limit to 4 atoms

        return acceptors, donors


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def get_available_reference_structures() -> Dict[str, List[str]]:
    """Get available reference PDB structures for each metal."""
    return REFERENCE_STRUCTURES.copy()


async def quick_scaffold_check(pdb_id: str) -> Dict[str, Any]:
    """
    Quick check if a PDB is suitable for scaffolding.

    Returns dict with:
    - valid: bool
    - metals_found: list
    - ligands_found: list
    - message: str
    """
    try:
        pdb_content = await asyncio.to_thread(_fetch_pdb_content, pdb_id)
        if not pdb_content:
            return {
                "valid": False,
                "metals_found": [],
                "ligands_found": [],
                "message": f"Could not fetch PDB: {pdb_id}",
            }

        metals = set()
        ligands = set()

        common_metals = {"ZN", "FE", "MG", "CA", "MN", "CU", "CO", "NI", "TB", "EU", "GD"}

        for line in pdb_content.split('\n'):
            if line.startswith('HETATM'):
                res_name = line[17:20].strip()
                if res_name in common_metals:
                    metals.add(res_name)
                elif len(res_name) == 3 and res_name not in {"HOH", "WAT"}:
                    ligands.add(res_name)

        return {
            "valid": len(metals) > 0 or len(ligands) > 0,
            "metals_found": list(metals),
            "ligands_found": list(ligands),
            "message": f"Found {len(metals)} metals, {len(ligands)} ligands",
        }

    except Exception as e:
        return {
            "valid": False,
            "metals_found": [],
            "ligands_found": [],
            "message": str(e),
        }
