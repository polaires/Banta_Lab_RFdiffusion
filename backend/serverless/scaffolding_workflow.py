"""
Scaffolding Workflow for AI Design Pipeline.

Wraps existing handler.py scaffolding functions for use with AI-driven design.
Enables queries like "scaffold the PQQ-Ca pocket of 4CVB".

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

logger = logging.getLogger(__name__)


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class ScaffoldResult:
    """Result from scaffolding workflow."""
    # Core outputs
    motif_pdb: str                         # Extracted motif ready for RFD3
    motif_residues: List[str] = field(default_factory=list)  # Residues to preserve (e.g., ["A10-15"])
    contig: str = ""                       # Full contig string for RFD3

    # Conditioning parameters for RFD3
    fixed_atoms: Dict[str, str] = field(default_factory=dict)      # {"X1": "all", "L1": "all"}
    rasa_targets: Dict[str, str] = field(default_factory=dict)     # Burial targets
    hbond_acceptors: Dict[str, str] = field(default_factory=dict)  # H-bond acceptor atoms
    hbond_donors: Dict[str, str] = field(default_factory=dict)     # H-bond donor atoms

    # Metadata
    source_info: Dict[str, Any] = field(default_factory=dict)      # Source structure info
    coordinating_residues: List[Dict] = field(default_factory=list)  # Coordination shell

    # Status
    success: bool = True
    error_message: str = ""


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
        chain_length: str = "80-120",
        coordination_cutoff: Optional[float] = None,
        include_second_shell: bool = False,
        include_all_ligand_contacts: bool = False,
        ligand_contact_cutoff: Optional[float] = None,
    ) -> ScaffoldResult:
        """
        Run scaffolding workflow.

        Steps:
        1. Fetch PDB structure from RCSB
        2. Find active site (metal + ligand + coordinating residues)
        3. Extract motif with fixed atoms
        4. Build scaffolding contig
        5. Determine conditioning parameters

        Args:
            pdb_id: 4-character PDB ID (e.g., "4CVB")
            metal: Metal element symbol (e.g., "CA", "ZN")
            ligand_code: Ligand 3-letter code (e.g., "PQQ", "CIT")
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

        logger.info(f"ScaffoldingWorkflow: PDB={pdb_id}, Metal={metal}, Ligand={ligand_code}, all_contacts={include_all_ligand_contacts}")

        try:
            # Step 1: Fetch structure from RCSB
            pdb_content = await asyncio.to_thread(_fetch_pdb_content, pdb_id)
            if not pdb_content:
                return ScaffoldResult(
                    motif_pdb="",
                    success=False,
                    error_message=f"Could not fetch PDB: {pdb_id}",
                )

            # Step 2: Find active site
            if metal and ligand_code:
                # Use metal-ligand active site finder
                active_site = await asyncio.to_thread(
                    find_metal_ligand_active_site,
                    pdb_id, metal, ligand_code, cutoff
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
                else:
                    active_site = None
            else:
                # Try to auto-detect metal from PDB
                active_site = self._auto_detect_active_site(pdb_content, cutoff)

            if not active_site:
                return ScaffoldResult(
                    motif_pdb="",
                    success=False,
                    error_message=f"No active site found in {pdb_id}",
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

            # Convert to residue ranges for contig
            motif_residues = self._residues_to_ranges(unique_residues)

            # Step 4: Build motif PDB (active site only)
            motif_pdb = self._extract_motif_pdb(
                pdb_content,
                active_site,
                unique_residues,
                metal,
                ligand_code,
            )

            # Step 5: Build scaffolding contig
            contig = self._build_scaffolding_contig(
                motif_residues=motif_residues,
                chain_length=chain_length,
            )

            # Step 6: Determine conditioning parameters
            fixed_atoms = self._determine_fixed_atoms(metal, ligand_code)
            rasa_targets = self._determine_rasa_targets(metal, ligand_code)
            hbond_acceptors, hbond_donors = self._determine_hbond_conditioning(
                pdb_content, ligand_code
            )

            return ScaffoldResult(
                motif_pdb=motif_pdb,
                motif_residues=motif_residues,
                contig=contig,
                fixed_atoms=fixed_atoms,
                rasa_targets=rasa_targets,
                hbond_acceptors=hbond_acceptors,
                hbond_donors=hbond_donors,
                coordinating_residues=unique_residues,
                source_info={
                    "pdb_id": pdb_id.upper(),
                    "metal": metal,
                    "ligand": ligand_code,
                    "metal_chain": active_site.get("metal_chain"),
                    "metal_resnum": active_site.get("metal_resnum"),
                    "coordination_number": active_site.get("coordination_number", 0),
                    "protein_donors": active_site.get("protein_donors", []),
                    "ligand_donors": active_site.get("ligand_donors", []),
                    "include_all_contacts": include_all_ligand_contacts,
                    "total_pocket_residues": len(unique_residues),
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
        gap_threshold: int = 10,
    ) -> List[str]:
        """
        Convert residue list to range strings like 'A10-15'.

        Groups residues that are within gap_threshold of each other.
        This prevents overly fragmented contigs.

        Args:
            residues: List of residue dicts
            gap_threshold: Max gap to merge into single range (default 10)

        Returns:
            List of range strings
        """
        if not residues:
            return []

        # Group by chain
        chains: Dict[str, List[int]] = {}
        for res in residues:
            chain = res.get("chain", "A")
            resnum = res.get("resnum", 0)
            if chain not in chains:
                chains[chain] = []
            chains[chain].append(resnum)

        ranges = []
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
                    # Start new range
                    ranges.append(f"{chain}{start}-{end}")
                    start = end = resnums[i]

            ranges.append(f"{chain}{start}-{end}")

        return ranges

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
            return f"0 {chain_length}"

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
            return f"0 {chain_length}"

        # Sort by start position
        protein_ranges.sort(key=lambda x: x[1])

        parts = []

        # N-terminal extension (short)
        parts.append("0 15-30")

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
                parts.append(f"0 {linker_min}-{linker_max}")

        # C-terminal extension (short)
        parts.append("0 15-30")

        # Add ligand/metal references at the end if present
        ligand_parts = []
        for r in motif_residues:
            match = re.match(r'([A-Z])(\d+)-(\d+)', r)
            if match:
                start = int(match.group(2))
                if start >= 1000:
                    # This is ligand/metal, add as fixed
                    ligand_parts.append(r)

        # Build final contig
        contig = "/".join(parts)
        if ligand_parts:
            contig += "/" + "/".join(ligand_parts)

        return contig

    def _extract_motif_pdb(
        self,
        pdb_content: str,
        active_site: Dict,
        coord_residues: List[Dict],
        metal: Optional[str],
        ligand_code: Optional[str],
    ) -> str:
        """
        Extract just the motif atoms from full PDB.

        Includes:
        - Metal ion (if present)
        - Ligand (if present)
        - Coordinating residues
        """
        lines = []

        # Build set of residues to keep
        coord_resnums = set()
        for res in coord_residues:
            key = (res.get("chain", "A"), res.get("resnum", 0))
            coord_resnums.add(key)

        metal_upper = metal.upper() if metal else None
        ligand_upper = ligand_code.upper() if ligand_code else None

        for line in pdb_content.split('\n'):
            if line.startswith('HETATM'):
                res_name = line[17:20].strip()

                # Keep metal
                if metal_upper and res_name == metal_upper:
                    lines.append(line)
                    continue

                # Keep ligand
                if ligand_upper and res_name == ligand_upper:
                    lines.append(line)
                    continue

            elif line.startswith('ATOM'):
                chain = line[21]
                try:
                    resnum = int(line[22:26].strip())
                except ValueError:
                    continue

                # Keep coordinating residues
                if (chain, resnum) in coord_resnums:
                    lines.append(line)

        lines.append("END")
        return '\n'.join(lines)

    def _determine_fixed_atoms(
        self,
        metal: Optional[str],
        ligand_code: Optional[str],
    ) -> Dict[str, str]:
        """
        Determine which atoms should be fixed during RFD3.
        """
        fixed = {}

        if metal:
            fixed["X1"] = "all"  # Metal always fixed

        if ligand_code:
            fixed["L1"] = "all"  # Ligand always fixed

        return fixed

    def _determine_rasa_targets(
        self,
        metal: Optional[str],
        ligand_code: Optional[str],
    ) -> Dict[str, str]:
        """
        Determine burial (RASA) targets for RFD3.
        """
        targets = {}

        if metal:
            targets["X1"] = "all"  # Bury metal

        # Optionally bury ligand hydrophobic portions
        # For now, default to burying metal only

        return targets

    def _determine_hbond_conditioning(
        self,
        pdb_content: str,
        ligand_code: Optional[str],
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Determine H-bond conditioning from ligand chemistry.

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

        if o_atoms:
            acceptors["L1"] = ",".join(o_atoms[:6])  # Limit to 6 atoms

        if n_atoms:
            donors["L1"] = ",".join(n_atoms[:4])  # Limit to 4 atoms

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
