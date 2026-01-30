# metal_site_fetcher.py
"""
PDB Metal Site Fetcher Module

Fetches metal binding site information from the RCSB PDB database.
Prioritizes experimental structures over calculated templates for accuracy.

This module provides:
- Query RCSB PDB for metal binding structures
- Extract coordination geometry from PDB content
- Generate templates from experimental structures
- Curated reference structures for common metals
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any, Tuple
import math
import logging

# Graceful import of requests - may not be available in all environments
try:
    import requests
    REQUESTS_AVAILABLE = True
except ImportError:
    REQUESTS_AVAILABLE = False

logger = logging.getLogger(__name__)

# =============================================================================
# CONSTANTS
# =============================================================================

RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_DATA_URL = "https://data.rcsb.org/rest/v1/core/entry"
RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download"

# Curated high-quality reference PDB structures for each metal
# These are well-characterized, high-resolution structures
REFERENCE_STRUCTURES: Dict[str, List[str]] = {
    # Zinc - catalytic and structural sites
    "ZN": ["1CA2", "1HET", "1A5T", "2JEZ", "1DSZ"],

    # Iron - heme and non-heme
    "FE": ["1RB9", "1MBS", "1WOC", "1HHO", "2CDV"],

    # Copper - type I, II, and binuclear
    "CU": ["1PLC", "1NWP", "2CUA", "1OAC", "1KDI"],

    # Calcium - EF-hand and non-EF-hand
    "CA": ["1CLL", "1W6S", "1PON", "3CLN", "1EXR"],

    # Magnesium - catalytic sites
    "MG": ["3PGK", "1ATP", "1PHK", "4ENL", "1CSA"],

    # Manganese - superoxide dismutase and others
    "MN": ["3MDS", "1N0J", "1VJ1", "1YAV", "1QNM"],

    # Cobalt - B12 and others
    "CO": ["1CBN", "3BCF", "1DUO", "2B2P", "1MJ4"],

    # Nickel - urease and others
    "NI": ["1FWJ", "1EJX", "2UBP", "3QGA", "1KVF"],

    # Lanthanides - lanmodulin and engineered sites
    "TB": ["6MI5", "7D4I", "5VNE", "6MQG", "6S0V"],
    "EU": ["6MI5", "7D4I", "5VNE", "6MQG", "6S0V"],
    "GD": ["6MI5", "7D4I", "5VNE", "6MQG", "6S0V"],
    "LA": ["6MI5", "7D4I", "5VNE", "6MQG", "6S0V"],
}


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class MetalSite:
    """Represents a metal binding site extracted from PDB."""
    pdb_id: str
    metal: str
    metal_chain: str = ""
    metal_resnum: int = 0
    metal_coords: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    coordination_number: int = 0
    coordinating_atoms: List[Dict[str, Any]] = field(default_factory=list)
    resolution: float = 0.0
    geometry: str = ""

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        return {
            "pdb_id": self.pdb_id,
            "metal": self.metal,
            "metal_chain": self.metal_chain,
            "metal_resnum": self.metal_resnum,
            "metal_coords": self.metal_coords,
            "coordination_number": self.coordination_number,
            "coordinating_atoms": self.coordinating_atoms,
            "resolution": self.resolution,
            "geometry": self.geometry,
        }


@dataclass
class MetalLigandSite(MetalSite):
    """Represents a metal-ligand complex site from PDB."""
    ligand_code: str = ""
    ligand_chain: str = ""
    ligand_resnum: int = 0
    ligand_atoms: List[Dict[str, Any]] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        base = super().to_dict()
        base.update({
            "ligand_code": self.ligand_code,
            "ligand_chain": self.ligand_chain,
            "ligand_resnum": self.ligand_resnum,
            "ligand_atoms": self.ligand_atoms,
        })
        return base


# =============================================================================
# RCSB API QUERY FUNCTIONS
# =============================================================================

def _build_metal_search_query(
    metal: str,
    resolution_max: float = 2.5,
    limit: int = 10,
) -> Dict[str, Any]:
    """
    Build RCSB search query for metal-containing structures.

    Args:
        metal: Metal element symbol (e.g., "ZN", "FE", "TB")
        resolution_max: Maximum resolution in Angstroms
        limit: Maximum number of results

    Returns:
        Query dict for RCSB API
    """
    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
                        "operator": "exact_match",
                        "value": metal.upper(),
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.resolution_combined",
                        "operator": "less_or_equal",
                        "value": resolution_max,
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "results_content_type": ["experimental"],
            "sort": [
                {
                    "sort_by": "rcsb_entry_info.resolution_combined",
                    "direction": "asc",
                }
            ],
            "paginate": {
                "start": 0,
                "rows": limit,
            },
        },
    }


def _build_metal_ligand_search_query(
    metal: str,
    ligand: str,
    resolution_max: float = 3.0,
    limit: int = 10,
) -> Dict[str, Any]:
    """
    Build RCSB search query for metal-ligand complex structures.

    Args:
        metal: Metal element symbol
        ligand: Ligand 3-letter code (e.g., "PQQ", "HEM")
        resolution_max: Maximum resolution
        limit: Maximum results

    Returns:
        Query dict for RCSB API
    """
    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
                        "operator": "exact_match",
                        "value": metal.upper(),
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
                        "operator": "exact_match",
                        "value": ligand.upper(),
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.resolution_combined",
                        "operator": "less_or_equal",
                        "value": resolution_max,
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "results_content_type": ["experimental"],
            "sort": [
                {
                    "sort_by": "rcsb_entry_info.resolution_combined",
                    "direction": "asc",
                }
            ],
            "paginate": {
                "start": 0,
                "rows": limit,
            },
        },
    }


def query_metal_sites(
    metal: str,
    resolution_max: float = 2.5,
    limit: int = 10,
) -> List[Dict[str, Any]]:
    """
    Query RCSB PDB for structures containing a specific metal.

    Args:
        metal: Metal element symbol (e.g., "ZN", "FE", "TB")
        resolution_max: Maximum resolution in Angstroms
        limit: Maximum number of results to return

    Returns:
        List of dicts with keys: pdb_id, metal, resolution
        Empty list if network unavailable or no results
    """
    if not REQUESTS_AVAILABLE:
        logger.warning("requests library not available - returning empty results")
        return []

    metal = metal.upper()
    query = _build_metal_search_query(metal, resolution_max, limit)

    try:
        response = requests.post(
            RCSB_SEARCH_URL,
            json=query,
            headers={"Content-Type": "application/json"},
            timeout=30,
        )
        response.raise_for_status()

        data = response.json()
        results = []

        for hit in data.get("result_set", []):
            pdb_id = hit.get("identifier", "")
            results.append({
                "pdb_id": pdb_id,
                "metal": metal,
                "resolution": hit.get("score", 0.0),  # Score is often resolution
            })

        return results

    except requests.RequestException as e:
        logger.warning(f"RCSB API query failed: {e}")
        return []
    except (KeyError, ValueError) as e:
        logger.warning(f"Error parsing RCSB response: {e}")
        return []


def query_metal_ligand_sites(
    metal: str,
    ligand: str,
    resolution_max: float = 3.0,
    limit: int = 10,
) -> List[Dict[str, Any]]:
    """
    Query RCSB PDB for structures containing a metal-ligand complex.

    Args:
        metal: Metal element symbol (e.g., "CA", "MG")
        ligand: Ligand 3-letter code (e.g., "PQQ", "ATP")
        resolution_max: Maximum resolution in Angstroms
        limit: Maximum number of results to return

    Returns:
        List of dicts with keys: pdb_id, metal, ligand, resolution
        Empty list if network unavailable or no results
    """
    if not REQUESTS_AVAILABLE:
        logger.warning("requests library not available - returning empty results")
        return []

    metal = metal.upper()
    ligand = ligand.upper()
    query = _build_metal_ligand_search_query(metal, ligand, resolution_max, limit)

    try:
        response = requests.post(
            RCSB_SEARCH_URL,
            json=query,
            headers={"Content-Type": "application/json"},
            timeout=30,
        )
        response.raise_for_status()

        data = response.json()
        results = []

        for hit in data.get("result_set", []):
            pdb_id = hit.get("identifier", "")
            results.append({
                "pdb_id": pdb_id,
                "metal": metal,
                "ligand": ligand,
                "resolution": hit.get("score", 0.0),
            })

        return results

    except requests.RequestException as e:
        logger.warning(f"RCSB API query failed: {e}")
        return []
    except (KeyError, ValueError) as e:
        logger.warning(f"Error parsing RCSB response: {e}")
        return []


# =============================================================================
# PDB CONTENT PARSING
# =============================================================================

def _parse_pdb_atom_line(line: str) -> Optional[Dict[str, Any]]:
    """
    Parse a single ATOM or HETATM line from PDB format.

    Args:
        line: Single line from PDB file

    Returns:
        Dict with atom info or None if not a valid atom line
    """
    if not line.startswith(("ATOM", "HETATM")):
        return None

    try:
        return {
            "record_type": line[0:6].strip(),
            "serial": int(line[6:11].strip()) if line[6:11].strip() else 0,
            "name": line[12:16].strip(),
            "alt_loc": line[16:17].strip(),
            "res_name": line[17:20].strip(),
            "chain_id": line[21:22].strip(),
            "res_seq": int(line[22:26].strip()) if line[22:26].strip() else 0,
            "x": float(line[30:38].strip()) if line[30:38].strip() else 0.0,
            "y": float(line[38:46].strip()) if line[38:46].strip() else 0.0,
            "z": float(line[46:54].strip()) if line[46:54].strip() else 0.0,
            "occupancy": float(line[54:60].strip()) if len(line) > 54 and line[54:60].strip() else 1.0,
            "temp_factor": float(line[60:66].strip()) if len(line) > 60 and line[60:66].strip() else 0.0,
            "element": line[76:78].strip() if len(line) > 76 else "",
        }
    except (ValueError, IndexError):
        return None


def _calculate_distance(
    coord1: Tuple[float, float, float],
    coord2: Tuple[float, float, float],
) -> float:
    """Calculate Euclidean distance between two 3D points."""
    return math.sqrt(
        (coord1[0] - coord2[0]) ** 2 +
        (coord1[1] - coord2[1]) ** 2 +
        (coord1[2] - coord2[2]) ** 2
    )


def _infer_geometry(coordination_number: int) -> str:
    """
    Infer coordination geometry from coordination number.

    Args:
        coordination_number: Number of coordinating atoms

    Returns:
        Geometry name string
    """
    geometry_map = {
        2: "linear",
        3: "trigonal_planar",
        4: "tetrahedral",  # or square_planar
        5: "trigonal_bipyramidal",  # or square_pyramidal
        6: "octahedral",
        7: "pentagonal_bipyramidal",
        8: "square_antiprismatic",
        9: "tricapped_trigonal_prismatic",
    }
    return geometry_map.get(coordination_number, "unknown")


def extract_metal_coordination(
    pdb_content: str,
    metal: str,
    cutoff: float = 3.0,
    metal_site_index: int = 0,
) -> Dict[str, Any]:
    """
    Extract coordinating atoms around a metal center from PDB content.

    Args:
        pdb_content: PDB file content as string
        metal: Metal element symbol to search for
        cutoff: Distance cutoff in Angstroms for coordination
        metal_site_index: Which metal site to extract (0 = first)

    Returns:
        Dict with keys:
            - metal: Metal symbol
            - metal_coords: (x, y, z) tuple
            - coordination_number: Number of coordinating atoms
            - coordinating_atoms: List of coordinating atom dicts
            - geometry: Inferred geometry
    """
    metal = metal.upper()
    atoms = []
    metal_atoms = []

    # Parse all atoms
    for line in pdb_content.split("\n"):
        atom = _parse_pdb_atom_line(line)
        if atom:
            atoms.append(atom)
            # Check if this is the target metal
            element = atom.get("element", "").upper()
            res_name = atom.get("res_name", "").upper()
            if element == metal or res_name == metal:
                metal_atoms.append(atom)

    # Return empty result if no metal found
    if not metal_atoms:
        return {
            "metal": metal,
            "metal_coords": (0.0, 0.0, 0.0),
            "coordination_number": 0,
            "coordinating_atoms": [],
            "geometry": "none",
        }

    # Select the target metal site
    if metal_site_index >= len(metal_atoms):
        metal_site_index = 0
    metal_atom = metal_atoms[metal_site_index]
    metal_coords = (metal_atom["x"], metal_atom["y"], metal_atom["z"])

    # Find coordinating atoms within cutoff distance
    # Coordinating atoms should be donor atoms (N, O, S)
    donor_elements = {"N", "O", "S"}
    coordinating = []

    for atom in atoms:
        # Skip the metal itself
        if atom["serial"] == metal_atom["serial"]:
            continue

        # Check if this is a potential donor atom
        element = atom.get("element", "").upper()
        atom_name = atom.get("name", "").upper()

        # Infer element from atom name if element field is empty
        if not element and atom_name:
            element = atom_name[0]

        if element not in donor_elements:
            continue

        # Calculate distance
        atom_coords = (atom["x"], atom["y"], atom["z"])
        distance = _calculate_distance(metal_coords, atom_coords)

        if distance <= cutoff:
            coordinating.append({
                "atom_name": atom["name"],
                "res_name": atom["res_name"],
                "res_seq": atom["res_seq"],
                "chain_id": atom["chain_id"],
                "element": element,
                "distance": round(distance, 2),
                "coords": atom_coords,
            })

    # Sort by distance
    coordinating.sort(key=lambda x: x["distance"])

    coordination_number = len(coordinating)
    geometry = _infer_geometry(coordination_number)

    return {
        "metal": metal,
        "metal_coords": metal_coords,
        "metal_chain": metal_atom.get("chain_id", ""),
        "metal_resnum": metal_atom.get("res_seq", 0),
        "coordination_number": coordination_number,
        "coordinating_atoms": coordinating,
        "geometry": geometry,
    }


# =============================================================================
# PDB FETCHING AND TEMPLATE GENERATION
# =============================================================================

def _fetch_pdb_content(pdb_id: str) -> Optional[str]:
    """
    Fetch PDB content from RCSB.

    Args:
        pdb_id: 4-character PDB ID

    Returns:
        PDB file content as string, or None if unavailable
    """
    if not REQUESTS_AVAILABLE:
        logger.warning("requests library not available - cannot fetch PDB")
        return None

    pdb_id = pdb_id.upper()
    url = f"{RCSB_DOWNLOAD_URL}/{pdb_id}.pdb"

    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        return response.text
    except requests.RequestException as e:
        logger.warning(f"Failed to fetch PDB {pdb_id}: {e}")
        return None


def generate_template_from_pdb(
    pdb_id: str,
    metal: str,
    cutoff: float = 3.0,
    metal_site_index: int = 0,
) -> Optional[Dict[str, Any]]:
    """
    Generate a metal coordination template from a PDB structure.

    Fetches the PDB file and extracts coordination information.

    Args:
        pdb_id: 4-character PDB ID
        metal: Metal element symbol
        cutoff: Distance cutoff for coordination
        metal_site_index: Which metal site to use

    Returns:
        Template dict with coordination info, or None if unavailable
    """
    pdb_content = _fetch_pdb_content(pdb_id)

    if not pdb_content:
        return None

    coord_info = extract_metal_coordination(
        pdb_content,
        metal=metal,
        cutoff=cutoff,
        metal_site_index=metal_site_index,
    )

    if coord_info["coordination_number"] == 0:
        logger.warning(f"No {metal} coordination found in PDB {pdb_id}")
        return None

    # Build template
    template = {
        "pdb_id": pdb_id.upper(),
        "metal": metal.upper(),
        "source": "pdb",
        "metal_coords": coord_info["metal_coords"],
        "coordination_number": coord_info["coordination_number"],
        "geometry": coord_info["geometry"],
        "coordinating_residues": [
            {
                "residue": atom["res_name"],
                "position": atom["res_seq"],
                "chain": atom["chain_id"],
                "atom": atom["atom_name"],
                "distance": atom["distance"],
            }
            for atom in coord_info["coordinating_atoms"]
        ],
    }

    return template


def get_reference_template(
    metal: str,
    ligand: Optional[str] = None,
    cutoff: float = 3.0,
) -> Optional[Dict[str, Any]]:
    """
    Get a template from curated reference structures.

    Tries reference structures in order until one succeeds.

    Args:
        metal: Metal element symbol
        ligand: Optional ligand code (for metal-ligand sites)
        cutoff: Distance cutoff for coordination

    Returns:
        Template dict or None if no reference available
    """
    metal = metal.upper()

    # Get reference PDB IDs for this metal
    pdb_ids = REFERENCE_STRUCTURES.get(metal, [])

    if not pdb_ids:
        logger.warning(f"No reference structures for metal {metal}")
        return None

    # Try each reference structure until one works
    for pdb_id in pdb_ids:
        template = generate_template_from_pdb(pdb_id, metal, cutoff)
        if template:
            return template

    logger.warning(f"Could not generate template from any reference structure for {metal}")
    return None


# =============================================================================
# BATCH OPERATIONS
# =============================================================================

def fetch_all_metal_sites_from_pdb(
    pdb_id: str,
    metals: Optional[List[str]] = None,
    cutoff: float = 3.0,
) -> List[Dict[str, Any]]:
    """
    Fetch all metal sites from a PDB structure.

    Args:
        pdb_id: 4-character PDB ID
        metals: List of metals to search for (None = all common metals)
        cutoff: Distance cutoff for coordination

    Returns:
        List of coordination info dicts for each metal site found
    """
    if metals is None:
        metals = ["ZN", "FE", "CU", "MN", "MG", "CA", "CO", "NI"]

    pdb_content = _fetch_pdb_content(pdb_id)

    if not pdb_content:
        return []

    results = []

    for metal in metals:
        # Extract all sites for this metal
        site_index = 0
        while True:
            coord_info = extract_metal_coordination(
                pdb_content,
                metal=metal,
                cutoff=cutoff,
                metal_site_index=site_index,
            )

            if coord_info["coordination_number"] == 0:
                break

            # Check if we're seeing the same site again (happens when index >= available)
            if site_index > 0:
                # Simple check: if coords are same as first site, we've wrapped
                first_site = extract_metal_coordination(
                    pdb_content, metal=metal, cutoff=cutoff, metal_site_index=0
                )
                if coord_info["metal_coords"] == first_site["metal_coords"]:
                    break

            coord_info["pdb_id"] = pdb_id.upper()
            results.append(coord_info)
            site_index += 1

            # Safety limit
            if site_index > 20:
                break

    return results


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def get_available_reference_metals() -> List[str]:
    """Get list of metals with curated reference structures."""
    return list(REFERENCE_STRUCTURES.keys())


def get_reference_pdb_ids(metal: str) -> List[str]:
    """Get curated reference PDB IDs for a metal."""
    return REFERENCE_STRUCTURES.get(metal.upper(), [])


# =============================================================================
# METAL-LIGAND ACTIVE SITE FUNCTIONS
# =============================================================================

def find_metal_ligand_active_site(
    pdb_id: str,
    metal: str,
    ligand: str,
    cutoff: float = 3.5,
) -> Optional[Dict[str, Any]]:
    """
    Find the metal site that directly coordinates a ligand in a PDB structure.

    This is the key function for finding active site metal coordination where
    the metal binds to an organic cofactor (e.g., PQQ-Ca, ATP-Mg, heme-Fe).

    Args:
        pdb_id: 4-character PDB ID
        metal: Metal element symbol (e.g., "CA", "MG", "FE")
        ligand: Ligand 3-letter code (e.g., "PQQ", "ATP", "HEM")
        cutoff: Distance cutoff for coordination (Angstroms)

    Returns:
        Dict with:
            - metal: Metal symbol
            - metal_chain: Chain ID of metal
            - metal_resnum: Residue number of metal
            - metal_coords: (x, y, z) coordinates
            - ligand_distance: Distance to closest ligand atom
            - coordination_number: Total coordinating atoms
            - coordinating_atoms: List of coordinating atom dicts
            - protein_donors: List of protein residues that coordinate
            - ligand_donors: List of ligand atoms that coordinate
        Returns None if not found or no coordination

    Example:
        >>> site = find_metal_ligand_active_site("1CQ1", "CA", "PQQ")
        >>> print(f"Active site Ca at {site['metal_chain']}{site['metal_resnum']}")
        >>> print(f"PQQ donors: {site['ligand_donors']}")
        >>> print(f"Protein donors: {site['protein_donors']}")
    """
    pdb_content = _fetch_pdb_content(pdb_id)
    if not pdb_content:
        logger.warning(f"Could not fetch PDB {pdb_id}")
        return None

    metal = metal.upper()
    ligand = ligand.upper()

    # Verify structure contains both metal and ligand
    hetatm_types = set()
    for line in pdb_content.split('\n'):
        if line.startswith('HETATM'):
            res_name = line[17:20].strip()
            hetatm_types.add(res_name)

    if metal not in hetatm_types:
        logger.warning(f"No {metal} found in {pdb_id}")
        return None

    if ligand not in hetatm_types:
        logger.warning(f"No {ligand} found in {pdb_id}")
        return None

    # Find all metal positions
    metal_sites = []
    for line in pdb_content.split('\n'):
        if line.startswith('HETATM') and line[17:20].strip() == metal:
            atom_name = line[12:16].strip()
            if atom_name == metal:
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    chain = line[21]
                    resnum = line[22:26].strip()
                    metal_sites.append({
                        'chain': chain, 'resnum': resnum, 'pos': (x, y, z)
                    })
                except (ValueError, IndexError):
                    pass

    if not metal_sites:
        return None

    # Find ligand donor atoms (O, N, S)
    ligand_atoms = []
    for line in pdb_content.split('\n'):
        if line.startswith('HETATM') and line[17:20].strip() == ligand:
            atom_name = line[12:16].strip()
            element = line[76:78].strip() if len(line) > 76 else atom_name[0]
            if element in ('O', 'N', 'S'):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    chain = line[21]
                    resnum = line[22:26].strip()
                    ligand_atoms.append({
                        'atom': atom_name, 'chain': chain,
                        'resnum': resnum, 'pos': (x, y, z)
                    })
                except (ValueError, IndexError):
                    pass

    if not ligand_atoms:
        return None

    # Find metal site closest to ligand
    best_site = None
    best_dist = float('inf')
    for site in metal_sites:
        for lig_atom in ligand_atoms:
            dist = _calculate_distance(site['pos'], lig_atom['pos'])
            if dist < best_dist:
                best_dist = dist
                best_site = site

    if not best_site or best_dist > cutoff:
        logger.warning(f"No {metal} within {cutoff}Å of {ligand}")
        return None

    # Extract full coordination sphere
    metal_pos = best_site['pos']
    coord_atoms = []
    protein_donors = []
    ligand_donors = []

    amino_acids = {
        "GLU", "ASP", "HIS", "CYS", "ASN", "TYR", "MET", "SER", "THR", "GLN",
        "LYS", "ARG", "GLY", "ALA", "VAL", "LEU", "ILE", "PRO", "PHE", "TRP"
    }

    for line in pdb_content.split('\n'):
        if not line.startswith(('ATOM', 'HETATM')):
            continue

        res_name = line[17:20].strip()
        atom_name = line[12:16].strip()
        chain = line[21]
        resnum = line[22:26].strip()

        if res_name == metal and atom_name == metal:
            continue
        if res_name == 'HOH':
            continue

        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            dist = _calculate_distance(metal_pos, (x, y, z))

            if dist <= cutoff:
                element = line[76:78].strip() if len(line) > 76 else atom_name[0]
                atom_info = {
                    'atom_name': atom_name,
                    'res_name': res_name,
                    'res_seq': int(resnum) if resnum.isdigit() else 0,
                    'chain_id': chain,
                    'element': element,
                    'distance': round(dist, 2),
                }
                coord_atoms.append(atom_info)

                if res_name == ligand:
                    ligand_donors.append(f"{atom_name}@{dist:.2f}Å")
                elif res_name in amino_acids:
                    protein_donors.append(f"{chain}{resnum}:{res_name} {atom_name}@{dist:.2f}Å")
        except (ValueError, IndexError):
            pass

    coord_atoms.sort(key=lambda x: x['distance'])

    return {
        'pdb_id': pdb_id.upper(),
        'metal': metal,
        'metal_chain': best_site['chain'],
        'metal_resnum': best_site['resnum'],
        'metal_coords': metal_pos,
        'ligand': ligand,
        'ligand_distance': round(best_dist, 2),
        'coordination_number': len(coord_atoms),
        'coordinating_atoms': coord_atoms,
        'protein_donors': protein_donors,
        'ligand_donors': ligand_donors,
    }


def search_metal_ligand_structures(
    metal: str,
    ligand: str,
    limit: int = 10,
) -> List[str]:
    """
    Search RCSB PDB for structures containing both a metal and a ligand.

    Args:
        metal: Metal element symbol
        ligand: Ligand 3-letter code
        limit: Maximum results

    Returns:
        List of PDB IDs that contain both the metal and ligand
    """
    if not REQUESTS_AVAILABLE:
        logger.warning("requests not available")
        return []

    # Search for ligand first (more specific)
    try:
        query = {
            'query': {
                'type': 'terminal',
                'service': 'text',
                'parameters': {
                    'attribute': 'rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id',
                    'operator': 'exact_match',
                    'value': ligand.upper()
                }
            },
            'return_type': 'entry',
            'request_options': {
                'results_content_type': ['experimental'],
                'paginate': {'start': 0, 'rows': limit * 3}
            }
        }
        response = requests.post(
            RCSB_SEARCH_URL,
            json=query,
            headers={"Content-Type": "application/json"},
            timeout=30,
        )
        response.raise_for_status()

        data = response.json()
        candidates = [h['identifier'] for h in data.get('result_set', [])]

        # Filter to those with the target metal
        results = []
        for pdb_id in candidates[:limit * 2]:
            pdb_content = _fetch_pdb_content(pdb_id)
            if pdb_content:
                for line in pdb_content.split('\n'):
                    if line.startswith('HETATM') and line[17:20].strip() == metal.upper():
                        results.append(pdb_id)
                        break
            if len(results) >= limit:
                break

        return results

    except Exception as e:
        logger.warning(f"Metal-ligand search failed: {e}")
        return []
