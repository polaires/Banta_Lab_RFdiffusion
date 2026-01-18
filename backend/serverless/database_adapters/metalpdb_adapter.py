# metalpdb_adapter.py
"""
MetalPDB Database Adapter

Provides a unified interface for querying metal coordination sites from
the MetalPDB database (https://metalpdb.cerm.unifi.it), with automatic
fallback to RCSB PDB queries when MetalPDB is unavailable.

This adapter is part of the AI-Driven Structure Discovery Infrastructure
for protein design workflows.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any, Tuple
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

METALPDB_BASE_URL = "https://metalpdb.cerm.unifi.it"
METALPDB_API_SEARCH = "/api/search"
METALPDB_API_SITE = "/api/site"

# Valid metal symbols for validation
VALID_METALS = {
    "ZN", "FE", "CU", "MN", "MG", "CA", "CO", "NI", "MO", "W",
    "TB", "EU", "GD", "LA", "CE", "PR", "ND", "SM", "DY", "HO",
    "ER", "TM", "YB", "LU", "Y", "SC", "V", "CR", "CD", "HG",
    "PB", "PT", "AU", "AG", "NA", "K", "LI", "BA", "SR",
}

# Common coordination geometries
GEOMETRY_MAP = {
    2: "linear",
    3: "trigonal_planar",
    4: "tetrahedral",
    5: "trigonal_bipyramidal",
    6: "octahedral",
    7: "pentagonal_bipyramidal",
    8: "square_antiprismatic",
    9: "tricapped_trigonal_prismatic",
}


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class MetalSiteResult:
    """
    Represents a metal binding site result from database queries.

    Attributes:
        pdb_id: 4-character PDB identifier
        metal: Metal element symbol (e.g., "ZN", "FE")
        coordination_number: Number of coordinating atoms
        geometry: Coordination geometry (e.g., "tetrahedral", "octahedral")
        resolution: Structure resolution in Angstroms
        coordinating_residues: List of coordinating residue strings (e.g., ["HIS94", "HIS96"])
        ligands: List of non-protein ligand codes
        source: Data source ("metalpdb" or "rcsb")
        site_index: Index of this site within the structure (for multi-metal)
        metal_coords: (x, y, z) coordinates of the metal center
        chain_id: Chain identifier where metal is located
        distances: Dict mapping residue to distance from metal
    """
    pdb_id: str
    metal: str
    coordination_number: int = 0
    geometry: str = ""
    resolution: float = 0.0
    coordinating_residues: List[str] = field(default_factory=list)
    ligands: List[str] = field(default_factory=list)
    source: str = "metalpdb"
    site_index: int = 0
    metal_coords: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    chain_id: str = ""
    distances: Dict[str, float] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        return {
            "pdb_id": self.pdb_id,
            "metal": self.metal,
            "coordination_number": self.coordination_number,
            "geometry": self.geometry,
            "resolution": self.resolution,
            "coordinating_residues": self.coordinating_residues,
            "ligands": self.ligands,
            "source": self.source,
            "site_index": self.site_index,
            "metal_coords": self.metal_coords,
            "chain_id": self.chain_id,
            "distances": self.distances,
        }


# =============================================================================
# METALPDB ADAPTER CLASS
# =============================================================================

class MetalPDBAdapter:
    """
    Adapter for querying metal coordination sites from MetalPDB database.

    Provides methods for searching metal binding sites by:
    - Metal type and coordination properties
    - Site-specific details including geometry
    - Coordination motif patterns (e.g., "His-His-Glu")

    Includes automatic fallback to RCSB PDB when MetalPDB is unavailable.

    Example:
        >>> adapter = MetalPDBAdapter()
        >>> results = adapter.search_by_metal("ZN", coordination_number=4, limit=10)
        >>> for site in results:
        ...     print(f"{site.pdb_id}: {site.geometry}")
    """

    def __init__(
        self,
        base_url: str = METALPDB_BASE_URL,
        timeout: int = 30,
        enable_fallback: bool = True,
    ):
        """
        Initialize the MetalPDB adapter.

        Args:
            base_url: Base URL for MetalPDB API (default: official MetalPDB)
            timeout: Request timeout in seconds
            enable_fallback: Whether to fallback to RCSB on MetalPDB failure
        """
        self.base_url = base_url
        self.timeout = timeout
        self.enable_fallback = enable_fallback

    def search_by_metal(
        self,
        metal: str,
        coordination_number: Optional[int] = None,
        geometry: Optional[str] = None,
        resolution_max: float = 3.0,
        limit: int = 20,
    ) -> List[MetalSiteResult]:
        """
        Search for metal binding sites by metal type.

        Args:
            metal: Metal element symbol (e.g., "ZN", "FE", "TB")
            coordination_number: Filter by specific coordination number
            geometry: Filter by geometry (e.g., "tetrahedral", "octahedral")
            resolution_max: Maximum structure resolution in Angstroms
            limit: Maximum number of results to return

        Returns:
            List of MetalSiteResult objects, empty list if no results or error
        """
        metal = metal.upper()

        # Validate metal symbol
        if metal not in VALID_METALS:
            logger.warning(f"Unknown metal symbol: {metal}")
            return []

        if not REQUESTS_AVAILABLE:
            logger.warning("requests library not available - returning empty results")
            return []

        # Try MetalPDB first
        results = self._search_metalpdb(
            metal=metal,
            coordination_number=coordination_number,
            geometry=geometry,
            resolution_max=resolution_max,
            limit=limit,
        )

        # Fallback to RCSB if MetalPDB failed and fallback is enabled
        if not results and self.enable_fallback:
            logger.info(f"MetalPDB returned no results, falling back to RCSB for {metal}")
            results = self._search_rcsb_fallback(
                metal=metal,
                resolution_max=resolution_max,
                limit=limit,
            )

        return results

    def get_site_details(
        self,
        pdb_id: str,
        metal: str,
        site_index: int = 0,
    ) -> Optional[MetalSiteResult]:
        """
        Get detailed coordination information for a specific metal site.

        Args:
            pdb_id: 4-character PDB identifier
            metal: Metal element symbol
            site_index: Index of the metal site (0 for first)

        Returns:
            MetalSiteResult with detailed coordination info, or None if not found
        """
        pdb_id = pdb_id.upper()
        metal = metal.upper()

        if not REQUESTS_AVAILABLE:
            logger.warning("requests library not available - cannot fetch site details")
            return None

        # Try MetalPDB first
        result = self._get_metalpdb_site_details(pdb_id, metal, site_index)

        # Fallback to RCSB if MetalPDB failed and fallback is enabled
        if result is None and self.enable_fallback:
            logger.info(f"MetalPDB site details failed, falling back to RCSB for {pdb_id}")
            result = self._get_rcsb_site_details(pdb_id, metal, site_index)

        return result

    def search_by_coordination_motif(
        self,
        motif: str,
        metal: Optional[str] = None,
        limit: int = 20,
    ) -> List[MetalSiteResult]:
        """
        Search for metal sites by coordination motif pattern.

        Args:
            motif: Coordination motif pattern (e.g., "His-His-Glu-Asp", "Cys-Cys-His")
            metal: Optional metal filter
            limit: Maximum number of results

        Returns:
            List of MetalSiteResult objects matching the motif
        """
        if not REQUESTS_AVAILABLE:
            logger.warning("requests library not available - returning empty results")
            return []

        # Parse motif into residue list
        residues = [r.strip().upper() for r in motif.split("-")]

        # Normalize residue names (3-letter codes)
        residue_map = {
            "HIS": "HIS", "H": "HIS",
            "GLU": "GLU", "E": "GLU",
            "ASP": "ASP", "D": "ASP",
            "CYS": "CYS", "C": "CYS",
            "ASN": "ASN", "N": "ASN",
            "GLN": "GLN", "Q": "GLN",
            "TYR": "TYR", "Y": "TYR",
            "SER": "SER", "S": "SER",
            "THR": "THR", "T": "THR",
            "LYS": "LYS", "K": "LYS",
            "ARG": "ARG", "R": "ARG",
            "MET": "MET", "M": "MET",
        }
        normalized_residues = [residue_map.get(r[:3].upper(), r[:3].upper()) for r in residues]

        # Try MetalPDB motif search
        results = self._search_metalpdb_motif(normalized_residues, metal, limit)

        # Fallback: search by metal and filter by motif post-hoc
        if not results and self.enable_fallback and metal:
            logger.info(f"MetalPDB motif search failed, falling back to RCSB + filtering")
            all_results = self._search_rcsb_fallback(
                metal=metal,
                resolution_max=3.0,
                limit=limit * 5,  # Fetch more to filter
            )
            # Filter results that might match motif (basic filtering)
            results = self._filter_by_motif(all_results, normalized_residues, limit)

        return results

    # =========================================================================
    # METALPDB API METHODS
    # =========================================================================

    def _search_metalpdb(
        self,
        metal: str,
        coordination_number: Optional[int],
        geometry: Optional[str],
        resolution_max: float,
        limit: int,
    ) -> List[MetalSiteResult]:
        """Query MetalPDB search API."""
        try:
            # Build query parameters
            params = {
                "metal": metal,
                "resolution_max": resolution_max,
                "limit": limit,
            }
            if coordination_number is not None:
                params["coordination_number"] = coordination_number
            if geometry is not None:
                params["geometry"] = geometry

            url = f"{self.base_url}{METALPDB_API_SEARCH}"
            response = requests.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()

            data = response.json()
            return self._parse_metalpdb_search_results(data, metal)

        except requests.RequestException as e:
            logger.warning(f"MetalPDB search failed: {e}")
            return []
        except (ValueError, KeyError) as e:
            logger.warning(f"Error parsing MetalPDB response: {e}")
            return []

    def _get_metalpdb_site_details(
        self,
        pdb_id: str,
        metal: str,
        site_index: int,
    ) -> Optional[MetalSiteResult]:
        """Get site details from MetalPDB API."""
        try:
            url = f"{self.base_url}{METALPDB_API_SITE}/{pdb_id}/{metal}/{site_index}"
            response = requests.get(url, timeout=self.timeout)
            response.raise_for_status()

            data = response.json()
            return self._parse_metalpdb_site_details(data, pdb_id, metal, site_index)

        except requests.RequestException as e:
            logger.warning(f"MetalPDB site details failed: {e}")
            return None
        except (ValueError, KeyError) as e:
            logger.warning(f"Error parsing MetalPDB site details: {e}")
            return None

    def _search_metalpdb_motif(
        self,
        residues: List[str],
        metal: Optional[str],
        limit: int,
    ) -> List[MetalSiteResult]:
        """Search MetalPDB by coordination motif."""
        try:
            params = {
                "motif": "-".join(residues),
                "limit": limit,
            }
            if metal:
                params["metal"] = metal

            url = f"{self.base_url}{METALPDB_API_SEARCH}"
            response = requests.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()

            data = response.json()
            return self._parse_metalpdb_search_results(data, metal or "")

        except requests.RequestException as e:
            logger.warning(f"MetalPDB motif search failed: {e}")
            return []
        except (ValueError, KeyError) as e:
            logger.warning(f"Error parsing MetalPDB motif response: {e}")
            return []

    def _parse_metalpdb_search_results(
        self,
        data: Dict[str, Any],
        metal: str,
    ) -> List[MetalSiteResult]:
        """Parse MetalPDB search API response."""
        results = []

        # Handle different response formats
        sites = data.get("results", data.get("sites", data.get("data", [])))

        if not isinstance(sites, list):
            return results

        for site in sites:
            try:
                result = MetalSiteResult(
                    pdb_id=site.get("pdb_id", site.get("pdbId", "")).upper(),
                    metal=site.get("metal", metal).upper(),
                    coordination_number=site.get("coordination_number", site.get("coordNum", 0)),
                    geometry=site.get("geometry", site.get("geom", "")),
                    resolution=site.get("resolution", site.get("res", 0.0)),
                    coordinating_residues=site.get("coordinating_residues", site.get("ligands", [])),
                    ligands=site.get("non_protein_ligands", []),
                    source="metalpdb",
                    site_index=site.get("site_index", site.get("siteIdx", 0)),
                )
                if result.pdb_id:  # Only add if valid PDB ID
                    results.append(result)
            except (KeyError, TypeError) as e:
                logger.debug(f"Skipping malformed site entry: {e}")
                continue

        return results

    def _parse_metalpdb_site_details(
        self,
        data: Dict[str, Any],
        pdb_id: str,
        metal: str,
        site_index: int,
    ) -> Optional[MetalSiteResult]:
        """Parse MetalPDB site details response."""
        try:
            # Extract coordination info
            coord_residues = []
            distances = {}

            ligands_data = data.get("ligands", data.get("coordinating_atoms", []))
            for lig in ligands_data:
                res_name = lig.get("res_name", lig.get("resName", ""))
                res_num = lig.get("res_seq", lig.get("resNum", 0))
                res_str = f"{res_name}{res_num}"
                coord_residues.append(res_str)

                distance = lig.get("distance", lig.get("dist", 0.0))
                if distance:
                    distances[res_str] = distance

            # Get metal coordinates
            metal_coords = (
                data.get("metal_x", data.get("x", 0.0)),
                data.get("metal_y", data.get("y", 0.0)),
                data.get("metal_z", data.get("z", 0.0)),
            )

            return MetalSiteResult(
                pdb_id=pdb_id,
                metal=metal,
                coordination_number=data.get("coordination_number", len(coord_residues)),
                geometry=data.get("geometry", GEOMETRY_MAP.get(len(coord_residues), "unknown")),
                resolution=data.get("resolution", 0.0),
                coordinating_residues=coord_residues,
                ligands=data.get("non_protein_ligands", []),
                source="metalpdb",
                site_index=site_index,
                metal_coords=metal_coords,
                chain_id=data.get("chain_id", data.get("chain", "")),
                distances=distances,
            )
        except (KeyError, TypeError) as e:
            logger.warning(f"Error parsing MetalPDB site details: {e}")
            return None

    # =========================================================================
    # RCSB FALLBACK METHODS
    # =========================================================================

    def _search_rcsb_fallback(
        self,
        metal: str,
        resolution_max: float,
        limit: int,
    ) -> List[MetalSiteResult]:
        """Fallback to RCSB PDB search when MetalPDB unavailable."""
        try:
            # Import the existing metal_site_fetcher module
            import sys
            import os
            # Add parent directory to path for imports
            parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            if parent_dir not in sys.path:
                sys.path.insert(0, parent_dir)

            from metal_site_fetcher import query_metal_sites

            rcsb_results = query_metal_sites(
                metal=metal,
                resolution_max=resolution_max,
                limit=limit,
            )

            # Convert to MetalSiteResult objects
            results = []
            for site in rcsb_results:
                result = MetalSiteResult(
                    pdb_id=site.get("pdb_id", "").upper(),
                    metal=site.get("metal", metal).upper(),
                    resolution=site.get("resolution", 0.0),
                    source="rcsb",
                )
                if result.pdb_id:
                    results.append(result)

            return results

        except ImportError as e:
            logger.warning(f"Could not import metal_site_fetcher for RCSB fallback: {e}")
            return []
        except Exception as e:
            logger.warning(f"RCSB fallback search failed: {e}")
            return []

    def _get_rcsb_site_details(
        self,
        pdb_id: str,
        metal: str,
        site_index: int,
    ) -> Optional[MetalSiteResult]:
        """Get site details from RCSB as fallback."""
        try:
            # Import the existing metal_site_fetcher module
            import sys
            import os
            parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            if parent_dir not in sys.path:
                sys.path.insert(0, parent_dir)

            from metal_site_fetcher import generate_template_from_pdb

            template = generate_template_from_pdb(
                pdb_id=pdb_id,
                metal=metal,
                metal_site_index=site_index,
            )

            if not template:
                return None

            # Convert coordinating_residues format
            coord_residues = []
            distances = {}
            for res in template.get("coordinating_residues", []):
                res_str = f"{res.get('residue', '')}{res.get('position', '')}"
                coord_residues.append(res_str)
                if res.get("distance"):
                    distances[res_str] = res.get("distance")

            return MetalSiteResult(
                pdb_id=template.get("pdb_id", pdb_id).upper(),
                metal=template.get("metal", metal).upper(),
                coordination_number=template.get("coordination_number", 0),
                geometry=template.get("geometry", ""),
                coordinating_residues=coord_residues,
                source="rcsb",
                site_index=site_index,
                metal_coords=template.get("metal_coords", (0.0, 0.0, 0.0)),
                distances=distances,
            )

        except ImportError as e:
            logger.warning(f"Could not import metal_site_fetcher for RCSB fallback: {e}")
            return None
        except Exception as e:
            logger.warning(f"RCSB fallback site details failed: {e}")
            return None

    def _filter_by_motif(
        self,
        results: List[MetalSiteResult],
        target_residues: List[str],
        limit: int,
    ) -> List[MetalSiteResult]:
        """Filter results by coordination motif (post-hoc filtering)."""
        filtered = []

        for result in results:
            if not result.coordinating_residues:
                continue

            # Extract residue names from coordinating_residues
            site_residues = []
            for res_str in result.coordinating_residues:
                # Extract 3-letter code (e.g., "HIS94" -> "HIS")
                res_name = ""
                for i, char in enumerate(res_str):
                    if char.isdigit():
                        res_name = res_str[:i]
                        break
                if not res_name:
                    res_name = res_str[:3]
                site_residues.append(res_name.upper())

            # Check if motif matches (subset matching)
            target_counts = {}
            for r in target_residues:
                target_counts[r] = target_counts.get(r, 0) + 1

            site_counts = {}
            for r in site_residues:
                site_counts[r] = site_counts.get(r, 0) + 1

            # Check if site has at least the required residues
            matches = True
            for res, count in target_counts.items():
                if site_counts.get(res, 0) < count:
                    matches = False
                    break

            if matches:
                filtered.append(result)
                if len(filtered) >= limit:
                    break

        return filtered
