# pubchem_adapter.py
"""
PubChem Database Adapter

Provides a unified interface for querying chemical compound information
from PubChem (https://pubchem.ncbi.nlm.nih.gov), including SMILES structures
for ligands used in protein design workflows.

This adapter is part of the AI-Driven Structure Discovery Infrastructure
for protein design workflows.
"""

from typing import Dict, List, Optional, Any
import logging
from urllib.parse import quote

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

PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

# Default configuration
DEFAULT_TIMEOUT = 30
DEFAULT_RESULT_LIMIT = 10

# Common ligands lookup table for fast retrieval
# Maps ligand name (uppercase) to PubChem CID
COMMON_LIGANDS = {
    "ATP": 5957,
    "ADP": 6022,
    "GTP": 6830,
    "NAD": 5893,
    "NADH": 928,
    "FAD": 643975,
    "HEME": 26945,
    "CITRATE": 311,
    "CITRIC ACID": 311,
    "PQQ": 1024,
    "AZOBENZENE": 2272,
}


# =============================================================================
# PUBCHEM ADAPTER CLASS
# =============================================================================

class PubChemAdapter:
    """
    Adapter for querying chemical compound data from the PubChem database.

    Provides methods for searching compounds by:
    - Name (common or IUPAC name)
    - PubChem CID (Compound ID)
    - Substructure SMILES pattern

    Also provides fast lookup for common ligands used in protein design.

    Example:
        >>> adapter = PubChemAdapter()
        >>> results = adapter.search_by_name("citric acid", limit=5)
        >>> for compound in results:
        ...     print(f"{compound['cid']}: {compound['smiles']}")
        >>> smiles = adapter.get_smiles("ATP")
        >>> print(smiles)
    """

    def __init__(
        self,
        base_url: str = PUBCHEM_BASE_URL,
        timeout: int = DEFAULT_TIMEOUT,
    ):
        """
        Initialize the PubChem adapter.

        Args:
            base_url: Base URL for PubChem PUG REST API (default: official PubChem)
            timeout: Request timeout in seconds
        """
        self.base_url = base_url
        self.timeout = timeout

    def search_by_name(
        self,
        name: str,
        limit: int = DEFAULT_RESULT_LIMIT,
    ) -> List[Dict[str, Any]]:
        """
        Search PubChem for compounds by name.

        Args:
            name: Compound name to search (common or IUPAC name)
            limit: Maximum number of results to return

        Returns:
            List of dictionaries with cid, smiles, and name
        """
        if not REQUESTS_AVAILABLE:
            logger.warning("requests library not available - returning empty results")
            return []

        if not name or not name.strip():
            logger.warning("Empty name provided to search_by_name")
            return []

        try:
            # First, search for CIDs by name
            name_clean = name.strip()
            name_encoded = quote(name_clean, safe='')
            url = f"{self.base_url}/compound/name/{name_encoded}/cids/JSON"

            response = requests.get(url, timeout=self.timeout)
            response.raise_for_status()

            data = response.json()
            cids = data.get("IdentifierList", {}).get("CID", [])

            if not cids:
                return []

            # Limit CIDs
            cids = cids[:limit]

            # Get properties for each CID
            results = []
            for cid in cids:
                compound_info = self._get_compound_properties(cid)
                if compound_info:
                    results.append(compound_info)

            return results

        except requests.RequestException as e:
            logger.warning(f"PubChem search by name failed: {e}")
            return []
        except (ValueError, KeyError) as e:
            logger.warning(f"Error parsing PubChem response: {e}")
            return []

    def get_compound(
        self,
        cid: int,
    ) -> Optional[Dict[str, Any]]:
        """
        Get compound information by PubChem CID.

        Args:
            cid: PubChem Compound ID

        Returns:
            Dictionary with smiles, molecular_formula, molecular_weight, iupac_name,
            or None if not found
        """
        if not REQUESTS_AVAILABLE:
            logger.warning("requests library not available - returning None")
            return None

        if cid <= 0:
            logger.warning(f"Invalid CID: {cid}")
            return None

        return self._get_compound_properties(cid)

    def get_smiles(
        self,
        ligand_name: str,
    ) -> Optional[str]:
        """
        Get SMILES structure for a ligand by name.

        Uses fast lookup table for common ligands (ATP, NAD, etc.) before
        falling back to PubChem API search.

        Args:
            ligand_name: Ligand name (case-insensitive)

        Returns:
            SMILES string, or None if not found
        """
        if not ligand_name or not ligand_name.strip():
            logger.warning("Empty ligand name provided to get_smiles")
            return None

        ligand_upper = ligand_name.strip().upper()

        # Check common ligands lookup table first
        if ligand_upper in COMMON_LIGANDS:
            cid = COMMON_LIGANDS[ligand_upper]
            compound = self.get_compound(cid)
            if compound:
                return compound.get("smiles")

        # Fall back to API search
        results = self.search_by_name(ligand_name, limit=1)
        if results:
            return results[0].get("smiles")

        return None

    def search_by_substructure(
        self,
        smiles: str,
        limit: int = DEFAULT_RESULT_LIMIT,
    ) -> List[Dict[str, Any]]:
        """
        Search PubChem for compounds containing a substructure.

        Args:
            smiles: SMILES pattern for substructure search
            limit: Maximum number of results to return

        Returns:
            List of dictionaries with cid, smiles, and name
        """
        if not REQUESTS_AVAILABLE:
            logger.warning("requests library not available - returning empty results")
            return []

        if not smiles or not smiles.strip():
            logger.warning("Empty SMILES provided to search_by_substructure")
            return []

        try:
            # PubChem substructure search is asynchronous
            # First, submit the search and get a list key
            smiles_clean = smiles.strip()
            smiles_encoded = quote(smiles_clean, safe='')
            url = f"{self.base_url}/compound/substructure/smiles/{smiles_encoded}/cids/JSON"

            response = requests.get(url, timeout=self.timeout)
            response.raise_for_status()

            data = response.json()

            # Check for waiting status (async search)
            if "Waiting" in data:
                # Poll for results
                list_key = data["Waiting"].get("ListKey")
                if list_key:
                    cids = self._poll_substructure_results(list_key, limit)
                else:
                    cids = []
            else:
                # Direct response
                cids = data.get("IdentifierList", {}).get("CID", [])

            if not cids:
                return []

            # Limit CIDs
            cids = cids[:limit]

            # Get properties for each CID
            results = []
            for cid in cids:
                compound_info = self._get_compound_properties(cid)
                if compound_info:
                    results.append(compound_info)

            return results

        except requests.RequestException as e:
            logger.warning(f"PubChem substructure search failed: {e}")
            return []
        except (ValueError, KeyError) as e:
            logger.warning(f"Error parsing PubChem substructure response: {e}")
            return []

    # =========================================================================
    # INTERNAL METHODS
    # =========================================================================

    def _get_compound_properties(
        self,
        cid: int,
    ) -> Optional[Dict[str, Any]]:
        """Get compound properties from PubChem by CID."""
        try:
            properties = "MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,IUPACName"
            url = f"{self.base_url}/compound/cid/{cid}/property/{properties}/JSON"

            response = requests.get(url, timeout=self.timeout)
            response.raise_for_status()

            data = response.json()
            props_list = data.get("PropertyTable", {}).get("Properties", [])

            if not props_list:
                return None

            props = props_list[0]

            # Prefer IsomericSMILES, fall back to CanonicalSMILES
            smiles = props.get("IsomericSMILES") or props.get("CanonicalSMILES", "")

            return {
                "cid": props.get("CID", cid),
                "smiles": smiles,
                "molecular_formula": props.get("MolecularFormula", ""),
                "molecular_weight": props.get("MolecularWeight", 0.0),
                "iupac_name": props.get("IUPACName", ""),
            }

        except requests.RequestException as e:
            logger.warning(f"PubChem get compound properties failed for CID {cid}: {e}")
            return None
        except (ValueError, KeyError) as e:
            logger.warning(f"Error parsing PubChem compound properties: {e}")
            return None

    def _poll_substructure_results(
        self,
        list_key: str,
        limit: int,
        max_attempts: int = 5,
        poll_interval: float = 1.0,
    ) -> List[int]:
        """Poll for asynchronous substructure search results."""
        import time

        url = f"{self.base_url}/compound/listkey/{list_key}/cids/JSON"

        for attempt in range(max_attempts):
            try:
                response = requests.get(url, timeout=self.timeout)
                response.raise_for_status()

                data = response.json()

                # Check if still waiting
                if "Waiting" in data:
                    # Exponential backoff: 1.0s, 1.5s, 2.25s, 3.375s, 5.0s
                    poll_interval = min(1.0 * (1.5 ** attempt), 5.0)
                    time.sleep(poll_interval)
                    continue

                # Got results
                cids = data.get("IdentifierList", {}).get("CID", [])
                return cids[:limit]

            except requests.RequestException as e:
                logger.warning(f"PubChem poll attempt {attempt + 1} failed: {e}")
                if attempt < max_attempts - 1:
                    # Exponential backoff: 1.0s, 1.5s, 2.25s, 3.375s, 5.0s
                    poll_interval = min(1.0 * (1.5 ** attempt), 5.0)
                    time.sleep(poll_interval)
                continue

        logger.warning(f"PubChem substructure search timed out after {max_attempts} attempts")
        return []
