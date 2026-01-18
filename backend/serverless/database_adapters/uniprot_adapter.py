# uniprot_adapter.py
"""
UniProt Database Adapter

Provides a unified interface for querying protein functions and metal binding
sites from the UniProt database (https://www.uniprot.org).

This adapter is part of the AI-Driven Structure Discovery Infrastructure
for protein design workflows.
"""

from typing import Dict, List, Optional, Any
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

UNIPROT_BASE_URL = "https://rest.uniprot.org"
UNIPROT_SEARCH_ENDPOINT = "/uniprotkb/search"
UNIPROT_ENTRY_ENDPOINT = "/uniprotkb"

# Query fields for UniProt REST API
# Note: ft_binding includes metal binding sites; xref_pdb for PDB cross-references
DEFAULT_QUERY_FIELDS = "accession,protein_name,organism_name,xref_pdb,sequence"

# Default configuration
DEFAULT_TIMEOUT = 30
DEFAULT_RESULT_LIMIT = 20


# =============================================================================
# UNIPROT ADAPTER CLASS
# =============================================================================

class UniProtAdapter:
    """
    Adapter for querying protein data from the UniProt database.

    Provides methods for searching proteins by:
    - Function keywords (e.g., "zinc finger", "dehydrogenase")
    - Metal binding properties
    - Getting PDB structure mappings
    - Extracting metal binding site annotations
    - Retrieving protein sequences

    Example:
        >>> adapter = UniProtAdapter()
        >>> results = adapter.search_by_function("zinc finger", organism="human", limit=10)
        >>> for protein in results:
        ...     print(f"{protein['uniprot_id']}: {protein['name']}")
    """

    def __init__(
        self,
        base_url: str = UNIPROT_BASE_URL,
        timeout: int = DEFAULT_TIMEOUT,
    ):
        """
        Initialize the UniProt adapter.

        Args:
            base_url: Base URL for UniProt REST API (default: official UniProt)
            timeout: Request timeout in seconds
        """
        self.base_url = base_url
        self.timeout = timeout

    def search_by_function(
        self,
        function_query: str,
        organism: Optional[str] = None,
        reviewed_only: bool = True,
        limit: int = DEFAULT_RESULT_LIMIT,
    ) -> List[Dict[str, Any]]:
        """
        Search UniProt for proteins by function.

        Args:
            function_query: Function keyword to search (e.g., "zinc finger", "dehydrogenase")
            organism: Optional organism filter (e.g., "human", "Escherichia coli")
            reviewed_only: If True, only return Swiss-Prot (reviewed) entries
            limit: Maximum number of results to return

        Returns:
            List of dictionaries with uniprot_id, name, function, organism
        """
        if not REQUESTS_AVAILABLE:
            logger.warning("requests library not available - returning empty results")
            return []

        try:
            # Build query string
            # Note: cc_function accepts text without quotes
            query_parts = [f"(cc_function:{function_query})"]

            if reviewed_only:
                query_parts.append("(reviewed:true)")

            if organism:
                query_parts.append(f"(organism_name:{organism})")

            query = " AND ".join(query_parts)

            # Make API request
            params = {
                "query": query,
                "format": "json",
                "fields": DEFAULT_QUERY_FIELDS,
                "size": limit,
            }

            url = f"{self.base_url}{UNIPROT_SEARCH_ENDPOINT}"
            response = requests.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()

            data = response.json()
            return self._parse_search_results(data)

        except requests.RequestException as e:
            logger.warning(f"UniProt search by function failed: {e}")
            return []
        except (ValueError, KeyError) as e:
            logger.warning(f"Error parsing UniProt response: {e}")
            return []

    def search_by_metal_binding(
        self,
        metal: str,
        organism: Optional[str] = None,
        limit: int = DEFAULT_RESULT_LIMIT,
    ) -> List[Dict[str, Any]]:
        """
        Search UniProt for proteins that bind a specific metal.

        Args:
            metal: Metal name to search (e.g., "Zinc", "Iron", "Copper")
            organism: Optional organism filter
            limit: Maximum number of results to return

        Returns:
            List of dictionaries with uniprot_id, name, organism
        """
        if not REQUESTS_AVAILABLE:
            logger.warning("requests library not available - returning empty results")
            return []

        try:
            # Build query using ft_binding feature type for metal binding
            # Note: ft_binding is the correct field for binding site features including metals
            query_parts = [f"(ft_binding:{metal})"]

            if organism:
                query_parts.append(f"(organism_name:{organism})")

            query = " AND ".join(query_parts)

            params = {
                "query": query,
                "format": "json",
                "fields": DEFAULT_QUERY_FIELDS,
                "size": limit,
            }

            url = f"{self.base_url}{UNIPROT_SEARCH_ENDPOINT}"
            response = requests.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()

            data = response.json()
            return self._parse_search_results(data)

        except requests.RequestException as e:
            logger.warning(f"UniProt search by metal binding failed: {e}")
            return []
        except (ValueError, KeyError) as e:
            logger.warning(f"Error parsing UniProt response: {e}")
            return []

    def get_pdb_mappings(
        self,
        uniprot_id: str,
    ) -> List[Dict[str, Any]]:
        """
        Get PDB structure mappings for a UniProt entry.

        Args:
            uniprot_id: UniProt accession ID (e.g., "P00918")

        Returns:
            List of dictionaries with pdb_id, method, resolution, chains
        """
        if not REQUESTS_AVAILABLE:
            logger.warning("requests library not available - returning empty results")
            return []

        try:
            # Fetch entry with cross-references
            url = f"{self.base_url}{UNIPROT_ENTRY_ENDPOINT}/{uniprot_id}"
            params = {
                "format": "json",
                "fields": "xref_pdb",
            }

            response = requests.get(url, params=params, timeout=self.timeout)

            if response.status_code == 404:
                logger.info(f"UniProt entry not found: {uniprot_id}")
                return []
            elif not response.ok:
                response.raise_for_status()

            data = response.json()
            return self._parse_pdb_mappings(data)

        except requests.RequestException as e:
            logger.warning(f"UniProt get PDB mappings failed: {e}")
            return []
        except (ValueError, KeyError) as e:
            logger.warning(f"Error parsing UniProt PDB mappings: {e}")
            return []

    def get_metal_binding_sites(
        self,
        uniprot_id: str,
    ) -> List[Dict[str, Any]]:
        """
        Get metal binding site annotations for a UniProt entry.

        Args:
            uniprot_id: UniProt accession ID (e.g., "P00918")

        Returns:
            List of dictionaries with position, metal, description
        """
        if not REQUESTS_AVAILABLE:
            logger.warning("requests library not available - returning empty results")
            return []

        try:
            # Fetch entry with feature annotations
            # Note: ft_binding includes all binding sites (metal and substrate)
            url = f"{self.base_url}{UNIPROT_ENTRY_ENDPOINT}/{uniprot_id}"
            params = {
                "format": "json",
                "fields": "ft_binding",
            }

            response = requests.get(url, params=params, timeout=self.timeout)

            if response.status_code == 404:
                logger.info(f"UniProt entry not found: {uniprot_id}")
                return []
            elif not response.ok:
                response.raise_for_status()

            data = response.json()
            return self._parse_metal_binding_sites(data)

        except requests.RequestException as e:
            logger.warning(f"UniProt get metal binding sites failed: {e}")
            return []
        except (ValueError, KeyError) as e:
            logger.warning(f"Error parsing UniProt metal binding sites: {e}")
            return []

    def get_sequence(
        self,
        uniprot_id: str,
    ) -> Optional[str]:
        """
        Get the protein sequence for a UniProt entry.

        Args:
            uniprot_id: UniProt accession ID (e.g., "P00918")

        Returns:
            Protein sequence string, or None if not found
        """
        if not REQUESTS_AVAILABLE:
            logger.warning("requests library not available - returning None")
            return None

        try:
            # Fetch entry with sequence
            url = f"{self.base_url}{UNIPROT_ENTRY_ENDPOINT}/{uniprot_id}"
            params = {
                "format": "json",
                "fields": "sequence",
            }

            response = requests.get(url, params=params, timeout=self.timeout)

            if response.status_code == 404:
                logger.info(f"UniProt entry not found: {uniprot_id}")
                return None
            elif not response.ok:
                response.raise_for_status()

            data = response.json()
            return self._parse_sequence(data)

        except requests.RequestException as e:
            logger.warning(f"UniProt get sequence failed: {e}")
            return None
        except (ValueError, KeyError) as e:
            logger.warning(f"Error parsing UniProt sequence: {e}")
            return None

    # =========================================================================
    # PARSING METHODS
    # =========================================================================

    def _parse_search_results(
        self,
        data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Parse UniProt search API response."""
        results = []

        # Handle search response format
        entries = data.get("results", [])

        for entry in entries:
            try:
                # Extract UniProt accession
                uniprot_id = entry.get("primaryAccession", "")

                # Extract protein name
                protein_desc = entry.get("proteinDescription", {})
                recommended_name = protein_desc.get("recommendedName", {})
                full_name = recommended_name.get("fullName", {})
                name = full_name.get("value", "Unknown")

                # Fallback to submitted name if no recommended name
                if name == "Unknown":
                    submitted_names = protein_desc.get("submissionNames", [])
                    if submitted_names:
                        name = submitted_names[0].get("fullName", {}).get("value", "Unknown")

                # Extract organism
                organism_info = entry.get("organism", {})
                organism = organism_info.get("scientificName", "")

                # Extract function from comments
                function = ""
                comments = entry.get("comments", [])
                for comment in comments:
                    if comment.get("commentType") == "FUNCTION":
                        texts = comment.get("texts", [])
                        if texts:
                            function = texts[0].get("value", "")
                        break

                result = {
                    "uniprot_id": uniprot_id,
                    "name": name,
                    "organism": organism,
                    "function": function,
                }

                if uniprot_id:
                    results.append(result)

            except (KeyError, TypeError) as e:
                logger.debug(f"Skipping malformed entry: {e}")
                continue

        return results

    def _parse_pdb_mappings(
        self,
        data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Parse PDB cross-references from UniProt entry."""
        results = []

        # Handle single entry format
        cross_refs = data.get("uniProtKBCrossReferences", [])

        for ref in cross_refs:
            try:
                if ref.get("database") == "PDB":
                    pdb_id = ref.get("id", "")

                    # Extract properties
                    properties = ref.get("properties", [])
                    method = ""
                    resolution = ""
                    chains = ""

                    for prop in properties:
                        key = prop.get("key", "")
                        value = prop.get("value", "")
                        if key == "Method":
                            method = value
                        elif key == "Resolution":
                            resolution = value
                        elif key == "Chains":
                            chains = value

                    result = {
                        "pdb_id": pdb_id.upper(),
                        "method": method,
                        "resolution": resolution,
                        "chains": chains,
                    }

                    if pdb_id:
                        results.append(result)

            except (KeyError, TypeError) as e:
                logger.debug(f"Skipping malformed PDB reference: {e}")
                continue

        return results

    def _parse_metal_binding_sites(
        self,
        data: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Parse metal binding site features from UniProt entry."""
        results = []

        # Extract features from entry
        features = data.get("features", [])

        for feature in features:
            try:
                feature_type = feature.get("type", "")

                # Look for metal binding features
                if feature_type in ("Metal binding", "Binding site"):
                    location = feature.get("location", {})
                    start = location.get("start", {}).get("value")
                    end = location.get("end", {}).get("value")

                    # Position can be single residue or range
                    if start == end:
                        position = str(start) if start else ""
                    else:
                        position = f"{start}-{end}" if start and end else ""

                    # Extract description/ligand
                    description = feature.get("description", "")
                    ligand = feature.get("ligand", {})
                    metal = ligand.get("name", "")

                    # Fallback to description if no explicit ligand
                    if not metal and description:
                        metal = description

                    result = {
                        "position": position,
                        "metal": metal,
                        "description": description,
                        "type": feature_type,
                    }

                    if position or metal:
                        results.append(result)

            except (KeyError, TypeError) as e:
                logger.debug(f"Skipping malformed feature: {e}")
                continue

        return results

    def _parse_sequence(
        self,
        data: Dict[str, Any],
    ) -> Optional[str]:
        """Parse sequence from UniProt entry."""
        try:
            sequence_info = data.get("sequence", {})
            sequence = sequence_info.get("value", "")
            return sequence if sequence else None
        except (KeyError, TypeError) as e:
            logger.debug(f"Error parsing sequence: {e}")
            return None
