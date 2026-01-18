# alphafold_adapter.py
"""
AlphaFold Database Adapter

Provides a unified interface for accessing predicted protein structures
from the AlphaFold Protein Structure Database (https://alphafold.ebi.ac.uk).

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

ALPHAFOLD_BASE_URL = "https://alphafold.ebi.ac.uk/api"
ALPHAFOLD_FILES_URL = "https://alphafold.ebi.ac.uk/files"

# Default configuration
DEFAULT_TIMEOUT = 30

# AlphaFold model version
DEFAULT_MODEL_VERSION = 4


# =============================================================================
# ALPHAFOLD ADAPTER CLASS
# =============================================================================

class AlphaFoldAdapter:
    """
    Adapter for accessing predicted protein structures from AlphaFold Database.

    Provides methods for:
    - Getting prediction metadata (pLDDT, model version, URLs)
    - Downloading predicted structures (PDB or CIF format)
    - Retrieving per-residue pLDDT confidence scores
    - Checking structure availability

    Example:
        >>> adapter = AlphaFoldAdapter()
        >>> prediction = adapter.get_prediction("P00918")
        >>> if prediction:
        ...     print(f"Model: {prediction['model_version']}, pLDDT: {prediction['plddt']}")
        >>> structure = adapter.download_structure("P00918", format="pdb")
        >>> if structure:
        ...     print(f"Downloaded {len(structure)} bytes")
    """

    def __init__(
        self,
        base_url: str = ALPHAFOLD_BASE_URL,
        files_url: str = ALPHAFOLD_FILES_URL,
        timeout: int = DEFAULT_TIMEOUT,
    ):
        """
        Initialize the AlphaFold adapter.

        Args:
            base_url: Base URL for AlphaFold API (default: official AlphaFold API)
            files_url: Base URL for AlphaFold structure files
            timeout: Request timeout in seconds
        """
        self.base_url = base_url
        self.files_url = files_url
        self.timeout = timeout

    def get_prediction(
        self,
        uniprot_id: str,
    ) -> Optional[Dict[str, Any]]:
        """
        Get prediction metadata for a UniProt accession.

        Args:
            uniprot_id: UniProt accession ID (e.g., "P00918")

        Returns:
            Dictionary with pdb_url, pae_url, plddt, model_version, etc.,
            or None if not found or error
        """
        if not REQUESTS_AVAILABLE:
            logger.warning("requests library not available - returning None")
            return None

        if not uniprot_id or not uniprot_id.strip():
            logger.warning("Empty UniProt ID provided to get_prediction")
            return None

        try:
            uniprot_id_clean = uniprot_id.strip().upper()
            url = f"{self.base_url}/prediction/{quote(uniprot_id_clean, safe='')}"

            response = requests.get(url, timeout=self.timeout)

            if response.status_code == 404:
                logger.info(f"AlphaFold prediction not found for: {uniprot_id_clean}")
                return None
            elif not response.ok:
                response.raise_for_status()

            data = response.json()
            return self._parse_prediction_response(data, uniprot_id_clean)

        except requests.RequestException as e:
            logger.warning(f"AlphaFold get prediction failed: {e}")
            return None
        except (ValueError, KeyError) as e:
            logger.warning(f"Error parsing AlphaFold response: {e}")
            return None

    def download_structure(
        self,
        uniprot_id: str,
        format: str = "pdb",
    ) -> Optional[str]:
        """
        Download predicted structure file.

        Args:
            uniprot_id: UniProt accession ID (e.g., "P00918")
            format: Structure format, either 'pdb' or 'cif' (default: 'pdb')

        Returns:
            Structure file content as string, or None if not found or error
        """
        if not REQUESTS_AVAILABLE:
            logger.warning("requests library not available - returning None")
            return None

        if not uniprot_id or not uniprot_id.strip():
            logger.warning("Empty UniProt ID provided to download_structure")
            return None

        format_lower = format.lower().strip()
        if format_lower not in ("pdb", "cif"):
            logger.warning(f"Invalid format: {format}. Must be 'pdb' or 'cif'")
            return None

        try:
            uniprot_id_clean = uniprot_id.strip().upper()

            # Construct file URL based on format
            if format_lower == "pdb":
                filename = f"AF-{uniprot_id_clean}-F1-model_v{DEFAULT_MODEL_VERSION}.pdb"
            else:
                filename = f"AF-{uniprot_id_clean}-F1-model_v{DEFAULT_MODEL_VERSION}.cif"

            url = f"{self.files_url}/{filename}"

            response = requests.get(url, timeout=self.timeout)

            if response.status_code == 404:
                logger.info(f"AlphaFold structure not found for: {uniprot_id_clean}")
                return None
            elif not response.ok:
                response.raise_for_status()

            return response.text

        except requests.RequestException as e:
            logger.warning(f"AlphaFold download structure failed: {e}")
            return None

    def get_plddt_scores(
        self,
        uniprot_id: str,
    ) -> Optional[List[float]]:
        """
        Get per-residue pLDDT confidence scores.

        Downloads the PDB file and extracts B-factor column which contains
        pLDDT scores for AlphaFold structures.

        Args:
            uniprot_id: UniProt accession ID (e.g., "P00918")

        Returns:
            List of pLDDT scores per residue (0-100 scale),
            or None if not found or error
        """
        if not REQUESTS_AVAILABLE:
            logger.warning("requests library not available - returning None")
            return None

        if not uniprot_id or not uniprot_id.strip():
            logger.warning("Empty UniProt ID provided to get_plddt_scores")
            return None

        try:
            # Download PDB structure
            pdb_content = self.download_structure(uniprot_id, format="pdb")

            if not pdb_content:
                return None

            # Parse PDB and extract B-factors (pLDDT scores)
            return self._extract_plddt_from_pdb(pdb_content)

        except Exception as e:
            logger.warning(f"Error extracting pLDDT scores: {e}")
            return None

    def check_availability(
        self,
        uniprot_id: str,
    ) -> bool:
        """
        Check if an AlphaFold prediction exists for a UniProt accession.

        Args:
            uniprot_id: UniProt accession ID (e.g., "P00918")

        Returns:
            True if prediction exists, False otherwise
        """
        if not REQUESTS_AVAILABLE:
            logger.warning("requests library not available - returning False")
            return False

        if not uniprot_id or not uniprot_id.strip():
            logger.warning("Empty UniProt ID provided to check_availability")
            return False

        try:
            uniprot_id_clean = uniprot_id.strip().upper()
            url = f"{self.base_url}/prediction/{quote(uniprot_id_clean, safe='')}"

            # Use HEAD request for efficiency
            response = requests.head(url, timeout=self.timeout)

            # AlphaFold API may not support HEAD, fallback to GET
            if response.status_code == 405:  # Method Not Allowed
                response = requests.get(url, timeout=self.timeout)

            return response.status_code == 200

        except requests.RequestException as e:
            logger.warning(f"AlphaFold availability check failed: {e}")
            return False

    # =========================================================================
    # INTERNAL METHODS
    # =========================================================================

    def _parse_prediction_response(
        self,
        data: Any,
        uniprot_id: str,
    ) -> Optional[Dict[str, Any]]:
        """Parse AlphaFold prediction API response."""
        try:
            # AlphaFold API returns a list of predictions for fragments
            # Most proteins have a single entry
            if isinstance(data, list):
                if not data:
                    return None
                # Use first entry (main prediction)
                entry = data[0]
            else:
                entry = data

            # Extract key fields
            model_url = entry.get("pdbUrl", entry.get("pdb_url", ""))
            cif_url = entry.get("cifUrl", entry.get("cif_url", ""))
            pae_url = entry.get("paeImageUrl", entry.get("pae_image_url", ""))
            pae_doc_url = entry.get("paeDocUrl", entry.get("pae_doc_url", ""))

            # Extract version from URL or use default
            model_version = DEFAULT_MODEL_VERSION
            if model_url:
                # Try to extract version from URL like "..._v4.pdb"
                if "_v" in model_url:
                    try:
                        version_part = model_url.split("_v")[-1]
                        version_str = version_part.split(".")[0]
                        model_version = int(version_str)
                    except (ValueError, IndexError):
                        pass

            # Get global pLDDT if available
            global_plddt = entry.get("globalMetricValue", entry.get("plddt", None))

            # Get sequence length
            seq_length = entry.get("uniprotEnd", 0) - entry.get("uniprotStart", 0)
            if seq_length <= 0:
                seq_length = entry.get("sequenceLength", None)

            result = {
                "uniprot_id": entry.get("uniprotAccession", uniprot_id),
                "pdb_url": model_url,
                "cif_url": cif_url,
                "pae_url": pae_url,
                "pae_doc_url": pae_doc_url,
                "model_version": model_version,
                "plddt": global_plddt,
                "sequence_length": seq_length,
                "gene": entry.get("gene", ""),
                "organism": entry.get("organismScientificName", ""),
                "entry_id": entry.get("entryId", ""),
            }

            return result

        except (KeyError, TypeError) as e:
            logger.warning(f"Error parsing prediction response: {e}")
            return None

    def _extract_plddt_from_pdb(
        self,
        pdb_content: str,
    ) -> List[float]:
        """
        Extract per-residue pLDDT scores from PDB B-factor column.

        AlphaFold stores pLDDT scores (0-100) in the B-factor column.
        We extract one score per residue (from CA atoms).
        """
        plddt_scores = []
        seen_residues = set()

        for line in pdb_content.split("\n"):
            # Only process ATOM records
            if not line.startswith("ATOM"):
                continue

            try:
                # PDB format: columns 13-16 = atom name, 23-26 = residue number
                # B-factor is in columns 61-66
                atom_name = line[12:16].strip()

                # Use CA atom for per-residue score
                if atom_name != "CA":
                    continue

                # Get residue identifier (chain + residue number)
                chain_id = line[21]
                res_num = line[22:26].strip()
                res_key = f"{chain_id}_{res_num}"

                # Skip if already seen this residue
                if res_key in seen_residues:
                    continue
                seen_residues.add(res_key)

                # Extract B-factor (pLDDT score)
                b_factor_str = line[60:66].strip()
                b_factor = float(b_factor_str)
                plddt_scores.append(b_factor)

            except (ValueError, IndexError) as e:
                # Skip malformed lines
                logger.debug(f"Skipping malformed PDB line: {e}")
                continue

        return plddt_scores
