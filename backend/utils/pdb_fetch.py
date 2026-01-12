"""
PDB Fetcher Module

Fetch PDB and mmCIF files from RCSB Protein Data Bank.
Includes validation, caching, and error handling.
"""

import os
import re
import tempfile
from typing import Optional, Dict, Any, Tuple
from pathlib import Path
import hashlib
from datetime import datetime, timedelta

# Try to import requests, fall back to urllib if not available
try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    import urllib.request
    import urllib.error
    HAS_REQUESTS = False


# RCSB PDB URLs
RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
RCSB_CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"
RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"

# Cache settings
CACHE_DIR = os.path.join(tempfile.gettempdir(), "rfdiffusion_pdb_cache")
CACHE_EXPIRY_HOURS = 24


def validate_pdb_id(pdb_id: str) -> Tuple[bool, str]:
    """
    Validate PDB ID format.

    Valid PDB IDs are 4 characters: 1 digit + 3 alphanumeric.

    Args:
        pdb_id: The PDB ID to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not pdb_id:
        return False, "PDB ID cannot be empty"

    pdb_id = pdb_id.strip().upper()

    if len(pdb_id) != 4:
        return False, f"PDB ID must be 4 characters, got {len(pdb_id)}"

    # First character must be a digit (1-9)
    if not pdb_id[0].isdigit() or pdb_id[0] == '0':
        return False, "PDB ID must start with a digit (1-9)"

    # Remaining 3 characters must be alphanumeric
    if not pdb_id[1:].isalnum():
        return False, "PDB ID must contain only alphanumeric characters"

    return True, ""


def _get_cache_path(pdb_id: str, format: str = "pdb") -> Path:
    """Get the cache file path for a PDB ID."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    return Path(CACHE_DIR) / f"{pdb_id.upper()}.{format}"


def _is_cache_valid(cache_path: Path) -> bool:
    """Check if cache file exists and is not expired."""
    if not cache_path.exists():
        return False

    mtime = datetime.fromtimestamp(cache_path.stat().st_mtime)
    expiry = datetime.now() - timedelta(hours=CACHE_EXPIRY_HOURS)

    return mtime > expiry


def _fetch_url(url: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Fetch content from URL.

    Returns:
        Tuple of (content, error_message)
    """
    if HAS_REQUESTS:
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                return response.text, None
            elif response.status_code == 404:
                return None, f"PDB entry not found at {url}"
            else:
                return None, f"HTTP error {response.status_code}: {response.reason}"
        except requests.exceptions.Timeout:
            return None, "Request timed out"
        except requests.exceptions.RequestException as e:
            return None, f"Request error: {str(e)}"
    else:
        try:
            with urllib.request.urlopen(url, timeout=30) as response:
                return response.read().decode('utf-8'), None
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None, f"PDB entry not found at {url}"
            return None, f"HTTP error {e.code}: {e.reason}"
        except urllib.error.URLError as e:
            return None, f"URL error: {str(e)}"
        except Exception as e:
            return None, f"Error: {str(e)}"


def fetch_pdb(pdb_id: str, use_cache: bool = True) -> Dict[str, Any]:
    """
    Fetch PDB file from RCSB.

    Args:
        pdb_id: 4-character PDB ID (e.g., "1BRF")
        use_cache: Whether to use local cache

    Returns:
        Dict with keys:
        - success: bool
        - content: PDB file content (if successful)
        - error: Error message (if failed)
        - pdb_id: Normalized PDB ID
        - source: "cache" or "rcsb"
    """
    # Validate PDB ID
    is_valid, error = validate_pdb_id(pdb_id)
    if not is_valid:
        return {
            "success": False,
            "error": error,
            "pdb_id": pdb_id,
            "content": None,
            "source": None,
        }

    pdb_id = pdb_id.strip().upper()

    # Check cache
    if use_cache:
        cache_path = _get_cache_path(pdb_id, "pdb")
        if _is_cache_valid(cache_path):
            with open(cache_path, 'r') as f:
                return {
                    "success": True,
                    "content": f.read(),
                    "pdb_id": pdb_id,
                    "source": "cache",
                    "error": None,
                }

    # Fetch from RCSB
    url = RCSB_PDB_URL.format(pdb_id=pdb_id)
    content, error = _fetch_url(url)

    if error:
        return {
            "success": False,
            "error": error,
            "pdb_id": pdb_id,
            "content": None,
            "source": None,
        }

    # Cache the result
    if use_cache and content:
        cache_path = _get_cache_path(pdb_id, "pdb")
        with open(cache_path, 'w') as f:
            f.write(content)

    return {
        "success": True,
        "content": content,
        "pdb_id": pdb_id,
        "source": "rcsb",
        "error": None,
    }


def fetch_cif(pdb_id: str, use_cache: bool = True) -> Dict[str, Any]:
    """
    Fetch mmCIF file from RCSB.

    mmCIF format contains more metadata than PDB format.

    Args:
        pdb_id: 4-character PDB ID (e.g., "1BRF")
        use_cache: Whether to use local cache

    Returns:
        Dict with keys:
        - success: bool
        - content: CIF file content (if successful)
        - error: Error message (if failed)
        - pdb_id: Normalized PDB ID
        - source: "cache" or "rcsb"
    """
    # Validate PDB ID
    is_valid, error = validate_pdb_id(pdb_id)
    if not is_valid:
        return {
            "success": False,
            "error": error,
            "pdb_id": pdb_id,
            "content": None,
            "source": None,
        }

    pdb_id = pdb_id.strip().upper()

    # Check cache
    if use_cache:
        cache_path = _get_cache_path(pdb_id, "cif")
        if _is_cache_valid(cache_path):
            with open(cache_path, 'r') as f:
                return {
                    "success": True,
                    "content": f.read(),
                    "pdb_id": pdb_id,
                    "source": "cache",
                    "error": None,
                }

    # Fetch from RCSB
    url = RCSB_CIF_URL.format(pdb_id=pdb_id)
    content, error = _fetch_url(url)

    if error:
        return {
            "success": False,
            "error": error,
            "pdb_id": pdb_id,
            "content": None,
            "source": None,
        }

    # Cache the result
    if use_cache and content:
        cache_path = _get_cache_path(pdb_id, "cif")
        with open(cache_path, 'w') as f:
            f.write(content)

    return {
        "success": True,
        "content": content,
        "pdb_id": pdb_id,
        "source": "rcsb",
        "error": None,
    }


def get_pdb_metadata(pdb_id: str) -> Dict[str, Any]:
    """
    Fetch PDB entry metadata from RCSB REST API.

    Args:
        pdb_id: 4-character PDB ID

    Returns:
        Dict with metadata or error
    """
    is_valid, error = validate_pdb_id(pdb_id)
    if not is_valid:
        return {"success": False, "error": error}

    pdb_id = pdb_id.strip().upper()
    url = RCSB_ENTRY_URL.format(pdb_id=pdb_id)

    if HAS_REQUESTS:
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                return {"success": True, "metadata": response.json()}
            else:
                return {"success": False, "error": f"HTTP {response.status_code}"}
        except Exception as e:
            return {"success": False, "error": str(e)}
    else:
        try:
            import json
            with urllib.request.urlopen(url, timeout=30) as response:
                return {"success": True, "metadata": json.loads(response.read().decode())}
        except Exception as e:
            return {"success": False, "error": str(e)}


def clear_cache(pdb_id: Optional[str] = None) -> Dict[str, Any]:
    """
    Clear PDB cache.

    Args:
        pdb_id: Specific PDB ID to clear, or None to clear all

    Returns:
        Dict with success status and count of cleared files
    """
    if not os.path.exists(CACHE_DIR):
        return {"success": True, "cleared": 0}

    cleared = 0

    if pdb_id:
        pdb_id = pdb_id.strip().upper()
        for ext in ["pdb", "cif"]:
            path = _get_cache_path(pdb_id, ext)
            if path.exists():
                os.remove(path)
                cleared += 1
    else:
        for f in os.listdir(CACHE_DIR):
            os.remove(os.path.join(CACHE_DIR, f))
            cleared += 1

    return {"success": True, "cleared": cleared}


# Convenience function for parsing PDB content
def parse_pdb_content(pdb_content: str) -> Dict[str, Any]:
    """
    Quick parse of PDB content to extract basic info.

    Args:
        pdb_content: Raw PDB file content

    Returns:
        Dict with basic structure information
    """
    info = {
        "title": None,
        "num_atoms": 0,
        "num_residues": 0,
        "chains": set(),
        "hetero_residues": [],  # HETATM residues (ligands, metals, etc.)
        "metals": [],
    }

    residues_seen = set()

    # Common metal ions in PDB
    METALS = {
        "FE", "ZN", "MG", "CA", "MN", "CO", "CU", "NI", "MO",
        "NA", "K", "LA", "TB", "GD", "EU", "CE", "SM", "YB",
        "CD", "HG", "PB", "W", "V", "CR", "PT", "AU", "AG",
    }

    for line in pdb_content.split('\n'):
        if line.startswith('TITLE'):
            title_text = line[10:].strip()
            if info["title"]:
                info["title"] += " " + title_text
            else:
                info["title"] = title_text

        elif line.startswith('ATOM') or line.startswith('HETATM'):
            info["num_atoms"] += 1

            # Extract chain, residue number, and residue name
            try:
                chain = line[21]
                res_num = line[22:26].strip()
                res_name = line[17:20].strip()

                info["chains"].add(chain)
                residues_seen.add((chain, res_num))

                if line.startswith('HETATM'):
                    het_id = f"{chain}:{res_name}{res_num}"
                    if het_id not in info["hetero_residues"]:
                        info["hetero_residues"].append(het_id)

                    # Check for metals
                    element = line[76:78].strip().upper() if len(line) > 76 else res_name[:2].upper()
                    if element in METALS or res_name in METALS:
                        metal_info = {
                            "chain": chain,
                            "residue": res_name,
                            "res_num": res_num,
                            "element": element if element in METALS else res_name,
                        }
                        if metal_info not in info["metals"]:
                            info["metals"].append(metal_info)
            except (IndexError, ValueError):
                pass

    info["num_residues"] = len(residues_seen)
    info["chains"] = sorted(list(info["chains"]))

    return info
