"""
Lightweight M-CSA (Mechanism and Catalytic Site Atlas) client for serverless.

No pydantic dependency - returns plain dicts.
Mirrors backend/app/services/mcsa_client.py but for the serverless environment.

M-CSA API structure:
- /api/entries/?format=json — paginated list of all enzyme entries
- /api/entries/<mcsa_id>/?format=json — single entry with residues
- Entries have 'residues' list, each with 'residue_chains' (LIST of chain mappings)
- residue_chains entries have: chain_name, pdb_id, resid, code, auth_resid

The API doesn't support direct PDB->residue queries. We search entries
and match by PDB ID within the residue_chains data.
"""

import logging
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

MCSA_API_BASE = "https://www.ebi.ac.uk/thornton-srv/m-csa/api"

# 3-letter to 3-letter normalization (M-CSA returns mixed case like "Asp")
_CODE_MAP = {
    "Ala": "ALA", "Arg": "ARG", "Asn": "ASN", "Asp": "ASP",
    "Cys": "CYS", "Gln": "GLN", "Glu": "GLU", "Gly": "GLY",
    "His": "HIS", "Ile": "ILE", "Leu": "LEU", "Lys": "LYS",
    "Met": "MET", "Phe": "PHE", "Pro": "PRO", "Ser": "SER",
    "Thr": "THR", "Trp": "TRP", "Tyr": "TYR", "Val": "VAL",
}


def _normalize_code(code: str) -> str:
    """Normalize residue code: 'Asp' -> 'ASP'."""
    return _CODE_MAP.get(code, code.upper() if code else "UNK")


async def query_mcsa_by_pdb(pdb_id: str) -> List[Dict]:
    """
    Query M-CSA for catalytic residues by PDB ID.

    Searches M-CSA entries and extracts residues mapped to the given PDB.
    The M-CSA API returns entries organized by enzyme mechanism, with
    residue_chains listing which PDB structures contain each catalytic residue.

    Args:
        pdb_id: 4-character PDB ID (e.g., "4CVB")

    Returns:
        List of dicts with keys: chain, resnum, resname, role, confidence, source
    """
    if not pdb_id or len(pdb_id) != 4:
        return []

    pdb_id_lower = pdb_id.lower()

    try:
        import httpx
        async with httpx.AsyncClient(timeout=15.0) as client:
            # Search entries — the pdb_id filter in the API isn't reliable,
            # so we fetch entries and check residue_chains for our PDB
            url = f"{MCSA_API_BASE}/entries/?format=json"
            response = await client.get(url)

            if response.status_code != 200:
                logger.warning(f"[M-CSA] HTTP {response.status_code}")
                return []

            data = response.json()

            # API returns paginated: {count, next, previous, results}
            entries = []
            if isinstance(data, dict):
                entries = data.get("results", [])
            elif isinstance(data, list):
                entries = data

            if not entries:
                return []

            residues = []
            seen = set()

            for entry in entries:
                role_summary = entry.get("roles_summary", "")

                # Each entry can have a 'residues' list
                entry_residues = entry.get("residues", [])
                if not entry_residues:
                    # Top-level entry might itself be a residue record
                    entry_residues = [entry]

                for res_record in entry_residues:
                    role = res_record.get("roles_summary", role_summary)

                    # residue_chains is a LIST of chain mappings
                    chains_list = res_record.get("residue_chains", [])
                    if isinstance(chains_list, dict):
                        chains_list = [chains_list]

                    for chain_info in chains_list:
                        # Check if this chain mapping is for our PDB
                        entry_pdb = (chain_info.get("pdb_id") or "").lower()
                        if entry_pdb != pdb_id_lower:
                            continue

                        chain_name = chain_info.get("chain_name", "A")
                        resid = chain_info.get("auth_resid") or chain_info.get("resid")
                        code = chain_info.get("code", "UNK")

                        if resid is None:
                            continue

                        key = f"{chain_name}{resid}"
                        if key in seen:
                            continue
                        seen.add(key)

                        residues.append({
                            "chain": chain_name,
                            "resnum": int(resid),
                            "resname": _normalize_code(code),
                            "role": role or None,
                            "confidence": 1.0,
                            "source": "mcsa",
                        })

            if residues:
                logger.info(f"[M-CSA] Found {len(residues)} catalytic residues for {pdb_id}")
            else:
                logger.info(f"[M-CSA] No entries found for PDB {pdb_id}")

            return residues

    except ImportError:
        logger.warning("[M-CSA] httpx not available, skipping M-CSA query")
        return []
    except Exception as e:
        logger.warning(f"[M-CSA] Query failed for {pdb_id}: {e}")
        return []
