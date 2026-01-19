"""
M-CSA (Mechanism and Catalytic Site Atlas) API client.
Queries curated catalytic residue data by PDB ID.
"""

import httpx
from typing import List, Optional
from pydantic import BaseModel


class CatalyticResidue(BaseModel):
    """A catalytic residue from M-CSA."""
    chain: str
    residue: int
    name: str
    role: Optional[str] = None
    confidence: float = 1.0
    source: str = "mcsa"


MCSA_API_BASE = "https://www.ebi.ac.uk/thornton-srv/m-csa/api"


async def query_mcsa_by_pdb(pdb_id: str) -> List[CatalyticResidue]:
    """
    Query M-CSA for catalytic residues by PDB ID.

    Args:
        pdb_id: 4-character PDB ID (e.g., "1TRZ")

    Returns:
        List of catalytic residues with chain, position, name, and role
    """
    if not pdb_id or len(pdb_id) != 4:
        return []

    pdb_id = pdb_id.upper()
    url = f"{MCSA_API_BASE}/residues/?pdb_id={pdb_id}&format=json"

    try:
        async with httpx.AsyncClient(timeout=10.0) as client:
            response = await client.get(url)

            if response.status_code != 200:
                return []

            data = response.json()

            if not data:
                return []

            residues = []
            seen = set()  # Deduplicate by chain+residue

            for entry in data:
                # Extract residue info from M-CSA response
                chains = entry.get("residue_chains", {})
                chain = chains.get("chain_name", "A")
                resid = chains.get("resid")
                code = chains.get("code", "UNK")
                role = entry.get("roles_summary", "")

                if resid is None:
                    continue

                key = f"{chain}{resid}"
                if key in seen:
                    continue
                seen.add(key)

                residues.append(CatalyticResidue(
                    chain=chain,
                    residue=int(resid),
                    name=code,
                    role=role if role else None,
                    confidence=1.0,
                    source="mcsa"
                ))

            return residues

    except Exception as e:
        print(f"[M-CSA] Query failed: {e}")
        return []
