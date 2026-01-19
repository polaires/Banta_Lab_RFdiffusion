"""
PrankWeb/P2Rank API client for binding site prediction.
Used as fallback when M-CSA has no data for a structure.
"""

import httpx
import asyncio
from typing import List, Optional
from pydantic import BaseModel


class PredictedResidue(BaseModel):
    """A predicted binding site residue from P2Rank."""
    chain: str
    residue: int
    name: str
    confidence: float
    source: str = "p2rank"


PRANKWEB_API = "https://prankweb.cz/api/v2"


async def query_prankweb(pdb_content: str, max_wait: int = 30) -> List[PredictedResidue]:
    """
    Submit structure to PrankWeb for binding site prediction.

    Args:
        pdb_content: PDB file content as string
        max_wait: Maximum seconds to wait for results

    Returns:
        List of predicted binding site residues from top pocket
    """
    if not pdb_content or len(pdb_content) < 100:
        return []

    try:
        async with httpx.AsyncClient(timeout=60.0) as client:
            # Submit prediction job
            submit_response = await client.post(
                f"{PRANKWEB_API}/predict/structure",
                json={"structure": pdb_content},
                headers={"Content-Type": "application/json"}
            )

            if submit_response.status_code not in (200, 201, 202):
                print(f"[PrankWeb] Submit failed: {submit_response.status_code}")
                return []

            job_data = submit_response.json()
            job_id = job_data.get("id") or job_data.get("taskId")

            if not job_id:
                # Some versions return results directly
                return _parse_prankweb_results(job_data)

            # Poll for results
            for _ in range(max_wait):
                await asyncio.sleep(1)

                status_response = await client.get(f"{PRANKWEB_API}/predict/{job_id}")

                if status_response.status_code != 200:
                    continue

                status_data = status_response.json()
                status = status_data.get("status", "").lower()

                if status in ("completed", "finished", "done"):
                    return _parse_prankweb_results(status_data)
                elif status in ("failed", "error"):
                    print(f"[PrankWeb] Job failed: {status_data.get('error')}")
                    return []

            print("[PrankWeb] Timeout waiting for results")
            return []

    except Exception as e:
        print(f"[PrankWeb] Query failed: {e}")
        return []


def _parse_prankweb_results(data: dict) -> List[PredictedResidue]:
    """Parse PrankWeb response into residue list."""
    residues = []
    seen = set()

    # Get pockets from response
    pockets = data.get("pockets", [])
    if not pockets:
        pockets = data.get("predictions", {}).get("pockets", [])

    if not pockets:
        return []

    # Use top pocket (highest score)
    top_pocket = pockets[0]
    pocket_residues = top_pocket.get("residues", [])
    pocket_score = top_pocket.get("score", 0.5)

    for res in pocket_residues:
        chain = res.get("chain", "A")
        resnum = res.get("residue_number") or res.get("resNum")
        resname = res.get("residue_name") or res.get("resName", "UNK")
        score = res.get("score", pocket_score)

        if resnum is None:
            continue

        key = f"{chain}{resnum}"
        if key in seen:
            continue
        seen.add(key)

        residues.append(PredictedResidue(
            chain=chain,
            residue=int(resnum),
            name=resname,
            confidence=float(score) if score else 0.5,
            source="p2rank"
        ))

    # Sort by confidence, take top 10
    residues.sort(key=lambda r: r.confidence, reverse=True)
    return residues[:10]
