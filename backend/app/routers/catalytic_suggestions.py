"""
Catalytic residue suggestions API endpoint.
Queries M-CSA first, falls back to P2Rank for unknown structures.
"""

import re
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Optional, Literal

from app.services.mcsa_client import query_mcsa_by_pdb, CatalyticResidue
from app.services.prankweb_client import query_prankweb, PredictedResidue


router = APIRouter(prefix="/api", tags=["catalytic"])


class SuggestionRequest(BaseModel):
    """Request for catalytic residue suggestions."""
    pdb_id: Optional[str] = None
    pdb_content: str


class SuggestionResidue(BaseModel):
    """A suggested catalytic residue."""
    chain: str
    residue: int
    name: str
    role: Optional[str] = None
    confidence: float
    source: Literal["mcsa", "p2rank"]


class SuggestionResponse(BaseModel):
    """Response containing catalytic residue suggestions."""
    source: Literal["mcsa", "p2rank", "none"]
    residues: List[SuggestionResidue]


def extract_pdb_id_from_content(pdb_content: str) -> Optional[str]:
    """Extract PDB ID from HEADER line of PDB file."""
    if not pdb_content:
        return None

    # Look for HEADER line with PDB ID
    # Format: HEADER    HYDROLASE                               01-JAN-00   1ABC
    for line in pdb_content.split('\n')[:20]:
        if line.startswith('HEADER'):
            # PDB ID is typically at positions 62-66
            if len(line) >= 66:
                pdb_id = line[62:66].strip()
                if len(pdb_id) == 4 and pdb_id.isalnum():
                    return pdb_id.upper()
            # Also try regex for flexibility
            match = re.search(r'\b([0-9][A-Za-z0-9]{3})\s*$', line)
            if match:
                return match.group(1).upper()

    return None


@router.post("/catalytic-suggestions", response_model=SuggestionResponse)
async def get_catalytic_suggestions(request: SuggestionRequest) -> SuggestionResponse:
    """
    Get catalytic residue suggestions for a protein structure.

    First queries M-CSA for curated catalytic residues.
    If no results, falls back to P2Rank binding site prediction.
    """
    # Determine PDB ID
    pdb_id = request.pdb_id
    if not pdb_id:
        pdb_id = extract_pdb_id_from_content(request.pdb_content)

    # Try M-CSA first if we have a PDB ID
    if pdb_id:
        mcsa_residues = await query_mcsa_by_pdb(pdb_id)

        if mcsa_residues:
            return SuggestionResponse(
                source="mcsa",
                residues=[
                    SuggestionResidue(
                        chain=r.chain,
                        residue=r.residue,
                        name=r.name,
                        role=r.role,
                        confidence=r.confidence,
                        source="mcsa"
                    )
                    for r in mcsa_residues
                ]
            )

    # Fallback to P2Rank
    p2rank_residues = await query_prankweb(request.pdb_content)

    if p2rank_residues:
        return SuggestionResponse(
            source="p2rank",
            residues=[
                SuggestionResidue(
                    chain=r.chain,
                    residue=r.residue,
                    name=r.name,
                    role=None,
                    confidence=r.confidence,
                    source="p2rank"
                )
                for r in p2rank_residues
            ]
        )

    # No suggestions found
    return SuggestionResponse(source="none", residues=[])
