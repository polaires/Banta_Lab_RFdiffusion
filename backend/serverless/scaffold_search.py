"""
Scaffold Search — Auto-discover PDB scaffolds for metal-ligand binding sites.

When a de novo design query specifies both a metal and ligand but no PDB ID,
this module searches RCSB for structures containing the ligand with the same
or HSAB-compatible metals, validates spatial proximity, ranks candidates with
a composite score, and recommends scaffolding if a good match is found.

Falls back to de novo gracefully if no suitable scaffold is found.
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from metal_site_fetcher import query_metal_ligand_sites, find_metal_ligand_active_site
from scaffolding_workflow import normalize_ligand_code, is_metal_substitution_compatible
from metal_chemistry import (
    METAL_DATABASE,
    get_coordination_number_range,
    validate_coordination_chemistry,
)

logger = logging.getLogger(__name__)

# Score threshold: candidates scoring >= this trigger automatic scaffolding
SCAFFOLD_SEARCH_THRESHOLD = 50.0


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class ScaffoldCandidate:
    """A validated PDB hit with scoring fields."""
    pdb_id: str
    source_metal: str
    target_metal: str
    ligand_code: str
    needs_substitution: bool = False

    # Coordination info from spatial validation
    coordination_number: int = 0
    protein_donors: List[str] = field(default_factory=list)
    ligand_donors: List[str] = field(default_factory=list)
    ligand_distance: float = 0.0
    coordinating_atoms: List[Dict[str, Any]] = field(default_factory=list)

    # Resolution from PDB search
    resolution: float = 0.0

    # Component scores (0-max for each dimension)
    score_cn: float = 0.0         # CN match (max 25)
    score_hsab: float = 0.0       # HSAB donor compat (max 20)
    score_lig_donors: float = 0.0  # Ligand donor count (max 20)
    score_prot_donors: float = 0.0  # Protein donor count (max 10)
    score_distance: float = 0.0   # Bond distance quality (max 15)
    score_resolution: float = 0.0  # Resolution (max 10)

    total_score: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "pdb_id": self.pdb_id,
            "source_metal": self.source_metal,
            "target_metal": self.target_metal,
            "ligand_code": self.ligand_code,
            "needs_substitution": self.needs_substitution,
            "coordination_number": self.coordination_number,
            "protein_donors": self.protein_donors,
            "ligand_donors": self.ligand_donors,
            "ligand_distance": self.ligand_distance,
            "resolution": self.resolution,
            "score_cn": round(self.score_cn, 1),
            "score_hsab": round(self.score_hsab, 1),
            "score_lig_donors": round(self.score_lig_donors, 1),
            "score_prot_donors": round(self.score_prot_donors, 1),
            "score_distance": round(self.score_distance, 1),
            "score_resolution": round(self.score_resolution, 1),
            "total_score": round(self.total_score, 1),
        }


@dataclass
class ScaffoldSearchResult:
    """Aggregated result from the scaffold search pipeline."""
    searched: bool = False
    query_metal: str = ""
    query_ligand: str = ""
    ligand_code: str = ""
    num_pdb_hits: int = 0
    num_validated: int = 0
    candidates: List[ScaffoldCandidate] = field(default_factory=list)
    best_candidate: Optional[ScaffoldCandidate] = None
    recommended_action: str = "de_novo"  # "scaffold" or "de_novo"
    reason: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "searched": self.searched,
            "query_metal": self.query_metal,
            "query_ligand": self.query_ligand,
            "ligand_code": self.ligand_code,
            "num_pdb_hits": self.num_pdb_hits,
            "num_validated": self.num_validated,
            "candidates": [c.to_dict() for c in self.candidates],
            "best_candidate": self.best_candidate.to_dict() if self.best_candidate else None,
            "recommended_action": self.recommended_action,
            "reason": self.reason,
        }


# =============================================================================
# Main Entry Point
# =============================================================================

def search_scaffold_candidates(
    metal: str,
    ligand_name: str,
    ligand_code: Optional[str] = None,
    ligand_smiles: Optional[str] = None,
    resolution_max: float = 3.0,
    limit: int = 10,
) -> ScaffoldSearchResult:
    """
    Search PDB for scaffold candidates containing a metal-ligand complex.

    Orchestrates the full search -> validate -> rank pipeline.

    Args:
        metal: Target metal symbol (e.g., "TB")
        ligand_name: Human-readable ligand name (e.g., "citrate")
        ligand_code: Optional pre-resolved PDB 3-letter code
        ligand_smiles: Optional SMILES for scoring
        resolution_max: Max resolution filter for PDB search
        limit: Max PDB hits to validate

    Returns:
        ScaffoldSearchResult with ranked candidates and recommendation
    """
    result = ScaffoldSearchResult(
        searched=True,
        query_metal=metal.upper(),
        query_ligand=ligand_name,
    )

    # Step 1: Resolve ligand code
    try:
        code = ligand_code or _resolve_ligand_code(ligand_name)
    except Exception as e:
        logger.warning(f"Ligand code resolution failed: {e}")
        result.reason = f"Could not resolve ligand code: {e}"
        return result

    result.ligand_code = code

    # Step 2: Search PDB for hits (exact + compatible metals)
    try:
        hits = _search_pdb_hits(
            ligand_code=code,
            target_metal=metal.upper(),
            resolution_max=resolution_max,
            limit=limit,
        )
    except Exception as e:
        logger.warning(f"PDB search failed: {e}")
        result.reason = f"PDB search error: {e}"
        return result

    result.num_pdb_hits = len(hits)

    if not hits:
        result.reason = f"No PDB structures found with {code} and compatible metals"
        result.recommended_action = "de_novo"
        return result

    # Step 3: Validate each hit (spatial proximity check)
    candidates = []
    for hit in hits:
        try:
            candidate = _validate_candidate(
                pdb_id=hit["pdb_id"],
                source_metal=hit["metal"],
                target_metal=metal.upper(),
                ligand_code=code,
            )
            if candidate:
                candidate.resolution = hit.get("resolution", 0.0)
                candidates.append(candidate)
        except Exception as e:
            logger.debug(f"Validation failed for {hit['pdb_id']}: {e}")
            continue

    result.num_validated = len(candidates)

    if not candidates:
        result.reason = "No PDB hits passed spatial proximity validation"
        result.recommended_action = "de_novo"
        return result

    # Step 4: Score and rank
    ranked = rank_candidates(candidates, metal.upper(), ligand_smiles)
    result.candidates = ranked
    result.best_candidate = ranked[0] if ranked else None

    # Step 5: Decide action
    if result.best_candidate and result.best_candidate.total_score >= SCAFFOLD_SEARCH_THRESHOLD:
        result.recommended_action = "scaffold"
        sub = " (metal substitution)" if result.best_candidate.needs_substitution else ""
        result.reason = (
            f"Found {result.best_candidate.pdb_id} with score "
            f"{result.best_candidate.total_score:.0f}/100{sub}"
        )
    else:
        result.recommended_action = "de_novo"
        best_score = result.best_candidate.total_score if result.best_candidate else 0
        result.reason = (
            f"Best candidate scored {best_score:.0f}/100, "
            f"below threshold {SCAFFOLD_SEARCH_THRESHOLD:.0f}"
        )

    return result


# =============================================================================
# Internal Functions
# =============================================================================

def _resolve_ligand_code(ligand_name: str) -> str:
    """Convert human-readable ligand name to PDB 3-letter code."""
    return normalize_ligand_code(ligand_name)


def _get_hsab_compatible_metals(target_metal: str) -> List[str]:
    """Return all METAL_DATABASE metals that are substitution-compatible with target."""
    target_metal = target_metal.upper()
    compatible = []
    for m in METAL_DATABASE:
        if m == target_metal:
            continue
        if is_metal_substitution_compatible(m, target_metal):
            compatible.append(m)
    return compatible


def _search_pdb_hits(
    ligand_code: str,
    target_metal: str,
    resolution_max: float = 3.0,
    limit: int = 10,
) -> List[Dict[str, Any]]:
    """
    Two-phase PDB search: exact metal+ligand, then compatible metal+ligand.

    Phase 1: Search for structures with the exact target metal + ligand.
    Phase 2: Search for structures with ALL HSAB-compatible metals + ligand.

    Collects hits from all compatible metals before applying the limit,
    so high-quality hits from less-common metals aren't crowded out by
    abundant hits from the first compatible metal queried.

    Returns combined, deduplicated list of hits (up to limit).
    """
    seen_pdb_ids = set()
    all_hits = []

    # Phase 1: Exact metal match
    exact_hits = query_metal_ligand_sites(
        metal=target_metal,
        ligand=ligand_code,
        resolution_max=resolution_max,
        limit=limit,
    )
    for hit in exact_hits:
        pdb_id = hit["pdb_id"]
        if pdb_id not in seen_pdb_ids:
            seen_pdb_ids.add(pdb_id)
            all_hits.append(hit)

    # Phase 2: Query ALL compatible metals (don't stop early at limit)
    # Each metal gets a per-metal cap to avoid one metal dominating,
    # but we collect from all metals before applying the final limit.
    compatible_metals = _get_hsab_compatible_metals(target_metal)
    per_metal_cap = max(limit * 3, 30)  # Cast wide net; scoring picks the best

    for compat_metal in compatible_metals:
        compat_hits = query_metal_ligand_sites(
            metal=compat_metal,
            ligand=ligand_code,
            resolution_max=resolution_max,
            limit=per_metal_cap,
        )
        for hit in compat_hits:
            pdb_id = hit["pdb_id"]
            if pdb_id not in seen_pdb_ids:
                seen_pdb_ids.add(pdb_id)
                all_hits.append(hit)

    # Apply final limit — exact metal hits come first (already at front),
    # then compatible metal hits in query order
    # Don't truncate here — let validation + scoring filter down to the best.
    # The caller applies its own limit after ranking.
    return all_hits


def _validate_candidate(
    pdb_id: str,
    source_metal: str,
    target_metal: str,
    ligand_code: str,
    cutoff: float = 3.5,
) -> Optional[ScaffoldCandidate]:
    """
    Spatial proximity check: verify metal directly coordinates ligand.

    Returns None if metal-ligand distance exceeds cutoff.
    """
    site = find_metal_ligand_active_site(
        pdb_id=pdb_id,
        metal=source_metal,
        ligand=ligand_code,
        cutoff=cutoff,
    )

    if site is None:
        return None

    return ScaffoldCandidate(
        pdb_id=pdb_id,
        source_metal=source_metal.upper(),
        target_metal=target_metal.upper(),
        ligand_code=ligand_code,
        needs_substitution=(source_metal.upper() != target_metal.upper()),
        coordination_number=site.get("coordination_number", 0),
        protein_donors=site.get("protein_donors", []),
        ligand_donors=site.get("ligand_donors", []),
        ligand_distance=site.get("ligand_distance", 0.0),
        coordinating_atoms=site.get("coordinating_atoms", []),
    )


def _score_candidate(
    candidate: ScaffoldCandidate,
    target_metal: str,
    ligand_smiles: Optional[str] = None,
) -> float:
    """
    Composite score (0-100) across 6 dimensions:

    1. CN match to target (25 pts) -- actual CN vs expected range
    2. HSAB donor compatibility (20 pts) -- validate_coordination_chemistry()
    3. Ligand donor count (20 pts) -- more ligand donors = deeper chelation
    4. Protein donor count (10 pts) -- protein coordination environment
    5. Bond distance quality (15 pts) -- avg distance vs optimal range
    6. Resolution (10 pts) -- lower resolution = better structure
    """
    target_metal = target_metal.upper()

    # --- 1. CN match (25 pts) ---
    try:
        ox_state = METAL_DATABASE.get(target_metal, {}).get("default_oxidation", 2)
        min_cn, max_cn = get_coordination_number_range(target_metal, ox_state)
        actual_cn = candidate.coordination_number

        if min_cn <= actual_cn <= max_cn:
            candidate.score_cn = 25.0
        elif actual_cn < min_cn:
            gap = min_cn - actual_cn
            candidate.score_cn = max(0.0, 25.0 - gap * 5.0)
        else:
            gap = actual_cn - max_cn
            candidate.score_cn = max(0.0, 25.0 - gap * 3.0)
    except (ValueError, KeyError):
        candidate.score_cn = 10.0  # Unknown metal, give partial credit

    # --- 2. HSAB donor compatibility (20 pts) ---
    try:
        # Extract residue names from protein donors (format: "A123:GLU OE1@2.35A")
        residue_names = []
        for donor_str in candidate.protein_donors:
            parts = donor_str.split(":")
            if len(parts) >= 2:
                res_part = parts[1].split()[0]  # "GLU"
                residue_names.append(res_part)

        if residue_names:
            ox_state = METAL_DATABASE.get(target_metal, {}).get("default_oxidation", 2)
            validation = validate_coordination_chemistry(
                target_metal, residue_names, ox_state
            )
            if validation.get("valid") and validation.get("hsab_compatible"):
                candidate.score_hsab = 20.0
            elif validation.get("hsab_compatible"):
                candidate.score_hsab = 15.0
            else:
                n_warnings = len(validation.get("warnings", []))
                candidate.score_hsab = max(0.0, 15.0 - n_warnings * 3.0)
        else:
            candidate.score_hsab = 5.0  # No protein donors to check
    except Exception:
        candidate.score_hsab = 5.0

    # --- 3. Ligand donor count (20 pts) -- rewards richer ligand coordination ---
    # More ligand donors = ligand is deeply chelating the metal (better template).
    # Scaled: 1->4, 2->8, 3->12, 4->15, 5->17, 6+->20
    n_lig = len(candidate.ligand_donors)
    if n_lig >= 6:
        candidate.score_lig_donors = 20.0
    elif n_lig >= 4:
        candidate.score_lig_donors = 15.0
    elif n_lig == 3:
        candidate.score_lig_donors = 12.0
    elif n_lig == 2:
        candidate.score_lig_donors = 8.0
    elif n_lig == 1:
        candidate.score_lig_donors = 4.0
    else:
        candidate.score_lig_donors = 0.0

    # --- 4. Protein donor count (10 pts) -- max at 4+ ---
    # Protein donors matter but less than ligand coordination for metal-ligand scaffolds.
    n_prot = len(candidate.protein_donors)
    if n_prot >= 4:
        candidate.score_prot_donors = 10.0
    elif n_prot == 3:
        candidate.score_prot_donors = 8.0
    elif n_prot == 2:
        candidate.score_prot_donors = 5.0
    elif n_prot == 1:
        candidate.score_prot_donors = 2.0
    else:
        candidate.score_prot_donors = 0.0

    # --- 5. Bond distance quality (15 pts) ---
    try:
        ox_state = METAL_DATABASE.get(target_metal, {}).get("default_oxidation", 2)
        bond_dists = METAL_DATABASE.get(target_metal, {}).get("bond_distances", {}).get(ox_state, {})

        if bond_dists and candidate.coordinating_atoms:
            distances = [a["distance"] for a in candidate.coordinating_atoms if "distance" in a]
            if distances:
                avg_dist = sum(distances) / len(distances)
                # Get optimal range (use O as reference since most common)
                opt_range = bond_dists.get("O", bond_dists.get("N", (2.0, 2.8)))
                opt_mid = (opt_range[0] + opt_range[1]) / 2
                deviation = abs(avg_dist - opt_mid)
                # 15 pts at 0 deviation, 0 at 1.5A deviation
                candidate.score_distance = max(0.0, 15.0 * (1.0 - deviation / 1.5))
            else:
                candidate.score_distance = 7.5
        else:
            candidate.score_distance = 7.5
    except Exception:
        candidate.score_distance = 7.5

    # --- 6. Resolution (10 pts) -- lower is better, clamped 0-5A ---
    res = candidate.resolution
    if res <= 0:
        candidate.score_resolution = 5.0  # Unknown resolution
    elif res <= 1.5:
        candidate.score_resolution = 10.0
    elif res <= 3.0:
        # Linear scale: 1.5A -> 10, 3.0A -> 0
        candidate.score_resolution = max(0.0, 10.0 * (3.0 - res) / 1.5)
    else:
        candidate.score_resolution = 0.0

    # Total
    candidate.total_score = (
        candidate.score_cn
        + candidate.score_hsab
        + candidate.score_lig_donors
        + candidate.score_prot_donors
        + candidate.score_distance
        + candidate.score_resolution
    )

    return candidate.total_score


def rank_candidates(
    candidates: List[ScaffoldCandidate],
    target_metal: str,
    ligand_smiles: Optional[str] = None,
) -> List[ScaffoldCandidate]:
    """Score and sort candidates descending by total_score."""
    for c in candidates:
        _score_candidate(c, target_metal, ligand_smiles)
    candidates.sort(key=lambda c: c.total_score, reverse=True)
    return candidates
