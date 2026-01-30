"""
Scaffold Searcher Pipeline Module

Auto-search PDB for scaffold candidates when a de novo query specifies
both a metal and ligand but no PDB ID. If a good scaffold is found,
reroutes the pipeline to scaffolding mode.

Reads: design_intent, resolved_ligand (optional)
Writes: scaffold_search_result to context
Side-effect: may mutate design_intent to scaffold mode
"""

import logging
from typing import Any, Dict, List, Optional

from pipeline_types import StepContext
from scaffold_search import (
    search_scaffold_candidates,
    SCAFFOLD_SEARCH_THRESHOLD,
)

logger = logging.getLogger(__name__)


class ScaffoldSearcher:
    """Auto-search PDB for scaffold candidates when metal+ligand specified."""

    name = "scaffold_searcher"
    description = "Auto-search PDB for scaffold candidates when metal+ligand specified"
    input_keys = ["design_intent"]
    output_keys = ["scaffold_search_result"]
    optional_keys = ["resolved_ligand"]

    def __init__(
        self,
        min_score: float = SCAFFOLD_SEARCH_THRESHOLD,
        resolution_max: float = 3.0,
        limit: int = 10,
        auto_route: bool = True,
        **kwargs,
    ):
        self.min_score = min_score
        self.resolution_max = resolution_max
        self.limit = limit
        self.auto_route = auto_route

    def run(self, context: StepContext) -> StepContext:
        """Search PDB for scaffolds and optionally reroute to scaffold mode."""
        intent = context.design_intent
        if not intent:
            logger.debug("scaffold_searcher: no design_intent, skipping")
            return context

        # Guard: skip if already in scaffold mode
        if getattr(intent, "design_mode", "de_novo") == "scaffold":
            logger.debug("scaffold_searcher: already in scaffold mode, skipping")
            return context

        # Guard: need both metal and ligand
        metal = getattr(intent, "metal_type", None)
        ligand = getattr(intent, "ligand_name", None)
        if not metal or not ligand:
            logger.debug("scaffold_searcher: no metal+ligand pair, skipping")
            return context

        # Get ligand code and SMILES from resolved_ligand if available
        ligand_code = None
        ligand_smiles = None
        resolved = context.resolved_ligand
        if resolved and getattr(resolved, "resolved", False):
            ligand_code = getattr(resolved, "ligand_code", None)
            ligand_smiles = getattr(resolved, "smiles", None)

        logger.info(f"scaffold_searcher: searching PDB for {metal}+{ligand} scaffolds")

        try:
            search_result = search_scaffold_candidates(
                metal=metal,
                ligand_name=ligand,
                ligand_code=ligand_code,
                ligand_smiles=ligand_smiles,
                resolution_max=self.resolution_max,
                limit=self.limit,
            )
        except Exception as e:
            logger.warning(f"scaffold_searcher: search failed: {e}")
            context.metadata["scaffold_searcher"] = {
                "status": "error",
                "error": str(e),
            }
            return context

        # Store result in context
        context.scaffold_search_result = search_result.to_dict()

        # Auto-route to scaffolding if score meets threshold
        best = search_result.best_candidate
        if (
            self.auto_route
            and best is not None
            and best.total_score >= self.min_score
        ):
            logger.info(
                f"scaffold_searcher: rerouting to scaffold mode "
                f"(PDB={best.pdb_id}, score={best.total_score:.0f})"
            )
            intent.design_mode = "scaffold"
            intent.source_pdb_id = best.pdb_id
            context.metadata["scaffold_searcher"] = {
                "status": "rerouted",
                "pdb_id": best.pdb_id,
                "score": round(best.total_score, 1),
                "needs_substitution": best.needs_substitution,
            }
        else:
            action = search_result.recommended_action
            context.metadata["scaffold_searcher"] = {
                "status": "searched",
                "recommended_action": action,
                "reason": search_result.reason,
                "num_candidates": search_result.num_validated,
            }

        return context
