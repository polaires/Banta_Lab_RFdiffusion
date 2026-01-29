"""
Scaffolder Module (M3)

Extracts motifs/pockets from existing PDB structures for scaffold-based design.
Wraps ScaffoldingWorkflow.run() — does NOT rewrite it.

Reads: params (pdb_id, ligand, metal), design_intent from context
Writes: scaffold_result, rfd3_config updates to context
"""

import logging
from typing import Any, Dict, List, Optional

from pipeline_types import StepContext

logger = logging.getLogger(__name__)


class Scaffolder:
    """Extract motifs and scaffold information from PDB structures."""

    name = "scaffolder"
    description = "Extract binding motifs and scaffold from PDB structures for motif-based design"
    input_keys = ["design_intent"]
    output_keys = ["scaffold_result"]
    optional_keys = ["params.pdb_id", "params.ligand", "params.metal"]

    def __init__(self, **kwargs):
        self.overrides = kwargs

    def run(self, context: StepContext) -> StepContext:
        """Run scaffolding workflow."""
        intent = context.design_intent

        # Gather scaffolding parameters
        pdb_id = self.overrides.get("pdb_id") or context.get_param("pdb_id")
        ligand = self.overrides.get("ligand") or context.get_param("ligand")
        metal = self.overrides.get("metal") or context.get_param("metal")
        enzyme_class = self.overrides.get("enzyme_class") or context.get_param("enzyme_class")

        # Fall back to design intent
        if intent:
            pdb_id = pdb_id or getattr(intent, "source_pdb_id", None)
            ligand = ligand or getattr(intent, "ligand_name", None)
            metal = metal or getattr(intent, "metal_type", None)
            enzyme_class = enzyme_class or getattr(intent, "enzyme_class", None)

        if not pdb_id:
            logger.warning("No PDB ID for scaffolding — skipping")
            return context

        logger.info(f"Running scaffolding: PDB={pdb_id}, ligand={ligand}, metal={metal}")

        try:
            from scaffolding_workflow import ScaffoldingWorkflow
            workflow = ScaffoldingWorkflow()
            scaffold_result = workflow.run(
                pdb_id=pdb_id,
                ligand_name=ligand,
                metal_type=metal,
                enzyme_class=enzyme_class,
            )

            # Store as dict for serialization
            if hasattr(scaffold_result, "__dict__"):
                context.scaffold_result = {
                    "pdb_content": getattr(scaffold_result, "motif_pdb", None)
                                   or getattr(scaffold_result, "pdb_content", None),
                    "pocket_description": getattr(scaffold_result, "pocket_description", ""),
                    "metal_type": getattr(scaffold_result, "metal_type", metal),
                    "ligand_code": getattr(scaffold_result, "ligand_code", ligand),
                    "residues": getattr(scaffold_result, "residues", []),
                    "substitution_info": getattr(scaffold_result, "substitution_info", None),
                }
            else:
                context.scaffold_result = scaffold_result

            # Update RFD3 config with scaffold params
            self._update_rfd3_config(context, scaffold_result)

            logger.info(f"Scaffolding complete: {context.scaffold_result.get('pocket_description', 'OK')}")

        except ImportError:
            logger.error("ScaffoldingWorkflow not available")
            raise
        except Exception as e:
            logger.error(f"Scaffolding failed: {e}")
            raise

        context.metadata["scaffolder"] = {
            "pdb_id": pdb_id,
            "ligand": ligand,
            "metal": metal,
            "enzyme_class": enzyme_class,
        }
        return context

    def _update_rfd3_config(self, context: StepContext, scaffold_result: Any) -> None:
        """Update RFD3 config with scaffold-derived parameters."""
        try:
            from rfd3_config_generator import scaffold_to_rfd3_params
            rfd3_updates = scaffold_to_rfd3_params(scaffold_result)
            if context.rfd3_config and isinstance(context.rfd3_config, dict):
                context.rfd3_config.update(rfd3_updates)
            elif context.rfd3_config and hasattr(context.rfd3_config, "__dict__"):
                for k, v in rfd3_updates.items():
                    if hasattr(context.rfd3_config, k):
                        setattr(context.rfd3_config, k, v)
            else:
                context.rfd3_config = rfd3_updates
        except (ImportError, AttributeError):
            logger.debug("scaffold_to_rfd3_params not available — skipping config update")
