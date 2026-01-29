"""
Backbone Generator Module (M4)

Generates protein backbones via RFD3.
Wraps inference_utils.run_rfd3_inference() — does NOT rewrite it.

Reads: rfd3_config (or params) from context
Writes: backbone_pdbs to context
"""

import logging
from typing import Any, Dict, List

from pipeline_types import StepContext

logger = logging.getLogger(__name__)

BATCH_SIZE = 10  # Process RFD3 in batches to manage memory


class BackboneGenerator:
    """Generate protein backbones using RFD3."""

    name = "backbone_generator"
    description = "Generate protein backbone structures using RFdiffusion3"
    input_keys = ["rfd3_config"]
    output_keys = ["backbone_pdbs"]
    optional_keys = ["scaffold_result"]

    def __init__(self, **kwargs):
        """Accept optional param overrides."""
        self.overrides = kwargs

    def run(self, context: StepContext) -> StepContext:
        """Execute RFD3 backbone generation.

        If SweepStrategy populated sweep_configs in metadata, iterates over
        each config. Otherwise runs a single config.
        """
        backend = context.backend
        if backend is None:
            raise RuntimeError("No inference backend configured on context")

        sweep_configs = context.metadata.get("sweep_configs")
        if sweep_configs:
            return self._run_sweep(context, backend, sweep_configs)
        return self._run_single(context, backend)

    def _run_single(self, context: StepContext, backend: Any) -> StepContext:
        """Run a single RFD3 config."""
        api_params = self._build_params(context)
        num_designs = api_params.get("num_designs", context.get_param("num_designs", 4))
        api_params["num_designs"] = num_designs

        logger.info(f"Generating {num_designs} backbone(s) via RFD3")
        all_pdbs = self._generate_batch(backend, api_params, num_designs)

        logger.info(f"Generated {len(all_pdbs)} backbone(s)")
        context.backbone_pdbs = all_pdbs
        context.metadata["backbone_generator"] = {
            "num_generated": len(all_pdbs),
            "num_requested": num_designs,
            "batch_size": BATCH_SIZE,
        }
        return context

    def _run_sweep(self, context: StepContext, backend: Any, sweep_configs: List[Dict[str, Any]]) -> StepContext:
        """Iterate over multiple RFD3 configs from SweepStrategy."""
        designs_per_config = context.metadata.get("sweep_designs_per_config", 10)
        base_params = self._build_params(context)

        logger.info(f"Sweep: running {len(sweep_configs)} config(s), "
                    f"{designs_per_config} designs each")

        all_pdbs: List[str] = []
        sweep_results: List[Dict[str, Any]] = []

        for i, cfg_overrides in enumerate(sweep_configs):
            params = {**base_params, **cfg_overrides, "num_designs": designs_per_config}
            logger.info(f"Sweep config {i + 1}/{len(sweep_configs)}: "
                       f"{cfg_overrides}")

            pdbs = self._generate_batch(backend, params, designs_per_config)
            all_pdbs.extend(pdbs)
            sweep_results.append({
                "config_index": i,
                "overrides": cfg_overrides,
                "num_generated": len(pdbs),
            })

        logger.info(f"Sweep complete: {len(all_pdbs)} total backbone(s) "
                    f"from {len(sweep_configs)} config(s)")
        context.backbone_pdbs = all_pdbs
        context.metadata["backbone_generator"] = {
            "num_generated": len(all_pdbs),
            "sweep_configs": len(sweep_configs),
            "designs_per_config": designs_per_config,
            "sweep_results": sweep_results,
        }
        return context

    def _generate_batch(self, backend: Any, api_params: Dict[str, Any], num_designs: int) -> List[str]:
        """Generate backbones in batches, return list of PDB strings."""
        all_pdbs: List[str] = []
        remaining = num_designs

        while remaining > 0:
            batch_size = min(remaining, BATCH_SIZE)
            batch_params = {**api_params, "num_designs": batch_size}

            result = backend.run_rfd3(batch_params)
            designs = result.get("designs", [])

            for design in designs:
                pdb = design.get("pdb_content", "")
                if pdb:
                    input_pdb = api_params.get("pdb_content")
                    if input_pdb:
                        pdb = _merge_hetatm(pdb, input_pdb)
                    all_pdbs.append(pdb)

            remaining -= batch_size

        return all_pdbs

    def _build_params(self, context: StepContext) -> Dict[str, Any]:
        """Build RFD3 API params from context."""
        params: Dict[str, Any] = {}

        # Prefer explicit rfd3_config on context
        if context.rfd3_config:
            if isinstance(context.rfd3_config, dict):
                params = dict(context.rfd3_config)
            elif hasattr(context.rfd3_config, "to_api_params"):
                params = context.rfd3_config.to_api_params()

        # Merge scaffold result into config if present
        if context.scaffold_result:
            scaffold = context.scaffold_result
            if isinstance(scaffold, dict):
                if "pdb_content" in scaffold and "pdb_content" not in params:
                    params["pdb_content"] = scaffold["pdb_content"]
                if "contig" in scaffold and "contig" not in params:
                    params["contig"] = scaffold["contig"]

        # Apply param-level overrides
        for key in ("num_designs", "contig", "pdb_content", "cfg_scale",
                     "step_scale", "num_timesteps", "seed"):
            val = context.get_param(key)
            if val is not None and key not in params:
                params[key] = val

        # Apply constructor overrides last
        params.update(self.overrides)

        return params


def _merge_hetatm(design_pdb: str, input_pdb: str) -> str:
    """
    Merge HETATM records from the input PDB into the design PDB.

    RFD3 outputs only ATOM records; ligands/metals need to be copied
    from the input PDB so downstream steps (MPNN, analysis) see them.
    """
    hetatm_lines = []
    for line in input_pdb.splitlines():
        if line.startswith("HETATM"):
            hetatm_lines.append(line)

    if not hetatm_lines:
        return design_pdb

    # Insert HETATM before END/TER at the end of the design PDB
    lines = design_pdb.splitlines()
    insert_idx = len(lines)
    for i in range(len(lines) - 1, -1, -1):
        if lines[i].startswith("END"):
            insert_idx = i
            break

    merged = lines[:insert_idx] + hetatm_lines + lines[insert_idx:]
    return "\n".join(merged)
