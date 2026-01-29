"""
Design Configurator Module (M2)

Converts DesignIntent + ResolvedLigand into RFD3 and MPNN configs.
Wraps design_types.infer_design_type(), design_rules.make_design_decisions(),
and RFD3ConfigGenerator.generate() — does NOT rewrite them.

Reads: design_intent, resolved_ligand from context
Writes: rfd3_config, mpnn_config to context
"""

import logging
from typing import Any, Dict, List, Optional

from pipeline_types import StepContext

logger = logging.getLogger(__name__)


class DesignConfigurator:
    """Generate RFD3 and MPNN configurations from design intent."""

    name = "design_configurator"
    description = "Generate RFD3 and MPNN configurations from parsed design intent"
    input_keys = ["design_intent"]
    output_keys = ["rfd3_config", "mpnn_config"]
    optional_keys = ["resolved_ligand"]

    def __init__(self, **kwargs):
        self.overrides = kwargs

    def run(self, context: StepContext) -> StepContext:
        """Generate configs from intent and ligand."""
        intent = context.design_intent
        if intent is None:
            raise ValueError("No design_intent in context — run intent_parser first")

        ligand = context.resolved_ligand

        # --- Import config generator ---
        from rfd3_config_generator import RFD3ConfigGenerator

        num_designs = self.overrides.get(
            "num_designs",
            context.get_param("num_designs", 4),
        )
        num_sequences = self.overrides.get(
            "num_sequences",
            context.get_param("num_sequences", 8),
        )

        generator = RFD3ConfigGenerator(
            default_num_designs=num_designs,
            default_num_sequences=num_sequences,
        )

        # --- Generate RFD3 config ---
        if ligand and hasattr(ligand, "resolved") and ligand.resolved:
            rfd3_config = generator.generate(intent, ligand, num_designs)
            mpnn_config = generator.generate_mpnn_config(intent, ligand, num_sequences)
            logger.info(f"Generated RFD3 config: cfg_scale={rfd3_config.cfg_scale}, "
                       f"step_scale={rfd3_config.step_scale}")
        else:
            # Fallback: generate basic config without resolved ligand
            logger.warning("No resolved ligand — generating basic config")
            rfd3_config = self._basic_rfd3_config(intent, num_designs)
            mpnn_config = self._basic_mpnn_config(num_sequences)

        # Apply overrides
        if "cfg_scale" in self.overrides:
            rfd3_config.cfg_scale = self.overrides["cfg_scale"]
        if "step_scale" in self.overrides:
            rfd3_config.step_scale = self.overrides["step_scale"]

        context.rfd3_config = rfd3_config
        context.mpnn_config = mpnn_config

        context.metadata["design_configurator"] = {
            "cfg_scale": rfd3_config.cfg_scale if hasattr(rfd3_config, "cfg_scale") else None,
            "num_designs": num_designs,
            "num_sequences": num_sequences,
            "has_ligand": ligand is not None and getattr(ligand, "resolved", False),
        }
        return context

    def _basic_rfd3_config(self, intent: Any, num_designs: int) -> Any:
        """Generate a basic RFD3 config without a resolved ligand."""
        from rfd3_config_generator import RFD3Config
        config = RFD3Config()
        config.num_designs = num_designs
        # Use intent chain length hints
        chain_min = getattr(intent, "chain_length_min", 80)
        chain_max = getattr(intent, "chain_length_max", 120)
        config.contig = f"/0,{chain_min}-{chain_max}"
        return config

    def _basic_mpnn_config(self, num_sequences: int) -> Any:
        """Generate a basic MPNN config."""
        from rfd3_config_generator import LigandMPNNConfig
        config = LigandMPNNConfig()
        config.num_sequences = num_sequences
        return config
