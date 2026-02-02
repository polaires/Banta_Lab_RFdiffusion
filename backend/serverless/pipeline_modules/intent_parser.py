"""
Intent Parser Module (M1)

Parses natural language design requests into structured DesignIntent.
Wraps NLDesignParser + LigandResolver — does NOT rewrite them.

Reads: params.query from context
Writes: design_intent, resolved_ligand to context
"""

import logging
from typing import Any, Dict, List, Optional

from pipeline_types import StepContext

logger = logging.getLogger(__name__)


class IntentParser:
    """Parse natural language into DesignIntent + resolve ligand."""

    name = "intent_parser"
    description = "Parse natural language design request into structured intent and resolve ligand"
    input_keys = ["params.query"]
    output_keys = ["design_intent", "resolved_ligand"]
    optional_keys = []

    def __init__(self, api_key: Optional[str] = None, api_url: Optional[str] = None, **kwargs):
        self.api_key = api_key
        self.api_url = api_url
        self.overrides = kwargs

    def run(self, context: StepContext) -> StepContext:
        """Parse NL query and resolve ligand."""
        query = context.get_param("query", "")
        if not query:
            raise ValueError("No 'query' parameter in context — set params.query")

        # --- Step 1: Parse NL query ---
        logger.info(f"Parsing design intent from: {query[:80]}...")
        from nl_design_parser import create_parser
        parser = create_parser(
            api_key=self.api_key or context.get_param("claude_api_key"),
            api_url=self.api_url or context.get_param("claude_api_url"),
            use_fallback=True,
        )
        intent = parser.parse(query)
        context.design_intent = intent
        logger.info(
            f"Parsed: metal={intent.metal_type}, ligand={intent.ligand_name}, "
            f"goal={intent.design_goal}, mode={'scaffold' if intent.is_scaffolding else 'de_novo'}"
        )

        # --- Step 2: Resolve ligand ---
        if intent.ligand_name or intent.metal_type:
            logger.info(f"Resolving ligand: {intent.ligand_name or 'metal-only'}")
            from ligand_resolver import LigandResolver
            resolver = LigandResolver()
            ligand = resolver.resolve(
                ligand_name=intent.ligand_name,
                metal_type=intent.metal_type,
                isomer_spec=getattr(intent, 'isomer_specification', None),
            )
            context.resolved_ligand = ligand
            if ligand.resolved:
                logger.info(f"Ligand resolved: {ligand.name} ({ligand.source})")
            else:
                logger.warning(f"Ligand not resolved: {ligand.warnings}")

        context.metadata["intent_parser"] = {
            "query": query,
            "metal": intent.metal_type,
            "ligand": intent.ligand_name,
            "confidence": intent.confidence,
            "is_scaffolding": intent.is_scaffolding,
        }
        return context
