"""Extracted handler modules for the RunPod serverless dispatcher.

Each module groups related API task handlers. handler.py imports from here
to keep the dispatcher routing thin while functions live in focused files.
"""

from handlers.history import handle_save_design_history, handle_check_lessons
from handlers.utility import handle_download_checkpoints, handle_delete_file
from handlers.esm3 import handle_esm3_score, handle_esm3_generate, handle_esm3_embed
from handlers.ligand import handle_resolve_ligand, handle_analyze_ligand_features
from handlers.analysis import (
    handle_validate_design,
    handle_binding_eval,
    handle_fastrelax,
    handle_detect_hotspots,
    handle_analyze_conservation,
    handle_interaction_analysis,
)
from handlers.orchestration import (
    handle_ai_parse,
    handle_scaffold_search,
    handle_scout_filter,
)
from handlers.design_analysis import handle_analyze_design

__all__ = [
    # History & learning
    "handle_save_design_history",
    "handle_check_lessons",
    # Utility
    "handle_download_checkpoints",
    "handle_delete_file",
    # ESM3
    "handle_esm3_score",
    "handle_esm3_generate",
    "handle_esm3_embed",
    # Ligand
    "handle_resolve_ligand",
    "handle_analyze_ligand_features",
    # Analysis & validation
    "handle_validate_design",
    "handle_binding_eval",
    "handle_fastrelax",
    "handle_detect_hotspots",
    "handle_analyze_conservation",
    "handle_interaction_analysis",
    # Orchestration
    "handle_ai_parse",
    "handle_scaffold_search",
    "handle_scout_filter",
    # Design analysis
    "handle_analyze_design",
]
