"""
AI Protein Engineering Module

This package contains the enhanced hybrid AI system for protein design:
- Chemistry rule base from research/literature
- LLM integration for parameter refinement
- AI-driven task detection and parameter generation
- Design evaluation and improvement suggestions
"""

from .recommender import (
    MetalChemistryRules,
    DesignStrategyRules,
    AIRecommender,
    generate_rfd3_parameters,
)

from .ai_engine import (
    AIEngine,
    TaskType,
    FormFieldConfig,
    AIResponse,
    HandbookLoader,
    analyze_user_request,
)

__all__ = [
    # Original recommender
    "MetalChemistryRules",
    "DesignStrategyRules",
    "AIRecommender",
    "generate_rfd3_parameters",
    # New AI engine
    "AIEngine",
    "TaskType",
    "FormFieldConfig",
    "AIResponse",
    "HandbookLoader",
    "analyze_user_request",
]
