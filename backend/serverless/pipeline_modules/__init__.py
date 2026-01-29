"""
Pipeline Modules — Composable steps for AI-driven protein design.

Each module implements the PipelineStep protocol and can be composed
by the AI assistant (via Python) or WorkflowRunner (via JSON spec).

Usage:
    from pipeline_modules import BackboneGenerator, MODULE_REGISTRY

    # Direct Python composition
    gen = BackboneGenerator()
    ctx = gen.run(ctx)

    # Registry lookup (for JSON workflow specs)
    module_cls = MODULE_REGISTRY["backbone_generator"]
    step = module_cls()
"""

from pipeline_modules.backbone_generator import BackboneGenerator
from pipeline_modules.sequence_designer import SequenceDesigner
from pipeline_modules.structure_predictor import StructurePredictor
from pipeline_modules.intent_parser import IntentParser
from pipeline_modules.design_configurator import DesignConfigurator
from pipeline_modules.scaffolder import Scaffolder
from pipeline_modules.analyzer import Analyzer
from pipeline_modules.reporter import Reporter

MODULE_REGISTRY = {
    "intent_parser": IntentParser,
    "design_configurator": DesignConfigurator,
    "scaffolder": Scaffolder,
    "backbone_generator": BackboneGenerator,
    "sequence_designer": SequenceDesigner,
    "structure_predictor": StructurePredictor,
    "analyzer": Analyzer,
    "reporter": Reporter,
}

__all__ = [
    "IntentParser",
    "DesignConfigurator",
    "Scaffolder",
    "BackboneGenerator",
    "SequenceDesigner",
    "StructurePredictor",
    "Analyzer",
    "Reporter",
    "MODULE_REGISTRY",
]
