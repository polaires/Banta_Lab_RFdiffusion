"""
Pipeline Types — Shared data structures for the modular pipeline architecture.

Defines the StepContext (shared state bag), PipelineStep protocol,
and result dataclasses used by all pipeline modules.
"""

import logging
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Protocol, runtime_checkable

logger = logging.getLogger(__name__)


# =============================================================================
# Result Dataclasses
# =============================================================================

@dataclass
class SequenceResult:
    """Result from a single sequence design."""
    sequence: str
    score: float = 0.0
    confidence: float = 0.0
    backbone_index: int = 0
    filename: str = ""
    pdb_content: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class PredictionResult:
    """Result from structure prediction (RF3 or ESMFold)."""
    sequence: str
    predicted_pdb: Optional[str] = None
    plddt: float = 0.0
    ptm: float = 0.0
    iptm: Optional[float] = None
    rmsd: Optional[float] = None
    passed: bool = False
    predictor: str = "rf3"  # "rf3" or "esmfold"
    backbone_index: int = 0
    sequence_index: int = 0
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class DesignCandidate:
    """A complete design candidate with backbone, sequence, and prediction."""
    backbone_pdb: str
    sequence: str
    predicted_pdb: Optional[str] = None
    plddt: float = 0.0
    ptm: float = 0.0
    rmsd: Optional[float] = None
    passed: bool = False
    analysis: Optional[Dict[str, Any]] = None
    rank: int = 0
    metadata: Dict[str, Any] = field(default_factory=dict)


# =============================================================================
# Inference Backend Protocol
# =============================================================================

@runtime_checkable
class InferenceBackend(Protocol):
    """Protocol for inference backends (in-process or HTTP)."""

    def run_rfd3(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Run RFD3 backbone generation."""
        ...

    def run_mpnn(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Run MPNN sequence design."""
        ...

    def run_rf3(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Run RF3 structure prediction."""
        ...


# =============================================================================
# Step Context — Shared state flowing through pipeline
# =============================================================================

@dataclass
class StepContext:
    """
    Mutable context that flows through pipeline steps.

    All data is string-based, matching actual inference patterns:
    - PDB content as strings
    - FASTA sequences as strings
    - Configs as dicts ready for API calls
    """

    # --- Core data ---
    # Design intent (from NL parsing)
    design_intent: Optional[Any] = None  # DesignIntent from nl_design_parser
    resolved_ligand: Optional[Any] = None  # ResolvedLigand from ligand_resolver

    # Configurations
    rfd3_config: Optional[Dict[str, Any]] = None
    mpnn_config: Optional[Dict[str, Any]] = None

    # Scaffold data
    scaffold_result: Optional[Dict[str, Any]] = None

    # Backbone PDB strings from RFD3
    backbone_pdbs: List[str] = field(default_factory=list)

    # Sequence design results
    sequences: List[SequenceResult] = field(default_factory=list)

    # Structure prediction results
    predictions: List[PredictionResult] = field(default_factory=list)

    # Analysis results
    analysis: Optional[Dict[str, Any]] = None

    # Report text
    report: Optional[str] = None

    # --- Infrastructure ---
    backend: Optional[Any] = None  # InferenceBackend instance
    working_dir: Optional[Path] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    # --- Params (any step can read these) ---
    params: Dict[str, Any] = field(default_factory=dict)

    def get_param(self, key: str, default: Any = None) -> Any:
        """Get a parameter from params dict with optional default."""
        return self.params.get(key, default)

    def to_dict(self) -> Dict[str, Any]:
        """Serialize context to dict for checkpointing / API response."""
        result: Dict[str, Any] = {
            "params": self.params,
            "metadata": self.metadata,
            "num_backbones": len(self.backbone_pdbs),
            "num_sequences": len(self.sequences),
            "num_predictions": len(self.predictions),
        }

        # Include design intent
        if self.design_intent and hasattr(self.design_intent, "to_dict"):
            result["design_intent"] = self.design_intent.to_dict()

        # Include configs
        if self.rfd3_config:
            result["rfd3_config"] = self.rfd3_config
        if self.mpnn_config:
            result["mpnn_config"] = self.mpnn_config

        # Include scaffold result
        if self.scaffold_result:
            result["scaffold_result"] = self.scaffold_result

        # Backbone PDbs — include count, truncate if large
        if self.backbone_pdbs:
            if len(self.backbone_pdbs) <= 10:
                result["backbone_pdbs"] = self.backbone_pdbs
            else:
                result["backbone_pdbs"] = self.backbone_pdbs[:10]
                result["backbone_pdbs_truncated"] = True

        # Sequences — lightweight summary
        if self.sequences:
            result["sequences"] = [
                {
                    "sequence": s.sequence,
                    "score": s.score,
                    "confidence": s.confidence,
                    "backbone_index": s.backbone_index,
                }
                for s in self.sequences
            ]

        # Predictions — summary
        if self.predictions:
            result["predictions"] = [
                {
                    "plddt": p.plddt,
                    "ptm": p.ptm,
                    "rmsd": p.rmsd,
                    "passed": p.passed,
                    "predictor": p.predictor,
                    "backbone_index": p.backbone_index,
                }
                for p in self.predictions
            ]

        # Analysis
        if self.analysis:
            result["analysis"] = self.analysis

        # Report
        if self.report:
            result["report"] = self.report

        return result


# =============================================================================
# Pipeline Step Protocol
# =============================================================================

@runtime_checkable
class PipelineStep(Protocol):
    """
    A composable pipeline step.

    Each module implements this protocol so the AI assistant or
    WorkflowRunner can compose them dynamically.
    """

    name: str
    description: str
    input_keys: List[str]
    output_keys: List[str]
    optional_keys: List[str]

    def run(self, context: StepContext) -> StepContext:
        """Execute step. Reads from context, writes results back."""
        ...
