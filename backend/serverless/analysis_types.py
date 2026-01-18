"""
Analysis Type Definitions

Provides structured types for design analysis results with clear status handling.
"""
from enum import Enum
from typing import Dict, Any, Optional
from dataclasses import dataclass


class AnalysisStatus(Enum):
    """Status of an analysis run."""
    SUCCESS = "success"
    NOT_APPLICABLE = "not_applicable"
    SKIPPED = "skipped"


class DesignType(Enum):
    """Types of protein designs."""
    MONOMER = "monomer"
    PROTEIN_DIMER = "protein_dimer"
    LIGAND_INTERFACE_DIMER = "ligand_interface_dimer"
    METAL_INTERFACE_DIMER = "metal_interface_dimer"
    METAL_LIGAND_INTERFACE_DIMER = "metal_ligand_interface_dimer"


@dataclass
class AnalysisResult:
    """Result of a single analysis module."""
    status: AnalysisStatus
    metrics: Optional[Dict[str, Any]] = None
    reason: Optional[str] = None

    @classmethod
    def success(cls, metrics: Dict[str, Any]) -> "AnalysisResult":
        """Create a successful analysis result."""
        return cls(status=AnalysisStatus.SUCCESS, metrics=metrics)

    @classmethod
    def not_applicable(cls, reason: str) -> "AnalysisResult":
        """Create a not-applicable result (analysis doesn't apply to this design)."""
        return cls(status=AnalysisStatus.NOT_APPLICABLE, reason=reason)

    @classmethod
    def skipped(cls, reason: str) -> "AnalysisResult":
        """Create a skipped result (analysis could apply but couldn't run)."""
        return cls(status=AnalysisStatus.SKIPPED, reason=reason)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        result = {"status": self.status.value}
        if self.metrics is not None:
            result["metrics"] = self.metrics
        if self.reason is not None:
            result["reason"] = self.reason
        return result


def detect_design_type(
    has_ligand: bool,
    has_metal: bool,
    chain_count: int,
) -> DesignType:
    """
    Detect design type from structural properties.

    Args:
        has_ligand: Whether a small molecule ligand is present
        has_metal: Whether a metal ion is present
        chain_count: Number of protein chains

    Returns:
        DesignType enum value
    """
    is_dimer = chain_count >= 2

    if has_ligand and has_metal and is_dimer:
        return DesignType.METAL_LIGAND_INTERFACE_DIMER
    elif has_ligand and is_dimer:
        return DesignType.LIGAND_INTERFACE_DIMER
    elif has_metal and is_dimer:
        return DesignType.METAL_INTERFACE_DIMER
    elif is_dimer:
        return DesignType.PROTEIN_DIMER
    else:
        return DesignType.MONOMER
