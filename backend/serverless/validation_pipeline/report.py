"""
Validation Report Schema

Structured output for validation results, issues, and suggested actions.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Any
import json


class IssueSeverity(Enum):
    """Severity levels for validation issues."""
    CRITICAL = "critical"  # Must fix before use
    WARNING = "warning"    # Should fix but usable
    INFO = "info"          # Informational only


class IssueType(Enum):
    """Types of validation issues."""
    CLASH = "clash"
    DISTANCE = "distance"
    COORDINATION = "coordination"
    GEOMETRY = "geometry"
    BOND_LENGTH = "bond_length"
    HSAB = "hsab_violation"


class ActionType(Enum):
    """Types of correction actions."""
    TRANSLATE = "translate"
    ROTATE = "rotate"
    REGENERATE = "regenerate"
    RELAX = "relax"
    REDESIGN = "redesign"


@dataclass
class ValidationIssue:
    """A single validation issue."""
    type: IssueType
    severity: IssueSeverity
    description: str
    atoms: List[str] = field(default_factory=list)  # e.g., ["CIT:C3", "VAL:A:45:CB"]
    measured_value: Optional[float] = None
    expected_range: Optional[tuple] = None
    suggested_action: Optional[ActionType] = None
    details: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "type": self.type.value,
            "severity": self.severity.value,
            "description": self.description,
            "atoms": self.atoms,
            "measured_value": self.measured_value,
            "expected_range": self.expected_range,
            "suggested_action": self.suggested_action.value if self.suggested_action else None,
            "details": self.details,
        }


@dataclass
class GeometryReport:
    """Geometry validation results."""
    metal_ligand_distances: Dict[str, float] = field(default_factory=dict)
    metal_carbon_distances: Dict[str, float] = field(default_factory=dict)
    ligand_bond_lengths: Dict[str, float] = field(default_factory=dict)
    all_carbons_far_from_metal: bool = True
    all_bonds_reasonable: bool = True
    min_carbon_metal_distance: float = 999.0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "metal_ligand_distances": self.metal_ligand_distances,
            "metal_carbon_distances": self.metal_carbon_distances,
            "ligand_bond_lengths": self.ligand_bond_lengths,
            "all_carbons_far_from_metal": self.all_carbons_far_from_metal,
            "all_bonds_reasonable": self.all_bonds_reasonable,
            "min_carbon_metal_distance": self.min_carbon_metal_distance,
        }


@dataclass
class ClashReport:
    """Clash detection results."""
    ligand_protein_clashes: List[Dict[str, Any]] = field(default_factory=list)
    metal_protein_clashes: List[Dict[str, Any]] = field(default_factory=list)
    internal_clashes: List[Dict[str, Any]] = field(default_factory=list)
    total_clash_count: int = 0
    worst_overlap: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "ligand_protein_clashes": self.ligand_protein_clashes,
            "metal_protein_clashes": self.metal_protein_clashes,
            "internal_clashes": self.internal_clashes,
            "total_clash_count": self.total_clash_count,
            "worst_overlap": self.worst_overlap,
        }


@dataclass
class CoordinationReport:
    """Coordination analysis results."""
    total_coordination: int = 0
    ligand_donors: int = 0
    protein_donors: int = 0
    target_coordination: int = 9
    coordination_complete: bool = False
    donor_residues: List[Dict[str, Any]] = field(default_factory=list)
    ligand_donor_atoms: List[Dict[str, Any]] = field(default_factory=list)
    geometry: str = "unknown"
    hsab_compliant: bool = True

    def to_dict(self) -> Dict[str, Any]:
        return {
            "total_coordination": self.total_coordination,
            "ligand_donors": self.ligand_donors,
            "protein_donors": self.protein_donors,
            "target_coordination": self.target_coordination,
            "coordination_complete": self.coordination_complete,
            "donor_residues": self.donor_residues,
            "ligand_donor_atoms": self.ligand_donor_atoms,
            "geometry": self.geometry,
            "hsab_compliant": self.hsab_compliant,
        }


@dataclass
class ValidationReport:
    """Complete validation report."""
    passed: bool = False
    quality_score: int = 0  # 0-100
    issues: List[ValidationIssue] = field(default_factory=list)
    geometry: Optional[GeometryReport] = None
    clashes: Optional[ClashReport] = None
    coordination: Optional[CoordinationReport] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def critical_issues(self) -> List[ValidationIssue]:
        """Get only critical issues."""
        return [i for i in self.issues if i.severity == IssueSeverity.CRITICAL]

    @property
    def warnings(self) -> List[ValidationIssue]:
        """Get warning-level issues."""
        return [i for i in self.issues if i.severity == IssueSeverity.WARNING]

    @property
    def has_clashes(self) -> bool:
        """Check if any clashes exist."""
        return self.clashes and self.clashes.total_clash_count > 0

    @property
    def has_geometry_issues(self) -> bool:
        """Check if geometry has issues."""
        if not self.geometry:
            return False
        return not (self.geometry.all_carbons_far_from_metal and
                    self.geometry.all_bonds_reasonable)

    def add_issue(self, issue: ValidationIssue):
        """Add an issue to the report."""
        self.issues.append(issue)
        # Recalculate pass status
        if issue.severity == IssueSeverity.CRITICAL:
            self.passed = False

    def calculate_score(self) -> int:
        """Calculate quality score based on issues."""
        score = 100

        for issue in self.issues:
            if issue.severity == IssueSeverity.CRITICAL:
                score -= 30
            elif issue.severity == IssueSeverity.WARNING:
                score -= 10
            elif issue.severity == IssueSeverity.INFO:
                score -= 2

        # Bonus for complete coordination
        if self.coordination and self.coordination.coordination_complete:
            score += 10

        # Penalty for clashes
        if self.clashes:
            score -= self.clashes.total_clash_count * 15

        self.quality_score = max(0, min(100, score))
        self.passed = self.quality_score >= 70 and len(self.critical_issues) == 0
        return self.quality_score

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "passed": self.passed,
            "quality_score": self.quality_score,
            "issues": [i.to_dict() for i in self.issues],
            "geometry": self.geometry.to_dict() if self.geometry else None,
            "clashes": self.clashes.to_dict() if self.clashes else None,
            "coordination": self.coordination.to_dict() if self.coordination else None,
            "metadata": self.metadata,
            "summary": {
                "critical_count": len(self.critical_issues),
                "warning_count": len(self.warnings),
                "has_clashes": self.has_clashes,
                "has_geometry_issues": self.has_geometry_issues,
            }
        }

    def to_json(self, indent: int = 2) -> str:
        """Convert to JSON string."""
        return json.dumps(self.to_dict(), indent=indent)

    def to_summary(self) -> str:
        """Generate human-readable summary for AI feedback."""
        lines = [
            f"Validation {'PASSED' if self.passed else 'FAILED'} (Score: {self.quality_score}/100)",
            "",
        ]

        if self.critical_issues:
            lines.append("CRITICAL ISSUES:")
            for issue in self.critical_issues:
                lines.append(f"  - {issue.description}")
            lines.append("")

        if self.warnings:
            lines.append("WARNINGS:")
            for issue in self.warnings:
                lines.append(f"  - {issue.description}")
            lines.append("")

        if self.geometry:
            lines.append("GEOMETRY:")
            lines.append(f"  Carbons >3Å from metal: {self.geometry.all_carbons_far_from_metal}")
            lines.append(f"  Min C-M distance: {self.geometry.min_carbon_metal_distance:.2f}Å")
            lines.append("")

        if self.coordination:
            lines.append("COORDINATION:")
            lines.append(f"  Total: {self.coordination.total_coordination}/{self.coordination.target_coordination}")
            lines.append(f"  Ligand donors: {self.coordination.ligand_donors}")
            lines.append(f"  Protein donors: {self.coordination.protein_donors}")
            lines.append("")

        if self.clashes and self.clashes.total_clash_count > 0:
            lines.append("CLASHES:")
            lines.append(f"  Total: {self.clashes.total_clash_count}")
            lines.append(f"  Worst overlap: {self.clashes.worst_overlap:.2f}Å")
            for clash in self.clashes.ligand_protein_clashes[:3]:  # Top 3
                lines.append(f"  - {clash.get('atom1', '?')} ↔ {clash.get('atom2', '?')}: {clash.get('distance', 0):.2f}Å")

        return "\n".join(lines)
