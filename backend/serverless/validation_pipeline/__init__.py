"""
Metal-Ligand Complex Validation Pipeline

Systematic validation and correction infrastructure for protein-metal-ligand designs.

Components:
- GeometryValidator: Check distances, bond lengths, coordination geometry
- ClashDetector: Find steric overlaps between atoms
- CoordinationChecker: Verify metal coordination sphere
- ValidationReport: Structured output for issues and actions

Usage:
    from validation_pipeline import ValidationPipeline

    pipeline = ValidationPipeline()
    report = pipeline.validate(pdb_content, metal="TB", ligand="CIT")

    if not report.passed:
        for issue in report.issues:
            print(f"{issue.severity}: {issue.description}")
"""

from .report import ValidationReport, ValidationIssue, IssueSeverity
from .geometry_validator import GeometryValidator
from .clash_detector import ClashDetector
from .pipeline import ValidationPipeline

__all__ = [
    "ValidationPipeline",
    "ValidationReport",
    "ValidationIssue",
    "IssueSeverity",
    "GeometryValidator",
    "ClashDetector",
]
