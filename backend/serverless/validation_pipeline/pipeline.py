"""
Validation Pipeline

Main orchestrator that runs all validators and produces unified reports.
"""

from typing import Optional, Dict, Any

from .report import (
    ValidationReport,
    ValidationIssue,
    IssueSeverity,
    IssueType,
    CoordinationReport,
)
from .geometry_validator import GeometryValidator
from .clash_detector import ClashDetector


class ValidationPipeline:
    """
    Main validation pipeline for metal-ligand-protein structures.

    Usage:
        pipeline = ValidationPipeline()
        report = pipeline.validate(pdb_content, metal="TB", ligand="CIT")

        if not report.passed:
            print(report.to_summary())
            for issue in report.critical_issues:
                print(f"FIX: {issue.description}")
    """

    def __init__(
        self,
        min_carbon_metal_dist: float = 3.0,
        clash_tolerance: float = 0.4,
        ignore_hydrogens: bool = True,
    ):
        """
        Initialize validation pipeline.

        Args:
            min_carbon_metal_dist: Minimum allowed C-M distance (Å)
            clash_tolerance: Tolerance for clash detection (Å)
            ignore_hydrogens: Skip hydrogen atoms in clash detection
        """
        self.geometry_validator = GeometryValidator(
            min_carbon_metal_dist=min_carbon_metal_dist,
        )
        self.clash_detector = ClashDetector(
            clash_tolerance=clash_tolerance,
            ignore_hydrogens=ignore_hydrogens,
        )

    def validate(
        self,
        pdb_content: str,
        metal: str,
        ligand: str,
        target_coordination: int = 9,
        expected_ligand_donors: int = 3,
    ) -> ValidationReport:
        """
        Run full validation pipeline.

        Args:
            pdb_content: PDB file content as string
            metal: Metal element code (e.g., "TB", "CA", "ZN")
            ligand: Ligand residue name (e.g., "CIT", "PQQ")
            target_coordination: Expected total coordination number
            expected_ligand_donors: Expected donors from ligand

        Returns:
            ValidationReport with all issues and scores
        """
        report = ValidationReport()
        report.metadata = {
            "metal": metal,
            "ligand": ligand,
            "target_coordination": target_coordination,
            "expected_ligand_donors": expected_ligand_donors,
        }

        # Run geometry validation
        geometry_report, geometry_issues = self.geometry_validator.validate(
            pdb_content, metal, ligand
        )
        report.geometry = geometry_report
        for issue in geometry_issues:
            report.add_issue(issue)

        # Run clash detection
        clash_report, clash_issues = self.clash_detector.detect(
            pdb_content, metal, ligand
        )
        report.clashes = clash_report
        for issue in clash_issues:
            report.add_issue(issue)

        # Analyze coordination
        coord_report = self._analyze_coordination(
            pdb_content, metal, ligand,
            target_coordination, expected_ligand_donors
        )
        report.coordination = coord_report

        # Check coordination issues
        if coord_report.total_coordination < target_coordination:
            deficit = target_coordination - coord_report.total_coordination
            severity = IssueSeverity.CRITICAL if deficit > 3 else IssueSeverity.WARNING
            report.add_issue(ValidationIssue(
                type=IssueType.COORDINATION,
                severity=severity,
                description=f"Incomplete coordination: {coord_report.total_coordination}/{target_coordination} (missing {deficit})",
                details={"deficit": deficit},
            ))

        if not coord_report.hsab_compliant:
            report.add_issue(ValidationIssue(
                type=IssueType.HSAB,
                severity=IssueSeverity.WARNING,
                description="HSAB theory violation: incompatible donor types",
            ))

        # Calculate final score
        report.calculate_score()

        return report

    def _analyze_coordination(
        self,
        pdb_content: str,
        metal: str,
        ligand: str,
        target_coordination: int,
        expected_ligand_donors: int,
    ) -> CoordinationReport:
        """Analyze metal coordination sphere."""
        report = CoordinationReport()
        report.target_coordination = target_coordination

        # Get coordinating atoms from geometry validator
        coordinating = self.geometry_validator.get_coordinating_atoms(max_distance=3.5)

        ligand_donors = [c for c in coordinating if c["source"] == "ligand"]
        protein_donors = [c for c in coordinating if c["source"] == "protein"]

        report.ligand_donors = len(ligand_donors)
        report.protein_donors = len(protein_donors)
        report.total_coordination = report.ligand_donors + report.protein_donors
        report.coordination_complete = report.total_coordination >= target_coordination

        report.ligand_donor_atoms = ligand_donors
        report.donor_residues = protein_donors

        # Check HSAB compliance
        report.hsab_compliant = self._check_hsab_compliance(
            metal, [d["element"] for d in coordinating]
        )

        return report

    def _check_hsab_compliance(self, metal: str, donor_elements: list) -> bool:
        """Check if donor elements are HSAB-compliant with metal."""
        hard_acids = {"TB", "EU", "GD", "LA", "CE", "SM", "YB", "CA", "MG"}
        soft_acids = {"CU", "AG", "AU", "HG"}

        metal_upper = metal.upper()

        if metal_upper in hard_acids:
            # Hard acids should not have S donors
            if "S" in donor_elements:
                return False
        elif metal_upper in soft_acids:
            # Soft acids prefer S, N over O
            pass  # Usually ok

        return True

    def validate_quick(
        self,
        pdb_content: str,
        metal: str,
        ligand: str,
    ) -> Dict[str, Any]:
        """
        Quick validation returning simple dict (for API responses).

        Returns:
            Dict with pass/fail status and key metrics
        """
        report = self.validate(pdb_content, metal, ligand)

        return {
            "passed": report.passed,
            "score": report.quality_score,
            "critical_issues": len(report.critical_issues),
            "warnings": len(report.warnings),
            "clashes": report.clashes.total_clash_count if report.clashes else 0,
            "coordination": report.coordination.total_coordination if report.coordination else 0,
            "summary": report.to_summary(),
        }

    def get_action_items(self, report: ValidationReport) -> list:
        """
        Get list of suggested actions to fix issues.

        Returns:
            List of action items with descriptions
        """
        actions = []

        for issue in report.issues:
            if issue.suggested_action:
                actions.append({
                    "action": issue.suggested_action.value,
                    "reason": issue.description,
                    "severity": issue.severity.value,
                    "atoms": issue.atoms,
                })

        # Sort by severity (critical first)
        severity_order = {"critical": 0, "warning": 1, "info": 2}
        actions.sort(key=lambda x: severity_order.get(x["severity"], 3))

        return actions
