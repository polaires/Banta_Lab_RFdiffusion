"""
Geometry Validator

Validates metal-ligand-protein geometry:
- Metal-ligand coordination distances
- Carbon-metal distances (must be >3Å)
- Ligand internal bond lengths
- Coordination geometry
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Any

from .report import (
    ValidationIssue,
    IssueSeverity,
    IssueType,
    ActionType,
    GeometryReport,
)


# Expected bond lengths (Å)
BOND_LENGTH_RANGES = {
    # C-C bonds
    ("C", "C"): (1.45, 1.60),  # Single bond
    # C-O bonds
    ("C", "O"): (1.20, 1.50),  # Ranges from C=O to C-O
    # C-N bonds
    ("C", "N"): (1.25, 1.55),
    # Metal-O coordination
    ("TB", "O"): (2.25, 2.70),
    ("EU", "O"): (2.30, 2.75),
    ("GD", "O"): (2.25, 2.70),
    ("CA", "O"): (2.30, 2.60),
    ("MG", "O"): (2.00, 2.30),
    ("ZN", "O"): (1.95, 2.25),
    ("ZN", "N"): (2.00, 2.15),
    ("ZN", "S"): (2.25, 2.40),
    ("FE", "O"): (1.90, 2.20),
    ("FE", "N"): (1.95, 2.15),
}

# Minimum distance from metal to carbon (should be >3Å)
MIN_CARBON_METAL_DISTANCE = 3.0


@dataclass
class Atom:
    """Parsed atom from PDB."""
    serial: int
    name: str
    res_name: str
    chain: str
    res_num: int
    coords: np.ndarray
    element: str

    @property
    def full_id(self) -> str:
        return f"{self.res_name}:{self.chain}:{self.res_num}:{self.name}"


class GeometryValidator:
    """Validate metal-ligand-protein geometry."""

    def __init__(
        self,
        min_carbon_metal_dist: float = MIN_CARBON_METAL_DISTANCE,
        custom_bond_ranges: Optional[Dict] = None,
    ):
        self.min_carbon_metal_dist = min_carbon_metal_dist
        self.bond_ranges = {**BOND_LENGTH_RANGES}
        if custom_bond_ranges:
            self.bond_ranges.update(custom_bond_ranges)

        self.atoms: List[Atom] = []
        self.metal_atom: Optional[Atom] = None
        self.ligand_atoms: List[Atom] = []
        self.protein_atoms: List[Atom] = []

    def validate(
        self,
        pdb_content: str,
        metal: str,
        ligand: str,
    ) -> Tuple[GeometryReport, List[ValidationIssue]]:
        """
        Run geometry validation.

        Args:
            pdb_content: PDB file content
            metal: Metal element code (e.g., "TB")
            ligand: Ligand residue name (e.g., "CIT")

        Returns:
            Tuple of (GeometryReport, list of issues)
        """
        self._parse_pdb(pdb_content, metal, ligand)

        issues = []
        report = GeometryReport()

        # Check metal-ligand distances
        ml_issues = self._check_metal_ligand_distances(metal)
        issues.extend(ml_issues)
        report.metal_ligand_distances = self._get_metal_ligand_distances()

        # Check carbon-metal distances
        cm_issues, cm_dists = self._check_carbon_metal_distances()
        issues.extend(cm_issues)
        report.metal_carbon_distances = cm_dists
        report.all_carbons_far_from_metal = len(cm_issues) == 0
        if cm_dists:
            report.min_carbon_metal_distance = min(cm_dists.values())

        # Check ligand bond lengths
        bl_issues = self._check_ligand_bond_lengths()
        issues.extend(bl_issues)
        report.all_bonds_reasonable = len(bl_issues) == 0

        return report, issues

    def _parse_pdb(self, pdb_content: str, metal: str, ligand: str):
        """Parse PDB content into atoms."""
        self.atoms = []
        self.metal_atom = None
        self.ligand_atoms = []
        self.protein_atoms = []

        metal_upper = metal.upper()
        ligand_upper = ligand.upper()

        for line in pdb_content.split("\n"):
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            try:
                serial = int(line[6:11].strip())
                name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21:22].strip() or "A"
                res_num = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                # Determine element
                if len(line) >= 78:
                    element = line[76:78].strip()
                else:
                    # Guess from atom name
                    element = name[0] if name else "X"

                atom = Atom(
                    serial=serial,
                    name=name,
                    res_name=res_name,
                    chain=chain,
                    res_num=res_num,
                    coords=np.array([x, y, z]),
                    element=element,
                )

                self.atoms.append(atom)

                # Categorize
                if res_name == metal_upper or element == metal_upper:
                    self.metal_atom = atom
                elif res_name == ligand_upper:
                    self.ligand_atoms.append(atom)
                elif line.startswith("ATOM"):
                    self.protein_atoms.append(atom)

            except (ValueError, IndexError):
                continue

    def _check_metal_ligand_distances(self, metal: str) -> List[ValidationIssue]:
        """Check metal to coordinating ligand atoms."""
        issues = []

        if not self.metal_atom or not self.ligand_atoms:
            return issues

        metal_upper = metal.upper()
        metal_coord = self.metal_atom.coords

        for atom in self.ligand_atoms:
            if atom.element not in ["O", "N", "S"]:
                continue  # Only check potential coordinating atoms

            dist = np.linalg.norm(atom.coords - metal_coord)

            # Get expected range
            key = (metal_upper, atom.element)
            expected = self.bond_ranges.get(key, (2.0, 3.0))

            # Check if coordinating (within ~3Å)
            if dist < 3.5:  # Potentially coordinating
                if dist < expected[0]:
                    issues.append(ValidationIssue(
                        type=IssueType.DISTANCE,
                        severity=IssueSeverity.WARNING,
                        description=f"Metal-ligand distance too short: {metal}-{atom.name} = {dist:.2f}Å (expected {expected[0]:.2f}-{expected[1]:.2f}Å)",
                        atoms=[self.metal_atom.full_id, atom.full_id],
                        measured_value=dist,
                        expected_range=expected,
                        suggested_action=ActionType.RELAX,
                    ))
                elif dist > expected[1]:
                    issues.append(ValidationIssue(
                        type=IssueType.DISTANCE,
                        severity=IssueSeverity.WARNING,
                        description=f"Metal-ligand distance too long: {metal}-{atom.name} = {dist:.2f}Å (expected {expected[0]:.2f}-{expected[1]:.2f}Å)",
                        atoms=[self.metal_atom.full_id, atom.full_id],
                        measured_value=dist,
                        expected_range=expected,
                        suggested_action=ActionType.TRANSLATE,
                    ))

        return issues

    def _get_metal_ligand_distances(self) -> Dict[str, float]:
        """Get distances from metal to all ligand atoms."""
        distances = {}

        if not self.metal_atom:
            return distances

        metal_coord = self.metal_atom.coords

        for atom in self.ligand_atoms:
            dist = np.linalg.norm(atom.coords - metal_coord)
            distances[atom.name] = round(dist, 2)

        return distances

    def _check_carbon_metal_distances(self) -> Tuple[List[ValidationIssue], Dict[str, float]]:
        """Check that all carbons are >3Å from metal."""
        issues = []
        distances = {}

        if not self.metal_atom:
            return issues, distances

        metal_coord = self.metal_atom.coords

        for atom in self.ligand_atoms:
            if atom.element != "C" and "C" not in atom.name:
                continue

            dist = np.linalg.norm(atom.coords - metal_coord)
            distances[atom.name] = round(dist, 2)

            if dist < self.min_carbon_metal_dist:
                issues.append(ValidationIssue(
                    type=IssueType.GEOMETRY,
                    severity=IssueSeverity.CRITICAL,
                    description=f"Carbon too close to metal: {atom.name} at {dist:.2f}Å (must be >{self.min_carbon_metal_dist}Å)",
                    atoms=[self.metal_atom.full_id, atom.full_id],
                    measured_value=dist,
                    expected_range=(self.min_carbon_metal_dist, None),
                    suggested_action=ActionType.REGENERATE,
                    details={"carbon_name": atom.name, "metal": self.metal_atom.element},
                ))

        return issues, distances

    def _check_ligand_bond_lengths(self) -> List[ValidationIssue]:
        """Check internal ligand bond lengths."""
        issues = []

        if len(self.ligand_atoms) < 2:
            return issues

        # Find likely bonds (atoms within 1.7Å)
        for i, atom1 in enumerate(self.ligand_atoms):
            for atom2 in self.ligand_atoms[i + 1:]:
                dist = np.linalg.norm(atom1.coords - atom2.coords)

                # Check if this looks like a bond (0.9-1.8Å typical)
                if 0.9 < dist < 1.8:
                    # Get expected range
                    key1 = (atom1.element, atom2.element)
                    key2 = (atom2.element, atom1.element)
                    expected = self.bond_ranges.get(key1) or self.bond_ranges.get(key2)

                    if expected:
                        if dist < expected[0] - 0.1:  # Allow small tolerance
                            issues.append(ValidationIssue(
                                type=IssueType.BOND_LENGTH,
                                severity=IssueSeverity.WARNING,
                                description=f"Bond too short: {atom1.name}-{atom2.name} = {dist:.2f}Å (expected {expected[0]:.2f}-{expected[1]:.2f}Å)",
                                atoms=[atom1.full_id, atom2.full_id],
                                measured_value=dist,
                                expected_range=expected,
                            ))
                        elif dist > expected[1] + 0.1:
                            issues.append(ValidationIssue(
                                type=IssueType.BOND_LENGTH,
                                severity=IssueSeverity.WARNING,
                                description=f"Bond too long: {atom1.name}-{atom2.name} = {dist:.2f}Å (expected {expected[0]:.2f}-{expected[1]:.2f}Å)",
                                atoms=[atom1.full_id, atom2.full_id],
                                measured_value=dist,
                                expected_range=expected,
                            ))

        return issues

    def get_coordinating_atoms(self, max_distance: float = 3.0) -> List[Dict[str, Any]]:
        """Get list of atoms coordinating the metal."""
        coordinating = []

        if not self.metal_atom:
            return coordinating

        metal_coord = self.metal_atom.coords

        # Check ligand atoms
        for atom in self.ligand_atoms:
            if atom.element in ["O", "N", "S"]:
                dist = np.linalg.norm(atom.coords - metal_coord)
                if dist < max_distance:
                    coordinating.append({
                        "source": "ligand",
                        "atom": atom.name,
                        "element": atom.element,
                        "distance": round(dist, 2),
                        "res_name": atom.res_name,
                    })

        # Check protein atoms
        for atom in self.protein_atoms:
            if atom.element in ["O", "N", "S"]:
                dist = np.linalg.norm(atom.coords - metal_coord)
                if dist < max_distance:
                    coordinating.append({
                        "source": "protein",
                        "atom": atom.name,
                        "element": atom.element,
                        "distance": round(dist, 2),
                        "res_name": atom.res_name,
                        "chain": atom.chain,
                        "res_num": atom.res_num,
                    })

        return sorted(coordinating, key=lambda x: x["distance"])
