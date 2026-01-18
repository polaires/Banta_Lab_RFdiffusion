"""
Clash Detector

Detects steric clashes between atoms based on van der Waals radii.
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Any

from .report import (
    ValidationIssue,
    IssueSeverity,
    IssueType,
    ActionType,
    ClashReport,
)


# Van der Waals radii in Angstroms
VDW_RADII = {
    # Common elements
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "P": 1.80,
    # Metals
    "CA": 1.74,
    "MG": 1.45,
    "ZN": 1.39,
    "FE": 1.40,
    "CU": 1.40,
    "MN": 1.40,
    "CO": 1.40,
    "NI": 1.40,
    # Lanthanides
    "TB": 1.75,
    "EU": 1.85,
    "GD": 1.80,
    "LA": 1.95,
    "CE": 1.85,
    "SM": 1.80,
    "YB": 1.70,
}

# Default radius for unknown elements
DEFAULT_VDW_RADIUS = 1.70

# Clash tolerance - atoms within (sum_radii - tolerance) are clashing
CLASH_TOLERANCE = 0.4


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

    def get_vdw_radius(self) -> float:
        """Get van der Waals radius for this atom."""
        return VDW_RADII.get(self.element, DEFAULT_VDW_RADIUS)


class ClashDetector:
    """Detect steric clashes between atoms."""

    def __init__(
        self,
        clash_tolerance: float = CLASH_TOLERANCE,
        ignore_hydrogens: bool = True,
        custom_radii: Optional[Dict[str, float]] = None,
    ):
        self.clash_tolerance = clash_tolerance
        self.ignore_hydrogens = ignore_hydrogens
        self.radii = {**VDW_RADII}
        if custom_radii:
            self.radii.update(custom_radii)

        self.atoms: List[Atom] = []
        self.metal_atom: Optional[Atom] = None
        self.ligand_atoms: List[Atom] = []
        self.protein_atoms: List[Atom] = []

    def detect(
        self,
        pdb_content: str,
        metal: str,
        ligand: str,
    ) -> Tuple[ClashReport, List[ValidationIssue]]:
        """
        Detect all clashes in the structure.

        Args:
            pdb_content: PDB file content
            metal: Metal element code (e.g., "TB")
            ligand: Ligand residue name (e.g., "CIT")

        Returns:
            Tuple of (ClashReport, list of issues)
        """
        self._parse_pdb(pdb_content, metal, ligand)

        issues = []
        report = ClashReport()

        # Ligand-protein clashes
        lp_clashes, lp_issues = self._check_ligand_protein_clashes()
        report.ligand_protein_clashes = lp_clashes
        issues.extend(lp_issues)

        # Metal-protein clashes
        mp_clashes, mp_issues = self._check_metal_protein_clashes()
        report.metal_protein_clashes = mp_clashes
        issues.extend(mp_issues)

        # Calculate totals
        report.total_clash_count = len(lp_clashes) + len(mp_clashes)

        all_clashes = lp_clashes + mp_clashes
        if all_clashes:
            overlaps = [c.get("overlap", 0) for c in all_clashes]
            report.worst_overlap = max(overlaps) if overlaps else 0.0

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

                # Skip hydrogens if configured
                if self.ignore_hydrogens and element == "H":
                    continue

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

    def _check_clash(self, atom1: Atom, atom2: Atom) -> Optional[Dict[str, Any]]:
        """Check if two atoms clash."""
        dist = np.linalg.norm(atom1.coords - atom2.coords)

        r1 = self.radii.get(atom1.element, DEFAULT_VDW_RADIUS)
        r2 = self.radii.get(atom2.element, DEFAULT_VDW_RADIUS)
        min_dist = r1 + r2 - self.clash_tolerance

        if dist < min_dist:
            overlap = min_dist - dist
            return {
                "atom1": atom1.full_id,
                "atom2": atom2.full_id,
                "distance": round(dist, 2),
                "min_allowed": round(min_dist, 2),
                "overlap": round(overlap, 2),
                "element1": atom1.element,
                "element2": atom2.element,
            }

        return None

    def _check_ligand_protein_clashes(self) -> Tuple[List[Dict], List[ValidationIssue]]:
        """Check for clashes between ligand and protein atoms."""
        clashes = []
        issues = []

        for lig_atom in self.ligand_atoms:
            for prot_atom in self.protein_atoms:
                clash = self._check_clash(lig_atom, prot_atom)
                if clash:
                    clashes.append(clash)

                    # Determine severity based on overlap
                    overlap = clash["overlap"]
                    if overlap > 1.0:
                        severity = IssueSeverity.CRITICAL
                    elif overlap > 0.5:
                        severity = IssueSeverity.WARNING
                    else:
                        severity = IssueSeverity.INFO

                    issues.append(ValidationIssue(
                        type=IssueType.CLASH,
                        severity=severity,
                        description=f"Ligand-protein clash: {lig_atom.name} ↔ {prot_atom.res_name}:{prot_atom.chain}:{prot_atom.res_num}:{prot_atom.name} ({clash['distance']:.2f}Å, overlap {overlap:.2f}Å)",
                        atoms=[lig_atom.full_id, prot_atom.full_id],
                        measured_value=clash["distance"],
                        expected_range=(clash["min_allowed"], None),
                        suggested_action=ActionType.TRANSLATE,
                        details=clash,
                    ))

        return clashes, issues

    def _check_metal_protein_clashes(self) -> Tuple[List[Dict], List[ValidationIssue]]:
        """Check for clashes between metal and protein backbone."""
        clashes = []
        issues = []

        if not self.metal_atom:
            return clashes, issues

        # Only check backbone atoms (N, CA, C, O)
        backbone_atoms = ["N", "CA", "C", "O"]

        for prot_atom in self.protein_atoms:
            # Skip if this is a coordinating sidechain
            if prot_atom.name not in backbone_atoms:
                continue

            clash = self._check_clash(self.metal_atom, prot_atom)
            if clash:
                clashes.append(clash)

                issues.append(ValidationIssue(
                    type=IssueType.CLASH,
                    severity=IssueSeverity.CRITICAL,
                    description=f"Metal-backbone clash: {self.metal_atom.element} ↔ {prot_atom.res_name}:{prot_atom.chain}:{prot_atom.res_num}:{prot_atom.name} ({clash['distance']:.2f}Å)",
                    atoms=[self.metal_atom.full_id, prot_atom.full_id],
                    measured_value=clash["distance"],
                    suggested_action=ActionType.TRANSLATE,
                    details=clash,
                ))

        return clashes, issues

    def get_closest_contacts(
        self,
        n: int = 10,
    ) -> List[Dict[str, Any]]:
        """Get N closest ligand-protein contacts (for debugging)."""
        contacts = []

        for lig_atom in self.ligand_atoms:
            for prot_atom in self.protein_atoms:
                dist = np.linalg.norm(lig_atom.coords - prot_atom.coords)
                contacts.append({
                    "ligand_atom": lig_atom.name,
                    "protein_atom": f"{prot_atom.res_name}:{prot_atom.chain}:{prot_atom.res_num}:{prot_atom.name}",
                    "distance": round(dist, 2),
                })

        return sorted(contacts, key=lambda x: x["distance"])[:n]
