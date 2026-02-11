"""
Analysis Type Definitions

Provides structured types for design analysis results with clear status handling.
Includes auto-detection utilities for metals and ligands in PDB files.
"""
from enum import Enum
from typing import Dict, Any, Optional, List, Tuple
from dataclasses import dataclass


# Common metal codes found in PDB files (including lanthanides)
METAL_CODES = {
    # Transition metals
    "FE", "ZN", "CA", "MG", "MN", "CO", "CU", "NI", "MO", "W",
    # Lanthanides (common in protein engineering)
    "TB", "EU", "GD", "DY", "SM", "ND", "LA", "CE", "PR", "PM",
    "HO", "ER", "TM", "YB", "LU",
    # Other metals
    "NA", "K", "CD", "HG", "PB", "AG", "AU", "PT", "RU", "RH",
    "PD", "IR", "BA", "SR", "AL", "GA", "IN", "TL", "CR", "V",
}

# Solvent/water residue names to exclude from ligand detection
SOLVENT_RESIDUES = {"HOH", "WAT", "DOD", "H2O", "SOL", "TIP", "TIP3", "SPC"}


@dataclass
class DetectedMetal:
    """Auto-detected metal from PDB content."""
    metal_type: str
    chain: str
    resnum: int
    x: float
    y: float
    z: float


@dataclass
class DetectedLigand:
    """Auto-detected ligand from PDB content."""
    ligand_name: str
    chain: str
    resnum: int
    atom_count: int


def detect_metal_from_pdb(pdb_content: str) -> Optional[DetectedMetal]:
    """
    Scan PDB HETATM records for metal codes.

    Checks residue name (cols 17-20), element symbol (cols 76-78), and
    atom name (cols 12-16) to handle RF3 outputs where metal may be in
    a ligand chain with non-standard residue names like "L:0".

    Args:
        pdb_content: PDB file content as string

    Returns:
        DetectedMetal with chain, resnum, coordinates if found, else None
    """
    for line in pdb_content.split("\n"):
        if not line.startswith("HETATM"):
            continue

        try:
            # PDB format columns:
            # 12-16: atom name
            # 17-20: residue name (right-justified in cols 17-20)
            # 21: chain ID
            # 22-26: residue sequence number
            # 30-38, 38-46, 46-54: x, y, z coordinates
            # 76-78: element symbol
            residue_name = line[17:20].strip().upper()

            # Check residue name first
            metal_type = None
            if residue_name in METAL_CODES:
                metal_type = residue_name
            else:
                # Check element symbol (cols 76-78) for RF3 outputs
                if len(line) > 76:
                    element = line[76:78].strip().upper()
                    if element in METAL_CODES:
                        metal_type = element
                # Check atom name (cols 12-16) as fallback
                if metal_type is None:
                    atom_name = line[12:16].strip().upper()
                    # Remove trailing numbers (e.g., TB0 -> TB)
                    atom_base = atom_name.rstrip("0123456789")
                    if atom_base in METAL_CODES:
                        metal_type = atom_base

            if metal_type:
                chain = line[21].strip() or "A"
                resnum = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())

                return DetectedMetal(
                    metal_type=metal_type,
                    chain=chain,
                    resnum=resnum,
                    x=x,
                    y=y,
                    z=z,
                )
        except (ValueError, IndexError):
            continue

    return None


def detect_ligand_from_pdb(pdb_content: str) -> Optional[DetectedLigand]:
    """
    Scan PDB HETATM records for non-metal, non-water residues.

    Args:
        pdb_content: PDB file content as string

    Returns:
        DetectedLigand with name, chain, resnum, atom count if found, else None
    """
    # Group HETATM records by (residue_name, chain, resnum)
    ligand_atoms: Dict[Tuple[str, str, int], int] = {}

    for line in pdb_content.split("\n"):
        if not line.startswith("HETATM"):
            continue

        try:
            residue_name = line[17:20].strip().upper()

            # Skip metals and solvents
            if residue_name in METAL_CODES or residue_name in SOLVENT_RESIDUES:
                continue

            chain = line[21].strip() or "A"
            resnum = int(line[22:26].strip())

            key = (residue_name, chain, resnum)
            ligand_atoms[key] = ligand_atoms.get(key, 0) + 1
        except (ValueError, IndexError):
            continue

    if not ligand_atoms:
        return None

    # Return the ligand with the most atoms (likely the main ligand)
    best_key = max(ligand_atoms, key=lambda k: ligand_atoms[k])
    return DetectedLigand(
        ligand_name=best_key[0],
        chain=best_key[1],
        resnum=best_key[2],
        atom_count=ligand_atoms[best_key],
    )


class AnalysisStatus(Enum):
    """Status of an analysis run."""
    SUCCESS = "success"
    NOT_APPLICABLE = "not_applicable"
    SKIPPED = "skipped"


class StructureType(Enum):
    """Structural composition of a protein design (monomer/dimer, metal/ligand).

    NOTE: This is different from design_types.DesignType which classifies
    the design INTENT (e.g., METAL_BINDING, ENZYME_ACTIVE_SITE).
    StructureType classifies the structure COMPOSITION after generation.
    """
    # Monomer designs (single chain)
    MONOMER = "monomer"
    METAL_MONOMER = "metal_monomer"
    LIGAND_MONOMER = "ligand_monomer"
    METAL_LIGAND_MONOMER = "metal_ligand_monomer"
    # Dimer designs (multiple chains)
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


@dataclass
class ConservationResult:
    """
    Conservation analysis result for a single position.

    Follows ConSurf GradesPE format for compatibility.
    Grade 1 = most conserved, Grade 9 = most variable.
    """
    position: int                      # 1-indexed position
    residue: str                       # Single-letter amino acid
    grade: int                         # 1-9 ConSurf grade
    score: float                       # Raw Rate4Site score
    confidence: Tuple[float, float]    # Bayesian confidence interval
    data_quality: str                  # "sufficient" or "insufficient"

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "position": self.position,
            "residue": self.residue,
            "grade": self.grade,
            "score": self.score,
            "confidence": list(self.confidence),
            "data_quality": self.data_quality,
        }


@dataclass
class ConservationSummary:
    """
    Summary of conservation analysis for design integration.

    Provides actionable information for scaffolding and partial diffusion:
    - highly_conserved_positions: Should be fixed during design
    - variable_positions: Can be freely redesigned
    """
    highly_conserved_positions: List[int]   # Grade 1-3 (most important)
    conserved_positions: List[int]          # Grade 4-5
    variable_positions: List[int]           # Grade 7-9 (can redesign)
    msa_depth: int                          # Number of sequences in alignment
    average_conservation: float             # Mean grade (1-9)
    method: str                             # "bayesian", "ml", or "entropy"
    reliable: bool                          # True if MSA depth >= 50

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "highly_conserved_positions": self.highly_conserved_positions,
            "conserved_positions": self.conserved_positions,
            "variable_positions": self.variable_positions,
            "msa_depth": self.msa_depth,
            "average_conservation": self.average_conservation,
            "method": self.method,
            "reliable": self.reliable,
        }

    def get_suggested_fixed_residues(self, chain: str = "A") -> List[str]:
        """
        Get residue IDs suggested for fixing during RFD3 design.

        Returns list in format ["A1", "A5", "A10"] for use with
        select_fixed_atoms in RFD3 config.
        """
        return [f"{chain}{pos}" for pos in self.highly_conserved_positions]


def detect_design_type(
    has_ligand: bool,
    has_metal: bool,
    chain_count: int,
) -> StructureType:
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

    # Dimer types (2+ protein chains)
    if has_ligand and has_metal and is_dimer:
        return StructureType.METAL_LIGAND_INTERFACE_DIMER
    elif has_ligand and is_dimer:
        return StructureType.LIGAND_INTERFACE_DIMER
    elif has_metal and is_dimer:
        return StructureType.METAL_INTERFACE_DIMER
    elif is_dimer:
        return StructureType.PROTEIN_DIMER
    # Monomer types (single protein chain)
    elif has_ligand and has_metal:
        return StructureType.METAL_LIGAND_MONOMER
    elif has_ligand:
        return StructureType.LIGAND_MONOMER
    elif has_metal:
        return StructureType.METAL_MONOMER
    else:
        return StructureType.MONOMER
