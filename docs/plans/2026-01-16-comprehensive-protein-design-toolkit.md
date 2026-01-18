# Comprehensive Protein Design Toolkit Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement a unified protein design toolkit that provides intelligent tool selection (RFD3 vs RFdiffusion, LigandMPNN vs ProteinMPNN, FastRelax decisions) based on design requirements, with optimized presets for metal-binding, small molecule, DNA/RNA, and protein-protein binder design.

**Architecture:** Create a `DesignOrchestrator` that implements the decision tree from our analysis, selecting backbone tools, sequence design tools, relaxation strategies, and validation pipelines based on design type. Integrate with existing handlers while adding new capabilities for automatic workflow selection.

**Tech Stack:** Python 3.12, RFD3/RFdiffusion (Foundry), LigandMPNN/ProteinMPNN, FastRelax (PyRosetta), ESMFold/AF2-Multimer, biotite

---

## Background: The Complete Decision Framework

From our analysis of the LigandMPNN paper (Nature Methods 2025) and RFD3 documentation:

### Tool Selection Decision Tree

```
Does design involve non-protein molecules?
                    │
    ┌───────────────┴───────────────┐
    │                               │
   NO                              YES
    │                               │
    ▼                               ▼
ProteinMPNN                    What type?
    │                               │
    │               ┌───────────────┼───────────────┐
    │               │               │               │
    │          Small molecule     Metal         Nucleotide
    │               │               │               │
    │               └───────────────┴───────────────┘
    │                               │
    ▼                               ▼
Need sidechain                 LigandMPNN
optimization?                       │
    │                               ▼
┌───┴───┐                   Need backbone
│       │                    relaxation?
NO     YES                          │
│       │                   ┌───────┴───────┐
▼       ▼                   │               │
MPNN   MPNN                NO              YES
alone  + FastRelax          │               │
                            ▼               ▼
                       LigandMPNN     LigandMPNN
                         alone        + FastRelax
```

### Quantitative Performance (From Papers)

**Sequence Recovery Near Ligands:**
| Method | Small Molecules | Nucleotides | Metals |
|--------|----------------|-------------|--------|
| Rosetta | 50.4% | 35.2% | 36.0% |
| ProteinMPNN | 50.5% | 34.0% | 40.6% |
| **LigandMPNN** | **63.3%** | **50.5%** | **77.5%** |

**Key Insight:** For metal coordination, LigandMPNN achieves 77.5% vs 40.6% for ProteinMPNN - a **90% relative improvement**.

### The Alanine Problem (Solved)

Without proper biasing, LigandMPNN generates 60%+ alanine ("AAAA") sequences. Our testing validated:
- `bias_AA: "A:-2.0"` reduces alanine from 60% to <10%
- `bias_AA: "H:2.0,E:1.0,D:1.0"` favors metal-coordinating residues
- `fixed_positions` preserves coordination geometry
- `temperature: 0.1` for conservative sampling

---

## Task 1: Create Design Type Enum and Configuration

**Files:**
- Create: `G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless\design_types.py`
- Test: `G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless\test_design_types.py`

**Step 1: Write the failing test**

```python
# test_design_types.py
"""Tests for design type configuration."""

import pytest
from design_types import (
    DesignType,
    DesignConfig,
    get_recommended_config,
    DESIGN_PRESETS,
)


def test_design_type_enum_has_all_types():
    """Should have all documented design types."""
    assert DesignType.METAL_BINDING in DesignType
    assert DesignType.SMALL_MOLECULE_BINDER in DesignType
    assert DesignType.DNA_RNA_BINDER in DesignType
    assert DesignType.PROTEIN_PROTEIN_BINDER in DesignType
    assert DesignType.ENZYME_ACTIVE_SITE in DesignType
    assert DesignType.SYMMETRIC_OLIGOMER in DesignType
    assert DesignType.GENERAL_SCAFFOLD in DesignType


def test_metal_binding_config_uses_ligandmpnn():
    """Metal binding should use LigandMPNN, not ProteinMPNN."""
    config = get_recommended_config(DesignType.METAL_BINDING)

    assert config.sequence_tool == "ligand_mpnn"
    assert config.backbone_tool == "rfd3"
    assert "A:-" in config.bias_AA  # Alanine penalty
    assert config.use_atom_context == True


def test_protein_binder_config_uses_proteinmpnn():
    """Protein-protein binder should use ProteinMPNN."""
    config = get_recommended_config(DesignType.PROTEIN_PROTEIN_BINDER)

    assert config.sequence_tool == "protein_mpnn"
    assert config.relaxation == "fastrelax"


def test_small_molecule_config_uses_ligandmpnn():
    """Small molecule binder should use LigandMPNN."""
    config = get_recommended_config(DesignType.SMALL_MOLECULE_BINDER)

    assert config.sequence_tool == "ligand_mpnn"
    assert config.ligand_cutoff == 6.0  # Tighter for small molecules


def test_config_allows_overrides():
    """Should allow overriding preset values."""
    config = get_recommended_config(
        DesignType.METAL_BINDING,
        temperature=0.2,
        num_sequences=16,
    )

    assert config.temperature == 0.2
    assert config.num_sequences == 16
    # But still uses correct sequence tool
    assert config.sequence_tool == "ligand_mpnn"
```

**Step 2: Run test to verify it fails**

Run: `wsl -d Ubuntu -e bash -c "cd /mnt/g/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless && python -m pytest test_design_types.py::test_design_type_enum_has_all_types -v"`
Expected: FAIL with "ModuleNotFoundError: No module named 'design_types'"

**Step 3: Write minimal implementation**

```python
# design_types.py
"""
Protein Design Type Configuration

Implements the decision tree for selecting appropriate tools based on design type.
Based on:
- LigandMPNN paper (Nature Methods 2025)
- RFD3 documentation
- BindCraft source code analysis
"""

from enum import Enum, auto
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any


class DesignType(Enum):
    """Types of protein design tasks."""
    METAL_BINDING = auto()          # Metal ion coordination
    SMALL_MOLECULE_BINDER = auto()  # Drug-like molecule binding
    DNA_RNA_BINDER = auto()         # Nucleic acid binding
    PROTEIN_PROTEIN_BINDER = auto() # Protein-protein interface
    ENZYME_ACTIVE_SITE = auto()     # Catalytic site scaffolding
    SYMMETRIC_OLIGOMER = auto()     # Homo-oligomeric assemblies
    SYMMETRIC_LIGAND_BINDER = auto() # Symmetric binder with ligand
    METAL_MEDIATED_DIMER = auto()   # Dimer with metal at interface
    GENERAL_SCAFFOLD = auto()       # Novel folds, no specific function
    PEPTIDE_BINDER = auto()         # Peptide binding (BindCraft-style)


@dataclass
class DesignConfig:
    """Configuration for a protein design task."""

    design_type: DesignType

    # Backbone generation
    backbone_tool: str = "rfd3"  # "rfd3", "rfdiffusion", "bindcraft"
    use_symmetry: bool = False
    symmetry_type: Optional[str] = None  # "C2", "C3", "D2", etc.

    # Sequence design
    sequence_tool: str = "ligand_mpnn"  # "ligand_mpnn", "protein_mpnn"
    temperature: float = 0.1
    num_sequences: int = 8
    bias_AA: Optional[str] = None
    omit_AA: Optional[str] = None
    mpnn_weights: str = "default"  # "default", "soluble"

    # LigandMPNN-specific
    use_atom_context: bool = True
    ligand_cutoff: float = 8.0
    pack_side_chains: bool = True
    number_of_packs: int = 4
    model_noise_level: str = "010"

    # Relaxation
    relaxation: str = "optional"  # "none", "optional", "fastrelax"

    # Validation
    validation: str = "af2"  # "af2", "af2_multimer", "esmfold", "af3"

    # Fixed positions (auto-detected or manual)
    auto_detect_fixed: bool = True
    fixed_positions: Optional[List[str]] = None
    coordination_cutoff: float = 3.5  # For metal coordination detection

    # Metadata
    description: str = ""
    rationale: str = ""


# =============================================================================
# DESIGN PRESETS
# =============================================================================

DESIGN_PRESETS: Dict[DesignType, DesignConfig] = {

    DesignType.METAL_BINDING: DesignConfig(
        design_type=DesignType.METAL_BINDING,
        backbone_tool="rfd3",
        sequence_tool="ligand_mpnn",
        temperature=0.1,
        bias_AA="A:-2.0,H:2.0,C:1.5,E:1.0,D:1.0",
        omit_AA=None,  # Keep cysteines for metal coordination
        use_atom_context=True,
        ligand_cutoff=6.0,  # Tighter for small metal ion
        pack_side_chains=True,
        model_noise_level="010",
        relaxation="optional",  # FastRelax may disrupt metal geometry
        validation="af2",
        auto_detect_fixed=True,
        description="Metal-binding protein design",
        rationale="LigandMPNN 77.5% vs ProteinMPNN 40.6% for metals. Bias toward H/C/E/D coordinating residues.",
    ),

    DesignType.SMALL_MOLECULE_BINDER: DesignConfig(
        design_type=DesignType.SMALL_MOLECULE_BINDER,
        backbone_tool="rfd3",
        sequence_tool="ligand_mpnn",
        temperature=0.1,
        bias_AA="A:-1.0,W:1.5,Y:1.5,F:1.0,L:0.5",
        omit_AA="C",  # Avoid disulfide complications
        use_atom_context=True,
        ligand_cutoff=6.0,
        pack_side_chains=True,
        relaxation="optional",
        validation="af2",  # + docking validation
        description="Small molecule binder design",
        rationale="LigandMPNN 63.3% vs 50.5% for small molecules. Favor aromatics for drug-like compounds.",
    ),

    DesignType.DNA_RNA_BINDER: DesignConfig(
        design_type=DesignType.DNA_RNA_BINDER,
        backbone_tool="rfd3",
        sequence_tool="ligand_mpnn",
        temperature=0.1,
        bias_AA="A:-1.0,R:2.0,K:1.5,H:1.0",  # Positive charge for DNA backbone
        omit_AA=None,
        use_atom_context=True,
        ligand_cutoff=10.0,  # Larger for nucleic acids
        pack_side_chains=True,
        relaxation="fastrelax",
        validation="af3",
        description="DNA/RNA binder design",
        rationale="LigandMPNN 50.5% vs 34.0% for nucleotides. Favor R/K/H for phosphate contacts.",
    ),

    DesignType.PROTEIN_PROTEIN_BINDER: DesignConfig(
        design_type=DesignType.PROTEIN_PROTEIN_BINDER,
        backbone_tool="rfdiffusion",  # Or bindcraft for end-to-end
        sequence_tool="protein_mpnn",  # NOT ligand_mpnn
        temperature=0.1,
        bias_AA=None,
        omit_AA="C",
        mpnn_weights="soluble",  # BindCraft default
        use_atom_context=False,  # No ligand
        pack_side_chains=False,  # FastRelax handles this
        relaxation="fastrelax",  # Critical for interface packing
        validation="af2_multimer",
        description="Protein-protein binder design",
        rationale="No non-protein molecules, so ProteinMPNN is appropriate. FastRelax for interface optimization.",
    ),

    DesignType.ENZYME_ACTIVE_SITE: DesignConfig(
        design_type=DesignType.ENZYME_ACTIVE_SITE,
        backbone_tool="rfd3",  # motif scaffolding mode
        sequence_tool="ligand_mpnn",
        temperature=0.1,
        bias_AA="A:-1.5",
        omit_AA=None,
        use_atom_context=True,
        pack_side_chains=True,
        relaxation="fastrelax",
        validation="af2",  # + docking for substrate
        auto_detect_fixed=True,  # Fix catalytic residues
        description="Enzyme active site scaffolding",
        rationale="LigandMPNN for substrate awareness. Fix catalytic residues via motif scaffolding.",
    ),

    DesignType.SYMMETRIC_OLIGOMER: DesignConfig(
        design_type=DesignType.SYMMETRIC_OLIGOMER,
        backbone_tool="rfdiffusion",
        use_symmetry=True,
        sequence_tool="protein_mpnn",
        temperature=0.1,
        bias_AA=None,
        omit_AA="C",
        pack_side_chains=False,
        relaxation="fastrelax",
        validation="af2_multimer",
        description="Symmetric homo-oligomer design",
        rationale="Pure protein design, symmetry at inference time.",
    ),

    DesignType.SYMMETRIC_LIGAND_BINDER: DesignConfig(
        design_type=DesignType.SYMMETRIC_LIGAND_BINDER,
        backbone_tool="rfd3",
        use_symmetry=True,
        sequence_tool="ligand_mpnn",
        temperature=0.1,
        bias_AA="A:-1.5,W:1.0,Y:1.0",
        omit_AA="C",
        use_atom_context=True,
        pack_side_chains=True,
        relaxation="optional",
        validation="af2_multimer",  # + docking
        description="Symmetric binder with ligand at interface",
        rationale="RFD3 symmetry + ligand features compose. LigandMPNN for ligand awareness.",
    ),

    DesignType.METAL_MEDIATED_DIMER: DesignConfig(
        design_type=DesignType.METAL_MEDIATED_DIMER,
        backbone_tool="rfd3",
        use_symmetry=True,
        symmetry_type="C2",
        sequence_tool="ligand_mpnn",
        temperature=0.1,
        bias_AA="A:-2.0,H:2.0,E:2.0,D:2.0",  # Heavy bias for coordinating residues
        omit_AA=None,
        use_atom_context=True,
        ligand_cutoff=6.0,
        pack_side_chains=True,
        relaxation="optional",  # May disrupt metal geometry
        validation="af2_multimer",
        auto_detect_fixed=True,
        description="Metal-mediated dimer design",
        rationale="Symmetry + ligand + RASA for zinc dimer case. LigandMPNN critical for metal context.",
    ),

    DesignType.GENERAL_SCAFFOLD: DesignConfig(
        design_type=DesignType.GENERAL_SCAFFOLD,
        backbone_tool="rfdiffusion",
        sequence_tool="protein_mpnn",
        temperature=0.1,
        bias_AA=None,
        omit_AA="C",
        pack_side_chains=False,
        relaxation="optional",
        validation="af2",
        description="General novel fold design",
        rationale="No specific function, ProteinMPNN is sufficient.",
    ),

    DesignType.PEPTIDE_BINDER: DesignConfig(
        design_type=DesignType.PEPTIDE_BINDER,
        backbone_tool="bindcraft",  # End-to-end optimized
        sequence_tool="protein_mpnn",
        temperature=0.1,
        bias_AA=None,
        omit_AA="C",
        mpnn_weights="soluble",
        pack_side_chains=False,  # BindCraft handles internally
        relaxation="built_in",  # BindCraft has built-in refinement
        validation="built_in",
        description="Peptide binder design (BindCraft-style)",
        rationale="BindCraft optimized for protein-protein, handles end-to-end.",
    ),
}


# =============================================================================
# METAL-SPECIFIC PRESETS
# =============================================================================

METAL_PRESETS: Dict[str, Dict[str, Any]] = {
    "zinc": {
        "bias_AA": "A:-2.0,H:2.0,C:1.5,E:1.0,D:1.0",
        "omit_AA": None,
        "coordinating_residues": ["HIS", "CYS", "GLU", "ASP"],
        "coordination_distance": 2.8,
        "coordination_numbers": [4, 5, 6],
    },
    "lanthanide": {
        "bias_AA": "A:-2.0,E:2.5,D:2.5,N:1.0,Q:1.0,C:-2.0",
        "omit_AA": "C",  # Soft donor incompatible
        "coordinating_residues": ["GLU", "ASP", "ASN", "GLN"],
        "coordination_distance": 2.4,
        "coordination_numbers": [8, 9],
    },
    "iron": {
        "bias_AA": "A:-2.0,H:2.0,E:1.5,D:1.5,C:1.0,M:0.5",
        "omit_AA": None,
        "coordinating_residues": ["HIS", "CYS", "GLU", "ASP", "MET"],
        "coordination_distance": 2.2,
        "coordination_numbers": [4, 5, 6],
    },
    "copper": {
        "bias_AA": "A:-2.0,H:2.0,C:1.5,M:1.5,E:1.0,D:1.0",
        "omit_AA": None,
        "coordinating_residues": ["HIS", "CYS", "MET", "GLU", "ASP"],
        "coordination_distance": 2.1,
        "coordination_numbers": [4, 5, 6],
    },
    "calcium": {
        "bias_AA": "A:-2.0,E:2.0,D:2.0,N:1.5,Q:1.0",
        "omit_AA": "C",
        "coordinating_residues": ["GLU", "ASP", "ASN"],
        "coordination_distance": 2.4,
        "coordination_numbers": [6, 7, 8],
    },
}


def get_recommended_config(
    design_type: DesignType,
    metal_type: Optional[str] = None,
    **overrides,
) -> DesignConfig:
    """
    Get recommended configuration for a design type.

    Args:
        design_type: Type of design task
        metal_type: For metal binding, specify metal (e.g., "zinc", "lanthanide")
        **overrides: Override specific config values

    Returns:
        DesignConfig with recommended settings
    """
    # Get base preset
    if design_type not in DESIGN_PRESETS:
        raise ValueError(f"Unknown design type: {design_type}")

    base = DESIGN_PRESETS[design_type]

    # Create copy with overrides
    config_dict = {
        "design_type": base.design_type,
        "backbone_tool": overrides.get("backbone_tool", base.backbone_tool),
        "use_symmetry": overrides.get("use_symmetry", base.use_symmetry),
        "symmetry_type": overrides.get("symmetry_type", base.symmetry_type),
        "sequence_tool": overrides.get("sequence_tool", base.sequence_tool),
        "temperature": overrides.get("temperature", base.temperature),
        "num_sequences": overrides.get("num_sequences", base.num_sequences),
        "bias_AA": overrides.get("bias_AA", base.bias_AA),
        "omit_AA": overrides.get("omit_AA", base.omit_AA),
        "mpnn_weights": overrides.get("mpnn_weights", base.mpnn_weights),
        "use_atom_context": overrides.get("use_atom_context", base.use_atom_context),
        "ligand_cutoff": overrides.get("ligand_cutoff", base.ligand_cutoff),
        "pack_side_chains": overrides.get("pack_side_chains", base.pack_side_chains),
        "number_of_packs": overrides.get("number_of_packs", base.number_of_packs),
        "model_noise_level": overrides.get("model_noise_level", base.model_noise_level),
        "relaxation": overrides.get("relaxation", base.relaxation),
        "validation": overrides.get("validation", base.validation),
        "auto_detect_fixed": overrides.get("auto_detect_fixed", base.auto_detect_fixed),
        "fixed_positions": overrides.get("fixed_positions", base.fixed_positions),
        "coordination_cutoff": overrides.get("coordination_cutoff", base.coordination_cutoff),
        "description": base.description,
        "rationale": base.rationale,
    }

    # Apply metal-specific overrides
    if metal_type and design_type in (DesignType.METAL_BINDING, DesignType.METAL_MEDIATED_DIMER):
        metal_key = metal_type.lower()
        if metal_key in METAL_PRESETS:
            metal = METAL_PRESETS[metal_key]
            if "bias_AA" not in overrides:
                config_dict["bias_AA"] = metal["bias_AA"]
            if "omit_AA" not in overrides:
                config_dict["omit_AA"] = metal["omit_AA"]
            config_dict["coordination_cutoff"] = metal.get("coordination_distance", 3.5)

    return DesignConfig(**config_dict)


def infer_design_type(
    has_ligand: bool = False,
    has_metal: bool = False,
    has_nucleotide: bool = False,
    has_target_protein: bool = False,
    is_symmetric: bool = False,
    has_motif: bool = False,
) -> DesignType:
    """
    Infer design type from characteristics.

    Args:
        has_ligand: Design involves small molecule
        has_metal: Design involves metal ion
        has_nucleotide: Design involves DNA/RNA
        has_target_protein: Design binds to a protein target
        is_symmetric: Design is symmetric oligomer
        has_motif: Design includes catalytic/functional motif

    Returns:
        Inferred DesignType
    """
    if has_metal:
        if is_symmetric:
            return DesignType.METAL_MEDIATED_DIMER
        return DesignType.METAL_BINDING

    if has_nucleotide:
        return DesignType.DNA_RNA_BINDER

    if has_ligand:
        if is_symmetric:
            return DesignType.SYMMETRIC_LIGAND_BINDER
        if has_motif:
            return DesignType.ENZYME_ACTIVE_SITE
        return DesignType.SMALL_MOLECULE_BINDER

    if has_target_protein:
        return DesignType.PROTEIN_PROTEIN_BINDER

    if is_symmetric:
        return DesignType.SYMMETRIC_OLIGOMER

    return DesignType.GENERAL_SCAFFOLD
```

**Step 4: Run test to verify it passes**

Run: `wsl -d Ubuntu -e bash -c "cd /mnt/g/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless && python -m pytest test_design_types.py -v"`
Expected: PASS

**Step 5: Commit**

```bash
git add backend/serverless/design_types.py backend/serverless/test_design_types.py
git commit -m "feat: add design type configuration with tool selection

Implement DesignType enum and DesignConfig dataclass that encapsulates
the decision tree for selecting appropriate backbone, sequence, and
relaxation tools based on design requirements.

Key presets:
- METAL_BINDING: LigandMPNN with H/C/E/D bias (77.5% recovery)
- SMALL_MOLECULE_BINDER: LigandMPNN with aromatic bias
- DNA_RNA_BINDER: LigandMPNN with R/K/H bias
- PROTEIN_PROTEIN_BINDER: ProteinMPNN + FastRelax
- METAL_MEDIATED_DIMER: Symmetry + LigandMPNN

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
"
```

---

## Task 2: Create Metal Coordination Site Detector

**Files:**
- Modify: `G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless\inference_utils.py`
- Test: `G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless\test_metal_coordination.py`

**Step 1: Write the failing test**

```python
# test_metal_coordination.py
"""Tests for metal coordination site detection."""

import pytest


ZINC_PDB = """
ATOM      1  N   HIS A  63      10.000  10.000  10.000  1.00 50.00           N
ATOM      2  CA  HIS A  63      11.000  10.000  10.000  1.00 50.00           C
ATOM      3  NE2 HIS A  63      12.500  10.000   7.230  1.00 50.00           N
ATOM     50  N   CYS B  39      10.000  10.000   0.000  1.00 50.00           N
ATOM     51  CA  CYS B  39      11.000  10.000   0.000  1.00 50.00           C
ATOM     52  SG  CYS B  39      12.300  10.000   3.440  1.00 50.00           S
HETATM  100  ZN  ZN  X   1      12.000  10.000   5.000  1.00 50.00          ZN
END
"""


def test_detect_coordinating_residues():
    """Should detect His and Cys coordinating Zn."""
    from inference_utils import detect_coordinating_residues

    sites = detect_coordinating_residues(ZINC_PDB, cutoff=4.0)

    assert len(sites) == 1
    site = sites[0]
    assert site.metal_code == "ZN"

    # Should find both coordinating residues
    chains_residues = [(r.chain, r.resnum) for r in site.coordinating_residues]
    assert ("A", 63) in chains_residues
    assert ("B", 39) in chains_residues


def test_get_fixed_positions_string():
    """Should generate fixed_positions for LigandMPNN."""
    from inference_utils import detect_coordinating_residues

    sites = detect_coordinating_residues(ZINC_PDB, cutoff=4.0)
    fixed = sites[0].get_fixed_positions_list()

    assert "A63" in fixed
    assert "B39" in fixed
```

**Step 2: Run test to verify it fails**

Run: `wsl -d Ubuntu -e bash -c "cd /mnt/g/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless && python -m pytest test_metal_coordination.py -v"`
Expected: FAIL with "cannot import name 'detect_coordinating_residues'"

**Step 3: Write minimal implementation**

Add to `inference_utils.py` after line 168 (after `detect_metal_in_pdb`):

```python
from dataclasses import dataclass
from typing import Tuple

@dataclass
class CoordinatingResidue:
    """A residue coordinating a metal ion."""
    chain: str
    resnum: int
    resname: str
    atom_name: str
    distance: float


@dataclass
class MetalCoordinationSite:
    """A metal coordination site with coordinating residues."""
    metal_code: str
    metal_coords: Tuple[float, float, float]
    coordinating_residues: List[CoordinatingResidue]

    def get_fixed_positions_string(self) -> str:
        """Generate fixed_positions string for LigandMPNN."""
        positions = sorted(set(
            f"{r.chain}{r.resnum}" for r in self.coordinating_residues
        ))
        return ",".join(positions)

    def get_fixed_positions_list(self) -> List[str]:
        """Generate fixed_positions list for LigandMPNN."""
        return sorted(set(
            f"{r.chain}{r.resnum}" for r in self.coordinating_residues
        ))


# Atoms that typically coordinate metals
COORDINATING_ATOMS = {
    "ND1", "NE2",  # Histidine imidazole
    "SG",          # Cysteine thiolate
    "OD1", "OD2",  # Aspartate carboxylate
    "OE1", "OE2",  # Glutamate carboxylate
    "OG", "OG1",   # Serine/Threonine hydroxyl
    "SD",          # Methionine thioether
    "NZ",          # Lysine amine (rare)
}


def detect_coordinating_residues(
    pdb_content: str,
    cutoff: float = 3.5,
    metal_codes: Optional[List[str]] = None,
) -> List[MetalCoordinationSite]:
    """
    Detect metal coordination sites and their coordinating residues.

    Args:
        pdb_content: PDB file content
        cutoff: Distance cutoff for coordination (Angstroms)
        metal_codes: Specific metals to look for (default: all known)

    Returns:
        List of MetalCoordinationSite objects
    """
    import math

    if metal_codes is None:
        metal_codes = list(METAL_CODES)

    # Parse metal positions
    metal_positions = []
    for line in pdb_content.split('\n'):
        if not line.startswith('HETATM'):
            continue
        try:
            res_name = line[17:20].strip()
            if res_name in metal_codes:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                metal_positions.append((res_name, (x, y, z)))
        except (ValueError, IndexError):
            continue

    if not metal_positions:
        return []

    # Parse protein atoms
    protein_atoms = []
    for line in pdb_content.split('\n'):
        if not line.startswith('ATOM'):
            continue
        try:
            atom_name = line[12:16].strip()
            if atom_name not in COORDINATING_ATOMS:
                continue

            res_name = line[17:20].strip()
            chain_id = line[21]
            res_num = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            protein_atoms.append({
                "atom_name": atom_name,
                "res_name": res_name,
                "chain": chain_id,
                "resnum": res_num,
                "coords": (x, y, z),
            })
        except (ValueError, IndexError):
            continue

    # Find coordinating residues for each metal
    sites = []
    for metal_code, metal_coords in metal_positions:
        coordinating = []

        for atom in protein_atoms:
            dx = atom["coords"][0] - metal_coords[0]
            dy = atom["coords"][1] - metal_coords[1]
            dz = atom["coords"][2] - metal_coords[2]
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)

            if dist <= cutoff:
                coordinating.append(CoordinatingResidue(
                    chain=atom["chain"],
                    resnum=atom["resnum"],
                    resname=atom["res_name"],
                    atom_name=atom["atom_name"],
                    distance=dist,
                ))

        if coordinating:
            sites.append(MetalCoordinationSite(
                metal_code=metal_code,
                metal_coords=metal_coords,
                coordinating_residues=coordinating,
            ))

    return sites
```

**Step 4: Run test to verify it passes**

Run: `wsl -d Ubuntu -e bash -c "cd /mnt/g/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless && python -m pytest test_metal_coordination.py -v"`
Expected: PASS

**Step 5: Commit**

```bash
git add backend/serverless/inference_utils.py backend/serverless/test_metal_coordination.py
git commit -m "feat: add metal coordination site detection

Implement detect_coordinating_residues() to automatically identify
metal-coordinating residues for LigandMPNN fixed_positions parameter.

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
"
```

---

## Task 3: Create Design Orchestrator

**Files:**
- Create: `G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless\design_orchestrator.py`
- Test: `G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless\test_design_orchestrator.py`

**Step 1: Write the failing test**

```python
# test_design_orchestrator.py
"""Tests for DesignOrchestrator."""

import pytest
from design_orchestrator import DesignOrchestrator
from design_types import DesignType


def test_orchestrator_initializes_with_design_type():
    """Should initialize with inferred config."""
    orchestrator = DesignOrchestrator(design_type=DesignType.METAL_BINDING)

    assert orchestrator.config.sequence_tool == "ligand_mpnn"
    assert orchestrator.config.use_atom_context == True


def test_orchestrator_builds_mpnn_request():
    """Should build valid MPNN request."""
    orchestrator = DesignOrchestrator(
        design_type=DesignType.METAL_BINDING,
        metal_type="zinc"
    )

    pdb = "ATOM...HETATM...ZN..."
    request = orchestrator.build_mpnn_request(pdb)

    assert request["model_type"] == "ligand_mpnn"
    assert "A:-" in request["bias_AA"]
    assert request["temperature"] == 0.1


def test_orchestrator_builds_rfd3_config():
    """Should build valid RFD3 config for backbone generation."""
    orchestrator = DesignOrchestrator(
        design_type=DesignType.METAL_MEDIATED_DIMER,
        metal_type="zinc"
    )

    rfd3_config = orchestrator.build_rfd3_config(
        ligand_code="ZN",
        chain_length="60-80",
    )

    assert "symmetry" in rfd3_config
    assert rfd3_config["symmetry"]["type"] == "C2"
    assert rfd3_config["ligand"] == "ZN"


def test_orchestrator_recommends_validation():
    """Should recommend appropriate validation method."""
    # Metal binding
    orchestrator = DesignOrchestrator(design_type=DesignType.METAL_BINDING)
    assert "af2" in orchestrator.get_validation_method()

    # Protein-protein
    orchestrator = DesignOrchestrator(design_type=DesignType.PROTEIN_PROTEIN_BINDER)
    assert "multimer" in orchestrator.get_validation_method()
```

**Step 2: Run test to verify it fails**

Run: `wsl -d Ubuntu -e bash -c "cd /mnt/g/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless && python -m pytest test_design_orchestrator.py -v"`
Expected: FAIL with "ModuleNotFoundError"

**Step 3: Write minimal implementation**

```python
# design_orchestrator.py
"""
Design Orchestrator

High-level orchestrator that selects and configures the appropriate
backbone, sequence, relaxation, and validation tools based on design type.

Implements the decision tree from our comprehensive analysis.
"""

from typing import Optional, List, Dict, Any
from design_types import (
    DesignType,
    DesignConfig,
    get_recommended_config,
    METAL_PRESETS,
)
from inference_utils import detect_coordinating_residues


class DesignOrchestrator:
    """
    Orchestrates protein design workflow with intelligent tool selection.

    Usage:
        orchestrator = DesignOrchestrator(
            design_type=DesignType.METAL_BINDING,
            metal_type="zinc"
        )

        # Build RFD3 config for backbone
        rfd3_config = orchestrator.build_rfd3_config(ligand_code="ZN")

        # Build MPNN request for sequence design
        mpnn_request = orchestrator.build_mpnn_request(pdb_content)

        # Get recommended validation
        validation = orchestrator.get_validation_method()
    """

    def __init__(
        self,
        design_type: DesignType,
        metal_type: Optional[str] = None,
        **config_overrides,
    ):
        """
        Initialize orchestrator with design type and optional overrides.

        Args:
            design_type: Type of design task
            metal_type: For metal designs, specify metal (zinc, lanthanide, etc.)
            **config_overrides: Override specific config values
        """
        self.design_type = design_type
        self.metal_type = metal_type
        self.config = get_recommended_config(
            design_type,
            metal_type=metal_type,
            **config_overrides,
        )

    def build_rfd3_config(
        self,
        ligand_code: Optional[str] = None,
        chain_length: str = "60-80",
        input_pdb_path: Optional[str] = None,
        hotspots: Optional[List[str]] = None,
        **extra_params,
    ) -> Dict[str, Any]:
        """
        Build RFD3 configuration for backbone generation.

        Args:
            ligand_code: Ligand/metal code (e.g., "ZN", "AZO")
            chain_length: Chain length range
            input_pdb_path: Path to input PDB with ligand
            hotspots: Target residues to contact
            **extra_params: Additional RFD3 parameters

        Returns:
            Dict suitable for RFD3 inference
        """
        config = {
            "contig": f"{chain_length}",
        }

        if input_pdb_path:
            config["input"] = input_pdb_path

        # Add ligand conditioning
        if ligand_code:
            config["ligand"] = ligand_code
            # Partially bury ligand at interface
            config["select_partially_buried"] = {ligand_code: "ALL"}

        # Add symmetry if needed
        if self.config.use_symmetry:
            sym_type = self.config.symmetry_type or "C2"
            config["symmetry"] = {"type": sym_type}
            # For symmetric designs, need chain break
            config["contig"] = f"{chain_length},/0,{chain_length}"
            # Center on ligand
            config["ori_token"] = [0.0, 0.0, 0.0]
            config["infer_ori_strategy"] = "com"

        # Add hotspots if specified
        if hotspots:
            config["select_hotspots"] = ",".join(hotspots)

        # Merge extra params
        config.update(extra_params)

        return config

    def build_mpnn_request(
        self,
        pdb_content: str,
        fixed_positions: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        """
        Build MPNN inference request from config.

        Args:
            pdb_content: PDB file content
            fixed_positions: Manual fixed positions (auto-detected if None)

        Returns:
            Dict suitable for run_mpnn_inference()
        """
        request = {
            "pdb_content": pdb_content,
            "num_sequences": self.config.num_sequences,
            "temperature": self.config.temperature,
            "model_type": self.config.sequence_tool,
            "model_noise_level": self.config.model_noise_level,
            "save_stats": True,
        }

        # Add bias/omit
        if self.config.bias_AA:
            request["bias_AA"] = self.config.bias_AA
        if self.config.omit_AA:
            request["omit_AA"] = self.config.omit_AA

        # LigandMPNN-specific parameters
        if self.config.sequence_tool == "ligand_mpnn":
            request["ligand_cutoff_for_score"] = self.config.ligand_cutoff
            request["pack_side_chains"] = self.config.pack_side_chains
            request["number_of_packs_per_design"] = self.config.number_of_packs

        # Handle fixed positions
        if fixed_positions:
            request["fixed_positions"] = fixed_positions
        elif self.config.auto_detect_fixed:
            # Auto-detect from metal coordination
            sites = detect_coordinating_residues(
                pdb_content,
                cutoff=self.config.coordination_cutoff
            )
            if sites:
                detected = sites[0].get_fixed_positions_list()
                if detected:
                    request["fixed_positions"] = detected
                    print(f"[Orchestrator] Auto-detected fixed positions: {detected}")

        return request

    def get_validation_method(self) -> str:
        """
        Get recommended validation method.

        Returns:
            Validation method string
        """
        return self.config.validation

    def should_use_fastrelax(self) -> bool:
        """
        Determine if FastRelax should be used.

        Returns:
            True if FastRelax is recommended
        """
        return self.config.relaxation == "fastrelax"

    def get_workflow_summary(self) -> Dict[str, Any]:
        """
        Get summary of recommended workflow.

        Returns:
            Dict with workflow steps and tools
        """
        return {
            "design_type": self.design_type.name,
            "backbone_tool": self.config.backbone_tool,
            "sequence_tool": self.config.sequence_tool,
            "relaxation": self.config.relaxation,
            "validation": self.config.validation,
            "description": self.config.description,
            "rationale": self.config.rationale,
            "parameters": {
                "temperature": self.config.temperature,
                "bias_AA": self.config.bias_AA,
                "omit_AA": self.config.omit_AA,
                "pack_side_chains": self.config.pack_side_chains,
                "use_symmetry": self.config.use_symmetry,
            }
        }


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def create_metal_design_orchestrator(
    metal_type: str,
    is_dimer: bool = False,
    **overrides,
) -> DesignOrchestrator:
    """
    Create orchestrator for metal-binding design.

    Args:
        metal_type: Metal type (zinc, lanthanide, iron, copper, calcium)
        is_dimer: True for metal-mediated dimer
        **overrides: Additional config overrides

    Returns:
        Configured DesignOrchestrator
    """
    design_type = (
        DesignType.METAL_MEDIATED_DIMER if is_dimer
        else DesignType.METAL_BINDING
    )
    return DesignOrchestrator(
        design_type=design_type,
        metal_type=metal_type,
        **overrides,
    )


def create_small_molecule_orchestrator(
    ligand_smiles: Optional[str] = None,
    favor_aromatics: bool = True,
    **overrides,
) -> DesignOrchestrator:
    """
    Create orchestrator for small molecule binder design.

    Args:
        ligand_smiles: SMILES string for bias optimization
        favor_aromatics: Bias toward aromatic residues
        **overrides: Additional config overrides

    Returns:
        Configured DesignOrchestrator
    """
    # Adjust bias based on ligand chemistry
    bias = "A:-1.0"
    if favor_aromatics:
        bias += ",W:1.5,Y:1.5,F:1.0"

    return DesignOrchestrator(
        design_type=DesignType.SMALL_MOLECULE_BINDER,
        bias_AA=overrides.pop("bias_AA", bias),
        **overrides,
    )


def create_protein_binder_orchestrator(**overrides) -> DesignOrchestrator:
    """
    Create orchestrator for protein-protein binder design.

    Args:
        **overrides: Config overrides

    Returns:
        Configured DesignOrchestrator
    """
    return DesignOrchestrator(
        design_type=DesignType.PROTEIN_PROTEIN_BINDER,
        **overrides,
    )
```

**Step 4: Run test to verify it passes**

Run: `wsl -d Ubuntu -e bash -c "cd /mnt/g/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless && python -m pytest test_design_orchestrator.py -v"`
Expected: PASS

**Step 5: Commit**

```bash
git add backend/serverless/design_orchestrator.py backend/serverless/test_design_orchestrator.py
git commit -m "feat: add DesignOrchestrator for intelligent workflow selection

Implement high-level orchestrator that builds RFD3 configs, MPNN requests,
and validation pipelines based on design type.

Features:
- Auto-selects LigandMPNN vs ProteinMPNN based on design type
- Auto-detects metal coordination sites for fixed_positions
- Configures symmetry, RASA, and ligand conditioning
- Recommends validation method (AF2, AF2-Multimer, ESMFold)

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
"
```

---

## Task 4: Add Unified Design API Endpoint

**Files:**
- Modify: `G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless\handler.py`
- Test: `G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless\test_unified_design.py`

**Step 1: Write the failing test**

```python
# test_unified_design.py
"""Tests for unified design API endpoint."""

import pytest


def test_unified_design_task_recognized():
    """Handler should recognize 'design' task."""
    # Import after handler module loads
    pass  # Placeholder - actual test requires handler import


def test_design_request_infers_type():
    """Should infer design type from request parameters."""
    from handler import infer_design_type_from_request

    # Metal binding
    request = {"metal": "ZN", "pdb_content": "..."}
    design_type = infer_design_type_from_request(request)
    assert design_type.name == "METAL_BINDING"

    # Protein binder
    request = {"target_pdb": "...", "pdb_content": "..."}
    design_type = infer_design_type_from_request(request)
    assert design_type.name == "PROTEIN_PROTEIN_BINDER"

    # Metal dimer
    request = {"metal": "ZN", "symmetry": "C2", "pdb_content": "..."}
    design_type = infer_design_type_from_request(request)
    assert design_type.name == "METAL_MEDIATED_DIMER"
```

**Step 2: Run test to verify it fails**

Run: `wsl -d Ubuntu -e bash -c "cd /mnt/g/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless && python -m pytest test_unified_design.py -v"`
Expected: FAIL

**Step 3: Write minimal implementation**

Add to `handler.py`:

```python
# Add to imports
from design_types import DesignType, infer_design_type
from design_orchestrator import DesignOrchestrator

# Add to task dispatch (around line 340)
elif task == "design":
    return handle_unified_design(job_input)

def infer_design_type_from_request(request: Dict[str, Any]) -> DesignType:
    """
    Infer design type from request parameters.

    Args:
        request: API request dict

    Returns:
        Inferred DesignType
    """
    has_metal = "metal" in request
    has_ligand = "ligand" in request or "ligand_smiles" in request
    has_nucleotide = "dna" in request or "rna" in request
    has_target = "target_pdb" in request
    is_symmetric = "symmetry" in request
    has_motif = "motif" in request or "catalytic_residues" in request

    return infer_design_type(
        has_ligand=has_ligand,
        has_metal=has_metal,
        has_nucleotide=has_nucleotide,
        has_target_protein=has_target,
        is_symmetric=is_symmetric,
        has_motif=has_motif,
    )


def handle_unified_design(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Unified design endpoint with intelligent tool selection.

    This endpoint automatically selects the appropriate backbone,
    sequence design, relaxation, and validation tools based on
    the design requirements.

    Input:
        # Core parameters
        pdb_content: str (required) - Input PDB content

        # Design type inference (any combination)
        metal: str (optional) - Metal ion code (ZN, FE, TB, etc.)
        ligand: str (optional) - Ligand CCD code
        target_pdb: str (optional) - Target protein for binder design
        symmetry: str (optional) - Symmetry type (C2, C3, D2, etc.)

        # Overrides
        design_type: str (optional) - Force specific design type
        sequence_tool: str (optional) - Force ligand_mpnn or protein_mpnn
        temperature: float (optional) - Override temperature
        num_sequences: int (optional) - Override sequence count
        bias_AA: str (optional) - Override AA bias

    Returns:
        {
            "status": "completed",
            "design_type": str,
            "workflow": {...},  # Tool selection summary
            "sequences": [...],
            "validation": {...}
        }
    """
    pdb_content = job_input.get("pdb_content")
    if not pdb_content:
        return {"status": "failed", "error": "Missing pdb_content"}

    # Infer or use explicit design type
    if "design_type" in job_input:
        try:
            design_type = DesignType[job_input["design_type"].upper()]
        except KeyError:
            return {
                "status": "failed",
                "error": f"Unknown design_type: {job_input['design_type']}"
            }
    else:
        design_type = infer_design_type_from_request(job_input)

    # Get metal type if applicable
    metal_type = job_input.get("metal", "").lower() or None
    if metal_type:
        # Map common metal codes to preset names
        metal_map = {
            "zn": "zinc", "fe": "iron", "cu": "copper", "ca": "calcium",
            "tb": "lanthanide", "gd": "lanthanide", "eu": "lanthanide",
        }
        metal_type = metal_map.get(metal_type, metal_type)

    # Build config overrides
    overrides = {}
    if "temperature" in job_input:
        overrides["temperature"] = job_input["temperature"]
    if "num_sequences" in job_input:
        overrides["num_sequences"] = job_input["num_sequences"]
    if "bias_AA" in job_input:
        overrides["bias_AA"] = job_input["bias_AA"]
    if "sequence_tool" in job_input:
        overrides["sequence_tool"] = job_input["sequence_tool"]

    # Create orchestrator
    orchestrator = DesignOrchestrator(
        design_type=design_type,
        metal_type=metal_type,
        **overrides,
    )

    # Build and run MPNN request
    mpnn_request = orchestrator.build_mpnn_request(
        pdb_content,
        fixed_positions=job_input.get("fixed_positions"),
    )

    # Run sequence design
    mpnn_result = run_mpnn_inference(**mpnn_request)

    if mpnn_result.get("status") != "completed":
        return mpnn_result

    # Build response
    result = {
        "status": "completed",
        "design_type": design_type.name,
        "workflow": orchestrator.get_workflow_summary(),
        "sequences": mpnn_result.get("result", {}).get("sequences", []),
        "config_used": {
            "sequence_tool": orchestrator.config.sequence_tool,
            "temperature": orchestrator.config.temperature,
            "bias_AA": orchestrator.config.bias_AA,
            "fixed_positions": mpnn_request.get("fixed_positions"),
        },
    }

    # Add packed structures if available
    packed = mpnn_result.get("result", {}).get("packed_structures")
    if packed:
        result["packed_structures"] = packed

    return result
```

**Step 4: Run test to verify it passes**

Run: `wsl -d Ubuntu -e bash -c "cd /mnt/g/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless && python -m pytest test_unified_design.py -v"`
Expected: PASS

**Step 5: Commit**

```bash
git add backend/serverless/handler.py backend/serverless/test_unified_design.py
git commit -m "feat: add unified 'design' endpoint with auto tool selection

Add handle_unified_design() that automatically selects appropriate tools
based on design requirements (metal, ligand, symmetry, target protein).

Usage:
  {\"task\": \"design\", \"metal\": \"ZN\", \"pdb_content\": \"...\"}
  -> Auto-selects LigandMPNN with zinc preset

  {\"task\": \"design\", \"target_pdb\": \"...\", \"pdb_content\": \"...\"}
  -> Auto-selects ProteinMPNN + FastRelax

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
"
```

---

## Task 5: Create Comprehensive Integration Tests

**Files:**
- Create: `G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless\test_design_integration.py`

**Step 1: Write integration tests**

```python
# test_design_integration.py
"""
Integration tests for the unified design system.

Run with: pytest test_design_integration.py -v -m integration
"""

import pytest
import requests
from pathlib import Path
from collections import Counter

API_URL = "http://localhost:8000/runsync"

pytestmark = pytest.mark.integration


def load_zinc_pdb():
    """Load zinc dimer test PDB."""
    paths = [
        Path(__file__).parent.parent.parent / "design_zn_dimer.pdb",
        Path("/mnt/g/Github_local_repo/Banta_Lab_RFdiffusion/design_zn_dimer.pdb"),
    ]
    for p in paths:
        if p.exists():
            return p.read_text()
    pytest.skip("Test PDB not found")


def parse_sequences(result):
    """Extract sequence strings from result."""
    sequences = []
    for seq in result.get("sequences", []):
        if isinstance(seq, dict):
            for line in seq.get("content", "").split("\n"):
                if not line.startswith(">"):
                    sequences.append(line.strip())
    return sequences


class TestMetalBindingDesign:
    """Tests for metal-binding design."""

    def test_unified_design_with_metal_uses_ligandmpnn(self):
        """Unified design with metal should use LigandMPNN."""
        pdb = load_zinc_pdb()

        payload = {
            "input": {
                "task": "design",
                "metal": "ZN",
                "pdb_content": pdb,
                "num_sequences": 2,
            }
        }

        response = requests.post(API_URL, json=payload, timeout=300)
        result = response.json()

        assert "error" not in result
        output = result.get("output", {})

        # Should have selected LigandMPNN
        assert output.get("workflow", {}).get("sequence_tool") == "ligand_mpnn"
        assert output.get("design_type") == "METAL_BINDING"

    def test_zinc_design_reduces_alanine(self):
        """Zinc design should have low alanine content."""
        pdb = load_zinc_pdb()

        payload = {
            "input": {
                "task": "design",
                "metal": "ZN",
                "pdb_content": pdb,
                "num_sequences": 4,
            }
        }

        response = requests.post(API_URL, json=payload, timeout=300)
        result = response.json()
        output = result.get("output", {})
        sequences = parse_sequences(output)

        combined = "".join(sequences)
        total = len(combined)
        ala_pct = (combined.count("A") / total * 100) if total > 0 else 0

        print(f"Alanine: {ala_pct:.1f}%")
        assert ala_pct < 25, f"Alanine too high: {ala_pct:.1f}%"

    def test_zinc_design_has_coordinating_residues(self):
        """Zinc design should have H/C/E/D residues."""
        pdb = load_zinc_pdb()

        payload = {
            "input": {
                "task": "design",
                "metal": "ZN",
                "pdb_content": pdb,
                "num_sequences": 4,
            }
        }

        response = requests.post(API_URL, json=payload, timeout=300)
        result = response.json()
        output = result.get("output", {})
        sequences = parse_sequences(output)

        combined = "".join(sequences)
        aa_counts = Counter(combined)

        coordinating = sum(aa_counts.get(aa, 0) for aa in "HCED")
        assert coordinating > 0, "Should have coordinating residues"


class TestProteinBinderDesign:
    """Tests for protein-protein binder design."""

    def test_protein_binder_uses_proteinmpnn(self):
        """Protein binder should use ProteinMPNN."""
        # Minimal test - just verify tool selection
        from handler import infer_design_type_from_request
        from design_types import DesignType

        request = {"target_pdb": "..."}
        design_type = infer_design_type_from_request(request)

        assert design_type == DesignType.PROTEIN_PROTEIN_BINDER


class TestSymmetricDesign:
    """Tests for symmetric design."""

    def test_metal_with_symmetry_creates_dimer(self):
        """Metal + symmetry should create metal-mediated dimer."""
        from handler import infer_design_type_from_request
        from design_types import DesignType

        request = {"metal": "ZN", "symmetry": "C2"}
        design_type = infer_design_type_from_request(request)

        assert design_type == DesignType.METAL_MEDIATED_DIMER


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
```

**Step 2: Commit**

```bash
git add backend/serverless/test_design_integration.py
git commit -m "test: add comprehensive integration tests for design system

Test coverage:
- Unified design auto-selects correct tools
- Metal binding uses LigandMPNN and reduces alanine
- Protein binder uses ProteinMPNN
- Symmetry + metal creates dimer design

Run with: pytest test_design_integration.py -v -m integration

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
"
```

---

## Task 6: Update Documentation

**Files:**
- Create: `G:\Github_local_repo\Banta_Lab_RFdiffusion\docs\DESIGN_TOOL_SELECTION.md`

**Step 1: Write comprehensive documentation**

```markdown
# Protein Design Tool Selection Guide

## Quick Reference Card

| I want to design... | Backbone Tool | Sequence Tool | Relaxation | Validation |
|---------------------|---------------|---------------|------------|------------|
| Metal-binding protein | RFD3 | LigandMPNN | Optional | AF2 + Metal validation |
| Small molecule binder | RFD3 | LigandMPNN | Optional | AF2 + Docking |
| DNA/RNA binder | RFD3 | LigandMPNN | FastRelax | AF3 |
| Enzyme with active site | RFD3 (motif) | LigandMPNN | FastRelax | AF2 + Docking |
| Protein-protein binder | RFDiffusion | ProteinMPNN | FastRelax | AF2-Multimer |
| Symmetric oligomer | RFDiffusion | ProteinMPNN | FastRelax | AF2-Multimer |
| Symmetric ligand binder | RFD3 + symmetry | LigandMPNN | Optional | AF2-Multimer |
| Metal-mediated dimer | RFD3 + symmetry | LigandMPNN | Optional | AF2-Multimer |
| General scaffold | RFDiffusion | ProteinMPNN | Optional | AF2 |
| Peptide binder | BindCraft | ProteinMPNN | Built-in | Built-in |

## The One Rule

**If your design involves ANY non-protein atoms (small molecules, metals, nucleotides), use LigandMPNN.**

The performance difference is substantial:
- Metals: 77.5% vs 40.6% sequence recovery
- Small molecules: 63.3% vs 50.5%
- Nucleotides: 50.5% vs 34.0%

## Using the Unified Design API

```python
# Automatic tool selection based on parameters
response = requests.post(API_URL, json={
    "input": {
        "task": "design",
        "metal": "ZN",           # Triggers LigandMPNN + metal presets
        "pdb_content": pdb,
        "num_sequences": 8,
    }
})

# Explicit design type
response = requests.post(API_URL, json={
    "input": {
        "task": "design",
        "design_type": "METAL_MEDIATED_DIMER",
        "metal": "ZN",
        "pdb_content": pdb,
    }
})
```

## Metal-Specific Presets

| Metal | bias_AA | omit_AA | Coordinating Residues |
|-------|---------|---------|----------------------|
| Zinc | A:-2.0,H:2.0,C:1.5,E:1.0,D:1.0 | - | HIS, CYS, GLU, ASP |
| Lanthanide | A:-2.0,E:2.5,D:2.5,N:1.0,Q:1.0,C:-2.0 | C | GLU, ASP, ASN, GLN |
| Iron | A:-2.0,H:2.0,E:1.5,D:1.5,C:1.0 | - | HIS, CYS, GLU, ASP, MET |
| Copper | A:-2.0,H:2.0,C:1.5,M:1.5,E:1.0,D:1.0 | - | HIS, CYS, MET, GLU, ASP |

## When FastRelax is Needed

**Use FastRelax:**
- After ProteinMPNN (no sidechain packing)
- Protein-protein interfaces
- RFDiffusion backbones (may have strain)

**Skip FastRelax:**
- Using LigandMPNN with pack_side_chains=True
- Metal coordination sites (may disrupt geometry)
- Quick screening

## Complete Workflow: Metal-Binding Dimer

```
1. BACKBONE GENERATION (RFD3)
   - symmetry: {type: C2}
   - ligand: "ZN"
   - select_partially_buried: {ZN: ALL}
   - ori_token: [0, 0, 0]

2. INTERFACE ANALYSIS
   - Auto-detect coordinating residues
   - Set fixed_positions for A63, B39

3. SEQUENCE DESIGN (LigandMPNN)
   - model_type: "ligand_mpnn"
   - bias_AA: "A:-2.0,H:2.0,E:1.0,D:1.0"
   - fixed_positions: ["A63", "B39"]
   - temperature: 0.1
   - pack_side_chains: True

4. VALIDATION
   - AF2-Multimer: Check dimer folds
   - Metal geometry: Verify coordination distances
   - pLDDT > 70, pAE < 10
```
```

**Step 2: Commit**

```bash
git add docs/DESIGN_TOOL_SELECTION.md
git commit -m "docs: add comprehensive tool selection guide

Document the decision tree for selecting backbone, sequence, relaxation,
and validation tools based on design requirements.

Includes:
- Quick reference card for all design types
- Metal-specific presets
- FastRelax usage guidelines
- Complete workflow examples

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
"
```

---

## Summary

This comprehensive plan implements:

1. **Design Type Configuration** (Task 1): Enum and presets for all design scenarios
2. **Metal Coordination Detection** (Task 2): Auto-detect fixed_positions for LigandMPNN
3. **Design Orchestrator** (Task 3): Intelligent tool selection based on design type
4. **Unified API Endpoint** (Task 4): Single `task: "design"` endpoint with auto-selection
5. **Integration Tests** (Task 5): Validate tool selection and alanine reduction
6. **Documentation** (Task 6): Complete guide for users

### Key Technical Decisions

| Decision | Rationale |
|----------|-----------|
| LigandMPNN for metals | 77.5% vs 40.6% sequence recovery |
| `bias_AA: "A:-2.0"` | Reduces alanine from 60% to <10% |
| Auto-detect fixed_positions | Preserves metal coordination geometry |
| Skip FastRelax for metals | May disrupt coordination geometry |
| RFD3 for ligand designs | Features compose (symmetry + ligand + RASA) |

### Execution Handoff

**Plan complete and saved to `docs/plans/2026-01-16-comprehensive-protein-design-toolkit.md`.**

**Two execution options:**

1. **Subagent-Driven (this session)** - I dispatch fresh subagent per task, review between tasks, fast iteration

2. **Parallel Session (separate)** - Open new session with executing-plans, batch execution with checkpoints

**Which approach?**
