# design_types.py
"""
Protein Design Type Configuration

Implements the decision tree for selecting appropriate tools based on design type.
Based on:
- LigandMPNN paper (Nature Methods 2025)
- RFD3 documentation
- BindCraft source code analysis
"""

import dataclasses
import re
import warnings
from enum import Enum, auto
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Tuple

from metal_chemistry import get_amino_acid_bias, get_coordination_number_range


# Valid single-letter amino acid codes
VALID_AA_CODES = set("ACDEFGHIKLMNPQRSTVWY")


def validate_bias_AA(bias_str: str) -> Tuple[bool, Optional[str]]:
    """
    Validate bias_AA format string.

    Expected format: "A:-2.0,H:2.0,E:1.0" (comma-separated AA:value pairs)

    Args:
        bias_str: The bias_AA string to validate

    Returns:
        Tuple of (is_valid, error_message)
        error_message is None if valid
    """
    if not bias_str:
        return True, None

    # Pattern: single letter AA code, colon, optional minus, number with optional decimal
    pattern = r'^[A-Z]:-?\d+\.?\d*$'

    parts = bias_str.split(',')
    for part in parts:
        part = part.strip()
        if not part:
            continue

        if not re.match(pattern, part):
            return False, f"Invalid bias_AA format: '{part}'. Expected format like 'A:-2.0' or 'H:2.0'"

        aa_code = part[0]
        if aa_code not in VALID_AA_CODES:
            return False, f"Invalid amino acid code: '{aa_code}'. Valid codes: {sorted(VALID_AA_CODES)}"

    return True, None


def validate_omit_AA(omit_str: str) -> Tuple[bool, Optional[str]]:
    """
    Validate omit_AA format string.

    Expected format: "C" or "CM" (single letter amino acid codes)

    Args:
        omit_str: The omit_AA string to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not omit_str:
        return True, None

    for aa in omit_str:
        if aa not in VALID_AA_CODES:
            return False, f"Invalid amino acid code in omit_AA: '{aa}'. Valid codes: {sorted(VALID_AA_CODES)}"

    return True, None


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

# Metal presets with HSAB-compliant chemistry
METAL_PRESETS: Dict[str, Dict[str, Any]] = {
    # Zinc (Borderline acid - prefers S, N donors)
    "zinc": {
        "description": "Zinc binding (general)",
        "bias_AA": "C:2.5,H:2.0,E:1.0,D:1.0,A:-2.0",
        "coordination_numbers": [4, 5, 6],
        "preferred_residues": ["C", "H", "D", "E"],
    },
    "zinc_structural": {
        "description": "Zinc structural site (Cys4 or Cys3His)",
        "bias_AA": get_amino_acid_bias("ZN", 2, site_type="structural"),
        "coordination_numbers": [4],
        "preferred_residues": ["C", "H"],
    },
    "zinc_catalytic": {
        "description": "Zinc catalytic site (His-rich + water)",
        "bias_AA": "H:3.0,E:2.0,D:2.0,C:1.0,A:-2.0",
        "coordination_numbers": [4, 5],
        "preferred_residues": ["H", "E", "D"],
    },
    # Lanthanides (Hard acids - O donors only, EXCLUDE Cys)
    "lanthanide": {
        "description": "Lanthanide binding (Tb, Eu, Gd, etc.)",
        "bias_AA": get_amino_acid_bias("TB", 3),
        "coordination_numbers": [8, 9],
        "preferred_residues": ["E", "D", "N", "Q"],
        "excluded_residues": ["C"],
    },
    "terbium": {
        "description": "Terbium binding site",
        "bias_AA": get_amino_acid_bias("TB", 3),
        "coordination_numbers": [8, 9],
        "geometry": "square_antiprism",
    },
    "europium": {
        "description": "Europium binding site",
        "bias_AA": get_amino_acid_bias("EU", 3),
        "coordination_numbers": [8, 9],
    },
    "gadolinium": {
        "description": "Gadolinium binding site (MRI)",
        "bias_AA": get_amino_acid_bias("GD", 3),
        "coordination_numbers": [8, 9],
    },
    # Iron (oxidation state dependent)
    "iron": {
        "description": "Iron binding (default Fe2+)",
        "bias_AA": get_amino_acid_bias("FE", 2),
        "coordination_numbers": [4, 5, 6],
    },
    "iron_ii": {
        "description": "Iron(II) - borderline, accepts S/N/O",
        "bias_AA": get_amino_acid_bias("FE", 2),
        "coordination_numbers": [4, 5, 6],
    },
    "iron_iii": {
        "description": "Iron(III) - harder, prefers O/N over S",
        "bias_AA": get_amino_acid_bias("FE", 3),
        "coordination_numbers": [5, 6],
    },
    "iron_sulfur": {
        "description": "Iron-sulfur cluster (FeS)",
        "bias_AA": "C:3.5,H:1.5,A:-2.0",
        "coordination_numbers": [4],
    },
    # Copper (oxidation state dependent)
    "copper": {
        "description": "Copper binding (default Cu2+)",
        "bias_AA": get_amino_acid_bias("CU", 2),
        "coordination_numbers": [4, 5, 6],
    },
    "copper_type1": {
        "description": "Type I blue copper (Cu+, 2His+Cys+Met)",
        "bias_AA": "C:3.0,H:2.5,M:2.0,A:-2.0",
        "coordination_numbers": [4],
        "histidine_tautomer": "ND1",
    },
    "copper_type2": {
        "description": "Type II copper (Cu2+, His-rich)",
        "bias_AA": "H:3.0,E:1.5,D:1.5,C:0.5,A:-2.0",
        "coordination_numbers": [4, 5, 6],
        "histidine_tautomer": "NE2",
    },
    # Other metals
    "calcium": {
        "description": "Calcium binding (EF-hand style)",
        "bias_AA": get_amino_acid_bias("CA", 2),
        "coordination_numbers": [6, 7, 8],
    },
    "magnesium": {
        "description": "Magnesium binding (octahedral)",
        "bias_AA": get_amino_acid_bias("MG", 2),
        "coordination_numbers": [6],
    },
    "manganese": {
        "description": "Manganese binding",
        "bias_AA": get_amino_acid_bias("MN", 2),
        "coordination_numbers": [5, 6],
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
    # Validate override keys
    valid_keys = {f.name for f in dataclasses.fields(DesignConfig)}
    invalid_keys = set(overrides.keys()) - valid_keys
    if invalid_keys:
        warnings.warn(f"Unknown config keys ignored: {invalid_keys}. Valid keys: {valid_keys}")

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
            if "omit_AA" not in overrides and "omit_AA" in metal:
                config_dict["omit_AA"] = metal["omit_AA"]
            config_dict["coordination_cutoff"] = metal.get("coordination_distance", 3.5)
        else:
            warnings.warn(f"Unknown metal type '{metal_type}'. Valid types: {list(METAL_PRESETS.keys())}")

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

    Priority order (first match wins):
    1. Metal (with symmetric -> METAL_MEDIATED_DIMER, else METAL_BINDING)
    2. Nucleotide -> DNA_RNA_BINDER
    3. Ligand (with symmetric -> SYMMETRIC_LIGAND_BINDER, with motif -> ENZYME_ACTIVE_SITE, else SMALL_MOLECULE_BINDER)
    4. Target protein -> PROTEIN_PROTEIN_BINDER
    5. Symmetric -> SYMMETRIC_OLIGOMER
    6. Default -> GENERAL_SCAFFOLD

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
