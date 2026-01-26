"""
Design Decision Rules for RFD3 + LigandMPNN Configuration.

This file provides structured decision tables for AI-driven protein design.
Rules are derived from successful experiment scripts (Dy_TriNOx v11, ln_citrate r7b).

USAGE: Import DESIGN_RULES and pass to Claude for config generation context.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional
from enum import Enum


# ═══════════════════════════════════════════════════════════════════
# 1. CHAIN LENGTH RULES
# ═══════════════════════════════════════════════════════════════════

class LigandSize(Enum):
    SMALL = "small"      # <10 heavy atoms (single metal ion, small molecules)
    MEDIUM = "medium"    # 10-30 heavy atoms (citrate, ATP, cofactors)
    LARGE = "large"      # >30 heavy atoms (PQQ, TriNOx, heme)


CHAIN_LENGTH_RULES: Dict[LigandSize, Dict[str, str]] = {
    LigandSize.SMALL: {
        "base_range": "80-100",
        "with_burial": "100-120",     # Need more residues to bury
        "symmetric": "60-80",          # Per subunit
    },
    LigandSize.MEDIUM: {
        "base_range": "100-130",
        "with_burial": "120-150",
        "symmetric": "80-100",
    },
    LigandSize.LARGE: {
        "base_range": "120-150",
        "with_burial": "140-180",
        "symmetric": "100-120",
    },
}


# ═══════════════════════════════════════════════════════════════════
# 2. HOTSPOT SELECTION RULES
# ═══════════════════════════════════════════════════════════════════

@dataclass
class HotspotStrategy:
    """Strategy for selecting hotspot atoms."""
    name: str
    description: str
    atom_selection: str  # "all", "donors", "hydrophobic", "coordinating"


HOTSPOT_STRATEGIES: Dict[str, HotspotStrategy] = {
    "hydrophobic_ligand": HotspotStrategy(
        name="hydrophobic_packing",
        description="Select aromatic C and aliphatic C for hydrophobic packing",
        atom_selection="hydrophobic",
    ),
    "polar_ligand": HotspotStrategy(
        name="donor_contacts",
        description="Select donor atoms (O, N) for H-bond network",
        atom_selection="donors",
    ),
    "metal_coordination": HotspotStrategy(
        name="coordination_sphere",
        description="Select metal-coordinating atoms only",
        atom_selection="coordinating",
    ),
    "full_coverage": HotspotStrategy(
        name="all_atoms",
        description="Hotspot all ligand atoms (simple default)",
        atom_selection="all",
    ),
}


def select_hotspot_strategy(ligand_character: str, has_metal: bool) -> str:
    """Select hotspot strategy based on ligand properties."""
    if ligand_character == "hydrophobic":
        return "hydrophobic_ligand"
    elif has_metal:
        return "metal_coordination"
    else:
        return "polar_ligand"


# ═══════════════════════════════════════════════════════════════════
# 3. BURIAL STRATEGY RULES
# ═══════════════════════════════════════════════════════════════════

@dataclass
class BurialStrategy:
    """Strategy for burial conditioning."""
    target: str           # "ligand", "metal", "both", "partial"
    atom_selection: str   # "all", "hydrophobic", "core"
    cfg_boost: float      # Additional CFG scale needed


BURIAL_STRATEGIES: Dict[str, BurialStrategy] = {
    "hydrophobic_full": BurialStrategy(
        target="ligand",
        atom_selection="all",
        cfg_boost=0.5,
    ),
    "polar_partial": BurialStrategy(
        target="ligand",
        atom_selection="hydrophobic",  # Only bury hydrophobic portions
        cfg_boost=0.0,
    ),
    "metal_focus": BurialStrategy(
        target="metal",
        atom_selection="all",
        cfg_boost=0.3,
    ),
    "coordination_pocket": BurialStrategy(
        target="both",
        atom_selection="all",
        cfg_boost=0.5,
    ),
}


def select_burial_strategy(
    ligand_character: str,
    metal_type: Optional[str],
    design_goal: str
) -> str:
    """Select burial strategy based on context."""
    if ligand_character == "hydrophobic":
        return "hydrophobic_full"
    elif metal_type and design_goal == "catalysis":
        return "coordination_pocket"
    elif metal_type:
        return "metal_focus"
    else:
        return "polar_partial"


# ═══════════════════════════════════════════════════════════════════
# 4. ORIENTATION STRATEGY RULES
# ═══════════════════════════════════════════════════════════════════

ORIENTATION_RULES: Dict[str, Optional[str]] = {
    # Design goal → infer_ori_strategy
    "binding": "com",           # Center on ligand center-of-mass
    "catalysis": "com",         # Center on active site
    "sensing": "hotspots",      # Approach from specific interface
    "scaffolding": None,        # Let RFD3 decide (motif scaffolding)
    "symmetric": None,          # Symmetry handles orientation
}


# ═══════════════════════════════════════════════════════════════════
# 5. MPNN AMINO ACID BIAS RULES
# ═══════════════════════════════════════════════════════════════════

@dataclass
class MPNNBiasProfile:
    """MPNN amino acid bias profile."""
    name: str
    bias_AA: str              # Comma-separated bias string
    omit_AA: Optional[str]    # Amino acids to exclude
    temperature: float        # Sampling temperature
    rationale: str


# Metal-specific profiles (from HSAB theory)
MPNN_METAL_PROFILES: Dict[str, MPNNBiasProfile] = {
    "hard": MPNNBiasProfile(
        name="lanthanide_coordination",
        bias_AA="E:3.0,D:3.0,N:2.4,Q:2.4",
        omit_AA="C",  # Cysteine incompatible with hard acids
        temperature=0.1,
        rationale="Hard acid: prefers O donors (carboxylates)",
    ),
    "soft": MPNNBiasProfile(
        name="soft_metal",
        bias_AA="C:4.0,H:3.0,M:2.0",
        omit_AA=None,
        temperature=0.1,
        rationale="Soft acid: prefers S/N donors",
    ),
    "borderline": MPNNBiasProfile(
        name="borderline_metal",
        bias_AA="H:3.0,D:2.0,E:2.0,C:1.5",
        omit_AA=None,
        temperature=0.1,
        rationale="Borderline: accepts O/N/S donors",
    ),
    "no_metal": MPNNBiasProfile(
        name="no_bias",
        bias_AA="",
        omit_AA=None,
        temperature=0.2,  # Higher diversity without metal constraints
        rationale="No metal coordination requirements",
    ),
}


# Context-aware bias adjustments
MPNN_CONTEXT_ADJUSTMENTS: Dict[str, Dict[str, float]] = {
    "hydrophobic_burial": {
        # Boost hydrophobic residues for packing around hydrophobic ligand
        "V": 1.5, "L": 1.5, "I": 1.5, "F": 2.0, "W": 1.5, "Y": 1.0,
    },
    "polar_interface": {
        # Boost polar residues for H-bond network
        "S": 1.0, "T": 1.0, "N": 1.0, "Q": 1.0,
    },
}


# ═══════════════════════════════════════════════════════════════════
# 6. CFG SCALE RULES
# ═══════════════════════════════════════════════════════════════════

CFG_SCALE_RULES: Dict[str, float] = {
    "base": 2.0,                    # Default CFG scale
    "with_hbond_conditioning": 2.5, # When using H-bond acceptor/donor
    "with_burial": 2.5,             # When using select_buried
    "complex_ligand": 2.8,          # Large ligand (>30 atoms)
    "symmetric": 3.0,               # Symmetric assemblies need stronger
    "max": 4.0,                     # Upper limit
}

STEP_SCALE_DEFAULT = 1.5  # From successful experiments


# ═══════════════════════════════════════════════════════════════════
# 7. DECISION ENGINE - Main Entry Point
# ═══════════════════════════════════════════════════════════════════

@dataclass
class DesignDecisions:
    """All decisions for a design run."""
    chain_length: str
    hotspot_strategy: str
    hotspot_atoms: List[str]
    burial_strategy: str
    burial_target: Dict[str, str]
    infer_ori_strategy: Optional[str]
    cfg_scale: float
    step_scale: float
    mpnn_bias_AA: str
    mpnn_omit_AA: Optional[str]
    mpnn_temperature: float
    rationale: Dict[str, str] = field(default_factory=dict)


def make_design_decisions(
    ligand_size: LigandSize,
    ligand_character: str,        # "hydrophobic" | "polar" | "amphipathic"
    ligand_atoms: Dict[str, List[str]],  # {"hydrophobic": ["C1"...], "donors": ["O1"...]}
    metal_type: Optional[str],
    metal_hsab: Optional[str],    # "hard" | "soft" | "borderline"
    design_goal: str,             # "binding" | "catalysis" | "sensing"
    target_topology: str,         # "monomer" | "dimer" | "symmetric"
    has_hbond_conditioning: bool,
) -> DesignDecisions:
    """
    Central decision function that applies all rules.
    Returns complete DesignDecisions for RFD3 + MPNN configuration.
    """
    # 1. Chain length
    length_rule = CHAIN_LENGTH_RULES[ligand_size]
    if target_topology == "symmetric":
        chain_length = length_rule["symmetric"]
    elif ligand_character == "hydrophobic":
        chain_length = length_rule["with_burial"]
    else:
        chain_length = length_rule["base_range"]

    # 2. Hotspot strategy
    hotspot_strat = select_hotspot_strategy(ligand_character, metal_type is not None)
    strategy = HOTSPOT_STRATEGIES[hotspot_strat]

    # Select atoms based on strategy
    if strategy.atom_selection == "hydrophobic":
        hotspot_atoms = ligand_atoms.get("hydrophobic", [])
    elif strategy.atom_selection == "donors":
        hotspot_atoms = ligand_atoms.get("donors", [])
    elif strategy.atom_selection == "coordinating":
        hotspot_atoms = ligand_atoms.get("coordinating", [])
    else:
        hotspot_atoms = ["all"]

    # Fallback if no atoms found for selected strategy
    if not hotspot_atoms:
        hotspot_atoms = ["all"]

    # 3. Burial strategy
    burial_strat = select_burial_strategy(ligand_character, metal_type, design_goal)
    burial = BURIAL_STRATEGIES[burial_strat]

    if burial.target == "metal":
        burial_target = {"X1": "all"}
    elif burial.target == "both":
        burial_target = {"X1": "all", "L1": "all"}
    elif burial.atom_selection == "hydrophobic":
        hydrophobic_atoms = ligand_atoms.get("hydrophobic", [])
        if hydrophobic_atoms:
            burial_target = {"L1": ",".join(hydrophobic_atoms[:10])}
        else:
            burial_target = {"L1": "all"}
    else:
        burial_target = {"L1": "all"}

    # 4. Orientation
    infer_ori = ORIENTATION_RULES.get(design_goal)

    # 5. CFG scale
    cfg = CFG_SCALE_RULES["base"]
    if has_hbond_conditioning:
        cfg = max(cfg, CFG_SCALE_RULES["with_hbond_conditioning"])
    if burial.cfg_boost > 0:
        cfg += burial.cfg_boost
    if ligand_size == LigandSize.LARGE:
        cfg = max(cfg, CFG_SCALE_RULES["complex_ligand"])
    if target_topology == "symmetric":
        cfg = max(cfg, CFG_SCALE_RULES["symmetric"])
    cfg = min(cfg, CFG_SCALE_RULES["max"])

    # 6. MPNN bias
    if metal_type and metal_hsab:
        profile = MPNN_METAL_PROFILES.get(metal_hsab, MPNN_METAL_PROFILES["borderline"])
    else:
        profile = MPNN_METAL_PROFILES["no_metal"]

    # Apply context adjustments
    bias_aa = profile.bias_AA
    if ligand_character == "hydrophobic":
        adjustments = MPNN_CONTEXT_ADJUSTMENTS.get("hydrophobic_burial", {})
        for aa, boost in adjustments.items():
            if boost > 0:
                bias_aa += f",{aa}:{boost}" if bias_aa else f"{aa}:{boost}"

    return DesignDecisions(
        chain_length=chain_length,
        hotspot_strategy=hotspot_strat,
        hotspot_atoms=hotspot_atoms,
        burial_strategy=burial_strat,
        burial_target=burial_target,
        infer_ori_strategy=infer_ori,
        cfg_scale=round(cfg, 1),
        step_scale=STEP_SCALE_DEFAULT,
        mpnn_bias_AA=bias_aa,
        mpnn_omit_AA=profile.omit_AA,
        mpnn_temperature=profile.temperature,
        rationale={
            "chain_length": f"Ligand size {ligand_size.value}, topology {target_topology}",
            "hotspots": f"Strategy: {strategy.description}",
            "burial": f"{burial_strat}: {burial.target} atoms",
            "orientation": f"Goal '{design_goal}' -> {infer_ori or 'auto'}",
            "mpnn": profile.rationale,
        }
    )


# ═══════════════════════════════════════════════════════════════════
# 8. COMPACT REFERENCE FOR CLAUDE PROMPTS
# ═══════════════════════════════════════════════════════════════════

DESIGN_RULES_SUMMARY = """
DESIGN RULES (apply based on ligand/metal):

CHAIN_LENGTH:
- Small (<10 atoms): 80-100, +20 if burial needed
- Medium (10-30): 100-130, +20 if burial
- Large (>30): 120-150, +20 if burial

HOTSPOTS (ligand character):
- Hydrophobic -> C atoms (aromatics, aliphatics)
- Polar -> O,N atoms (donors)
- Metal-coord -> coordinating atoms only

BURIAL:
- Hydrophobic ligand -> bury all ligand
- Polar ligand -> bury metal only
- Catalysis -> bury both

ORIENTATION:
- binding/catalysis -> "com"
- sensing -> "hotspots"
- symmetric -> None

MPNN_BIAS (metal HSAB):
- Hard (Ln,Ca,Mg): E:3,D:3,N:2.4,Q:2.4; omit C
- Soft (Cu,Ag): C:4,H:3,M:2
- Borderline (Zn,Fe): H:3,D:2,E:2,C:1.5

CFG_SCALE:
- Base: 2.0
- +0.5 for H-bond conditioning
- +0.3-0.5 for burial
- Max: 4.0

STEP_SCALE: 1.5 (always)
"""
