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
# 6B. STABILITY PROFILE RULES
# ═══════════════════════════════════════════════════════════════════

@dataclass
class StabilityProfile:
    """
    RFD3 parameter profile for stability optimization.

    Based on RFD3 documentation:
    - step_scale: 0.5-3.0 (higher = more designable structures)
    - gamma_0: 0.3-1.0 (lower = more designable, less diverse)
    - cfg_scale: 0.5-4.0 (higher = stronger constraint adherence)
    - num_timesteps: 50-500 (more = higher quality)
    """
    name: str
    cfg_scale: float          # Classifier-free guidance strength
    step_scale: float         # Step scale for designability
    gamma_0: float            # Gamma parameter (lower = more designable)
    num_timesteps: int        # Denoising steps (more = higher quality)
    plddt_threshold: float    # ESMFold validation threshold
    description: str


STABILITY_PROFILES: Dict[str, StabilityProfile] = {
    "balanced": StabilityProfile(
        name="balanced",
        cfg_scale=2.5,
        step_scale=1.5,
        gamma_0=0.5,
        num_timesteps=200,
        plddt_threshold=0.75,
        description="Default balance of diversity and stability",
    ),
    "stability_focused": StabilityProfile(
        name="stability_focused",
        cfg_scale=3.0,
        step_scale=1.2,
        gamma_0=0.4,
        num_timesteps=250,
        plddt_threshold=0.82,
        description="Prioritize fold stability over diversity",
    ),
    "ultra_stable": StabilityProfile(
        name="ultra_stable",
        cfg_scale=3.5,
        step_scale=1.0,
        gamma_0=0.35,
        num_timesteps=300,
        plddt_threshold=0.88,
        description="Maximum stability, reduced sequence diversity",
    ),
}


def select_stability_profile(
    stability_focus: bool,
    design_goal: str,
    enzyme_class: Optional[str] = None,
) -> StabilityProfile:
    """
    Select stability profile based on design context.

    Args:
        stability_focus: Whether user requested stability optimization
        design_goal: Design goal (binding, catalysis, sensing, structural)
        enzyme_class: Detected enzyme class (if any)

    Returns:
        StabilityProfile with appropriate parameters
    """
    if stability_focus:
        # Explicit stability request
        if design_goal == "structural":
            return STABILITY_PROFILES["ultra_stable"]
        elif enzyme_class:
            # Enzymes need reliable fold but also some flexibility
            return STABILITY_PROFILES["stability_focused"]
        else:
            return STABILITY_PROFILES["stability_focused"]

    # Default behavior based on design goal
    if design_goal == "structural":
        return STABILITY_PROFILES["stability_focused"]

    return STABILITY_PROFILES["balanced"]


# ═══════════════════════════════════════════════════════════════════
# 7. DECISION ENGINE - Main Entry Point
# ═══════════════════════════════════════════════════════════════════

@dataclass
class DesignDecisions:
    """All decisions for a design run."""
    # Core design parameters
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

    # Stability optimization parameters
    stability_profile: str = "balanced"
    gamma_0: Optional[float] = None
    num_timesteps: int = 200
    plddt_threshold: float = 0.75

    # Ligand fixing strategy for RFD3 select_fixed_atoms
    ligand_fixing_strategy: str = "fix_all"           # "fix_coordination" | "diffuse_all" | "fix_all"
    coordination_fixed_atoms: List[str] = field(default_factory=list)  # e.g. ["O2","O5","O7"]

    # Enzyme activity preservation parameters
    select_exposed: Dict[str, str] = field(default_factory=dict)
    catalytic_fixed_residues: List[str] = field(default_factory=list)
    preserve_hbond_network: bool = False
    enzyme_class: Optional[str] = None


def make_design_decisions(
    ligand_size: LigandSize,
    ligand_character: str,        # "hydrophobic" | "polar" | "amphipathic"
    ligand_atoms: Dict[str, List[str]],  # {"hydrophobic": ["C1"...], "donors": ["O1"...]}
    metal_type: Optional[str],
    metal_hsab: Optional[str],    # "hard" | "soft" | "borderline"
    design_goal: str,             # "binding" | "catalysis" | "sensing"
    target_topology: str,         # "monomer" | "dimer" | "symmetric"
    has_hbond_conditioning: bool,
    stability_focus: bool = False,        # User requested stability optimization
    enzyme_class: Optional[str] = None,   # Detected enzyme class
    preserve_function: bool = False,      # Preserve enzymatic function
    ligand_coordination_donors: Optional[List[str]] = None,  # Known metal-coordinating atoms on ligand
    is_scaffolding: bool = False,         # Whether this is a scaffolding (motif) design
    has_ligand: bool = False,             # Whether a ligand is present
    recommended_fixing_strategy: Optional[str] = None,  # Strategy from ligand resolver
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

    # 7. Stability profile selection (NEW)
    stability = select_stability_profile(stability_focus, design_goal, enzyme_class)

    # Override base cfg_scale with stability profile if stability is requested
    if stability_focus:
        cfg = max(cfg, stability.cfg_scale)

    # 8. Enzyme activity preservation (NEW)
    select_exposed: Dict[str, str] = {}
    catalytic_fixed_residues: List[str] = []
    preserve_hbond_network = False

    if enzyme_class and preserve_function:
        # Boost cfg_scale for enzyme preservation
        cfg = max(cfg, 2.8)
        preserve_hbond_network = True
        # Note: select_exposed and catalytic_fixed_residues will be populated
        # by enzyme_chemistry module based on PDB analysis

    # 9. Ligand fixing strategy
    # Rules:
    # - Scaffolding (is_scaffolding=True) → fix_all (preserve extracted motif geometry)
    # - If resolver already set a strategy (e.g. PubChem → diffuse_all), use it
    # - Metal + ligand with known coordination donors → fix_coordination
    # - Metal + ligand without coordination info → fix_all (safe fallback)
    # - Organic ligand only (no metal) → diffuse_all (let RFD3 place it)
    # - No ligand → fix_all (only metal, if present)
    ligand_fixing_strategy = "fix_all"
    coordination_fixed_atoms: List[str] = []

    if is_scaffolding:
        ligand_fixing_strategy = "fix_all"
    elif recommended_fixing_strategy:
        ligand_fixing_strategy = recommended_fixing_strategy
    elif metal_type and has_ligand:
        if ligand_coordination_donors and len(ligand_coordination_donors) > 0:
            ligand_fixing_strategy = "fix_coordination"
            coordination_fixed_atoms = list(ligand_coordination_donors)
        else:
            ligand_fixing_strategy = "fix_all"
    elif has_ligand and not metal_type:
        ligand_fixing_strategy = "diffuse_all"

    return DesignDecisions(
        chain_length=chain_length,
        hotspot_strategy=hotspot_strat,
        hotspot_atoms=hotspot_atoms,
        burial_strategy=burial_strat,
        burial_target=burial_target,
        infer_ori_strategy=infer_ori,
        cfg_scale=round(cfg, 1),
        step_scale=stability.step_scale,  # Use stability profile
        mpnn_bias_AA=bias_aa,
        mpnn_omit_AA=profile.omit_AA,
        mpnn_temperature=profile.temperature,
        ligand_fixing_strategy=ligand_fixing_strategy,
        coordination_fixed_atoms=coordination_fixed_atoms,
        rationale={
            "chain_length": f"Ligand size {ligand_size.value}, topology {target_topology}",
            "hotspots": f"Strategy: {strategy.description}",
            "burial": f"{burial_strat}: {burial.target} atoms",
            "orientation": f"Goal '{design_goal}' -> {infer_ori or 'auto'}",
            "mpnn": profile.rationale,
            "stability": f"Profile: {stability.name} - {stability.description}",
        },
        # Stability parameters
        stability_profile=stability.name,
        gamma_0=stability.gamma_0,
        num_timesteps=stability.num_timesteps,
        plddt_threshold=stability.plddt_threshold,
        # Enzyme parameters
        select_exposed=select_exposed,
        catalytic_fixed_residues=catalytic_fixed_residues,
        preserve_hbond_network=preserve_hbond_network,
        enzyme_class=enzyme_class,
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

STABILITY_PROFILES:
- balanced: cfg=2.5, step=1.5, gamma=0.5, steps=200, plddt=0.75
- stability_focused: cfg=3.0, step=1.2, gamma=0.4, steps=250, plddt=0.82
- ultra_stable: cfg=3.5, step=1.0, gamma=0.35, steps=300, plddt=0.88

ENZYME_PRESERVATION:
- preserve_function=True -> boost cfg to 2.8+, preserve H-bond network
- select_exposed for substrate channel atoms
- catalytic_fixed_residues for active site geometry
"""
