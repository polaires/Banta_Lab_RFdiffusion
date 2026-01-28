"""
RFD3 Configuration Generator

Generates optimal RFdiffusion3 configuration from DesignIntent and ResolvedLigand.

Combines:
- NL-parsed design intent (metal, ligand, goal, topology)
- Ligand chemistry analysis (hotspots, H-bonds, burial)
- HSAB metal chemistry (amino acid bias)
- Learned lessons (cfg_scale, conditioning strategies)
- Decision engine rules (design_rules.py)

Outputs complete RFD3 API parameters ready for inference.
"""

import logging
from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional, Tuple

# Import from our modules
from nl_design_parser import DesignIntent
from ligand_resolver import ResolvedLigand, LigandChemistry, get_recommended_rfd3_params
from metal_chemistry import (
    METAL_DATABASE,
    get_amino_acid_bias,
    get_coordination_number_range,
    is_lanthanide,
    get_hsab_class_simple,
)

# Import decision engine
from design_rules import (
    LigandSize,
    make_design_decisions,
    DesignDecisions,
    STEP_SCALE_DEFAULT,
)

logger = logging.getLogger(__name__)


# =============================================================================
# Configuration Constants
# =============================================================================

# Default RFD3 parameters (from lessons learned)
DEFAULT_RFD3_PARAMS = {
    "task": "rfd3",
    "use_classifier_free_guidance": True,  # CRITICAL: Always use CFG
    "cfg_scale": 2.5,                       # Moderate conditioning
    "step_scale": 1.5,                      # Default step scale
    "num_timesteps": 200,                   # Full denoising
    "num_designs": 4,                       # Default batch size
}

# CFG scale recommendations by design complexity
CFG_SCALE_BY_COMPLEXITY = {
    "simple": 2.0,        # Simple binding, few constraints
    "moderate": 2.5,      # Standard metal-ligand
    "complex": 3.0,       # Multiple constraints, dimer
    "very_complex": 3.5,  # Symmetric, many constraints
}

# Chain length by topology
CHAIN_LENGTH_BY_TOPOLOGY = {
    "monomer": (80, 120),
    "dimer": (100, 150),
    "symmetric": (60, 100),  # Per subunit
    "custom": (80, 140),
}

# Design goal specific settings
GOAL_SETTINGS = {
    "binding": {
        "description": "Optimize for tight ligand binding",
        "prefer_burial": True,
        "cfg_boost": 0.0,
    },
    "catalysis": {
        "description": "Optimize for catalytic activity",
        "prefer_burial": False,  # Need some solvent access
        "cfg_boost": 0.5,  # Stronger constraints
    },
    "sensing": {
        "description": "Optimize for luminescence/sensing",
        "prefer_burial": False,  # Need emission access
        "cfg_boost": 0.0,
    },
    "structural": {
        "description": "Optimize for scaffold stability",
        "prefer_burial": True,
        "cfg_boost": -0.5,  # Rely more on natural folding
    },
}


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class RFD3Config:
    """
    Complete RFD3 configuration for API call.

    Note: Parameter names must match run_rfd3_inference() signature in inference_utils.py.
    """
    # Core parameters
    task: str = "rfd3"
    contig: str = "X1,L1,/0,80-120"
    num_designs: int = 4
    num_timesteps: int = 200

    # CFG parameters
    use_classifier_free_guidance: bool = True
    cfg_scale: float = 2.5
    step_scale: float = 1.5  # From decision engine

    # Stability optimization parameters (Phase 4: AI Infrastructure Enhancement)
    gamma_0: Optional[float] = None  # Noise schedule parameter (lower = more stable)

    # Orientation strategy
    infer_ori_strategy: Optional[str] = None  # "com" | "hotspots" | None

    # Atom selections (Dict format for RFD3 API)
    select_fixed_atoms: Dict[str, str] = field(default_factory=dict)
    select_buried: Dict[str, str] = field(default_factory=dict)
    select_exposed: Dict[str, str] = field(default_factory=dict)  # Enzyme activity: substrate channel

    # Hotspots - can be List[str] or Dict[str, str] depending on mode
    hotspots: List[str] = field(default_factory=list)
    select_hotspots: Dict[str, str] = field(default_factory=dict)  # Atom-level hotspots

    # H-bond conditioning
    select_hbond_acceptor: Dict[str, str] = field(default_factory=dict)
    select_hbond_donor: Dict[str, str] = field(default_factory=dict)

    # Optional parameters
    ori_token: Optional[List[float]] = None

    # Input PDB
    pdb_content: Optional[str] = None

    # Metadata (not passed to API)
    config_source: str = "ai_generated"
    warnings: List[str] = field(default_factory=list)
    notes: Dict[str, Any] = field(default_factory=dict)  # Changed to Dict for decision rationale

    def to_api_params(self) -> Dict[str, Any]:
        """
        Convert to API-ready parameters dict.

        Excludes None values and metadata fields.
        Maps to run_rfd3_inference() parameter names.
        """
        params = {
            "task": self.task,
            "contig": self.contig,
            "num_designs": self.num_designs,
            "num_timesteps": self.num_timesteps,
            "use_classifier_free_guidance": self.use_classifier_free_guidance,
            "cfg_scale": self.cfg_scale,
            "step_scale": self.step_scale,
        }

        # Stability optimization: gamma_0 (noise schedule)
        if self.gamma_0 is not None:
            params["gamma_0"] = self.gamma_0

        # Add orientation strategy if specified
        if self.infer_ori_strategy:
            params["infer_ori_strategy"] = self.infer_ori_strategy

        # Add selections if not empty
        if self.select_fixed_atoms:
            params["select_fixed_atoms"] = self.select_fixed_atoms

        # Hotspots: prefer select_hotspots (atom-level) over hotspots (residue-level)
        if self.select_hotspots:
            params["select_hotspots"] = self.select_hotspots
        elif self.hotspots:
            params["hotspots"] = self.hotspots

        if self.select_buried:
            params["select_buried"] = self.select_buried

        # Enzyme activity preservation: substrate channel exposure
        if self.select_exposed:
            params["select_exposed"] = self.select_exposed

        if self.select_hbond_acceptor:
            params["select_hbond_acceptor"] = self.select_hbond_acceptor
        if self.select_hbond_donor:
            params["select_hbond_donor"] = self.select_hbond_donor

        # Add optional params
        if self.ori_token:
            params["ori_token"] = self.ori_token
        if self.pdb_content:
            params["pdb_content"] = self.pdb_content

        return params

    def to_dict(self) -> Dict[str, Any]:
        """Convert to full dictionary including metadata."""
        return asdict(self)


@dataclass
class LigandMPNNConfig:
    """
    LigandMPNN configuration for sequence design.

    Note: Parameter names must match run_mpnn_inference() signature in inference_utils.py.
    """
    model_type: str = "ligand_mpnn"
    temperature: float = 0.1
    num_sequences: int = 8

    # Residue biasing (uses bias_AA and omit_AA in inference)
    bias_AA: Optional[str] = None
    omit_AA: Optional[str] = None

    # Context options (actual param is use_side_chain_context)
    use_side_chain_context: bool = True

    # Fixed residues (will be populated from structure)
    fixed_positions: Optional[List[str]] = None

    def to_api_params(self) -> Dict[str, Any]:
        """
        Convert to API parameters.

        Maps to run_mpnn_inference() parameter names.
        """
        params = {
            "task": "mpnn",
            "model_type": self.model_type,
            "temperature": self.temperature,
            "num_sequences": self.num_sequences,
            "use_side_chain_context": self.use_side_chain_context,
        }

        if self.bias_AA:
            params["bias_AA"] = self.bias_AA
        if self.omit_AA:
            params["omit_AA"] = self.omit_AA
        if self.fixed_positions:
            params["fixed_positions"] = self.fixed_positions

        return params


# =============================================================================
# Helper Functions for Decision Engine
# =============================================================================

def classify_ligand_size(chemistry: LigandChemistry) -> LigandSize:
    """Classify ligand by heavy atom count."""
    count = chemistry.heavy_atom_count if chemistry.heavy_atom_count > 0 else chemistry.total_atoms
    if count < 10:
        return LigandSize.SMALL
    elif count <= 30:
        return LigandSize.MEDIUM
    else:
        return LigandSize.LARGE


def extract_ligand_atoms(chemistry: LigandChemistry) -> Dict[str, List[str]]:
    """Extract categorized atoms from ligand chemistry for decision engine."""
    return {
        "hydrophobic": chemistry.hydrophobic_atoms or [],
        "donors": chemistry.donor_atoms if chemistry.donor_atoms else [],
        "coordinating": chemistry.coordinating_atoms or chemistry.hbond_acceptors or [],
        "acceptors": chemistry.hbond_acceptors or [],
    }


# =============================================================================
# Config Generator Class
# =============================================================================

class RFD3ConfigGenerator:
    """
    Generates RFD3 and LigandMPNN configurations from design intent.
    """

    def __init__(
        self,
        default_num_designs: int = 4,
        default_num_sequences: int = 8,
    ):
        """
        Initialize generator.

        Args:
            default_num_designs: Default number of backbone designs
            default_num_sequences: Default number of sequences per backbone
        """
        self.default_num_designs = default_num_designs
        self.default_num_sequences = default_num_sequences

    def generate(
        self,
        intent: DesignIntent,
        ligand: ResolvedLigand,
        num_designs: Optional[int] = None,
    ) -> RFD3Config:
        """
        Generate RFD3 configuration using the decision engine.

        Args:
            intent: Parsed design intent
            ligand: Resolved ligand with chemistry
            num_designs: Override for number of designs

        Returns:
            RFD3Config ready for API call
        """
        config = RFD3Config()
        config.num_designs = num_designs or self.default_num_designs

        # === Use Decision Engine ===
        # 1. Classify ligand and extract chemistry data
        ligand_size = LigandSize.MEDIUM  # Default
        ligand_atoms: Dict[str, List[str]] = {}
        ligand_character = "polar"  # Default
        has_hbond_conditioning = False

        if ligand.chemistry:
            ligand_size = classify_ligand_size(ligand.chemistry)
            ligand_atoms = extract_ligand_atoms(ligand.chemistry)
            ligand_character = ligand.chemistry.ligand_character or "polar"
            has_hbond_conditioning = bool(ligand.chemistry.hbond_acceptors)

        # 2. Get HSAB classification for metal
        metal_hsab = None
        if intent.metal_type:
            metal_hsab = get_hsab_class_simple(intent.metal_type)

        # 3. Call decision engine
        decisions = make_design_decisions(
            ligand_size=ligand_size,
            ligand_character=ligand_character,
            ligand_atoms=ligand_atoms,
            metal_type=intent.metal_type,
            metal_hsab=metal_hsab,
            design_goal=intent.design_goal or "binding",
            target_topology=intent.target_topology or "monomer",
            has_hbond_conditioning=has_hbond_conditioning,
        )

        # === Apply Decisions to Config ===
        # 4. Build contig string with decision-engine chain length
        chain_length = decisions.chain_length
        config.contig = self._build_contig_from_decisions(intent, ligand, chain_length)

        # 5. Set CFG and step scale from decisions
        config.cfg_scale = decisions.cfg_scale
        config.step_scale = decisions.step_scale

        # 6. Set orientation strategy from decisions
        config.infer_ori_strategy = decisions.infer_ori_strategy

        # 7. Set hotspots from decisions
        if decisions.hotspot_atoms and decisions.hotspot_atoms != ["all"]:
            # Atom-level hotspots (e.g., {"L1": "C4,C5,C6"})
            config.select_hotspots = {"L1": ",".join(decisions.hotspot_atoms)}
        else:
            # Residue-level hotspots (e.g., ["L1"])
            if ligand.resolved:
                config.hotspots = ["L1"]

        # 8. Set burial from decisions
        config.select_buried = decisions.burial_target

        # 9. Set fixed atoms (metal + ligand)
        config.select_fixed_atoms = self._get_fixed_atoms(intent, ligand)

        # 10. Set PDB content
        if ligand.pdb_content:
            config.pdb_content = ligand.pdb_content

        # 11. Get H-bond conditioning from ligand chemistry
        if ligand.chemistry and ligand.chemistry.hbond_acceptors:
            # H-bond acceptors for polar ligands
            config.select_hbond_acceptor = {"L1": ",".join(ligand.chemistry.hbond_acceptors)}

        # === Store Decision Rationale ===
        config.notes = decisions.rationale.copy()
        config.notes["contig"] = config.contig
        config.notes["ligand_size"] = ligand_size.value
        config.notes["metal_hsab"] = metal_hsab or "none"

        # 12. Validate configuration
        self._validate_config(config, intent, ligand)

        return config

    def _build_contig_from_decisions(
        self,
        intent: DesignIntent,
        ligand: ResolvedLigand,
        chain_length: str,
    ) -> str:
        """Build contig string using decision-engine chain length."""
        parts = []

        # Metal chain (X1)
        if intent.metal_type:
            parts.append("X1")

        # Ligand chain (L1)
        if ligand.resolved:
            parts.append("L1")

        # Protein chain (new protein to design)
        # chain_length from decision engine is like "100-130"
        if intent.target_topology == "dimer":
            # Two chains for dimer
            parts.append(f"/0,{chain_length}")
            parts.append(f"/0,{chain_length}")
        elif intent.target_topology == "symmetric":
            # Single chain, symmetry handled by sym parameter
            parts.append(f"/0,{chain_length}")
        else:
            # Monomer: single chain
            parts.append(f"/0,{chain_length}")

        return ",".join(parts)

    def generate_mpnn_config(
        self,
        intent: DesignIntent,
        ligand: ResolvedLigand,
        num_sequences: Optional[int] = None,
    ) -> LigandMPNNConfig:
        """
        Generate LigandMPNN configuration using decision engine.

        Args:
            intent: Parsed design intent
            ligand: Resolved ligand
            num_sequences: Override for number of sequences

        Returns:
            LigandMPNNConfig for sequence design
        """
        config = LigandMPNNConfig()
        config.num_sequences = num_sequences or self.default_num_sequences

        # === Use Decision Engine for MPNN biasing ===
        # 1. Classify ligand
        ligand_size = LigandSize.MEDIUM
        ligand_atoms: Dict[str, List[str]] = {}
        ligand_character = "polar"
        has_hbond = False

        if ligand.chemistry:
            ligand_size = classify_ligand_size(ligand.chemistry)
            ligand_atoms = extract_ligand_atoms(ligand.chemistry)
            ligand_character = ligand.chemistry.ligand_character or "polar"
            has_hbond = bool(ligand.chemistry.hbond_acceptors)

        # 2. Get HSAB classification for metal
        metal_hsab = None
        if intent.metal_type:
            metal_hsab = get_hsab_class_simple(intent.metal_type)

        # 3. Call decision engine
        decisions = make_design_decisions(
            ligand_size=ligand_size,
            ligand_character=ligand_character,
            ligand_atoms=ligand_atoms,
            metal_type=intent.metal_type,
            metal_hsab=metal_hsab,
            design_goal=intent.design_goal or "binding",
            target_topology=intent.target_topology or "monomer",
            has_hbond_conditioning=has_hbond,
        )

        # 4. Apply decision engine's MPNN parameters
        config.bias_AA = decisions.mpnn_bias_AA if decisions.mpnn_bias_AA else None
        config.omit_AA = decisions.mpnn_omit_AA
        config.temperature = decisions.mpnn_temperature

        return config

    def _build_contig(
        self,
        intent: DesignIntent,
        ligand: ResolvedLigand,
    ) -> str:
        """Build contig string for RFD3."""
        parts = []

        # Metal chain (X1)
        if intent.metal_type:
            parts.append("X1")

        # Ligand chain (L1)
        if ligand.resolved:
            parts.append("L1")

        # Protein chain (new protein to design)
        length_range = f"{intent.chain_length_min}-{intent.chain_length_max}"

        # Adjust for topology
        if intent.target_topology == "dimer":
            # Two chains for dimer
            parts.append(f"/0,{length_range}")
            parts.append(f"/0,{length_range}")
        elif intent.target_topology == "symmetric":
            # Single chain, symmetry handled by sym parameter
            parts.append(f"/0,{length_range}")
        else:
            # Monomer: single chain
            parts.append(f"/0,{length_range}")

        return ",".join(parts)

    def _get_fixed_atoms(
        self,
        intent: DesignIntent,
        ligand: ResolvedLigand,
    ) -> Dict[str, str]:
        """Get fixed atoms selection."""
        fixed = {}

        # Fix metal
        if intent.metal_type:
            fixed["X1"] = "all"

        # Fix ligand
        if ligand.resolved:
            fixed["L1"] = "all"

        return fixed

    def _apply_chemistry_settings(
        self,
        config: RFD3Config,
        intent: DesignIntent,
        ligand: ResolvedLigand,
    ) -> None:
        """Apply settings based on ligand chemistry."""
        if not ligand.chemistry:
            return

        chem = ligand.chemistry

        # For hydrophobic ligands: focus on burial
        if chem.ligand_character == "hydrophobic":
            if not config.select_buried:
                config.select_buried = {"L1": "all"}
            config.notes.append("Hydrophobic ligand: emphasizing burial")

        # For polar ligands: use H-bond conditioning
        elif chem.ligand_character == "polar":
            # H-bond acceptors already set from recommendations
            config.notes.append("Polar ligand: using H-bond conditioning")

        # Set hotspots if not already set
        if not config.hotspots and chem.hotspot_atoms:
            # For ligand binding, hotspot is the residue reference (e.g., "L1")
            # NOT individual atoms (RFD3 expects ChainResID format like "A20")
            config.hotspots.append("L1")

    def _apply_goal_settings(
        self,
        config: RFD3Config,
        intent: DesignIntent,
    ) -> None:
        """Apply settings based on design goal."""
        goal = intent.design_goal
        settings = GOAL_SETTINGS.get(goal, GOAL_SETTINGS["binding"])

        # Apply burial preference
        if settings.get("prefer_burial"):
            if "L1" not in config.select_buried:
                config.select_buried["L1"] = "all"
                config.notes.append(f"Goal '{goal}': enabling full burial")

    def _calculate_cfg_scale(
        self,
        intent: DesignIntent,
        ligand: ResolvedLigand,
    ) -> float:
        """Calculate appropriate CFG scale based on complexity."""
        # Start with ligand recommendation
        if ligand.chemistry:
            base_cfg = ligand.chemistry.cfg_scale_suggestion
        else:
            base_cfg = 2.5

        # Adjust for topology complexity
        if intent.target_topology == "dimer":
            base_cfg += 0.5
        elif intent.target_topology == "symmetric":
            base_cfg += 0.5

        # Adjust for goal
        goal_settings = GOAL_SETTINGS.get(intent.design_goal, {})
        base_cfg += goal_settings.get("cfg_boost", 0.0)

        # Adjust for metal coordination requirements
        if intent.metal_type:
            if is_lanthanide(intent.metal_type):
                # Lanthanides need higher CN, more constraint
                base_cfg += 0.3

        # Clamp to reasonable range
        return max(1.5, min(4.0, base_cfg))

    def _validate_config(
        self,
        config: RFD3Config,
        intent: DesignIntent,
        ligand: ResolvedLigand,
    ) -> None:
        """Validate configuration and add warnings."""
        # Check for common issues

        # CFG without conditioning
        if config.use_classifier_free_guidance:
            has_conditioning = (
                config.hotspots or
                config.select_buried or
                config.select_hbond_acceptor or
                config.select_hbond_donor
            )
            if not has_conditioning:
                config.warnings.append(
                    "CFG enabled but no conditioning specified"
                )

        # Missing PDB content
        if not config.pdb_content:
            config.warnings.append("No PDB content - will need to provide separately")

        # Lanthanide without high CN
        if intent.metal_type and is_lanthanide(intent.metal_type):
            expected_cn = get_coordination_number_range(intent.metal_type, 3)
            if intent.chain_length_max < 80:
                config.warnings.append(
                    f"Lanthanide {intent.metal_type} needs {expected_cn[0]}-{expected_cn[1]} "
                    f"coordination - ensure sufficient protein size"
                )


# =============================================================================
# Convenience Functions
# =============================================================================

def generate_rfd3_config(
    intent: DesignIntent,
    ligand: ResolvedLigand,
    num_designs: int = 4,
) -> RFD3Config:
    """
    Convenience function to generate RFD3 config.

    Args:
        intent: Design intent
        ligand: Resolved ligand
        num_designs: Number of designs to generate

    Returns:
        RFD3Config
    """
    generator = RFD3ConfigGenerator()
    return generator.generate(intent, ligand, num_designs)


def generate_mpnn_config(
    intent: DesignIntent,
    ligand: ResolvedLigand,
    num_sequences: int = 8,
) -> LigandMPNNConfig:
    """
    Convenience function to generate MPNN config.

    Args:
        intent: Design intent
        ligand: Resolved ligand
        num_sequences: Number of sequences per design

    Returns:
        LigandMPNNConfig
    """
    generator = RFD3ConfigGenerator()
    return generator.generate_mpnn_config(intent, ligand, num_sequences)


# =============================================================================
# Test Function
# =============================================================================

def test_config_generator():
    """Test the config generator with decision engine."""
    from nl_design_parser import DesignIntent
    from ligand_resolver import ResolvedLigand, LigandChemistry

    print("=" * 60)
    print("TEST 1: Polar ligand (Citrate-like) with Terbium")
    print("=" * 60)

    # Test citrate-like polar ligand with lanthanide
    intent1 = DesignIntent(
        metal_type="TB",
        ligand_name="citrate",
        design_goal="binding",
        target_topology="monomer",
        chain_length_min=80,
        chain_length_max=120,
        confidence=0.9,
    )

    ligand1 = ResolvedLigand(
        name="citrate",
        metal_type="TB",
        residue_code="CIT",
        resolved=True,
        source="template",
        template_name="citrate_tb",
        pdb_content="HETATM    1  C1  CIT L   1      50.000  50.000  50.000  1.00  0.00           C\nEND",
        chemistry=LigandChemistry(
            donor_atoms=["O2", "O5", "O7"],
            hbond_acceptors=["O1", "O3", "O4", "O6"],
            hydrophobic_atoms=["C1", "C2", "C3", "C4", "C5", "C6"],
            total_atoms=13,
            polar_atoms=7,
            nonpolar_atoms=6,
            ligand_character="polar",
            recommended_burial="partial",
            hotspot_atoms=["O1", "O3", "O4", "O6"],
            cfg_scale_suggestion=2.5,
            heavy_atom_count=13,
            coordinating_atoms=["O2", "O5", "O7"],
        ),
    )

    generator = RFD3ConfigGenerator()
    rfd3_config = generator.generate(intent1, ligand1)
    api_params = rfd3_config.to_api_params()

    print("\nRFD3 Configuration:")
    for key, value in api_params.items():
        if key != "pdb_content":
            print(f"  {key}: {value}")
    print(f"\nDecision Rationale:")
    for key, value in rfd3_config.notes.items():
        print(f"  {key}: {value}")

    # Verify expected values for citrate/TB
    assert "100-130" in api_params["contig"], f"Expected 100-130 chain length, got {api_params['contig']}"
    assert api_params.get("infer_ori_strategy") == "com", f"Expected infer_ori_strategy='com'"
    assert api_params.get("step_scale") == 1.5, f"Expected step_scale=1.5"
    print("\n✓ Citrate/TB test passed!")

    print("\n" + "=" * 60)
    print("TEST 2: Hydrophobic ligand (TriNOx-like) with Dysprosium")
    print("=" * 60)

    intent2 = DesignIntent(
        metal_type="DY",
        ligand_name="trinox",
        design_goal="binding",
        target_topology="monomer",
        chain_length_min=80,
        chain_length_max=120,
        confidence=0.9,
    )

    ligand2 = ResolvedLigand(
        name="trinox",
        metal_type="DY",
        residue_code="TNX",
        resolved=True,
        source="template",
        template_name="trinox_dy",
        pdb_content="HETATM    1  C1  TNX L   1      50.000  50.000  50.000  1.00  0.00           C\nEND",
        chemistry=LigandChemistry(
            donor_atoms=[],
            hbond_acceptors=["O1", "O2"],
            hydrophobic_atoms=["C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12",
                               "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20",
                               "C21", "C22", "C23", "C24", "C25", "C26", "C27", "C28",
                               "C29", "C30", "C31", "C32", "C33", "C34", "C35", "C36"],
            total_atoms=40,
            polar_atoms=5,
            nonpolar_atoms=35,
            ligand_character="hydrophobic",
            recommended_burial="full",
            hotspot_atoms=["C4", "C5", "C6", "C7", "C8"],
            cfg_scale_suggestion=2.8,
            heavy_atom_count=40,
            coordinating_atoms=["O1", "O2", "N1"],
        ),
    )

    rfd3_config2 = generator.generate(intent2, ligand2)
    api_params2 = rfd3_config2.to_api_params()

    print("\nRFD3 Configuration:")
    for key, value in api_params2.items():
        if key != "pdb_content":
            print(f"  {key}: {value}")
    print(f"\nDecision Rationale:")
    for key, value in rfd3_config2.notes.items():
        print(f"  {key}: {value}")

    # Verify expected values for TriNOx/DY
    assert "140-180" in api_params2["contig"], f"Expected 140-180 chain length, got {api_params2['contig']}"
    assert "select_hotspots" in api_params2, "Expected atom-level hotspots for hydrophobic ligand"
    print("\n✓ TriNOx/DY test passed!")

    print("\n" + "=" * 60)
    print("TEST 3: MPNN Biasing with Decision Engine")
    print("=" * 60)

    mpnn_config1 = generator.generate_mpnn_config(intent1, ligand1)
    mpnn_params1 = mpnn_config1.to_api_params()
    print("\nMPNN Config for Citrate/TB (polar ligand, hard acid):")
    for key, value in mpnn_params1.items():
        print(f"  {key}: {value}")

    # Verify hard acid (lanthanide) biasing
    assert mpnn_config1.omit_AA == "C", f"Expected omit_AA='C' for lanthanide, got {mpnn_config1.omit_AA}"
    assert "E:" in (mpnn_config1.bias_AA or ""), f"Expected Glu bias for lanthanide"
    print("\n✓ MPNN Lanthanide biasing test passed!")

    print("\n" + "=" * 60)
    print("All decision engine tests passed!")
    print("=" * 60)


if __name__ == "__main__":
    test_config_generator()
