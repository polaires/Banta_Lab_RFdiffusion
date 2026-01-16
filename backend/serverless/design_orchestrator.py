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
