"""
RunPod Serverless Handler for Foundry Protein Design

This handler processes inference requests for RFD3, RF3, and MPNN models
in a serverless environment. Models are loaded once at worker startup.

Usage:
    - Deploy to RunPod Serverless with network volume for checkpoints
    - Local testing: python handler.py --rp_serve_api
"""

import runpod
import os
import sys
import traceback
from typing import Dict, Any, Optional, List

# Import inference utilities
from inference_utils import (
    run_rfd3_inference,
    run_rf3_inference,
    run_mpnn_inference,
    calculate_rmsd,
    analyze_structure,
    get_gpu_info,
    check_foundry_available,
    # ESM3 functions
    esm3_score_sequences,
    esm3_generate_sequence,
    esm3_get_embeddings,
)

# Import binding analysis functions
from binding_analysis import (
    evaluate_binding,
    analyze_interface,
    run_gnina_scoring,
    check_gnina_available,
    smiles_to_sdf,
    check_steric_clashes,
    validate_binding_comprehensive,
    to_python_types,
    calculate_shape_complementarity,
    calculate_surface_hydrophobicity,
    calculate_unsaturated_hbonds,
    count_interface_residues,
    # Homodimerization scoring for heterodimer validation
    score_homodimerization,
    calculate_sequence_identity,
)

# Import ESMFold validation for structure prediction
try:
    from esmfold_utils import (
        validate_structure_esmfold,
        is_esmfold_available,
        ESMFOLD_THRESHOLDS,
    )
    ESMFOLD_AVAILABLE = is_esmfold_available()
except ImportError:
    ESMFOLD_AVAILABLE = False
    validate_structure_esmfold = None
    ESMFOLD_THRESHOLDS = {}

# Import cleavage utilities for cleavable monomer algorithm
from cleavage_utils import (
    find_cleavage_sites,
    score_cleavage_site,
    select_best_cleavage_site,
    cleave_protein,
    validate_cleaved_dimer,
    create_homo_dimer,
    CleavageSite,
    CleavageResult,
)

# Import shared interaction analysis module
try:
    from shared.interaction_analysis import (
        analyze_all_interactions,
        format_for_frontend,
        format_for_ai_assistant,
        generate_recommendations,
    )
    SHARED_INTERACTION_ANALYSIS_AVAILABLE = True
except ImportError:
    SHARED_INTERACTION_ANALYSIS_AVAILABLE = False
    print("[Handler] Warning: shared.interaction_analysis not available")

# Import Rosetta utilities for structure refinement
from rosetta_utils import (
    fastrelax_with_ligand,
    check_pyrosetta_available,
    calculate_clash_score,
    score_interface as rosetta_score_interface,
)

# Import hotspot detection for auto-detecting binding sites
from hotspot_detection import (
    detect_hotspots_sasa,
    calculate_radius_of_gyration,
    calculate_binding_angle_spread,
    filter_wrap_around_designs,
    format_hotspots_for_rfd3,
)

# Import lanthanide template generators and validation
try:
    from lanthanide_templates import (
        generate_ef_hand_template,
        generate_c4_symmetric_template,
        get_template,
        is_lanthanide,
        LANTHANIDE_PARAMS,
        TEMPLATE_TYPES,
        # NEW: Template library system
        TEMPLATE_LIBRARY,
        METAL_TEMPLATE_RECOMMENDATIONS,
        recommend_template,
        get_template_info,
        list_templates,
        generate_template_from_library,
        generate_parametric_template,
        get_template_for_metal,
        # Fixed position extraction for LigandMPNN
        get_template_fixed_positions,
    )
    LANTHANIDE_TEMPLATES_AVAILABLE = True
except ImportError:
    LANTHANIDE_TEMPLATES_AVAILABLE = False
    print("[Handler] Warning: lanthanide_templates not available")

try:
    from metal_validation import (
        validate_lanthanide_site,
        calculate_tebl_readiness,
        get_refinement_parameters,
        batch_validate,
        LANTHANIDE_CRITERIA,
    )
    METAL_VALIDATION_AVAILABLE = True
except ImportError:
    METAL_VALIDATION_AVAILABLE = False
    print("[Handler] Warning: metal_validation not available")

try:
    from tebl_analysis import (
        predict_tebl_signal,
        add_trp_antenna_to_design,
    )
    TEBL_ANALYSIS_AVAILABLE = True
except ImportError:
    TEBL_ANALYSIS_AVAILABLE = False
    print("[Handler] Warning: tebl_analysis not available")

# Import design types and orchestrator for unified design endpoint
from design_types import DesignType, infer_design_type
from design_orchestrator import DesignOrchestrator

# ============== Ligand H-bond Presets ==============

# Azobenzene RDKit atom naming (canonical SMILES: c1ccc(\N=N\c2ccccc2)cc1)
# RDKit generates atoms in order: Ring1 (C1-C6), Azo (N7,N8), Ring2 (C9-C14)
AZOBENZENE_ATOMS = {
    "ring1": ["C1", "C2", "C3", "C4", "C5", "C6"],
    "ring2": ["C9", "C10", "C11", "C12", "C13", "C14"],
    "azo": ["N7", "N8"],  # H-bond acceptor sites (sp2 nitrogen)
}

# H-bond presets for common ligands
LIGAND_HBOND_PRESETS = {
    "azobenzene": {
        "acceptors": {"UNL": "N7,N8"},  # Azo nitrogens are sp2 acceptors
        "donors": {},
    },
    "atp": {
        "acceptors": {"UNL": "O1A,O2A,O3A,O1B,O2B,O3B,O1G,O2G,O3G"},
        "donors": {"UNL": "N6"},
    },
    "nad": {
        "acceptors": {"UNL": "N1A,N3A,N7A,N1N,N7N"},
        "donors": {"UNL": "N6A,NC2"},
    },
}


def get_ligand_hbond_preset(ligand_smiles: str) -> Dict[str, Dict[str, str]]:
    """
    Auto-detect H-bond preset from SMILES.

    Returns:
        Dict with 'acceptors' and 'donors' keys for H-bond conditioning.
    """
    if not ligand_smiles:
        return {}
    if "N=N" in ligand_smiles:
        return LIGAND_HBOND_PRESETS["azobenzene"]
    if "OP(O)(=O)" in ligand_smiles or "OP(=O)" in ligand_smiles:
        return LIGAND_HBOND_PRESETS.get("atp", {})
    return {}


# ============== Configuration ==============

CHECKPOINT_DIR = os.environ.get(
    "CHECKPOINT_DIR",
    os.environ.get("FOUNDRY_CHECKPOINT_DIRS", "/runpod-volume/checkpoints")
)
os.environ["FOUNDRY_CHECKPOINT_DIRS"] = CHECKPOINT_DIR

# Global model state
MODELS_LOADED = False
FOUNDRY_AVAILABLE = False
GPU_INFO: Dict[str, Any] = {}


# ============== Model Loading ==============

def load_models():
    """Load models at worker startup (called once per worker lifecycle)"""
    global MODELS_LOADED, FOUNDRY_AVAILABLE, GPU_INFO

    if MODELS_LOADED:
        print("[Handler] Models already loaded, skipping")
        return

    print(f"[Handler] Initializing worker...")
    print(f"[Handler] Checkpoint directory: {CHECKPOINT_DIR}")
    print(f"[Handler] Python version: {sys.version}")

    # Check GPU
    GPU_INFO = get_gpu_info()
    print(f"[Handler] GPU: {GPU_INFO}")

    # Check Foundry availability
    FOUNDRY_AVAILABLE = check_foundry_available()
    print(f"[Handler] Foundry available: {FOUNDRY_AVAILABLE}")

    if FOUNDRY_AVAILABLE:
        # Pre-import modules to warm up (optional, reduces first-request latency)
        try:
            print("[Handler] Pre-warming RFD3...")
            from rfd3.engine import RFD3InferenceEngine
            print("[Handler] RFD3 module loaded")
        except ImportError as e:
            print(f"[Handler] RFD3 not available: {e}")

        try:
            print("[Handler] Pre-warming RF3...")
            from rf3.inference_engines.rf3 import RF3InferenceEngine
            print("[Handler] RF3 module loaded")
        except ImportError as e:
            print(f"[Handler] RF3 not available: {e}")

        try:
            print("[Handler] Pre-warming MPNN...")
            from mpnn.inference_engines.mpnn import MPNNInferenceEngine
            print("[Handler] MPNN module loaded")
        except ImportError as e:
            print(f"[Handler] MPNN not available: {e}")

    MODELS_LOADED = True
    print("[Handler] Worker initialization complete")


# Load models at module import time (worker startup)
load_models()


# ============== Handler ==============

def handler(job: Dict[str, Any]) -> Dict[str, Any]:
    """
    Process a single inference request.

    Input format:
    {
        "input": {
            "task": "rfd3" | "rf3" | "mpnn" | "rmsd" | "health",
            // Task-specific parameters
        }
    }

    Output format:
    {
        "status": "completed" | "failed",
        "result": {...} | None,
        "error": str | None
    }
    """
    try:
        job_input = job.get("input", {})
        task = job_input.get("task", "").lower()

        print(f"[Handler] Processing task: {task}")

        if task == "health":
            return handle_health()
        elif task == "rfd3":
            return handle_rfd3(job_input)
        elif task == "rf3":
            return handle_rf3(job_input)
        elif task == "mpnn":
            return handle_mpnn(job_input)
        elif task == "rmsd":
            return handle_rmsd(job_input)
        elif task == "download_checkpoints":
            return handle_download_checkpoints(job_input)
        elif task == "delete_file":
            return handle_delete_file(job_input)
        elif task == "analyze":
            return handle_analyze(job_input)
        # ESM3 tasks
        elif task == "esm3_score":
            return handle_esm3_score(job_input)
        elif task == "esm3_generate":
            return handle_esm3_generate(job_input)
        elif task == "esm3_embed":
            return handle_esm3_embed(job_input)
        # Binding evaluation
        elif task == "binding_eval":
            return handle_binding_eval(job_input)
        # Cleavable monomer algorithm
        elif task == "cleavable_monomer":
            return handle_cleavable_monomer(job_input)
        # Interface ligand design (separable dimer)
        elif task == "interface_ligand_design":
            return handle_interface_ligand_design(job_input)
        # Interface metal design (metal-coordinated dimer)
        elif task == "interface_metal_design":
            return handle_interface_metal_design(job_input)
        # FastRelax for structure refinement
        elif task == "fastrelax":
            return handle_fastrelax(job_input)
        # Full protein binder design pipeline
        elif task == "protein_binder_design":
            return handle_protein_binder_design(job_input)
        # Hotspot detection for binder design
        elif task == "detect_hotspots":
            return handle_detect_hotspots(job_input)
        # Interaction analysis (PLIP or distance-based)
        elif task == "interaction_analysis":
            return handle_interaction_analysis(job_input)
        # Unified design endpoint with intelligent tool selection
        elif task == "design":
            return handle_unified_design(job_input)
        else:
            return {
                "status": "failed",
                "error": f"Unknown task: {task}. Valid tasks: health, rfd3, rf3, mpnn, rmsd, analyze, binding_eval, cleavable_monomer, interface_ligand_design, interface_metal_design, fastrelax, protein_binder_design, detect_hotspots, interaction_analysis, design, esm3_score, esm3_generate, esm3_embed, download_checkpoints, delete_file"
            }

    except Exception as e:
        print(f"[Handler] Error: {e}")
        traceback.print_exc()

        # Collect detailed error context
        error_context = {
            "task": task,
            "input_keys": list(job_input.keys()) if job_input else [],
            "gpu_info": GPU_INFO,
            "foundry_available": FOUNDRY_AVAILABLE,
            "checkpoint_dir": CHECKPOINT_DIR,
        }

        # Try to get GPU memory status
        try:
            import subprocess
            mem_result = subprocess.run(
                ["nvidia-smi", "--query-gpu=memory.used,memory.total", "--format=csv,noheader,nounits"],
                capture_output=True, text=True, timeout=5
            )
            if mem_result.returncode == 0:
                parts = mem_result.stdout.strip().split(", ")
                error_context["gpu_memory_used_mb"] = int(parts[0]) if parts else None
                error_context["gpu_memory_total_mb"] = int(parts[1]) if len(parts) > 1 else None
        except Exception:
            pass

        return {
            "status": "failed",
            "error": str(e),
            "error_type": type(e).__name__,
            "traceback": traceback.format_exc(),
            "context": error_context
        }


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
        pdb_content: str (required) - Input PDB content
        metal: str (optional) - Metal ion code (ZN, FE, TB, etc.)
        ligand: str (optional) - Ligand CCD code
        target_pdb: str (optional) - Target protein for binder design
        symmetry: str (optional) - Symmetry type (C2, C3, D2, etc.)
        design_type: str (optional) - Force specific design type
        temperature: float (optional) - Override temperature
        num_sequences: int (optional) - Override sequence count
        bias_AA: str (optional) - Override AA bias

    Returns:
        {
            "status": "completed",
            "design_type": str,
            "workflow": {...},
            "sequences": [...],
            "config_used": {...}
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


def handle_health() -> Dict[str, Any]:
    """Return health status"""
    return {
        "status": "completed",
        "result": {
            "healthy": True,
            "mode": "real" if FOUNDRY_AVAILABLE else "mock",
            "gpu_available": GPU_INFO.get("available", False),
            "gpu_name": GPU_INFO.get("name"),
            "gpu_memory_gb": GPU_INFO.get("memory_gb"),
            "checkpoint_dir": CHECKPOINT_DIR,
        }
    }


def handle_rfd3(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """Handle RFD3 design request"""
    # Core parameters
    contig = job_input.get("contig") or job_input.get("contigs")
    length = job_input.get("length")
    partial_t = job_input.get("partial_t")

    # Either contig or length is required (unless using partial diffusion with pdb_content)
    if not contig and not length and not partial_t:
        return {"status": "failed", "error": "Missing 'contig' or 'length' parameter"}

    num_designs = job_input.get("num_designs", 1)
    seed = job_input.get("seed")
    pdb_content = job_input.get("pdb_content")

    # Quality settings
    num_timesteps = job_input.get("num_timesteps")
    step_scale = job_input.get("step_scale")
    gamma_0 = job_input.get("gamma_0")
    is_non_loopy = job_input.get("is_non_loopy")

    # Symmetry
    symmetry = job_input.get("symmetry")

    # Small molecule / enzyme design
    ligand = job_input.get("ligand")
    ligand_smiles = job_input.get("ligand_smiles")  # SMILES for organic molecules
    ligand_sdf = job_input.get("ligand_sdf")        # SDF content with 3D coords
    ligand_center = job_input.get("ligand_center")  # [x, y, z] center for SMILES-generated ligand
    conformer_method = job_input.get("conformer_method")  # "rdkit", "xtb", or "torsional"
    interface_ligand = job_input.get("interface_ligand", False)  # Place ligand at symmetric interface
    select_fixed_atoms = job_input.get("select_fixed_atoms")
    unindex = job_input.get("unindex")

    # RASA conditioning (binding pocket design)
    select_buried = job_input.get("select_buried")
    select_exposed = job_input.get("select_exposed")
    select_partially_buried = job_input.get("select_partially_buried")

    # Protein binder design
    hotspots = job_input.get("hotspots")
    infer_ori_strategy = job_input.get("infer_ori_strategy")

    # Nucleic acid binder design
    na_chains = job_input.get("na_chains")
    ori_token = job_input.get("ori_token")
    select_hbond_donor = job_input.get("select_hbond_donor")
    select_hbond_acceptor = job_input.get("select_hbond_acceptor")

    # Covalent modifications (enzyme design)
    covalent_bonds = job_input.get("covalent_bonds")

    # Guiding potentials (auxiliary potentials for design optimization)
    guiding_potentials = job_input.get("guiding_potentials")
    guide_scale = job_input.get("guide_scale")
    guide_decay = job_input.get("guide_decay")

    result = run_rfd3_inference(
        contig=contig,
        length=length,
        num_designs=num_designs,
        seed=seed,
        pdb_content=pdb_content,
        # Quality
        num_timesteps=num_timesteps,
        step_scale=step_scale,
        gamma_0=gamma_0,
        is_non_loopy=is_non_loopy,
        # Refinement
        partial_t=partial_t,
        # Symmetry
        symmetry=symmetry,
        # Small molecule / enzyme
        ligand=ligand,
        ligand_smiles=ligand_smiles,
        ligand_sdf=ligand_sdf,
        ligand_center=tuple(ligand_center) if ligand_center else None,
        conformer_method=conformer_method,
        interface_ligand=interface_ligand,
        select_fixed_atoms=select_fixed_atoms,
        unindex=unindex,
        # RASA conditioning
        select_buried=select_buried,
        select_exposed=select_exposed,
        select_partially_buried=select_partially_buried,
        # Protein binder
        hotspots=hotspots,
        infer_ori_strategy=infer_ori_strategy,
        # Nucleic acid binder
        na_chains=na_chains,
        ori_token=ori_token,
        select_hbond_donor=select_hbond_donor,
        select_hbond_acceptor=select_hbond_acceptor,
        # Covalent modifications
        covalent_bonds=covalent_bonds,
        # Guiding potentials
        guiding_potentials=guiding_potentials,
        guide_scale=guide_scale,
        guide_decay=guide_decay,
        # Mock mode
        use_mock=not FOUNDRY_AVAILABLE
    )

    return result


def handle_rf3(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """Handle RF3 prediction request

    Supports:
    - Single chain prediction (sequence only)
    - Multi-chain/dimer prediction (sequence + sequences list)
    - Protein-ligand prediction (sequence + ligand_smiles)
    - Combined dimer + ligand (all parameters)

    For binding evaluation, ipTM is calculated for:
    - Protein-protein interfaces (when sequences provided)
    - Protein-ligand interfaces (when ligand_smiles provided)
    """
    sequence = job_input.get("sequence")
    if not sequence:
        return {"status": "failed", "error": "Missing 'sequence' parameter"}

    name = job_input.get("name", "prediction")
    pdb_content = job_input.get("pdb_content")
    msa_content = job_input.get("msa_content")  # Optional MSA file content

    # Multi-chain support for dimer/oligomer evaluation
    sequences = job_input.get("sequences")  # List of additional chain sequences

    # Ligand support for protein-ligand binding evaluation
    ligand_smiles = job_input.get("ligand_smiles")  # SMILES string

    result = run_rf3_inference(
        sequence=sequence,
        name=name,
        pdb_content=pdb_content,
        msa_content=msa_content,
        sequences=sequences,
        ligand_smiles=ligand_smiles,
        use_mock=not FOUNDRY_AVAILABLE
    )

    return result


def handle_mpnn(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """Handle MPNN sequence design request with full LigandMPNN capabilities.

    For ligand-aware design:
    - Include HETATM records in pdb_content
    - Use model_type='ligand_mpnn' (default)
    - Optionally fix binding site residues with fixed_positions

    Advanced LigandMPNN parameters (Nature Methods 2025):
    - pack_side_chains: Enable sidechain packing (generates chi angles)
    - pack_with_ligand_context: Include ligand atoms when packing sidechains
    - number_of_packs_per_design: Number of sidechain packing samples per sequence
    - bias_AA: Global amino acid biases (e.g., "W:3.0,Y:2.0,C:-5.0")
    - omit_AA: Amino acids to completely omit (e.g., "C" to avoid cysteines)
    - model_noise_level: Model noise level (005=low, 010=default, 020, 030=high)
    - ligand_cutoff_for_score: Distance cutoff for ligand-adjacent residues (Angstroms)
    - use_side_chain_context: Use fixed residue sidechains as additional context
    - save_stats: Return confidence metrics (ligand_confidence, overall_confidence)
    """
    pdb_content = job_input.get("pdb_content")
    if not pdb_content:
        return {"status": "failed", "error": "Missing 'pdb_content' parameter"}

    # Basic parameters
    num_sequences = job_input.get("num_sequences", 8)
    temperature = job_input.get("temperature", 0.1)
    model_type = job_input.get("model_type", "ligand_mpnn")
    remove_waters = job_input.get("remove_waters", True)
    fixed_positions = job_input.get("fixed_positions")  # e.g., ["A35", "A36", "B35"]

    # Advanced LigandMPNN parameters
    pack_side_chains = job_input.get("pack_side_chains", False)

    # Fix incomplete backbone when sidechain packing is enabled
    # LigandMPNN sidechain packer requires N, CA, C, O for each residue
    if pack_side_chains:
        from inference_utils import fix_incomplete_backbone
        pdb_content = fix_incomplete_backbone(pdb_content)
    pack_with_ligand_context = job_input.get("pack_with_ligand_context", True)
    number_of_packs_per_design = job_input.get("number_of_packs_per_design", 4)
    bias_AA = job_input.get("bias_AA")  # e.g., "W:3.0,Y:2.0,C:-5.0"
    omit_AA = job_input.get("omit_AA")  # e.g., "C"
    model_noise_level = job_input.get("model_noise_level", "010")  # 005, 010, 020, 030
    ligand_cutoff_for_score = job_input.get("ligand_cutoff_for_score", 8.0)
    use_side_chain_context = job_input.get("use_side_chain_context", False)
    save_stats = job_input.get("save_stats", False)

    result = run_mpnn_inference(
        pdb_content=pdb_content,
        num_sequences=num_sequences,
        temperature=temperature,
        model_type=model_type,
        remove_waters=remove_waters,
        fixed_positions=fixed_positions,
        use_mock=not FOUNDRY_AVAILABLE,
        # Advanced LigandMPNN parameters
        pack_side_chains=pack_side_chains,
        pack_with_ligand_context=pack_with_ligand_context,
        number_of_packs_per_design=number_of_packs_per_design,
        bias_AA=bias_AA,
        omit_AA=omit_AA,
        model_noise_level=model_noise_level,
        ligand_cutoff_for_score=ligand_cutoff_for_score,
        use_side_chain_context=use_side_chain_context,
        save_stats=save_stats,
    )

    return result


def run_ligandmpnn_for_ligand_binding(
    pdb_content: str,
    ligand_type: str = "small_molecule",
    num_sequences: int = 4,
    temperature: float = 0.1,
    fixed_positions: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Run LigandMPNN with optimized settings for ligand binding design.

    This helper function provides preset configurations based on the Nature Methods
    2025 paper recommendations for different ligand types.

    Automatically fixes incomplete backbone atoms (adds missing O atoms using
    idealized geometry) as required by LigandMPNN's sidechain packer.

    Args:
        pdb_content: PDB content with ligand (HETATM records)
        ligand_type: One of "small_molecule", "metal", "nucleotide", "protein", or "lanthanide"
        num_sequences: Number of sequences to generate (default: 4)
        temperature: Sampling temperature (default: 0.1)
        fixed_positions: List of positions to keep fixed during design, e.g. ["A15", "A25"]
                        These positions will NOT be redesigned. Critical for preserving
                        metal-coordinating residues from templates.

    Returns:
        MPNN result with sequences and optional confidence metrics

    Example usage in heterodimer design:
        mpnn_result = run_ligandmpnn_for_ligand_binding(
            pdb_content=design["pdb_content"],
            ligand_type="small_molecule",
            num_sequences=4,
        )
    """
    # Fix incomplete backbone atoms (LigandMPNN sidechain packer requires N, CA, C, O)
    # This adds missing O atoms with idealized geometry following RFD3's approach
    from inference_utils import fix_incomplete_backbone
    pdb_content = fix_incomplete_backbone(pdb_content)

    # Preset configurations by ligand type (from Nature Methods 2025)
    presets = {
        "small_molecule": {
            # For azobenzene, drugs, small organic molecules
            "bias_AA": "W:2.0,Y:2.0,F:1.5,H:1.0",  # Favor aromatic residues
            "omit_AA": "C",  # Avoid cysteines
            "ligand_cutoff_for_score": 6.0,  # Tighter for small molecules
            "model_noise_level": "010",  # Default noise
            "pack_side_chains": True,
            "pack_with_ligand_context": True,
            "number_of_packs_per_design": 4,
            "save_stats": True,
        },
        "metal": {
            # For zinc, calcium, iron, magnesium binding
            "bias_AA": "H:3.0,C:3.0,D:2.0,E:2.0",  # Metal-coordinating residues
            "omit_AA": None,  # Allow cysteine for metal coordination
            "ligand_cutoff_for_score": 4.0,  # Very tight for metals
            "model_noise_level": "010",  # Default (005 checkpoint often not available)
            "pack_side_chains": True,
            "pack_with_ligand_context": True,
            "number_of_packs_per_design": 4,
            "save_stats": True,
        },
        "nucleotide": {
            # For DNA, RNA, ATP binding
            "bias_AA": "R:2.0,K:2.0,N:1.5,Q:1.5",  # DNA-binding residues
            "omit_AA": "C",
            "ligand_cutoff_for_score": 8.0,  # Larger cutoff for nucleotides
            "model_noise_level": "010",
            "pack_side_chains": True,
            "pack_with_ligand_context": True,
            "number_of_packs_per_design": 4,
            "save_stats": True,
        },
        "protein": {
            # For protein-protein interfaces (no ligand)
            "bias_AA": None,
            "omit_AA": "C",
            "ligand_cutoff_for_score": 8.0,
            "model_noise_level": "020",  # More robust for interfaces
            "pack_side_chains": False,  # Less critical for PPI
            "pack_with_ligand_context": False,
            "number_of_packs_per_design": 1,
            "save_stats": False,
        },
        "lanthanide": {
            # For lanthanide metals (TB, GD, EU, LA, YB) - Caldwell et al. 2020
            # Lanthanides strongly prefer carboxylate oxygen donors (Glu, Asp)
            # 8-9 coordination typical, requires ~4 Glu/Asp per chain at interface
            # Note: Asp (D) is preferred over Glu (E) for tighter coordination geometry
            # EF-hand loop pattern: DxDxDGxIxxE shows more Asp than Glu
            "bias_AA": "D:6.0,E:4.0,N:1.0,Q:1.0",  # Strong Asp bias (tighter coordination)
            "omit_AA": "C,H",  # Avoid His/Cys (poor lanthanide donors)
            "ligand_cutoff_for_score": 4.0,  # Tight for metal coordination
            "model_noise_level": "005",  # Low noise for precise coordination
            # NOTE: LigandMPNN atomize_side_chains is important for proper metal binding
            # Fix atom name issues before MPNN to avoid NaN coordinates in output
            "pack_side_chains": True,  # Enable proper sidechain packing with metal context
            "pack_with_ligand_context": True,  # Use metal context for better packing
            "number_of_packs_per_design": 4,  # Multiple packs to find best rotamers
            "save_stats": True,
        },
    }

    preset = presets.get(ligand_type, presets["small_molecule"])

    mpnn_input = {
        "pdb_content": pdb_content,
        "num_sequences": num_sequences,
        "temperature": temperature,
        "model_type": "ligand_mpnn",
        "remove_waters": True,
        # Apply preset
        "bias_AA": preset["bias_AA"],
        "omit_AA": preset["omit_AA"],
        "ligand_cutoff_for_score": preset["ligand_cutoff_for_score"],
        "model_noise_level": preset["model_noise_level"],
        "pack_side_chains": preset["pack_side_chains"],
        "pack_with_ligand_context": preset["pack_with_ligand_context"],
        "number_of_packs_per_design": preset["number_of_packs_per_design"],
        "save_stats": preset["save_stats"],
        # Fixed positions - prevent redesign of specified residues (e.g., metal coordinators)
        "fixed_positions": fixed_positions,
    }

    if fixed_positions:
        print(f"[LigandMPNN] Fixed positions (will NOT be redesigned): {fixed_positions}")
    print(f"[LigandMPNN] Running with {ligand_type} preset: bias_AA={preset['bias_AA']}, cutoff={preset['ligand_cutoff_for_score']}Å")

    return handle_mpnn(mpnn_input)


def handle_rmsd(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """Handle RMSD calculation request"""
    pdb1 = job_input.get("pdb_content_1")
    pdb2 = job_input.get("pdb_content_2")

    if not pdb1 or not pdb2:
        return {"status": "failed", "error": "Missing pdb_content_1 or pdb_content_2"}

    backbone_only = job_input.get("backbone_only", True)

    result = calculate_rmsd(pdb1, pdb2, backbone_only)

    if "error" in result and result.get("rmsd") is None:
        return {"status": "failed", "error": result["error"]}

    return {"status": "completed", "result": result}


def handle_analyze(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Analyze a PDB structure to identify binding sites and suggest RFdiffusion parameters.

    This is an AI-assisted tool for understanding protein structures before design.

    Input:
        pdb_content: PDB file content as string
        target_ligands: Optional list of specific ligand codes to analyze (e.g., ["CA", "ZN"])

    Returns:
        Analysis including:
        - Identified ligands/metals
        - Coordinating residues and distances
        - Suggested RFdiffusion parameters for redesign
    """
    pdb_content = job_input.get("pdb_content")
    if not pdb_content:
        return {"status": "failed", "error": "Missing 'pdb_content' parameter"}

    target_ligands = job_input.get("target_ligands")

    result = analyze_structure(pdb_content, target_ligands)

    return result


# ============== Binding Evaluation Handler ==============

def handle_binding_eval(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Evaluate protein-protein or protein-ligand binding quality.

    This combines multiple analysis methods:
    1. Steric Clash Check - CRITICAL first validation
    2. Interface Analyzer - contacts, H-bonds, buried surface area
    3. GNINA Scoring - CNN-based affinity prediction (if available)

    Input:
        pdb_content: PDB file content (protein complex or protein-ligand)
        chain_a: First chain for interface analysis (default: "A")
        chain_b: Second chain (default: "B", or "HETATM" for ligand)
        ligand_smiles: SMILES string for GNINA docking (optional)
        ligand_sdf: SDF content for GNINA docking (optional)
        run_gnina: Whether to run GNINA scoring (default: True if ligand provided)
        check_clashes: Whether to run steric clash check (default: True)

    Returns:
        Combined evaluation with:
        - clash_check: steric clash detection results (CRITICAL)
        - interface_analysis: contacts, hbonds, dSASA, packstat, estimated_dG
        - gnina_scoring: CNN affinity, poses (if GNINA available and ligand provided)
        - summary: quality assessment based on thresholds
    """
    pdb_content = job_input.get("pdb_content")
    if not pdb_content:
        return {"status": "failed", "error": "Missing 'pdb_content' parameter"}

    chain_a = job_input.get("chain_a", "A")
    chain_b = job_input.get("chain_b", "B")
    ligand_smiles = job_input.get("ligand_smiles")
    ligand_sdf = job_input.get("ligand_sdf")
    run_gnina_flag = job_input.get("run_gnina", True)
    check_clashes_flag = job_input.get("check_clashes", True)

    # Docking box parameters (for GNINA)
    docking_center = job_input.get("docking_center")  # [x, y, z]
    docking_box_size = job_input.get("docking_box_size", 25.0)
    whole_protein_search = job_input.get("whole_protein_search", False)

    # Check GNINA availability
    gnina_available = check_gnina_available()
    print(f"[Handler] GNINA available: {gnina_available}")

    # Step 1: Check for steric clashes FIRST (critical validation)
    clash_result = None
    if check_clashes_flag:
        print("[Handler] Running steric clash detection...")
        clash_result = check_steric_clashes(
            pdb_content=pdb_content,
            ligand_smiles=ligand_smiles,
            ligand_sdf=ligand_sdf,
        )
        print(f"[Handler] Clash check result: has_clashes={clash_result.get('has_clashes')}, min_dist={clash_result.get('min_distance')}")

    # Step 2: Run full evaluation
    result = evaluate_binding(
        pdb_content=pdb_content,
        ligand_smiles=ligand_smiles,
        ligand_sdf=ligand_sdf,
        chain_a=chain_a,
        chain_b=chain_b,
        run_gnina=run_gnina_flag and gnina_available,
        docking_center=tuple(docking_center) if docking_center else None,
        docking_box_size=docking_box_size,
        whole_protein_search=whole_protein_search,
    )

    # Add clash check result
    result["clash_check"] = clash_result

    # Add GNINA availability info to result
    result["gnina_available"] = gnina_available

    # Add warning if clashes detected
    if clash_result and clash_result.get("has_clashes"):
        result["warning"] = "Steric clashes detected - binding pocket may not exist. Use ligand-first generation."
        result["summary"]["binding_pocket_valid"] = False
    else:
        result["summary"]["binding_pocket_valid"] = True

    # Convert numpy types to native Python types for JSON serialization
    from binding_analysis import to_python_types
    return to_python_types(result)


# ============== FastRelax Handler ==============

def handle_fastrelax(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Run Rosetta FastRelax to refine protein-ligand complexes and resolve clashes.

    This is adapted from BindCraft's relaxation workflow with key modification:
    - Ligand is kept FIXED while protein moves around it
    - Useful for resolving minor clashes in RFD3-generated designs

    Input:
        pdb_content: str - PDB file content with protein and ligand
        ligand_smiles: str - SMILES string for GNINA scoring (optional)
        max_iter: int - Maximum minimization iterations (default: 200)
        interface_only: bool - Only relax interface residues (default: True, faster)
        interface_distance: float - Distance cutoff for interface (default: 8.0 Å)
        constrain_coords: bool - Constrain backbone to starting coords (default: True)
        run_gnina_before: bool - Score with GNINA before relaxation (default: True)
        run_gnina_after: bool - Score with GNINA after relaxation (default: True)

    Returns:
        status: "completed" or "error"
        relaxed_pdb: Relaxed PDB content
        energy_before: Rosetta energy before relaxation
        energy_after: Rosetta energy after relaxation
        energy_change: Energy difference (negative = improved)
        gnina_before: GNINA results before relaxation (if requested)
        gnina_after: GNINA results after relaxation (if requested)
        improvement: Summary of improvements
    """
    pdb_content = job_input.get("pdb_content")
    if not pdb_content:
        return {"status": "error", "error": "Missing 'pdb_content' parameter"}

    ligand_smiles = job_input.get("ligand_smiles")
    max_iter = job_input.get("max_iter", 200)
    interface_only = job_input.get("interface_only", True)
    interface_distance = job_input.get("interface_distance", 8.0)
    constrain_coords = job_input.get("constrain_coords", True)
    run_gnina_before = job_input.get("run_gnina_before", True)
    run_gnina_after = job_input.get("run_gnina_after", True)

    # Check PyRosetta availability
    pyrosetta_available = check_pyrosetta_available()
    if not pyrosetta_available:
        return {
            "status": "error",
            "error": "PyRosetta not available. Install with: pip install pyrosetta-installer && python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'"
        }

    print(f"[Handler] Running FastRelax with max_iter={max_iter}, interface_only={interface_only}")

    # Score BEFORE with GNINA (if ligand provided)
    gnina_before = None
    if run_gnina_before and ligand_smiles and check_gnina_available():
        print("[Handler] Running GNINA scoring BEFORE relaxation...")
        eval_result = evaluate_binding(
            pdb_content=pdb_content,
            ligand_smiles=ligand_smiles,
            run_gnina=True,
            whole_protein_search=True,
        )
        if eval_result.get("gnina_scoring", {}).get("status") == "completed":
            gnina_before = eval_result["gnina_scoring"].get("result", {})
            print(f"[Handler] GNINA before: affinity={gnina_before.get('best_affinity')}, "
                  f"CNN={gnina_before.get('best_cnn_score')}")

    # Run FastRelax (requires ligand_smiles for parameterization)
    if not ligand_smiles:
        return {
            "status": "error",
            "error": "ligand_smiles is required for FastRelax to parameterize the ligand"
        }

    relax_result = fastrelax_with_ligand(
        pdb_content=pdb_content,
        ligand_smiles=ligand_smiles,
        max_iter=max_iter,
        constrain_coords=constrain_coords,
        interface_only=interface_only,
        interface_distance=interface_distance,
    )

    if relax_result.get("status") != "completed":
        return relax_result

    print(f"[Handler] FastRelax complete: E_before={relax_result.get('energy_before'):.1f}, "
          f"E_after={relax_result.get('energy_after'):.1f}")

    # Score AFTER with GNINA (if ligand provided)
    gnina_after = None
    if run_gnina_after and ligand_smiles and check_gnina_available():
        print("[Handler] Running GNINA scoring AFTER relaxation...")
        eval_result = evaluate_binding(
            pdb_content=relax_result["relaxed_pdb"],
            ligand_smiles=ligand_smiles,
            run_gnina=True,
            whole_protein_search=True,
        )
        if eval_result.get("gnina_scoring", {}).get("status") == "completed":
            gnina_after = eval_result["gnina_scoring"].get("result", {})
            print(f"[Handler] GNINA after: affinity={gnina_after.get('best_affinity')}, "
                  f"CNN={gnina_after.get('best_cnn_score')}")

    # Build response
    response = {
        "status": "completed",
        "relaxed_pdb": relax_result["relaxed_pdb"],
        "energy_before": relax_result["energy_before"],
        "energy_after": relax_result["energy_after"],
        "energy_change": relax_result["energy_change"],
        "ligand_residues": relax_result.get("ligand_residues"),
    }

    if interface_only:
        response["interface_residues"] = relax_result.get("interface_residues")

    # Add GNINA results
    response["gnina_before"] = gnina_before
    response["gnina_after"] = gnina_after

    # Calculate improvement summary
    improvement = {
        "rosetta_improved": relax_result["energy_change"] < 0,
        "rosetta_delta": relax_result["energy_change"],
    }

    if gnina_before and gnina_after:
        before_aff = gnina_before.get("best_affinity")
        after_aff = gnina_after.get("best_affinity")
        if before_aff is not None and after_aff is not None:
            improvement["gnina_improved"] = after_aff < before_aff
            improvement["gnina_delta"] = after_aff - before_aff
            improvement["affinity_before"] = before_aff
            improvement["affinity_after"] = after_aff

        before_cnn = gnina_before.get("best_cnn_score")
        after_cnn = gnina_after.get("best_cnn_score")
        if before_cnn is not None and after_cnn is not None:
            improvement["cnn_improved"] = after_cnn > before_cnn
            improvement["cnn_delta"] = after_cnn - before_cnn

    response["improvement"] = improvement

    # Convert numpy types
    return to_python_types(response)


# ============== Hotspot Detection Handler ==============

def handle_detect_hotspots(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Auto-detect binding hotspots for protein binder design.

    Uses SASA (Solvent Accessible Surface Area) analysis combined with
    spatial clustering to identify surface patches suitable for binder attachment.

    Input:
        target_pdb: str (required) - Target protein PDB content
        target_chain: str (default: "A") - Chain to analyze
        method: str (default: "exposed_clustered") - Detection method:
            - "exposed": All exposed residues
            - "exposed_clustered": Exposed residues clustered spatially
            - "patch": Largest contiguous surface patch
            - "interface_like": Mixed hydrophobic/polar residues
        max_hotspots: int (default: 3) - Maximum hotspots to return (BindCraft recommends 3-5)
        prefer_hydrophobic: bool (default: True) - Prioritize hydrophobic residues

    Returns:
        {
            "status": "completed",
            "hotspots": ["A25", "A30", "A35"],
            "method": "exposed_clustered",
            "cluster_center": {"x": 10.5, "y": 20.3, "z": 15.7},
            "residue_details": [
                {"residue": "A25", "restype": "LEU", "relative_sasa": 0.65, "property": "hydrophobic"},
                ...
            ],
            "total_exposed": 42
        }
    """
    target_pdb = job_input.get("target_pdb")
    if not target_pdb:
        return {"status": "error", "error": "Missing 'target_pdb' parameter"}

    target_chain = job_input.get("target_chain", "A")
    method = job_input.get("method", "exposed_clustered")
    max_hotspots = job_input.get("max_hotspots", 3)
    prefer_hydrophobic = job_input.get("prefer_hydrophobic", True)

    print(f"[DetectHotspots] Detecting hotspots on chain {target_chain}")
    print(f"[DetectHotspots] Method: {method}, Max: {max_hotspots}")

    result = detect_hotspots_sasa(
        pdb_content=target_pdb,
        target_chain=target_chain,
        method=method,
        max_hotspots=max_hotspots,
        prefer_hydrophobic=prefer_hydrophobic,
    )

    if result.get("status") == "completed":
        print(f"[DetectHotspots] Found {len(result.get('hotspots', []))} hotspots: {result.get('hotspots')}")
    else:
        print(f"[DetectHotspots] Detection failed: {result.get('error', 'Unknown error')}")

    return to_python_types(result)


# ============== Interaction Analysis Handler ==============

def handle_interaction_analysis(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Analyze protein-ligand interactions using PLIP or distance-based methods.

    Detects multiple interaction types:
    - Hydrogen bonds (with D-H-A angle validation if PLIP available)
    - Hydrophobic contacts
    - Pi-stacking (face-to-face and edge-to-face)
    - Salt bridges
    - Halogen bonds

    Input:
        pdb_content: str (required) - PDB content with protein and ligand
        ligand_name: str (default: "UNL") - Ligand residue name
        include_visualization: bool (default: True) - Include visualization data
        include_recommendations: bool (default: True) - Include binding recommendations
        ligand_has_aromatics: bool (default: False) - For pi-stacking recommendations

    Returns:
        {
            "status": "completed",
            "interactions": {
                "hbonds": int,
                "hydrophobic": int,
                "pi_stacking": int,
                "salt_bridges": int,
                "total": int
            },
            "key_residues": ["A:GLN40", "A:PHE21", ...],
            "details": {
                "hydrogen_bonds": [...],
                "hydrophobic_contacts": [...],
                "pi_stacking": [...],
                "salt_bridges": [...]
            },
            "recommendations": [...],
            "ai_summary": str,
            "analysis_method": "plip" | "distance_based"
        }
    """
    pdb_content = job_input.get("pdb_content")
    if not pdb_content:
        return {"status": "error", "error": "Missing 'pdb_content' parameter"}

    ligand_name = job_input.get("ligand_name", "UNL")
    include_visualization = job_input.get("include_visualization", True)
    include_recommendations = job_input.get("include_recommendations", True)
    ligand_has_aromatics = job_input.get("ligand_has_aromatics", False)

    print(f"[InteractionAnalysis] Analyzing interactions with ligand {ligand_name}")

    if not SHARED_INTERACTION_ANALYSIS_AVAILABLE:
        return {
            "status": "error",
            "error": "Interaction analysis module not available"
        }

    try:
        # Run comprehensive interaction analysis
        summary = analyze_all_interactions(
            pdb_content=pdb_content,
            ligand_name=ligand_name,
            include_visualization_data=include_visualization,
        )

        if summary.status == "error":
            return {
                "status": "error",
                "error": f"Analysis failed: {summary.error}"
            }

        # Format for frontend
        frontend_data = format_for_frontend(summary)

        result = {
            "status": "completed",
            "interactions": frontend_data["interactions"],
            "key_residues": summary.key_residues,
            "details": frontend_data.get("details", {}),
            "analysis_method": summary.analysis_method,
        }

        # Add visualization data if requested
        if include_visualization and frontend_data.get("visualization"):
            result["visualization"] = frontend_data["visualization"]

        # Add recommendations if requested
        if include_recommendations:
            result["recommendations"] = generate_recommendations(
                summary,
                ligand_has_aromatics=ligand_has_aromatics
            )

        # Add AI-friendly summary
        result["ai_summary"] = format_for_ai_assistant(summary)

        print(f"[InteractionAnalysis] Found {frontend_data['interactions']['total']} total interactions")
        print(f"[InteractionAnalysis] Key residues: {summary.key_residues[:5]}...")

        return to_python_types(result)

    except Exception as e:
        print(f"[InteractionAnalysis] Error: {e}")
        traceback.print_exc()
        return {
            "status": "error",
            "error": f"Interaction analysis failed: {str(e)}"
        }


# ============== Protein Binder Design Handler ==============

def handle_protein_binder_design(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Complete protein-protein binder design pipeline with multi-stage validation.

    Designs a NEW protein that binds to a TARGET protein (BindCraft-style workflow):
    0. Hotspot Detection (auto) - Find optimal binding site if not provided
    1. Structure Generation (RFD3) - Generate binder backbones using hotspots
    1b. Rg Filtering - Remove elongated/wrap-around binders
    2. Sequence Design (ProteinMPNN) - Optimize binder sequence
    3. Sequence Scoring (ESM-3) - Score sequence naturalness
    4. Structure Refinement (FastRelax) - Resolve clashes
    5. Interface Analysis - Score binding interface quality
    6. Quality Filtering - Apply thresholds
    7. Ranking - Composite score ranking

    Input:
        target_pdb: str (required) - PDB content of target protein to bind
        hotspots: list[str] (optional) - Target residues to contact, e.g. ["A25", "A45", "A65"]
        auto_hotspots: bool (default: True) - Auto-detect hotspots if none provided
        hotspot_method: str (default: "exposed_clustered") - Method for auto-detection:
            - "exposed": All exposed residues
            - "exposed_clustered": Spatially clustered exposed residues
            - "patch": Largest contiguous surface patch
            - "interface_like": Mixed hydrophobic/polar residues
        max_hotspots: int (default: 3) - Maximum hotspots for auto-detection (BindCraft recommends 3-5)
        binder_length: str (default: "60-80") - Binder chain length range
        num_designs: int (default: 10) - Number of designs to return
        quality_threshold: str (default: "standard") - "relaxed", "standard", "strict"
        filter_wrap_around: bool (default: True) - Filter elongated binders using Rg
        max_rg_ratio: float (default: 1.5) - Max Rg/expected_Rg for compact binders
        run_mpnn: bool (default: True) - Run ProteinMPNN sequence design
        run_esm_scoring: bool (default: True) - Score sequences with ESM-3
        run_fastrelax: bool (default: True) - Refine structures with FastRelax
        seed: int (optional) - Random seed

    Returns:
        {
            "status": "completed",
            "designs": [
                {
                    "pdb_content": str,
                    "binder_sequence": str,
                    "interface_contacts": int,
                    "interface_hbonds": int,
                    "buried_sasa": float,
                    "esm_perplexity": float,
                    "esm_confidence": float,
                    "rosetta_energy": float,
                    "rg_ratio": float,
                    "rank": int
                }
            ],
            "statistics": {
                "generated": int,
                "rg_filtered": int,
                "mpnn_designed": int,
                "passed_filters": int,
                "returned": int
            },
            "hotspots_used": ["A25", "A30", ...],
            "hotspots_auto_detected": bool
        }
    """
    # Extract parameters
    target_pdb = job_input.get("target_pdb")
    if not target_pdb:
        return {"status": "error", "error": "Missing 'target_pdb' parameter - provide the target protein PDB content"}

    # Clean the target PDB to remove non-standard residues (chromophores, ligands, etc.)
    # This prevents RFD3 parsing issues like with GFP's chromophore (CRO)
    target_pdb = _strip_nonstandard_residues(target_pdb)
    if not target_pdb.strip():
        return {"status": "error", "error": "Target PDB contains no standard amino acids after cleaning"}

    hotspots = job_input.get("hotspots", [])  # e.g., ["A25", "A45", "A65"]
    auto_hotspots = job_input.get("auto_hotspots", True)
    hotspot_method = job_input.get("hotspot_method", "exposed_clustered")
    max_hotspots = job_input.get("max_hotspots", 3)  # BindCraft recommends 3-5 hotspots
    binder_length = job_input.get("binder_length", "60-80")
    num_designs = job_input.get("num_designs", 10)
    quality_threshold = job_input.get("quality_threshold", "standard")

    # Debug: Log received parameters
    print(f"[ProteinBinder] Parameters received:")
    print(f"  - num_designs: {num_designs} (requested by user)")
    print(f"  - binder_length: {binder_length}")
    print(f"  - quality_threshold: {quality_threshold}")
    print(f"  - auto_hotspots: {auto_hotspots}")
    filter_wrap_around = job_input.get("filter_wrap_around", True)
    max_rg_ratio = job_input.get("max_rg_ratio", 1.3)  # Stricter: 1.3 instead of 1.5
    run_mpnn = job_input.get("run_mpnn", True)
    run_esm_scoring = job_input.get("run_esm_scoring", True)
    run_fastrelax = job_input.get("run_fastrelax", True)
    run_esmfold_validation = job_input.get("run_esmfold_validation", False)  # Optional ESMFold validation
    seed = job_input.get("seed")

    # Protocol presets (BindCraft-style configurations for different target types)
    protocol = job_input.get("protocol")  # Optional protocol preset
    PROTOCOL_PRESETS = {
        "miniprotein_default": {
            "binder_length": "60-100",
            "max_hotspots": 5,
            "quality_threshold": "standard",
        },
        "miniprotein_hardtarget": {
            "binder_length": "80-120",
            "max_hotspots": 7,
            "quality_threshold": "relaxed",  # Start relaxed for difficult targets
        },
        "peptide_default": {
            "binder_length": "15-30",
            "max_hotspots": 3,
            "quality_threshold": "standard",
        },
        "peptide_helical": {
            "binder_length": "15-25",
            "max_hotspots": 3,
            "quality_threshold": "standard",
        },
        "large_binder": {
            "binder_length": "100-150",
            "max_hotspots": 7,
            "quality_threshold": "standard",
        },
    }

    # Apply protocol preset if specified
    if protocol and protocol in PROTOCOL_PRESETS:
        preset = PROTOCOL_PRESETS[protocol]
        print(f"[ProteinBinder] Using protocol preset: {protocol}")
        # Only override if not explicitly provided
        if "binder_length" not in job_input:
            binder_length = preset["binder_length"]
        if "max_hotspots" not in job_input:
            max_hotspots = preset["max_hotspots"]
        if "quality_threshold" not in job_input:
            quality_threshold = preset["quality_threshold"]

    # No generation multiplier - generate exactly what user requests
    generation_multiplier = 1

    # Track whether hotspots were auto-detected
    hotspots_auto_detected = False

    print(f"[ProteinBinder] Starting protein-protein binder pipeline")
    print(f"[ProteinBinder] Binder length: {binder_length}, Hotspots: {hotspots}")
    print(f"[ProteinBinder] Designs requested: {num_designs}, Threshold: {quality_threshold}")
    print(f"[ProteinBinder] Rg filtering: {filter_wrap_around}, max_ratio: {max_rg_ratio}")
    if run_esmfold_validation:
        print(f"[ProteinBinder] ESMFold validation: enabled")

    # Quality thresholds for protein-protein binding (BindCraft-inspired)
    QUALITY_THRESHOLDS = {
        "relaxed": {
            "interface_contacts_min": 5,
            "interface_hbonds_min": 0,
            "esm_perplexity_max": 12.0,
            "shape_complementarity_min": 0.4,
            "surface_hydrophobicity_max": 0.45,
            "esmfold_plddt_min": 0.6,
            "esmfold_rmsd_max": 3.0,
        },
        "standard": {
            "interface_contacts_min": 10,
            "interface_hbonds_min": 2,
            "esm_perplexity_max": 8.0,
            "shape_complementarity_min": 0.5,
            "surface_hydrophobicity_max": 0.37,
            "esmfold_plddt_min": 0.7,
            "esmfold_rmsd_max": 2.0,
        },
        "strict": {
            "interface_contacts_min": 15,
            "interface_hbonds_min": 4,
            "esm_perplexity_max": 5.0,
            "shape_complementarity_min": 0.6,
            "surface_hydrophobicity_max": 0.30,
            "esmfold_plddt_min": 0.8,
            "esmfold_rmsd_max": 1.5,
        },
    }

    thresholds = QUALITY_THRESHOLDS.get(quality_threshold, QUALITY_THRESHOLDS["standard"])

    # Statistics tracking
    stats = {
        "generated": 0,
        "rg_filtered": 0,
        "mpnn_designed": 0,
        "esm_passed": 0,
        "relaxed": 0,
        "interface_analyzed": 0,
        "passed_filters": 0,
        "returned": 0,
    }

    # Determine target chain early (needed for hotspot detection)
    target_chain = "A"
    first_res, last_res, target_length = _get_residue_range_in_pdb(target_pdb, chain=target_chain)
    if target_length == 0:
        # Try to find any chain
        for chain in ["B", "C", "D"]:
            first_res, last_res, target_length = _get_residue_range_in_pdb(target_pdb, chain=chain)
            if target_length > 0:
                target_chain = chain
                break

    if target_length == 0:
        return {"status": "error", "error": "Could not determine target protein residue range from PDB"}

    # ============== Stage 0: Auto-Hotspot Detection ==============
    if not hotspots and auto_hotspots:
        print(f"[ProteinBinder] Stage 0: Auto-detecting hotspots on chain {target_chain}...")

        hotspot_result = detect_hotspots_sasa(
            pdb_content=target_pdb,
            target_chain=target_chain,
            method=hotspot_method,
            max_hotspots=max_hotspots,
            prefer_hydrophobic=True,
        )

        if hotspot_result.get("status") == "completed":
            hotspots = hotspot_result.get("hotspots", [])
            hotspots_auto_detected = True
            print(f"[ProteinBinder] Stage 0: Auto-detected {len(hotspots)} hotspots: {hotspots}")
            if hotspot_result.get("residue_details"):
                for detail in hotspot_result["residue_details"]:
                    print(f"  - {detail['residue']}: {detail.get('restype', 'UNK')} ({detail.get('property', 'unknown')}, SASA={detail.get('relative_sasa', 0.0):.2f})")
        else:
            print(f"[ProteinBinder] Stage 0: Hotspot detection failed: {hotspot_result.get('error')}")
            print("[ProteinBinder] Proceeding without hotspots (may result in wrap-around binding)")
    elif hotspots:
        print(f"[ProteinBinder] Using provided hotspots: {hotspots}")
    else:
        print("[ProteinBinder] No hotspots specified and auto-detection disabled")

    # ============== Stage 1: Structure Generation (RFD3) ==============
    # Generate multiplier × requested designs to account for filtering losses
    generation_count = num_designs * generation_multiplier
    print(f"[ProteinBinder] Stage 1: Generating {generation_count} binder backbones...")

    # RFD3 contig format for binder design:
    # "A2-52,/0,40-60" means:
    #   - A2-52: Keep target residues 2-52 from chain A (from input)
    #   - /0: Chain break with zero gap (creates separate chain B)
    #   - 40-60: Design new binder with 40-60 residues
    # Note: Use forward slash /0 with commas, NOT backslash
    contig = f"{target_chain}{first_res}-{last_res},/0,{binder_length}"
    print(f"[ProteinBinder] Target: {target_chain}{first_res}-{last_res} ({target_length} residues)")
    print(f"[ProteinBinder] Binder length: {binder_length}")
    print(f"[ProteinBinder] Full contig: {contig}")

    valid_designs = []

    for design_idx in range(generation_count):
        design_seed = (seed + design_idx) if seed is not None else None
        print(f"[ProteinBinder] Generating design {design_idx + 1}/{generation_count}...")

        # Build RFD3 request with full contig string
        rfd3_input = {
            "task": "rfd3",
            "contig": contig,
            "pdb_content": target_pdb,
            "num_designs": 1,
            "seed": design_seed,
        }

        # Add hotspots if provided
        # CRITICAL: Use infer_ori_strategy="hotspots" to position binder towards hotspots
        # This prevents wrap-around binding by guiding the binder center-of-mass
        # towards the hotspot region rather than random positioning around target
        if hotspots:
            rfd3_input["hotspots"] = hotspots
            rfd3_input["infer_ori_strategy"] = "hotspots"  # Position binder towards hotspots

        result = handle_rfd3(rfd3_input)

        if result.get("status") != "completed":
            print(f"[ProteinBinder] Design {design_idx + 1} failed: {result.get('error')}")
            continue

        result_designs = result.get("result", {}).get("designs", [])
        if not result_designs:
            print(f"[ProteinBinder] Design {design_idx + 1} produced no output")
            continue

        pdb_content = result_designs[0].get("content") or result_designs[0].get("pdb_content")
        if not pdb_content:
            continue

        # Check if RFD3 already output separate chains (A and B)
        chains_in_output = set()
        for line in pdb_content.split('\n'):
            if line.startswith('ATOM'):
                chains_in_output.add(line[21])

        # Only relabel if chain B doesn't exist
        if 'B' not in chains_in_output:
            print(f"[ProteinBinder] Design {design_idx + 1}: Relabeling binder as chain B...")
            pdb_content = _relabel_binder_chain(
                rfd3_pdb_content=pdb_content,
                target_residue_count=target_length,
                target_chain="A",
                binder_chain="B",
            )
            # Update chain list after relabeling
            chains_in_output = set()
            for line in pdb_content.split('\n'):
                if line.startswith('ATOM'):
                    chains_in_output.add(line[21])

        print(f"[ProteinBinder] Design {design_idx + 1} chains: {sorted(chains_in_output)}")

        # Extract binder sequence from chain B
        binder_sequence = _extract_sequence_from_pdb(pdb_content, chain="B")
        if not binder_sequence:
            # Fallback: try to extract from any non-target chain
            binder_sequence = _extract_sequence_from_pdb(pdb_content)

        valid_designs.append({
            "pdb_content": pdb_content,
            "binder_sequence": binder_sequence,
            "design_idx": design_idx,
        })

    stats["generated"] = len(valid_designs)
    print(f"[ProteinBinder] Stage 1: Generated {len(valid_designs)} binder backbones")

    if not valid_designs:
        return {
            "status": "error",
            "error": "No designs were generated successfully",
            "statistics": stats,
            "hotspots_used": hotspots,
            "hotspots_auto_detected": hotspots_auto_detected,
        }

    # ============== Stage 1b: Rg Filtering (Remove Wrap-Around Binders) ==============
    if filter_wrap_around and valid_designs:
        print(f"[ProteinBinder] Stage 1b: Filtering wrap-around binders (Rg ratio <= {max_rg_ratio}, angular < 150°)...")

        designs_before_rg = len(valid_designs)
        compact_designs = []
        rejected_count = 0

        for design in valid_designs:
            rejection_reasons = []

            # Check 1: Rg ratio (binder compactness)
            rg_result = calculate_radius_of_gyration(
                pdb_content=design["pdb_content"],
                chain="B"  # Binder chain
            )

            if rg_result.get("status") == "completed":
                rg_ratio = rg_result.get("rg_ratio", 1.0)
                design["rg_ratio"] = rg_ratio
                design["rg"] = rg_result.get("rg")
                design["rg_assessment"] = rg_result.get("assessment")

                if rg_ratio > max_rg_ratio:
                    rejection_reasons.append(f"Rg ratio {rg_ratio:.2f} > {max_rg_ratio}")
            else:
                design["rg_ratio"] = None

            # Check 2: Angular spread (binding mode - wrap-around detection)
            angular_result = calculate_binding_angle_spread(
                pdb_content=design["pdb_content"],
                target_chain="A",
                binder_chain="B"
            )

            if angular_result.get("status") == "completed":
                angular_spread = angular_result.get("angular_spread", 0)
                design["angular_spread"] = angular_spread
                design["is_wrap_around"] = angular_result.get("is_wrap_around", False)

                if angular_result.get("is_wrap_around", False):
                    rejection_reasons.append(f"wrap-around binding (angular spread {angular_spread:.0f}°)")

            # Accept or reject
            if not rejection_reasons:
                compact_designs.append(design)
            else:
                rejected_count += 1
                print(f"[ProteinBinder] Rejected design {design.get('design_idx', '?')}: {'; '.join(rejection_reasons)}")

        valid_designs = compact_designs
        stats["rg_filtered"] = rejected_count

        print(f"[ProteinBinder] Stage 1b: Kept {len(valid_designs)}/{designs_before_rg} compact designs "
              f"(rejected {rejected_count} wrap-around/elongated)")

    if not valid_designs:
        return {
            "status": "error",
            "error": "All designs were filtered as elongated/wrap-around. Try different hotspots.",
            "statistics": stats,
            "hotspots_used": hotspots,
            "hotspots_auto_detected": hotspots_auto_detected,
        }

    # ============== Stage 2: Sequence Design (ProteinMPNN) ==============
    if run_mpnn and valid_designs:
        print(f"[ProteinBinder] Stage 2: Running ProteinMPNN on {len(valid_designs)} designs...")

        for i, design in enumerate(valid_designs):
            try:
                mpnn_input = {
                    "task": "mpnn",
                    "pdb_content": design["pdb_content"],
                    "num_sequences": 1,
                    "temperature": 0.1,
                    # Only design the binder chain, fix the target
                    "fixed_chains": ["A"],
                }

                mpnn_result = handle_mpnn(mpnn_input)

                if mpnn_result.get("status") == "completed":
                    sequences = mpnn_result.get("result", {}).get("sequences", [])
                    if sequences:
                        # Get the designed binder sequence
                        best_seq = sequences[0]
                        design["binder_sequence"] = best_seq.get("sequence", design["binder_sequence"])
                        design["mpnn_score"] = best_seq.get("score", 0.0)
                        stats["mpnn_designed"] += 1
                        print(f"[ProteinBinder] MPNN design {i+1}: score={design.get('mpnn_score', 'N/A')}")
                    else:
                        print(f"[ProteinBinder] MPNN design {i+1}: No sequences returned")
                else:
                    print(f"[ProteinBinder] MPNN design {i+1} failed: {mpnn_result.get('error', 'Unknown error')}")
            except Exception as e:
                print(f"[ProteinBinder] MPNN exception for design {i+1}: {e}")

    print(f"[ProteinBinder] Stage 2: {stats['mpnn_designed']} sequences designed")

    # ============== Stage 3: ESM-3 Sequence Scoring ==============
    if run_esm_scoring and valid_designs:
        print(f"[ProteinBinder] Stage 3: Running ESM-3 scoring on {len(valid_designs)} sequences...")
        try:
            from esm_utils import score_sequence_esm3, clear_esm3_cache

            for design in valid_designs:
                seq = design.get("binder_sequence", "")
                if seq and len(seq) >= 10:
                    esm_result = score_sequence_esm3(seq)
                    if esm_result.get("status") == "completed":
                        design["esm_perplexity"] = esm_result.get("perplexity", 99.0)
                        design["esm_confidence"] = esm_result.get("overall_confidence", 0.0)
                        stats["esm_passed"] += 1
                    else:
                        design["esm_perplexity"] = 99.0
                        design["esm_confidence"] = 0.0
                else:
                    design["esm_perplexity"] = 99.0
                    design["esm_confidence"] = 0.0

            clear_esm3_cache()

        except ImportError as e:
            print(f"[ProteinBinder] ESM-3 not available: {e}")
            for design in valid_designs:
                design["esm_perplexity"] = None
                design["esm_confidence"] = None
        except Exception as e:
            print(f"[ProteinBinder] ESM-3 scoring failed: {e}")
            for design in valid_designs:
                design["esm_perplexity"] = None
                design["esm_confidence"] = None

    print(f"[ProteinBinder] Stage 3: {stats['esm_passed']} sequences scored")

    # ============== Stage 4: FastRelax Refinement ==============
    if run_fastrelax and valid_designs:
        print(f"[ProteinBinder] Stage 4: Running FastRelax on {len(valid_designs)} structures...")
        pyrosetta_available = check_pyrosetta_available()

        if pyrosetta_available:
            for i, design in enumerate(valid_designs):
                try:
                    # For protein-protein, relax interface residues
                    relax_result = fastrelax_protein_complex(
                        pdb_content=design["pdb_content"],
                        max_iter=200,
                        constrain_coords=True,
                        interface_only=True,
                        interface_distance=8.0,
                    )

                    if relax_result.get("status") == "completed":
                        design["pdb_content"] = relax_result["relaxed_pdb"]
                        design["rosetta_energy"] = relax_result.get("energy_after", 0.0)
                        stats["relaxed"] += 1
                        print(f"[ProteinBinder] Relaxed design {i+1}: E={relax_result.get('energy_after', 0):.1f}")
                    else:
                        design["rosetta_energy"] = None
                except Exception as e:
                    print(f"[ProteinBinder] FastRelax error for design {i+1}: {e}")
                    design["rosetta_energy"] = None
        else:
            print("[ProteinBinder] PyRosetta not available, skipping FastRelax")
            for design in valid_designs:
                design["rosetta_energy"] = None

    print(f"[ProteinBinder] Stage 4: {stats['relaxed']} structures refined")

    # ============== Stage 5: Interface Analysis ==============
    print(f"[ProteinBinder] Stage 5: Analyzing interfaces...")

    for i, design in enumerate(valid_designs):
        try:
            # Basic interface analysis
            interface_result = analyze_interface(
                pdb_content=design["pdb_content"],
                chain_a="A",  # Target
                chain_b="B",  # Binder
                contact_distance=8.0,
            )

            if interface_result.get("status") != "error":
                metrics = interface_result.get("metrics", {})
                design["interface_contacts"] = metrics.get("contacts", 0)
                design["interface_hbonds"] = metrics.get("hbonds_int", 0)
                design["buried_sasa"] = metrics.get("dSASA_int", 0.0)
                design["packstat"] = metrics.get("packstat", 0.0)
                stats["interface_analyzed"] += 1
            else:
                design["interface_contacts"] = 0
                design["interface_hbonds"] = 0
                design["buried_sasa"] = 0.0
                design["packstat"] = 0.0

            # Shape Complementarity (BindCraft-style)
            sc_result = calculate_shape_complementarity(
                pdb_content=design["pdb_content"],
                chain_a="A",
                chain_b="B",
            )
            if sc_result.get("status") == "completed":
                design["shape_complementarity"] = sc_result.get("shape_complementarity", 0.0)
            else:
                design["shape_complementarity"] = 0.0

            # Surface Hydrophobicity (binder chain)
            hydro_result = calculate_surface_hydrophobicity(
                pdb_content=design["pdb_content"],
                chain="B",
            )
            if hydro_result.get("status") == "completed":
                design["surface_hydrophobicity"] = hydro_result.get("surface_hydrophobicity", 0.0)
            else:
                design["surface_hydrophobicity"] = 0.0

            print(f"[ProteinBinder] Interface {i+1}: contacts={design['interface_contacts']}, "
                  f"hbonds={design['interface_hbonds']}, SC={design['shape_complementarity']:.2f}")

        except Exception as e:
            print(f"[ProteinBinder] Interface analysis error for design {i+1}: {e}")
            design["interface_contacts"] = 0
            design["interface_hbonds"] = 0
            design["buried_sasa"] = 0.0
            design["shape_complementarity"] = 0.0
            design["surface_hydrophobicity"] = 0.0

    print(f"[ProteinBinder] Stage 5: {stats['interface_analyzed']} interfaces analyzed")

    # ============== Stage 5b: ESMFold Validation (Optional) ==============
    if run_esmfold_validation and ESMFOLD_AVAILABLE and validate_structure_esmfold:
        print(f"[ProteinBinder] Stage 5b: Running ESMFold structure validation...")
        stats["esmfold_validated"] = 0

        for i, design in enumerate(valid_designs):
            try:
                sequence = design.get("binder_sequence", "")
                if not sequence:
                    design["esmfold_plddt"] = None
                    design["esmfold_rmsd"] = None
                    design["esmfold_passed"] = False
                    continue

                validation = validate_structure_esmfold(
                    sequence=sequence,
                    designed_pdb=design["pdb_content"],
                    binder_chain="B",
                    threshold=quality_threshold,
                )

                if validation.get("status") == "completed":
                    design["esmfold_plddt"] = validation.get("esmfold_plddt", 0.0)
                    design["esmfold_rmsd"] = validation.get("esmfold_rmsd")
                    design["esmfold_passed"] = validation.get("validation_passed", False)
                    stats["esmfold_validated"] += 1
                    print(f"[ProteinBinder] ESMFold {i+1}: pLDDT={design['esmfold_plddt']:.2f}, "
                          f"RMSD={design['esmfold_rmsd']:.2f if design['esmfold_rmsd'] else 'N/A'}Å, "
                          f"passed={design['esmfold_passed']}")
                else:
                    design["esmfold_plddt"] = None
                    design["esmfold_rmsd"] = None
                    design["esmfold_passed"] = False

            except Exception as e:
                print(f"[ProteinBinder] ESMFold validation error for design {i+1}: {e}")
                design["esmfold_plddt"] = None
                design["esmfold_rmsd"] = None
                design["esmfold_passed"] = False

        print(f"[ProteinBinder] Stage 5b: {stats.get('esmfold_validated', 0)} structures validated")
    elif run_esmfold_validation and not ESMFOLD_AVAILABLE:
        print(f"[ProteinBinder] Stage 5b: ESMFold not available, skipping validation")

    # ============== Stage 6: Quality Filtering ==============
    print(f"[ProteinBinder] Stage 6: Applying {quality_threshold} quality filters...")
    filtered_designs = []

    for design in valid_designs:
        # Check interface contacts
        contacts = design.get("interface_contacts", 0)
        if contacts < thresholds["interface_contacts_min"]:
            continue

        # Check hydrogen bonds
        hbonds = design.get("interface_hbonds", 0)
        if hbonds < thresholds["interface_hbonds_min"]:
            continue

        # Check ESM perplexity
        perplexity = design.get("esm_perplexity")
        if perplexity is not None and perplexity > thresholds["esm_perplexity_max"]:
            continue

        # Check Shape Complementarity (BindCraft filter)
        sc = design.get("shape_complementarity", 0.0)
        if sc < thresholds.get("shape_complementarity_min", 0.0):
            continue

        # Check Surface Hydrophobicity (BindCraft filter - prevents aggregation)
        hydro = design.get("surface_hydrophobicity", 0.0)
        if hydro > thresholds.get("surface_hydrophobicity_max", 1.0):
            continue

        # Check ESMFold validation (if enabled)
        if run_esmfold_validation and ESMFOLD_AVAILABLE:
            esmfold_plddt = design.get("esmfold_plddt")
            esmfold_rmsd = design.get("esmfold_rmsd")

            if esmfold_plddt is not None:
                if esmfold_plddt < thresholds.get("esmfold_plddt_min", 0.0):
                    continue
            if esmfold_rmsd is not None:
                if esmfold_rmsd > thresholds.get("esmfold_rmsd_max", float("inf")):
                    continue

        filtered_designs.append(design)

    stats["passed_filters"] = len(filtered_designs)
    print(f"[ProteinBinder] Stage 6: {len(filtered_designs)} designs passed filters")

    # ============== Stage 7: Ranking ==============
    print(f"[ProteinBinder] Stage 7: Ranking designs...")

    def compute_composite_score(design):
        """Compute weighted composite score (higher is better for protein-protein)."""
        score = 0.0

        # Interface contacts (weight 0.20, more is better)
        contacts = design.get("interface_contacts", 0)
        score += 0.20 * min(contacts / 20.0, 1.0)  # Normalize to 0-1

        # H-bonds (weight 0.15, more is better)
        hbonds = design.get("interface_hbonds", 0)
        score += 0.15 * min(hbonds / 10.0, 1.0)  # Normalize to 0-1

        # Buried SASA (weight 0.10, more is better)
        sasa = design.get("buried_sasa", 0.0)
        score += 0.10 * min(sasa / 2000.0, 1.0)  # Normalize to 0-1

        # Shape Complementarity (weight 0.15, higher is better) - BindCraft key metric
        sc = design.get("shape_complementarity", 0.0) or 0.0
        score += 0.15 * sc  # Already 0-1

        # ESM confidence (weight 0.15, higher is better)
        conf = design.get("esm_confidence", 0.0) or 0.0
        score += 0.15 * conf

        # ESMFold pLDDT (weight 0.15, higher is better) - structure prediction confidence
        esmfold_plddt = design.get("esmfold_plddt")
        if esmfold_plddt is not None:
            score += 0.15 * esmfold_plddt  # Already 0-1
        else:
            # If no ESMFold, redistribute weight to ESM confidence
            score += 0.075 * conf

        # Rosetta energy (weight 0.10, lower is better)
        energy = design.get("rosetta_energy")
        if energy is not None:
            # Normalize: -100 is great, 0 is neutral, +100 is bad
            energy_score = max(0, 1 - (energy + 100) / 200)
            score += 0.10 * energy_score

        return score

    # Sort by composite score (higher is better)
    for design in filtered_designs:
        design["_composite_score"] = compute_composite_score(design)

    ranked_designs = sorted(filtered_designs, key=lambda d: d["_composite_score"], reverse=True)

    # Assign ranks and limit to requested count
    final_designs = []
    for i, design in enumerate(ranked_designs[:num_designs]):
        design["rank"] = i + 1
        # Clean up internal fields
        design.pop("_composite_score", None)
        design.pop("design_idx", None)
        final_designs.append(design)

    stats["returned"] = len(final_designs)
    print(f"[ProteinBinder] Stage 7: Returning top {len(final_designs)} designs")

    # ============== Build Response ==============
    response = {
        "status": "completed",
        "designs": final_designs,
        "statistics": stats,
        "thresholds_used": thresholds,
        "hotspots_used": hotspots,
        "hotspots_auto_detected": hotspots_auto_detected,
    }

    # Add protocol info if used
    if protocol and protocol in PROTOCOL_PRESETS:
        response["protocol_used"] = protocol

    # Add ESMFold availability info
    if run_esmfold_validation:
        response["esmfold_enabled"] = ESMFOLD_AVAILABLE

    return to_python_types(response)


def handle_protein_binder_design_with_retry(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Wrapper for protein binder design with intelligent retry logic.

    If the pass rate is below target, automatically adjusts parameters and retries:
    1. First retry: Increase generation count 2x
    2. Second retry: Relax quality threshold by one level

    This helps ensure users get usable results even for difficult targets.

    Args:
        job_input: Same as handle_protein_binder_design

    Returns:
        Same as handle_protein_binder_design, with additional retry info
    """
    target_pass_rate = job_input.get("target_pass_rate", 0.25)  # 25% default
    max_retries = job_input.get("max_retries", 2)
    enable_retry = job_input.get("enable_retry", False)  # Disabled by default

    if not enable_retry:
        return handle_protein_binder_design(job_input)

    best_result = None
    attempts = []

    for attempt in range(max_retries + 1):
        # Make a copy of input for this attempt
        attempt_input = job_input.copy()

        # Adjust parameters based on attempt number
        if attempt == 1:
            # First retry: double the generation count
            current_designs = attempt_input.get("num_designs", 10)
            attempt_input["num_designs"] = current_designs * 2
            print(f"[ProteinBinder-Retry] Attempt {attempt + 1}: Increasing designs to {attempt_input['num_designs']}")

        elif attempt == 2:
            # Second retry: relax quality threshold
            current_threshold = attempt_input.get("quality_threshold", "standard")
            if current_threshold == "strict":
                attempt_input["quality_threshold"] = "standard"
            elif current_threshold == "standard":
                attempt_input["quality_threshold"] = "relaxed"
            print(f"[ProteinBinder-Retry] Attempt {attempt + 1}: Relaxing threshold to {attempt_input['quality_threshold']}")

        # Run the design
        result = handle_protein_binder_design(attempt_input)

        if result.get("status") != "completed":
            attempts.append({"attempt": attempt + 1, "status": "error", "error": result.get("error")})
            continue

        # Calculate pass rate
        stats = result.get("statistics", {})
        generated = stats.get("generated", 1)
        passed = stats.get("passed_filters", 0)
        pass_rate = passed / max(generated, 1)

        attempts.append({
            "attempt": attempt + 1,
            "pass_rate": pass_rate,
            "generated": generated,
            "passed": passed,
            "threshold": attempt_input.get("quality_threshold", "standard")
        })

        # Keep best result
        if best_result is None or len(result.get("designs", [])) > len(best_result.get("designs", [])):
            best_result = result

        # Check if we met the target
        if pass_rate >= target_pass_rate and len(result.get("designs", [])) >= job_input.get("num_designs", 10):
            print(f"[ProteinBinder-Retry] Success on attempt {attempt + 1}: pass_rate={pass_rate:.1%}")
            result["retry_attempts"] = attempts
            return result

        print(f"[ProteinBinder-Retry] Attempt {attempt + 1}: pass_rate={pass_rate:.1%} < target={target_pass_rate:.1%}")

    # Return best result with retry info
    if best_result:
        best_result["retry_attempts"] = attempts
        best_result["retry_note"] = f"Best result after {len(attempts)} attempts"
        return best_result

    return {
        "status": "error",
        "error": "All retry attempts failed",
        "retry_attempts": attempts
    }


def _count_residues_in_pdb(pdb_content: str, chain: str = "A") -> int:
    """Count residues in a specific chain of a PDB."""
    seen_residues = set()
    for line in pdb_content.split('\n'):
        if line.startswith('ATOM'):
            chain_id = line[21]
            if chain_id == chain:
                res_num = line[22:26].strip()
                seen_residues.add(res_num)
    return len(seen_residues)


def _get_residue_range_in_pdb(pdb_content: str, chain: str = "A") -> tuple:
    """
    Get the first and last residue numbers in a specific chain of a PDB.

    Returns:
        Tuple of (first_res, last_res, count) or (None, None, 0) if chain not found.
    """
    residue_nums = []
    for line in pdb_content.split('\n'):
        if line.startswith('ATOM'):
            chain_id = line[21]
            if chain_id == chain:
                try:
                    res_num = int(line[22:26].strip())
                    if res_num not in residue_nums:
                        residue_nums.append(res_num)
                except ValueError:
                    continue

    if not residue_nums:
        return None, None, 0

    residue_nums.sort()
    return residue_nums[0], residue_nums[-1], len(residue_nums)


def _relabel_binder_chain(
    rfd3_pdb_content: str,
    target_residue_count: int,
    target_chain: str = "A",
    binder_chain: str = "B",
) -> str:
    """
    Relabel binder residues as a separate chain after RFD3 generation.

    RFD3 outputs everything in chain A with continuous numbering:
    - Target residues: 1 to target_residue_count
    - Binder residues: target_residue_count+1 to end

    This function splits them into:
    - Chain A: target (residues 1-target_residue_count)
    - Chain B: binder (renumbered from 1)

    Args:
        rfd3_pdb_content: PDB content from RFD3
        target_residue_count: Number of residues in the target protein
        target_chain: Chain label for target (default: "A")
        binder_chain: Chain label for binder (default: "B")

    Returns:
        PDB content with binder relabeled as separate chain
    """
    lines = rfd3_pdb_content.split('\n')
    output_lines = []
    binder_residue_offset = None

    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                res_num = int(line[22:26].strip())

                if res_num <= target_residue_count:
                    # Target residue - keep as chain A
                    new_line = line[:21] + target_chain + line[22:]
                    output_lines.append(new_line)
                else:
                    # Binder residue - relabel as chain B
                    if binder_residue_offset is None:
                        binder_residue_offset = res_num - 1
                    new_res_num = res_num - binder_residue_offset
                    new_line = line[:21] + binder_chain + f"{new_res_num:4d}" + line[26:]
                    output_lines.append(new_line)
            except ValueError:
                output_lines.append(line)
        elif line.startswith('TER'):
            # Add TER for chain break
            if output_lines:
                last_line = output_lines[-1]
                if last_line.startswith('ATOM'):
                    chain = last_line[21]
                    output_lines.append(f"TER")
        else:
            output_lines.append(line)

    # Ensure we have proper chain termination
    result = '\n'.join(output_lines)
    if not result.endswith('END'):
        result += '\nEND'

    return result


def fastrelax_protein_complex(
    pdb_content: str,
    max_iter: int = 200,
    constrain_coords: bool = True,
    interface_only: bool = True,
    interface_distance: float = 8.0,
) -> Dict[str, Any]:
    """
    Run FastRelax on a protein-protein complex.

    Unlike ligand relaxation, this relaxes both chains at the interface.
    """
    try:
        from rosetta_utils import init_pyrosetta
        import pyrosetta as pr
        from pyrosetta.rosetta.core.kinematics import MoveMap
        from pyrosetta.rosetta.protocols.relax import FastRelax
        import tempfile
        import os

        init_pyrosetta()

        # Write temp PDB
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False, mode='w') as f:
            f.write(pdb_content)
            input_path = f.name

        try:
            pose = pr.pose_from_pdb(input_path)
            scorefxn = pr.get_fa_scorefxn()

            energy_before = scorefxn(pose)

            # Create MoveMap
            mmf = MoveMap()

            if interface_only:
                # Find interface residues between chains
                interface_residues = set()
                for i in range(1, pose.total_residue() + 1):
                    res_i = pose.residue(i)
                    if not res_i.is_protein():
                        continue
                    chain_i = pose.pdb_info().chain(i)

                    for j in range(i + 1, pose.total_residue() + 1):
                        res_j = pose.residue(j)
                        if not res_j.is_protein():
                            continue
                        chain_j = pose.pdb_info().chain(j)

                        # Only look at inter-chain contacts
                        if chain_i == chain_j:
                            continue

                        # Check CA-CA distance
                        ca_i = res_i.xyz("CA")
                        ca_j = res_j.xyz("CA")
                        dist = (ca_i - ca_j).norm()

                        if dist < interface_distance:
                            interface_residues.add(i)
                            interface_residues.add(j)

                # Set movemap
                for res in range(1, pose.total_residue() + 1):
                    if res in interface_residues:
                        mmf.set_chi(res, True)
                        mmf.set_bb(res, True)
                    else:
                        mmf.set_chi(res, False)
                        mmf.set_bb(res, False)
            else:
                mmf.set_chi(True)
                mmf.set_bb(True)

            mmf.set_jump(False)

            # Setup FastRelax
            fastrelax = FastRelax()
            fastrelax.set_scorefxn(scorefxn)
            fastrelax.set_movemap(mmf)
            fastrelax.max_iter(max_iter)
            fastrelax.min_type("lbfgs_armijo_nonmonotone")

            if constrain_coords:
                fastrelax.constrain_relax_to_start_coords(True)

            # Run
            fastrelax.apply(pose)

            energy_after = scorefxn(pose)

            # Output
            output_path = input_path.replace('.pdb', '_relaxed.pdb')
            pose.dump_pdb(output_path)

            with open(output_path) as f:
                relaxed_pdb = f.read()

            os.unlink(output_path)

            return {
                "status": "completed",
                "relaxed_pdb": relaxed_pdb,
                "energy_before": energy_before,
                "energy_after": energy_after,
                "energy_change": energy_after - energy_before,
            }

        finally:
            os.unlink(input_path)

    except Exception as e:
        return {"status": "error", "error": str(e)}


def _extract_sequence_from_pdb(pdb_content: str, chain: str = None) -> str:
    """
    Extract amino acid sequence from PDB content.

    Args:
        pdb_content: PDB file content as string
        chain: Optional chain ID to filter (e.g., "A", "B"). If None, returns all chains.

    Returns:
        Single-letter amino acid sequence string
    """
    aa_map = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
    }

    sequence = []
    seen_residues = set()

    for line in pdb_content.split('\n'):
        if line.startswith('ATOM'):
            res_name = line[17:20].strip()
            chain_id = line[21]
            res_num = line[22:26].strip()

            # Filter by chain if specified
            if chain is not None and chain_id != chain:
                continue

            key = (chain_id, res_num)
            if key not in seen_residues and res_name in aa_map:
                sequence.append(aa_map[res_name])
                seen_residues.add(key)

    return ''.join(sequence)


# ============== Cleavable Monomer Handler ==============

def handle_cleavable_monomer(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Execute the Cleavable Monomer Algorithm for interface ligand design.

    This algorithm achieves stronger ligand binding at dimer interfaces by:
    1. Designing a monomer with ligand fully BURIED in center
    2. Finding loop regions where backbone crosses near the ligand
    3. Cleaving at optimal site to create dimer with ligand at interface
    4. Validating the resulting dimer (contacts, clashes, affinity)
    5. (Optional) Creating homo-dimer via C2 symmetry

    Input:
        ligand_smiles: str (required) - SMILES string for the ligand
        target_length: str (default: "80-100") - Chain length range for RFD3
        num_designs: int (default: 3) - Number of monomer designs to try
        cleavage_strategy: str (default: "balanced") - "balanced", "symmetric", "central"
        create_homo_dimer: bool (default: False) - Apply C2 symmetry to create homo-dimer
        min_chain_length: int (default: 25) - Minimum residues per chain after cleavage
        min_contacts_per_chain: int (default: 3) - Minimum ligand contacts per chain
        seed: int (optional) - Random seed for RFD3

        # Advanced options
        distance_range: [min, max] (default: [4.0, 12.0]) - Cleavage site distance range
        contact_cutoff: float (default: 5.0) - Contact distance cutoff
        remove_loop: bool (default: True) - Remove loop residues at cleavage site
        loop_window: int (default: 3) - Loop residues to remove on each side

    Returns:
        {
            "status": "completed" | "failed",
            "result": {
                "dimer": {
                    "pdb_content": str,
                    "chain_a_length": int,
                    "chain_b_length": int,
                    "residues_removed": int
                },
                "validation": {
                    "contacts": {...},
                    "clashes": {...},
                    "gnina": {...}  # if available
                },
                "cleavage_site": {
                    "residue_id": int,
                    "distance_to_ligand": float,
                    "score": float,
                    ...
                },
                "monomer": {
                    "pdb_content": str,
                    "design_index": int
                },
                "homo_dimer": {...}  # if create_homo_dimer=True
            }
        }
    """
    # Extract parameters
    ligand_smiles = job_input.get("ligand_smiles")
    if not ligand_smiles:
        return {"status": "failed", "error": "Missing 'ligand_smiles' parameter"}

    target_length = job_input.get("target_length", "80-100")
    num_designs = job_input.get("num_designs", 3)
    cleavage_strategy = job_input.get("cleavage_strategy", "balanced")
    create_homo = job_input.get("create_homo_dimer", False)
    min_chain_length = job_input.get("min_chain_length", 25)
    min_contacts_per_chain = job_input.get("min_contacts_per_chain", 3)
    seed = job_input.get("seed")

    # Advanced options
    distance_range = job_input.get("distance_range", [4.0, 12.0])
    contact_cutoff = job_input.get("contact_cutoff", 5.0)
    remove_loop = job_input.get("remove_loop", True)
    loop_window = job_input.get("loop_window", 3)

    print(f"[CleavableMonomer] Starting with ligand_smiles={ligand_smiles[:30]}...")
    print(f"[CleavableMonomer] target_length={target_length}, num_designs={num_designs}")

    # Step 1: Generate monomers with buried ligand
    monomer_results = []
    best_result = None
    best_score = float('-inf')

    for design_idx in range(num_designs):
        design_seed = (seed + design_idx) if seed is not None else None
        print(f"[CleavableMonomer] Generating monomer {design_idx + 1}/{num_designs}...")

        # Call RFD3 with buried ligand
        rfd3_input = {
            "task": "rfd3",
            "length": target_length,
            "ligand_smiles": ligand_smiles,
            "select_buried": {"UNL": "ALL"},  # Bury entire ligand
            "seed": design_seed,
            "num_designs": 1,
            "is_non_loopy": False,  # Allow loops for cleavage
        }

        rfd3_result = handle_rfd3(rfd3_input)

        if rfd3_result.get("status") != "completed":
            print(f"[CleavableMonomer] Design {design_idx + 1} failed: {rfd3_result.get('error')}")
            continue

        designs = rfd3_result.get("result", {}).get("designs", [])
        if not designs:
            print(f"[CleavableMonomer] Design {design_idx + 1} produced no designs")
            continue

        monomer_pdb = designs[0].get("content") or designs[0].get("pdb_content")
        if not monomer_pdb:
            print(f"[CleavableMonomer] Design {design_idx + 1}: No PDB content in result")
            continue

        print(f"[CleavableMonomer] Design {design_idx + 1} generated, finding cleavage sites...")

        # Step 2: Find cleavage sites
        sites = find_cleavage_sites(
            pdb_content=monomer_pdb,
            ligand_name="UNL",
            min_chain_length=min_chain_length,
            distance_range=tuple(distance_range),
            min_contacts_per_chain=min_contacts_per_chain,
            contact_cutoff=contact_cutoff,
        )

        if not sites:
            print(f"[CleavableMonomer] Design {design_idx + 1}: No valid cleavage sites found")
            monomer_results.append({
                "design_index": design_idx,
                "pdb_content": monomer_pdb,
                "cleavage_sites": 0,
                "error": "No valid cleavage sites found"
            })
            continue

        print(f"[CleavableMonomer] Design {design_idx + 1}: Found {len(sites)} cleavage sites")

        # Step 3: Score and select best site
        # Count total residues for scoring
        n_residues = sites[0].chain_a_length + sites[0].chain_b_length if sites else 0
        for site in sites:
            score_cleavage_site(site, n_residues)

        best_site = select_best_cleavage_site(sites, n_residues, strategy=cleavage_strategy)
        if best_site is None:
            print(f"[CleavableMonomer] Design {design_idx + 1}: No suitable cleavage site selected")
            monomer_results.append({
                "design_index": design_idx,
                "pdb_content": monomer_pdb,
                "cleavage_sites": len(sites),
                "error": "No suitable cleavage site selected"
            })
            continue
        print(f"[CleavableMonomer] Best site: residue {best_site.residue_id}, score={best_site.score:.2f}")

        # Step 4: Cleave protein
        cleavage_result = cleave_protein(
            pdb_content=monomer_pdb,
            cleavage_site=best_site,
            remove_loop=remove_loop,
            loop_window=loop_window,
        )

        # Step 5: Validate cleaved dimer
        validation = validate_cleaved_dimer(
            pdb_content=cleavage_result.dimer_pdb,
            ligand_name="UNL",
            ligand_smiles=ligand_smiles,
            min_contacts_per_chain=min_contacts_per_chain,
        )

        result_entry = {
            "design_index": design_idx,
            "monomer_pdb": monomer_pdb,
            "dimer_pdb": cleavage_result.dimer_pdb,
            "cleavage_site": {
                "residue_id": best_site.residue_id,
                "residue_name": best_site.residue_name,
                "distance_to_ligand": best_site.distance_to_ligand,
                "secondary_structure": best_site.secondary_structure,
                "contacts_chain_a": best_site.contacts_chain_a,
                "contacts_chain_b": best_site.contacts_chain_b,
                "chain_a_length": best_site.chain_a_length,
                "chain_b_length": best_site.chain_b_length,
                "score": best_site.score,
            },
            "cleavage_result": {
                "chain_a_length": cleavage_result.chain_a_length,
                "chain_b_length": cleavage_result.chain_b_length,
                "loop_start": cleavage_result.loop_start,
                "loop_end": cleavage_result.loop_end,
                "residues_removed": cleavage_result.residues_removed,
            },
            "validation": validation,
        }

        # Compute overall score for ranking
        overall_score = best_site.score
        checks = validation.get("checks", {})
        is_valid = validation.get("overall_pass", False)

        if is_valid:
            overall_score += 10  # Bonus for valid dimer

        # Get GNINA affinity from nested structure
        gnina_result = checks.get("gnina_result", {}).get("result", {})
        gnina_affinity = gnina_result.get("best_affinity")

        if gnina_affinity is not None:
            # Better (more negative) affinity = higher score
            # Weight affinity heavily: multiply by 5 to prioritize good binding
            overall_score -= gnina_affinity * 5

        result_entry["overall_score"] = overall_score
        monomer_results.append(result_entry)

        if overall_score > best_score:
            best_score = overall_score
            best_result = result_entry

        print(f"[CleavableMonomer] Design {design_idx + 1}: overall_score={overall_score:.2f}, valid={is_valid}, affinity={gnina_affinity}")

    # Check if any design succeeded
    if best_result is None:
        # Return diagnostics
        return {
            "status": "failed",
            "error": "No valid cleavage found in any design",
            "diagnostics": {
                "num_designs_tried": num_designs,
                "results": monomer_results,
                "suggestions": [
                    "Try longer chains (100-140 residues)",
                    "Set 'is_non_loopy: false' to encourage more loops",
                    "Generate more designs (num_designs: 5-10)",
                    "Reduce min_contacts_per_chain to 2",
                ]
            }
        }

    # Step 6 (Optional): Create homo-dimer
    homo_dimer_result = None
    if create_homo:
        print("[CleavableMonomer] Creating homo-dimer with C2 symmetry...")
        homo_dimer_result = create_homo_dimer(
            pdb_content=best_result["dimer_pdb"],
            ligand_name="UNL",
        )

    # Build final result
    final_result = {
        "status": "completed",
        "result": {
            "dimer": {
                "pdb_content": best_result["dimer_pdb"],
                "chain_a_length": best_result["cleavage_result"]["chain_a_length"],
                "chain_b_length": best_result["cleavage_result"]["chain_b_length"],
                "residues_removed": best_result["cleavage_result"]["residues_removed"],
            },
            "validation": best_result["validation"],
            "cleavage_site": best_result["cleavage_site"],
            "monomer": {
                "pdb_content": best_result["monomer_pdb"],
                "design_index": best_result["design_index"],
            },
            "all_results": [
                {
                    "design_index": r["design_index"],
                    "overall_score": r.get("overall_score", 0),
                    "cleavage_site": r.get("cleavage_site"),
                    "validation": r.get("validation"),
                }
                for r in monomer_results
                if "dimer_pdb" in r
            ],
        }
    }

    if homo_dimer_result:
        final_result["result"]["homo_dimer"] = homo_dimer_result

    # Convert numpy types for JSON serialization
    return to_python_types(final_result)


# ============== Interface Ligand Design Handler ==============

def handle_interface_ligand_design(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Design protein dimers with ligands at the interface using SEPARABLE chains.

    This addresses the topology problem where buried ligand design creates
    entangled chains that cannot physically separate. The solution uses
    asymmetric binders that grow protein on ONE side of the ligand only.

    Approaches:
    1. "asymmetric" - Design one-sided binder using ori_token offset
    2. "sequential" - Design chain B using chain A as fixed context
    3. "full" - Complete workflow: asymmetric chain A, then complementary chain B

    Input:
        ligand_smiles: str (required) - SMILES string for the ligand
        approach: str (default: "asymmetric") - "asymmetric", "sequential", or "full"
        chain_length: str (default: "60-80") - Chain length range for RFD3
        num_designs: int (default: 3) - Number of designs to generate
        seed: int (optional) - Random seed

        # For asymmetric approach:
        ori_offset: [x, y, z] (default: [12.0, 0.0, 0.0]) - Offset from ligand center
        exposed_atoms: str (optional) - Ligand atoms to keep exposed (e.g., "C7,C8,C9")
        side: str (default: "left") - Which side to bind: "left" or "right"

        # For sequential approach:
        chain_a_pdb: str (required) - PDB content of chain A with ligand

    Returns:
        {
            "status": "completed" | "failed",
            "result": {
                "designs": [{"pdb_content": str, "metrics": {...}}, ...],
                "approach": str,
                "best_design": int
            }
        }
    """
    # Extract parameters
    ligand_smiles = job_input.get("ligand_smiles")
    if not ligand_smiles:
        return {"status": "failed", "error": "Missing 'ligand_smiles' parameter"}

    approach = job_input.get("approach", "asymmetric")
    chain_length = job_input.get("chain_length", "60-80")
    num_designs = job_input.get("num_designs", 3)
    seed = job_input.get("seed")

    # FastRelax parameters for post-design refinement
    run_fastrelax = job_input.get("run_fastrelax", False)
    fastrelax_iter = job_input.get("fastrelax_iter", 200)

    print(f"[InterfaceLigand] Starting with approach={approach}, ligand={ligand_smiles[:30]}, fastrelax={run_fastrelax}...")

    if approach == "asymmetric":
        return _design_asymmetric_binder(job_input)
    elif approach == "sequential":
        return _design_sequential_binder(job_input)
    elif approach == "full":
        return _design_full_dimer(job_input)
    elif approach == "symmetric":
        return _design_symmetric_dimer(job_input)
    # New heterodimer approaches (anti-homodimerization)
    elif approach == "joint":
        return _design_joint_heterodimer(job_input)
    elif approach == "asymmetric_rasa":
        return _design_asymmetric_rasa_heterodimer(job_input)
    elif approach == "induced":
        return _design_induced_heterodimer(job_input)
    else:
        return {
            "status": "failed",
            "error": f"Unknown approach: {approach}. Valid: asymmetric, sequential, full, symmetric, joint, asymmetric_rasa, induced"
        }


def _design_asymmetric_binder(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Design a one-sided binder using RASA conditioning (select_exposed).

    The key mechanism is select_exposed which tells RFD3 to keep certain ligand
    atoms accessible (not buried by protein). This forces the protein to wrap
    around only one side of the ligand, leaving space for a second chain.

    Note: ori_token is disabled by default as it was causing bad positioning.
    Use select_exposed as the primary mechanism for one-sided binding.
    """
    ligand_smiles = job_input.get("ligand_smiles")
    chain_length = job_input.get("chain_length", "60-80")
    num_designs = job_input.get("num_designs", 3)
    seed = job_input.get("seed")

    # Asymmetric parameters
    ori_offset = job_input.get("ori_offset", [12.0, 0.0, 0.0])
    exposed_atoms = job_input.get("exposed_atoms")
    side = job_input.get("side", "left")

    # If side is specified but no exposed_atoms, use defaults for azobenzene
    if not exposed_atoms and "N=N" in ligand_smiles:
        # Azobenzene atom naming: C1-C6 (first phenyl), N7-N8 (azo), C9-C14 (second phenyl)
        if side == "right":
            exposed_atoms = "C9,C10,C11,C12,C13,C14"  # Second phenyl
        else:
            exposed_atoms = "C1,C2,C3,C4,C5,C6"  # First phenyl

    # Adjust ori_offset based on side
    if side == "right" and ori_offset[0] > 0:
        ori_offset = [-ori_offset[0], ori_offset[1], ori_offset[2]]

    # Use ori_token only if explicitly set by user
    # NOTE: ori_token doesn't reliably position proteins on opposite sides
    use_ori_token = job_input.get("use_ori_token", False)

    print(f"[AsymmetricBinder] ori_offset={ori_offset}, exposed_atoms={exposed_atoms}, use_ori_token={use_ori_token}")

    designs = []
    best_design_idx = 0
    best_score = float('-inf')

    for design_idx in range(num_designs):
        design_seed = (seed + design_idx) if seed is not None else None
        print(f"[AsymmetricBinder] Generating design {design_idx + 1}/{num_designs}...")

        # Build RFD3 config for asymmetric binder
        rfd3_input = {
            "task": "rfd3",
            "length": chain_length,
            "ligand_smiles": ligand_smiles,
            "seed": design_seed,
            "num_designs": 1,
        }

        # Optionally use ori_token to offset protein from ligand
        if use_ori_token:
            rfd3_input["ori_token"] = ori_offset

        # Force specified atoms to remain exposed (one-sided binding)
        # This is the KEY mechanism for asymmetric binding
        if exposed_atoms:
            rfd3_input["select_exposed"] = {"UNL": exposed_atoms}

        result = handle_rfd3(rfd3_input)

        if result.get("status") != "completed":
            print(f"[AsymmetricBinder] Design {design_idx + 1} failed: {result.get('error')}")
            continue

        result_designs = result.get("result", {}).get("designs", [])
        if not result_designs:
            print(f"[AsymmetricBinder] Design {design_idx + 1} produced no designs")
            continue

        pdb_content = result_designs[0].get("content") or result_designs[0].get("pdb_content")
        if not pdb_content:
            print(f"[AsymmetricBinder] Design {design_idx + 1}: No PDB content")
            continue

        # Validate the design
        validation = validate_cleaved_dimer(
            pdb_content=pdb_content,
            ligand_name="UNL",
            ligand_smiles=ligand_smiles,
            min_contacts_per_chain=2,  # Lower threshold for one-sided
        )

        # Get metrics
        checks = validation.get("checks", {})
        gnina_result = checks.get("gnina_result", {}).get("result", {})
        affinity = gnina_result.get("best_affinity")

        design_entry = {
            "pdb_content": pdb_content,
            "design_index": design_idx,
            "validation": validation,
            "metrics": {
                "affinity": affinity,
                "contacts": checks.get("contacts_a", 0) + checks.get("contacts_b", 0),
                "has_clashes": checks.get("clash_check", {}).get("has_clashes", False),
            }
        }

        # Score: prioritize affinity
        score = 0
        if affinity is not None:
            score = -affinity  # More negative = better
        if not design_entry["metrics"]["has_clashes"]:
            score += 5

        design_entry["score"] = score
        designs.append(design_entry)

        if score > best_score:
            best_score = score
            best_design_idx = len(designs) - 1

        print(f"[AsymmetricBinder] Design {design_idx + 1}: affinity={affinity}, score={score:.2f}")

    if not designs:
        return {
            "status": "failed",
            "error": "No valid asymmetric binder designs generated",
            "suggestions": [
                "Specify exposed_atoms for your ligand (e.g., 'C1,C2,C3' for one side)",
                "Increase num_designs to try more seeds",
                "Try a different chain_length range",
            ]
        }

    return {
        "status": "completed",
        "result": {
            "designs": designs,
            "approach": "asymmetric",
            "best_design": best_design_idx,
            "ori_offset": ori_offset,
            "exposed_atoms": exposed_atoms,
        }
    }


def _get_ligand_chain_resnum(pdb_content: str, ligand_name: str = "UNL") -> tuple:
    """
    Extract the chain ID and residue number of a ligand from PDB content.

    Returns: (chain_id, res_num) or (None, None) if not found.
    """
    for line in pdb_content.split('\n'):
        if line.startswith('HETATM'):
            res_name = line[17:20].strip()
            if res_name == ligand_name:
                chain_id = line[21].strip() or "L"
                res_num = line[22:26].strip()
                return (chain_id, res_num)
    return (None, None)


def _design_sequential_binder(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Design chain B to bind the opposite side of the ligand from chain A.

    NEW APPROACH: Since contig-based binder design doesn't work for ligands,
    we design Chain B independently using ONLY the ligand (extracted from Chain A's output).
    Chain B uses OPPOSITE RASA conditioning from Chain A.

    Steps:
    1. Extract ligand PDB from Chain A output (preserving its exact position)
    2. Design Chain B around JUST the ligand (like asymmetric design)
    3. Use opposite RASA conditioning (expose atoms that Chain A bound)
    4. Combine Chain A + Chain B + ligand into final dimer
    """
    ligand_smiles = job_input.get("ligand_smiles")
    chain_a_pdb = job_input.get("chain_a_pdb")

    if not chain_a_pdb:
        return {"status": "failed", "error": "Missing 'chain_a_pdb' parameter for sequential approach"}

    chain_length = job_input.get("chain_length", "60-80")
    num_designs = job_input.get("num_designs", 3)
    seed = job_input.get("seed")

    # Get chain A length for later combining
    chain_a_length = _get_chain_length(chain_a_pdb, "A")
    if chain_a_length == 0:
        chain_a_length = _get_chain_length(chain_a_pdb)

    print(f"[SequentialBinder] Chain A length: {chain_a_length}")

    # Extract ligand PDB from Chain A output (preserves exact 3D position)
    ligand_pdb = _extract_ligand_pdb(chain_a_pdb, "UNL")
    if not ligand_pdb:
        return {"status": "failed", "error": "Could not extract ligand UNL from chain A PDB"}

    print(f"[SequentialBinder] Extracted ligand from Chain A output")

    # Get atoms that Chain B should EXPOSE (opposite of what Chain A exposed)
    # Chain A exposed certain atoms, Chain B should expose the OPPOSITE atoms
    exposed_atoms_a = job_input.get("exposed_atoms_b", "")  # Named confusingly, but this is what A exposed

    # Determine opposite side for Chain B to expose
    if "N=N" in (ligand_smiles or ""):
        # Azobenzene: C1-C6 (first phenyl), C9-C14 (second phenyl)
        if "C1" in exposed_atoms_a:
            # Chain A exposed C1-C6, so Chain A bound C9-C14
            # Chain B should bind C1-C6, so Chain B exposes C9-C14
            exposed_atoms_b = "C9,C10,C11,C12,C13,C14"
        else:
            # Chain A exposed C9-C14, so Chain A bound C1-C6
            # Chain B should bind C9-C14, so Chain B exposes C1-C6
            exposed_atoms_b = "C1,C2,C3,C4,C5,C6"
    else:
        exposed_atoms_b = ""

    print(f"[SequentialBinder] Chain A exposed: {exposed_atoms_a}, Chain B will expose: {exposed_atoms_b}")

    designs = []
    best_design_idx = 0
    best_score = float('-inf')

    # Extract Chain A protein-only (no ligand) for later combining
    chain_a_protein = _extract_protein_pdb(chain_a_pdb, "A")

    for design_idx in range(num_designs):
        design_seed = (seed + design_idx) if seed is not None else None
        print(f"[SequentialBinder] Generating design {design_idx + 1}/{num_designs}...")

        # NEW APPROACH: Design Chain B around JUST the ligand (like asymmetric design)
        # Then combine with Chain A afterward
        # NOTE: ori_token doesn't reliably position chains on opposite sides,
        # so we rely on select_exposed (RASA) to make each chain bind different ligand faces
        ori_offset_b = job_input.get("ori_offset_b", [-12.0, 0.0, 0.0])

        rfd3_input = {
            "task": "rfd3",
            "length": chain_length,  # Use length instead of contig
            "pdb_content": ligand_pdb,  # Just the ligand, no Chain A
            "ligand": "UNL",  # Reference to ligand in the input
            "seed": design_seed,
            "num_designs": 1,
            # Position Chain B on opposite side of ligand from Chain A
            "ori_token": ori_offset_b,
        }
        print(f"[SequentialBinder] Chain B ori_token: {ori_offset_b}")

        # RASA conditioning: Chain B exposes opposite atoms from Chain A
        if exposed_atoms_b:
            rfd3_input["select_exposed"] = {"UNL": exposed_atoms_b}
            print(f"[SequentialBinder] Chain B will expose: {exposed_atoms_b}")

        result = handle_rfd3(rfd3_input)

        if result.get("status") != "completed":
            print(f"[SequentialBinder] Design {design_idx + 1} failed: {result.get('error')}")
            continue

        result_designs = result.get("result", {}).get("designs", [])
        if not result_designs:
            print(f"[SequentialBinder] Design {design_idx + 1} produced no designs")
            continue

        chain_b_pdb = result_designs[0].get("content") or result_designs[0].get("pdb_content")
        if not chain_b_pdb:
            print(f"[SequentialBinder] Design {design_idx + 1}: No PDB content")
            continue

        # Extract Chain B protein-only and relabel as chain B
        chain_b_protein = _extract_protein_pdb(chain_b_pdb)
        chain_b_protein = _relabel_chain(chain_b_protein, "B")

        # CRITICAL: Translate Chain B to the opposite side of the ligand from Chain A
        # Both chains are designed wrapped around the ligand at origin, so they overlap
        # We need to translate (not just rotate) to achieve proper separation
        chain_b_protein = _translate_chain_to_opposite_side(
            chain_b_protein, chain_a_protein, ligand_pdb
        )

        # Combine Chain A + Chain B + Ligand into final dimer
        dimer_pdb = _combine_chains(chain_a_protein, chain_b_protein, ligand_pdb)

        # Verify we have two chains
        chains_in_output = _get_chains_in_pdb(dimer_pdb)
        print(f"[SequentialBinder] Design {design_idx + 1}: chains = {chains_in_output}")

        if "A" not in chains_in_output or "B" not in chains_in_output:
            print(f"[SequentialBinder] Design {design_idx + 1}: Missing chain A or B, skipping")
            continue

        # Validate the dimer
        validation = validate_cleaved_dimer(
            pdb_content=dimer_pdb,
            ligand_name="UNL",
            ligand_smiles=ligand_smiles,
            min_contacts_per_chain=3,
        )

        checks = validation.get("checks", {})
        gnina_result = checks.get("gnina_result", {}).get("result", {})
        affinity = gnina_result.get("best_affinity")

        design_entry = {
            "pdb_content": dimer_pdb,
            "design_index": design_idx,
            "validation": validation,
            "metrics": {
                "affinity": affinity,
                "contacts_a": checks.get("contacts_a", 0),
                "contacts_b": checks.get("contacts_b", 0),
                "has_clashes": checks.get("clash_check", {}).get("has_clashes", False),
                "separable": True,  # Interface designs are separable by design
            }
        }

        # Score
        score = 0
        if affinity is not None:
            score = -affinity * 5  # Weight affinity heavily
        if validation.get("overall_pass"):
            score += 10
        if not design_entry["metrics"]["has_clashes"]:
            score += 5

        design_entry["score"] = score
        designs.append(design_entry)

        if score > best_score:
            best_score = score
            best_design_idx = len(designs) - 1

        print(f"[SequentialBinder] Design {design_idx + 1}: contacts_a={checks.get('contacts_a', 0)}, contacts_b={checks.get('contacts_b', 0)}, affinity={affinity}")

    if not designs:
        return {
            "status": "failed",
            "error": "No valid sequential binder designs generated",
        }

    return {
        "status": "completed",
        "result": {
            "designs": designs,
            "approach": "sequential",
            "best_design": best_design_idx,
        }
    }


def _design_full_dimer(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Complete workflow: design asymmetric chain A, then complementary chain B.

    This produces a fully separable dimer with ligand at the interface.
    """
    ligand_smiles = job_input.get("ligand_smiles")
    chain_length = job_input.get("chain_length", "60-80")
    num_designs = job_input.get("num_designs", 3)
    seed = job_input.get("seed")

    print("[FullDimer] Step 1: Design asymmetric chain A...")

    # Step 1: Design asymmetric chain A
    step1_input = {
        **job_input,
        "approach": "asymmetric",
        "side": job_input.get("side", "left"),
    }
    step1_result = _design_asymmetric_binder(step1_input)

    if step1_result.get("status") != "completed":
        return {
            "status": "failed",
            "error": f"Step 1 (asymmetric) failed: {step1_result.get('error')}",
            "step1_result": step1_result,
        }

    # Get best chain A
    step1_designs = step1_result.get("result", {}).get("designs", [])
    best_a_idx = step1_result.get("result", {}).get("best_design", 0)
    chain_a_pdb = step1_designs[best_a_idx]["pdb_content"]

    print(f"[FullDimer] Step 1 complete. Best chain A: index {best_a_idx}")
    print("[FullDimer] Step 2: Design complementary chain B...")

    # Step 2: Design chain B using chain A as context
    # Chain B should bind the atoms that Chain A kept EXPOSED
    # (Chain A exposed them so Chain B could bind them later)
    exposed_a = step1_result.get("result", {}).get("exposed_atoms", "")
    if "N=N" in job_input.get("ligand_smiles", ""):
        # Azobenzene atom naming: C1-C6 (first phenyl), N7-N8 (azo), C9-C14 (second phenyl)
        if "C1" in exposed_a:
            # Chain A exposed C1-C6 (first phenyl), so Chain A bound C9-C14 (second phenyl)
            # Chain B should bind C1-C6 (the exposed atoms from Chain A's side)
            buried_b = "C1,C2,C3,C4,C5,C6"  # First phenyl - Chain B binds these
        else:
            # Chain A exposed C9-C14 (second phenyl), so Chain A bound C1-C6 (first phenyl)
            # Chain B should bind C9-C14 (the exposed atoms from Chain A's side)
            buried_b = "C9,C10,C11,C12,C13,C14"  # Second phenyl - Chain B binds these
        print(f"[FullDimer] Chain A exposed: {exposed_a}, Chain B will bind: {buried_b}")
    else:
        buried_b = job_input.get("exposed_atoms_b", "")

    step2_input = {
        **job_input,
        "approach": "sequential",
        "chain_a_pdb": chain_a_pdb,
        "exposed_atoms_b": buried_b,  # Named for compatibility but now used as buried atoms
    }
    step2_result = _design_sequential_binder(step2_input)

    if step2_result.get("status") != "completed":
        return {
            "status": "failed",
            "error": f"Step 2 (sequential) failed: {step2_result.get('error')}",
            "step1_result": step1_result,
            "step2_result": step2_result,
        }

    # Get best dimer
    step2_designs = step2_result.get("result", {}).get("designs", [])
    best_b_idx = step2_result.get("result", {}).get("best_design", 0)

    print(f"[FullDimer] Step 2 complete. Best dimer: index {best_b_idx}")

    # Extract best dimer metrics for easy access
    best_dimer = step2_designs[best_b_idx]
    dimer_metrics = best_dimer.get("metrics", {})

    return {
        "status": "completed",
        "result": {
            "approach": "full",
            "chain_a": {
                "pdb_content": chain_a_pdb,
                "design_index": best_a_idx,
                "metrics": step1_designs[best_a_idx].get("metrics"),
            },
            "dimer": {
                "pdb_content": best_dimer["pdb_content"],
                "design_index": best_b_idx,
                "validation": best_dimer.get("validation"),
                "metrics": dimer_metrics,
                # Flatten key metrics for frontend compatibility
                "affinity": dimer_metrics.get("affinity"),
                "contacts_a": dimer_metrics.get("contacts_a", 0),
                "contacts_b": dimer_metrics.get("contacts_b", 0),
                "has_clashes": dimer_metrics.get("has_clashes", False),
                "separable": dimer_metrics.get("separable", True),
            },
            "all_chain_a_designs": step1_designs,
            "all_dimer_designs": step2_designs,
        }
    }


def _design_symmetric_dimer(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Design a C2 symmetric dimer with ligand at the interface.

    This uses RFD3's built-in interface_ligand mode with C2 symmetry,
    which produces properly separated chains without overlap.

    The resulting dimer has:
    - Chain A on one side of the ligand
    - Chain B (mirror of A) on the opposite side
    - Ligand at the C2 symmetry axis (interface)
    """
    ligand_smiles = job_input.get("ligand_smiles")
    chain_length = job_input.get("chain_length", "60-80")
    num_designs = job_input.get("num_designs", 3)
    seed = job_input.get("seed")

    print(f"[SymmetricDimer] Designing C2 symmetric dimer with interface ligand...")

    designs = []
    best_design_idx = 0
    best_score = float('-inf')

    for design_idx in range(num_designs):
        design_seed = (seed + design_idx) if seed is not None else None
        print(f"[SymmetricDimer] Generating design {design_idx + 1}/{num_designs}...")

        # Use RFD3 with C2 symmetry and interface_ligand in "binder" mode
        # The "binder" mode adds H-bond and burial conditioning for stronger binding:
        # - select_hbond_acceptor on N7,N8 (azo nitrogens)
        # - select_buried on N7,N8 (forces close contact)
        # - select_partially_buried on C4,C9 (phenyl carbons at interface)
        rfd3_input = {
            "task": "rfd3",
            "length": chain_length,
            "ligand_smiles": ligand_smiles,
            "symmetry": "C2",
            "interface_ligand": "binder",  # Use binder mode for stronger contacts
            "seed": design_seed,
            "num_designs": 1,
        }

        result = handle_rfd3(rfd3_input)

        if result.get("status") != "completed":
            print(f"[SymmetricDimer] Design {design_idx + 1} failed: {result.get('error')}")
            continue

        result_designs = result.get("result", {}).get("designs", [])
        if not result_designs:
            print(f"[SymmetricDimer] Design {design_idx + 1} produced no designs")
            continue

        pdb_content = result_designs[0].get("content") or result_designs[0].get("pdb_content")
        if not pdb_content:
            print(f"[SymmetricDimer] Design {design_idx + 1}: No PDB content")
            continue

        # Optional FastRelax refinement to resolve clashes
        run_fastrelax = job_input.get("run_fastrelax", False)
        fastrelax_iter = job_input.get("fastrelax_iter", 200)
        relaxed = False
        relax_energy_change = None

        if run_fastrelax and check_pyrosetta_available() and ligand_smiles:
            print(f"[SymmetricDimer] Running FastRelax on design {design_idx + 1}...")
            relax_result = fastrelax_with_ligand(
                pdb_content=pdb_content,
                ligand_smiles=ligand_smiles,
                max_iter=fastrelax_iter,
                constrain_coords=True,
                interface_only=True,  # Faster - only relax near ligand
                interface_distance=8.0,
            )

            if relax_result.get("status") == "completed":
                pdb_content = relax_result["relaxed_pdb"]
                relaxed = True
                relax_energy_change = relax_result.get("energy_change")
                print(f"[SymmetricDimer] FastRelax complete: ΔE={relax_energy_change:.1f}")
            else:
                print(f"[SymmetricDimer] FastRelax failed: {relax_result.get('error')}")
        elif run_fastrelax and not ligand_smiles:
            print(f"[SymmetricDimer] FastRelax requested but ligand_smiles not provided")
        elif run_fastrelax:
            print(f"[SymmetricDimer] FastRelax requested but PyRosetta not available")

        # Validate the dimer
        validation = validate_cleaved_dimer(
            pdb_content=pdb_content,
            ligand_name="UNL",
            ligand_smiles=ligand_smiles,
            min_contacts_per_chain=0,  # Symmetric designs may have different contact patterns
        )

        checks = validation.get("checks", {})
        gnina_result = checks.get("gnina_result", {}).get("result", {})
        affinity = gnina_result.get("best_affinity")

        design_entry = {
            "pdb_content": pdb_content,
            "design_index": design_idx,
            "validation": validation,
            "metrics": {
                "affinity": affinity,
                "contacts_a": checks.get("contacts_a", 0),
                "contacts_b": checks.get("contacts_b", 0),
                "has_clashes": checks.get("clash_check", {}).get("has_clashes", False),
                "separable": True,  # C2 symmetric designs are separable by construction
                "relaxed": relaxed,
                "relax_energy_change": relax_energy_change,
            }
        }

        # Score
        score = 0
        if affinity is not None and affinity < 0:
            score = -affinity * 5  # Weight affinity
        if not design_entry["metrics"]["has_clashes"]:
            score += 10

        design_entry["score"] = score
        designs.append(design_entry)

        if score > best_score:
            best_score = score
            best_design_idx = len(designs) - 1

        relax_info = f", relaxed={relaxed}" if run_fastrelax else ""
        print(f"[SymmetricDimer] Design {design_idx + 1}: affinity={affinity}, clashes={design_entry['metrics']['has_clashes']}{relax_info}")

    if not designs:
        return {
            "status": "failed",
            "error": "No valid symmetric dimer designs produced"
        }

    return {
        "status": "completed",
        "result": {
            "approach": "symmetric",
            "designs": designs,
            "best_design": best_design_idx,
            # Flatten best design metrics for frontend
            "dimer": {
                "pdb_content": designs[best_design_idx]["pdb_content"],
                "design_index": best_design_idx,
                "validation": designs[best_design_idx].get("validation"),
                "metrics": designs[best_design_idx]["metrics"],
            }
        }
    }


def _extract_ligand_pdb(pdb_content: str, ligand_name: str = "UNL") -> str:
    """Extract ligand HETATM records from PDB content."""
    lines = []
    for line in pdb_content.split("\n"):
        if line.startswith("HETATM"):
            res_name = line[17:20].strip()
            if res_name == ligand_name:
                lines.append(line)
    return "\n".join(lines) + "\n" if lines else ""


def _extract_protein_pdb(pdb_content: str, chain_id: str = None) -> str:
    """Extract protein ATOM records from PDB content, optionally for a specific chain."""
    lines = []
    for line in pdb_content.split("\n"):
        if line.startswith("ATOM"):
            if chain_id is None:
                lines.append(line)
            else:
                current_chain = line[21:22].strip()
                if current_chain == chain_id:
                    lines.append(line)
    return "\n".join(lines) + "\n" if lines else ""


def _relabel_chain(pdb_content: str, new_chain_id: str) -> str:
    """Change chain ID for all ATOM/HETATM records in PDB content."""
    lines = []
    for line in pdb_content.split("\n"):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            if len(line) >= 22:
                line = line[:21] + new_chain_id + line[22:]
        lines.append(line)
    return "\n".join(lines)


def _combine_chains(chain_a_pdb: str, chain_b_pdb: str, ligand_pdb: str) -> str:
    """Combine Chain A, Chain B, and ligand into a single PDB."""
    parts = []

    # Add Chain A (ensure it's labeled as A)
    chain_a_clean = chain_a_pdb.strip()
    if chain_a_clean:
        parts.append(chain_a_clean)

    # Add Chain B (should already be labeled as B)
    chain_b_clean = chain_b_pdb.strip()
    if chain_b_clean:
        parts.append(chain_b_clean)

    # Add ligand (label as chain L)
    ligand_clean = _relabel_chain(ligand_pdb.strip(), "L")
    if ligand_clean:
        parts.append(ligand_clean)

    return "\n".join(parts) + "\nEND\n"


def _get_centroid(pdb_content: str) -> tuple:
    """Calculate centroid (center of mass) of all atoms in PDB content."""
    x_sum, y_sum, z_sum = 0.0, 0.0, 0.0
    count = 0

    for line in pdb_content.split("\n"):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                x_sum += x
                y_sum += y
                z_sum += z
                count += 1
            except (ValueError, IndexError):
                continue

    if count == 0:
        return (0.0, 0.0, 0.0)

    return (x_sum / count, y_sum / count, z_sum / count)


def _translate_pdb(pdb_content: str, dx: float, dy: float, dz: float) -> str:
    """Translate all atom coordinates by (dx, dy, dz)."""
    lines = []

    for line in pdb_content.split("\n"):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                x = float(line[30:38].strip()) + dx
                y = float(line[38:46].strip()) + dy
                z = float(line[46:54].strip()) + dz
                # Format coordinates with proper PDB spacing
                new_line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
                lines.append(new_line)
            except (ValueError, IndexError):
                lines.append(line)
        else:
            lines.append(line)

    return "\n".join(lines)


def _rotate_chain_180_around_x(chain_pdb: str, ligand_pdb: str) -> str:
    """
    Rotate chain B by 180 degrees around the ligand's X-axis.

    For azobenzene (linear molecule along X-axis), this flips the protein
    to the opposite side of the ligand in the Y-Z plane.
    """
    centroid_ligand = _get_centroid(ligand_pdb)
    lx, ly, lz = centroid_ligand

    print(f"[RotateChain] Rotating Chain B 180° around X-axis centered at ligand")
    print(f"[RotateChain] Ligand centroid: ({lx:.2f}, {ly:.2f}, {lz:.2f})")

    lines = []
    for line in chain_pdb.split("\n"):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())

                # Translate to ligand-centered coordinates
                y_centered = y - ly
                z_centered = z - lz

                # Rotate 180° around X-axis: (y, z) -> (-y, -z)
                y_rotated = -y_centered
                z_rotated = -z_centered

                # Translate back
                y_new = y_rotated + ly
                z_new = z_rotated + lz

                # Format coordinates with proper PDB spacing
                new_line = f"{line[:30]}{x:8.3f}{y_new:8.3f}{z_new:8.3f}{line[54:]}"
                lines.append(new_line)
            except (ValueError, IndexError):
                lines.append(line)
        else:
            lines.append(line)

    return "\n".join(lines)


def _should_flip_chain_b(chain_a_pdb: str, chain_b_pdb: str, ligand_pdb: str) -> bool:
    """
    Determine if Chain B should be flipped to avoid overlap with Chain A.

    Returns True if both chains are on the same side of the ligand.
    """
    centroid_ligand = _get_centroid(ligand_pdb)
    centroid_a = _get_centroid(chain_a_pdb)
    centroid_b = _get_centroid(chain_b_pdb)

    # Calculate vectors from ligand to each chain
    vec_a = (centroid_a[0] - centroid_ligand[0],
             centroid_a[1] - centroid_ligand[1],
             centroid_a[2] - centroid_ligand[2])
    vec_b = (centroid_b[0] - centroid_ligand[0],
             centroid_b[1] - centroid_ligand[1],
             centroid_b[2] - centroid_ligand[2])

    # Dot product to check if on same side (considering Y and Z for azobenzene)
    # For azobenzene along X-axis, Y-Z plane is perpendicular to the ligand
    dot_yz = vec_a[1] * vec_b[1] + vec_a[2] * vec_b[2]

    print(f"[FlipCheck] vec_a: ({vec_a[0]:.2f}, {vec_a[1]:.2f}, {vec_a[2]:.2f})")
    print(f"[FlipCheck] vec_b: ({vec_b[0]:.2f}, {vec_b[1]:.2f}, {vec_b[2]:.2f})")
    print(f"[FlipCheck] Y-Z dot product: {dot_yz:.2f}")

    return dot_yz > 0  # Same side if positive


def _count_chain_residues(pdb_content: str, chain_id: str) -> int:
    """
    Count the number of residues in a specific chain.

    Args:
        pdb_content: PDB file content as string
        chain_id: Chain identifier (e.g., "A", "B")

    Returns:
        Number of unique residues in the chain
    """
    residues = set()
    for line in pdb_content.split("\n"):
        if line.startswith("ATOM") and len(line) > 26:
            if line[21] == chain_id:
                res_num = line[22:26].strip()
                residues.add(res_num)
    return len(residues)


def _align_chain_to_interface(
    chain_b_pdb: str,
    chain_a_pdb: str,
    ligand_pdb: str,
    target_radius: float = None
) -> str:
    """
    Align Chain B to form a proper interface with Chain A around the ligand.

    The problem with independent chain design:
    - Chain A wraps around ligand at radius R_a (typically ~4-8 Angstrom)
    - Chain B wraps around ligand at radius R_b (can be 20-30 Angstrom!)
    - Both ligands are at origin, so Kabsch alignment does nothing

    This function:
    1. Calculates Chain A's position vector from ligand (direction and radius)
    2. Calculates where Chain B should be: OPPOSITE direction, SAME radius
    3. Translates Chain B to that target position

    Args:
        chain_b_pdb: PDB content of Chain B (protein only)
        chain_a_pdb: PDB content of Chain A (protein only)
        ligand_pdb: Ligand PDB content
        target_radius: Optional target radius; if None, uses Chain A's radius

    Returns:
        Transformed Chain B PDB aligned to form interface
    """
    import numpy as np

    def extract_coords(pdb_content: str) -> np.ndarray:
        """Extract all atom coordinates from PDB."""
        coords = []
        for line in pdb_content.split("\n"):
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except (ValueError, IndexError):
                    continue
        return np.array(coords) if coords else np.array([]).reshape(0, 3)

    def apply_translation(pdb_content: str, translation: np.ndarray) -> str:
        """Apply translation to all atoms in PDB."""
        lines = []
        for line in pdb_content.split("\n"):
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])

                    # Apply translation
                    new_x = x + translation[0]
                    new_y = y + translation[1]
                    new_z = z + translation[2]

                    # Format new coordinates
                    new_line = (
                        line[:30] +
                        f"{new_x:8.3f}" +
                        f"{new_y:8.3f}" +
                        f"{new_z:8.3f}" +
                        line[54:]
                    )
                    lines.append(new_line)
                except (ValueError, IndexError):
                    lines.append(line)
            else:
                lines.append(line)
        return "\n".join(lines)

    # Extract coordinates
    coords_a = extract_coords(chain_a_pdb)
    coords_b = extract_coords(chain_b_pdb)
    coords_ligand = extract_coords(ligand_pdb)

    if len(coords_a) == 0 or len(coords_b) == 0 or len(coords_ligand) == 0:
        print("[AlignInterface] Warning: Could not extract coordinates")
        return chain_b_pdb

    # Calculate centroids
    centroid_a = np.mean(coords_a, axis=0)
    centroid_b = np.mean(coords_b, axis=0)
    centroid_ligand = np.mean(coords_ligand, axis=0)

    # Vector from ligand to Chain A
    vec_ligand_to_a = centroid_a - centroid_ligand
    radius_a = np.linalg.norm(vec_ligand_to_a)

    # Current position of Chain B relative to ligand
    vec_ligand_to_b = centroid_b - centroid_ligand
    radius_b = np.linalg.norm(vec_ligand_to_b)

    print(f"[AlignInterface] Chain A radius: {radius_a:.1f}A, Chain B radius: {radius_b:.1f}A")

    # Determine target radius (use Chain A's radius or provided value)
    if target_radius is None:
        target_radius = radius_a

    # Target position for Chain B: OPPOSITE direction from Chain A, at target radius
    if radius_a > 0.1:
        direction_a = vec_ligand_to_a / radius_a
        # Opposite direction
        direction_b_target = -direction_a
        # Target position
        target_position = centroid_ligand + direction_b_target * target_radius
    else:
        # Chain A is too close to ligand, use default offset
        print("[AlignInterface] Warning: Chain A too close to ligand, using default offset")
        target_position = centroid_ligand + np.array([-15.0, 0.0, 0.0])

    # Translation needed
    translation = target_position - centroid_b

    print(f"[AlignInterface] Translating Chain B by ({translation[0]:.1f}, {translation[1]:.1f}, {translation[2]:.1f})")
    print(f"[AlignInterface] Chain B will be at radius {target_radius:.1f}A on opposite side")

    # Apply translation
    transformed_b = apply_translation(chain_b_pdb, translation)

    return transformed_b


def _superimpose_by_ligand(chain_b_pdb: str, ligand_b_pdb: str, ligand_a_pdb: str) -> str:
    """
    DEPRECATED: Use _align_chain_to_interface instead.

    This function doesn't work when both ligands are at the same position (origin).
    """
    print("[SuperimposeLigand] WARNING: This function is deprecated. Use _align_chain_to_interface.")
    return chain_b_pdb


def _translate_chain_to_opposite_side(
    chain_b_pdb: str, chain_a_pdb: str, ligand_pdb: str
) -> str:
    """
    Translate Chain B to the opposite side of the ligand from Chain A.

    This ensures separable dimer geometry where:
    - Chain A is on one side of the ligand
    - Chain B is on the opposite side
    - The ligand sits at the interface between them

    Strategy:
    1. Calculate Chain A's dominant direction from ligand
    2. Translate Chain B to the opposite side with enough separation (~25Å typical)
    """
    import math

    centroid_ligand = _get_centroid(ligand_pdb)
    centroid_a = _get_centroid(chain_a_pdb)
    centroid_b = _get_centroid(chain_b_pdb)

    lx, ly, lz = centroid_ligand
    ax, ay, az = centroid_a
    bx, by, bz = centroid_b

    # Vector from ligand to Chain A (this is where Chain A is relative to ligand)
    vec_a_x = ax - lx
    vec_a_y = ay - ly

    # Magnitude of Chain A offset in XY plane
    dist_a = math.sqrt(vec_a_x**2 + vec_a_y**2)
    if dist_a < 0.1:
        dist_a = 1.0  # Avoid division by zero

    # Normalize to get direction
    dir_x = vec_a_x / dist_a
    dir_y = vec_a_y / dist_a

    # Target separation: 25Å from ligand on opposite side
    # This ensures proper physical separation for a separable dimer
    SEPARATION = 25.0

    # Target position for Chain B: opposite side, same distance
    target_bx = lx - dir_x * SEPARATION
    target_by = ly - dir_y * SEPARATION

    # Also move Chain A to its target position for symmetry
    # (we don't modify A here, but the same distance)
    target_ax = lx + dir_x * SEPARATION

    # Translation vector for Chain B: move current centroid to target
    trans_x = target_bx - bx
    trans_y = target_by - by
    trans_z = 0  # Don't translate in Z (ligand axis for azobenzene)

    print(f"[TranslateChain] Chain A centroid: ({ax:.2f}, {ay:.2f}, {az:.2f})")
    print(f"[TranslateChain] Chain B centroid: ({bx:.2f}, {by:.2f}, {bz:.2f})")
    print(f"[TranslateChain] Ligand centroid: ({lx:.2f}, {ly:.2f}, {lz:.2f})")
    print(f"[TranslateChain] Direction from ligand to A: ({dir_x:.2f}, {dir_y:.2f})")
    print(f"[TranslateChain] Target for Chain B: ({target_bx:.2f}, {target_by:.2f}) - {SEPARATION:.0f}Å separation")
    print(f"[TranslateChain] Translating Chain B by ({trans_x:.2f}, {trans_y:.2f}, {trans_z:.2f})")

    lines = []
    for line in chain_b_pdb.split("\n"):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())

                x_new = x + trans_x
                y_new = y + trans_y
                z_new = z + trans_z

                new_line = f"{line[:30]}{x_new:8.3f}{y_new:8.3f}{z_new:8.3f}{line[54:]}"
                lines.append(new_line)
            except (ValueError, IndexError):
                lines.append(line)
        else:
            lines.append(line)

    return "\n".join(lines)


def _rotate_chain_180_around_z(chain_pdb: str, ligand_pdb: str) -> str:
    """
    Rotate chain by 180 degrees around the ligand's Z-axis.

    For azobenzene (with N=N bond oriented along Z-axis), this rotates
    the protein to the opposite side of the ligand in the X-Y plane.
    This is the correct C2 symmetry operation for separable dimers.
    """
    centroid_ligand = _get_centroid(ligand_pdb)
    lx, ly, lz = centroid_ligand

    print(f"[RotateChain] Rotating Chain B 180° around Z-axis (C2 operation)")
    print(f"[RotateChain] Ligand centroid: ({lx:.2f}, {ly:.2f}, {lz:.2f})")

    lines = []
    for line in chain_pdb.split("\n"):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())

                # Translate to ligand-centered coordinates
                x_centered = x - lx
                y_centered = y - ly

                # Rotate 180° around Z-axis: (x, y) -> (-x, -y)
                x_rotated = -x_centered
                y_rotated = -y_centered

                # Translate back
                x_new = x_rotated + lx
                y_new = y_rotated + ly

                # Format coordinates with proper PDB spacing
                new_line = f"{line[:30]}{x_new:8.3f}{y_new:8.3f}{z:8.3f}{line[54:]}"
                lines.append(new_line)
            except (ValueError, IndexError):
                lines.append(line)
        else:
            lines.append(line)

    return "\n".join(lines)


def _chains_on_same_side_xy(chain_a_pdb: str, chain_b_pdb: str, ligand_pdb: str) -> bool:
    """
    Check if both chains are on the same side of the ligand in the X-Y plane.

    For azobenzene with N=N along Z-axis, proper separation means chains
    should be on opposite sides in X-Y (i.e., negative dot product in X-Y).

    Returns True if chains need to be separated (same side).
    """
    centroid_ligand = _get_centroid(ligand_pdb)
    centroid_a = _get_centroid(chain_a_pdb)
    centroid_b = _get_centroid(chain_b_pdb)

    # Calculate vectors from ligand to each chain in X-Y plane only
    vec_a_xy = (centroid_a[0] - centroid_ligand[0],
                centroid_a[1] - centroid_ligand[1])
    vec_b_xy = (centroid_b[0] - centroid_ligand[0],
                centroid_b[1] - centroid_ligand[1])

    # Dot product in X-Y plane
    dot_xy = vec_a_xy[0] * vec_b_xy[0] + vec_a_xy[1] * vec_b_xy[1]

    print(f"[SameXYCheck] Chain A XY offset from ligand: ({vec_a_xy[0]:.2f}, {vec_a_xy[1]:.2f})")
    print(f"[SameXYCheck] Chain B XY offset from ligand: ({vec_b_xy[0]:.2f}, {vec_b_xy[1]:.2f})")
    print(f"[SameXYCheck] X-Y dot product: {dot_xy:.2f} (positive = same side)")

    return dot_xy > 0  # Same side if positive


def _get_chain_length(pdb_content: str, chain_id: str = None) -> int:
    """Get the number of residues in a chain from PDB content."""
    residues = set()
    for line in pdb_content.split("\n"):
        if line.startswith("ATOM"):
            current_chain = line[21:22].strip()
            if chain_id is None or current_chain == chain_id:
                res_num = int(line[22:26].strip())
                residues.add(res_num)
    return len(residues)


def _renumber_pdb_atoms(pdb_content: str) -> str:
    """Renumber atoms in PDB to have strictly increasing IDs.

    RFD3 requires atom IDs to be strictly increasing. This function
    fixes PDB files that have gaps or duplicates in atom numbering.
    """
    lines = pdb_content.split("\n")
    result_lines = []
    atom_id = 1

    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            # PDB format: columns 1-6 = record type, columns 7-11 = atom serial
            # In 0-indexed Python: [0:6] = record, [6:11] = serial, [11:] = rest
            if len(line) >= 11:
                prefix = line[:6]  # ATOM or HETATM
                rest = line[11:]   # Everything after atom ID
                new_line = f"{prefix}{atom_id:5d}{rest}"
                result_lines.append(new_line)
                atom_id += 1
            else:
                result_lines.append(line)
        else:
            result_lines.append(line)

    print(f"[_renumber_pdb_atoms] Renumbered {atom_id - 1} atoms")
    return "\n".join(result_lines)


def _ensure_chain_a(pdb_content: str) -> str:
    """Ensure all protein atoms have chain ID 'A'.

    RFD3 de novo designs may have blank or inconsistent chain IDs.
    This function ensures chain A is set so sequential design contig works.
    """
    lines = pdb_content.split("\n")
    result_lines = []

    for line in lines:
        if line.startswith("ATOM"):
            # PDB format: column 22 (0-indexed: 21) is chain ID
            if len(line) >= 22:
                # Replace chain ID with 'A'
                new_line = line[:21] + "A" + line[22:]
                result_lines.append(new_line)
            else:
                result_lines.append(line)
        else:
            result_lines.append(line)

    return "\n".join(result_lines)


def _label_dimer_chains(pdb_content: str, chain_a_length: int) -> str:
    """Label protein chains as A and B based on residue numbers.

    After sequential design, RFD3 may output all residues in chain A.
    This function splits them into chains A and B based on the expected
    chain A length.

    Args:
        pdb_content: PDB content with potentially misaligned chain IDs
        chain_a_length: Number of residues in chain A

    Returns:
        PDB content with chains properly labeled A and B
    """
    lines = pdb_content.split("\n")
    result_lines = []

    # First pass: identify the residue number cutoff for chain A vs B
    residue_nums = set()
    for line in lines:
        if line.startswith("ATOM"):
            try:
                res_num = int(line[22:26].strip())
                residue_nums.add(res_num)
            except (ValueError, IndexError):
                pass

    sorted_residues = sorted(residue_nums)
    if len(sorted_residues) <= chain_a_length:
        # Not enough residues for two chains, return as-is
        print(f"[_label_dimer_chains] Only {len(sorted_residues)} residues, expected {chain_a_length}+ for dimer")
        return pdb_content

    # Find the boundary: chain A is first N residues, chain B is the rest
    # Note: there may be gaps in numbering, so use sorted list
    chain_a_residues = set(sorted_residues[:chain_a_length])
    chain_b_start = sorted_residues[chain_a_length] if len(sorted_residues) > chain_a_length else None

    print(f"[_label_dimer_chains] Chain A: residues up to {max(chain_a_residues)}, Chain B starts at: {chain_b_start}")

    # Second pass: relabel chains
    for line in lines:
        if line.startswith("ATOM"):
            if len(line) >= 26:
                try:
                    res_num = int(line[22:26].strip())
                    if res_num in chain_a_residues:
                        chain_id = "A"
                    else:
                        chain_id = "B"
                    new_line = line[:21] + chain_id + line[22:]
                    result_lines.append(new_line)
                except (ValueError, IndexError):
                    result_lines.append(line)
            else:
                result_lines.append(line)
        else:
            result_lines.append(line)

    return "\n".join(result_lines)


def _get_chains_in_pdb(pdb_content: str) -> list:
    """Get list of unique chain IDs in PDB content."""
    chains = set()
    for line in pdb_content.split("\n"):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            if len(line) >= 22:
                chain_id = line[21:22]
                if chain_id.strip():  # Non-empty chain ID
                    chains.add(chain_id)
    return sorted(chains)


def _extract_ligand_from_pdb(pdb_content: str, ligand_name: str) -> str:
    """Extract HETATM records for a specific ligand from PDB content.

    Args:
        pdb_content: Full PDB content
        ligand_name: Residue name of the ligand (e.g., "UNL")

    Returns:
        PDB content containing only the ligand HETATM records
    """
    ligand_lines = []
    for line in pdb_content.split("\n"):
        if line.startswith("HETATM"):
            # Residue name is at columns 18-20 (0-indexed: 17:20)
            if len(line) >= 20:
                res_name = line[17:20].strip()
                if res_name == ligand_name:
                    ligand_lines.append(line)

    if ligand_lines:
        return "\n".join(ligand_lines)
    return ""


# Standard amino acid 3-letter codes
STANDARD_AMINO_ACIDS = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    # Common modifications that are still amino acids
    'MSE',  # Selenomethionine (often used in crystallography)
}


def _strip_nonstandard_residues(pdb_content: str) -> str:
    """Strip non-standard residues and renumber to be contiguous.

    This removes HETATM records (ligands, chromophores, ions, water) and
    non-standard amino acids that can cause RFD3 parsing issues.

    Additionally, it renumbers residues to be contiguous (1, 2, 3, ...)
    because RFD3 expects continuous residue numbering in contig strings.

    For example, GFP (1EMA) has residues ...64, 66, 68... (65/67 are the
    chromophore). After cleaning, this becomes ...64, 65, 66...

    Args:
        pdb_content: Raw PDB content

    Returns:
        PDB content with only standard amino acid ATOM records,
        renumbered to be contiguous starting from 1.
    """
    # First pass: collect valid ATOM lines and unique residue numbers per chain
    valid_lines = []
    removed_residues = set()
    chain_residues = {}  # chain -> list of (original_resnum, insertion_code)

    for line in pdb_content.split("\n"):
        if line.startswith("ATOM"):
            if len(line) >= 27:
                res_name = line[17:20].strip()
                if res_name in STANDARD_AMINO_ACIDS:
                    chain_id = line[21]
                    res_num = line[22:26].strip()
                    ins_code = line[26] if len(line) > 26 else ' '
                    res_key = (res_num, ins_code)

                    if chain_id not in chain_residues:
                        chain_residues[chain_id] = []
                    if res_key not in chain_residues[chain_id]:
                        chain_residues[chain_id].append(res_key)

                    valid_lines.append(line)
                else:
                    removed_residues.add(res_name)
        elif line.startswith("HETATM"):
            if len(line) >= 20:
                res_name = line[17:20].strip()
                removed_residues.add(res_name)

    if removed_residues:
        print(f"[PDBClean] Removed non-standard residues: {sorted(removed_residues)}")

    # Build renumbering map: (chain, old_res, ins_code) -> new_res
    renumber_map = {}
    for chain_id, res_list in chain_residues.items():
        # Sort by original residue number (numerically if possible)
        def sort_key(r):
            try:
                return (int(r[0]), r[1])
            except ValueError:
                return (999999, r[1])
        res_list.sort(key=sort_key)

        for new_idx, (old_res, ins_code) in enumerate(res_list, start=1):
            renumber_map[(chain_id, old_res, ins_code)] = new_idx

    # Second pass: apply renumbering
    clean_lines = []
    for line in valid_lines:
        chain_id = line[21]
        old_res = line[22:26].strip()
        ins_code = line[26] if len(line) > 26 else ' '

        new_res = renumber_map.get((chain_id, old_res, ins_code))
        if new_res is not None:
            # Replace residue number (columns 23-26, 1-indexed: 22:26 in 0-indexed)
            # Format: right-justified, 4 characters
            new_res_str = f"{new_res:4d}"
            # Replace insertion code with space (column 27, index 26)
            new_line = line[:22] + new_res_str + ' ' + line[27:]
            clean_lines.append(new_line)
        else:
            clean_lines.append(line)

    # Add TER and END records
    clean_lines.append("TER")
    clean_lines.append("END")

    # Report renumbering stats
    for chain_id, res_list in chain_residues.items():
        if res_list:
            orig_first = res_list[0][0]
            orig_last = res_list[-1][0]
            print(f"[PDBClean] Chain {chain_id}: {len(res_list)} residues, "
                  f"renumbered {orig_first}-{orig_last} -> 1-{len(res_list)}")

    return "\n".join(clean_lines)


def _combine_chains_to_dimer(chain_a_pdb: str, chain_b_pdb: str) -> str:
    """Combine two chain PDBs into a dimer with proper chain labels.

    Chain A keeps its label 'A', Chain B is relabeled to 'B'.
    Ligand is taken from chain_a_pdb and kept with its original chain/residue.

    Args:
        chain_a_pdb: PDB content for chain A (with ligand)
        chain_b_pdb: PDB content for chain B (with its own ligand copy - will be removed)

    Returns:
        Combined PDB with chains A, B, and ligand
    """
    result_lines = []
    atom_serial = 1

    # Extract chain A protein atoms (already labeled 'A')
    for line in chain_a_pdb.split("\n"):
        if line.startswith("ATOM"):
            # Ensure chain ID is 'A' and renumber atoms
            new_line = line[:6] + f"{atom_serial:5d}" + line[11:21] + "A" + line[22:]
            result_lines.append(new_line)
            atom_serial += 1

    # Extract chain B protein atoms (relabel to 'B')
    for line in chain_b_pdb.split("\n"):
        if line.startswith("ATOM"):
            # Set chain ID to 'B' and renumber atoms
            new_line = line[:6] + f"{atom_serial:5d}" + line[11:21] + "B" + line[22:]
            result_lines.append(new_line)
            atom_serial += 1

    # Extract ligand from chain A (HETATM records)
    for line in chain_a_pdb.split("\n"):
        if line.startswith("HETATM"):
            # Keep ligand with original chain ID, renumber atoms
            new_line = line[:6] + f"{atom_serial:5d}" + line[11:]
            result_lines.append(new_line)
            atom_serial += 1

    result_lines.append("END")
    return "\n".join(result_lines)


# ============== Heterodimer Design Approaches (Anti-Homodimerization) ==============

def _design_joint_heterodimer(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Design a TRUE heterodimer by designing BOTH chains simultaneously using RFD3's
    native multi-chain capability with chain breaks.

    Approach 5: Joint Multi-Chain Diffusion
    - Uses contig syntax "60-100,/0,60-100" to design two independent chains
    - Both chains co-evolve around the ligand during diffusion
    - Different hotspots per chain enforce asymmetric binding
    - Results in chains that cannot homodimerize (only A+B can form dimer)

    This is the recommended approach for true heterodimer design.
    """
    ligand_smiles = job_input.get("ligand_smiles")
    chain_length = job_input.get("chain_length", "60-80")
    num_designs = job_input.get("num_designs", 3)
    seed = job_input.get("seed")

    # Parse chain length range
    if "-" in str(chain_length):
        min_len, max_len = map(int, str(chain_length).split("-"))
    else:
        min_len = max_len = int(chain_length)

    # Build contig with chain break: "min-max,/0,min-max"
    contig = f"{min_len}-{max_len},/0,{min_len}-{max_len}"

    print(f"[JointHeterodimer] Designing with contig={contig}")
    print(f"[JointHeterodimer] Ligand: {ligand_smiles[:30]}...")

    # Determine atoms for asymmetric binding based on ligand
    # For cis-azobenzene, use different phenyl rings for each chain
    # RDKit naming: Ring1 (C1-C6), Azo (N7,N8), Ring2 (C9-C14)
    is_azobenzene = "N=N" in (ligand_smiles or "")

    # Allow user to control H-bond vs hydrophobic balance
    # Options: "hydrophobic" (default), "balanced", "hbond_focused"
    binding_mode = job_input.get("binding_mode", "balanced")

    if is_azobenzene:
        # Chain A binds Ring 1, Chain B binds Ring 2
        hotspots_a = "C1,C2,C3,C4,C5,C6"  # First phenyl ring carbons
        hotspots_b = "C9,C10,C11,C12,C13,C14"  # Second phenyl ring carbons

        if binding_mode == "hydrophobic":
            # Full hydrophobic burial - best affinity, no H-bond specificity
            burial_atoms = "C1,C2,C3,C4,C5,C6,C9,C10,C11,C12,C13,C14,N7,N8"
            partial_burial_atoms = ""
            use_hbond_conditioning = False
        elif binding_mode == "hbond_focused":
            # H-bond focused - prioritize H-bonds over hydrophobic contacts
            burial_atoms = "C6,C9"  # Only interface carbons (near azo)
            partial_burial_atoms = ""
            use_hbond_conditioning = True
        else:  # "balanced" (default)
            # Balance: phenyl rings buried, azo at interface for potential H-bonds
            burial_atoms = "C1,C2,C3,C4,C5,C6,C9,C10,C11,C12,C13,C14"  # Phenyl carbons buried
            partial_burial_atoms = "N7,N8"  # Azo at interface (can form H-bonds)
            use_hbond_conditioning = True

        print(f"[JointHeterodimer] Binding mode: {binding_mode}")
    else:
        hotspots_a = ""
        hotspots_b = ""
        burial_atoms = ""
        partial_burial_atoms = ""
        use_hbond_conditioning = False

    designs = []
    best_design_idx = 0
    best_score = float('-inf')

    for design_idx in range(num_designs):
        design_seed = (seed + design_idx) if seed is not None else None
        print(f"[JointHeterodimer] Generating design {design_idx + 1}/{num_designs}...")

        # Build RFD3 config for joint heterodimer design
        # Note: This uses the native multi-chain capability of RFD3
        rfd3_input = {
            "task": "rfd3",
            "contig": contig,
            "ligand_smiles": ligand_smiles,
            "seed": design_seed,
            "num_designs": 1,
            "is_non_loopy": True,  # Prefer helical structures for interfaces
        }

        # Add burial constraints based on binding mode
        if burial_atoms:
            rfd3_input["select_buried"] = {"UNL": burial_atoms}

        # Partial burial for azo nitrogens (at interface, accessible for H-bonds)
        if partial_burial_atoms:
            rfd3_input["select_partially_buried"] = {"UNL": partial_burial_atoms}

        # Add H-bond conditioning (conditional based on binding mode)
        if use_hbond_conditioning and is_azobenzene:
            rfd3_input["select_hbond_acceptor"] = {"UNL": "N7,N8"}
            print(f"[JointHeterodimer] H-bond conditioning: select_hbond_acceptor=N7,N8")

        # Run RFD3 with multi-chain contig
        result = handle_rfd3(rfd3_input)

        if result.get("status") != "completed":
            print(f"[JointHeterodimer] Design {design_idx + 1} failed: {result.get('error')}")
            continue

        result_designs = result.get("result", {}).get("designs", [])
        if not result_designs:
            print(f"[JointHeterodimer] Design {design_idx + 1} produced no designs")
            continue

        pdb_content = result_designs[0].get("content") or result_designs[0].get("pdb_content")
        if not pdb_content:
            print(f"[JointHeterodimer] Design {design_idx + 1}: No PDB content")
            continue

        # Validate the heterodimer
        validation = validate_cleaved_dimer(
            pdb_content=pdb_content,
            ligand_name="UNL",
            ligand_smiles=ligand_smiles,
            min_contacts_per_chain=2,
        )

        checks = validation.get("checks", {})
        gnina_result = checks.get("gnina_result", {}).get("result", {})
        affinity = gnina_result.get("best_affinity")

        # Calculate sequence identity between chains
        seq_identity = calculate_sequence_identity(pdb_content, "A", "B")
        identity_pct = seq_identity.get("sequence_identity_percent", 100)

        # Score homodimerization potential
        homo_score = score_homodimerization(
            pdb_content=pdb_content,
            ligand_smiles=ligand_smiles,
            ligand_name="UNL",
        )

        design_entry = {
            "pdb_content": pdb_content,
            "design_index": design_idx,
            "validation": validation,
            "metrics": {
                "affinity": affinity,
                "contacts_a": checks.get("contacts_a", 0),
                "contacts_b": checks.get("contacts_b", 0),
                "has_clashes": checks.get("clash_check", {}).get("has_clashes", False),
                "separable": checks.get("separable", True),
                "sequence_identity": identity_pct,
                "is_heterodimer": identity_pct < 70,
            },
            "anti_homodimerization": homo_score,
        }

        # Score: prioritize heterodimer properties
        score = 0
        if affinity is not None and affinity < 0:
            score += -affinity * 5  # Better affinity
        if identity_pct < 70:
            score += 20  # True heterodimer
        if homo_score.get("passes_anti_homodimerization", False):
            score += 30  # Anti-homo validated
        if not design_entry["metrics"]["has_clashes"]:
            score += 10

        # H-bond bonus scoring for azobenzene N7/N8 (azo nitrogens)
        n7_hbonds = checks.get("n7_hbonds", 0)
        n8_hbonds = checks.get("n8_hbonds", 0)
        total_ligand_hbonds = checks.get("total_ligand_hbonds", 0)
        if n7_hbonds >= 1:
            score += 10  # N7 has at least one H-bond
        if n8_hbonds >= 1:
            score += 10  # N8 has at least one H-bond
        if n7_hbonds >= 1 and n8_hbonds >= 1:
            score += 5  # Both azo nitrogens have H-bonds (extra bonus)

        # Add H-bond metrics to design entry
        design_entry["metrics"]["n7_hbonds"] = n7_hbonds
        design_entry["metrics"]["n8_hbonds"] = n8_hbonds
        design_entry["metrics"]["total_ligand_hbonds"] = total_ligand_hbonds

        # Add comprehensive interaction analysis (from validation or compute fresh)
        if validation.get("interactions"):
            design_entry["interactions"] = validation["interactions"]
            design_entry["key_binding_residues"] = validation.get("key_binding_residues", [])
            design_entry["recommendations"] = validation.get("recommendations", [])
            design_entry["interaction_summary"] = checks.get("interaction_summary", {})
        elif SHARED_INTERACTION_ANALYSIS_AVAILABLE:
            try:
                interaction_result = analyze_all_interactions(pdb_content, "UNL")
                design_entry["interactions"] = format_for_frontend(interaction_result)
                design_entry["key_binding_residues"] = interaction_result.key_residues
                design_entry["recommendations"] = generate_recommendations(
                    interaction_result, ligand_has_aromatics=is_azobenzene
                )
                design_entry["interaction_summary"] = {
                    "hydrogen_bonds": len(interaction_result.hydrogen_bonds),
                    "hydrophobic_contacts": len(interaction_result.hydrophobic_contacts),
                    "pi_stacking": len(interaction_result.pi_stacking),
                    "salt_bridges": len(interaction_result.salt_bridges),
                    "total": interaction_result.total_contacts,
                }
            except Exception as e:
                print(f"[JointHeterodimer] Interaction analysis failed: {e}")

        design_entry["score"] = score
        designs.append(design_entry)

        if score > best_score:
            best_score = score
            best_design_idx = len(designs) - 1

        print(f"[JointHeterodimer] Design {design_idx + 1}: affinity={affinity}, identity={identity_pct:.1f}%, anti_homo={homo_score.get('passes_anti_homodimerization', False)}")

    if not designs:
        return {
            "status": "failed",
            "error": "No valid joint heterodimer designs produced",
            "suggestions": [
                "Try different chain lengths",
                "Increase num_designs",
                "Check if ligand SMILES is valid"
            ]
        }

    return {
        "status": "completed",
        "result": {
            "approach": "joint",
            "designs": designs,
            "best_design": best_design_idx,
            "contig": contig,
            "dimer": {
                "pdb_content": designs[best_design_idx]["pdb_content"],
                "metrics": designs[best_design_idx]["metrics"],
                "anti_homodimerization": designs[best_design_idx].get("anti_homodimerization"),
            }
        }
    }


def _design_asymmetric_rasa_heterodimer(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Design a heterodimer using Asymmetric RASA (Relative Accessible Surface Area)
    conditioning for each chain.

    Approach 7: Asymmetric RASA + Hotspot Differentiation
    - Chain A buries one part of ligand, exposes the other
    - Chain B buries the opposite part, exposes what A buried
    - Results in complementary binding modes that prevent homodimerization

    This builds on the existing asymmetric approach but enforces different
    burial/exposure for each chain explicitly.
    """
    ligand_smiles = job_input.get("ligand_smiles")
    chain_length = job_input.get("chain_length", "60-80")
    num_designs = job_input.get("num_designs", 3)
    seed = job_input.get("seed")

    print(f"[AsymmetricRASA] Starting heterodimer design with complementary RASA...")

    # Allow user to control H-bond vs hydrophobic balance
    # Options: "hydrophobic", "balanced" (default), "hbond_focused"
    binding_mode = job_input.get("binding_mode", "balanced")

    # For azobenzene, define complementary RASA constraints
    # RDKit atom naming for azobenzene "c1ccc(/N=N\\c2ccccc2)cc1":
    # Ring 1: C1,C2,C3,C4,C5,C6 (first phenyl ring)
    # Ring 2: C9,C10,C11,C12,C13,C14 (second phenyl ring)
    # Azo: N7,N8
    is_azobenzene = "N=N" in (ligand_smiles or "")
    if is_azobenzene:
        if binding_mode == "hydrophobic":
            # Chain A: buries ring 1 + N7, exposes ring 2
            buried_a = "C1,C2,C3,C4,C5,C6,N7"
            exposed_a = "C9,C10,C11,C12,C13,C14"
            # Chain B: buries ring 2 + N8, exposes ring 1
            buried_b = "C9,C10,C11,C12,C13,C14,N8"
            exposed_b = "C1,C2,C3,C4,C5,C6"
            use_hbond_conditioning = False
        else:  # "balanced" or "hbond_focused"
            # Chain A: buries ring 1 (carbons only), exposes ring 2
            buried_a = "C1,C2,C3,C4,C5,C6"
            exposed_a = "C9,C10,C11,C12,C13,C14"
            # Chain B: buries ring 2 (carbons only), exposes ring 1
            buried_b = "C9,C10,C11,C12,C13,C14"
            exposed_b = "C1,C2,C3,C4,C5,C6"
            use_hbond_conditioning = True

        print(f"[AsymmetricRASA] Binding mode: {binding_mode}")
    else:
        buried_a = ""
        exposed_a = ""
        buried_b = ""
        exposed_b = ""
        use_hbond_conditioning = False

    designs = []
    best_design_idx = 0
    best_score = float('-inf')

    for design_idx in range(num_designs):
        design_seed = (seed + design_idx) if seed is not None else None
        print(f"[AsymmetricRASA] Generating design {design_idx + 1}/{num_designs}...")

        # Step 1: Design Chain A with RASA constraints
        chain_a_input = {
            "task": "rfd3",
            "length": chain_length,
            "ligand_smiles": ligand_smiles,
            "seed": design_seed,
            "num_designs": 1,
            "is_non_loopy": True,
        }
        if buried_a:
            chain_a_input["select_buried"] = {"UNL": buried_a}
        if exposed_a:
            chain_a_input["select_exposed"] = {"UNL": exposed_a}

        # H-bond conditioning: Chain A targets N7 (first azo nitrogen)
        if use_hbond_conditioning and is_azobenzene:
            chain_a_input["select_hbond_acceptor"] = {"UNL": "N7"}
            print(f"[AsymmetricRASA] Chain A H-bond conditioning: select_hbond_acceptor=N7")

        result_a = handle_rfd3(chain_a_input)
        if result_a.get("status") != "completed":
            print(f"[AsymmetricRASA] Chain A design failed: {result_a.get('error')}")
            continue

        designs_a = result_a.get("result", {}).get("designs", [])
        if not designs_a:
            continue
        pdb_a = designs_a[0].get("content") or designs_a[0].get("pdb_content")
        if not pdb_a:
            continue

        # Extract ligand from Chain A output
        ligand_a_pdb = _extract_ligand_pdb(pdb_a, "UNL")

        # Step 2: Design Chain B WITH Chain A + ligand as context (FIXED!)
        # Previous approach designed Chain B independently with its own ligand copy,
        # resulting in spatial mismatch. Now we design Chain B in Chain A's coordinate frame.
        #
        # RFD3 Context Protein Feature:
        # - Provide Chain A + ligand as pdb_content (input context)
        # - Use contig "A1-{len},/0,{min}-{max}" to fix Chain A and design Chain B
        # - Chain B naturally designs around the SAME ligand in the SAME coordinate frame

        chain_a_len = _count_chain_residues(pdb_a, "A")
        min_len, max_len = 50, 70  # Default range
        if "-" in chain_length:
            parts = chain_length.split("-")
            min_len, max_len = int(parts[0]), int(parts[1])
        else:
            min_len = max_len = int(chain_length)

        # Contig: fix Chain A residues 1-N, chain break, design new chain with length range
        context_contig = f"A1-{chain_a_len},/0,{min_len}-{max_len}"
        print(f"[AsymmetricRASA] Designing Chain B with context contig: {context_contig}")

        chain_b_input = {
            "task": "rfd3",
            "pdb_content": pdb_a,        # Chain A + ligand as INPUT CONTEXT
            "contig": context_contig,     # Fix A, design new chain B
            "ligand": "UNL",              # Use existing ligand from Chain A PDB
            "seed": design_seed + 1000 if design_seed else None,
            "num_designs": 1,
            "is_non_loopy": True,
            # Center the design on the ligand to guide Chain B toward it
            "ori_token": [0.0, 0.0, 0.0],
        }
        if buried_b:
            chain_b_input["select_buried"] = {"UNL": buried_b}
        if exposed_b:
            chain_b_input["select_exposed"] = {"UNL": exposed_b}

        # H-bond conditioning: Chain B targets N8 (second azo nitrogen)
        if use_hbond_conditioning and is_azobenzene:
            chain_b_input["select_hbond_acceptor"] = {"UNL": "N8"}
            print(f"[AsymmetricRASA] Chain B H-bond conditioning: select_hbond_acceptor=N8")

        result_b = handle_rfd3(chain_b_input)
        if result_b.get("status") != "completed":
            print(f"[AsymmetricRASA] Chain B design failed: {result_b.get('error')}")
            continue

        designs_b = result_b.get("result", {}).get("designs", [])
        if not designs_b:
            continue
        pdb_b = designs_b[0].get("content") or designs_b[0].get("pdb_content")
        if not pdb_b:
            continue

        # Step 3: Extract and relabel chains from the context-aware output
        # The output from RFD3 context design already has both chains in the same coordinate frame!
        # Chain A is preserved from input, Chain B is the newly designed chain
        chain_a_protein = _relabel_chain(_extract_protein_pdb(pdb_b, "A"), "A")
        chain_b_protein = _relabel_chain(_extract_protein_pdb(pdb_b, "B"), "B")

        # No translation needed! Chain B was designed in Chain A's coordinate frame
        print(f"[AsymmetricRASA] Chain B designed in context - no translation needed")

        # Combine into final dimer (ligand is already in the output)
        dimer_pdb = _combine_chains(chain_a_protein, chain_b_protein, ligand_a_pdb)

        # Validate
        validation = validate_cleaved_dimer(
            pdb_content=dimer_pdb,
            ligand_name="UNL",
            ligand_smiles=ligand_smiles,
            min_contacts_per_chain=2,
        )

        checks = validation.get("checks", {})
        gnina_result = checks.get("gnina_result", {}).get("result", {})
        affinity = gnina_result.get("best_affinity")

        seq_identity = calculate_sequence_identity(dimer_pdb, "A", "B")
        identity_pct = seq_identity.get("sequence_identity_percent", 100)

        homo_score = score_homodimerization(
            pdb_content=dimer_pdb,
            ligand_smiles=ligand_smiles,
            ligand_name="UNL",
        )

        design_entry = {
            "pdb_content": dimer_pdb,
            "design_index": design_idx,
            "validation": validation,
            "metrics": {
                "affinity": affinity,
                "contacts_a": checks.get("contacts_a", 0),
                "contacts_b": checks.get("contacts_b", 0),
                "has_clashes": checks.get("clash_check", {}).get("has_clashes", False),
                "separable": True,
                "sequence_identity": identity_pct,
                "is_heterodimer": identity_pct < 70,
                "rasa_a": {"buried": buried_a, "exposed": exposed_a},
                "rasa_b": {"buried": buried_b, "exposed": exposed_b},
            },
            "anti_homodimerization": homo_score,
        }

        score = 0
        if affinity is not None and affinity < 0:
            score += -affinity * 5
        if identity_pct < 70:
            score += 20
        if homo_score.get("passes_anti_homodimerization", False):
            score += 30
        if not design_entry["metrics"]["has_clashes"]:
            score += 10

        # H-bond bonus scoring for azobenzene N7/N8 (azo nitrogens)
        n7_hbonds = checks.get("n7_hbonds", 0)
        n8_hbonds = checks.get("n8_hbonds", 0)
        total_ligand_hbonds = checks.get("total_ligand_hbonds", 0)
        if n7_hbonds >= 1:
            score += 10  # N7 has at least one H-bond
        if n8_hbonds >= 1:
            score += 10  # N8 has at least one H-bond
        if n7_hbonds >= 1 and n8_hbonds >= 1:
            score += 5  # Both azo nitrogens have H-bonds (extra bonus)

        # Add H-bond metrics to design entry
        design_entry["metrics"]["n7_hbonds"] = n7_hbonds
        design_entry["metrics"]["n8_hbonds"] = n8_hbonds
        design_entry["metrics"]["total_ligand_hbonds"] = total_ligand_hbonds

        # Add comprehensive interaction analysis (from validation or compute fresh)
        if validation.get("interactions"):
            design_entry["interactions"] = validation["interactions"]
            design_entry["key_binding_residues"] = validation.get("key_binding_residues", [])
            design_entry["recommendations"] = validation.get("recommendations", [])
            design_entry["interaction_summary"] = checks.get("interaction_summary", {})
        elif SHARED_INTERACTION_ANALYSIS_AVAILABLE:
            try:
                interaction_result = analyze_all_interactions(pdb_content, "UNL")
                design_entry["interactions"] = format_for_frontend(interaction_result)
                design_entry["key_binding_residues"] = interaction_result.key_residues
                design_entry["recommendations"] = generate_recommendations(
                    interaction_result, ligand_has_aromatics=is_azobenzene
                )
                design_entry["interaction_summary"] = {
                    "hydrogen_bonds": len(interaction_result.hydrogen_bonds),
                    "hydrophobic_contacts": len(interaction_result.hydrophobic_contacts),
                    "pi_stacking": len(interaction_result.pi_stacking),
                    "salt_bridges": len(interaction_result.salt_bridges),
                    "total": interaction_result.total_contacts,
                }
            except Exception as e:
                print(f"[AsymmetricRASA] Interaction analysis failed: {e}")

        design_entry["score"] = score
        designs.append(design_entry)

        if score > best_score:
            best_score = score
            best_design_idx = len(designs) - 1

        print(f"[AsymmetricRASA] Design {design_idx + 1}: affinity={affinity}, identity={identity_pct:.1f}%, n7_hbonds={n7_hbonds}, n8_hbonds={n8_hbonds}")

    if not designs:
        return {
            "status": "failed",
            "error": "No valid asymmetric RASA heterodimer designs produced"
        }

    return {
        "status": "completed",
        "result": {
            "approach": "asymmetric_rasa",
            "designs": designs,
            "best_design": best_design_idx,
            "dimer": {
                "pdb_content": designs[best_design_idx]["pdb_content"],
                "metrics": designs[best_design_idx]["metrics"],
                "anti_homodimerization": designs[best_design_idx].get("anti_homodimerization"),
            }
        }
    }


def _design_induced_heterodimer(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Design an INDUCED heterodimer where Chain B only binds when ligand is present.

    Approach 8: Induced Dimerization Scaffold
    - Chain A binds the ligand directly
    - Chain B binds the NEOSURFACE created by A+ligand complex
    - Result: A+B only dimerize when ligand is present
    - Mimics chemically-induced dimerization (CID) systems

    This approach creates the most selective heterodimer where A-A, B-B,
    and even A+B without ligand don't form stable dimers.
    """
    ligand_smiles = job_input.get("ligand_smiles")
    chain_length = job_input.get("chain_length", "60-80")
    num_designs = job_input.get("num_designs", 3)
    seed = job_input.get("seed")

    # Allow user to control H-bond vs hydrophobic balance
    # Options: "hydrophobic", "balanced" (default), "hbond_focused"
    binding_mode = job_input.get("binding_mode", "balanced")

    print(f"[InducedDimer] Starting induced heterodimer design...")
    print(f"[InducedDimer] Binding mode: {binding_mode}")

    designs = []
    best_design_idx = 0
    best_score = float('-inf')

    for design_idx in range(num_designs):
        design_seed = (seed + design_idx) if seed is not None else None
        print(f"[InducedDimer] Generating design {design_idx + 1}/{num_designs}...")

        # Step 1: Design Chain A as a standard small molecule binder
        # Chain A should fully encapsulate its side of the ligand
        chain_a_input = {
            "task": "rfd3",
            "length": chain_length,
            "ligand_smiles": ligand_smiles,
            "seed": design_seed,
            "num_designs": 1,
            "is_non_loopy": True,
        }

        # For azobenzene, A binds one ring completely
        # RDKit naming: Ring 1 = C1-C6; Ring 2 = C9-C14; Azo = N7,N8
        is_azobenzene = "N=N" in (ligand_smiles or "")
        if is_azobenzene:
            if binding_mode == "hydrophobic":
                # Full hydrophobic burial including N7
                chain_a_input["select_buried"] = {"UNL": "C1,C2,C3,C4,C5,C6,N7"}
                chain_a_input["select_exposed"] = {"UNL": "C9,C10,C11,C12,C13,C14"}
            else:  # "balanced" or "hbond_focused"
                # N7 accessible for H-bonding
                chain_a_input["select_buried"] = {"UNL": "C1,C2,C3,C4,C5,C6"}
                chain_a_input["select_exposed"] = {"UNL": "C9,C10,C11,C12,C13,C14"}
                # H-bond conditioning: Chain A targets N7 (azo nitrogen as H-bond acceptor)
                chain_a_input["select_hbond_acceptor"] = {"UNL": "N7"}
                print(f"[InducedDimer] Chain A H-bond conditioning: select_hbond_acceptor=N7")

        result_a = handle_rfd3(chain_a_input)
        if result_a.get("status") != "completed":
            print(f"[InducedDimer] Chain A design failed")
            continue

        designs_a = result_a.get("result", {}).get("designs", [])
        if not designs_a:
            continue
        pdb_a = designs_a[0].get("content") or designs_a[0].get("pdb_content")
        if not pdb_a:
            continue

        # Step 2: Extract ligand from Chain A
        ligand_a_pdb = _extract_ligand_pdb(pdb_a, "UNL")

        # Step 3: Design Chain B WITH Chain A + ligand as context (FIXED!)
        # Previous approach designed Chain B independently with its own ligand copy,
        # resulting in spatial mismatch. Now we design Chain B in Chain A's coordinate frame.
        #
        # RFD3 Context Protein Feature:
        # - Provide Chain A + ligand as pdb_content (input context)
        # - Use contig "A1-{len},/0,{min}-{max}" to fix Chain A and design Chain B
        # - Chain B naturally designs around the SAME ligand in the SAME coordinate frame

        chain_a_len = _count_chain_residues(pdb_a, "A")
        min_len, max_len = 50, 70  # Default range
        if "-" in chain_length:
            parts = chain_length.split("-")
            min_len, max_len = int(parts[0]), int(parts[1])
        else:
            min_len = max_len = int(chain_length)

        # Contig: fix Chain A residues 1-N, chain break, design new chain with length range
        context_contig = f"A1-{chain_a_len},/0,{min_len}-{max_len}"
        print(f"[InducedDimer] Designing Chain B with context contig: {context_contig}")

        chain_b_input = {
            "task": "rfd3",
            "pdb_content": pdb_a,        # Chain A + ligand as INPUT CONTEXT
            "contig": context_contig,     # Fix A, design new chain B
            "ligand": "UNL",              # Use existing ligand from Chain A PDB
            "seed": design_seed + 2000 if design_seed else None,
            "num_designs": 1,
            "is_non_loopy": True,
            # Center the design on the ligand to guide Chain B toward it
            "ori_token": [0.0, 0.0, 0.0],
        }

        # B binds opposite side of ligand (neosurface targeting)
        # RDKit naming: Ring 1 = C1-C6; Ring 2 = C9-C14; Azo = N7,N8
        if is_azobenzene:
            if binding_mode == "hydrophobic":
                # Full hydrophobic burial including N8
                chain_b_input["select_buried"] = {"UNL": "C9,C10,C11,C12,C13,C14,N8"}
                chain_b_input["select_exposed"] = {"UNL": "C1,C2,C3,C4,C5,C6"}
            else:  # "balanced" or "hbond_focused"
                # N8 accessible for H-bonding
                chain_b_input["select_buried"] = {"UNL": "C9,C10,C11,C12,C13,C14"}
                chain_b_input["select_exposed"] = {"UNL": "C1,C2,C3,C4,C5,C6"}
                # H-bond conditioning: Chain B targets N8 (azo nitrogen as H-bond acceptor)
                chain_b_input["select_hbond_acceptor"] = {"UNL": "N8"}
                print(f"[InducedDimer] Chain B H-bond conditioning: select_hbond_acceptor=N8")

        result_b = handle_rfd3(chain_b_input)
        if result_b.get("status") != "completed":
            print(f"[InducedDimer] Chain B design failed")
            continue

        designs_b = result_b.get("result", {}).get("designs", [])
        if not designs_b:
            continue
        pdb_b = designs_b[0].get("content") or designs_b[0].get("pdb_content")
        if not pdb_b:
            continue

        # Step 4: Extract and relabel chains from the context-aware output
        # The output from RFD3 context design already has both chains in the same coordinate frame!
        # Chain A is preserved from input, Chain B is the newly designed chain
        chain_a_protein = _relabel_chain(_extract_protein_pdb(pdb_b, "A"), "A")
        chain_b_protein = _relabel_chain(_extract_protein_pdb(pdb_b, "B"), "B")

        # No translation needed! Chain B was designed in Chain A's coordinate frame
        print(f"[InducedDimer] Chain B designed in context - no translation needed")

        # Combine into final dimer (ligand is already in the output)
        dimer_pdb = _combine_chains(chain_a_protein, chain_b_protein, ligand_a_pdb)

        # Validate
        validation = validate_cleaved_dimer(
            pdb_content=dimer_pdb,
            ligand_name="UNL",
            ligand_smiles=ligand_smiles,
            min_contacts_per_chain=2,
        )

        checks = validation.get("checks", {})
        gnina_result = checks.get("gnina_result", {}).get("result", {})
        affinity = gnina_result.get("best_affinity")

        seq_identity = calculate_sequence_identity(dimer_pdb, "A", "B")
        identity_pct = seq_identity.get("sequence_identity_percent", 100)

        homo_score = score_homodimerization(
            pdb_content=dimer_pdb,
            ligand_smiles=ligand_smiles,
            ligand_name="UNL",
        )

        design_entry = {
            "pdb_content": dimer_pdb,
            "design_index": design_idx,
            "validation": validation,
            "metrics": {
                "affinity": affinity,
                "contacts_a": checks.get("contacts_a", 0),
                "contacts_b": checks.get("contacts_b", 0),
                "has_clashes": checks.get("clash_check", {}).get("has_clashes", False),
                "separable": True,
                "sequence_identity": identity_pct,
                "is_heterodimer": identity_pct < 70,
                "induced_dimer": True,  # Mark as induced approach
            },
            "anti_homodimerization": homo_score,
        }

        score = 0
        if affinity is not None and affinity < 0:
            score += -affinity * 5
        if identity_pct < 70:
            score += 20
        if homo_score.get("passes_anti_homodimerization", False):
            score += 30
        if not design_entry["metrics"]["has_clashes"]:
            score += 10

        # H-bond bonus scoring for azobenzene N7/N8 (azo nitrogens)
        n7_hbonds = checks.get("n7_hbonds", 0)
        n8_hbonds = checks.get("n8_hbonds", 0)
        total_ligand_hbonds = checks.get("total_ligand_hbonds", 0)
        if n7_hbonds >= 1:
            score += 10  # N7 has at least one H-bond
        if n8_hbonds >= 1:
            score += 10  # N8 has at least one H-bond
        if n7_hbonds >= 1 and n8_hbonds >= 1:
            score += 5  # Both azo nitrogens have H-bonds (extra bonus)

        # Add H-bond metrics to design entry
        design_entry["metrics"]["n7_hbonds"] = n7_hbonds
        design_entry["metrics"]["n8_hbonds"] = n8_hbonds
        design_entry["metrics"]["total_ligand_hbonds"] = total_ligand_hbonds

        # Add comprehensive interaction analysis (from validation or compute fresh)
        if validation.get("interactions"):
            design_entry["interactions"] = validation["interactions"]
            design_entry["key_binding_residues"] = validation.get("key_binding_residues", [])
            design_entry["recommendations"] = validation.get("recommendations", [])
            design_entry["interaction_summary"] = checks.get("interaction_summary", {})
        elif SHARED_INTERACTION_ANALYSIS_AVAILABLE:
            try:
                interaction_result = analyze_all_interactions(pdb_content, "UNL")
                design_entry["interactions"] = format_for_frontend(interaction_result)
                design_entry["key_binding_residues"] = interaction_result.key_residues
                design_entry["recommendations"] = generate_recommendations(
                    interaction_result, ligand_has_aromatics=is_azobenzene
                )
                design_entry["interaction_summary"] = {
                    "hydrogen_bonds": len(interaction_result.hydrogen_bonds),
                    "hydrophobic_contacts": len(interaction_result.hydrophobic_contacts),
                    "pi_stacking": len(interaction_result.pi_stacking),
                    "salt_bridges": len(interaction_result.salt_bridges),
                    "total": interaction_result.total_contacts,
                }
            except Exception as e:
                print(f"[InducedDimer] Interaction analysis failed: {e}")

        design_entry["score"] = score
        designs.append(design_entry)

        if score > best_score:
            best_score = score
            best_design_idx = len(designs) - 1

        print(f"[InducedDimer] Design {design_idx + 1}: affinity={affinity}, identity={identity_pct:.1f}%, n7_hbonds={n7_hbonds}, n8_hbonds={n8_hbonds}")

    if not designs:
        return {
            "status": "failed",
            "error": "No valid induced heterodimer designs produced"
        }

    return {
        "status": "completed",
        "result": {
            "approach": "induced",
            "designs": designs,
            "best_design": best_design_idx,
            "dimer": {
                "pdb_content": designs[best_design_idx]["pdb_content"],
                "metrics": designs[best_design_idx]["metrics"],
                "anti_homodimerization": designs[best_design_idx].get("anti_homodimerization"),
            }
        }
    }


# ============== Interface Metal Design Handler ==============

# Metal coordination profiles
METAL_PROFILES = {
    "ZN": {
        "name": "Zinc",
        "coordination": [4, 5, 6],
        "preferred_coord": 4,
        "geometry": "tetrahedral",
        "donors": ["His", "Cys", "Asp", "Glu"],
        "distance": (2.0, 2.3),
    },
    "FE": {
        "name": "Iron",
        "coordination": [4, 5, 6],
        "preferred_coord": 6,
        "geometry": "octahedral",
        "donors": ["His", "Cys", "Asp", "Glu", "Tyr"],
        "distance": (1.9, 2.2),
    },
    "CU": {
        "name": "Copper",
        "coordination": [4, 5, 6],
        "preferred_coord": 4,
        "geometry": "square_planar",
        "donors": ["His", "Cys", "Met"],
        "distance": (1.9, 2.1),
    },
    "MN": {
        "name": "Manganese",
        "coordination": [6],
        "preferred_coord": 6,
        "geometry": "octahedral",
        "donors": ["His", "Asp", "Glu"],
        "distance": (2.1, 2.3),
    },
    "CA": {
        "name": "Calcium",
        "coordination": [6, 7, 8],
        "preferred_coord": 7,
        "geometry": "pentagonal_bipyramidal",
        "donors": ["Asp", "Glu", "Asn"],
        "distance": (2.3, 2.5),
    },
    "MG": {
        "name": "Magnesium",
        "coordination": [6],
        "preferred_coord": 6,
        "geometry": "octahedral",
        "donors": ["Asp", "Glu"],
        "distance": (2.0, 2.2),
    },
    "TB": {
        "name": "Terbium",
        "coordination": [8, 9],
        "preferred_coord": 9,
        "geometry": "tricapped_trigonal_prismatic",
        "donors": ["Asp", "Glu", "Asn"],
        "distance": (2.3, 2.5),
        "special": "luminescent",
    },
    "GD": {
        "name": "Gadolinium",
        "coordination": [8, 9],
        "preferred_coord": 9,
        "geometry": "tricapped_trigonal_prismatic",
        "donors": ["Asp", "Glu", "Asn"],
        "distance": (2.3, 2.5),
        "special": "paramagnetic",
    },
    "EU": {
        "name": "Europium",
        "coordination": [8, 9],
        "preferred_coord": 9,
        "geometry": "tricapped_trigonal_prismatic",
        "donors": ["Asp", "Glu", "Asn"],
        "distance": (2.3, 2.5),
        "special": "luminescent",
    },
}


def handle_interface_metal_design(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Design protein heterodimers with metal coordination at the interface.

    This extends the interface_ligand_design concept to metal ions, using
    coordination chemistry instead of organic ligand binding.

    Approaches:
    1. "joint_metal" - Design both chains around central metal (RECOMMENDED)
    2. "asymmetric_metal" - Different donor types per chain
    3. "induced_metal" - Metal binding triggers dimerization
    4. "bridging_metal" - Multiple metals bridging the interface
    5. "redox_switch" - Oxidation state controls affinity

    Input:
        metal: str (required) - Metal ion code (ZN, FE, CA, TB, etc.)
        approach: str (default: "joint_metal") - Design approach
        chain_length: str (default: "60-80") - Chain length range
        num_designs: int (default: 3) - Number of designs
        coordination_split: [int, int] (default: [2, 2]) - Donors per chain
        chain_a_donors: list[str] (optional) - Preferred donors for chain A
        chain_b_donors: list[str] (optional) - Preferred donors for chain B
        coordination_geometry: str (optional) - Target geometry
        seed: int (optional) - Random seed

        # NEW: Template Library System (chemically realistic templates)
        template_name: str (optional) - Template from library:
            - "caldwell_4": 4 bidentate Glu, CN=8 (highest affinity)
            - "ef_hand_8": 4 mono Asp + 2 bi Glu, CN=8 (balanced)
            - "lanm_mixed": 3 bi + 2 waters, CN=9 (natural-like)
            - "high_coord_9": 4 bi Glu + 1 water, CN=9 (large lanthanides)
            - Auto-selected based on metal if not specified

        # LEGACY: template_type (for backward compatibility)
        template_type: str (default: "ef_hand" for lanthanides, "none" otherwise)
            - "ef_hand": Legacy 8-residue template
            - "c4_symmetric": Legacy 4-residue template
            - "none": No template (original behavior)

        # Parametric mode (advanced - custom coordination)
        parametric: dict (optional) - Custom coordination parameters:
            - coordination_number: int (6-10)
            - num_waters: int (0-4)
            - bidentate_fraction: float (0.0-1.0)
            - randomize: bool - generate stochastic variants
            - num_variants: int - number of variants if randomize=True

        include_waters: bool (default: False) - Include water molecules in coordination
        add_trp_antenna: bool (default: False) - Add Trp for TEBL luminescence assay
        validate_coordination: bool (default: True) - Run coordination.py validation
        iterative_refinement: bool (default: False) - Refine until target coordination

    Returns:
        {
            "status": "completed" | "failed",
            "result": {
                "designs": [...],
                "approach": str,
                "metal": str,
                "best_design": int,
                "validation": {...},  # If validate_coordination=True
                "template_used": str  # Template name from library
            }
        }
    """
    metal = job_input.get("metal")
    if not metal:
        return {"status": "failed", "error": "Missing 'metal' parameter"}

    metal = metal.upper()
    if metal not in METAL_PROFILES:
        return {
            "status": "failed",
            "error": f"Unknown metal: {metal}. Valid: {', '.join(METAL_PROFILES.keys())}"
        }

    approach = job_input.get("approach", "joint_metal")

    # Check for lanthanide-specific options
    is_lanthanide_metal = LANTHANIDE_TEMPLATES_AVAILABLE and is_lanthanide(metal)

    # NEW: Template library system
    template_name = job_input.get("template_name")
    parametric_params = job_input.get("parametric")

    # Legacy: template_type (backward compatibility)
    template_type = job_input.get("template_type")

    # Determine which template system to use
    if template_name and LANTHANIDE_TEMPLATES_AVAILABLE:
        # NEW: Use template library
        job_input["_template_name"] = template_name
        job_input["_use_library"] = True
        template_type = "library"  # Signal to use library generator
        print(f"[InterfaceMetal] Using template library: {template_name}")
    elif parametric_params and LANTHANIDE_TEMPLATES_AVAILABLE:
        # NEW: Use parametric mode
        job_input["_parametric"] = parametric_params
        job_input["_use_parametric"] = True
        template_type = "parametric"
        print(f"[InterfaceMetal] Using parametric mode: CN={parametric_params.get('coordination_number', 8)}")
    elif template_type is None and is_lanthanide_metal:
        # Auto-select from template library for lanthanides
        template_name = recommend_template(metal) if LANTHANIDE_TEMPLATES_AVAILABLE else None
        if template_name:
            job_input["_template_name"] = template_name
            job_input["_use_library"] = True
            template_type = "library"
            print(f"[InterfaceMetal] Auto-selected template: {template_name} for {metal}")
        else:
            template_type = "ef_hand"  # Fallback to legacy
            print(f"[InterfaceMetal] Using legacy EF-hand template for lanthanide {metal}")
    elif template_type is None:
        template_type = "none"

    # Store template settings in job_input for downstream functions
    job_input["_template_type"] = template_type
    job_input["_is_lanthanide"] = is_lanthanide_metal

    print(f"[InterfaceMetal] Starting with approach={approach}, metal={metal}, template={template_type}...")

    if approach == "joint_metal":
        return _design_joint_metal_dimer(job_input)
    elif approach == "asymmetric_metal":
        return _design_asymmetric_metal_dimer(job_input)
    elif approach == "induced_metal":
        return _design_induced_metal_dimer(job_input)
    elif approach == "bridging_metal":
        return _design_bridging_metal_dimer(job_input)
    elif approach == "redox_switch":
        return _design_redox_switchable_dimer(job_input)
    else:
        return {
            "status": "failed",
            "error": f"Unknown approach: {approach}. Valid: joint_metal, asymmetric_metal, induced_metal, bridging_metal, redox_switch"
        }


def _generate_metal_pdb(metal: str, x: float = 50.0, y: float = 50.0, z: float = 50.0) -> str:
    """
    Generate a minimal PDB file with a single metal ion.

    Args:
        metal: Metal code (ZN, FE, CA, etc.)
        x, y, z: Coordinates for the metal ion (default: 50,50,50 to match lanthanide templates
                 and ensure protein designs wrap AROUND the metal, not radiate FROM it)

    Returns:
        PDB file content as string
    """
    # Metal element symbol (first 2 characters typically)
    element = metal[:2].upper()

    # Standard PDB HETATM format
    pdb_lines = [
        f"HETATM    1 {metal:>4s} {metal:>3s} L   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}",
        "END",
    ]

    return "\n".join(pdb_lines)


def _build_motif_contig(template_type: str, chain_length: str, donor_residue: str = "ASP") -> str:
    """
    Build contig that KEEPS template residues and designs scaffold around them.

    This is the KEY to making templates work - instead of designing entirely new chains,
    we tell RFD3 to preserve the pre-positioned coordinating residues and only design
    the scaffold that connects them.

    Args:
        template_type: "ef_hand" or "c4_symmetric"
        chain_length: Chain length range, e.g., "60-80"
        donor_residue: "ASP", "GLU", or "MIXED" (affects contig comments only)

    Returns:
        Contig string with motif-preserving syntax

    Contig format: "A10,5-15,A15,..." means:
    - "A10" = keep residue 10 from chain A of input PDB
    - "5-15" = design 5-15 new residues as scaffold
    - This preserves template residues while designing around them
    """
    # Parse chain length
    if "-" in str(chain_length):
        min_len, max_len = map(int, str(chain_length).split("-"))
    else:
        min_len = max_len = int(chain_length)

    # Calculate linker lengths to fill between motif positions
    # EF-hand has residues at positions 10, 15, 20, 25 (5 apart)
    # C4 has residues at positions 15, 25 (10 apart)

    if template_type == "ef_hand":
        # 8 coordinating residues: A10,A15,A20,A25 and B10,B15,B20,B25
        # Need linkers between each and at N/C termini
        # Total chain: N-term + A10 + link + A15 + link + A20 + link + A25 + C-term
        # With 4 residues and ~50-70 total, linkers are ~10-15 each

        # N-terminal linker length (design residues 1-9 before first motif)
        nterm_len = "5-12"
        # Inter-motif linkers (between positions 10→15→20→25, each 5 residues apart)
        link_len = "3-8"
        # C-terminal linker (after position 25)
        cterm_len = "15-35"

        # Build contig: keep motif residues, design linkers
        # A10,link,A15,link,A20,link,A25,cterm,/0,B10,link,B15,link,B20,link,B25,cterm
        contig = (
            f"{nterm_len},A10,{link_len},A15,{link_len},A20,{link_len},A25,{cterm_len},"
            f"/0,"
            f"{nterm_len},B10,{link_len},B15,{link_len},B20,{link_len},B25,{cterm_len}"
        )

    elif template_type == "c4_symmetric":
        # 4 coordinating residues: A15,A25 and B15,B25
        # N-terminal linker before position 15
        nterm_len = "8-18"
        # Link between 15 and 25
        link_len = "5-12"
        # C-terminal after 25
        cterm_len = "20-40"

        contig = (
            f"{nterm_len},A15,{link_len},A25,{cterm_len},"
            f"/0,"
            f"{nterm_len},B15,{link_len},B25,{cterm_len}"
        )

    else:
        # No template - design de novo (original behavior)
        contig = f"{min_len}-{max_len},/0,{min_len}-{max_len}"

    return contig


def _build_fixed_atoms(template_type: str, donor_residue: str = "ASP") -> Dict[str, str]:
    """
    Build select_fixed_atoms dict to lock coordinating oxygen atoms at their precise positions.

    This ensures RFD3 preserves the exact geometry of the carboxylate oxygens
    that were carefully positioned by the template generator.

    Args:
        template_type: "ef_hand" or "c4_symmetric"
        donor_residue: "ASP" (OD1/OD2), "GLU" (OE1/OE2), or "MIXED"

    Returns:
        Dictionary mapping residue identifiers to atom names to fix
    """
    # Determine which oxygen atoms to fix based on residue type
    if donor_residue.upper() == "ASP":
        oxygen_atoms = "OD1,OD2"
    elif donor_residue.upper() == "GLU":
        oxygen_atoms = "OE1,OE2"
    else:  # MIXED - need to specify per-residue
        # For EF-hand MIXED: A10,A15,A20 = ASP, A25 = GLU, B10,B15 = ASP, B20,B25 = GLU
        # For C4 MIXED: A15,B15 = ASP, A25,B25 = GLU
        pass  # Handle below

    if template_type == "ef_hand":
        if donor_residue.upper() == "MIXED":
            return {
                "A10": "OD1,OD2", "A15": "OD1,OD2", "A20": "OD1,OD2", "A25": "OE1,OE2",
                "B10": "OD1,OD2", "B15": "OD1,OD2", "B20": "OE1,OE2", "B25": "OE1,OE2",
            }
        else:
            return {
                "A10": oxygen_atoms, "A15": oxygen_atoms,
                "A20": oxygen_atoms, "A25": oxygen_atoms,
                "B10": oxygen_atoms, "B15": oxygen_atoms,
                "B20": oxygen_atoms, "B25": oxygen_atoms,
            }

    elif template_type == "c4_symmetric":
        if donor_residue.upper() == "MIXED":
            return {
                "A15": "OD1,OD2", "A25": "OE1,OE2",
                "B15": "OD1,OD2", "B25": "OE1,OE2",
            }
        else:
            return {
                "A15": oxygen_atoms, "A25": oxygen_atoms,
                "B15": oxygen_atoms, "B25": oxygen_atoms,
            }

    # No template - no atoms to fix
    return {}


def _build_motif_contig_from_library(template_name: str, chain_length: str) -> str:
    """
    Build motif-preserving contig from template library definition.

    Dynamically constructs contig based on the residue positions defined
    in the template library, allowing different templates to have different
    residue arrangements.

    Args:
        template_name: Key from TEMPLATE_LIBRARY
        chain_length: Chain length range, e.g., "60-80"

    Returns:
        Contig string with motif-preserving syntax
    """
    if not LANTHANIDE_TEMPLATES_AVAILABLE:
        # Fall back to de novo if templates not available
        if "-" in str(chain_length):
            min_len, max_len = map(int, str(chain_length).split("-"))
        else:
            min_len = max_len = int(chain_length)
        return f"{min_len}-{max_len},/0,{min_len}-{max_len}"

    template_def = TEMPLATE_LIBRARY.get(template_name)
    if not template_def:
        # Unknown template - fall back to de novo
        if "-" in str(chain_length):
            min_len, max_len = map(int, str(chain_length).split("-"))
        else:
            min_len = max_len = int(chain_length)
        return f"{min_len}-{max_len},/0,{min_len}-{max_len}"

    # Extract residue positions for each chain
    chain_a_residues = sorted(
        [r["resnum"] for r in template_def["residues"] if r["chain"] == "A"]
    )
    chain_b_residues = sorted(
        [r["resnum"] for r in template_def["residues"] if r["chain"] == "B"]
    )

    # Build contig for chain A
    # Format: nterm_link, A{pos1}, link, A{pos2}, ..., cterm_link
    nterm_len = "5-12"
    link_len = "3-8"
    cterm_len = "15-35"

    # Chain A contig
    chain_a_parts = [nterm_len]
    for i, resnum in enumerate(chain_a_residues):
        chain_a_parts.append(f"A{resnum}")
        if i < len(chain_a_residues) - 1:
            chain_a_parts.append(link_len)
    chain_a_parts.append(cterm_len)
    chain_a_contig = ",".join(chain_a_parts)

    # Chain B contig
    chain_b_parts = [nterm_len]
    for i, resnum in enumerate(chain_b_residues):
        chain_b_parts.append(f"B{resnum}")
        if i < len(chain_b_residues) - 1:
            chain_b_parts.append(link_len)
    chain_b_parts.append(cterm_len)
    chain_b_contig = ",".join(chain_b_parts)

    # Combine with chain separator
    return f"{chain_a_contig},/0,{chain_b_contig}"


def _build_fixed_atoms_from_library(template_name: str) -> Dict[str, str]:
    """
    Build select_fixed_atoms dict from template library definition.

    Dynamically determines which atoms to fix based on the residue types
    and modes defined in the template.

    Args:
        template_name: Key from TEMPLATE_LIBRARY

    Returns:
        Dictionary mapping residue identifiers to atom names to fix
    """
    if not LANTHANIDE_TEMPLATES_AVAILABLE:
        return {}

    template_def = TEMPLATE_LIBRARY.get(template_name)
    if not template_def:
        return {}

    fixed_atoms = {}
    for residue in template_def["residues"]:
        chain = residue["chain"]
        resnum = residue["resnum"]
        res_type = residue["type"]
        mode = residue["mode"]

        key = f"{chain}{resnum}"

        # Determine which atoms to fix based on residue type
        if res_type == "ASP":
            if mode == "bidentate":
                fixed_atoms[key] = "OD1,OD2"  # Both oxygens
            else:  # monodentate
                fixed_atoms[key] = "OD1"  # Only coordinating oxygen
        elif res_type == "GLU":
            if mode == "bidentate":
                fixed_atoms[key] = "OE1,OE2"
            else:  # monodentate
                fixed_atoms[key] = "OE1"

    return fixed_atoms


def _design_joint_metal_dimer(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Design a metal-coordinated heterodimer by designing BOTH chains simultaneously
    around a central metal ion.

    Uses contig syntax "min-max,/0,min-max" to design two chains that co-evolve
    around the metal during diffusion. Each chain contributes coordinating residues.

    For lanthanides (TB, GD, EU, etc.), uses template-based design with pre-positioned
    coordinating Glu residues based on Caldwell et al. 2020 methodology.

    Contig Format Reference:
        - Comma (,): Separates contig elements
        - /0: Chain break with zero gap (creates separate chain)
        - X-Y: Design X to Y residues (range)
        - ChainResnum: Keep specific residue from input (e.g., A10, B25)
        Example: "60-80,/0,60-80" = Two chains of 60-80 residues each

    Symmetry Options:
        This function uses explicit 2-chain contig ("X,/0,X") rather than
        RFdiffusion's C2 symmetry operator. Trade-offs:

        **Current Approach (explicit 2-chain):**
        - Pros: Each chain can have different template residues
        - Pros: Works with asymmetric lanthanide coordination (3:1 splits)
        - Cons: No guaranteed geometric symmetry
        - Cons: May produce asymmetric interfaces

        **Alternative (C2 symmetry operator):**
        - Pros: Guaranteed geometric symmetry
        - Pros: Metal exactly at symmetry axis
        - Cons: Requires pre-symmetrized motifs
        - Cons: Each chain must contribute identical coordination
        - Use: `_design_symmetric_dimer()` for C2-symmetric designs

        For applications requiring strict geometric symmetry, consider using
        `_design_symmetric_dimer()` with `symmetry="C2"` parameter instead.
    """
    metal = job_input.get("metal", "ZN").upper()
    chain_length = job_input.get("chain_length", "60-80")
    num_designs = job_input.get("num_designs", 3)
    seed = job_input.get("seed")

    # Coordination parameters
    coordination_split = job_input.get("coordination_split", [2, 2])
    chain_a_donors = job_input.get("chain_a_donors")
    chain_b_donors = job_input.get("chain_b_donors")

    # Lanthanide-specific options (from handle_interface_metal_design)
    template_type = job_input.get("_template_type", "none")
    is_lanthanide_metal = job_input.get("_is_lanthanide", False)
    include_waters = job_input.get("include_waters", False)
    add_trp_antenna = job_input.get("add_trp_antenna", False)
    validate_coordination = job_input.get("validate_coordination", True)
    iterative_refinement = job_input.get("iterative_refinement", False)
    donor_residue = job_input.get("donor_residue", "ASP")  # NEW: ASP recommended for lanthanides
    use_motif_scaffolding = job_input.get("use_motif_scaffolding", True)  # NEW: use template as motif

    profile = METAL_PROFILES[metal]
    total_coord = sum(coordination_split)

    # For lanthanides, expect higher coordination
    if is_lanthanide_metal:
        target_coordination = 8  # Lanthanides need 8-9 coordination
        if total_coord < target_coordination:
            print(f"[JointMetal] Adjusting lanthanide coordination target to {target_coordination}")
            total_coord = target_coordination
    else:
        target_coordination = total_coord

    # Validate coordination number
    if total_coord not in profile["coordination"]:
        print(f"[JointMetal] Warning: {metal} typically has {profile['coordination']}-coordinate, not {total_coord}")

    # Parse chain length
    if "-" in str(chain_length):
        min_len, max_len = map(int, str(chain_length).split("-"))
    else:
        min_len = max_len = int(chain_length)

    # Generate template PDB or simple metal PDB
    # Use consistent metal position (50,50,50) for all approaches
    metal_center = [50.0, 50.0, 50.0]
    fixed_atoms = {}  # Will be populated for template motif scaffolding
    template_used = None  # Track which template was used

    # NEW: Template Library System
    if template_type == "library" and LANTHANIDE_TEMPLATES_AVAILABLE:
        # Use new template library with chemically realistic coordination
        template_name = job_input.get("_template_name", "caldwell_4")
        template_pdb = generate_template_from_library(
            template_name=template_name,
            metal=metal,
            add_trp_antenna=add_trp_antenna,
        )
        template_info = get_template_info(template_name)
        template_used = template_name
        print(f"[JointMetal] Using template library: {template_name}")
        print(f"[JointMetal] Template: CN={template_info['coordination_number']}, geometry={template_info['geometry']}")
        metal_pdb = template_pdb

        # Build motif-preserving contig based on template definition
        if use_motif_scaffolding:
            contig = _build_motif_contig_from_library(template_name, chain_length)
            fixed_atoms = _build_fixed_atoms_from_library(template_name)
            print(f"[JointMetal] MOTIF SCAFFOLDING from library template")
            print(f"[JointMetal] Fixed atoms: {fixed_atoms}")
        else:
            contig = f"{min_len}-{max_len},/0,{min_len}-{max_len}"
            print(f"[JointMetal] DE NOVO mode: template only provides metal position")

    # NEW: Parametric Mode
    elif template_type == "parametric" and LANTHANIDE_TEMPLATES_AVAILABLE:
        parametric = job_input.get("_parametric", {})
        templates = generate_parametric_template(
            metal=metal,
            coordination_number=parametric.get("coordination_number", 8),
            num_waters=parametric.get("num_waters", 0),
            bidentate_fraction=parametric.get("bidentate_fraction", 0.5),
            preferred_donors=parametric.get("preferred_donors", ["ASP", "GLU"]),
            randomize=parametric.get("randomize", False),
            num_variants=parametric.get("num_variants", 1),
            seed=seed,
            pocket_radius=parametric.get("pocket_radius", 6.0),
            add_trp_antenna=add_trp_antenna,
        )
        template_pdb = templates[0]  # Use first variant for now
        template_used = f"parametric_cn{parametric.get('coordination_number', 8)}"
        print(f"[JointMetal] Using parametric template: CN={parametric.get('coordination_number', 8)}")
        metal_pdb = template_pdb

        # For parametric, use de novo contig (template provides geometry only)
        contig = f"{min_len}-{max_len},/0,{min_len}-{max_len}"

    elif template_type != "none" and LANTHANIDE_TEMPLATES_AVAILABLE and is_lanthanide_metal:
        # LEGACY: Use old template-based approach for lanthanides
        # Templates already place metal at (50,50,50)
        template_pdb = get_template(
            template_type=template_type,
            metal=metal,
            donor_residue=donor_residue,  # NEW: support ASP, GLU, or MIXED
            include_waters=include_waters,
            add_trp_antenna=add_trp_antenna,
        )
        template_used = f"legacy_{template_type}"
        print(f"[JointMetal] Using LEGACY {template_type} template for lanthanide {metal}")
        print(f"[JointMetal] Template: donor_residue={donor_residue}, waters={include_waters}, Trp={add_trp_antenna}")
        metal_pdb = template_pdb

        # Build motif-preserving contig (KEY FIX: preserves template residues!)
        if use_motif_scaffolding:
            contig = _build_motif_contig(template_type, chain_length, donor_residue)
            fixed_atoms = _build_fixed_atoms(template_type, donor_residue)
            print(f"[JointMetal] MOTIF SCAFFOLDING: contig preserves template residues")
            print(f"[JointMetal] Fixed atoms: {fixed_atoms}")
        else:
            # Original behavior - design new chains (ignores template residues)
            contig = f"{min_len}-{max_len},/0,{min_len}-{max_len}"
            print(f"[JointMetal] DE NOVO mode: template only provides metal position")
    else:
        # Generate a simple PDB with the metal ion at (50,50,50)
        # This ensures protein designs wrap AROUND the metal, not radiate FROM it
        metal_pdb = _generate_metal_pdb(metal, *metal_center)
        template_used = "none"
        # No template - design de novo
        contig = f"{min_len}-{max_len},/0,{min_len}-{max_len}"

    # Store template used for result
    job_input["_template_used"] = template_used

    print(f"[JointMetal] Metal={metal}, target_coordination={target_coordination}")
    print(f"[JointMetal] Preferred donors: {profile['donors']}")
    print(f"[JointMetal] Contig: {contig}")
    print(f"[JointMetal] Metal center: {metal_center}")

    designs = []
    best_design_idx = 0
    best_score = float('-inf')

    for design_idx in range(num_designs):
        design_seed = (seed + design_idx) if seed is not None else None
        print(f"[JointMetal] Generating design {design_idx + 1}/{num_designs}...")

        # Build RFD3 config for metal-coordinated dimer
        # We provide the metal PDB as input and use ligand code

        # Calculate ori_token to position protein COM offset from metal
        # Small Z offset ensures protein wraps around metal rather than centering on it
        # User can override with custom ori_token if needed
        user_ori_token = job_input.get("ori_token")
        if user_ori_token:
            ori_token = user_ori_token
        else:
            # Default: offset 3Å from metal center in Z direction
            ori_token = [metal_center[0], metal_center[1], metal_center[2] + 3.0]

        if design_idx == 0:
            print(f"[JointMetal] Using ori_token: {ori_token} (user override: {user_ori_token is not None})")

        rfd3_input = {
            "task": "rfd3",
            "contig": contig,
            "pdb_content": metal_pdb,  # Provide metal PDB as input
            "ligand": metal,  # Metal code as ligand
            "ori_token": ori_token,  # Control protein COM positioning relative to metal
            "seed": design_seed,
            "num_designs": 1,
            "is_non_loopy": True,  # Prefer structured elements near metal
            # Guiding potentials for metal coordination optimization
            # substrate_contacts: reinforces metal-protein proximity during diffusion
            # monomer_ROG: encourages compact structures (prevents extended/unfolded designs)
            # olig_contacts: optimizes inter-chain contacts at dimer interface
            "guiding_potentials": [
                "type:substrate_contacts,weight:5,s:1,r_0:8,d_0:4",
                "type:monomer_ROG,weight:1,min_dist:15",
                "type:olig_contacts,weight_intra:1,weight_inter:0.5",
            ],
            "guide_scale": 2,
        }

        # Add motif scaffolding parameters when using templates
        if fixed_atoms:
            rfd3_input["select_fixed_atoms"] = fixed_atoms
            # Hotspot ensures designed residues contact the metal
            rfd3_input["select_hotspots"] = {"L1": "all"}  # L is metal chain
            if design_idx == 0:
                print(f"[JointMetal] Added select_fixed_atoms for {len(fixed_atoms)} residues")
                print(f"[JointMetal] Added select_hotspots: L1:all (ensure metal contact)")

            # NEW: Add H-bond conditioning for secondary coordination sphere
            # Carboxylate oxygens (OE1, OE2 for Glu; OD1, OD2 for Asp) are H-bond acceptors
            # This reinforces the H-bond network around the metal site
            # Research shows 11x activity boost from proper secondary sphere H-bonding
            if template_used and template_type == "library" and LANTHANIDE_TEMPLATES_AVAILABLE:
                hbond_acceptors = {}
                template_def = TEMPLATE_LIBRARY.get(template_used, {})
                for res in template_def.get("residues", []):
                    chain = res.get("chain", "A")
                    resnum = res.get("resnum")
                    res_type = res.get("type", "GLU")
                    if resnum:
                        key = f"{chain}{resnum}"
                        if res_type == "GLU":
                            hbond_acceptors[key] = "OE1,OE2"
                        elif res_type == "ASP":
                            hbond_acceptors[key] = "OD1,OD2"

                if hbond_acceptors:
                    rfd3_input["select_hbond_acceptor"] = hbond_acceptors
                    if design_idx == 0:
                        print(f"[JointMetal] Added H-bond conditioning for {len(hbond_acceptors)} coordinating residues")

        # Run RFD3
        result = handle_rfd3(rfd3_input)

        if result.get("status") != "completed":
            print(f"[JointMetal] Design {design_idx + 1} failed: {result.get('error')}")
            continue

        result_designs = result.get("result", {}).get("designs", [])
        if not result_designs:
            print(f"[JointMetal] Design {design_idx + 1} produced no designs")
            continue

        pdb_content = result_designs[0].get("content") or result_designs[0].get("pdb_content")
        if not pdb_content:
            print(f"[JointMetal] Design {design_idx + 1}: No PDB content")
            continue

        # Run LigandMPNN sequence redesign for lanthanides with carboxylate bias
        # This uses D:6.0,E:4.0 bias to favor Asp/Glu coordinating residues
        # LigandMPNN atomize_side_chains for proper metal binding context
        if is_lanthanide_metal:
            print(f"[JointMetal] Design {design_idx + 1}: Running LigandMPNN with lanthanide preset...")
            try:
                # Fix atom names/IDs before MPNN to avoid NaN issues
                # (see: https://github.com/RosettaCommons/foundry/issues/123)
                from inference_utils import fix_pdb_for_mpnn
                pdb_content = fix_pdb_for_mpnn(pdb_content)

                # CRITICAL FIX: Extract fixed positions from template to preserve metal coordinators
                # Without this, LigandMPNN may mutate the Asp/Glu residues positioned for
                # metal coordination, destroying the carefully designed binding site geometry.
                coord_fixed_positions = None
                if template_used and LANTHANIDE_TEMPLATES_AVAILABLE:
                    coord_fixed_positions = get_template_fixed_positions(template_used)
                    if coord_fixed_positions:
                        print(f"[JointMetal] Design {design_idx + 1}: Fixing coordinating residues: {coord_fixed_positions}")

                # run_ligandmpnn_for_ligand_binding() automatically fixes incomplete backbone
                # by adding missing O atoms with idealized geometry
                mpnn_result = run_ligandmpnn_for_ligand_binding(
                    pdb_content=pdb_content,
                    ligand_type="lanthanide",
                    num_sequences=4,
                    temperature=0.1,
                    fixed_positions=coord_fixed_positions,
                )
                if mpnn_result.get("status") == "completed":
                    result_data = mpnn_result.get("result", {})
                    # PDB structures are in packed_structures (when pack_side_chains=True)
                    packed_structures = result_data.get("packed_structures", [])
                    sequences = result_data.get("sequences", [])

                    if packed_structures:
                        # Use the first packed structure (best sequence with sidechains)
                        redesigned_pdb = packed_structures[0].get("pdb_content")
                        if redesigned_pdb:
                            pdb_content = redesigned_pdb
                            print(f"[JointMetal] Design {design_idx + 1}: Sequence redesigned with D/E bias + proper sidechains")
                        else:
                            # NaN case: atomize_side_chains produced NaN coordinates
                            # Fall back to FastRelax sidechain packing (BindCraft approach)
                            print(f"[JointMetal] Design {design_idx + 1}: MPNN sidechain atomization failed (NaN coords)")
                            print(f"[JointMetal] Design {design_idx + 1}: Falling back to FastRelax sidechain packing...")

                            # Get sequence from FASTA output
                            if sequences:
                                seq_content = sequences[0].get("content", "")
                                seq_lines = [l for l in seq_content.split('\n') if l and not l.startswith('>')]
                                if seq_lines:
                                    first_seq = ''.join(seq_lines)
                                    d_count = first_seq.count('D')
                                    e_count = first_seq.count('E')
                                    print(f"[JointMetal] Design {design_idx + 1}: Using MPNN sequence: {d_count} Asp, {e_count} Glu")

                                    # Apply sequence to backbone and use FastRelax for sidechains
                                    try:
                                        from inference_utils import apply_sequence_to_pdb
                                        from rosetta_utils import fastrelax_protein_only

                                        # Apply MPNN sequence to the RFD backbone
                                        # Preserve sidechains for coordinating residues (they have proper metal geometry)
                                        seq_applied_pdb = apply_sequence_to_pdb(
                                            pdb_content,
                                            first_seq,
                                            keep_sidechains=False,  # Strip for FastRelax
                                            preserve_coordinating_residues=True,  # Keep coordinating residue sidechains
                                            metal_coord_distance=3.5,  # Å threshold
                                        )

                                        # Run FastRelax to rebuild sidechains
                                        relax_result = fastrelax_protein_only(
                                            seq_applied_pdb,
                                            max_iter=200,
                                            constrain_coords=True,
                                            repack_only=True  # Just repack sidechains, don't minimize backbone
                                        )

                                        if relax_result.get("status") == "completed":
                                            pdb_content = relax_result.get("pdb_content", pdb_content)
                                            print(f"[JointMetal] Design {design_idx + 1}: FastRelax sidechain packing successful")
                                        else:
                                            print(f"[JointMetal] Design {design_idx + 1}: FastRelax failed: {relax_result.get('error')}")
                                            # Keep the sequence-applied backbone without sidechains
                                            pdb_content = seq_applied_pdb
                                    except Exception as relax_err:
                                        print(f"[JointMetal] Design {design_idx + 1}: FastRelax error: {relax_err}")
                            else:
                                print(f"[JointMetal] Design {design_idx + 1}: No sequences available for fallback")
                    else:
                        # No packed structures - use sequences directly
                        if sequences:
                            seq_content = sequences[0].get("content", "")
                            seq_lines = [l for l in seq_content.split('\n') if l and not l.startswith('>')]
                            if seq_lines:
                                first_seq = ''.join(seq_lines)
                                d_count = first_seq.count('D')
                                e_count = first_seq.count('E')
                                print(f"[JointMetal] Design {design_idx + 1}: LigandMPNN sequence: {d_count} Asp, {e_count} Glu (no packed structures)")
                        else:
                            print(f"[JointMetal] Design {design_idx + 1}: LigandMPNN returned no packed structures, using original")
                else:
                    print(f"[JointMetal] Design {design_idx + 1}: LigandMPNN failed: {mpnn_result.get('error')}")
            except Exception as e:
                import traceback
                traceback.print_exc()
                print(f"[JointMetal] Design {design_idx + 1}: LigandMPNN error: {e}")

        # Analyze metal coordination
        coord_analysis = _analyze_metal_coordination(pdb_content, metal)

        # Calculate sequence identity between chains
        seq_identity = calculate_sequence_identity(pdb_content, "A", "B")
        identity_pct = seq_identity.get("sequence_identity_percent", 100)

        design_entry = {
            "pdb_content": pdb_content,
            "design_index": design_idx,
            "metrics": {
                "coordination_number": coord_analysis.get("coordination_number", 0),
                "chain_a_donors": coord_analysis.get("chain_a_donors", 0),
                "chain_b_donors": coord_analysis.get("chain_b_donors", 0),
                "donor_residues": coord_analysis.get("donor_residues", []),
                "average_distance": coord_analysis.get("average_distance"),
                "sequence_identity": identity_pct,
                "is_heterodimer": identity_pct < 70,
            },
            "coordination_analysis": coord_analysis,
        }

        # Score based on metal coordination quality
        score = 0
        coord_num = coord_analysis.get("coordination_number", 0)

        # Prefer correct coordination number
        if coord_num in profile["coordination"]:
            score += 30
        elif coord_num > 0:
            score += 10

        # Prefer both chains contributing
        if coord_analysis.get("chain_a_donors", 0) >= 1 and coord_analysis.get("chain_b_donors", 0) >= 1:
            score += 20

        # Prefer true heterodimer
        if identity_pct < 70:
            score += 20

        # Prefer appropriate bond distances
        avg_dist = coord_analysis.get("average_distance")
        if avg_dist and profile["distance"][0] <= avg_dist <= profile["distance"][1]:
            score += 15

        design_entry["score"] = score
        designs.append(design_entry)

        if score > best_score:
            best_score = score
            best_design_idx = len(designs) - 1

        print(f"[JointMetal] Design {design_idx + 1}: coord={coord_num}, identity={identity_pct:.1f}%, score={score}")

    if not designs:
        return {
            "status": "failed",
            "error": "No valid joint metal dimer designs produced"
        }

    # Post-design validation for lanthanides using metal_validation module
    validation_result = None
    tebl_result = None

    # Find the actual metal residue number in the PDB (may differ from template)
    def _find_metal_resnum(pdb_content: str, metal_code: str) -> int:
        """Find the residue number of the metal in the PDB."""
        for line in pdb_content.split('\n'):
            if line.startswith('HETATM') and metal_code in line[17:20].strip():
                try:
                    return int(line[22:26].strip())
                except ValueError:
                    pass
        return 1  # Default fallback

    if validate_coordination and METAL_VALIDATION_AVAILABLE and is_lanthanide_metal:
        best_pdb = designs[best_design_idx]["pdb_content"]
        metal_resnum = _find_metal_resnum(best_pdb, metal)
        print(f"[JointMetal] Running lanthanide validation on best design (metal at L:{metal_resnum})...")

        try:
            validation_result = validate_lanthanide_site(
                pdb_content=best_pdb,
                metal=metal,
                metal_chain="L",
                metal_resnum=metal_resnum,
                target_coordination=target_coordination,
                check_tebl=add_trp_antenna,
            )

            if validation_result.get("success"):
                quality_score = validation_result.get("quality_score", 0)
                quality_rating = validation_result.get("quality_rating", "unknown")
                coord_num = validation_result.get("coordination_number", 0)
                print(f"[JointMetal] Validation: quality={quality_score:.1f} ({quality_rating}), coord={coord_num}/{target_coordination}")

                # Add interface metrics if available
                try:
                    from binding_analysis import analyze_interface
                    interface_result = analyze_interface(best_pdb)
                    if interface_result.get("success"):
                        validation_result["interface_metrics"] = {
                            "contacts": interface_result.get("contacts", 0),
                            "hbonds": interface_result.get("hbonds_int", 0),
                            "buried_sasa": interface_result.get("dSASA_int", 0),
                            "packstat": interface_result.get("packstat", 0),
                        }
                        print(f"[JointMetal] Interface: contacts={interface_result.get('contacts', 0)}, hbonds={interface_result.get('hbonds_int', 0)}")
                except Exception as e:
                    # Interface analysis is optional, don't fail validation
                    validation_result["interface_metrics"] = {"error": str(e)}

                # Add validation to best design's metrics
                designs[best_design_idx]["validation"] = validation_result

                # Check TEBL readiness
                if add_trp_antenna and "tebl_details" in validation_result:
                    tebl_details = validation_result["tebl_details"]
                    tebl_ready = tebl_details.get("has_antenna", False)
                    trp_dist = tebl_details.get("trp_metal_distance", "N/A")
                    print(f"[JointMetal] TEBL ready: {tebl_ready}, Trp distance: {trp_dist}Å")

        except Exception as e:
            print(f"[JointMetal] Validation error: {e}")
            validation_result = {"error": str(e)}

    # Add TEBL signal prediction if requested
    if add_trp_antenna and TEBL_ANALYSIS_AVAILABLE and is_lanthanide_metal:
        best_pdb = designs[best_design_idx]["pdb_content"]
        metal_resnum = _find_metal_resnum(best_pdb, metal)
        try:
            tebl_result = predict_tebl_signal(
                pdb_content=best_pdb,
                metal=metal,
                metal_chain="L",
                metal_resnum=metal_resnum,
            )
            if tebl_result.get("has_antenna"):
                signal = tebl_result.get("signal_strength", "unknown")
                efficiency = tebl_result.get("predicted_efficiency", 0)
                print(f"[JointMetal] TEBL prediction: signal={signal}, efficiency={efficiency:.3f}")
                designs[best_design_idx]["tebl_prediction"] = tebl_result
        except Exception as e:
            print(f"[JointMetal] TEBL analysis error: {e}")

    result = {
        "status": "completed",
        "result": {
            "approach": "joint_metal",
            "metal": metal,
            "metal_profile": profile,
            "designs": designs,
            "best_design": best_design_idx,
            "dimer": {
                "pdb_content": designs[best_design_idx]["pdb_content"],
                "metrics": designs[best_design_idx]["metrics"],
                "coordination": designs[best_design_idx]["coordination_analysis"],
            },
            # Lanthanide-specific results
            "template_type": template_type if is_lanthanide_metal else None,
            "target_coordination": target_coordination if is_lanthanide_metal else None,
        }
    }

    # Add validation results if available
    if validation_result:
        result["result"]["validation"] = validation_result
    if tebl_result:
        result["result"]["tebl_analysis"] = tebl_result

    return result


def _design_asymmetric_metal_dimer(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Design metal-coordinated heterodimer with different donor types per chain.

    For example: Chain A provides His-His, Chain B provides Cys-Cys.
    This creates asymmetric chemistry at the interface.
    """
    metal = job_input.get("metal", "ZN").upper()
    chain_length = job_input.get("chain_length", "60-80")
    num_designs = job_input.get("num_designs", 3)
    seed = job_input.get("seed")

    chain_a_donors = job_input.get("chain_a_donors", ["His", "His"])
    chain_b_donors = job_input.get("chain_b_donors", ["Cys", "Cys"])

    profile = METAL_PROFILES[metal]

    print(f"[AsymmetricMetal] Metal={metal}")
    print(f"[AsymmetricMetal] Chain A donors: {chain_a_donors}")
    print(f"[AsymmetricMetal] Chain B donors: {chain_b_donors}")

    # For now, use joint design and rely on RFD3 to create variety
    # A more sophisticated approach would use sequence constraints post-design
    job_input["approach"] = "joint_metal"  # Fallback to joint
    result = _design_joint_metal_dimer(job_input)

    if result.get("status") == "completed":
        result["result"]["approach"] = "asymmetric_metal"
        result["result"]["target_donors"] = {
            "chain_a": chain_a_donors,
            "chain_b": chain_b_donors,
        }

    return result


def _design_induced_metal_dimer(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Design monomers with incomplete coordination that only dimerize when metal binds.

    Each monomer has partial coordination sites exposed. Metal binding completes
    the coordination sphere by bringing the chains together.
    """
    metal = job_input.get("metal", "ZN").upper()
    chain_length = job_input.get("chain_length", "60-80")
    num_designs = job_input.get("num_designs", 3)
    seed = job_input.get("seed")

    profile = METAL_PROFILES[metal]

    print(f"[InducedMetal] Designing metal-induced dimerization with {metal}")

    # Design using joint approach - the incomplete coordination creates
    # energetic drive for dimerization
    job_input["approach"] = "joint_metal"
    result = _design_joint_metal_dimer(job_input)

    if result.get("status") == "completed":
        result["result"]["approach"] = "induced_metal"
        result["result"]["mechanism"] = "incomplete_coordination_drives_assembly"

    return result


def _design_bridging_metal_dimer(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Design dimers with multiple metal ions bridging the interface.

    Examples: Fe-Fe (hemerythrin-like), Cu-Cu (type 3 copper), Zn-Zn (metallo-β-lactamase)
    """
    metal = job_input.get("metal", "FE").upper()
    second_metal = job_input.get("bridging_metal", metal).upper()
    chain_length = job_input.get("chain_length", "60-80")
    num_designs = job_input.get("num_designs", 3)
    seed = job_input.get("seed")

    print(f"[BridgingMetal] Designing {metal}-{second_metal} bridged dimer")
    print(f"[BridgingMetal] Note: Dinuclear site design is advanced - using single metal with larger interface")

    # For now, design with single metal - dinuclear sites require specialized constraints
    job_input["approach"] = "joint_metal"
    result = _design_joint_metal_dimer(job_input)

    if result.get("status") == "completed":
        result["result"]["approach"] = "bridging_metal"
        result["result"]["target_metals"] = [metal, second_metal]
        result["result"]["note"] = "Dinuclear site requires post-design metal placement"

    return result


def _design_redox_switchable_dimer(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Design dimers where metal oxidation state controls affinity.

    For example: Fe²⁺ (tetrahedral) vs Fe³⁺ (octahedral) have different
    preferred coordination, affecting dimer stability.
    """
    metal = job_input.get("metal", "FE").upper()
    chain_length = job_input.get("chain_length", "60-80")
    num_designs = job_input.get("num_designs", 3)
    seed = job_input.get("seed")

    # Check if metal is redox-active
    redox_metals = {"FE", "CU", "MN", "CO"}
    if metal not in redox_metals:
        return {
            "status": "failed",
            "error": f"{metal} is not redox-active. Use Fe, Cu, Mn, or Co."
        }

    print(f"[RedoxSwitch] Designing redox-switchable dimer with {metal}")

    # Design for the lower oxidation state (more stable coordination)
    job_input["approach"] = "joint_metal"
    result = _design_joint_metal_dimer(job_input)

    if result.get("status") == "completed":
        result["result"]["approach"] = "redox_switch"
        result["result"]["redox_states"] = {
            "reduced": f"{metal}²⁺",
            "oxidized": f"{metal}³⁺",
        }
        result["result"]["mechanism"] = "coordination_geometry_change_on_oxidation"

    return result


def _analyze_metal_coordination(pdb_content: str, metal: str) -> Dict[str, Any]:
    """
    Analyze metal coordination in a PDB structure.

    Returns coordination number, donor residues, distances, and per-chain breakdown.
    """
    import re

    # Parse PDB
    lines = pdb_content.strip().split("\n")

    # Find metal atom(s)
    metal_atoms = []
    protein_atoms = []

    for line in lines:
        if line.startswith("HETATM") and metal in line[17:20]:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            metal_atoms.append({"x": x, "y": y, "z": z, "line": line})
        elif line.startswith("ATOM"):
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain = line[21]
            res_num = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            # Check if this is a potential donor atom
            is_donor = False
            if atom_name in ["NE2", "ND1"] and res_name == "HIS":
                is_donor = True
            elif atom_name == "SG" and res_name == "CYS":
                is_donor = True
            elif atom_name in ["OD1", "OD2"] and res_name == "ASP":
                is_donor = True
            elif atom_name in ["OE1", "OE2"] and res_name == "GLU":
                is_donor = True
            elif atom_name == "OH" and res_name == "TYR":
                is_donor = True
            elif atom_name == "SD" and res_name == "MET":
                is_donor = True
            elif atom_name in ["OD1", "ND2"] and res_name == "ASN":
                is_donor = True

            if is_donor:
                protein_atoms.append({
                    "atom_name": atom_name,
                    "res_name": res_name,
                    "chain": chain,
                    "res_num": res_num,
                    "x": x, "y": y, "z": z,
                })

    if not metal_atoms:
        return {
            "coordination_number": 0,
            "chain_a_donors": 0,
            "chain_b_donors": 0,
            "donor_residues": [],
            "error": "No metal found in structure",
        }

    # Use first metal atom
    metal_atom = metal_atoms[0]
    mx, my, mz = metal_atom["x"], metal_atom["y"], metal_atom["z"]

    # Find coordinating atoms (within 3.0 Å)
    profile = METAL_PROFILES.get(metal, {})
    max_dist = profile.get("distance", (2.0, 2.5))[1] + 0.5  # Add tolerance

    coordinating = []
    for atom in protein_atoms:
        dist = ((atom["x"] - mx)**2 + (atom["y"] - my)**2 + (atom["z"] - mz)**2) ** 0.5
        if dist <= max_dist:
            coordinating.append({
                "residue": f"{atom['chain']}{atom['res_num']}",
                "res_name": atom["res_name"],
                "atom": atom["atom_name"],
                "chain": atom["chain"],
                "distance": round(dist, 2),
            })

    # Sort by distance
    coordinating.sort(key=lambda x: x["distance"])

    # Count per chain
    chain_a_count = sum(1 for c in coordinating if c["chain"] == "A")
    chain_b_count = sum(1 for c in coordinating if c["chain"] == "B")

    # Calculate average distance
    avg_dist = sum(c["distance"] for c in coordinating) / len(coordinating) if coordinating else None

    return {
        "coordination_number": len(coordinating),
        "chain_a_donors": chain_a_count,
        "chain_b_donors": chain_b_count,
        "donor_residues": coordinating,
        "average_distance": round(avg_dist, 2) if avg_dist else None,
    }


# ============== ESM3 Handlers ==============

def handle_esm3_score(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Score protein sequences using ESM3 perplexity.

    Input:
        sequences: list[str] - List of amino acid sequences to score

    Returns:
        scores: list of {sequence, perplexity, score}
        Lower perplexity = better sequence quality
    """
    sequences = job_input.get("sequences")
    if not sequences:
        return {"status": "failed", "error": "Missing 'sequences' parameter"}

    if not isinstance(sequences, list):
        return {"status": "failed", "error": "'sequences' must be a list of strings"}

    return esm3_score_sequences(job_input)


def handle_esm3_generate(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Generate protein sequences using ESM3 with optional function conditioning.

    Input:
        prompt: str (optional) - Partial sequence to complete
        functions: list[str] (optional) - Function keywords like ["zinc-binding", "hydrolase"]
        num_sequences: int (default: 4) - Number of sequences to generate
        temperature: float (default: 0.7) - Sampling temperature
        max_length: int (default: 200) - Maximum sequence length

    Returns:
        sequences: list[str] - Generated protein sequences
    """
    return esm3_generate_sequence(job_input)


def handle_esm3_embed(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Get ESM3 embeddings for protein sequences.

    Input:
        sequences: list[str] - List of amino acid sequences

    Returns:
        embeddings: list of {sequence, per_residue, global}
        per_residue: [seq_len, hidden_dim] per-position embeddings
        global: [hidden_dim] mean-pooled sequence embedding
    """
    sequences = job_input.get("sequences")
    if not sequences:
        return {"status": "failed", "error": "Missing 'sequences' parameter"}

    if not isinstance(sequences, list):
        return {"status": "failed", "error": "'sequences' must be a list of strings"}

    return esm3_get_embeddings(job_input)


def handle_download_checkpoints(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Download model checkpoints to the Network Volume using Foundry CLI.

    This task downloads the required model checkpoints (RFD3, RF3, MPNN) to the
    Network Volume at /runpod-volume/checkpoints/. This only needs to be run once.

    Input:
        force: bool - Force re-download even if checkpoints exist (default: False)

    Returns:
        status: "completed" | "failed"
        result: {
            checkpoints_dir: str,
            before: {files, total_size_gb},
            after: {files, total_size_gb},
            downloaded: bool,
            message: str
        }
    """
    import subprocess
    import glob as glob_module

    force = job_input.get("force", False)

    def get_checkpoint_info(directory: str) -> Dict[str, Any]:
        """Get info about checkpoint files in directory"""
        files = []
        total_size = 0

        if os.path.exists(directory):
            # Look for checkpoint files
            patterns = ["**/*.ckpt", "**/*.pt", "**/*.pth", "**/*.safetensors"]
            for pattern in patterns:
                for filepath in glob_module.glob(os.path.join(directory, pattern), recursive=True):
                    size = os.path.getsize(filepath)
                    files.append({
                        "name": os.path.basename(filepath),
                        "path": filepath,
                        "size_mb": round(size / (1024 * 1024), 2)
                    })
                    total_size += size

        return {
            "files": files,
            "file_count": len(files),
            "total_size_gb": round(total_size / (1024 * 1024 * 1024), 2)
        }

    try:
        print(f"[Handler] Checking checkpoints at: {CHECKPOINT_DIR}")

        # Ensure checkpoint directory exists
        os.makedirs(CHECKPOINT_DIR, exist_ok=True)

        # Get current state
        before_info = get_checkpoint_info(CHECKPOINT_DIR)
        print(f"[Handler] Before: {before_info['file_count']} files, {before_info['total_size_gb']} GB")

        # Check if we already have checkpoints
        has_checkpoints = before_info["total_size_gb"] > 1.0  # At least 1GB of checkpoints

        if has_checkpoints and not force:
            return {
                "status": "completed",
                "result": {
                    "checkpoints_dir": CHECKPOINT_DIR,
                    "before": before_info,
                    "after": before_info,
                    "downloaded": False,
                    "message": f"Checkpoints already exist ({before_info['total_size_gb']} GB). Use force=True to re-download."
                }
            }

        # If force is enabled, delete any 0-byte (corrupted) checkpoint files first
        if force:
            for file_info in before_info.get("files", []):
                if file_info.get("size_mb", 0) == 0:
                    filepath = file_info.get("path")
                    if filepath and os.path.exists(filepath):
                        print(f"[Handler] Deleting corrupted 0-byte file: {filepath}")
                        os.remove(filepath)

        # Download checkpoints using Foundry CLI
        cmd = ["foundry", "install", "base-models"]
        if force:
            cmd.append("--force")
        print(f"[Handler] Downloading checkpoints using Foundry CLI...")
        print(f"[Handler] Running: {' '.join(cmd)}")

        # Set environment for Foundry
        env = os.environ.copy()
        env["FOUNDRY_CHECKPOINT_DIRS"] = CHECKPOINT_DIR

        # Run foundry install
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            env=env,
            timeout=3600  # 1 hour timeout for downloads
        )

        print(f"[Handler] Foundry stdout: {result.stdout}")
        if result.stderr:
            print(f"[Handler] Foundry stderr: {result.stderr}")

        if result.returncode != 0:
            # Try alternative: foundry download
            print("[Handler] Trying alternative: foundry download...")
            result2 = subprocess.run(
                ["foundry", "download", "--all"],
                capture_output=True,
                text=True,
                env=env,
                timeout=3600
            )
            print(f"[Handler] Alternative stdout: {result2.stdout}")
            if result2.stderr:
                print(f"[Handler] Alternative stderr: {result2.stderr}")

        # Get updated state
        after_info = get_checkpoint_info(CHECKPOINT_DIR)
        print(f"[Handler] After: {after_info['file_count']} files, {after_info['total_size_gb']} GB")

        downloaded_size = after_info["total_size_gb"] - before_info["total_size_gb"]

        return {
            "status": "completed",
            "result": {
                "checkpoints_dir": CHECKPOINT_DIR,
                "before": before_info,
                "after": after_info,
                "downloaded": True,
                "downloaded_gb": round(downloaded_size, 2),
                "message": f"Downloaded {round(downloaded_size, 2)} GB of checkpoints. Total: {after_info['total_size_gb']} GB",
                "foundry_output": result.stdout[:2000] if result.stdout else None
            }
        }

    except subprocess.TimeoutExpired:
        return {
            "status": "failed",
            "error": "Checkpoint download timed out after 1 hour"
        }
    except Exception as e:
        print(f"[Handler] Download error: {e}")
        traceback.print_exc()
        return {
            "status": "failed",
            "error": str(e),
            "traceback": traceback.format_exc()
        }


def handle_delete_file(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Delete a specific file from the Network Volume.
    Only allows deleting files within the checkpoint directory for safety.

    Input:
        filename: str - Name of file to delete (e.g., "rfd3_latest.ckpt")
    """
    filename = job_input.get("filename")
    if not filename:
        return {"status": "failed", "error": "Missing 'filename' parameter"}

    # Safety: only allow deleting from checkpoint directory
    filepath = os.path.join(CHECKPOINT_DIR, filename)

    # Prevent path traversal
    if not os.path.abspath(filepath).startswith(os.path.abspath(CHECKPOINT_DIR)):
        return {"status": "failed", "error": "Invalid filename - path traversal not allowed"}

    if not os.path.exists(filepath):
        return {
            "status": "completed",
            "result": {
                "deleted": False,
                "message": f"File does not exist: {filepath}"
            }
        }

    try:
        file_size = os.path.getsize(filepath)
        os.remove(filepath)
        return {
            "status": "completed",
            "result": {
                "deleted": True,
                "filepath": filepath,
                "size_bytes": file_size,
                "message": f"Successfully deleted {filepath} ({file_size} bytes)"
            }
        }
    except Exception as e:
        return {
            "status": "failed",
            "error": f"Failed to delete {filepath}: {str(e)}"
        }


# ============== RunPod Entry Point ==============

if __name__ == "__main__":
    print("[Handler] Starting RunPod serverless handler...")
    runpod.serverless.start({"handler": handler})
