"""Validation and analysis handlers.

Tasks: validate_design, binding_eval, fastrelax, detect_hotspots,
       analyze_conservation, interaction_analysis
"""

import traceback
from typing import Dict, Any

from binding_analysis import (
    evaluate_binding,
    check_gnina_available,
    check_steric_clashes,
    to_python_types,
)
from rosetta_utils import (
    fastrelax_with_ligand,
    check_pyrosetta_available,
)
from hotspot_detection import detect_hotspots_sasa


def handle_validate_design(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Run full design validation pipeline: MPNN sequence design + ESMFold structure prediction.

    This validates if an RFD3 backbone design is "designable" - whether MPNN-designed
    sequences fold back to the intended structure.

    Input:
        pdb_content: PDB file content as string (required unless pdb_path provided)
        pdb_path: Path to PDB file (alternative to pdb_content)
        num_sequences: Number of MPNN sequences to generate (default: 8)
        use_ligandmpnn: Use LigandMPNN instead of ProteinMPNN (default: True)
        ligand_name: Ligand residue name, e.g., "UNL" (optional, auto-detected)
        metal_type: Metal type code, e.g., "DY", "TB" (optional, auto-detected)
        max_backbone_rmsd: Maximum acceptable backbone RMSD (default: 1.5)
        min_plddt: Minimum acceptable mean pLDDT (default: 0.7)
        min_metal_coordination: Minimum metal coordination number (default: 6)
        min_protein_donors: Minimum protein donor atoms for metal (default: 1)
        run_pyrosetta: Whether to run PyRosetta refinement (default: True)
        temperature: MPNN sampling temperature (default: 0.1)
        use_rf3: Whether to also validate with RF3 (default: True)
        ligand_smiles: SMILES string for ligand-aware RF3 prediction (optional)

    Returns:
        Validation results including candidates, metrics, and pass/fail info.
    """
    try:
        from design_validation_pipeline import DesignValidationPipeline
    except ImportError:
        return {
            "status": "failed",
            "error": "Validation pipeline not available. Check design_validation_pipeline.py."
        }

    # Get PDB content from either direct content or file path
    pdb_content = job_input.get("pdb_content")
    pdb_path = job_input.get("pdb_path")

    if not pdb_content and not pdb_path:
        return {"status": "failed", "error": "Missing 'pdb_content' or 'pdb_path' parameter"}

    # If path provided, read content
    if pdb_path and not pdb_content:
        try:
            with open(pdb_path, 'r') as f:
                pdb_content = f.read()
        except Exception as e:
            return {"status": "failed", "error": f"Failed to read PDB file: {e}"}

    # Initialize pipeline
    pipeline = DesignValidationPipeline()

    # Run validation
    try:
        result = pipeline.validate(
            pdb_content=pdb_content,
            pdb_path=pdb_path,
            num_sequences=job_input.get("num_sequences", 8),
            use_ligandmpnn=job_input.get("use_ligandmpnn", True),
            ligand_name=job_input.get("ligand_name"),
            metal_type=job_input.get("metal_type"),
            max_backbone_rmsd=job_input.get("max_backbone_rmsd", 1.5),
            min_plddt=job_input.get("min_plddt", 0.7),
            min_metal_coordination=job_input.get("min_metal_coordination", 6),
            min_protein_donors=job_input.get("min_protein_donors", 1),
            run_pyrosetta=job_input.get("run_pyrosetta", True),
            temperature=job_input.get("temperature", 0.1),
            use_rf3=job_input.get("use_rf3", True),
            ligand_smiles=job_input.get("ligand_smiles"),
        )

        # Convert ValidationResult to dict
        return {
            "status": "completed",
            **result.to_dict()
        }
    except Exception as e:
        return {
            "status": "failed",
            "error": str(e),
            "traceback": traceback.format_exc()
        }


def handle_binding_eval(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Evaluate protein-protein or protein-ligand binding quality.

    Combines: steric clash check, interface analysis, GNINA scoring.

    Input:
        pdb_content: PDB file content (protein complex or protein-ligand)
        chain_a: First chain for interface analysis (default: "A")
        chain_b: Second chain (default: "B", or "HETATM" for ligand)
        ligand_smiles: SMILES string for GNINA docking (optional)
        ligand_sdf: SDF content for GNINA docking (optional)
        run_gnina: Whether to run GNINA scoring (default: True if ligand provided)
        check_clashes: Whether to run steric clash check (default: True)

    Returns:
        Combined evaluation with clash_check, interface_analysis, gnina_scoring, summary.
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
    return to_python_types(result)


def handle_fastrelax(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Run Rosetta FastRelax to refine protein-ligand complexes and resolve clashes.

    Ligand is kept FIXED while protein moves around it.

    Input:
        pdb_content: str - PDB file content with protein and ligand
        ligand_smiles: str - SMILES string for GNINA scoring (optional)
        max_iter: int - Maximum minimization iterations (default: 200)
        interface_only: bool - Only relax interface residues (default: True)
        interface_distance: float - Distance cutoff for interface (default: 8.0)
        constrain_coords: bool - Constrain backbone to starting coords (default: True)
        run_gnina_before: bool - Score with GNINA before relaxation (default: True)
        run_gnina_after: bool - Score with GNINA after relaxation (default: True)

    Returns:
        Relaxed PDB, energy changes, optional GNINA before/after comparison.
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


def handle_detect_hotspots(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Auto-detect binding hotspots for protein binder design.

    Uses SASA analysis combined with spatial clustering.

    Input:
        target_pdb: str (required) - Target protein PDB content
        target_chain: str (default: "A") - Chain to analyze
        method: str (default: "exposed_clustered") - Detection method
        max_hotspots: int (default: 3) - Maximum hotspots to return
        prefer_hydrophobic: bool (default: True) - Prioritize hydrophobic residues

    Returns:
        Hotspot positions, method used, cluster details, residue properties.
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


def handle_analyze_conservation(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Analyze evolutionary conservation using ConSurf methodology.

    Input:
        pdb_content: str (required) - PDB file content
        chain: str (default: "A") - Chain to analyze
        method: str (default: "bayesian") - Rate4Site method

    Returns:
        Per-residue conservation grades, conserved/variable positions, MSA depth.
    """
    import asyncio

    try:
        from conservation_analyzer import ConservationAnalyzer
    except ImportError:
        return {"status": "error", "error": "Conservation analyzer not available"}

    pdb_content = job_input.get("pdb_content")
    if not pdb_content:
        return {"status": "error", "error": "Missing 'pdb_content' parameter"}

    chain = job_input.get("chain", "A")
    method = job_input.get("method", "bayesian")

    print(f"[Conservation] Analyzing conservation for chain {chain}, method={method}")

    try:
        analyzer = ConservationAnalyzer()
        result = asyncio.run(analyzer.analyze(
            pdb_content=pdb_content,
            chain=chain,
            method=method,
        ))

        response = result.to_dict()
        response["status"] = "success"
        response["highly_conserved_positions"] = response.pop("highly_conserved", [])
        response["conserved_positions"] = response.pop("conserved", [])
        response["variable_positions"] = response.pop("variable", [])
        response["reliable"] = result.msa_depth >= 30

        print(f"[Conservation] Completed: {len(response['grades'])} residues, "
              f"MSA depth={result.msa_depth}, "
              f"conserved={len(response['highly_conserved_positions'])}, "
              f"variable={len(response['variable_positions'])}")

        return to_python_types(response)

    except Exception as e:
        error_msg = str(e)
        print(f"[Conservation] Analysis failed: {error_msg}")
        traceback.print_exc()
        return {
            "status": "error",
            "error": f"Conservation analysis failed: {error_msg}",
        }


def handle_interaction_analysis(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Analyze protein-ligand interactions using PLIP or distance-based methods.

    Input:
        pdb_content: str (required) - PDB content with protein and ligand
        ligand_name: str (default: "UNL") - Ligand residue name
        include_visualization: bool (default: True) - Include visualization data
        include_recommendations: bool (default: True) - Include binding recommendations
        ligand_has_aromatics: bool (default: False) - For pi-stacking recommendations

    Returns:
        Interaction counts, key residues, details, recommendations, AI summary.
    """
    pdb_content = job_input.get("pdb_content")
    if not pdb_content:
        return {"status": "error", "error": "Missing 'pdb_content' parameter"}

    try:
        from shared.interaction_analysis import (
            analyze_all_interactions,
            format_for_frontend,
            format_for_ai_assistant,
            generate_recommendations,
        )
    except ImportError:
        return {"status": "error", "error": "Interaction analysis module not available"}

    ligand_name = job_input.get("ligand_name", "UNL")
    include_visualization = job_input.get("include_visualization", True)
    include_recommendations = job_input.get("include_recommendations", True)
    ligand_has_aromatics = job_input.get("ligand_has_aromatics", False)

    print(f"[InteractionAnalysis] Analyzing interactions with ligand {ligand_name}")

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
