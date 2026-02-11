"""Unified design analysis handler with filter evaluation.

Tasks: analyze_design
"""

import traceback
from typing import Dict, Any, Optional, List, Tuple


def _flatten_analysis_metrics(analysis: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extract flat metric dict from UnifiedDesignAnalyzer's nested output.

    Maps nested analysis results into a flat dict suitable for filter evaluation.
    """
    flat: Dict[str, Any] = {}

    analyses = analysis.get("analyses", {})

    # Structure confidence -> plddt
    sc = analyses.get("structure_confidence", {})
    if sc.get("status") == "success":
        metrics = sc.get("metrics", {})
        if "plddt" in metrics:
            flat["plddt"] = metrics["plddt"]

    # Metal coordination -> coordination_number, geometry_rmsd
    mc = analyses.get("metal_coordination", {})
    if mc.get("status") == "success":
        metrics = mc.get("metrics", {})
        if "coordination_number" in metrics:
            flat["coordination_number"] = metrics["coordination_number"]
        if "geometry_rmsd" in metrics:
            flat["geometry_rmsd"] = metrics["geometry_rmsd"]
        if "ligand_coordination" in metrics:
            flat["ligand_coordination"] = metrics["ligand_coordination"]
        if "protein_coordination" in metrics:
            flat["protein_coordination"] = metrics["protein_coordination"]

    # Ligand interactions
    li = analyses.get("ligand_interactions", {})
    if li.get("status") == "success":
        metrics = li.get("metrics", {})
        for key in ["total_contacts", "hydrogen_bonds", "hydrophobic_contacts", "pi_stacking", "salt_bridges"]:
            if key in metrics:
                flat[key] = metrics[key]
        # Alias: ligand_contacts = total_contacts (used by chemistry-aware presets)
        if "total_contacts" in metrics:
            flat["ligand_contacts"] = metrics["total_contacts"]

    # Sequence composition
    seq = analyses.get("sequence_composition", {})
    if seq.get("status") == "success":
        metrics = seq.get("metrics", {})
        if "total_residues" in metrics:
            flat["total_residues"] = metrics["total_residues"]
        ala = metrics.get("alanine", {})
        if "percentage" in ala:
            flat["alanine_pct"] = ala["percentage"]
        arom = metrics.get("aromatics", {})
        if "percentage" in arom:
            flat["aromatic_pct"] = arom["percentage"]

    # Interface quality
    iq = analyses.get("interface_quality", {})
    if iq.get("status") == "success":
        metrics = iq.get("metrics", {})
        for key in ["dSASA", "contacts", "interface_residues"]:
            if key in metrics:
                flat[key] = metrics[key]

    # Topology
    topo = analyses.get("topology", {})
    if topo.get("status") == "success":
        metrics = topo.get("metrics", {})
        if "valid" in metrics:
            flat["topology_valid"] = metrics["valid"]

    return flat


def _pick_filter_preset(
    design_type: Optional[str],
    metal: Optional[str] = None,
    ligand_name: Optional[str] = None,
    tier: str = "standard",
) -> Tuple[str, Dict[str, Any]]:
    """Pick or generate the appropriate filter preset.

    If metal is provided and recognized, generates chemistry-aware thresholds
    from metal_chemistry.py data (CN ranges, HSAB class). Falls back to
    hardcoded FILTER_PRESETS otherwise.

    Returns:
        Tuple of (preset_name, preset_dict).
    """
    # Chemistry-aware path: generate thresholds from metal database
    if metal:
        try:
            from filter_evaluator import generate_chemistry_aware_preset
            name, preset = generate_chemistry_aware_preset(
                metal=metal,
                ligand_name=ligand_name,
                tier=tier,
                design_type=design_type or "metal_binding",
            )
            return name, preset
        except Exception:
            pass  # Fall through to hardcoded presets

    # Hardcoded fallback
    from filter_evaluator import FILTER_PRESETS
    if design_type in ("metal", "metal_ligand", "metal_monomer", "metal_ligand_monomer"):
        name = "metal_binding"
    elif design_type in ("enzyme", "scaffold", "enzyme_scaffold"):
        name = "enzyme_scaffold"
    else:
        name = "scout_relaxed"
    return name, dict(FILTER_PRESETS.get(name, {}))


def _evaluate_flat_metrics(
    metrics: Dict[str, Any], preset: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Evaluate flat metrics dict against a filter preset.

    Returns dict with 'passed' bool and 'failed_filters' list.
    """
    failed_filters: List[Dict[str, Any]] = []

    for metric_name, threshold in preset.items():
        value = metrics.get(metric_name)
        if value is None:
            continue

        passed = True
        if "min" in threshold and value < threshold["min"]:
            passed = False
        if "max" in threshold and value > threshold["max"]:
            passed = False

        if not passed:
            failed_filters.append({
                "metric": metric_name,
                "value": value,
                "threshold": threshold,
            })

    return {
        "passed": len(failed_filters) == 0,
        "failed_filters": failed_filters,
    }


def handle_analyze_design(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Run UnifiedDesignAnalyzer + FILTER_PRESETS filtering on a design PDB.

    Input:
        pdb_content: PDB file content (required)
        metal_type: Metal element code (auto-detected if omitted)
        ligand_name: Ligand name for context
        design_type: Design type for filter preset selection
        design_params: Optional params dict passed to analyzer
        filter_tier: str - Filter stringency tier (default: "standard")

    Returns:
        Enriched analysis with flat metrics, filter pass/fail, and full nested data.
    """
    try:
        from unified_analyzer import UnifiedDesignAnalyzer

        pdb_content = job_input.get("pdb_content")
        if not pdb_content:
            return {"status": "failed", "error": "pdb_content is required"}

        # Auto-convert CIF to PDB format (RF3 outputs CIF despite .pdb extension)
        if pdb_content.lstrip().startswith("data_"):
            try:
                from biotite.structure.io import pdbx, pdb as pdb_io
                import io as _io
                cif_file = pdbx.CIFFile.read(_io.StringIO(pdb_content))
                # RF3 CIF may lack pdbx_PDB_model_num - catch any exception from model=1
                try:
                    atom_array = pdbx.get_structure(cif_file, model=1)
                except Exception:
                    # No model field or other issue - get first structure without model specification
                    atom_array = pdbx.get_structure(cif_file)
                pdb_file = pdb_io.PDBFile()
                pdb_io.set_structure(pdb_file, atom_array)
                pdb_str = _io.StringIO()
                pdb_file.write(pdb_str)
                pdb_content = pdb_str.getvalue()
                # Fix CIF residue names (L:0, L:1) -> standard codes using element/atom_name
                from analysis_types import METAL_CODES
                fixed_lines = []
                metal_found = False
                for line in pdb_content.split("\n"):
                    if line.startswith("HETATM"):
                        res_name = line[17:20].strip()
                        if ":" in res_name:
                            # Check element column (cols 76-78) or atom_name (cols 12-16)
                            element = line[76:78].strip().upper() if len(line) > 76 else ""
                            atom_name = line[12:16].strip().rstrip("0123456789").upper()
                            new_res = None
                            if element and element in METAL_CODES:
                                new_res = element
                                metal_found = True
                            elif atom_name and atom_name in METAL_CODES:
                                new_res = atom_name
                                metal_found = True
                            if new_res:
                                line = line[:17] + f"{new_res:>3}" + line[20:]
                    fixed_lines.append(line)
                pdb_content = "\n".join(fixed_lines)
                print(f"[analyze_design] Converted CIF input to PDB format (metal_found={metal_found})")
            except Exception as e:
                print(f"[analyze_design] Warning: CIF-to-PDB conversion failed: {e}")

        metal_type = job_input.get("metal_type")
        ligand_name = job_input.get("ligand_name")
        design_type = job_input.get("design_type")
        design_params = job_input.get("design_params", {})
        filter_tier = job_input.get("filter_tier", "standard")

        # 1. Run UnifiedDesignAnalyzer
        analyzer = UnifiedDesignAnalyzer()
        analysis = analyzer.analyze(
            pdb_content,
            design_params,
            metal_type=metal_type,
        )

        # 2. Flatten nested analysis metrics
        flat_metrics = _flatten_analysis_metrics(analysis)

        # 3. Pick filter preset — chemistry-aware if metal provided
        preset_name, preset = _pick_filter_preset(
            design_type,
            metal=metal_type,
            ligand_name=ligand_name,
            tier=filter_tier,
        )

        # 4. Evaluate against preset
        filter_result = _evaluate_flat_metrics(flat_metrics, preset)

        # 5. Build chemistry context for frontend display
        chemistry_context = None
        if metal_type:
            try:
                from metal_chemistry import get_hsab_class_simple, get_coordination_number_range, METAL_DATABASE
                metal_upper = metal_type.upper()
                if metal_upper in METAL_DATABASE:
                    default_ox = METAL_DATABASE[metal_upper].get("default_oxidation", 2)
                    cn_min, cn_max = get_coordination_number_range(metal_upper, default_ox)
                    chemistry_context = {
                        "metal": metal_upper,
                        "hsab_class": get_hsab_class_simple(metal_upper),
                        "formal_cn_range": [cn_min, cn_max],
                        "tier": filter_tier,
                    }
            except Exception:
                pass

        return {
            "status": "completed",
            "result": {
                "design_type": analysis.get("design_type", "unknown"),
                "metrics": flat_metrics,
                "filter_preset": preset_name,
                "filter_passed": filter_result["passed"],
                "failed_filters": filter_result["failed_filters"],
                "filter_tier": filter_tier,
                "chemistry_context": chemistry_context,
                "auto_detected": analysis.get("auto_detected", {}),
                "analyses": analysis.get("analyses", {}),
            },
        }
    except Exception as e:
        traceback.print_exc()
        return {
            "status": "failed",
            "error": f"Design analysis failed: {str(e)}",
        }
