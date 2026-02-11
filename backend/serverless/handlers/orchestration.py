"""Orchestration handlers for NL parsing, scaffold search, and scout filtering.

Tasks: ai_parse, scaffold_search, scout_filter
"""

import os
import traceback
from typing import Dict, Any

from inference_utils import run_mpnn_inference, run_rf3_inference


def handle_ai_parse(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Parse a natural language design query into structured DesignIntent using Claude AI.

    Falls back to keyword-based parsing if no Claude API key is available.

    Input:
        query: str - Natural language design request
        claude_api_key: str - Optional Claude API key (falls back to env vars)

    Returns:
        Parsed DesignIntent fields with parser_type indicator.
    """
    query = job_input.get("query", "")
    if not query:
        return {"status": "failed", "error": "Must provide 'query' parameter"}

    api_key = (
        job_input.get("claude_api_key")
        or os.environ.get("CLAUDE_API_KEY")
        or os.environ.get("ANTHROPIC_API_KEY")
    )

    try:
        from nl_design_parser import create_parser
        parser = create_parser(api_key=api_key, use_fallback=True)
        intent = parser.parse(query)
        result = intent.to_dict()
        result["parser_type"] = "ai" if api_key else "fallback"
        return {"status": "completed", "result": result}
    except Exception as e:
        print(f"[AI Parse] Error: {e}")
        traceback.print_exc()
        return {"status": "failed", "error": f"Parse error: {str(e)}"}


def handle_scaffold_search(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Search RCSB PDB for scaffold candidates containing a metal-ligand complex.

    Input:
        metal: str - Metal symbol (e.g., "TB", "ZN", "CA")
        ligand_name: str - Human-readable ligand name (e.g., "citrate")
        ligand_code: str - Optional PDB 3-letter code (e.g., "CIT")
        ligand_smiles: str - Optional SMILES for scoring
        resolution_max: float - Max resolution filter (default 3.0)
        limit: int - Max PDB hits to validate (default 10)
        fetch_pdb: bool - If true, fetch PDB content for candidates (default true)

    Returns:
        Search results with candidates, scores, optional PDB content.
    """
    try:
        from scaffold_search import search_scaffold_candidates
    except ImportError as e:
        return {"status": "failed", "error": f"Scaffold search not available: {e}"}

    metal = job_input.get("metal")
    ligand_name = job_input.get("ligand_name", "")
    ligand_code = job_input.get("ligand_code")
    ligand_smiles = job_input.get("ligand_smiles")
    resolution_max = job_input.get("resolution_max", 3.0)
    limit = job_input.get("limit", 10)
    fetch_pdb = job_input.get("fetch_pdb", True)

    if not metal:
        return {"status": "failed", "error": "Must provide metal symbol (e.g., TB, ZN, CA)"}
    if not ligand_name and not ligand_code:
        return {"status": "failed", "error": "Must provide ligand_name or ligand_code"}

    print(f"[ScaffoldSearch] Searching for {metal} + {ligand_name or ligand_code}...")

    try:
        search_result = search_scaffold_candidates(
            metal=metal,
            ligand_name=ligand_name,
            ligand_code=ligand_code,
            ligand_smiles=ligand_smiles,
            resolution_max=resolution_max,
            limit=limit,
        )
    except Exception as e:
        print(f"[ScaffoldSearch] Search failed: {e}")
        return {"status": "failed", "error": f"Scaffold search error: {e}"}

    result_dict = search_result.to_dict()

    # Optionally fetch PDB content for candidates so they can be visualized
    if fetch_pdb and search_result.candidates:
        import requests as _req
        candidate_pdbs = {}  # pdb_id -> pdb_content
        for candidate in search_result.candidates:
            pdb_id = candidate.pdb_id
            try:
                pdb_url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
                resp = _req.get(pdb_url, timeout=30)
                resp.raise_for_status()
                candidate_pdbs[pdb_id] = resp.text
                print(f"[ScaffoldSearch] Fetched PDB {pdb_id} ({len(resp.text)} bytes)")
            except Exception as e:
                print(f"[ScaffoldSearch] Failed to fetch PDB {pdb_id}: {e}")
                candidate_pdbs[pdb_id] = None
        result_dict["candidate_pdbs"] = candidate_pdbs
        # Also set best_pdb_content for backwards compatibility
        if search_result.best_candidate:
            result_dict["best_pdb_content"] = candidate_pdbs.get(search_result.best_candidate.pdb_id)

    print(f"[ScaffoldSearch] Done: {search_result.num_pdb_hits} hits, "
          f"{search_result.num_validated} validated, "
          f"action={search_result.recommended_action}")

    return {"status": "completed", "result": result_dict}


def handle_scout_filter(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Filter backbones by generating 1 scout sequence each and validating.

    For each backbone PDB: run MPNN (1 seq) -> RF3 -> check pTM/pLDDT thresholds.
    Includes CPU-only pre-filter to catch broken backbones before GPU work.
    Returns only passing backbones with scout metrics.
    """
    backbone_pdbs = job_input.get("backbone_pdbs", [])
    if not backbone_pdbs:
        return {
            "status": "completed",
            "result": {
                "filtered_pdbs": [],
                "original_count": 0,
                "passing_count": 0,
                "scout_results": [],
                "ptm_threshold": job_input.get("ptm_threshold", 0.6),
                "plddt_threshold": job_input.get("plddt_threshold", 0.65),
            },
        }

    ptm_threshold = job_input.get("ptm_threshold", 0.6)
    plddt_threshold = job_input.get("plddt_threshold", 0.65)
    target_metal = job_input.get("target_metal")
    ligand_smiles = job_input.get("ligand_smiles")
    ligand_name = job_input.get("ligand_name")
    pre_filter_only = job_input.get("pre_filter_only", False)

    # --- CPU Pre-filter: catch obviously bad backbones before GPU work ---
    pre_filter_results = []
    pre_filter_passed_indices = set()
    pre_filter_failed_count = 0

    try:
        from backbone_pre_filter import run_backbone_pre_filter

        for i, backbone_pdb in enumerate(backbone_pdbs):
            try:
                pf_result = run_backbone_pre_filter(
                    backbone_pdb,
                    target_metal=target_metal,
                    ligand_name=ligand_name,
                )
                pre_filter_results.append(pf_result)
                if pf_result["passed"]:
                    pre_filter_passed_indices.add(i)
                    print(f"[Scout Pre-filter] Backbone {i}: PASSED")
                else:
                    pre_filter_failed_count += 1
                    print(f"[Scout Pre-filter] Backbone {i}: FAILED — {pf_result['failed_checks']}")
            except Exception as e:
                # Non-fatal: allow backbone through
                pre_filter_results.append({"passed": True, "error": str(e), "checks": {}, "failed_checks": [], "skipped_checks": []})
                pre_filter_passed_indices.add(i)
                print(f"[Scout Pre-filter] Backbone {i}: ERROR (allowing through) — {e}")

        print(f"[Scout Pre-filter] {len(pre_filter_passed_indices)}/{len(backbone_pdbs)} passed pre-filter")
    except ImportError:
        # Module not available — allow all backbones through
        print("[Scout Pre-filter] backbone_pre_filter module not available, skipping")
        pre_filter_passed_indices = set(range(len(backbone_pdbs)))

    # Pre-filter only mode: return results without running MPNN+RF3 GPU loop
    if pre_filter_only:
        print(f"[Scout] Pre-filter only mode: returning {len(pre_filter_passed_indices)}/{len(backbone_pdbs)} passed")
        scout_results_pf = []
        filtered_pdbs_pf = []
        for i, backbone_pdb in enumerate(backbone_pdbs):
            passed = i in pre_filter_passed_indices
            scout_results_pf.append({
                "backbone_index": i,
                "ptm": 0.0,
                "plddt": 0.0,
                "passed": passed,
                "sequence": "",
                "pre_filter_failed": not passed,
            })
            if passed:
                filtered_pdbs_pf.append(backbone_pdb)

        return {
            "status": "completed",
            "result": {
                "filtered_pdbs": filtered_pdbs_pf,
                "original_count": len(backbone_pdbs),
                "passing_count": len(filtered_pdbs_pf),
                "scout_results": scout_results_pf,
                "ptm_threshold": ptm_threshold,
                "plddt_threshold": plddt_threshold,
                "pre_filter_results": pre_filter_results,
                "pre_filter_passed": len(pre_filter_passed_indices),
                "pre_filter_failed": pre_filter_failed_count,
                "pre_filter_only": True,
            },
        }

    filtered_pdbs = []
    scout_results = []

    for i, backbone_pdb in enumerate(backbone_pdbs):
        scout_entry = {
            "backbone_index": i,
            "ptm": 0.0,
            "plddt": 0.0,
            "passed": False,
            "sequence": "",
            "pre_filter_failed": i not in pre_filter_passed_indices,
        }

        # Skip GPU work for pre-filter failures
        if i not in pre_filter_passed_indices:
            print(f"[Scout] Backbone {i}: skipped (pre-filter failed)")
            scout_results.append(scout_entry)
            continue

        try:
            # 1. Generate 1 scout sequence via MPNN
            mpnn_params = {
                "pdb_content": backbone_pdb,
                "num_sequences": 1,
                "temperature": 0.1,
                "model_type": "ligand_mpnn",
            }
            mpnn_result = run_mpnn_inference(mpnn_params)
            sequences = mpnn_result.get("sequences", [])

            # Parse sequence from FASTA if needed
            sequence = ""
            if sequences:
                first = sequences[0]
                content = first.get("content", first.get("sequence", "")) if isinstance(first, dict) else str(first)
                if content.startswith(">"):
                    for line in content.splitlines():
                        if not line.startswith(">") and line.strip():
                            sequence = line.strip()
                            break
                else:
                    sequence = content.strip()

            if not sequence:
                print(f"[Scout] Backbone {i}: no sequence generated, skipping")
                scout_results.append(scout_entry)
                continue

            scout_entry["sequence"] = sequence

            # 2. Validate with RF3 (with ligand + metal context)
            rf3_params: Dict[str, Any] = {"sequence": sequence, "name": f"scout_{i:03d}"}
            if ligand_smiles:
                rf3_params["ligand_smiles"] = ligand_smiles
            if target_metal:
                rf3_params["metal"] = target_metal

            rf3_result = run_rf3_inference(**rf3_params)
            predictions = rf3_result.get("predictions", [{}])
            pred = predictions[0] if predictions else rf3_result
            confidences = pred.get("confidences", rf3_result.get("confidences", {}))
            summary = confidences.get("summary_confidences", confidences) if isinstance(confidences, dict) else {}

            ptm = float(summary.get("ptm", pred.get("ptm", 0.0)))
            plddt = float(summary.get("overall_plddt", pred.get("mean_plddt", 0.0)))

            scout_entry["ptm"] = ptm
            scout_entry["plddt"] = plddt

            # 3. Check thresholds
            if ptm >= ptm_threshold and plddt >= plddt_threshold:
                scout_entry["passed"] = True
                filtered_pdbs.append(backbone_pdb)
                print(f"[Scout] Backbone {i}: PASSED (ptm={ptm:.3f}, plddt={plddt:.3f})")
            else:
                print(f"[Scout] Backbone {i}: FAILED (ptm={ptm:.3f}, plddt={plddt:.3f})")

        except Exception as e:
            print(f"[Scout] Backbone {i}: error — {e}")

        scout_results.append(scout_entry)

    print(f"[Scout] Complete: {len(filtered_pdbs)}/{len(backbone_pdbs)} backbones passed")

    return {
        "status": "completed",
        "result": {
            "filtered_pdbs": filtered_pdbs,
            "original_count": len(backbone_pdbs),
            "passing_count": len(filtered_pdbs),
            "scout_results": scout_results,
            "ptm_threshold": ptm_threshold,
            "plddt_threshold": plddt_threshold,
            "pre_filter_results": pre_filter_results,
            "pre_filter_passed": len(pre_filter_passed_indices),
            "pre_filter_failed": pre_filter_failed_count,
        },
    }
