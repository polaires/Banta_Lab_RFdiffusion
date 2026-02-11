"""Unified design analysis handler with filter evaluation.

Tasks: analyze_design
"""

import traceback
from typing import Dict, Any, Optional, List, Tuple


def _compute_ligand_rmsd(
    predicted_pdb: str,
    backbone_pdb: str,
    metal_type: Optional[str] = None,
) -> Optional[float]:
    """
    Compute backbone-aligned ligand RMSD between predicted and template structures.

    The RFD3 paper (bioRxiv 2025.09.18.676967v2) uses backbone-aligned ligand RMSD < 5 Å
    as an AF3 success criterion for ligand binding.

    Process:
      1. Extract CA atoms from both PDBs for backbone alignment
      2. Compute Kabsch rotation/translation on CA atoms
      3. Apply rotation to predicted ligand HETATM coordinates
      4. Compute RMSD between aligned predicted and template ligand atoms

    Returns None if either PDB lacks HETATM or backbone alignment fails.
    """
    if not predicted_pdb or not backbone_pdb:
        return None

    try:
        import numpy as np
    except ImportError:
        return None

    try:
        from analysis_types import METAL_CODES
    except ImportError:
        METAL_CODES = {
            "FE", "ZN", "CA", "MG", "MN", "CO", "CU", "NI", "MO", "W",
            "TB", "EU", "GD", "DY", "SM", "ND", "LA", "CE", "PR",
        }

    def _parse_hetatm_coords(pdb_str: str, exclude_metals: bool = True) -> List[List[float]]:
        """Extract HETATM non-metal coordinates."""
        coords = []
        for line in pdb_str.split("\n"):
            if not line.startswith("HETATM"):
                continue
            res_name = line[17:20].strip().upper()
            # Skip metals — we want ligand atoms only
            if exclude_metals and res_name in METAL_CODES:
                continue
            if metal_type and res_name.upper() == metal_type.upper():
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
            except (ValueError, IndexError):
                continue
        return coords

    def _parse_ca_coords(pdb_str: str) -> List[List[float]]:
        """Extract CA atom coordinates."""
        coords = []
        for line in pdb_str.split("\n"):
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
            except (ValueError, IndexError):
                continue
        return coords

    # Get ligand coords from both structures
    pred_lig = _parse_hetatm_coords(predicted_pdb)
    ref_lig = _parse_hetatm_coords(backbone_pdb)
    if not pred_lig or not ref_lig:
        return None
    if len(pred_lig) != len(ref_lig):
        return None  # Atom count mismatch

    # Get CA coords for backbone alignment
    pred_ca = _parse_ca_coords(predicted_pdb)
    ref_ca = _parse_ca_coords(backbone_pdb)
    if not pred_ca or not ref_ca:
        return None
    n_ca = min(len(pred_ca), len(ref_ca))
    if n_ca < 10:
        return None  # Too few CA atoms for reliable alignment

    pred_ca_arr = np.array(pred_ca[:n_ca])
    ref_ca_arr = np.array(ref_ca[:n_ca])

    # Kabsch alignment on CA atoms
    pred_center = pred_ca_arr.mean(axis=0)
    ref_center = ref_ca_arr.mean(axis=0)
    p = pred_ca_arr - pred_center
    q = ref_ca_arr - ref_center
    h = p.T @ q
    u, s, vt = np.linalg.svd(h)
    d = np.linalg.det(vt.T @ u.T)
    sign_matrix = np.eye(3)
    sign_matrix[2, 2] = np.sign(d)
    rotation = vt.T @ sign_matrix @ u.T

    # Apply rotation to predicted ligand coords
    pred_lig_arr = np.array(pred_lig) - pred_center
    aligned_pred_lig = (rotation @ pred_lig_arr.T).T + ref_center
    ref_lig_arr = np.array(ref_lig)

    # Compute RMSD
    diff = aligned_pred_lig - ref_lig_arr
    rmsd = float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))
    return round(rmsd, 3)


def _fix_rf3_hetatm(
    protein_pdb: str,
    backbone_pdb: str,
    metal_type: Optional[str] = None,
) -> str:
    """Fix RF3 HETATM naming to standard CCD conventions.

    RF3 outputs HETATM with non-standard names:
      - Metal: residue 'L:0', atom 'CA0' (should be 'CA', 'CA')
      - Ligand: residue 'L:1', atoms 'O0,C0,...' (should be 'PQQ', 'O5,N6,...')

    This function renames RF3 HETATM records using the backbone PDB as reference,
    while keeping RF3's predicted 3D positions (which are self-consistent with
    the RF3 protein structure). Falls back to grafting backbone HETATM when RF3
    has no HETATM at all.

    Handles both raw RF3 PDB output and post-CIF-conversion PDB (where the CIF
    converter may have partially fixed metal names but left ligand as L:1).
    """
    if not backbone_pdb:
        return protein_pdb

    # Import METAL_CODES for robust metal identification
    try:
        from analysis_types import METAL_CODES
    except ImportError:
        METAL_CODES = {
            "FE", "ZN", "CA", "MG", "MN", "CO", "CU", "NI", "MO", "W",
            "TB", "EU", "GD", "DY", "SM", "ND", "LA", "CE", "PR", "PM",
            "HO", "ER", "TM", "YB", "LU", "NA", "K", "CD", "HG",
        }

    # Extract metal and ligand residue names from backbone
    bb_metal_res = None
    bb_ligand_res = None
    for line in backbone_pdb.split("\n"):
        if not line.startswith("HETATM"):
            continue
        res_name = line[17:20].strip()
        if not res_name:
            continue
        if metal_type and res_name.upper() == metal_type.upper():
            bb_metal_res = res_name
        elif res_name.upper() in METAL_CODES and not bb_metal_res:
            bb_metal_res = res_name
        elif not bb_ligand_res and res_name != bb_metal_res:
            bb_ligand_res = res_name

    if not bb_metal_res and not bb_ligand_res:
        return protein_pdb

    # Collect HETATM lines from RF3 output
    hetatm_lines = [
        l for l in protein_pdb.split("\n") if l.startswith("HETATM")
    ]

    if not hetatm_lines:
        # RF3 output is protein-only — graft HETATM from backbone
        bb_hetatm = [l for l in backbone_pdb.split("\n") if l.startswith("HETATM")]
        bb_conect = [l for l in backbone_pdb.split("\n") if l.startswith("CONECT")]
        if bb_hetatm:
            merged = [l for l in protein_pdb.split("\n") if not l.startswith("END")]
            merged.append("TER")
            merged.extend(bb_hetatm)
            if bb_conect:
                merged.extend(bb_conect)
            merged.append("END")
            print(f"[_fix_rf3_hetatm] Grafted {len(bb_hetatm)} HETATM from backbone")
            return "\n".join(merged)
        return protein_pdb

    # Check if any HETATM has non-standard naming (colon in residue name)
    has_nonstandard = any(":" in l[17:20] for l in hetatm_lines)
    if not has_nonstandard:
        return protein_pdb  # Already standard naming

    # Rename non-standard HETATM in-place (keep RF3 positions)
    fixed_lines = []
    metal_count = 0
    ligand_count = 0
    for line in protein_pdb.split("\n"):
        if line.startswith("HETATM"):
            res_name = line[17:20].strip()
            if ":" in res_name:
                # Identify if this HETATM is metal or ligand
                is_metal = False

                # Method 1: element symbol column (76-78) — most reliable
                if len(line) > 76:
                    element = line[76:78].strip().upper()
                    if element in METAL_CODES:
                        is_metal = True

                # Method 2: atom name (12-16), strip trailing digits
                if not is_metal:
                    atom_name = line[12:16].strip().rstrip("0123456789").upper()
                    if atom_name in METAL_CODES:
                        is_metal = True

                # Method 3: L:0 convention (RF3 uses L:0 for metal)
                if not is_metal and res_name == "L:0" and bb_metal_res:
                    is_metal = True

                if is_metal and bb_metal_res:
                    metal_upper = bb_metal_res.upper()
                    # PDB atom name: 2-char elements left-justify at col 13
                    if len(metal_upper) >= 2:
                        atom_field = f"{metal_upper:<4}"
                    else:
                        atom_field = f" {metal_upper:<3}"
                    line = (
                        line[:12] + atom_field + line[16:17]
                        + f"{metal_upper:>3}" + line[20:]
                    )
                    metal_count += 1
                elif bb_ligand_res:
                    # Rename ligand residue (keep RF3 atom names)
                    line = line[:17] + f"{bb_ligand_res:>3}" + line[20:]
                    ligand_count += 1
        fixed_lines.append(line)

    print(f"[_fix_rf3_hetatm] Renamed {metal_count} metal + {ligand_count} ligand HETATM "
          f"(metal={bb_metal_res}, ligand={bb_ligand_res})")
    return "\n".join(fixed_lines)


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


def _merge_rf3_confidences(
    flat_metrics: Dict[str, Any],
    rf3_confidences: Optional[Dict[str, Any]] = None,
    iptm: Optional[float] = None,
    min_chain_pair_pae: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Merge RF3/AF3 ligand-interface confidence metrics into flat metrics.

    These metrics are critical for filtering metal-ligand designs — the RFD3 paper
    (bioRxiv 2025.09.18.676967v2) defines AF3 success as:
      iPTM > 0.8, min chain-pair PAE < 1.5
    Without these, designs with well-folded protein cores but disordered
    metal-ligand pockets pass filters incorrectly.

    Accepts either an rf3_confidences dict (from RF3 inference output) or
    explicit iptm/min_chain_pair_pae values. Explicit values take precedence.
    """
    if rf3_confidences:
        if "iptm" in rf3_confidences and rf3_confidences["iptm"] is not None:
            flat_metrics.setdefault("iptm", rf3_confidences["iptm"])
        # Compute min chain-pair PAE from PAE matrix if available
        # (fallback — preferred path is pre-computed in extract_rf3_confidences)
        if "pae_matrix" in rf3_confidences and rf3_confidences["pae_matrix"]:
            prot_len = rf3_confidences.get("protein_length", 0)
            pae_val = _compute_min_chain_pair_pae(rf3_confidences["pae_matrix"], prot_len)
            if pae_val is not None:
                flat_metrics.setdefault("min_chain_pair_pae", pae_val)
        # Direct min_chain_pair_pae from caller (e.g., AF3 server result)
        if "min_chain_pair_pae" in rf3_confidences and rf3_confidences["min_chain_pair_pae"] is not None:
            flat_metrics.setdefault("min_chain_pair_pae", rf3_confidences["min_chain_pair_pae"])

    # Explicit values override rf3_confidences
    if iptm is not None:
        flat_metrics["iptm"] = iptm
    if min_chain_pair_pae is not None:
        flat_metrics["min_chain_pair_pae"] = min_chain_pair_pae

    return flat_metrics


def _compute_min_chain_pair_pae(
    pae_matrix: List[List[float]],
    protein_length: int = 0,
) -> Optional[float]:
    """
    Compute mean inter-chain PAE from the PAE matrix.

    For protein-ligand complexes, the PAE matrix is (N_prot + N_lig) × (N_prot + N_lig).
    The inter-chain blocks (protein→ligand + ligand→protein) give the protein↔ligand PAE.
    Returns the mean of the inter-chain block values.

    The RFD3 paper (bioRxiv 2025.09.18.676967v2) uses min chain-pair PAE < 1.5
    as a hard gate for ligand binding success.

    Note: This is a fallback for when min_chain_pair_pae wasn't pre-computed
    in extract_rf3_confidences (e.g., from a truncated PAE matrix). The
    preferred path is to compute it from the full matrix in inference_utils.py.

    Args:
        pae_matrix: PAE matrix (may be truncated to 100×100).
        protein_length: Number of protein residues. When > 0, extracts the
            inter-chain off-diagonal block. When 0, returns overall mean PAE
            as a rough estimate.
    """
    if not pae_matrix:
        return None
    try:
        n_total = len(pae_matrix)
        if protein_length > 0 and n_total > protein_length:
            # Compute mean of inter-chain block (protein↔ligand)
            inter_chain_values = []
            for i in range(min(protein_length, n_total)):
                row = pae_matrix[i]
                for j in range(protein_length, min(len(row), n_total)):
                    inter_chain_values.append(float(row[j]))
            for i in range(protein_length, n_total):
                row = pae_matrix[i]
                for j in range(min(protein_length, len(row))):
                    inter_chain_values.append(float(row[j]))
            if inter_chain_values:
                return round(sum(inter_chain_values) / len(inter_chain_values), 3)
        else:
            # No chain boundary info — return overall mean PAE as rough estimate.
            # This is less accurate but better than nothing.
            all_values = []
            for row in pae_matrix:
                for val in row:
                    all_values.append(float(val))
            if all_values:
                return round(sum(all_values) / len(all_values), 3)
        return None
    except (TypeError, ValueError):
        return None


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
        backbone_pdb: Optional backbone PDB with HETATM records for metal/ligand.
                      When provided, HETATM records are grafted onto pdb_content
                      so that coordination analysis can detect the metal site.
                      Required when pdb_content comes from RF3 (protein-only output).
        rf3_confidences: Optional dict with RF3/AF3 confidence outputs (iptm, pae_matrix,
                         min_chain_pair_pae). Used for ligand-interface quality filtering.
        iptm: Optional float - explicit iPTM value (overrides rf3_confidences).
        min_chain_pair_pae: Optional float - explicit min chain-pair PAE (overrides rf3_confidences).

    Note on AF3 paper thresholds (bioRxiv 2025.09.18.676967v2):
        The RFD3 paper defines AF3 success for ligand-binding designs as:
        iPTM > 0.8 AND min chain-pair PAE < 1.5. These are now checked by
        metal_binding filter presets. Callers should pass rf3_confidences or
        explicit iptm/min_chain_pair_pae values for proper filtering.

    Returns:
        Enriched analysis with flat metrics, filter pass/fail, and full nested data.
    """
    try:
        from unified_analyzer import UnifiedDesignAnalyzer

        pdb_content = job_input.get("pdb_content")
        if not pdb_content:
            return {"status": "failed", "error": "pdb_content is required"}

        backbone_pdb = job_input.get("backbone_pdb")

        # Step 1: Auto-convert CIF to PDB format (RF3 outputs CIF despite .pdb extension)
        if pdb_content.lstrip().startswith("data_"):
            try:
                from biotite.structure.io import pdbx, pdb as pdb_io
                import io as _io
                cif_file = pdbx.CIFFile.read(_io.StringIO(pdb_content))
                try:
                    atom_array = pdbx.get_structure(cif_file, model=1)
                except Exception:
                    atom_array = pdbx.get_structure(cif_file)
                pdb_file = pdb_io.PDBFile()
                pdb_io.set_structure(pdb_file, atom_array)
                pdb_str = _io.StringIO()
                pdb_file.write(pdb_str)
                pdb_content = pdb_str.getvalue()
                print("[analyze_design] Converted CIF input to PDB format")
            except Exception as e:
                print(f"[analyze_design] Warning: CIF-to-PDB conversion failed: {e}")

        # Step 2a: Compute ligand RMSD BEFORE HETATM grafting.
        # After grafting, ligand positions come from the backbone template (RMSD = 0).
        # Must compare RF3-predicted HETATM (if any) to template HETATM.
        # RFD3 paper (bioRxiv 2025.09.18.676967v2): ligand RMSD < 5 Å for AF3 success.
        ligand_rmsd = None
        if backbone_pdb:
            ligand_rmsd = _compute_ligand_rmsd(
                pdb_content, backbone_pdb, metal_type=job_input.get("metal_type")
            )

        # Step 2b: Fix RF3 HETATM naming (must run AFTER CIF conversion)
        # RF3 uses non-standard names (L:0 for metal, L:1 for ligand).
        # The backbone PDB has the correct CCD residue names as reference.
        if backbone_pdb:
            pdb_content = _fix_rf3_hetatm(
                pdb_content, backbone_pdb, metal_type=job_input.get("metal_type")
            )

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

        # 2b. Add ligand RMSD if computed (RFD3 paper: < 5 Å for AF3 success)
        if ligand_rmsd is not None:
            flat_metrics["ligand_rmsd"] = ligand_rmsd

        # 2c. Merge RF3/AF3 ligand-interface confidence (iPTM, chain-pair PAE)
        # These are critical for metal-ligand designs — without them, designs with
        # well-folded cores but disordered ligand pockets pass incorrectly.
        # See AF3 paper thresholds: iPTM > 0.8, min chain-pair PAE < 1.5
        flat_metrics = _merge_rf3_confidences(
            flat_metrics,
            rf3_confidences=job_input.get("rf3_confidences"),
            iptm=job_input.get("iptm"),
            min_chain_pair_pae=job_input.get("min_chain_pair_pae"),
        )

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
