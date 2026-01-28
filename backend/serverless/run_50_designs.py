"""
Run 50 designs with the enhanced AI infrastructure and check quality.

This script runs the full AI design pipeline with the original query
that prompted the infrastructure enhancement.

Can run in two modes:
1. API mode (default): Calls Docker API via HTTP (for running from Windows)
2. Mock mode: Uses mock data for testing without API
"""
import sys
import time
import json
import requests
from datetime import datetime
from pathlib import Path

# Test query
TEST_QUERY = "scaffold the active pocket pqq-ca of 4cvb and keep all interacting residues and make protein more stable, and i want it able to preserve its ADHD functionality."


def _merge_hetatm_to_backbone(
    backbone_pdb: str,
    scaffold_pdb: str,
    blocker_residues: list = None,
) -> str:
    """
    Merge HETATM records from scaffold into backbone PDB.

    RFD3 only outputs protein backbone (ATOM records). The metal and ligand
    (HETATM) from the input scaffold need to be merged back.

    Also removes blocker ALA residues that were only needed during RFD3.

    Args:
        backbone_pdb: RFD3 output (ATOM records only)
        scaffold_pdb: Original scaffold with HETATM (metal M chain, ligand L chain)
        blocker_residues: List of blocker residue IDs to remove (e.g., ["A22", "A23"])

    Returns:
        PDB with both ATOM and HETATM records (blockers removed)
    """
    # Parse blocker residue numbers
    blocker_resnums = set()
    if blocker_residues:
        for res_id in blocker_residues:
            try:
                resnum = int(res_id[1:])  # Strip chain letter, get number
                blocker_resnums.add(resnum)
            except (ValueError, IndexError):
                pass

    # Extract ATOM lines from backbone, filtering out blockers
    atom_lines = []
    for line in backbone_pdb.split('\n'):
        if line.startswith('ATOM'):
            # Check if this is a blocker residue
            if blocker_resnums:
                try:
                    resnum = int(line[22:26].strip())
                    if resnum in blocker_resnums:
                        continue  # Skip blocker
                except ValueError:
                    pass
            atom_lines.append(line)
        elif line.startswith('TER'):
            atom_lines.append(line)

    # Extract HETATM from scaffold
    hetatm_lines = []
    for line in scaffold_pdb.split('\n'):
        if line.startswith('HETATM'):
            hetatm_lines.append(line)

    if not hetatm_lines:
        print(f"  [Warning] No HETATM found in scaffold")
        return backbone_pdb

    # Find max atom serial in backbone
    max_atom_num = 0
    for line in atom_lines:
        if line.startswith('ATOM'):
            try:
                atom_num = int(line[6:11].strip())
                max_atom_num = max(max_atom_num, atom_num)
            except ValueError:
                pass

    # Renumber HETATM to continue from max atom serial
    renumbered_hetatm = []
    for i, line in enumerate(hetatm_lines):
        new_atom_num = max_atom_num + i + 1
        new_line = f"{line[:6]}{new_atom_num:5d}{line[11:]}"
        renumbered_hetatm.append(new_line)

    # Combine ATOM + renumbered HETATM
    result_pdb = '\n'.join(atom_lines) + '\n' + '\n'.join(renumbered_hetatm) + '\nEND\n'

    return result_pdb


def _check_ligand_clashes(
    pdb_content: str,
    ligand_name: str,
    clash_threshold: float = 2.0,
) -> dict:
    """
    Check for clashes between designed protein and ligand atoms.

    Uses a tiered approach:
    - Hard clashes (< 1.2 Å): Always flagged — atoms overlap
    - Polar-polar contacts (1.2-2.0 Å between O/N/S atoms): Tolerated as
      potential H-bond contacts that FastRelax can optimize
    - Non-polar clashes (1.2-threshold Å): Flagged — no physical basis

    Args:
        pdb_content: PDB content with protein and ligand HETATM
        ligand_name: 3-letter ligand code (e.g., "PQQ")
        clash_threshold: Distance threshold for non-polar clashes (default 2.0 Å)

    Returns:
        Dict with has_clashes, clash_count, min_distance, soft_contacts
    """
    import math

    POLAR_ELEMENTS = {'O', 'N', 'S'}
    HARD_CLASH_THRESHOLD = 1.2  # Atoms overlapping — reject always

    def _is_polar(atom_name: str) -> bool:
        """Check if atom is polar based on first non-digit character."""
        for ch in atom_name:
            if ch.isalpha():
                return ch in POLAR_ELEMENTS
        return False

    # Parse ligand atoms
    ligand_atoms = []
    for line in pdb_content.split('\n'):
        if line.startswith('HETATM') and line[17:20].strip().upper() == ligand_name.upper():
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atomname = line[12:16].strip()
                ligand_atoms.append({'name': atomname, 'coord': (x, y, z)})
            except (ValueError, IndexError):
                pass

    if not ligand_atoms:
        return {'has_clashes': False, 'clash_count': 0, 'min_distance': float('inf'),
                'soft_contacts': 0}

    # Parse protein atoms (all chains, including A)
    protein_atoms = []
    for line in pdb_content.split('\n'):
        if line.startswith('ATOM'):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atomname = line[12:16].strip()
                resname = line[17:20].strip()
                resnum = line[22:26].strip()
                protein_atoms.append({'name': atomname, 'resname': resname, 'resnum': resnum, 'coord': (x, y, z)})
            except (ValueError, IndexError):
                pass

    if not protein_atoms:
        return {'has_clashes': False, 'clash_count': 0, 'min_distance': float('inf'),
                'soft_contacts': 0}

    # Check distances
    clashes = []
    soft_contacts = []
    min_distance = float('inf')

    for pa in protein_atoms:
        px, py, pz = pa['coord']
        for la in ligand_atoms:
            lx, ly, lz = la['coord']
            dist = math.sqrt((px - lx)**2 + (py - ly)**2 + (pz - lz)**2)
            min_distance = min(min_distance, dist)

            if dist >= clash_threshold:
                continue

            contact = {
                'protein': f"{pa['resnum']}:{pa['resname']}.{pa['name']}",
                'ligand': f"{ligand_name}.{la['name']}",
                'distance': dist,
            }

            # Hard clash: atoms overlapping
            if dist < HARD_CLASH_THRESHOLD:
                clashes.append(contact)
            # Polar-polar: potential H-bond contact, tolerate
            elif _is_polar(pa['name']) and _is_polar(la['name']):
                soft_contacts.append(contact)
            # Non-polar or mixed: steric clash
            else:
                clashes.append(contact)

    return {
        'has_clashes': len(clashes) > 0,
        'clash_count': len(clashes),
        'clashes': clashes[:5],  # Keep first 5 for debugging
        'soft_contacts': len(soft_contacts),
        'soft_contact_details': soft_contacts[:5],
        'min_distance': min_distance,
    }


###############################################################################
# VALIDATION HELPERS
###############################################################################

# PQQ canonical SMILES (PubChem CID 1024)
PQQ_SMILES = "C1=C(C2=C(C(=O)C(=O)C3=C2NC(=C3)C(=O)O)N=C1C(=O)O)C(=O)O"

BASE_URL = "http://localhost:8000"


# Import shared utilities from esmfold_utils (avoids code duplication)
try:
    from esmfold_utils import (
        _kabsch_rmsd as _kabsch_rmsd_shared,
        _extract_ca_coords_from_pdb as _extract_ca_pdb_shared,
        _extract_ca_coords_from_cif as _extract_ca_cif_shared,
    )

    def _kabsch_rmsd_numpy(coords1, coords2):
        """Wrapper: delegates to esmfold_utils._kabsch_rmsd."""
        import numpy as np
        c1 = np.array(coords1, dtype=np.float64)
        c2 = np.array(coords2, dtype=np.float64)
        n = min(len(c1), len(c2))
        if n == 0:
            return float("inf")
        return float(_kabsch_rmsd_shared(c1[:n], c2[:n]))

    def _parse_ca_from_pdb(pdb_text: str):
        """Wrapper: delegates to esmfold_utils._extract_ca_coords_from_pdb."""
        coords = _extract_ca_pdb_shared(pdb_text)
        return [tuple(row) for row in coords.tolist()] if len(coords) > 0 else []

    def _parse_ca_from_cif(cif_text: str):
        """Wrapper: delegates to esmfold_utils._extract_ca_coords_from_cif."""
        coords = _extract_ca_cif_shared(cif_text)
        return [tuple(row) for row in coords.tolist()] if len(coords) > 0 else []

except ImportError:
    # Standalone fallback — keep local definitions for out-of-container use
    def _kabsch_rmsd_numpy(coords1, coords2):
        """Proper Kabsch RMSD using numpy SVD."""
        import numpy as np
        n = min(len(coords1), len(coords2))
        if n == 0:
            return float("inf")
        P = np.array(coords1[:n], dtype=np.float64)
        Q = np.array(coords2[:n], dtype=np.float64)
        P_centered = P - P.mean(axis=0)
        Q_centered = Q - Q.mean(axis=0)
        H = P_centered.T @ Q_centered
        U, S, Vt = np.linalg.svd(H)
        d = np.sign(np.linalg.det(Vt.T @ U.T))
        D = np.diag([1.0, 1.0, d])
        R = Vt.T @ D @ U.T
        Q_aligned = Q_centered @ R
        diff = P_centered - Q_aligned
        return float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))

    def _parse_ca_from_pdb(pdb_text: str):
        """Extract CA atom coordinates from PDB text."""
        cas = []
        for line in pdb_text.splitlines():
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    cas.append((x, y, z))
                except (ValueError, IndexError):
                    pass
        return cas

    def _parse_ca_from_cif(cif_text: str):
        """Extract CA atom coordinates from mmCIF text."""
        cas = []
        in_atom_site = False
        col_map = {}
        col_idx = 0
        for line in cif_text.splitlines():
            if line.startswith("_atom_site."):
                in_atom_site = True
                field = line.strip().split(".")[1].strip()
                col_map[field] = col_idx
                col_idx += 1
                continue
            if in_atom_site and (line.startswith("_") or line.startswith("#")):
                in_atom_site = False
                col_map = {}
                col_idx = 0
                continue
            if in_atom_site and line.strip() and not line.startswith("_") and not line.startswith("#"):
                if not col_map:
                    continue
                parts = line.split()
                if len(parts) < max(col_map.values()) + 1:
                    continue
                label_atom_id = col_map.get("label_atom_id") or col_map.get("auth_atom_id")
                if label_atom_id is None:
                    continue
                atom_name = parts[label_atom_id].strip('"')
                if atom_name != "CA":
                    continue
                group_idx = col_map.get("group_PDB")
                if group_idx is not None and parts[group_idx] != "ATOM":
                    continue
                x_idx = col_map.get("Cartn_x")
                y_idx = col_map.get("Cartn_y")
                z_idx = col_map.get("Cartn_z")
                if x_idx is None or y_idx is None or z_idx is None:
                    continue
                try:
                    cas.append((float(parts[x_idx]), float(parts[y_idx]), float(parts[z_idx])))
                except ValueError:
                    continue
        return cas


def run_validate_design_api(pdb_content, sequences=None, max_backbone_rmsd=2.0,
                            min_plddt=0.70, num_sequences=8):
    """Call validate_design handler for proper Kabsch RMSD validation (ESMFold).

    Args:
        pdb_content: Designed backbone PDB content
        sequences: Optional list of sequences (if None, MPNN generates them)
        max_backbone_rmsd: Pass threshold for backbone RMSD
        min_plddt: Pass threshold for pLDDT
        num_sequences: Number of sequences to generate if sequences not provided

    Returns:
        Dict with validation results including per-sequence RMSD and pLDDT
    """
    import requests

    payload = {
        "task": "validate_design",
        "pdb_content": pdb_content,
        "num_sequences": num_sequences,
        "max_backbone_rmsd": max_backbone_rmsd,
        "min_plddt": min_plddt,
        "use_ligandmpnn": True,
        "temperature": 0.1,
    }

    try:
        resp = requests.post(
            f"{BASE_URL}/runsync",
            json={"input": payload},
            timeout=600,
        )
        raw = resp.json()
        if raw.get("status") == "COMPLETED":
            return raw.get("output", {})
        return {"error": f"status={raw.get('status')}", "raw": raw}
    except Exception as e:
        return {"error": str(e)}


def run_rf3_prediction(sequence, name, ligand_smiles=None):
    """Call RF3 structure prediction API.

    Args:
        sequence: Amino acid sequence
        name: Name for the prediction
        ligand_smiles: Optional SMILES string for ligand-aware prediction

    Returns:
        Dict with mean_plddt, ptm, cif_content, or error
    """
    import requests

    payload = {"task": "rf3", "sequence": sequence, "name": name}
    if ligand_smiles:
        payload["ligand_smiles"] = ligand_smiles

    try:
        resp = requests.post(
            f"{BASE_URL}/runsync",
            json={"input": payload},
            timeout=600,
        )
        raw = resp.json()
    except Exception as e:
        return {"error": str(e)}

    if raw.get("status") != "COMPLETED":
        return {"error": f"status={raw.get('status')}"}

    output = raw.get("output", {})
    if output.get("status") not in ("success", "completed"):
        return {"error": output.get("error", f"status={output.get('status')}")}

    result = output.get("result", {})
    confidences = result.get("confidences", {})
    predictions = result.get("predictions", [])
    cif_content = predictions[0].get("content", "") if predictions else ""

    return {
        "mean_plddt": confidences.get("mean_plddt"),
        "ptm": confidences.get("ptm"),
        "overall_plddt": confidences.get("overall_plddt"),
        "overall_pae": confidences.get("overall_pae"),
        "ranking_score": confidences.get("ranking_score"),
        "cif_content": cif_content,
    }


def run_rf3_validation(sequence, backbone_pdb, name, ligand_smiles=None):
    """RF3 structure prediction + Kabsch RMSD vs designed backbone.

    Args:
        sequence: Amino acid sequence
        backbone_pdb: Designed backbone PDB content (reference)
        name: Name for the prediction
        ligand_smiles: Optional SMILES for ligand-aware prediction

    Returns:
        Dict with rf3_plddt, rf3_ptm, rf3_rmsd, and pass/fail status
    """
    pred = run_rf3_prediction(sequence, name, ligand_smiles=ligand_smiles)
    if "error" in pred:
        return {"error": pred["error"], "rf3_plddt": None, "rf3_ptm": None, "rf3_rmsd": None}

    plddt = pred.get("mean_plddt")
    ptm = pred.get("ptm")

    # Compute Kabsch RMSD vs designed backbone
    rmsd = None
    if pred.get("cif_content"):
        pred_cas = _parse_ca_from_cif(pred["cif_content"])
        ref_cas = _parse_ca_from_pdb(backbone_pdb)
        if pred_cas and ref_cas:
            rmsd = _kabsch_rmsd_numpy(ref_cas, pred_cas)

    return {
        "rf3_plddt": plddt,
        "rf3_ptm": ptm,
        "rf3_rmsd": rmsd,
        "cif_content": pred.get("cif_content", ""),
    }


def run_validation_stage(backbone_pdbs, all_sequences, num_sequences_per_design,
                         ligand_smiles=None, output_dir=None):
    """Run ESMFold + RF3 validation on all sequences.

    Args:
        backbone_pdbs: List of backbone PDB strings
        all_sequences: List of amino acid sequences (flat list)
        num_sequences_per_design: Number of sequences per backbone
        ligand_smiles: Optional SMILES for ligand-aware RF3 prediction
        output_dir: Directory to save results

    Returns:
        List of validation result dicts
    """
    print("\n" + "=" * 70)
    print("VALIDATION STAGE")
    print("=" * 70)

    validation_results = []

    # Test PQQ SMILES with RF3 if provided
    pqq_works = False
    if ligand_smiles and all_sequences:
        print("\nTesting ligand SMILES with RF3...")
        test_pred = run_rf3_prediction(all_sequences[0], "smiles_test", ligand_smiles=ligand_smiles)
        if "error" not in test_pred:
            print(f"  Ligand SMILES works! pLDDT={test_pred['mean_plddt']:.3f}, pTM={test_pred['ptm']:.3f}")
            pqq_works = True
        else:
            print(f"  Ligand SMILES failed: {test_pred['error']} — running protein-only RF3")

    # Validate each sequence
    total = len(all_sequences)
    for idx, seq in enumerate(all_sequences):
        backbone_idx = idx // num_sequences_per_design
        seq_idx = idx % num_sequences_per_design
        name = f"bb{backbone_idx+1}_seq{seq_idx+1}"

        if backbone_idx >= len(backbone_pdbs):
            print(f"  [WARN] No backbone for sequence {idx}, skipping")
            continue

        backbone_pdb = backbone_pdbs[backbone_idx]

        print(f"\n  [{idx+1}/{total}] {name} ({len(seq)} res)...", end=" ")

        result = {
            "backbone": backbone_idx + 1,
            "seq": seq_idx + 1,
            "sequence_length": len(seq),
            "sequence": seq,
        }

        # RF3 protein-only validation
        rf3_result = run_rf3_validation(seq, backbone_pdb, name)
        if "error" in rf3_result:
            print(f"RF3 ERROR: {rf3_result['error']}")
            result.update({"rf3_plddt": None, "rf3_ptm": None, "rf3_rmsd": None, "rf3_error": rf3_result["error"]})
        else:
            plddt = rf3_result["rf3_plddt"]
            ptm = rf3_result["rf3_ptm"]
            rmsd = rf3_result["rf3_rmsd"]
            print(f"pLDDT={plddt:.3f}, pTM={ptm:.3f}", end="")
            if rmsd is not None:
                print(f", RMSD={rmsd:.2f}A", end="")
            print()
            result.update({
                "rf3_plddt": plddt,
                "rf3_ptm": ptm,
                "rf3_rmsd": rmsd,
            })

        # RF3 + ligand validation (if SMILES works)
        if pqq_works:
            rf3_lig = run_rf3_validation(seq, backbone_pdb, f"{name}_lig", ligand_smiles=ligand_smiles)
            if "error" not in rf3_lig:
                result["rf3_lig_plddt"] = rf3_lig["rf3_plddt"]
                result["rf3_lig_ptm"] = rf3_lig["rf3_ptm"]
                result["rf3_lig_rmsd"] = rf3_lig["rf3_rmsd"]
                print(f"    +Ligand: pLDDT={rf3_lig['rf3_plddt']:.3f}, pTM={rf3_lig['rf3_ptm']:.3f}", end="")
                if rf3_lig["rf3_rmsd"] is not None:
                    print(f", RMSD={rf3_lig['rf3_rmsd']:.2f}A", end="")
                print()
            else:
                result["rf3_lig_plddt"] = None
                result["rf3_lig_ptm"] = None
                result["rf3_lig_rmsd"] = None

        # Determine pass/fail
        passed_rf3 = (
            result.get("rf3_rmsd") is not None
            and result["rf3_rmsd"] < 2.0
            and result.get("rf3_plddt") is not None
            and result["rf3_plddt"] > 0.70
        )
        passed_rf3_lig = (
            result.get("rf3_lig_ptm") is not None
            and result["rf3_lig_ptm"] > 0.80
            and result.get("rf3_lig_plddt") is not None
            and result["rf3_lig_plddt"] > 0.75
        )
        result["passed_rf3"] = passed_rf3
        result["passed_rf3_ligand"] = passed_rf3_lig
        result["passed_any"] = passed_rf3 or passed_rf3_lig

        validation_results.append(result)

    # Summary
    print("\n" + "-" * 70)
    print("VALIDATION SUMMARY")
    print("-" * 70)
    n_total = len(validation_results)
    n_passed_rf3 = sum(1 for r in validation_results if r.get("passed_rf3"))
    n_passed_lig = sum(1 for r in validation_results if r.get("passed_rf3_ligand"))
    n_passed_any = sum(1 for r in validation_results if r.get("passed_any"))

    rmsds = [r["rf3_rmsd"] for r in validation_results if r.get("rf3_rmsd") is not None]
    plddts = [r["rf3_plddt"] for r in validation_results if r.get("rf3_plddt") is not None]
    ptms = [r["rf3_ptm"] for r in validation_results if r.get("rf3_ptm") is not None]

    print(f"\nTotal sequences validated: {n_total}")
    print(f"Passed RF3 (RMSD<2A + pLDDT>0.7): {n_passed_rf3}/{n_total} ({100*n_passed_rf3/max(n_total,1):.1f}%)")
    if pqq_works:
        print(f"Passed RF3+Ligand (pTM>0.8 + pLDDT>0.75): {n_passed_lig}/{n_total} ({100*n_passed_lig/max(n_total,1):.1f}%)")
    print(f"Passed ANY criterion: {n_passed_any}/{n_total} ({100*n_passed_any/max(n_total,1):.1f}%)")

    if rmsds:
        print(f"\nRF3 RMSD: min={min(rmsds):.2f}A, max={max(rmsds):.2f}A, mean={sum(rmsds)/len(rmsds):.2f}A")
    if plddts:
        print(f"RF3 pLDDT: min={min(plddts):.3f}, max={max(plddts):.3f}, mean={sum(plddts)/len(plddts):.3f}")
    if ptms:
        print(f"RF3 pTM: min={min(ptms):.3f}, max={max(ptms):.3f}, mean={sum(ptms)/len(ptms):.3f}")

    # Per-sequence table
    print(f"\n{'BB':>3} {'Seq':>3} | {'pLDDT':>6} {'pTM':>5} {'RMSD':>7} {'Pass':>4}", end="")
    if pqq_works:
        print(f" | {'Lig_pLDDT':>9} {'Lig_pTM':>7} {'Lig_RMSD':>8}")
    else:
        print()
    for r in validation_results:
        plddt_s = f"{r['rf3_plddt']:.3f}" if r.get("rf3_plddt") is not None else "ERR"
        ptm_s = f"{r['rf3_ptm']:.3f}" if r.get("rf3_ptm") is not None else "ERR"
        rmsd_s = f"{r['rf3_rmsd']:.2f}A" if r.get("rf3_rmsd") is not None else "N/A"
        pass_s = "Y" if r.get("passed_any") else "N"
        print(f"{r['backbone']:>3} {r['seq']:>3} | {plddt_s:>6} {ptm_s:>5} {rmsd_s:>7} {pass_s:>4}", end="")
        if pqq_works:
            lp = f"{r['rf3_lig_plddt']:.3f}" if r.get("rf3_lig_plddt") is not None else "N/A"
            lt = f"{r['rf3_lig_ptm']:.3f}" if r.get("rf3_lig_ptm") is not None else "N/A"
            lr = f"{r['rf3_lig_rmsd']:.2f}A" if r.get("rf3_lig_rmsd") is not None else "N/A"
            print(f" | {lp:>9} {lt:>7} {lr:>8}")
        else:
            print()

    # Save results
    if output_dir:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        val_file = output_path / "validation_results.json"
        # Strip sequence text for smaller file (keep first 20 chars)
        save_results = []
        for r in validation_results:
            sr = dict(r)
            if "sequence" in sr:
                sr["sequence"] = sr["sequence"][:20] + "..."
            if "cif_content" in sr:
                del sr["cif_content"]
            save_results.append(sr)
        with open(val_file, 'w') as f:
            json.dump(save_results, f, indent=2)
        print(f"\nValidation results saved to {val_file}")

    return validation_results


def run_designs_via_api(num_designs: int = 50, num_sequences_per_design: int = 4):
    """
    Run designs using the HTTP API (for calling from outside Docker).
    """
    import asyncio

    async def _run():
        print("=" * 70)
        print(f"AI DESIGN PIPELINE - {num_designs} DESIGNS [API MODE]")
        print("=" * 70)
        print(f"\nQuery: {TEST_QUERY}")
        print(f"\nTarget: {num_designs} backbone designs, {num_sequences_per_design} sequences each")
        print(f"Total sequences: {num_designs * num_sequences_per_design}")
        print("=" * 70)

        # Import components
        from nl_design_parser import SimpleFallbackParser
        from design_rules import select_stability_profile
        from scaffolding_workflow import ScaffoldingWorkflow, check_backbone_continuity
        from api_client import run_rfd3_api, run_mpnn_api, check_health

        # Check API health
        print("\nChecking API health...")
        if not check_health():
            print("[FAIL] API is not available. Start Docker container first.")
            return None

        print("[OK] API is healthy")

        # Parse query
        print("\nParsing query...")
        parser = SimpleFallbackParser()
        intent = parser.parse(TEST_QUERY)

        print(f"  - PDB: {intent.source_pdb_id}")
        print(f"  - Ligand: {intent.ligand_name}")
        print(f"  - Metal: {intent.metal_type}")
        print(f"  - Enzyme: {intent.enzyme_class}")
        print(f"  - Stability: {intent.stability_focus}")

        # Get stability profile
        stability = select_stability_profile(intent, intent.design_goal or "binding")
        print(f"  - Profile: {stability.name} (cfg={stability.cfg_scale}, step={stability.step_scale})")

        # Decide fixed_atom_type based on user intent:
        #   preserve_function=True → ALL (preserve exact catalytic geometry)
        #   include_all_contacts only → TIP (keep functional groups, allow adaptation)
        #   stability_focus only     → BKBN (allow sidechain redesign for packing)
        #   default                  → BKBN
        if intent.preserve_function:
            fixed_atom_type = "ALL"
            print(f"\n  [Intent] preserve_function=True -> fixing ALL atoms on catalytic residues")
        elif intent.include_all_contacts and not intent.stability_focus:
            fixed_atom_type = "TIP"
            print(f"\n  [Intent] include_all_contacts=True -> fixing TIP atoms on catalytic residues")
        else:
            fixed_atom_type = "BKBN"
            print(f"\n  [Intent] default -> fixing BKBN atoms on catalytic residues")

        # Extract scaffold (async)
        print("\nExtracting scaffold from PDB...")
        workflow = ScaffoldingWorkflow()
        scaffold_result = await workflow.run(
            pdb_id=intent.source_pdb_id,
            ligand_code=intent.ligand_name,
            metal=intent.metal_type,
            include_all_ligand_contacts=intent.include_all_contacts,
            fixed_atom_type=fixed_atom_type,
            use_minimal_motif=True,          # Evidence-based residue selection
            enzyme_class=intent.enzyme_class, # For enzyme chemistry evidence
        )

        if not scaffold_result or not scaffold_result.success:
            print(f"[FAIL] Could not extract scaffold: {getattr(scaffold_result, 'error', 'unknown')}")
            return None

        print(f"  - Design approach: {scaffold_result.source_info.get('design_approach', 'unknown')}")
        print(f"  - Catalytic residues: {scaffold_result.source_info.get('catalytic_residue_ids', [])[:5]}...")
        print(f"  - Fixed atoms: {len(scaffold_result.fixed_atoms)} entries")
        print(f"  - Coordination: {scaffold_result.source_info.get('coordination_number', 'N/A')}")
        print(f"  - Minimal motif: {scaffold_result.source_info.get('use_minimal_motif', False)}")

        # Print minimal motif details if available
        if scaffold_result.motif_result:
            mr = scaffold_result.motif_result
            tier_counts = {}
            for r in mr.motif_residues:
                tier_counts[r.tier] = tier_counts.get(r.tier, 0) + 1
            print(f"\n  Minimal Motif Selection:")
            print(f"    - Motif residues: {len(mr.motif_residues)} (Tier 1-3)")
            print(f"    - Context residues: {len(mr.context_residues)} (Tier 4, excluded)")
            print(f"    - Total fixed atoms: {mr.total_fixed_atoms}")
            print(f"    - Tier breakdown: {tier_counts}")
            print(f"    - Evidence sources: {mr.evidence_summary}")
            for r in mr.motif_residues:
                atoms_brief = ",".join(r.functional_atoms[:4])
                if len(r.functional_atoms) > 4:
                    atoms_brief += f"..({len(r.functional_atoms)})"
                print(f"      T{r.tier} {r.chain}{r.resnum:>4d} {r.resname:3s} "
                      f"fix={r.fix_type:4s} atoms=[{atoms_brief}]")
            if mr.warnings:
                for w in mr.warnings:
                    print(f"    [WARN] {w}")

        # Build RFD3 params using ENZYME SCAFFOLD DESIGN approach
        # This uses length + unindex + ligand (NOT contig with blockers)
        # This is the approach that works correctly in the UI

        print(f"\n  Using Enzyme Scaffold Design approach:")
        print(f"    - Length: {scaffold_result.length}")
        print(f"    - Unindex: {scaffold_result.unindex}")
        print(f"    - Ligand codes: {scaffold_result.ligand_codes}")

        # NOTE: For Enzyme Scaffold Design, use balanced parameters
        # The "ultra_stable" profile with step_scale=1.0, gamma_0=0.35, etc.
        # can cause more clashes. Use defaults that work reliably.
        rfd3_params = {
            # Enzyme Scaffold Design parameters (correct approach)
            "length": scaffold_result.length,
            "unindex": scaffold_result.unindex,
            "ligand": scaffold_result.ligand_codes,
            "pdb_content": scaffold_result.motif_pdb,
            "select_fixed_atoms": scaffold_result.fixed_atoms,
            "select_buried": scaffold_result.rasa_targets if scaffold_result.rasa_targets else None,
            # Design settings - use balanced defaults for scaffold design
            "num_designs": num_designs,
            "is_non_loopy": True,
            "cfg_scale": 2.5,  # Balanced, not ultra_stable
            # Don't override step_scale, num_timesteps, gamma_0 - use RFD3 defaults
            "use_classifier_free_guidance": True,
        }
        print(f"    - Using balanced parameters (cfg_scale=2.5)")

        # H-bond conditioning: tell RFD3 that ligand O atoms are H-bond acceptors
        # and N atoms are donors. This guides backbone placement away from clashes.
        if scaffold_result.hbond_acceptors:
            rfd3_params["select_hbond_acceptor"] = scaffold_result.hbond_acceptors
        if scaffold_result.hbond_donors:
            rfd3_params["select_hbond_donor"] = scaffold_result.hbond_donors

        # Run RFD3
        print(f"\nRunning RFD3 ({num_designs} designs)...")
        start_time = time.time()
        rfd3_result = run_rfd3_api(**rfd3_params, timeout=600)
        rfd3_time = time.time() - start_time

        if rfd3_result["status"] != "completed":
            print(f"[FAIL] RFD3 failed: {rfd3_result.get('error')}")
            return None

        designs = rfd3_result["result"].get("designs", [])
        raw_backbone_pdbs = [d.get("content", "") for d in designs if d.get("content")]
        print(f"  - Generated: {len(raw_backbone_pdbs)} backbones in {rfd3_time:.1f}s")

        # Check if HETATM records (metal + ligand) need to be merged
        # With Enzyme Scaffold Design approach, RFD3 may include HETATM in output
        print("\nChecking HETATM records in designs...")
        merged_pdbs = []
        needs_merge_count = 0
        for pdb in raw_backbone_pdbs:
            has_hetatm = any(line.startswith('HETATM') for line in pdb.split('\n'))
            if has_hetatm:
                merged_pdbs.append(pdb)  # Already has HETATM
            else:
                # Merge HETATM from scaffold (no blockers in new approach)
                merged_pdb = _merge_hetatm_to_backbone(pdb, scaffold_result.motif_pdb, [])
                merged_pdbs.append(merged_pdb)
                needs_merge_count += 1
        if needs_merge_count > 0:
            print(f"  - Merged HETATM records to {needs_merge_count} backbones")
        else:
            print(f"  - All {len(merged_pdbs)} designs already have HETATM records")

        # CRITICAL: Filter designs with ligand clashes
        # RFD3 does NOT perform collision detection - designs may pass through ligand
        print("\nFiltering ligand clashes...")
        ligand_code = intent.ligand_name.upper() if intent.ligand_name else "PQQ"
        backbone_pdbs = []
        clash_count = 0
        for i, pdb in enumerate(merged_pdbs):
            clash_result = _check_ligand_clashes(pdb, ligand_code, clash_threshold=2.0)
            if clash_result['has_clashes']:
                clash_count += 1
                if clash_count <= 3:  # Show first 3 failures
                    soft = clash_result.get('soft_contacts', 0)
                    print(f"  [CLASH] Design {i+1}: {clash_result['clash_count']} hard clashes, "
                          f"{soft} soft contacts, min dist {clash_result['min_distance']:.2f} Å")
            else:
                backbone_pdbs.append(pdb)

        print(f"  - Passed: {len(backbone_pdbs)}/{len(merged_pdbs)} ({100*len(backbone_pdbs)/len(merged_pdbs):.0f}%)")
        print(f"  - Rejected: {clash_count} designs with hard ligand clashes")

        if not backbone_pdbs:
            print("[FAIL] All designs had ligand clashes. Try different parameters.")
            return None

        # Filter backbone continuity (chain breaks)
        print("\nChecking backbone continuity...")
        continuous_pdbs = []
        for i, pdb in enumerate(backbone_pdbs):
            cont = check_backbone_continuity(pdb)
            if cont["continuous"]:
                continuous_pdbs.append(pdb)
            else:
                print(f"  [BREAK] Design {i+1}: {cont['num_breaks']} breaks "
                      f"(max CA-CA={cont['max_observed']}A)")
        if continuous_pdbs:
            print(f"  - Continuous: {len(continuous_pdbs)}/{len(backbone_pdbs)}")
        else:
            print(f"  [WARN] No continuous backbones! Using all {len(backbone_pdbs)} designs anyway.")
            continuous_pdbs = backbone_pdbs  # Fallback: use all
        backbone_pdbs = continuous_pdbs

        # Run MPNN with proper enzyme-aware settings
        # Fix catalytic residues so MPNN does NOT redesign them
        catalytic_residue_ids = scaffold_result.source_info.get("catalytic_residue_ids", [])

        # Determine bias_AA based on ligand types present
        # NOTE: bias_AA is GLOBAL - applies to ALL positions, not just binding site.
        # LigandMPNN already uses atom context to design good binding residues.
        # Keep biases MILD to avoid destroying the hydrophobic core needed for folding.
        has_metal = bool(intent.metal_type)
        has_small_molecule = bool(intent.ligand_name)
        if has_metal and has_small_molecule:
            # Mild bias: slight preference for coordinating/aromatic residues
            mpnn_bias = "H:0.5,D:0.5,E:0.5,W:0.3,Y:0.3"
            mpnn_omit = None  # Allow cysteine for metal coordination
        elif has_metal:
            mpnn_bias = "H:1.0,C:1.0,D:0.5,E:0.5"
            mpnn_omit = None
        elif has_small_molecule:
            mpnn_bias = "W:0.5,Y:0.5,F:0.3"
            mpnn_omit = "C"
        else:
            mpnn_bias = None
            mpnn_omit = "C"

        print(f"\nRunning LigandMPNN ({num_sequences_per_design} seqs/backbone)...")
        print(f"  - Fixed positions (catalytic): {catalytic_residue_ids[:5]}{'...' if len(catalytic_residue_ids) > 5 else ''} ({len(catalytic_residue_ids)} residues)")
        print(f"  - bias_AA: {mpnn_bias}")
        print(f"  - omit_AA: {mpnn_omit}")
        print(f"  - pack_side_chains: True (ligand-aware)")
        start_time = time.time()
        all_sequences = []

        for i, pdb in enumerate(backbone_pdbs):
            mpnn_result = run_mpnn_api(
                pdb_content=pdb,
                num_seqs=num_sequences_per_design,
                temperature=0.1,
                fixed_positions=catalytic_residue_ids if catalytic_residue_ids else None,
                bias_AA=mpnn_bias,
                omit_AA=mpnn_omit,
                pack_side_chains=True,
                pack_with_ligand_context=True,
                ligand_cutoff_for_score=5.0,
                use_side_chain_context=bool(catalytic_residue_ids),
                save_stats=True,
            )
            if mpnn_result["status"] == "completed":
                seqs = mpnn_result["result"].get("sequences", [])
                for seq_data in seqs:
                    content = seq_data.get("content", "")
                    for line in content.split('\n'):
                        if line and not line.startswith('>'):
                            all_sequences.append(line.strip())
            if (i + 1) % 10 == 0:
                print(f"    Processed {i + 1}/{len(backbone_pdbs)} backbones...")

        mpnn_time = time.time() - start_time
        print(f"  - Generated: {len(all_sequences)} sequences in {mpnn_time:.1f}s")

        # Step 7: Validation (RF3 + Kabsch RMSD)
        ligand_smiles = PQQ_SMILES if intent.ligand_name and intent.ligand_name.upper() == "PQQ" else None
        validation_start = time.time()
        validation_results = run_validation_stage(
            backbone_pdbs=backbone_pdbs,
            all_sequences=all_sequences,
            num_sequences_per_design=num_sequences_per_design,
            ligand_smiles=ligand_smiles,
            output_dir=output_dir if 'output_dir' in dir() else Path("design_results"),
        )
        validation_time = time.time() - validation_start

        # Results summary
        print("\n" + "=" * 70)
        print("RESULTS SUMMARY")
        print("=" * 70)
        print(f"\nBackbones: {len(backbone_pdbs)}")
        print(f"Sequences: {len(all_sequences)}")
        print(f"\nStability Profile: {stability.name}")
        print(f"  - cfg_scale: {stability.cfg_scale}")
        print(f"  - step_scale: {stability.step_scale}")
        print(f"  - gamma_0: {stability.gamma_0}")
        print(f"  - num_timesteps: {stability.num_timesteps}")
        print(f"  - plddt_threshold: {stability.plddt_threshold}")

        if intent.enzyme_class:
            print(f"\nEnzyme Preservation:")
            print(f"  - Class: {intent.enzyme_class}")
            print(f"  - Preserve Function: {intent.preserve_function}")

        print(f"\nTiming:")
        print(f"  - RFD3: {rfd3_time:.1f}s")
        print(f"  - MPNN: {mpnn_time:.1f}s")
        print(f"  - Validation: {validation_time:.1f}s")
        print(f"  - Total: {rfd3_time + mpnn_time + validation_time:.1f}s")

        # Save results
        output_dir = Path("design_results")
        output_dir.mkdir(exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # Save backbones
        for i, pdb in enumerate(backbone_pdbs):
            pdb_file = output_dir / f"backbone_{timestamp}_{i+1}.pdb"
            with open(pdb_file, 'w') as f:
                f.write(pdb)

        # Save sequences
        seq_file = output_dir / f"sequences_{timestamp}.fasta"
        with open(seq_file, 'w') as f:
            for i, seq in enumerate(all_sequences):
                f.write(f">design_{i+1}\n{seq}\n")

        print(f"\nSaved {len(backbone_pdbs)} PDBs and {len(all_sequences)} sequences to {output_dir}")

        # Create result object for compatibility
        n_passed = sum(1 for v in validation_results if v.get("passed_any")) if validation_results else 0
        n_total = len(validation_results) if validation_results else len(all_sequences)
        rmsds = [v["rf3_rmsd"] for v in validation_results if v.get("rf3_rmsd") is not None] if validation_results else []
        plddts_val = [v["rf3_plddt"] for v in validation_results if v.get("rf3_plddt") is not None] if validation_results else []

        class Result:
            success = True
            num_backbones = len(backbone_pdbs)
            num_sequences = len(all_sequences)
            pass_rate = n_passed / max(n_total, 1)
            best_rmsd = min(rmsds) if rmsds else None
            best_plddt = max(plddts_val) if plddts_val else 0.0
            error = None

        return Result()

    # Run the async function
    return asyncio.run(_run())


def run_designs(num_designs: int = 50, num_sequences_per_design: int = 4, use_mock: bool = False):
    """
    Run designs through the AI pipeline.

    Args:
        num_designs: Total number of backbone designs to generate
        num_sequences_per_design: Number of sequences per backbone
        use_mock: Use mock inference (for testing without API)
    """
    print("=" * 70)
    print(f"AI DESIGN PIPELINE - {num_designs} DESIGNS" + (" [MOCK MODE]" if use_mock else ""))
    print("=" * 70)
    print(f"\nQuery: {TEST_QUERY}")
    print(f"\nTarget: {num_designs} backbone designs, {num_sequences_per_design} sequences each")
    print(f"Total sequences: {num_designs * num_sequences_per_design}")
    print("=" * 70)

    # Import pipeline
    from ai_design_pipeline import AIDesignPipeline

    # Check for API key
    import os
    claude_api_key = os.environ.get("CLAUDE_API_KEY") or os.environ.get("ANTHROPIC_API_KEY")

    # Initialize pipeline
    pipeline = AIDesignPipeline(
        claude_api_key=claude_api_key,
        use_mock=use_mock,
    )

    print(f"\nPipeline initialized")
    print(f"  - Scaffolding available: {pipeline.scaffolding_workflow is not None}")
    print(f"  - Analyzer available: {pipeline.analyzer is not None}")

    # Run pipeline
    print(f"\nStarting design run at {datetime.now().strftime('%H:%M:%S')}...")
    start_time = time.time()

    result = pipeline.run(
        query=TEST_QUERY,
        num_designs=num_designs,
        num_sequences=num_sequences_per_design,
        validate_sequences=True,
    )

    elapsed = time.time() - start_time

    # Print results
    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)

    print(f"\nStatus: {'SUCCESS' if result.success else 'FAILED'}")
    if result.error:
        print(f"Error: {result.error}")

    print(f"\nDesign Statistics:")
    print(f"  - Backbones generated: {result.num_backbones}")
    print(f"  - Sequences designed: {result.num_sequences}")
    print(f"  - Validation pass rate: {result.pass_rate:.1%}")

    if result.best_rmsd and result.best_rmsd < float('inf'):
        print(f"\nBest Result:")
        print(f"  - RMSD: {result.best_rmsd:.2f} A")
        print(f"  - pLDDT: {result.best_plddt:.2f}")

    if result.best_sequence:
        print(f"\nBest Sequence ({len(result.best_sequence.sequence)} residues):")
        seq = result.best_sequence.sequence
        # Print sequence in rows of 60
        for i in range(0, len(seq), 60):
            print(f"  {seq[i:i+60]}")

    # Print timing
    print(f"\nTiming:")
    for stage, duration in result.timings.items():
        print(f"  - {stage}: {duration:.1f}s")
    print(f"  - TOTAL: {elapsed:.1f}s ({elapsed/60:.1f} min)")

    # Print analysis results if available
    if result.analysis_results:
        print(f"\nAnalysis Results:")

        # Stability profile
        if "stability_profile" in result.analysis_results:
            sp = result.analysis_results["stability_profile"]
            print(f"  Stability Profile: {sp.get('name', 'N/A')}")
            print(f"    - cfg_scale: {sp.get('cfg_scale', 'N/A')}")
            print(f"    - step_scale: {sp.get('step_scale', 'N/A')}")
            print(f"    - gamma_0: {sp.get('gamma_0', 'N/A')}")
            print(f"    - plddt_threshold: {sp.get('plddt_threshold', 'N/A')}")

        # Enzyme preservation
        if "enzyme_preservation" in result.analysis_results:
            ep = result.analysis_results["enzyme_preservation"]
            print(f"  Enzyme Preservation:")
            print(f"    - enzyme_class: {ep.get('enzyme_class', 'N/A')}")
            print(f"    - preserve_function: {ep.get('preserve_function', 'N/A')}")

        # Scaffolding info
        if "scaffolding_info" in result.analysis_results:
            si = result.analysis_results["scaffolding_info"]
            print(f"  Scaffolding Source:")
            print(f"    - PDB: {si.get('pdb_id', 'N/A')}")
            print(f"    - Metal: {si.get('metal', 'N/A')}")
            print(f"    - Ligand: {si.get('ligand', 'N/A')}")
            print(f"    - Coordination: {si.get('coordination_number', 'N/A')}")

    # Validation breakdown
    if result.validation_results:
        print(f"\nValidation Breakdown:")
        passed = sum(1 for v in result.validation_results if v.passed)
        failed = len(result.validation_results) - passed
        print(f"  - Passed: {passed}")
        print(f"  - Failed: {failed}")

        # RMSD distribution
        rmsds = [v.rmsd for v in result.validation_results if v.rmsd < float('inf')]
        if rmsds:
            print(f"  - RMSD range: {min(rmsds):.2f} - {max(rmsds):.2f} A")
            print(f"  - RMSD mean: {sum(rmsds)/len(rmsds):.2f} A")

        # pLDDT distribution
        plddts = [v.plddt for v in result.validation_results if v.plddt > 0]
        if plddts:
            print(f"  - pLDDT range: {min(plddts):.2f} - {max(plddts):.2f}")
            print(f"  - pLDDT mean: {sum(plddts)/len(plddts):.2f}")

    # Recommendations
    if result.recommendations:
        print(f"\nRecommendations:")
        for rec in result.recommendations[:5]:
            print(f"  - {rec}")

    # Save results
    output_dir = Path("design_results")
    output_dir.mkdir(exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = output_dir / f"design_run_{timestamp}.json"

    with open(output_file, 'w') as f:
        json.dump(result.to_dict(), f, indent=2, default=str)

    print(f"\nResults saved to: {output_file}")

    # Quality assessment
    print("\n" + "=" * 70)
    print("QUALITY ASSESSMENT")
    print("=" * 70)

    quality_score = 0
    quality_notes = []

    # Check pass rate
    if result.pass_rate >= 0.5:
        quality_score += 25
        quality_notes.append("[GOOD] Pass rate >= 50%")
    elif result.pass_rate >= 0.2:
        quality_score += 10
        quality_notes.append("[OK] Pass rate >= 20%")
    else:
        quality_notes.append("[POOR] Pass rate < 20%")

    # Check best RMSD
    if result.best_rmsd and result.best_rmsd < float('inf'):
        if result.best_rmsd < 1.5:
            quality_score += 25
            quality_notes.append("[EXCELLENT] Best RMSD < 1.5 A")
        elif result.best_rmsd < 2.5:
            quality_score += 15
            quality_notes.append("[GOOD] Best RMSD < 2.5 A")
        elif result.best_rmsd < 4.0:
            quality_score += 5
            quality_notes.append("[OK] Best RMSD < 4.0 A")
        else:
            quality_notes.append(f"[POOR] Best RMSD = {result.best_rmsd:.2f} A")

    # Check pLDDT
    if result.best_plddt >= 0.8:
        quality_score += 25
        quality_notes.append("[EXCELLENT] Best pLDDT >= 0.8")
    elif result.best_plddt >= 0.7:
        quality_score += 15
        quality_notes.append("[GOOD] Best pLDDT >= 0.7")
    elif result.best_plddt >= 0.5:
        quality_score += 5
        quality_notes.append("[OK] Best pLDDT >= 0.5")
    else:
        quality_notes.append(f"[POOR] Best pLDDT = {result.best_plddt:.2f}")

    # Check number of designs
    if result.num_backbones >= num_designs * 0.9:
        quality_score += 25
        quality_notes.append(f"[GOOD] Generated {result.num_backbones}/{num_designs} backbones")
    elif result.num_backbones >= num_designs * 0.5:
        quality_score += 10
        quality_notes.append(f"[OK] Generated {result.num_backbones}/{num_designs} backbones")
    else:
        quality_notes.append(f"[POOR] Only {result.num_backbones}/{num_designs} backbones")

    print(f"\nQuality Score: {quality_score}/100")
    print("\nAssessment:")
    for note in quality_notes:
        print(f"  {note}")

    # Overall verdict
    print("\n" + "-" * 40)
    if quality_score >= 75:
        print("VERDICT: EXCELLENT - Ready for experimental validation")
    elif quality_score >= 50:
        print("VERDICT: GOOD - Consider additional optimization")
    elif quality_score >= 25:
        print("VERDICT: FAIR - Needs parameter tuning")
    else:
        print("VERDICT: POOR - Review design strategy")

    return result


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(description="Run AI design pipeline")
    parser.add_argument("--num-designs", type=int, default=50, help="Number of backbone designs")
    parser.add_argument("--num-seqs", type=int, default=4, help="Sequences per backbone")
    parser.add_argument("--mock", action="store_true", help="Use mock inference for testing")
    parser.add_argument("--api", action="store_true", help="Use HTTP API mode (for running from Windows)")

    args = parser.parse_args()

    try:
        if args.api:
            # Use HTTP API mode
            result = run_designs_via_api(
                num_designs=args.num_designs,
                num_sequences_per_design=args.num_seqs,
            )
        else:
            # Use direct mode (inside Docker) or mock
            result = run_designs(
                num_designs=args.num_designs,
                num_sequences_per_design=args.num_seqs,
                use_mock=args.mock,
            )
        return 0 if result and result.success else 1
    except Exception as e:
        print(f"\n[ERROR] Pipeline failed: {e}")
        import traceback
        traceback.print_exc()
        return 1


NUM_DESIGNS = 5  # Set to 50 for full run, 5 for quick test


if __name__ == "__main__":
    sys.exit(main())
