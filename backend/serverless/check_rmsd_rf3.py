"""
RF3 validation of all designed sequences.
- Protein-only RF3 prediction
- RF3 with PQQ ligand (PubChem canonical SMILES)
- RMSD computation from mmCIF output vs designed backbones
"""
import sys
import json
import math
import re
import requests
import builtins
from pathlib import Path

_orig_print = builtins.print
def print(*args, **kwargs):
    kwargs.setdefault('flush', True)
    _orig_print(*args, **kwargs)

BASE_URL = "http://localhost:8000"
OUTPUT_DIR = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless\output_minimal_motif_25")

# PubChem canonical SMILES for PQQ (CID 1024)
PQQ_SMILES = "C1=C(C2=C(C(=O)C(=O)C3=C2NC(=C3)C(=O)O)N=C1C(=O)O)C(=O)O"


# Import shared utilities from esmfold_utils (avoids code duplication)
try:
    from esmfold_utils import (
        _kabsch_rmsd as _kabsch_rmsd_shared,
        _extract_ca_coords_from_pdb as _extract_ca_pdb_shared,
        _extract_ca_coords_from_cif as _extract_ca_cif_shared,
    )

    def parse_ca_from_pdb(pdb_text):
        """Wrapper: delegates to esmfold_utils._extract_ca_coords_from_pdb."""
        coords = _extract_ca_pdb_shared(pdb_text)
        return [tuple(row) for row in coords.tolist()] if len(coords) > 0 else []

    def parse_ca_from_cif(cif_text):
        """Wrapper: delegates to esmfold_utils._extract_ca_coords_from_cif."""
        coords = _extract_ca_cif_shared(cif_text)
        return [tuple(row) for row in coords.tolist()] if len(coords) > 0 else []

    def kabsch_rmsd(coords1, coords2):
        """Wrapper: delegates to esmfold_utils._kabsch_rmsd."""
        import numpy as np
        n = min(len(coords1), len(coords2))
        if n == 0:
            return float("inf")
        c1 = np.array(coords1[:n], dtype=np.float64)
        c2 = np.array(coords2[:n], dtype=np.float64)
        return float(_kabsch_rmsd_shared(c1, c2))

    def _kabsch_rmsd_numpy(coords1, coords2):
        """Wrapper: delegates to esmfold_utils._kabsch_rmsd."""
        import numpy as np
        c1 = np.array(coords1, dtype=np.float64)
        c2 = np.array(coords2, dtype=np.float64)
        return float(_kabsch_rmsd_shared(c1, c2))

except ImportError:
    # Standalone fallback — keep local definitions for out-of-container use
    def parse_ca_from_pdb(pdb_text):
        """Extract CA atom coordinates from PDB text."""
        cas = []
        for line in pdb_text.splitlines():
            if (line.startswith("ATOM") or line.startswith("HETATM")):
                atom_name = line[12:16].strip()
                if atom_name == "CA":
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    cas.append((x, y, z))
        return cas

    def parse_ca_from_cif(cif_text):
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
            if in_atom_site and line.startswith("_"):
                in_atom_site = False
                col_map = {}
                col_idx = 0
                continue
            if in_atom_site and line.startswith("#"):
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
                label_atom_id = col_map.get("label_atom_id")
                if label_atom_id is None:
                    label_atom_id = col_map.get("auth_atom_id")
                if label_atom_id is None:
                    continue
                atom_name = parts[label_atom_id].strip('"')
                if atom_name != "CA":
                    continue
                group_idx = col_map.get("group_PDB")
                if group_idx is not None and parts[group_idx] not in ("ATOM",):
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

    def kabsch_rmsd(coords1, coords2):
        """Compute RMSD after optimal Kabsch superposition."""
        n = min(len(coords1), len(coords2))
        if n == 0:
            return float("inf")
        try:
            import numpy as np
            return _kabsch_rmsd_numpy(coords1[:n], coords2[:n])
        except ImportError:
            pass
        # Fallback: unaligned RMSD (upper bound)
        c1 = coords1[:n]
        c2 = coords2[:n]
        cx1 = sum(x for x, y, z in c1) / n
        cy1 = sum(y for x, y, z in c1) / n
        cz1 = sum(z for x, y, z in c1) / n
        cx2 = sum(x for x, y, z in c2) / n
        cy2 = sum(y for x, y, z in c2) / n
        cz2 = sum(z for x, y, z in c2) / n
        p = [(x - cx1, y - cy1, z - cz1) for x, y, z in c1]
        q = [(x - cx2, y - cy2, z - cz2) for x, y, z in c2]
        rmsd_sq = 0.0
        for i in range(n):
            dx = p[i][0] - q[i][0]
            dy = p[i][1] - q[i][1]
            dz = p[i][2] - q[i][2]
            rmsd_sq += dx * dx + dy * dy + dz * dz
        return math.sqrt(rmsd_sq / n)

    def _kabsch_rmsd_numpy(coords1, coords2):
        """Proper Kabsch RMSD using numpy SVD."""
        import numpy as np
        P = np.array(coords1, dtype=np.float64)
        Q = np.array(coords2, dtype=np.float64)
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


def predict_rf3(sequence, name, ligand_smiles=None):
    """Call RF3 prediction API."""
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
        return {"error": f"status={raw.get('status')}", "raw": raw}

    output = raw.get("output", {})
    if output.get("status") not in ("success", "completed"):
        return {"error": output.get("error", f"status={output.get('status')}"), "raw": raw}

    result = output.get("result", {})
    confidences = result.get("confidences", {})
    predictions = result.get("predictions", [])

    cif_content = ""
    if predictions and isinstance(predictions, list):
        cif_content = predictions[0].get("content", "")

    return {
        "mean_plddt": confidences.get("mean_plddt"),
        "ptm": confidences.get("ptm"),
        "overall_plddt": confidences.get("overall_plddt"),
        "overall_pae": confidences.get("overall_pae"),
        "ranking_score": confidences.get("ranking_score"),
        "cif_content": cif_content,
    }


# ---- Load sequences ----
seqs = []
seq_names = []
with open(OUTPUT_DIR / "sequences.fasta") as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            seq_names.append(line[1:])
        elif line:
            seqs.append(line)

print(f"Loaded {len(seqs)} sequences from {len(set(seq_names))} designs")

# ---- Load designed backbone CA coordinates ----
backbone_cas = {}
for i in range(1, 4):
    pdb_path = OUTPUT_DIR / f"backbone_{i}.pdb"
    if pdb_path.exists():
        backbone_cas[i] = parse_ca_from_pdb(pdb_path.read_text())
        print(f"Backbone {i}: {len(backbone_cas[i])} CA atoms")

# ---- Test PQQ SMILES first ----
print("\n=== Testing PQQ SMILES with RF3 ===")
test_result = predict_rf3(seqs[0], "smiles_test", ligand_smiles=PQQ_SMILES)
if "error" in test_result:
    print(f"PQQ SMILES FAILED: {test_result['error']}")
    pqq_works = False
else:
    print(f"PQQ SMILES works! pLDDT={test_result['mean_plddt']:.3f}, pTM={test_result['ptm']:.3f}")
    pqq_works = True

# ---- Run protein-only RF3 on all sequences ----
print("\n=== RF3 Protein-Only Validation ===")
results = []
for idx, seq in enumerate(seqs):
    backbone_num = (idx // 4) + 1
    seq_num = (idx % 4) + 1
    name = f"bb{backbone_num}_seq{seq_num}"

    print(f"\n[{idx+1}/{len(seqs)}] {name} ({len(seq)} res)...", end=" ")
    pred = predict_rf3(seq, name)

    if "error" in pred:
        print(f"ERROR: {pred['error']}")
        results.append({
            "backbone": backbone_num, "seq": seq_num,
            "rf3_plddt": None, "rf3_ptm": None, "rf3_rmsd": None,
            "error": pred["error"],
        })
        continue

    plddt = pred["mean_plddt"]
    ptm = pred["ptm"]

    # Compute RMSD if we have CIF and backbone
    rmsd = None
    if pred["cif_content"] and backbone_num in backbone_cas:
        pred_cas = parse_ca_from_cif(pred["cif_content"])
        ref_cas = backbone_cas[backbone_num]
        if pred_cas:
            rmsd = kabsch_rmsd(ref_cas, pred_cas)

    print(f"pLDDT={plddt:.3f}, pTM={ptm:.3f}", end="")
    if rmsd is not None:
        print(f", RMSD={rmsd:.2f}A ({len(parse_ca_from_cif(pred['cif_content']))} vs {len(backbone_cas.get(backbone_num, []))} CAs)")
    else:
        print()

    results.append({
        "backbone": backbone_num, "seq": seq_num,
        "rf3_plddt": plddt, "rf3_ptm": ptm, "rf3_rmsd": rmsd,
    })

# ---- Run RF3 with PQQ if SMILES works ----
if pqq_works:
    print("\n=== RF3 + PQQ Ligand Validation ===")
    for idx, seq in enumerate(seqs):
        backbone_num = (idx // 4) + 1
        seq_num = (idx % 4) + 1
        name = f"bb{backbone_num}_seq{seq_num}_pqq"

        print(f"\n[{idx+1}/{len(seqs)}] {name}...", end=" ")
        pred = predict_rf3(seq, name, ligand_smiles=PQQ_SMILES)

        if "error" in pred:
            print(f"ERROR: {pred['error']}")
            results[idx]["rf3_pqq_plddt"] = None
            results[idx]["rf3_pqq_ptm"] = None
            continue

        plddt = pred["mean_plddt"]
        ptm = pred["ptm"]
        print(f"pLDDT={plddt:.3f}, pTM={ptm:.3f}")
        results[idx]["rf3_pqq_plddt"] = plddt
        results[idx]["rf3_pqq_ptm"] = ptm

# ---- Summary ----
print("\n\n=== SUMMARY ===")
print(f"{'BB':>3} {'Seq':>3} | {'pLDDT':>6} {'pTM':>5} {'RMSD':>7}", end="")
if pqq_works:
    print(f" | {'PQQ_pLDDT':>9} {'PQQ_pTM':>7}")
else:
    print()

for r in results:
    plddt_str = f"{r['rf3_plddt']:.3f}" if r.get('rf3_plddt') is not None else "ERROR"
    ptm_str = f"{r['rf3_ptm']:.3f}" if r.get('rf3_ptm') is not None else "ERROR"
    rmsd_str = f"{r['rf3_rmsd']:.2f}A" if r.get('rf3_rmsd') is not None else "N/A"
    print(f"{r['backbone']:>3} {r['seq']:>3} | {plddt_str:>6} {ptm_str:>5} {rmsd_str:>7}", end="")
    if pqq_works:
        pqq_plddt = f"{r['rf3_pqq_plddt']:.3f}" if r.get('rf3_pqq_plddt') is not None else "N/A"
        pqq_ptm = f"{r['rf3_pqq_ptm']:.3f}" if r.get('rf3_pqq_ptm') is not None else "N/A"
        print(f" | {pqq_plddt:>9} {pqq_ptm:>7}")
    else:
        print()

# Save results
out_path = OUTPUT_DIR / "rf3_validation_results.json"
with open(out_path, "w") as f:
    json.dump(results, f, indent=2)
print(f"\nResults saved to {out_path}")
