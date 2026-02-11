#!/usr/bin/env python3
"""Run v10 test: Focus on ligand N atoms with correct H-bond acceptor conditioning."""

import json
import subprocess
import os
from pathlib import Path
import numpy as np

# Configuration
BASE_DIR = Path("G:/Github_local_repo/Banta_Lab_RFdiffusion/experiments/Dy_TriNOx_scaffold")
OUTPUT_DIR = BASE_DIR / "outputs_v10"
INPUT_PDB = BASE_DIR / "inputs" / "Dy_TriNOx_UNL.pdb"
API_URL = "http://localhost:8000/runsync"

# V10 test config
TEST_CONFIG = {
    "name": "v10_N_focus",
    "contig": "X1,L1,/0,120-140",
    "num_designs": 10
}

# V10 base params - correct H-bond conditioning
BASE_PARAMS = {
    "task": "rfd3",
    "use_classifier_free_guidance": True,
    "cfg_scale": 3.0,  # Higher for stronger conditioning
    "step_scale": 1.5,
    "num_timesteps": 200,
    "select_fixed_atoms": {"X1": "all", "L1": "all"},
    # Hotspots ONLY on N atoms (no aromatic C's)
    "select_hotspots": {"L1": "N1,N2,N3,N4"},
    "select_buried": {"L1": "all"},
    # CORRECT: N atoms are H-bond ACCEPTORS (have lone pairs)
    # Protein provides H-bond DONORS (NH, OH groups)
    "select_hbond_acceptor": {"L1": "N1,N2,N3,N4"},
    # Smaller negative Z - closer to where N atoms actually are
    "ori_token": [0.0, 8.0, -2.0]
}


def run_rfd3_batch(name: str, contig: str, num_designs: int, pdb_content: str) -> list:
    """Run a single RFD3 batch and return designs."""
    params = {**BASE_PARAMS, "contig": contig, "num_designs": num_designs, "pdb_content": pdb_content}

    payload = json.dumps({"input": params})

    print(f"Running {name} with params:")
    print(f"  contig: {contig}")
    print(f"  num_designs: {num_designs}")
    print(f"  select_hotspots: {params['select_hotspots']}")
    print(f"  select_hbond_acceptor: {params['select_hbond_acceptor']}")
    print(f"  ori_token: {params['ori_token']}")
    print(f"  cfg_scale: {params['cfg_scale']}")

    result = subprocess.run(
        ["curl", "-s", "-X", "POST", API_URL, "-H", "Content-Type: application/json", "-d", payload],
        capture_output=True, text=True, timeout=600
    )

    try:
        response = json.loads(result.stdout)
        if response.get("status") == "COMPLETED":
            return response.get("output", {}).get("result", {}).get("designs", [])
        else:
            error = response.get("output", {}).get("error", response.get("error", "Unknown error"))
            print(f"  ERROR: {str(error)[:200]}")
            return []
    except json.JSONDecodeError as e:
        print(f"  ERROR: Failed to parse response: {e}")
        print(f"  Response: {result.stdout[:500]}")
        return []


def check_ligand_hbonds(pdb_content: str) -> dict:
    """Check H-bonds from protein to ligand N atoms."""
    ligand_N_coords = {}  # atom_name -> coord
    protein_donors = []  # (coord, atom_name, res_name, res_num)
    protein_coords = []
    metal_coord = None

    for line in pdb_content.strip().split('\n'):
        if line.startswith('HETATM') or line.startswith('ATOM'):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coord = np.array([x, y, z])
                res_name = line[17:20].strip()
                atom_name = line[12:16].strip()
                res_num = int(line[22:26])

                if line.startswith('HETATM'):
                    if res_name == 'DY':
                        metal_coord = coord
                    elif res_name == 'UNL' and atom_name in ['N1', 'N2', 'N3', 'N4']:
                        ligand_N_coords[atom_name] = coord
                elif line.startswith('ATOM'):
                    protein_coords.append(coord)
                    # H-bond donors (NH, OH groups that can donate H)
                    if atom_name in ['N', 'NE', 'NH1', 'NH2', 'ND1', 'ND2', 'NE1', 'NE2', 'NZ', 'OG', 'OG1', 'OH']:
                        protein_donors.append((coord, atom_name, res_name, res_num))
            except:
                continue

    result = {
        'hbonds': [],
        'n_distances': {},
        'z_offset': 0,
    }

    if not ligand_N_coords or not protein_coords:
        return result

    protein_coords = np.array(protein_coords)

    # Check H-bonds to each ligand N
    for n_atom, n_coord in ligand_N_coords.items():
        min_dist = 999
        closest = None
        for (pc, pa, pr, pn) in protein_donors:
            d = np.linalg.norm(n_coord - pc)
            if d < min_dist:
                min_dist = d
                closest = f"{pr}{pn}.{pa}"
            if d < 3.5:  # H-bond distance
                result['hbonds'].append({
                    'ligand_atom': n_atom,
                    'protein': f"{pr}{pn}.{pa}",
                    'distance': d
                })
        result['n_distances'][n_atom] = (min_dist, closest)

    # Z offset
    if metal_coord is not None:
        result['z_offset'] = np.mean(protein_coords[:, 2]) - metal_coord[2]

    return result


def main():
    import time

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Read input PDB
    with open(INPUT_PDB) as f:
        pdb_content = f.read()

    name = TEST_CONFIG["name"]
    contig = TEST_CONFIG["contig"]
    num_designs = TEST_CONFIG["num_designs"]

    print(f"\n{'='*60}")
    print(f"V10 TEST: {name}")
    print(f"Key changes from v9:")
    print(f"  - CORRECT: select_hbond_acceptor (not donor)")
    print(f"  - Hotspots ONLY on N atoms")
    print(f"  - Smaller Z offset (-2 vs -5)")
    print(f"  - Higher cfg_scale (3.0)")
    print(f"{'='*60}")

    start_time = time.time()
    designs = run_rfd3_batch(name, contig, num_designs, pdb_content)
    elapsed = time.time() - start_time

    # Save designs
    for i, design in enumerate(designs):
        filename = OUTPUT_DIR / f"{name}_{i:03d}.pdb"
        with open(filename, "w") as f:
            f.write(design["content"])

    print(f"\n  Generated {len(designs)} designs in {elapsed:.2f} seconds")
    print(f"  Saved to {OUTPUT_DIR}")

    # H-bond analysis
    if designs:
        print(f"\n{'='*60}")
        print("H-BOND ANALYSIS (protein donors -> ligand N acceptors):")
        print(f"{'='*60}")

        total_hbonds = 0
        designs_with_hbonds = 0

        for i, design in enumerate(designs):
            result = check_ligand_hbonds(design["content"])
            n_hbonds = len(result['hbonds'])
            total_hbonds += n_hbonds
            if n_hbonds > 0:
                designs_with_hbonds += 1

            status = 'EXCELLENT' if n_hbonds >= 3 else 'GOOD' if n_hbonds >= 2 else 'OK' if n_hbonds >= 1 else 'NONE'
            print(f"\nDesign {i:03d}: {n_hbonds} H-bonds, Z-offset={result['z_offset']:+.1f}A [{status}]")

            # Show distances to each N
            for n_atom, (dist, closest) in sorted(result['n_distances'].items()):
                hb_mark = " <-- H-BOND" if dist < 3.5 else " (close)" if dist < 4.5 else ""
                print(f"    {n_atom}: {dist:.1f}A to {closest}{hb_mark}")

        print(f"\n{'='*60}")
        print(f"SUMMARY:")
        print(f"  Designs with H-bonds: {designs_with_hbonds}/{len(designs)}")
        print(f"  Total H-bonds: {total_hbonds}")
        print(f"  Avg H-bonds/design: {total_hbonds/len(designs):.1f}")
        print(f"{'='*60}")


if __name__ == "__main__":
    main()
