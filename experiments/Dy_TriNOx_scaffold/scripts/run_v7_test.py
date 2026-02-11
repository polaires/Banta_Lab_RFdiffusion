#!/usr/bin/env python3
"""Run v7 test: 10 large scaffold designs with inner hotspots (near metal)."""

import json
import subprocess
import os
from pathlib import Path

# Configuration
BASE_DIR = Path("/mnt/g/Github_local_repo/Banta_Lab_RFdiffusion/experiments/Dy_TriNOx_scaffold")
OUTPUT_DIR = BASE_DIR / "outputs_v7"
INPUT_PDB = BASE_DIR / "inputs" / "Dy_TriNOx_UNL.pdb"
API_URL = "http://localhost:8000/runsync"

# V7 test config - inner hotspots near metal
TEST_CONFIG = {
    "name": "v7_large_test",
    "contig": "X1,L1,/0,120-140",
    "num_designs": 10
}

# V7 base params - hotspots on metal + phenoxide oxygens
BASE_PARAMS = {
    "task": "rfd3",
    "use_classifier_free_guidance": True,
    "cfg_scale": 2.5,
    "step_scale": 1.5,
    "num_timesteps": 200,
    "select_fixed_atoms": {"X1": "all", "L1": "all"},
    # Inner hotspots: metal + phenoxide oxygens + cap nitrogen
    "select_hotspots": {"X1": "1", "L1": "O1,O2,O3,N2"},
    "select_buried": {"L1": "all"},
    # H-bond acceptors: phenoxide oxygens directly coordinate Dy
    "select_hbond_acceptor": {"L1": "O1,O2,O3"},
    # Position protein above the complex (toward open coordination site)
    "ori_token": [0.0, 10.0, 5.0]
}


def run_rfd3_batch(name: str, contig: str, num_designs: int, pdb_content: str) -> list:
    """Run a single RFD3 batch and return designs."""
    params = {**BASE_PARAMS, "contig": contig, "num_designs": num_designs, "pdb_content": pdb_content}

    payload = json.dumps({"input": params})

    print(f"Running {name} with params:")
    print(f"  contig: {contig}")
    print(f"  num_designs: {num_designs}")
    print(f"  select_hotspots: {params['select_hotspots']}")
    print(f"  select_buried: {params['select_buried']}")
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
    print(f"V7 TEST: {name}")
    print(f"Key change: Inner hotspots (phenoxide O + metal + N2)")
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

    # Distance check on all designs
    if designs:
        print(f"\n{'='*60}")
        print("Distance analysis for all designs:")
        print(f"{'='*60}")

        successes = 0
        for i, design in enumerate(designs):
            print(f"\nDesign {i:03d}:")
            if check_design_distances(design["content"]):
                successes += 1

        print(f"\n{'='*60}")
        print(f"SUMMARY: {successes}/{len(designs)} designs with protein near metal (<5A)")
        print(f"{'='*60}")


def check_design_distances(pdb_content: str) -> bool:
    """Check distances from protein to metal AND ligand in a design."""
    import numpy as np

    metal_coord = None
    ligand_coords = []
    inner_ligand_coords = []  # O1, O2, O3, N2
    protein_coords = []

    inner_atoms = {'O1', 'O2', 'O3', 'N2'}

    for line in pdb_content.strip().split('\n'):
        if line.startswith('HETATM') or line.startswith('ATOM'):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coord = np.array([x, y, z])

                res_name = line[17:20].strip()
                atom_name = line[12:16].strip()

                if line.startswith('HETATM'):
                    if res_name == 'DY' or atom_name == 'DY':
                        metal_coord = coord
                    elif res_name == 'UNL':
                        ligand_coords.append(coord)
                        if atom_name in inner_atoms:
                            inner_ligand_coords.append(coord)
                elif line.startswith('ATOM'):
                    protein_coords.append(coord)
            except:
                continue

    if metal_coord is None:
        print("  No metal found!")
        return False

    if not protein_coords:
        print("  No protein atoms found!")
        return False

    protein_coords = np.array(protein_coords)
    ligand_coords = np.array(ligand_coords) if ligand_coords else None
    inner_ligand_coords = np.array(inner_ligand_coords) if inner_ligand_coords else None

    # Distance to metal
    metal_distances = np.linalg.norm(protein_coords - metal_coord, axis=1)
    metal_min = np.min(metal_distances)
    metal_close_3A = np.sum(metal_distances < 3.0)
    metal_close_5A = np.sum(metal_distances < 5.0)

    # Distance to inner ligand atoms (O1, O2, O3, N2)
    if inner_ligand_coords is not None and len(inner_ligand_coords) > 0:
        inner_min_dists = []
        for pc in protein_coords:
            dists = np.linalg.norm(inner_ligand_coords - pc, axis=1)
            inner_min_dists.append(np.min(dists))
        inner_min_dists = np.array(inner_min_dists)
        inner_min = np.min(inner_min_dists)
        inner_close_3A = np.sum(inner_min_dists < 3.0)
        inner_close_5A = np.sum(inner_min_dists < 5.0)
    else:
        inner_min = None

    # Distance to all ligand atoms
    if ligand_coords is not None:
        ligand_min_dists = []
        for pc in protein_coords:
            dists = np.linalg.norm(ligand_coords - pc, axis=1)
            ligand_min_dists.append(np.min(dists))
        ligand_min_dists = np.array(ligand_min_dists)
        ligand_min = np.min(ligand_min_dists)
        ligand_close_5A = np.sum(ligand_min_dists < 5.0)
    else:
        ligand_min = None

    print(f"  Protein atoms: {len(protein_coords)}")
    print(f"  Metal (Dy) - Min dist: {metal_min:.2f}A, <3A: {metal_close_3A}, <5A: {metal_close_5A}")
    if inner_min is not None:
        print(f"  Inner (O1,O2,O3,N2) - Min dist: {inner_min:.2f}A, <3A: {inner_close_3A}, <5A: {inner_close_5A}")
    if ligand_min is not None:
        print(f"  Ligand (all) - Min dist: {ligand_min:.2f}A, <5A: {ligand_close_5A}")

    # Success criteria
    success = False
    if metal_close_3A > 0:
        print(f"  -> EXCELLENT: Protein within coordination distance of metal!")
        success = True
    elif metal_min < 5.0:
        print(f"  -> GOOD: Protein close to metal (<5A)")
        success = True
    elif inner_min is not None and inner_min < 3.0:
        print(f"  -> OK: Protein contacts inner atoms")
        success = True
    else:
        print(f"  -> FAIL: Protein not close enough to metal")

    return success


if __name__ == "__main__":
    main()
