#!/usr/bin/env python3
"""Run v6 test: 10 large scaffold designs with hotspots on TriNOx ligand."""

import json
import subprocess
import os
from pathlib import Path

# Configuration
BASE_DIR = Path("/mnt/g/Github_local_repo/Banta_Lab_RFdiffusion/experiments/Dy_TriNOx_scaffold")
OUTPUT_DIR = BASE_DIR / "outputs_v6"
INPUT_PDB = BASE_DIR / "inputs" / "Dy_TriNOx_UNL.pdb"
API_URL = "http://localhost:8000/runsync"

# V6 test config - hotspots on TriNOx ligand
TEST_CONFIG = {
    "name": "v6_large_test",
    "contig": "X1,L1,/0,120-140",
    "num_designs": 10
}

# V6 base params - priority is binding TriNOx, not metal directly
BASE_PARAMS = {
    "task": "rfd3",
    "use_classifier_free_guidance": True,
    "cfg_scale": 2.5,
    "step_scale": 1.5,
    "num_timesteps": 200,
    "select_fixed_atoms": {"X1": "all", "L1": "all"},
    # Hotspots on TriNOx peripheral atoms - aromatic carbons + amide nitrogens
    "select_hotspots": {"L1": "C4,C5,C15,C16,C26,C27,N1,N3,N4"},
    "select_buried": {"L1": "all"},
    # H-bond acceptors: phenoxide oxygens + amide nitrogens
    "select_hbond_acceptor": {"L1": "O1,O2,O3,N1,N3,N4"},
    "ori_token": [0.0, 5.0, 0.0]
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
    print(f"V6 TEST: {name}")
    print(f"Key change: Hotspots on TriNOx ligand (not metal)")
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

        for i, design in enumerate(designs):
            print(f"\nDesign {i:03d}:")
            check_design_distances(design["content"])


def check_design_distances(pdb_content: str):
    """Check distances from protein to metal AND ligand in a design."""
    import numpy as np

    metal_coord = None
    ligand_coords = []
    protein_coords = []

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
                elif line.startswith('ATOM'):
                    protein_coords.append(coord)
            except:
                continue

    if metal_coord is None:
        print("  No metal found!")
        return

    if not protein_coords:
        print("  No protein atoms found!")
        return

    protein_coords = np.array(protein_coords)
    ligand_coords = np.array(ligand_coords) if ligand_coords else None

    # Distance to metal
    metal_distances = np.linalg.norm(protein_coords - metal_coord, axis=1)
    metal_min = np.min(metal_distances)
    metal_close_3A = np.sum(metal_distances < 3.0)
    metal_close_5A = np.sum(metal_distances < 5.0)

    # Distance to ligand (nearest ligand atom)
    if ligand_coords is not None:
        # For each protein atom, find distance to nearest ligand atom
        ligand_min_dists = []
        for pc in protein_coords:
            dists = np.linalg.norm(ligand_coords - pc, axis=1)
            ligand_min_dists.append(np.min(dists))
        ligand_min_dists = np.array(ligand_min_dists)
        ligand_min = np.min(ligand_min_dists)
        ligand_close_3A = np.sum(ligand_min_dists < 3.0)
        ligand_close_5A = np.sum(ligand_min_dists < 5.0)
    else:
        ligand_min = None

    print(f"  Protein atoms: {len(protein_coords)}")
    print(f"  Metal - Min dist: {metal_min:.2f}A, <3A: {metal_close_3A}, <5A: {metal_close_5A}")
    if ligand_min is not None:
        print(f"  Ligand - Min dist: {ligand_min:.2f}A, <3A: {ligand_close_3A}, <5A: {ligand_close_5A}")

    # Success criteria
    success = False
    if ligand_min is not None and ligand_min < 4.0:
        print(f"  -> GOOD: Protein contacts ligand!")
        success = True
    if metal_min < 4.0:
        print(f"  -> GREAT: Protein close to metal!")
        success = True
    if not success:
        print(f"  -> FAIL: Protein too far from ligand/metal")


if __name__ == "__main__":
    main()
