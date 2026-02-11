#!/usr/bin/env python3
"""Run v5 test: 10 large scaffold designs with hotspot + explicit origin."""

import json
import subprocess
import os
from pathlib import Path

# Configuration
BASE_DIR = Path("/mnt/g/Github_local_repo/Banta_Lab_RFdiffusion/experiments/Dy_TriNOx_scaffold")
OUTPUT_DIR = BASE_DIR / "outputs_v5"
INPUT_PDB = BASE_DIR / "inputs" / "Dy_TriNOx_UNL.pdb"
API_URL = "http://localhost:8000/runsync"

# V5 test config
TEST_CONFIG = {
    "name": "v5_large_test",
    "contig": "X1,L1,/0,120-140",
    "num_designs": 10
}

BASE_PARAMS = {
    "task": "rfd3",
    "use_classifier_free_guidance": True,
    "cfg_scale": 2.5,
    "step_scale": 1.5,
    "num_timesteps": 200,
    "select_fixed_atoms": {"X1": "all", "L1": "all"},
    "select_hotspots": {"X1": "1"},
    "select_buried": {"L1": "all"},
    "select_hbond_acceptor": {"L1": "O1,O2,O3"},
    "ori_token": [0.0, 7.8, 5.5]
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
    print(f"V5 TEST: {name}")
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

    # Quick distance check on first design if any
    if designs:
        print(f"\n{'='*60}")
        print("Quick distance check on first design:")
        print(f"{'='*60}")
        check_design_distances(designs[0]["content"])


def check_design_distances(pdb_content: str):
    """Quick check of distances in a design."""
    import numpy as np

    metal_coord = None
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
    distances = np.linalg.norm(protein_coords - metal_coord, axis=1)

    min_dist = np.min(distances)
    close_3A = np.sum(distances < 3.0)
    close_5A = np.sum(distances < 5.0)

    print(f"  Metal position: {metal_coord}")
    print(f"  Protein atoms: {len(protein_coords)}")
    print(f"  Min distance to metal: {min_dist:.2f} Å")
    print(f"  Atoms within 3Å: {close_3A}")
    print(f"  Atoms within 5Å: {close_5A}")

    if min_dist < 4.0:
        print(f"\n  ✓ SUCCESS: Protein is close to metal!")
    else:
        print(f"\n  ✗ FAIL: Protein still too far from metal")


if __name__ == "__main__":
    main()
