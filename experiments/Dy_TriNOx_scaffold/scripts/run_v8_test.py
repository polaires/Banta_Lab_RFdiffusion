#!/usr/bin/env python3
"""Run v8 test: Force protein approach from +Z (above) through open coordination site."""

import json
import subprocess
import os
from pathlib import Path

# Configuration
BASE_DIR = Path("G:/Github_local_repo/Banta_Lab_RFdiffusion/experiments/Dy_TriNOx_scaffold")
OUTPUT_DIR = BASE_DIR / "outputs_v8"
INPUT_PDB = BASE_DIR / "inputs" / "Dy_TriNOx_UNL.pdb"
API_URL = "http://localhost:8000/runsync"

# V8 test config
TEST_CONFIG = {
    "name": "v8_large_test",
    "contig": "X1,L1,/0,120-140",
    "num_designs": 10
}

# V8 base params - force approach from above (+Z)
BASE_PARAMS = {
    "task": "rfd3",
    "use_classifier_free_guidance": True,
    "cfg_scale": 3.0,  # Stronger guidance
    "step_scale": 1.5,
    "num_timesteps": 200,
    "select_fixed_atoms": {"X1": "all", "L1": "all"},
    # Peripheral hotspots (v6 style - proven to work)
    "select_hotspots": {"L1": "C4,C5,C15,C16,C26,C27,N1,N3,N4"},
    "select_buried": {"L1": "all"},
    "select_hbond_acceptor": {"L1": "O1,O2,O3"},
    # Strong +Z ori_token - force protein ABOVE the complex
    "ori_token": [0.0, 8.0, 15.0]
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
    print(f"V8 TEST: {name}")
    print(f"Key change: Strong +Z ori_token to force approach from above")
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
        close_to_metal = []
        for i, design in enumerate(designs):
            print(f"\nDesign {i:03d}:")
            result = check_design_distances(design["content"])
            if result['metal_min'] < 5.0:
                successes += 1
            close_to_metal.append(result['metal_min'])

        print(f"\n{'='*60}")
        print(f"SUMMARY:")
        print(f"  Designs with protein <5A from metal: {successes}/{len(designs)}")
        print(f"  Metal distances: min={min(close_to_metal):.2f}A, max={max(close_to_metal):.2f}A, avg={sum(close_to_metal)/len(close_to_metal):.2f}A")
        print(f"{'='*60}")


def check_design_distances(pdb_content: str) -> dict:
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

    result = {'metal_min': 999, 'ligand_min': 999}

    if metal_coord is None:
        print("  No metal found!")
        return result

    if not protein_coords:
        print("  No protein atoms found!")
        return result

    protein_coords = np.array(protein_coords)
    ligand_coords = np.array(ligand_coords) if ligand_coords else None

    # Distance to metal
    metal_distances = np.linalg.norm(protein_coords - metal_coord, axis=1)
    metal_min = np.min(metal_distances)
    metal_close_3A = np.sum(metal_distances < 3.0)
    metal_close_5A = np.sum(metal_distances < 5.0)
    result['metal_min'] = metal_min

    # Distance to ligand
    if ligand_coords is not None:
        ligand_min_dists = []
        for pc in protein_coords:
            dists = np.linalg.norm(ligand_coords - pc, axis=1)
            ligand_min_dists.append(np.min(dists))
        ligand_min_dists = np.array(ligand_min_dists)
        ligand_min = np.min(ligand_min_dists)
        ligand_close_5A = np.sum(ligand_min_dists < 5.0)
        result['ligand_min'] = ligand_min
    else:
        ligand_min = None
        ligand_close_5A = 0

    # Calculate protein Z centroid relative to metal
    protein_z_mean = np.mean(protein_coords[:, 2])
    metal_z = metal_coord[2]
    z_offset = protein_z_mean - metal_z

    print(f"  Protein atoms: {len(protein_coords)}")
    print(f"  Metal (Dy) - Min dist: {metal_min:.2f}A, <3A: {metal_close_3A}, <5A: {metal_close_5A}")
    if ligand_min is not None:
        print(f"  Ligand (all) - Min dist: {ligand_min:.2f}A, <5A: {ligand_close_5A}")
    print(f"  Protein Z offset from metal: {z_offset:+.1f}A (positive=above)")

    # Success criteria
    if metal_close_3A > 0:
        print(f"  -> EXCELLENT: Protein within coordination distance!")
    elif metal_min < 5.0:
        print(f"  -> GOOD: Protein close to metal (<5A)")
    elif ligand_min is not None and ligand_min < 3.0:
        print(f"  -> OK: Contacts ligand")
    else:
        print(f"  -> FAIL: Not close enough")

    return result


if __name__ == "__main__":
    main()
