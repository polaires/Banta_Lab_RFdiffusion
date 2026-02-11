#!/usr/bin/env python3
"""Run v9 test: Maximize protein-TriNOx H-bonding by targeting accessible N atoms."""

import json
import subprocess
import os
from pathlib import Path
import numpy as np

# Configuration
BASE_DIR = Path("G:/Github_local_repo/Banta_Lab_RFdiffusion/experiments/Dy_TriNOx_scaffold")
OUTPUT_DIR = BASE_DIR / "outputs_v9"
INPUT_PDB = BASE_DIR / "inputs" / "Dy_TriNOx_UNL.pdb"
API_URL = "http://localhost:8000/runsync"

# V9 test config
TEST_CONFIG = {
    "name": "v9_hbond_test",
    "contig": "X1,L1,/0,120-140",
    "num_designs": 10
}

# V9 base params - target ligand N atoms for H-bonding
BASE_PARAMS = {
    "task": "rfd3",
    "use_classifier_free_guidance": True,
    "cfg_scale": 2.5,
    "step_scale": 1.5,
    "num_timesteps": 200,
    "select_fixed_atoms": {"X1": "all", "L1": "all"},
    # Hotspots on N atoms + aromatic C's
    "select_hotspots": {"L1": "N1,N2,N3,N4,C4,C5,C15,C16,C26,C27"},
    "select_buried": {"L1": "all"},
    # Key: Tell RFD3 these N atoms should RECEIVE H-bonds from protein
    "select_hbond_donor": {"L1": "N1,N2,N3,N4"},
    # Negative Z to push protein BELOW the complex (where ligand N atoms are)
    "ori_token": [0.0, 8.0, -5.0]
}


def run_rfd3_batch(name: str, contig: str, num_designs: int, pdb_content: str) -> list:
    """Run a single RFD3 batch and return designs."""
    params = {**BASE_PARAMS, "contig": contig, "num_designs": num_designs, "pdb_content": pdb_content}

    payload = json.dumps({"input": params})

    print(f"Running {name} with params:")
    print(f"  contig: {contig}")
    print(f"  num_designs: {num_designs}")
    print(f"  select_hotspots: {params['select_hotspots']}")
    print(f"  select_hbond_donor: {params['select_hbond_donor']}")
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
                    # H-bond donors (NH, OH groups)
                    if atom_name in ['N', 'NE', 'NH1', 'NH2', 'ND1', 'ND2', 'NE1', 'NE2', 'NZ', 'OG', 'OG1', 'OH']:
                        protein_donors.append((coord, atom_name, res_name, res_num))
            except:
                continue

    result = {
        'hbonds': [],
        'ligand_min': 999,
        'z_offset': 0,
        'prot_atoms_4A': 0
    }

    if not ligand_N_coords or not protein_coords:
        return result

    protein_coords = np.array(protein_coords)

    # Check H-bonds to each ligand N
    for n_atom, n_coord in ligand_N_coords.items():
        for (pc, pa, pr, pn) in protein_donors:
            d = np.linalg.norm(n_coord - pc)
            if d < 3.5:  # H-bond distance
                result['hbonds'].append({
                    'ligand_atom': n_atom,
                    'protein': f"{pr}{pn}.{pa}",
                    'distance': d
                })

    # Overall stats
    if metal_coord is not None:
        result['z_offset'] = np.mean(protein_coords[:, 2]) - metal_coord[2]

    # Ligand contact stats (all UNL atoms)
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
    print(f"V9 TEST: {name}")
    print(f"Key change: Target ligand N atoms for H-bonding")
    print(f"  - Hotspots on N1,N2,N3,N4 + aromatic C's")
    print(f"  - select_hbond_donor on N atoms")
    print(f"  - Negative Z ori_token to approach from below")
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
        print("H-BOND ANALYSIS (protein -> ligand N atoms):")
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

            for hb in result['hbonds']:
                print(f"    {hb['ligand_atom']} <- {hb['protein']} ({hb['distance']:.2f}A)")

        print(f"\n{'='*60}")
        print(f"SUMMARY:")
        print(f"  Designs with H-bonds: {designs_with_hbonds}/{len(designs)}")
        print(f"  Total H-bonds: {total_hbonds}")
        print(f"  Avg H-bonds/design: {total_hbonds/len(designs):.1f}")
        print(f"{'='*60}")


if __name__ == "__main__":
    main()
