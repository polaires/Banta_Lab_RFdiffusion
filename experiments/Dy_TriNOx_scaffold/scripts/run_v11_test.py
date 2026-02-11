#!/usr/bin/env python3
"""Run v11 test: Optimize hydrophobic packing around TriNOx aromatic rings and tert-butyls."""

import json
import subprocess
import os
from pathlib import Path
import numpy as np

# Configuration
BASE_DIR = Path("G:/Github_local_repo/Banta_Lab_RFdiffusion/experiments/Dy_TriNOx_scaffold")
OUTPUT_DIR = BASE_DIR / "outputs_v11"
INPUT_PDB = BASE_DIR / "inputs" / "Dy_TriNOx_UNL.pdb"
API_URL = "http://localhost:8000/runsync"

# V11 test config
TEST_CONFIG = {
    "name": "v11_hydrophobic",
    "contig": "X1,L1,/0,120-140",
    "num_designs": 10
}

# V11 base params - hydrophobic packing focus
# TriNOx is mostly hydrophobic: aromatic rings + tert-butyls
# NO H-bond conditioning - the N-OH groups face the metal, not outward
BASE_PARAMS = {
    "task": "rfd3",
    "use_classifier_free_guidance": True,
    "cfg_scale": 2.5,
    "step_scale": 1.5,
    "num_timesteps": 200,
    "select_fixed_atoms": {"X1": "all", "L1": "all"},
    # Hotspots on ALL hydrophobic carbons (aromatics + tert-butyls)
    "select_hotspots": {"L1": "C4,C5,C6,C7,C8,C9,C10,C15,C16,C17,C18,C19,C20,C21,C26,C27,C28,C29,C30,C31,C32"},
    # Bury the entire ligand for tight packing
    "select_buried": {"L1": "all"}
    # NO hbond conditioning - this is a hydrophobic ligand
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
    print(f"  cfg_scale: {params['cfg_scale']}")
    print(f"  (NO H-bond conditioning - TriNOx is hydrophobic)")

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


def analyze_hydrophobic_contacts(pdb_content: str) -> dict:
    """Analyze hydrophobic contacts to TriNOx aromatic rings and tert-butyls."""
    # Aromatic carbons (3 rings)
    aromatic_c = ['C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C26', 'C27', 'C28', 'C29', 'C30', 'C31']
    # Tert-butyl carbons
    tbutyl_c = ['C10', 'C21', 'C32']

    ligand_coords = {}
    protein_hydrophobic = []  # (coord, atom_name, res_name, res_num)
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
                    elif res_name == 'UNL':
                        ligand_coords[atom_name] = coord
                elif line.startswith('ATOM'):
                    protein_coords.append(coord)
                    # Hydrophobic atoms (C atoms from hydrophobic residues)
                    if res_name in ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'TYR', 'PRO'] and atom_name.startswith('C'):
                        protein_hydrophobic.append((coord, atom_name, res_name, res_num))
            except:
                continue

    result = {
        'aromatic_contacts': 0,
        'tbutyl_contacts': 0,
        'total_contacts_4A': 0,
        'ligand_min': 999,
        'z_offset': 0,
        'contact_residues': set()
    }

    if not ligand_coords or not protein_coords:
        return result

    protein_coords = np.array(protein_coords)

    # Count contacts to aromatic carbons
    for c_name in aromatic_c:
        if c_name in ligand_coords:
            lc = ligand_coords[c_name]
            for (pc, pa, pr, pn) in protein_hydrophobic:
                d = np.linalg.norm(lc - pc)
                if d < 4.5:  # Hydrophobic contact distance
                    result['aromatic_contacts'] += 1
                    result['contact_residues'].add(f"{pr}{pn}")

    # Count contacts to tert-butyl carbons
    for c_name in tbutyl_c:
        if c_name in ligand_coords:
            lc = ligand_coords[c_name]
            for (pc, pa, pr, pn) in protein_hydrophobic:
                d = np.linalg.norm(lc - pc)
                if d < 4.5:
                    result['tbutyl_contacts'] += 1
                    result['contact_residues'].add(f"{pr}{pn}")

    # Overall ligand contacts
    all_lig = np.array(list(ligand_coords.values()))
    min_dists = []
    for pc in protein_coords:
        d = np.linalg.norm(all_lig - pc, axis=1)
        min_dists.append(np.min(d))
    min_dists = np.array(min_dists)

    result['ligand_min'] = np.min(min_dists)
    result['total_contacts_4A'] = np.sum(min_dists < 4.0)

    if metal_coord is not None:
        result['z_offset'] = np.mean(protein_coords[:, 2]) - metal_coord[2]

    return result


def main():
    import time

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    with open(INPUT_PDB) as f:
        pdb_content = f.read()

    name = TEST_CONFIG["name"]
    contig = TEST_CONFIG["contig"]
    num_designs = TEST_CONFIG["num_designs"]

    print(f"\n{'='*60}")
    print(f"V11 TEST: {name}")
    print(f"Strategy: HYDROPHOBIC PACKING")
    print(f"  - TriNOx is mostly hydrophobic (aromatics + tert-butyls)")
    print(f"  - Hotspots on ALL hydrophobic carbons")
    print(f"  - select_buried for tight packing")
    print(f"  - NO H-bond conditioning")
    print(f"{'='*60}")

    start_time = time.time()
    designs = run_rfd3_batch(name, contig, num_designs, pdb_content)
    elapsed = time.time() - start_time

    for i, design in enumerate(designs):
        filename = OUTPUT_DIR / f"{name}_{i:03d}.pdb"
        with open(filename, "w") as f:
            f.write(design["content"])

    print(f"\n  Generated {len(designs)} designs in {elapsed:.2f} seconds")
    print(f"  Saved to {OUTPUT_DIR}")

    if designs:
        print(f"\n{'='*60}")
        print("HYDROPHOBIC CONTACT ANALYSIS:")
        print(f"{'='*60}")

        for i, design in enumerate(designs):
            r = analyze_hydrophobic_contacts(design["content"])

            total_hydro = r['aromatic_contacts'] + r['tbutyl_contacts']
            status = 'EXCELLENT' if total_hydro >= 20 else 'GOOD' if total_hydro >= 10 else 'OK' if total_hydro >= 5 else 'LOW'

            print(f"\nDesign {i:03d}: [{status}]")
            print(f"  Ligand min dist: {r['ligand_min']:.1f}A")
            print(f"  Total <4A contacts: {r['total_contacts_4A']}")
            print(f"  Aromatic contacts: {r['aromatic_contacts']}")
            print(f"  Tert-butyl contacts: {r['tbutyl_contacts']}")
            print(f"  Hydrophobic residues: {len(r['contact_residues'])}")
            print(f"  Z-offset: {r['z_offset']:+.1f}A")

        print(f"\n{'='*60}")
        print(f"SUMMARY:")
        avg_aromatic = np.mean([analyze_hydrophobic_contacts(d["content"])['aromatic_contacts'] for d in designs])
        avg_tbutyl = np.mean([analyze_hydrophobic_contacts(d["content"])['tbutyl_contacts'] for d in designs])
        print(f"  Avg aromatic contacts: {avg_aromatic:.1f}")
        print(f"  Avg tert-butyl contacts: {avg_tbutyl:.1f}")
        print(f"{'='*60}")


if __name__ == "__main__":
    main()
