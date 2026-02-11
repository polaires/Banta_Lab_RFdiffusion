#!/usr/bin/env python3
"""Check distances from protein to Dy metal in v4 designs."""

import numpy as np
from pathlib import Path

OUTPUT_DIR = Path("/mnt/g/Github_local_repo/Banta_Lab_RFdiffusion/experiments/Dy_TriNOx_scaffold/outputs_v4")

def parse_pdb_coords(pdb_content):
    """Extract coordinates from PDB content."""
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

    return metal_coord, np.array(ligand_coords) if ligand_coords else None, np.array(protein_coords) if protein_coords else None


def main():
    design_files = sorted(OUTPUT_DIR.glob("v4_small_*.pdb"))[:5]  # Check first 5

    print("Distance analysis for v4_small designs:")
    print("=" * 60)

    for pdb_file in design_files:
        with open(pdb_file) as f:
            pdb_content = f.read()

        metal_coord, ligand_coords, protein_coords = parse_pdb_coords(pdb_content)

        if metal_coord is None:
            print(f"{pdb_file.stem}: No metal found!")
            continue

        if protein_coords is None or len(protein_coords) == 0:
            print(f"{pdb_file.stem}: No protein atoms found!")
            continue

        # Calculate distances from metal to protein atoms
        metal_to_protein = np.linalg.norm(protein_coords - metal_coord, axis=1)
        min_dist = np.min(metal_to_protein)
        mean_dist = np.mean(metal_to_protein)

        # Calculate distances from ligand center to protein
        if ligand_coords is not None:
            ligand_center = np.mean(ligand_coords, axis=0)
            ligand_to_protein = np.linalg.norm(protein_coords - ligand_center, axis=1)
            ligand_min_dist = np.min(ligand_to_protein)
        else:
            ligand_min_dist = None

        # Count atoms within coordination distance
        close_atoms = np.sum(metal_to_protein < 3.0)  # Typical Dy coordination distance
        semi_close = np.sum(metal_to_protein < 5.0)

        print(f"{pdb_file.stem}:")
        print(f"  Metal at: {metal_coord}")
        print(f"  Protein atoms: {len(protein_coords)}")
        print(f"  Min distance to metal: {min_dist:.2f} Å")
        print(f"  Atoms within 3Å: {close_atoms}, within 5Å: {semi_close}")
        if ligand_min_dist:
            print(f"  Min distance to ligand center: {ligand_min_dist:.2f} Å")
        print()


if __name__ == "__main__":
    main()
