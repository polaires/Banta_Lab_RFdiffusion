"""Analyze motif-scaffolded designs for metal coordination."""

import math
from pathlib import Path

OUTPUT_DIR = Path("G:/Github_local_repo/Banta_Lab_RFdiffusion/experiments/Dy_TriNOx_scaffold/outputs_v2")

def parse_pdb(pdb_content):
    atoms = []
    for line in pdb_content.split('\n'):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21]
                res_num = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms.append({
                    "name": atom_name,
                    "res_name": res_name,
                    "chain": chain,
                    "res_num": res_num,
                    "x": x, "y": y, "z": z
                })
            except:
                pass
    return atoms

def distance(a1, a2):
    return math.sqrt(
        (a1["x"] - a2["x"])**2 +
        (a1["y"] - a2["y"])**2 +
        (a1["z"] - a2["z"])**2
    )

# Analyze each design
for pdb_file in sorted(OUTPUT_DIR.glob("test_motif_scaffold_*.pdb")):
    print(f"\n{'='*60}")
    print(f"Analyzing: {pdb_file.name}")
    print('='*60)

    with open(pdb_file, 'r') as f:
        content = f.read()

    atoms = parse_pdb(content)

    # Find Dy
    dy = None
    for a in atoms:
        if a["res_name"] == "DY" or a["name"] == "DY":
            dy = a
            break

    if not dy:
        print("  No Dy found!")
        continue

    print(f"  Dy position: ({dy['x']:.2f}, {dy['y']:.2f}, {dy['z']:.2f})")

    # Categorize atoms
    ligand_atoms = [a for a in atoms if a["res_name"] == "UNL"]
    asp_atoms = [a for a in atoms if a["res_name"] == "ASP"]
    other_protein = [a for a in atoms if a["res_name"] not in ["DY", "UNL", "ASP"]
                     and a["res_name"] in ["ALA", "ARG", "ASN", "CYS", "GLN", "GLU",
                                           "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
                                           "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]]

    print(f"  Ligand atoms: {len(ligand_atoms)}")
    print(f"  Asp atoms: {len(asp_atoms)}")
    print(f"  Other protein atoms: {len(other_protein)}")

    # Find Asp OD1/OD2 distances to Dy
    print(f"\n  Asp coordination:")
    for a in asp_atoms:
        if a["name"] in ["OD1", "OD2"]:
            d = distance(dy, a)
            status = "COORDINATING!" if d < 3.0 else "too far"
            print(f"    {a['res_name']}{a['res_num']} {a['name']}: {d:.2f} A - {status}")

    # Find closest protein atom to Dy
    all_protein = asp_atoms + other_protein
    if all_protein:
        closest = min(all_protein, key=lambda a: distance(dy, a))
        closest_dist = distance(dy, closest)
        print(f"\n  Closest protein atom to Dy:")
        print(f"    {closest['res_name']}{closest['res_num']} {closest['name']}: {closest_dist:.2f} A")

    # Count atoms within coordination shell
    within_3 = [a for a in all_protein if distance(dy, a) <= 3.0]
    within_5 = [a for a in all_protein if distance(dy, a) <= 5.0]
    print(f"\n  Atoms within 3A: {len(within_3)}")
    print(f"  Atoms within 5A: {len(within_5)}")

    if within_3:
        print(f"  Coordinating atoms:")
        for a in within_3:
            print(f"    {a['res_name']}{a['res_num']} {a['name']}: {distance(dy, a):.2f} A")
