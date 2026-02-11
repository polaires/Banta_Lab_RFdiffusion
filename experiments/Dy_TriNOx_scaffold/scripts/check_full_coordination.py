"""Check FULL coordination: TriNOx + Asp to Dy."""

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


# Analyze each v2 design
for pdb_file in sorted(OUTPUT_DIR.glob("v2_*.pdb"))[:3]:  # Just first 3
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

    # Count protein residues
    protein_residues = set()
    for a in atoms:
        if a["res_name"] in ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                             "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
                             "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]:
            protein_residues.add((a["chain"], a["res_num"], a["res_name"]))

    print(f"  Protein residues: {len(protein_residues)}")

    # Find TriNOx oxygens (O1, O2, O3) - these should coordinate Dy
    trinox_oxygens = [a for a in atoms if a["res_name"] == "UNL" and a["name"] in ["O1", "O2", "O3"]]
    print(f"\n  TriNOx coordination (O1, O2, O3 to Dy):")
    for a in trinox_oxygens[:3]:  # Only first 3 (avoid duplicates)
        d = distance(dy, a)
        status = "OK" if d < 3.0 else "BAD"
        print(f"    {a['name']}: {d:.2f} A [{status}]")

    # Find Asp coordination
    asp_atoms = [a for a in atoms if a["res_name"] == "ASP" and a["name"] in ["OD1", "OD2"]]
    print(f"\n  Asp coordination (OD1/OD2 to Dy):")
    for a in asp_atoms[:2]:  # First Asp only
        d = distance(dy, a)
        status = "OK" if d < 3.0 else "BAD"
        print(f"    {a['res_name']}{a['res_num']} {a['name']}: {d:.2f} A [{status}]")
