"""Verify FULL coordination in v3 designs - find the ACTUAL coordinating Asp."""

import math
from pathlib import Path

OUTPUT_DIR = Path("G:/Github_local_repo/Banta_Lab_RFdiffusion/experiments/Dy_TriNOx_scaffold/outputs_v3")


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


print("=" * 80)
print("V3 DESIGN COORDINATION VERIFICATION")
print("=" * 80)

# Analyze all v3 designs
for pdb_file in sorted(OUTPUT_DIR.glob("v3_*.pdb")):
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
        print(f"{pdb_file.stem}: NO DY FOUND")
        continue

    # TriNOx oxygens
    trinox_o = [a for a in atoms if a["res_name"] == "UNL" and a["name"] in ["O1", "O2", "O3"]]
    trinox_dists = []
    seen = set()
    for a in trinox_o:
        if a["name"] not in seen:
            trinox_dists.append(distance(dy, a))
            seen.add(a["name"])

    trinox_ok = all(d < 3.0 for d in trinox_dists[:3])

    # Find ALL Asp OD1/OD2 and check which one coordinates
    asp_atoms = [a for a in atoms if a["res_name"] == "ASP" and a["name"] in ["OD1", "OD2"]]
    coord_asp = None
    min_asp_dist = 999
    for a in asp_atoms:
        d = distance(dy, a)
        if d < min_asp_dist:
            min_asp_dist = d
            coord_asp = a

    asp_ok = min_asp_dist < 3.0

    # Count protein residues
    protein_res = set()
    for a in atoms:
        if a["res_name"] in ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                             "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
                             "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]:
            protein_res.add((a["chain"], a["res_num"]))

    # Protein contacts to TriNOx
    trinox_all = [a for a in atoms if a["res_name"] == "UNL"]
    protein_atoms = [a for a in atoms if a["res_name"] in
                    ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                     "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
                     "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]]

    contact_residues = set()
    for ta in trinox_all:
        for pa in protein_atoms:
            if distance(ta, pa) <= 4.0:
                contact_residues.add((pa["chain"], pa["res_num"]))

    status = "OK" if (trinox_ok and asp_ok and len(contact_residues) >= 5) else "CHECK"

    print(f"{pdb_file.stem}: {len(protein_res)} res, "
          f"TriNOx-Dy={trinox_dists[0]:.2f}A [{['BAD','OK'][trinox_ok]}], "
          f"Asp{coord_asp['res_num'] if coord_asp else '?'}-Dy={min_asp_dist:.2f}A [{['BAD','OK'][asp_ok]}], "
          f"pocket={len(contact_residues)} res [{status}]")

print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
