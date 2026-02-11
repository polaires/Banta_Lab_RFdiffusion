"""Check FULL coordination in v3 designs."""

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


# Analyze one design from each config
for pdb_file in [
    OUTPUT_DIR / "v3_medium_000.pdb",
    OUTPUT_DIR / "v3_large_000.pdb",
]:
    print(f"\n{'='*70}")
    print(f"Analyzing: {pdb_file.name}")
    print('='*70)

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
    protein_residues = {}
    for a in atoms:
        if a["res_name"] in ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                             "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
                             "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]:
            key = (a["chain"], a["res_num"])
            protein_residues[key] = a["res_name"]

    print(f"  Protein residues: {len(protein_residues)}")

    # TriNOx coordination to Dy
    trinox_oxygens = [a for a in atoms if a["res_name"] == "UNL" and a["name"] in ["O1", "O2", "O3"]]
    print(f"\n  TriNOx-Dy coordination:")
    seen = set()
    for a in trinox_oxygens:
        if a["name"] not in seen:
            d = distance(dy, a)
            status = "OK" if d < 3.0 else "BAD"
            print(f"    {a['name']}: {d:.2f} A [{status}]")
            seen.add(a["name"])

    # Asp coordination to Dy
    asp_od = [a for a in atoms if a["res_name"] == "ASP" and a["name"] in ["OD1", "OD2"]]
    print(f"\n  Asp-Dy coordination:")
    for a in asp_od[:2]:
        d = distance(dy, a)
        status = "OK" if d < 3.0 else "far"
        print(f"    ASP{a['res_num']} {a['name']}: {d:.2f} A [{status}]")

    # Protein atoms near TriNOx (pocket formation)
    trinox_all = [a for a in atoms if a["res_name"] == "UNL"]
    protein_atoms = [a for a in atoms if a["res_name"] in
                    ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                     "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
                     "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]]

    # Find residues that contact TriNOx
    contacting_residues = set()
    for ta in trinox_all:
        for pa in protein_atoms:
            d = distance(ta, pa)
            if d <= 4.0:
                contacting_residues.add((pa["chain"], pa["res_num"], pa["res_name"]))

    print(f"\n  Protein residues contacting TriNOx (<4A):")
    print(f"    {len(contacting_residues)} residues")
    for chain, res_num, res_name in sorted(contacting_residues)[:10]:
        print(f"      {res_name}{res_num}")
    if len(contacting_residues) > 10:
        print(f"      ... and {len(contacting_residues) - 10} more")
