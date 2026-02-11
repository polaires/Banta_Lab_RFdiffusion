"""
Analyze Dy-TriNOx scaffold designs for metal coordination.
Check which residues are near the Dy center and could coordinate.
"""

import json
from pathlib import Path
import math

OUTPUT_DIR = Path("G:/Github_local_repo/Banta_Lab_RFdiffusion/experiments/Dy_TriNOx_scaffold/outputs")
INPUT_DIR = Path("G:/Github_local_repo/Banta_Lab_RFdiffusion/experiments/Dy_TriNOx_scaffold/inputs")

# Coordinating atoms by residue type
COORD_ATOMS = {
    "ASP": ["OD1", "OD2"],      # Carboxylate oxygens
    "GLU": ["OE1", "OE2"],      # Carboxylate oxygens
    "HIS": ["ND1", "NE2"],      # Imidazole nitrogens
    "CYS": ["SG"],              # Thiolate sulfur
    "ASN": ["OD1"],             # Amide oxygen
    "GLN": ["OE1"],             # Amide oxygen
    "SER": ["OG"],              # Hydroxyl oxygen
    "THR": ["OG1"],             # Hydroxyl oxygen
    "TYR": ["OH"],              # Phenolic oxygen
    "MET": ["SD"],              # Thioether sulfur
}

# Typical coordination distances for Ln3+ (in Angstroms)
COORD_DISTANCE_MAX = 3.0  # Max distance for direct coordination
NEAR_DISTANCE_MAX = 5.0   # Max distance for potential secondary shell


def parse_pdb(pdb_path):
    """Parse PDB file and extract atom coordinates."""
    atoms = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
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
    return atoms


def distance(a1, a2):
    """Calculate distance between two atoms."""
    return math.sqrt(
        (a1["x"] - a2["x"])**2 +
        (a1["y"] - a2["y"])**2 +
        (a1["z"] - a2["z"])**2
    )


def find_dy_atom(atoms):
    """Find the Dy atom."""
    for a in atoms:
        if a["res_name"] == "DY" or a["name"] == "DY":
            return a
    return None


def analyze_coordination(pdb_path):
    """Analyze coordination environment around Dy."""
    atoms = parse_pdb(pdb_path)
    dy = find_dy_atom(atoms)

    if not dy:
        return {"error": "No Dy found"}

    # Find protein atoms near Dy
    coordinating = []
    near_shell = []

    for a in atoms:
        # Skip non-protein atoms
        if a["chain"] in ["X", "L"]:
            continue
        if a["res_name"] in ["DY", "UNL"]:
            continue

        d = distance(dy, a)

        # Check if this is a potential coordinating atom
        is_coord_atom = (a["res_name"] in COORD_ATOMS and
                        a["name"] in COORD_ATOMS.get(a["res_name"], []))

        if d <= COORD_DISTANCE_MAX:
            coordinating.append({
                "res": f"{a['res_name']}{a['res_num']}",
                "atom": a["name"],
                "distance": round(d, 2),
                "is_coord_atom": is_coord_atom
            })
        elif d <= NEAR_DISTANCE_MAX and is_coord_atom:
            near_shell.append({
                "res": f"{a['res_name']}{a['res_num']}",
                "atom": a["name"],
                "distance": round(d, 2)
            })

    # Also check backbone atoms (N, O) within coordination distance
    backbone_near = []
    for a in atoms:
        if a["chain"] in ["X", "L"]:
            continue
        if a["name"] in ["N", "O"] and a["res_name"] not in ["DY", "UNL"]:
            d = distance(dy, a)
            if d <= COORD_DISTANCE_MAX:
                backbone_near.append({
                    "res": f"{a['res_name']}{a['res_num']}",
                    "atom": a["name"],
                    "distance": round(d, 2)
                })

    return {
        "coordinating": coordinating,
        "near_shell": near_shell,
        "backbone_near": backbone_near,
        "dy_coords": [dy["x"], dy["y"], dy["z"]]
    }


def main():
    print("=" * 70)
    print("DY-TRINOX COORDINATION ANALYSIS")
    print("=" * 70)

    # Find all design PDBs
    pdb_files = sorted([f for f in OUTPUT_DIR.glob("r1*.pdb")])
    print(f"\nAnalyzing {len(pdb_files)} backbone designs...")

    results = []

    for pdb_file in pdb_files:
        print(f"\n{pdb_file.stem}:")

        analysis = analyze_coordination(pdb_file)

        if "error" in analysis:
            print(f"  ERROR: {analysis['error']}")
            continue

        coord = analysis["coordinating"]
        near = analysis["near_shell"]
        backbone = analysis["backbone_near"]

        # Count true coordinating atoms
        true_coord = [c for c in coord if c["is_coord_atom"]]

        print(f"  Dy position: ({analysis['dy_coords'][0]:.1f}, {analysis['dy_coords'][1]:.1f}, {analysis['dy_coords'][2]:.1f})")
        print(f"  Atoms within 3A of Dy: {len(coord)}")
        print(f"  Potential coord atoms within 3A: {len(true_coord)}")
        print(f"  Coord atoms in 3-5A shell: {len(near)}")
        print(f"  Backbone atoms within 3A: {len(backbone)}")

        if true_coord:
            print(f"  COORDINATING RESIDUES:")
            for c in true_coord:
                print(f"    {c['res']} {c['atom']}: {c['distance']} A")

        if near:
            print(f"  NEAR SHELL (3-5A):")
            for n in near:
                print(f"    {n['res']} {n['atom']}: {n['distance']} A")

        results.append({
            "name": pdb_file.stem,
            "coordinating_atoms": len(true_coord),
            "near_shell_atoms": len(near),
            "backbone_near": len(backbone),
            "details": analysis
        })

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    total_coord = sum(r["coordinating_atoms"] for r in results)
    designs_with_coord = sum(1 for r in results if r["coordinating_atoms"] > 0)

    print(f"Total designs: {len(results)}")
    print(f"Designs with coordinating atoms (<3A): {designs_with_coord}")
    print(f"Total coordinating atoms found: {total_coord}")

    if designs_with_coord == 0:
        print("\nWARNING: No designs have direct metal coordination!")
        print("The protein scaffolds may not be properly positioned around the Dy center.")
        print("\nPossible issues:")
        print("  1. RFD3 may not be placing the protein close enough to the metal")
        print("  2. The burial conditioning may not be forcing close contact")
        print("  3. Need to add explicit hotspot or H-bond conditioning")

    # Save analysis
    analysis_file = OUTPUT_DIR / "coordination_analysis.json"
    with open(analysis_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nFull analysis saved to {analysis_file}")


if __name__ == "__main__":
    main()
