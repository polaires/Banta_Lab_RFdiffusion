#!/usr/bin/env python3
"""
Analyze Round 5 designs for TB-citrate coordination preservation
and protein coordination around the metal-citrate complex.
"""

import json
import math
from pathlib import Path
from collections import defaultdict

EXPERIMENT_DIR = Path(__file__).parent.parent
R5_DIR = EXPERIMENT_DIR / "outputs" / "round_05"


def parse_pdb(filepath: Path) -> tuple:
    """Parse PDB file to extract atoms."""
    protein_atoms = []
    ligand_atoms = []
    metal_atoms = []

    with open(filepath) as f:
        for line in f:
            if line.startswith("ATOM"):
                atom = parse_atom_line(line)
                if atom:
                    protein_atoms.append(atom)
            elif line.startswith("HETATM"):
                atom = parse_atom_line(line)
                if atom:
                    res_name = atom.get("res_name", "").strip()
                    if res_name == "TB":
                        metal_atoms.append(atom)
                    else:
                        ligand_atoms.append(atom)

    return protein_atoms, ligand_atoms, metal_atoms


def parse_atom_line(line: str) -> dict:
    """Parse a single ATOM/HETATM line."""
    try:
        return {
            "record": line[0:6].strip(),
            "serial": int(line[6:11].strip()),
            "name": line[12:16].strip(),
            "res_name": line[17:20].strip(),
            "chain": line[21].strip(),
            "res_seq": line[22:26].strip(),
            "x": float(line[30:38].strip()),
            "y": float(line[38:46].strip()),
            "z": float(line[46:54].strip()),
        }
    except (ValueError, IndexError):
        return None


def distance(a1: dict, a2: dict) -> float:
    """Calculate 3D distance between two atoms."""
    return math.sqrt(
        (a1["x"] - a2["x"])**2 +
        (a1["y"] - a2["y"])**2 +
        (a1["z"] - a2["z"])**2
    )


def find_coordinating_residues(protein_atoms: list, center: dict, cutoff: float = 5.0) -> list:
    """Find residues with coordinating atoms near a center point."""
    coordinating = []

    # Group atoms by residue
    residues = defaultdict(list)
    for atom in protein_atoms:
        key = (atom["res_name"], atom["res_seq"])
        residues[key].append(atom)

    # Check each residue for coordinating atoms
    for (res_name, res_seq), atoms in residues.items():
        # Define coordinating atoms for each residue type
        coord_atoms = {
            "ASP": ["OD1", "OD2"],
            "GLU": ["OE1", "OE2"],
            "ASN": ["OD1"],
            "GLN": ["OE1"],
            "HIS": ["ND1", "NE2"],
            "SER": ["OG"],
            "THR": ["OG1"],
        }

        if res_name in coord_atoms:
            for atom in atoms:
                if atom["name"] in coord_atoms[res_name]:
                    dist = distance(atom, center)
                    if dist <= cutoff:
                        coordinating.append({
                            "res_name": res_name,
                            "res_seq": res_seq,
                            "atom": atom["name"],
                            "dist": round(dist, 2)
                        })

    return sorted(coordinating, key=lambda x: x["dist"])


def analyze_design(filepath: Path) -> dict:
    """Analyze a single design."""
    result = {"name": filepath.stem}

    protein_atoms, ligand_atoms, metal_atoms = parse_pdb(filepath)

    if not metal_atoms:
        result["error"] = "No TB metal found"
        return result

    if not ligand_atoms:
        result["error"] = "No citrate found"
        return result

    # Get TB position
    tb = metal_atoms[0]
    result["tb_position"] = (round(tb["x"], 2), round(tb["y"], 2), round(tb["z"], 2))

    # Get citrate oxygens and calculate distances to TB
    cit_oxygens = [a for a in ligand_atoms if a["name"].startswith("O")]
    tb_o_distances = []
    for o in cit_oxygens:
        dist = distance(tb, o)
        tb_o_distances.append({
            "atom": o["name"],
            "dist": round(dist, 2)
        })

    tb_o_distances.sort(key=lambda x: x["dist"])
    result["tb_citrate_distances"] = tb_o_distances
    result["min_tb_o_dist"] = tb_o_distances[0]["dist"] if tb_o_distances else None
    result["max_tb_o_dist"] = tb_o_distances[-1]["dist"] if tb_o_distances else None

    # Count residues
    res_ids = set()
    for atom in protein_atoms:
        res_ids.add(atom["res_seq"])
    result["n_residues"] = len(res_ids)

    # Find protein residues coordinating TB
    tb_coord = find_coordinating_residues(protein_atoms, tb, cutoff=4.0)
    result["tb_coordination"] = tb_coord
    result["tb_CN"] = len(tb_coord)

    # Find protein residues near citrate center
    cit_center = {
        "x": sum(a["x"] for a in ligand_atoms) / len(ligand_atoms),
        "y": sum(a["y"] for a in ligand_atoms) / len(ligand_atoms),
        "z": sum(a["z"] for a in ligand_atoms) / len(ligand_atoms),
    }
    cit_coord = find_coordinating_residues(protein_atoms, cit_center, cutoff=5.0)
    result["cit_nearby"] = cit_coord
    result["cit_CN"] = len(cit_coord)

    # Classify by coordinating residue types
    coord_types = defaultdict(int)
    for res in tb_coord:
        coord_types[res["res_name"]] += 1
    result["coord_types"] = dict(coord_types)

    return result


def main():
    print("=" * 60)
    print("Round 5 Analysis: TB-Citrate Coordination Verification")
    print("=" * 60)

    # Find all designs
    pdb_files = sorted(R5_DIR.glob("*.pdb"))
    if not pdb_files:
        print(f"No PDB files found in {R5_DIR}")
        return

    print(f"\nFound {len(pdb_files)} designs\n")

    all_results = []
    tb_o_distances = []
    cn_values = []

    for pdb_file in pdb_files:
        result = analyze_design(pdb_file)
        all_results.append(result)

        if "error" not in result:
            tb_o_distances.append(result["min_tb_o_dist"])
            cn_values.append(result["tb_CN"])

    # Summary table
    print(f"{'Design':<30} {'TB-O min':>8} {'TB-O max':>8} {'TB CN':>6} {'CIT CN':>7}")
    print("-" * 65)

    for r in all_results:
        name = r["name"][:28]
        if "error" in r:
            print(f"{name:<30} ERROR: {r['error']}")
        else:
            tb_min = f"{r['min_tb_o_dist']:.2f}" if r.get("min_tb_o_dist") else "N/A"
            tb_max = f"{r['max_tb_o_dist']:.2f}" if r.get("max_tb_o_dist") else "N/A"
            tb_cn = str(r.get("tb_CN", 0))
            cit_cn = str(r.get("cit_CN", 0))
            print(f"{name:<30} {tb_min:>8} {tb_max:>8} {tb_cn:>6} {cit_cn:>7}")

    # Statistics
    if tb_o_distances:
        print(f"\n{'=' * 60}")
        print("Summary Statistics")
        print(f"{'=' * 60}")
        print(f"\nTB-Citrate O distances:")
        print(f"  Min: {min(tb_o_distances):.2f} A")
        print(f"  Max: {max(tb_o_distances):.2f} A")
        print(f"  Mean: {sum(tb_o_distances)/len(tb_o_distances):.2f} A")

        print(f"\nProtein coordination (TB within 4A):")
        print(f"  Min CN: {min(cn_values)}")
        print(f"  Max CN: {max(cn_values)}")
        print(f"  Mean CN: {sum(cn_values)/len(cn_values):.1f}")

        # Target: TB-O should be 2.3-2.7 A (preserved from input)
        good_dist = sum(1 for d in tb_o_distances if 2.0 <= d <= 3.0)
        print(f"\n  Designs with TB-O in 2-3 A range: {good_dist}/{len(tb_o_distances)}")

        # Count designs with protein coordination
        good_cn = sum(1 for cn in cn_values if cn >= 1)
        print(f"  Designs with protein CN >= 1: {good_cn}/{len(cn_values)}")

    # Best designs (highest protein CN)
    print(f"\n{'=' * 60}")
    print("Top Designs (by protein coordination number)")
    print(f"{'=' * 60}")

    sorted_results = sorted(
        [r for r in all_results if "error" not in r],
        key=lambda x: x.get("tb_CN", 0),
        reverse=True
    )

    for r in sorted_results[:10]:
        print(f"\n{r['name']}:")
        print(f"  TB-O range: {r['min_tb_o_dist']:.2f} - {r['max_tb_o_dist']:.2f} A")
        print(f"  Protein TB CN: {r['tb_CN']}")
        if r.get("tb_coordination"):
            coord_str = ", ".join(f"{c['res_name']}{c['res_seq']}({c['dist']:.1f}A)"
                                  for c in r["tb_coordination"][:5])
            print(f"  Coordinating: {coord_str}")
        if r.get("coord_types"):
            print(f"  Types: {r['coord_types']}")

    # Save results
    summary = {
        "round": 5,
        "n_designs": len(all_results),
        "n_successful": len([r for r in all_results if "error" not in r]),
        "tb_o_distance_stats": {
            "min": min(tb_o_distances) if tb_o_distances else None,
            "max": max(tb_o_distances) if tb_o_distances else None,
            "mean": sum(tb_o_distances)/len(tb_o_distances) if tb_o_distances else None,
        },
        "protein_cn_stats": {
            "min": min(cn_values) if cn_values else None,
            "max": max(cn_values) if cn_values else None,
            "mean": sum(cn_values)/len(cn_values) if cn_values else None,
        },
        "designs": all_results
    }

    summary_path = R5_DIR / "r5_analysis.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"\n\nAnalysis saved to: {summary_path}")


if __name__ == "__main__":
    main()
