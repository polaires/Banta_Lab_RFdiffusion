#!/usr/bin/env python3
"""
Analyze Round 5 RF3 predictions for coordination and fold quality.
"""

import json
import math
from pathlib import Path
from collections import defaultdict

EXPERIMENT_DIR = Path(__file__).parent.parent
RF3_DIR = EXPERIMENT_DIR / "outputs" / "round_05_rf3"


def parse_mmcif_atoms(filepath: Path) -> tuple:
    """Parse mmCIF format file."""
    protein_atoms = []
    ligand_atoms = []
    metal_atoms = []

    in_atom_site = False
    header_fields = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()

            if line.startswith("_atom_site."):
                in_atom_site = True
                field = line.split(".")[1].strip()
                header_fields.append(field)
                continue

            if in_atom_site and (line.startswith("#") or line.startswith("_") or line.startswith("loop_")):
                in_atom_site = False
                continue

            if in_atom_site and (line.startswith("ATOM") or line.startswith("HETATM")):
                parts = line.split()
                if len(parts) < len(header_fields):
                    continue

                atom = {}
                for i, field in enumerate(header_fields):
                    if i < len(parts):
                        atom[field] = parts[i]

                try:
                    atom["x"] = float(atom.get("Cartn_x", 0))
                    atom["y"] = float(atom.get("Cartn_y", 0))
                    atom["z"] = float(atom.get("Cartn_z", 0))
                    atom["plddt"] = float(atom.get("B_iso_or_equiv", 0))
                except ValueError:
                    continue

                group = atom.get("group_PDB", "")
                comp_id = atom.get("label_comp_id", "")

                metal_names = ["TB", "GD", "EU", "LA", "CE", "MG", "CA", "ZN", "FE"]
                if comp_id.upper() in metal_names or atom.get("type_symbol", "").upper() in metal_names:
                    metal_atoms.append(atom)
                elif group == "HETATM":
                    ligand_atoms.append(atom)
                else:
                    protein_atoms.append(atom)

    return protein_atoms, ligand_atoms, metal_atoms


def distance(a1: dict, a2: dict) -> float:
    return math.sqrt((a1["x"] - a2["x"])**2 + (a1["y"] - a2["y"])**2 + (a1["z"] - a2["z"])**2)


def find_coordinating_residues(protein_atoms: list, center: dict, cutoff: float = 5.0) -> list:
    """Find residues with coordinating atoms near center."""
    coordinating = []
    residues = defaultdict(list)

    for atom in protein_atoms:
        key = (atom.get("label_comp_id", ""), atom.get("label_seq_id", ""))
        residues[key].append(atom)

    coord_atoms = {
        "ASP": ["OD1", "OD2"],
        "GLU": ["OE1", "OE2"],
        "ASN": ["OD1"],
        "GLN": ["OE1"],
        "HIS": ["ND1", "NE2"],
        "SER": ["OG"],
        "THR": ["OG1"],
    }

    for (res_name, res_seq), atoms in residues.items():
        if res_name in coord_atoms:
            for atom in atoms:
                if atom.get("label_atom_id", "") in coord_atoms[res_name]:
                    dist = distance(atom, center)
                    if dist <= cutoff:
                        coordinating.append({
                            "res": f"{res_name}{res_seq}",
                            "atom": atom.get("label_atom_id"),
                            "dist": round(dist, 2)
                        })

    return sorted(coordinating, key=lambda x: x["dist"])


def analyze_rf3(filepath: Path) -> dict:
    """Analyze a single RF3 prediction."""
    result = {"name": filepath.stem}

    protein_atoms, ligand_atoms, metal_atoms = parse_mmcif_atoms(filepath)

    if not protein_atoms:
        result["error"] = "No protein atoms"
        return result

    # Protein metrics
    plddt_vals = [a["plddt"] for a in protein_atoms if a["plddt"] > 0]
    if plddt_vals:
        result["mean_plddt"] = round(sum(plddt_vals) / len(plddt_vals) * 100, 1)
        result["min_plddt"] = round(min(plddt_vals) * 100, 1)

    res_ids = set(a.get("label_seq_id") for a in protein_atoms)
    result["n_residues"] = len(res_ids)

    # Ligand metrics
    result["has_ligand"] = len(ligand_atoms) > 0
    if ligand_atoms:
        lig_plddt = [a["plddt"] for a in ligand_atoms if a["plddt"] > 0]
        if lig_plddt:
            result["ligand_plddt"] = round(sum(lig_plddt) / len(lig_plddt) * 100, 1)

        # Ligand center
        lig_center = {
            "x": sum(a["x"] for a in ligand_atoms) / len(ligand_atoms),
            "y": sum(a["y"] for a in ligand_atoms) / len(ligand_atoms),
            "z": sum(a["z"] for a in ligand_atoms) / len(ligand_atoms),
        }

        # Find coordinating residues
        coord_4A = find_coordinating_residues(protein_atoms, lig_center, 4.0)
        coord_6A = find_coordinating_residues(protein_atoms, lig_center, 6.0)

        result["CN_4A"] = len(coord_4A)
        result["CN_6A"] = len(coord_6A)
        result["coord_residues"] = coord_4A[:10]

        # Count by type
        types = defaultdict(int)
        for c in coord_4A:
            res_name = c["res"][:3]
            types[res_name] += 1
        result["coord_types"] = dict(types)

    return result


def main():
    print("=" * 60)
    print("Round 5 RF3 Analysis: Fold Quality & Coordination")
    print("=" * 60)

    pdb_files = sorted(RF3_DIR.glob("*_rf3.pdb"))
    if not pdb_files:
        print(f"No PDB files found in {RF3_DIR}")
        return

    print(f"\nFound {len(pdb_files)} predictions\n")

    all_results = []

    print(f"{'Name':<35} {'pLDDT':>7} {'Lig':>5} {'CN4':>4} {'CN6':>4} {'Coord Types'}")
    print("-" * 75)

    for pdb_file in pdb_files:
        r = analyze_rf3(pdb_file)
        all_results.append(r)

        if "error" in r:
            print(f"{r['name']:<35} ERROR: {r['error']}")
        else:
            name = r["name"][:33]
            plddt = f"{r.get('mean_plddt', 0):.1f}%"
            lig = f"{r.get('ligand_plddt', 0):.0f}%" if r.get("ligand_plddt") else "No"
            cn4 = str(r.get("CN_4A", 0))
            cn6 = str(r.get("CN_6A", 0))
            types = str(r.get("coord_types", {})) if r.get("coord_types") else "-"
            print(f"{name:<35} {plddt:>7} {lig:>5} {cn4:>4} {cn6:>4} {types}")

    # Summary stats
    successful = [r for r in all_results if "error" not in r]
    if successful:
        plddts = [r["mean_plddt"] for r in successful if r.get("mean_plddt")]
        cn4s = [r.get("CN_4A", 0) for r in successful]
        cn6s = [r.get("CN_6A", 0) for r in successful]

        print(f"\n{'=' * 60}")
        print("Summary Statistics")
        print(f"{'=' * 60}")
        print(f"\nProtein fold quality (pLDDT):")
        print(f"  Mean: {sum(plddts)/len(plddts):.1f}%")
        print(f"  Min:  {min(plddts):.1f}%")
        print(f"  Max:  {max(plddts):.1f}%")

        print(f"\nCoordination potential:")
        print(f"  Mean CN (4A): {sum(cn4s)/len(cn4s):.1f}")
        print(f"  Max CN (4A):  {max(cn4s)}")
        print(f"  Designs with CN >= 3: {sum(1 for c in cn4s if c >= 3)}/{len(cn4s)}")

    # Top designs
    print(f"\n{'=' * 60}")
    print("Top Designs (by CN)")
    print(f"{'=' * 60}")

    sorted_results = sorted(successful, key=lambda x: x.get("CN_4A", 0), reverse=True)
    for r in sorted_results[:5]:
        print(f"\n{r['name']}:")
        print(f"  pLDDT: {r.get('mean_plddt', 0):.1f}%")
        print(f"  Ligand pLDDT: {r.get('ligand_plddt', 'N/A')}")
        print(f"  Coordination (4A): {r.get('CN_4A', 0)}")
        if r.get("coord_residues"):
            coords = ", ".join(f"{c['res']}({c['dist']:.1f}A)" for c in r["coord_residues"][:5])
            print(f"  Near ligand: {coords}")

    # Save results
    summary_path = RF3_DIR / "rf3_analysis.json"
    with open(summary_path, "w") as f:
        json.dump({"predictions": all_results}, f, indent=2)
    print(f"\n\nSaved: {summary_path}")


if __name__ == "__main__":
    main()
