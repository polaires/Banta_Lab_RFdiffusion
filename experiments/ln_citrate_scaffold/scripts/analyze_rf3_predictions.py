#!/usr/bin/env python3
"""
Analyze RF3 predicted structures for metal coordination potential.

Goals:
1. Extract protein fold quality (pLDDT)
2. Find citrate position
3. Identify D/E residues near citrate
4. Evaluate potential metal coordination sites
"""

import json
import math
from pathlib import Path
from collections import defaultdict

EXPERIMENT_DIR = Path(__file__).parent.parent
RF3_DIR = EXPERIMENT_DIR / "outputs" / "round_04_rf3"


def parse_mmcif_atoms(filepath: Path) -> tuple:
    """
    Parse mmCIF format file to extract atom coordinates and pLDDT.

    Returns:
        protein_atoms: list of protein atom dicts
        ligand_atoms: list of ligand atom dicts
        metal_atoms: list of metal atom dicts
    """
    protein_atoms = []
    ligand_atoms = []
    metal_atoms = []

    in_atom_site = False
    header_fields = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()

            # Detect start of atom_site loop
            if line.startswith("_atom_site."):
                in_atom_site = True
                field = line.split(".")[1].strip()
                header_fields.append(field)
                continue

            # End of loop - when we hit # or another section
            if in_atom_site and (line.startswith("#") or line.startswith("_") or line.startswith("loop_")):
                in_atom_site = False
                continue

            # Parse atom records
            if in_atom_site and (line.startswith("ATOM") or line.startswith("HETATM")):
                parts = line.split()
                if len(parts) < len(header_fields):
                    continue

                atom = {}
                for i, field in enumerate(header_fields):
                    if i < len(parts):
                        atom[field] = parts[i]

                # Extract coordinates
                try:
                    atom["x"] = float(atom.get("Cartn_x", 0))
                    atom["y"] = float(atom.get("Cartn_y", 0))
                    atom["z"] = float(atom.get("Cartn_z", 0))
                    atom["plddt"] = float(atom.get("B_iso_or_equiv", 0))
                except ValueError:
                    continue

                # Categorize atoms
                group = atom.get("group_PDB", "")
                comp_id = atom.get("label_comp_id", "")

                # Check for metals
                metal_names = ["TB", "GD", "EU", "LA", "CE", "MG", "CA", "ZN", "FE", "CU", "MN"]
                if comp_id.upper() in metal_names or atom.get("type_symbol", "").upper() in metal_names:
                    metal_atoms.append(atom)
                elif group == "HETATM":
                    ligand_atoms.append(atom)
                else:
                    protein_atoms.append(atom)

    return protein_atoms, ligand_atoms, metal_atoms


def distance(atom1: dict, atom2: dict) -> float:
    """Calculate 3D distance between two atoms."""
    return math.sqrt(
        (atom1["x"] - atom2["x"])**2 +
        (atom1["y"] - atom2["y"])**2 +
        (atom1["z"] - atom2["z"])**2
    )


def get_residue_atoms(atoms: list, res_type: str) -> dict:
    """Group atoms by residue for given residue type (D, E, etc.)."""
    residues = defaultdict(list)
    for atom in atoms:
        comp = atom.get("label_comp_id", "")
        if comp == res_type:
            seq_id = atom.get("label_seq_id", "")
            residues[seq_id].append(atom)
    return dict(residues)


def find_coordinating_residues(protein_atoms: list, ligand_atoms: list,
                                cutoff: float = 5.0) -> list:
    """
    Find D/E residues with carboxylate oxygens near the ligand.

    Args:
        protein_atoms: list of protein atom dicts
        ligand_atoms: list of ligand atom dicts
        cutoff: distance cutoff in Angstroms

    Returns:
        list of (residue_type, residue_id, min_distance, coordinating_atoms)
    """
    # Get ligand center and oxygen positions
    ligand_oxygens = [a for a in ligand_atoms if a.get("type_symbol", "").upper() == "O"]
    if not ligand_oxygens:
        # No oxygens found by symbol, try atom name
        ligand_oxygens = [a for a in ligand_atoms if a.get("label_atom_id", "").upper().startswith("O")]

    if not ligand_oxygens:
        return []

    # Calculate ligand center
    lig_center = {
        "x": sum(a["x"] for a in ligand_atoms) / len(ligand_atoms),
        "y": sum(a["y"] for a in ligand_atoms) / len(ligand_atoms),
        "z": sum(a["z"] for a in ligand_atoms) / len(ligand_atoms)
    }

    coordinating = []

    # Check ASP (D) - carboxylate is OD1, OD2
    asp_residues = get_residue_atoms(protein_atoms, "ASP")
    for res_id, atoms in asp_residues.items():
        carboxyl_o = [a for a in atoms if a.get("label_atom_id", "") in ["OD1", "OD2"]]
        if carboxyl_o:
            min_dist = min(distance(o, lig_center) for o in carboxyl_o)
            if min_dist <= cutoff:
                coordinating.append(("ASP", res_id, round(min_dist, 2), [a.get("label_atom_id") for a in carboxyl_o]))

    # Check GLU (E) - carboxylate is OE1, OE2
    glu_residues = get_residue_atoms(protein_atoms, "GLU")
    for res_id, atoms in glu_residues.items():
        carboxyl_o = [a for a in atoms if a.get("label_atom_id", "") in ["OE1", "OE2"]]
        if carboxyl_o:
            min_dist = min(distance(o, lig_center) for o in carboxyl_o)
            if min_dist <= cutoff:
                coordinating.append(("GLU", res_id, round(min_dist, 2), [a.get("label_atom_id") for a in carboxyl_o]))

    # Check ASN (N) - amide oxygen is OD1
    asn_residues = get_residue_atoms(protein_atoms, "ASN")
    for res_id, atoms in asn_residues.items():
        amide_o = [a for a in atoms if a.get("label_atom_id", "") == "OD1"]
        if amide_o:
            min_dist = min(distance(o, lig_center) for o in amide_o)
            if min_dist <= cutoff:
                coordinating.append(("ASN", res_id, round(min_dist, 2), ["OD1"]))

    # Check GLN (Q) - amide oxygen is OE1
    gln_residues = get_residue_atoms(protein_atoms, "GLN")
    for res_id, atoms in gln_residues.items():
        amide_o = [a for a in atoms if a.get("label_atom_id", "") == "OE1"]
        if amide_o:
            min_dist = min(distance(o, lig_center) for o in amide_o)
            if min_dist <= cutoff:
                coordinating.append(("GLN", res_id, round(min_dist, 2), ["OE1"]))

    # Check HIS (H) - imidazole nitrogens
    his_residues = get_residue_atoms(protein_atoms, "HIS")
    for res_id, atoms in his_residues.items():
        imid_n = [a for a in atoms if a.get("label_atom_id", "") in ["ND1", "NE2"]]
        if imid_n:
            min_dist = min(distance(n, lig_center) for n in imid_n)
            if min_dist <= cutoff:
                coordinating.append(("HIS", res_id, round(min_dist, 2), [a.get("label_atom_id") for a in imid_n]))

    return sorted(coordinating, key=lambda x: x[2])


def analyze_structure(filepath: Path) -> dict:
    """Analyze a single RF3 prediction."""
    result = {"name": filepath.stem, "file": filepath.name}

    protein_atoms, ligand_atoms, metal_atoms = parse_mmcif_atoms(filepath)

    if not protein_atoms:
        result["error"] = "No protein atoms found"
        return result

    # Protein metrics
    plddt_values = [a["plddt"] for a in protein_atoms if a["plddt"] > 0]
    if plddt_values:
        result["mean_plddt"] = round(sum(plddt_values) / len(plddt_values) * 100, 1)
        result["min_plddt"] = round(min(plddt_values) * 100, 1)
        result["max_plddt"] = round(max(plddt_values) * 100, 1)

    # Count residues
    residue_ids = set()
    for atom in protein_atoms:
        seq_id = atom.get("label_seq_id", "")
        if seq_id:
            residue_ids.add(seq_id)
    result["n_residues"] = len(residue_ids)

    # Ligand metrics
    result["has_ligand"] = len(ligand_atoms) > 0
    if ligand_atoms:
        lig_plddt = [a["plddt"] for a in ligand_atoms if a["plddt"] > 0]
        if lig_plddt:
            result["ligand_plddt"] = round(sum(lig_plddt) / len(lig_plddt) * 100, 1)
        result["n_ligand_atoms"] = len(ligand_atoms)

    # Metal metrics
    result["has_metal"] = len(metal_atoms) > 0
    if metal_atoms:
        result["metals"] = list(set(a.get("label_comp_id", "") for a in metal_atoms))

    # Coordination analysis
    if ligand_atoms:
        coordinating = find_coordinating_residues(protein_atoms, ligand_atoms, cutoff=6.0)
        result["nearby_residues"] = len(coordinating)
        result["coordinating_residues"] = [
            {"type": r[0], "id": r[1], "dist": r[2]}
            for r in coordinating[:10]  # Top 10 closest
        ]

        # Count by type
        type_counts = defaultdict(int)
        for res in coordinating:
            type_counts[res[0]] += 1
        result["coord_counts"] = dict(type_counts)

        # Potential coordination number (residues within 4.0 A of ligand center)
        close_residues = find_coordinating_residues(protein_atoms, ligand_atoms, cutoff=4.0)
        result["potential_CN"] = len(close_residues)

    return result


def main():
    print("=" * 60)
    print("RF3 Prediction Analysis - Metal Coordination Potential")
    print("=" * 60)

    # Find all RF3 predictions
    pdb_files = sorted(RF3_DIR.glob("*_rf3.pdb"))

    if not pdb_files:
        print(f"No RF3 predictions found in {RF3_DIR}")
        return

    print(f"\nFound {len(pdb_files)} predictions\n")

    all_results = []

    for pdb_file in pdb_files:
        print(f"\n{'='*40}")
        print(f"Analyzing: {pdb_file.name}")

        result = analyze_structure(pdb_file)
        all_results.append(result)

        if "error" in result:
            print(f"  ERROR: {result['error']}")
            continue

        print(f"  Residues: {result.get('n_residues', 'N/A')}")
        print(f"  Mean pLDDT: {result.get('mean_plddt', 'N/A')}%")
        print(f"  Ligand present: {result.get('has_ligand', False)}")
        if result.get("has_ligand"):
            print(f"  Ligand pLDDT: {result.get('ligand_plddt', 'N/A')}%")
        print(f"  Metal present: {result.get('has_metal', False)}")

        if result.get("coordinating_residues"):
            print(f"\n  Potential coordinating residues (within 6A of ligand):")
            for res in result["coordinating_residues"][:5]:
                print(f"    {res['type']}{res['id']}: {res['dist']} A")
            print(f"\n  Coordination potential:")
            print(f"    Residues within 4A: {result.get('potential_CN', 0)}")
            print(f"    Residues within 6A: {result.get('nearby_residues', 0)}")
            if result.get("coord_counts"):
                print(f"    By type: {result['coord_counts']}")

    # Save results
    summary = {
        "n_predictions": len(all_results),
        "predictions": all_results
    }

    summary_path = RF3_DIR / "rf3_analysis.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    # Summary table
    print(f"\n{'='*60}")
    print("Summary")
    print("="*60)
    print(f"\n{'Name':<40} {'pLDDT':>8} {'Lig':>6} {'CN4':>5} {'CN6':>5}")
    print("-"*70)

    for r in all_results:
        name = r["name"][:38]
        plddt = f"{r.get('mean_plddt', 0):.1f}%" if r.get('mean_plddt') else "N/A"
        lig = f"{r.get('ligand_plddt', 0):.0f}%" if r.get('ligand_plddt') else "No"
        cn4 = str(r.get('potential_CN', 0))
        cn6 = str(r.get('nearby_residues', 0))
        print(f"{name:<40} {plddt:>8} {lig:>6} {cn4:>5} {cn6:>5}")

    print(f"\nAnalysis saved to: {summary_path}")

    # Best candidates
    candidates = [r for r in all_results if r.get("potential_CN", 0) >= 3]
    if candidates:
        print(f"\n{'='*60}")
        print(f"Top candidates (CN >= 3 within 4A):")
        for c in sorted(candidates, key=lambda x: x.get("potential_CN", 0), reverse=True):
            print(f"  - {c['name']}: CN={c['potential_CN']}, pLDDT={c.get('mean_plddt', 0):.1f}%")
            if c.get("coordinating_residues"):
                types = [f"{r['type']}{r['id']}" for r in c["coordinating_residues"][:5]]
                print(f"    Near ligand: {', '.join(types)}")


if __name__ == "__main__":
    main()
