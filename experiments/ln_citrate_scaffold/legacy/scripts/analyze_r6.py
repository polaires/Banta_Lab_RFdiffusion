#!/usr/bin/env python3
"""
Round 6 Analysis Script - Comprehensive Coordination Analysis

Analyzes designs from:
- r6_a: Partial diffusion from R5 best (partial_t=0.5)
- r6_b: Aggressive partial diffusion (partial_t=0.7)
- r6_c: Homodimer with C2 symmetry
- r6_d: Extended de novo with stronger H-bond conditioning

Key metrics:
- Coordination number at 3.5A, 4.0A, 5.0A cutoffs
- Donor atom types (O vs N vs S)
- Per-chain contribution for dimers
- Best designs ranked by total CN
"""

import json
import math
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Any


def parse_pdb(pdb_path: Path) -> Tuple[Dict[str, List], List[Dict], List[Dict]]:
    """
    Parse PDB file into atoms, residues, and heteroatoms.
    """
    atoms_by_residue = {}
    metal_atoms = []
    ligand_atoms = []

    metal_res_names = {'TB', 'LA', 'CE', 'PR', 'ND', 'SM', 'EU', 'GD', 'DY',
                       'HO', 'ER', 'TM', 'YB', 'LU', 'MG', 'ZN', 'FE', 'CU', 'MN', 'CO', 'NI'}

    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                record_type = line[0:6].strip()
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21]
                res_num = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                atom = {
                    "name": atom_name,
                    "res_name": res_name,
                    "chain": chain,
                    "res_num": res_num,
                    "x": x, "y": y, "z": z
                }

                if record_type == "HETATM" and res_name in metal_res_names:
                    metal_atoms.append(atom)
                elif res_name == "CIT":
                    ligand_atoms.append(atom)
                elif record_type == "ATOM":
                    resid = f"{chain}{res_num}"
                    if resid not in atoms_by_residue:
                        atoms_by_residue[resid] = []
                    atoms_by_residue[resid].append(atom)

    return atoms_by_residue, metal_atoms, ligand_atoms


def distance(atom1: Dict, atom2: Dict) -> float:
    """Euclidean distance between two atoms."""
    return math.sqrt(
        (atom1["x"] - atom2["x"])**2 +
        (atom1["y"] - atom2["y"])**2 +
        (atom1["z"] - atom2["z"])**2
    )


def find_coordinating_atoms(metal: Dict, atoms_by_residue: Dict, cutoff: float = 3.5) -> List[Dict]:
    """
    Find sidechain O/N atoms within cutoff distance of metal.
    Returns detailed information about each coordinating atom.
    """
    donors = []

    # Sidechain atoms that can coordinate (by element)
    oxygen_atoms = ["OE1", "OE2", "OD1", "OD2", "OG", "OG1", "OH"]
    nitrogen_atoms = ["ND1", "ND2", "NE", "NE1", "NE2", "NH1", "NH2", "NZ"]
    sulfur_atoms = ["SG"]

    for resid, atoms in atoms_by_residue.items():
        for atom in atoms:
            atom_name = atom["name"]
            element = None

            if atom_name in oxygen_atoms:
                element = "O"
            elif atom_name in nitrogen_atoms:
                element = "N"
            elif atom_name in sulfur_atoms:
                element = "S"

            if element:
                dist = distance(metal, atom)
                if dist <= cutoff:
                    donors.append({
                        "resid": resid,
                        "res_name": atom["res_name"],
                        "atom": atom_name,
                        "element": element,
                        "distance": round(dist, 2),
                        "chain": atom["chain"]
                    })

    return donors


def calculate_coordination_metrics(metal: Dict, atoms_by_residue: Dict) -> Dict:
    """
    Calculate comprehensive coordination metrics at multiple cutoffs.
    """
    metrics = {}

    for cutoff in [3.0, 3.5, 4.0, 5.0, 6.0]:
        donors = find_coordinating_atoms(metal, atoms_by_residue, cutoff)

        # Count by element
        o_count = sum(1 for d in donors if d["element"] == "O")
        n_count = sum(1 for d in donors if d["element"] == "N")
        s_count = sum(1 for d in donors if d["element"] == "S")

        # Count by chain (for dimers)
        chains = set(d["chain"] for d in donors)
        by_chain = {c: sum(1 for d in donors if d["chain"] == c) for c in chains}

        # Carboxylate-specific (Glu/Asp O atoms)
        carboxylate_count = sum(1 for d in donors
                                if d["res_name"] in ["GLU", "ASP"]
                                and d["element"] == "O")

        metrics[f"cn_{cutoff}A"] = len(donors)
        metrics[f"o_{cutoff}A"] = o_count
        metrics[f"n_{cutoff}A"] = n_count
        metrics[f"s_{cutoff}A"] = s_count
        metrics[f"carboxylate_{cutoff}A"] = carboxylate_count

        if len(by_chain) > 1:
            metrics[f"by_chain_{cutoff}A"] = by_chain

    return metrics


def analyze_design(pdb_path: Path) -> Dict:
    """
    Analyze a single design PDB.
    """
    atoms_by_residue, metal_atoms, ligand_atoms = parse_pdb(pdb_path)

    if not metal_atoms:
        return {"error": "No metal found", "file": pdb_path.name}

    metal = metal_atoms[0]

    # Basic info
    chains = set(a["chain"] for atoms in atoms_by_residue.values() for a in atoms)
    n_residues = len(atoms_by_residue)

    # Coordination metrics
    coord_metrics = calculate_coordination_metrics(metal, atoms_by_residue)

    # Get detailed donors at 4A
    donors_4a = find_coordinating_atoms(metal, atoms_by_residue, 4.0)

    # Classify quality
    cn_4a = coord_metrics.get("cn_4.0A", 0)
    carboxylate_4a = coord_metrics.get("carboxylate_4.0A", 0)

    if cn_4a >= 6 and carboxylate_4a >= 4:
        quality = "excellent"
    elif cn_4a >= 4 and carboxylate_4a >= 3:
        quality = "good"
    elif cn_4a >= 2:
        quality = "moderate"
    else:
        quality = "poor"

    # Score (higher is better)
    score = (
        cn_4a * 10 +
        carboxylate_4a * 5 +
        coord_metrics.get("cn_3.5A", 0) * 3 -
        coord_metrics.get("s_4.0A", 0) * 5  # Penalize sulfur for Ln
    )

    return {
        "file": pdb_path.name,
        "n_residues": n_residues,
        "n_chains": len(chains),
        "chains": list(chains),
        "metal": metal["res_name"],
        "metal_pos": [round(metal["x"], 2), round(metal["y"], 2), round(metal["z"], 2)],
        "quality": quality,
        "score": score,
        **coord_metrics,
        "donors_4a": donors_4a[:10] if donors_4a else []  # Top 10 closest
    }


def analyze_round6(output_dir: Path) -> Dict:
    """
    Analyze all Round 6 designs.
    """
    results = {
        "round": "6",
        "timestamp": datetime.now().isoformat(),
        "configs": {},
        "best_designs": [],
        "summary": {}
    }

    # Group by config prefix
    pdbs = sorted(output_dir.glob("*.pdb"))

    if not pdbs:
        print(f"No PDB files found in {output_dir}")
        return results

    print(f"Found {len(pdbs)} PDB files")

    # Analyze each design
    all_analyses = []
    for pdb_path in pdbs:
        analysis = analyze_design(pdb_path)
        all_analyses.append(analysis)

        # Group by config
        config_name = "_".join(pdb_path.stem.split("_")[:-1])
        if config_name not in results["configs"]:
            results["configs"][config_name] = []
        results["configs"][config_name].append(analysis)

    # Summary stats
    print(f"\n{'='*70}")
    print("ROUND 6 ANALYSIS SUMMARY")
    print(f"{'='*70}")

    for config_name, analyses in results["configs"].items():
        cn_values = [a.get("cn_4.0A", 0) for a in analyses if "error" not in a]
        carb_values = [a.get("carboxylate_4.0A", 0) for a in analyses if "error" not in a]

        if cn_values:
            print(f"\n{config_name}:")
            print(f"  Designs: {len(cn_values)}")
            print(f"  CN@4A:   min={min(cn_values)}, max={max(cn_values)}, mean={sum(cn_values)/len(cn_values):.1f}")
            print(f"  Carbox:  min={min(carb_values)}, max={max(carb_values)}, mean={sum(carb_values)/len(carb_values):.1f}")

            results["summary"][config_name] = {
                "n_designs": len(cn_values),
                "cn_4a_range": [min(cn_values), max(cn_values)],
                "cn_4a_mean": round(sum(cn_values)/len(cn_values), 2),
                "carboxylate_mean": round(sum(carb_values)/len(carb_values), 2),
            }

    # Top designs overall
    successful = [a for a in all_analyses if "error" not in a]
    sorted_designs = sorted(successful, key=lambda x: x.get("score", 0), reverse=True)

    print(f"\n{'='*70}")
    print("TOP 10 DESIGNS BY SCORE")
    print(f"{'='*70}")

    for i, d in enumerate(sorted_designs[:10]):
        print(f"\n{i+1}. {d['file']}")
        print(f"   Score: {d['score']}, Quality: {d['quality']}")
        print(f"   CN@4A: {d.get('cn_4.0A', 0)} (Carboxylate: {d.get('carboxylate_4.0A', 0)})")
        print(f"   CN@3.5A: {d.get('cn_3.5A', 0)}")
        if d.get("donors_4a"):
            donors_str = ", ".join(f"{d['res_name']}{d['resid'][1:]}({d['distance']}A)"
                                  for d in d["donors_4a"][:5])
            print(f"   Donors: {donors_str}")

    results["best_designs"] = sorted_designs[:20]

    return results


def main():
    experiment_dir = Path(__file__).parent.parent
    output_dir = experiment_dir / "outputs" / "round_06"
    analysis_dir = experiment_dir / "analysis"

    analysis_dir.mkdir(exist_ok=True)

    print(f"Analyzing Round 6 outputs from: {output_dir}")

    if not output_dir.exists():
        print(f"Output directory does not exist: {output_dir}")
        return

    results = analyze_round6(output_dir)

    # Save results
    output_path = analysis_dir / "round_06_analysis.json"
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)

    print(f"\n\nAnalysis saved to: {output_path}")

    # Generate markdown summary
    md_path = analysis_dir / "ROUND_6_SUMMARY.md"
    with open(md_path, "w") as f:
        f.write("# Round 6 Analysis Summary\n\n")
        f.write(f"**Date:** {results['timestamp'][:10]}\n\n")

        f.write("## Configuration Results\n\n")
        for config, stats in results.get("summary", {}).items():
            f.write(f"### {config}\n")
            f.write(f"- Designs: {stats.get('n_designs', 0)}\n")
            f.write(f"- CN@4A range: {stats.get('cn_4a_range', [0,0])}\n")
            f.write(f"- CN@4A mean: {stats.get('cn_4a_mean', 0)}\n")
            f.write(f"- Carboxylate mean: {stats.get('carboxylate_mean', 0)}\n\n")

        f.write("## Top Designs\n\n")
        f.write("| Rank | File | Score | CN@4A | Carboxylate | Quality |\n")
        f.write("|------|------|-------|-------|-------------|----------|\n")
        for i, d in enumerate(results.get("best_designs", [])[:10]):
            f.write(f"| {i+1} | {d['file']} | {d.get('score', 0)} | {d.get('cn_4.0A', 0)} | {d.get('carboxylate_4.0A', 0)} | {d.get('quality', 'N/A')} |\n")

    print(f"Markdown summary saved to: {md_path}")


if __name__ == "__main__":
    main()
