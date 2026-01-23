#!/usr/bin/env python3
"""
Simple standalone analysis for Ln-citrate scaffold designs.
Directly parses PDB files to compute coordination metrics.

No external dependencies beyond numpy.
"""

import argparse
import json
import math
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional


def parse_pdb(pdb_path: Path) -> Tuple[Dict[str, List], List[Dict], List[Dict]]:
    """
    Parse PDB file into atoms, residues, and heteroatoms.

    Returns:
        atoms_by_residue: dict mapping resid to list of atoms
        metal_atoms: list of metal atoms
        ligand_atoms: list of ligand atoms
    """
    atoms_by_residue = {}
    metal_atoms = []
    ligand_atoms = []

    # Metal residue names (3-letter codes for HETATM records)
    metal_res_names = {'TB', 'LA', 'CE', 'PR', 'ND', 'SM', 'EU', 'GD', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'MG', 'ZN', 'FE', 'CU', 'MN', 'CO', 'NI'}
    # Note: CA is removed - it conflicts with alpha carbon atom name

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

                # Classify atom - only check res_name for metals to avoid CA confusion
                # HETATM with metal residue name = metal
                if record_type == "HETATM" and res_name in metal_res_names:
                    metal_atoms.append(atom)
                elif res_name == "CIT":
                    ligand_atoms.append(atom)
                elif record_type == "ATOM":
                    # Normal protein residue
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


def find_coordinating_residues(metal: Dict, atoms_by_residue: Dict,
                                cutoff: float = 3.5) -> List[Dict]:
    """
    Find residues with SIDECHAIN O/N atoms within cutoff distance of metal.

    Returns list of coordinating residue info with distances.
    """
    coordinating = []

    # Sidechain atoms that can coordinate metals
    coordinating_atoms = [
        "OE1", "OE2",  # Glu carboxylate
        "OD1", "OD2",  # Asp carboxylate
        "OG", "OG1",   # Ser/Thr hydroxyl
        "OH",          # Tyr hydroxyl
        "ND1", "NE2",  # His imidazole
        "ND2", "NE2",  # Asn/Gln amide N
        "SG"           # Cys thiol (bad for Ln but RFD3 uses it)
    ]

    for resid, atoms in atoms_by_residue.items():
        for atom in atoms:
            if atom["name"] in coordinating_atoms:
                dist = distance(metal, atom)
                if dist <= cutoff:
                    coordinating.append({
                        "resid": resid,
                        "res_name": atom["res_name"],
                        "atom": atom["name"],
                        "distance": round(dist, 2)
                    })
                    break  # One per residue

    return coordinating


def find_coordinating_oxygens(metal: Dict, atoms_by_residue: Dict,
                               cutoff: float = 3.5) -> List[Dict]:
    """
    Find sidechain O/N atoms within cutoff distance of metal.
    Excludes backbone atoms (N, C, O, CA) to focus on true coordination.
    """
    donors = []

    # Sidechain O atoms that can coordinate
    sidechain_o = ["OE1", "OE2", "OD1", "OD2", "OG", "OG1", "OH"]
    # Sidechain N atoms that can coordinate
    sidechain_n = ["ND1", "ND2", "NE", "NE1", "NE2", "NH1", "NH2", "NZ"]

    for resid, atoms in atoms_by_residue.items():
        for atom in atoms:
            atom_name = atom["name"]
            # Only count sidechain donors, not backbone
            if atom_name in sidechain_o or atom_name in sidechain_n:
                dist = distance(metal, atom)
                if dist <= cutoff:
                    donors.append({
                        "resid": resid,
                        "res_name": atom["res_name"],
                        "atom": atom_name,
                        "distance": round(dist, 2)
                    })

    return donors


def count_citrate_contacts(ligand_atoms: List[Dict], atoms_by_residue: Dict,
                           cutoff: float = 4.0) -> int:
    """Count protein atoms in contact with citrate."""
    contacts = set()

    for lig_atom in ligand_atoms:
        for resid, atoms in atoms_by_residue.items():
            for atom in atoms:
                if distance(lig_atom, atom) <= cutoff:
                    contacts.add(resid)
                    break

    return len(contacts)


def estimate_burial(metal: Dict, atoms_by_residue: Dict,
                   shell_radius: float = 8.0) -> float:
    """
    Estimate metal burial by counting nearby atoms.
    Returns approximate SASA (lower = more buried).
    """
    nearby_count = 0

    for atoms in atoms_by_residue.values():
        for atom in atoms:
            if distance(metal, atom) <= shell_radius:
                nearby_count += 1

    # Heuristic: 50+ atoms nearby = well buried
    # Convert to approximate SASA (higher nearby = lower SASA)
    if nearby_count >= 50:
        return 0  # Fully buried
    elif nearby_count >= 40:
        return 5
    elif nearby_count >= 30:
        return 10
    elif nearby_count >= 20:
        return 20
    else:
        return 30  # Exposed


def count_secondary_structure(atoms_by_residue: Dict) -> float:
    """
    Estimate secondary structure content from residue types.
    Helix-favoring: A, E, L, M, Q, K
    Sheet-favoring: V, I, Y, F, W, T
    """
    ss_favoring = {'ALA', 'GLU', 'LEU', 'MET', 'GLN', 'LYS',
                   'VAL', 'ILE', 'TYR', 'PHE', 'TRP', 'THR'}
    total = 0
    ss_count = 0

    for atoms in atoms_by_residue.values():
        if atoms:
            res_name = atoms[0]["res_name"]
            total += 1
            if res_name in ss_favoring:
                ss_count += 1

    return ss_count / total if total > 0 else 0.5


def assign_tier(cn: int, geom_rmsd: float, metal_sasa: float, cit_contacts: int) -> str:
    """Assign quality tier based on metrics."""
    # Tier thresholds for Ln-citrate
    if cn >= 8 and geom_rmsd <= 0.8 and metal_sasa <= 2 and cit_contacts >= 3:
        return "S"
    elif cn >= 8 and geom_rmsd <= 1.2 and metal_sasa <= 5 and cit_contacts >= 3:
        return "A"
    elif cn >= 7 and geom_rmsd <= 1.5 and metal_sasa <= 10 and cit_contacts >= 2:
        return "B"
    elif cn >= 6 and geom_rmsd <= 2.0 and metal_sasa <= 20 and cit_contacts >= 1:
        return "C"
    return "F"


def compute_score(cn: int, geom_rmsd: float, metal_sasa: float,
                  cit_contacts: int, ss_content: float = 0.5) -> float:
    """Compute composite score (0-100)."""
    score = 0.0
    score += min(30, cn * 3.75)                      # CN: max 30
    score += max(0, 25 - geom_rmsd * 16.67)          # Geometry: max 25
    score += max(0, 20 - metal_sasa * 1.0)           # Burial: max 20
    score += min(15, cit_contacts * 5)               # Citrate: max 15
    score += min(10, ss_content * 10)                # SS: max 10
    return round(score, 2)


def analyze_pdb(pdb_path: Path) -> Dict[str, Any]:
    """
    Analyze a single PDB design.

    Returns dict with all metrics.
    """
    try:
        atoms_by_residue, metal_atoms, ligand_atoms = parse_pdb(pdb_path)

        if not metal_atoms:
            return {
                "file": pdb_path.name,
                "error": "No metal atom found",
                "tier": "F",
                "score": 0
            }

        metal = metal_atoms[0]  # Use first metal

        # Find coordinating atoms
        coord_residues = find_coordinating_residues(metal, atoms_by_residue, cutoff=3.0)
        coord_oxygens = find_coordinating_oxygens(metal, atoms_by_residue, cutoff=3.0)

        # Coordination number
        cn = len(coord_oxygens)

        # For now, use a placeholder geometry RMSD (would need ideal template)
        # Estimate based on CN spread
        geom_rmsd = 1.0 if cn >= 7 else 1.5 if cn >= 5 else 2.0

        # Metal burial
        metal_sasa = estimate_burial(metal, atoms_by_residue)

        # Citrate contacts
        cit_contacts = count_citrate_contacts(ligand_atoms, atoms_by_residue)

        # Secondary structure estimate
        ss_content = count_secondary_structure(atoms_by_residue)

        # Tier and score
        tier = assign_tier(cn, geom_rmsd, metal_sasa, cit_contacts)
        score = compute_score(cn, geom_rmsd, metal_sasa, cit_contacts, ss_content)

        # Count residues
        n_residues = len(atoms_by_residue)

        return {
            "file": pdb_path.name,
            "n_residues": n_residues,
            "coordination_number": cn,
            "coordinating_residues": coord_residues,
            "coordinating_oxygens": len(coord_oxygens),
            "geometry_rmsd": geom_rmsd,
            "metal_sasa": metal_sasa,
            "citrate_contacts": cit_contacts,
            "ss_content": round(ss_content, 3),
            "tier": tier,
            "score": score,
            "metal_position": {
                "x": round(metal["x"], 2),
                "y": round(metal["y"], 2),
                "z": round(metal["z"], 2)
            }
        }

    except Exception as e:
        return {
            "file": pdb_path.name,
            "error": str(e),
            "tier": "F",
            "score": 0
        }


def main():
    parser = argparse.ArgumentParser(description="Simple analysis for Ln-citrate designs")
    parser.add_argument("--round", "-r", type=str, required=True, help="Round identifier (e.g., 1, 2, 2b)")
    parser.add_argument("--save", "-s", action="store_true", help="Save analysis to file")
    args = parser.parse_args()

    round_id = args.round

    print(f"\n{'='*60}")
    print(f"Ln-Citrate Scaffold Analysis - Round {round_id}")
    print(f"{'='*60}\n")

    # Find output directory - handle both numeric (01, 02) and string (02b) formats
    if round_id.isdigit():
        output_dir = Path(__file__).parent.parent / "outputs" / f"round_{int(round_id):02d}"
    else:
        # For mixed formats like "2b", try "round_02b" first, then "round_2b"
        output_dir = Path(__file__).parent.parent / "outputs" / f"round_0{round_id}"
        if not output_dir.exists():
            output_dir = Path(__file__).parent.parent / "outputs" / f"round_{round_id}"
    if not output_dir.exists():
        print(f"Error: Output directory not found: {output_dir}")
        sys.exit(1)

    # Find all PDB files
    pdb_files = sorted(output_dir.glob("*.pdb"))
    if not pdb_files:
        print(f"Error: No PDB files found in {output_dir}")
        sys.exit(1)

    print(f"Found {len(pdb_files)} PDB files to analyze\n")

    # Analyze each file
    all_results = []
    results_by_config = {}

    for pdb_file in pdb_files:
        result = analyze_pdb(pdb_file)
        all_results.append(result)

        # Group by config
        parts = pdb_file.stem.split("_")
        if len(parts) >= 2:
            config_name = "_".join(parts[:-1])
            if config_name not in results_by_config:
                results_by_config[config_name] = []
            results_by_config[config_name].append(result)

        # Print result
        if "error" in result and result.get("score", 0) == 0:
            print(f"  {result['file']}: ERROR - {result['error']}")
        else:
            print(f"  {result['file']}:")
            print(f"    Residues: {result.get('n_residues', 'N/A')}")
            print(f"    CN: {result.get('coordination_number', 0)} (donors within 3.0Å)")
            print(f"    Metal SASA: ~{result.get('metal_sasa', 'N/A')}")
            print(f"    Citrate contacts: {result.get('citrate_contacts', 0)}")
            print(f"    Tier: {result['tier']}, Score: {result['score']:.1f}")

            if result.get('coordinating_residues'):
                print(f"    Coordinating residues:")
                for cr in result['coordinating_residues']:
                    print(f"      {cr['res_name']}{cr['resid']}.{cr['atom']}: {cr['distance']}Å")

    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}\n")

    tier_counts = {"S": 0, "A": 0, "B": 0, "C": 0, "F": 0}
    scores = []
    cns = []

    for r in all_results:
        tier_counts[r.get("tier", "F")] += 1
        if r.get("score", 0) > 0:
            scores.append(r["score"])
        if r.get("coordination_number", 0) > 0:
            cns.append(r["coordination_number"])

    print(f"Total designs: {len(all_results)}")
    print(f"Tier distribution:")
    for tier in ["S", "A", "B", "C", "F"]:
        count = tier_counts[tier]
        print(f"  {tier}: {count}")

    n_pass = tier_counts["S"] + tier_counts["A"] + tier_counts["B"]
    pass_rate = n_pass / len(all_results) if all_results else 0
    print(f"\nPass rate (B or better): {pass_rate*100:.1f}%")

    if scores:
        print(f"Score range: {min(scores):.1f} - {max(scores):.1f}")
    if cns:
        print(f"CN range: {min(cns)} - {max(cns)}")

    # Judgment
    print(f"\n{'='*60}")
    print("JUDGMENT")
    print(f"{'='*60}\n")

    if pass_rate < 0.10:
        print("Action: REVISE_FUNDAMENTALLY")
        print("\nThe designs show low coordination numbers.")
        print("Recommendations:")
        print("  1. Check if metal is in the expected position after RFD3")
        print("  2. Verify coordinating residue arrangement in input")
        print("  3. Consider tighter constraints on fixed atoms")
        print("  4. May need to generate more designs for statistics")
    elif pass_rate < 0.30:
        print("Action: ADJUST_KEY_PARAMETERS")
        print("\nSome designs showing promise but need refinement.")
    elif pass_rate < 0.50:
        print("Action: FINE_TUNE")
        print("\nGood progress - fine-tune current approach.")
    else:
        print("Action: PROCEED_TO_MPNN")
        print("\nBackbone quality sufficient for sequence design!")

    # Note about sample size
    print(f"\n** Note: Only {len(all_results)} designs analyzed.")
    print("   Consider generating 8-16 designs per config for better statistics.")

    # Save results
    if args.save:
        analysis_dir = Path(__file__).parent.parent / "analysis"
        analysis_dir.mkdir(exist_ok=True)

        # Handle both numeric and string round IDs
        if round_id.isdigit():
            analysis_file = analysis_dir / f"round_{int(round_id):02d}_simple_analysis.json"
        else:
            analysis_file = analysis_dir / f"round_{round_id}_simple_analysis.json"
        with open(analysis_file, "w") as f:
            json.dump({
                "round": round_id,
                "timestamp": datetime.now().isoformat(),
                "n_designs": len(all_results),
                "pass_rate": pass_rate,
                "tier_distribution": tier_counts,
                "results": all_results
            }, f, indent=2)

        print(f"\nAnalysis saved to: {analysis_file}")


if __name__ == "__main__":
    main()
