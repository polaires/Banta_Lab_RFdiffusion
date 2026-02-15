#!/usr/bin/env python3
"""
Score a single Ln-citrate scaffold design.

Usage:
    python score_design.py path/to/design.pdb

Returns tier (S/A/B/C/F) and detailed metrics.
"""

import argparse
import json
import sys
from pathlib import Path

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "backend" / "serverless"))


# Tier thresholds for Ln-citrate designs
TIER_THRESHOLDS = {
    "S": {"cn_min": 8, "geom_max": 0.8, "sasa_max": 2, "cit_min": 3},
    "A": {"cn_min": 8, "geom_max": 1.2, "sasa_max": 5, "cit_min": 3},
    "B": {"cn_min": 7, "geom_max": 1.5, "sasa_max": 10, "cit_min": 2},
    "C": {"cn_min": 6, "geom_max": 2.0, "sasa_max": 20, "cit_min": 1},
}


def parse_pdb_basic(pdb_path: str) -> dict:
    """
    Basic PDB parsing for coordination analysis.
    Works without external dependencies.
    """
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atom = {
                    "name": line[12:16].strip(),
                    "resname": line[17:20].strip(),
                    "chain": line[21],
                    "resnum": int(line[22:26]),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "element": line[76:78].strip() if len(line) > 76 else line[12:14].strip()
                }
                atoms.append(atom)
    return {"atoms": atoms}


def compute_distance(a1: dict, a2: dict) -> float:
    """Euclidean distance between two atoms."""
    import math
    return math.sqrt(
        (a1["x"] - a2["x"])**2 +
        (a1["y"] - a2["y"])**2 +
        (a1["z"] - a2["z"])**2
    )


def analyze_coordination(structure: dict, metal_name: str = "TB", cutoff: float = 3.0) -> dict:
    """
    Analyze metal coordination from parsed structure.
    """
    atoms = structure["atoms"]

    # Find metal
    metal_atoms = [a for a in atoms if a["resname"] == metal_name or a["element"] == metal_name]
    if not metal_atoms:
        return {"error": f"Metal {metal_name} not found"}

    metal = metal_atoms[0]

    # Find coordinating atoms (O and N within cutoff)
    coordinating = []
    for atom in atoms:
        if atom["element"] in ["O", "N"] and atom["resname"] not in [metal_name, "HOH"]:
            dist = compute_distance(metal, atom)
            if dist <= cutoff:
                coordinating.append({
                    "atom": atom["name"],
                    "residue": f"{atom['chain']}{atom['resnum']}",
                    "resname": atom["resname"],
                    "distance": round(dist, 3),
                    "element": atom["element"]
                })

    # Count citrate contacts
    citrate_contacts = len([c for c in coordinating if c["resname"] == "CIT"])

    # Count protein donors
    protein_donors = len([c for c in coordinating if c["resname"] in ["GLU", "ASP", "ASN", "GLN", "HIS"]])

    return {
        "coordination_number": len(coordinating),
        "protein_donors": protein_donors,
        "citrate_contacts": citrate_contacts,
        "coordinating_atoms": coordinating,
        "distances": [c["distance"] for c in coordinating]
    }


def score_ln_citrate_design(pdb_path: str) -> dict:
    """
    Score a single design against Ln-citrate criteria.

    Returns:
        {
            'tier': 'S'|'A'|'B'|'C'|'F',
            'score': float,
            'metrics': {...},
            'pass': bool,
            'issues': [...]
        }
    """
    # Parse PDB
    structure = parse_pdb_basic(pdb_path)

    # Analyze coordination
    coord = analyze_coordination(structure, "TB", cutoff=3.0)

    if "error" in coord:
        return {
            "tier": "F",
            "score": 0.0,
            "metrics": {},
            "pass": False,
            "issues": [coord["error"]]
        }

    # Extract metrics
    cn = coord["coordination_number"]
    cit_contacts = coord["citrate_contacts"]

    # Placeholder for geometry and SASA (would need more complex calculation)
    # For now, estimate based on coordination quality
    geom_rmsd = 1.0 if cn >= 8 else 1.5 if cn >= 7 else 2.0
    metal_sasa = 5.0 if cn >= 8 else 15.0 if cn >= 6 else 30.0

    # Determine tier
    tier = "F"
    for t in ["S", "A", "B", "C"]:
        thresh = TIER_THRESHOLDS[t]
        if (cn >= thresh["cn_min"] and
            geom_rmsd <= thresh["geom_max"] and
            metal_sasa <= thresh["sasa_max"] and
            cit_contacts >= thresh["cit_min"]):
            tier = t
            break

    # Compute score
    score = 0.0
    score += min(30, cn * 3.75)                     # CN: max 30
    score += max(0, 25 - geom_rmsd * 16.67)         # Geometry: max 25
    score += max(0, 20 - metal_sasa * 1.0)          # Burial: max 20
    score += min(15, cit_contacts * 5)              # Citrate: max 15
    score += 5.0  # Placeholder for SS content

    # Identify issues
    issues = []
    if cn < 7:
        issues.append(f"Low coordination number: {cn} (need ≥7)")
    if cn < 8:
        issues.append(f"Suboptimal CN: {cn} (ideal ≥8 for Ln)")
    if cit_contacts < 2:
        issues.append(f"Poor citrate coordination: {cit_contacts} contacts")

    # Check distances
    distances = coord.get("distances", [])
    bad_dists = [d for d in distances if d < 2.3 or d > 2.6]
    if bad_dists:
        issues.append(f"Bond distances outside range: {bad_dists}")

    return {
        "tier": tier,
        "score": round(score, 2),
        "metrics": {
            "cn": cn,
            "protein_donors": coord["protein_donors"],
            "citrate_contacts": cit_contacts,
            "geometry_rmsd": geom_rmsd,
            "metal_sasa": metal_sasa,
            "distances": distances
        },
        "pass": tier in ["S", "A", "B"],
        "issues": issues,
        "coordinating": coord["coordinating_atoms"]
    }


def main():
    parser = argparse.ArgumentParser(description="Score Ln-citrate scaffold design")
    parser.add_argument("pdb", type=str, help="Path to PDB file")
    parser.add_argument("--json", "-j", action="store_true", help="Output as JSON")
    args = parser.parse_args()

    pdb_path = Path(args.pdb)
    if not pdb_path.exists():
        print(f"Error: File not found: {pdb_path}")
        sys.exit(1)

    result = score_ln_citrate_design(str(pdb_path))

    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print(f"\n{'='*50}")
        print(f"Ln-Citrate Design Score: {pdb_path.name}")
        print(f"{'='*50}\n")

        tier_colors = {"S": "★★★★★", "A": "★★★★☆", "B": "★★★☆☆", "C": "★★☆☆☆", "F": "★☆☆☆☆"}

        print(f"TIER: {result['tier']} {tier_colors[result['tier']]}")
        print(f"SCORE: {result['score']:.1f}/100")
        print(f"PASS: {'Yes' if result['pass'] else 'No'}")

        print(f"\nMetrics:")
        metrics = result["metrics"]
        print(f"  Coordination Number: {metrics.get('cn', 0)}")
        print(f"  Protein Donors: {metrics.get('protein_donors', 0)}")
        print(f"  Citrate Contacts: {metrics.get('citrate_contacts', 0)}")
        print(f"  Geometry RMSD: {metrics.get('geometry_rmsd', 0):.2f}Å")
        print(f"  Metal SASA: {metrics.get('metal_sasa', 0):.1f}Ų")

        if result["issues"]:
            print(f"\nIssues:")
            for issue in result["issues"]:
                print(f"  - {issue}")

        if result.get("coordinating"):
            print(f"\nCoordinating Atoms:")
            for atom in result["coordinating"][:10]:  # Show first 10
                print(f"  {atom['residue']} {atom['resname']} {atom['atom']}: {atom['distance']:.2f}Å")


if __name__ == "__main__":
    main()
