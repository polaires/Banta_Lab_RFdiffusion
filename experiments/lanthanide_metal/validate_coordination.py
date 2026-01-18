#!/usr/bin/env python3
"""
Validation script for lanthanide coordination geometry.

Analyzes PDB files from test runs and validates:
1. Coordination number (target: 8-9 for Tb)
2. Carboxylate donor count (target: >= 6)
3. Bond distances (target: 2.3-2.5 A)
4. Chain contribution (target: 4+4 for EF-hand)
5. Geometry quality

Usage:
    python validate_coordination.py design.pdb
    python validate_coordination.py outputs/ef_hand_asp_*/ --batch
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional
import re


def parse_pdb_atoms(pdb_content: str) -> Tuple[Dict, List[Dict]]:
    """Parse PDB content and extract metal and coordinating atoms."""
    metal_atom = None
    coord_atoms = []

    # Known lanthanides and metals
    metals = {"TB", "GD", "EU", "LA", "YB", "ZN", "FE", "CA", "MG", "MN", "CO", "NI", "CU"}

    # Coordinating atom names
    coord_atom_names = {
        "OD1", "OD2",  # Asp carboxylate
        "OE1", "OE2",  # Glu carboxylate
        "NE2", "ND1",  # His imidazole
        "SG",          # Cys thiol
        "OG", "OG1",   # Ser/Thr hydroxyl
        "NZ",          # Lys amino
        "O",           # Backbone carbonyl
        "OH",          # Tyr hydroxyl
    }

    for line in pdb_content.split("\n"):
        if line.startswith("HETATM") or line.startswith("ATOM"):
            try:
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21]
                res_num = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                # Check if this is the metal
                if line.startswith("HETATM") and res_name in metals:
                    metal_atom = {
                        "name": atom_name,
                        "resname": res_name,
                        "chain": chain,
                        "resnum": res_num,
                        "coords": (x, y, z),
                    }
                # Check if this is a potential coordinating atom
                elif atom_name in coord_atom_names:
                    coord_atoms.append({
                        "atom": atom_name,
                        "resname": res_name,
                        "chain": chain,
                        "resnum": res_num,
                        "coords": (x, y, z),
                    })
            except (ValueError, IndexError):
                continue

    return metal_atom, coord_atoms


def calculate_distance(p1: Tuple[float, float, float], p2: Tuple[float, float, float]) -> float:
    """Calculate Euclidean distance between two points."""
    return ((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2) ** 0.5


def analyze_coordination(
    pdb_content: str,
    max_coord_distance: float = 3.5
) -> Dict[str, Any]:
    """Analyze metal coordination geometry from PDB content."""
    metal_atom, coord_atoms = parse_pdb_atoms(pdb_content)

    if not metal_atom:
        return {"error": "No metal atom found in PDB"}

    metal_pos = metal_atom["coords"]

    # Find atoms within coordination distance
    coordinating = []
    for atom in coord_atoms:
        dist = calculate_distance(metal_pos, atom["coords"])
        if dist <= max_coord_distance:
            atom["distance"] = round(dist, 3)
            coordinating.append(atom)

    # Sort by distance
    coordinating.sort(key=lambda x: x["distance"])

    # Analyze by type
    carboxylate_donors = [a for a in coordinating if a["atom"] in ["OD1", "OD2", "OE1", "OE2"]]
    asp_donors = [a for a in coordinating if a["resname"] == "ASP" and a["atom"] in ["OD1", "OD2"]]
    glu_donors = [a for a in coordinating if a["resname"] == "GLU" and a["atom"] in ["OE1", "OE2"]]
    his_donors = [a for a in coordinating if a["atom"] in ["NE2", "ND1"]]
    cys_donors = [a for a in coordinating if a["atom"] == "SG"]

    # Analyze by chain
    chain_a = [a for a in coordinating if a["chain"] == "A"]
    chain_b = [a for a in coordinating if a["chain"] == "B"]

    # Calculate average distance
    if coordinating:
        avg_distance = sum(a["distance"] for a in coordinating) / len(coordinating)
    else:
        avg_distance = None

    # Unique residues contributing
    unique_residues = set((a["chain"], a["resnum"]) for a in coordinating)

    return {
        "metal": metal_atom,
        "coordination_number": len(coordinating),
        "unique_residues": len(unique_residues),
        "average_distance": round(avg_distance, 3) if avg_distance else None,
        "carboxylate_count": len(carboxylate_donors),
        "asp_count": len(asp_donors),
        "glu_count": len(glu_donors),
        "his_count": len(his_donors),
        "cys_count": len(cys_donors),
        "chain_a_donors": len(chain_a),
        "chain_b_donors": len(chain_b),
        "coordinating_atoms": coordinating,
        "by_residue": sorted(list(unique_residues)),
    }


def validate_lanthanide(analysis: Dict[str, Any]) -> Dict[str, Any]:
    """Validate coordination against lanthanide requirements."""
    issues = []
    score = 100

    # Check coordination number (target: 8-9)
    coord_num = analysis.get("coordination_number", 0)
    if coord_num < 6:
        issues.append(f"Low coordination: {coord_num} (target: 8-9)")
        score -= 30
    elif coord_num < 8:
        issues.append(f"Suboptimal coordination: {coord_num} (target: 8-9)")
        score -= 10

    # Check carboxylate donors (target: >= 6)
    carboxylate = analysis.get("carboxylate_count", 0)
    if carboxylate < 4:
        issues.append(f"Few carboxylate donors: {carboxylate} (target: >= 6)")
        score -= 25
    elif carboxylate < 6:
        issues.append(f"Moderate carboxylate donors: {carboxylate} (target: >= 6)")
        score -= 10

    # Check average distance (target: 2.3-2.5 A)
    avg_dist = analysis.get("average_distance")
    if avg_dist:
        if avg_dist < 2.0:
            issues.append(f"Bonds too short: {avg_dist} A (target: 2.3-2.5)")
            score -= 15
        elif avg_dist > 3.0:
            issues.append(f"Bonds too long: {avg_dist} A (target: 2.3-2.5)")
            score -= 15
        elif not (2.2 <= avg_dist <= 2.6):
            issues.append(f"Suboptimal bond distance: {avg_dist} A (target: 2.3-2.5)")
            score -= 5

    # Check chain balance (target: 4+4 for EF-hand)
    chain_a = analysis.get("chain_a_donors", 0)
    chain_b = analysis.get("chain_b_donors", 0)
    if chain_a == 0 or chain_b == 0:
        issues.append(f"Asymmetric contribution: A={chain_a}, B={chain_b}")
        score -= 20
    elif abs(chain_a - chain_b) > 3:
        issues.append(f"Unbalanced chains: A={chain_a}, B={chain_b}")
        score -= 10

    return {
        "score": max(0, score),
        "pass": score >= 70,
        "issues": issues,
    }


def format_report(analysis: Dict[str, Any], validation: Dict[str, Any]) -> str:
    """Format a human-readable report."""
    lines = []
    lines.append("="*60)
    lines.append("LANTHANIDE COORDINATION VALIDATION REPORT")
    lines.append("="*60)

    # Metal info
    metal = analysis.get("metal", {})
    lines.append(f"\nMetal: {metal.get('resname', 'N/A')}")
    lines.append(f"Position: {metal.get('coords', 'N/A')}")

    # Coordination summary
    lines.append(f"\n--- Coordination Summary ---")
    lines.append(f"Coordination number: {analysis.get('coordination_number', 0)}")
    lines.append(f"Unique residues: {analysis.get('unique_residues', 0)}")
    lines.append(f"Average distance: {analysis.get('average_distance', 'N/A')} A")

    # Donor breakdown
    lines.append(f"\n--- Donor Breakdown ---")
    lines.append(f"Carboxylate (OD/OE): {analysis.get('carboxylate_count', 0)}")
    lines.append(f"  - Asp: {analysis.get('asp_count', 0)}")
    lines.append(f"  - Glu: {analysis.get('glu_count', 0)}")
    lines.append(f"Histidine (NE2/ND1): {analysis.get('his_count', 0)}")
    lines.append(f"Cysteine (SG): {analysis.get('cys_count', 0)}")

    # Chain contribution
    lines.append(f"\n--- Chain Contribution ---")
    lines.append(f"Chain A: {analysis.get('chain_a_donors', 0)} donors")
    lines.append(f"Chain B: {analysis.get('chain_b_donors', 0)} donors")

    # Coordinating residues
    lines.append(f"\n--- Coordinating Residues ---")
    for chain, resnum in sorted(analysis.get("by_residue", [])):
        lines.append(f"  {chain}{resnum}")

    # Validation results
    lines.append(f"\n--- Validation ---")
    lines.append(f"Score: {validation.get('score', 0)}/100")
    lines.append(f"Pass: {'YES' if validation.get('pass') else 'NO'}")

    if validation.get("issues"):
        lines.append(f"\nIssues:")
        for issue in validation["issues"]:
            lines.append(f"  - {issue}")

    lines.append("="*60)
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Validate lanthanide coordination")
    parser.add_argument("path", help="PDB file or directory to analyze")
    parser.add_argument("--batch", action="store_true", help="Process all PDB files in directory")
    parser.add_argument("--max-distance", type=float, default=3.5, help="Max coordination distance")
    parser.add_argument("--json", action="store_true", help="Output as JSON")

    args = parser.parse_args()

    path = Path(args.path)

    if args.batch or path.is_dir():
        # Process all PDB files
        pdb_files = list(path.glob("**/*.pdb"))
        results = []

        for pdb_file in pdb_files:
            with open(pdb_file) as f:
                pdb_content = f.read()

            analysis = analyze_coordination(pdb_content, args.max_distance)
            validation = validate_lanthanide(analysis)

            results.append({
                "file": str(pdb_file),
                "analysis": analysis,
                "validation": validation,
            })

        if args.json:
            print(json.dumps(results, indent=2))
        else:
            for r in results:
                print(f"\nFile: {r['file']}")
                print(format_report(r["analysis"], r["validation"]))
    else:
        # Single file
        with open(path) as f:
            pdb_content = f.read()

        analysis = analyze_coordination(pdb_content, args.max_distance)
        validation = validate_lanthanide(analysis)

        if args.json:
            print(json.dumps({"analysis": analysis, "validation": validation}, indent=2))
        else:
            print(format_report(analysis, validation))

    return 0


if __name__ == "__main__":
    sys.exit(main())
