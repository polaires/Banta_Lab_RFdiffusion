"""
Unified Analysis for Dy-TriNOx Scaffold Designs.

Uses metal_validation.py for comprehensive coordination validation.

Chemistry of Dy-TriNOx:
- Dy³⁺ coordination number: 8-9 (lanthanide)
- TriNOx ligand (UNL) provides:
  - 3 phenolate O donors (O1, O2, O3) coordinating Dy
  - 3 amine N atoms (N1, N3, N4) that may also coordinate
  - Total ligand coordination: 3-6 donors
- Protein must provide remaining donors to reach CN=8-9
  - Minimum: 1 Asp/Glu carboxylate
  - Ideal: 2-3 additional donors
"""

import sys
import os
import math
from pathlib import Path
from typing import Dict, Any, List

# Add serverless dir to path
serverless_path = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless")
sys.path.insert(0, str(serverless_path))

from metal_validation import validate_metal_ligand_complex_site

# Paths
SCRIPT_DIR = Path(__file__).parent
OUTPUTS_V3 = SCRIPT_DIR.parent / "outputs_v3"

# Dy-TriNOx specific parameters
DY_TRINOX_PARAMS = {
    "metal": "DY",
    "ligand_name": "UNL",  # TriNOx is stored as UNL
    "expected_ligand_donors": 3,  # O1, O2, O3 phenolates (conservative)
    "expected_protein_donors": 2,  # At least 2 from protein
    "metal_coordination_target": 8,  # Dy³⁺ typical CN
    "distance_cutoff": 3.0,  # Lanthanide coordination sphere
    "ligand_hbond_cutoff": 3.5,
    "ligand_hbond_acceptors": ["N1", "N2", "N3", "N4"],  # Amine N for H-bonds
}


def analyze_trinox_burial(pdb_content: str) -> Dict[str, Any]:
    """
    Analyze how well the TriNOx ligand is buried in the protein pocket.

    Returns:
        - burial_ratio: fraction of TriNOx atoms with protein contacts
        - contact_residues: list of protein residues near TriNOx
        - pocket_quality: assessment of pocket formation
    """
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

    # Separate TriNOx and protein atoms
    trinox_atoms = [a for a in atoms if a["res_name"] == "UNL"]
    protein_aa = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
                  "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    protein_atoms = [a for a in atoms if a["res_name"] in protein_aa]

    if not trinox_atoms:
        return {"error": "No TriNOx atoms found"}

    # Count TriNOx atoms with protein contacts
    contact_cutoff = 4.5  # Å
    trinox_with_contacts = 0
    contact_residues = set()

    for ta in trinox_atoms:
        has_contact = False
        for pa in protein_atoms:
            dx = ta["x"] - pa["x"]
            dy = ta["y"] - pa["y"]
            dz = ta["z"] - pa["z"]
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)
            if dist <= contact_cutoff:
                has_contact = True
                contact_residues.add((pa["chain"], pa["res_num"], pa["res_name"]))
        if has_contact:
            trinox_with_contacts += 1

    burial_ratio = trinox_with_contacts / len(trinox_atoms)

    # Assess pocket quality
    if len(contact_residues) >= 10 and burial_ratio >= 0.8:
        pocket_quality = "excellent"
    elif len(contact_residues) >= 6 and burial_ratio >= 0.6:
        pocket_quality = "good"
    elif len(contact_residues) >= 3 and burial_ratio >= 0.3:
        pocket_quality = "partial"
    else:
        pocket_quality = "poor"

    return {
        "total_trinox_atoms": len(trinox_atoms),
        "trinox_with_contacts": trinox_with_contacts,
        "burial_ratio": round(burial_ratio, 2),
        "contact_residues": len(contact_residues),
        "contact_residue_list": sorted(contact_residues),
        "pocket_quality": pocket_quality,
    }


def analyze_sequence_composition(pdb_content: str) -> Dict[str, Any]:
    """Analyze sequence composition for metal-binding residues."""
    # Extract sequence from PDB
    residues = {}
    for line in pdb_content.split('\n'):
        if line.startswith("ATOM"):
            try:
                chain = line[21]
                res_num = int(line[22:26].strip())
                res_name = line[17:20].strip()
                residues[(chain, res_num)] = res_name
            except:
                pass

    # Count amino acid types
    aa_3to1 = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }

    sequence = "".join([aa_3to1.get(res_name, "X") for (chain, res_num), res_name in sorted(residues.items())])
    length = len(sequence)

    if length == 0:
        return {"error": "No protein sequence found"}

    # Count key residues
    counts = {
        "D": sequence.count("D"),  # Asp
        "E": sequence.count("E"),  # Glu
        "H": sequence.count("H"),  # His
        "A": sequence.count("A"),  # Ala
    }

    # Check for problematic composition
    ala_ratio = counts["A"] / length
    carbox_count = counts["D"] + counts["E"]

    issues = []
    if ala_ratio > 0.30:
        issues.append(f"High alanine content ({ala_ratio:.0%}) - may indicate poor designability")
    if carbox_count < 2:
        issues.append(f"Low carboxylate content ({carbox_count}) - may lack coordination capacity")

    return {
        "length": length,
        "sequence": sequence,
        "asp_count": counts["D"],
        "glu_count": counts["E"],
        "carboxylate_total": carbox_count,
        "his_count": counts["H"],
        "ala_count": counts["A"],
        "ala_ratio": round(ala_ratio, 2),
        "issues": issues,
    }


def run_unified_analysis(pdb_path: Path) -> Dict[str, Any]:
    """Run comprehensive unified analysis on a design."""
    with open(pdb_path, 'r') as f:
        pdb_content = f.read()

    result = {
        "design": pdb_path.stem,
        "path": str(pdb_path),
    }

    # 1. Metal-ligand complex validation
    coordination = validate_metal_ligand_complex_site(
        pdb_content=pdb_content,
        **DY_TRINOX_PARAMS
    )
    result["coordination"] = coordination

    # 2. TriNOx burial analysis
    burial = analyze_trinox_burial(pdb_content)
    result["burial"] = burial

    # 3. Sequence composition
    composition = analyze_sequence_composition(pdb_content)
    result["composition"] = composition

    # 4. Overall quality assessment
    quality_issues = []
    quality_score = 0

    # Coordination checks
    if coordination.get("success"):
        quality_score += 40
    else:
        if coordination.get("issues"):
            quality_issues.extend(coordination["issues"])

    if coordination.get("ligand_coordination", 0) >= 3:
        quality_score += 20  # TriNOx coordinating
    else:
        quality_issues.append(f"TriNOx-Dy coordination broken: {coordination.get('ligand_coordination', 0)}/3")

    if coordination.get("protein_coordination", 0) >= 1:
        quality_score += 20  # At least 1 protein donor
    else:
        quality_issues.append("No protein coordination to Dy")

    # Burial checks
    if burial.get("pocket_quality") in ["excellent", "good"]:
        quality_score += 15
    elif burial.get("pocket_quality") == "partial":
        quality_score += 5
        quality_issues.append("TriNOx only partially buried")
    else:
        quality_issues.append("TriNOx not buried in protein pocket")

    # Composition checks
    if composition.get("issues"):
        quality_issues.extend(composition["issues"])
    else:
        quality_score += 5

    result["quality_score"] = quality_score
    result["quality_issues"] = quality_issues
    result["quality_rating"] = (
        "excellent" if quality_score >= 90 else
        "good" if quality_score >= 70 else
        "acceptable" if quality_score >= 50 else
        "poor"
    )

    return result


def main():
    print("=" * 80)
    print("DY-TRINOX UNIFIED ANALYSIS")
    print("=" * 80)

    # Find all v3 designs
    pdb_files = sorted(OUTPUTS_V3.glob("v3_*.pdb"))
    print(f"\nFound {len(pdb_files)} v3 designs")

    if not pdb_files:
        print("ERROR: No v3 designs found!")
        return

    all_results = []

    for pdb_file in pdb_files:
        print(f"\n{'='*60}")
        print(f"Analyzing: {pdb_file.stem}")
        print('='*60)

        result = run_unified_analysis(pdb_file)
        all_results.append(result)

        # Print summary
        coord = result["coordination"]
        burial = result["burial"]
        comp = result["composition"]

        print(f"\n[COORDINATION]")
        print(f"  TriNOx -> Dy: {coord.get('ligand_coordination', 0)}/{DY_TRINOX_PARAMS['expected_ligand_donors']} donors")
        if coord.get("ligand_donors"):
            for d in coord["ligand_donors"]:
                print(f"    {d['atom']}: {d['distance']}A")

        print(f"  Protein -> Dy: {coord.get('protein_coordination', 0)}/{DY_TRINOX_PARAMS['expected_protein_donors']} donors")
        if coord.get("donor_residues"):
            for d in coord["donor_residues"]:
                print(f"    {d['chain']}{d['resnum']} {d['resname']} {d['atom']}: {d['distance']}A")

        print(f"  Total coordination: {coord.get('coordination_number', 0)}/{DY_TRINOX_PARAMS['metal_coordination_target']}")

        print(f"\n[BURIAL]")
        print(f"  Burial ratio: {burial.get('burial_ratio', 0):.0%}")
        print(f"  Contact residues: {burial.get('contact_residues', 0)}")
        print(f"  Pocket quality: {burial.get('pocket_quality', 'unknown')}")

        print(f"\n[COMPOSITION]")
        print(f"  Length: {comp.get('length', 0)} aa")
        print(f"  Carboxylate (D+E): {comp.get('carboxylate_total', 0)}")
        print(f"  Alanine: {comp.get('ala_count', 0)} ({comp.get('ala_ratio', 0):.0%})")

        print(f"\n[QUALITY]")
        print(f"  Score: {result['quality_score']}/100")
        print(f"  Rating: {result['quality_rating']}")
        if result["quality_issues"]:
            print("  Issues:")
            for issue in result["quality_issues"]:
                print(f"    - {issue}")

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    # Rank by quality score
    ranked = sorted(all_results, key=lambda x: x["quality_score"], reverse=True)

    print(f"\n{'Design':<20} {'Score':>6} {'Rating':<12} {'Coord':>6} {'Burial':>8} {'Issues'}")
    print("-" * 80)
    for r in ranked:
        coord = r["coordination"]
        burial = r["burial"]
        print(f"{r['design']:<20} {r['quality_score']:>6} {r['quality_rating']:<12} "
              f"{coord.get('coordination_number', 0):>6} {burial.get('burial_ratio', 0):>7.0%} "
              f"{len(r['quality_issues'])}")

    # Best design
    best = ranked[0]
    print(f"\nBest design: {best['design']}")
    print(f"  Quality score: {best['quality_score']}/100")
    print(f"  Quality rating: {best['quality_rating']}")

    # Count by quality
    excellent = sum(1 for r in all_results if r["quality_rating"] == "excellent")
    good = sum(1 for r in all_results if r["quality_rating"] == "good")
    acceptable = sum(1 for r in all_results if r["quality_rating"] == "acceptable")
    poor = sum(1 for r in all_results if r["quality_rating"] == "poor")

    print(f"\nQuality distribution:")
    print(f"  Excellent: {excellent}")
    print(f"  Good: {good}")
    print(f"  Acceptable: {acceptable}")
    print(f"  Poor: {poor}")


if __name__ == "__main__":
    main()
