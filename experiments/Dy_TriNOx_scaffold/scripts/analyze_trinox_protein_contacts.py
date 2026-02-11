"""
Analyze TriNOx-Protein interactions (more important than Dy-protein coordination).

TriNOx structure:
- O1, O2, O3: phenolate oxygens (coordinate Dy, could accept H-bonds from protein)
- N1, N3, N4: amine nitrogens (coordinate Dy, could donate/accept H-bonds)
- N2: central tertiary amine (could accept H-bonds)
- C1-C33: aromatic carbons (need hydrophobic packing)

Key interactions needed:
1. H-bonds from protein backbone/sidechains to TriNOx N/O
2. Hydrophobic contacts to TriNOx aromatic rings
3. Proper encapsulation without steric clashes
"""

import math
from pathlib import Path
from collections import defaultdict

OUTPUT_DIR = Path("G:/Github_local_repo/Banta_Lab_RFdiffusion/experiments/Dy_TriNOx_scaffold/outputs_v3")

# H-bond donor atoms in protein
PROTEIN_HBOND_DONORS = {
    "N",      # backbone
    "NE", "NH1", "NH2",  # Arg
    "ND2", "NE2",  # Asn, Gln, His
    "NZ",     # Lys
    "OG", "OG1",  # Ser, Thr
    "OH",     # Tyr
    "NE1",    # Trp
}

# H-bond acceptor atoms in protein
PROTEIN_HBOND_ACCEPTORS = {
    "O",      # backbone
    "OD1", "OD2",  # Asp
    "OE1", "OE2",  # Glu
    "OG", "OG1",  # Ser, Thr
    "OH",     # Tyr
    "ND1", "NE2",  # His (can also accept)
}

# Hydrophobic sidechains
HYDROPHOBIC_RESIDUES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "TYR", "PRO"}

# TriNOx atoms that can accept H-bonds
TRINOX_HBOND_ACCEPTORS = {"O1", "O2", "O3", "N1", "N2", "N3", "N4"}

# TriNOx aromatic carbons (need hydrophobic packing)
TRINOX_AROMATIC = {f"C{i}" for i in range(1, 34)}


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


def analyze_design(pdb_path: Path):
    with open(pdb_path, 'r') as f:
        content = f.read()

    atoms = parse_pdb(content)

    # Separate atoms
    trinox_atoms = [a for a in atoms if a["res_name"] == "UNL"]
    protein_aa = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
                  "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    protein_atoms = [a for a in atoms if a["res_name"] in protein_aa]

    # Remove duplicate TriNOx atoms (keep first occurrence of each atom name)
    seen_trinox = set()
    unique_trinox = []
    for a in trinox_atoms:
        if a["name"] not in seen_trinox:
            unique_trinox.append(a)
            seen_trinox.add(a["name"])
    trinox_atoms = unique_trinox

    results = {
        "design": pdb_path.stem,
        "trinox_atoms": len(trinox_atoms),
        "protein_residues": len(set((a["chain"], a["res_num"]) for a in protein_atoms)),
    }

    # 1. Find H-bonds from protein to TriNOx acceptors
    hbond_cutoff = 3.5  # Å
    hbonds_to_trinox = []

    trinox_acceptors = [a for a in trinox_atoms if a["name"] in TRINOX_HBOND_ACCEPTORS]
    protein_donors = [a for a in protein_atoms if a["name"] in PROTEIN_HBOND_DONORS]

    for ta in trinox_acceptors:
        for pa in protein_donors:
            d = distance(ta, pa)
            if d <= hbond_cutoff:
                hbonds_to_trinox.append({
                    "trinox_atom": ta["name"],
                    "protein": f"{pa['res_name']}{pa['res_num']}",
                    "protein_atom": pa["name"],
                    "distance": round(d, 2),
                })

    results["hbonds_to_trinox"] = hbonds_to_trinox
    results["hbond_count"] = len(hbonds_to_trinox)

    # Count unique H-bond acceptors on TriNOx that are satisfied
    satisfied_acceptors = set(hb["trinox_atom"] for hb in hbonds_to_trinox)
    results["trinox_acceptors_satisfied"] = len(satisfied_acceptors)
    results["trinox_acceptors_total"] = len(TRINOX_HBOND_ACCEPTORS)

    # 2. Find hydrophobic contacts to TriNOx aromatic carbons
    hydrophobic_cutoff = 4.5  # Å
    hydrophobic_contacts = []

    trinox_carbons = [a for a in trinox_atoms if a["name"].startswith("C")]
    hydrophobic_atoms = [a for a in protein_atoms if a["res_name"] in HYDROPHOBIC_RESIDUES]

    contacted_carbons = set()
    for tc in trinox_carbons:
        for pa in hydrophobic_atoms:
            d = distance(tc, pa)
            if d <= hydrophobic_cutoff:
                contacted_carbons.add(tc["name"])
                hydrophobic_contacts.append({
                    "trinox_carbon": tc["name"],
                    "protein": f"{pa['res_name']}{pa['res_num']}",
                    "distance": round(d, 2),
                })

    results["hydrophobic_contacts"] = len(hydrophobic_contacts)
    results["carbons_contacted"] = len(contacted_carbons)
    results["carbons_total"] = len(trinox_carbons)
    results["carbon_coverage"] = round(len(contacted_carbons) / max(len(trinox_carbons), 1), 2)

    # 3. Identify contacting residue types
    contact_cutoff = 4.5
    contacting_residues = defaultdict(int)
    for ta in trinox_atoms:
        for pa in protein_atoms:
            if distance(ta, pa) <= contact_cutoff:
                contacting_residues[pa["res_name"]] += 1

    results["contacting_residue_types"] = dict(contacting_residues)

    # 4. Check for aromatic stacking (Phe, Tyr, Trp near TriNOx rings)
    aromatic_residues = [a for a in protein_atoms if a["res_name"] in ["PHE", "TYR", "TRP"]]
    aromatic_contacts = 0
    for ta in trinox_carbons:
        for pa in aromatic_residues:
            if distance(ta, pa) <= 5.0:
                aromatic_contacts += 1
                break

    results["aromatic_stacking_contacts"] = aromatic_contacts

    # 5. Overall quality assessment
    quality_score = 0
    issues = []

    # H-bonds (max 30 points)
    if len(hbonds_to_trinox) >= 4:
        quality_score += 30
    elif len(hbonds_to_trinox) >= 2:
        quality_score += 15
    elif len(hbonds_to_trinox) >= 1:
        quality_score += 5
    else:
        issues.append("No H-bonds to TriNOx")

    # Hydrophobic coverage (max 30 points)
    coverage = results["carbon_coverage"]
    if coverage >= 0.7:
        quality_score += 30
    elif coverage >= 0.5:
        quality_score += 20
    elif coverage >= 0.3:
        quality_score += 10
    else:
        issues.append(f"Poor hydrophobic coverage ({coverage:.0%})")

    # Aromatic stacking (max 20 points)
    if aromatic_contacts >= 5:
        quality_score += 20
    elif aromatic_contacts >= 2:
        quality_score += 10
    else:
        issues.append("Few aromatic stacking interactions")

    # Residue diversity (max 20 points) - penalize all-Ala
    ala_contacts = contacting_residues.get("ALA", 0)
    total_contacts = sum(contacting_residues.values())
    if total_contacts > 0:
        ala_ratio = ala_contacts / total_contacts
        if ala_ratio < 0.3:
            quality_score += 20
        elif ala_ratio < 0.5:
            quality_score += 10
        else:
            issues.append(f"High Ala in contacts ({ala_ratio:.0%})")

    results["quality_score"] = quality_score
    results["quality_issues"] = issues
    results["quality_rating"] = (
        "excellent" if quality_score >= 80 else
        "good" if quality_score >= 60 else
        "acceptable" if quality_score >= 40 else
        "poor"
    )

    return results


def main():
    print("=" * 80)
    print("TRINOX-PROTEIN INTERACTION ANALYSIS")
    print("Focus: H-bonds and hydrophobic contacts to TriNOx ligand")
    print("=" * 80)

    pdb_files = sorted(OUTPUT_DIR.glob("v3_*.pdb"))
    print(f"\nAnalyzing {len(pdb_files)} v3 designs")

    all_results = []

    for pdb_file in pdb_files:
        result = analyze_design(pdb_file)
        all_results.append(result)

        print(f"\n{'-'*60}")
        print(f"{result['design']}")
        print(f"{'-'*60}")

        print(f"\n  H-bonds to TriNOx: {result['hbond_count']}")
        print(f"    Acceptors satisfied: {result['trinox_acceptors_satisfied']}/{result['trinox_acceptors_total']}")
        if result['hbonds_to_trinox']:
            for hb in result['hbonds_to_trinox'][:5]:
                print(f"      {hb['trinox_atom']} <- {hb['protein']} {hb['protein_atom']}: {hb['distance']}A")
            if len(result['hbonds_to_trinox']) > 5:
                print(f"      ... and {len(result['hbonds_to_trinox']) - 5} more")

        print(f"\n  Hydrophobic coverage: {result['carbon_coverage']:.0%}")
        print(f"    Carbons contacted: {result['carbons_contacted']}/{result['carbons_total']}")
        print(f"    Aromatic stacking: {result['aromatic_stacking_contacts']} contacts")

        print(f"\n  Contacting residues:")
        sorted_contacts = sorted(result['contacting_residue_types'].items(), key=lambda x: -x[1])
        for res, count in sorted_contacts[:5]:
            print(f"    {res}: {count}")

        print(f"\n  Quality: {result['quality_score']}/100 ({result['quality_rating']})")
        if result['quality_issues']:
            for issue in result['quality_issues']:
                print(f"    - {issue}")

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    # Rank by quality
    ranked = sorted(all_results, key=lambda x: x["quality_score"], reverse=True)

    print(f"\n{'Design':<20} {'Score':>6} {'Rating':<10} {'H-bonds':>8} {'HphobCov':>9} {'AromStack':>10}")
    print("-" * 75)
    for r in ranked:
        print(f"{r['design']:<20} {r['quality_score']:>6} {r['quality_rating']:<10} "
              f"{r['hbond_count']:>8} {r['carbon_coverage']:>8.0%} {r['aromatic_stacking_contacts']:>10}")

    # Quality distribution
    excellent = sum(1 for r in all_results if r["quality_rating"] == "excellent")
    good = sum(1 for r in all_results if r["quality_rating"] == "good")
    acceptable = sum(1 for r in all_results if r["quality_rating"] == "acceptable")
    poor = sum(1 for r in all_results if r["quality_rating"] == "poor")

    print(f"\nQuality distribution:")
    print(f"  Excellent: {excellent}")
    print(f"  Good: {good}")
    print(f"  Acceptable: {acceptable}")
    print(f"  Poor: {poor}")

    # Key finding
    best = ranked[0]
    print(f"\nBest design: {best['design']}")
    print(f"  H-bonds to TriNOx: {best['hbond_count']}")
    print(f"  Hydrophobic coverage: {best['carbon_coverage']:.0%}")


if __name__ == "__main__":
    main()
