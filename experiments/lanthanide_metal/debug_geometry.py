#!/usr/bin/env python3
"""
Debug script to analyze the geometry issues in template generation.

Examines bond distances, angles, and atom positions to identify
why the structures look unrealistic.
"""

import sys
import math
from pathlib import Path

# Add backend to path
SCRIPT_DIR = Path(__file__).parent
BACKEND_DIR = SCRIPT_DIR.parent.parent / "backend" / "serverless"
sys.path.insert(0, str(BACKEND_DIR))

from lanthanide_templates import generate_template_from_library


def parse_pdb_atoms(pdb_content):
    """Parse all atoms from PDB content."""
    atoms = {}
    for line in pdb_content.split("\n"):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21]
                res_num = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                key = f"{chain}{res_num}_{atom_name}"
                atoms[key] = {
                    "name": atom_name,
                    "resname": res_name,
                    "chain": chain,
                    "resnum": res_num,
                    "coords": (x, y, z),
                }
            except (ValueError, IndexError):
                continue
    return atoms


def distance(p1, p2):
    """Calculate Euclidean distance."""
    return math.sqrt(sum((a - b)**2 for a, b in zip(p1, p2)))


def angle(p1, p2, p3):
    """Calculate angle at p2 between p1-p2-p3."""
    v1 = tuple(a - b for a, b in zip(p1, p2))
    v2 = tuple(a - b for a, b in zip(p3, p2))

    dot = sum(a * b for a, b in zip(v1, v2))
    mag1 = math.sqrt(sum(a**2 for a in v1))
    mag2 = math.sqrt(sum(a**2 for a in v2))

    if mag1 < 0.01 or mag2 < 0.01:
        return 0.0

    cos_angle = dot / (mag1 * mag2)
    cos_angle = max(-1, min(1, cos_angle))
    return math.degrees(math.acos(cos_angle))


def analyze_asp_geometry(atoms, chain, resnum):
    """Analyze geometry of an ASP residue."""
    prefix = f"{chain}{resnum}_"

    required = ["N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"]
    for atom in required:
        if prefix + atom not in atoms:
            print(f"  Missing atom: {prefix}{atom}")
            return None

    # Get coordinates
    n = atoms[prefix + "N"]["coords"]
    ca = atoms[prefix + "CA"]["coords"]
    c = atoms[prefix + "C"]["coords"]
    cb = atoms[prefix + "CB"]["coords"]
    cg = atoms[prefix + "CG"]["coords"]
    od1 = atoms[prefix + "OD1"]["coords"]
    od2 = atoms[prefix + "OD2"]["coords"]

    # Expected bond lengths
    expected = {
        "CA-CB": 1.54,
        "CB-CG": 1.52,
        "CG-OD1": 1.25,
        "CG-OD2": 1.25,
        "CA-N": 1.47,
        "CA-C": 1.52,
    }

    # Actual bond lengths
    actual = {
        "CA-CB": distance(ca, cb),
        "CB-CG": distance(cb, cg),
        "CG-OD1": distance(cg, od1),
        "CG-OD2": distance(cg, od2),
        "CA-N": distance(ca, n),
        "CA-C": distance(ca, c),
    }

    print(f"\n  ASP {chain}{resnum} Bond Lengths:")
    print(f"  {'Bond':<12} {'Expected':>10} {'Actual':>10} {'Delta':>10}")
    print(f"  {'-'*44}")

    issues = []
    for bond, exp in expected.items():
        act = actual[bond]
        delta = act - exp
        flag = "[!]" if abs(delta) > 0.3 else ""
        print(f"  {bond:<12} {exp:>10.3f} {act:>10.3f} {delta:>+10.3f} {flag}")
        if abs(delta) > 0.3:
            issues.append(f"{bond}: expected {exp:.2f}, got {act:.2f}")

    # Check angles
    print(f"\n  ASP {chain}{resnum} Bond Angles:")

    angle_ca_cb_cg = angle(ca, cb, cg)
    angle_cb_cg_od1 = angle(cb, cg, od1)
    angle_cb_cg_od2 = angle(cb, cg, od2)
    angle_od1_cg_od2 = angle(od1, cg, od2)

    print(f"  CA-CB-CG:    {angle_ca_cb_cg:>6.1f} deg (expected ~111)")
    print(f"  CB-CG-OD1:   {angle_cb_cg_od1:>6.1f} deg (expected ~118)")
    print(f"  CB-CG-OD2:   {angle_cb_cg_od2:>6.1f} deg (expected ~118)")
    print(f"  OD1-CG-OD2:  {angle_od1_cg_od2:>6.1f} deg (expected ~124)")

    if abs(angle_ca_cb_cg - 111) > 30:
        issues.append(f"CA-CB-CG angle: {angle_ca_cb_cg:.1f} (expected ~111)")
    if abs(angle_od1_cg_od2 - 124) > 30:
        issues.append(f"OD1-CG-OD2 angle: {angle_od1_cg_od2:.1f} (expected ~124)")

    return issues


def analyze_glu_geometry(atoms, chain, resnum):
    """Analyze geometry of a GLU residue."""
    prefix = f"{chain}{resnum}_"

    required = ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"]
    for atom in required:
        if prefix + atom not in atoms:
            print(f"  Missing atom: {prefix}{atom}")
            return None

    # Get coordinates
    n = atoms[prefix + "N"]["coords"]
    ca = atoms[prefix + "CA"]["coords"]
    c = atoms[prefix + "C"]["coords"]
    cb = atoms[prefix + "CB"]["coords"]
    cg = atoms[prefix + "CG"]["coords"]
    cd = atoms[prefix + "CD"]["coords"]
    oe1 = atoms[prefix + "OE1"]["coords"]
    oe2 = atoms[prefix + "OE2"]["coords"]

    # Expected bond lengths
    expected = {
        "CA-CB": 1.54,
        "CB-CG": 1.52,
        "CG-CD": 1.52,
        "CD-OE1": 1.25,
        "CD-OE2": 1.25,
        "CA-N": 1.47,
        "CA-C": 1.52,
    }

    # Actual bond lengths
    actual = {
        "CA-CB": distance(ca, cb),
        "CB-CG": distance(cb, cg),
        "CG-CD": distance(cg, cd),
        "CD-OE1": distance(cd, oe1),
        "CD-OE2": distance(cd, oe2),
        "CA-N": distance(ca, n),
        "CA-C": distance(ca, c),
    }

    print(f"\n  GLU {chain}{resnum} Bond Lengths:")
    print(f"  {'Bond':<12} {'Expected':>10} {'Actual':>10} {'Delta':>10}")
    print(f"  {'-'*44}")

    issues = []
    for bond, exp in expected.items():
        act = actual[bond]
        delta = act - exp
        flag = "[!]" if abs(delta) > 0.3 else ""
        print(f"  {bond:<12} {exp:>10.3f} {act:>10.3f} {delta:>+10.3f} {flag}")
        if abs(delta) > 0.3:
            issues.append(f"{bond}: expected {exp:.2f}, got {act:.2f}")

    # Check angles
    print(f"\n  GLU {chain}{resnum} Bond Angles:")

    angle_ca_cb_cg = angle(ca, cb, cg)
    angle_cb_cg_cd = angle(cb, cg, cd)
    angle_cg_cd_oe1 = angle(cg, cd, oe1)
    angle_cg_cd_oe2 = angle(cg, cd, oe2)
    angle_oe1_cd_oe2 = angle(oe1, cd, oe2)

    print(f"  CA-CB-CG:    {angle_ca_cb_cg:>6.1f} deg (expected ~111)")
    print(f"  CB-CG-CD:    {angle_cb_cg_cd:>6.1f} deg (expected ~111)")
    print(f"  CG-CD-OE1:   {angle_cg_cd_oe1:>6.1f} deg (expected ~118)")
    print(f"  CG-CD-OE2:   {angle_cg_cd_oe2:>6.1f} deg (expected ~118)")
    print(f"  OE1-CD-OE2:  {angle_oe1_cd_oe2:>6.1f} deg (expected ~124)")

    if abs(angle_ca_cb_cg - 111) > 30:
        issues.append(f"CA-CB-CG angle: {angle_ca_cb_cg:.1f} (expected ~111)")
    if abs(angle_oe1_cd_oe2 - 124) > 30:
        issues.append(f"OE1-CD-OE2 angle: {angle_oe1_cd_oe2:.1f} (expected ~124)")

    return issues


def main():
    print("=" * 70)
    print("TEMPLATE GEOMETRY ANALYSIS")
    print("=" * 70)

    # Generate caldwell_4 template
    print("\n>>> Analyzing caldwell_4 template (4 bidentate Glu)...")
    pdb = generate_template_from_library("caldwell_4", metal="TB")

    # Save for inspection
    output_path = SCRIPT_DIR / "outputs" / "debug_caldwell_4.pdb"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(pdb)
    print(f"Saved to: {output_path}")

    atoms = parse_pdb_atoms(pdb)

    # Find all residues
    residues = set()
    for key, atom in atoms.items():
        if atom["resname"] in ["GLU", "ASP"]:
            residues.add((atom["chain"], atom["resnum"], atom["resname"]))

    all_issues = []
    for chain, resnum, resname in sorted(residues):
        print(f"\n{'='*50}")
        print(f"Residue: {resname} {chain}{resnum}")
        print("=" * 50)

        if resname == "ASP":
            issues = analyze_asp_geometry(atoms, chain, resnum)
        else:
            issues = analyze_glu_geometry(atoms, chain, resnum)

        if issues:
            all_issues.extend(issues)

    # Check metal distances
    print("\n" + "=" * 50)
    print("METAL COORDINATION DISTANCES")
    print("=" * 50)

    metal_key = None
    metal_pos = None
    for key, atom in atoms.items():
        if atom["resname"] in ["TB", "GD", "EU", "LA"]:
            metal_key = key
            metal_pos = atom["coords"]
            break

    if metal_pos:
        print(f"\nMetal position: {metal_pos}")
        print(f"\nDistances to coordinating oxygens:")

        for key, atom in atoms.items():
            if atom["name"] in ["OD1", "OD2", "OE1", "OE2"]:
                dist = distance(metal_pos, atom["coords"])
                flag = "[OK]" if 2.0 <= dist <= 3.0 else "[!]"
                print(f"  {key}: {dist:.3f} A {flag}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    if all_issues:
        print(f"\nFound {len(all_issues)} geometry issues:")
        for issue in all_issues:
            print(f"  - {issue}")
    else:
        print("\nNo major geometry issues found.")

    # Now analyze ef_hand_8 which has mixed mono/bidentate
    print("\n\n>>> Analyzing ef_hand_8 template (mixed mono/bidentate)...")
    pdb2 = generate_template_from_library("ef_hand_8", metal="TB")

    output_path2 = SCRIPT_DIR / "outputs" / "debug_ef_hand_8.pdb"
    with open(output_path2, "w") as f:
        f.write(pdb2)
    print(f"Saved to: {output_path2}")

    atoms2 = parse_pdb_atoms(pdb2)

    residues2 = set()
    for key, atom in atoms2.items():
        if atom["resname"] in ["GLU", "ASP"]:
            residues2.add((atom["chain"], atom["resnum"], atom["resname"]))

    for chain, resnum, resname in sorted(residues2):
        print(f"\n{'='*50}")
        print(f"Residue: {resname} {chain}{resnum}")
        print("=" * 50)

        if resname == "ASP":
            analyze_asp_geometry(atoms2, chain, resnum)
        else:
            analyze_glu_geometry(atoms2, chain, resnum)

    return 0


if __name__ == "__main__":
    sys.exit(main())
