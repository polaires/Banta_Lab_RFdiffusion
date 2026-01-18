#!/usr/bin/env python3
"""
Deep analysis of coordination geometry.

Checks:
1. Metal-Oxygen distances (should be ~2.3-2.6 A)
2. O-O distances between coordinating oxygens (should be ~2.8-3.5 A, not clashing)
3. O-M-O angles (should match coordination geometry)
4. CA-CA distances (should be realistic for protein scaffold)
5. Overall coordination polyhedron shape
"""

import sys
import math
import numpy as np
from pathlib import Path
from itertools import combinations

SCRIPT_DIR = Path(__file__).parent
BACKEND_DIR = SCRIPT_DIR.parent.parent / "backend" / "serverless"
sys.path.insert(0, str(BACKEND_DIR))

from lanthanide_templates import (
    generate_template_from_library,
    generate_parametric_template,
    COORDINATION_GEOMETRIES,
)


def parse_pdb(pdb_content):
    """Parse PDB to extract atoms."""
    atoms = []
    metal_pos = None

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

                if res_name in ["TB", "GD", "EU", "LA", "CE", "YB"]:
                    metal_pos = np.array([x, y, z])
                else:
                    atoms.append({
                        "name": atom_name,
                        "resname": res_name,
                        "chain": chain,
                        "resnum": res_num,
                        "coords": np.array([x, y, z]),
                    })
            except (ValueError, IndexError):
                continue

    return atoms, metal_pos


def distance(p1, p2):
    """Euclidean distance."""
    return np.linalg.norm(p1 - p2)


def angle_at_vertex(p1, vertex, p2):
    """Angle at vertex between p1-vertex-p2."""
    v1 = p1 - vertex
    v2 = p2 - vertex

    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    cos_angle = np.clip(cos_angle, -1, 1)
    return math.degrees(math.acos(cos_angle))


def analyze_coordination_geometry(pdb_content, template_name, cutoff=3.0):
    """Comprehensive analysis of coordination geometry."""
    atoms, metal_pos = parse_pdb(pdb_content)

    if metal_pos is None:
        print(f"  ERROR: No metal found in {template_name}")
        return

    # Find coordinating oxygens
    coord_oxygens = []
    for atom in atoms:
        if atom["name"] in ["OD1", "OD2", "OE1", "OE2"]:
            dist = distance(atom["coords"], metal_pos)
            if dist <= cutoff:
                coord_oxygens.append({
                    **atom,
                    "metal_dist": dist,
                })

    # Find CA atoms
    ca_atoms = [a for a in atoms if a["name"] == "CA"]

    print(f"\n{'='*70}")
    print(f"TEMPLATE: {template_name}")
    print(f"{'='*70}")

    # 1. Metal-Oxygen distances
    print(f"\n1. METAL-OXYGEN DISTANCES (should be 2.3-2.6 A)")
    print(f"   {'Residue':<15} {'Oxygen':<6} {'Distance':<10} {'Status'}")
    print(f"   {'-'*45}")

    for o in coord_oxygens:
        status = "[OK]" if 2.0 <= o["metal_dist"] <= 2.8 else "[!] OUT OF RANGE"
        print(f"   {o['resname']} {o['chain']}{o['resnum']:<10} {o['name']:<6} {o['metal_dist']:<10.2f} {status}")

    avg_dist = np.mean([o["metal_dist"] for o in coord_oxygens]) if coord_oxygens else 0
    print(f"\n   Average: {avg_dist:.2f} A, Count: {len(coord_oxygens)}")

    # 2. O-O distances (inter-oxygen)
    print(f"\n2. OXYGEN-OXYGEN DISTANCES (should be 2.8-3.5 A, no clashes < 2.4 A)")
    print(f"   {'Pair':<30} {'Distance':<10} {'Status'}")
    print(f"   {'-'*55}")

    o_o_distances = []
    clashes = 0
    too_far = 0

    for (o1, o2) in combinations(coord_oxygens, 2):
        dist = distance(o1["coords"], o2["coords"])
        o_o_distances.append(dist)

        label1 = f"{o1['resname']}{o1['chain']}{o1['resnum']}.{o1['name']}"
        label2 = f"{o2['resname']}{o2['chain']}{o2['resnum']}.{o2['name']}"

        if dist < 2.4:
            status = "[!] CLASH!"
            clashes += 1
        elif dist > 5.0:
            status = "[!] FAR"
            too_far += 1
        else:
            status = "[OK]"

        print(f"   {label1} - {label2:<15} {dist:<10.2f} {status}")

    if o_o_distances:
        print(f"\n   Min: {min(o_o_distances):.2f} A, Max: {max(o_o_distances):.2f} A, Avg: {np.mean(o_o_distances):.2f} A")
        print(f"   Clashes (<2.4A): {clashes}, Too far (>5.0A): {too_far}")

    # 3. O-M-O angles
    print(f"\n3. O-METAL-O ANGLES (SAP CN=8: expect 70.5°, 109.5°, 141.5°)")
    print(f"   {'Pair':<30} {'Angle':<10}")
    print(f"   {'-'*45}")

    angles = []
    for (o1, o2) in combinations(coord_oxygens, 2):
        ang = angle_at_vertex(o1["coords"], metal_pos, o2["coords"])
        angles.append(ang)

        label1 = f"{o1['resname']}{o1['chain']}{o1['resnum']}.{o1['name']}"
        label2 = f"{o2['resname']}{o2['chain']}{o2['resnum']}.{o2['name']}"
        print(f"   {label1} - {label2:<15} {ang:<10.1f}°")

    if angles:
        print(f"\n   Min: {min(angles):.1f}°, Max: {max(angles):.1f}°")

        # Ideal SAP angles
        ideal_angles = [70.5, 109.5, 141.5, 180.0]
        print(f"   Ideal SAP angles: {ideal_angles}")

    # 4. CA-CA distances
    print(f"\n4. CA-CA DISTANCES (should be 8-15 A for interface design)")
    print(f"   {'Pair':<25} {'Distance':<10}")
    print(f"   {'-'*40}")

    ca_distances = []
    for (ca1, ca2) in combinations(ca_atoms, 2):
        dist = distance(ca1["coords"], ca2["coords"])
        ca_distances.append(dist)

        label1 = f"{ca1['resname']}{ca1['chain']}{ca1['resnum']}"
        label2 = f"{ca2['resname']}{ca2['chain']}{ca2['resnum']}"
        print(f"   {label1} - {label2:<12} {dist:<10.2f}")

    if ca_distances:
        print(f"\n   Min: {min(ca_distances):.2f} A, Max: {max(ca_distances):.2f} A")

    # 5. Coordination polyhedron analysis
    print(f"\n5. COORDINATION POLYHEDRON ANALYSIS")

    if len(coord_oxygens) >= 4:
        # Calculate centroid
        centroid = np.mean([o["coords"] for o in coord_oxygens], axis=0)
        dist_to_metal = distance(centroid, metal_pos)
        print(f"   Centroid of O atoms: {centroid}")
        print(f"   Distance from centroid to metal: {dist_to_metal:.2f} A (should be ~0)")

        # Check symmetry - distances from centroid to each O
        radii = [distance(o["coords"], centroid) for o in coord_oxygens]
        print(f"   Radii from centroid: min={min(radii):.2f}, max={max(radii):.2f}, std={np.std(radii):.2f}")
        print(f"   (Low std = more symmetric)")

    return {
        "coord_count": len(coord_oxygens),
        "avg_metal_dist": avg_dist,
        "o_o_distances": o_o_distances,
        "clashes": clashes,
        "angles": angles,
        "ca_distances": ca_distances,
    }


def main():
    print("=" * 70)
    print("DEEP COORDINATION GEOMETRY ANALYSIS")
    print("=" * 70)

    # Analyze library templates
    templates_to_test = [
        ("caldwell_4", 3.0),
        ("ef_hand_8", 3.0),
        ("high_coord_9", 3.5),
    ]

    for template_name, cutoff in templates_to_test:
        pdb = generate_template_from_library(template_name, metal="TB")
        analyze_coordination_geometry(pdb, template_name, cutoff=cutoff)

    # Analyze parametric template
    print("\n" + "=" * 70)
    print("PARAMETRIC TEMPLATE (CN=8, 50% bidentate)")
    print("=" * 70)

    templates = generate_parametric_template(
        metal="TB",
        coordination_number=8,
        num_waters=0,
        bidentate_fraction=0.5,
        seed=42,
    )
    analyze_coordination_geometry(templates[0], "parametric_8", cutoff=3.0)


if __name__ == "__main__":
    main()
