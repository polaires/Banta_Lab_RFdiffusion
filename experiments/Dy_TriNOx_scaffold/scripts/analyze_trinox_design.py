"""
Analyze Dy-TriNOx protein designs for coordination quality.

Key metrics:
1. Dy-protein coordination distance (target: 2.2-2.7 Å)
2. Coordination completeness (aim for 8-9 total coordination)
3. Metal burial (>70% SASA reduction)
4. Coordination geometry angles

Usage:
    python analyze_trinox_design.py design.pdb
    python analyze_trinox_design.py outputs/*.pdb --batch
"""

import argparse
import json
import math
import os
import sys
from pathlib import Path
from dataclasses import dataclass, asdict
from typing import List, Tuple, Optional


@dataclass
class Atom:
    """Simplified atom representation."""
    serial: int
    name: str
    resname: str
    chain: str
    resnum: int
    x: float
    y: float
    z: float
    element: str


@dataclass
class CoordinationAnalysis:
    """Results of coordination analysis."""
    design_name: str
    dy_position: Tuple[float, float, float]
    trinox_coordination: int  # Should be 7
    protein_coordination: int  # Target: 1-2
    total_coordination: int  # Target: 8-9
    closest_protein_distance: float
    closest_protein_atom: str
    closest_protein_residue: str
    coordinating_residues: List[str]
    coordination_distances: List[float]
    coordination_angles: List[float]
    burial_estimate: float
    score: float
    warnings: List[str]


def parse_pdb(pdb_path: str) -> List[Atom]:
    """Parse PDB file and return list of atoms."""
    atoms = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    serial = int(line[6:11])
                    name = line[12:16].strip()
                    resname = line[17:20].strip()
                    chain = line[21:22]
                    resnum = int(line[22:26])
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    element = line[76:78].strip() if len(line) > 76 else name[0]

                    atoms.append(Atom(
                        serial=serial,
                        name=name,
                        resname=resname,
                        chain=chain,
                        resnum=resnum,
                        x=x, y=y, z=z,
                        element=element
                    ))
                except (ValueError, IndexError):
                    continue
    return atoms


def distance(a1: Atom, a2: Atom) -> float:
    """Calculate distance between two atoms."""
    return math.sqrt(
        (a1.x - a2.x)**2 +
        (a1.y - a2.y)**2 +
        (a1.z - a2.z)**2
    )


def angle(a1: Atom, center: Atom, a2: Atom) -> float:
    """Calculate angle a1-center-a2 in degrees."""
    # Vectors from center to a1 and a2
    v1 = (a1.x - center.x, a1.y - center.y, a1.z - center.z)
    v2 = (a2.x - center.x, a2.y - center.y, a2.z - center.z)

    # Dot product
    dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

    # Magnitudes
    mag1 = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
    mag2 = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)

    if mag1 < 0.001 or mag2 < 0.001:
        return 0.0

    cos_angle = max(-1.0, min(1.0, dot / (mag1 * mag2)))
    return math.degrees(math.acos(cos_angle))


def find_dy(atoms: List[Atom]) -> Optional[Atom]:
    """Find the Dysprosium atom."""
    for atom in atoms:
        if atom.element == 'DY' or atom.name.upper().startswith('DY'):
            return atom
    return None


def find_trinox_donors(atoms: List[Atom]) -> List[Atom]:
    """Find the 7 donor atoms from TriNOx ligand."""
    donors = []
    for atom in atoms:
        if atom.chain == 'L':  # Ligand chain
            # Phenolate oxygens (O1, O2, O3)
            if atom.name in ['O1', 'O2', 'O3']:
                donors.append(atom)
            # Amine nitrogens (N1, N3, N4)
            if atom.name in ['N1', 'N3', 'N4']:
                donors.append(atom)
            # Apical nitrogen (N2)
            if atom.name == 'N2':
                donors.append(atom)
    return donors


def find_protein_donors(atoms: List[Atom], dy: Atom, cutoff: float = 4.0) -> List[Tuple[Atom, float]]:
    """
    Find potential coordinating atoms from the protein.

    Returns list of (atom, distance) tuples for O/N atoms within cutoff.
    """
    candidates = []

    for atom in atoms:
        # Skip non-protein atoms (ligand/metal chains)
        if atom.chain in ['X', 'L']:
            continue

        # Only consider O and N as potential donors
        if atom.element not in ['O', 'N']:
            continue

        # Skip backbone N (usually not good coordinators)
        if atom.name == 'N' and atom.resname != 'PRO':
            continue

        dist = distance(atom, dy)
        if dist <= cutoff:
            candidates.append((atom, dist))

    # Sort by distance
    candidates.sort(key=lambda x: x[1])
    return candidates


def estimate_burial(atoms: List[Atom], dy: Atom) -> float:
    """
    Estimate metal burial based on nearby protein atoms.

    Simple heuristic: count C atoms within 6Å shell.
    More C atoms = more buried.
    """
    shell_atoms = 0
    for atom in atoms:
        if atom.chain in ['X', 'L']:
            continue
        dist = distance(atom, dy)
        if 3.0 <= dist <= 6.0 and atom.element == 'C':
            shell_atoms += 1

    # Normalize: ~20 C atoms in shell = 100% burial estimate
    return min(1.0, shell_atoms / 20.0)


def analyze_design(pdb_path: str) -> CoordinationAnalysis:
    """
    Analyze a single design for coordination quality.
    """
    design_name = Path(pdb_path).stem
    atoms = parse_pdb(pdb_path)
    warnings = []

    # Find Dy
    dy = find_dy(atoms)
    if dy is None:
        return CoordinationAnalysis(
            design_name=design_name,
            dy_position=(0, 0, 0),
            trinox_coordination=0,
            protein_coordination=0,
            total_coordination=0,
            closest_protein_distance=999.0,
            closest_protein_atom="N/A",
            closest_protein_residue="N/A",
            coordinating_residues=[],
            coordination_distances=[],
            coordination_angles=[],
            burial_estimate=0.0,
            score=0.0,
            warnings=["Dy atom not found"]
        )

    dy_pos = (dy.x, dy.y, dy.z)

    # Find TriNOx donors
    trinox_donors = find_trinox_donors(atoms)
    trinox_coord = len(trinox_donors)

    if trinox_coord != 7:
        warnings.append(f"Expected 7 TriNOx donors, found {trinox_coord}")

    # Find protein donors
    protein_candidates = find_protein_donors(atoms, dy, cutoff=4.0)

    # Count coordinating atoms (distance < 3.0 Å)
    coordinating = [(a, d) for a, d in protein_candidates if d < 3.0]
    protein_coord = len(coordinating)

    # Extract info about coordinating residues
    coord_residues = []
    coord_distances = []
    for atom, dist in coordinating:
        res_id = f"{atom.resname}{atom.resnum}:{atom.name}"
        coord_residues.append(res_id)
        coord_distances.append(round(dist, 2))

    # Closest protein atom
    if protein_candidates:
        closest_atom, closest_dist = protein_candidates[0]
        closest_name = closest_atom.name
        closest_res = f"{closest_atom.resname}{closest_atom.resnum}"
    else:
        closest_dist = 999.0
        closest_name = "N/A"
        closest_res = "N/A"

    # Calculate some coordination angles
    angles = []
    if len(coordinating) >= 2:
        for i in range(min(3, len(coordinating))):
            for j in range(i + 1, min(3, len(coordinating))):
                ang = angle(coordinating[i][0], dy, coordinating[j][0])
                angles.append(round(ang, 1))

    # Estimate burial
    burial = estimate_burial(atoms, dy)

    # Calculate score
    total_coord = trinox_coord + protein_coord

    # Scoring:
    # - Distance score: optimal at 2.2-2.5 Å
    # - Coordination score: target 8-9
    # - Burial score: higher is better

    dist_score = 0.0
    if closest_dist < 2.0:
        dist_score = 0.5  # Too close
    elif closest_dist < 2.7:
        dist_score = 1.0  # Optimal
    elif closest_dist < 3.0:
        dist_score = 0.7  # Acceptable
    elif closest_dist < 4.0:
        dist_score = 0.3  # Marginal
    else:
        dist_score = 0.0  # No coordination

    coord_score = 0.0
    if total_coord == 8:
        coord_score = 1.0
    elif total_coord == 9:
        coord_score = 0.95
    elif total_coord == 7 and protein_coord >= 1:
        coord_score = 0.8
    elif total_coord == 10:
        coord_score = 0.7

    burial_score = burial

    # Weighted average
    final_score = 0.4 * dist_score + 0.3 * coord_score + 0.3 * burial_score

    # Add warnings
    if protein_coord == 0:
        warnings.append("No protein coordination detected")
    elif protein_coord > 2:
        warnings.append(f"High protein coordination ({protein_coord}) - may be crowded")

    if burial < 0.3:
        warnings.append("Low burial - metal may be too exposed")

    return CoordinationAnalysis(
        design_name=design_name,
        dy_position=dy_pos,
        trinox_coordination=trinox_coord,
        protein_coordination=protein_coord,
        total_coordination=total_coord,
        closest_protein_distance=round(closest_dist, 2),
        closest_protein_atom=closest_name,
        closest_protein_residue=closest_res,
        coordinating_residues=coord_residues,
        coordination_distances=coord_distances,
        coordination_angles=angles,
        burial_estimate=round(burial, 2),
        score=round(final_score, 3),
        warnings=warnings
    )


def print_analysis(analysis: CoordinationAnalysis):
    """Print analysis results in a readable format."""
    print(f"\n{'='*60}")
    print(f"Design: {analysis.design_name}")
    print(f"{'='*60}")

    print(f"\nCoordination Summary:")
    print(f"  TriNOx donors:    {analysis.trinox_coordination}/7")
    print(f"  Protein donors:   {analysis.protein_coordination} (target: 1-2)")
    print(f"  Total:            {analysis.total_coordination} (target: 8-9)")

    print(f"\nClosest Protein Contact:")
    print(f"  Distance:  {analysis.closest_protein_distance} Å")
    print(f"  Atom:      {analysis.closest_protein_atom}")
    print(f"  Residue:   {analysis.closest_protein_residue}")

    if analysis.coordinating_residues:
        print(f"\nCoordinating Residues:")
        for res, dist in zip(analysis.coordinating_residues, analysis.coordination_distances):
            status = "✓" if dist < 2.7 else "~"
            print(f"  {status} {res}: {dist} Å")

    print(f"\nBurial Estimate: {analysis.burial_estimate*100:.0f}%")
    print(f"\nOverall Score: {analysis.score:.3f}")

    if analysis.warnings:
        print(f"\nWarnings:")
        for w in analysis.warnings:
            print(f"  ⚠ {w}")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze Dy-TriNOx protein designs"
    )
    parser.add_argument(
        "pdbs",
        nargs="+",
        help="PDB file(s) to analyze"
    )
    parser.add_argument(
        "--batch", "-b",
        action="store_true",
        help="Batch mode - output JSON summary"
    )
    parser.add_argument(
        "--output", "-o",
        help="Output JSON file (batch mode)"
    )
    parser.add_argument(
        "--top", "-t",
        type=int,
        default=5,
        help="Show top N designs (batch mode)"
    )

    args = parser.parse_args()

    # Analyze all designs
    results = []
    for pdb_path in args.pdbs:
        if not os.path.exists(pdb_path):
            print(f"Warning: File not found: {pdb_path}")
            continue

        analysis = analyze_design(pdb_path)
        results.append(analysis)

        if not args.batch:
            print_analysis(analysis)

    if args.batch and results:
        # Sort by score
        results.sort(key=lambda x: x.score, reverse=True)

        print(f"\n{'='*60}")
        print(f"BATCH ANALYSIS SUMMARY")
        print(f"{'='*60}")
        print(f"Total designs analyzed: {len(results)}")

        print(f"\nTop {args.top} designs:")
        print("-" * 60)
        for i, r in enumerate(results[:args.top], 1):
            print(f"{i}. {r.design_name}")
            print(f"   Score: {r.score:.3f} | Coord: {r.total_coordination} | "
                  f"Dist: {r.closest_protein_distance}Å | Burial: {r.burial_estimate*100:.0f}%")

        # Save JSON if requested
        if args.output:
            output_data = {
                "summary": {
                    "total_designs": len(results),
                    "avg_score": sum(r.score for r in results) / len(results),
                    "best_score": results[0].score if results else 0,
                    "best_design": results[0].design_name if results else None,
                },
                "designs": [asdict(r) for r in results]
            }
            with open(args.output, 'w') as f:
                json.dump(output_data, f, indent=2)
            print(f"\nSaved to: {args.output}")


if __name__ == "__main__":
    main()
