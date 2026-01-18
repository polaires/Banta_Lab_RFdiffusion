#!/usr/bin/env python3
"""
Analyze heterodimer structures to understand why approaches failed.

Checks:
1. Ligand contacts per chain
2. Chain positions relative to ligand
3. Whether chains are on opposite sides or same side
"""

import sys
import math
import numpy as np
from pathlib import Path


def extract_coords_by_chain(pdb_content: str):
    """Extract atom coordinates grouped by chain."""
    chains = {}
    ligand_coords = []

    for line in pdb_content.split("\n"):
        if line.startswith("ATOM"):
            chain = line[21]
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                if chain not in chains:
                    chains[chain] = []
                chains[chain].append([x, y, z])
            except (ValueError, IndexError):
                continue
        elif line.startswith("HETATM"):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                ligand_coords.append([x, y, z])
            except (ValueError, IndexError):
                continue

    return chains, ligand_coords


def calculate_centroid(coords):
    """Calculate centroid of coordinates."""
    if not coords:
        return None
    coords = np.array(coords)
    return np.mean(coords, axis=0)


def calculate_contacts(chain_coords, ligand_coords, cutoff=5.0):
    """Count contacts within cutoff distance."""
    contacts = 0
    chain_arr = np.array(chain_coords)
    ligand_arr = np.array(ligand_coords)

    for lc in ligand_arr:
        distances = np.linalg.norm(chain_arr - lc, axis=1)
        contacts += np.sum(distances < cutoff)

    return contacts


def analyze_heterodimer(pdb_path: str):
    """Analyze a heterodimer PDB file."""
    print(f"\n{'='*60}")
    print(f"Analyzing: {Path(pdb_path).name}")
    print(f"{'='*60}")

    with open(pdb_path, 'r') as f:
        pdb_content = f.read()

    chains, ligand_coords = extract_coords_by_chain(pdb_content)

    if not ligand_coords:
        print("ERROR: No ligand (HETATM) found!")
        return

    ligand_centroid = calculate_centroid(ligand_coords)
    print(f"\nLigand centroid: ({ligand_centroid[0]:.1f}, {ligand_centroid[1]:.1f}, {ligand_centroid[2]:.1f})")

    chain_centroids = {}
    chain_contacts = {}
    chain_vectors = {}

    for chain_id, coords in sorted(chains.items()):
        centroid = calculate_centroid(coords)
        chain_centroids[chain_id] = centroid

        # Calculate contacts
        contacts = calculate_contacts(coords, ligand_coords, cutoff=5.0)
        chain_contacts[chain_id] = contacts

        # Vector from ligand to chain
        vector = centroid - ligand_centroid
        chain_vectors[chain_id] = vector
        dist = np.linalg.norm(vector)

        print(f"\nChain {chain_id}:")
        print(f"  Centroid: ({centroid[0]:.1f}, {centroid[1]:.1f}, {centroid[2]:.1f})")
        print(f"  Distance to ligand: {dist:.1f} A")
        print(f"  Contacts with ligand (5A cutoff): {contacts}")
        print(f"  Vector from ligand: ({vector[0]:.1f}, {vector[1]:.1f}, {vector[2]:.1f})")

    # Check if chains are on opposite sides
    if 'A' in chain_vectors and 'B' in chain_vectors:
        vec_a = chain_vectors['A']
        vec_b = chain_vectors['B']

        # Normalize vectors
        norm_a = vec_a / np.linalg.norm(vec_a)
        norm_b = vec_b / np.linalg.norm(vec_b)

        # Dot product: negative = opposite sides, positive = same side
        dot_product = np.dot(norm_a, norm_b)
        angle = np.arccos(np.clip(dot_product, -1, 1)) * 180 / np.pi

        print(f"\n{'='*40}")
        print(f"CHAIN POSITIONING ANALYSIS:")
        print(f"{'='*40}")
        print(f"Dot product of direction vectors: {dot_product:.3f}")
        print(f"Angle between chains (from ligand): {angle:.1f} deg")

        if dot_product < -0.5:
            print(f"[OK] GOOD: Chains on OPPOSITE sides (angle > 120 deg)")
        elif dot_product < 0:
            print(f"[WARN] MARGINAL: Chains at ~90 deg angle")
        else:
            print(f"[FAIL] BAD: Chains on SAME SIDE (angle < 90 deg)")

        # Check contacts
        print(f"\n{'='*40}")
        print(f"CONTACT ANALYSIS:")
        print(f"{'='*40}")
        ca = chain_contacts.get('A', 0)
        cb = chain_contacts.get('B', 0)

        if ca >= 5 and cb >= 5:
            print(f"[OK] GOOD: Both chains have >=5 contacts (A={ca}, B={cb})")
        elif ca >= 5 or cb >= 5:
            print(f"[WARN] PARTIAL: Only one chain has contacts (A={ca}, B={cb})")
        else:
            print(f"[FAIL] BAD: Neither chain has adequate contacts (A={ca}, B={cb})")

    return {
        'chain_contacts': chain_contacts,
        'chain_centroids': chain_centroids,
        'ligand_centroid': ligand_centroid,
    }


def main():
    output_dir = Path("G:/Github_local_repo/Banta_Lab_RFdiffusion/experiments/azobenzene_dimer/outputs/heterodimer_test_20260113_020810")

    # Analyze all designs
    files = [
        "induced_design_00.pdb",
        "induced_design_01.pdb",
        "induced_design_02.pdb",
        "induced_design_03.pdb",
        "induced_design_04.pdb",
        "asymmetric_rasa_design_00.pdb",
        "asymmetric_rasa_design_01.pdb",
        "asymmetric_rasa_design_02.pdb",
        "asymmetric_rasa_design_03.pdb",
        "asymmetric_rasa_design_04.pdb",
        "joint_design_00.pdb",
        "joint_design_01.pdb",
        "joint_design_02.pdb",
        "joint_design_03.pdb",
        "joint_design_04.pdb",
    ]

    for f in files:
        path = output_dir / f
        if path.exists():
            analyze_heterodimer(str(path))
        else:
            print(f"File not found: {path}")


if __name__ == "__main__":
    main()
