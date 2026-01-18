#!/usr/bin/env python3
"""
Detailed coordination analysis of designed metal dimer
"""

import math
from typing import Dict, List, Tuple

def parse_pdb_detailed(pdb_content: str) -> Dict:
    """Parse PDB with detailed atom information"""
    atoms = []
    metal_atom = None

    for line in pdb_content.split('\n'):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21]
                res_num = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                atom = {
                    'name': atom_name,
                    'res_name': res_name,
                    'chain': chain,
                    'res_num': res_num,
                    'coords': (x, y, z),
                    'line': line
                }
                atoms.append(atom)

                # Check if this is metal (chain L or metal residue names)
                if chain == 'L' or res_name in ['ZN', 'FE', 'CU', 'MN', 'CO', 'NI', 'MG', 'CA', 'TB', 'EU', 'GD']:
                    metal_atom = atom

            except (ValueError, IndexError):
                continue

    return {'atoms': atoms, 'metal': metal_atom}


def calculate_distance(coord1: Tuple, coord2: Tuple) -> float:
    """Calculate Euclidean distance"""
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(coord1, coord2)))


def find_coordinating_atoms(pdb_data: Dict, cutoff: float = 3.5) -> List[Dict]:
    """Find atoms within coordination distance of metal"""
    metal = pdb_data['metal']
    if not metal:
        print("No metal atom found!")
        return []

    metal_coords = metal['coords']
    coordinating = []

    # Donor atom types
    donor_atoms = {
        'HIS': ['NE2', 'ND1'],
        'CYS': ['SG'],
        'ASP': ['OD1', 'OD2'],
        'GLU': ['OE1', 'OE2'],
        'ASN': ['OD1', 'ND2'],
        'GLN': ['OE1', 'NE2'],
        'SER': ['OG'],
        'THR': ['OG1'],
        'TYR': ['OH'],
        'MET': ['SD'],
    }

    for atom in pdb_data['atoms']:
        if atom['chain'] == metal['chain']:
            continue  # Skip metal itself

        # Check if this is a potential donor atom
        res_name = atom['res_name']
        atom_name = atom['name']

        is_donor = False
        if res_name in donor_atoms:
            if atom_name in donor_atoms[res_name]:
                is_donor = True

        if is_donor:
            dist = calculate_distance(atom['coords'], metal_coords)
            if dist <= cutoff:
                coordinating.append({
                    'chain': atom['chain'],
                    'res_num': atom['res_num'],
                    'res_name': res_name,
                    'atom_name': atom_name,
                    'distance': dist,
                    'coords': atom['coords']
                })

    return sorted(coordinating, key=lambda x: x['distance'])


def analyze_geometry(coordinating_atoms: List[Dict], metal_coords: Tuple) -> Dict:
    """Analyze coordination geometry"""
    n = len(coordinating_atoms)

    if n < 2:
        return {'geometry': 'insufficient', 'coordination_number': n}

    # Calculate angles between donor atoms through metal center
    angles = []
    for i in range(n):
        for j in range(i + 1, n):
            c1 = coordinating_atoms[i]['coords']
            c2 = coordinating_atoms[j]['coords']

            # Vectors from metal to donors
            v1 = tuple(c1[k] - metal_coords[k] for k in range(3))
            v2 = tuple(c2[k] - metal_coords[k] for k in range(3))

            # Dot product and magnitudes
            dot = sum(v1[k] * v2[k] for k in range(3))
            mag1 = math.sqrt(sum(v1[k]**2 for k in range(3)))
            mag2 = math.sqrt(sum(v2[k]**2 for k in range(3)))

            if mag1 > 0 and mag2 > 0:
                cos_angle = max(-1, min(1, dot / (mag1 * mag2)))
                angle = math.degrees(math.acos(cos_angle))
                angles.append({
                    'donors': (i, j),
                    'angle': angle
                })

    # Determine geometry based on coordination number and angles
    avg_angle = sum(a['angle'] for a in angles) / len(angles) if angles else 0

    geometry = 'unknown'
    ideal_angles = {}

    if n == 4:
        # Tetrahedral: ideal angle ~109.5°
        # Square planar: ideal angles 90° and 180°
        if 100 < avg_angle < 120:
            geometry = 'tetrahedral'
            ideal_angles = {'target': 109.5, 'deviation': abs(avg_angle - 109.5)}
        elif 85 < avg_angle < 95:
            geometry = 'square_planar'
            ideal_angles = {'target': 90.0, 'deviation': abs(avg_angle - 90.0)}
    elif n == 5:
        geometry = 'trigonal_bipyramidal' if avg_angle > 100 else 'square_pyramidal'
    elif n == 6:
        geometry = 'octahedral'
        ideal_angles = {'target': 90.0, 'deviation': abs(avg_angle - 90.0)}
    elif n == 8:
        geometry = 'square_antiprismatic'
    elif n == 9:
        geometry = 'tricapped_trigonal_prismatic'

    return {
        'coordination_number': n,
        'geometry': geometry,
        'average_angle': round(avg_angle, 1),
        'angles': angles,
        'ideal_angles': ideal_angles
    }


def main():
    # Read the designed PDB
    try:
        with open('design_zn_dimer.pdb', 'r') as f:
            pdb_content = f.read()
    except FileNotFoundError:
        print("PDB file not found. Run analyze_metal_dimer.py first.")
        return

    print("="*70)
    print("DETAILED COORDINATION ANALYSIS")
    print("="*70)

    pdb_data = parse_pdb_detailed(pdb_content)

    # Find metal
    metal = pdb_data['metal']
    if metal:
        print(f"\nMetal found:")
        print(f"  Chain: {metal['chain']}")
        print(f"  Residue: {metal['res_name']} {metal['res_num']}")
        print(f"  Coordinates: ({metal['coords'][0]:.2f}, {metal['coords'][1]:.2f}, {metal['coords'][2]:.2f})")
    else:
        print("\nNo metal atom found in structure!")
        return

    # Find coordinating atoms (3.5 Å cutoff for metals)
    coordinating = find_coordinating_atoms(pdb_data, cutoff=3.5)

    print(f"\nCoordinating atoms (within 3.5 Å):")
    print("-"*60)

    if not coordinating:
        print("  No coordinating atoms found within cutoff!")
        # Try larger cutoff
        print("\n  Trying 5.0 Å cutoff...")
        coordinating = find_coordinating_atoms(pdb_data, cutoff=5.0)

    for atom in coordinating:
        print(f"  {atom['chain']}:{atom['res_name']}{atom['res_num']}.{atom['atom_name']} - {atom['distance']:.2f} Å")

    # Analyze geometry
    if coordinating and metal:
        geom = analyze_geometry(coordinating, metal['coords'])

        print(f"\nCoordination Geometry Analysis:")
        print("-"*60)
        print(f"  Coordination number: {geom['coordination_number']}")
        print(f"  Predicted geometry: {geom['geometry']}")
        print(f"  Average donor-metal-donor angle: {geom['average_angle']}°")

        if geom.get('ideal_angles'):
            print(f"  Ideal angle: {geom['ideal_angles']['target']}°")
            print(f"  Deviation from ideal: {geom['ideal_angles']['deviation']:.1f}°")

        # Quality assessment
        print(f"\nQuality Assessment:")
        print("-"*60)

        # Check coordination number
        cn = geom['coordination_number']
        cn_ok = 3 <= cn <= 6  # Reasonable for Zn
        print(f"  [{'PASS' if cn_ok else 'WARN'}] Coordination number ({cn}) - Zn typically 4-6")

        # Check distances
        distances = [a['distance'] for a in coordinating]
        avg_dist = sum(distances) / len(distances) if distances else 0
        dist_ok = 1.8 <= avg_dist <= 2.8  # Typical Zn-donor distances
        print(f"  [{'PASS' if dist_ok else 'WARN'}] Average distance ({avg_dist:.2f} Å) - ideal 2.0-2.5 Å")

        # Check both chains contribute
        chains = set(a['chain'] for a in coordinating)
        interface_ok = len(chains) >= 2
        print(f"  [{'PASS' if interface_ok else 'WARN'}] Interface formation - {len(chains)} chain(s) coordinating")

        # Check donor diversity
        donor_types = set(a['res_name'] for a in coordinating)
        print(f"  [INFO] Donor residue types: {', '.join(donor_types)}")

    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)


if __name__ == "__main__":
    main()
