#!/usr/bin/env python3
"""
Detailed coordination analysis for Terbium (lanthanide) metal dimer
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

                # Check if this is metal (chain L or lanthanide residue names)
                if chain == 'L' or res_name in ['TB', 'EU', 'GD', 'YB', 'DY', 'SM', 'ND', 'LA', 'CE', 'PR']:
                    metal_atom = atom

            except (ValueError, IndexError):
                continue

    return {'atoms': atoms, 'metal': metal_atom}


def calculate_distance(coord1: Tuple, coord2: Tuple) -> float:
    """Calculate Euclidean distance"""
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(coord1, coord2)))


def find_coordinating_atoms(pdb_data: Dict, cutoff: float = 4.0) -> List[Dict]:
    """Find atoms within coordination distance of metal (larger cutoff for lanthanides)"""
    metal = pdb_data['metal']
    if not metal:
        print("No metal atom found!")
        return []

    metal_coords = metal['coords']
    coordinating = []

    # Donor atom types - lanthanides prefer oxygen donors
    donor_atoms = {
        'ASP': ['OD1', 'OD2'],  # Bidentate carboxylate
        'GLU': ['OE1', 'OE2'],  # Bidentate carboxylate
        'ASN': ['OD1', 'ND2'],  # Amide oxygen
        'GLN': ['OE1', 'NE2'],  # Amide oxygen
        'SER': ['OG'],
        'THR': ['OG1'],
        'TYR': ['OH'],
        'HIS': ['NE2', 'ND1'],  # Less common for lanthanides
        'HOH': ['O'],  # Water
        'WAT': ['O'],
    }

    for atom in pdb_data['atoms']:
        if atom['chain'] == metal['chain']:
            continue  # Skip metal itself

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


def analyze_lanthanide_geometry(coordinating_atoms: List[Dict], metal_coords: Tuple) -> Dict:
    """Analyze coordination geometry for lanthanides (CN 8-9)"""
    n = len(coordinating_atoms)

    if n < 2:
        return {'geometry': 'insufficient', 'coordination_number': n}

    # Calculate average distance
    distances = [a['distance'] for a in coordinating_atoms]
    avg_dist = sum(distances) / len(distances)

    # Determine geometry based on coordination number
    geometry = 'unknown'
    ideal_cn = 8  # Target for Tb

    if n == 8:
        geometry = 'square_antiprismatic'
    elif n == 9:
        geometry = 'tricapped_trigonal_prismatic'
    elif n == 7:
        geometry = 'capped_octahedral'
    elif n == 6:
        geometry = 'octahedral'
    elif 4 <= n <= 5:
        geometry = 'distorted_tetrahedral'

    # Check for bidentate carboxylates (both oxygens coordinating)
    bidentate_residues = {}
    for atom in coordinating_atoms:
        key = (atom['chain'], atom['res_num'], atom['res_name'])
        if key not in bidentate_residues:
            bidentate_residues[key] = []
        bidentate_residues[key].append(atom['atom_name'])

    bidentate_count = sum(1 for atoms in bidentate_residues.values() if len(atoms) >= 2)

    return {
        'coordination_number': n,
        'geometry': geometry,
        'average_distance': round(avg_dist, 2),
        'min_distance': round(min(distances), 2),
        'max_distance': round(max(distances), 2),
        'bidentate_carboxylates': bidentate_count,
        'unique_residues': len(bidentate_residues),
        'ideal_cn': ideal_cn
    }


def main():
    # Read the designed PDB
    try:
        with open('design_tb_dimer.pdb', 'r') as f:
            pdb_content = f.read()
    except FileNotFoundError:
        print("PDB file not found. Run analyze_metal_dimer.py TB first.")
        return

    print("="*70)
    print("TERBIUM COORDINATION ANALYSIS")
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

    # Find coordinating atoms (4.0 Å cutoff for lanthanides - longer bonds)
    print(f"\nSearching for coordinating atoms (cutoff: 4.0 Å)...")
    coordinating = find_coordinating_atoms(pdb_data, cutoff=4.0)

    print(f"\nCoordinating atoms found: {len(coordinating)}")
    print("-"*60)

    if not coordinating:
        print("  No coordinating atoms found within cutoff!")
        print("\n  Trying 5.0 Å cutoff...")
        coordinating = find_coordinating_atoms(pdb_data, cutoff=5.0)

    # Group by residue
    by_residue = {}
    for atom in coordinating:
        key = f"{atom['chain']}:{atom['res_name']}{atom['res_num']}"
        if key not in by_residue:
            by_residue[key] = []
        by_residue[key].append(atom)

    for res_key, atoms in sorted(by_residue.items()):
        print(f"\n  {res_key}:")
        for atom in atoms:
            print(f"    {atom['atom_name']}: {atom['distance']:.2f} Å")

    # Analyze geometry
    if coordinating and metal:
        geom = analyze_lanthanide_geometry(coordinating, metal['coords'])

        print(f"\n" + "="*70)
        print("COORDINATION GEOMETRY ANALYSIS")
        print("="*70)
        print(f"\n  Coordination number: {geom['coordination_number']} (ideal for Tb: {geom['ideal_cn']}-9)")
        print(f"  Predicted geometry: {geom['geometry']}")
        print(f"  Average Tb-O distance: {geom['average_distance']} Å (ideal: 2.3-2.5 Å)")
        print(f"  Distance range: {geom['min_distance']} - {geom['max_distance']} Å")
        print(f"  Bidentate carboxylates: {geom['bidentate_carboxylates']}")
        print(f"  Unique coordinating residues: {geom['unique_residues']}")

        # Chain distribution
        chain_dist = {}
        for atom in coordinating:
            chain_dist[atom['chain']] = chain_dist.get(atom['chain'], 0) + 1

        print(f"\n  Chain distribution:")
        for chain, count in sorted(chain_dist.items()):
            print(f"    Chain {chain}: {count} donor atoms")

        # Donor type distribution
        donor_types = {}
        for atom in coordinating:
            donor_types[atom['res_name']] = donor_types.get(atom['res_name'], 0) + 1

        print(f"\n  Donor residue types:")
        for dtype, count in sorted(donor_types.items(), key=lambda x: -x[1]):
            print(f"    {dtype}: {count}")

        # Quality assessment
        print(f"\n" + "="*70)
        print("QUALITY ASSESSMENT")
        print("="*70)

        # Check coordination number
        cn = geom['coordination_number']
        cn_ok = 6 <= cn <= 10  # Reasonable for lanthanides
        cn_ideal = 8 <= cn <= 9
        print(f"\n  [{'PASS' if cn_ideal else ('OK' if cn_ok else 'WARN')}] Coordination number ({cn}) - Tb typically 8-9")

        # Check distances
        dist_ok = 2.0 <= geom['average_distance'] <= 3.5  # Lanthanide-O distances
        dist_ideal = 2.2 <= geom['average_distance'] <= 2.8
        print(f"  [{'PASS' if dist_ideal else ('OK' if dist_ok else 'WARN')}] Average distance ({geom['average_distance']} Å) - ideal 2.3-2.5 Å")

        # Check both chains contribute
        interface_ok = len(chain_dist) >= 2
        print(f"  [{'PASS' if interface_ok else 'WARN'}] Interface formation - {len(chain_dist)} chain(s) coordinating")

        # Check for correct donor types (Asp/Glu for lanthanides)
        correct_donors = set(donor_types.keys()).issubset({'ASP', 'GLU', 'ASN', 'GLN', 'SER', 'THR'})
        print(f"  [{'PASS' if correct_donors else 'OK'}] Donor types - {', '.join(donor_types.keys())}")

        # Check for bidentate coordination (typical for lanthanides)
        bidentate_ok = geom['bidentate_carboxylates'] >= 2
        print(f"  [{'PASS' if bidentate_ok else 'OK'}] Bidentate carboxylates - {geom['bidentate_carboxylates']} found")

        # Lanthanide-specific features
        print(f"\n  Lanthanide-specific assessment:")
        print(f"    - Oxygen-rich coordination sphere: {'YES' if 'ASP' in donor_types or 'GLU' in donor_types else 'NO'}")
        print(f"    - High coordination number (8-9): {'YES' if cn >= 8 else 'NO'}")
        print(f"    - Bidentate carboxylate binding: {'YES' if geom['bidentate_carboxylates'] > 0 else 'NO'}")

    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)


if __name__ == "__main__":
    main()
