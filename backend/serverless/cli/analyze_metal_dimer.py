#!/usr/bin/env python3
"""
Analyze Interface Metal Dimer Design Quality

This script designs interface metal dimers and analyzes:
1. Metal coordination geometry
2. Donor residue distances
3. Secondary structure content
4. Chain interface properties
5. Sequence composition
"""

import requests
import json
import sys
import re
from collections import Counter
from typing import Dict, List, Tuple
import math

BASE_URL = "http://localhost:8000"

# Amino acid categories
HYDROPHOBIC = set(['A', 'V', 'L', 'I', 'M', 'F', 'W', 'P'])
POLAR = set(['S', 'T', 'N', 'Q', 'Y'])
CHARGED_POS = set(['K', 'R', 'H'])
CHARGED_NEG = set(['D', 'E'])
SPECIAL = set(['G', 'C'])

# Metal coordination donors
DONOR_ATOMS = {
    'HIS': ['NE2', 'ND1'],  # Histidine - nitrogen donors
    'CYS': ['SG'],          # Cysteine - sulfur donor
    'ASP': ['OD1', 'OD2'],  # Aspartate - oxygen donors
    'GLU': ['OE1', 'OE2'],  # Glutamate - oxygen donors
}

def parse_pdb(pdb_content: str) -> Dict:
    """Parse PDB content and extract structural information"""
    atoms = []
    chains = set()
    residues = {}

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

                atoms.append({
                    'name': atom_name,
                    'res_name': res_name,
                    'chain': chain,
                    'res_num': res_num,
                    'coords': (x, y, z)
                })
                chains.add(chain)

                key = (chain, res_num)
                if key not in residues:
                    residues[key] = {'name': res_name, 'atoms': []}
                residues[key]['atoms'].append({
                    'name': atom_name,
                    'coords': (x, y, z)
                })
            except (ValueError, IndexError):
                continue

    return {
        'atoms': atoms,
        'chains': sorted(chains),
        'residues': residues,
        'num_residues': len(residues)
    }


def get_sequence(pdb_data: Dict, chain: str) -> str:
    """Extract sequence for a chain"""
    aa_map = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }

    residues = [(k[1], v['name']) for k, v in pdb_data['residues'].items() if k[0] == chain]
    residues.sort(key=lambda x: x[0])

    return ''.join([aa_map.get(r[1], 'X') for r in residues])


def calculate_distance(coord1: Tuple, coord2: Tuple) -> float:
    """Calculate Euclidean distance between two coordinates"""
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(coord1, coord2)))


def find_potential_donors(pdb_data: Dict) -> List[Dict]:
    """Find potential metal donor residues (His, Cys, Asp, Glu)"""
    donors = []
    donor_residues = {'HIS', 'CYS', 'ASP', 'GLU'}

    for (chain, res_num), res_data in pdb_data['residues'].items():
        if res_data['name'] in donor_residues:
            # Find donor atoms
            donor_atoms_for_res = DONOR_ATOMS.get(res_data['name'], [])
            for atom in res_data['atoms']:
                if atom['name'] in donor_atoms_for_res:
                    donors.append({
                        'chain': chain,
                        'res_num': res_num,
                        'res_name': res_data['name'],
                        'atom_name': atom['name'],
                        'coords': atom['coords']
                    })

    return donors


def analyze_coordination(donors: List[Dict], target_metal: str = 'ZN') -> Dict:
    """Analyze potential metal coordination geometry"""
    if len(donors) < 2:
        return {'error': 'Insufficient donors found'}

    # Calculate centroid of donor atoms (potential metal site)
    coords = [d['coords'] for d in donors]
    centroid = tuple(sum(c[i] for c in coords) / len(coords) for i in range(3))

    # Calculate distances from centroid
    distances = [calculate_distance(d['coords'], centroid) for d in donors]

    # Analyze donor composition
    donor_types = Counter([d['res_name'] for d in donors])
    chain_distribution = Counter([d['chain'] for d in donors])

    return {
        'num_potential_donors': len(donors),
        'centroid': centroid,
        'avg_distance_to_centroid': sum(distances) / len(distances),
        'min_distance': min(distances),
        'max_distance': max(distances),
        'donor_types': dict(donor_types),
        'chain_distribution': dict(chain_distribution),
        'donors': donors
    }


def analyze_sequence_composition(sequence: str) -> Dict:
    """Analyze amino acid composition"""
    total = len(sequence)
    if total == 0:
        return {}

    counts = Counter(sequence)

    hydrophobic_count = sum(counts.get(aa, 0) for aa in HYDROPHOBIC)
    polar_count = sum(counts.get(aa, 0) for aa in POLAR)
    pos_count = sum(counts.get(aa, 0) for aa in CHARGED_POS)
    neg_count = sum(counts.get(aa, 0) for aa in CHARGED_NEG)

    return {
        'length': total,
        'hydrophobic_pct': round(100 * hydrophobic_count / total, 1),
        'polar_pct': round(100 * polar_count / total, 1),
        'positive_pct': round(100 * pos_count / total, 1),
        'negative_pct': round(100 * neg_count / total, 1),
        'his_count': counts.get('H', 0),
        'cys_count': counts.get('C', 0),
        'asp_count': counts.get('D', 0),
        'glu_count': counts.get('E', 0),
    }


def estimate_secondary_structure(sequence: str) -> Dict:
    """Estimate secondary structure propensity based on sequence"""
    # Simple propensity-based estimation
    helix_formers = set(['A', 'E', 'L', 'M', 'Q', 'K', 'R'])
    sheet_formers = set(['V', 'I', 'Y', 'F', 'W', 'T'])

    total = len(sequence)
    if total == 0:
        return {}

    helix_count = sum(1 for aa in sequence if aa in helix_formers)
    sheet_count = sum(1 for aa in sequence if aa in sheet_formers)

    return {
        'helix_propensity_pct': round(100 * helix_count / total, 1),
        'sheet_propensity_pct': round(100 * sheet_count / total, 1),
        'coil_propensity_pct': round(100 * (total - helix_count - sheet_count) / total, 1)
    }


def design_and_analyze(metal: str = 'ZN', approach: str = 'joint_metal',
                       num_designs: int = 3) -> Dict:
    """Design interface metal dimers and analyze quality"""

    # Metal-specific parameters
    metal_configs = {
        'ZN': {
            'coordination_split': [2, 2],
            'geometry': 'tetrahedral',
            'chain_a_donors': ['His', 'His'],
            'chain_b_donors': ['Cys', 'Cys'],
        },
        'FE': {
            'coordination_split': [3, 3],
            'geometry': 'octahedral',
            'chain_a_donors': ['His', 'His', 'Cys'],
            'chain_b_donors': ['His', 'Cys', 'Cys'],
        },
        'TB': {
            'coordination_split': [4, 4],
            'geometry': 'square_antiprismatic',
            'chain_a_donors': ['Asp', 'Asp', 'Glu', 'Glu'],
            'chain_b_donors': ['Asp', 'Asp', 'Glu', 'Glu'],
        }
    }

    config = metal_configs.get(metal, metal_configs['ZN'])

    # Build request
    request = {
        'task': 'interface_metal_design',
        'metal': metal,
        'approach': approach,
        'chain_length': '60-80',
        'num_designs': num_designs,
        'num_timesteps': 50,
        'step_scale': 1.0,
        'gamma_0': 0.2,
        'metal_config': {
            **config,
            'include_waters': False,
        }
    }

    print(f"\n{'='*70}")
    print(f"Designing {metal} Interface Metal Dimer ({approach})")
    print(f"{'='*70}")
    print(f"Configuration: {json.dumps(config, indent=2)}")
    print(f"Generating {num_designs} designs...")

    # Call API
    try:
        response = requests.post(
            f"{BASE_URL}/runsync",
            json={"input": request},
            headers={"Content-Type": "application/json"},
            timeout=900
        )
        result = response.json()
    except Exception as e:
        return {'error': str(e)}

    if result.get('status') != 'COMPLETED':
        return {'error': f"Design failed: {result}"}

    output = result['output']['result']
    designs = output.get('designs', [])

    print(f"\nGenerated {len(designs)} designs successfully!")

    # Analyze each design
    analyses = []
    for i, design in enumerate(designs):
        print(f"\n{'-'*50}")
        print(f"Analyzing Design {i+1}")
        print(f"{'-'*50}")

        pdb_content = design.get('pdb_content', '')
        pdb_data = parse_pdb(pdb_content)

        # Basic stats
        print(f"  Chains: {pdb_data['chains']}")
        print(f"  Total residues: {pdb_data['num_residues']}")

        # Per-chain analysis
        chain_analyses = {}
        for chain in pdb_data['chains']:
            seq = get_sequence(pdb_data, chain)
            comp = analyze_sequence_composition(seq)
            ss = estimate_secondary_structure(seq)

            chain_analyses[chain] = {
                'sequence': seq,
                'composition': comp,
                'secondary_structure': ss
            }

            print(f"\n  Chain {chain}:")
            print(f"    Length: {comp.get('length', 0)} residues")
            print(f"    Sequence: {seq[:30]}..." if len(seq) > 30 else f"    Sequence: {seq}")
            print(f"    Hydrophobic: {comp.get('hydrophobic_pct', 0)}%")
            print(f"    Charged (+/-): {comp.get('positive_pct', 0)}% / {comp.get('negative_pct', 0)}%")
            print(f"    His/Cys/Asp/Glu: {comp.get('his_count', 0)}/{comp.get('cys_count', 0)}/{comp.get('asp_count', 0)}/{comp.get('glu_count', 0)}")
            print(f"    Helix propensity: {ss.get('helix_propensity_pct', 0)}%")

        # Metal coordination analysis
        donors = find_potential_donors(pdb_data)
        coord_analysis = analyze_coordination(donors, metal)

        print(f"\n  Metal Coordination Analysis:")
        print(f"    Potential donors found: {coord_analysis.get('num_potential_donors', 0)}")
        print(f"    Donor types: {coord_analysis.get('donor_types', {})}")
        print(f"    Chain distribution: {coord_analysis.get('chain_distribution', {})}")
        if 'avg_distance_to_centroid' in coord_analysis:
            print(f"    Avg distance to centroid: {coord_analysis['avg_distance_to_centroid']:.2f} A")

        analyses.append({
            'design_index': i,
            'pdb_data': pdb_data,
            'chain_analyses': chain_analyses,
            'coordination': coord_analysis,
            'pdb_content': pdb_content
        })

    return {
        'metal': metal,
        'approach': approach,
        'num_designs': len(designs),
        'metal_profile': output.get('metal_profile', {}),
        'analyses': analyses
    }


def main():
    metal = sys.argv[1].upper() if len(sys.argv) > 1 else 'ZN'
    num_designs = int(sys.argv[2]) if len(sys.argv) > 2 else 3

    result = design_and_analyze(metal=metal, num_designs=num_designs)

    if 'error' in result:
        print(f"\n[ERROR] {result['error']}")
        sys.exit(1)

    # Summary
    print(f"\n{'='*70}")
    print("QUALITY SUMMARY")
    print(f"{'='*70}")

    print(f"\nMetal: {result['metal']}")
    print(f"Approach: {result['approach']}")
    print(f"Designs generated: {result['num_designs']}")

    if result.get('metal_profile'):
        mp = result['metal_profile']
        print(f"\nMetal Profile:")
        print(f"  Name: {mp.get('name')}")
        print(f"  Preferred coordination: {mp.get('preferred_coord')}")
        print(f"  Geometry: {mp.get('geometry')}")
        print(f"  Preferred donors: {mp.get('donors')}")

    # Aggregate statistics
    total_donors = 0
    chain_lengths = []
    for analysis in result['analyses']:
        total_donors += analysis['coordination'].get('num_potential_donors', 0)
        for chain, ca in analysis['chain_analyses'].items():
            chain_lengths.append(ca['composition'].get('length', 0))

    avg_donors = total_donors / len(result['analyses']) if result['analyses'] else 0
    avg_length = sum(chain_lengths) / len(chain_lengths) if chain_lengths else 0

    print(f"\nAggregate Statistics:")
    print(f"  Average donors per design: {avg_donors:.1f}")
    print(f"  Average chain length: {avg_length:.1f} residues")

    # Quality assessment
    print(f"\nQuality Assessment:")

    # Check donor count matches expected
    expected_donors = 4  # For Zn tetrahedral
    if result['metal'] == 'FE':
        expected_donors = 6
    elif result['metal'] == 'TB':
        expected_donors = 8

    donor_match = abs(avg_donors - expected_donors) <= 2
    print(f"  [{'PASS' if donor_match else 'WARN'}] Donor count (expected ~{expected_donors}, got {avg_donors:.1f})")

    # Check chain lengths in range
    length_ok = all(50 <= l <= 90 for l in chain_lengths)
    print(f"  [{'PASS' if length_ok else 'WARN'}] Chain lengths in range (60-80 +/- 10)")

    # Check both chains have donors
    chain_dist_ok = all(
        len(a['coordination'].get('chain_distribution', {})) >= 2
        for a in result['analyses']
    )
    print(f"  [{'PASS' if chain_dist_ok else 'WARN'}] Donors distributed across chains")

    print(f"\n[COMPLETE] Analysis finished successfully")

    # Save first design PDB for further analysis
    if result['analyses']:
        pdb_file = f"design_{result['metal'].lower()}_dimer.pdb"
        with open(pdb_file, 'w') as f:
            f.write(result['analyses'][0]['pdb_content'])
        print(f"\nFirst design saved to: {pdb_file}")

    return result


if __name__ == "__main__":
    main()
