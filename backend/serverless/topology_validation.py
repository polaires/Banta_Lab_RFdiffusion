"""
Topology Validation for Dimer Designs

This module provides functions to validate that designed dimers have
proper topology - i.e., the two chains can physically separate without
passing through each other.

Key Functions:
- compute_convex_hull_overlap: Check if chains are entangled
- translation_separation_test: Test if chains can move apart
- validate_dimer_topology: Main validation function
"""

import numpy as np
from typing import Tuple, Optional, Dict, List, Any
from pathlib import Path

try:
    from scipy.spatial import ConvexHull, Delaunay
    from scipy.spatial.distance import cdist
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    print("Warning: scipy not available. Topology validation will use fallback methods.")

try:
    from Bio.PDB import PDBParser
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    print("Warning: BioPython not available. Will use simple PDB parsing.")


def parse_pdb_simple(pdb_path: str) -> Dict[str, np.ndarray]:
    """
    Simple PDB parser that extracts coordinates by chain.
    Fallback when BioPython is not available.

    Args:
        pdb_path: Path to PDB file

    Returns:
        Dictionary mapping chain IDs to coordinate arrays
    """
    chains = {}
    ligand_coords = []

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain_id = line[21]
                res_name = line[17:20].strip()
                element = line[76:78].strip() if len(line) > 76 else line[12:14].strip()

                # Skip hydrogens
                if element == 'H':
                    continue

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                coord = [x, y, z]

                # Check if ligand (common ligand residue names)
                if res_name in ['AZO', 'UNL', 'LIG', 'DRG'] or line.startswith('HETATM'):
                    ligand_coords.append(coord)
                else:
                    if chain_id not in chains:
                        chains[chain_id] = []
                    chains[chain_id].append(coord)

    # Convert to numpy arrays
    result = {}
    for chain_id, coords in chains.items():
        if coords:
            result[chain_id] = np.array(coords)

    if ligand_coords:
        result['LIGAND'] = np.array(ligand_coords)

    return result


def parse_pdb_biopython(pdb_path: str) -> Dict[str, np.ndarray]:
    """
    Parse PDB using BioPython for more robust handling.

    Args:
        pdb_path: Path to PDB file

    Returns:
        Dictionary mapping chain IDs to coordinate arrays
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_path)

    chains = {}
    ligand_coords = []

    for model in structure:
        for chain in model:
            chain_id = chain.id
            coords = []

            for residue in chain:
                res_name = residue.resname

                # Check if ligand
                is_ligand = res_name in ['AZO', 'UNL', 'LIG', 'DRG'] or residue.id[0] != ' '

                for atom in residue:
                    # Skip hydrogens
                    if atom.element == 'H':
                        continue

                    if is_ligand:
                        ligand_coords.append(atom.coord)
                    else:
                        coords.append(atom.coord)

            if coords:
                chains[chain_id] = np.array(coords)

    if ligand_coords:
        chains['LIGAND'] = np.array(ligand_coords)

    return chains


def parse_pdb(pdb_path: str) -> Dict[str, np.ndarray]:
    """
    Parse PDB file and extract coordinates by chain.
    Uses BioPython if available, otherwise falls back to simple parser.

    Args:
        pdb_path: Path to PDB file

    Returns:
        Dictionary mapping chain IDs to coordinate arrays
    """
    if BIOPYTHON_AVAILABLE:
        return parse_pdb_biopython(pdb_path)
    else:
        return parse_pdb_simple(pdb_path)


def compute_convex_hull_overlap(
    chain_a_coords: np.ndarray,
    chain_b_coords: np.ndarray
) -> float:
    """
    Compute the overlap between convex hulls of two chains.

    High overlap (>0.3) indicates the chains are entangled and
    cannot physically separate.

    Args:
        chain_a_coords: Nx3 array of chain A atom coordinates
        chain_b_coords: Mx3 array of chain B atom coordinates

    Returns:
        Overlap fraction (0-1), where <0.30 is acceptable
    """
    if not SCIPY_AVAILABLE:
        # Fallback: use bounding box overlap
        return compute_bounding_box_overlap(chain_a_coords, chain_b_coords)

    try:
        # Compute convex hulls
        hull_a = ConvexHull(chain_a_coords)
        hull_b = ConvexHull(chain_b_coords)

        # Create Delaunay triangulations for point-in-hull tests
        del_a = Delaunay(chain_a_coords[hull_a.vertices])
        del_b = Delaunay(chain_b_coords[hull_b.vertices])

        # Count points of B inside hull A
        inside_a = np.sum(del_a.find_simplex(chain_b_coords) >= 0)

        # Count points of A inside hull B
        inside_b = np.sum(del_b.find_simplex(chain_a_coords) >= 0)

        # Calculate overlap fraction
        total_points = len(chain_a_coords) + len(chain_b_coords)
        overlap_fraction = (inside_a + inside_b) / total_points

        return overlap_fraction

    except Exception as e:
        print(f"Warning: ConvexHull computation failed: {e}")
        return compute_bounding_box_overlap(chain_a_coords, chain_b_coords)


def compute_bounding_box_overlap(
    chain_a_coords: np.ndarray,
    chain_b_coords: np.ndarray
) -> float:
    """
    Fallback method: compute overlap using bounding boxes.

    Args:
        chain_a_coords: Nx3 array of chain A atom coordinates
        chain_b_coords: Mx3 array of chain B atom coordinates

    Returns:
        Overlap fraction estimate
    """
    # Get bounding boxes
    min_a = np.min(chain_a_coords, axis=0)
    max_a = np.max(chain_a_coords, axis=0)
    min_b = np.min(chain_b_coords, axis=0)
    max_b = np.max(chain_b_coords, axis=0)

    # Calculate intersection
    intersection_min = np.maximum(min_a, min_b)
    intersection_max = np.minimum(max_a, max_b)

    # Check if there's any intersection
    if np.any(intersection_min >= intersection_max):
        return 0.0

    # Calculate volumes
    intersection_vol = np.prod(intersection_max - intersection_min)
    vol_a = np.prod(max_a - min_a)
    vol_b = np.prod(max_b - min_b)

    # Overlap as fraction of smaller volume
    min_vol = min(vol_a, vol_b)
    if min_vol == 0:
        return 0.0

    return intersection_vol / min_vol


def translation_separation_test(
    chain_a_coords: np.ndarray,
    chain_b_coords: np.ndarray,
    ligand_coords: Optional[np.ndarray] = None,
    directions: Optional[List[Tuple[float, float, float]]] = None,
    max_distance: float = 50.0,
    step_size: float = 2.0,
    clash_threshold: float = 3.0
) -> Tuple[bool, Optional[np.ndarray]]:
    """
    Test if chains can be separated by translation without clashing.

    This simulates pulling chain B away from chain A along various
    directions to see if it can escape without passing through chain A.

    Args:
        chain_a_coords: Nx3 array of chain A atom coordinates
        chain_b_coords: Mx3 array of chain B atom coordinates
        ligand_coords: Optional ligand coordinates (chain B moves away from ligand too)
        directions: List of direction vectors to test
        max_distance: Maximum translation distance to test
        step_size: Step size for translation
        clash_threshold: Minimum distance to avoid clashes

    Returns:
        Tuple of (can_separate, separation_direction)
    """
    if directions is None:
        # Test multiple directions
        directions = [
            (1, 0, 0), (-1, 0, 0),
            (0, 1, 0), (0, -1, 0),
            (0, 0, 1), (0, 0, -1),
            (1, 1, 0), (1, -1, 0), (-1, 1, 0), (-1, -1, 0),
            (1, 0, 1), (1, 0, -1), (-1, 0, 1), (-1, 0, -1),
            (0, 1, 1), (0, 1, -1), (0, -1, 1), (0, -1, -1),
            (1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1),
            (-1, -1, 1), (-1, 1, -1), (1, -1, -1), (-1, -1, -1)
        ]

    for direction in directions:
        direction = np.array(direction, dtype=float)
        direction = direction / np.linalg.norm(direction)

        can_separate = True

        # Translate chain B along direction
        for distance in np.arange(step_size, max_distance + step_size, step_size):
            translated_b = chain_b_coords + direction * distance

            # Check minimum distance to chain A
            if SCIPY_AVAILABLE:
                distances = cdist(chain_a_coords, translated_b)
                min_dist = np.min(distances)
            else:
                # Fallback: sample distances
                min_dist = float('inf')
                for a_coord in chain_a_coords[::10]:  # Sample every 10th atom
                    for b_coord in translated_b[::10]:
                        dist = np.linalg.norm(a_coord - b_coord)
                        min_dist = min(min_dist, dist)

            # If we have good clearance at this distance
            if min_dist > clash_threshold:
                # Continue checking if clearance increases
                continue
            else:
                # Still clashing, break inner loop
                can_separate = False
                break

        # If we made it through without permanent clashes, chains can separate
        if can_separate:
            return True, direction

    return False, None


def count_contacts(
    coords_a: np.ndarray,
    coords_b: np.ndarray,
    distance_threshold: float = 4.0
) -> int:
    """
    Count number of contacts between two sets of coordinates.

    Args:
        coords_a: Nx3 array
        coords_b: Mx3 array
        distance_threshold: Maximum distance for a contact

    Returns:
        Number of atom pairs within threshold
    """
    if SCIPY_AVAILABLE:
        distances = cdist(coords_a, coords_b)
        return np.sum(distances < distance_threshold)
    else:
        count = 0
        for a_coord in coords_a:
            for b_coord in coords_b:
                if np.linalg.norm(a_coord - b_coord) < distance_threshold:
                    count += 1
        return count


def validate_dimer_topology(
    pdb_path: str,
    ligand_resname: str = 'AZO',
    hull_overlap_threshold: float = 0.30,
    contact_threshold: int = 3,
    contact_distance: float = 4.0
) -> Dict[str, Any]:
    """
    Main validation function for dimer topology.

    Checks:
    1. Convex hull overlap (< threshold means separable)
    2. Translation separation test (can chains physically separate?)
    3. Contact counts per chain (both should contact ligand)

    Args:
        pdb_path: Path to PDB file
        ligand_resname: Residue name of ligand
        hull_overlap_threshold: Maximum acceptable hull overlap
        contact_threshold: Minimum contacts per chain
        contact_distance: Distance threshold for contacts

    Returns:
        Dictionary with validation results
    """
    result = {
        'pdb_path': str(pdb_path),
        'valid': False,
        'error': None
    }

    # Parse PDB
    try:
        chains = parse_pdb(pdb_path)
    except Exception as e:
        result['error'] = f"Failed to parse PDB: {e}"
        return result

    # Get chain coordinates (excluding ligand)
    chain_ids = [k for k in chains.keys() if k != 'LIGAND']

    if len(chain_ids) < 2:
        result['error'] = f"Less than 2 chains found (found: {chain_ids})"
        result['num_chains'] = len(chain_ids)
        return result

    # Use first two chains
    chain_a_coords = chains[chain_ids[0]]
    chain_b_coords = chains[chain_ids[1]]
    ligand_coords = chains.get('LIGAND', None)

    result['chain_a_id'] = chain_ids[0]
    result['chain_b_id'] = chain_ids[1]
    result['chain_a_atoms'] = len(chain_a_coords)
    result['chain_b_atoms'] = len(chain_b_coords)

    # 1. Convex hull overlap test
    hull_overlap = compute_convex_hull_overlap(chain_a_coords, chain_b_coords)
    result['hull_overlap'] = hull_overlap
    result['hull_overlap_pass'] = hull_overlap < hull_overlap_threshold

    # 2. Translation separation test
    separable, sep_direction = translation_separation_test(
        chain_a_coords, chain_b_coords, ligand_coords
    )
    result['separable'] = separable
    result['separation_direction'] = sep_direction.tolist() if sep_direction is not None else None

    # 3. Contact analysis (if ligand present)
    if ligand_coords is not None:
        contacts_a = count_contacts(chain_a_coords, ligand_coords, contact_distance)
        contacts_b = count_contacts(chain_b_coords, ligand_coords, contact_distance)

        result['contacts_chain_a'] = contacts_a
        result['contacts_chain_b'] = contacts_b
        result['contacts_chain_a_pass'] = contacts_a >= contact_threshold
        result['contacts_chain_b_pass'] = contacts_b >= contact_threshold
        result['both_chains_contact'] = contacts_a >= contact_threshold and contacts_b >= contact_threshold
        result['ligand_atoms'] = len(ligand_coords)
    else:
        result['contacts_chain_a'] = None
        result['contacts_chain_b'] = None
        result['both_chains_contact'] = None
        result['ligand_atoms'] = 0

    # 4. Overall topology validity
    result['topology_valid'] = result['hull_overlap_pass'] and result['separable']

    # 5. Overall validity (topology + contacts)
    if ligand_coords is not None:
        result['valid'] = (
            result['topology_valid'] and
            result.get('both_chains_contact', False)
        )
    else:
        result['valid'] = result['topology_valid']

    return result


def validate_batch(
    pdb_paths: List[str],
    **kwargs
) -> List[Dict[str, Any]]:
    """
    Validate multiple PDB files.

    Args:
        pdb_paths: List of paths to PDB files
        **kwargs: Arguments passed to validate_dimer_topology

    Returns:
        List of validation results
    """
    results = []
    for pdb_path in pdb_paths:
        result = validate_dimer_topology(pdb_path, **kwargs)
        results.append(result)
    return results


def filter_valid_designs(
    pdb_paths: List[str],
    require_contacts: bool = True,
    **kwargs
) -> Tuple[List[str], List[Dict[str, Any]]]:
    """
    Filter PDB files to only those with valid topology.

    Args:
        pdb_paths: List of paths to PDB files
        require_contacts: Whether to require ligand contacts
        **kwargs: Arguments passed to validate_dimer_topology

    Returns:
        Tuple of (valid_paths, all_results)
    """
    results = validate_batch(pdb_paths, **kwargs)

    valid_paths = []
    for pdb_path, result in zip(pdb_paths, results):
        if require_contacts:
            if result.get('valid', False):
                valid_paths.append(pdb_path)
        else:
            if result.get('topology_valid', False):
                valid_paths.append(pdb_path)

    return valid_paths, results


# CLI interface
if __name__ == '__main__':
    import sys
    import json

    if len(sys.argv) < 2:
        print("Usage: python topology_validation.py <pdb_file> [ligand_resname]")
        print("       python topology_validation.py --batch <pdb_file1> <pdb_file2> ...")
        sys.exit(1)

    if sys.argv[1] == '--batch':
        pdb_paths = sys.argv[2:]
        results = validate_batch(pdb_paths)

        valid_count = sum(1 for r in results if r.get('valid', False))
        print(f"\nResults: {valid_count}/{len(results)} designs have valid topology")

        for r in results:
            status = "VALID" if r.get('valid', False) else "INVALID"
            print(f"  {Path(r['pdb_path']).name}: {status}")
            if not r.get('valid', False):
                if r.get('error'):
                    print(f"    Error: {r['error']}")
                elif not r.get('hull_overlap_pass', True):
                    print(f"    Hull overlap: {r.get('hull_overlap', 'N/A'):.3f} (max 0.30)")
                elif not r.get('separable', True):
                    print(f"    Chains cannot separate")
                elif not r.get('both_chains_contact', True):
                    print(f"    Contacts: A={r.get('contacts_chain_a')}, B={r.get('contacts_chain_b')}")
    else:
        pdb_path = sys.argv[1]
        ligand_resname = sys.argv[2] if len(sys.argv) > 2 else 'AZO'

        result = validate_dimer_topology(pdb_path, ligand_resname=ligand_resname)

        print(json.dumps(result, indent=2, default=str))

        if result.get('valid', False):
            print("\n✓ VALID: Design has correct topology")
        else:
            print("\n✗ INVALID: Design has topology issues")
            if result.get('error'):
                print(f"  Error: {result['error']}")
