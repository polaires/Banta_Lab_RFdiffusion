"""
Hotspot Detection Module

Auto-detect optimal binding hotspots for protein binder design.
Uses SASA analysis combined with spatial clustering to identify
surface patches suitable for binder attachment.

Based on research from:
- BindCraft: Recommends clustered hydrophobic patches (F,Y,W,I,L,M)
- ECMIS: 7 Å spatial clustering, combined scoring
- PPI-hotspotID: SASA + hydrophobicity + conservation
- FTMap: Probe-based consensus site detection

Methods:
1. SASA-based: Find exposed residues suitable for binding
2. Clustering: Group nearby surface residues into binding patches
3. Property-based: Prioritize hydrophobic/aromatic residues at surface
"""

import sys
import os
from typing import List, Dict, Any, Optional, Tuple
import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.sasa import calculate_sasa

# BindCraft-recommended hotspot residues (best for binder design)
# "Ideal target patches contain either a phenylalanine, tyrosine, tryptophan,
# isoleucine, leucine, or methionine at the binding site"
BINDCRAFT_PREFERRED = {"PHE", "TYR", "TRP", "ILE", "LEU", "MET"}

# Amino acid properties for scoring
HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}
AROMATIC = {"PHE", "TRP", "TYR", "HIS"}
POLAR = {"SER", "THR", "ASN", "GLN", "TYR", "CYS"}
CHARGED_POSITIVE = {"ARG", "LYS", "HIS"}
CHARGED_NEGATIVE = {"ASP", "GLU"}
CHARGED = CHARGED_POSITIVE | CHARGED_NEGATIVE

# BindCraft notes that lysine is "poorly targetable"
POORLY_TARGETABLE = {"LYS"}

# Spatial clustering parameters (from ECMIS paper)
CLUSTER_DISTANCE_THRESHOLD = 7.0  # Angstroms - residues within this form a patch
PATCH_RADIUS = 12.0  # Angstroms - maximum patch size


def detect_hotspots_sasa(
    pdb_content: str,
    target_chain: str = "A",
    method: str = "exposed_clustered",
    max_hotspots: int = 5,
    prefer_hydrophobic: bool = True,
    min_relative_sasa: float = 0.25,
) -> Dict[str, Any]:
    """
    Detect binding hotspots using SASA analysis.

    Methods:
        - "exposed": All exposed residues (relative_sasa > min_relative_sasa)
        - "exposed_clustered": Exposed residues clustered spatially
        - "patch": Largest contiguous surface patch
        - "interface_like": Residues with mixed polar/hydrophobic character

    Args:
        pdb_content: PDB file content as string
        target_chain: Chain to analyze (default: "A")
        method: Detection method (default: "exposed_clustered")
        max_hotspots: Maximum hotspots to return (default: 5)
        prefer_hydrophobic: Prioritize hydrophobic residues (default: True)
        min_relative_sasa: Minimum relative SASA for exposed (default: 0.25)

    Returns:
        Dict with:
            - hotspots: List of hotspot residue IDs ["A25", "A30", ...]
            - method: Detection method used
            - residue_details: Details for each hotspot
            - total_exposed: Total exposed residues found
    """
    # Step 1: Calculate per-residue SASA
    sasa_result = calculate_sasa(pdb_content)

    if not sasa_result.get("success"):
        return {
            "status": "error",
            "error": f"SASA calculation failed: {sasa_result.get('error', 'Unknown error')}",
            "hotspots": [],
        }

    # Step 2: Extract CA coordinates for clustering
    ca_coords = _extract_ca_coords(pdb_content, target_chain)

    # Step 3: Filter to target chain and exposed residues
    exposed_residues = []
    per_residue = sasa_result.get("per_residue", {})

    for res_key, data in per_residue.items():
        # Parse residue key format "A:ALA15"
        if ":" not in res_key:
            continue

        chain_part, res_part = res_key.split(":", 1)
        if chain_part != target_chain:
            continue

        # Check if exposed
        relative_sasa = data.get("relative_sasa", 0)
        burial = data.get("burial_classification", "buried")

        if burial == "exposed" or relative_sasa >= min_relative_sasa:
            # Extract residue number
            res_name = res_part[:3] if len(res_part) >= 3 else res_part
            res_num_str = res_part[3:] if len(res_part) > 3 else ""

            try:
                res_num = int(res_num_str)
            except ValueError:
                continue

            # Get CA coordinates if available
            coords = ca_coords.get(res_num)

            exposed_residues.append({
                "residue_id": f"{target_chain}{res_num}",
                "residue_name": res_name,
                "residue_number": res_num,
                "relative_sasa": relative_sasa,
                "coords": coords,
                "property": _classify_residue_property(res_name),
            })

    if not exposed_residues:
        return {
            "status": "warning",
            "hotspots": [],
            "method": method,
            "total_exposed": 0,
            "message": f"No exposed residues found in chain {target_chain}",
        }

    # Step 4: Apply method-specific selection
    if method == "exposed":
        selected = exposed_residues[:max_hotspots]

    elif method == "exposed_clustered":
        selected = _cluster_residues(exposed_residues, max_hotspots)

    elif method == "patch":
        selected = _find_surface_patch(exposed_residues, max_hotspots)

    elif method == "interface_like":
        selected = _select_interface_like(exposed_residues, max_hotspots)

    else:
        # Default to exposed_clustered
        selected = _cluster_residues(exposed_residues, max_hotspots)

    # Step 5: Optionally prioritize hydrophobic residues
    if prefer_hydrophobic and len(selected) > max_hotspots:
        selected = _prioritize_hydrophobic(selected, max_hotspots)

    # Step 6: Ensure we have the right number
    selected = selected[:max_hotspots]

    # Format output
    hotspot_list = [r["residue_id"] for r in selected]

    # Calculate cluster center if we have coordinates
    cluster_center = None
    coords_list = [r["coords"] for r in selected if r.get("coords") is not None]
    if coords_list:
        coords_array = np.array(coords_list)
        cluster_center = {
            "x": float(coords_array[:, 0].mean()),
            "y": float(coords_array[:, 1].mean()),
            "z": float(coords_array[:, 2].mean()),
        }

    return {
        "status": "completed",
        "hotspots": hotspot_list,
        "method": method,
        "cluster_center": cluster_center,
        "residue_details": [
            {
                "residue": r["residue_id"],
                "restype": r["residue_name"],  # Match frontend interface
                "relative_sasa": round(r["relative_sasa"], 3),  # Match frontend interface
                "property": r["property"],
            }
            for r in selected
        ],
        "total_exposed": len(exposed_residues),
    }


def _extract_ca_coords(pdb_content: str, chain: str) -> Dict[int, Tuple[float, float, float]]:
    """
    Extract CA atom coordinates for each residue in a chain.

    Args:
        pdb_content: PDB file content
        chain: Chain identifier

    Returns:
        Dict mapping residue_number -> (x, y, z) coordinates
    """
    ca_coords = {}

    for line in pdb_content.split('\n'):
        if not line.startswith('ATOM'):
            continue

        # Parse atom name (columns 13-16)
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue

        # Parse chain (column 22)
        atom_chain = line[21].strip()
        if atom_chain != chain:
            continue

        # Parse residue number (columns 23-26)
        try:
            res_num = int(line[22:26].strip())
        except ValueError:
            continue

        # Parse coordinates (columns 31-54)
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            ca_coords[res_num] = (x, y, z)
        except ValueError:
            continue

    return ca_coords


def _classify_residue_property(res_name: str) -> str:
    """Classify residue by chemical property."""
    res = res_name.upper()[:3]

    if res in HYDROPHOBIC:
        return "hydrophobic"
    elif res in AROMATIC:
        return "aromatic"
    elif res in POLAR:
        return "polar"
    elif res in CHARGED_POSITIVE:
        return "positive"
    elif res in CHARGED_NEGATIVE:
        return "negative"
    else:
        return "other"


def _calculate_hotspot_score(residue: Dict, all_residues: List[Dict]) -> float:
    """
    Calculate comprehensive hotspot score based on research.

    Scoring based on:
    - BindCraft: F,Y,W,I,L,M are ideal (+3 points)
    - ECMIS: Spatial clustering within 7 Å enhances score
    - PPI-hotspotID: SASA contributes to score
    - BindCraft: K (lysine) is poorly targetable (-2 points)

    Args:
        residue: Residue dict with coords, residue_name, relative_sasa
        all_residues: All candidate residues for calculating clustering bonus

    Returns:
        Combined hotspot score (higher = better)
    """
    restype = residue.get("residue_name", "")[:3].upper()
    sasa = residue.get("relative_sasa", 0)
    coords = residue.get("coords")

    # Base score from SASA (0-1 normalized)
    score = sasa * 1.0

    # BindCraft-preferred residues get strong bonus
    if restype in BINDCRAFT_PREFERRED:
        score += 3.0
    elif restype in AROMATIC:
        score += 2.5  # Aromatics good for binding
    elif restype in HYDROPHOBIC:
        score += 2.0  # Other hydrophobics also good
    elif restype in POLAR:
        score += 0.5  # Polar OK for specificity
    elif restype in CHARGED_NEGATIVE:
        score -= 0.5  # Charged less ideal
    elif restype in CHARGED_POSITIVE:
        if restype in POORLY_TARGETABLE:  # LYS
            score -= 2.0  # BindCraft: lysine poorly targetable
        else:
            score -= 0.5

    # ECMIS-inspired: Spatial clustering bonus (7 Å proximity)
    if coords is not None:
        coords_arr = np.array(coords)
        nearby_count = 0
        nearby_preferred = 0

        for other in all_residues:
            if other is residue:
                continue
            other_coords = other.get("coords")
            if other_coords is None:
                continue

            distance = np.linalg.norm(coords_arr - np.array(other_coords))
            if distance < CLUSTER_DISTANCE_THRESHOLD:  # 7 Å
                nearby_count += 1
                other_restype = other.get("residue_name", "")[:3].upper()
                if other_restype in BINDCRAFT_PREFERRED or other_restype in HYDROPHOBIC:
                    nearby_preferred += 1

        # Bonus for having nearby residues (forms a patch)
        score += nearby_count * 0.3
        # Extra bonus for nearby hydrophobic/preferred residues
        score += nearby_preferred * 0.2

    return score


def _cluster_residues(residues: List[Dict], max_clusters: int) -> List[Dict]:
    """
    Find the best compact patch of hotspot residues.

    Based on BindCraft recommendation: "several surface residues within a radial patch"
    and ECMIS: residues within 7 Å form efficient binding patches.

    Strategy:
    1. Score all residues using comprehensive hotspot scoring
    2. Find the best "seed" residue (highest score)
    3. Expand patch by adding nearby high-scoring residues
    4. Ensure all selected residues are within PATCH_RADIUS (12 Å) of each other

    Args:
        residues: List of residue dicts with 'coords' field
        max_clusters: Maximum residues to return

    Returns:
        List of residues forming a compact binding patch
    """
    if len(residues) <= max_clusters:
        return residues

    # Filter to residues with coordinates
    with_coords = [r for r in residues if r.get("coords") is not None]

    if len(with_coords) < 2:
        # Not enough coordinates, fall back to score-based selection
        for r in residues:
            r["hotspot_score"] = _calculate_hotspot_score(r, residues)
        return sorted(residues, key=lambda x: x.get("hotspot_score", 0), reverse=True)[:max_clusters]

    # Calculate hotspot score for each residue
    for r in with_coords:
        r["hotspot_score"] = _calculate_hotspot_score(r, with_coords)

    # Sort by hotspot score (best first)
    sorted_residues = sorted(with_coords, key=lambda x: x["hotspot_score"], reverse=True)

    # Find best compact patch using greedy expansion from best seed
    best_patch = _find_best_patch(sorted_residues, max_clusters, PATCH_RADIUS)

    return best_patch


def _find_best_patch(
    sorted_residues: List[Dict],
    max_size: int,
    max_radius: float = 12.0
) -> List[Dict]:
    """
    Find the best compact patch starting from the highest-scoring residue.

    Ensures all residues in the patch are within max_radius of each other,
    forming a focused binding site rather than scattered hotspots.

    Args:
        sorted_residues: Residues sorted by hotspot_score (descending)
        max_size: Maximum patch size
        max_radius: Maximum distance between any two residues in patch (Angstroms)

    Returns:
        List of residues forming the best compact patch
    """
    if len(sorted_residues) <= max_size:
        return sorted_residues

    best_patch = []
    best_score = 0

    # Try starting from each of the top candidates as seed
    num_seeds_to_try = min(5, len(sorted_residues))

    for seed_idx in range(num_seeds_to_try):
        seed = sorted_residues[seed_idx]
        seed_coords = np.array(seed["coords"])

        # Start patch with seed
        patch = [seed]
        patch_coords = [seed_coords]

        # Add nearby high-scoring residues
        for candidate in sorted_residues:
            if candidate in patch:
                continue
            if len(patch) >= max_size:
                break

            cand_coords = np.array(candidate["coords"])

            # Check if candidate is within max_radius of ALL residues in patch
            is_within_radius = True
            for existing_coords in patch_coords:
                distance = np.linalg.norm(cand_coords - existing_coords)
                if distance > max_radius:
                    is_within_radius = False
                    break

            if is_within_radius:
                patch.append(candidate)
                patch_coords.append(cand_coords)

        # Calculate total patch score
        patch_score = sum(r.get("hotspot_score", 0) for r in patch)

        # Bonus for having preferred residues in patch
        preferred_count = sum(
            1 for r in patch
            if r.get("residue_name", "")[:3].upper() in BINDCRAFT_PREFERRED
        )
        patch_score += preferred_count * 0.5

        if patch_score > best_score:
            best_score = patch_score
            best_patch = patch

    return best_patch


def _simple_cluster(residues: List[Dict], max_clusters: int, min_distance: float = 10.0) -> List[Dict]:
    """
    Simple greedy clustering without scipy.

    Selects residues that are at least min_distance apart, preferring
    those with higher SASA values.

    Args:
        residues: List of residue dicts with 'coords' field
        max_clusters: Maximum clusters to return
        min_distance: Minimum distance between selected residues (Angstroms)

    Returns:
        List of selected residues
    """
    # Sort by SASA (most exposed first)
    sorted_residues = sorted(residues, key=lambda x: x["relative_sasa"], reverse=True)

    selected = []
    for residue in sorted_residues:
        if len(selected) >= max_clusters:
            break

        coords = np.array(residue["coords"])

        # Check if far enough from all already-selected residues
        is_far_enough = True
        for existing in selected:
            existing_coords = np.array(existing["coords"])
            distance = np.linalg.norm(coords - existing_coords)
            if distance < min_distance:
                is_far_enough = False
                break

        if is_far_enough:
            selected.append(residue)

    return selected


def _find_surface_patch(residues: List[Dict], max_count: int) -> List[Dict]:
    """
    Find the largest contiguous surface patch.

    Uses a simplified approach: find the residue with most neighbors
    and select nearby residues.

    Args:
        residues: List of residue dicts with 'coords' field
        max_count: Maximum residues to return

    Returns:
        List of residues forming the largest patch
    """
    if len(residues) <= max_count:
        return residues

    with_coords = [r for r in residues if r.get("coords") is not None]

    if len(with_coords) < 2:
        return residues[:max_count]

    # For each residue, count neighbors within 12 Angstroms
    neighbor_distance = 12.0
    neighbor_counts = []

    for i, res in enumerate(with_coords):
        coords_i = np.array(res["coords"])
        count = 0
        for j, other in enumerate(with_coords):
            if i != j:
                coords_j = np.array(other["coords"])
                distance = np.linalg.norm(coords_i - coords_j)
                if distance < neighbor_distance:
                    count += 1
        neighbor_counts.append(count)

    # Find residue with most neighbors (center of largest patch)
    max_neighbors_idx = np.argmax(neighbor_counts)
    center_residue = with_coords[max_neighbors_idx]
    center_coords = np.array(center_residue["coords"])

    # Select residues closest to the center
    distances_to_center = []
    for res in with_coords:
        coords = np.array(res["coords"])
        dist = np.linalg.norm(coords - center_coords)
        distances_to_center.append((dist, res))

    # Sort by distance and return closest ones
    distances_to_center.sort(key=lambda x: x[0])
    selected = [r for _, r in distances_to_center[:max_count]]

    return selected


def _select_interface_like(residues: List[Dict], max_count: int) -> List[Dict]:
    """
    Select residues with interface-like properties.

    Interfaces typically have a mix of hydrophobic and polar residues,
    with some aromatic and charged residues for specificity.

    Args:
        residues: List of residue dicts
        max_count: Maximum residues to return

    Returns:
        List of selected residues
    """
    # Score each residue for interface-like character
    scored = []

    for r in residues:
        restype = r.get("residue_name", "")[:3].upper()
        sasa = r.get("relative_sasa", 0)

        # Base score from SASA (more exposed = higher base score)
        score = sasa * 2

        # Property-based bonuses
        if restype in AROMATIC:
            score += 3  # Aromatics often at interfaces (pi-stacking, cation-pi)
        elif restype in HYDROPHOBIC:
            score += 2  # Hydrophobic core of interfaces
        elif restype in POLAR:
            score += 1.5  # Polar residues for specificity
        elif restype in CHARGED:
            score += 1  # Charged for salt bridges

        # Prefer medium SASA (0.3-0.6) over fully exposed
        if 0.3 <= sasa <= 0.6:
            score += 1

        r["interface_score"] = score
        scored.append(r)

    # Sort by interface score
    scored.sort(key=lambda x: x.get("interface_score", 0), reverse=True)

    return scored[:max_count]


def _prioritize_hydrophobic(residues: List[Dict], max_count: int) -> List[Dict]:
    """
    Prioritize hydrophobic residues among the candidates.

    Args:
        residues: List of residue dicts
        max_count: Maximum residues to return

    Returns:
        List of selected residues, prioritizing hydrophobic
    """
    # Separate hydrophobic and non-hydrophobic
    hydrophobic = []
    other = []

    for r in residues:
        restype = r.get("residue_name", "")[:3].upper()
        if restype in HYDROPHOBIC or restype in AROMATIC:
            hydrophobic.append(r)
        else:
            other.append(r)

    # Sort each group by SASA
    hydrophobic.sort(key=lambda x: x.get("relative_sasa", 0), reverse=True)
    other.sort(key=lambda x: x.get("relative_sasa", 0), reverse=True)

    # Take hydrophobic first, then fill with others
    selected = hydrophobic[:max_count]
    remaining = max_count - len(selected)
    if remaining > 0:
        selected.extend(other[:remaining])

    return selected


def calculate_radius_of_gyration(pdb_content: str, chain: str = "B") -> Dict[str, Any]:
    """
    Calculate radius of gyration for a protein chain.

    Lower Rg = more compact. Higher Rg = more elongated/wrap-around.

    Typical values for globular proteins:
        - Rg ≈ 2.2 * N^0.38 (Flory scaling for compact proteins)
        - Compact binder (50 res): ~12-15 Å
        - Elongated binder (50 res): ~20-25 Å

    Args:
        pdb_content: PDB file content
        chain: Chain to analyze

    Returns:
        Dict with:
            - rg: Radius of gyration (Angstroms)
            - n_residues: Number of residues
            - expected_rg: Expected Rg for globular protein of this size
            - rg_ratio: actual_rg / expected_rg (>1.5 suggests elongated)
            - is_compact: Whether design is compact (rg_ratio < 1.5)
    """
    try:
        # Extract CA coordinates
        ca_coords = _extract_ca_coords(pdb_content, chain)

        if len(ca_coords) < 3:
            return {
                "status": "error",
                "error": f"Not enough CA atoms found in chain {chain}",
            }

        # Convert to array
        coords = np.array(list(ca_coords.values()))
        n_residues = len(coords)

        # Calculate center of mass
        center = coords.mean(axis=0)

        # Calculate Rg
        squared_distances = np.sum((coords - center) ** 2, axis=1)
        rg = np.sqrt(np.mean(squared_distances))

        # Expected Rg for globular protein (Flory scaling)
        expected_rg = 2.2 * (n_residues ** 0.38)

        # Rg ratio
        rg_ratio = rg / expected_rg

        return {
            "status": "completed",
            "rg": round(float(rg), 2),
            "n_residues": n_residues,
            "expected_rg": round(float(expected_rg), 2),
            "rg_ratio": round(float(rg_ratio), 3),
            "is_compact": rg_ratio < 1.5,
            "assessment": "compact" if rg_ratio < 1.3 else ("acceptable" if rg_ratio < 1.5 else "elongated"),
        }

    except Exception as e:
        return {
            "status": "error",
            "error": str(e),
        }


def calculate_binding_angle_spread(
    pdb_content: str,
    target_chain: str = "A",
    binder_chain: str = "B",
    contact_distance: float = 8.0
) -> Dict[str, Any]:
    """
    Calculate the angular spread of binder contacts around the target.

    A wrap-around binder will have contacts spread over a wide angle (>180°),
    while a well-positioned binder contacts from one side (<120°).

    Args:
        pdb_content: PDB file content
        target_chain: Target protein chain
        binder_chain: Binder chain
        contact_distance: Distance threshold for contacts

    Returns:
        Dict with angular_spread (degrees), is_wrap_around (bool)
    """
    try:
        # Extract CA coordinates for both chains
        target_cas = _extract_ca_coords(pdb_content, target_chain)
        binder_cas = _extract_ca_coords(pdb_content, binder_chain)

        if len(target_cas) < 3 or len(binder_cas) < 3:
            return {"status": "error", "error": "Not enough atoms"}

        target_coords = np.array(list(target_cas.values()))
        binder_coords = np.array(list(binder_cas.values()))

        # Find target center
        target_center = target_coords.mean(axis=0)

        # Find binder residues that contact the target
        contact_vectors = []
        for binder_coord in binder_coords:
            # Check if this binder residue contacts any target residue
            distances = np.linalg.norm(target_coords - binder_coord, axis=1)
            if np.min(distances) <= contact_distance:
                # Vector from target center to contact point
                vec = binder_coord - target_center
                vec = vec / (np.linalg.norm(vec) + 1e-8)  # Normalize
                contact_vectors.append(vec)

        if len(contact_vectors) < 2:
            return {
                "status": "completed",
                "angular_spread": 0.0,
                "n_contacts": len(contact_vectors),
                "is_wrap_around": False,
                "assessment": "insufficient_contacts"
            }

        contact_vectors = np.array(contact_vectors)

        # Calculate pairwise angles between contact vectors
        max_angle = 0.0
        for i in range(len(contact_vectors)):
            for j in range(i + 1, len(contact_vectors)):
                dot = np.clip(np.dot(contact_vectors[i], contact_vectors[j]), -1, 1)
                angle = np.degrees(np.arccos(dot))
                max_angle = max(max_angle, angle)

        # Also check the "cone angle" - spread from mean direction
        mean_direction = contact_vectors.mean(axis=0)
        mean_direction = mean_direction / (np.linalg.norm(mean_direction) + 1e-8)

        angles_from_mean = []
        for vec in contact_vectors:
            dot = np.clip(np.dot(vec, mean_direction), -1, 1)
            angles_from_mean.append(np.degrees(np.arccos(dot)))

        cone_half_angle = max(angles_from_mean) if angles_from_mean else 0

        # Wrap-around if max pairwise angle > 150° or cone angle > 90°
        is_wrap_around = max_angle > 150 or cone_half_angle > 90

        return {
            "status": "completed",
            "angular_spread": round(float(max_angle), 1),
            "cone_half_angle": round(float(cone_half_angle), 1),
            "n_contacts": len(contact_vectors),
            "is_wrap_around": is_wrap_around,
            "assessment": "wrap_around" if is_wrap_around else "single_sided"
        }

    except Exception as e:
        return {"status": "error", "error": str(e)}


def filter_wrap_around_designs(
    designs: List[Dict],
    binder_chain: str = "B",
    target_chain: str = "A",
    max_rg_ratio: float = 1.3,  # Tighter threshold
    check_angular_spread: bool = True,
    max_angular_spread: float = 150.0,
) -> Tuple[List[Dict], List[Dict]]:
    """
    Filter out designs with elongated/wrap-around binders.

    Uses two metrics:
    1. Rg ratio: Measures binder compactness (>1.3 = elongated)
    2. Angular spread: Measures how binder contacts spread around target (>150° = wrap-around)

    Args:
        designs: List of design dicts with 'pdb_content' field
        binder_chain: Chain ID for the binder (default: "B")
        target_chain: Chain ID for the target (default: "A")
        max_rg_ratio: Maximum allowed Rg/expected_Rg ratio (default: 1.3, stricter)
        check_angular_spread: Also check angular spread of contacts
        max_angular_spread: Maximum allowed angular spread (default: 150°)

    Returns:
        Tuple of (accepted_designs, rejected_designs)
    """
    accepted = []
    rejected = []

    for design in designs:
        pdb_content = design.get("pdb_content", "")

        if not pdb_content:
            rejected.append(design)
            continue

        # Check 1: Rg ratio (binder compactness)
        rg_result = calculate_radius_of_gyration(pdb_content, chain=binder_chain)
        design["rg_analysis"] = rg_result

        rg_passed = True
        rg_ratio = None
        if rg_result.get("status") == "completed":
            rg_ratio = rg_result.get("rg_ratio", 0)
            design["rg_ratio"] = rg_ratio
            design["rg"] = rg_result.get("rg")
            design["rg_assessment"] = rg_result.get("assessment")
            if rg_ratio > max_rg_ratio:
                rg_passed = False

        # Check 2: Angular spread (binding mode)
        angular_passed = True
        if check_angular_spread:
            angular_result = calculate_binding_angle_spread(
                pdb_content, target_chain=target_chain, binder_chain=binder_chain
            )
            design["angular_analysis"] = angular_result

            if angular_result.get("status") == "completed":
                angular_spread = angular_result.get("angular_spread", 0)
                design["angular_spread"] = angular_spread
                design["is_wrap_around"] = angular_result.get("is_wrap_around", False)

                if angular_result.get("is_wrap_around", False):
                    angular_passed = False

        # Decide: reject if either check fails
        if rg_passed and angular_passed:
            accepted.append(design)
        else:
            reasons = []
            if not rg_passed:
                reasons.append(f"Rg ratio {rg_ratio:.2f} > {max_rg_ratio} (elongated)")
            if not angular_passed:
                reasons.append(f"Angular spread > 150° (wrap-around binding)")
            design["rejection_reason"] = "; ".join(reasons)
            rejected.append(design)

    return accepted, rejected


def format_hotspots_for_rfd3(hotspots: List[str]) -> Dict[str, str]:
    """
    Convert hotspot list to RFD3 format.

    RFD3 expects hotspots as: {"A25": "ALL", "A30": "ALL", ...}

    Args:
        hotspots: List of residue IDs like ["A25", "A30", "A35"]

    Returns:
        Dict mapping residue ID to atom selection ("ALL" for all atoms)
    """
    return {res_id: "ALL" for res_id in hotspots}


if __name__ == "__main__":
    # Test with sample PDB
    import argparse

    parser = argparse.ArgumentParser(description="Detect binding hotspots in a protein")
    parser.add_argument("pdb_file", help="Path to PDB file")
    parser.add_argument("--chain", default="A", help="Target chain (default: A)")
    parser.add_argument("--method", default="exposed_clustered",
                        choices=["exposed", "exposed_clustered", "patch", "interface_like"],
                        help="Detection method")
    parser.add_argument("--max-hotspots", type=int, default=5, help="Maximum hotspots to return")

    args = parser.parse_args()

    with open(args.pdb_file) as f:
        pdb_content = f.read()

    result = detect_hotspots_sasa(
        pdb_content=pdb_content,
        target_chain=args.chain,
        method=args.method,
        max_hotspots=args.max_hotspots,
    )

    import json
    print(json.dumps(result, indent=2))
