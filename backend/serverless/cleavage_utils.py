"""
Cleavage Utilities for Interface Ligand Design

Implements the Cleavable Monomer Algorithm:
1. Design a monomer with ligand buried in center
2. Find loop regions where backbone crosses near ligand
3. Cleave at optimal site to create dimer with ligand at interface

This approach achieves stronger binding than direct symmetric design
because the binding pocket is designed as a complete unit, then split.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Any, Optional, Tuple
import sys
import os

# Add parent paths for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'utils'))

from binding_analysis import (
    count_ligand_contacts_by_residue,
    get_ligand_centroid,
    get_ca_positions,
    analyze_interface,
    check_steric_clashes,
    run_gnina_scoring,
    smiles_to_sdf,
    to_python_types,
)


@dataclass
class CleavageSite:
    """Represents a candidate cleavage site in the protein backbone."""
    residue_index: int  # 0-based index in residue list
    residue_id: int  # Actual PDB residue number
    residue_name: str  # 3-letter residue code
    chain: str  # Chain ID
    distance_to_ligand: float  # Distance from CA to ligand centroid
    secondary_structure: str  # 'H', 'E', or 'C'
    contacts_chain_a: int  # Ligand contacts for residues before this
    contacts_chain_b: int  # Ligand contacts for residues after this
    chain_a_length: int  # Number of residues before cleavage
    chain_b_length: int  # Number of residues after cleavage
    score: float = 0.0  # Computed score (higher = better)


@dataclass
class CleavageResult:
    """Result of cleaving a protein at a specific site."""
    dimer_pdb: str  # PDB content of resulting dimer
    chain_a_length: int  # Residues in chain A
    chain_b_length: int  # Residues in chain B
    loop_start: int  # First residue removed
    loop_end: int  # Last residue removed
    residues_removed: int  # Total residues removed
    cleavage_site: CleavageSite  # The site used


def find_cleavage_sites(
    pdb_content: str,
    ligand_name: str = "UNL",
    min_chain_length: int = 25,
    distance_range: Tuple[float, float] = (4.0, 12.0),
    min_contacts_per_chain: int = 3,
    contact_cutoff: float = 5.0,
) -> List[CleavageSite]:
    """
    Find candidate cleavage sites where backbone crosses near ligand.

    A good cleavage site:
    1. Is in a loop (not helix/sheet) - preserves secondary structure
    2. Is close to ligand (backbone crosses binding pocket)
    3. Results in both chains having ligand contacts
    4. Results in both chains having sufficient length

    Args:
        pdb_content: PDB file content as string
        ligand_name: Residue name of ligand (default "UNL")
        min_chain_length: Minimum residues per chain after cleavage (default 25)
        distance_range: (min, max) distance from CA to ligand centroid
        min_contacts_per_chain: Minimum contacts required per chain (default 3)
        contact_cutoff: Distance cutoff for contact counting (default 5.0)

    Returns:
        List of CleavageSite objects, sorted by residue index
    """
    try:
        from dssp import assign_secondary_structure, get_ss_codes_by_residue
    except ImportError:
        # Try relative import
        try:
            from backend.utils.dssp import assign_secondary_structure, get_ss_codes_by_residue
        except ImportError:
            print("[CleavageUtils] Warning: Could not import DSSP, using fallback")
            get_ss_codes_by_residue = lambda x: {}
            assign_secondary_structure = lambda x: {"success": False}

    # Get ligand centroid
    ligand_centroid = get_ligand_centroid(pdb_content, ligand_name)
    if ligand_centroid is None:
        print(f"[CleavageUtils] Warning: Ligand '{ligand_name}' not found")
        return []

    ligand_centroid = np.array(ligand_centroid)

    # Get CA positions
    ca_positions = get_ca_positions(pdb_content)
    if not ca_positions:
        print("[CleavageUtils] Warning: No CA atoms found")
        return []

    n_residues = len(ca_positions)

    # Get secondary structure
    ss_result = assign_secondary_structure(pdb_content)
    ss_codes = get_ss_codes_by_residue(ss_result)

    # Get per-residue ligand contacts
    residue_contacts = count_ligand_contacts_by_residue(
        pdb_content, ligand_name, contact_cutoff
    )

    # Build cumulative contact counts for efficient range queries
    # cumulative[i] = total contacts for residues 0..i-1
    res_ids = [r[0] for r in ca_positions]
    contact_counts = [residue_contacts.get(r, 0) for r in res_ids]
    cumulative = [0]
    for c in contact_counts:
        cumulative.append(cumulative[-1] + c)

    candidates = []
    min_dist, max_dist = distance_range

    for i in range(min_chain_length, n_residues - min_chain_length):
        res_id = res_ids[i]
        pos = ca_positions[i][1]

        # Check 1: Is this residue in a loop?
        ss = ss_codes.get(res_id, 'C')
        if ss != 'C':  # Not coil/loop
            continue

        # Check 2: Is this residue in distance range to ligand?
        dist_to_ligand = float(np.linalg.norm(pos - ligand_centroid))
        if dist_to_ligand < min_dist or dist_to_ligand > max_dist:
            continue

        # Check 3: Count contacts for each resulting chain
        # Chain A: residues 0..i-1 (cumulative[i] - cumulative[0])
        # Chain B: residues i+1..N-1 (cumulative[N] - cumulative[i+1])
        contacts_a = cumulative[i] - cumulative[0]
        contacts_b = cumulative[n_residues] - cumulative[i + 1]

        if contacts_a < min_contacts_per_chain or contacts_b < min_contacts_per_chain:
            continue

        # Get residue name from PDB
        res_name = get_residue_name(pdb_content, res_id)

        # Valid candidate
        candidates.append(CleavageSite(
            residue_index=i,
            residue_id=res_id,
            residue_name=res_name,
            chain="A",  # Monomer is typically chain A
            distance_to_ligand=dist_to_ligand,
            secondary_structure=ss,
            contacts_chain_a=contacts_a,
            contacts_chain_b=contacts_b,
            chain_a_length=i,
            chain_b_length=n_residues - i - 1
        ))

    return candidates


def get_residue_name(pdb_content: str, res_id: int) -> str:
    """Extract residue name from PDB for a given residue ID."""
    for line in pdb_content.split('\n'):
        if line.startswith('ATOM'):
            try:
                line_res_id = int(line[22:26].strip())
                if line_res_id == res_id:
                    return line[17:20].strip()
            except ValueError:
                continue
    return "UNK"


def score_cleavage_site(
    site: CleavageSite,
    n_residues: int,
    weights: Optional[Dict[str, float]] = None
) -> float:
    """
    Score a cleavage site. Higher score = better site.

    Scoring factors:
    1. Contact score: total contacts with ligand
    2. Symmetry: balanced contacts between chains
    3. Length balance: prefer cuts near middle
    4. Distance penalty: very close to ligand may cause issues

    Args:
        site: CleavageSite to score
        n_residues: Total residues in protein
        weights: Optional custom weights for scoring factors

    Returns:
        Composite score (higher = better)
    """
    # Default weights
    w = weights or {
        "contacts": 1.0,
        "symmetry": 2.0,
        "length_balance": 0.5,
        "distance_penalty": 2.0
    }

    # Contact score: more contacts = better
    contact_score = (site.contacts_chain_a + site.contacts_chain_b) * w["contacts"]

    # Symmetry score: balanced contacts = better (0-1)
    max_contacts = max(site.contacts_chain_a, site.contacts_chain_b)
    min_contacts = min(site.contacts_chain_a, site.contacts_chain_b)
    symmetry = min_contacts / max_contacts if max_contacts > 0 else 0
    symmetry_score = symmetry * w["symmetry"] * 10  # Scale up

    # Length balance: prefer cuts near middle (0-1)
    length_balance = 1 - abs(site.residue_index - n_residues / 2) / (n_residues / 2)
    length_score = length_balance * w["length_balance"] * 10

    # Distance penalty: very close to ligand might cause issues
    if site.distance_to_ligand < 5.0:
        distance_penalty = (5.0 - site.distance_to_ligand) * w["distance_penalty"]
    else:
        distance_penalty = 0

    return contact_score + symmetry_score + length_score - distance_penalty


def select_best_cleavage_site(
    candidates: List[CleavageSite],
    n_residues: int,
    strategy: str = "balanced"
) -> Optional[CleavageSite]:
    """
    Select the best cleavage site from candidates.

    Args:
        candidates: List of candidate CleavageSite objects
        n_residues: Total residues in protein
        strategy: Selection strategy:
            - "balanced": Balance contacts, symmetry, and length
            - "symmetric": Prioritize contact symmetry
            - "central": Prioritize central cuts

    Returns:
        Best CleavageSite or None if no candidates
    """
    if not candidates:
        return None

    # Score all candidates
    for site in candidates:
        if strategy == "symmetric":
            weights = {"contacts": 0.5, "symmetry": 4.0, "length_balance": 0.3, "distance_penalty": 1.0}
        elif strategy == "central":
            weights = {"contacts": 0.5, "symmetry": 1.0, "length_balance": 3.0, "distance_penalty": 1.0}
        else:  # balanced
            weights = None  # Use defaults

        site.score = score_cleavage_site(site, n_residues, weights)

    # Return highest scoring
    return max(candidates, key=lambda s: s.score)


def cleave_protein(
    pdb_content: str,
    cleavage_site: CleavageSite,
    ligand_name: str = "UNL",
    remove_loop: bool = True,
    loop_window: int = 3,
) -> CleavageResult:
    """
    Cleave protein at specified site to create two chains.

    Args:
        pdb_content: PDB file content
        cleavage_site: CleavageSite specifying where to cut
        ligand_name: Ligand residue name
        remove_loop: If True, remove entire loop around cleavage site
        loop_window: Residues to check on each side for loop removal

    Returns:
        CleavageResult with dimer PDB and metadata
    """
    try:
        from dssp import assign_secondary_structure, get_ss_codes_by_residue
    except ImportError:
        try:
            from backend.utils.dssp import assign_secondary_structure, get_ss_codes_by_residue
        except ImportError:
            get_ss_codes_by_residue = lambda x: {}
            assign_secondary_structure = lambda x: {"success": False}

    cleavage_res_id = cleavage_site.residue_id

    if remove_loop:
        # Get secondary structure to find loop boundaries
        ss_result = assign_secondary_structure(pdb_content)
        ss_codes = get_ss_codes_by_residue(ss_result)

        # Expand cleavage to full loop
        loop_start = cleavage_res_id
        loop_end = cleavage_res_id

        # Get all residue IDs
        all_res_ids = sorted(ss_codes.keys()) if ss_codes else [cleavage_res_id]

        # Expand backward
        for j in range(cleavage_res_id - 1,
                       max(min(all_res_ids) if all_res_ids else 1, cleavage_res_id - loop_window) - 1, -1):
            if ss_codes.get(j, 'C') == 'C':
                loop_start = j
            else:
                break

        # Expand forward
        max_res = max(all_res_ids) if all_res_ids else cleavage_res_id + loop_window
        for j in range(cleavage_res_id + 1,
                       min(max_res + 1, cleavage_res_id + loop_window + 1)):
            if ss_codes.get(j, 'C') == 'C':
                loop_end = j
            else:
                break
    else:
        loop_start = cleavage_res_id
        loop_end = cleavage_res_id

    # Parse PDB and split into chains
    chain_a_lines = []
    chain_b_lines = []
    ligand_lines = []
    other_lines = []

    min_res_b = None  # Track minimum residue ID for chain B renumbering

    for line in pdb_content.split('\n'):
        if line.startswith('ATOM'):
            try:
                res_id = int(line[22:26].strip())

                if res_id < loop_start:
                    # Chain A: before loop
                    new_line = line[:21] + 'A' + line[22:]
                    chain_a_lines.append(new_line)
                elif res_id > loop_end:
                    # Chain B: after loop
                    if min_res_b is None:
                        min_res_b = res_id
                    # Renumber residue starting from 1
                    new_res_id = res_id - min_res_b + 1
                    new_line = line[:21] + 'B' + f"{new_res_id:4d}" + line[26:]
                    chain_b_lines.append(new_line)
                # else: residue is in the loop - skip it

            except ValueError:
                other_lines.append(line)

        elif line.startswith('HETATM'):
            # Keep ligand with chain L
            res_name = line[17:20].strip()
            if res_name == ligand_name:
                new_line = line[:21] + 'L' + line[22:]
                ligand_lines.append(new_line)
            else:
                other_lines.append(line)

        elif line.startswith('TER') or line.startswith('END'):
            pass  # Will add our own TER/END records

        else:
            if line.strip():
                other_lines.append(line)

    # Count chain lengths
    chain_a_res_ids = set()
    chain_b_res_ids = set()

    for line in chain_a_lines:
        try:
            res_id = int(line[22:26].strip())
            chain_a_res_ids.add(res_id)
        except ValueError:
            pass

    for line in chain_b_lines:
        try:
            res_id = int(line[22:26].strip())
            chain_b_res_ids.add(res_id)
        except ValueError:
            pass

    chain_a_length = len(chain_a_res_ids)
    chain_b_length = len(chain_b_res_ids)

    # Assemble dimer PDB
    output_lines = []

    # Header
    output_lines.append("REMARK   Cleaved dimer from cleavable monomer algorithm")
    output_lines.append(f"REMARK   Cleavage site: residue {cleavage_res_id}")
    output_lines.append(f"REMARK   Loop removed: {loop_start}-{loop_end}")

    # Chain A
    output_lines.extend(chain_a_lines)
    if chain_a_lines:
        output_lines.append("TER")

    # Chain B
    output_lines.extend(chain_b_lines)
    if chain_b_lines:
        output_lines.append("TER")

    # Ligand
    output_lines.extend(ligand_lines)
    if ligand_lines:
        output_lines.append("TER")

    output_lines.append("END")

    dimer_pdb = '\n'.join(output_lines)

    return CleavageResult(
        dimer_pdb=dimer_pdb,
        chain_a_length=chain_a_length,
        chain_b_length=chain_b_length,
        loop_start=loop_start,
        loop_end=loop_end,
        residues_removed=loop_end - loop_start + 1,
        cleavage_site=cleavage_site
    )


def validate_cleaved_dimer(
    pdb_content: str,
    ligand_name: str = "UNL",
    ligand_smiles: Optional[str] = None,
    min_contacts_per_chain: int = 3,
    max_affinity: float = -5.0,
    run_gnina: bool = True,
) -> Dict[str, Any]:
    """
    Validate a cleaved dimer for structural and binding quality.

    Checks:
    1. Chain lengths are sufficient
    2. No interchain clashes
    3. Both chains contact ligand
    4. Contact symmetry is acceptable
    5. GNINA affinity meets threshold (if enabled)

    Args:
        pdb_content: Dimer PDB content
        ligand_name: Ligand residue name
        ligand_smiles: SMILES string for GNINA scoring
        min_contacts_per_chain: Minimum contacts per chain
        max_affinity: Maximum (least negative) acceptable GNINA affinity
        run_gnina: Whether to run GNINA scoring

    Returns:
        Validation results with pass/fail for each criterion
    """
    results = {
        "status": "completed",
        "checks": {},
        "overall_pass": False,
    }

    # Parse structure to count chain lengths
    chain_a_res = set()
    chain_b_res = set()
    has_ligand = False

    for line in pdb_content.split('\n'):
        if line.startswith('ATOM'):
            chain = line[21]
            res_id = int(line[22:26].strip())
            if chain == 'A':
                chain_a_res.add(res_id)
            elif chain == 'B':
                chain_b_res.add(res_id)
        elif line.startswith('HETATM'):
            res_name = line[17:20].strip()
            if res_name == ligand_name:
                has_ligand = True

    chain_a_length = len(chain_a_res)
    chain_b_length = len(chain_b_res)

    # Check 1: Chain lengths
    results["checks"]["chain_a_length"] = chain_a_length
    results["checks"]["chain_b_length"] = chain_b_length
    results["checks"]["length_pass"] = min(chain_a_length, chain_b_length) >= 20

    # Check 2: Ligand present
    results["checks"]["ligand_present"] = has_ligand

    if not has_ligand:
        results["checks"]["contact_pass"] = False
        results["checks"]["symmetry_pass"] = False
        results["overall_pass"] = False
        results["error"] = f"Ligand '{ligand_name}' not found in structure"
        return to_python_types(results)

    # Check 3: Contacts per chain
    # Count contacts from chain A and B separately to ligand
    contacts_a = 0
    contacts_b = 0
    contact_cutoff = 4.0

    # Parse ligand coordinates
    ligand_coords = []
    for line in pdb_content.split('\n'):
        if line.startswith('HETATM'):
            res_name = line[17:20].strip()
            if res_name == ligand_name:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                ligand_coords.append(np.array([x, y, z]))

    # Parse chain A and B coordinates
    for line in pdb_content.split('\n'):
        if line.startswith('ATOM'):
            chain = line[21]
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            pos = np.array([x, y, z])

            for lig_pos in ligand_coords:
                if np.linalg.norm(pos - lig_pos) < contact_cutoff:
                    if chain == 'A':
                        contacts_a += 1
                    elif chain == 'B':
                        contacts_b += 1
                    break

    results["checks"]["contacts_a"] = contacts_a
    results["checks"]["contacts_b"] = contacts_b
    results["checks"]["contact_pass"] = (
        contacts_a >= min_contacts_per_chain and
        contacts_b >= min_contacts_per_chain
    )

    # Check 4: Symmetry
    total_contacts = contacts_a + contacts_b
    if total_contacts > 0:
        symmetry = 1.0 - abs(contacts_a - contacts_b) / total_contacts
    else:
        symmetry = 0.0
    results["checks"]["symmetry_score"] = symmetry
    results["checks"]["symmetry_pass"] = symmetry > 0.4

    # Check 5: No clashes between chains
    clash_result = check_steric_clashes(pdb_content)
    results["checks"]["clash_check"] = clash_result
    results["checks"]["clash_pass"] = not clash_result.get("has_clashes", True)

    # Check 6: GNINA affinity (optional)
    if run_gnina and ligand_smiles:
        ligand_sdf = smiles_to_sdf(ligand_smiles)
        if ligand_sdf:
            # Extract protein only for receptor
            receptor_lines = []
            for line in pdb_content.split('\n'):
                if line.startswith('ATOM') or line.startswith('TER') or line.startswith('END'):
                    receptor_lines.append(line)
            receptor_pdb = '\n'.join(receptor_lines)

            gnina_result = run_gnina_scoring(
                receptor_pdb=receptor_pdb,
                ligand_sdf=ligand_sdf,
                whole_protein=True,
            )
            results["checks"]["gnina_result"] = gnina_result

            if gnina_result.get("status") == "completed":
                affinity = gnina_result.get("result", {}).get("best_affinity")
                results["checks"]["gnina_affinity"] = affinity
                results["checks"]["affinity_pass"] = affinity is not None and affinity < max_affinity
            else:
                results["checks"]["affinity_pass"] = None  # Could not determine
        else:
            results["checks"]["affinity_pass"] = None
    else:
        results["checks"]["affinity_pass"] = None

    # Overall pass
    required_checks = ["length_pass", "contact_pass", "clash_pass"]
    optional_checks = ["symmetry_pass", "affinity_pass"]

    all_required = all(results["checks"].get(c, False) for c in required_checks)
    results["overall_pass"] = all_required

    return to_python_types(results)


def create_homo_dimer(
    pdb_content: str,
    ligand_name: str = "UNL",
    chain_to_use: str = "larger",
) -> Tuple[Optional[str], Dict[str, Any]]:
    """
    Create a homo-dimer from hetero-dimer by applying C2 symmetry.

    Takes one chain and creates its C2-symmetric copy. The ligand
    must be on or near the C2 axis for this to work properly.

    Args:
        pdb_content: Hetero-dimer PDB content
        ligand_name: Ligand residue name
        chain_to_use: "larger", "A", or "B" - which chain to duplicate

    Returns:
        (homo_dimer_pdb, metadata) or (None, error_dict)
    """
    # Parse chains
    chain_a_lines = []
    chain_b_lines = []
    ligand_lines = []

    chain_a_coords = []
    chain_b_coords = []
    ligand_coords = []

    for line in pdb_content.split('\n'):
        if line.startswith('ATOM'):
            chain = line[21]
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            if chain == 'A':
                chain_a_lines.append(line)
                chain_a_coords.append([x, y, z])
            elif chain == 'B':
                chain_b_lines.append(line)
                chain_b_coords.append([x, y, z])

        elif line.startswith('HETATM'):
            res_name = line[17:20].strip()
            if res_name == ligand_name:
                ligand_lines.append(line)
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                ligand_coords.append([x, y, z])

    if not ligand_coords:
        return None, {"error": f"Ligand '{ligand_name}' not found"}

    # Calculate ligand center (C2 axis passes through here)
    ligand_center = np.mean(ligand_coords, axis=0)

    # Select template chain
    if chain_to_use == "larger":
        template_lines = chain_a_lines if len(chain_a_lines) > len(chain_b_lines) else chain_b_lines
        template_coords = chain_a_coords if len(chain_a_lines) > len(chain_b_lines) else chain_b_coords
    elif chain_to_use == "A":
        template_lines = chain_a_lines
        template_coords = chain_a_coords
    else:
        template_lines = chain_b_lines
        template_coords = chain_b_coords

    if not template_lines:
        return None, {"error": f"Chain {chain_to_use} not found or empty"}

    # Apply C2 rotation (180° around Z axis through ligand center)
    # C2: (x, y, z) -> (2*cx - x, 2*cy - y, z)
    rotated_lines = []
    for i, line in enumerate(template_lines):
        x, y, z = template_coords[i]

        # Apply C2 transformation
        new_x = 2 * ligand_center[0] - x
        new_y = 2 * ligand_center[1] - y
        new_z = z  # Z unchanged for C2 around Z

        # Create new line with chain B and new coordinates
        new_line = (
            line[:21] + 'B' +
            f"{new_x:8.3f}{new_y:8.3f}{new_z:8.3f}" +
            line[54:]
        )
        rotated_lines.append(new_line)

    # Check for clashes between template and rotated
    min_dist = float('inf')
    for coord in template_coords:
        for i, line in enumerate(template_lines):
            rot_coord = [
                2 * ligand_center[0] - template_coords[i][0],
                2 * ligand_center[1] - template_coords[i][1],
                template_coords[i][2]
            ]
            dist = np.linalg.norm(np.array(coord) - np.array(rot_coord))
            if dist < min_dist:
                min_dist = dist

    if min_dist < 2.5:
        return None, {
            "error": "Chains clash after C2 symmetry",
            "min_distance": min_dist
        }

    # Assemble homo-dimer PDB
    output_lines = []
    output_lines.append("REMARK   Homo-dimer from C2 symmetry")
    output_lines.append(f"REMARK   C2 axis at: {ligand_center[0]:.2f}, {ligand_center[1]:.2f}, Z")

    # Chain A (original template, relabeled)
    for line in template_lines:
        output_lines.append(line[:21] + 'A' + line[22:])
    output_lines.append("TER")

    # Chain B (rotated copy)
    output_lines.extend(rotated_lines)
    output_lines.append("TER")

    # Ligand
    output_lines.extend(ligand_lines)
    output_lines.append("TER")

    output_lines.append("END")

    homo_pdb = '\n'.join(output_lines)

    metadata = {
        "min_interchain_distance": min_dist,
        "chain_length": len(set(int(line[22:26].strip()) for line in template_lines if line.startswith('ATOM'))),
        "c2_center": ligand_center.tolist(),
    }

    return homo_pdb, to_python_types(metadata)


def count_residues(pdb_content: str) -> int:
    """Count unique residues in a PDB structure."""
    res_ids = set()
    for line in pdb_content.split('\n'):
        if line.startswith('ATOM'):
            try:
                res_id = int(line[22:26].strip())
                res_ids.add(res_id)
            except ValueError:
                pass
    return len(res_ids)
