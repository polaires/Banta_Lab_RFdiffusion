"""
DSSP Secondary Structure Assignment Module

Assign secondary structure to protein residues using DSSP algorithm.
Uses mdtraj or biotite if available, with fallback to simple H-bond analysis.
"""

import tempfile
import os
from typing import Dict, List, Any, Optional, Tuple
import numpy as np


# Try to import DSSP libraries
MDTRAJ_AVAILABLE = False
BIOTITE_DSSP_AVAILABLE = False

try:
    import mdtraj
    MDTRAJ_AVAILABLE = True
except ImportError:
    pass

try:
    from biotite.structure import annotate_sse
    from biotite.structure.io.pdb import PDBFile
    import io
    BIOTITE_DSSP_AVAILABLE = True
except ImportError:
    pass


# Secondary structure codes
SS_CODES = {
    "H": "alpha_helix",
    "G": "310_helix",
    "I": "pi_helix",
    "E": "beta_sheet",
    "B": "beta_bridge",
    "T": "turn",
    "S": "bend",
    "C": "coil",
    " ": "coil",
    "": "coil",
    "a": "alpha_helix",  # Biotite uses lowercase
    "b": "beta_sheet",
    "c": "coil",
}


def _assign_ss_mdtraj(pdb_content: str) -> Dict[str, Any]:
    """
    Assign secondary structure using MDTraj's DSSP implementation.

    Args:
        pdb_content: PDB file content

    Returns:
        Dict with secondary structure assignments
    """
    # Write to temp file
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode='w') as f:
        f.write(pdb_content)
        temp_path = f.name

    try:
        # Load structure
        traj = mdtraj.load(temp_path)

        # Compute DSSP
        dssp = mdtraj.compute_dssp(traj, simplified=False)[0]

        # Get topology info
        topology = traj.topology

        per_residue = {}
        for i, residue in enumerate(topology.residues):
            key = f"{residue.chain.chain_id}:{residue.name}{residue.resSeq}"
            ss_code = dssp[i] if i < len(dssp) else "C"

            per_residue[key] = {
                "chain": str(residue.chain.chain_id),
                "residue_name": residue.name,
                "residue_number": residue.resSeq,
                "ss_code": ss_code,
                "ss_type": SS_CODES.get(ss_code, "coil"),
            }

        return {
            "success": True,
            "method": "mdtraj_dssp",
            "per_residue": per_residue,
        }

    finally:
        os.unlink(temp_path)


def _assign_ss_biotite(pdb_content: str) -> Dict[str, Any]:
    """
    Assign secondary structure using Biotite's SSE annotation.

    Args:
        pdb_content: PDB file content

    Returns:
        Dict with secondary structure assignments
    """
    # Parse PDB
    pdb_file = PDBFile.read(io.StringIO(pdb_content))
    structure = pdb_file.get_structure(model=1)

    # Filter to protein atoms only
    from biotite.structure import filter_amino_acids
    protein = structure[filter_amino_acids(structure)]

    if len(protein) == 0:
        return {
            "success": False,
            "error": "No protein atoms found",
        }

    # Annotate secondary structure
    sse = annotate_sse(protein)

    # Get unique residues
    from biotite.structure import get_residue_starts
    res_starts = get_residue_starts(protein)

    per_residue = {}
    for i, start_idx in enumerate(res_starts):
        chain = protein.chain_id[start_idx]
        res_name = protein.res_name[start_idx]
        res_num = protein.res_id[start_idx]
        key = f"{chain}:{res_name}{res_num}"

        ss_code = sse[i] if i < len(sse) else "c"
        per_residue[key] = {
            "chain": chain,
            "residue_name": res_name,
            "residue_number": int(res_num),
            "ss_code": ss_code,
            "ss_type": SS_CODES.get(ss_code, "coil"),
        }

    return {
        "success": True,
        "method": "biotite_sse",
        "per_residue": per_residue,
    }


def _assign_ss_simple(pdb_content: str) -> Dict[str, Any]:
    """
    Simple secondary structure estimation based on backbone geometry.
    Fallback when no DSSP library is available.

    Uses backbone dihedral angles to estimate structure:
    - Alpha helix: phi ~ -60, psi ~ -45
    - Beta sheet: phi ~ -120, psi ~ 120

    Args:
        pdb_content: PDB file content

    Returns:
        Dict with estimated secondary structure
    """
    # Parse backbone atoms
    backbone_atoms = {}  # residue_key -> {N, CA, C}

    for line in pdb_content.split('\n'):
        if not line.startswith('ATOM'):
            continue

        try:
            atom_name = line[12:16].strip()
            if atom_name not in ['N', 'CA', 'C']:
                continue

            chain = line[21]
            res_name = line[17:20].strip()
            res_num = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            key = (chain, res_num)
            if key not in backbone_atoms:
                backbone_atoms[key] = {
                    "chain": chain,
                    "res_name": res_name,
                    "res_num": res_num,
                }
            backbone_atoms[key][atom_name] = np.array([x, y, z])

        except (ValueError, IndexError):
            continue

    # Sort by residue number
    sorted_residues = sorted(backbone_atoms.keys(), key=lambda x: (x[0], x[1]))

    per_residue = {}
    for i, key in enumerate(sorted_residues):
        data = backbone_atoms[key]
        res_key = f"{data['chain']}:{data['res_name']}{data['res_num']}"

        # Default to coil
        ss_code = "C"

        # Check if we have enough atoms for dihedral calculation
        if all(atom in data for atom in ['N', 'CA', 'C']):
            # Check for helix pattern (H-bond pattern would be i+4)
            if i + 4 < len(sorted_residues):
                next_key = sorted_residues[i + 4]
                if next_key[0] == key[0]:  # Same chain
                    next_data = backbone_atoms[next_key]
                    if 'N' in next_data and 'C' in data:
                        # Check C=O...H-N distance (simplified)
                        c_pos = data['C']
                        n_pos = next_data['N']
                        distance = np.linalg.norm(c_pos - n_pos)
                        if distance < 4.5:  # Approximate H-bond distance
                            ss_code = "H"

            # Check for sheet pattern (parallel or antiparallel)
            # This is very simplified - real DSSP is much more complex
            if ss_code == "C" and i > 0 and i < len(sorted_residues) - 1:
                prev_data = backbone_atoms[sorted_residues[i - 1]]
                next_data = backbone_atoms[sorted_residues[i + 1]]

                if all(a in prev_data for a in ['CA']) and all(a in next_data for a in ['CA']):
                    # Check for extended conformation
                    ca_prev = prev_data['CA']
                    ca_curr = data['CA']
                    ca_next = next_data['CA']

                    v1 = ca_curr - ca_prev
                    v2 = ca_next - ca_curr

                    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-8)
                    if cos_angle > 0.7:  # Extended conformation
                        ss_code = "E"

        per_residue[res_key] = {
            "chain": data["chain"],
            "residue_name": data["res_name"],
            "residue_number": data["res_num"],
            "ss_code": ss_code,
            "ss_type": SS_CODES.get(ss_code, "coil"),
        }

    return {
        "success": True,
        "method": "simple_estimation",
        "per_residue": per_residue,
        "warning": "Using simple estimation - install mdtraj or biotite for accurate DSSP",
    }


def assign_secondary_structure(pdb_content: str) -> Dict[str, Any]:
    """
    Assign secondary structure to all residues.

    Uses MDTraj DSSP if available, falls back to Biotite, then simple estimation.

    Args:
        pdb_content: PDB file content

    Returns:
        Dict with secondary structure assignments including:
        - per_residue: Dict mapping residue IDs to SS assignments
        - ss_segments: List of contiguous SS segments
        - content: Helix/sheet/coil percentages
    """
    # Try different methods
    result = None

    if MDTRAJ_AVAILABLE:
        try:
            result = _assign_ss_mdtraj(pdb_content)
        except Exception as e:
            print(f"MDTraj DSSP failed: {e}, trying biotite...")

    if result is None and BIOTITE_DSSP_AVAILABLE:
        try:
            result = _assign_ss_biotite(pdb_content)
        except Exception as e:
            print(f"Biotite SSE failed: {e}, using simple estimation...")

    if result is None:
        result = _assign_ss_simple(pdb_content)

    if not result.get("success"):
        return result

    # Analyze secondary structure content
    per_residue = result["per_residue"]
    ss_counts = {"helix": 0, "sheet": 0, "coil": 0}

    for data in per_residue.values():
        ss_type = data["ss_type"]
        if "helix" in ss_type:
            ss_counts["helix"] += 1
        elif "sheet" in ss_type or "bridge" in ss_type:
            ss_counts["sheet"] += 1
        else:
            ss_counts["coil"] += 1

    total = sum(ss_counts.values())
    if total > 0:
        content = {k: round(v / total * 100, 1) for k, v in ss_counts.items()}
    else:
        content = {"helix": 0.0, "sheet": 0.0, "coil": 0.0}

    # Extract contiguous segments
    segments = extract_ss_segments(per_residue)

    result["content"] = content
    result["ss_segments"] = segments
    result["residue_count"] = total

    return result


def extract_ss_segments(per_residue: Dict[str, Dict]) -> List[Dict]:
    """
    Extract contiguous secondary structure segments.

    Args:
        per_residue: Dict mapping residue IDs to SS data

    Returns:
        List of segment dictionaries
    """
    # Sort residues by chain and number
    sorted_residues = sorted(
        per_residue.items(),
        key=lambda x: (x[1]["chain"], x[1]["residue_number"])
    )

    segments = []
    current_segment = None

    for res_id, data in sorted_residues:
        ss_type = data["ss_type"]

        # Check if we need to start a new segment
        if current_segment is None or current_segment["ss_type"] != ss_type:
            if current_segment is not None:
                segments.append(current_segment)

            current_segment = {
                "ss_type": ss_type,
                "chain": data["chain"],
                "start_residue": data["residue_number"],
                "end_residue": data["residue_number"],
                "length": 1,
                "residues": [res_id],
            }
        else:
            # Continue current segment
            current_segment["end_residue"] = data["residue_number"]
            current_segment["length"] += 1
            current_segment["residues"].append(res_id)

    if current_segment is not None:
        segments.append(current_segment)

    return segments


def get_ss_context(
    ss_result: Dict[str, Any],
    residue_ids: List[str]
) -> Dict[str, Any]:
    """
    Get secondary structure context for specific residues.

    Useful for understanding structural context of binding site residues.

    Args:
        ss_result: Output from assign_secondary_structure()
        residue_ids: List of residue identifiers

    Returns:
        Dict with SS context for the residues
    """
    if not ss_result.get("success"):
        return ss_result

    per_residue = ss_result["per_residue"]
    context = []

    for res_id in residue_ids:
        # Try to find matching residue
        matched = None
        for key in per_residue:
            if res_id in key or key.endswith(res_id.split(":")[-1] if ":" in res_id else res_id):
                matched = per_residue[key]
                matched["residue_id"] = key
                break

        if matched:
            context.append(matched)
        else:
            context.append({
                "residue_id": res_id,
                "ss_type": "unknown",
                "ss_code": "?",
                "warning": "Residue not found",
            })

    # Summarize context
    ss_summary = {}
    for ctx in context:
        ss_type = ctx.get("ss_type", "unknown")
        ss_summary[ss_type] = ss_summary.get(ss_type, 0) + 1

    return {
        "success": True,
        "context": context,
        "summary": ss_summary,
    }


def suggest_ss_constraints(
    ss_result: Dict[str, Any],
    preserve_types: List[str] = ["alpha_helix", "beta_sheet"]
) -> Dict[str, Any]:
    """
    Suggest RFD3 parameters to preserve secondary structure.

    Args:
        ss_result: Output from assign_secondary_structure()
        preserve_types: SS types to preserve

    Returns:
        Dict with RFD3 constraint suggestions
    """
    if not ss_result.get("success"):
        return ss_result

    segments = ss_result.get("ss_segments", [])

    # Find segments to preserve
    preserve_segments = [
        seg for seg in segments
        if seg["ss_type"] in preserve_types and seg["length"] >= 3
    ]

    # Generate contig constraints
    # Format: fixed regions from SS elements, design regions in between
    constraints = []

    for seg in preserve_segments:
        chain = seg["chain"]
        start = seg["start_residue"]
        end = seg["end_residue"]

        constraints.append({
            "chain": chain,
            "start": start,
            "end": end,
            "ss_type": seg["ss_type"],
            "length": seg["length"],
            "suggested_contig": f"{chain}{start}-{end}",
        })

    return {
        "success": True,
        "preserved_segments": len(preserve_segments),
        "constraints": constraints,
        "suggestion": "Use these residue ranges in 'contig' to preserve secondary structure",
    }


def get_ss_codes_by_residue(ss_result: Dict[str, Any]) -> Dict[int, str]:
    """
    Get secondary structure codes indexed by residue number.

    Normalizes all SS codes to simple H/E/C format:
    - H: helix (includes alpha, 3-10, pi)
    - E: sheet (includes extended, bridge)
    - C: coil (everything else)

    Args:
        ss_result: Output from assign_secondary_structure()

    Returns:
        Dict mapping residue_number -> ss_code ('H', 'E', 'C')
        Returns empty dict if ss_result is not successful.
    """
    if not ss_result.get("success"):
        return {}

    codes = {}
    per_residue = ss_result.get("per_residue", {})

    for key, data in per_residue.items():
        res_num = data.get("residue_number")
        ss_code = data.get("ss_code", "C")

        if res_num is not None:
            # Normalize to H/E/C
            if ss_code in ['H', 'G', 'I', 'a']:  # Helices
                codes[int(res_num)] = 'H'
            elif ss_code in ['E', 'B', 'b']:  # Sheets/bridges
                codes[int(res_num)] = 'E'
            else:  # Coil, turn, bend, or unknown
                codes[int(res_num)] = 'C'

    return codes
