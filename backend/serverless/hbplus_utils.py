#!/usr/bin/env python3
"""
HBPLUS Utilities for Hydrogen Bond Analysis.

HBPLUS (Hydrogen Bond PLot Utility Suite) provides accurate H-bond detection
with geometry validation beyond simple distance cutoffs.

HBPLUS Criteria:
- Donor-Acceptor distance: <= 3.5 A
- Hydrogen-Acceptor distance: <= 2.5 A
- Donor-Hydrogen-Acceptor angle: >= 90 degrees

This module provides:
- check_hbplus_available(): Check if HBPLUS executable is installed
- run_hbplus(): Run HBPLUS and parse output
- get_ligand_hbonds(): Filter for ligand-specific H-bonds
- fallback_distance_hbonds(): Simple distance-based detection (when HBPLUS unavailable)

Installation:
    Download from: https://www.ebi.ac.uk/thornton-srv/software/HBPLUS/
    Set environment: export HBPLUS_PATH=/path/to/hbplus
"""

import os
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from dataclasses import dataclass

# Try to get HBPLUS path from environment
HBPLUS_PATH = os.environ.get("HBPLUS_PATH", "/usr/local/bin/hbplus")


@dataclass
class HydrogenBond:
    """Represents a hydrogen bond detected by HBPLUS."""
    donor_chain: str
    donor_residue: int
    donor_resname: str
    donor_atom: str
    acceptor_chain: str
    acceptor_residue: int
    acceptor_resname: str
    acceptor_atom: str
    da_distance: float  # Donor-Acceptor distance
    ha_distance: float  # Hydrogen-Acceptor distance
    dha_angle: float    # Donor-Hydrogen-Acceptor angle
    category: str       # MM (main-main), MS (main-side), SS (side-side)

    def involves_ligand(self, ligand_name: str = "UNL") -> bool:
        """Check if this H-bond involves the specified ligand."""
        return self.donor_resname == ligand_name or self.acceptor_resname == ligand_name

    def involves_atom(self, ligand_name: str, atom_name: str) -> bool:
        """Check if this H-bond involves a specific ligand atom."""
        if self.donor_resname == ligand_name and self.donor_atom == atom_name:
            return True
        if self.acceptor_resname == ligand_name and self.acceptor_atom == atom_name:
            return True
        return False


def check_hbplus_available() -> Tuple[bool, str]:
    """
    Check if HBPLUS executable is available.

    Returns:
        (available, message) - Boolean and status message
    """
    if not os.path.exists(HBPLUS_PATH):
        # Try common paths
        common_paths = [
            "/usr/local/bin/hbplus",
            "/usr/bin/hbplus",
            "/opt/hbplus/hbplus",
            os.path.expanduser("~/bin/hbplus"),
        ]
        for path in common_paths:
            if os.path.exists(path):
                global HBPLUS_PATH
                HBPLUS_PATH = path
                return True, f"Found HBPLUS at {path}"

        return False, f"HBPLUS not found at {HBPLUS_PATH} or common paths"

    # Check if executable
    try:
        result = subprocess.run(
            [HBPLUS_PATH, "-h"],
            capture_output=True,
            timeout=5
        )
        if result.returncode in [0, 1]:  # HBPLUS returns 1 for -h
            return True, f"HBPLUS available at {HBPLUS_PATH}"
        else:
            return False, f"HBPLUS at {HBPLUS_PATH} failed: {result.stderr.decode()}"
    except subprocess.TimeoutExpired:
        return False, "HBPLUS timed out"
    except Exception as e:
        return False, f"Error running HBPLUS: {e}"


def run_hbplus(
    pdb_content: str,
    timeout: int = 30,
) -> Dict[str, Any]:
    """
    Run HBPLUS on a PDB structure.

    Args:
        pdb_content: PDB file content as string
        timeout: Maximum execution time in seconds

    Returns:
        Dictionary with:
        - status: "completed" or "failed"
        - hbonds: List[HydrogenBond] of detected H-bonds
        - raw_output: Raw .hb2 file content
        - error: Error message if failed
    """
    result = {
        "status": "pending",
        "hbonds": [],
        "raw_output": "",
        "error": None,
    }

    available, msg = check_hbplus_available()
    if not available:
        result["status"] = "failed"
        result["error"] = msg
        return result

    # Create temporary directory and files
    with tempfile.TemporaryDirectory() as tmpdir:
        pdb_path = Path(tmpdir) / "input.pdb"
        hb2_path = Path(tmpdir) / "input.hb2"

        # Write PDB
        pdb_path.write_text(pdb_content)

        # Run HBPLUS
        try:
            proc = subprocess.run(
                [HBPLUS_PATH, "-R", str(pdb_path)],  # -R for reduced output
                cwd=tmpdir,
                capture_output=True,
                timeout=timeout,
            )

            # HBPLUS creates <basename>.hb2 in the same directory
            if hb2_path.exists():
                hb2_content = hb2_path.read_text()
                result["raw_output"] = hb2_content
                result["hbonds"] = parse_hb2_file(hb2_content)
                result["status"] = "completed"
            else:
                # Check for error
                stderr = proc.stderr.decode() if proc.stderr else ""
                stdout = proc.stdout.decode() if proc.stdout else ""
                result["status"] = "failed"
                result["error"] = f"HBPLUS did not produce output. stderr: {stderr}, stdout: {stdout}"

        except subprocess.TimeoutExpired:
            result["status"] = "failed"
            result["error"] = f"HBPLUS timed out after {timeout}s"
        except Exception as e:
            result["status"] = "failed"
            result["error"] = str(e)

    return result


def parse_hb2_file(hb2_content: str) -> List[HydrogenBond]:
    """
    Parse HBPLUS .hb2 output file.

    HB2 Format (columns are fixed-width):
    Columns:
    1-5: Donor chain, residue number
    6-9: Donor residue name
    10-13: Donor atom name
    15-19: Acceptor chain, residue number
    20-23: Acceptor residue name
    24-27: Acceptor atom name
    28-32: D-A distance
    33-37: Category (MM, MS, SS, HS, etc.)
    38-42: H-A distance
    43-47: D-H-A angle

    Args:
        hb2_content: Content of .hb2 file

    Returns:
        List of HydrogenBond objects
    """
    hbonds = []

    for line in hb2_content.split('\n'):
        # Skip header and empty lines
        if not line.strip() or line.startswith('-') or 'DONOR' in line:
            continue

        # HB2 format varies slightly, use regex for flexibility
        # Example line: "A0001-ALA N   A0002-ALA O   2.89MM 1.95180.5"
        match = re.match(
            r'([A-Z]?)(\d+)-(\w+)\s+(\w+)\s+'  # Donor
            r'([A-Z]?)(\d+)-(\w+)\s+(\w+)\s+'  # Acceptor
            r'(\d+\.\d+)(\w+)\s*'              # D-A distance and category
            r'(\d+\.\d+)?(\d+\.\d+)?',         # H-A distance and angle (optional)
            line
        )

        if match:
            try:
                hbond = HydrogenBond(
                    donor_chain=match.group(1) or "A",
                    donor_residue=int(match.group(2)),
                    donor_resname=match.group(3),
                    donor_atom=match.group(4),
                    acceptor_chain=match.group(5) or "A",
                    acceptor_residue=int(match.group(6)),
                    acceptor_resname=match.group(7),
                    acceptor_atom=match.group(8),
                    da_distance=float(match.group(9)),
                    category=match.group(10),
                    ha_distance=float(match.group(11)) if match.group(11) else 0.0,
                    dha_angle=float(match.group(12)) if match.group(12) else 0.0,
                )
                hbonds.append(hbond)
            except (ValueError, TypeError):
                # Skip malformed lines
                continue

    return hbonds


def get_ligand_hbonds(
    pdb_content: str,
    ligand_name: str = "UNL",
    ligand_atoms: Optional[List[str]] = None,
    use_hbplus: bool = True,
) -> Dict[str, Any]:
    """
    Get H-bonds involving a specific ligand.

    Args:
        pdb_content: PDB file content
        ligand_name: Residue name of ligand (default "UNL")
        ligand_atoms: Optional list of specific atoms to track (e.g., ["N5", "N6"])
        use_hbplus: Whether to use HBPLUS (falls back to distance if unavailable)

    Returns:
        Dictionary with:
        - total_hbonds: Total H-bonds to ligand
        - hbonds: List of H-bond details
        - by_atom: Dict mapping atom names to their H-bonds
        - method: "hbplus" or "distance_fallback"
    """
    result = {
        "total_hbonds": 0,
        "hbonds": [],
        "by_atom": {},
        "method": "unknown",
    }

    if use_hbplus:
        hbplus_result = run_hbplus(pdb_content)
        if hbplus_result["status"] == "completed":
            result["method"] = "hbplus"

            # Filter for ligand H-bonds
            ligand_hbonds = [
                hb for hb in hbplus_result["hbonds"]
                if hb.involves_ligand(ligand_name)
            ]

            result["total_hbonds"] = len(ligand_hbonds)

            # Convert to dicts
            for hb in ligand_hbonds:
                hb_dict = {
                    "donor_chain": hb.donor_chain,
                    "donor_residue": hb.donor_residue,
                    "donor_resname": hb.donor_resname,
                    "donor_atom": hb.donor_atom,
                    "acceptor_chain": hb.acceptor_chain,
                    "acceptor_residue": hb.acceptor_residue,
                    "acceptor_resname": hb.acceptor_resname,
                    "acceptor_atom": hb.acceptor_atom,
                    "da_distance": hb.da_distance,
                    "ha_distance": hb.ha_distance,
                    "dha_angle": hb.dha_angle,
                    "category": hb.category,
                }
                result["hbonds"].append(hb_dict)

            # Group by atom if requested
            if ligand_atoms:
                for atom in ligand_atoms:
                    atom_hbonds = [
                        hb for hb in ligand_hbonds
                        if hb.involves_atom(ligand_name, atom)
                    ]
                    result["by_atom"][atom] = [
                        {
                            "partner_chain": hb.donor_chain if hb.acceptor_atom == atom else hb.acceptor_chain,
                            "partner_residue": hb.donor_residue if hb.acceptor_atom == atom else hb.acceptor_residue,
                            "partner_resname": hb.donor_resname if hb.acceptor_atom == atom else hb.acceptor_resname,
                            "partner_atom": hb.donor_atom if hb.acceptor_atom == atom else hb.acceptor_atom,
                            "distance": hb.da_distance,
                            "angle": hb.dha_angle,
                        }
                        for hb in atom_hbonds
                    ]

            return result

    # Fallback to distance-based detection
    print("[HBPlus] HBPLUS unavailable, using distance-based fallback")
    from binding_analysis import analyze_ligand_hbonds as distance_analyze

    fallback_result = distance_analyze(
        pdb_content=pdb_content,
        ligand_name=ligand_name,
        ligand_atoms=ligand_atoms,
        hbond_distance=3.5,
    )

    fallback_result["method"] = "distance_fallback"
    return fallback_result


def compare_methods(
    pdb_content: str,
    ligand_name: str = "UNL",
    ligand_atoms: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Compare HBPLUS vs distance-based H-bond detection.

    Useful for validation and benchmarking.

    Returns:
        Dictionary comparing both methods' results
    """
    result = {
        "hbplus": None,
        "distance": None,
        "agreement": None,
    }

    # Run HBPLUS
    hbplus_result = get_ligand_hbonds(
        pdb_content=pdb_content,
        ligand_name=ligand_name,
        ligand_atoms=ligand_atoms,
        use_hbplus=True,
    )

    if hbplus_result["method"] == "hbplus":
        result["hbplus"] = {
            "total": hbplus_result["total_hbonds"],
            "by_atom": {k: len(v) for k, v in hbplus_result.get("by_atom", {}).items()},
        }

    # Run distance-based
    from binding_analysis import analyze_ligand_hbonds as distance_analyze
    dist_result = distance_analyze(
        pdb_content=pdb_content,
        ligand_name=ligand_name,
        ligand_atoms=ligand_atoms,
        hbond_distance=3.5,
    )

    result["distance"] = {
        "total": dist_result["total_hbonds"],
        "by_atom": {k: len(v) for k, v in dist_result.get("by_atom", {}).items()},
    }

    # Calculate agreement
    if result["hbplus"] and result["distance"]:
        hbplus_total = result["hbplus"]["total"]
        dist_total = result["distance"]["total"]
        if max(hbplus_total, dist_total) > 0:
            agreement = min(hbplus_total, dist_total) / max(hbplus_total, dist_total)
        else:
            agreement = 1.0  # Both found 0
        result["agreement"] = agreement

    return result


# Quick test if run directly
if __name__ == "__main__":
    available, msg = check_hbplus_available()
    print(f"HBPLUS available: {available}")
    print(f"Message: {msg}")

    # Example usage
    example_pdb = """HEADER    TEST
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.251   2.400   0.000  1.00  0.00           O
HETATM    5  N5  UNL X   1       0.000   3.000   0.000  1.00  0.00           N
HETATM    6  N6  UNL X   1       1.000   3.000   0.000  1.00  0.00           N
END
"""

    if available:
        result = get_ligand_hbonds(
            pdb_content=example_pdb,
            ligand_name="UNL",
            ligand_atoms=["N5", "N6"],
        )
        print(f"\nH-bond result: {result}")
    else:
        print("\nSkipping H-bond test (HBPLUS not available)")
