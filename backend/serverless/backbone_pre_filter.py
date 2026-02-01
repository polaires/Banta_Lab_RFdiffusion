"""
Backbone pre-filter for scout step.

CPU-only checks that run BEFORE the MPNN+RF3 GPU loop to catch obviously
bad backbones (broken chains, missing metal coordination, missing ligands)
without wasting GPU compute.

Generalizable: runs only applicable checks based on what's detected in the PDB.
Errors in any check = allow backbone through (non-fatal).
"""

from typing import Any, Dict, List, Optional
import math


def run_backbone_pre_filter(
    pdb_content: str,
    target_metal: Optional[str] = None,
    ligand_name: Optional[str] = None,
    thresholds: Optional[Dict[str, float]] = None,
) -> Dict[str, Any]:
    """
    Run applicable pre-filter checks on a backbone PDB.

    Args:
        pdb_content: PDB file content as string
        target_metal: Expected metal code (TB, CA, ZN, etc.) — enables metal check
        ligand_name: Expected ligand residue name (CIT, PQQ, etc.) — enables ligand check
        thresholds: Override default thresholds (ca_gap_max, cn_min, geom_rmsd_max, etc.)

    Returns:
        {
            "passed": bool,
            "checks": {"chain_continuity": {...}, "metal_coordination": {...}, ...},
            "failed_checks": ["chain_continuity"],
            "skipped_checks": ["ligand_presence"],
        }
    """
    defaults = {
        "ca_gap_max": 4.5,       # Max CA-CA distance before flagging chain break
        "cn_min": 4,             # Minimum coordination number for metal
        "geom_rmsd_max": 1.5,   # Max geometry RMSD for metal coordination
        "oo_clash_max": 2,       # Max inter-residue O-O clashes allowed
        "ligand_min_atoms": 3,   # Min atoms for ligand to be considered present
        "pocket_ca_min": 2,      # Min CA atoms within pocket_radius of metal
        "pocket_radius": 8.0,    # Radius (A) for pocket density check
        "ligand_metal_max_dist": 5.0,  # Max ligand-metal distance (A)
    }

    # Chemistry-aware defaults: adjust thresholds based on metal type
    if target_metal:
        try:
            from metal_chemistry import METAL_DATABASE, get_coordination_number_range
            metal_upper = target_metal.upper()
            if metal_upper in METAL_DATABASE:
                metal_data = METAL_DATABASE[metal_upper]
                default_ox = metal_data.get("default_oxidation", 2)
                cn_min_formal, cn_max_formal = get_coordination_number_range(
                    metal_upper, default_ox
                )
                # Pre-filter is lenient: formal_cn_min - 2, floor at 3
                defaults["cn_min"] = max(3, cn_min_formal - 2)
        except (ImportError, ValueError, KeyError):
            pass  # Use hardcoded defaults

    if thresholds:
        defaults.update(thresholds)

    checks: Dict[str, Dict[str, Any]] = {}
    failed_checks: List[str] = []
    skipped_checks: List[str] = []

    # --- Check 1: Chain Continuity (always runs) ---
    try:
        chain_result = _check_chain_continuity(pdb_content, defaults["ca_gap_max"])
        checks["chain_continuity"] = chain_result
        if not chain_result["passed"]:
            failed_checks.append("chain_continuity")
    except Exception as e:
        # Non-fatal: allow through
        checks["chain_continuity"] = {"passed": True, "error": str(e)}

    # --- Detect metal from PDB (for Check 2) ---
    detected_metal = _detect_metal_hetatm(pdb_content)
    has_metal = detected_metal is not None

    # --- Check 2: Metal Coordination (if metal detected or expected) ---
    if has_metal:
        try:
            # If both metal AND ligand, use complex site validation
            detected_ligand = _detect_ligand_hetatm(pdb_content) if ligand_name else None
            if detected_ligand and ligand_name:
                metal_result = _check_metal_ligand_complex(
                    pdb_content,
                    detected_metal["res_name"],
                    ligand_name,
                )
                checks["metal_coordination"] = metal_result
            else:
                metal_result = _check_metal_coordination(
                    pdb_content,
                    defaults["cn_min"],
                    defaults["geom_rmsd_max"],
                    defaults["oo_clash_max"],
                )
                checks["metal_coordination"] = metal_result

            if not metal_result["passed"]:
                failed_checks.append("metal_coordination")
        except Exception as e:
            checks["metal_coordination"] = {"passed": True, "error": str(e)}
    else:
        skipped_checks.append("metal_coordination")

    # --- Check 3: Ligand Presence (if ligand_name provided) ---
    if ligand_name:
        try:
            ligand_result = _check_ligand_presence(
                pdb_content, ligand_name, int(defaults["ligand_min_atoms"])
            )
            checks["ligand_presence"] = ligand_result
            if not ligand_result["passed"]:
                failed_checks.append("ligand_presence")
        except Exception as e:
            checks["ligand_presence"] = {"passed": True, "error": str(e)}
    else:
        skipped_checks.append("ligand_presence")

    # --- Check 4: Coordination Pocket Density (if metal detected) ---
    # At backbone stage, sidechains are poly-Ala — only CA positions are meaningful.
    # Count CA atoms near the metal to assess whether MPNN has enough residues
    # to place coordinating donors (Glu/Asp/His).
    if has_metal and detected_metal:
        try:
            pocket_result = _check_coordination_pocket(
                pdb_content,
                detected_metal,
                int(defaults["pocket_ca_min"]),
                float(defaults["pocket_radius"]),
            )
            checks["coordination_pocket"] = pocket_result
            if not pocket_result["passed"]:
                failed_checks.append("coordination_pocket")
        except Exception as e:
            checks["coordination_pocket"] = {"passed": True, "error": str(e)}
    else:
        skipped_checks.append("coordination_pocket")

    # --- Check 5: Ligand-Metal Distance (if both metal and ligand present) ---
    # Verifies the pre-formed metal-ligand complex is intact (both are HETATM
    # records placed by RFD3, so positions are valid at backbone stage).
    if has_metal and detected_metal and ligand_name:
        try:
            lig_dist_result = _check_ligand_metal_distance(
                pdb_content,
                detected_metal,
                ligand_name,
                float(defaults["ligand_metal_max_dist"]),
            )
            checks["ligand_metal_distance"] = lig_dist_result
            if not lig_dist_result["passed"]:
                failed_checks.append("ligand_metal_distance")
        except Exception as e:
            checks["ligand_metal_distance"] = {"passed": True, "error": str(e)}
    else:
        skipped_checks.append("ligand_metal_distance")

    return {
        "passed": len(failed_checks) == 0,
        "checks": checks,
        "failed_checks": failed_checks,
        "skipped_checks": skipped_checks,
    }


# ============== Check Implementations ==============


def _check_chain_continuity(
    pdb_content: str, ca_gap_max: float
) -> Dict[str, Any]:
    """
    Parse CA atoms, check that consecutive residues have CA-CA distance ~3.8A.
    Fail if any gap > ca_gap_max.
    """
    # Parse CA atoms by chain
    ca_atoms: Dict[str, List[tuple]] = {}  # chain -> [(resnum, x, y, z)]

    for line in pdb_content.split("\n"):
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue

        try:
            chain = line[21].strip() or "A"
            resnum = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            if chain not in ca_atoms:
                ca_atoms[chain] = []
            ca_atoms[chain].append((resnum, x, y, z))
        except (ValueError, IndexError):
            continue

    if not ca_atoms:
        return {"passed": True, "total_residues": 0, "chain_breaks": [],
                "note": "No CA atoms found"}

    # Sort by residue number within each chain
    chain_breaks: List[Dict[str, Any]] = []
    total_residues = 0

    for chain, atoms in ca_atoms.items():
        atoms.sort(key=lambda a: a[0])
        total_residues += len(atoms)

        for i in range(len(atoms) - 1):
            resnum_i, xi, yi, zi = atoms[i]
            resnum_j, xj, yj, zj = atoms[i + 1]

            # Only check consecutive residue numbers
            if resnum_j - resnum_i != 1:
                continue

            dist = math.sqrt(
                (xj - xi) ** 2 + (yj - yi) ** 2 + (zj - zi) ** 2
            )

            if dist > ca_gap_max:
                chain_breaks.append({
                    "chain": chain,
                    "res_i": resnum_i,
                    "res_j": resnum_j,
                    "distance": round(dist, 2),
                })

    return {
        "passed": len(chain_breaks) == 0,
        "chain_breaks": chain_breaks,
        "total_residues": total_residues,
    }


def _check_metal_coordination(
    pdb_content: str,
    cn_min: int,
    geom_rmsd_max: float,
    oo_clash_max: int,
) -> Dict[str, Any]:
    """
    Validate metal coordination geometry using existing validate_coordination_geometry.
    """
    try:
        from geometry_validation import validate_coordination_geometry
    except ImportError:
        return {"passed": True, "error": "geometry_validation not available"}

    result = validate_coordination_geometry(pdb_content)

    if not result.get("valid") and result.get("error"):
        # No metal found — not a failure, just skip
        if "No lanthanide" in str(result.get("error", "")):
            return {"passed": True, "note": "No lanthanide detected (non-lanthanide metal)"}
        return {"passed": True, "error": result.get("error")}

    cn = result.get("coordination_number", 0)
    rmsd = result.get("geometry_rmsd")
    clashes = result.get("inter_residue_clashes", result.get("clash_count", 0))

    reasons = []
    if cn < cn_min:
        reasons.append(f"CN={cn} < {cn_min}")
    if rmsd is not None and rmsd > geom_rmsd_max:
        reasons.append(f"geometry RMSD={rmsd:.2f} > {geom_rmsd_max}")
    if isinstance(clashes, int) and clashes > oo_clash_max:
        reasons.append(f"O-O clashes={clashes} > {oo_clash_max}")

    return {
        "passed": len(reasons) == 0,
        "coordination_number": cn,
        "geometry_rmsd": round(rmsd, 3) if rmsd is not None else None,
        "clash_count": clashes,
        "reasons": reasons,
    }


def _check_metal_ligand_complex(
    pdb_content: str,
    metal: str,
    ligand_name: str,
) -> Dict[str, Any]:
    """
    Validate metal-ligand complex site using existing validate_metal_ligand_complex_site.
    Uses relaxed thresholds for pre-filter (just catching obviously broken sites).
    """
    try:
        from metal_validation import validate_metal_ligand_complex_site
    except ImportError:
        return {"passed": True, "error": "metal_validation not available"}

    result = validate_metal_ligand_complex_site(
        pdb_content, metal, ligand_name
    )

    if not result.get("success", True):
        error = result.get("error", "")
        # Metal not found = definite failure
        if "not found" in error.lower():
            return {
                "passed": False,
                "reasons": [error],
                "coordination_number": 0,
            }
        return {"passed": True, "error": error}

    cn = result.get("coordination_number", 0)
    ligand_cn = result.get("ligand_coordination", 0)

    reasons = []
    # Very lenient for pre-filter: just catch total absence
    if cn < 3:
        reasons.append(f"Total CN={cn} < 3 (metal barely coordinated)")
    if ligand_cn == 0:
        reasons.append("Ligand provides 0 donors to metal")

    return {
        "passed": len(reasons) == 0,
        "coordination_number": cn,
        "ligand_coordination": ligand_cn,
        "protein_coordination": result.get("protein_coordination", 0),
        "reasons": reasons,
    }


def _check_ligand_presence(
    pdb_content: str, ligand_name: str, min_atoms: int
) -> Dict[str, Any]:
    """
    Verify expected ligand HETATM exists with sufficient atoms.
    Uses detect_ligand_from_pdb if available, falls back to simple HETATM scan.
    """
    try:
        from analysis_types import detect_ligand_from_pdb
        detected = detect_ligand_from_pdb(pdb_content)
        if detected and detected.ligand_name.upper() == ligand_name.upper():
            if detected.atom_count >= min_atoms:
                return {
                    "passed": True,
                    "ligand_name": detected.ligand_name,
                    "atom_count": detected.atom_count,
                }
            else:
                return {
                    "passed": False,
                    "ligand_name": detected.ligand_name,
                    "atom_count": detected.atom_count,
                    "reasons": [
                        f"Ligand {ligand_name} has only {detected.atom_count} atoms "
                        f"(need >= {min_atoms})"
                    ],
                }
        # Ligand detected but different name — check if expected one exists at all
        # Fall through to manual scan
    except ImportError:
        pass

    # Fallback: manual HETATM scan for the specific ligand
    atom_count = 0
    for line in pdb_content.split("\n"):
        if not line.startswith("HETATM"):
            continue
        res_name = line[17:20].strip().upper()
        if res_name == ligand_name.upper():
            atom_count += 1

    if atom_count == 0:
        return {
            "passed": False,
            "ligand_name": ligand_name,
            "atom_count": 0,
            "reasons": [f"Ligand {ligand_name} not found in PDB"],
        }

    if atom_count < min_atoms:
        return {
            "passed": False,
            "ligand_name": ligand_name,
            "atom_count": atom_count,
            "reasons": [
                f"Ligand {ligand_name} has only {atom_count} atoms "
                f"(need >= {min_atoms})"
            ],
        }

    return {
        "passed": True,
        "ligand_name": ligand_name,
        "atom_count": atom_count,
    }


def _check_coordination_pocket(
    pdb_content: str,
    detected_metal: Dict[str, Any],
    pocket_ca_min: int,
    pocket_radius: float,
) -> Dict[str, Any]:
    """
    Count backbone CA atoms within pocket_radius of the metal center.

    At the backbone/pre-filter stage, sidechains are poly-alanine placeholders.
    This check verifies the backbone forms a pocket around the metal where
    MPNN can later place coordinating residues (each nearby CA position can
    contribute ~2 sidechain donor atoms from Glu/Asp).

    Empirically validated: passing TB designs average 8.6 CAs within 8.0A.
    Threshold of 2 catches only completely surface-exposed metals.
    """
    mx, my, mz = detected_metal["x"], detected_metal["y"], detected_metal["z"]

    ca_near = 0
    ca_total = 0
    for line in pdb_content.split("\n"):
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        ca_total += 1
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except (ValueError, IndexError):
            continue
        dist = math.sqrt((x - mx) ** 2 + (y - my) ** 2 + (z - mz) ** 2)
        if dist <= pocket_radius:
            ca_near += 1

    passed = ca_near >= pocket_ca_min
    reasons = []
    if not passed:
        reasons.append(
            f"Only {ca_near} CA atoms within {pocket_radius}A of metal "
            f"(need >= {pocket_ca_min}) — metal too exposed for coordination shell"
        )

    return {
        "passed": passed,
        "ca_within_radius": ca_near,
        "pocket_radius": pocket_radius,
        "ca_total": ca_total,
        "reasons": reasons,
    }


def _check_ligand_metal_distance(
    pdb_content: str,
    detected_metal: Dict[str, Any],
    ligand_name: str,
    max_distance: float,
) -> Dict[str, Any]:
    """
    Verify the closest ligand atom is within max_distance of the metal center.

    Both metal and ligand are HETATM records placed by RFD3, so their positions
    are valid at the backbone stage. This catches cases where the metal-ligand
    complex has dissociated (e.g., ligand drifted away during diffusion).

    Empirically: in passing designs, ligand O atoms are within 0.86-2.74A of TB.
    Threshold of 5.0A is very generous, only catching truly dissociated complexes.
    """
    mx, my, mz = detected_metal["x"], detected_metal["y"], detected_metal["z"]
    ligand_upper = ligand_name.upper()

    min_dist = float("inf")
    ligand_atom_count = 0

    for line in pdb_content.split("\n"):
        if not line.startswith("HETATM"):
            continue
        res_name = line[17:20].strip().upper()
        if res_name != ligand_upper:
            continue
        ligand_atom_count += 1
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except (ValueError, IndexError):
            continue
        dist = math.sqrt((x - mx) ** 2 + (y - my) ** 2 + (z - mz) ** 2)
        if dist < min_dist:
            min_dist = dist

    if ligand_atom_count == 0:
        # No ligand found — don't fail here (ligand_presence check handles that)
        return {"passed": True, "note": f"No {ligand_name} atoms found for distance check"}

    passed = min_dist <= max_distance
    reasons = []
    if not passed:
        reasons.append(
            f"Nearest ligand atom is {min_dist:.1f}A from metal "
            f"(max {max_distance}A) — complex may be dissociated"
        )

    return {
        "passed": passed,
        "min_distance": round(min_dist, 2),
        "max_allowed": max_distance,
        "ligand_atoms_found": ligand_atom_count,
        "reasons": reasons,
    }


# ============== PDB Helpers ==============


def _detect_metal_hetatm(pdb_content: str) -> Optional[Dict[str, Any]]:
    """Quick scan for any metal HETATM record."""
    METAL_CODES = {
        "TB", "CA", "ZN", "MG", "FE", "MN", "CO", "NI", "CU",
        "CD", "HG", "PB", "EU", "GD", "DY", "HO", "ER", "YB",
        "LU", "LA", "CE", "PR", "ND", "SM", "TM",
    }

    for line in pdb_content.split("\n"):
        if not line.startswith("HETATM"):
            continue
        res_name = line[17:20].strip().upper()
        if res_name in METAL_CODES:
            try:
                return {
                    "res_name": res_name,
                    "chain": line[21].strip() or "A",
                    "resnum": int(line[22:26].strip()),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                }
            except (ValueError, IndexError):
                continue
    return None


def _detect_ligand_hetatm(pdb_content: str) -> Optional[Dict[str, Any]]:
    """Quick scan for any non-metal, non-water HETATM residue."""
    METAL_CODES = {
        "TB", "CA", "ZN", "MG", "FE", "MN", "CO", "NI", "CU",
        "CD", "HG", "PB", "EU", "GD", "DY", "HO", "ER", "YB",
        "LU", "LA", "CE", "PR", "ND", "SM", "TM",
    }
    SOLVENT = {"HOH", "WAT", "DOD", "SOL"}

    ligand_atoms: Dict[str, int] = {}

    for line in pdb_content.split("\n"):
        if not line.startswith("HETATM"):
            continue
        res_name = line[17:20].strip().upper()
        if res_name in METAL_CODES or res_name in SOLVENT:
            continue
        ligand_atoms[res_name] = ligand_atoms.get(res_name, 0) + 1

    if not ligand_atoms:
        return None

    # Return the one with most atoms
    best = max(ligand_atoms.items(), key=lambda x: x[1])
    return {"ligand_name": best[0], "atom_count": best[1]}
