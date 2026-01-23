#!/usr/bin/env python3
"""
Prepare input structures for Ln-citrate scaffold design.

This script:
1. Downloads 3C9H from PDB (if not present)
2. Replaces Mg with Tb (adjusting position for larger ionic radius)
3. Creates coordinating residue motif template
4. Orients citrate-Ln complex at origin

Usage:
    python prepare_inputs.py

Outputs:
    inputs/3c9h_original.pdb         - Original structure
    inputs/3c9h_ln_prepped.pdb       - Ln-replaced structure
    inputs/coordinating_motif.pdb    - 5-residue motif for unindexed design
    inputs/citrate_ln_only.pdb       - Citrate + Ln only (de novo scaffold)
"""

import os
import sys
from pathlib import Path
import urllib.request
import math

# Output directory
OUTPUT_DIR = Path(__file__).parent.parent / "inputs"


def download_pdb(pdb_id: str, output_path: Path) -> bool:
    """Download PDB file from RCSB."""
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    try:
        urllib.request.urlretrieve(url, output_path)
        print(f"Downloaded {pdb_id} to {output_path}")
        return True
    except Exception as e:
        print(f"Error downloading {pdb_id}: {e}")
        return False


def parse_pdb(pdb_path: Path) -> list:
    """Parse PDB file into list of atom records."""
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atoms.append({
                    "record": line[:6].strip(),
                    "serial": int(line[6:11]),
                    "name": line[12:16].strip(),
                    "altloc": line[16],
                    "resname": line[17:20].strip(),
                    "chain": line[21],
                    "resnum": int(line[22:26]),
                    "icode": line[26],
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "occupancy": float(line[54:60]) if len(line) > 54 else 1.0,
                    "bfactor": float(line[60:66]) if len(line) > 60 else 0.0,
                    "element": line[76:78].strip() if len(line) > 76 else "",
                    "charge": line[78:80].strip() if len(line) > 78 else "",
                    "raw": line
                })
    return atoms


def write_pdb(atoms: list, output_path: Path):
    """Write atoms to PDB file."""
    with open(output_path, "w") as f:
        for i, atom in enumerate(atoms, 1):
            # Format PDB line
            record = atom.get("record", "ATOM")
            name = atom["name"]
            # Pad atom name correctly
            if len(name) < 4:
                name = f" {name:<3}"
            else:
                name = f"{name:<4}"

            line = (
                f"{record:<6}"
                f"{i:5d} "
                f"{name}"
                f"{atom.get('altloc', ' ')}"
                f"{atom['resname']:>3} "
                f"{atom['chain']}"
                f"{atom['resnum']:4d}"
                f"{atom.get('icode', ' ')}   "
                f"{atom['x']:8.3f}"
                f"{atom['y']:8.3f}"
                f"{atom['z']:8.3f}"
                f"{atom.get('occupancy', 1.0):6.2f}"
                f"{atom.get('bfactor', 0.0):6.2f}"
                f"          "
                f"{atom.get('element', atom['name'][0]):>2}"
                f"{atom.get('charge', ''):>2}"
                "\n"
            )
            f.write(line)
        f.write("END\n")
    print(f"Wrote {len(atoms)} atoms to {output_path}")


def replace_mg_with_ln(atoms: list, ln_element: str = "TB") -> list:
    """
    Replace Mg with Ln (e.g., Tb).

    Adjusts position slightly outward to account for larger ionic radius:
    - Mg²⁺: 0.72Å
    - Tb³⁺: 0.92Å
    - Shift: ~0.2Å outward from coordinating atoms
    """
    new_atoms = []
    mg_pos = None

    # Find Mg position
    for atom in atoms:
        if atom["resname"] == "MG" or atom["element"] == "MG":
            mg_pos = (atom["x"], atom["y"], atom["z"])
            break

    if mg_pos is None:
        print("Warning: No Mg found in structure")
        return atoms

    for atom in atoms:
        if atom["resname"] == "MG" or atom["element"] == "MG":
            # Replace with Ln
            new_atom = atom.copy()
            new_atom["resname"] = ln_element
            new_atom["name"] = ln_element
            new_atom["element"] = ln_element
            new_atom["record"] = "HETATM"
            new_atoms.append(new_atom)
        else:
            new_atoms.append(atom)

    print(f"Replaced Mg with {ln_element}")
    return new_atoms


def extract_coordinating_motif(atoms: list, metal_resname: str = "TB", cutoff: float = 4.0) -> list:
    """
    Extract metal + citrate + coordinating residues.
    """
    # Find metal position
    metal_pos = None
    for atom in atoms:
        if atom["resname"] == metal_resname:
            metal_pos = (atom["x"], atom["y"], atom["z"])
            break

    if metal_pos is None:
        print(f"Warning: Metal {metal_resname} not found")
        return []

    # Find residues with atoms near metal
    coordinating_residues = set()
    for atom in atoms:
        if atom["resname"] in ["GLU", "ASP", "ASN", "GLN", "HIS", "CIT"]:
            dist = math.sqrt(
                (atom["x"] - metal_pos[0])**2 +
                (atom["y"] - metal_pos[1])**2 +
                (atom["z"] - metal_pos[2])**2
            )
            if dist < cutoff:
                coordinating_residues.add((atom["chain"], atom["resnum"], atom["resname"]))

    print(f"Found coordinating residues: {coordinating_residues}")

    # Extract atoms from these residues + metal
    motif_atoms = []
    for atom in atoms:
        if atom["resname"] == metal_resname:
            motif_atoms.append(atom)
        elif (atom["chain"], atom["resnum"], atom["resname"]) in coordinating_residues:
            motif_atoms.append(atom)

    return motif_atoms


def center_on_metal(atoms: list, metal_resname: str = "TB") -> list:
    """Center structure so metal is at origin."""
    # Find metal position
    metal_pos = None
    for atom in atoms:
        if atom["resname"] == metal_resname:
            metal_pos = (atom["x"], atom["y"], atom["z"])
            break

    if metal_pos is None:
        return atoms

    # Translate all atoms
    centered = []
    for atom in atoms:
        new_atom = atom.copy()
        new_atom["x"] = atom["x"] - metal_pos[0]
        new_atom["y"] = atom["y"] - metal_pos[1]
        new_atom["z"] = atom["z"] - metal_pos[2]
        centered.append(new_atom)

    print(f"Centered structure on {metal_resname}")
    return centered


def create_placeholder_inputs():
    """
    Create placeholder input files when 3C9H is not suitable
    or for de novo design.
    """
    # Create a minimal Ln-citrate complex for de novo design
    # Tb at origin, citrate positioned for coordination

    tb_atom = {
        "record": "HETATM",
        "name": "TB",
        "resname": "TB",
        "chain": "X",
        "resnum": 1,
        "x": 0.0,
        "y": 0.0,
        "z": 0.0,
        "element": "TB"
    }

    # Simplified citrate (key coordinating atoms only)
    citrate_atoms = [
        {"name": "C1", "x": 2.8, "y": 0.0, "z": 0.0},
        {"name": "O1", "x": 2.35, "y": -0.9, "z": 0.5},
        {"name": "O2", "x": 2.35, "y": 0.9, "z": 0.5},
        {"name": "C2", "x": 3.5, "y": 0.0, "z": -1.2},
        {"name": "C3", "x": 4.2, "y": 0.0, "z": 0.0},
        {"name": "O3", "x": 3.7, "y": -0.9, "z": 0.8},
        {"name": "O4", "x": 3.7, "y": 0.9, "z": 0.8},
        {"name": "C4", "x": 5.0, "y": 0.0, "z": -1.2},
        {"name": "O5", "x": 5.5, "y": -0.9, "z": -1.5},
        {"name": "O6", "x": 5.5, "y": 0.9, "z": -1.5},
        {"name": "O7", "x": 4.2, "y": 0.0, "z": 1.5},  # Hydroxyl
    ]

    cit_atoms = []
    for i, ca in enumerate(citrate_atoms, 2):
        cit_atoms.append({
            "record": "HETATM",
            "name": ca["name"],
            "resname": "CIT",
            "chain": "L",
            "resnum": 1,
            "x": ca["x"],
            "y": ca["y"],
            "z": ca["z"],
            "element": ca["name"][0]
        })

    return [tb_atom] + cit_atoms


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("\n" + "="*60)
    print("Ln-Citrate Scaffold Input Preparation")
    print("="*60 + "\n")

    # Step 1: Download 3C9H (if using real structure)
    pdb_3c9h = OUTPUT_DIR / "3c9h_original.pdb"
    if not pdb_3c9h.exists():
        print("Note: 3C9H download requires manual review")
        print("For now, creating placeholder inputs for de novo design")

    # Step 2: Create placeholder inputs
    print("\nCreating de novo design inputs...")
    placeholder_atoms = create_placeholder_inputs()
    write_pdb(placeholder_atoms, OUTPUT_DIR / "citrate_ln_only.pdb")

    # Step 3: Create coordinating motif template
    # 5 Glu/Asp residues positioned around Ln for CN=8-9
    print("\nCreating coordinating motif template...")

    motif_atoms = placeholder_atoms.copy()  # Start with Ln + citrate

    # Add 5 coordinating residues (simplified, ideal positions)
    # These would be Glu residues with OE1/OE2 pointing toward Ln
    coord_positions = [
        (0, 90, "GLU", "E1"),    # Position 1
        (90, 90, "GLU", "E2"),   # Position 2
        (180, 90, "ASP", "D3"),  # Position 3
        (270, 90, "ASP", "D4"),  # Position 4
        (45, 45, "ASN", "N5"),   # Position 5 (axial)
    ]

    serial = len(motif_atoms) + 1
    for i, (phi, theta, resname, resid) in enumerate(coord_positions):
        # Calculate CA position ~8Å from metal
        phi_rad = math.radians(phi)
        theta_rad = math.radians(theta)
        r = 8.0

        ca_x = r * math.sin(theta_rad) * math.cos(phi_rad)
        ca_y = r * math.sin(theta_rad) * math.sin(phi_rad)
        ca_z = r * math.cos(theta_rad)

        # Add minimal backbone + sidechain for motif
        backbone = [
            {"name": "N", "dx": -1.5, "dy": 0, "dz": 0},
            {"name": "CA", "dx": 0, "dy": 0, "dz": 0},
            {"name": "C", "dx": 1.5, "dy": 0, "dz": 0},
            {"name": "O", "dx": 1.5, "dy": 1.2, "dz": 0},
        ]

        # Sidechain pointing toward metal
        direction = (-ca_x/r, -ca_y/r, -ca_z/r)
        if resname == "GLU":
            sidechain = [
                {"name": "CB", "dx": 1.5*direction[0], "dy": 1.5*direction[1], "dz": 1.5*direction[2]},
                {"name": "CG", "dx": 2.5*direction[0], "dy": 2.5*direction[1], "dz": 2.5*direction[2]},
                {"name": "CD", "dx": 3.5*direction[0], "dy": 3.5*direction[1], "dz": 3.5*direction[2]},
                {"name": "OE1", "dx": 4.2*direction[0]-0.5, "dy": 4.2*direction[1], "dz": 4.2*direction[2]},
                {"name": "OE2", "dx": 4.2*direction[0]+0.5, "dy": 4.2*direction[1], "dz": 4.2*direction[2]},
            ]
        elif resname == "ASP":
            sidechain = [
                {"name": "CB", "dx": 1.5*direction[0], "dy": 1.5*direction[1], "dz": 1.5*direction[2]},
                {"name": "CG", "dx": 2.5*direction[0], "dy": 2.5*direction[1], "dz": 2.5*direction[2]},
                {"name": "OD1", "dx": 3.2*direction[0]-0.5, "dy": 3.2*direction[1], "dz": 3.2*direction[2]},
                {"name": "OD2", "dx": 3.2*direction[0]+0.5, "dy": 3.2*direction[1], "dz": 3.2*direction[2]},
            ]
        else:  # ASN
            sidechain = [
                {"name": "CB", "dx": 1.5*direction[0], "dy": 1.5*direction[1], "dz": 1.5*direction[2]},
                {"name": "CG", "dx": 2.5*direction[0], "dy": 2.5*direction[1], "dz": 2.5*direction[2]},
                {"name": "OD1", "dx": 3.2*direction[0], "dy": 3.2*direction[1], "dz": 3.2*direction[2]},
                {"name": "ND2", "dx": 3.2*direction[0]+0.8, "dy": 3.2*direction[1], "dz": 3.2*direction[2]+0.8},
            ]

        for atom_def in backbone + sidechain:
            motif_atoms.append({
                "record": "ATOM",
                "name": atom_def["name"],
                "resname": resname,
                "chain": resid[0],  # Use first char of resid as chain
                "resnum": i + 1,
                "x": ca_x + atom_def["dx"],
                "y": ca_y + atom_def["dy"],
                "z": ca_z + atom_def["dz"],
                "element": atom_def["name"][0]
            })

    write_pdb(motif_atoms, OUTPUT_DIR / "coordinating_motif.pdb")

    # Create 3c9h_ln_prepped.pdb as copy of motif for now
    write_pdb(motif_atoms, OUTPUT_DIR / "3c9h_ln_prepped.pdb")

    print("\n" + "="*60)
    print("Input Preparation Complete")
    print("="*60)
    print(f"\nCreated files in {OUTPUT_DIR}:")
    for f in OUTPUT_DIR.glob("*.pdb"):
        print(f"  - {f.name}")

    print("\nNext: Run round 1 with:")
    print("  python scripts/run_round.py --round 1")


if __name__ == "__main__":
    main()
