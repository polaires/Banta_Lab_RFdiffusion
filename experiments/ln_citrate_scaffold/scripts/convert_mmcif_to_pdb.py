#!/usr/bin/env python3
"""
Convert RF3 mmCIF outputs to standard PDB format.
"""

import re
from pathlib import Path

EXPERIMENT_DIR = Path(__file__).parent.parent
RF3_DIR = EXPERIMENT_DIR / "outputs" / "round_05_rf3"


def parse_mmcif_to_atoms(filepath: Path) -> list:
    """Parse mmCIF file and extract atom records."""
    atoms = []
    in_atom_site = False
    header_fields = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()

            if line.startswith("_atom_site."):
                in_atom_site = True
                field = line.split(".")[1].strip()
                header_fields.append(field)
                continue

            if in_atom_site and (line.startswith("#") or line.startswith("_") or line.startswith("loop_")):
                in_atom_site = False
                continue

            if in_atom_site and (line.startswith("ATOM") or line.startswith("HETATM")):
                parts = line.split()
                if len(parts) < len(header_fields):
                    continue

                atom = {}
                for i, field in enumerate(header_fields):
                    if i < len(parts):
                        atom[field] = parts[i]

                atoms.append(atom)

    return atoms


def format_pdb_line(atom: dict, serial: int) -> str:
    """Format atom dictionary as PDB ATOM/HETATM line."""
    record = atom.get("group_PDB", "ATOM")
    name = atom.get("label_atom_id", "")
    res_name = atom.get("label_comp_id", "")
    chain = atom.get("label_asym_id", "A")
    res_seq = atom.get("label_seq_id", "1")

    try:
        x = float(atom.get("Cartn_x", 0))
        y = float(atom.get("Cartn_y", 0))
        z = float(atom.get("Cartn_z", 0))
        occupancy = float(atom.get("occupancy", 1.0))
        b_factor = float(atom.get("B_iso_or_equiv", 0))
    except ValueError:
        return None

    element = atom.get("type_symbol", name[0] if name else "X")

    # Format atom name (4 chars, left-justified if 4 chars, otherwise pad)
    if len(name) < 4:
        atom_name = f" {name:<3}"
    else:
        atom_name = name[:4]

    # Format residue name (3 chars, right-justified)
    res_name_fmt = f"{res_name:>3}"[:3]

    # Handle numeric chain IDs (convert to letters A-Z)
    if chain.isdigit():
        chain_letter = chr(ord('A') + int(chain) % 26)
    else:
        chain_letter = chain[0] if chain else "A"

    # Format residue sequence (4 chars, right-justified)
    try:
        res_seq_int = int(res_seq) if res_seq else 1
    except ValueError:
        res_seq_int = 1

    # PDB format:
    # ATOM   1234  N   ALA A  12      1.234   5.678   9.012  1.00 20.00           N
    line = (
        f"{record:<6}"       # 1-6: ATOM/HETATM
        f"{serial:>5} "      # 7-11: serial + space
        f"{atom_name}"       # 13-16: atom name
        f" "                 # 17: alt loc
        f"{res_name_fmt} "   # 18-20: res name + space
        f"{chain_letter}"    # 22: chain
        f"{res_seq_int:>4}"  # 23-26: res seq
        f"    "              # 27-30: insertion code + spaces
        f"{x:>8.3f}"         # 31-38: x
        f"{y:>8.3f}"         # 39-46: y
        f"{z:>8.3f}"         # 47-54: z
        f"{occupancy:>6.2f}" # 55-60: occupancy
        f"{b_factor:>6.2f}"  # 61-66: temp factor
        f"          "        # 67-76: spaces
        f"{element:>2}"      # 77-78: element
    )

    return line


def convert_mmcif_to_pdb(input_path: Path, output_path: Path):
    """Convert mmCIF file to PDB format."""
    atoms = parse_mmcif_to_atoms(input_path)

    if not atoms:
        print(f"  No atoms found in {input_path.name}")
        return False

    with open(output_path, 'w') as f:
        for i, atom in enumerate(atoms):
            line = format_pdb_line(atom, i)
            if line:
                f.write(line + "\n")
        f.write("END\n")

    return True


def main():
    print("=" * 60)
    print("Converting RF3 mmCIF to PDB Format")
    print("=" * 60)

    # Find all mmCIF files (named *.pdb but in mmCIF format)
    mmcif_files = sorted(RF3_DIR.glob("*_rf3.pdb"))

    if not mmcif_files:
        print(f"No RF3 files found in {RF3_DIR}")
        return

    print(f"\nFound {len(mmcif_files)} RF3 predictions\n")

    converted = 0
    for mmcif_file in mmcif_files:
        # Skip if already has a converted version
        output_name = mmcif_file.stem + "_converted.pdb"
        output_path = RF3_DIR / output_name

        if output_path.exists():
            print(f"  {mmcif_file.name} -> already converted")
            converted += 1
            continue

        print(f"  {mmcif_file.name} -> {output_name}... ", end="")
        if convert_mmcif_to_pdb(mmcif_file, output_path):
            print("OK")
            converted += 1
        else:
            print("FAILED")

    print(f"\nConverted: {converted}/{len(mmcif_files)} files")
    print(f"Output: {RF3_DIR}")


if __name__ == "__main__":
    main()
