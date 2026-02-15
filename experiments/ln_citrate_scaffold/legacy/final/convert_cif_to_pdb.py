"""Convert CIF-format .pdb files to proper PDB format.

Scans multiple directories for .pdb files that are actually mmCIF format
(start with 'data_') and converts them using biotite.

Usage (inside Docker container):
    python /experiments/ln_citrate_scaffold/final/convert_cif_to_pdb.py
"""
import io
import sys
from pathlib import Path

from biotite.structure.io.pdbx import CIFFile, get_structure
from biotite.structure.io.pdb import PDBFile


def convert_directory(directory: Path, recursive: bool = False) -> tuple:
    """Convert all CIF-in-PDB files in a directory. Returns (converted, skipped, failed)."""
    if not directory.exists():
        print(f"  Directory not found: {directory}")
        return 0, 0, 0

    pattern = "**/*.pdb" if recursive else "*.pdb"
    pdb_files = sorted(directory.glob(pattern))
    if not pdb_files:
        return 0, 0, 0

    converted = skipped = failed = 0

    for pdb_path in pdb_files:
        with open(pdb_path) as f:
            first_line = f.readline()

        # Skip files already in PDB format
        if not first_line.startswith("data_"):
            skipped += 1
            continue

        try:
            cif_file = CIFFile.read(str(pdb_path))
            block = list(cif_file.values())[0]
            atom_array = get_structure(block, model=1)

            pdb_file = PDBFile()
            pdb_file.set_structure(atom_array)
            buf = io.StringIO()
            pdb_file.write(buf)
            pdb_content = buf.getvalue()

            if not pdb_content.strip():
                print(f"  WARN: Empty PDB output for {pdb_path}")
                failed += 1
                continue

            # Back up original CIF
            cif_backup = pdb_path.with_suffix(".cif")
            if not cif_backup.exists():
                pdb_path.rename(cif_backup)
            else:
                # CIF backup already exists (re-run), just overwrite PDB
                pass

            with open(pdb_path, "w") as f:
                f.write(pdb_content)

            converted += 1

        except Exception as e:
            print(f"  FAIL: {pdb_path}: {e}")
            failed += 1

    return converted, skipped, failed


# Resolve base path (container vs WSL)
base = Path("/experiments/ln_citrate_scaffold/final")
if not base.exists():
    base = Path("/mnt/g/Github_local_repo/Banta_Lab_RFdiffusion/experiments/ln_citrate_scaffold/final")

dirs_to_convert = [
    ("promising/", base / "promising", False),
    ("top_candidates/ (recursive)", base / "top_candidates", True),
]

total_converted = total_skipped = total_failed = 0

for label, directory, recursive in dirs_to_convert:
    print(f"\n--- {label} ---")
    c, s, f = convert_directory(directory, recursive=recursive)
    print(f"  {c} converted, {s} already PDB, {f} failed")
    total_converted += c
    total_skipped += s
    total_failed += f

print(f"\n=== TOTAL: {total_converted} converted, {total_skipped} already PDB, {total_failed} failed ===")
