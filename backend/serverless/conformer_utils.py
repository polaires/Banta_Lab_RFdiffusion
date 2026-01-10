"""
Multi-method conformer generation utilities.

Supports three methods for generating 3D conformers from SMILES:
1. RDKit ETKDG + MMFF (fast, default)
2. xTB/GFN2 (high-quality quantum-based)
3. Torsional Diffusion (ML state-of-the-art)
"""

from typing import Optional, Tuple
from enum import Enum
import tempfile
import subprocess
import os
import shutil


class ConformerMethod(str, Enum):
    """Available conformer generation methods."""
    RDKIT = "rdkit"          # Fast, default
    XTB = "xtb"              # High-quality quantum
    TORSIONAL = "torsional"  # ML state-of-the-art


def generate_conformer(
    smiles: str,
    method: ConformerMethod = ConformerMethod.RDKIT,
    name: str = "LIG",
    center: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    optimize: bool = True,
    fallback: bool = True,
) -> Optional[str]:
    """
    Generate 3D conformer from SMILES using specified method.

    Args:
        smiles: SMILES string for the molecule
        method: Conformer generation method (rdkit, xtb, or torsional)
        name: 3-letter residue name for PDB output (default: "LIG")
        center: (x, y, z) coordinates for ligand center
        optimize: Whether to optimize the geometry (default: True)
        fallback: If True, fall back to RDKit if requested method fails

    Returns:
        PDB string with HETATM records, or None if all methods fail
    """
    methods_to_try = [method]
    if fallback and method != ConformerMethod.RDKIT:
        methods_to_try.append(ConformerMethod.RDKIT)

    for m in methods_to_try:
        try:
            if m == ConformerMethod.RDKIT:
                result = _generate_rdkit(smiles, name, center, optimize)
            elif m == ConformerMethod.XTB:
                result = _generate_xtb(smiles, name, center)
            elif m == ConformerMethod.TORSIONAL:
                result = _generate_torsional_diffusion(smiles, name, center)
            else:
                raise ValueError(f"Unknown method: {m}")

            if result:
                print(f"[Conformer] Generated 3D structure using {m.value} method")
                return result

        except Exception as e:
            print(f"[Conformer] {m.value} method failed: {e}")
            if not fallback or m == methods_to_try[-1]:
                raise

    return None


def _generate_rdkit(
    smiles: str,
    name: str,
    center: Tuple[float, float, float],
    optimize: bool = True
) -> Optional[str]:
    """
    RDKit ETKDG + MMFF optimization.

    Fast baseline method using distance geometry embedding
    followed by MMFF94 force field optimization.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        import numpy as np

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Failed to parse SMILES: {smiles}")

        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)

        # Use ETKDGv3 for better conformer generation
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        result = AllChem.EmbedMolecule(mol, params)

        if result == -1:
            # Fallback to basic embedding if ETKDGv3 fails
            AllChem.EmbedMolecule(mol, randomSeed=42)

        if optimize:
            # MMFF optimization
            try:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            except Exception:
                # If MMFF fails, try UFF
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)

        # Get conformer
        conf = mol.GetConformer()

        # Calculate current center
        coords = []
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            coords.append([pos.x, pos.y, pos.z])
        coords = np.array(coords)
        current_center = np.mean(coords, axis=0)

        # Calculate translation to target center
        translation = np.array(center) - current_center

        # Count heavy atoms
        heavy_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
        print(f"[Conformer] RDKit: {mol.GetNumAtoms()} atoms, {heavy_atoms} heavy atoms")

        # Generate PDB content
        return _coords_to_pdb(mol, conf, translation, name)

    except ImportError:
        raise ImportError("RDKit not available. Install with: pip install rdkit")


def _generate_xtb(
    smiles: str,
    name: str,
    center: Tuple[float, float, float],
    timeout: int = 60
) -> Optional[str]:
    """
    xTB GFN2 optimization for high-quality conformers.

    Uses semi-empirical quantum chemistry for accurate geometries,
    especially important for photoswitches and metal complexes.
    """
    # Check if xtb is available
    xtb_path = shutil.which("xtb")
    if not xtb_path:
        raise RuntimeError("xTB not found. Install with: conda install -c conda-forge xtb")

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        import numpy as np

        # Generate initial structure with RDKit
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Failed to parse SMILES: {smiles}")

        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)

        # Get conformer
        conf = mol.GetConformer()

        # Create temp directory for xTB
        with tempfile.TemporaryDirectory() as tmpdir:
            xyz_input = os.path.join(tmpdir, "input.xyz")
            xyz_output = os.path.join(tmpdir, "xtbopt.xyz")

            # Write XYZ file
            _write_xyz(mol, conf, xyz_input)

            # Run xTB optimization
            print("[Conformer] Running GFN2-xTB optimization...")
            result = subprocess.run(
                ["xtb", "input.xyz", "--opt", "tight", "--gfn", "2"],
                cwd=tmpdir,
                capture_output=True,
                text=True,
                timeout=timeout
            )

            if result.returncode != 0:
                print(f"[Conformer] xTB warning: {result.stderr[:200]}")

            # Read optimized coordinates
            if os.path.exists(xyz_output):
                optimized_coords = _read_xyz(xyz_output)
            else:
                # xTB might output to different file
                xyz_output = os.path.join(tmpdir, "xtbopt.coord")
                if os.path.exists(xyz_output):
                    optimized_coords = _read_coord(xyz_output)
                else:
                    raise RuntimeError("xTB did not produce output file")

            # Update conformer with optimized coordinates
            for i, (x, y, z) in enumerate(optimized_coords):
                conf.SetAtomPosition(i, (x, y, z))

            # Calculate translation to center
            coords = np.array(optimized_coords)
            current_center = np.mean(coords, axis=0)
            translation = np.array(center) - current_center

            # Count heavy atoms
            heavy_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
            print(f"[Conformer] xTB: {mol.GetNumAtoms()} atoms, {heavy_atoms} heavy atoms")

            return _coords_to_pdb(mol, conf, translation, name)

    except subprocess.TimeoutExpired:
        raise RuntimeError(f"xTB optimization timed out after {timeout}s")


def _generate_torsional_diffusion(
    smiles: str,
    name: str,
    center: Tuple[float, float, float],
    n_conformers: int = 1
) -> Optional[str]:
    """
    Torsional Diffusion ML-based conformer generation.

    State-of-the-art ML method for drug-like molecules.
    Best for novel molecules with complex torsional profiles.
    """
    try:
        # Lazy import to avoid loading model until needed
        from torsional_diffusion import TorsionalDiffusion
        import numpy as np

        print("[Conformer] Loading Torsional Diffusion model...")
        model = TorsionalDiffusion.from_pretrained("drugs")

        print("[Conformer] Generating conformers with Torsional Diffusion...")
        conformers = model.generate(smiles, n_conformers=n_conformers)

        if not conformers:
            raise RuntimeError("Torsional Diffusion generated no conformers")

        # Get best conformer (first one, ranked by likelihood)
        best = conformers[0]

        # Convert to RDKit mol for PDB output
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)

        # Set coordinates from Torsional Diffusion output
        conf = Chem.Conformer(mol.GetNumAtoms())
        for i, (x, y, z) in enumerate(best.coordinates):
            conf.SetAtomPosition(i, (x, y, z))
        mol.AddConformer(conf, assignId=True)

        # Calculate translation to center
        coords = np.array(best.coordinates)
        current_center = np.mean(coords, axis=0)
        translation = np.array(center) - current_center

        return _coords_to_pdb(mol, conf, translation, name)

    except ImportError:
        raise ImportError(
            "Torsional Diffusion not available. "
            "Install with: pip install torsional-diffusion"
        )


def _coords_to_pdb(mol, conf, translation, name: str) -> str:
    """Convert RDKit molecule with conformer to PDB HETATM records."""
    import numpy as np

    pdb_lines = []
    atom_serial = 1

    for i, atom in enumerate(mol.GetAtoms()):
        pos = conf.GetAtomPosition(i)
        coord = np.array([pos.x, pos.y, pos.z]) + translation
        symbol = atom.GetSymbol()
        atom_name = f"{symbol}{i+1}"[:4]

        line = (
            f"HETATM{atom_serial:5d}  {atom_name:<3s} {name:3s} L   1    "
            f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}  1.00  0.00          {symbol:>2s}"
        )
        pdb_lines.append(line)
        atom_serial += 1

    pdb_lines.append("END")
    return '\n'.join(pdb_lines)


def _mol_to_sdf(mol, conf, translation, name: str) -> str:
    """Convert RDKit molecule with conformer to SDF format.

    SDF is preferred over PDB for RFD3 ligand input as it preserves
    full bond information and avoids CCD naming conflicts.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import numpy as np

    # Apply translation to conformer
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        new_pos = (pos.x + translation[0], pos.y + translation[1], pos.z + translation[2])
        conf.SetAtomPosition(i, new_pos)

    # Set molecule name
    mol.SetProp("_Name", name)

    # Generate SDF string
    sdf_string = Chem.MolToMolBlock(mol)
    sdf_string += "$$$$\n"

    return sdf_string


def generate_conformer_sdf(
    smiles: str,
    method: ConformerMethod = ConformerMethod.RDKIT,
    name: str = "L:LIG",  # Use L: prefix to avoid CCD conflicts
    center: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    optimize: bool = True,
    fallback: bool = True,
) -> Optional[str]:
    """
    Generate 3D conformer from SMILES and return as SDF string.

    This is the PREFERRED method for RFD3 ligand input.
    SDF preserves bond information and the L: prefix avoids CCD conflicts.

    Args:
        smiles: SMILES string for the molecule
        method: Conformer generation method (rdkit, xtb, or torsional)
        name: Molecule name (use L: prefix for custom ligands)
        center: (x, y, z) coordinates for ligand center
        optimize: Whether to optimize the geometry (default: True)
        fallback: If True, fall back to RDKit if requested method fails

    Returns:
        SDF string, or None if all methods fail
    """
    methods_to_try = [method]
    if fallback and method != ConformerMethod.RDKIT:
        methods_to_try.append(ConformerMethod.RDKIT)

    for m in methods_to_try:
        try:
            mol, conf, translation = _generate_mol_and_conformer(smiles, m, center, optimize)
            if mol and conf:
                print(f"[Conformer] Generated 3D structure using {m.value} method")
                heavy_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
                print(f"[Conformer] {mol.GetNumAtoms()} atoms, {heavy_atoms} heavy atoms")
                return _mol_to_sdf(mol, conf, translation, name)

        except Exception as e:
            print(f"[Conformer] {m.value} method failed: {e}")
            if not fallback or m == methods_to_try[-1]:
                raise

    return None


def _generate_mol_and_conformer(
    smiles: str,
    method: ConformerMethod,
    center: Tuple[float, float, float],
    optimize: bool
) -> tuple:
    """Generate molecule and conformer, return (mol, conf, translation)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import numpy as np

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Failed to parse SMILES: {smiles}")

    mol = Chem.AddHs(mol)

    if method == ConformerMethod.RDKIT:
        # Use ETKDGv3 for better conformer generation
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        result = AllChem.EmbedMolecule(mol, params)
        if result == -1:
            AllChem.EmbedMolecule(mol, randomSeed=42)
        if optimize:
            try:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            except Exception:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    elif method == ConformerMethod.XTB:
        # xTB optimization - first embed with RDKit, then optimize
        AllChem.EmbedMolecule(mol, randomSeed=42)
        mol = _optimize_with_xtb(mol)
    elif method == ConformerMethod.TORSIONAL:
        raise NotImplementedError("Torsional Diffusion not yet integrated for SDF output")

    conf = mol.GetConformer()

    # Calculate translation to center
    coords = []
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        coords.append([pos.x, pos.y, pos.z])
    coords = np.array(coords)
    current_center = np.mean(coords, axis=0)
    translation = np.array(center) - current_center

    return mol, conf, translation


def _optimize_with_xtb(mol):
    """Optimize molecule geometry with xTB."""
    import subprocess
    import os

    xtb_path = shutil.which("xtb")
    if not xtb_path:
        print("[Conformer] xTB not found, falling back to RDKit")
        from rdkit.Chem import AllChem
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        return mol

    conf = mol.GetConformer()

    with tempfile.TemporaryDirectory() as tmpdir:
        xyz_input = os.path.join(tmpdir, "input.xyz")
        xyz_output = os.path.join(tmpdir, "xtbopt.xyz")

        _write_xyz(mol, conf, xyz_input)

        print("[Conformer] Running GFN2-xTB optimization...")
        result = subprocess.run(
            ["xtb", "input.xyz", "--opt", "tight", "--gfn", "2"],
            cwd=tmpdir,
            capture_output=True,
            text=True,
            timeout=60
        )

        if os.path.exists(xyz_output):
            optimized_coords = _read_xyz(xyz_output)
            for i, (x, y, z) in enumerate(optimized_coords):
                conf.SetAtomPosition(i, (x, y, z))

    return mol


def _write_xyz(mol, conf, filepath: str):
    """Write RDKit molecule to XYZ format for xTB."""
    with open(filepath, 'w') as f:
        n_atoms = mol.GetNumAtoms()
        f.write(f"{n_atoms}\n")
        f.write("Generated by conformer_utils\n")

        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            symbol = atom.GetSymbol()
            f.write(f"{symbol:2s} {pos.x:15.8f} {pos.y:15.8f} {pos.z:15.8f}\n")


def _read_xyz(filepath: str) -> list:
    """Read XYZ file and return list of (x, y, z) coordinates."""
    coords = []
    with open(filepath, 'r') as f:
        lines = f.readlines()
        n_atoms = int(lines[0].strip())
        for i in range(2, 2 + n_atoms):
            parts = lines[i].split()
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            coords.append((x, y, z))
    return coords


def _read_coord(filepath: str) -> list:
    """Read xTB coord file format and return coordinates in Angstroms."""
    coords = []
    bohr_to_angstrom = 0.529177

    with open(filepath, 'r') as f:
        lines = f.readlines()

    in_coords = False
    for line in lines:
        line = line.strip()
        if line.startswith('$coord'):
            in_coords = True
            continue
        if line.startswith('$'):
            in_coords = False
            continue
        if in_coords and line:
            parts = line.split()
            # Coord file uses Bohr units
            x = float(parts[0]) * bohr_to_angstrom
            y = float(parts[1]) * bohr_to_angstrom
            z = float(parts[2]) * bohr_to_angstrom
            coords.append((x, y, z))

    return coords


# Utility function to check available methods
def get_available_methods() -> list:
    """Return list of available conformer generation methods."""
    available = [ConformerMethod.RDKIT]  # Always available if RDKit is installed

    # Check xTB
    if shutil.which("xtb"):
        available.append(ConformerMethod.XTB)

    # Check Torsional Diffusion
    try:
        import torsional_diffusion
        available.append(ConformerMethod.TORSIONAL)
    except ImportError:
        pass

    return available
