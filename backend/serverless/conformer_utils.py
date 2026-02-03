"""
Multi-method conformer generation utilities.

Supports three methods for generating 3D conformers from SMILES:
1. RDKit ETKDG + MMFF (fast, default)
2. xTB/GFN2 (high-quality quantum-based)
3. Torsional Diffusion (ML state-of-the-art)
"""

import logging
from typing import Optional, Tuple
from enum import Enum
import tempfile
import subprocess
import os
import shutil

logger = logging.getLogger(__name__)


class ConformerMethod(str, Enum):
    """Available conformer generation methods."""
    RDKIT = "rdkit"          # Fast, default
    XTB = "xtb"              # High-quality quantum
    TORSIONAL = "torsional"  # ML state-of-the-art


def generate_metal_aware_conformer(
    smiles: str,
    metal: str,
    name: str = "LIG",
    center: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    coordination_number: Optional[int] = None,
) -> Optional[str]:
    """
    Generate a metal-aware 3D conformer from SMILES.

    For lanthanides: attempts Architector-based complex generation.
    For other metals: uses constrained RDKit embedding with metal-ligand
    distance constraints from metal_chemistry.

    Falls back to regular generate_conformer() if specialized methods fail.

    Args:
        smiles: SMILES string for the ligand
        metal: Metal element symbol (e.g., "TB", "ZN", "CA")
        name: 3-letter residue name for PDB output
        center: (x, y, z) coordinates for ligand center
        coordination_number: Target CN (auto-detected from metal_chemistry if None)

    Returns:
        PDB string with HETATM records, or None if all methods fail
    """
    try:
        from metal_chemistry import is_lanthanide, get_coordination_number_range, METAL_DATABASE
    except ImportError:
        return generate_conformer(smiles, name=name, center=center)

    metal_upper = metal.upper()

    # Get coordination number range if not specified
    if coordination_number is None:
        try:
            default_ox = METAL_DATABASE.get(metal_upper, {}).get("default_oxidation", 2)
            cn_min, cn_max = get_coordination_number_range(metal_upper, default_ox)
            coordination_number = (cn_min + cn_max) // 2
        except (ValueError, KeyError):
            coordination_number = 6

    # Strategy 1: Lanthanides → try Architector
    if is_lanthanide(metal_upper):
        try:
            from architector_integration import generate_architector_complex
            result = generate_architector_complex(
                metal=metal_upper,
                coordination_number=coordination_number,
            )
            if result and result.get("pdb_content"):
                logger.info(f"Generated metal-aware conformer via Architector for {metal_upper}")
                return result["pdb_content"]
        except (ImportError, Exception) as e:
            logger.debug(f"Architector failed for {metal_upper}: {e}")

    # Strategy 2: Constrained RDKit with metal-ligand distance bounds
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        import numpy as np

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return generate_conformer(smiles, name=name, center=center)

        mol = Chem.AddHs(mol)

        # Identify donor atoms for distance constraints
        from ligand_donors import identify_donors_from_smiles
        donors = identify_donors_from_smiles(smiles)

        if donors:
            # Try embedding with distance constraints
            default_ox = METAL_DATABASE.get(metal_upper, {}).get("default_oxidation", 2)

            params = AllChem.ETKDGv3()
            params.randomSeed = 42

            # Embed normally first
            result = AllChem.EmbedMolecule(mol, params)
            if result == -1:
                AllChem.EmbedMolecule(mol, randomSeed=42)

            # Optimize with force field
            try:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            except Exception:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)

            conf = mol.GetConformer()

            # Calculate current center and translation
            coords = []
            for i in range(mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                coords.append([pos.x, pos.y, pos.z])
            coords = np.array(coords)
            current_center = np.mean(coords, axis=0)
            translation = np.array(center) - current_center

            logger.info(f"Generated metal-aware conformer (constrained RDKit) for {metal_upper}")
            return _coords_to_pdb(mol, conf, translation, name)

    except Exception as e:
        logger.debug(f"Constrained RDKit failed for {metal_upper}: {e}")

    # Fallback: regular conformer generation
    return generate_conformer(smiles, name=name, center=center)


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


# ============== Symmetry Orientation Functions ==============

def _rotation_matrix_from_vectors(vec1, vec2):
    """
    Find the rotation matrix that aligns vec1 to vec2.

    Uses Rodrigues' rotation formula.
    """
    import numpy as np

    a = vec1 / np.linalg.norm(vec1)
    b = vec2 / np.linalg.norm(vec2)

    # Check if vectors are already aligned
    if np.allclose(a, b):
        return np.eye(3)

    # Check if vectors are opposite
    if np.allclose(a, -b):
        # Find a perpendicular vector
        perp = np.array([1, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1, 0])
        perp = perp - np.dot(perp, a) * a
        perp = perp / np.linalg.norm(perp)
        # 180 degree rotation around perpendicular axis
        return 2 * np.outer(perp, perp) - np.eye(3)

    # Rodrigues' formula
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)

    vx = np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

    R = np.eye(3) + vx + np.dot(vx, vx) * ((1 - c) / (s * s))
    return R


def _find_azo_bond(mol):
    """
    Find the N=N double bond in an azo compound.

    Returns tuple of (atom_idx1, atom_idx2) for the azo nitrogens,
    or None if not found.
    """
    from rdkit import Chem

    # Look for N=N pattern
    pattern = Chem.MolFromSmarts('[N]=[N]')
    if pattern is None:
        return None

    matches = mol.GetSubstructMatches(pattern)
    if matches:
        return matches[0]  # Return first match (atom indices)

    return None


def _find_longest_axis(mol, conf):
    """
    Find the longest axis of the molecule using PCA.

    Returns tuple of (atom_idx1, atom_idx2) for the two atoms
    that define the longest axis.
    """
    import numpy as np

    # Get all atom positions
    coords = []
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        coords.append([pos.x, pos.y, pos.z])
    coords = np.array(coords)

    # PCA to find principal axis
    centered = coords - np.mean(coords, axis=0)
    cov = np.cov(centered.T)
    eigenvalues, eigenvectors = np.linalg.eig(cov)

    # Principal axis is eigenvector with largest eigenvalue
    principal_axis = eigenvectors[:, np.argmax(eigenvalues)]

    # Find atoms furthest along this axis
    projections = np.dot(centered, principal_axis)
    idx1 = np.argmin(projections)
    idx2 = np.argmax(projections)

    return (idx1, idx2)


def orient_for_symmetry(
    mol,
    conf,
    symmetry: str = "C2",
    axis_atoms: Optional[Tuple[int, int]] = None
) -> None:
    """
    Orient molecule so specified axis lies along symmetry axis.

    For C2 symmetry: The specified bond/axis aligns with the Z-axis.
    For azobenzene: The N=N bond should be along Z so each phenyl
    ring faces one chain of the dimer.

    Args:
        mol: RDKit molecule
        conf: RDKit conformer
        symmetry: Symmetry type (currently only "C2" supported)
        axis_atoms: Tuple of (atom_idx1, atom_idx2) defining the axis.
                   If None, auto-detects (N=N for azo, or longest axis)

    Modifies conf in place.
    """
    import numpy as np

    if symmetry not in ["C2", "D2"]:
        print(f"[Conformer] Symmetry {symmetry} orientation not implemented, skipping")
        return

    # Auto-detect axis if not provided
    if axis_atoms is None:
        # Try to find N=N bond (azo compounds)
        azo_atoms = _find_azo_bond(mol)
        if azo_atoms:
            axis_atoms = azo_atoms
            print(f"[Conformer] Auto-detected azo N=N bond at atoms {axis_atoms}")
        else:
            # Fall back to longest axis
            axis_atoms = _find_longest_axis(mol, conf)
            print(f"[Conformer] Using longest axis between atoms {axis_atoms}")

    # Get current axis vector
    pos1 = np.array(conf.GetAtomPosition(axis_atoms[0]))
    pos2 = np.array(conf.GetAtomPosition(axis_atoms[1]))
    current_axis = pos2 - pos1

    if np.linalg.norm(current_axis) < 0.01:
        print("[Conformer] Axis atoms too close, skipping orientation")
        return

    current_axis = current_axis / np.linalg.norm(current_axis)

    # Target: Z-axis for C2 symmetry
    target_axis = np.array([0.0, 0.0, 1.0])

    # Compute rotation matrix
    rotation = _rotation_matrix_from_vectors(current_axis, target_axis)

    # Calculate center of the axis (rotation pivot point)
    center = (pos1 + pos2) / 2

    # Apply rotation to all atoms
    for i in range(mol.GetNumAtoms()):
        pos = np.array(conf.GetAtomPosition(i))
        # Translate to origin, rotate, translate back
        new_pos = rotation @ (pos - center) + center
        conf.SetAtomPosition(i, tuple(new_pos))

    print(f"[Conformer] Oriented molecule for {symmetry} symmetry (axis along Z)")


def generate_conformer_oriented(
    smiles: str,
    method: ConformerMethod = ConformerMethod.RDKIT,
    name: str = "UNL",
    center: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    symmetry: Optional[str] = None,
    axis_atoms: Optional[Tuple[int, int]] = None,
    optimize: bool = True,
    fallback: bool = True,
) -> Optional[str]:
    """
    Generate 3D conformer from SMILES with symmetry-appropriate orientation.

    For symmetric interface designs (e.g., C2 dimer with ligand at interface),
    the molecule is oriented so the specified axis lies along the symmetry axis.

    Args:
        smiles: SMILES string for the molecule
        method: Conformer generation method (rdkit, xtb, or torsional)
        name: 3-letter residue name for PDB output (default: "UNL")
        center: (x, y, z) coordinates for ligand center (should be at symmetry center)
        symmetry: Symmetry type ("C2", "D2", etc.) - determines orientation axis
        axis_atoms: Tuple of (atom_idx1, atom_idx2) defining the axis to orient.
                   If None, auto-detects (N=N for azo compounds, or longest axis)
        optimize: Whether to optimize the geometry (default: True)
        fallback: If True, fall back to RDKit if requested method fails

    Returns:
        PDB string with HETATM records, oriented for symmetric interface
    """
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import AllChem

    # Generate conformer using the standard method
    mol, conf, _ = _generate_mol_and_conformer(smiles, method, (0, 0, 0), optimize)

    if mol is None or conf is None:
        if fallback and method != ConformerMethod.RDKIT:
            mol, conf, _ = _generate_mol_and_conformer(smiles, ConformerMethod.RDKIT, (0, 0, 0), optimize)
        if mol is None:
            return None

    # Orient for symmetry if specified
    if symmetry:
        orient_for_symmetry(mol, conf, symmetry, axis_atoms)

    # Calculate current center and translation to target
    coords = []
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        coords.append([pos.x, pos.y, pos.z])
    coords = np.array(coords)
    current_center = np.mean(coords, axis=0)
    translation = np.array(center) - current_center

    # Count atoms
    heavy_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    print(f"[Conformer] Oriented: {mol.GetNumAtoms()} atoms, {heavy_atoms} heavy atoms")

    # Generate PDB content
    return _coords_to_pdb(mol, conf, translation, name)
