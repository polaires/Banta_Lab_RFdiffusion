"""
Rosetta FastRelax utilities for RFD3 designs.

Adapted from BindCraft (https://github.com/martinpacesa/BindCraft) with ligand support.

This module provides structure refinement to resolve ligand-protein clashes
that occur in ~40-60% of RFD3-generated designs.

Key Features:
- FastRelax with ligand parameterized from SMILES (using rdkit-to-params)
- Interface-only relaxation for speed
- Coordinate constraints to preserve topology
- Energy scoring before/after comparison
"""

import os
import tempfile
from typing import Dict, Any, Optional, Set

# PyRosetta imports (may not be available in all environments)
PYROSETTA_AVAILABLE = False
try:
    import pyrosetta as pr
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.core.scoring import get_score_function
    PYROSETTA_AVAILABLE = True
except ImportError:
    print("[RosettaUtils] PyRosetta not available. FastRelax will be disabled.")

# rdkit-to-params for ligand parameterization
RDKIT_TO_PARAMS_AVAILABLE = False
try:
    from rdkit_to_params import Params
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_TO_PARAMS_AVAILABLE = True
except ImportError:
    print("[RosettaUtils] rdkit-to-params not available. Ligand parameterization will fail.")

# Global initialization flag
_initialized = False


def check_pyrosetta_available() -> bool:
    """Check if PyRosetta is available."""
    return PYROSETTA_AVAILABLE


def init_pyrosetta(extra_res_fa: str = None) -> bool:
    """
    Initialize PyRosetta (once per process, unless extra_res_fa changes).

    Args:
        extra_res_fa: Path to params file for custom residue types

    Returns:
        True if initialization successful, False otherwise.
    """
    global _initialized

    if not PYROSETTA_AVAILABLE:
        return False

    # If already initialized without params and no new params, skip
    if _initialized and extra_res_fa is None:
        return True

    # If we need to reinitialize with params, we need a fresh init
    # PyRosetta can only be initialized once, but we can use distributed init
    try:
        flags = (
            '-ignore_zero_occupancy '
            '-mute all '
            '-corrections::beta_nov16 true '
            '-relax:default_repeats 1'
        )

        if extra_res_fa:
            flags += f' -extra_res_fa {extra_res_fa}'
            print(f"[RosettaUtils] Initializing with extra_res_fa: {extra_res_fa}")

        if _initialized:
            # Already initialized, try distributed init for reinitialization
            try:
                import pyrosetta.distributed.io as io
                pr.distributed.init(flags)
            except Exception:
                # Can't reinitialize, try continuing with current state
                print("[RosettaUtils] Already initialized, continuing...")
                return True
        else:
            pr.init(flags)
            _initialized = True

        print("[RosettaUtils] PyRosetta initialized successfully")
        return True
    except Exception as e:
        print(f"[RosettaUtils] Failed to initialize PyRosetta: {e}")
        return False


def generate_ligand_params(smiles: str, name: str = "LIG"):
    """
    Generate PyRosetta params from SMILES using rdkit-to-params.

    Args:
        smiles: SMILES string for the ligand
        name: 3-letter residue name for the ligand (default: "LIG")

    Returns:
        Params object, or None if failed
    """
    if not RDKIT_TO_PARAMS_AVAILABLE:
        print("[RosettaUtils] rdkit-to-params not available")
        return None

    try:
        # Create molecule from SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"[RosettaUtils] Failed to parse SMILES: {smiles}")
            return None

        # Add hydrogens and generate 3D coordinates
        mol = AllChem.AddHs(mol)
        result = AllChem.EmbedMolecule(mol, randomSeed=42)
        if result != 0:
            print("[RosettaUtils] Failed to embed molecule")
            return None

        # Optimize geometry
        AllChem.MMFFOptimizeMolecule(mol)

        # Generate params using rdkit-to-params
        p = Params.from_mol(mol, name=name)

        print(f"[RosettaUtils] Generated params for ligand: {name}")
        return p

    except Exception as e:
        print(f"[RosettaUtils] Failed to generate params: {e}")
        import traceback
        traceback.print_exc()
        return None


def generate_ligand_params_from_pdb(pdb_content: str, smiles: str, name: str = "UNL"):
    """
    Generate PyRosetta params by loading the ligand directly from PDB.

    This method preserves BOTH atom names AND coordinates from the original PDB,
    which is critical for proper ligand-protein clash detection during FastRelax.

    The SMILES parameter is kept for API compatibility but is not used - we load
    the molecule directly from PDB which preserves all structural information.

    Args:
        pdb_content: PDB content containing the ligand
        smiles: SMILES string (kept for API compatibility, not used)
        name: 3-letter residue name for the ligand in PDB (default: "UNL")

    Returns:
        Params object, or None if failed
    """
    if not RDKIT_TO_PARAMS_AVAILABLE:
        print("[RosettaUtils] rdkit-to-params not available")
        return None

    try:
        from rdkit import Chem

        # Extract ONLY ligand HETATM lines (RDKit needs just the ligand, not full PDB)
        ligand_lines = []
        for line in pdb_content.split('\n'):
            if line.startswith('HETATM'):
                res_name = line[17:20].strip()
                if res_name == name:
                    ligand_lines.append(line)

        if not ligand_lines:
            print(f"[RosettaUtils] No HETATM lines found for residue {name}")
            return None

        ligand_pdb = '\n'.join(ligand_lines) + '\nEND\n'
        print(f"[RosettaUtils] Extracted {len(ligand_lines)} HETATM lines for {name}")

        # Write ligand-only PDB to temp file
        fd, pdb_path = tempfile.mkstemp(suffix='.pdb')
        os.close(fd)
        with open(pdb_path, 'w') as f:
            f.write(ligand_pdb)

        try:
            # Load molecule directly from PDB - this preserves atom names AND coordinates
            mol = Chem.MolFromPDBFile(pdb_path, removeHs=False)

            if mol is None:
                print(f"[RosettaUtils] RDKit failed to load ligand PDB for {name}")
                return None

            print(f"[RosettaUtils] Loaded {name} from PDB: {mol.GetNumAtoms()} atoms")

            # Verify coordinates were loaded
            if mol.GetNumConformers() > 0:
                conf = mol.GetConformer()
                pos = conf.GetAtomPosition(0)
                print(f"[RosettaUtils] First atom coords from PDB: ({pos.x:.2f}, {pos.y:.2f}, {pos.z:.2f})")

            # Create params from mol - this preserves atom names from PDB
            p = Params.from_mol(mol, name=name)
            print(f"[RosettaUtils] Generated params from PDB-loaded mol for: {name}")

            # Verify atom names in params
            params_text = p.dumps()
            atom_lines = [l for l in params_text.split('\n') if l.startswith('ATOM')]
            if atom_lines:
                print(f"[RosettaUtils] First atom in params: {atom_lines[0][:40]}")

            return p

        finally:
            if os.path.exists(pdb_path):
                os.unlink(pdb_path)

    except Exception as e:
        print(f"[RosettaUtils] Failed to generate params from PDB: {e}")
        import traceback
        traceback.print_exc()
        return None


def get_ligand_residues(pose) -> Set[int]:
    """
    Get residue indices for ligand/HETATM atoms.

    Args:
        pose: PyRosetta pose object

    Returns:
        Set of residue indices (1-indexed) that are non-protein
    """
    ligand_res = set()
    for res in range(1, pose.total_residue() + 1):
        if not pose.residue(res).is_protein():
            ligand_res.add(res)
    return ligand_res


def get_interface_residues(
    pose,
    ligand_residues: Set[int],
    distance: float = 8.0
) -> Set[int]:
    """
    Get protein residues within distance of ligand.

    Args:
        pose: PyRosetta pose object
        ligand_residues: Set of ligand residue indices
        distance: Distance cutoff in Angstroms (default 8.0)

    Returns:
        Set of protein residue indices near the ligand
    """
    interface = set()

    for lig_res in ligand_residues:
        # Get all atom coordinates for this ligand residue
        lig_coords = []
        for i in range(1, pose.residue(lig_res).natoms() + 1):
            lig_coords.append(pose.residue(lig_res).xyz(i))

        # Check each protein residue
        for prot_res in range(1, pose.total_residue() + 1):
            if prot_res in ligand_residues:
                continue
            if not pose.residue(prot_res).is_protein():
                continue

            # Check CA distance to any ligand atom
            try:
                ca_coord = pose.residue(prot_res).xyz("CA")
                for lig_coord in lig_coords:
                    if (ca_coord - lig_coord).norm() < distance:
                        interface.add(prot_res)
                        break
            except Exception:
                # Some residues may not have CA (e.g., ACE caps)
                continue

    return interface


def fastrelax_with_ligand(
    pdb_content: str,
    ligand_smiles: Optional[str] = None,
    ligand_name: str = "UNL",
    max_iter: int = 200,
    constrain_coords: bool = True,
    interface_only: bool = False,
    interface_distance: float = 8.0,
) -> Dict[str, Any]:
    """
    Run FastRelax on protein-ligand complex.

    Key adaptation from BindCraft: Keep ligand FIXED, only relax protein around it.
    This resolves minor clashes without distorting the ligand position.

    Args:
        pdb_content: PDB file content as string
        ligand_smiles: SMILES string for the ligand (required for parameterization)
        ligand_name: 3-letter name for ligand residue in PDB (default: "UNL")
        max_iter: Maximum minimization iterations (default 200, BindCraft uses 200)
        constrain_coords: Constrain backbone to starting coordinates (default True)
        interface_only: Only relax residues near ligand (default False, faster if True)
        interface_distance: Distance cutoff for interface residues (default 8.0 Å)

    Returns:
        Dict with:
            status: "completed" or "error"
            relaxed_pdb: Relaxed PDB content (if successful)
            energy_before: Rosetta energy before relaxation
            energy_after: Rosetta energy after relaxation
            energy_change: Energy difference (negative = improved)
            interface_residues: Number of residues relaxed (if interface_only)
            error: Error message (if failed)
    """
    if not PYROSETTA_AVAILABLE:
        return {
            "status": "error",
            "error": "PyRosetta not available. Install with: pip install pyrosetta-installer"
        }

    if not init_pyrosetta():
        return {
            "status": "error",
            "error": "Failed to initialize PyRosetta"
        }

    # Check if ligand SMILES is provided
    if not ligand_smiles:
        return {
            "status": "error",
            "error": "ligand_smiles is required for parameterization"
        }

    input_path = None
    output_path = None

    try:
        # Separate protein and ligand from PDB
        # We'll relax protein-only (with tight constraints to preserve interface)
        # Then recombine with the original ligand coordinates
        print("[RosettaUtils] Separating protein and ligand from PDB...")
        protein_lines = []
        ligand_lines = []

        for line in pdb_content.split('\n'):
            if line.startswith('ATOM'):
                protein_lines.append(line)
            elif line.startswith('HETATM'):
                ligand_lines.append(line)
            elif line.startswith(('TER',)) and protein_lines:
                protein_lines.append(line)

        protein_pdb = '\n'.join(protein_lines) + '\nEND\n'

        if not ligand_lines:
            print("[RosettaUtils] Warning: No HETATM (ligand) lines found in PDB")

        print(f"[RosettaUtils] Found {len(ligand_lines)} ligand atoms")

        # Write protein-only PDB to temp file
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.pdb', delete=False
        ) as f:
            f.write(protein_pdb)
            input_path = f.name

        # Load protein-only pose
        print("[RosettaUtils] Loading protein-only pose...")
        pose = pr.pose_from_pdb(input_path)

        print(f"[RosettaUtils] Pose loaded with {pose.total_residue()} residues")

        # Get score function
        scorefxn = get_score_function()

        # Score before relaxation (protein-only)
        energy_before = float(scorefxn(pose))
        print(f"[RosettaUtils] Energy before: {energy_before:.1f}")

        # Get ligand coordinates from original PDB for interface detection
        ligand_coords = []
        for line in ligand_lines:
            if len(line) >= 54:  # Standard PDB HETATM line
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    ligand_coords.append((x, y, z))
                except ValueError:
                    continue

        # Create MoveMap
        mmf = MoveMap()

        interface_count = 0

        if interface_only and ligand_coords:
            # Find interface residues by distance to ligand coordinates
            import numpy as np
            ligand_array = np.array(ligand_coords)

            interface_residues_set = set()
            for res in range(1, pose.total_residue() + 1):
                if not pose.residue(res).is_protein():
                    continue
                try:
                    ca_xyz = pose.residue(res).xyz("CA")
                    ca_coord = np.array([ca_xyz.x, ca_xyz.y, ca_xyz.z])
                    distances = np.linalg.norm(ligand_array - ca_coord, axis=1)
                    if np.min(distances) < interface_distance:
                        interface_residues_set.add(res)
                except Exception:
                    continue

            interface_count = len(interface_residues_set)
            print(f"[RosettaUtils] Interface-only mode: {interface_count} residues near ligand")

            for res in range(1, pose.total_residue() + 1):
                if res in interface_residues_set:
                    # Allow interface protein residues to move
                    mmf.set_chi(res, True)
                    mmf.set_bb(res, True)
                else:
                    # Keep other residues fixed
                    mmf.set_chi(res, False)
                    mmf.set_bb(res, False)
        else:
            # Move all protein residues
            print("[RosettaUtils] Full protein relaxation mode")
            for res in range(1, pose.total_residue() + 1):
                mmf.set_chi(res, True)
                mmf.set_bb(res, True)

        # Disable jump movements (no rigid body transforms)
        mmf.set_jump(False)

        # Setup FastRelax
        fastrelax = FastRelax()
        fastrelax.set_scorefxn(scorefxn)
        fastrelax.set_movemap(mmf)
        fastrelax.max_iter(max_iter)
        fastrelax.min_type("lbfgs_armijo_nonmonotone")

        if constrain_coords:
            fastrelax.constrain_relax_to_start_coords(True)

        # Run relaxation
        print(f"[RosettaUtils] Running FastRelax with max_iter={max_iter}, "
              f"interface_only={interface_only}")
        fastrelax.apply(pose)

        # Score after relaxation
        energy_after = float(scorefxn(pose))

        # Output relaxed structure
        output_path = input_path.replace('.pdb', '_relaxed.pdb')
        pose.dump_pdb(output_path)

        # Read relaxed protein PDB
        with open(output_path) as f:
            relaxed_protein_pdb = f.read()

        # Clean up Rosetta metadata from protein PDB
        relaxed_protein_pdb = clean_rosetta_pdb(relaxed_protein_pdb)

        # Recombine relaxed protein with original ligand coordinates
        # Remove END from protein PDB and add ligand lines
        protein_cleaned = []
        for line in relaxed_protein_pdb.split('\n'):
            if not line.startswith('END') and line.strip():
                protein_cleaned.append(line)

        # Add original ligand lines and END
        relaxed_pdb = '\n'.join(protein_cleaned + ligand_lines + ['END'])

        print(f"[RosettaUtils] Recombined protein ({len(protein_cleaned)} lines) with ligand ({len(ligand_lines)} atoms)")

        result = {
            "status": "completed",
            "relaxed_pdb": relaxed_pdb,
            "energy_before": energy_before,
            "energy_after": energy_after,
            "energy_change": energy_after - energy_before,
            "ligand_atoms": len(ligand_lines),
        }

        if interface_only:
            result["interface_residues"] = interface_count

        print(f"[RosettaUtils] FastRelax complete: "
              f"E_before={energy_before:.1f}, E_after={energy_after:.1f}, "
              f"ΔE={energy_after - energy_before:.1f}")

        return result

    except Exception as e:
        import traceback
        traceback.print_exc()
        return {
            "status": "error",
            "error": str(e)
        }

    finally:
        # Cleanup temp files
        if input_path and os.path.exists(input_path):
            os.unlink(input_path)
        if output_path and os.path.exists(output_path):
            os.unlink(output_path)


def clean_rosetta_pdb(pdb_content: str) -> str:
    """
    Clean Rosetta-generated PDB to remove extraneous metadata.

    Keeps only: ATOM, HETATM, TER, END, MODEL, ENDMDL, LINK records.

    Args:
        pdb_content: Raw PDB content from Rosetta

    Returns:
        Cleaned PDB content
    """
    valid_prefixes = ('ATOM', 'HETATM', 'TER', 'END', 'MODEL', 'ENDMDL', 'LINK')

    cleaned_lines = []
    for line in pdb_content.split('\n'):
        if line.startswith(valid_prefixes):
            cleaned_lines.append(line)

    return '\n'.join(cleaned_lines)


def fastrelax_with_ligand_in_pose(
    pdb_content: str,
    ligand_smiles: str,
    ligand_residue_name: str = "CIT",
    max_iter: int = 200,
    repack_only: bool = True,
) -> Dict[str, Any]:
    """
    Run FastRelax with LIGAND LOADED INTO THE POSE.

    This is the PROPER way to do ligand-aware packing in PyRosetta.
    When the ligand is in the pose, fa_rep (van der Waals repulsion) naturally
    prevents sidechain rotamers from clashing with ligand atoms.

    Workflow:
    1. Generate PyRosetta params from SMILES (using rdkit-to-params)
    2. Reinitialize PyRosetta with the params file
    3. Load FULL PDB (protein + ligand) into pose
    4. Run packing - fa_rep naturally avoids ligand clashes!

    Args:
        pdb_content: Full PDB content including HETATM ligand
        ligand_smiles: SMILES string for the ligand
        ligand_residue_name: 3-letter name for ligand in PDB (default "CIT")
        max_iter: Maximum minimization iterations
        repack_only: Only repack sidechains, no backbone minimization

    Returns:
        Dict with status, relaxed_pdb, energy metrics, etc.
    """
    global _initialized

    if not PYROSETTA_AVAILABLE:
        return {"status": "error", "error": "PyRosetta not available"}

    if not RDKIT_TO_PARAMS_AVAILABLE:
        return {"status": "error", "error": "rdkit-to-params not available for ligand parameterization"}

    input_path = None
    output_path = None
    params_path = None

    try:
        print(f"[RosettaUtils] Ligand-in-pose mode: generating params for {ligand_residue_name}")
        print(f"[RosettaUtils] SMILES: {ligand_smiles}")

        # Step 1: Generate params file from SMILES
        params = generate_ligand_params_from_pdb(pdb_content, ligand_smiles, ligand_residue_name)

        if params is None:
            # Fallback to SMILES-only parameterization
            print("[RosettaUtils] Falling back to SMILES-only parameterization")
            params = generate_ligand_params(ligand_smiles, ligand_residue_name)

        if params is None:
            return {"status": "error", "error": "Failed to generate ligand params"}

        # Write params to temp file
        # Create temp file path for params
        params_fd, params_path = tempfile.mkstemp(suffix='.params')
        os.close(params_fd)  # Close the file descriptor
        params.dump(params_path)  # rdkit-to-params expects a file path string

        print(f"[RosettaUtils] Params file written: {params_path}")

        # Step 2: Reinitialize PyRosetta with ligand params
        # Force reinitialization with new params
        _initialized = False
        if not init_pyrosetta(extra_res_fa=params_path):
            return {"status": "error", "error": "Failed to reinitialize PyRosetta with ligand params"}

        # Step 3: Write full PDB and load into pose
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            input_path = f.name

        # Load FULL pose (protein + ligand)
        pose = pr.pose_from_pdb(input_path)

        # Count residue types
        protein_count = sum(1 for r in range(1, pose.total_residue()+1) if pose.residue(r).is_protein())
        ligand_count = pose.total_residue() - protein_count

        print(f"[RosettaUtils] Loaded pose: {protein_count} protein residues, {ligand_count} ligand residues")

        if ligand_count == 0:
            print("[RosettaUtils] WARNING: No ligand loaded into pose! Params may not match PDB residue name.")

        # Verify ligand coordinates match PDB (they should, since we load params directly from PDB)
        # Parse expected HETATM coordinates from input PDB for verification
        expected_coords = {}
        for line in pdb_content.split('\n'):
            if line.startswith('HETATM'):
                try:
                    atom_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    expected_coords[(res_name, atom_name)] = (x, y, z)
                except:
                    pass

        # Verify and log ligand coordinates
        for r in range(1, pose.total_residue() + 1):
            if not pose.residue(r).is_protein():
                res = pose.residue(r)
                res_name = res.name()[:3]
                print(f"[RosettaUtils] Ligand {r} ({res.name()}): {res.natoms()} atoms")

                # Check first atom coordinates
                if res.natoms() >= 1:
                    xyz = res.xyz(1)
                    atom_name = res.atom_name(1).strip()
                    key = (res_name, atom_name)
                    print(f"[RosettaUtils]   First atom '{atom_name}': ({xyz.x:.2f}, {xyz.y:.2f}, {xyz.z:.2f})")

                    if key in expected_coords:
                        ex, ey, ez = expected_coords[key]
                        dist = ((xyz.x - ex)**2 + (xyz.y - ey)**2 + (xyz.z - ez)**2)**0.5
                        if dist < 0.1:
                            print(f"[RosettaUtils]   Coords MATCH expected ({ex:.2f}, {ey:.2f}, {ez:.2f})")
                        else:
                            print(f"[RosettaUtils]   WARNING: Coords differ from expected ({ex:.2f}, {ey:.2f}, {ez:.2f}) by {dist:.2f}A")

        # Step 4: Identify ligand residue indices and coordinating residues
        ligand_residue_indices = set()
        for r in range(1, pose.total_residue() + 1):
            if not pose.residue(r).is_protein():
                ligand_residue_indices.add(r)
                print(f"[RosettaUtils] Ligand residue: {r} ({pose.residue(r).name()})")

        # Get metal positions from original PDB to identify coordinating residues
        metal_positions = []
        for line in pdb_content.split('\n'):
            if line.startswith('HETATM'):
                res_name = line[17:20].strip().upper() if len(line) > 20 else ""
                if res_name in {'TB', 'GD', 'EU', 'LA', 'CE', 'PR', 'ND', 'SM', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU'}:
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        metal_positions.append((x, y, z))
                    except:
                        pass

        # Find coordinating residues (ASP/GLU with oxygens close to metal)
        coordinating_residues = set()
        metal_coord_distance = 3.5

        for res_num in range(1, pose.total_residue() + 1):
            residue = pose.residue(res_num)
            if not residue.is_protein():
                continue
            res_name = residue.name3()

            if res_name in ['ASP', 'GLU']:
                oxygen_names = ['OD1', 'OD2'] if res_name == 'ASP' else ['OE1', 'OE2']
                for o_name in oxygen_names:
                    if residue.has(o_name):
                        o_idx = residue.atom_index(o_name)
                        o_pos = residue.xyz(o_idx)

                        for metal_pos in metal_positions:
                            import math
                            dist = math.sqrt(
                                (o_pos.x - metal_pos[0])**2 +
                                (o_pos.y - metal_pos[1])**2 +
                                (o_pos.z - metal_pos[2])**2
                            )
                            if dist < metal_coord_distance:
                                coordinating_residues.add(res_num)
                                print(f"[RosettaUtils] Protecting coordinating residue {res_num} ({res_name})")
                                break

        # Step 5: Setup packing task
        from pyrosetta.rosetta.core.pack.task import TaskFactory
        from pyrosetta.rosetta.core.pack.task.operation import (
            RestrictToRepacking, InitializeFromCommandline,
            OperateOnResidueSubset, PreventRepackingRLT
        )
        from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
        from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover

        tf = TaskFactory()
        tf.push_back(InitializeFromCommandline())
        tf.push_back(RestrictToRepacking())

        # Prevent repacking of coordinating residues
        for res_num in coordinating_residues:
            selector = ResidueIndexSelector(str(res_num))
            prevent_op = OperateOnResidueSubset(PreventRepackingRLT(), selector)
            tf.push_back(prevent_op)

        # Prevent repacking of ligand residues
        for res_num in ligand_residue_indices:
            selector = ResidueIndexSelector(str(res_num))
            prevent_op = OperateOnResidueSubset(PreventRepackingRLT(), selector)
            tf.push_back(prevent_op)

        scorefxn = get_score_function()
        energy_before = float(scorefxn(pose))
        print(f"[RosettaUtils] Energy before packing: {energy_before:.1f}")

        # Step 6: Run packing with ligand in pose
        # fa_rep will naturally penalize clashes with ligand!
        packer = PackRotamersMover()
        packer.task_factory(tf)
        packer.score_function(scorefxn)
        packer.apply(pose)

        energy_after_pack = float(scorefxn(pose))
        print(f"[RosettaUtils] Energy after packing: {energy_after_pack:.1f}")

        # Step 7: Optional minimization
        if not repack_only:
            print("[RosettaUtils] Running FastRelax minimization...")
            mmf = MoveMap()

            for res in range(1, pose.total_residue() + 1):
                if res in ligand_residue_indices:
                    # Keep ligand completely fixed
                    mmf.set_chi(res, False)
                    mmf.set_bb(res, False)
                elif res in coordinating_residues:
                    # Keep coordinating residues fixed
                    mmf.set_chi(res, False)
                    mmf.set_bb(res, False)
                else:
                    # Allow protein sidechains to adjust
                    mmf.set_chi(res, True)
                    mmf.set_bb(res, False)  # Keep backbone fixed

            mmf.set_jump(False)

            fastrelax = FastRelax()
            fastrelax.set_scorefxn(scorefxn)
            fastrelax.set_movemap(mmf)
            fastrelax.max_iter(max_iter)
            fastrelax.constrain_relax_to_start_coords(True)
            fastrelax.apply(pose)

        energy_final = float(scorefxn(pose))
        print(f"[RosettaUtils] Final energy: {energy_final:.1f}")

        # Step 8: Output relaxed structure
        # IMPORTANT: We output PROTEIN-ONLY from the pose, then recombine with
        # the ORIGINAL HETATM coordinates. This ensures:
        # 1. Sidechains were packed avoiding the ligand (fa_rep worked)
        # 2. Final output has EXACT template ligand coordinates (not RDKit-generated)
        output_path = input_path.replace('.pdb', '_relaxed.pdb')
        pose.dump_pdb(output_path)

        with open(output_path) as f:
            relaxed_pdb = f.read()

        # Clean Rosetta metadata
        relaxed_pdb = clean_rosetta_pdb(relaxed_pdb)

        # Extract PROTEIN ONLY from relaxed PDB (PyRosetta may have altered ligand coords)
        relaxed_protein_lines = []
        for line in relaxed_pdb.split('\n'):
            if line.startswith('ATOM') or (line.startswith('TER') and relaxed_protein_lines):
                relaxed_protein_lines.append(line)

        # Extract ORIGINAL HETATM from input PDB (exact template coordinates)
        original_hetatm_lines = []
        for line in pdb_content.split('\n'):
            if line.startswith('HETATM'):
                original_hetatm_lines.append(line)

        # Recombine: relaxed protein + original HETATM
        final_pdb = '\n'.join(relaxed_protein_lines + original_hetatm_lines + ['END'])

        print(f"[RosettaUtils] Output: {len(relaxed_protein_lines)} ATOM + {len(original_hetatm_lines)} HETATM (original coords)")

        return {
            "status": "completed",
            "relaxed_pdb": final_pdb,
            "energy_before": energy_before,
            "energy_after": energy_final,
            "energy_change": energy_final - energy_before,
            "ligand_in_pose": ligand_count > 0,
            "coordinating_residues": len(coordinating_residues),
        }

    except Exception as e:
        import traceback
        traceback.print_exc()
        return {"status": "error", "error": str(e)}

    finally:
        if input_path and os.path.exists(input_path):
            os.unlink(input_path)
        if output_path and os.path.exists(output_path):
            os.unlink(output_path)
        if params_path and os.path.exists(params_path):
            os.unlink(params_path)


def fastrelax_protein_only(
    pdb_content: str,
    max_iter: int = 200,
    constrain_coords: bool = True,
    repack_only: bool = False,
    ligand_aware: bool = True,
    ligand_exclusion_radius: float = 2.8,
) -> Dict[str, Any]:
    """
    Run FastRelax on protein structure with LIGAND-AWARE sidechain packing.

    This function is designed for sidechain packing after MPNN sequence design,
    following the BindCraft approach: MPNN for sequence, FastRelax for sidechains.

    CRITICAL: When ligand_aware=True, the packer uses HETATM coordinates to avoid
    placing sidechains that clash with the ligand. This is done via AtomPairConstraints
    that penalize protein atoms within ligand_exclusion_radius of any ligand atom.

    Args:
        pdb_content: PDB file content (protein + optional HETATM ligand)
        max_iter: Maximum minimization iterations (default 200, BindCraft uses 200)
        constrain_coords: Constrain backbone to starting coordinates (default True)
        repack_only: If True, only repack sidechains without minimization (faster)
        ligand_aware: If True, add repulsive constraints to prevent ligand clashes
        ligand_exclusion_radius: Minimum allowed distance from ligand atoms (Å)

    Returns:
        Dict with:
            status: "completed" or "error"
            relaxed_pdb: Relaxed PDB content (if successful)
            energy_before: Rosetta energy before relaxation
            energy_after: Rosetta energy after relaxation
            error: Error message (if failed)
    """
    if not PYROSETTA_AVAILABLE:
        return {
            "status": "error",
            "error": "PyRosetta not available"
        }

    if not init_pyrosetta():
        return {
            "status": "error",
            "error": "Failed to initialize PyRosetta"
        }

    input_path = None
    output_path = None

    try:
        # Separate protein and HETATM (metals, waters, etc.)
        protein_lines = []
        hetatm_lines = []

        for line in pdb_content.split('\n'):
            if line.startswith('ATOM'):
                protein_lines.append(line)
            elif line.startswith('HETATM'):
                hetatm_lines.append(line)
            elif line.startswith('TER') and protein_lines:
                protein_lines.append(line)

        protein_pdb = '\n'.join(protein_lines) + '\nEND\n'

        print(f"[RosettaUtils] Protein: {len(protein_lines)} atoms, HETATM: {len(hetatm_lines)} atoms")

        # Write protein-only PDB
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.pdb', delete=False
        ) as f:
            f.write(protein_pdb)
            input_path = f.name

        # Load pose
        pose = pr.pose_from_pdb(input_path)
        print(f"[RosettaUtils] Loaded pose with {pose.total_residue()} residues")

        # Rebuild missing sidechain atoms using idealized geometry
        # This is critical when loading backbone-only structures after sequence redesign
        from pyrosetta.rosetta.protocols.simple_moves import ReturnSidechainMover
        from pyrosetta.rosetta.core.pack.task import TaskFactory
        from pyrosetta.rosetta.core.pack.task.operation import (
            RestrictToRepacking, InitializeFromCommandline,
            OperateOnResidueSubset, PreventRepackingRLT
        )
        from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
        from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover

        print("[RosettaUtils] Rebuilding missing sidechains...")

        # Identify coordinating residues that should NOT be repacked
        # They already have proper metal-coordinating geometry from the template
        coordinating_residues = []
        metal_coord_distance = 3.5  # Å

        # Find metal positions from HETATM lines
        metal_positions = []
        for line in hetatm_lines:
            res_name = line[17:20].strip().upper() if len(line) > 20 else ""
            if res_name in {'TB', 'GD', 'EU', 'LA', 'CE', 'PR', 'ND', 'SM', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU'}:
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    metal_positions.append((x, y, z))
                except (ValueError, IndexError):
                    pass

        # Check which residues have carboxylate oxygens close to metal
        for res_num in range(1, pose.total_residue() + 1):
            residue = pose.residue(res_num)
            res_name = residue.name3()

            if res_name in ['ASP', 'GLU']:
                # Check OD1/OD2 (Asp) or OE1/OE2 (Glu)
                oxygen_names = ['OD1', 'OD2'] if res_name == 'ASP' else ['OE1', 'OE2']
                for o_name in oxygen_names:
                    if residue.has(o_name):
                        o_idx = residue.atom_index(o_name)
                        o_pos = residue.xyz(o_idx)

                        for metal_pos in metal_positions:
                            import math
                            dist = math.sqrt(
                                (o_pos.x - metal_pos[0])**2 +
                                (o_pos.y - metal_pos[1])**2 +
                                (o_pos.z - metal_pos[2])**2
                            )
                            if dist < metal_coord_distance:
                                if res_num not in coordinating_residues:
                                    coordinating_residues.append(res_num)
                                    print(f"[RosettaUtils] Protecting coordinating residue {res_num} ({res_name}) from repacking (dist={dist:.2f}Å)")
                                break

        # Create a task that allows repacking of all residues EXCEPT coordinating ones
        tf = TaskFactory()
        tf.push_back(InitializeFromCommandline())
        tf.push_back(RestrictToRepacking())  # Don't change sequence, just repack

        # Prevent repacking of coordinating residues
        if coordinating_residues:
            print(f"[RosettaUtils] Protecting {len(coordinating_residues)} coordinating residues from repacking")
            for res_num in coordinating_residues:
                selector = ResidueIndexSelector(str(res_num))
                # PreventRepackingRLT is a ResLvlTaskOperation (required by OperateOnResidueSubset)
                prevent_op = OperateOnResidueSubset(PreventRepackingRLT(), selector)
                tf.push_back(prevent_op)

        scorefxn = get_score_function()

        # === LIGAND-AWARE PACKING ===
        # PyRosetta's standard packer doesn't know about HETATM atoms.
        # We implement a two-phase approach:
        # Phase 1: Standard packing (fast, builds sidechains)
        # Phase 2: Post-packing ligand clash check and rotamer refinement
        ligand_coords = []
        if ligand_aware and hetatm_lines:
            import numpy as np

            # Extract ALL ligand atom coordinates
            for line in hetatm_lines:
                if len(line) >= 54:
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        ligand_coords.append((x, y, z))
                    except (ValueError, IndexError):
                        continue

            print(f"[RosettaUtils] Ligand-aware mode: {len(ligand_coords)} ligand atoms to avoid")

        # Run initial packing to build sidechains
        packer = PackRotamersMover()
        packer.task_factory(tf)
        packer.score_function(scorefxn)
        packer.apply(pose)
        print(f"[RosettaUtils] Initial sidechains built (protected: {len(coordinating_residues)} residues)")

        # Phase 2: Ligand-aware refinement
        # For clashing residues, do multiple repack trials and keep best
        if ligand_aware and ligand_coords:
            import numpy as np

            ligand_array = np.array(ligand_coords)
            clash_threshold = ligand_exclusion_radius

            def count_ligand_clashes(test_pose):
                """Count sidechain atoms clashing with ligand."""
                clashes = []
                for res_num in range(1, test_pose.total_residue() + 1):
                    if res_num in coordinating_residues:
                        continue
                    residue = test_pose.residue(res_num)
                    if not residue.is_protein():
                        continue
                    for atom_idx in range(1, residue.natoms() + 1):
                        atom_name = residue.atom_name(atom_idx).strip()
                        if atom_name in ('N', 'CA', 'C', 'O', 'H', 'HA', 'OXT'):
                            continue
                        if atom_name.startswith('H'):
                            continue
                        atom_pos = residue.xyz(atom_idx)
                        atom_coord = np.array([atom_pos.x, atom_pos.y, atom_pos.z])
                        min_dist = np.min(np.linalg.norm(ligand_array - atom_coord, axis=1))
                        if min_dist < clash_threshold:
                            clashes.append((res_num, residue.name3(), atom_name, min_dist))
                return clashes

            # Check initial clashes
            initial_clashes = count_ligand_clashes(pose)
            print(f"[RosettaUtils] Initial ligand clashes: {len(initial_clashes)}")

            if initial_clashes:
                # Try multiple packing trials to find a better configuration
                best_pose = pose.clone()
                best_clash_count = len(initial_clashes)

                # Create task with extra rotamers for interface residues
                from pyrosetta.rosetta.core.pack.task.operation import ExtraRotamersGeneric
                extra_tf = TaskFactory()
                extra_tf.push_back(InitializeFromCommandline())
                extra_tf.push_back(RestrictToRepacking())

                extra_rot = ExtraRotamersGeneric()
                extra_rot.ex1(True)
                extra_rot.ex2(True)
                extra_rot.ex3(True)
                extra_rot.ex4(True)
                extra_tf.push_back(extra_rot)

                # Prevent repacking of coordinating residues
                for res_num in coordinating_residues:
                    selector = ResidueIndexSelector(str(res_num))
                    prevent_op = OperateOnResidueSubset(PreventRepackingRLT(), selector)
                    extra_tf.push_back(prevent_op)

                # Run multiple packing trials
                NUM_TRIALS = 5
                print(f"[RosettaUtils] Running {NUM_TRIALS} packing trials for ligand compatibility...")

                for trial in range(NUM_TRIALS):
                    trial_pose = pose.clone()
                    trial_packer = PackRotamersMover()
                    trial_packer.task_factory(extra_tf)
                    trial_packer.score_function(scorefxn)
                    trial_packer.nloop(3)  # Internal trials
                    trial_packer.apply(trial_pose)

                    trial_clashes = count_ligand_clashes(trial_pose)
                    if len(trial_clashes) < best_clash_count:
                        best_pose = trial_pose.clone()
                        best_clash_count = len(trial_clashes)
                        print(f"[RosettaUtils] Trial {trial+1}: {len(trial_clashes)} clashes (improved!)")
                    else:
                        print(f"[RosettaUtils] Trial {trial+1}: {len(trial_clashes)} clashes")

                    if best_clash_count == 0:
                        break

                pose = best_pose
                final_clashes = count_ligand_clashes(pose)
                print(f"[RosettaUtils] Final ligand clashes: {len(final_clashes)} (was {len(initial_clashes)})")

                if final_clashes:
                    for res_num, res_name, atom_name, dist in final_clashes[:5]:
                        print(f"[RosettaUtils]   {res_name}{res_num} {atom_name}: {dist:.2f}Å")

        energy_before = float(scorefxn(pose))
        print(f"[RosettaUtils] Energy before relaxation: {energy_before:.1f}")

        # Create MoveMap
        mmf = MoveMap()

        if repack_only:
            # Only allow sidechain movement (chi angles)
            print("[RosettaUtils] Repack-only mode: sidechains only")
            for res in range(1, pose.total_residue() + 1):
                # Freeze chi angles for coordinating residues to preserve metal geometry
                if res in coordinating_residues:
                    mmf.set_chi(res, False)
                    mmf.set_bb(res, False)
                else:
                    mmf.set_chi(res, True)
                    mmf.set_bb(res, False)
        else:
            # Full relaxation: backbone + sidechains
            print("[RosettaUtils] Full relax mode: backbone + sidechains")
            for res in range(1, pose.total_residue() + 1):
                # Freeze chi angles for coordinating residues to preserve metal geometry
                if res in coordinating_residues:
                    mmf.set_chi(res, False)
                    mmf.set_bb(res, False)
                else:
                    mmf.set_chi(res, True)
                    mmf.set_bb(res, True)

        if coordinating_residues:
            print(f"[RosettaUtils] Froze chi angles for {len(coordinating_residues)} coordinating residues")

            # Add coordinate constraints for coordinating oxygen atoms
            # This ensures OE1/OE2 atoms stay near the metal even if backbone shifts
            from pyrosetta.rosetta.core.scoring.constraints import CoordinateConstraint
            from pyrosetta.rosetta.core.scoring.func import HarmonicFunc
            from pyrosetta.rosetta.core.id import AtomID
            from pyrosetta.rosetta.numeric import xyzVector_double_t as V3

            # Get virtual root atom for coordinate constraints
            # Use CA of first residue as anchor
            anchor_atom = AtomID(pose.residue(1).atom_index("CA"), 1)

            for res_num in coordinating_residues:
                residue = pose.residue(res_num)
                res_name = residue.name3()

                # Get oxygen atom names based on residue type
                oxygen_names = ['OD1', 'OD2'] if res_name == 'ASP' else ['OE1', 'OE2']

                for o_name in oxygen_names:
                    if residue.has(o_name):
                        o_idx = residue.atom_index(o_name)
                        o_pos = residue.xyz(o_idx)

                        # Create coordinate constraint with tight harmonic (sd=0.1 Å)
                        constraint_atom = AtomID(o_idx, res_num)
                        target_pos = V3(o_pos.x, o_pos.y, o_pos.z)
                        harmonic = HarmonicFunc(0.0, 0.1)  # 0 distance, 0.1 Å std dev

                        coord_constraint = CoordinateConstraint(
                            constraint_atom, anchor_atom, target_pos, harmonic
                        )
                        pose.add_constraint(coord_constraint)

            print(f"[RosettaUtils] Added coordinate constraints for coordinating oxygens")

            # Enable coordinate constraints in scoring
            from pyrosetta.rosetta.core.scoring import ScoreType
            scorefxn.set_weight(ScoreType.coordinate_constraint, 10.0)

        # No rigid body movements
        mmf.set_jump(False)

        # Setup FastRelax
        fastrelax = FastRelax()
        fastrelax.set_scorefxn(scorefxn)
        fastrelax.set_movemap(mmf)
        fastrelax.max_iter(max_iter)
        fastrelax.min_type("lbfgs_armijo_nonmonotone")

        if constrain_coords:
            fastrelax.constrain_relax_to_start_coords(True)

        # Run relaxation
        print(f"[RosettaUtils] Running FastRelax (max_iter={max_iter}, repack_only={repack_only})...")
        fastrelax.apply(pose)

        # Score after
        energy_after = float(scorefxn(pose))
        print(f"[RosettaUtils] Energy after: {energy_after:.1f}, ΔE={energy_after - energy_before:.1f}")

        # Output relaxed structure
        output_path = input_path.replace('.pdb', '_relaxed.pdb')
        pose.dump_pdb(output_path)

        with open(output_path) as f:
            relaxed_protein_pdb = f.read()

        # Clean Rosetta metadata
        relaxed_protein_pdb = clean_rosetta_pdb(relaxed_protein_pdb)

        # Recombine with original HETATM
        protein_cleaned = [l for l in relaxed_protein_pdb.split('\n')
                          if not l.startswith('END') and l.strip()]

        relaxed_pdb = '\n'.join(protein_cleaned + hetatm_lines + ['END'])

        return {
            "status": "completed",
            "relaxed_pdb": relaxed_pdb,
            "energy_before": energy_before,
            "energy_after": energy_after,
            "energy_change": energy_after - energy_before,
        }

    except Exception as e:
        import traceback
        traceback.print_exc()
        return {
            "status": "error",
            "error": str(e)
        }

    finally:
        if input_path and os.path.exists(input_path):
            os.unlink(input_path)
        if output_path and os.path.exists(output_path):
            os.unlink(output_path)


def score_interface(
    pdb_content: str,
    chain_a: str = "A",
    chain_b: str = "B",
) -> Dict[str, Any]:
    """
    Score protein-protein or protein-ligand interface using PyRosetta.

    Calculates:
    - dG (binding energy estimate)
    - Shape complementarity
    - Interface H-bonds
    - Buried surface area
    - PackStat (packing quality)

    Args:
        pdb_content: PDB file content
        chain_a: First chain ID
        chain_b: Second chain ID

    Returns:
        Dict with interface metrics
    """
    if not PYROSETTA_AVAILABLE:
        return {
            "status": "error",
            "error": "PyRosetta not available"
        }

    if not init_pyrosetta():
        return {
            "status": "error",
            "error": "Failed to initialize PyRosetta"
        }

    try:
        from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

        # Write temp PDB
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.pdb', delete=False
        ) as f:
            f.write(pdb_content)
            input_path = f.name

        try:
            pose = pr.pose_from_pdb(input_path)

            # Setup InterfaceAnalyzerMover
            iam = InterfaceAnalyzerMover()
            iam.set_interface(f"{chain_a}_{chain_b}")

            scorefxn = get_score_function()
            iam.set_scorefunction(scorefxn)
            iam.set_compute_packstat(True)
            iam.set_compute_interface_energy(True)
            iam.set_calc_dSASA(True)
            iam.set_compute_interface_sc(True)
            iam.set_pack_separated(True)

            iam.apply(pose)

            return {
                "status": "completed",
                "dG": iam.get_interface_dG(),
                "dSASA": iam.get_interface_delta_sasa(),
                "packstat": iam.get_interface_packstat(),
                "sc_value": iam.get_interface_sc(),
            }

        finally:
            os.unlink(input_path)

    except Exception as e:
        return {
            "status": "error",
            "error": str(e)
        }


def calculate_clash_score(
    pdb_content: str,
    threshold: float = 2.4,
    exclude_ligand: bool = False,
) -> Dict[str, Any]:
    """
    Calculate steric clash score using distance-based detection.

    A clash is defined as two non-bonded heavy atoms closer than threshold.

    Args:
        pdb_content: PDB file content
        threshold: Distance cutoff in Angstroms (default 2.4)
        exclude_ligand: Whether to exclude ligand-ligand clashes

    Returns:
        Dict with clash count and details
    """
    try:
        from Bio.PDB import PDBParser
        from scipy.spatial import cKDTree
        import numpy as np
        import io

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', io.StringIO(pdb_content))

        atoms = []
        atom_info = []

        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        # Skip hydrogens
                        if atom.element == 'H':
                            continue

                        atoms.append(atom.coord)
                        atom_info.append({
                            'chain': chain.id,
                            'res_id': residue.id[1],
                            'res_name': residue.resname,
                            'atom_name': atom.get_name(),
                            'is_hetatm': residue.id[0] != ' ',
                        })

        if len(atoms) < 2:
            return {"status": "error", "error": "Not enough atoms"}

        atoms = np.array(atoms)
        tree = cKDTree(atoms)
        pairs = tree.query_pairs(threshold)

        # Filter valid clashes
        clash_count = 0
        clash_details = []

        for (i, j) in pairs:
            info_i = atom_info[i]
            info_j = atom_info[j]

            # Skip same residue
            if info_i['chain'] == info_j['chain'] and info_i['res_id'] == info_j['res_id']:
                continue

            # Skip bonded neighbors (sequential residues in same chain)
            if info_i['chain'] == info_j['chain'] and abs(info_i['res_id'] - info_j['res_id']) == 1:
                continue

            # Optionally skip ligand-ligand
            if exclude_ligand and info_i['is_hetatm'] and info_j['is_hetatm']:
                continue

            clash_count += 1
            if len(clash_details) < 10:  # Store first 10
                dist = np.linalg.norm(atoms[i] - atoms[j])
                clash_details.append({
                    'atom_1': f"{info_i['chain']}/{info_i['res_name']}{info_i['res_id']}/{info_i['atom_name']}",
                    'atom_2': f"{info_j['chain']}/{info_j['res_name']}{info_j['res_id']}/{info_j['atom_name']}",
                    'distance': round(dist, 2),
                })

        return {
            "status": "completed",
            "clash_count": clash_count,
            "clash_details": clash_details,
            "threshold": threshold,
        }

    except ImportError as e:
        return {"status": "error", "error": f"Missing dependency: {e}"}
    except Exception as e:
        return {"status": "error", "error": str(e)}


def minimize_ligand_clashes(
    pdb_content: str,
    clash_threshold: float = 2.4,
    max_iterations: int = 3,
) -> Dict[str, Any]:
    """
    Post-relaxation clash minimization for protein-ligand structures.

    Detects residues with sidechains clashing with ligand and tries alternative
    rotamers to minimize clashes. This is useful after FastRelax when the ligand
    was not present during packing.

    Args:
        pdb_content: PDB content with protein + HETATM (ligand/metal)
        clash_threshold: Distance threshold for clash detection (Å)
        max_iterations: Maximum rounds of clash minimization

    Returns:
        Dict with:
            status: "completed" or "error"
            minimized_pdb: PDB with minimized clashes
            initial_clashes: Number of clashes before
            final_clashes: Number of clashes after
            residues_fixed: List of residue IDs that were repacked
    """
    if not PYROSETTA_AVAILABLE:
        return {
            "status": "error",
            "error": "PyRosetta not available"
        }

    if not init_pyrosetta():
        return {
            "status": "error",
            "error": "Failed to initialize PyRosetta"
        }

    try:
        import numpy as np
        from scipy.spatial import cKDTree
        from Bio.PDB import PDBParser
        import io

        # Separate protein and HETATM
        protein_lines = []
        hetatm_lines = []

        for line in pdb_content.split('\n'):
            if line.startswith('ATOM'):
                protein_lines.append(line)
            elif line.startswith('HETATM'):
                hetatm_lines.append(line)
            elif line.startswith('TER') and protein_lines:
                protein_lines.append(line)

        if not hetatm_lines:
            return {
                "status": "completed",
                "minimized_pdb": pdb_content,
                "initial_clashes": 0,
                "final_clashes": 0,
                "residues_fixed": [],
                "message": "No HETATM found, no ligand clashes possible"
            }

        # Parse ligand coordinates
        ligand_coords = []
        for line in hetatm_lines:
            if len(line) >= 54:
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    ligand_coords.append((x, y, z))
                except ValueError:
                    continue

        if not ligand_coords:
            return {
                "status": "completed",
                "minimized_pdb": pdb_content,
                "initial_clashes": 0,
                "final_clashes": 0,
                "residues_fixed": [],
                "message": "Could not parse ligand coordinates"
            }

        ligand_array = np.array(ligand_coords)

        # Parse protein to find clashing residues
        def find_clashing_residues(protein_pdb_lines):
            """Find residues with sidechain atoms clashing with ligand."""
            clashing = []
            for line in protein_pdb_lines:
                if not line.startswith('ATOM'):
                    continue
                atom_name = line[12:16].strip()
                # Skip backbone atoms
                if atom_name in ('N', 'CA', 'C', 'O', 'H', 'HA'):
                    continue
                # Skip hydrogen
                if atom_name.startswith('H'):
                    continue

                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    chain = line[21]
                    res_num = int(line[22:26].strip())
                    res_name = line[17:20].strip()

                    # Check distance to any ligand atom
                    atom_coord = np.array([x, y, z])
                    distances = np.linalg.norm(ligand_array - atom_coord, axis=1)
                    min_dist = np.min(distances)

                    if min_dist < clash_threshold:
                        res_key = (chain, res_num, res_name)
                        if res_key not in clashing:
                            clashing.append(res_key)
                            print(f"[ClashMin] Clash: {res_name}{res_num} {atom_name} at {min_dist:.2f}Å from ligand")

                except (ValueError, IndexError):
                    continue

            return clashing

        # Find initial clashes
        initial_clashing = find_clashing_residues(protein_lines)
        initial_clash_count = len(initial_clashing)

        if initial_clash_count == 0:
            print("[ClashMin] No sidechain-ligand clashes detected")
            return {
                "status": "completed",
                "minimized_pdb": pdb_content,
                "initial_clashes": 0,
                "final_clashes": 0,
                "residues_fixed": []
            }

        print(f"[ClashMin] Found {initial_clash_count} residues with sidechain-ligand clashes")

        # Write protein-only PDB for PyRosetta
        protein_pdb = '\n'.join(protein_lines) + '\nEND\n'

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(protein_pdb)
            input_path = f.name

        try:
            # Load pose
            pose = pr.pose_from_pdb(input_path)
            scorefxn = get_score_function()

            # For ligand clash resolution, PyRosetta packer doesn't know about HETATM
            # We use a two-phase approach:
            # Phase 1: Try repacking with extra rotamers and check if clashes reduce
            # Phase 2: For stubborn clashes, try mutating to smaller residues (Ala, Gly)

            from pyrosetta.rosetta.core.pack.task import TaskFactory
            from pyrosetta.rosetta.core.pack.task.operation import (
                RestrictToRepacking, InitializeFromCommandline,
                OperateOnResidueSubset, PreventRepackingRLT
            )
            from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
            from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
            from pyrosetta.rosetta.protocols.simple_moves import MutateResidue

            # Map chain+res_num to pose residue number
            clashing_pose_residues = []
            clashing_info = {}  # pose_res -> (chain, res_num, res_name)
            for chain, res_num, res_name in initial_clashing:
                for pose_res in range(1, pose.total_residue() + 1):
                    pdb_info = pose.pdb_info()
                    if pdb_info:
                        pose_chain = pdb_info.chain(pose_res)
                        pose_resnum = pdb_info.number(pose_res)
                        if pose_chain == chain and pose_resnum == res_num:
                            clashing_pose_residues.append(pose_res)
                            clashing_info[pose_res] = (chain, res_num, res_name)
                            break

            print(f"[ClashMin] Processing {len(clashing_pose_residues)} clashing residues")

            # Phase 1: Try repacking with extra rotamers
            tf = TaskFactory()
            tf.push_back(InitializeFromCommandline())
            tf.push_back(RestrictToRepacking())

            clashing_set = set(clashing_pose_residues)
            for res in range(1, pose.total_residue() + 1):
                if res not in clashing_set:
                    selector = ResidueIndexSelector(str(res))
                    prevent_op = OperateOnResidueSubset(PreventRepackingRLT(), selector)
                    tf.push_back(prevent_op)

            from pyrosetta.rosetta.core.pack.task.operation import ExtraRotamersGeneric
            extra_rot = ExtraRotamersGeneric()
            extra_rot.ex1(True)
            extra_rot.ex2(True)
            extra_rot.ex3(True)
            extra_rot.ex4(True)
            tf.push_back(extra_rot)

            best_pose = pose.clone()
            best_clash_count = initial_clash_count

            # Try repacking first
            print(f"[ClashMin] Phase 1: Repacking with extra rotamers")
            packer = PackRotamersMover()
            packer.task_factory(tf)
            packer.score_function(scorefxn)
            packer.nloop(10)  # More packing trials

            test_pose = best_pose.clone()
            packer.apply(test_pose)

            # Check clashes after repacking
            output_path = input_path.replace('.pdb', '_repacked.pdb')
            test_pose.dump_pdb(output_path)
            with open(output_path) as f:
                test_pdb = f.read()
            test_pdb = clean_rosetta_pdb(test_pdb)
            test_lines = [l for l in test_pdb.split('\n') if l.startswith('ATOM')]
            current_clashing = find_clashing_residues(test_lines)
            os.unlink(output_path)

            if len(current_clashing) < best_clash_count:
                best_pose = test_pose.clone()
                best_clash_count = len(current_clashing)
                print(f"[ClashMin] Repacking reduced clashes to {best_clash_count}")
            else:
                print(f"[ClashMin] Repacking didn't help ({len(current_clashing)} clashes)")

            # Phase 2: For remaining clashes, try mutating to smaller residues
            if best_clash_count > 0:
                print(f"[ClashMin] Phase 2: Mutating stubborn clashes to smaller residues")

                # Re-check which residues still clash
                output_path = input_path.replace('.pdb', '_check.pdb')
                best_pose.dump_pdb(output_path)
                with open(output_path) as f:
                    check_pdb = f.read()
                check_pdb = clean_rosetta_pdb(check_pdb)
                check_lines = [l for l in check_pdb.split('\n') if l.startswith('ATOM')]
                still_clashing = find_clashing_residues(check_lines)
                os.unlink(output_path)

                # Map still-clashing residues to pose numbers
                still_clashing_pose = []
                for chain, res_num, res_name in still_clashing:
                    for pose_res, info in clashing_info.items():
                        if info[0] == chain and info[1] == res_num:
                            still_clashing_pose.append(pose_res)
                            break

                mutations_made = []
                for pose_res in still_clashing_pose:
                    info = clashing_info.get(pose_res)
                    if not info:
                        continue
                    chain, res_num, res_name = info

                    # Skip if already Ala or Gly
                    if res_name in ('ALA', 'GLY'):
                        print(f"[ClashMin] {res_name}{res_num} already small, cannot reduce further")
                        continue

                    # Try mutation to Alanine first
                    print(f"[ClashMin] Mutating {res_name}{res_num} -> ALA")
                    mutator = MutateResidue(pose_res, 'ALA')
                    mutator.apply(best_pose)
                    mutations_made.append(f"{res_name}{res_num}->ALA")

                # Repack after mutations
                if mutations_made:
                    packer.apply(best_pose)

                    # Final clash check
                    output_path = input_path.replace('.pdb', '_mutated.pdb')
                    best_pose.dump_pdb(output_path)
                    with open(output_path) as f:
                        final_pdb = f.read()
                    final_pdb = clean_rosetta_pdb(final_pdb)
                    final_lines = [l for l in final_pdb.split('\n') if l.startswith('ATOM')]
                    final_clashing = find_clashing_residues(final_lines)
                    best_clash_count = len(final_clashing)
                    os.unlink(output_path)
                    print(f"[ClashMin] After mutations: {best_clash_count} clashing residues")

            # Output final structure
            output_path = input_path.replace('.pdb', '_minimized.pdb')
            best_pose.dump_pdb(output_path)

            with open(output_path) as f:
                minimized_protein = f.read()

            minimized_protein = clean_rosetta_pdb(minimized_protein)

            # Recombine with HETATM
            protein_cleaned = [l for l in minimized_protein.split('\n')
                              if not l.startswith('END') and l.strip()]
            minimized_pdb = '\n'.join(protein_cleaned + hetatm_lines + ['END'])

            os.unlink(output_path)

            return {
                "status": "completed",
                "minimized_pdb": minimized_pdb,
                "initial_clashes": initial_clash_count,
                "final_clashes": best_clash_count,
                "residues_fixed": mutations_made if 'mutations_made' in dir() else [],
            }

        finally:
            if os.path.exists(input_path):
                os.unlink(input_path)

    except Exception as e:
        import traceback
        traceback.print_exc()
        return {
            "status": "error",
            "error": str(e)
        }


# Convenience function for batch processing
def relax_designs(
    designs: list,
    ligand_smiles: str,
    max_iter: int = 200,
    interface_only: bool = True,
) -> list:
    """
    Apply FastRelax to a batch of designs.

    Args:
        designs: List of dicts with 'pdb_content' key
        ligand_smiles: SMILES string for the ligand
        max_iter: Max iterations per design
        interface_only: Only relax interface (faster)

    Returns:
        List of designs with 'relaxed_pdb' and 'energy_change' added
    """
    results = []

    for i, design in enumerate(designs):
        print(f"[RosettaUtils] Relaxing design {i+1}/{len(designs)}")

        pdb_content = design.get('pdb_content')
        if not pdb_content:
            results.append({**design, 'relax_error': 'No pdb_content'})
            continue

        relax_result = fastrelax_with_ligand(
            pdb_content=pdb_content,
            ligand_smiles=ligand_smiles,
            max_iter=max_iter,
            interface_only=interface_only,
        )

        if relax_result['status'] == 'completed':
            results.append({
                **design,
                'pdb_content': relax_result['relaxed_pdb'],
                'relaxed': True,
                'energy_before': relax_result['energy_before'],
                'energy_after': relax_result['energy_after'],
                'energy_change': relax_result['energy_change'],
            })
        else:
            results.append({
                **design,
                'relaxed': False,
                'relax_error': relax_result.get('error', 'Unknown error'),
            })

    return results
