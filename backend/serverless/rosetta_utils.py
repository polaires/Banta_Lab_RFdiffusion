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
    Generate PyRosetta params from SMILES, matching atom names from PDB.

    This is the preferred method when you have an existing PDB with ligand coordinates,
    as it ensures atom names match between the params file and PDB.

    Args:
        pdb_content: PDB content containing the ligand
        smiles: SMILES string for the ligand
        name: 3-letter residue name for the ligand in PDB (default: "UNL")

    Returns:
        Params object, or None if failed
    """
    if not RDKIT_TO_PARAMS_AVAILABLE:
        print("[RosettaUtils] rdkit-to-params not available")
        return None

    try:
        # Write PDB to temp file
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.pdb', delete=False
        ) as f:
            f.write(pdb_content)
            pdb_path = f.name

        try:
            # Use from_smiles_w_pdbfile to match atom names
            p = Params.from_smiles_w_pdbfile(
                pdb_file=pdb_path,
                smiles=smiles,
                name=name
            )
            print(f"[RosettaUtils] Generated params from PDB for ligand: {name}")
            return p
        finally:
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
