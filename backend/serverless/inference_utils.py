"""
Inference Utilities for RunPod Serverless Handler

Shared functions for running RFD3, RF3, MPNN inference and utilities.
Extracted from main.py for use in serverless environment.
"""

import os
import io
import json
import random
import subprocess
import tempfile
from typing import Dict, Any, Optional, List, Tuple


# ============== Utility Functions ==============

def get_gpu_info() -> Dict[str, Any]:
    """Get GPU information via nvidia-smi"""
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=name,memory.total", "--format=csv,noheader,nounits"],
            capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0:
            parts = result.stdout.strip().split(", ")
            return {
                "available": True,
                "name": parts[0] if parts else "Unknown",
                "memory_gb": float(parts[1]) / 1024 if len(parts) > 1 else 0,
            }
    except Exception:
        pass
    return {"available": False, "name": None, "memory_gb": None}


def check_foundry_available() -> bool:
    """Check if Foundry Python modules are importable"""
    try:
        from rfd3.engine import RFD3InferenceEngine
        return True
    except ImportError:
        pass

    try:
        from rf3.inference_engines.rf3 import RF3InferenceEngine
        return True
    except ImportError:
        pass

    return False


# Standard amino acid 3-letter codes
STANDARD_AMINO_ACIDS = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    # Modified/alternate names
    'MSE', 'SEC', 'PYL',  # Selenomethionine, Selenocysteine, Pyrrolysine
}


def extract_protein_contig_from_pdb(pdb_content: str, ligand_code: str = None) -> str:
    """
    Extract a proper contig string from PDB content, excluding metals/ligands.

    This is critical for partial diffusion (partial_t) because RFD3 separates
    non-polymer residues into a different chain, causing residue number gaps.

    Args:
        pdb_content: PDB file content
        ligand_code: Optional ligand code to exclude (e.g., "TB", "ZN")

    Returns:
        Contig string like "A1-50,A52-100" with proper gaps for non-protein residues
    """
    chain_residues: Dict[str, List[int]] = {}

    for line in pdb_content.split('\n'):
        if not line.startswith('ATOM'):
            continue

        try:
            res_name = line[17:20].strip()
            chain_id = line[21] if len(line) > 21 else 'A'
            res_num = int(line[22:26].strip())

            # Skip non-standard amino acids (metals, ligands, etc.)
            if res_name not in STANDARD_AMINO_ACIDS:
                continue

            # Skip specified ligand code
            if ligand_code and res_name == ligand_code:
                continue

            if chain_id not in chain_residues:
                chain_residues[chain_id] = []

            if res_num not in chain_residues[chain_id]:
                chain_residues[chain_id].append(res_num)
        except (ValueError, IndexError):
            continue

    # Build contig segments for each chain
    all_contigs = []

    for chain_id in sorted(chain_residues.keys()):
        residues = sorted(chain_residues[chain_id])
        if not residues:
            continue

        # Find contiguous segments
        segments = []
        start = residues[0]
        end = residues[0]

        for i in range(1, len(residues)):
            if residues[i] == end + 1:
                # Continue current segment
                end = residues[i]
            else:
                # Gap found - save current segment and start new one
                segments.append((start, end))
                start = residues[i]
                end = residues[i]

        # Don't forget the last segment
        segments.append((start, end))

        # Build contig string for this chain
        for start, end in segments:
            if start == end:
                all_contigs.append(f"{chain_id}{start}")
            else:
                all_contigs.append(f"{chain_id}{start}-{end}")

    return ",".join(all_contigs)


# Common metal ion codes (single-atom HETATM entries)
METAL_CODES = {
    'FE', 'ZN', 'CA', 'MG', 'MN', 'CO', 'CU', 'NI',  # Common transition metals
    'TB', 'EU', 'GD', 'SM', 'LA', 'CE', 'PR', 'ND',  # Lanthanides
    'NA', 'K', 'CD', 'HG', 'PB', 'AG', 'AU',          # Other metals
}


def detect_metal_in_pdb(pdb_content: str) -> Optional[str]:
    """
    Detect the first metal ion code found in PDB content.

    Args:
        pdb_content: PDB file content

    Returns:
        Metal code (e.g., "FE", "ZN") or None if no metal found
    """
    for line in pdb_content.split('\n'):
        if not line.startswith('HETATM'):
            continue
        try:
            res_name = line[17:20].strip()
            if res_name in METAL_CODES:
                return res_name
        except (ValueError, IndexError):
            continue
    return None


def replace_metal_in_pdb(pdb_content: str, target_metal: str, source_metal: str = None) -> str:
    """
    Replace metal ion in PDB content with a different metal.

    This is needed for metal binding redesign where the original structure
    has one metal (e.g., FE) but we want to redesign for a different metal (e.g., TB).

    Args:
        pdb_content: PDB file content
        target_metal: Target metal code (e.g., "TB")
        source_metal: Source metal to replace (auto-detected if None)

    Returns:
        Modified PDB content with metal replaced
    """
    if not source_metal:
        source_metal = detect_metal_in_pdb(pdb_content)
        if not source_metal:
            print(f"[RFD3] No metal found in PDB to replace")
            return pdb_content

    if source_metal == target_metal:
        print(f"[RFD3] Metal already matches target: {target_metal}")
        return pdb_content

    print(f"[RFD3] Replacing metal {source_metal} -> {target_metal}")

    # Prepare replacement strings with proper padding
    # PDB format: columns 17-20 for residue name (3 chars, right-padded)
    source_res = source_metal.ljust(3)  # e.g., "FE "
    target_res = target_metal.ljust(3)  # e.g., "TB "

    # Element symbol is in columns 76-78 (2 chars, right-justified)
    source_elem = source_metal.rjust(2)  # e.g., "FE"
    target_elem = target_metal[:2].rjust(2)  # e.g., "TB" (max 2 chars for element)

    modified_lines = []
    for line in pdb_content.split('\n'):
        if line.startswith('HETATM'):
            try:
                res_name = line[17:20]
                if res_name.strip() == source_metal:
                    # Replace residue name at columns 17-20
                    line = line[:17] + target_res + line[20:]

                    # Replace atom name at columns 12-16 (usually same as element for metals)
                    atom_name = line[12:16].strip()
                    if atom_name == source_metal:
                        target_atom = target_metal.ljust(4)[:4]  # 4 chars, left-padded
                        line = line[:12] + target_atom + line[16:]

                    # Replace element symbol at columns 76-78 if present
                    if len(line) >= 78:
                        elem = line[76:78]
                        if elem.strip() == source_metal or elem.strip() == source_metal[:2]:
                            line = line[:76] + target_elem + line[78:]
            except (ValueError, IndexError):
                pass
        modified_lines.append(line)

    return '\n'.join(modified_lines)


def atom_array_to_pdb(atom_array) -> str:
    """Convert biotite AtomArray to PDB string"""
    from biotite.structure.io.pdb import PDBFile

    pdb_file = PDBFile()
    pdb_file.set_structure(atom_array)
    buf = io.StringIO()
    pdb_file.write(buf)
    return buf.getvalue()


def atom_array_to_cif(atom_array) -> str:
    """Convert biotite AtomArray to CIF string"""
    try:
        from atomworks.io.utils.io_utils import to_cif_file

        with tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as f:
            temp_path = f.name

        to_cif_file(atom_array, temp_path)
        with open(temp_path) as f:
            content = f.read()
        os.unlink(temp_path)
        return content
    except ImportError:
        # Fallback to biotite CIF writer
        from biotite.structure.io.pdbx import PDBxFile, set_structure

        pdbx_file = PDBxFile()
        set_structure(pdbx_file, atom_array)
        buf = io.StringIO()
        pdbx_file.write(buf)
        return buf.getvalue()


def extract_plddt_from_cif(cif_content: str) -> Optional[List[float]]:
    """
    Extract pLDDT values from CIF file.
    RF3 stores pLDDT in B_iso_or_equiv column for CA atoms.
    """
    import re

    plddt_values = []
    ca_bfactors = {}

    # Parse atom_site loop for CA atoms
    in_atom_site = False
    columns = []

    for line in cif_content.split('\n'):
        line = line.strip()

        if line.startswith('loop_'):
            in_atom_site = False
            columns = []
        elif line.startswith('_atom_site.'):
            in_atom_site = True
            columns.append(line.replace('_atom_site.', ''))
        elif in_atom_site and line and not line.startswith('_') and not line.startswith('#'):
            # Data row
            if not columns:
                continue

            parts = line.split()
            if len(parts) >= len(columns):
                try:
                    row = dict(zip(columns, parts))
                    atom_id = row.get('label_atom_id', row.get('auth_atom_id', ''))
                    if atom_id == 'CA':
                        res_seq = int(row.get('label_seq_id', row.get('auth_seq_id', 0)))
                        bfactor = float(row.get('B_iso_or_equiv', 0))
                        if res_seq > 0:
                            ca_bfactors[res_seq] = bfactor
                except (ValueError, KeyError):
                    continue

    if ca_bfactors:
        # Sort by residue number and return
        for i in sorted(ca_bfactors.keys()):
            plddt_values.append(ca_bfactors[i])

    return plddt_values if plddt_values else None


def pdb_to_atom_array(content: str):
    """Convert PDB or CIF string to biotite AtomArray"""
    content_stripped = content.strip()

    if content_stripped.startswith('data_'):
        # CIF format
        from biotite.structure.io.pdbx import CIFFile, get_structure

        cif_file = CIFFile.read(io.StringIO(content))
        block = list(cif_file.values())[0] if cif_file else None
        if block is None:
            raise ValueError("Empty CIF file")
        return get_structure(block, model=1)
    else:
        # PDB format
        from biotite.structure.io.pdb import PDBFile

        pdb_file = PDBFile.read(io.StringIO(content))
        return pdb_file.get_structure(model=1)


def _generate_seed_for_symmetry(symmetry_id: str, offset: float = 10.0) -> str:
    """
    Generate seed residue(s) positioned for symmetric interface design.

    For C2 symmetry: Single residue at (+offset, 0, 0)
    RFD3 will replicate it to (-offset, 0, 0) for chain B.

    The seed provides protein atoms for RFD3's symmetry handler,
    allowing native symmetry to work with ligand-only designs.
    """
    if symmetry_id == "C2":
        # Single seed at +X offset
        return _generate_seed_residue(offset=(offset, 0.0, 0.0))
    elif symmetry_id == "C3":
        # Seed at 120° spacing
        return _generate_seed_residue(offset=(offset, 0.0, 0.0))
    else:
        # Default: single seed at offset
        return _generate_seed_residue(offset=(offset, 0.0, 0.0))


def _generate_seed_residue(offset: tuple = (8.0, 0.0, 0.0)) -> str:
    """
    Generate a minimal glycine residue at the specified offset from origin.

    This seed gives RFD3's symmetry handler something to work with for de novo
    symmetric designs with ligand.

    Args:
        offset: (x, y, z) offset from origin in Angstroms

    Returns:
        PDB string for a single glycine residue
    """
    x, y, z = offset
    # Glycine backbone coordinates (relative to CA at origin)
    # Standard geometry: N-CA=1.47Å, CA-C=1.52Å, C-O=1.23Å
    pdb_lines = [
        f"ATOM      1  N   GLY A   1    {x-1.47:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           N",
        f"ATOM      2  CA  GLY A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C",
        f"ATOM      3  C   GLY A   1    {x+1.52:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C",
        f"ATOM      4  O   GLY A   1    {x+1.52:8.3f}{y+1.23:8.3f}{z:8.3f}  1.00  0.00           O",
    ]
    return "\n".join(pdb_lines)


def _apply_post_symmetry(atom_array, symmetry_id: str):
    """
    Apply symmetry transformation to an atom array in post-processing.

    For C2 symmetry with ligand at interface:
    1. Iteratively translate protein away from ligand center
    2. Apply 180° rotation around Z-axis for chain B
    3. Check for clashes - if found, increase translation
    4. Keep ligand at center between both chains

    This creates a proper dimer with ligand at the interface, not overlapping chains.
    """
    import numpy as np
    from biotite.structure import AtomArray

    if symmetry_id != "C2":
        print(f"[RFD3] Warning: Post-processing symmetry only supports C2, got {symmetry_id}")
        return atom_array

    # Separate protein and ligand
    protein_mask = atom_array.chain_id != "L"
    ligand_mask = atom_array.chain_id == "L"

    protein_atoms = atom_array[protein_mask]
    ligand_atoms = atom_array[ligand_mask] if ligand_mask.any() else None

    if len(protein_atoms) == 0:
        print("[RFD3] Warning: No protein atoms found for symmetry")
        return atom_array

    # Get ligand center (should be near origin since we oriented it there)
    if ligand_atoms is not None and len(ligand_atoms) > 0:
        ligand_center = ligand_atoms.coord.mean(axis=0)
    else:
        ligand_center = np.array([0.0, 0.0, 0.0])

    # Store original coordinates for iterative testing
    original_coords = protein_atoms.coord.copy()

    # Calculate the protein's extent in X direction
    x_min = original_coords[:, 0].min()
    x_max = original_coords[:, 0].max()
    protein_x_extent = x_max - x_min
    print(f"[RFD3] Protein X extent: {protein_x_extent:.1f}Å (from {x_min:.1f} to {x_max:.1f})")

    # For interface design, we need to balance:
    # 1. Minimize clashes between chain A and chain B
    # 2. Maximize contacts between protein chains and ligand
    #
    # Strategy: Test translations and pick the one with best combined score
    # Key insight: To avoid clashes, protein atoms after translation must not
    # extend past origin (where chain B's atoms will be after rotation).

    clash_threshold = 2.0  # Angstroms - actual steric clash (was 2.5)
    contact_threshold = 4.5  # Angstroms - binding contacts
    ligand_coords = ligand_atoms.coord if ligand_atoms is not None else np.array([[0.0, 0.0, 0.0]])

    # VERIFY ONE-SIDED BINDING before applying symmetry
    # For azobenzene: atoms 0-5 = first phenyl (should be buried), atoms 6-13 = rest (should be exposed)
    # Check that protein primarily contacts the first half of the ligand
    if ligand_atoms is not None and len(ligand_atoms) > 10:
        # Split ligand into two halves
        half_idx = len(ligand_atoms) // 2
        first_half_coords = ligand_atoms.coord[:half_idx]
        second_half_coords = ligand_atoms.coord[half_idx:]

        # Count contacts with each half
        from scipy.spatial import distance_matrix
        dist_first = distance_matrix(original_coords, first_half_coords)
        dist_second = distance_matrix(original_coords, second_half_coords)

        contacts_first = int(np.sum(dist_first < contact_threshold))
        contacts_second = int(np.sum(dist_second < contact_threshold))

        print(f"[RFD3] One-sided binding check: contacts with first half={contacts_first}, second half={contacts_second}")

        # Ratio check: we want > 2:1 asymmetry for good half-binding
        if contacts_first > 0 and contacts_second > 0:
            ratio = contacts_first / contacts_second if contacts_second > 0 else float('inf')
            if ratio < 2.0 and ratio > 0.5:
                print(f"[RFD3] WARNING: Monomer contacts both halves of ligand (ratio={ratio:.2f})")
                print(f"[RFD3] This design may have clashes after C2 rotation")
            elif ratio >= 2.0:
                print(f"[RFD3] Good: Monomer primarily contacts first half (ratio={ratio:.2f})")
            else:
                print(f"[RFD3] Good: Monomer primarily contacts second half (ratio={1/ratio:.2f})")

    # For interface design, we want chains CLOSE to ligand but not clashing with each other
    # Test a range of translations starting from smaller values
    # The key is finding the sweet spot where both chains can bind without interchain clashes

    # Start from a smaller distance - chains CAN overlap slightly in space as long as
    # they don't have severe steric clashes
    min_translation = 6  # Start from 6Å
    max_translation = 20  # Up to 20Å

    # Test more distances with finer granularity
    translation_distances = list(range(min_translation, max_translation + 1, 1))  # Test every 1Å

    print(f"[RFD3] Testing translation distances from {min_translation}Å to {max_translation}Å")

    best_distance = 10.0  # Default
    best_score = float('-inf')
    results = []

    for push_distance in translation_distances:
        # Test this translation distance
        test_coords = original_coords.copy()
        test_coords[:, 0] += push_distance

        # Apply C2 rotation to get chain B coordinates
        rotated_coords = test_coords.copy()
        rotated_coords[:, 0] = 2 * ligand_center[0] - rotated_coords[:, 0]
        rotated_coords[:, 1] = 2 * ligand_center[1] - rotated_coords[:, 1]

        # Count clashes between chain A and chain B
        try:
            from scipy.spatial import distance_matrix
            clash_matrix = distance_matrix(test_coords, rotated_coords)
            clash_count = int(np.sum(clash_matrix < clash_threshold))

            # Count contacts: how many protein atoms are within contact_threshold of ligand
            contacts_a_matrix = distance_matrix(test_coords, ligand_coords)
            contacts_b_matrix = distance_matrix(rotated_coords, ligand_coords)
            contacts_a = int(np.sum(contacts_a_matrix < contact_threshold))
            contacts_b = int(np.sum(contacts_b_matrix < contact_threshold))
        except ImportError:
            # Fallback to manual calculation
            clash_count = 0
            contacts_a = 0
            contacts_b = 0
            for i in range(len(test_coords)):
                # Clashes
                diffs = rotated_coords - test_coords[i]
                distances = np.sqrt(np.sum(diffs * diffs, axis=1))
                clash_count += np.sum(distances < clash_threshold)
                # Contacts with ligand
                for lc in ligand_coords:
                    dist_to_lig = np.sqrt(np.sum((test_coords[i] - lc)**2))
                    if dist_to_lig < contact_threshold:
                        contacts_a += 1
                    dist_to_lig_b = np.sqrt(np.sum((rotated_coords[i] - lc)**2))
                    if dist_to_lig_b < contact_threshold:
                        contacts_b += 1

        # Calculate combined score:
        # Goal: Find BALANCE between close enough for binding and no clashes
        #
        # Key insights:
        # - Too close (10Å): 40+ clashes, 200+ contacts → severe steric problems
        # - Too far (18Å): 0 clashes, <10 contacts → weak binding
        # - Sweet spot: ~12-14Å with moderate contacts (20-50) and few clashes (<10)
        #
        # Scoring strategy:
        # - Reward contacts but CAP the benefit (diminishing returns above ~50)
        # - SEVERE penalty for interchain clashes - these cause +60 kcal/mol affinity
        # - Bonus for balanced contacts between chains

        total_contacts = contacts_a + contacts_b
        min_contacts = min(contacts_a, contacts_b)

        # Reward contacts with diminishing returns
        # First 50 contacts are valuable, above that the benefit plateaus
        if total_contacts <= 50:
            contact_score = total_contacts * 5
        else:
            # Diminishing returns above 50 contacts
            contact_score = 50 * 5 + (total_contacts - 50) * 1

        # Symmetry bonus: reward balanced contacts between chains
        symmetry_bonus = min_contacts * 3

        # Clash penalty: SEVERE - each clash contributes to positive affinity
        # We want ZERO clashes ideally
        # Any clashes > 5 means chains are overlapping significantly
        if clash_count > 20:
            # Severe overlap - reject this translation
            clash_penalty = 10000  # Effectively disqualify
        elif clash_count > 5:
            # Moderate clashing - strong penalty
            clash_penalty = clash_count * 50
        else:
            # Minimal clashing - mild penalty
            clash_penalty = clash_count * 10

        score = contact_score + symmetry_bonus - clash_penalty

        results.append({
            'distance': push_distance,
            'clashes': clash_count,
            'contacts_a': contacts_a,
            'contacts_b': contacts_b,
            'score': score
        })

        if score > best_score:
            best_score = score
            best_distance = push_distance

    # Print summary of best options
    sorted_results = sorted(results, key=lambda x: x['score'], reverse=True)[:3]
    print(f"[RFD3] Top translations by score:")
    for r in sorted_results:
        print(f"  {r['distance']}Å: score={r['score']:.0f} (clashes={r['clashes']}, contacts A={r['contacts_a']}, B={r['contacts_b']})")

    # Apply the best translation distance
    protein_atoms.coord = original_coords.copy()
    protein_atoms.coord[:, 0] += best_distance
    print(f"[RFD3] Applied translation +{best_distance:.0f}Å along X (score={best_score:.0f})")

    # Rename protein to chain A
    protein_atoms.chain_id[:] = "A"

    # Create C2 rotated copy (180° around Z-axis through ligand center)
    rotated_coords = protein_atoms.coord.copy()
    rotated_coords[:, 0] = 2 * ligand_center[0] - rotated_coords[:, 0]
    rotated_coords[:, 1] = 2 * ligand_center[1] - rotated_coords[:, 1]

    # Create chain B atoms
    chain_b = protein_atoms.copy()
    chain_b.coord = rotated_coords
    chain_b.chain_id[:] = "B"

    # Renumber residues for chain B to avoid overlap
    max_res_a = protein_atoms.res_id.max() if len(protein_atoms) > 0 else 0
    chain_b.res_id = chain_b.res_id + max_res_a

    # Combine: chain A + chain B + ligand (if present)
    if ligand_atoms is not None and len(ligand_atoms) > 0:
        result = protein_atoms + chain_b + ligand_atoms
    else:
        result = protein_atoms + chain_b

    print(f"[RFD3] Applied C2 symmetry: chain A ({len(protein_atoms)} atoms) + chain B")
    return result


# ============== Mock Implementations ==============

def generate_mock_pdb(length: int = 100) -> str:
    """Generate mock PDB content"""
    lines = ["HEADER    MOCK PROTEIN STRUCTURE"]
    for i in range(1, length + 1):
        x, y, z = i * 0.5, i * 0.3, i * 0.2
        lines.append(
            f"ATOM  {i:5d}  CA  ALA A{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
        )
    lines.append("END")
    return "\n".join(lines)


def generate_mock_fasta(num_seq: int = 8, length: int = 100) -> str:
    """Generate mock FASTA sequences"""
    aa = "ACDEFGHIKLMNPQRSTVWY"
    lines = []
    for i in range(num_seq):
        seq = "".join(random.choices(aa, k=length))
        lines.append(f">design_{i+1}")
        lines.append(seq)
    return "\n".join(lines)


def generate_mock_confidences(length: int = 100) -> Dict[str, Any]:
    """Generate mock confidence metrics"""
    return {
        "summary_confidences": {
            "overall_plddt": round(random.uniform(0.7, 0.95), 3),
            "overall_pae": round(random.uniform(2.0, 8.0), 2),
            "overall_pde": round(random.uniform(0.1, 0.3), 3),
            "ptm": round(random.uniform(0.6, 0.9), 3),
            "iptm": None,
            "ranking_score": round(random.uniform(0.5, 0.85), 3),
            "has_clash": False,
        },
        "per_residue_plddt": [round(random.uniform(0.6, 1.0), 2) for _ in range(length)],
    }


# ============== RFD3 Inference ==============

def run_rfd3_inference(
    contig: Optional[str] = None,
    length: Optional[str] = None,
    num_designs: int = 1,
    seed: Optional[int] = None,
    pdb_content: Optional[str] = None,
    # Quality settings
    num_timesteps: Optional[int] = None,
    step_scale: Optional[float] = None,
    gamma_0: Optional[float] = None,
    is_non_loopy: Optional[bool] = None,
    # Partial diffusion (refinement)
    partial_t: Optional[int] = None,
    # Symmetry
    symmetry: Optional[Dict[str, Any]] = None,
    # Small molecule / enzyme design
    ligand: Optional[str] = None,
    ligand_smiles: Optional[str] = None,  # SMILES string for organic molecules
    ligand_sdf: Optional[str] = None,     # SDF content with 3D coordinates
    ligand_center: Optional[Tuple[float, float, float]] = None,  # Center for SMILES-generated ligand
    conformer_method: Optional[str] = None,  # "rdkit", "xtb", or "torsional"
    interface_ligand: bool = False,  # Place ligand at symmetric interface (requires symmetry)
    select_fixed_atoms: Optional[Dict[str, str]] = None,
    unindex: Optional[str] = None,
    # RASA conditioning (binding pocket design)
    select_buried: Optional[Dict[str, str]] = None,
    select_exposed: Optional[Dict[str, str]] = None,
    select_partially_buried: Optional[Dict[str, str]] = None,
    # Protein binder design
    hotspots: Optional[List[str]] = None,
    infer_ori_strategy: Optional[str] = None,
    # Nucleic acid binder design
    na_chains: Optional[str] = None,
    ori_token: Optional[List[float]] = None,
    select_hbond_donor: Optional[Dict[str, str]] = None,
    select_hbond_acceptor: Optional[Dict[str, str]] = None,
    # Covalent modifications (enzyme design)
    covalent_bonds: Optional[List[Dict[str, Any]]] = None,
    # Mock mode
    use_mock: bool = False
) -> Dict[str, Any]:
    """
    Run RFD3 inference with full parameter support.

    Supports multiple design tasks:
    - De novo protein design (length only)
    - Protein binder design (contig + hotspots)
    - Small molecule binder design (ligand + RASA conditioning)
    - Nucleic acid binder design (na_chains + ori_token + H-bond conditioning)
    - Enzyme scaffold design (ligand + unindex + fixed atoms)
    - Symmetric oligomer design (symmetry config)
    - Structure refinement (partial_t)

    Ligand Input Options:
    - ligand: Code for metal ions/small molecules already in PDB (e.g., "ZN", "CA", "ATP")
    - ligand_smiles: SMILES string for organic molecules (RFD3 generates 3D coords)
    - ligand_sdf: SDF content with explicit 3D coordinates for the ligand

    For best binding pocket generation, use RASA conditioning:
    - select_buried: Force residues near ligand to be buried (creates enclosed pocket)
    - select_exposed: Force surface residues to be solvent-accessible
    - unindex: Mark coordinating residues as flexible during diffusion
    """

    if use_mock:
        mock_length = length or contig or "100"
        return run_rfd3_mock(mock_length, num_designs)

    try:
        from rfd3.engine import RFD3InferenceConfig, RFD3InferenceEngine

        # Set seed for reproducibility
        if seed is not None:
            try:
                from lightning.fabric import seed_everything
                seed_everything(seed)
                print(f"[RFD3] Set seed to {seed}")
            except ImportError:
                pass

        # Build specification dictionary
        spec: Dict[str, Any] = {}

        # Handle contig parameter (residues to fix from input structure)
        # For partial diffusion with metals/ligands, we need to auto-generate a proper contig
        # because RFD3 separates non-polymer residues into a different chain
        contig_to_use = contig
        if partial_t is not None and pdb_content and not contig:
            # Auto-generate contig from PDB, excluding metals/ligands
            auto_contig = extract_protein_contig_from_pdb(pdb_content, ligand)
            if auto_contig:
                contig_to_use = auto_contig
                print(f"[RFD3] Auto-generated contig for partial diffusion: {auto_contig}")

        if contig_to_use:
            contig_str = contig_to_use
            is_simple_length = contig_str.replace("-", "").isdigit()
            if is_simple_length:
                # If contig looks like a length (e.g., "100" or "40-60"), treat as length
                if "-" in contig_str:
                    spec["length"] = contig_str
                else:
                    spec["length"] = int(contig_str)
            else:
                # Real contig with chain info (e.g., "A2-52")
                # For partial diffusion, validate that the contig doesn't have gaps from metals
                if partial_t is not None and pdb_content:
                    # Re-extract proper contig to handle any gaps from non-protein residues
                    proper_contig = extract_protein_contig_from_pdb(pdb_content, ligand)
                    if proper_contig and proper_contig != contig_str:
                        print(f"[RFD3] Corrected contig for partial diffusion: {contig_str} -> {proper_contig}")
                        contig_str = proper_contig
                spec["contig"] = contig_str

        # Handle length parameter (de novo chain to design)
        # Can be used TOGETHER with contig for binder design
        if length:
            length_str = str(length)
            if "-" in length_str:
                spec["length"] = length_str
            elif length_str.isdigit():
                spec["length"] = int(length_str)
            else:
                spec["length"] = length_str

        # Ligand for small molecule / enzyme design
        # Priority: ligand_sdf > ligand_smiles > ligand (code)
        temp_ligand_path = None
        apply_symmetry_post = None  # For post-processing symmetry (deprecated - use native symmetry)
        force_symmetry_cli = False  # Force symmetry CLI flag for interface designs
        if ligand_sdf:
            # Write SDF to temp file for RFD3
            with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False, mode='w') as f:
                f.write(ligand_sdf)
                temp_ligand_path = f.name
            spec["ligand"] = temp_ligand_path
            print(f"[RFD3] Using SDF ligand from temp file: {temp_ligand_path}")
        elif ligand_smiles:
            # Generate 3D conformer from SMILES and embed in input PDB
            # Use L:XXX naming convention to avoid CCD conflicts (per RFD3 docs)
            from conformer_utils import generate_conformer, generate_conformer_oriented, ConformerMethod

            # Determine conformer method
            method = ConformerMethod.RDKIT
            if conformer_method:
                try:
                    method = ConformerMethod(conformer_method.lower())
                except ValueError:
                    print(f"[RFD3] Unknown conformer method '{conformer_method}', using rdkit")

            # Use provided center or default to origin
            center = ligand_center or (0.0, 0.0, 0.0)

            # Use L:X residue name pattern to avoid CCD conflicts
            # RFD3 docs: "We suggest renaming all custom ligands to begin with L: to avoid all clashes with the CCD"
            ligand_res_name = "UNL"  # "Unknown Ligand" - a valid 3-char code

            # Check if interface_ligand mode is requested with symmetry
            sym_id = None
            if symmetry:
                sym_id = symmetry.get("id") if isinstance(symmetry, dict) else symmetry

            if interface_ligand and sym_id:
                # Interface ligand mode: Orient ligand for symmetric interface
                print(f"[RFD3] Interface ligand mode: orienting for {sym_id} symmetry...")
                ligand_pdb = generate_conformer_oriented(
                    smiles=ligand_smiles,
                    method=method,
                    name=ligand_res_name,
                    center=center,
                    symmetry=sym_id,
                    axis_atoms=None,  # Auto-detect (N=N for azo, or longest axis)
                    optimize=True,
                    fallback=True
                )

                # Add RASA conditioning to keep ligand at interface (partially buried)
                # BUT only if we have an existing input structure - for de novo design,
                # there are no residues to select from, so we rely on ligand orientation only
                if pdb_content:
                    if select_partially_buried is None:
                        select_partially_buried = {}
                    # L1 = ligand residue 1 in chain L
                    select_partially_buried["L1"] = "ALL"
                    print(f"[RFD3] Added RASA conditioning: select_partially_buried={{L1: ALL}}")
                else:
                    print(f"[RFD3] De novo design: using ligand orientation only (no RASA for de novo)")
            else:
                # Standard mode: Generate conformer without orientation
                print(f"[RFD3] Generating 3D structure from SMILES using {method.value} method...")
                ligand_pdb = generate_conformer(
                    smiles=ligand_smiles,
                    method=method,
                    name=ligand_res_name,
                    center=center,
                    optimize=True,
                    fallback=True
                )

            if ligand_pdb:
                if pdb_content:
                    # Existing structure - embed ligand into it
                    pdb_content = pdb_content.replace("END\n", "").replace("END", "").rstrip()
                    pdb_content = pdb_content + "\n" + ligand_pdb
                    spec["ligand"] = ligand_res_name
                    print(f"[RFD3] Embedded ligand as residue {ligand_res_name}")
                elif interface_ligand and symmetry:
                    # De novo symmetric design with ligand at interface
                    #
                    # IMPORTANT: Native RFD3 symmetry does NOT work with ligand-only input
                    # because there are no protein entities to symmetrize.
                    # Instead, use POST-PROCESSING symmetry: design a monomer then apply C2 rotation.
                    #
                    # Check if binder mode is requested (uses hotspots + H-bond conditioning)
                    # interface_ligand can be True (default) or "binder" for binder mode
                    binder_mode = interface_ligand == "binder" if isinstance(interface_ligand, str) else False

                    sym_id = symmetry.get("id") if isinstance(symmetry, dict) else symmetry

                    # Ligand becomes the input structure
                    pdb_content = ligand_pdb
                    spec["ligand"] = ligand_res_name

                    # Center on ligand for proper geometry
                    spec["ori_token"] = [0.0, 0.0, 0.0]

                    # POST-PROCESSING SYMMETRY MODE:
                    # Design a ONE-SIDED binder, then apply C2 rotation for the other chain
                    #
                    # CRITICAL: The protein must only bind ONE SIDE of the ligand (180°, not 360°)
                    # Otherwise both chains wrap around the ligand and trap it inside.
                    #
                    # For azobenzene oriented along Z-axis:
                    # - C1-C6: first phenyl ring (negative Z side)
                    # - N7,N8: azo nitrogens (center, on C2 axis)
                    # - C9-C14: second phenyl ring (positive Z side)
                    #
                    # Strategy: Use select_exposed on one phenyl ring so protein only binds the other side
                    # Then C2 rotation creates a chain that binds the exposed side

                    # EXPOSE MORE than half the ligand to STRICTLY limit binding to one side
                    # The protein can ONLY contact C1-C6 (first phenyl ring)
                    # Everything else (N7,N8,C9-C14) must be exposed for the C2 partner
                    #
                    # This is critical: if we only expose C9-C14, the protein can still
                    # approach from the "top" or "bottom" of the ligand plane and wrap around.
                    # By also exposing N7,N8 (the azo nitrogens on the C2 axis), we ensure
                    # the protein only approaches from ONE side of the ligand.
                    spec["select_exposed"] = {ligand_res_name: "N7,N8,C9,C10,C11,C12,C13,C14"}

                    # Bury the binding side phenyl carbons to ensure strong binding
                    spec["select_buried"] = {ligand_res_name: "C1,C2,C3,C4,C5,C6"}

                    if binder_mode:
                        # BINDER MODE: Add H-bond conditioning for the azo nitrogens
                        # After C2 rotation, the partner chain will H-bond with N7,N8
                        # Since N7,N8 are exposed for THIS chain, H-bond goes to partner
                        print(f"[RFD3] Interface ligand mode: STRICT HALF binder + POST-PROCESSING {sym_id}")
                        print(f"[RFD3] Exposed: N7,N8,C9-C14 (8 atoms); Buried: C1-C6 (6 atoms)")
                    else:
                        print(f"[RFD3] Interface ligand mode: STRICT HALF binder + POST-PROCESSING {sym_id}")
                        print(f"[RFD3] Exposed: N7,N8,C9-C14; Buried: C1-C6")

                    # Use POST-PROCESSING symmetry - native doesn't work with ligand-only input
                    apply_symmetry_post = sym_id
                    force_symmetry_cli = False
                else:
                    # Standard de novo design - ligand becomes input structure
                    pdb_content = ligand_pdb
                    spec["ligand"] = ligand_res_name
                    print(f"[RFD3] Embedded ligand as residue {ligand_res_name}")
            else:
                raise ValueError(f"Failed to generate 3D structure from SMILES: {ligand_smiles}")
        elif ligand:
            spec["ligand"] = ligand
            # For metal binding redesign: replace original metal with target metal in PDB
            if pdb_content and ligand in METAL_CODES:
                pdb_content = replace_metal_in_pdb(pdb_content, target_metal=ligand)

        # Unindex for enzyme design (residues with inferred positions)
        if unindex:
            spec["unindex"] = unindex

        # Symmetry configuration - RFD3 expects SymmetryConfig dict with id, is_unsym_motif, is_symmetric_motif
        # Skip if using post-processing symmetry (apply_symmetry_post is set)
        if symmetry and apply_symmetry_post is None:
            sym_id = symmetry.get("id") if isinstance(symmetry, dict) else symmetry
            if sym_id:
                # Build SymmetryConfig dict with all fields
                sym_config = {"id": sym_id}
                if isinstance(symmetry, dict):
                    if symmetry.get("is_symmetric_motif") is not None:
                        sym_config["is_symmetric_motif"] = symmetry.get("is_symmetric_motif")
                    if symmetry.get("is_unsym_motif"):
                        sym_config["is_unsym_motif"] = symmetry["is_unsym_motif"]
                spec["symmetry"] = sym_config
        elif apply_symmetry_post:
            print(f"[RFD3] Skipping native symmetry - will apply {apply_symmetry_post} in post-processing")

        # is_non_loopy goes in SPEC (not top-level)
        if is_non_loopy is not None:
            spec["is_non_loopy"] = is_non_loopy

        # partial_t goes in SPEC (not top-level)
        if partial_t is not None:
            spec["partial_t"] = partial_t

        # Hotspots -> select_hotspots in SPEC
        if hotspots:
            # Convert list to dict format: ["A15", "A20"] -> {"A15": "ALL", "A20": "ALL"}
            spec["select_hotspots"] = {h: "ALL" for h in hotspots}

        # Origin and strategy go in SPEC
        if ori_token:
            spec["ori_token"] = ori_token
        if infer_ori_strategy:
            spec["infer_ori_strategy"] = infer_ori_strategy

        # RASA conditioning goes in SPEC
        if select_buried:
            spec["select_buried"] = select_buried
        if select_exposed:
            spec["select_exposed"] = select_exposed
        if select_partially_buried:
            spec["select_partially_buried"] = select_partially_buried

        # Fixed atoms go in SPEC
        if select_fixed_atoms:
            spec["select_fixed_atoms"] = select_fixed_atoms

        # H-bond conditioning goes in SPEC
        if select_hbond_donor:
            spec["select_hbond_donor"] = select_hbond_donor
        if select_hbond_acceptor:
            spec["select_hbond_acceptor"] = select_hbond_acceptor

        # Covalent bonds for enzyme design (struct_conn format)
        # Format: [{"protein": {"chain": "A", "res_name": "CYS", "res_num": 145, "atom_name": "SG"},
        #          "ligand": {"chain": "L", "res_name": "LIG", "res_num": 1, "atom_name": "C1"}}]
        if covalent_bonds and len(covalent_bonds) > 0:
            # Convert to struct_conn format for RFD3
            struct_conn = []
            for bond in covalent_bonds:
                protein = bond.get("protein", {})
                ligand = bond.get("ligand", {})
                # Format: "chain/resName/resNum/atomName"
                protein_spec = f"{protein.get('chain', 'A')}/{protein.get('res_name', '')}/{protein.get('res_num', 0)}/{protein.get('atom_name', '')}"
                ligand_spec = f"{ligand.get('chain', 'L')}/{ligand.get('res_name', '')}/{ligand.get('res_num', 1)}/{ligand.get('atom_name', '')}"
                struct_conn.append({
                    "partner_1": protein_spec,
                    "partner_2": ligand_spec,
                    "conn_type": "covale"  # Covalent bond type
                })
            spec["struct_conn"] = struct_conn
            print(f"[RFD3] Added {len(struct_conn)} covalent bond(s) to spec")

        # Build inference sampler config for DIFFUSION parameters ONLY
        sampler_config: Dict[str, Any] = {}
        if num_timesteps is not None:
            sampler_config["num_timesteps"] = num_timesteps
        if step_scale is not None:
            sampler_config["step_scale"] = step_scale
        if gamma_0 is not None:
            sampler_config["gamma_0"] = gamma_0

        # Symmetry requires special sampler kind
        # But NOT when using post-processing symmetry (apply_symmetry_post is set)
        if (symmetry or force_symmetry_cli) and apply_symmetry_post is None:
            sampler_config["kind"] = "symmetry"
            print(f"[RFD3] Using symmetry sampler (kind=symmetry)")
        elif apply_symmetry_post:
            print(f"[RFD3] Using standard sampler (post-processing {apply_symmetry_post} symmetry)")

        # Build config kwargs - ONLY these top-level params are allowed!
        config_kwargs: Dict[str, Any] = {
            "specification": spec,
            "diffusion_batch_size": num_designs,
        }

        # Add sampler config if any diffusion parameters were set
        if sampler_config:
            config_kwargs["inference_sampler"] = sampler_config

        # High-order symmetry requires low memory mode (only for native symmetry)
        if symmetry and apply_symmetry_post is None:
            sym_id = symmetry.get("id") if isinstance(symmetry, dict) else symmetry
            if sym_id in ["T", "O", "I"]:
                config_kwargs["low_memory_mode"] = True

        # Handle PDB input - write to temp file and pass path via spec["input"]
        temp_pdb_path = None
        if pdb_content:
            # Write PDB to temp file - RFD3 expects file path, not AtomArray
            with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode='w') as f:
                f.write(pdb_content)
                temp_pdb_path = f.name
            spec["input"] = temp_pdb_path
            print(f"[RFD3] Wrote input structure to {temp_pdb_path}")

        print(f"[RFD3] Spec: {json.dumps({k: str(v)[:100] for k, v in spec.items()}, indent=2)}")
        print(f"[RFD3] Config: {json.dumps({k: str(v)[:100] for k, v in config_kwargs.items()}, indent=2)}")

        # Create config and run
        config = RFD3InferenceConfig(**config_kwargs)
        model = RFD3InferenceEngine(**config)

        # Run inference - inputs are passed via spec["input"], not model.run()
        outputs_dict = model.run(inputs=None, out_dir=None, n_batches=1)

        # Process outputs
        designs = []
        if outputs_dict:
            for key, result_list in outputs_dict.items():
                items = result_list if isinstance(result_list, list) else [result_list]

                for idx, item in enumerate(items):
                    if hasattr(item, 'atom_array'):
                        try:
                            atom_array = item.atom_array

                            # Apply post-processing symmetry if needed
                            if apply_symmetry_post:
                                atom_array = _apply_post_symmetry(atom_array, apply_symmetry_post)
                                print(f"[RFD3] Applied {apply_symmetry_post} symmetry in post-processing")

                            output_pdb = atom_array_to_pdb(atom_array)
                            filename = f"{key}_{idx}.pdb" if len(items) > 1 else f"{key}.pdb"

                            design: Dict[str, Any] = {"filename": filename, "content": output_pdb}

                            try:
                                design["cif_content"] = atom_array_to_cif(atom_array)
                            except Exception:
                                pass

                            designs.append(design)
                        except Exception as e:
                            print(f"[RFD3] Error converting output: {e}")

        # Cleanup temp files
        if temp_pdb_path and os.path.exists(temp_pdb_path):
            os.unlink(temp_pdb_path)
        if temp_ligand_path and os.path.exists(temp_ligand_path):
            os.unlink(temp_ligand_path)

        return {
            "status": "completed",
            "result": {
                "designs": designs,
                "mode": "real",
                "seed": seed,
                "ligand_used": ligand_smiles or ligand or None,
            }
        }

    except Exception as e:
        import traceback
        # Cleanup temp files on error
        if 'temp_pdb_path' in locals() and temp_pdb_path and os.path.exists(temp_pdb_path):
            os.unlink(temp_pdb_path)
        if 'temp_ligand_path' in locals() and temp_ligand_path and os.path.exists(temp_ligand_path):
            os.unlink(temp_ligand_path)
        return {
            "status": "failed",
            "error": str(e),
            "traceback": traceback.format_exc()
        }


def run_rfd3_mock(contig: str, num_designs: int) -> Dict[str, Any]:
    """Mock RFD3 inference"""
    try:
        length = int(contig) if contig.isdigit() else 100
    except Exception:
        length = 100

    designs = []
    for i in range(num_designs):
        designs.append({
            "filename": f"design_{i+1}.pdb",
            "content": generate_mock_pdb(length)
        })

    return {
        "status": "completed",
        "result": {
            "designs": designs,
            "mode": "mock",
        }
    }


# ============== RF3 Inference ==============

def ensure_project_root():
    """Ensure .project-root file exists for rootutils (needed by rf3 CLI)"""
    locations = [
        "/workspace/.project-root",
        "/.project-root",
        os.path.expanduser("~/.project-root"),
    ]

    # Try to find rf3 package location
    try:
        import rf3
        rf3_pkg_dir = os.path.dirname(rf3.__file__)
        locations.extend([
            os.path.join(rf3_pkg_dir, ".project-root"),
            os.path.join(os.path.dirname(rf3_pkg_dir), ".project-root"),
            os.path.join(os.path.dirname(os.path.dirname(rf3_pkg_dir)), ".project-root"),
        ])
    except ImportError:
        pass

    for loc in locations:
        try:
            parent_dir = os.path.dirname(loc)
            if parent_dir and not os.path.exists(parent_dir):
                os.makedirs(parent_dir, exist_ok=True)
            if not os.path.exists(loc):
                with open(loc, "w") as f:
                    f.write("")
                print(f"[RF3] Created .project-root at {loc}")
        except Exception:
            pass  # Silently continue


def run_rf3_inference(
    sequence: str,
    name: str = "prediction",
    pdb_content: Optional[str] = None,
    msa_content: Optional[str] = None,
    sequences: Optional[List[str]] = None,  # Multi-chain: additional chain sequences
    ligand_smiles: Optional[str] = None,     # Small molecule SMILES for protein-ligand prediction
    use_mock: bool = False
) -> Dict[str, Any]:
    """Run RF3 structure prediction

    Args:
        sequence: Protein sequence for chain A (required)
        name: Name for the prediction output
        pdb_content: Optional PDB content for structure-based input
        msa_content: Optional MSA content (.a3m or .fasta format) for improved predictions
        sequences: Optional list of additional chain sequences [chain_B, chain_C, ...]
                   For dimer evaluation, pass the second chain here
        ligand_smiles: Optional ligand SMILES for protein-ligand binding evaluation
                       This enables ipTM scoring for the protein-ligand interface
        use_mock: Whether to use mock mode

    Returns:
        Dict with predictions and confidences including:
        - plddt_per_residue: Per-residue confidence
        - mean_plddt: Average pLDDT
        - ptm: Predicted TM score
        - ipTM: Interface pTM (only for multi-chain or protein-ligand)
        - chain_pair_pae: PAE between chain pairs
    """

    if use_mock:
        return run_rf3_mock(sequence, name)

    # Ensure .project-root exists before any rf3 imports (fixes rootutils issue)
    ensure_project_root()

    try:
        # Try Python API first (supports MSA, multi-chain, ligand)
        return run_rf3_python_api(sequence, name, pdb_content, msa_content, sequences, ligand_smiles)
    except ImportError as e:
        print(f"[RF3] Python API not available: {e}")
        # Fallback to CLI with multi-chain and ligand support
        return run_rf3_cli(sequence, name, msa_content, sequences, ligand_smiles)
    except Exception as e:
        import traceback
        return {
            "status": "failed",
            "error": str(e),
            "traceback": traceback.format_exc()
        }


def run_rf3_python_api(
    sequence: str,
    name: str,
    pdb_content: Optional[str],
    msa_content: Optional[str] = None,
    sequences: Optional[List[str]] = None,
    ligand_smiles: Optional[str] = None
) -> Dict[str, Any]:
    """RF3 via Python API with optional MSA support"""
    from rf3.inference_engines.rf3 import RF3InferenceEngine
    from rf3.utils.inference import InferenceInput

    # For sequence-only without MSA, use CLI approach (supports multi-chain/ligand)
    if not pdb_content and not msa_content:
        return run_rf3_cli(sequence, name, msa_content, sequences, ligand_smiles)

    # Initialize engine
    inference_engine = RF3InferenceEngine(ckpt_path='rf3', verbose=False)

    # Handle MSA file if provided
    msa_path = None
    if msa_content:
        # Write MSA to temp file
        with tempfile.NamedTemporaryFile(suffix=".a3m", delete=False, mode='w') as f:
            f.write(msa_content)
            msa_path = f.name
        print(f"[RF3] Wrote MSA to {msa_path}")

    try:
        if pdb_content:
            # Structure-based input
            atom_array = pdb_to_atom_array(pdb_content)
            input_data = InferenceInput.from_atom_array(atom_array, example_id=name)
            if msa_path:
                input_data.msa_path = msa_path
        else:
            # Sequence-only with MSA - use JSON config approach
            return run_rf3_cli(sequence, name, msa_content)

        # Run inference
        rf3_outputs = inference_engine.run(inputs=input_data)

        # Process outputs
        predictions = []
        confidences = None

        if name in rf3_outputs:
            results = rf3_outputs[name]

            for idx, rf3_output in enumerate(results):
                if hasattr(rf3_output, 'atom_array'):
                    pdb_out = atom_array_to_pdb(rf3_output.atom_array)

                    pred = {
                        "filename": f"{name}_{idx}.pdb" if len(results) > 1 else f"{name}.pdb",
                        "content": pdb_out,
                    }

                    try:
                        pred["cif_content"] = atom_array_to_cif(rf3_output.atom_array)
                    except Exception:
                        pass

                    predictions.append(pred)

                # Extract confidences from first output
                if idx == 0 and confidences is None:
                    confidences = extract_rf3_confidences(rf3_output)

        return {
            "status": "completed",
            "result": {
                "predictions": predictions,
                "confidences": confidences,
                "mode": "real",
            }
        }
    finally:
        # Cleanup temp MSA file
        if msa_path and os.path.exists(msa_path):
            os.unlink(msa_path)


def run_rf3_cli(
    sequence: str,
    name: str,
    msa_content: Optional[str] = None,
    sequences: Optional[List[str]] = None,  # Multi-chain support
    ligand_smiles: Optional[str] = None,     # Small molecule support
) -> Dict[str, Any]:
    """
    RF3 via CLI for sequence-only input with optional MSA support.

    Args:
        sequence: Primary sequence (chain A)
        name: Job name
        msa_content: Optional MSA content
        sequences: Optional list of additional chain sequences [chain_B, chain_C, ...]
        ligand_smiles: Optional ligand SMILES string for protein-ligand prediction

    Returns:
        Prediction results with confidences including ipTM for interfaces
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # Write MSA to temp file if provided
        msa_path = None
        if msa_content:
            msa_path = os.path.join(tmpdir, "msa.a3m")
            with open(msa_path, "w") as f:
                f.write(msa_content)
            print(f"[RF3] Wrote MSA to {msa_path}")

        # Create JSON config with components
        config_path = os.path.join(tmpdir, "input.json")
        components = []

        # Add primary chain A (RF3 uses "seq" not "sequence")
        component_a = {"seq": sequence, "chain_id": "A"}
        if msa_path:
            component_a["msa_path"] = msa_path
        components.append(component_a)

        # Add additional chains if provided (for dimer/complex evaluation)
        if sequences:
            chain_ids = "BCDEFGHIJKLMNOPQRSTUVWXYZ"
            for i, seq in enumerate(sequences):
                if i < len(chain_ids):
                    components.append({
                        "seq": seq,
                        "chain_id": chain_ids[i]
                    })

        # Add ligand if provided (RF3 uses smiles directly)
        if ligand_smiles:
            components.append({
                "smiles": ligand_smiles
            })

        # RF3 expects input as a list of prediction jobs
        json_content = [{
            "name": name,
            "components": components
        }]
        with open(config_path, "w") as f:
            json.dump(json_content, f)
        print(f"[RF3] JSON config: {json.dumps(json_content, indent=2)}")

        out_dir = os.path.join(tmpdir, "output")
        os.makedirs(out_dir, exist_ok=True)

        # Create .project-root in temp directory (rootutils searches upward from cwd)
        project_root_file = os.path.join(tmpdir, ".project-root")
        with open(project_root_file, "w") as f:
            f.write("")

        # Set environment variable as alternative for rootutils
        env = os.environ.copy()
        env["PROJECT_ROOT"] = tmpdir

        # Get checkpoint directory from environment
        checkpoint_dir = os.environ.get("FOUNDRY_CHECKPOINT_DIRS", "/runpod-volume/checkpoints")
        rf3_ckpt = os.path.join(checkpoint_dir, "rf3_foundry_01_24_latest.ckpt")

        # Run CLI from temp directory where .project-root exists
        # Pass checkpoint path explicitly
        cmd = f"rf3 predict out_dir={out_dir} inputs={config_path} ckpt_path={rf3_ckpt}"

        result = subprocess.run(
            cmd, shell=True, capture_output=True, text=True,
            timeout=3600, cwd=tmpdir, env=env
        )

        if result.returncode == 0:
            predictions = []
            confidences = {}
            # Find output files
            for root, _, files in os.walk(out_dir):
                for filename in files:
                    filepath = os.path.join(root, filename)
                    if filename.endswith((".pdb", ".cif")):
                        with open(filepath) as f:
                            content = f.read()
                            predictions.append({"filename": filename, "content": content})
                            # Extract pLDDT from CIF B-factors
                            if filename.endswith(".cif"):
                                plddt_values = extract_plddt_from_cif(content)
                                if plddt_values:
                                    confidences["plddt_per_residue"] = plddt_values
                                    confidences["mean_plddt"] = sum(plddt_values) / len(plddt_values)
                    elif filename.endswith(".json"):
                        # Check for confidence JSON files
                        try:
                            with open(filepath) as f:
                                conf_data = json.load(f)
                                if "plddt" in conf_data or "ptm" in conf_data:
                                    confidences.update(conf_data)
                        except:
                            pass

            return {
                "status": "completed",
                "result": {
                    "predictions": predictions,
                    "confidences": confidences if confidences else None,
                    "mode": "real (cli)",
                }
            }
        else:
            return {
                "status": "failed",
                "error": result.stderr or result.stdout or "Unknown error"
            }


def run_rf3_mock(sequence: str, name: str) -> Dict[str, Any]:
    """Mock RF3 inference"""
    length = len(sequence)
    pdb_content = generate_mock_pdb(length)
    confidences = generate_mock_confidences(length)

    return {
        "status": "completed",
        "result": {
            "predictions": [{
                "filename": f"{name}.pdb",
                "content": pdb_content,
            }],
            "confidences": confidences,
            "mode": "mock",
        }
    }


def extract_rf3_confidences(rf3_output) -> Optional[Dict[str, Any]]:
    """Extract confidence metrics from RF3 output"""
    confidences = {}

    if hasattr(rf3_output, 'summary_confidences') and rf3_output.summary_confidences:
        sc = rf3_output.summary_confidences
        confidences["summary_confidences"] = {
            "overall_plddt": float(sc.get('overall_plddt', 0)),
            "overall_pae": float(sc.get('overall_pae', 0)),
            "overall_pde": float(sc.get('overall_pde', 0)) if 'overall_pde' in sc else None,
            "ptm": float(sc.get('ptm', 0)),
            "iptm": float(sc.get('iptm')) if sc.get('iptm') is not None else None,
            "ranking_score": float(sc.get('ranking_score', 0)),
            "has_clash": bool(sc.get('has_clash', False)),
        }

    if hasattr(rf3_output, 'confidences') and rf3_output.confidences:
        conf = rf3_output.confidences
        if 'atom_plddts' in conf:
            atom_plddts = list(conf['atom_plddts'])
            confidences["per_residue_plddt"] = [round(float(p), 3) for p in atom_plddts[:500]]

        if 'pae' in conf:
            pae = conf['pae']
            pae_list = [[round(float(x), 2) for x in row[:100]] for row in pae[:100]]
            confidences["pae_matrix"] = pae_list

    return confidences if confidences else None


# ============== MPNN Inference ==============

def run_mpnn_inference(
    pdb_content: str,
    num_sequences: int = 8,
    temperature: float = 0.1,
    model_type: str = "ligand_mpnn",
    remove_waters: bool = True,
    fixed_positions: Optional[List[str]] = None,  # e.g., ["A35", "A36", "B35"]
    use_mock: bool = False,
    # Advanced LigandMPNN parameters (from Nature Methods 2025 paper)
    pack_side_chains: bool = False,
    pack_with_ligand_context: bool = True,
    number_of_packs_per_design: int = 4,
    bias_AA: Optional[str] = None,  # e.g., "W:3.0,P:-3.0,C:-5.0"
    omit_AA: Optional[str] = None,  # e.g., "C" to omit cysteine
    model_noise_level: str = "010",  # 005, 010, 020, 030
    ligand_cutoff_for_score: float = 8.0,  # Angstroms
    use_side_chain_context: bool = False,
    save_stats: bool = False,
) -> Dict[str, Any]:
    """Run MPNN sequence design with full LigandMPNN capabilities.

    Args:
        pdb_content: PDB content including HETATM records for ligand awareness
        num_sequences: Number of sequences to design
        temperature: Sampling temperature (lower = more conservative, 0.05-0.1 recommended)
        model_type: 'ligand_mpnn' for ligand-aware design, 'proteinmpnn' otherwise
        remove_waters: Whether to remove water molecules
        fixed_positions: List of residue positions to keep fixed (e.g., ["A35", "A36"])
        use_mock: Use mock mode for testing

        Advanced LigandMPNN parameters:
        pack_side_chains: Enable sidechain packing (generates chi angles)
        pack_with_ligand_context: Include ligand atoms when packing sidechains
        number_of_packs_per_design: Number of sidechain packing samples per sequence
        bias_AA: Global amino acid biases (e.g., "W:3.0,Y:2.0,C:-5.0")
        omit_AA: Amino acids to completely omit (e.g., "C" to avoid cysteines)
        model_noise_level: Model noise level (005=low, 010=default, 020, 030=high)
        ligand_cutoff_for_score: Distance cutoff for ligand-adjacent residues (Angstroms)
        use_side_chain_context: Use fixed residue sidechains as additional context
        save_stats: Return confidence metrics (ligand_confidence, overall_confidence)

    Returns:
        Dict with designed sequences and optional confidence metrics

    Note: For best ligand binding design:
    - Use model_type='ligand_mpnn' with pack_side_chains=True
    - Bias toward aromatic residues for hydrophobic ligands: bias_AA="W:2.0,Y:2.0,F:1.5"
    - Omit cysteines to avoid disulfide complications: omit_AA="C"
    - Use ligand_cutoff_for_score=6.0 for small molecules (tighter contact)
    """

    if use_mock:
        return run_mpnn_mock(pdb_content, num_sequences)

    # Build kwargs dict for the advanced parameters
    advanced_params = {
        "pack_side_chains": pack_side_chains,
        "pack_with_ligand_context": pack_with_ligand_context,
        "number_of_packs_per_design": number_of_packs_per_design,
        "bias_AA": bias_AA,
        "omit_AA": omit_AA,
        "model_noise_level": model_noise_level,
        "ligand_cutoff_for_score": ligand_cutoff_for_score,
        "use_side_chain_context": use_side_chain_context,
        "save_stats": save_stats,
    }

    try:
        return run_mpnn_python_api(
            pdb_content, num_sequences, temperature, model_type, remove_waters, fixed_positions,
            **advanced_params
        )
    except ImportError as e:
        print(f"[MPNN] Python API not available: {e}")
        return run_mpnn_cli(pdb_content, num_sequences, temperature, model_type, remove_waters, fixed_positions,
                          **advanced_params)
    except Exception as e:
        import traceback
        # Try CLI fallback
        try:
            return run_mpnn_cli(pdb_content, num_sequences, temperature, model_type, remove_waters, fixed_positions,
                              **advanced_params)
        except Exception as cli_error:
            return {
                "status": "failed",
                "error": f"Python API: {str(e)}\nCLI fallback: {str(cli_error)}",
                "traceback": traceback.format_exc()
            }


def run_mpnn_python_api(
    pdb_content: str,
    num_sequences: int,
    temperature: float,
    model_type: str,
    remove_waters: bool,
    fixed_positions: Optional[List[str]] = None,
    # Advanced LigandMPNN parameters
    pack_side_chains: bool = False,
    pack_with_ligand_context: bool = True,
    number_of_packs_per_design: int = 4,
    bias_AA: Optional[str] = None,
    omit_AA: Optional[str] = None,
    model_noise_level: str = "010",
    ligand_cutoff_for_score: float = 8.0,
    use_side_chain_context: bool = False,
    save_stats: bool = False,
) -> Dict[str, Any]:
    """MPNN via Python API with full LigandMPNN capabilities."""
    from mpnn.inference_engines.mpnn import MPNNInferenceEngine
    from biotite.structure import get_residue_starts
    from biotite.sequence import ProteinSequence

    # Parse input
    atom_array = pdb_to_atom_array(pdb_content)

    # Configure engine with model type
    # Note: model_type must be "ligand_mpnn" or "protein_mpnn"
    # Noise level is controlled via checkpoint_path (e.g., ligandmpnn_v_32_010_25.pt)
    checkpoint_path = None
    if model_noise_level != "010":
        # Non-default noise level - construct checkpoint path
        # Format: /runpod-volume/checkpoints/ligandmpnn_v_32_{noise}_25.pt
        if model_type == "ligand_mpnn":
            checkpoint_path = f"/runpod-volume/checkpoints/ligandmpnn_v_32_{model_noise_level}_25.pt"
        else:
            checkpoint_path = f"/runpod-volume/checkpoints/proteinmpnn_v_48_{model_noise_level}_25.pt"
        print(f"[MPNN] Using checkpoint: {checkpoint_path}")

    engine_config = {
        "model_type": model_type,  # Must be "ligand_mpnn" or "protein_mpnn"
        "is_legacy_weights": True,
        "out_directory": None,
        "write_structures": pack_side_chains,  # Write structures if packing sidechains
        "write_fasta": False,
    }

    # Add checkpoint path if non-default noise level
    if checkpoint_path:
        engine_config["checkpoint_path"] = checkpoint_path

    input_configs = [{
        "batch_size": num_sequences,
        "remove_waters": remove_waters,
        "temperature": temperature,
    }]

    # Add fixed positions if specified
    if fixed_positions:
        input_configs[0]["fixed_positions"] = fixed_positions
        print(f"[MPNN] Fixed positions: {fixed_positions}")

    # Single-letter to 3-letter amino acid code conversion
    AA_1TO3 = {
        "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
        "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
        "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
        "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
    }

    # Add amino acid biasing (uses 3-letter codes)
    if bias_AA:
        # Parse bias string like "W:3.0,Y:2.0,C:-5.0" into dict with 3-letter codes
        bias_dict = {}
        for item in bias_AA.split(","):
            aa, val = item.strip().split(":")
            aa1 = aa.strip().upper()
            aa3 = AA_1TO3.get(aa1, aa1)  # Convert to 3-letter if single-letter
            bias_dict[aa3] = float(val.strip())
        input_configs[0]["bias"] = bias_dict  # Note: key is "bias", not "bias_AA"
        print(f"[MPNN] AA bias: {bias_dict}")

    # Add amino acid omission (uses 3-letter codes)
    if omit_AA:
        # Convert single-letter codes to 3-letter codes
        omit_list = []
        for aa in omit_AA.split(","):
            aa1 = aa.strip().upper()
            aa3 = AA_1TO3.get(aa1, aa1)
            omit_list.append(aa3)
        input_configs[0]["omit"] = omit_list  # Note: key is "omit", not "omit_AA"
        print(f"[MPNN] Omitting AAs: {omit_list}")

    # LigandMPNN-specific parameters
    if model_type == "ligand_mpnn":
        input_configs[0]["ligand_mpnn_use_atom_context"] = 1  # Enable ligand awareness
        input_configs[0]["ligand_mpnn_cutoff_for_score"] = ligand_cutoff_for_score
        if use_side_chain_context:
            input_configs[0]["ligand_mpnn_use_side_chain_context"] = 1
        print(f"[MPNN] Ligand cutoff for score: {ligand_cutoff_for_score}Å")

    # Sidechain packing parameters
    if pack_side_chains:
        input_configs[0]["pack_side_chains"] = 1
        input_configs[0]["number_of_packs_per_design"] = number_of_packs_per_design
        if pack_with_ligand_context and model_type == "ligand_mpnn":
            input_configs[0]["pack_with_ligand_context"] = 1
        print(f"[MPNN] Sidechain packing enabled: {number_of_packs_per_design} packs/design")

    # Run inference
    model = MPNNInferenceEngine(**engine_config)
    mpnn_outputs = model.run(input_dicts=input_configs, atom_arrays=[atom_array])

    # Extract sequences and stats
    sequences = []
    stats_list = []
    pdb_structures = []

    for i, item in enumerate(mpnn_outputs):
        if hasattr(item, 'atom_array'):
            res_starts = get_residue_starts(item.atom_array)
            seq = ''.join(
                ProteinSequence.convert_letter_3to1(res_name)
                for res_name in item.atom_array.res_name[res_starts]
                if res_name in ProteinSequence._dict_3to1
            )
            sequences.append(f">design_{i+1}\n{seq}")

            # Extract confidence/stats if available
            if save_stats and hasattr(item, 'stats'):
                stats_list.append({
                    "design": i + 1,
                    "overall_confidence": getattr(item.stats, 'overall_confidence', None),
                    "ligand_confidence": getattr(item.stats, 'ligand_confidence', None),
                    "score": getattr(item.stats, 'score', None),
                })

            # Extract packed structure if available
            if pack_side_chains and hasattr(item, 'pdb_content'):
                pdb_structures.append({
                    "design": i + 1,
                    "pdb_content": item.pdb_content,
                })

    fasta_content = "\n".join(sequences)

    result = {
        "status": "completed",
        "result": {
            "sequences": [{"filename": "sequences.fasta", "content": fasta_content}],
            "num_sequences": len(sequences),
            "model_type": model_type,
            "model_noise_level": model_noise_level,
            "mode": "real",
        }
    }

    # Add stats if requested
    if save_stats and stats_list:
        result["result"]["stats"] = stats_list
        # Compute averages
        confidences = [s["overall_confidence"] for s in stats_list if s["overall_confidence"] is not None]
        ligand_confs = [s["ligand_confidence"] for s in stats_list if s["ligand_confidence"] is not None]
        if confidences:
            result["result"]["avg_overall_confidence"] = sum(confidences) / len(confidences)
        if ligand_confs:
            result["result"]["avg_ligand_confidence"] = sum(ligand_confs) / len(ligand_confs)

    # Add packed structures if generated
    if pack_side_chains and pdb_structures:
        result["result"]["packed_structures"] = pdb_structures

    return result


def run_mpnn_cli(
    pdb_content: str,
    num_sequences: int,
    temperature: float,
    model_type: str,
    remove_waters: bool,
    fixed_positions: Optional[List[str]] = None,
    # Advanced LigandMPNN parameters
    pack_side_chains: bool = False,
    pack_with_ligand_context: bool = True,
    number_of_packs_per_design: int = 4,
    bias_AA: Optional[str] = None,
    omit_AA: Optional[str] = None,
    model_noise_level: str = "010",
    ligand_cutoff_for_score: float = 8.0,
    use_side_chain_context: bool = False,
    save_stats: bool = False,
) -> Dict[str, Any]:
    """MPNN via CLI with full LigandMPNN capabilities."""
    with tempfile.TemporaryDirectory() as tmpdir:
        pdb_path = os.path.join(tmpdir, "input.pdb")
        with open(pdb_path, "w") as f:
            f.write(pdb_content)

        out_dir = os.path.join(tmpdir, "output")
        os.makedirs(out_dir, exist_ok=True)

        # Build command - model_type must be "ligand_mpnn" or "protein_mpnn"
        # Checkpoint path is required by the CLI
        if model_type == "ligand_mpnn":
            checkpoint = f"/runpod-volume/checkpoints/ligandmpnn_v_32_{model_noise_level}_25.pt"
        else:
            checkpoint = f"/runpod-volume/checkpoints/proteinmpnn_v_48_{model_noise_level}_25.pt"

        cmd_parts = [
            "mpnn",
            f"--structure_path {pdb_path}",
            f"--out_directory {out_dir}",
            "--name design",
            f"--model_type {model_type}",  # Must be "ligand_mpnn" or "protein_mpnn"
            f"--checkpoint_path {checkpoint}",  # Required by CLI
            f"--batch_size {num_sequences}",
            f"--temperature {temperature}",
            "--is_legacy_weights True",
            "--write_fasta True",
            f"--write_structures {'True' if pack_side_chains else 'False'}",
            f"--remove_waters {'True' if remove_waters else 'None'}",
        ]
        print(f"[MPNN] Using checkpoint: {checkpoint}")

        # Add fixed positions if specified (atomworks CLI uses --fixed_residues)
        if fixed_positions:
            fixed_str = ",".join(fixed_positions)
            cmd_parts.append(f"--fixed_residues {fixed_str}")
            print(f"[MPNN] CLI fixed_residues: {fixed_str}")

        # Add amino acid biasing (atomworks CLI uses --bias with JSON format)
        # Format: {"TRP": 3.0, "TYR": 2.0} - uses 3-letter codes
        # Single-letter to 3-letter conversion
        AA_1TO3 = {
            "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
            "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
            "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
            "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
        }
        if bias_AA:
            # Convert "W:3.0,Y:2.0" to JSON format with 3-letter codes
            bias_dict = {}
            for item in bias_AA.split(","):
                aa, val = item.strip().split(":")
                aa1 = aa.strip().upper()
                aa3 = AA_1TO3.get(aa1, aa1)  # Convert to 3-letter if single-letter
                bias_dict[aa3] = float(val.strip())
            import json
            bias_json = json.dumps(bias_dict)
            cmd_parts.append(f"--bias '{bias_json}'")
            print(f"[MPNN] CLI bias: {bias_json}")

        # Add amino acid omission (atomworks CLI uses --omit with JSON list format)
        if omit_AA:
            # Convert single-letter codes to 3-letter codes as JSON list
            omit_list = []
            for aa in omit_AA.split(","):
                aa1 = aa.strip().upper()
                aa3 = AA_1TO3.get(aa1, aa1)  # Convert to 3-letter if single-letter
                omit_list.append(aa3)
            import json
            omit_json = json.dumps(omit_list)
            cmd_parts.append(f"--omit '{omit_json}'")
            print(f"[MPNN] CLI omit: {omit_json}")

        # Note: LigandMPNN-specific parameters (cutoff, atom context, side chain context)
        # are handled automatically by the atomworks CLI when model_type=ligand_mpnn
        # The ligand awareness is built into the model itself
        if model_type == "ligand_mpnn":
            print(f"[MPNN] Using ligand_mpnn model (ligand awareness automatic)")

        # Sidechain atomization and structure output
        if pack_side_chains:
            cmd_parts.append("--atomize_side_chains True")
            # write_structures is already set above
            print(f"[MPNN] Sidechain atomization enabled")

        cmd = " ".join(cmd_parts)
        print(f"[MPNN] Running CLI: {cmd}")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=1800)

        if result.returncode == 0:
            sequences = []
            packed_structures = []
            stats_data = None

            for root, _, files in os.walk(out_dir):
                for filename in files:
                    filepath = os.path.join(root, filename)
                    if filename.endswith((".fa", ".fasta")):
                        with open(filepath) as f:
                            sequences.append({"filename": filename, "content": f.read()})
                    elif filename.endswith(".pdb") and pack_side_chains:
                        with open(filepath) as f:
                            packed_structures.append({"filename": filename, "pdb_content": f.read()})
                    elif filename.endswith("_stats.json") and save_stats:
                        import json
                        with open(filepath) as f:
                            stats_data = json.load(f)

            if sequences:
                response = {
                    "status": "completed",
                    "result": {
                        "sequences": sequences,
                        "model_type": model_type,
                        "model_noise_level": model_noise_level,
                        "mode": "real (cli)",
                    }
                }
                if packed_structures:
                    response["result"]["packed_structures"] = packed_structures
                if stats_data:
                    response["result"]["stats"] = stats_data
                return response
            else:
                return {
                    "status": "failed",
                    "error": f"No output sequences found. stdout: {result.stdout}"
                }
        else:
            return {
                "status": "failed",
                "error": result.stderr or result.stdout or "Unknown error"
            }


def run_mpnn_mock(pdb_content: str, num_sequences: int) -> Dict[str, Any]:
    """Mock MPNN inference"""
    lines = pdb_content.split('\n')
    ca_count = sum(1 for l in lines if l.startswith('ATOM') and ' CA ' in l)
    length = ca_count if ca_count > 0 else 100

    return {
        "status": "completed",
        "result": {
            "sequences": [{
                "filename": "sequences.fasta",
                "content": generate_mock_fasta(num_sequences, length),
            }],
            "model_type": "mock",
            "mode": "mock",
        }
    }


# ============== RMSD Calculation ==============

def calculate_rmsd(pdb1: str, pdb2: str, backbone_only: bool = True) -> Dict[str, Any]:
    """Calculate RMSD between two structures"""
    try:
        from biotite.structure import rmsd, superimpose
        import numpy as np

        BACKBONE_ATOMS = ['N', 'CA', 'C', 'O']

        aa1 = pdb_to_atom_array(pdb1)
        aa2 = pdb_to_atom_array(pdb2)

        if backbone_only:
            mask1 = np.isin(aa1.atom_name, BACKBONE_ATOMS)
            mask2 = np.isin(aa2.atom_name, BACKBONE_ATOMS)
            aa1 = aa1[mask1]
            aa2 = aa2[mask2]

        if len(aa1) != len(aa2):
            min_len = min(len(aa1), len(aa2))
            aa1 = aa1[:min_len]
            aa2 = aa2[:min_len]

        aa2_fitted, _ = superimpose(aa1, aa2)
        rmsd_value = float(rmsd(aa1, aa2_fitted))

        if rmsd_value < 1.0:
            interpretation = "Excellent"
            description = "Very high designability - structure is highly likely to fold as designed"
        elif rmsd_value < 2.0:
            interpretation = "Good"
            description = "Good designability - structure will likely fold correctly"
        elif rmsd_value < 3.0:
            interpretation = "Moderate"
            description = "Moderate designability - some structural deviation expected"
        else:
            interpretation = "Poor"
            description = "Low designability - significant structural deviation"

        return {
            "rmsd": round(rmsd_value, 3),
            "interpretation": interpretation,
            "description": description,
            "backbone_only": backbone_only,
            "num_atoms_compared": len(aa1),
        }

    except Exception as e:
        return {
            "error": str(e),
            "rmsd": None,
        }


# ============== Structure Analysis ==============

def analyze_structure(pdb_content: str, target_ligands: Optional[List[str]] = None) -> Dict[str, Any]:
    """
    Analyze a PDB structure for binding site information.

    This is an AI-assisted tool that:
    1. Identifies metal ions and ligands (HETATM records)
    2. Finds coordinating residues within cutoff distance
    3. Calculates binding pocket properties
    4. Suggests RFdiffusion parameters for redesign

    Args:
        pdb_content: PDB file content as string
        target_ligands: Optional list of specific ligand codes to analyze

    Returns:
        Analysis results including binding sites, coordinating residues, and design suggestions
    """
    try:
        import numpy as np
        from biotite.structure import get_residues

        atom_array = pdb_to_atom_array(pdb_content)

        # Identify HETATM (ligands/metals)
        hetatm_mask = atom_array.hetero
        hetatms = atom_array[hetatm_mask]

        # Get unique ligand residues
        ligand_info = []
        seen_ligands = set()

        for i, atom in enumerate(hetatms):
            res_name = atom.res_name
            chain_id = atom.chain_id
            res_id = atom.res_id
            key = (res_name, chain_id, res_id)

            if key not in seen_ligands:
                seen_ligands.add(key)

                # Skip water
                if res_name in ['HOH', 'WAT', 'DOD', 'H2O']:
                    continue

                # Filter by target ligands if specified
                if target_ligands and res_name not in target_ligands:
                    continue

                # Get atom coordinates for this ligand
                lig_mask = (hetatms.res_name == res_name) & (hetatms.chain_id == chain_id) & (hetatms.res_id == res_id)
                lig_atoms = hetatms[lig_mask]
                lig_center = np.mean(lig_atoms.coord, axis=0)

                ligand_info.append({
                    "name": res_name,
                    "chain": chain_id,
                    "res_id": int(res_id),
                    "num_atoms": len(lig_atoms),
                    "center": lig_center.tolist(),
                    "atom_names": list(lig_atoms.atom_name),
                })

        # Find coordinating residues for each ligand
        binding_sites = []
        protein_atoms = atom_array[~atom_array.hetero]

        for lig in ligand_info:
            lig_center = np.array(lig["center"])

            # Calculate distances from all protein atoms to ligand center
            distances = np.linalg.norm(protein_atoms.coord - lig_center, axis=1)

            # Find atoms within 4.5 Angstroms (typical coordination distance)
            coord_mask = distances < 4.5
            coord_atoms = protein_atoms[coord_mask]
            coord_distances = distances[coord_mask]

            # Group by residue
            coord_residues = {}
            for i, atom in enumerate(coord_atoms):
                res_key = f"{atom.chain_id}{atom.res_id}"
                if res_key not in coord_residues:
                    coord_residues[res_key] = {
                        "chain": atom.chain_id,
                        "res_id": int(atom.res_id),
                        "res_name": atom.res_name,
                        "atoms": [],
                        "min_distance": float('inf'),
                    }
                coord_residues[res_key]["atoms"].append({
                    "name": atom.atom_name,
                    "distance": round(float(coord_distances[i]), 2),
                })
                coord_residues[res_key]["min_distance"] = min(
                    coord_residues[res_key]["min_distance"],
                    float(coord_distances[i])
                )

            # Sort by distance
            sorted_residues = sorted(
                coord_residues.values(),
                key=lambda x: x["min_distance"]
            )

            # Identify likely coordinating atoms (carboxyl oxygens for metals)
            likely_coordinators = []
            for res in sorted_residues:
                res_name = res["res_name"]
                # Acidic residues often coordinate metals
                if res_name in ['ASP', 'GLU']:
                    for atom in res["atoms"]:
                        if atom["name"] in ['OD1', 'OD2', 'OE1', 'OE2'] and atom["distance"] < 3.0:
                            likely_coordinators.append(f"{res['chain']}{res['res_id']}")
                            break
                # Histidine can also coordinate
                elif res_name == 'HIS':
                    for atom in res["atoms"]:
                        if atom["name"] in ['ND1', 'NE2'] and atom["distance"] < 3.0:
                            likely_coordinators.append(f"{res['chain']}{res['res_id']}")
                            break
                # Backbone carbonyl can coordinate
                elif any(a["name"] == 'O' and a["distance"] < 3.0 for a in res["atoms"]):
                    likely_coordinators.append(f"{res['chain']}{res['res_id']}")

            binding_sites.append({
                "ligand": lig,
                "coordinating_residues": sorted_residues[:10],  # Top 10 closest
                "likely_coordinators": list(set(likely_coordinators)),
                "coordination_number": len(likely_coordinators),
            })

        # Get protein info
        protein_residues = get_residues(protein_atoms)
        num_residues = len(protein_residues[0])

        # Generate RFdiffusion suggestions
        suggestions = []
        for site in binding_sites:
            lig_name = site["ligand"]["name"]
            coordinators = site["likely_coordinators"]

            # Metal ions that can be replaced with lanthanides
            METAL_IONS = ['CA', 'MG', 'ZN', 'FE', 'MN', 'CO', 'NI', 'CU']
            LANTHANIDES = ['LA', 'CE', 'PR', 'ND', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU']

            if lig_name in METAL_IONS:
                suggestions.append({
                    "type": "metal_replacement",
                    "description": f"Replace {lig_name} with a lanthanide (e.g., TB, GD, EU)",
                    "rfd3_params": {
                        "ligand": "TB",  # Example: terbium
                        "partial_t": 15,  # Moderate noise for pocket refinement
                        "unindex": ",".join(coordinators) if coordinators else None,
                        "is_non_loopy": True,
                        "num_timesteps": 50,
                    },
                    "notes": [
                        f"Current coordination number: {site['coordination_number']}",
                        f"Lanthanides prefer 8-9 coordination (vs {lig_name}'s typical 6-7)",
                        "Consider adding more coordinating residues for optimal lanthanide binding",
                        f"Coordinating residues to keep: {', '.join(coordinators)}",
                    ]
                })
            elif lig_name in LANTHANIDES:
                suggestions.append({
                    "type": "lanthanide_optimization",
                    "description": f"Optimize binding pocket for {lig_name}",
                    "rfd3_params": {
                        "ligand": lig_name,
                        "partial_t": 10,  # Lower noise for fine-tuning
                        "unindex": ",".join(coordinators) if coordinators else None,
                        "is_non_loopy": True,
                    },
                    "notes": [
                        f"Current coordination number: {site['coordination_number']}",
                        "Lanthanides prefer higher coordination numbers (8-9)",
                    ]
                })

        return {
            "status": "completed",
            "result": {
                "num_residues": num_residues,
                "ligands": ligand_info,
                "binding_sites": binding_sites,
                "suggestions": suggestions,
            }
        }

    except Exception as e:
        import traceback
        return {
            "status": "failed",
            "error": str(e),
            "traceback": traceback.format_exc(),
        }


# ============== Theozyme Creation ==============

def create_theozyme(
    ligand_pdb: str,
    binding_residues: List[Dict[str, Any]],
    chain_id: str = "A"
) -> str:
    """
    Create a theozyme PDB with ligand + key binding residues.

    A theozyme is an idealized 3D arrangement of catalytic/binding functional
    groups positioned around a target molecule. RFD3 can then scaffold around
    this theozyme to create a complete protein structure.

    This is the RECOMMENDED approach for ligand-conditioned protein design:
    1. Position ligand at desired location (usually origin)
    2. Define key binding residues (aromatic for pi-stacking, polar for H-bonds)
    3. Run RFD3 with theozyme as input for motif scaffolding

    Args:
        ligand_pdb: PDB content for the ligand (HETATM records)
        binding_residues: List of dicts defining key binding residues:
            [
                {"type": "PHE", "position": [3.5, 0.0, 2.0], "interaction": "pi_stack"},
                {"type": "TYR", "position": [-3.5, 0.0, -2.0], "interaction": "pi_stack"},
                {"type": "SER", "position": [0.0, 4.0, 0.0], "interaction": "hbond"},
            ]
        chain_id: Chain ID for the binding residues (default: "A")

    Returns:
        Combined PDB content with ligand + binding residue atoms

    Example usage:
        theozyme = create_theozyme(
            ligand_pdb=azobenzene_pdb,
            binding_residues=[
                {"type": "PHE", "position": [3.5, 0.0, 2.0], "interaction": "pi_stack"},
                {"type": "PHE", "position": [-3.5, 0.0, -2.0], "interaction": "pi_stack"},
            ]
        )
        # Then run RFD3 with theozyme for motif scaffolding
    """
    # Standard atom templates for common binding residues
    # These are idealized positions relative to the CA atom
    RESIDUE_TEMPLATES = {
        "PHE": {  # Phenylalanine - for pi-stacking
            "atoms": [
                ("N", [-1.458, 0.0, 0.0]),
                ("CA", [0.0, 0.0, 0.0]),
                ("C", [0.524, 1.420, 0.0]),
                ("O", [-0.216, 2.409, 0.0]),
                ("CB", [0.529, -0.777, -1.209]),
                ("CG", [0.011, -0.368, -2.575]),
                ("CD1", [-1.331, -0.147, -2.848]),
                ("CD2", [0.896, -0.195, -3.635]),
                ("CE1", [-1.785, 0.237, -4.104]),
                ("CE2", [0.447, 0.189, -4.893]),
                ("CZ", [-0.888, 0.405, -5.140]),
            ]
        },
        "TYR": {  # Tyrosine - for pi-stacking + H-bonding
            "atoms": [
                ("N", [-1.458, 0.0, 0.0]),
                ("CA", [0.0, 0.0, 0.0]),
                ("C", [0.524, 1.420, 0.0]),
                ("O", [-0.216, 2.409, 0.0]),
                ("CB", [0.529, -0.777, -1.209]),
                ("CG", [0.011, -0.368, -2.575]),
                ("CD1", [-1.331, -0.147, -2.848]),
                ("CD2", [0.896, -0.195, -3.635]),
                ("CE1", [-1.785, 0.237, -4.104]),
                ("CE2", [0.447, 0.189, -4.893]),
                ("CZ", [-0.888, 0.405, -5.140]),
                ("OH", [-1.330, 0.784, -6.363]),
            ]
        },
        "TRP": {  # Tryptophan - for pi-stacking
            "atoms": [
                ("N", [-1.458, 0.0, 0.0]),
                ("CA", [0.0, 0.0, 0.0]),
                ("C", [0.524, 1.420, 0.0]),
                ("O", [-0.216, 2.409, 0.0]),
                ("CB", [0.529, -0.777, -1.209]),
                ("CG", [0.011, -0.368, -2.575]),
                ("CD1", [-1.260, -0.024, -2.992]),
                ("CD2", [0.776, -0.195, -3.775]),
                ("NE1", [-1.322, 0.298, -4.330]),
                ("CE2", [-0.047, 0.180, -4.881]),
                ("CE3", [2.157, -0.409, -3.964]),
                ("CZ2", [0.383, 0.427, -6.192]),
                ("CZ3", [2.576, -0.165, -5.271]),
                ("CH2", [1.672, 0.279, -6.361]),
            ]
        },
        "SER": {  # Serine - for H-bonding
            "atoms": [
                ("N", [-1.458, 0.0, 0.0]),
                ("CA", [0.0, 0.0, 0.0]),
                ("C", [0.524, 1.420, 0.0]),
                ("O", [-0.216, 2.409, 0.0]),
                ("CB", [0.529, -0.777, -1.209]),
                ("OG", [0.011, -0.368, -2.438]),
            ]
        },
        "THR": {  # Threonine - for H-bonding
            "atoms": [
                ("N", [-1.458, 0.0, 0.0]),
                ("CA", [0.0, 0.0, 0.0]),
                ("C", [0.524, 1.420, 0.0]),
                ("O", [-0.216, 2.409, 0.0]),
                ("CB", [0.529, -0.777, -1.209]),
                ("OG1", [0.011, -0.368, -2.438]),
                ("CG2", [2.054, -0.777, -1.209]),
            ]
        },
        "ASN": {  # Asparagine - for H-bonding
            "atoms": [
                ("N", [-1.458, 0.0, 0.0]),
                ("CA", [0.0, 0.0, 0.0]),
                ("C", [0.524, 1.420, 0.0]),
                ("O", [-0.216, 2.409, 0.0]),
                ("CB", [0.529, -0.777, -1.209]),
                ("CG", [0.011, -0.368, -2.575]),
                ("OD1", [-1.175, -0.139, -2.819]),
                ("ND2", [0.918, -0.170, -3.520]),
            ]
        },
        "GLN": {  # Glutamine - for H-bonding
            "atoms": [
                ("N", [-1.458, 0.0, 0.0]),
                ("CA", [0.0, 0.0, 0.0]),
                ("C", [0.524, 1.420, 0.0]),
                ("O", [-0.216, 2.409, 0.0]),
                ("CB", [0.529, -0.777, -1.209]),
                ("CG", [0.011, -0.368, -2.575]),
                ("CD", [0.529, -1.145, -3.784]),
                ("OE1", [-0.102, -1.208, -4.842]),
                ("NE2", [1.719, -1.748, -3.634]),
            ]
        },
    }

    import numpy as np

    output_lines = []
    atom_serial = 1
    res_num = 1

    # First add binding residues as ATOM records
    for residue in binding_residues:
        res_type = residue.get("type", "PHE").upper()
        position = residue.get("position", [0.0, 0.0, 0.0])
        res_chain = residue.get("chain", chain_id)

        if res_type not in RESIDUE_TEMPLATES:
            print(f"[Theozyme] Warning: Unknown residue type {res_type}, skipping")
            continue

        template = RESIDUE_TEMPLATES[res_type]
        pos = np.array(position)

        for atom_name, rel_coord in template["atoms"]:
            coord = pos + np.array(rel_coord)
            # Format PDB ATOM record
            line = f"ATOM  {atom_serial:5d}  {atom_name:<3s} {res_type:3s} {res_chain}{res_num:4d}    {coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}  1.00  0.00           {atom_name[0]:>1s}"
            output_lines.append(line)
            atom_serial += 1

        res_num += 1

    # Add TER record between protein and ligand
    output_lines.append("TER")

    # Add ligand (HETATM records)
    ligand_lines = [l for l in ligand_pdb.split('\n') if l.startswith('HETATM')]
    output_lines.extend(ligand_lines)

    # Add END
    output_lines.append("END")

    return '\n'.join(output_lines)


def generate_ligand_pdb(
    smiles: str,
    name: str = "LIG",
    center: Tuple[float, float, float] = (0.0, 0.0, 0.0)
) -> Optional[str]:
    """
    Generate a 3D PDB structure from a SMILES string.

    Uses RDKit to generate 3D coordinates and outputs as PDB HETATM records.
    The ligand is centered at the specified position.

    Args:
        smiles: SMILES string for the molecule
        name: 3-letter residue name for the ligand (default: "LIG")
        center: (x, y, z) coordinates for ligand center (default: origin)

    Returns:
        PDB content with HETATM records, or None if generation fails
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        import numpy as np

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)

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

        # Generate PDB content
        pdb_lines = []
        atom_serial = 1

        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            coord = np.array([pos.x, pos.y, pos.z]) + translation
            symbol = atom.GetSymbol()
            atom_name = f"{symbol}{i+1}"[:4]

            line = f"HETATM{atom_serial:5d}  {atom_name:<3s} {name:3s} L   1    {coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}  1.00  0.00          {symbol:>2s}"
            pdb_lines.append(line)
            atom_serial += 1

        pdb_lines.append("END")
        return '\n'.join(pdb_lines)

    except ImportError:
        print("[Theozyme] RDKit not available for SMILES to PDB conversion")
        return None
    except Exception as e:
        print(f"[Theozyme] Error generating PDB from SMILES: {e}")
        return None


# ============== ESM3 Inference ==============

# Global lazy-loaded ESM3 model
_esm3_model = None
_esm3_tokenizers = None


def get_esm3_model():
    """Lazy load ESM3 model on first use.

    Note: User must set HF_TOKEN environment variable and accept
    ESM3 license on HuggingFace before first use.
    """
    global _esm3_model, _esm3_tokenizers

    if _esm3_model is None:
        print("[ESM3] Loading model esm3-sm-open-v1...")
        from esm.models.esm3 import ESM3
        from esm.tokenization import get_model_tokenizers

        # Model: EvolutionaryScale/esm3-sm-open-v1 (~1.4B params, ~10GB VRAM)
        _esm3_model = ESM3.from_pretrained("esm3-sm-open-v1").to("cuda")
        _esm3_tokenizers = get_model_tokenizers("esm3-sm-open-v1")
        print("[ESM3] Model loaded successfully")

    return _esm3_model, _esm3_tokenizers


def unload_esm3():
    """Free ESM3 model memory when not needed.

    Call this before loading memory-intensive models like RF3.
    """
    global _esm3_model, _esm3_tokenizers
    import torch
    import gc

    if _esm3_model is not None:
        print("[ESM3] Unloading model to free memory...")
        del _esm3_model
        del _esm3_tokenizers
        _esm3_model = None
        _esm3_tokenizers = None
        torch.cuda.empty_cache()
        gc.collect()
        print("[ESM3] Model unloaded")


def esm3_score_sequences(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """Score sequences using ESM3 perplexity.

    Lower perplexity = higher quality/more natural sequence.
    Score = 1/perplexity (higher = better).

    Args:
        job_input: Dict with 'sequences' key containing list of amino acid sequences

    Returns:
        Dict with scores list containing {sequence, perplexity, score} for each input
    """
    import torch
    import traceback

    sequences = job_input.get("sequences", [])
    if not sequences:
        return {"status": "failed", "error": "No sequences provided"}

    try:
        from esm.sdk.api import ESMProtein, SamplingConfig, SamplingTrackConfig

        model, tokenizers = get_esm3_model()
        results = []

        print(f"[ESM3] Scoring {len(sequences)} sequences")

        for seq in sequences:
            try:
                # Create protein object
                protein = ESMProtein(sequence=seq)

                # Forward pass to get logits
                with torch.no_grad():
                    output = model.forward_and_sample(
                        protein,
                        SamplingTrackConfig(sequence=SamplingConfig(temperature=0)),
                    )

                    # Compute perplexity from logits
                    perplexity = 1.0  # Default
                    if hasattr(output, 'sequence_logits') and output.sequence_logits is not None:
                        logits = output.sequence_logits
                        tokens = tokenizers.sequence.encode(seq)
                        # Get log probabilities
                        log_probs = torch.nn.functional.log_softmax(logits, dim=-1)
                        # Get the log prob of the actual tokens
                        token_log_probs = log_probs.gather(-1, tokens.unsqueeze(-1)).squeeze(-1)
                        # Average negative log likelihood = perplexity
                        perplexity = torch.exp(-token_log_probs.mean()).item()

                results.append({
                    "sequence": seq[:50] + "..." if len(seq) > 50 else seq,
                    "length": len(seq),
                    "perplexity": round(perplexity, 3),
                    "score": round(1.0 / max(perplexity, 0.001), 3),
                })
            except Exception as e:
                results.append({
                    "sequence": seq[:50] + "..." if len(seq) > 50 else seq,
                    "length": len(seq),
                    "error": str(e),
                    "perplexity": None,
                    "score": None,
                })

        return {
            "status": "completed",
            "result": {"scores": results}
        }

    except ImportError as e:
        return {
            "status": "failed",
            "error": f"ESM3 not available: {e}. Make sure esm package is installed and HF_TOKEN is set."
        }
    except Exception as e:
        return {
            "status": "failed",
            "error": str(e),
            "traceback": traceback.format_exc()
        }


def esm3_generate_sequence(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """Generate protein sequences using ESM3 with optional function conditioning.

    Args:
        job_input: Dict with:
            - prompt: Optional partial sequence to complete
            - functions: Optional list of function keywords (e.g., ["zinc-binding", "hydrolase"])
            - num_sequences: Number of sequences to generate (default: 4)
            - temperature: Sampling temperature (default: 0.7)
            - max_length: Maximum sequence length (default: 200)

    Returns:
        Dict with 'sequences' list of generated amino acid sequences
    """
    import torch
    import traceback

    prompt = job_input.get("prompt", "")
    function_keywords = job_input.get("functions", [])
    num_sequences = job_input.get("num_sequences", 4)
    temperature = job_input.get("temperature", 0.7)
    max_length = job_input.get("max_length", 200)

    try:
        from esm.sdk.api import ESMProtein, GenerationConfig

        model, tokenizers = get_esm3_model()
        generated = []

        print(f"[ESM3] Generating {num_sequences} sequences, temp={temperature}, max_len={max_length}")
        if function_keywords:
            print(f"[ESM3] Function conditioning: {function_keywords}")

        for i in range(num_sequences):
            try:
                # Build protein with optional function annotations
                protein_kwargs = {}
                if prompt:
                    protein_kwargs["sequence"] = prompt

                # Add function annotations if provided
                if function_keywords and hasattr(tokenizers, 'function'):
                    try:
                        function_tokens = tokenizers.function.encode_keywords(function_keywords)
                        protein_kwargs["function_annotations"] = function_tokens
                    except Exception as e:
                        print(f"[ESM3] Warning: Could not encode function keywords: {e}")

                protein = ESMProtein(**protein_kwargs) if protein_kwargs else ESMProtein()

                # Configure generation
                gen_config = GenerationConfig(
                    track="sequence",
                    num_steps=max_length,
                    temperature=temperature,
                )

                # Generate
                with torch.no_grad():
                    output = model.generate(protein, gen_config)

                if hasattr(output, 'sequence') and output.sequence:
                    generated.append(output.sequence)
                    print(f"[ESM3] Generated sequence {i+1}: {len(output.sequence)} residues")
            except Exception as e:
                print(f"[ESM3] Error generating sequence {i+1}: {e}")
                continue

        if not generated:
            return {
                "status": "failed",
                "error": "Failed to generate any sequences"
            }

        return {
            "status": "completed",
            "result": {
                "sequences": generated,
                "num_generated": len(generated),
                "temperature": temperature,
                "max_length": max_length,
                "function_keywords": function_keywords if function_keywords else None,
            }
        }

    except ImportError as e:
        return {
            "status": "failed",
            "error": f"ESM3 not available: {e}. Make sure esm package is installed and HF_TOKEN is set."
        }
    except Exception as e:
        return {
            "status": "failed",
            "error": str(e),
            "traceback": traceback.format_exc()
        }


def esm3_get_embeddings(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """Get ESM3 embeddings for protein sequences.

    Embeddings can be used for:
    - Sequence similarity comparison
    - Clustering proteins by function
    - Feature input for downstream ML models

    Args:
        job_input: Dict with 'sequences' key containing list of amino acid sequences

    Returns:
        Dict with embeddings list containing {sequence, per_residue, global} for each input
    """
    import torch
    import traceback

    sequences = job_input.get("sequences", [])
    if not sequences:
        return {"status": "failed", "error": "No sequences provided"}

    # Limit embedding size for response
    max_residues_per_seq = 100

    try:
        from esm.sdk.api import ESMProtein

        model, tokenizers = get_esm3_model()
        results = []

        print(f"[ESM3] Computing embeddings for {len(sequences)} sequences")

        for seq in sequences:
            try:
                protein = ESMProtein(sequence=seq)

                with torch.no_grad():
                    output = model.encode(protein)

                if hasattr(output, 'embeddings') and output.embeddings is not None:
                    embeddings = output.embeddings

                    # Per-residue embeddings (truncate for response size)
                    per_residue = embeddings[:max_residues_per_seq].cpu().numpy().tolist()

                    # Global embedding (mean pooled)
                    global_emb = embeddings.mean(dim=0).cpu().numpy().tolist()

                    results.append({
                        "sequence": seq[:50] + "..." if len(seq) > 50 else seq,
                        "length": len(seq),
                        "embedding_dim": len(global_emb),
                        "per_residue": per_residue,
                        "global": global_emb,
                    })
                    print(f"[ESM3] Computed embedding for sequence of length {len(seq)}")
                else:
                    results.append({
                        "sequence": seq[:50] + "..." if len(seq) > 50 else seq,
                        "length": len(seq),
                        "error": "No embeddings returned",
                    })

            except Exception as e:
                results.append({
                    "sequence": seq[:50] + "..." if len(seq) > 50 else seq,
                    "length": len(seq),
                    "error": str(e),
                })

        return {
            "status": "completed",
            "result": {
                "embeddings": results,
                "num_processed": len(results),
            }
        }

    except ImportError as e:
        return {
            "status": "failed",
            "error": f"ESM3 not available: {e}. Make sure esm package is installed and HF_TOKEN is set."
        }
    except Exception as e:
        return {
            "status": "failed",
            "error": str(e),
            "traceback": traceback.format_exc()
        }
