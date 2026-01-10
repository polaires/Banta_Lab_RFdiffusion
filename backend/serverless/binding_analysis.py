"""
Binding Analysis Module

Provides tools for evaluating protein-protein and protein-ligand binding quality:
1. Interface Analyzer - Python-based interface metrics (no Rosetta license needed)
2. GNINA Integration - CNN-based docking and affinity scoring

These tools complement RF3's ipTM by providing:
- Detailed interaction profiling
- Binding energy estimates
- Shape complementarity metrics
"""

import os
import json
import tempfile
import subprocess
from typing import Dict, Any, Optional, List, Tuple
import numpy as np


def to_python_types(obj: Any) -> Any:
    """Recursively convert numpy types to native Python types for JSON serialization."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, (np.float32, np.float64)):
        return float(obj)
    elif isinstance(obj, (np.int32, np.int64)):
        return int(obj)
    elif isinstance(obj, np.str_):
        return str(obj)
    elif isinstance(obj, dict):
        return {k: to_python_types(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [to_python_types(item) for item in obj]
    elif isinstance(obj, tuple):
        return tuple(to_python_types(item) for item in obj)
    elif isinstance(obj, set):
        return {to_python_types(item) for item in obj}
    return obj


# ============== Interface Analyzer (Python-based) ==============

def analyze_interface(
    pdb_content: str,
    chain_a: str = "A",
    chain_b: str = "B",
    contact_distance: float = 8.0,
    hbond_distance: float = 3.5,
) -> Dict[str, Any]:
    """
    Analyze protein-protein or protein-ligand interface.

    Calculates metrics similar to Rosetta InterfaceAnalyzer:
    - dSASA_int: Buried surface area at interface
    - nres_int: Number of interface residues
    - hbonds_int: Cross-interface hydrogen bonds
    - contacts: Number of atomic contacts
    - packstat: Packing quality estimate (0-1)

    Args:
        pdb_content: PDB file content as string
        chain_a: First chain identifier
        chain_b: Second chain identifier (or 'HETATM' for ligand)
        contact_distance: Distance cutoff for contacts (default 8 Angstrom)
        hbond_distance: Distance cutoff for H-bonds (default 3.5 Angstrom)

    Returns:
        Dict with interface metrics and per-residue details
    """
    try:
        from biotite.structure.io.pdb import PDBFile
        from biotite.structure import sasa, hbond
        import io

        # Parse PDB
        pdb_file = PDBFile.read(io.StringIO(pdb_content))
        structure = pdb_file.get_structure(model=1)

        # Separate chains
        if chain_b == "HETATM":
            # Protein-ligand interface
            chain_a_atoms = structure[structure.chain_id == chain_a]
            chain_b_atoms = structure[structure.hetero]
        else:
            # Protein-protein interface
            chain_a_atoms = structure[structure.chain_id == chain_a]
            chain_b_atoms = structure[structure.chain_id == chain_b]

        if len(chain_a_atoms) == 0 or len(chain_b_atoms) == 0:
            return {
                "status": "error",
                "error": f"Could not find atoms for chains {chain_a} and {chain_b}"
            }

        # Calculate contacts
        contacts = []
        interface_residues_a = set()
        interface_residues_b = set()

        for i, atom_a in enumerate(chain_a_atoms):
            for j, atom_b in enumerate(chain_b_atoms):
                dist = np.linalg.norm(atom_a.coord - atom_b.coord)
                if dist < contact_distance:
                    contacts.append({
                        "atom_a": atom_a.atom_name,
                        "res_a": f"{atom_a.res_name}{atom_a.res_id}",
                        "atom_b": atom_b.atom_name,
                        "res_b": f"{atom_b.res_name}{atom_b.res_id}",
                        "distance": round(dist, 2)
                    })
                    interface_residues_a.add((atom_a.res_name, atom_a.res_id))
                    interface_residues_b.add((atom_b.res_name, atom_b.res_id))

        # Identify hydrogen bonds (donor-acceptor pairs)
        hbonds_int = []
        donor_atoms = ['N', 'NE', 'NH1', 'NH2', 'ND1', 'ND2', 'NE1', 'NE2', 'NZ', 'OG', 'OG1', 'OH']
        acceptor_atoms = ['O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH', 'ND1', 'NE2']

        for contact in contacts:
            if contact["distance"] < hbond_distance:
                atom_a = contact["atom_a"]
                atom_b = contact["atom_b"]
                # Check if this could be an H-bond
                if (atom_a in donor_atoms and atom_b in acceptor_atoms) or \
                   (atom_a in acceptor_atoms and atom_b in donor_atoms):
                    hbonds_int.append(contact)

        # Calculate SASA for buried surface area
        try:
            # SASA of complex
            complex_sasa = sasa(structure)
            total_complex_sasa = np.sum(complex_sasa)

            # SASA of chain A alone
            chain_a_sasa = sasa(chain_a_atoms)
            total_a_sasa = np.sum(chain_a_sasa)

            # SASA of chain B alone
            chain_b_sasa = sasa(chain_b_atoms)
            total_b_sasa = np.sum(chain_b_sasa)

            # Buried surface area = (SASA_A + SASA_B - SASA_complex) / 2
            dSASA_int = (total_a_sasa + total_b_sasa - total_complex_sasa) / 2
        except Exception:
            dSASA_int = None

        # Estimate packing quality (simplified packstat)
        # Based on contact density at interface
        if len(interface_residues_a) > 0 and len(interface_residues_b) > 0:
            contacts_per_residue = len(contacts) / (len(interface_residues_a) + len(interface_residues_b))
            # Normalize to 0-1 range (typical good interfaces have 5-15 contacts per residue)
            packstat = min(1.0, contacts_per_residue / 10.0)
        else:
            packstat = 0.0

        # Classify residues by type at interface
        hydrophobic_residues = ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO']
        polar_residues = ['SER', 'THR', 'ASN', 'GLN', 'TYR', 'CYS']
        charged_residues = ['LYS', 'ARG', 'HIS', 'ASP', 'GLU']

        interface_composition = {
            "hydrophobic": 0,
            "polar": 0,
            "charged": 0
        }
        all_interface_residues = list(interface_residues_a) + list(interface_residues_b)
        for res_name, _ in all_interface_residues:
            if res_name in hydrophobic_residues:
                interface_composition["hydrophobic"] += 1
            elif res_name in polar_residues:
                interface_composition["polar"] += 1
            elif res_name in charged_residues:
                interface_composition["charged"] += 1

        # Estimate binding energy (simplified)
        # Based on: contacts, H-bonds, and buried surface area
        # Rough approximation: -0.5 kcal/mol per H-bond, -0.02 kcal/mol per contact
        estimated_dG = -0.5 * len(hbonds_int) - 0.02 * len(contacts)
        if dSASA_int:
            # Add contribution from hydrophobic burial (~-25 cal/mol per A^2)
            estimated_dG -= 0.025 * dSASA_int * (interface_composition["hydrophobic"] / max(1, len(all_interface_residues)))

        result = {
            "status": "completed",
            "metrics": {
                "nres_int": len(interface_residues_a) + len(interface_residues_b),
                "nres_chain_a": len(interface_residues_a),
                "nres_chain_b": len(interface_residues_b),
                "contacts": len(contacts),
                "hbonds_int": len(hbonds_int),
                "dSASA_int": round(float(dSASA_int), 1) if dSASA_int is not None else None,
                "packstat": round(float(packstat), 3),
                "estimated_dG": round(float(estimated_dG), 2),
                "interface_composition": interface_composition
            },
            "hbond_details": hbonds_int[:20],  # Top 20 H-bonds
            "contact_details": sorted(contacts, key=lambda x: x["distance"])[:30],  # Top 30 contacts
            "interface_residues": {
                "chain_a": [f"{r[0]}{r[1]}" for r in sorted(interface_residues_a, key=lambda x: x[1])],
                "chain_b": [f"{r[0]}{r[1]}" for r in sorted(interface_residues_b, key=lambda x: x[1])]
            }
        }
        # Convert numpy types to native Python types for JSON serialization
        return to_python_types(result)

    except ImportError as e:
        return {
            "status": "error",
            "error": f"Missing dependency: {e}"
        }
    except Exception as e:
        return {
            "status": "error",
            "error": str(e)
        }


# ============== GNINA Integration ==============

def check_gnina_available() -> bool:
    """Check if GNINA is available in the environment."""
    try:
        result = subprocess.run(
            ["gnina", "--version"],
            capture_output=True,
            text=True,
            timeout=10
        )
        return result.returncode == 0
    except Exception:
        return False


def run_gnina_scoring(
    receptor_pdb: str,
    ligand_sdf: str,
    autobox_ligand: Optional[str] = None,
    center: Optional[Tuple[float, float, float]] = None,  # Explicit box center
    box_size: float = 25.0,  # Larger search box for better coverage
    exhaustiveness: int = 32,  # Higher for more thorough search
    num_modes: int = 20,  # More poses to explore
    cnn_scoring: str = "all",  # Full CNN scoring for best results
    minimize: bool = True,
    seed: int = 42,  # Reproducibility
    whole_protein: bool = False,  # Search entire protein surface
) -> Dict[str, Any]:
    """
    Run GNINA for protein-ligand docking/scoring.

    GNINA provides CNN-based scoring that outperforms Vina by 20-25%.

    Args:
        receptor_pdb: Receptor PDB content as string
        ligand_sdf: Ligand SDF/MOL content as string (or SMILES)
        autobox_ligand: Reference ligand for defining search box (optional)
        center: Explicit (x, y, z) box center coordinates. If None, uses autobox_ligand.
        box_size: Search box size in Angstroms (default 25.0 for thorough coverage)
        exhaustiveness: Search thoroughness (default 16 for quality results)
        num_modes: Maximum binding modes to generate (default 20)
        cnn_scoring: CNN scoring mode: 'none', 'rescore', 'refinement', 'all' (default 'all')
        minimize: Whether to minimize poses (default True)
        seed: Random seed for reproducibility (default 42)

    Returns:
        Dict with docking results including CNN scores and affinities
    """
    if not check_gnina_available():
        return {
            "status": "error",
            "error": "GNINA not available. Install with: apt-get install gnina or use Docker image."
        }

    with tempfile.TemporaryDirectory() as tmpdir:
        # Write receptor
        receptor_path = os.path.join(tmpdir, "receptor.pdb")
        with open(receptor_path, "w") as f:
            f.write(receptor_pdb)

        # Write ligand
        ligand_path = os.path.join(tmpdir, "ligand.sdf")
        with open(ligand_path, "w") as f:
            f.write(ligand_sdf)

        # Output path
        output_path = os.path.join(tmpdir, "docked.sdf")

        # Build command
        cmd = [
            "gnina",
            "-r", receptor_path,
            "-l", ligand_path,
            "-o", output_path,
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes", str(num_modes),
            "--cnn_scoring", cnn_scoring,
            "--seed", str(seed),
        ]

        # Define search box
        if whole_protein:
            # Calculate protein bounding box to search entire surface
            protein_center, protein_size = calculate_protein_bounding_box(receptor_pdb)
            if protein_center:
                cmd.extend([
                    "--center_x", str(protein_center[0]),
                    "--center_y", str(protein_center[1]),
                    "--center_z", str(protein_center[2]),
                    "--size_x", str(protein_size + 10),  # Add padding
                    "--size_y", str(protein_size + 10),
                    "--size_z", str(protein_size + 10),
                ])
                print(f"[GNINA] Searching whole protein: center={protein_center}, size={protein_size}")
            else:
                cmd.extend(["--autobox_ligand", ligand_path, "--autobox_add", "20.0"])
        elif center is not None:
            # Use explicit box center and size
            cmd.extend([
                "--center_x", str(center[0]),
                "--center_y", str(center[1]),
                "--center_z", str(center[2]),
                "--size_x", str(box_size),
                "--size_y", str(box_size),
                "--size_z", str(box_size),
            ])
        elif autobox_ligand:
            autobox_path = os.path.join(tmpdir, "autobox.sdf")
            with open(autobox_path, "w") as f:
                f.write(autobox_ligand)
            cmd.extend(["--autobox_ligand", autobox_path, "--autobox_add", "4.0"])
        else:
            # Use ligand itself for autobox with padding
            cmd.extend(["--autobox_ligand", ligand_path, "--autobox_add", "4.0"])

        if minimize:
            cmd.append("--minimize")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout
            )

            if result.returncode != 0:
                return {
                    "status": "error",
                    "error": result.stderr or result.stdout or "GNINA failed"
                }

            # Parse output SDF for scores
            if os.path.exists(output_path):
                with open(output_path) as f:
                    docked_content = f.read()

                # Parse SDF properties
                poses = parse_gnina_output(docked_content)

                return {
                    "status": "completed",
                    "result": {
                        "num_poses": len(poses),
                        "poses": poses,
                        "best_affinity": poses[0]["affinity"] if poses else None,
                        "best_cnn_score": poses[0]["cnn_score"] if poses else None,
                        "docked_sdf": docked_content
                    }
                }
            else:
                return {
                    "status": "error",
                    "error": "No output generated"
                }

        except subprocess.TimeoutExpired:
            return {
                "status": "error",
                "error": "GNINA timed out after 5 minutes"
            }
        except Exception as e:
            return {
                "status": "error",
                "error": str(e)
            }


def parse_gnina_output(sdf_content: str) -> List[Dict[str, Any]]:
    """Parse GNINA SDF output to extract scores."""
    poses = []
    current_pose = {}

    lines = sdf_content.split('\n')
    for i, line in enumerate(lines):
        # Handle both '>  <tag>' and '> <tag>' formats
        line_stripped = line.strip()
        if '<minimizedAffinity>' in line_stripped:
            if i + 1 < len(lines):
                try:
                    current_pose["affinity"] = float(lines[i + 1].strip())
                except ValueError:
                    pass
        elif '<CNNscore>' in line_stripped:
            if i + 1 < len(lines):
                try:
                    current_pose["cnn_score"] = float(lines[i + 1].strip())
                except ValueError:
                    pass
        elif '<CNNaffinity>' in line_stripped:
            if i + 1 < len(lines):
                try:
                    current_pose["cnn_affinity"] = float(lines[i + 1].strip())
                except ValueError:
                    pass
        elif line.strip() == '$$$$':
            if current_pose:
                poses.append(current_pose)
                current_pose = {}

    # Sort by CNN score (higher is better)
    poses.sort(key=lambda x: x.get("cnn_score", 0), reverse=True)
    return poses


def smiles_to_sdf(smiles: str, name: str = "ligand") -> Optional[str]:
    """Convert SMILES to SDF format using RDKit."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)

        # Convert to SDF
        mol.SetProp("_Name", name)
        sdf = Chem.MolToMolBlock(mol)
        return sdf + "\n$$$$\n"

    except ImportError:
        return None
    except Exception:
        return None


# ============== Combined Binding Evaluation ==============

def evaluate_binding(
    pdb_content: str,
    ligand_smiles: Optional[str] = None,
    ligand_sdf: Optional[str] = None,
    chain_a: str = "A",
    chain_b: str = "B",
    run_gnina: bool = True,
    docking_center: Optional[Tuple[float, float, float]] = None,
    docking_box_size: float = 25.0,
    whole_protein_search: bool = False,
) -> Dict[str, Any]:
    """
    Comprehensive binding evaluation combining multiple methods.

    Runs:
    1. Interface Analyzer - contacts, H-bonds, buried surface area
    2. GNINA Scoring - CNN-based affinity prediction (if available)

    Args:
        pdb_content: PDB content (protein complex or protein-ligand)
        ligand_smiles: Ligand SMILES (for GNINA, if no ligand in PDB)
        ligand_sdf: Ligand SDF content (alternative to SMILES)
        chain_a: First chain for interface analysis
        chain_b: Second chain (or 'HETATM' for ligand)
        run_gnina: Whether to run GNINA scoring
        docking_center: Explicit (x, y, z) center for docking box
        docking_box_size: Size of docking box in Angstroms

    Returns:
        Combined binding evaluation results
    """
    results = {
        "status": "completed",
        "interface_analysis": None,
        "gnina_scoring": None,
        "summary": {}
    }

    # Run interface analysis
    interface_result = analyze_interface(
        pdb_content=pdb_content,
        chain_a=chain_a,
        chain_b=chain_b
    )
    results["interface_analysis"] = interface_result

    # Run GNINA if requested and ligand is available
    if run_gnina and (ligand_smiles or ligand_sdf):
        # Prepare ligand SDF
        if ligand_sdf is None and ligand_smiles:
            ligand_sdf = smiles_to_sdf(ligand_smiles)

        if ligand_sdf:
            # Extract protein-only PDB for receptor
            receptor_pdb = extract_protein_only(pdb_content)

            # Calculate interface center if not provided and not doing whole protein search
            center = docking_center
            if center is None and not whole_protein_search:
                center = calculate_interface_center(pdb_content, chain_a, chain_b)

            gnina_result = run_gnina_scoring(
                receptor_pdb=receptor_pdb,
                ligand_sdf=ligand_sdf,
                center=center,
                box_size=docking_box_size,
                whole_protein=whole_protein_search,
            )
            results["gnina_scoring"] = gnina_result

    # Generate summary
    summary = {}

    if interface_result.get("status") == "completed":
        metrics = interface_result.get("metrics", {})
        summary["interface"] = {
            "contacts": metrics.get("contacts", 0),
            "hbonds": metrics.get("hbonds_int", 0),
            "buried_area": metrics.get("dSASA_int"),
            "packstat": metrics.get("packstat", 0),
            "estimated_dG": metrics.get("estimated_dG")
        }

        # Quality assessment
        contacts = metrics.get("contacts", 0)
        hbonds = metrics.get("hbonds_int", 0)
        packstat = metrics.get("packstat", 0)

        if contacts >= 50 and hbonds >= 4 and packstat >= 0.5:
            summary["interface_quality"] = "good"
        elif contacts >= 20 and hbonds >= 2:
            summary["interface_quality"] = "moderate"
        else:
            summary["interface_quality"] = "poor"

    if results.get("gnina_scoring", {}).get("status") == "completed":
        gnina_metrics = results["gnina_scoring"].get("result", {})
        summary["gnina"] = {
            "best_affinity": gnina_metrics.get("best_affinity"),
            "best_cnn_score": gnina_metrics.get("best_cnn_score"),
            "num_poses": gnina_metrics.get("num_poses", 0)
        }

        # Affinity assessment
        affinity = gnina_metrics.get("best_affinity")
        if affinity is not None:
            if affinity < -8:
                summary["affinity_quality"] = "strong"
            elif affinity < -6:
                summary["affinity_quality"] = "moderate"
            else:
                summary["affinity_quality"] = "weak"

    # Calculate composite binding score (0-100)
    summary["composite_score"] = calculate_composite_score(
        interface_result,
        results.get("gnina_scoring")
    )

    results["summary"] = summary
    # Convert numpy types to native Python types for JSON serialization
    return to_python_types(results)


def calculate_composite_score(
    interface_result: Dict[str, Any],
    gnina_result: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    Calculate a composite binding score from multiple metrics.

    Scoring rubric (0-100 scale):
    - Interface quality (40 points max):
      - Contacts density: 0-15 points
      - H-bonds: 0-10 points
      - Buried surface area: 0-10 points
      - Packing quality: 0-5 points
    - GNINA affinity (40 points max):
      - CNN affinity: 0-20 points
      - CNN score: 0-20 points
    - Composition bonus (20 points max):
      - Hydrophobic burial: 0-10 points
      - Aromatic contacts: 0-10 points

    Returns:
        Dict with composite score and breakdown
    """
    score_breakdown = {
        "interface_contacts": 0,
        "interface_hbonds": 0,
        "interface_burial": 0,
        "interface_packing": 0,
        "gnina_affinity": 0,
        "gnina_cnn_score": 0,
        "composition_bonus": 0,
    }
    max_scores = {
        "interface_contacts": 15,
        "interface_hbonds": 10,
        "interface_burial": 10,
        "interface_packing": 5,
        "gnina_affinity": 20,
        "gnina_cnn_score": 20,
        "composition_bonus": 20,
    }

    # Interface scoring
    if interface_result and interface_result.get("status") == "completed":
        metrics = interface_result.get("metrics", {})

        # Contacts (scaled: 100+ contacts = max score)
        contacts = metrics.get("contacts", 0)
        score_breakdown["interface_contacts"] = min(15, contacts / 100 * 15)

        # H-bonds (scaled: 10+ H-bonds = max score)
        hbonds = metrics.get("hbonds_int", 0)
        score_breakdown["interface_hbonds"] = min(10, hbonds / 10 * 10)

        # Buried surface area (scaled: 1500+ Å² = max score)
        dsasa = metrics.get("dSASA_int")
        if dsasa:
            score_breakdown["interface_burial"] = min(10, dsasa / 1500 * 10)

        # Packing quality (already 0-1)
        packstat = metrics.get("packstat", 0)
        score_breakdown["interface_packing"] = packstat * 5

        # Composition bonus
        composition = metrics.get("interface_composition", {})
        hydrophobic = composition.get("hydrophobic", 0)
        # Hydrophobic burial bonus (up to 10 hydrophobic residues)
        score_breakdown["composition_bonus"] = min(10, hydrophobic)

    # GNINA scoring
    if gnina_result and gnina_result.get("status") == "completed":
        result = gnina_result.get("result", {})

        # Affinity (scaled: -8 or better = max score)
        affinity = result.get("best_affinity")
        if affinity is not None:
            # More negative = better
            if affinity <= -8:
                score_breakdown["gnina_affinity"] = 20
            elif affinity <= 0:
                # Linear scale from 0 to -8
                score_breakdown["gnina_affinity"] = abs(affinity) / 8 * 20
            else:
                # Positive affinity = poor binding
                score_breakdown["gnina_affinity"] = 0

        # CNN score (scaled: 0.7+ = max score)
        cnn_score = result.get("best_cnn_score")
        if cnn_score is not None:
            score_breakdown["gnina_cnn_score"] = min(20, cnn_score / 0.7 * 20)

    # Calculate total
    total = sum(score_breakdown.values())
    max_total = sum(max_scores.values())

    # Quality assessment
    if total >= 70:
        quality = "excellent"
    elif total >= 50:
        quality = "good"
    elif total >= 30:
        quality = "moderate"
    else:
        quality = "poor"

    return {
        "total": round(total, 1),
        "max_possible": max_total,
        "percentage": round(total / max_total * 100, 1),
        "quality": quality,
        "breakdown": {k: round(v, 1) for k, v in score_breakdown.items()},
        "max_scores": max_scores,
    }


def extract_protein_only(pdb_content: str) -> str:
    """Extract only ATOM records (remove HETATM) from PDB."""
    lines = pdb_content.split('\n')
    protein_lines = []
    for line in lines:
        if line.startswith('ATOM') or line.startswith('TER') or line.startswith('END'):
            protein_lines.append(line)
    return '\n'.join(protein_lines)


def calculate_protein_bounding_box(
    pdb_content: str
) -> Tuple[Optional[Tuple[float, float, float]], float]:
    """
    Calculate the bounding box of a protein structure.

    Returns:
        (center, max_dimension): Center coordinates and largest dimension
    """
    try:
        lines = pdb_content.split('\n')
        atom_lines = [l for l in lines if l.startswith('ATOM')]

        if not atom_lines:
            return None, 0

        coords = []
        for line in atom_lines:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords.append([x, y, z])

        coords = np.array(coords)
        min_coords = coords.min(axis=0)
        max_coords = coords.max(axis=0)

        center = (min_coords + max_coords) / 2
        dimensions = max_coords - min_coords
        max_dim = float(dimensions.max())

        return tuple(center.tolist()), max_dim

    except Exception:
        return None, 0


def calculate_interface_center(
    pdb_content: str,
    chain_a: str = "A",
    chain_b: str = "B"
) -> Optional[Tuple[float, float, float]]:
    """
    Calculate the center of the interface between two chains.

    For symmetric dimers, this is typically near the C2 axis (origin).

    Args:
        pdb_content: PDB content as string
        chain_a: First chain identifier
        chain_b: Second chain identifier

    Returns:
        (x, y, z) center coordinates, or None if calculation fails
    """
    try:
        lines = pdb_content.split('\n')
        atom_lines = [l for l in lines if l.startswith('ATOM')]

        coords_a = []
        coords_b = []

        for line in atom_lines:
            chain = line[21]
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            if chain == chain_a:
                coords_a.append((x, y, z))
            elif chain == chain_b:
                coords_b.append((x, y, z))

        if not coords_a or not coords_b:
            return None

        # Calculate centroids
        centroid_a = np.array(coords_a).mean(axis=0)
        centroid_b = np.array(coords_b).mean(axis=0)

        # Interface center is midpoint between chain centroids
        center = (centroid_a + centroid_b) / 2

        return tuple(center.tolist())

    except Exception:
        return None


# ============== Steric Clash Detection ==============

def check_steric_clashes(
    pdb_content: str,
    ligand_smiles: Optional[str] = None,
    ligand_sdf: Optional[str] = None,
    clash_threshold: float = 2.0,
) -> Dict[str, Any]:
    """
    Check for steric clashes between protein and ligand.

    A clash is defined as any protein-ligand atom pair with distance < threshold.
    This is a CRITICAL validation step - positive GNINA affinity typically means
    the ligand was placed inside the protein with no real binding pocket.

    Args:
        pdb_content: PDB content with protein (may include HETATM ligand)
        ligand_smiles: SMILES string for ligand (if not in PDB)
        ligand_sdf: SDF content for ligand (if not in PDB)
        clash_threshold: Distance below which atoms are considered clashing (default 2.0 Å)

    Returns:
        Dict with clash analysis:
        - has_clashes: bool
        - min_distance: float
        - clash_count: int
        - clash_details: list of clashing atom pairs
        - recommendation: str
    """
    try:
        from biotite.structure.io.pdb import PDBFile
        import io

        # Parse PDB
        pdb_file = PDBFile.read(io.StringIO(pdb_content))
        structure = pdb_file.get_structure(model=1)

        # Get protein atoms (ATOM records)
        protein_atoms = structure[~structure.hetero]

        # Get ligand atoms (HETATM records)
        ligand_atoms = structure[structure.hetero]

        # If no HETATM in PDB, try to parse from SMILES/SDF
        if len(ligand_atoms) == 0:
            if ligand_smiles or ligand_sdf:
                # Need to generate ligand coordinates and check against protein
                # For now, return warning that ligand wasn't found
                return {
                    "status": "warning",
                    "has_clashes": None,
                    "message": "Ligand not found in PDB - provide ligand in PDB format for clash detection"
                }
            else:
                return {
                    "status": "error",
                    "error": "No ligand found in PDB and no SMILES/SDF provided"
                }

        if len(protein_atoms) == 0:
            return {
                "status": "error",
                "error": "No protein atoms found in PDB"
            }

        # Calculate all protein-ligand distances
        clashes = []
        min_distance = float('inf')

        for i, ligand_atom in enumerate(ligand_atoms):
            for j, protein_atom in enumerate(protein_atoms):
                dist = np.linalg.norm(ligand_atom.coord - protein_atom.coord)

                if dist < min_distance:
                    min_distance = dist

                if dist < clash_threshold:
                    clashes.append({
                        "ligand_atom": ligand_atom.atom_name,
                        "ligand_res": f"{ligand_atom.res_name}{ligand_atom.res_id}",
                        "protein_atom": protein_atom.atom_name,
                        "protein_res": f"{protein_atom.res_name}{protein_atom.res_id}",
                        "chain": protein_atom.chain_id,
                        "distance": round(dist, 3)
                    })

        has_clashes = len(clashes) > 0

        # Generate recommendation
        if has_clashes:
            if min_distance < 1.0:
                recommendation = "CRITICAL: Severe atomic overlaps detected. Ligand is inside protein - no binding pocket exists. Use ligand-first generation."
            elif min_distance < 1.5:
                recommendation = "SEVERE: Major steric clashes. Binding pocket is blocked. Consider redesigning pocket or using ligand-conditioned generation."
            else:
                recommendation = "WARNING: Minor clashes detected. May need pocket refinement or ligand repositioning."
        else:
            if min_distance < 3.0:
                recommendation = "Good: No clashes, but tight fit. Check for favorable interactions."
            elif min_distance < 5.0:
                recommendation = "Acceptable: Some space for ligand mobility. Verify binding affinity with GNINA."
            else:
                recommendation = "Warning: Ligand may be too far from protein for effective binding."

        return {
            "status": "completed",
            "has_clashes": has_clashes,
            "clash_count": len(clashes),
            "min_distance": round(min_distance, 3),
            "avg_min_distance": round(np.mean([c["distance"] for c in clashes[:10]]), 3) if clashes else None,
            "clash_details": clashes[:20],  # Top 20 clashes
            "recommendation": recommendation,
            "binding_pocket_exists": not has_clashes and min_distance > 2.0 and min_distance < 6.0
        }

    except ImportError as e:
        return {
            "status": "error",
            "error": f"Missing dependency: {e}"
        }
    except Exception as e:
        return {
            "status": "error",
            "error": str(e)
        }


def validate_binding_comprehensive(
    pdb_content: str,
    ligand_smiles: str,
    chain_a: str = "A",
    chain_b: str = "B",
) -> Dict[str, Any]:
    """
    Comprehensive binding validation that checks for clashes FIRST.

    This is the recommended validation function that:
    1. Checks for steric clashes (fatal if present)
    2. Runs interface analysis
    3. Runs GNINA docking (if no clashes)
    4. Calculates composite score

    Args:
        pdb_content: PDB content with protein and ligand
        ligand_smiles: SMILES string for the ligand
        chain_a: First chain identifier
        chain_b: Second chain identifier

    Returns:
        Dict with comprehensive validation results
    """
    results = {
        "status": "completed",
        "clash_check": None,
        "interface_analysis": None,
        "gnina_scoring": None,
        "composite_score": None,
        "recommendation": None
    }

    # Step 1: Check for steric clashes FIRST
    clash_result = check_steric_clashes(pdb_content)
    results["clash_check"] = clash_result

    # If severe clashes, abort further analysis
    if clash_result.get("status") == "completed" and clash_result.get("has_clashes"):
        min_dist = clash_result.get("min_distance", 0)
        if min_dist < 1.5:
            results["status"] = "failed"
            results["fatal_error"] = "Severe steric clashes detected - no binding pocket exists"
            results["recommendation"] = clash_result.get("recommendation")
            return results

    # Step 2: Run interface analysis
    interface_result = analyze_interface(
        pdb_content=pdb_content,
        chain_a=chain_a,
        chain_b=chain_b
    )
    results["interface_analysis"] = interface_result

    # Step 3: Run GNINA scoring (only if pocket likely exists)
    pocket_exists = clash_result.get("binding_pocket_exists", False)
    if pocket_exists or not clash_result.get("has_clashes", True):
        ligand_sdf = smiles_to_sdf(ligand_smiles)
        if ligand_sdf:
            receptor_pdb = extract_protein_only(pdb_content)
            gnina_result = run_gnina_scoring(
                receptor_pdb=receptor_pdb,
                ligand_sdf=ligand_sdf,
                whole_protein=True,  # Search whole surface
            )
            results["gnina_scoring"] = gnina_result
    else:
        results["gnina_scoring"] = {
            "status": "skipped",
            "reason": "Clashes detected - GNINA scoring would be invalid"
        }

    # Step 4: Calculate composite score
    results["composite_score"] = calculate_composite_score(
        interface_result,
        results.get("gnina_scoring") if results.get("gnina_scoring", {}).get("status") == "completed" else None
    )

    # Generate final recommendation
    if clash_result.get("has_clashes"):
        results["recommendation"] = "Use ligand-first generation approach to create proper binding pocket"
    elif results.get("gnina_scoring", {}).get("result", {}).get("best_affinity", 0) > 0:
        results["recommendation"] = "Positive GNINA affinity suggests poor binding. Redesign pocket using theozyme approach."
    else:
        results["recommendation"] = "Design passes basic validation. Proceed to experimental testing."

    return to_python_types(results)


# ============== Quality Thresholds ==============

BINDING_QUALITY_THRESHOLDS = {
    "interface": {
        "contacts": {"good": 50, "moderate": 20, "poor": 0},
        "hbonds": {"good": 4, "moderate": 2, "poor": 0},
        "packstat": {"good": 0.5, "moderate": 0.3, "poor": 0},
        "dSASA": {"good": 1000, "moderate": 500, "poor": 0}  # Angstrom^2
    },
    "gnina": {
        "affinity": {"strong": -8, "moderate": -6, "weak": -4},  # kcal/mol
        "cnn_score": {"good": 0.7, "moderate": 0.5, "poor": 0.3}
    },
    "rf3": {
        "ipTM": {"good": 0.6, "moderate": 0.4, "poor": 0.2},
        "pLDDT": {"good": 0.8, "moderate": 0.7, "poor": 0.5}
    }
}
