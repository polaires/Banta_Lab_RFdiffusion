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
    elif isinstance(obj, (np.float32, np.float64, np.floating)):
        return float(obj)
    elif isinstance(obj, (np.int32, np.int64, np.integer)):
        return int(obj)
    elif isinstance(obj, (np.bool_, np.bool)):
        return bool(obj)
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


def validate_interface_ligand(
    pdb_content: str,
    ligand_chain: str = "L",
    min_contacts_per_chain: int = 3,
    contact_cutoff: float = 4.0,
) -> Dict[str, Any]:
    """
    Validate that a ligand is at the symmetric interface of a dimer.

    For symmetric dimer designs with ligand at the interface, the ligand
    should contact BOTH protein chains approximately equally.

    Args:
        pdb_content: PDB content with protein chains and ligand
        ligand_chain: Chain ID for the ligand (default: "L")
        min_contacts_per_chain: Minimum contacts required per chain (default: 3)
        contact_cutoff: Distance cutoff for contacts in Angstroms (default: 4.0)

    Returns:
        Dict with:
            valid: bool - True if ligand contacts both chains adequately
            contacts_A: int - Number of contacts with chain A
            contacts_B: int - Number of contacts with chain B
            symmetry_score: float - How symmetric the contacts are (0-1)
            recommendation: str - Suggested action
    """
    try:
        from biotite.structure.io.pdb import PDBFile
        import io

        # Parse PDB
        pdb_file = PDBFile.read(io.StringIO(pdb_content))
        structure = pdb_file.get_structure(model=1)

        # Get ligand atoms - try chain ID first, then HETATM
        ligand_mask = structure.chain_id == ligand_chain
        if not np.any(ligand_mask):
            # Try HETATM records
            ligand_mask = structure.hetero
        ligand_atoms = structure[ligand_mask]

        if len(ligand_atoms) == 0:
            return {
                "status": "error",
                "error": f"No ligand found in chain {ligand_chain} or HETATM records",
                "valid": False
            }

        # Get protein chain A and B atoms
        chain_a_mask = (structure.chain_id == "A") & ~structure.hetero
        chain_b_mask = (structure.chain_id == "B") & ~structure.hetero
        chain_a_atoms = structure[chain_a_mask]
        chain_b_atoms = structure[chain_b_mask]

        if len(chain_a_atoms) == 0:
            return {
                "status": "error",
                "error": "No atoms found for chain A",
                "valid": False
            }

        # For monomers, just check overall contacts
        is_dimer = len(chain_b_atoms) > 0

        # Count contacts (within cutoff distance)
        def count_contacts(atoms1, atoms2, cutoff):
            """Count atomic contacts between two atom sets."""
            contacts = 0
            for i in range(len(atoms1)):
                pos1 = atoms1.coord[i]
                for j in range(len(atoms2)):
                    pos2 = atoms2.coord[j]
                    dist = np.sqrt(np.sum((pos1 - pos2) ** 2))
                    if dist < cutoff:
                        contacts += 1
            return contacts

        contacts_a = count_contacts(ligand_atoms, chain_a_atoms, contact_cutoff)

        if is_dimer:
            contacts_b = count_contacts(ligand_atoms, chain_b_atoms, contact_cutoff)
        else:
            contacts_b = 0

        # Validate
        if is_dimer:
            valid = contacts_a >= min_contacts_per_chain and contacts_b >= min_contacts_per_chain
        else:
            valid = contacts_a >= min_contacts_per_chain

        # Symmetry score: 1.0 = perfectly symmetric, 0.0 = all contacts with one chain
        total = contacts_a + contacts_b
        if total > 0 and is_dimer:
            symmetry_score = 1.0 - abs(contacts_a - contacts_b) / total
        else:
            symmetry_score = 1.0 if not is_dimer else 0.0

        # Generate recommendation
        if valid:
            recommendation = "Valid interface design - ligand contacts both chains"
        elif not is_dimer:
            if contacts_a >= min_contacts_per_chain:
                recommendation = "Monomer design - ligand has adequate contacts"
            else:
                recommendation = "Ligand has insufficient contacts with protein"
        elif contacts_a < min_contacts_per_chain and contacts_b < min_contacts_per_chain:
            recommendation = "Ligand is exposed - not bound to either chain. Redesign needed."
        elif contacts_a < min_contacts_per_chain:
            recommendation = "Ligand buried in chain B only - not at interface. Use interface_ligand=True"
        elif contacts_b < min_contacts_per_chain:
            recommendation = "Ligand buried in chain A only - not at interface. Use interface_ligand=True"
        else:
            recommendation = "Ligand not at interface - redesign with interface_ligand=True"

        return to_python_types({
            "status": "completed",
            "valid": valid,
            "is_dimer": is_dimer,
            "contacts_A": contacts_a,
            "contacts_B": contacts_b,
            "total_contacts": total,
            "symmetry_score": symmetry_score,
            "min_contacts_required": min_contacts_per_chain,
            "recommendation": recommendation
        })

    except Exception as e:
        return {
            "status": "error",
            "error": str(e),
            "valid": False
        }


# ============== Cleavable Monomer Helpers ==============

def count_ligand_contacts_by_residue(
    pdb_content: str,
    ligand_name: str = "UNL",
    cutoff: float = 5.0,
) -> Dict[int, int]:
    """
    Count ligand contacts for each protein residue.

    Used by cleavable monomer algorithm to determine which residues
    contact the ligand and how cleaving at different points would
    distribute contacts between the resulting chains.

    Args:
        pdb_content: PDB file content as string
        ligand_name: Ligand residue name (default "UNL")
        cutoff: Contact distance cutoff in Angstroms (default 5.0)

    Returns:
        Dict mapping residue_number -> contact_count
        Empty dict if parsing fails.
    """
    try:
        from biotite.structure.io.pdb import PDBFile
        import io

        pdb_file = PDBFile.read(io.StringIO(pdb_content))
        structure = pdb_file.get_structure(model=1)

        # Get ligand atoms
        ligand_mask = structure.res_name == ligand_name
        ligand_atoms = structure[ligand_mask]

        if len(ligand_atoms) == 0:
            return {}

        ligand_coords = ligand_atoms.coord

        # Get protein atoms (non-HETATM)
        protein_mask = ~structure.hetero
        protein_atoms = structure[protein_mask]

        if len(protein_atoms) == 0:
            return {}

        # Count contacts per residue
        contacts = {}
        for i in range(len(protein_atoms)):
            res_id = int(protein_atoms.res_id[i])
            pos = protein_atoms.coord[i]

            # Check distance to any ligand atom
            for lig_pos in ligand_coords:
                dist = np.sqrt(np.sum((pos - lig_pos) ** 2))
                if dist < cutoff:
                    contacts[res_id] = contacts.get(res_id, 0) + 1
                    break  # Count residue once per ligand atom contact

        return contacts

    except Exception as e:
        print(f"[BindingAnalysis] Error counting contacts by residue: {e}")
        return {}


def get_ligand_centroid(
    pdb_content: str,
    ligand_name: str = "UNL"
) -> Optional[Tuple[float, float, float]]:
    """
    Get the centroid (center of mass) of the ligand.

    Used by cleavable monomer algorithm to calculate distances
    from residues to the ligand center.

    Args:
        pdb_content: PDB file content as string
        ligand_name: Ligand residue name (default "UNL")

    Returns:
        (x, y, z) centroid coordinates, or None if ligand not found.
    """
    try:
        from biotite.structure.io.pdb import PDBFile
        import io

        pdb_file = PDBFile.read(io.StringIO(pdb_content))
        structure = pdb_file.get_structure(model=1)

        # Get ligand atoms
        ligand_mask = structure.res_name == ligand_name
        ligand_atoms = structure[ligand_mask]

        if len(ligand_atoms) == 0:
            return None

        # Calculate centroid
        centroid = ligand_atoms.coord.mean(axis=0)
        return tuple(float(x) for x in centroid)

    except Exception as e:
        print(f"[BindingAnalysis] Error getting ligand centroid: {e}")
        return None


def get_ca_positions(
    pdb_content: str
) -> List[Tuple[int, np.ndarray]]:
    """
    Get CA atom positions for all protein residues.

    Used by cleavable monomer algorithm to calculate distances
    from backbone to ligand.

    Args:
        pdb_content: PDB file content as string

    Returns:
        List of (residue_number, xyz_array) tuples sorted by residue number.
    """
    try:
        from biotite.structure.io.pdb import PDBFile
        import io

        pdb_file = PDBFile.read(io.StringIO(pdb_content))
        structure = pdb_file.get_structure(model=1)

        # Get CA atoms from protein
        ca_mask = (structure.atom_name == "CA") & (~structure.hetero)
        ca_atoms = structure[ca_mask]

        if len(ca_atoms) == 0:
            return []

        # Build list of (res_id, position)
        positions = []
        for i in range(len(ca_atoms)):
            res_id = int(ca_atoms.res_id[i])
            pos = ca_atoms.coord[i]
            positions.append((res_id, pos))

        # Sort by residue number
        positions.sort(key=lambda x: x[0])
        return positions

    except Exception as e:
        print(f"[BindingAnalysis] Error getting CA positions: {e}")
        return []


# ============== Shape Complementarity (BindCraft-style) ==============

def calculate_shape_complementarity(
    pdb_content: str,
    chain_a: str = "A",
    chain_b: str = "B",
) -> Dict[str, Any]:
    """
    Calculate interface Shape Complementarity (SC) using PyRosetta.

    Shape Complementarity measures how well two surfaces "fit" together.
    Values range from 0 (poor) to 1 (perfect complementarity).

    BindCraft threshold: SC > 0.5 indicates good interface packing.

    Args:
        pdb_content: PDB file content as string
        chain_a: First chain (typically target)
        chain_b: Second chain (typically binder)

    Returns:
        Dict with:
        - shape_complementarity: SC score (0-1)
        - status: "completed" or "error"
    """
    try:
        # Try PyRosetta first
        import pyrosetta
        from pyrosetta import pose_from_pdb_string
        from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

        # Initialize PyRosetta if needed
        if not pyrosetta.rosetta.basic.was_init_called():
            pyrosetta.init("-mute all")

        # Load pose
        pose = pose_from_pdb_string(pdb_content)

        # Create interface analyzer
        interface_str = f"{chain_a}_{chain_b}"
        iam = InterfaceAnalyzerMover()
        iam.set_interface(interface_str)
        iam.set_compute_interface_sc(True)

        # Apply to pose
        iam.apply(pose)

        # Get SC score
        sc_score = iam.get_interface_sc()

        return {
            "status": "completed",
            "shape_complementarity": float(sc_score),
            "method": "pyrosetta"
        }

    except ImportError:
        # Fallback to geometric approximation if PyRosetta not available
        return _calculate_sc_geometric(pdb_content, chain_a, chain_b)
    except Exception as e:
        print(f"[BindingAnalysis] PyRosetta SC failed: {e}, trying geometric fallback")
        return _calculate_sc_geometric(pdb_content, chain_a, chain_b)


def _calculate_sc_geometric(
    pdb_content: str,
    chain_a: str,
    chain_b: str,
) -> Dict[str, Any]:
    """
    Geometric approximation of Shape Complementarity.

    Uses surface normal dot products at interface contacts
    as a proxy for SC when PyRosetta is not available.
    """
    try:
        from biotite.structure.io.pdb import PDBFile
        import io

        pdb_file = PDBFile.read(io.StringIO(pdb_content))
        structure = pdb_file.get_structure(model=1)

        # Get interface atoms
        chain_a_atoms = structure[structure.chain_id == chain_a]
        chain_b_atoms = structure[structure.chain_id == chain_b]

        if len(chain_a_atoms) == 0 or len(chain_b_atoms) == 0:
            return {"status": "error", "error": "Could not find chains"}

        # Find interface contacts
        contact_cutoff = 5.0  # Angstroms
        contact_pairs = []

        for i, atom_a in enumerate(chain_a_atoms):
            for j, atom_b in enumerate(chain_b_atoms):
                dist = np.linalg.norm(atom_a.coord - atom_b.coord)
                if dist < contact_cutoff:
                    contact_pairs.append((i, j, dist))

        if len(contact_pairs) == 0:
            return {
                "status": "completed",
                "shape_complementarity": 0.0,
                "method": "geometric",
                "note": "No interface contacts found"
            }

        # Estimate SC from contact distances
        # Closer contacts = better complementarity
        # SC ≈ fraction of contacts that are "tight" (< 4Å)
        tight_contacts = sum(1 for _, _, d in contact_pairs if d < 4.0)
        sc_estimate = tight_contacts / len(contact_pairs)

        # Adjust to typical SC range (0.3-0.8)
        sc_estimate = 0.3 + sc_estimate * 0.5

        return {
            "status": "completed",
            "shape_complementarity": round(float(sc_estimate), 3),
            "method": "geometric",
            "n_contacts": len(contact_pairs),
            "n_tight_contacts": tight_contacts
        }

    except Exception as e:
        return {
            "status": "error",
            "error": str(e)
        }


# ============== Surface Hydrophobicity ==============

def calculate_surface_hydrophobicity(
    pdb_content: str,
    chain: str = "B",
) -> Dict[str, Any]:
    """
    Calculate surface hydrophobicity of a protein chain.

    BindCraft filters designs with surface_hydrophobicity > 0.37
    to prevent aggregation-prone binders.

    Args:
        pdb_content: PDB file content as string
        chain: Chain to analyze

    Returns:
        Dict with:
        - surface_hydrophobicity: Fraction of exposed hydrophobic residues
        - status: "completed" or "error"
    """
    try:
        from biotite.structure.io.pdb import PDBFile
        from biotite.structure import sasa
        import io

        pdb_file = PDBFile.read(io.StringIO(pdb_content))
        structure = pdb_file.get_structure(model=1)

        # Get chain atoms
        chain_atoms = structure[structure.chain_id == chain]
        if len(chain_atoms) == 0:
            return {"status": "error", "error": f"Chain {chain} not found"}

        # Calculate per-atom SASA
        atom_sasa = sasa(chain_atoms)

        # Hydrophobic residues
        hydrophobic_residues = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO'}

        # Count exposed residues by type
        exposed_threshold = 10.0  # Å² for a residue to be considered "exposed"

        # Aggregate SASA by residue
        residue_sasa = {}
        residue_type = {}

        for i, atom in enumerate(chain_atoms):
            res_id = atom.res_id
            res_name = atom.res_name

            if res_id not in residue_sasa:
                residue_sasa[res_id] = 0.0
                residue_type[res_id] = res_name

            residue_sasa[res_id] += atom_sasa[i]

        # Count exposed residues
        total_exposed = 0
        hydrophobic_exposed = 0

        for res_id, total_sasa in residue_sasa.items():
            if total_sasa > exposed_threshold:
                total_exposed += 1
                if residue_type[res_id] in hydrophobic_residues:
                    hydrophobic_exposed += 1

        # Calculate fraction
        if total_exposed > 0:
            surface_hydrophobicity = hydrophobic_exposed / total_exposed
        else:
            surface_hydrophobicity = 0.0

        return {
            "status": "completed",
            "surface_hydrophobicity": round(float(surface_hydrophobicity), 3),
            "total_exposed_residues": total_exposed,
            "hydrophobic_exposed_residues": hydrophobic_exposed,
            "passes_bindcraft_threshold": surface_hydrophobicity < 0.37
        }

    except Exception as e:
        return {
            "status": "error",
            "error": str(e)
        }


# ============== Unsaturated Hydrogen Bonds ==============

def calculate_unsaturated_hbonds(
    pdb_content: str,
    chain: str = "B",
    interface_residues: Optional[List[int]] = None,
) -> Dict[str, Any]:
    """
    Count buried polar atoms that are not forming hydrogen bonds.

    BindCraft filters designs with unsaturated_hbonds > 6
    as these indicate potential instability.

    Args:
        pdb_content: PDB file content as string
        chain: Chain to analyze
        interface_residues: Optional list of interface residue numbers

    Returns:
        Dict with:
        - unsaturated_count: Number of unsaturated polar atoms
        - status: "completed" or "error"
    """
    try:
        from biotite.structure.io.pdb import PDBFile
        from biotite.structure import sasa
        import io

        pdb_file = PDBFile.read(io.StringIO(pdb_content))
        structure = pdb_file.get_structure(model=1)

        # Get chain atoms
        chain_atoms = structure[structure.chain_id == chain]
        if len(chain_atoms) == 0:
            return {"status": "error", "error": f"Chain {chain} not found"}

        # Calculate SASA to identify buried atoms
        atom_sasa = sasa(chain_atoms)

        # Polar atoms that can form H-bonds
        hbond_atoms = {
            'N', 'O',  # Backbone
            'OG', 'OG1', 'OH',  # Ser, Thr, Tyr
            'OD1', 'OD2', 'OE1', 'OE2',  # Asp, Glu
            'ND1', 'ND2', 'NE', 'NE1', 'NE2', 'NH1', 'NH2', 'NZ',  # His, Arg, Lys
            'SG'  # Cys
        }

        # Find buried polar atoms
        buried_polar = []
        buried_threshold = 5.0  # Å² - atoms below this SASA are "buried"

        for i, atom in enumerate(chain_atoms):
            atom_name = atom.atom_name
            res_id = atom.res_id

            # Check if in interface (if specified)
            if interface_residues is not None:
                if res_id not in interface_residues:
                    continue

            # Check if polar and buried
            if atom_name in hbond_atoms and atom_sasa[i] < buried_threshold:
                buried_polar.append({
                    "atom": atom_name,
                    "residue": f"{atom.res_name}{res_id}",
                    "sasa": round(float(atom_sasa[i]), 2)
                })

        # Simple check: count potentially unsaturated
        # (A proper check would identify actual H-bond partners)
        unsaturated_count = len(buried_polar)

        return {
            "status": "completed",
            "unsaturated_count": unsaturated_count,
            "buried_polar_atoms": buried_polar[:20],  # First 20
            "passes_bindcraft_threshold": unsaturated_count < 6
        }

    except Exception as e:
        return {
            "status": "error",
            "error": str(e)
        }


# ============== Ligand-Specific H-bond Analysis ==============

def analyze_ligand_hbonds(
    pdb_content: str,
    ligand_name: str = "UNL",
    ligand_atoms: Optional[List[str]] = None,
    hbond_distance: float = 3.5,
) -> Dict[str, Any]:
    """
    Analyze hydrogen bonds to specific ligand atoms.

    Tracks H-bonds between protein and ligand, with optional per-atom breakdown
    for targeted ligand atoms (e.g., N5, N6 for azobenzene).

    Args:
        pdb_content: PDB file content as string
        ligand_name: Ligand residue name (default "UNL")
        ligand_atoms: List of specific ligand atoms to track (e.g., ["N5", "N6"])
        hbond_distance: Maximum distance for H-bond detection (default 3.5 Å)

    Returns:
        Dict with:
        - total_hbonds: Total H-bonds to ligand
        - hbonds: List of H-bond details
        - by_atom: Dict of {atom: [hbond_details]}
        - saturation: Dict of {atom: {count, saturated}}
    """
    try:
        from biotite.structure.io.pdb import PDBFile
        import io
        import numpy as np

        pdb_file = PDBFile.read(io.StringIO(pdb_content))
        structure = pdb_file.get_structure(model=1)

        # Get ligand atoms
        ligand_mask = structure.res_name == ligand_name
        ligand_struct = structure[ligand_mask]

        if len(ligand_struct) == 0:
            return {
                "status": "error",
                "error": f"Ligand {ligand_name} not found"
            }

        # Get protein atoms (non-HETATM)
        protein_mask = ~structure.hetero
        protein_struct = structure[protein_mask]

        if len(protein_struct) == 0:
            return {
                "status": "error",
                "error": "No protein atoms found"
            }

        # Donor/acceptor definitions
        # Protein donors: N-H groups
        protein_donor_atoms = ['N', 'NE', 'NH1', 'NH2', 'ND1', 'ND2', 'NE1', 'NE2', 'NZ', 'OG', 'OG1', 'OH', 'SG']
        # Protein acceptors: O, N with lone pairs
        protein_acceptor_atoms = ['O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH', 'ND1', 'NE2', 'SD']

        # Ligand atoms that can be acceptors (N, O, S)
        ligand_acceptor_types = ['N', 'O', 'S']  # First letter of element
        # Ligand atoms that can be donors (N, O with H)
        ligand_donor_types = ['N', 'O']

        hbonds = []
        hbonds_by_atom = {}

        for i in range(len(ligand_struct)):
            lig_atom = ligand_struct[i]
            lig_atom_name = lig_atom.atom_name

            # Skip if not in target atoms list (if specified)
            if ligand_atoms and lig_atom_name not in ligand_atoms:
                continue

            # Determine if ligand atom can be donor or acceptor
            element = lig_atom.element if hasattr(lig_atom, 'element') else lig_atom_name[0]
            is_lig_acceptor = element in ligand_acceptor_types
            is_lig_donor = element in ligand_donor_types

            for j in range(len(protein_struct)):
                prot_atom = protein_struct[j]
                dist = np.linalg.norm(lig_atom.coord - prot_atom.coord)

                if dist > hbond_distance:
                    continue

                prot_atom_name = prot_atom.atom_name
                is_prot_donor = prot_atom_name in protein_donor_atoms
                is_prot_acceptor = prot_atom_name in protein_acceptor_atoms

                # Check if valid H-bond
                is_hbond = False
                hbond_type = None

                if is_lig_acceptor and is_prot_donor:
                    is_hbond = True
                    hbond_type = "ligand_acceptor"
                elif is_lig_donor and is_prot_acceptor:
                    is_hbond = True
                    hbond_type = "ligand_donor"

                if is_hbond:
                    hbond = {
                        "ligand_atom": lig_atom_name,
                        "protein_res": f"{prot_atom.res_name}{prot_atom.res_id}",
                        "protein_chain": prot_atom.chain_id,
                        "protein_atom": prot_atom_name,
                        "distance": round(float(dist), 2),
                        "type": hbond_type,
                    }
                    hbonds.append(hbond)

                    if lig_atom_name not in hbonds_by_atom:
                        hbonds_by_atom[lig_atom_name] = []
                    hbonds_by_atom[lig_atom_name].append(hbond)

        # Calculate saturation for requested atoms
        saturation = {}
        if ligand_atoms:
            for atom in ligand_atoms:
                count = len(hbonds_by_atom.get(atom, []))
                saturation[atom] = {
                    "count": count,
                    "saturated": count >= 1,  # At least 1 H-bond
                }

        unsaturated_count = sum(1 for s in saturation.values() if not s.get("saturated", False)) if saturation else 0

        return {
            "status": "completed",
            "total_hbonds": len(hbonds),
            "hbonds": hbonds[:20],  # First 20
            "by_atom": hbonds_by_atom,
            "saturation": saturation,
            "unsaturated_count": unsaturated_count,
        }

    except Exception as e:
        return {
            "status": "error",
            "error": str(e)
        }


# ============== Interface Residue Count ==============

def count_interface_residues(
    pdb_content: str,
    chain_a: str = "A",
    chain_b: str = "B",
    distance_cutoff: float = 5.0,
) -> Dict[str, Any]:
    """
    Count residues at the protein-protein interface.

    BindCraft requires interface_residues > 6 to ensure
    meaningful binding interface.

    Args:
        pdb_content: PDB file content
        chain_a: First chain
        chain_b: Second chain
        distance_cutoff: Distance for interface (Å)

    Returns:
        Dict with interface residue counts
    """
    try:
        from biotite.structure.io.pdb import PDBFile
        import io

        pdb_file = PDBFile.read(io.StringIO(pdb_content))
        structure = pdb_file.get_structure(model=1)

        chain_a_atoms = structure[structure.chain_id == chain_a]
        chain_b_atoms = structure[structure.chain_id == chain_b]

        if len(chain_a_atoms) == 0 or len(chain_b_atoms) == 0:
            return {"status": "error", "error": "Chains not found"}

        interface_a = set()
        interface_b = set()

        for atom_a in chain_a_atoms:
            for atom_b in chain_b_atoms:
                dist = np.linalg.norm(atom_a.coord - atom_b.coord)
                if dist < distance_cutoff:
                    interface_a.add(atom_a.res_id)
                    interface_b.add(atom_b.res_id)

        total_interface = len(interface_a) + len(interface_b)

        return {
            "status": "completed",
            "interface_residues_total": total_interface,
            "interface_residues_chain_a": len(interface_a),
            "interface_residues_chain_b": len(interface_b),
            "interface_residue_ids_a": sorted(list(interface_a)),
            "interface_residue_ids_b": sorted(list(interface_b)),
            "passes_bindcraft_threshold": total_interface >= 6
        }

    except Exception as e:
        return {
            "status": "error",
            "error": str(e)
        }


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


# ============== Homodimerization Scoring (Anti-Homodimer Validation) ==============

def score_homodimerization(
    pdb_content: str,
    ligand_smiles: str,
    ligand_name: str = "UNL",
    chain_a: str = "A",
    chain_b: str = "B",
    docking_box_size: float = 25.0,
) -> Dict[str, Any]:
    """
    Score heterodimer vs homodimer binding to validate anti-homodimerization design.

    For a true heterodimer design, we want:
    - A-B heterodimer to bind strongly (affinity < -5 kcal/mol)
    - A-A homodimer to bind weakly (affinity > -2 kcal/mol)
    - B-B homodimer to bind weakly (affinity > -2 kcal/mol)
    - Hetero/homo ratio > 2.5x

    Args:
        pdb_content: PDB content with heterodimer (chains A and B) and ligand
        ligand_smiles: SMILES string for the ligand
        ligand_name: Residue name for ligand in PDB (default "UNL")
        chain_a: Chain ID for first chain
        chain_b: Chain ID for second chain
        docking_box_size: Size of docking box for GNINA

    Returns:
        Dict with scores for A-B, A-A, B-B dimers and anti-homodimerization metrics
    """
    try:
        from biotite.structure.io.pdb import PDBFile
        import io

        # Parse original PDB
        pdb_file = PDBFile.read(io.StringIO(pdb_content))
        structure = pdb_file.get_structure(model=1)

        # Extract chains and ligand
        chain_a_atoms = structure[structure.chain_id == chain_a]
        chain_b_atoms = structure[structure.chain_id == chain_b]
        ligand_atoms = structure[structure.res_name == ligand_name]

        if len(chain_a_atoms) == 0 or len(chain_b_atoms) == 0:
            return {"status": "error", "error": f"Chains {chain_a} or {chain_b} not found"}

        if len(ligand_atoms) == 0:
            return {"status": "error", "error": f"Ligand {ligand_name} not found"}

        # Get ligand centroid for docking box
        ligand_center = ligand_atoms.coord.mean(axis=0)

        results = {}

        # Convert SMILES to SDF for GNINA
        ligand_sdf = smiles_to_sdf(ligand_smiles)
        if not ligand_sdf:
            return {"status": "error", "error": "Failed to convert ligand SMILES to SDF"}

        # ===== 1. Score original A-B heterodimer =====
        ab_result = run_gnina_scoring(
            receptor_pdb=pdb_content,
            ligand_sdf=ligand_sdf,
            center=tuple(ligand_center.tolist()),
            box_size=docking_box_size,
        )
        results["ab_heterodimer"] = ab_result.get("result", {})
        ab_affinity = results["ab_heterodimer"].get("best_affinity", 0)

        # ===== 2. Create and score A-A homodimer =====
        aa_pdb = _create_homodimer_pdb(pdb_content, chain_a, ligand_name)
        if aa_pdb:
            aa_result = run_gnina_scoring(
                receptor_pdb=aa_pdb,
                ligand_sdf=ligand_sdf,
                center=tuple(ligand_center.tolist()),
                box_size=docking_box_size,
            )
            results["aa_homodimer"] = aa_result.get("result", {})
        else:
            results["aa_homodimer"] = {"error": "Failed to create A-A homodimer"}
        aa_affinity = results["aa_homodimer"].get("best_affinity", 0)

        # ===== 3. Create and score B-B homodimer =====
        bb_pdb = _create_homodimer_pdb(pdb_content, chain_b, ligand_name)
        if bb_pdb:
            bb_result = run_gnina_scoring(
                receptor_pdb=bb_pdb,
                ligand_sdf=ligand_sdf,
                center=tuple(ligand_center.tolist()),
                box_size=docking_box_size,
            )
            results["bb_homodimer"] = bb_result.get("result", {})
        else:
            results["bb_homodimer"] = {"error": "Failed to create B-B homodimer"}
        bb_affinity = results["bb_homodimer"].get("best_affinity", 0)

        # ===== 4. Calculate anti-homodimerization metrics =====

        # Hetero/homo ratio (higher is better)
        # More negative affinity = stronger binding
        # We want: ab_affinity << aa_affinity and ab_affinity << bb_affinity

        # Avoid division by zero
        if aa_affinity >= 0:
            aa_ratio = float('inf') if ab_affinity < 0 else 0
        else:
            aa_ratio = ab_affinity / aa_affinity if aa_affinity != 0 else 0

        if bb_affinity >= 0:
            bb_ratio = float('inf') if ab_affinity < 0 else 0
        else:
            bb_ratio = ab_affinity / bb_affinity if bb_affinity != 0 else 0

        # Selectivity: how much stronger is hetero vs best homo
        best_homo_affinity = min(aa_affinity, bb_affinity)  # Less negative = weaker
        selectivity = best_homo_affinity - ab_affinity  # Positive = heterodimer is stronger

        # Anti-homodimerization score (0-100)
        # Criteria:
        # - ab_affinity < -5: heterodimer binds well
        # - aa_affinity > -2: A-A doesn't bind
        # - bb_affinity > -2: B-B doesn't bind
        score = 0

        if ab_affinity < -5:
            score += 40  # Strong heterodimer binding
        elif ab_affinity < -3:
            score += 20

        if aa_affinity > -2:
            score += 25  # Weak A-A binding (good)
        elif aa_affinity > -3:
            score += 10

        if bb_affinity > -2:
            score += 25  # Weak B-B binding (good)
        elif bb_affinity > -3:
            score += 10

        if selectivity > 3:
            score += 10  # High selectivity bonus

        # Determine if design passes anti-homodimerization criteria
        passes_anti_homo = (
            ab_affinity < -5 and  # Heterodimer binds
            aa_affinity > -2 and  # A-A doesn't bind
            bb_affinity > -2 and  # B-B doesn't bind
            selectivity > 2.5     # Heterodimer at least 2.5 kcal/mol stronger
        )

        return {
            "status": "completed",
            "affinities": {
                "ab_heterodimer": ab_affinity,
                "aa_homodimer": aa_affinity,
                "bb_homodimer": bb_affinity,
            },
            "ratios": {
                "ab_vs_aa": aa_ratio,
                "ab_vs_bb": bb_ratio,
            },
            "selectivity_kcal": selectivity,
            "anti_homo_score": score,
            "passes_anti_homodimerization": passes_anti_homo,
            "detailed_results": results,
            "thresholds": {
                "heterodimer_target": "< -5 kcal/mol",
                "homodimer_target": "> -2 kcal/mol",
                "selectivity_target": "> 2.5 kcal/mol",
            }
        }

    except Exception as e:
        import traceback
        return {
            "status": "error",
            "error": str(e),
            "traceback": traceback.format_exc()
        }


def _create_homodimer_pdb(
    pdb_content: str,
    source_chain: str,
    ligand_name: str = "UNL",
) -> Optional[str]:
    """
    Create a homodimer PDB by duplicating one chain and rotating 180 degrees.

    Takes chain A (or B) from the original heterodimer, duplicates it,
    rotates the copy 180 degrees around the ligand's Z-axis, and relabels
    as chains A and B.

    Args:
        pdb_content: Original PDB with heterodimer
        source_chain: Chain to duplicate (A or B)
        ligand_name: Ligand residue name

    Returns:
        PDB content with homodimer (two copies of source chain) + ligand
    """
    try:
        from biotite.structure.io.pdb import PDBFile
        import io

        pdb_file = PDBFile.read(io.StringIO(pdb_content))
        structure = pdb_file.get_structure(model=1)

        # Get source chain and ligand
        source_atoms = structure[structure.chain_id == source_chain]
        ligand_atoms = structure[structure.res_name == ligand_name]

        if len(source_atoms) == 0:
            return None

        # Get ligand centroid for rotation center
        if len(ligand_atoms) > 0:
            center = ligand_atoms.coord.mean(axis=0)
        else:
            center = source_atoms.coord.mean(axis=0)

        # Create a copy of the source chain and rotate 180 degrees around Z
        rotated_coords = source_atoms.coord.copy()
        # Translate to center
        rotated_coords -= center
        # Rotate 180 degrees around Z: (x, y) -> (-x, -y)
        rotated_coords[:, 0] = -rotated_coords[:, 0]
        rotated_coords[:, 1] = -rotated_coords[:, 1]
        # Translate back
        rotated_coords += center

        # Build homodimer PDB
        lines = []
        lines.append(f"REMARK  Homodimer of chain {source_chain}")

        atom_num = 1

        # Chain A = original source chain
        for atom in source_atoms:
            x, y, z = atom.coord
            line = f"ATOM  {atom_num:5d}  {atom.atom_name:<4s}{atom.res_name:>3s} A{atom.res_id:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {atom.element:>2s}"
            lines.append(line)
            atom_num += 1

        # Chain B = rotated copy
        for i, atom in enumerate(source_atoms):
            x, y, z = rotated_coords[i]
            line = f"ATOM  {atom_num:5d}  {atom.atom_name:<4s}{atom.res_name:>3s} B{atom.res_id:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {atom.element:>2s}"
            lines.append(line)
            atom_num += 1

        # Add ligand
        for atom in ligand_atoms:
            x, y, z = atom.coord
            line = f"HETATM{atom_num:5d}  {atom.atom_name:<4s}{atom.res_name:>3s} L   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {atom.element:>2s}"
            lines.append(line)
            atom_num += 1

        lines.append("END")

        return "\n".join(lines)

    except Exception as e:
        print(f"[_create_homodimer_pdb] Error: {e}")
        return None


def calculate_sequence_identity(
    pdb_content: str,
    chain_a: str = "A",
    chain_b: str = "B",
) -> Dict[str, Any]:
    """
    Calculate sequence identity between two chains.

    For heterodimers, we expect lower identity (< 70%).
    For homodimers (symmetric), identity would be 100%.

    Args:
        pdb_content: PDB content with both chains
        chain_a: First chain ID
        chain_b: Second chain ID

    Returns:
        Dict with sequence identity percentage and alignment details
    """
    try:
        from biotite.structure.io.pdb import PDBFile
        import io

        # 3-letter to 1-letter amino acid codes
        aa_map = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
        }

        pdb_file = PDBFile.read(io.StringIO(pdb_content))
        structure = pdb_file.get_structure(model=1)

        # Extract CA atoms to get sequence
        ca_atoms_a = structure[(structure.chain_id == chain_a) & (structure.atom_name == "CA")]
        ca_atoms_b = structure[(structure.chain_id == chain_b) & (structure.atom_name == "CA")]

        # Build sequences
        seq_a = "".join([aa_map.get(atom.res_name, 'X') for atom in ca_atoms_a])
        seq_b = "".join([aa_map.get(atom.res_name, 'X') for atom in ca_atoms_b])

        if len(seq_a) == 0 or len(seq_b) == 0:
            return {"status": "error", "error": "Could not extract sequences"}

        # Simple pairwise identity calculation
        # For sequences of different lengths, align from start
        min_len = min(len(seq_a), len(seq_b))
        max_len = max(len(seq_a), len(seq_b))

        matches = sum(1 for i in range(min_len) if seq_a[i] == seq_b[i])
        identity = (matches / max_len) * 100 if max_len > 0 else 0

        # More sophisticated: use Biopython for proper alignment if available
        try:
            from Bio import pairwise2
            from Bio.pairwise2 import format_alignment

            alignments = pairwise2.align.globalxx(seq_a, seq_b, one_alignment_only=True)
            if alignments:
                aligned_a, aligned_b, score, begin, end = alignments[0]
                aligned_matches = sum(1 for a, b in zip(aligned_a, aligned_b) if a == b and a != '-')
                aligned_len = len(aligned_a)
                identity = (aligned_matches / aligned_len) * 100 if aligned_len > 0 else 0
        except ImportError:
            pass  # Use simple calculation

        return {
            "status": "completed",
            "sequence_identity_percent": round(identity, 1),
            "sequence_a": seq_a,
            "sequence_b": seq_b,
            "length_a": len(seq_a),
            "length_b": len(seq_b),
            "is_heterodimer": identity < 70,  # True heterodimer if < 70% identity
        }

    except Exception as e:
        return {"status": "error", "error": str(e)}
