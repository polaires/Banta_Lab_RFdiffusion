"""
Shared Interaction Analysis Module

Provides unified protein-ligand interaction analysis used by:
- Interface ligand dimer workflow (handler.py)
- Validation pipelines (cleavage_utils.py)
- AI design assistant (ai_engine.py)
- API endpoints (main.py)

Analysis methods (in priority order):
1. ProLIF - Modern, pip-installable, good for fingerprinting
2. PLIP - Comprehensive, requires OpenBabel
3. Biotite distance-based - Fallback when others unavailable

RFdiffusion3 note: The RFD3 paper uses Rosetta DDGnorepack for binding
energy and AlphaFold3 for interface validation, not detailed profiling.
"""

import os
import tempfile
import math
from typing import Dict, Any, List, Optional, Tuple
from dataclasses import dataclass, field
import numpy as np


@dataclass
class HydrogenBond:
    """Hydrogen bond interaction details."""
    ligand_atom: str
    protein_residue: str
    protein_chain: str
    protein_atom: str
    distance: float
    angle: Optional[float] = None  # D-H-A angle in degrees
    donor_coords: Optional[Tuple[float, float, float]] = None
    acceptor_coords: Optional[Tuple[float, float, float]] = None
    type: str = "unknown"  # "ligand_donor" or "ligand_acceptor"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "ligand_atom": self.ligand_atom,
            "protein_residue": self.protein_residue,
            "protein_chain": self.protein_chain,
            "protein_atom": self.protein_atom,
            "distance": self.distance,
            "angle": self.angle,
            "donor_coords": self.donor_coords,
            "acceptor_coords": self.acceptor_coords,
            "type": self.type,
        }


@dataclass
class HydrophobicContact:
    """Hydrophobic contact interaction details."""
    ligand_atom: str
    protein_residue: str
    protein_chain: str
    protein_atom: str
    distance: float
    coords: Optional[Tuple[float, float, float]] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "ligand_atom": self.ligand_atom,
            "protein_residue": self.protein_residue,
            "protein_chain": self.protein_chain,
            "protein_atom": self.protein_atom,
            "distance": self.distance,
            "coords": self.coords,
        }


@dataclass
class PiStacking:
    """Pi-stacking interaction details."""
    ligand_ring: str
    protein_residue: str
    protein_chain: str
    type: str  # "face-to-face" or "edge-to-face"
    distance: float
    angle: Optional[float] = None
    ligand_center: Optional[Tuple[float, float, float]] = None
    protein_center: Optional[Tuple[float, float, float]] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "ligand_ring": self.ligand_ring,
            "protein_residue": self.protein_residue,
            "protein_chain": self.protein_chain,
            "type": self.type,
            "distance": self.distance,
            "angle": self.angle,
            "ligand_center": self.ligand_center,
            "protein_center": self.protein_center,
        }


@dataclass
class SaltBridge:
    """Salt bridge interaction details."""
    ligand_group: str
    protein_residue: str
    protein_chain: str
    distance: float
    ligand_charge: str  # "positive" or "negative"
    protein_charge: str
    ligand_coords: Optional[Tuple[float, float, float]] = None
    protein_coords: Optional[Tuple[float, float, float]] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "ligand_group": self.ligand_group,
            "protein_residue": self.protein_residue,
            "protein_chain": self.protein_chain,
            "distance": self.distance,
            "ligand_charge": self.ligand_charge,
            "protein_charge": self.protein_charge,
            "ligand_coords": self.ligand_coords,
            "protein_coords": self.protein_coords,
        }


@dataclass
class HalogenBond:
    """Halogen bond interaction details."""
    ligand_atom: str
    protein_residue: str
    protein_chain: str
    protein_atom: str
    distance: float
    angle: Optional[float] = None
    halogen_type: str = ""  # Cl, Br, I, F

    def to_dict(self) -> Dict[str, Any]:
        return {
            "ligand_atom": self.ligand_atom,
            "protein_residue": self.protein_residue,
            "protein_chain": self.protein_chain,
            "protein_atom": self.protein_atom,
            "distance": self.distance,
            "angle": self.angle,
            "halogen_type": self.halogen_type,
        }


@dataclass
class InteractionSummary:
    """Standardized interaction data structure for all workflows."""
    hydrogen_bonds: List[HydrogenBond] = field(default_factory=list)
    hydrophobic_contacts: List[HydrophobicContact] = field(default_factory=list)
    pi_stacking: List[PiStacking] = field(default_factory=list)
    salt_bridges: List[SaltBridge] = field(default_factory=list)
    halogen_bonds: List[HalogenBond] = field(default_factory=list)
    total_contacts: int = 0
    key_residues: List[str] = field(default_factory=list)
    visualization_data: Optional[Dict[str, Any]] = None
    analysis_method: str = "unknown"  # "plip" or "distance_based"
    status: str = "completed"
    error: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "hydrogen_bonds": [hb.to_dict() for hb in self.hydrogen_bonds],
            "hydrophobic_contacts": [hc.to_dict() for hc in self.hydrophobic_contacts],
            "pi_stacking": [ps.to_dict() for ps in self.pi_stacking],
            "salt_bridges": [sb.to_dict() for sb in self.salt_bridges],
            "halogen_bonds": [hb.to_dict() for hb in self.halogen_bonds],
            "total_contacts": self.total_contacts,
            "key_residues": self.key_residues,
            "visualization_data": self.visualization_data,
            "analysis_method": self.analysis_method,
            "status": self.status,
            "error": self.error,
        }


def analyze_interactions_plip(
    pdb_content: str,
    ligand_name: str = "UNL",
) -> InteractionSummary:
    """
    Comprehensive interaction analysis using PLIP.

    PLIP detects 8 interaction types:
    - Hydrogen bonds (with angle validation)
    - Hydrophobic contacts
    - Pi-stacking (face-to-face and edge-to-face)
    - Salt bridges
    - Pi-cation interactions
    - Halogen bonds
    - Water bridges
    - Metal coordination

    Args:
        pdb_content: PDB file content as string
        ligand_name: Ligand residue name (default "UNL")

    Returns:
        InteractionSummary with all detected interactions
    """
    try:
        from plip.structure.preparation import PDBComplex
        from plip.exchange.report import BindingSiteReport
    except ImportError:
        raise ImportError("PLIP not installed. Install with: pip install plip")

    summary = InteractionSummary(analysis_method="plip")

    # Write PDB to temp file (PLIP requires file path)
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        f.write(pdb_content)
        tmp_path = f.name

    try:
        # Run PLIP analysis
        mol = PDBComplex()
        mol.load_pdb(tmp_path)
        mol.analyze()

        # Get binding site for our ligand
        binding_sites = mol.interaction_sets

        for key, site in binding_sites.items():
            # Key format is like "UNL:A:1" (name:chain:resnum)
            if ligand_name in key:
                # Extract hydrogen bonds
                for hbond in site.hbonds_ldon + site.hbonds_pdon:
                    donor_coords = (hbond.d.x, hbond.d.y, hbond.d.z) if hasattr(hbond, 'd') else None
                    acceptor_coords = (hbond.a.x, hbond.a.y, hbond.a.z) if hasattr(hbond, 'a') else None

                    hb = HydrogenBond(
                        ligand_atom=getattr(hbond, 'atype', 'UNK'),
                        protein_residue=f"{hbond.restype}{hbond.resnr}",
                        protein_chain=hbond.reschain,
                        protein_atom=getattr(hbond, 'dtype', 'UNK'),
                        distance=round(hbond.dist_h_a if hasattr(hbond, 'dist_h_a') else hbond.dist_d_a, 2),
                        angle=round(hbond.angle, 1) if hasattr(hbond, 'angle') else None,
                        donor_coords=donor_coords,
                        acceptor_coords=acceptor_coords,
                        type="ligand_donor" if hbond in site.hbonds_ldon else "ligand_acceptor",
                    )
                    summary.hydrogen_bonds.append(hb)

                # Extract hydrophobic contacts
                for hydrophobic in site.hydrophobic_contacts:
                    hc = HydrophobicContact(
                        ligand_atom=getattr(hydrophobic, 'ligatom_orig_name', 'UNK'),
                        protein_residue=f"{hydrophobic.restype}{hydrophobic.resnr}",
                        protein_chain=hydrophobic.reschain,
                        protein_atom=getattr(hydrophobic, 'bsatom_orig_name', 'UNK'),
                        distance=round(hydrophobic.distance, 2),
                        coords=(hydrophobic.ligatom.x, hydrophobic.ligatom.y, hydrophobic.ligatom.z)
                               if hasattr(hydrophobic, 'ligatom') else None,
                    )
                    summary.hydrophobic_contacts.append(hc)

                # Extract pi-stacking
                for pistacking in site.pistacking:
                    ps = PiStacking(
                        ligand_ring=str(getattr(pistacking, 'ligandring', 'ring')),
                        protein_residue=f"{pistacking.restype}{pistacking.resnr}",
                        protein_chain=pistacking.reschain,
                        type="face-to-face" if pistacking.type == 'P' else "edge-to-face",
                        distance=round(pistacking.distance, 2),
                        angle=round(pistacking.angle, 1) if hasattr(pistacking, 'angle') else None,
                    )
                    summary.pi_stacking.append(ps)

                # Extract salt bridges
                for saltbridge in site.saltbridge_lneg + site.saltbridge_pneg:
                    sb = SaltBridge(
                        ligand_group=str(getattr(saltbridge, 'lig_group', 'charged')),
                        protein_residue=f"{saltbridge.restype}{saltbridge.resnr}",
                        protein_chain=saltbridge.reschain,
                        distance=round(saltbridge.distance, 2),
                        ligand_charge="negative" if saltbridge in site.saltbridge_lneg else "positive",
                        protein_charge="positive" if saltbridge in site.saltbridge_lneg else "negative",
                    )
                    summary.salt_bridges.append(sb)

                # Extract halogen bonds
                for halogenbond in site.halogen_bonds:
                    hb = HalogenBond(
                        ligand_atom=getattr(halogenbond, 'don_atom_orig_name', 'X'),
                        protein_residue=f"{halogenbond.restype}{halogenbond.resnr}",
                        protein_chain=halogenbond.reschain,
                        protein_atom=getattr(halogenbond, 'acc_atom_orig_name', 'O'),
                        distance=round(halogenbond.distance, 2),
                        angle=round(halogenbond.don_angle, 1) if hasattr(halogenbond, 'don_angle') else None,
                        halogen_type=getattr(halogenbond, 'halogen_type', 'X'),
                    )
                    summary.halogen_bonds.append(hb)

        # Calculate total contacts
        summary.total_contacts = (
            len(summary.hydrogen_bonds) +
            len(summary.hydrophobic_contacts) +
            len(summary.pi_stacking) +
            len(summary.salt_bridges) +
            len(summary.halogen_bonds)
        )

        # Extract key residues (most contacts)
        residue_counts = {}
        for hb in summary.hydrogen_bonds:
            res = f"{hb.protein_chain}:{hb.protein_residue}"
            residue_counts[res] = residue_counts.get(res, 0) + 2  # H-bonds weighted
        for hc in summary.hydrophobic_contacts:
            res = f"{hc.protein_chain}:{hc.protein_residue}"
            residue_counts[res] = residue_counts.get(res, 0) + 1
        for ps in summary.pi_stacking:
            res = f"{ps.protein_chain}:{ps.protein_residue}"
            residue_counts[res] = residue_counts.get(res, 0) + 3  # Pi-stacking weighted
        for sb in summary.salt_bridges:
            res = f"{sb.protein_chain}:{sb.protein_residue}"
            residue_counts[res] = residue_counts.get(res, 0) + 2

        # Top 10 key residues
        sorted_residues = sorted(residue_counts.items(), key=lambda x: x[1], reverse=True)
        summary.key_residues = [res for res, _ in sorted_residues[:10]]

    except Exception as e:
        summary.status = "error"
        summary.error = str(e)
    finally:
        # Clean up temp file
        try:
            os.unlink(tmp_path)
        except:
            pass

    return summary


def analyze_interactions_prolif(
    pdb_content: str,
    ligand_name: str = "UNL",
) -> InteractionSummary:
    """
    Protein-ligand interaction analysis using ProLIF.

    ProLIF is a modern, pip-installable alternative to PLIP that provides:
    - H-bond detection (donor/acceptor)
    - Hydrophobic contacts
    - Pi-stacking (face-to-face and edge-to-face)
    - Salt bridges
    - Cation-pi interactions
    - Halogen bonds

    Args:
        pdb_content: PDB file content as string
        ligand_name: Ligand residue name (default "UNL")

    Returns:
        InteractionSummary with all detected interactions
    """
    try:
        import prolif as plf
    except ImportError:
        raise ImportError("ProLIF not installed. Install with: pip install prolif")

    summary = InteractionSummary(analysis_method="prolif")

    # Write PDB to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        f.write(pdb_content)
        tmp_path = f.name

    try:
        # Use MDAnalysis to load PDB (more robust than RDKit for metal complexes)
        import MDAnalysis as mda

        u = mda.Universe(tmp_path)
        protein_sel = u.select_atoms("protein")
        ligand_sel = u.select_atoms(f"resname {ligand_name}")

        if len(ligand_sel) == 0:
            raise ValueError(f"Ligand '{ligand_name}' not found in PDB")

        if len(protein_sel) == 0:
            raise ValueError("No protein atoms found in PDB")

        # Convert to ProLIF molecules
        # NoImplicit=False allows structures without explicit hydrogens (like RFD3 output)
        protein_mol = plf.Molecule.from_mda(protein_sel, NoImplicit=False)
        ligand_mol = plf.Molecule.from_mda(ligand_sel, NoImplicit=False)

        # Run fingerprint analysis
        fp = plf.Fingerprint()
        fp.run_from_iterable([ligand_mol], protein_mol)

        # Extract interactions from fingerprint
        df = fp.to_dataframe()

        # Parse interactions
        for col in df.columns:
            if df[col].any():
                # Column format: (LigandResname, ProteinResname, InteractionType)
                if len(col) >= 3:
                    lig_res, prot_res, interaction_type = col[0], col[1], col[2]

                    # Map to our data structures
                    if 'HBDonor' in interaction_type or 'HBAcceptor' in interaction_type:
                        hb = HydrogenBond(
                            ligand_atom=str(lig_res),
                            protein_residue=str(prot_res),
                            protein_chain=str(prot_res).split('.')[0] if '.' in str(prot_res) else 'A',
                            protein_atom='',
                            distance=0.0,  # ProLIF doesn't provide distances directly
                            type="ligand_donor" if 'Donor' in interaction_type else "ligand_acceptor",
                        )
                        summary.hydrogen_bonds.append(hb)

                    elif 'Hydrophobic' in interaction_type:
                        hc = HydrophobicContact(
                            ligand_atom=str(lig_res),
                            protein_residue=str(prot_res),
                            protein_chain=str(prot_res).split('.')[0] if '.' in str(prot_res) else 'A',
                            protein_atom='',
                            distance=0.0,
                        )
                        summary.hydrophobic_contacts.append(hc)

                    elif 'PiStack' in interaction_type or 'EdgeToFace' in interaction_type:
                        ps = PiStacking(
                            ligand_ring=str(lig_res),
                            protein_residue=str(prot_res),
                            protein_chain=str(prot_res).split('.')[0] if '.' in str(prot_res) else 'A',
                            type="edge-to-face" if 'Edge' in interaction_type else "face-to-face",
                            distance=0.0,
                        )
                        summary.pi_stacking.append(ps)

                    elif 'Ionic' in interaction_type or 'Salt' in interaction_type:
                        sb = SaltBridge(
                            ligand_group=str(lig_res),
                            protein_residue=str(prot_res),
                            protein_chain=str(prot_res).split('.')[0] if '.' in str(prot_res) else 'A',
                            distance=0.0,
                            ligand_charge="charged",
                            protein_charge="charged",
                        )
                        summary.salt_bridges.append(sb)

                    elif 'Halogen' in interaction_type or 'XBond' in interaction_type:
                        xb = HalogenBond(
                            ligand_atom=str(lig_res),
                            protein_residue=str(prot_res),
                            protein_chain=str(prot_res).split('.')[0] if '.' in str(prot_res) else 'A',
                            protein_atom='',
                            distance=0.0,
                        )
                        summary.halogen_bonds.append(xb)

        # Calculate totals
        summary.total_contacts = (
            len(summary.hydrogen_bonds) +
            len(summary.hydrophobic_contacts) +
            len(summary.pi_stacking) +
            len(summary.salt_bridges) +
            len(summary.halogen_bonds)
        )

        # Key residues
        residue_counts = {}
        for hb in summary.hydrogen_bonds:
            res = f"{hb.protein_chain}:{hb.protein_residue}"
            residue_counts[res] = residue_counts.get(res, 0) + 2
        for hc in summary.hydrophobic_contacts:
            res = f"{hc.protein_chain}:{hc.protein_residue}"
            residue_counts[res] = residue_counts.get(res, 0) + 1
        for ps in summary.pi_stacking:
            res = f"{ps.protein_chain}:{ps.protein_residue}"
            residue_counts[res] = residue_counts.get(res, 0) + 3

        sorted_residues = sorted(residue_counts.items(), key=lambda x: x[1], reverse=True)
        summary.key_residues = [res for res, _ in sorted_residues[:10]]

    except Exception as e:
        summary.status = "error"
        summary.error = str(e)
    finally:
        try:
            os.unlink(tmp_path)
        except:
            pass

    return summary


def _calculate_angle(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    """Calculate angle in degrees at point p2 (D-H-A angle)."""
    v1 = p1 - p2
    v2 = p3 - p2

    cosine = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    cosine = np.clip(cosine, -1.0, 1.0)  # Avoid numerical errors

    return math.degrees(math.acos(cosine))


def _analyze_hbonds_distance_based(
    pdb_content: str,
    ligand_name: str = "UNL",
    hbond_distance: float = 3.5,
    min_angle: float = 120.0,
) -> List[HydrogenBond]:
    """
    Distance and angle-based H-bond detection as fallback.

    Criteria:
    - D-A distance: <= 3.5 Angstrom
    - D-H-A angle: >= 120 degrees (when H position available)
    """
    from biotite.structure.io.pdb import PDBFile
    import io

    pdb_file = PDBFile.read(io.StringIO(pdb_content))
    structure = pdb_file.get_structure(model=1)

    # Get ligand and protein atoms
    ligand_mask = structure.res_name == ligand_name
    ligand_struct = structure[ligand_mask]

    protein_mask = ~structure.hetero
    protein_struct = structure[protein_mask]

    if len(ligand_struct) == 0 or len(protein_struct) == 0:
        return []

    # Donor/acceptor definitions
    protein_donor_atoms = ['N', 'NE', 'NH1', 'NH2', 'ND1', 'ND2', 'NE1', 'NE2', 'NZ', 'OG', 'OG1', 'OH', 'SG']
    protein_acceptor_atoms = ['O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH', 'ND1', 'NE2', 'SD']
    ligand_acceptor_types = ['N', 'O', 'S']
    ligand_donor_types = ['N', 'O']

    # Build atom lookup for angle calculation
    protein_atoms_by_res = {}
    for atom in protein_struct:
        key = (atom.chain_id, atom.res_id)
        if key not in protein_atoms_by_res:
            protein_atoms_by_res[key] = {}
        protein_atoms_by_res[key][atom.atom_name] = atom.coord

    hbonds = []

    for i in range(len(ligand_struct)):
        lig_atom = ligand_struct[i]
        lig_atom_name = lig_atom.atom_name
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

            is_hbond = False
            hbond_type = None

            if is_lig_acceptor and is_prot_donor:
                is_hbond = True
                hbond_type = "ligand_acceptor"
            elif is_lig_donor and is_prot_acceptor:
                is_hbond = True
                hbond_type = "ligand_donor"

            if is_hbond:
                # Try to calculate angle
                angle = None
                res_key = (prot_atom.chain_id, prot_atom.res_id)

                # For backbone N as donor, H is approximately along N-CA vector
                if hbond_type == "ligand_acceptor" and prot_atom_name == 'N':
                    res_atoms = protein_atoms_by_res.get(res_key, {})
                    if 'CA' in res_atoms:
                        # Estimate H position (along N-CA direction)
                        n_coord = prot_atom.coord
                        ca_coord = res_atoms['CA']
                        h_direction = n_coord - ca_coord
                        h_direction = h_direction / np.linalg.norm(h_direction)
                        h_coord = n_coord + h_direction * 1.0  # ~1 Angstrom

                        # D-H-A angle (N-H-Ligand)
                        angle = _calculate_angle(n_coord, h_coord, lig_atom.coord)

                # Skip if angle is too small (not a real H-bond)
                if angle is not None and angle < min_angle:
                    continue

                hb = HydrogenBond(
                    ligand_atom=lig_atom_name,
                    protein_residue=f"{prot_atom.res_name}{prot_atom.res_id}",
                    protein_chain=prot_atom.chain_id,
                    protein_atom=prot_atom_name,
                    distance=round(float(dist), 2),
                    angle=round(angle, 1) if angle else None,
                    donor_coords=tuple(prot_atom.coord) if hbond_type == "ligand_acceptor" else tuple(lig_atom.coord),
                    acceptor_coords=tuple(lig_atom.coord) if hbond_type == "ligand_acceptor" else tuple(prot_atom.coord),
                    type=hbond_type,
                )
                hbonds.append(hb)

    return hbonds


def _analyze_hydrophobic_distance_based(
    pdb_content: str,
    ligand_name: str = "UNL",
    contact_distance: float = 4.0,
) -> List[HydrophobicContact]:
    """Distance-based hydrophobic contact detection."""
    from biotite.structure.io.pdb import PDBFile
    import io

    pdb_file = PDBFile.read(io.StringIO(pdb_content))
    structure = pdb_file.get_structure(model=1)

    ligand_mask = structure.res_name == ligand_name
    ligand_struct = structure[ligand_mask]

    protein_mask = ~structure.hetero
    protein_struct = structure[protein_mask]

    if len(ligand_struct) == 0 or len(protein_struct) == 0:
        return []

    # Hydrophobic residues and their carbon atoms
    hydrophobic_residues = ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO']
    carbon_atoms = ['C', 'CA', 'CB', 'CG', 'CG1', 'CG2', 'CD', 'CD1', 'CD2', 'CE', 'CE1', 'CE2', 'CZ']

    contacts = []

    for i in range(len(ligand_struct)):
        lig_atom = ligand_struct[i]
        element = lig_atom.element if hasattr(lig_atom, 'element') else lig_atom.atom_name[0]

        # Skip non-carbon ligand atoms
        if element != 'C':
            continue

        for j in range(len(protein_struct)):
            prot_atom = protein_struct[j]

            # Check if hydrophobic residue and carbon atom
            if prot_atom.res_name not in hydrophobic_residues:
                continue
            if prot_atom.atom_name not in carbon_atoms:
                continue

            dist = np.linalg.norm(lig_atom.coord - prot_atom.coord)

            if dist <= contact_distance:
                hc = HydrophobicContact(
                    ligand_atom=lig_atom.atom_name,
                    protein_residue=f"{prot_atom.res_name}{prot_atom.res_id}",
                    protein_chain=prot_atom.chain_id,
                    protein_atom=prot_atom.atom_name,
                    distance=round(float(dist), 2),
                    coords=tuple(lig_atom.coord),
                )
                contacts.append(hc)

    return contacts


def analyze_all_interactions(
    pdb_content: str,
    ligand_name: str = "UNL",
    include_visualization_data: bool = True,
    use_plip: bool = True,
    use_prolif: bool = True,
) -> InteractionSummary:
    """
    Unified interaction analysis used by all workflows.

    Priority order:
    1. ProLIF - Modern, pip-installable (requires rdkit, MDAnalysis)
    2. PLIP - Comprehensive (requires OpenBabel)
    3. Biotite distance-based - Fallback

    Args:
        pdb_content: PDB file content as string
        ligand_name: Ligand residue name (default "UNL")
        include_visualization_data: Include data for 3D visualization
        use_plip: Try to use PLIP (default True)
        use_prolif: Try to use ProLIF (default True)

    Returns:
        InteractionSummary with all detected interactions
    """
    # Check if PDB has explicit hydrogens (needed for ProLIF/PLIP to work properly)
    has_hydrogens = ' H ' in pdb_content or ' H\n' in pdb_content

    # For structures with explicit hydrogens, try advanced tools
    if has_hydrogens:
        # Try ProLIF first (modern, pip-installable)
        if use_prolif:
            try:
                summary = analyze_interactions_prolif(pdb_content, ligand_name)
                if summary.status == "completed" and summary.total_contacts > 0:
                    if include_visualization_data:
                        summary.visualization_data = _generate_visualization_data(summary)
                    return summary
            except ImportError:
                pass  # ProLIF not installed, fall back
            except Exception as e:
                pass  # ProLIF failed, fall back

        # Try PLIP second (comprehensive but requires OpenBabel)
        if use_plip:
            try:
                summary = analyze_interactions_plip(pdb_content, ligand_name)
                if summary.status == "completed" and summary.total_contacts > 0:
                    if include_visualization_data:
                        summary.visualization_data = _generate_visualization_data(summary)
                    return summary
            except ImportError:
                pass  # PLIP not installed, fall back
            except Exception as e:
                pass  # PLIP failed, fall back

    # Use biotite distance-based analysis
    # This works well for RFD3 outputs which lack explicit hydrogens
    summary = InteractionSummary(analysis_method="distance_based")

    try:
        # H-bonds with angle validation
        summary.hydrogen_bonds = _analyze_hbonds_distance_based(
            pdb_content, ligand_name, hbond_distance=3.5, min_angle=120.0
        )

        # Hydrophobic contacts
        summary.hydrophobic_contacts = _analyze_hydrophobic_distance_based(
            pdb_content, ligand_name, contact_distance=4.0
        )

        # Calculate totals
        summary.total_contacts = (
            len(summary.hydrogen_bonds) +
            len(summary.hydrophobic_contacts)
        )

        # Extract key residues
        residue_counts = {}
        for hb in summary.hydrogen_bonds:
            res = f"{hb.protein_chain}:{hb.protein_residue}"
            residue_counts[res] = residue_counts.get(res, 0) + 2
        for hc in summary.hydrophobic_contacts:
            res = f"{hc.protein_chain}:{hc.protein_residue}"
            residue_counts[res] = residue_counts.get(res, 0) + 1

        sorted_residues = sorted(residue_counts.items(), key=lambda x: x[1], reverse=True)
        summary.key_residues = [res for res, _ in sorted_residues[:10]]

        if include_visualization_data:
            summary.visualization_data = _generate_visualization_data(summary)

        summary.status = "completed"

    except Exception as e:
        summary.status = "error"
        summary.error = str(e)

    return summary


def _generate_visualization_data(summary: InteractionSummary) -> Dict[str, Any]:
    """Generate data for 3D visualization (Molstar/py3Dmol compatible)."""
    viz_data = {
        "lines": [],      # For H-bonds (dashed lines)
        "spheres": [],    # For hydrophobic contacts
        "planes": [],     # For pi-stacking
        "highlights": [], # For key residues
    }

    # H-bond lines (blue dashed)
    for hb in summary.hydrogen_bonds:
        if hb.donor_coords and hb.acceptor_coords:
            viz_data["lines"].append({
                "start": hb.donor_coords,
                "end": hb.acceptor_coords,
                "color": "#3B82F6",  # Blue
                "dashed": True,
                "label": f"{hb.distance} A",
            })

    # Hydrophobic contact spheres (yellow)
    for hc in summary.hydrophobic_contacts:
        if hc.coords:
            viz_data["spheres"].append({
                "center": hc.coords,
                "radius": 0.5,
                "color": "#FCD34D",  # Yellow
                "opacity": 0.6,
            })

    # Pi-stacking visualization
    for ps in summary.pi_stacking:
        if ps.ligand_center and ps.protein_center:
            viz_data["planes"].append({
                "center1": ps.ligand_center,
                "center2": ps.protein_center,
                "color": "#F97316",  # Orange
                "type": ps.type,
            })

    # Key residue highlights
    viz_data["highlights"] = summary.key_residues

    return viz_data


def format_for_frontend(summary: InteractionSummary) -> Dict[str, Any]:
    """
    Format interaction data for frontend display.

    Returns a structure suitable for React components.
    """
    return {
        "interactions": {
            "hbonds": len(summary.hydrogen_bonds),
            "hydrophobic": len(summary.hydrophobic_contacts),
            "pi_stacking": len(summary.pi_stacking),
            "salt_bridges": len(summary.salt_bridges),
            "halogen_bonds": len(summary.halogen_bonds),
            "total": summary.total_contacts,
        },
        "key_residues": summary.key_residues,
        "details": {
            "hydrogen_bonds": [hb.to_dict() for hb in summary.hydrogen_bonds[:20]],
            "hydrophobic_contacts": [hc.to_dict() for hc in summary.hydrophobic_contacts[:30]],
            "pi_stacking": [ps.to_dict() for ps in summary.pi_stacking],
            "salt_bridges": [sb.to_dict() for sb in summary.salt_bridges],
        },
        "visualization": summary.visualization_data,
        "analysis_method": summary.analysis_method,
        "status": summary.status,
        "error": summary.error,
    }


def format_for_ai_assistant(summary: InteractionSummary) -> str:
    """
    Format interaction data as natural language for AI assistant.

    Returns human-readable summary text.
    """
    lines = []

    # Summary counts
    lines.append(f"Found {summary.total_contacts} total protein-ligand interactions:")
    lines.append(f"  - {len(summary.hydrogen_bonds)} hydrogen bonds")
    lines.append(f"  - {len(summary.hydrophobic_contacts)} hydrophobic contacts")

    if summary.pi_stacking:
        lines.append(f"  - {len(summary.pi_stacking)} pi-stacking interactions")
    if summary.salt_bridges:
        lines.append(f"  - {len(summary.salt_bridges)} salt bridges")
    if summary.halogen_bonds:
        lines.append(f"  - {len(summary.halogen_bonds)} halogen bonds")

    # Key binding residues
    if summary.key_residues:
        lines.append(f"\nKey binding residues: {', '.join(summary.key_residues[:5])}")

    # Notable H-bonds
    if summary.hydrogen_bonds:
        lines.append("\nNotable hydrogen bonds:")
        for hb in summary.hydrogen_bonds[:3]:
            angle_str = f", {hb.angle}deg" if hb.angle else ""
            lines.append(f"  - {hb.protein_residue}.{hb.protein_atom} - {hb.ligand_atom} ({hb.distance}A{angle_str})")

    # Analysis method note
    if summary.analysis_method == "distance_based":
        lines.append("\n(Analysis used distance-based method; PLIP unavailable for detailed profiling)")

    return "\n".join(lines)


def generate_recommendations(summary: InteractionSummary, ligand_has_aromatics: bool = False) -> List[str]:
    """
    Generate actionable recommendations based on interaction analysis.

    Args:
        summary: InteractionSummary from analysis
        ligand_has_aromatics: Whether ligand contains aromatic rings

    Returns:
        List of recommendation strings
    """
    recs = []

    # H-bond recommendations
    if len(summary.hydrogen_bonds) < 2:
        recs.append("Consider adding polar residues (Gln, Asn, Ser) near ligand polar groups to increase H-bonding")
    elif len(summary.hydrogen_bonds) >= 4:
        recs.append("Strong H-bond network detected - good for specificity")

    # Hydrophobic recommendations
    if len(summary.hydrophobic_contacts) < 5:
        recs.append("Binding pocket may benefit from more hydrophobic residues (Leu, Ile, Val, Phe)")
    elif len(summary.hydrophobic_contacts) >= 15:
        recs.append("Excellent hydrophobic burial - contributes to binding affinity")

    # Pi-stacking recommendations
    if ligand_has_aromatics and not summary.pi_stacking:
        recs.append("Ligand has aromatic rings - aromatic residues (Phe, Tyr, Trp) could form pi-stacking")

    # Salt bridge recommendations
    if not summary.salt_bridges:
        recs.append("No salt bridges detected - consider charged residues if ligand has ionizable groups")

    # Overall assessment
    if summary.total_contacts < 10:
        recs.append("Low total contacts - binding may be weak; consider redesigning binding pocket")
    elif summary.total_contacts >= 25:
        recs.append("High contact count suggests strong binding potential")

    return recs
