"""
Unified Design Analyzer

Orchestrates all analysis tools into a single entry point producing structured JSON.
Auto-detects design type and runs appropriate analyses.
"""
import os
import json
import shutil
from datetime import datetime
from typing import Dict, Any, Optional, List

from analysis_types import (
    AnalysisResult,
    AnalysisStatus,
    DesignType,
    detect_design_type,
    detect_metal_from_pdb,
    detect_ligand_from_pdb,
    DetectedMetal,
    DetectedLigand,
    METAL_CODES,
    SOLVENT_RESIDUES,
)

# Try to import PLIP-based interaction analysis
try:
    from shared.interaction_analysis import analyze_all_interactions, InteractionSummary
    PLIP_AVAILABLE = True
except ImportError:
    PLIP_AVAILABLE = False


def sanitize_pdb(pdb_content: str, remove_duplicates: bool = True) -> tuple:
    """
    Sanitize PDB content to fix common issues from RFD3 output.

    RFdiffusion3 generates PDBs with known format issues:
    1. Duplicate atom serial numbers (sidechain atoms share backbone number)
    2. Atom numbering starting at 0 instead of 1
    3. Missing TER records between chains
    4. Missing END record
    5. Duplicate HETATM records (ligand/metal appears twice)
    6. Negative B-factors (minor)

    Args:
        pdb_content: Raw PDB file content
        remove_duplicates: If True, remove duplicate HETATM records

    Returns:
        Tuple of (sanitized_pdb_content, list_of_issues_fixed)
    """
    import math

    # Guard: detect CIF content that was accidentally passed as PDB
    if pdb_content and pdb_content.strip().startswith("data_"):
        import logging as _log
        _log.getLogger(__name__).warning("sanitize_pdb received CIF content instead of PDB — skipping sanitization")
        return pdb_content, ["warning: content is CIF format, not PDB"]

    lines = pdb_content.split('\n')
    result_lines = []
    issues_fixed = []
    atom_serial = 1

    # Track seen HETATM records to detect duplicates
    # Key: (res_name, chain, resnum, atom_name, x, y, z) -> first occurrence line
    seen_hetatm = {}

    # Track chains to add TER records
    current_chain = None
    last_record_type = None  # 'ATOM' or 'HETATM'

    # Track atom names per residue for duplicate detection
    residue_atom_names: Dict[tuple, set] = {}

    for line in lines:
        # Skip empty lines
        if not line.strip():
            continue

        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                record_type = line[:6].strip()
                atom_name = line[12:16].strip()
                alt_loc = line[16] if len(line) > 16 else ' '
                res_name = line[17:20].strip()
                chain_id = line[21] if len(line) > 21 else 'A'
                res_num = int(line[22:26]) if len(line) > 25 else 1
                ins_code = line[26] if len(line) > 26 else ' '
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                # Skip atoms with NaN coordinates
                if math.isnan(x) or math.isnan(y) or math.isnan(z):
                    issues_fixed.append(f"removed NaN coord atom: {atom_name} {res_name}{res_num}")
                    continue

                occupancy = float(line[54:60]) if len(line) > 59 else 1.0
                b_factor = float(line[60:66]) if len(line) > 65 else 0.0

                # Fix negative B-factors
                if b_factor < 0:
                    b_factor = 0.0
                    if "fixed negative B-factor" not in issues_fixed:
                        issues_fixed.append("fixed negative B-factor(s)")

                element = line[76:78].strip() if len(line) > 77 else atom_name[0]

            except (ValueError, IndexError) as e:
                # Malformed line - skip
                issues_fixed.append(f"skipped malformed line: {line[:30]}...")
                continue

            # Check for duplicate HETATM records
            # RFD3 may output same ligand/metal twice with different chain IDs
            # Compare by (res_name, atom_name, coordinates) ignoring chain/resnum
            if record_type == 'HETATM' and remove_duplicates:
                # Round coordinates to detect duplicates with minor float differences
                hetatm_key = (res_name, atom_name,
                              round(x, 2), round(y, 2), round(z, 2))
                if hetatm_key in seen_hetatm:
                    orig_chain = seen_hetatm[hetatm_key]
                    if "removed duplicate HETATM" not in [i.split(':')[0] for i in issues_fixed]:
                        issues_fixed.append(f"removed duplicate HETATM: {res_name} (kept chain {orig_chain}, removed chain {chain_id})")
                    continue
                seen_hetatm[hetatm_key] = chain_id  # Track which chain we kept

            # Check for duplicate atom names within residue
            residue_key = (chain_id, res_num, ins_code)
            if residue_key not in residue_atom_names:
                residue_atom_names[residue_key] = set()

            original_atom_name = atom_name
            if atom_name in residue_atom_names[residue_key]:
                # Duplicate atom name - this is expected for RFD3 sidechain atoms
                # Don't rename, just track for reporting
                pass  # We'll still add it - biotite/analysis tools handle this
            residue_atom_names[residue_key].add(atom_name)

            # Add TER record when switching from ATOM to HETATM or between chains
            if record_type == 'ATOM':
                if last_record_type == 'ATOM' and current_chain and current_chain != chain_id:
                    result_lines.append(f"TER   {atom_serial:5d}")
                    atom_serial += 1
                    if "added TER records" not in issues_fixed:
                        issues_fixed.append("added TER records")
                current_chain = chain_id
            elif record_type == 'HETATM' and last_record_type == 'ATOM':
                result_lines.append(f"TER   {atom_serial:5d}")
                atom_serial += 1
                if "added TER records" not in issues_fixed:
                    issues_fixed.append("added TER records")

            last_record_type = record_type

            # Format atom name for PDB (4 characters, proper alignment)
            if len(atom_name) < 4:
                if len(element) == 1:
                    atom_name_fmt = f" {atom_name:<3}"
                else:
                    atom_name_fmt = f"{atom_name:<4}"
            else:
                atom_name_fmt = atom_name[:4]

            # Reconstruct line with sequential atom serial
            if record_type == 'HETATM':
                new_line = (
                    f"HETATM{atom_serial:5d} {atom_name_fmt}{alt_loc}{res_name:>3} "
                    f"{chain_id}{res_num:4d}{ins_code}   "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{b_factor:6.2f}          {element:>2}  "
                )
            else:
                new_line = (
                    f"ATOM  {atom_serial:5d} {atom_name_fmt}{alt_loc}{res_name:>3} "
                    f"{chain_id}{res_num:4d}{ins_code}   "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{b_factor:6.2f}          {element:>2}  "
                )

            result_lines.append(new_line)
            atom_serial += 1

        elif line.startswith('TER'):
            # Keep existing TER but renumber
            result_lines.append(f"TER   {atom_serial:5d}")
            atom_serial += 1

        elif line.startswith('END'):
            # Skip - we'll add at the end
            continue

        elif line.startswith('CONECT') or line.startswith('MASTER'):
            # Skip connectivity records
            continue

        else:
            # Keep other lines (REMARK, CRYST1, HEADER, etc.)
            result_lines.append(line)

    # Add END record if not present
    if not any(l.startswith('END') for l in result_lines):
        result_lines.append("END")
        if "added END record" not in issues_fixed:
            issues_fixed.append("added END record")

    # Track if we renumbered atoms
    if atom_serial > 1:
        issues_fixed.insert(0, f"renumbered {atom_serial - 1} atoms sequentially")

    return '\n'.join(result_lines), issues_fixed


def validate_metal_coordination(
    pdb_content: str,
    metal: str,
    expected_cn: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Validate metal coordination geometry in a predicted structure.

    Parses PDB for metal HETATM position, finds coordinating atoms within 3.0A,
    and calls metal_chemistry.validate_coordination() to check bond distances,
    coordination number, and geometry.

    Args:
        pdb_content: PDB file content (e.g., RF3 prediction output)
        metal: Metal element symbol (e.g., "TB", "ZN", "CA")
        expected_cn: Expected coordination number (auto-detected if None)

    Returns:
        Dict with keys:
            - valid: bool - whether coordination geometry is acceptable
            - coordination_number: int - actual CN found
            - distances: List[float] - metal-ligand distances
            - geometry_match: str or None - best matching ideal geometry
            - issues: List[str] - list of detected problems
    """
    import math

    metal_upper = metal.upper()

    # Find metal position
    metal_pos = None
    for line in pdb_content.split('\n'):
        if not line.startswith('HETATM'):
            continue
        try:
            res_name = line[17:20].strip()
            element = line[76:78].strip() if len(line) >= 78 else line[12:16].strip()
            if res_name.upper() == metal_upper or element.upper() == metal_upper:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                metal_pos = (x, y, z)
                break
        except (ValueError, IndexError):
            continue

    if metal_pos is None:
        return {
            "valid": False,
            "coordination_number": 0,
            "distances": [],
            "geometry_match": None,
            "issues": [f"Metal {metal_upper} not found in PDB"],
        }

    # Find coordinating atoms within 3.0A (O, N, S atoms)
    coord_positions = []
    coord_atoms = []
    mx, my, mz = metal_pos

    for line in pdb_content.split('\n'):
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            continue
        try:
            element = line[76:78].strip() if len(line) >= 78 else line[12:16].strip()[0]
            if element.upper() not in ('O', 'N', 'S'):
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            dist = math.sqrt((x - mx)**2 + (y - my)**2 + (z - mz)**2)
            if dist < 3.0 and dist > 0.5:  # Within coordination distance
                coord_positions.append((x, y, z))
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                res_num = line[22:26].strip()
                coord_atoms.append(f"{res_name}{res_num}:{atom_name}")
        except (ValueError, IndexError):
            continue

    # Call metal_chemistry.validate_coordination if available
    try:
        from metal_chemistry import validate_coordination
        result = validate_coordination(coord_positions, metal_pos, metal_upper, expected_cn)
        result["coordinating_atoms"] = coord_atoms
        return result
    except ImportError:
        # Fallback: basic validation
        actual_cn = len(coord_positions)
        distances = [
            math.sqrt((x - mx)**2 + (y - my)**2 + (z - mz)**2)
            for x, y, z in coord_positions
        ]
        issues = []
        if expected_cn and actual_cn < expected_cn:
            issues.append(f"Low CN: found {actual_cn}, expected {expected_cn}")

        return {
            "valid": len(issues) == 0 and actual_cn >= 3,
            "coordination_number": actual_cn,
            "distances": distances,
            "geometry_match": None,
            "issues": issues,
            "coordinating_atoms": coord_atoms,
        }


class UnifiedDesignAnalyzer:
    """
    Unified analyzer that orchestrates all design analysis tools.

    Usage:
        analyzer = UnifiedDesignAnalyzer()
        result = analyzer.analyze(pdb_content, design_params, metal_type="TB")
    """

    def __init__(self):
        """Initialize analyzer with available analysis modules."""
        self._analysis_modules = self._discover_modules()

    def _discover_modules(self) -> Dict[str, bool]:
        """Discover which analysis modules are available."""
        modules = {
            "binding_analysis": False,
            "metal_validation": False,
            "metal_chemistry": False,
            "tebl_analysis": False,
            "hotspot_detection": False,
            "topology_validation": False,
            "gnina": False,
            "interaction_analysis": PLIP_AVAILABLE,
            "pyrosetta": False,
            "conservation": False,
        }

        # Check for Python modules
        try:
            import binding_analysis
            modules["binding_analysis"] = True
        except ImportError:
            pass

        try:
            import metal_validation
            modules["metal_validation"] = True
        except ImportError:
            pass

        try:
            import metal_chemistry
            modules["metal_chemistry"] = True
        except ImportError:
            pass

        try:
            import tebl_analysis
            modules["tebl_analysis"] = True
        except ImportError:
            pass

        try:
            import hotspot_detection
            modules["hotspot_detection"] = True
        except ImportError:
            pass

        try:
            import topology_validation
            modules["topology_validation"] = True
        except ImportError:
            pass

        # Check for GNINA binary (cross-platform)
        modules["gnina"] = shutil.which("gnina") is not None

        # Check for PyRosetta (interface scoring, FastRelax)
        try:
            import rosetta_utils
            modules["pyrosetta"] = rosetta_utils.check_pyrosetta_available()
        except ImportError:
            pass

        # Check for conservation analyzer (ConSurf methodology)
        try:
            from conservation_analyzer import check_conservation_available
            availability = check_conservation_available()
            modules["conservation"] = availability.get("basic_pipeline", False)
        except ImportError:
            pass

        return modules

    def _detect_design_type(
        self,
        pdb_content: str,
        metal_type: Optional[str] = None,
        ligand_sdf: Optional[str] = None,
    ) -> tuple:
        """
        Detect design type from PDB content with auto-detection of metal/ligand.

        Returns:
            Tuple of (DesignType, DetectedMetal or None, DetectedLigand or None)
        """
        # Count protein chains
        chains = set()
        for line in pdb_content.split("\n"):
            if line.startswith("ATOM"):
                chain_id = line[21] if len(line) > 21 else ""
                if chain_id.strip():
                    chains.add(chain_id)

        chain_count = len(chains)

        # Auto-detect metal from PDB
        detected_metal = detect_metal_from_pdb(pdb_content)

        # Auto-detect ligand from PDB
        detected_ligand = detect_ligand_from_pdb(pdb_content)

        # Determine presence based on explicit params OR auto-detection
        has_metal = metal_type is not None or detected_metal is not None
        has_ligand = ligand_sdf is not None or detected_ligand is not None

        design_type = detect_design_type(
            has_ligand=has_ligand,
            has_metal=has_metal,
            chain_count=chain_count,
        )

        return design_type, detected_metal, detected_ligand

    def _generate_design_id(self, design_type: DesignType) -> str:
        """Generate unique design ID."""
        timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
        type_short = design_type.value.replace("_", "-")
        return f"{timestamp}_{type_short}"

    def analyze(
        self,
        pdb_content: str,
        design_params: Dict[str, Any],
        pdb_path: Optional[str] = None,
        ligand_sdf: Optional[str] = None,
        metal_type: Optional[str] = None,
        metal_chain: Optional[str] = None,
        metal_resnum: Optional[int] = None,
    ) -> Dict[str, Any]:
        """
        Run comprehensive analysis on a design.

        Args:
            pdb_content: PDB file content as string
            design_params: Parameters used to generate the design
            pdb_path: Optional path to PDB file (for tools that need file path)
            ligand_sdf: Optional SDF content for ligand
            metal_type: Optional metal type code (e.g., "TB", "ZN"). Auto-detected if None.
            metal_chain: Chain ID of metal. Auto-detected if None.
            metal_resnum: Residue number of metal. Auto-detected if None.

        Returns:
            Structured metrics JSON with all analysis results
        """
        # Sanitize PDB content to fix RFD3 output issues
        sanitized_pdb, sanitization_issues = sanitize_pdb(pdb_content)

        # Use sanitized content for all analysis
        pdb_content = sanitized_pdb

        # Detect design type and auto-detect metal/ligand
        design_type, detected_metal, detected_ligand = self._detect_design_type(
            pdb_content, metal_type, ligand_sdf
        )
        design_id = self._generate_design_id(design_type)

        # Use auto-detected values if explicit values not provided
        effective_metal_type = metal_type
        effective_metal_chain = metal_chain
        effective_metal_resnum = metal_resnum
        metal_auto_detected = False

        if detected_metal:
            if effective_metal_type is None:
                effective_metal_type = detected_metal.metal_type
                metal_auto_detected = True
            if effective_metal_chain is None:
                effective_metal_chain = detected_metal.chain
                metal_auto_detected = True
            if effective_metal_resnum is None:
                effective_metal_resnum = detected_metal.resnum
                metal_auto_detected = True

        # Initialize result structure
        result = {
            "design_id": design_id,
            "design_type": design_type.value,
            "timestamp": datetime.now().isoformat(),
            "design_params": design_params,
            "analyses": {},
        }

        # Add sanitization info if any issues were fixed
        if sanitization_issues:
            result["pdb_sanitization"] = {
                "issues_fixed": sanitization_issues,
                "count": len(sanitization_issues),
            }

        # Add auto-detected information
        auto_detected = {}
        if detected_metal:
            auto_detected["metal"] = {
                "type": detected_metal.metal_type,
                "chain": detected_metal.chain,
                "resnum": detected_metal.resnum,
                "coords": [detected_metal.x, detected_metal.y, detected_metal.z],
            }
        if detected_ligand:
            auto_detected["ligand"] = {
                "name": detected_ligand.ligand_name,
                "chain": detected_ligand.chain,
                "resnum": detected_ligand.resnum,
                "atom_count": detected_ligand.atom_count,
            }
        if auto_detected:
            result["auto_detected"] = auto_detected

        # Run structure confidence analysis
        result["analyses"]["structure_confidence"] = self._analyze_structure_confidence(
            pdb_content
        ).to_dict()

        # Run interface analysis (for dimers only)
        monomer_types = {
            DesignType.MONOMER, DesignType.METAL_MONOMER,
            DesignType.LIGAND_MONOMER, DesignType.METAL_LIGAND_MONOMER
        }
        if design_type not in monomer_types:
            result["analyses"]["interface_quality"] = self._analyze_interface(
                pdb_content
            ).to_dict()
        else:
            result["analyses"]["interface_quality"] = AnalysisResult.not_applicable(
                "monomer design - no interface"
            ).to_dict()

        # Run metal coordination analysis (uses effective values from auto-detection)
        # Pass ligand info for metal-ligand complex analysis
        effective_ligand_name = detected_ligand.ligand_name if detected_ligand else None
        if effective_metal_type:
            result["analyses"]["metal_coordination"] = self._analyze_metal_coordination(
                pdb_content, effective_metal_type, effective_metal_chain, effective_metal_resnum,
                ligand_name=effective_ligand_name
            ).to_dict()
        else:
            result["analyses"]["metal_coordination"] = AnalysisResult.not_applicable(
                "no metal detected or specified"
            ).to_dict()

        # Run ligand binding analysis (GNINA scoring with SDF)
        if ligand_sdf:
            result["analyses"]["ligand_binding"] = self._analyze_ligand_binding(
                pdb_content, pdb_path, ligand_sdf
            ).to_dict()
        else:
            result["analyses"]["ligand_binding"] = AnalysisResult.not_applicable(
                "no ligand SDF provided"
            ).to_dict()

        # Run PLIP-based ligand interaction analysis (works with PDB-embedded ligands)
        if detected_ligand:
            result["analyses"]["ligand_interactions"] = self._analyze_ligand_interactions(
                pdb_content, detected_ligand.ligand_name
            ).to_dict()
        else:
            result["analyses"]["ligand_interactions"] = AnalysisResult.not_applicable(
                "no ligand detected in PDB"
            ).to_dict()

        # Run sequence composition analysis (Ala check, aromatics, etc.)
        result["analyses"]["sequence_composition"] = self._analyze_sequence_composition(
            pdb_content
        ).to_dict()

        # Run topology validation
        result["analyses"]["topology"] = self._analyze_topology(pdb_content).to_dict()

        # Run symmetry analysis (for dimers only)
        if design_type not in monomer_types:
            result["analyses"]["symmetry"] = self._analyze_symmetry(pdb_content).to_dict()
        else:
            result["analyses"]["symmetry"] = AnalysisResult.not_applicable(
                "monomer design - no symmetry analysis"
            ).to_dict()

        # Run PyRosetta interface scoring (BindCraft-style metrics for dimers)
        if design_type not in monomer_types:
            # Get chain IDs from protein atoms
            protein_chains = []
            for line in pdb_content.split("\n"):
                if line.startswith("ATOM"):
                    chain_id = line[21] if len(line) > 21 else ""
                    if chain_id.strip() and chain_id not in protein_chains:
                        protein_chains.append(chain_id)

            if len(protein_chains) >= 2:
                result["analyses"]["pyrosetta_interface"] = self._analyze_pyrosetta_interface(
                    pdb_content, chain_a=protein_chains[0], chain_b=protein_chains[1]
                ).to_dict()
            else:
                result["analyses"]["pyrosetta_interface"] = AnalysisResult.not_applicable(
                    "need at least 2 protein chains for interface analysis"
                ).to_dict()
        else:
            result["analyses"]["pyrosetta_interface"] = AnalysisResult.not_applicable(
                "monomer design - no interface analysis"
            ).to_dict()

        # Run enzyme activity preservation analysis (Phase 6: AI Infrastructure)
        # Triggered when enzyme_class is specified in design_params
        enzyme_class = design_params.get("enzyme_class")
        if enzyme_class:
            result["analyses"]["enzyme_preservation"] = self._analyze_enzyme_preservation(
                designed_pdb=pdb_content,
                enzyme_class=enzyme_class,
                reference_pdb=design_params.get("reference_pdb"),
                catalytic_residues=design_params.get("catalytic_residues", []),
                exposed_atoms=design_params.get("select_exposed", {}),
            ).to_dict()
        else:
            result["analyses"]["enzyme_preservation"] = AnalysisResult.not_applicable(
                "no enzyme class specified"
            ).to_dict()

        # Run evolutionary conservation analysis (ConSurf methodology)
        # Note: This is an async operation that may take time for BLAST search
        # Skip by default in quick analysis, enable with design_params["analyze_conservation"]
        if design_params.get("analyze_conservation", False):
            result["analyses"]["conservation"] = self._analyze_conservation(
                pdb_content,
                chain=design_params.get("conservation_chain", "A"),
            ).to_dict()
        else:
            result["analyses"]["conservation"] = AnalysisResult.skipped(
                "conservation analysis disabled (set analyze_conservation=True to enable)"
            ).to_dict()

        return result

    def _analyze_structure_confidence(self, pdb_content: str) -> AnalysisResult:
        """
        Extract structure confidence metrics (pLDDT, pAE).

        Note: RFD3 PDBs do NOT store pLDDT in B-factors (values are ~0.00).
        pLDDT is only available from structure prediction (AF2/ESMFold/RF3 predict).
        This method checks if B-factors contain valid pLDDT values.
        """
        # Try to extract pLDDT from B-factor column (common convention)
        plddt_values = []

        for line in pdb_content.split("\n"):
            if line.startswith("ATOM"):
                try:
                    bfactor = float(line[60:66].strip())
                    # pLDDT is typically stored as 0-100 in B-factor
                    if 0 <= bfactor <= 100:
                        plddt_values.append(bfactor / 100.0)
                except (ValueError, IndexError):
                    pass

        if plddt_values:
            mean_plddt = sum(plddt_values) / len(plddt_values)
            max_plddt = max(plddt_values)

            # Check if B-factors are valid pLDDT values
            # RFD3 PDBs have B-factors near 0 (not pLDDT)
            # Valid pLDDT from structure prediction should have higher values
            if max_plddt < 0.05:  # All values < 5% suggests not real pLDDT
                return AnalysisResult.skipped(
                    "B-factors appear to be from RFD3 (not pLDDT). "
                    "Run structure prediction for confidence scores."
                )

            return AnalysisResult.success({
                "plddt": mean_plddt,
                "plddt_min": min(plddt_values),
                "plddt_max": max_plddt,
            })
        else:
            return AnalysisResult.skipped("could not extract pLDDT from B-factors")

    def _analyze_interface(self, pdb_content: str) -> AnalysisResult:
        """Analyze protein-protein interface."""
        if not self._analysis_modules.get("binding_analysis"):
            return AnalysisResult.skipped("binding_analysis module not available")

        try:
            from binding_analysis import analyze_interface
            result = analyze_interface(pdb_content)

            if result.get("status") == "error":
                return AnalysisResult.skipped(result.get("error", "unknown error"))

            return AnalysisResult.success({
                "dSASA": result.get("dSASA_int", 0),
                "contacts": result.get("contacts", 0),
                "interface_residues": result.get("nres_int", 0),
                "interface_hbonds": result.get("hbonds_int", 0),
                "packstat": result.get("packstat", 0),
            })
        except Exception as e:
            return AnalysisResult.skipped(f"interface analysis failed: {str(e)}")

    def _analyze_metal_coordination(
        self,
        pdb_content: str,
        metal_type: str,
        metal_chain: str,
        metal_resnum: int,
        ligand_name: Optional[str] = None,
    ) -> AnalysisResult:
        """
        Analyze metal coordination geometry.

        Uses validate_metal_ligand_complex_site when a ligand is present
        (better for metal-ligand complexes like Dy-TriNOx).
        Falls back to validate_lanthanide_site for pure metal sites.
        """
        if not self._analysis_modules.get("metal_validation"):
            return AnalysisResult.skipped("metal_validation module not available")

        try:
            # Use metal-ligand complex validation when ligand is present
            if ligand_name:
                from metal_validation import validate_metal_ligand_complex_site

                result = validate_metal_ligand_complex_site(
                    pdb_content=pdb_content,
                    metal=metal_type,
                    ligand_name=ligand_name,
                    distance_cutoff=3.0,  # Tighter for lanthanides
                )

                return AnalysisResult.success({
                    "coordination_number": result.get("coordination_number", 0),
                    "ligand_coordination": result.get("ligand_coordination", 0),
                    "protein_coordination": result.get("protein_coordination", 0),
                    "ligand_donors": result.get("ligand_donors", []),
                    "protein_donors": [
                        f"{d['chain']}{d['resnum']} {d['resname']}.{d['atom']}"
                        for d in result.get("donor_residues", [])
                    ],
                    "issues": result.get("issues", []),
                    "success": result.get("success", False),
                })
            else:
                # Fall back to simple lanthanide validation for pure metal sites
                from metal_validation import validate_lanthanide_site

                result = validate_lanthanide_site(
                    pdb_content=pdb_content,
                    metal=metal_type,
                    metal_chain=metal_chain,
                    metal_resnum=metal_resnum,
                )

                return AnalysisResult.success({
                    "coordination_number": result.get("coordination_number", 0),
                    "geometry_type": result.get("geometry_type", "unknown"),
                    "geometry_rmsd": result.get("geometry_rmsd", 0),
                    "coordination_distance": result.get("mean_distance", 0),
                    "hsab_compatible": result.get("hsab_compatible", False),
                    "coordinating_residues": result.get("coordinating_residues", []),
                })
        except Exception as e:
            return AnalysisResult.skipped(f"metal coordination analysis failed: {str(e)}")

    def _analyze_ligand_binding(
        self,
        pdb_content: str,
        pdb_path: Optional[str],
        ligand_sdf: str,
    ) -> AnalysisResult:
        """Analyze protein-ligand binding with GNINA."""
        if not self._analysis_modules.get("gnina"):
            return AnalysisResult.skipped("GNINA binary not available")

        if not self._analysis_modules.get("binding_analysis"):
            return AnalysisResult.skipped("binding_analysis module not available")

        try:
            from binding_analysis import run_gnina_scoring

            result = run_gnina_scoring(
                receptor_pdb=pdb_content,
                ligand_sdf=ligand_sdf,
            )

            if result.get("status") == "error":
                return AnalysisResult.skipped(result.get("error", "GNINA failed"))

            return AnalysisResult.success({
                "gnina_affinity": result.get("affinity", 0),
                "gnina_cnn_score": result.get("cnn_score", 0),
                "gnina_cnn_affinity": result.get("cnn_affinity", 0),
            })
        except Exception as e:
            return AnalysisResult.skipped(f"ligand binding analysis failed: {str(e)}")

    def _analyze_ligand_interactions(
        self,
        pdb_content: str,
        ligand_name: str,
    ) -> AnalysisResult:
        """
        Analyze protein-ligand interactions using PLIP.

        Detects:
        - Hydrogen bonds (with angle validation)
        - Hydrophobic contacts
        - Pi-stacking (face-to-face and edge-to-face)
        - Salt bridges
        - Halogen bonds

        Args:
            pdb_content: PDB file content
            ligand_name: Ligand residue name (e.g., "UNL")

        Returns:
            AnalysisResult with detailed interaction metrics
        """
        if not self._analysis_modules.get("interaction_analysis"):
            return AnalysisResult.skipped("interaction_analysis module not available")

        try:
            # Use PLIP-based analysis from shared module
            summary = analyze_all_interactions(
                pdb_content=pdb_content,
                ligand_name=ligand_name,
                include_visualization_data=False,
                use_plip=True,
            )

            if summary.status == "error":
                return AnalysisResult.skipped(f"PLIP analysis failed: {summary.error}")

            # Format H-bond details
            hbond_details = []
            for hb in summary.hydrogen_bonds[:10]:  # Top 10
                hbond_details.append({
                    "protein": f"{hb.protein_chain}:{hb.protein_residue}.{hb.protein_atom}",
                    "ligand_atom": hb.ligand_atom,
                    "distance": hb.distance,
                    "angle": hb.angle,
                    "type": hb.type,
                })

            # Format hydrophobic contact details
            hydrophobic_details = []
            for hc in summary.hydrophobic_contacts[:10]:  # Top 10
                hydrophobic_details.append({
                    "protein": f"{hc.protein_chain}:{hc.protein_residue}.{hc.protein_atom}",
                    "ligand_atom": hc.ligand_atom,
                    "distance": hc.distance,
                })

            # Format pi-stacking details
            pi_stacking_details = []
            for ps in summary.pi_stacking:
                pi_stacking_details.append({
                    "protein": f"{ps.protein_chain}:{ps.protein_residue}",
                    "type": ps.type,
                    "distance": ps.distance,
                })

            # Format salt bridge details
            salt_bridge_details = []
            for sb in summary.salt_bridges:
                salt_bridge_details.append({
                    "protein": f"{sb.protein_chain}:{sb.protein_residue}",
                    "distance": sb.distance,
                    "ligand_charge": sb.ligand_charge,
                })

            return AnalysisResult.success({
                "total_contacts": summary.total_contacts,
                "hydrogen_bonds": len(summary.hydrogen_bonds),
                "hydrophobic_contacts": len(summary.hydrophobic_contacts),
                "pi_stacking": len(summary.pi_stacking),
                "salt_bridges": len(summary.salt_bridges),
                "halogen_bonds": len(summary.halogen_bonds),
                "key_residues": summary.key_residues,
                "analysis_method": summary.analysis_method,
                "hbond_details": hbond_details,
                "hydrophobic_details": hydrophobic_details,
                "pi_stacking_details": pi_stacking_details,
                "salt_bridge_details": salt_bridge_details,
            })
        except Exception as e:
            return AnalysisResult.skipped(f"ligand interaction analysis failed: {str(e)}")

    def _analyze_topology(self, pdb_content: str) -> AnalysisResult:
        """Validate topology (ring closure, etc.)."""
        if not self._analysis_modules.get("topology_validation"):
            return AnalysisResult.skipped("topology_validation module not available")

        try:
            from topology_validation import validate_topology

            result = validate_topology(pdb_content)

            return AnalysisResult.success({
                "valid": result.get("valid", False),
                "issues": result.get("issues", []),
            })
        except Exception as e:
            return AnalysisResult.skipped(f"topology validation failed: {str(e)}")

    def _analyze_symmetry(self, pdb_content: str) -> AnalysisResult:
        """Analyze chain similarity for dimers (basic check via length comparison)."""
        # Basic symmetry analysis - compare chain lengths as a simple heuristic
        try:
            chain_coords = {}
            for line in pdb_content.split("\n"):
                if line.startswith("ATOM") and " CA " in line:
                    chain_id = line[21]
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    if chain_id not in chain_coords:
                        chain_coords[chain_id] = []
                    chain_coords[chain_id].append((x, y, z))

            chains = list(chain_coords.keys())
            if len(chains) < 2:
                return AnalysisResult.not_applicable("need at least 2 chains for symmetry")

            # Simple symmetry score based on chain length similarity
            len_a = len(chain_coords[chains[0]])
            len_b = len(chain_coords[chains[1]])
            length_ratio = min(len_a, len_b) / max(len_a, len_b) if max(len_a, len_b) > 0 else 0

            return AnalysisResult.success({
                "c2_score": length_ratio,
                "chain_lengths": {chains[0]: len_a, chains[1]: len_b},
            })
        except Exception as e:
            return AnalysisResult.skipped(f"symmetry analysis failed: {str(e)}")

    def _analyze_pyrosetta_interface(
        self,
        pdb_content: str,
        chain_a: str = "A",
        chain_b: str = "B",
    ) -> AnalysisResult:
        """
        Analyze protein interface using PyRosetta (BindCraft-style metrics).

        Calculates:
        - dG: Binding energy (negative = favorable)
        - dSASA: Buried surface area
        - sc_value: Shape complementarity (>0.6 good)
        - packstat: Packing quality (>0.6 good)
        - n_InterfaceHbonds: Hydrogen bond count (≥3 good)

        Args:
            pdb_content: PDB file content
            chain_a: First chain ID
            chain_b: Second chain ID

        Returns:
            AnalysisResult with interface metrics
        """
        if not self._analysis_modules.get("pyrosetta"):
            return AnalysisResult.skipped("PyRosetta not available")

        try:
            import rosetta_utils

            result = rosetta_utils.score_interface(
                pdb_content=pdb_content,
                chain_a=chain_a,
                chain_b=chain_b,
            )

            if result.get("status") == "error":
                return AnalysisResult.skipped(f"PyRosetta error: {result.get('error')}")

            # Extract metrics
            metrics = {
                "dG": result.get("dG"),
                "dSASA": result.get("dSASA"),
                "sc_value": result.get("sc_value"),
                "packstat": result.get("packstat"),
            }

            # Add optional H-bond metrics if available
            if "n_InterfaceHbonds" in result:
                metrics["n_InterfaceHbonds"] = result["n_InterfaceHbonds"]
            if "n_InterfaceUnsatHbonds" in result:
                metrics["n_InterfaceUnsatHbonds"] = result["n_InterfaceUnsatHbonds"]
            if "interface_energy" in result:
                metrics["interface_energy"] = result["interface_energy"]

            # Evaluate quality based on BindCraft thresholds
            passes_filters = True
            issues = []

            if metrics.get("dG") is not None and metrics["dG"] > 0:
                issues.append(f"Unfavorable binding energy: dG={metrics['dG']:.1f} (want ≤0)")
                passes_filters = False

            if metrics.get("sc_value") is not None and metrics["sc_value"] < 0.6:
                issues.append(f"Low shape complementarity: {metrics['sc_value']:.2f} (want ≥0.6)")
                passes_filters = False

            if metrics.get("n_InterfaceHbonds") is not None and metrics["n_InterfaceHbonds"] < 3:
                issues.append(f"Few interface H-bonds: {metrics['n_InterfaceHbonds']} (want ≥3)")

            metrics["passes_bindcraft_filters"] = passes_filters
            metrics["issues"] = issues

            return AnalysisResult.success(metrics)

        except Exception as e:
            return AnalysisResult.skipped(f"PyRosetta interface analysis failed: {str(e)}")

    def _analyze_sequence_composition(self, pdb_content: str) -> AnalysisResult:
        """
        Analyze sequence composition for potential design issues.

        Checks:
        - High Alanine content (>25% = warning, >35% = concern)
        - Aromatic residue count (W, Y, F) for binding pocket quality
        - Overall amino acid distribution

        Args:
            pdb_content: PDB file content

        Returns:
            AnalysisResult with composition stats and warnings
        """
        try:
            # Extract sequence from CA atoms (one per residue)
            residues = []
            seen_residues = set()  # Track (chain, resnum) to avoid duplicates

            # Standard amino acid 3-letter to 1-letter mapping
            aa_map = {
                'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
                'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
                'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
                'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
            }

            for line in pdb_content.split("\n"):
                if line.startswith("ATOM") and " CA " in line:
                    try:
                        res_name = line[17:20].strip()
                        chain_id = line[21] if len(line) > 21 else 'A'
                        res_num = int(line[22:26])

                        res_key = (chain_id, res_num)
                        if res_key not in seen_residues and res_name in aa_map:
                            seen_residues.add(res_key)
                            residues.append(res_name)
                    except (ValueError, IndexError):
                        continue

            if not residues:
                return AnalysisResult.skipped("no protein residues found")

            total = len(residues)

            # Count each amino acid
            aa_counts = {}
            for res in residues:
                aa_counts[res] = aa_counts.get(res, 0) + 1

            # Calculate percentages
            aa_percentages = {aa: (count / total) * 100 for aa, count in aa_counts.items()}

            # Specific checks
            ala_count = aa_counts.get('ALA', 0)
            ala_pct = (ala_count / total) * 100

            # Aromatic residues (important for ligand binding pockets)
            trp_count = aa_counts.get('TRP', 0)
            tyr_count = aa_counts.get('TYR', 0)
            phe_count = aa_counts.get('PHE', 0)
            aromatic_count = trp_count + tyr_count + phe_count
            aromatic_pct = (aromatic_count / total) * 100

            # Hydrophobic residues
            hydrophobic_aas = ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'TYR', 'PRO']
            hydrophobic_count = sum(aa_counts.get(aa, 0) for aa in hydrophobic_aas)
            hydrophobic_pct = (hydrophobic_count / total) * 100

            # Charged residues
            charged_aas = ['ASP', 'GLU', 'LYS', 'ARG', 'HIS']
            charged_count = sum(aa_counts.get(aa, 0) for aa in charged_aas)
            charged_pct = (charged_count / total) * 100

            # Build warnings list
            warnings = []

            # Alanine check (main concern from v11 designs)
            if ala_pct > 35:
                warnings.append(f"HIGH Alanine content: {ala_pct:.1f}% ({ala_count}/{total}) - may indicate poor sequence design")
            elif ala_pct > 25:
                warnings.append(f"Elevated Alanine content: {ala_pct:.1f}% ({ala_count}/{total}) - consider aromatic bias for binding pocket")

            # Low aromatic check (important for hydrophobic ligand binding)
            if aromatic_pct < 3:
                warnings.append(f"Low aromatic content: {aromatic_pct:.1f}% ({aromatic_count} W/Y/F) - may weaken ligand binding pocket")

            # Typical protein composition reference (for context)
            # Ala ~8%, Trp ~1%, Tyr ~3%, Phe ~4% in natural proteins

            return AnalysisResult.success({
                "total_residues": total,
                "alanine": {
                    "count": ala_count,
                    "percentage": round(ala_pct, 1),
                },
                "aromatics": {
                    "total": aromatic_count,
                    "percentage": round(aromatic_pct, 1),
                    "W": trp_count,
                    "Y": tyr_count,
                    "F": phe_count,
                },
                "hydrophobic": {
                    "count": hydrophobic_count,
                    "percentage": round(hydrophobic_pct, 1),
                },
                "charged": {
                    "count": charged_count,
                    "percentage": round(charged_pct, 1),
                },
                "warnings": warnings,
                "composition": {aa_map.get(aa, aa): count for aa, count in sorted(aa_counts.items())},
            })

        except Exception as e:
            return AnalysisResult.skipped(f"sequence composition analysis failed: {str(e)}")

    # =========================================================================
    # Enzyme Activity Preservation Analysis (Phase 6: AI Infrastructure)
    # =========================================================================

    def _analyze_active_site_rmsd(
        self,
        designed_pdb: str,
        reference_pdb: str,
        catalytic_residues: List[str],
    ) -> AnalysisResult:
        """
        Calculate RMSD of catalytic residues between designed and reference structures.

        Args:
            designed_pdb: Designed protein PDB content
            reference_pdb: Reference (template) PDB content
            catalytic_residues: List of residue identifiers (e.g., ["A45", "A78", "A123"])

        Returns:
            AnalysisResult with RMSD metrics. Pass threshold: <1.5Å
        """
        try:
            import math

            def parse_ca_coords(pdb_content: str, residues: List[str]) -> Dict[str, tuple]:
                """Extract CA coordinates for specified residues."""
                coords = {}
                for line in pdb_content.split("\n"):
                    if line.startswith("ATOM") and " CA " in line:
                        try:
                            chain = line[21] if len(line) > 21 else 'A'
                            resnum = int(line[22:26])
                            res_id = f"{chain}{resnum}"
                            if res_id in residues:
                                x = float(line[30:38])
                                y = float(line[38:46])
                                z = float(line[46:54])
                                coords[res_id] = (x, y, z)
                        except (ValueError, IndexError):
                            continue
                return coords

            # Get CA coordinates for catalytic residues
            designed_coords = parse_ca_coords(designed_pdb, catalytic_residues)
            reference_coords = parse_ca_coords(reference_pdb, catalytic_residues)

            if not designed_coords or not reference_coords:
                return AnalysisResult.skipped(
                    f"Could not find catalytic residues in structures. "
                    f"Designed: {len(designed_coords)}, Reference: {len(reference_coords)}"
                )

            # Calculate RMSD for matching residues
            matched_residues = []
            squared_distances = []

            for res_id in catalytic_residues:
                if res_id in designed_coords and res_id in reference_coords:
                    d_coords = designed_coords[res_id]
                    r_coords = reference_coords[res_id]
                    sq_dist = sum((d - r) ** 2 for d, r in zip(d_coords, r_coords))
                    squared_distances.append(sq_dist)
                    matched_residues.append(res_id)

            if not squared_distances:
                return AnalysisResult.skipped("No matching catalytic residues found")

            rmsd = math.sqrt(sum(squared_distances) / len(squared_distances))

            # Evaluate against threshold
            passed = rmsd < 1.5
            status = "good" if rmsd < 1.0 else ("acceptable" if rmsd < 1.5 else "poor")

            return AnalysisResult.success({
                "active_site_rmsd": round(rmsd, 3),
                "matched_residues": matched_residues,
                "total_catalytic": len(catalytic_residues),
                "matched_count": len(matched_residues),
                "passed": passed,
                "status": status,
                "threshold": 1.5,
            })

        except Exception as e:
            return AnalysisResult.skipped(f"active site RMSD analysis failed: {str(e)}")

    def _analyze_substrate_access(
        self,
        pdb_content: str,
        exposed_atoms: Dict[str, str],
    ) -> AnalysisResult:
        """
        Validate substrate channel accessibility via RASA calculation.

        Args:
            pdb_content: PDB content to analyze
            exposed_atoms: Dict of chain->atoms that should remain exposed (e.g., {"L1": "O1,O2,C3"})

        Returns:
            AnalysisResult with accessibility metrics. Pass threshold: >30% exposed
        """
        try:
            # Simple accessibility check based on atom presence and burial
            # For full RASA calculation, would need external tool like FreeSASA

            # Parse atoms from PDB
            atom_coords = {}
            for line in pdb_content.split("\n"):
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    try:
                        chain = line[21] if len(line) > 21 else 'A'
                        atom_name = line[12:16].strip()
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        atom_coords[(chain, atom_name)] = (x, y, z)
                    except (ValueError, IndexError):
                        continue

            # Check which exposed atoms are present
            found_exposed = []
            missing_exposed = []

            for chain, atoms_str in exposed_atoms.items():
                for atom_name in atoms_str.split(","):
                    atom_name = atom_name.strip()
                    if (chain, atom_name) in atom_coords:
                        found_exposed.append(f"{chain}:{atom_name}")
                    else:
                        missing_exposed.append(f"{chain}:{atom_name}")

            total_expected = len(found_exposed) + len(missing_exposed)
            if total_expected == 0:
                return AnalysisResult.skipped("No exposed atoms specified")

            # Simple accessibility metric: presence ratio
            # Full RASA calculation would require external library
            presence_ratio = len(found_exposed) / total_expected if total_expected > 0 else 0

            # For now, report presence check
            # TODO: Integrate FreeSASA for proper RASA calculation
            passed = presence_ratio > 0.7  # At least 70% of expected atoms present

            return AnalysisResult.success({
                "substrate_channel": {
                    "atoms_found": len(found_exposed),
                    "atoms_expected": total_expected,
                    "presence_ratio": round(presence_ratio, 2),
                    "found": found_exposed[:10],  # First 10 for brevity
                    "missing": missing_exposed[:10],
                },
                "passed": passed,
                "note": "Basic presence check. Full RASA analysis requires FreeSASA integration.",
            })

        except Exception as e:
            return AnalysisResult.skipped(f"substrate access analysis failed: {str(e)}")

    def _analyze_enzyme_preservation(
        self,
        designed_pdb: str,
        enzyme_class: str,
        reference_pdb: Optional[str] = None,
        catalytic_residues: Optional[List[str]] = None,
        exposed_atoms: Optional[Dict[str, str]] = None,
    ) -> AnalysisResult:
        """
        Combined enzyme activity preservation analysis.

        Evaluates:
        1. Active site RMSD (if reference provided)
        2. Substrate channel accessibility (if exposed atoms specified)
        3. Enzyme class-specific checks

        Args:
            designed_pdb: Designed protein PDB content
            enzyme_class: Enzyme class name (e.g., "dehydrogenase", "quinoprotein")
            reference_pdb: Optional reference PDB for RMSD comparison
            catalytic_residues: Optional list of catalytic residue IDs
            exposed_atoms: Optional dict of atoms that should remain exposed

        Returns:
            AnalysisResult with combined enzyme preservation metrics
        """
        try:
            results = {
                "enzyme_class": enzyme_class,
                "checks_performed": [],
                "issues": [],
            }

            overall_passed = True

            # Check 1: Active site RMSD
            if reference_pdb and catalytic_residues:
                rmsd_result = self._analyze_active_site_rmsd(
                    designed_pdb, reference_pdb, catalytic_residues
                )
                if rmsd_result.status == AnalysisStatus.SUCCESS:
                    results["active_site_rmsd"] = rmsd_result.metrics
                    results["checks_performed"].append("active_site_rmsd")
                    if not rmsd_result.metrics.get("passed", False):
                        overall_passed = False
                        results["issues"].append(
                            f"Active site RMSD {rmsd_result.metrics.get('active_site_rmsd', '?')}Å > 1.5Å threshold"
                        )
                else:
                    results["active_site_rmsd_skipped"] = rmsd_result.message

            # Check 2: Substrate channel accessibility
            if exposed_atoms:
                access_result = self._analyze_substrate_access(designed_pdb, exposed_atoms)
                if access_result.status == AnalysisStatus.SUCCESS:
                    results["substrate_access"] = access_result.metrics
                    results["checks_performed"].append("substrate_access")
                    if not access_result.metrics.get("passed", False):
                        overall_passed = False
                        results["issues"].append("Substrate channel atoms may be blocked")
                else:
                    results["substrate_access_skipped"] = access_result.message

            # Check 3: Enzyme class-specific validation
            # Import enzyme chemistry for class-specific checks
            try:
                from enzyme_chemistry import get_preservation_requirements
                preservation_reqs = get_preservation_requirements(enzyme_class)

                if preservation_reqs:
                    results["enzyme_requirements"] = {
                        "cofactor_types": preservation_reqs.get("cofactor_types", []),
                        "hbond_network_critical": preservation_reqs.get("hbond_network_critical", False),
                        "metal_role": preservation_reqs.get("metal_role", "none"),
                        "substrate_access_required": preservation_reqs.get("substrate_access_required", False),
                    }
                    results["checks_performed"].append("enzyme_requirements")
            except ImportError:
                results["enzyme_requirements_skipped"] = "enzyme_chemistry module not available"

            # Overall assessment
            results["passed"] = overall_passed and len(results["checks_performed"]) > 0
            results["preservation_score"] = (
                1.0 if overall_passed else 0.5
            ) if results["checks_performed"] else 0.0

            # Status message
            if not results["checks_performed"]:
                results["status"] = "no_checks_available"
            elif overall_passed:
                results["status"] = "preserved"
            else:
                results["status"] = "concerns"

            return AnalysisResult.success(results)

        except Exception as e:
            return AnalysisResult.skipped(f"enzyme preservation analysis failed: {str(e)}")

    # =========================================================================
    # Evolutionary Conservation Analysis (ConSurf Methodology)
    # =========================================================================

    def _analyze_conservation(
        self,
        pdb_content: str,
        chain: str = "A",
    ) -> AnalysisResult:
        """
        Analyze evolutionary conservation using ConSurf methodology.

        Uses:
        - NCBI BLAST API for homolog search
        - MMseqs2 for redundancy filtering
        - MUSCLE for alignment
        - Rate4Site for conservation scoring (or Shannon entropy fallback)

        Conservation grades follow ConSurf 1-9 scale:
        - Grade 1-3: Highly conserved (functionally important)
        - Grade 4-6: Average conservation
        - Grade 7-9: Variable (can be redesigned)

        Args:
            pdb_content: PDB file content
            chain: Chain ID to analyze

        Returns:
            AnalysisResult with conservation grades and suggestions for design
        """
        if not self._analysis_modules.get("conservation"):
            return AnalysisResult.skipped("conservation analyzer module not available")

        try:
            import asyncio
            from conservation_analyzer import (
                ConservationAnalyzer,
                get_conservation_summary,
            )

            analyzer = ConservationAnalyzer()

            # Run async analysis (may take time for BLAST)
            try:
                loop = asyncio.get_event_loop()
            except RuntimeError:
                loop = asyncio.new_event_loop()
                asyncio.set_event_loop(loop)

            analysis = loop.run_until_complete(
                analyzer.analyze(pdb_content, chain=chain)
            )

            # Get summary for design integration
            summary = get_conservation_summary(analysis)

            return AnalysisResult.success({
                "grades": [g.to_dict() for g in analysis.grades],
                "highly_conserved_positions": analysis.highly_conserved,
                "conserved_positions": analysis.conserved,
                "variable_positions": analysis.variable,
                "msa_depth": analysis.msa_depth,
                "average_conservation": round(analysis.average_conservation, 2),
                "method": analysis.method,
                "reliable": summary["reliable_analysis"],
                "suggested_fixed_residues": summary["suggested_fixed_residues"],
                "notes": analysis.notes,
            })

        except ValueError as e:
            # Insufficient homologs or sequence extraction failed
            return AnalysisResult.skipped(f"conservation analysis failed: {str(e)}")

        except Exception as e:
            return AnalysisResult.skipped(f"conservation analysis error: {str(e)}")
