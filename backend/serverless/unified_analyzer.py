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

from analysis_types import AnalysisResult, AnalysisStatus, DesignType, detect_design_type


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

        return modules

    def _detect_design_type(
        self,
        pdb_content: str,
        metal_type: Optional[str] = None,
        ligand_sdf: Optional[str] = None,
    ) -> DesignType:
        """Detect design type from PDB content and provided metadata."""
        # Count chains
        chains = set()
        has_hetatm = False

        for line in pdb_content.split("\n"):
            if line.startswith("ATOM"):
                chain_id = line[21] if len(line) > 21 else ""
                if chain_id.strip():
                    chains.add(chain_id)
            elif line.startswith("HETATM"):
                has_hetatm = True

        chain_count = len(chains)
        has_ligand = ligand_sdf is not None or has_hetatm
        has_metal = metal_type is not None

        return detect_design_type(
            has_ligand=has_ligand,
            has_metal=has_metal,
            chain_count=chain_count,
        )

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
        metal_chain: str = "L",
        metal_resnum: int = 1,
    ) -> Dict[str, Any]:
        """
        Run comprehensive analysis on a design.

        Args:
            pdb_content: PDB file content as string
            design_params: Parameters used to generate the design
            pdb_path: Optional path to PDB file (for tools that need file path)
            ligand_sdf: Optional SDF content for ligand
            metal_type: Optional metal type code (e.g., "TB", "ZN")
            metal_chain: Chain ID of metal (default "L")
            metal_resnum: Residue number of metal (default 1)

        Returns:
            Structured metrics JSON with all analysis results
        """
        # Detect design type
        design_type = self._detect_design_type(pdb_content, metal_type, ligand_sdf)
        design_id = self._generate_design_id(design_type)

        # Initialize result structure
        result = {
            "design_id": design_id,
            "design_type": design_type.value,
            "timestamp": datetime.now().isoformat(),
            "design_params": design_params,
            "analyses": {},
        }

        # Run structure confidence analysis
        result["analyses"]["structure_confidence"] = self._analyze_structure_confidence(
            pdb_content
        ).to_dict()

        # Run interface analysis (for dimers)
        if design_type != DesignType.MONOMER:
            result["analyses"]["interface_quality"] = self._analyze_interface(
                pdb_content
            ).to_dict()
        else:
            result["analyses"]["interface_quality"] = AnalysisResult.not_applicable(
                "monomer design - no interface"
            ).to_dict()

        # Run metal coordination analysis
        if metal_type:
            result["analyses"]["metal_coordination"] = self._analyze_metal_coordination(
                pdb_content, metal_type, metal_chain, metal_resnum
            ).to_dict()
        else:
            result["analyses"]["metal_coordination"] = AnalysisResult.not_applicable(
                "no metal type specified"
            ).to_dict()

        # Run ligand binding analysis
        if ligand_sdf:
            result["analyses"]["ligand_binding"] = self._analyze_ligand_binding(
                pdb_content, pdb_path, ligand_sdf
            ).to_dict()
        else:
            result["analyses"]["ligand_binding"] = AnalysisResult.not_applicable(
                "no ligand provided"
            ).to_dict()

        # Run topology validation
        result["analyses"]["topology"] = self._analyze_topology(pdb_content).to_dict()

        # Run symmetry analysis (for dimers)
        if design_type != DesignType.MONOMER:
            result["analyses"]["symmetry"] = self._analyze_symmetry(pdb_content).to_dict()
        else:
            result["analyses"]["symmetry"] = AnalysisResult.not_applicable(
                "monomer design - no symmetry analysis"
            ).to_dict()

        return result

    def _analyze_structure_confidence(self, pdb_content: str) -> AnalysisResult:
        """Extract structure confidence metrics (pLDDT, pAE)."""
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
            return AnalysisResult.success({
                "plddt": sum(plddt_values) / len(plddt_values),
                "plddt_min": min(plddt_values),
                "plddt_max": max(plddt_values),
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
    ) -> AnalysisResult:
        """Analyze metal coordination geometry."""
        if not self._analysis_modules.get("metal_validation"):
            return AnalysisResult.skipped("metal_validation module not available")

        try:
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
