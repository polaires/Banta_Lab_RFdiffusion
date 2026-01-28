"""
AI-Powered Protein Design Pipeline

Orchestrates the full design pipeline from natural language input to final report:
1. NL Parse → DesignIntent
2. Ligand Resolve → ResolvedLigand
3. Config Generate → RFD3Config
4. RFD3 Backbone → PDB designs
5. LigandMPNN Sequence → Designed sequences
6. ESMFold Validation → Structure verification
7. Unified Analysis → Comprehensive metrics
8. Report Generation → Human-readable summary

Usage:
    pipeline = AIDesignPipeline(claude_api_key="...", history_dir="/path/to/history")
    result = pipeline.run("Design a protein to bind citrate with terbium")
"""

import logging
import time
import json
from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional, Tuple
from pathlib import Path

# Import our modules
from nl_design_parser import NLDesignParser, DesignIntent, create_parser
from ligand_resolver import LigandResolver, ResolvedLigand, resolve_ligand
from rfd3_config_generator import (
    RFD3ConfigGenerator,
    RFD3Config,
    LigandMPNNConfig,
    generate_rfd3_config,
    generate_mpnn_config,
)

# Import stability profiles from decision engine
from design_rules import select_stability_profile, StabilityProfile

# Import enzyme chemistry for activity preservation
try:
    from enzyme_chemistry import (
        detect_enzyme_class,
        get_preservation_requirements,
        get_substrate_channel_atoms,
        get_catalytic_hbond_network,
    )
    ENZYME_CHEMISTRY_AVAILABLE = True
except ImportError:
    ENZYME_CHEMISTRY_AVAILABLE = False
    detect_enzyme_class = None
    get_preservation_requirements = None
    get_substrate_channel_atoms = None
    get_catalytic_hbond_network = None

# Import inference utilities
from inference_utils import run_rfd3_inference, run_mpnn_inference

# Import scaffolding workflow for motif-based design
try:
    from scaffolding_workflow import ScaffoldingWorkflow, ScaffoldResult
    SCAFFOLDING_AVAILABLE = True
except ImportError:
    SCAFFOLDING_AVAILABLE = False
    ScaffoldingWorkflow = None
    ScaffoldResult = None

# Import ESMFold validation
try:
    from esmfold_utils import validate_structure_esmfold, is_esmfold_available
    ESMFOLD_AVAILABLE = is_esmfold_available()
except ImportError:
    ESMFOLD_AVAILABLE = False
    validate_structure_esmfold = None

# Import RF3 validation
try:
    from esmfold_utils import validate_structure_rf3, is_rf3_available
    RF3_AVAILABLE = is_rf3_available()
except (ImportError, AttributeError):
    RF3_AVAILABLE = False
    validate_structure_rf3 = None

# Import unified analyzer
try:
    from unified_analyzer import UnifiedDesignAnalyzer
    ANALYZER_AVAILABLE = True
except ImportError:
    ANALYZER_AVAILABLE = False
    UnifiedDesignAnalyzer = None

# Import design history
try:
    from design_history import DesignHistoryManager
    HISTORY_AVAILABLE = True
except ImportError:
    HISTORY_AVAILABLE = False
    DesignHistoryManager = None

# Import lesson detector
try:
    from lesson_detector import LessonDetector
    LESSON_DETECTOR_AVAILABLE = True
except ImportError:
    LESSON_DETECTOR_AVAILABLE = False
    LessonDetector = None

logger = logging.getLogger(__name__)


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class SequenceResult:
    """Result from sequence design."""
    sequence: str
    confidence: float = 0.0
    filename: str = ""
    pdb_content: Optional[str] = None


@dataclass
class ValidationResult:
    """Result from structure validation."""
    passed: bool = False
    rmsd: float = float('inf')
    plddt: float = 0.0
    predicted_pdb: Optional[str] = None
    # RF3 metrics
    rf3_rmsd: Optional[float] = None
    rf3_plddt: Optional[float] = None
    rf3_ptm: Optional[float] = None
    # Best-of metrics (across ESMFold and RF3)
    best_rmsd: Optional[float] = None
    best_plddt: Optional[float] = None
    warnings: List[str] = field(default_factory=list)


@dataclass
class PipelineResult:
    """Complete result from AI design pipeline."""
    # Status
    success: bool = False
    error: Optional[str] = None

    # Input parsing
    design_intent: Optional[DesignIntent] = None
    resolved_ligand: Optional[ResolvedLigand] = None

    # Configurations used
    rfd3_config: Optional[RFD3Config] = None
    mpnn_config: Optional[LigandMPNNConfig] = None

    # Backbone designs
    backbone_pdbs: List[str] = field(default_factory=list)
    num_backbones: int = 0

    # Sequence designs
    sequences: List[SequenceResult] = field(default_factory=list)
    num_sequences: int = 0

    # Best results
    best_sequence: Optional[SequenceResult] = None
    best_sequence_pdb: Optional[str] = None
    best_rmsd: float = float('inf')
    best_plddt: float = 0.0

    # Validation
    validation_results: List[ValidationResult] = field(default_factory=list)
    pass_rate: float = 0.0

    # Analysis
    analysis_results: Dict[str, Any] = field(default_factory=dict)

    # Report
    report: str = ""
    recommendations: List[str] = field(default_factory=list)

    # Lessons
    lessons_triggered: List[str] = field(default_factory=list)

    # Timing
    timings: Dict[str, float] = field(default_factory=dict)
    total_time: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for API response."""
        result = {
            "success": self.success,
            "error": self.error,
            "num_backbones": self.num_backbones,
            "num_sequences": self.num_sequences,
            "pass_rate": self.pass_rate,
            "best_rmsd": self.best_rmsd if self.best_rmsd < float('inf') else None,
            "best_plddt": self.best_plddt,
            "report": self.report,
            "recommendations": self.recommendations,
            "timings": self.timings,
            "total_time": self.total_time,
        }

        # Include design intent if available
        if self.design_intent:
            result["design_intent"] = self.design_intent.to_dict()

        # Include best sequence if available
        if self.best_sequence:
            result["best_sequence"] = {
                "sequence": self.best_sequence.sequence,
                "confidence": self.best_sequence.confidence,
            }

        # Include best PDB if available
        if self.best_sequence_pdb:
            result["best_sequence_pdb"] = self.best_sequence_pdb

        # Include backbone PDBs
        if self.backbone_pdbs:
            result["backbone_pdbs"] = self.backbone_pdbs

        # Include all sequences
        if self.sequences:
            result["sequences"] = [
                {"sequence": s.sequence, "confidence": s.confidence}
                for s in self.sequences
            ]

        # Include analysis if available
        if self.analysis_results:
            result["analysis"] = self.analysis_results

        return result


# =============================================================================
# Main Pipeline Class
# =============================================================================

class AIDesignPipeline:
    """
    AI-powered protein design pipeline.

    Orchestrates the full design workflow from natural language input to
    validated protein structures with comprehensive analysis.
    """

    def __init__(
        self,
        claude_api_key: Optional[str] = None,
        claude_api_url: str = "https://yinli.one/v1/messages",
        history_dir: Optional[str] = None,
        use_mock: bool = False,
    ):
        """
        Initialize the pipeline.

        Args:
            claude_api_key: API key for Claude NL parsing
            claude_api_url: API endpoint for Claude
            history_dir: Directory for design history (optional)
            use_mock: Use mock inference for testing
        """
        self.claude_api_key = claude_api_key
        self.claude_api_url = claude_api_url
        self.history_dir = history_dir
        self.use_mock = use_mock

        # Initialize components
        self.parser = create_parser(
            api_key=claude_api_key,
            api_url=claude_api_url,
            use_fallback=True,
        )
        self.resolver = LigandResolver()
        self.config_generator = RFD3ConfigGenerator()

        # Initialize optional components
        self.history_manager = None
        if HISTORY_AVAILABLE and history_dir:
            self.history_manager = DesignHistoryManager(history_dir)

        self.lesson_detector = None
        if LESSON_DETECTOR_AVAILABLE and history_dir:
            lessons_path = Path(history_dir) / "lessons"
            if lessons_path.exists():
                self.lesson_detector = LessonDetector(str(lessons_path))

        self.analyzer = None
        if ANALYZER_AVAILABLE:
            self.analyzer = UnifiedDesignAnalyzer()

        # Initialize scaffolding workflow
        self.scaffolding_workflow = None
        if SCAFFOLDING_AVAILABLE:
            self.scaffolding_workflow = ScaffoldingWorkflow()

    def run(
        self,
        query: str,
        num_designs: int = 4,
        num_sequences: int = 8,
        validate_sequences: bool = True,
        session_name: Optional[str] = None,
    ) -> PipelineResult:
        """
        Run the full design pipeline from natural language query.

        Args:
            query: Natural language design request
            num_designs: Number of backbone designs to generate
            num_sequences: Number of sequences per backbone
            validate_sequences: Whether to validate with ESMFold
            session_name: Optional session name for history tracking

        Returns:
            PipelineResult with all outputs
        """
        result = PipelineResult()
        start_time = time.time()

        try:
            # Stage 1: Parse natural language
            stage_start = time.time()
            intent = self._parse_query(query)
            result.design_intent = intent
            result.timings["parse"] = time.time() - stage_start

            if intent.confidence < 0.3:
                result.error = f"Low confidence parsing: {intent.warnings}"
                logger.warning(f"Low confidence parsing: {intent.confidence}")

            # === DESIGN MODE ROUTER ===
            # Route to scaffolding workflow if PDB ID and scaffolding intent detected
            if intent.is_scaffolding and self.scaffolding_workflow:
                logger.info(f"Routing to scaffolding workflow: PDB={intent.source_pdb_id}")
                return self._run_scaffolding_pipeline(
                    intent=intent,
                    result=result,
                    start_time=start_time,
                    num_designs=num_designs,
                    num_sequences=num_sequences,
                    validate_sequences=validate_sequences,
                )
            # Otherwise continue with de novo design workflow

            # Stage 2: Resolve ligand
            stage_start = time.time()
            ligand = self._resolve_ligand(intent)
            result.resolved_ligand = ligand
            result.timings["resolve"] = time.time() - stage_start

            if not ligand.resolved:
                result.error = f"Could not resolve ligand: {ligand.warnings}"
                logger.error(f"Ligand resolution failed: {ligand.warnings}")
                return self._finalize_result(result, start_time)

            # Stage 3: Generate configurations
            stage_start = time.time()
            rfd3_config = self._generate_rfd3_config(intent, ligand, num_designs)
            mpnn_config = self._generate_mpnn_config(intent, ligand, num_sequences)
            result.rfd3_config = rfd3_config
            result.mpnn_config = mpnn_config
            result.timings["config"] = time.time() - stage_start

            # Stage 4: Run RFD3 backbone generation
            stage_start = time.time()
            backbone_pdbs = self._run_rfd3(rfd3_config)
            result.backbone_pdbs = backbone_pdbs
            result.num_backbones = len(backbone_pdbs)
            result.timings["rfd3"] = time.time() - stage_start

            if not backbone_pdbs:
                result.error = "RFD3 backbone generation failed"
                return self._finalize_result(result, start_time)

            # Stage 5: Run LigandMPNN sequence design
            stage_start = time.time()
            all_sequences = []
            for backbone_pdb in backbone_pdbs:
                sequences = self._run_mpnn(backbone_pdb, mpnn_config)
                all_sequences.extend(sequences)
            result.sequences = all_sequences
            result.num_sequences = len(all_sequences)
            result.timings["mpnn"] = time.time() - stage_start

            if not all_sequences:
                result.error = "LigandMPNN sequence design failed"
                return self._finalize_result(result, start_time)

            # Stage 6: Validate with ESMFold (+ RF3 if available)
            if validate_sequences and ESMFOLD_AVAILABLE:
                stage_start = time.time()
                # Pass ligand SMILES for RF3 ligand-aware prediction
                _ligand_smiles = ligand.smiles if (ligand and ligand.resolved and hasattr(ligand, 'smiles')) else None
                validation_results = self._validate_sequences(
                    all_sequences,
                    backbone_pdbs,
                    intent,
                    ligand_smiles=_ligand_smiles,
                )
                result.validation_results = validation_results
                result.timings["validation"] = time.time() - stage_start

                # Calculate pass rate
                passed = sum(1 for v in validation_results if v.passed)
                result.pass_rate = passed / len(validation_results) if validation_results else 0.0

                # Find best result
                self._find_best_result(result, validation_results, all_sequences)

            # Stage 7: Run unified analysis on best result
            if self.analyzer and result.best_sequence_pdb:
                stage_start = time.time()
                analysis = self._run_analysis(result.best_sequence_pdb, intent)
                result.analysis_results = analysis
                result.timings["analysis"] = time.time() - stage_start

            # Stage 8: Check for lessons
            if self.lesson_detector and result.analysis_results:
                lessons = self._check_lessons(result.analysis_results)
                result.lessons_triggered = lessons

            # Stage 9: Generate report
            result.report = self._generate_report(result)
            result.recommendations = self._generate_recommendations(result)

            # Mark success
            result.success = True

        except Exception as e:
            logger.exception(f"Pipeline error: {e}")
            result.error = str(e)
            result.success = False

        return self._finalize_result(result, start_time)

    def _parse_query(self, query: str) -> DesignIntent:
        """Parse natural language query into design intent."""
        logger.info(f"Parsing query: {query}")
        return self.parser.parse(query)

    def _resolve_ligand(self, intent: DesignIntent) -> ResolvedLigand:
        """Resolve ligand from design intent."""
        logger.info(f"Resolving ligand: {intent.ligand_name} with {intent.metal_type}")
        return self.resolver.resolve(
            intent.ligand_name or "",
            intent.metal_type,
        )

    def _generate_rfd3_config(
        self,
        intent: DesignIntent,
        ligand: ResolvedLigand,
        num_designs: int,
    ) -> RFD3Config:
        """Generate RFD3 configuration."""
        logger.info("Generating RFD3 configuration")
        return self.config_generator.generate(intent, ligand, num_designs)

    def _generate_mpnn_config(
        self,
        intent: DesignIntent,
        ligand: ResolvedLigand,
        num_sequences: int,
    ) -> LigandMPNNConfig:
        """Generate LigandMPNN configuration."""
        logger.info("Generating LigandMPNN configuration")
        config = self.config_generator.generate_mpnn_config(intent, ligand)
        config.num_sequences = num_sequences
        return config

    def _run_scaffolding_pipeline(
        self,
        intent: DesignIntent,
        result: PipelineResult,
        start_time: float,
        num_designs: int = 4,
        num_sequences: int = 8,
        validate_sequences: bool = True,
    ) -> PipelineResult:
        """
        Run scaffolding workflow for motif-based design.

        This is called when the query requests scaffolding from an existing PDB
        (e.g., "scaffold the PQQ-Ca pocket of 4CVB").

        Args:
            intent: Parsed design intent with source_pdb_id
            result: PipelineResult to populate
            start_time: Pipeline start time
            num_designs: Number of scaffold designs to generate
            num_sequences: Number of sequences per design
            validate_sequences: Whether to validate with ESMFold

        Returns:
            Populated PipelineResult
        """
        import asyncio

        try:
            # Stage 2: Fetch and extract scaffold
            stage_start = time.time()
            logger.info(f"Fetching scaffold from PDB: {intent.source_pdb_id}")

            # Run async scaffolding workflow
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)
            try:
                scaffold_result = loop.run_until_complete(
                    self.scaffolding_workflow.run(
                        pdb_id=intent.source_pdb_id,
                        metal=intent.metal_type,
                        ligand_code=intent.ligand_name,
                        include_all_ligand_contacts=intent.include_all_contacts,  # Full pocket
                        use_minimal_motif=True,          # Evidence-based residue selection
                        enzyme_class=intent.enzyme_class, # For enzyme chemistry evidence
                    )
                )
            finally:
                loop.close()

            result.timings["scaffold_fetch"] = time.time() - stage_start

            if not scaffold_result.success:
                result.error = f"Scaffolding failed: {scaffold_result.error_message}"
                return self._finalize_result(result, start_time)

            # Log pocket size
            pocket_mode = "full pocket" if intent.include_all_contacts else "metal coordination"
            logger.info(f"Scaffold extracted ({pocket_mode}): {len(scaffold_result.coordinating_residues)} residues")
            logger.info(f"Contig: {scaffold_result.contig}")

            # Stage 3: Build RFD3 config from scaffold
            stage_start = time.time()

            # === Stability Profile Selection (Phase 5: AI Infrastructure Enhancement) ===
            # Use decision engine's stability profile instead of hardcoded values
            stability = select_stability_profile(intent, intent.design_goal or "binding")
            logger.info(f"Selected stability profile: {stability.name} (cfg_scale={stability.cfg_scale}, step_scale={stability.step_scale})")

            rfd3_params = {
                "pdb_content": scaffold_result.motif_pdb,
                "contig": scaffold_result.contig,
                "num_designs": num_designs,
                "select_fixed_atoms": scaffold_result.fixed_atoms,
                "select_buried": scaffold_result.rasa_targets,
                "use_classifier_free_guidance": True,
                "cfg_scale": stability.cfg_scale,
                "step_scale": stability.step_scale,
                "num_timesteps": stability.num_timesteps,
            }

            # Add gamma_0 if stability profile specifies it
            if stability.gamma_0 is not None:
                rfd3_params["gamma_0"] = stability.gamma_0

            # Add H-bond conditioning if available
            if scaffold_result.hbond_acceptors:
                rfd3_params["select_hbond_acceptor"] = scaffold_result.hbond_acceptors
            if scaffold_result.hbond_donors:
                rfd3_params["select_hbond_donor"] = scaffold_result.hbond_donors

            # === Enzyme Activity Preservation (Phase 5: AI Infrastructure Enhancement) ===
            # Apply substrate channel exposure if enzyme class is detected and function preservation requested
            if ENZYME_CHEMISTRY_AVAILABLE and intent.enzyme_class and intent.preserve_function:
                logger.info(f"Enzyme activity preservation enabled: class={intent.enzyme_class}")

                # Get substrate channel atoms for RASA conditioning
                try:
                    exposed_atoms = get_substrate_channel_atoms(
                        pdb_content=scaffold_result.motif_pdb,
                        enzyme_class=intent.enzyme_class,
                        ligand_name=intent.ligand_name,
                    )
                    if exposed_atoms:
                        rfd3_params["select_exposed"] = exposed_atoms
                        logger.info(f"Substrate channel exposure: {exposed_atoms}")
                except Exception as e:
                    logger.warning(f"Could not get substrate channel atoms: {e}")

                # Get catalytic H-bond network preservation
                # Note: This requires identified catalytic residues from the scaffold
                # For now, rely on scaffold_result.hbond_acceptors/donors which are already added above
                # Future enhancement: use identify_catalytic_residues_from_pdb() with enzyme_class
                if hasattr(scaffold_result, 'coordinating_residues') and scaffold_result.coordinating_residues:
                    try:
                        # Convert coordinating residues to the format expected by get_catalytic_hbond_network
                        catalytic_residue_list = []
                        for res_str in scaffold_result.coordinating_residues:
                            # Parse residue strings like "A45" or "A 45"
                            import re
                            match = re.match(r'([A-Z])(\d+)', res_str.replace(' ', ''))
                            if match:
                                catalytic_residue_list.append({
                                    'chain': match.group(1),
                                    'resnum': int(match.group(2)),
                                })

                        if catalytic_residue_list:
                            catalytic_hbonds = get_catalytic_hbond_network(
                                scaffold_result.motif_pdb,
                                catalytic_residue_list,
                            )
                            if catalytic_hbonds.get("acceptors"):
                                existing = rfd3_params.get("select_hbond_acceptor", {})
                                for chain, atoms in catalytic_hbonds["acceptors"].items():
                                    if chain in existing:
                                        existing[chain] = existing[chain] + "," + atoms
                                    else:
                                        existing[chain] = atoms
                                rfd3_params["select_hbond_acceptor"] = existing
                            if catalytic_hbonds.get("donors"):
                                existing = rfd3_params.get("select_hbond_donor", {})
                                for chain, atoms in catalytic_hbonds["donors"].items():
                                    if chain in existing:
                                        existing[chain] = existing[chain] + "," + atoms
                                    else:
                                        existing[chain] = atoms
                                rfd3_params["select_hbond_donor"] = existing
                    except Exception as e:
                        logger.warning(f"Could not get catalytic H-bond network: {e}")

            # Log stability and enzyme settings
            logger.info(f"RFD3 params: cfg_scale={stability.cfg_scale}, step_scale={stability.step_scale}, "
                       f"num_timesteps={stability.num_timesteps}, gamma_0={stability.gamma_0}")

            result.timings["config"] = time.time() - stage_start

            # Stage 4: Run RFD3 backbone generation
            stage_start = time.time()
            logger.info(f"Running RFD3 scaffolding with contig: {scaffold_result.contig}")

            if self.use_mock:
                backbone_pdbs = [scaffold_result.motif_pdb] * num_designs
            else:
                rfd3_result = run_rfd3_inference(**rfd3_params)
                if rfd3_result.get("status") == "completed":
                    designs = rfd3_result.get("result", {}).get("designs", [])
                    backbone_pdbs = [d.get("content", "") for d in designs if d.get("content")]
                else:
                    backbone_pdbs = []

            result.backbone_pdbs = backbone_pdbs
            result.num_backbones = len(backbone_pdbs)
            result.timings["rfd3"] = time.time() - stage_start

            if not backbone_pdbs:
                result.error = "RFD3 scaffolding failed"
                return self._finalize_result(result, start_time)

            # Filter backbone continuity (chain breaks)
            from scaffolding_workflow import check_backbone_continuity
            continuous_pdbs = []
            for pdb in backbone_pdbs:
                cont = check_backbone_continuity(pdb)
                if cont["continuous"]:
                    continuous_pdbs.append(pdb)
                else:
                    logger.warning(f"Backbone has {cont['num_breaks']} chain breaks "
                                   f"(max CA-CA={cont['max_observed']}A), skipping")
            if continuous_pdbs:
                logger.info(f"Continuity filter: {len(continuous_pdbs)}/{len(backbone_pdbs)} passed")
                backbone_pdbs = continuous_pdbs
                result.backbone_pdbs = backbone_pdbs
                result.num_backbones = len(backbone_pdbs)
            else:
                logger.warning("No continuous backbones, using all designs as fallback")

            # Stage 5: Run LigandMPNN sequence design
            stage_start = time.time()
            all_sequences = []

            # Build MPNN config using decision engine (not hardcoded)
            mpnn_config = self._generate_mpnn_config(intent, ligand, num_sequences)
            mpnn_params = mpnn_config.to_api_params()
            mpnn_params.pop("task", None)  # Remove task key, run_mpnn_inference doesn't need it

            # Fix catalytic residues so MPNN does NOT redesign them
            catalytic_ids = scaffold_result.source_info.get("catalytic_residue_ids", [])
            if catalytic_ids:
                mpnn_params["fixed_positions"] = catalytic_ids
                mpnn_params["use_side_chain_context"] = True
                logger.info(f"MPNN: Fixing {len(catalytic_ids)} catalytic residues: {catalytic_ids[:5]}...")

            # Enable ligand-aware sidechain packing for enzyme designs
            if intent.ligand_name or intent.metal_type:
                mpnn_params["pack_side_chains"] = True
                mpnn_params["pack_with_ligand_context"] = True
                mpnn_params["number_of_packs_per_design"] = 4

            # Combined bias for metal + small molecule (override decision engine if both present)
            # NOTE: bias_AA is GLOBAL - keep mild to preserve hydrophobic core for folding.
            # LigandMPNN atom context already handles binding site design.
            if intent.metal_type and intent.ligand_name:
                mpnn_params["bias_AA"] = "H:0.5,D:0.5,E:0.5,W:0.3,Y:0.3"
                mpnn_params.pop("omit_AA", None)  # Allow cysteine for metal coordination

            logger.info(f"MPNN params: bias_AA={mpnn_params.get('bias_AA')}, "
                       f"fixed={len(catalytic_ids)} residues, "
                       f"pack_side_chains={mpnn_params.get('pack_side_chains', False)}")

            for backbone_pdb in backbone_pdbs:
                if self.use_mock:
                    all_sequences.append(SequenceResult(
                        sequence="MOCK" * 25,
                        confidence=0.9,
                    ))
                    continue

                mpnn_params["pdb_content"] = backbone_pdb
                mpnn_result = run_mpnn_inference(**mpnn_params)

                if mpnn_result.get("status") == "completed":
                    sequences_data = mpnn_result.get("result", {}).get("sequences", [])
                    for seq_data in sequences_data:
                        content = seq_data.get("content", "")
                        for line in content.split('\n'):
                            if line and not line.startswith('>'):
                                all_sequences.append(SequenceResult(
                                    sequence=line.strip(),
                                    confidence=seq_data.get("confidence", 0.8),
                                    filename=seq_data.get("filename", ""),
                                ))

            result.sequences = all_sequences
            result.num_sequences = len(all_sequences)
            result.timings["mpnn"] = time.time() - stage_start

            if not all_sequences:
                result.error = "LigandMPNN sequence design failed"
                return self._finalize_result(result, start_time)

            # Stage 6: Validate with ESMFold (+ RF3 if available)
            if validate_sequences and ESMFOLD_AVAILABLE:
                stage_start = time.time()
                # Pass ligand SMILES for RF3 ligand-aware prediction (scaffolding path)
                _ligand_smiles = intent.ligand_smiles if hasattr(intent, 'ligand_smiles') else None
                validation_results = self._validate_sequences(
                    all_sequences,
                    backbone_pdbs,
                    intent,
                    ligand_smiles=_ligand_smiles,
                )
                result.validation_results = validation_results
                result.timings["validation"] = time.time() - stage_start

                passed = sum(1 for v in validation_results if v.passed)
                result.pass_rate = passed / len(validation_results) if validation_results else 0.0

                self._find_best_result(result, validation_results, all_sequences)

            # Stage 7: Run unified analysis
            if self.analyzer and result.best_sequence_pdb:
                stage_start = time.time()
                analysis = self._run_analysis(result.best_sequence_pdb, intent)
                result.analysis_results = analysis
                result.timings["analysis"] = time.time() - stage_start

            # Stage 8: Generate report
            result.report = self._generate_scaffolding_report(result, scaffold_result)
            result.recommendations = self._generate_recommendations(result)

            # Add scaffolding info to analysis
            result.analysis_results["scaffolding_info"] = scaffold_result.source_info

            # Add stability profile info (Phase 5: for post-design filtering)
            result.analysis_results["stability_profile"] = {
                "name": stability.name,
                "plddt_threshold": stability.plddt_threshold,
                "cfg_scale": stability.cfg_scale,
                "step_scale": stability.step_scale,
                "gamma_0": stability.gamma_0,
                "num_timesteps": stability.num_timesteps,
            }

            # Add enzyme preservation info if applicable
            if intent.enzyme_class:
                result.analysis_results["enzyme_preservation"] = {
                    "enzyme_class": intent.enzyme_class,
                    "preserve_function": intent.preserve_function,
                    "enzyme_class_confidence": getattr(intent, 'enzyme_class_confidence', 0.0),
                }

            result.success = True

        except Exception as e:
            logger.exception(f"Scaffolding pipeline error: {e}")
            result.error = str(e)
            result.success = False

        return self._finalize_result(result, start_time)

    def _generate_scaffolding_report(
        self,
        result: PipelineResult,
        scaffold_result: 'ScaffoldResult',
    ) -> str:
        """Generate report for scaffolding workflow."""
        lines = ["# AI Design Pipeline Report (Scaffolding Mode)\n"]

        # Scaffold source
        lines.append("## Scaffold Source")
        info = scaffold_result.source_info
        lines.append(f"- PDB ID: {info.get('pdb_id', 'N/A')}")
        lines.append(f"- Metal: {info.get('metal', 'N/A')}")
        lines.append(f"- Ligand: {info.get('ligand', 'N/A')}")
        lines.append(f"- Coordination number: {info.get('coordination_number', 0)}")

        if info.get('protein_donors'):
            lines.append(f"- Protein donors: {', '.join(info['protein_donors'][:5])}")
        if info.get('ligand_donors'):
            lines.append(f"- Ligand donors: {', '.join(info['ligand_donors'][:5])}")
        lines.append("")

        # Contig
        lines.append("## Scaffolding Configuration")
        lines.append(f"- Contig: `{scaffold_result.contig}`")
        lines.append(f"- Motif residues: {', '.join(scaffold_result.motif_residues)}")
        lines.append("")

        # Results (same as standard report)
        lines.append("## Results Summary")
        lines.append(f"- Backbones generated: {result.num_backbones}")
        lines.append(f"- Sequences designed: {result.num_sequences}")

        if result.validation_results:
            lines.append(f"- Validation pass rate: {result.pass_rate:.1%}")

        if result.best_rmsd and result.best_rmsd < float('inf'):
            lines.append(f"- Best RMSD: {result.best_rmsd:.2f}Å")
            lines.append(f"- Best pLDDT: {result.best_plddt:.2f}")
        lines.append("")

        # Best sequence
        if result.best_sequence:
            lines.append("## Best Sequence")
            lines.append("```")
            lines.append(result.best_sequence.sequence)
            lines.append("```")
            lines.append("")

        # Timing
        lines.append("## Timing")
        for stage, duration in result.timings.items():
            lines.append(f"- {stage}: {duration:.1f}s")
        lines.append(f"- Total: {result.total_time:.1f}s")

        return "\n".join(lines)

    def _run_rfd3(self, config: RFD3Config) -> List[str]:
        """Run RFD3 backbone generation."""
        logger.info(f"Running RFD3 with {config.num_designs} designs")

        try:
            params = config.to_api_params()
            # Remove 'task' key - that's for the handler router, not the inference function
            params.pop("task", None)

            if self.use_mock:
                # Return mock PDBs for testing
                return [config.pdb_content or "MOCK PDB"] * config.num_designs

            result = run_rfd3_inference(**params)

            if result.get("status") == "completed":
                designs = result.get("result", {}).get("designs", [])
                return [d.get("content", "") for d in designs if d.get("content")]

        except Exception as e:
            logger.error(f"RFD3 inference failed: {e}")

        return []

    def _run_mpnn(
        self,
        backbone_pdb: str,
        config: LigandMPNNConfig,
    ) -> List[SequenceResult]:
        """Run LigandMPNN sequence design."""
        logger.info(f"Running LigandMPNN with {config.num_sequences} sequences")

        try:
            params = config.to_api_params()
            # Remove 'task' key - that's for the handler router, not the inference function
            params.pop("task", None)
            params["pdb_content"] = backbone_pdb

            if self.use_mock:
                # Return mock sequences for testing
                return [
                    SequenceResult(sequence="MOCK" * 25, confidence=0.9)
                    for _ in range(config.num_sequences)
                ]

            result = run_mpnn_inference(**params)

            if result.get("status") == "completed":
                sequences_data = result.get("result", {}).get("sequences", [])
                sequences = []

                for seq_data in sequences_data:
                    content = seq_data.get("content", "")
                    # Parse FASTA content
                    for line in content.split('\n'):
                        if line and not line.startswith('>'):
                            sequences.append(SequenceResult(
                                sequence=line.strip(),
                                confidence=seq_data.get("confidence", 0.8),
                                filename=seq_data.get("filename", ""),
                            ))

                return sequences

        except Exception as e:
            logger.error(f"MPNN inference failed: {e}")

        return []

    def _detect_protein_chain(self, pdb_content: str) -> str:
        """
        Detect the protein chain from backbone PDB content.

        For metal-ligand designs, the backbone typically has:
        - Metal atoms (HETATM with element like TB, ZN)
        - Ligand atoms (HETATM with residue like CIT)
        - Protein residues (ATOM with standard amino acid names)

        Returns the chain ID containing the most ATOM records.
        """
        if not pdb_content:
            return "A"  # Default fallback

        chain_counts = {}
        for line in pdb_content.split('\n'):
            if line.startswith('ATOM'):
                chain = line[21] if len(line) > 21 else 'A'
                chain_counts[chain] = chain_counts.get(chain, 0) + 1

        if chain_counts:
            # Return chain with most ATOM records (protein chain)
            return max(chain_counts, key=chain_counts.get)

        # If no ATOM records, look for any chain in HETATM (all atoms may be HETATM)
        for line in pdb_content.split('\n'):
            if line.startswith('HETATM'):
                chain = line[21] if len(line) > 21 else 'A'
                return chain

        return "A"  # Default fallback

    def _validate_sequences(
        self,
        sequences: List[SequenceResult],
        backbone_pdbs: List[str],
        intent: DesignIntent,
        validation_threshold: str = "exploratory",
        ligand_smiles: Optional[str] = None,
    ) -> List[ValidationResult]:
        """Validate sequences with ESMFold and optionally RF3.

        Args:
            sequences: List of sequences to validate
            backbone_pdbs: List of backbone PDBs to compare against
            intent: Design intent for context
            validation_threshold: Threshold level ("strict", "standard", "relaxed", "exploratory")
                Default is "exploratory" for early-stage metal-ligand designs where
                high RMSD is expected due to the challenging nature of de novo design.
            ligand_smiles: Optional SMILES for ligand-aware RF3 prediction
        """
        methods = ["ESMFold"]
        if RF3_AVAILABLE:
            methods.append("RF3")
        logger.info(f"Validating {len(sequences)} sequences with {'+'.join(methods)} (threshold={validation_threshold})")

        results = []

        # Detect protein chain from backbone PDB
        protein_chain = self._detect_protein_chain(backbone_pdbs[0] if backbone_pdbs else "")
        logger.info(f"Detected protein chain: {protein_chain}")

        for seq in sequences[:8]:  # Limit to 8 for speed
            try:
                if self.use_mock:
                    results.append(ValidationResult(
                        passed=True,
                        rmsd=1.5,
                        plddt=0.85,
                        best_rmsd=1.5,
                        best_plddt=0.85,
                    ))
                    continue

                # --- ESMFold validation ---
                validation = validate_structure_esmfold(
                    seq.sequence,
                    backbone_pdbs[0] if backbone_pdbs else None,
                    binder_chain=protein_chain,
                    threshold=validation_threshold,
                )

                esm_rmsd = validation.get("esmfold_rmsd")
                esm_plddt = validation.get("esmfold_plddt", 0.0)
                esm_passed = validation.get("validation_passed", False)

                logger.info(f"ESMFold: pLDDT={esm_plddt:.2f}, RMSD={esm_rmsd if esm_rmsd else 'N/A'}Å, passed={esm_passed}")

                val_result = ValidationResult(
                    passed=esm_passed,
                    rmsd=esm_rmsd if esm_rmsd is not None else float('inf'),
                    plddt=esm_plddt,
                    predicted_pdb=validation.get("predicted_pdb"),
                )

                # --- RF3 validation (if available) ---
                if RF3_AVAILABLE and validate_structure_rf3 is not None:
                    try:
                        rf3_val = validate_structure_rf3(
                            sequence=seq.sequence,
                            designed_pdb=backbone_pdbs[0] if backbone_pdbs else None,
                            binder_chain=protein_chain,
                            ligand_smiles=ligand_smiles,
                            threshold=validation_threshold,
                        )

                        if rf3_val.get("status") == "completed":
                            val_result.rf3_plddt = rf3_val.get("rf3_plddt")
                            val_result.rf3_ptm = rf3_val.get("rf3_ptm")
                            val_result.rf3_rmsd = rf3_val.get("rf3_rmsd")

                            rf3_plddt_str = f"{val_result.rf3_plddt:.2f}" if val_result.rf3_plddt else "N/A"
                            rf3_rmsd_str = f"{val_result.rf3_rmsd:.2f}" if val_result.rf3_rmsd else "N/A"
                            logger.info(f"RF3: pLDDT={rf3_plddt_str}, RMSD={rf3_rmsd_str}Å")
                        else:
                            logger.warning(f"RF3 validation failed: {rf3_val.get('error')}")

                    except Exception as e:
                        logger.error(f"RF3 validation error: {e}")

                # --- Compute best-of metrics ---
                plddt_values = [v for v in [esm_plddt, val_result.rf3_plddt] if v is not None and v > 0]
                rmsd_values = [v for v in [esm_rmsd, val_result.rf3_rmsd] if v is not None]
                val_result.best_plddt = max(plddt_values) if plddt_values else esm_plddt
                val_result.best_rmsd = min(rmsd_values) if rmsd_values else (esm_rmsd if esm_rmsd is not None else float('inf'))

                # Re-evaluate pass/fail using best-of metrics
                from esmfold_utils import ESMFOLD_THRESHOLDS
                thresholds = ESMFOLD_THRESHOLDS.get(validation_threshold, ESMFOLD_THRESHOLDS["standard"])
                best_passes_plddt = val_result.best_plddt is not None and val_result.best_plddt >= thresholds["plddt_min"]
                best_passes_rmsd = val_result.best_rmsd is None or val_result.best_rmsd <= thresholds["rmsd_max"]
                val_result.passed = best_passes_plddt and best_passes_rmsd

                results.append(val_result)

            except Exception as e:
                logger.error(f"Validation failed for sequence: {e}")
                results.append(ValidationResult(
                    passed=False,
                    warnings=[str(e)],
                ))

        return results

    def _find_best_result(
        self,
        result: PipelineResult,
        validations: List[ValidationResult],
        sequences: List[SequenceResult],
    ) -> None:
        """Find the best result from validation.

        First looks for passing designs (sorted by RMSD).
        If none pass, still reports the best non-passing result (highest pLDDT).
        """
        best_passing_idx = -1
        best_passing_rmsd = float('inf')

        best_overall_idx = -1
        best_overall_plddt = 0.0

        for i, val in enumerate(validations):
            # Use best_rmsd/best_plddt if available (best-of across ESMFold + RF3)
            effective_rmsd = val.best_rmsd if val.best_rmsd is not None else val.rmsd
            effective_plddt = val.best_plddt if val.best_plddt is not None else val.plddt

            # Track best passing design (by RMSD)
            if val.passed and effective_rmsd < best_passing_rmsd:
                best_passing_rmsd = effective_rmsd
                best_passing_idx = i

            # Track best overall design (by pLDDT) - for reporting even if none pass
            if effective_plddt > best_overall_plddt:
                best_overall_plddt = effective_plddt
                best_overall_idx = i

        # Use passing design if available, otherwise use best overall
        best_idx = best_passing_idx if best_passing_idx >= 0 else best_overall_idx

        if best_idx >= 0 and best_idx < len(sequences):
            result.best_sequence = sequences[best_idx]
            best_val = validations[best_idx]
            effective_rmsd = best_val.best_rmsd if best_val.best_rmsd is not None else best_val.rmsd
            effective_plddt = best_val.best_plddt if best_val.best_plddt is not None else best_val.plddt
            result.best_rmsd = effective_rmsd if effective_rmsd != float('inf') else None
            result.best_plddt = effective_plddt
            result.best_sequence_pdb = best_val.predicted_pdb

            passed_str = "passed" if best_val.passed else "not passed"
            rmsd_str = f"{result.best_rmsd:.2f}Å" if result.best_rmsd else "N/A"
            logger.info(f"Best result: idx={best_idx}, pLDDT={result.best_plddt:.2f}, RMSD={rmsd_str}, {passed_str}")

    def _run_analysis(
        self,
        pdb_content: str,
        intent: DesignIntent,
    ) -> Dict[str, Any]:
        """Run unified analysis on the best result."""
        logger.info("Running unified analysis")

        if not self.analyzer:
            return {}

        try:
            results = self.analyzer.analyze(pdb_content)
            return results.to_dict() if hasattr(results, 'to_dict') else {}
        except Exception as e:
            logger.error(f"Analysis failed: {e}")
            return {}

    def _check_lessons(self, analysis: Dict[str, Any]) -> List[str]:
        """Check if any lessons are triggered."""
        if not self.lesson_detector:
            return []

        try:
            return self.lesson_detector.check_triggers(analysis)
        except Exception as e:
            logger.error(f"Lesson detection failed: {e}")
            return []

    def _generate_report(self, result: PipelineResult) -> str:
        """Generate human-readable report."""
        lines = ["# AI Design Pipeline Report\n"]

        # Intent
        if result.design_intent:
            intent = result.design_intent
            lines.append("## Design Intent")
            lines.append(f"- Metal: {intent.metal_type or 'None'}")
            lines.append(f"- Ligand: {intent.ligand_name or 'None'}")
            lines.append(f"- Goal: {intent.design_goal}")
            lines.append(f"- Topology: {intent.target_topology}")
            lines.append(f"- Chain length: {intent.chain_length_range}")
            lines.append("")

        # Results summary
        lines.append("## Results Summary")
        lines.append(f"- Backbones generated: {result.num_backbones}")
        lines.append(f"- Sequences designed: {result.num_sequences}")

        if result.validation_results:
            lines.append(f"- Validation pass rate: {result.pass_rate:.1%}")

        if result.best_rmsd < float('inf'):
            lines.append(f"- Best RMSD: {result.best_rmsd:.2f}Å")
            lines.append(f"- Best pLDDT: {result.best_plddt:.2f}")

        lines.append("")

        # Best sequence
        if result.best_sequence:
            lines.append("## Best Sequence")
            lines.append(f"```")
            lines.append(result.best_sequence.sequence)
            lines.append(f"```")
            lines.append("")

        # Timing
        lines.append("## Timing")
        for stage, duration in result.timings.items():
            lines.append(f"- {stage}: {duration:.1f}s")
        lines.append(f"- Total: {result.total_time:.1f}s")

        return "\n".join(lines)

    def _generate_recommendations(self, result: PipelineResult) -> List[str]:
        """Generate actionable recommendations."""
        recommendations = []

        # Based on pass rate
        if result.pass_rate < 0.2:
            recommendations.append(
                "Low pass rate - consider adjusting chain length or conditioning"
            )
        elif result.pass_rate > 0.8:
            recommendations.append(
                "High pass rate - design strategy is working well"
            )

        # Based on RMSD
        if result.best_rmsd < 1.5:
            recommendations.append(
                "Excellent fold prediction - proceed to experimental validation"
            )
        elif result.best_rmsd < 2.0:
            recommendations.append(
                "Good fold prediction - consider AF3 validation with ligand"
            )
        elif result.best_rmsd < float('inf'):
            recommendations.append(
                "Moderate RMSD - may need additional refinement"
            )

        # From intent suggestions
        if result.design_intent:
            recommendations.extend(result.design_intent.suggestions[:3])

        # From lessons
        if result.lessons_triggered:
            recommendations.append(
                f"Triggered lessons: {', '.join(result.lessons_triggered)}"
            )

        return recommendations

    def _finalize_result(
        self,
        result: PipelineResult,
        start_time: float,
    ) -> PipelineResult:
        """Finalize result with timing."""
        result.total_time = time.time() - start_time
        return result


# =============================================================================
# Handler Function for API
# =============================================================================

def handle_ai_design(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Handler function for AI design task.

    Input:
        query: str - Natural language design request
        claude_api_key: str - API key for NL parsing
        num_designs: int (optional, default 4)
        num_sequences: int (optional, default 8)
        session_name: str (optional) - For history tracking
        validate: bool (optional, default True) - Run ESMFold validation

    Output:
        status: "completed" | "failed"
        result: PipelineResult as dict
    """
    try:
        # Extract parameters
        query = job_input.get("query")
        if not query:
            return {"status": "failed", "error": "Missing 'query' parameter"}

        claude_api_key = job_input.get("claude_api_key")
        claude_api_url = job_input.get("claude_api_url", "https://yinli.one/v1/messages")
        num_designs = job_input.get("num_designs", 4)
        num_sequences = job_input.get("num_sequences", 8)
        session_name = job_input.get("session_name")
        validate = job_input.get("validate", True)
        use_mock = job_input.get("use_mock", False)

        # Initialize pipeline
        pipeline = AIDesignPipeline(
            claude_api_key=claude_api_key,
            claude_api_url=claude_api_url,
            use_mock=use_mock,
        )

        # Run pipeline
        result = pipeline.run(
            query=query,
            num_designs=num_designs,
            num_sequences=num_sequences,
            validate_sequences=validate,
            session_name=session_name,
        )

        return {
            "status": "completed" if result.success else "failed",
            "result": result.to_dict(),
        }

    except Exception as e:
        logger.exception(f"AI design handler error: {e}")
        return {
            "status": "failed",
            "error": str(e),
        }


# =============================================================================
# Test Function
# =============================================================================

def test_pipeline():
    """Test the pipeline with mock inference."""
    pipeline = AIDesignPipeline(
        claude_api_key=None,  # Use fallback parser
        use_mock=True,
    )

    # Test 1: De novo design query
    print("\n" + "=" * 60)
    print("Test 1: De Novo Design")
    print("=" * 60)

    result = pipeline.run(
        query="Design a protein to bind citrate with terbium",
        num_designs=2,
        num_sequences=4,
        validate_sequences=True,
    )

    print(f"Success: {result.success}")
    print(f"Design Mode: {result.design_intent.design_mode if result.design_intent else 'N/A'}")
    print(f"Error: {result.error}")
    print(f"Backbones: {result.num_backbones}")
    print(f"Sequences: {result.num_sequences}")
    print(f"Total time: {result.total_time:.1f}s")

    # Test 2: Scaffolding query
    print("\n" + "=" * 60)
    print("Test 2: Scaffolding Query")
    print("=" * 60)

    result2 = pipeline.run(
        query="Scaffold the PQQ-Ca pocket of 4CVB",
        num_designs=2,
        num_sequences=4,
        validate_sequences=False,  # Skip validation for speed
    )

    print(f"Success: {result2.success}")
    print(f"Design Mode: {result2.design_intent.design_mode if result2.design_intent else 'N/A'}")
    print(f"Source PDB: {result2.design_intent.source_pdb_id if result2.design_intent else 'N/A'}")
    print(f"Error: {result2.error}")
    print(f"Total time: {result2.total_time:.1f}s")


if __name__ == "__main__":
    test_pipeline()
