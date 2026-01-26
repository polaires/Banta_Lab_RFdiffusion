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

# Import inference utilities
from inference_utils import run_rfd3_inference, run_mpnn_inference

# Import ESMFold validation
try:
    from esmfold_utils import validate_structure_esmfold, is_esmfold_available
    ESMFOLD_AVAILABLE = is_esmfold_available()
except ImportError:
    ESMFOLD_AVAILABLE = False
    validate_structure_esmfold = None

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

            # Stage 6: Validate with ESMFold
            if validate_sequences and ESMFOLD_AVAILABLE:
                stage_start = time.time()
                validation_results = self._validate_sequences(
                    all_sequences,
                    backbone_pdbs,
                    intent,
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
    ) -> List[ValidationResult]:
        """Validate sequences with ESMFold.

        Args:
            sequences: List of sequences to validate
            backbone_pdbs: List of backbone PDBs to compare against
            intent: Design intent for context
            validation_threshold: Threshold level ("strict", "standard", "relaxed", "exploratory")
                Default is "exploratory" for early-stage metal-ligand designs where
                high RMSD is expected due to the challenging nature of de novo design.
        """
        logger.info(f"Validating {len(sequences)} sequences with ESMFold (threshold={validation_threshold})")

        results = []

        # Detect protein chain from backbone PDB
        # For metal-ligand designs, protein is often on chain A or B
        protein_chain = self._detect_protein_chain(backbone_pdbs[0] if backbone_pdbs else "")
        logger.info(f"Detected protein chain: {protein_chain}")

        for seq in sequences[:8]:  # Limit to 8 for speed
            try:
                if self.use_mock:
                    results.append(ValidationResult(
                        passed=True,
                        rmsd=1.5,
                        plddt=0.85,
                    ))
                    continue

                validation = validate_structure_esmfold(
                    seq.sequence,
                    backbone_pdbs[0] if backbone_pdbs else None,
                    binder_chain=protein_chain,
                    threshold=validation_threshold,
                )

                # Note: validate_structure_esmfold returns "esmfold_plddt" and "esmfold_rmsd"
                rmsd = validation.get("esmfold_rmsd")  # Can be None if calculation failed
                plddt = validation.get("esmfold_plddt", 0.0)

                # Use validation_passed from the function directly
                passed = validation.get("validation_passed", False)

                logger.info(f"Seq validation: pLDDT={plddt:.2f}, RMSD={rmsd if rmsd else 'N/A'}Å, passed={passed}")

                results.append(ValidationResult(
                    passed=passed,
                    rmsd=rmsd if rmsd is not None else float('inf'),
                    plddt=plddt,
                    predicted_pdb=validation.get("predicted_pdb"),
                ))

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
            # Track best passing design (by RMSD)
            if val.passed and val.rmsd < best_passing_rmsd:
                best_passing_rmsd = val.rmsd
                best_passing_idx = i

            # Track best overall design (by pLDDT) - for reporting even if none pass
            if val.plddt > best_overall_plddt:
                best_overall_plddt = val.plddt
                best_overall_idx = i

        # Use passing design if available, otherwise use best overall
        best_idx = best_passing_idx if best_passing_idx >= 0 else best_overall_idx

        if best_idx >= 0 and best_idx < len(sequences):
            result.best_sequence = sequences[best_idx]
            result.best_rmsd = validations[best_idx].rmsd if validations[best_idx].rmsd != float('inf') else None
            result.best_plddt = validations[best_idx].plddt
            result.best_sequence_pdb = validations[best_idx].predicted_pdb

            passed_str = "passed" if validations[best_idx].passed else "not passed"
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

    result = pipeline.run(
        query="Design a protein to bind citrate with terbium",
        num_designs=2,
        num_sequences=4,
        validate_sequences=True,
    )

    print("\n" + "=" * 60)
    print("Pipeline Test Results")
    print("=" * 60)
    print(f"Success: {result.success}")
    print(f"Error: {result.error}")
    print(f"Backbones: {result.num_backbones}")
    print(f"Sequences: {result.num_sequences}")
    print(f"Pass rate: {result.pass_rate:.1%}")
    print(f"Total time: {result.total_time:.1f}s")
    print("\n" + result.report)


if __name__ == "__main__":
    test_pipeline()
