"""
Conservation Analyzer - ConSurf Methodology Implementation

Calculates evolutionary conservation scores using Rate4Site algorithm.
Implements ConSurf's MSA preparation pipeline with lightweight alternatives:
- NCBI BLAST API (instead of local SwissProt)
- MMseqs2 (instead of CD-HIT)
- MUSCLE (alignment)
- Rate4Site (conservation calculation)

ConSurf Parameters Preserved:
- BLAST_MAX_HOMOLOGUES_TO_DISPLAY = 500
- FRAGMENT_REDUNDANCY_RATE = 95% (MMseqs2 --min-seq-id 0.95)
- MINIMUM_FRAGMENTS_FOR_MSA = 5
- BAYES_INTERVAL = 3 (for confidence intervals)

Usage:
    analyzer = ConservationAnalyzer()
    result = await analyzer.analyze(pdb_content, chain="A")

    # Get highly conserved positions for scaffolding
    for pos in result.highly_conserved:
        print(f"Position {pos}: highly conserved (grade 1-3)")

Output matches ConSurf GradesPE format:
    Grade 1-3: Highly conserved (maroon in ConSurf)
    Grade 4-6: Average conservation
    Grade 7-9: Variable (turquoise in ConSurf)
"""

import asyncio
import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any, Tuple
from datetime import datetime

logger = logging.getLogger(__name__)


# Import utilities
try:
    from utils.blast_client import (
        NCBIBlastClient,
        BlastResult,
        extract_sequence_from_pdb,
    )
    BLAST_CLIENT_AVAILABLE = True
except ImportError:
    BLAST_CLIENT_AVAILABLE = False

try:
    from utils.msa_utils import (
        prepare_msa_from_blast,
        filter_blast_hits,
        cluster_sequences_mmseqs,
        run_muscle_alignment,
        MSAResult,
    )
    MSA_UTILS_AVAILABLE = True
except ImportError:
    MSA_UTILS_AVAILABLE = False

try:
    from utils.rate4site_wrapper import (
        run_rate4site,
        Rate4SiteResult,
        Rate4SiteScore,
        assign_consurf_grades,
        calculate_simple_conservation,
        scores_to_consurf_format,
        create_conservation_pdb,
        RATE4SITE_AVAILABLE,
    )
except ImportError:
    RATE4SITE_AVAILABLE = False


@dataclass
class ConservationGrade:
    """
    Single residue conservation result matching ConSurf format.

    Attributes:
        position: 1-indexed position in sequence
        residue: Single-letter amino acid code
        score: Rate4Site score (lower = more conserved)
        grade: ConSurf grade 1-9 (1 = most conserved)
        confidence_interval: Bayesian confidence interval (low, high)
        data_quality: "sufficient" or "insufficient" based on MSA depth
    """
    position: int
    residue: str
    score: float
    grade: int
    confidence_interval: Tuple[float, float]
    data_quality: str

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "position": self.position,
            "residue": self.residue,
            "score": self.score,
            "grade": self.grade,
            "confidence_interval": list(self.confidence_interval),
            "data_quality": self.data_quality,
        }


@dataclass
class ConservationAnalysis:
    """
    Complete conservation analysis result.

    Provides ConSurf-compatible conservation data for protein design:
    - Per-residue grades (1-9 scale)
    - Identification of conserved/variable regions
    - MSA quality metrics

    For scaffolding/partial diffusion:
    - Use highly_conserved positions as fixed residues
    - Allow variable positions to be redesigned
    """
    grades: List[ConservationGrade]
    query_sequence: str
    query_name: str
    msa_depth: int              # Number of sequences in alignment
    method: str                 # "bayesian", "ml", or "entropy"
    highly_conserved: List[int]    # Positions with grade 1-3
    conserved: List[int]           # Positions with grade 4-5
    variable: List[int]            # Positions with grade 7-9
    average_conservation: float    # Mean conservation grade
    notes: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "grades": [g.to_dict() for g in self.grades],
            "query_sequence": self.query_sequence,
            "query_name": self.query_name,
            "msa_depth": self.msa_depth,
            "method": self.method,
            "highly_conserved": self.highly_conserved,
            "conserved": self.conserved,
            "variable": self.variable,
            "average_conservation": self.average_conservation,
            "notes": self.notes,
        }

    def get_grade(self, position: int) -> Optional[int]:
        """Get conservation grade for a specific position."""
        for g in self.grades:
            if g.position == position:
                return g.grade
        return None

    def get_grades_dict(self) -> Dict[int, int]:
        """Get position -> grade mapping."""
        return {g.position: g.grade for g in self.grades}

    def get_fixed_residue_suggestions(
        self,
        threshold: int = 3,
    ) -> List[str]:
        """
        Get residue IDs suggested for fixing during design.

        Args:
            threshold: Maximum grade to consider "highly conserved"

        Returns:
            List of residue IDs in format "A{position}"
        """
        fixed = []
        for g in self.grades:
            if g.grade <= threshold and g.data_quality == "sufficient":
                fixed.append(f"A{g.position}")
        return fixed


class ConservationAnalyzer:
    """
    Main analyzer class implementing ConSurf methodology.

    Performs evolutionary conservation analysis by:
    1. Extracting sequence from PDB
    2. Running BLAST to find homologs
    3. Filtering and clustering results
    4. Creating multiple sequence alignment
    5. Calculating conservation scores with Rate4Site
    6. Assigning ConSurf grades (1-9)

    ConSurf Parameters (preserved exactly):
    - BLAST_MAX_HOMOLOGS = 500
    - REDUNDANCY_CUTOFF = 0.95 (95% identity clustering)
    - MIN_SEQUENCES_FOR_MSA = 5
    - BAYES_INTERVAL = 3
    """

    # ConSurf parameters
    BLAST_MAX_HOMOLOGS = 500
    REDUNDANCY_CUTOFF = 0.95
    MIN_SEQUENCES_FOR_MSA = 5
    BAYES_INTERVAL = 3

    # Grade thresholds
    HIGHLY_CONSERVED_THRESHOLD = 3  # Grades 1-3
    VARIABLE_THRESHOLD = 7           # Grades 7-9

    def __init__(
        self,
        blast_email: Optional[str] = None,
        use_cache: bool = True,
    ):
        """
        Initialize conservation analyzer.

        Args:
            blast_email: Email for NCBI BLAST API (recommended for heavy usage)
            use_cache: Whether to cache BLAST results
        """
        self.blast_email = blast_email
        self.use_cache = use_cache

        if BLAST_CLIENT_AVAILABLE:
            self.blast_client = NCBIBlastClient(email=blast_email)
        else:
            self.blast_client = None

    async def analyze(
        self,
        pdb_content: str,
        chain: str = "A",
        method: str = "bayesian",
        precomputed_msa: Optional[str] = None,
        skip_blast: bool = False,
    ) -> ConservationAnalysis:
        """
        Analyze evolutionary conservation of a protein structure.

        Args:
            pdb_content: PDB file content
            chain: Chain ID to analyze
            method: Conservation calculation method ("bayesian" or "ml")
            precomputed_msa: Optional pre-computed MSA in FASTA format
            skip_blast: Skip BLAST search (requires precomputed_msa)

        Returns:
            ConservationAnalysis with per-residue grades

        Raises:
            RuntimeError: If analysis fails
            ValueError: If insufficient homologs found
        """
        notes = {}
        start_time = datetime.now()

        # Step 1: Extract sequence from PDB
        logger.info(f"Extracting sequence from chain {chain}")
        query_sequence = extract_sequence_from_pdb(pdb_content, chain)

        if not query_sequence:
            raise ValueError(f"Could not extract sequence from chain {chain}")

        query_name = f"query_{chain}"
        notes["query_length"] = len(query_sequence)
        logger.info(f"Query sequence: {len(query_sequence)} residues")

        # Step 2: Get MSA (either provided or via BLAST)
        if precomputed_msa:
            logger.info("Using precomputed MSA")
            msa_fasta = precomputed_msa
            msa_depth = precomputed_msa.count(">")
            notes["msa_source"] = "precomputed"

        elif skip_blast:
            raise ValueError("skip_blast=True requires precomputed_msa")

        else:
            # Run BLAST pipeline
            msa_fasta, msa_depth, blast_notes = await self._run_blast_pipeline(
                query_sequence,
                query_name,
            )
            notes.update(blast_notes)

        # Step 3: Calculate conservation scores
        if RATE4SITE_AVAILABLE:
            logger.info(f"Running Rate4Site ({method} method)")
            try:
                r4s_result = run_rate4site(
                    msa_fasta,
                    method=method,
                    query_name=query_name,
                )
                scores = r4s_result.scores
                notes["rate4site_method"] = method
                notes["rate4site_model"] = r4s_result.model

            except Exception as e:
                logger.warning(f"Rate4Site failed: {e}, falling back to entropy")
                scores = calculate_simple_conservation(msa_fasta)
                method = "entropy"
                notes["rate4site_error"] = str(e)

        else:
            logger.info("Rate4Site not available, using Shannon entropy")
            scores = calculate_simple_conservation(msa_fasta)
            method = "entropy"
            notes["rate4site_available"] = False

        # Step 4: Assign ConSurf grades
        logger.info("Assigning ConSurf grades (1-9)")
        grades = assign_consurf_grades(scores)

        # Step 5: Build result
        conservation_grades = []
        highly_conserved = []
        conserved = []
        variable = []

        for score, grade in zip(scores, grades):
            cg = ConservationGrade(
                position=score.position,
                residue=score.residue,
                score=score.score,
                grade=grade,
                confidence_interval=score.qq_interval,
                data_quality=score.data_quality,
            )
            conservation_grades.append(cg)

            if grade <= self.HIGHLY_CONSERVED_THRESHOLD:
                highly_conserved.append(score.position)
            elif grade <= 5:
                conserved.append(score.position)
            elif grade >= self.VARIABLE_THRESHOLD:
                variable.append(score.position)

        # Calculate average conservation
        avg_grade = sum(g.grade for g in conservation_grades) / len(conservation_grades) if conservation_grades else 5.0

        # Timing
        elapsed = (datetime.now() - start_time).total_seconds()
        notes["analysis_time_seconds"] = elapsed
        logger.info(f"Conservation analysis completed in {elapsed:.1f}s")

        return ConservationAnalysis(
            grades=conservation_grades,
            query_sequence=query_sequence,
            query_name=query_name,
            msa_depth=msa_depth,
            method=method,
            highly_conserved=highly_conserved,
            conserved=conserved,
            variable=variable,
            average_conservation=avg_grade,
            notes=notes,
        )

    async def _run_blast_pipeline(
        self,
        query_sequence: str,
        query_name: str,
    ) -> Tuple[str, int, Dict[str, Any]]:
        """
        Run BLAST search and prepare MSA.

        Returns:
            Tuple of (MSA FASTA, sequence count, notes dict)
        """
        notes = {}

        if not BLAST_CLIENT_AVAILABLE:
            raise RuntimeError(
                "BLAST client not available. Install aiohttp or provide precomputed_msa."
            )

        # Run BLAST search
        logger.info(f"Running NCBI BLAST search (max {self.BLAST_MAX_HOMOLOGS} hits)")
        blast_result = await self.blast_client.search(
            query_sequence,
            database="nr",
            max_hits=self.BLAST_MAX_HOMOLOGS,
            use_cache=self.use_cache,
        )

        notes["blast_hits"] = len(blast_result.hits)
        notes["blast_rid"] = blast_result.rid
        notes["blast_cached"] = blast_result.cached
        logger.info(f"BLAST found {len(blast_result.hits)} hits (cached={blast_result.cached})")

        # Prepare MSA
        if not MSA_UTILS_AVAILABLE:
            raise RuntimeError("MSA utilities not available")

        msa_result = prepare_msa_from_blast(
            query_sequence=query_sequence,
            query_name=query_name,
            blast_hits=blast_result.hits,
            max_sequences=self.BLAST_MAX_HOMOLOGS,
            min_sequences=self.MIN_SEQUENCES_FOR_MSA,
            identity_cutoff=self.REDUNDANCY_CUTOFF,
        )

        notes["msa_sequences"] = msa_result.n_sequences
        notes["msa_length"] = msa_result.alignment_length
        notes["msa_method"] = msa_result.method

        return msa_result.alignment, msa_result.n_sequences, notes

    def analyze_sync(
        self,
        pdb_content: str,
        chain: str = "A",
        **kwargs,
    ) -> ConservationAnalysis:
        """
        Synchronous wrapper for analyze().

        Use this when calling from non-async code.
        """
        return asyncio.run(self.analyze(pdb_content, chain, **kwargs))


# Convenience functions for integration

async def analyze_conservation(
    pdb_content: str,
    chain: str = "A",
    **kwargs,
) -> ConservationAnalysis:
    """
    Convenience function for conservation analysis.

    Args:
        pdb_content: PDB file content
        chain: Chain to analyze
        **kwargs: Additional arguments for ConservationAnalyzer.analyze()

    Returns:
        ConservationAnalysis result
    """
    analyzer = ConservationAnalyzer()
    return await analyzer.analyze(pdb_content, chain, **kwargs)


def get_conservation_summary(analysis: ConservationAnalysis) -> Dict[str, Any]:
    """
    Get summary of conservation analysis for design decisions.

    Args:
        analysis: ConservationAnalysis result

    Returns:
        Dict with key metrics for scaffolding/design
    """
    return {
        "msa_depth": analysis.msa_depth,
        "method": analysis.method,
        "average_grade": analysis.average_conservation,
        "n_highly_conserved": len(analysis.highly_conserved),
        "n_variable": len(analysis.variable),
        "highly_conserved_positions": analysis.highly_conserved,
        "suggested_fixed_residues": analysis.get_fixed_residue_suggestions(),
        "reliable_analysis": analysis.msa_depth >= 50,
    }


def apply_conservation_to_fixed_atoms(
    fixed_atoms: Dict[str, str],
    analysis: ConservationAnalysis,
    threshold: int = 3,
    fix_type: str = "CA",
) -> Dict[str, str]:
    """
    Add conserved residues to fixed atoms dict for RFD3.

    Args:
        fixed_atoms: Existing fixed atoms dict
        analysis: Conservation analysis result
        threshold: Max grade to fix (1-3 = highly conserved)
        fix_type: What to fix ("CA" for backbone, "all" for everything)

    Returns:
        Updated fixed_atoms dict
    """
    updated = fixed_atoms.copy()

    for grade in analysis.grades:
        if grade.grade <= threshold and grade.data_quality == "sufficient":
            residue_id = f"A{grade.position}"
            if residue_id not in updated:
                updated[residue_id] = fix_type

    return updated


# Module availability check
def check_conservation_available() -> Dict[str, bool]:
    """Check which conservation analysis components are available."""
    return {
        "blast_client": BLAST_CLIENT_AVAILABLE,
        "msa_utils": MSA_UTILS_AVAILABLE,
        "rate4site": RATE4SITE_AVAILABLE,
        "full_pipeline": (
            BLAST_CLIENT_AVAILABLE and
            MSA_UTILS_AVAILABLE and
            RATE4SITE_AVAILABLE
        ),
        "basic_pipeline": BLAST_CLIENT_AVAILABLE and MSA_UTILS_AVAILABLE,
    }
