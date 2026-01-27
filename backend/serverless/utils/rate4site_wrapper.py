"""
Rate4Site Binary Wrapper

Python wrapper for the Rate4Site C++ binary that calculates evolutionary
conservation scores from multiple sequence alignments.

Implements ConSurf methodology:
- Bayesian empirical method (default)
- Maximum likelihood method (optional)
- JTT substitution model (default)
- 1-9 grade assignment (ConSurf color scale)

Rate4Site Output Format (GradesPE):
    POS  SEQ  SCORE    QQ-INTERVAL     STD    MSA_DATA
    1    M   -0.982   [-1.281,-0.663]  0.316     50/50
"""

import subprocess
import tempfile
import shutil
import os
import re
import logging
from dataclasses import dataclass
from typing import List, Optional, Tuple, Dict, Any
from pathlib import Path

logger = logging.getLogger(__name__)


# Check if Rate4Site is available
RATE4SITE_AVAILABLE = False
RATE4SITE_PATH = None

# Check multiple possible locations
_RATE4SITE_PATHS = [
    os.environ.get("RATE4SITE_PATH"),
    shutil.which("rate4site"),
    "/usr/local/bin/rate4site",
    "/opt/rate4site/rate4site",
]

for path in _RATE4SITE_PATHS:
    if path and os.path.isfile(path) and os.access(path, os.X_OK):
        RATE4SITE_PATH = path
        RATE4SITE_AVAILABLE = True
        break


@dataclass
class Rate4SiteScore:
    """
    Rate4Site conservation score for a single position.

    Matches ConSurf GradesPE output format.
    """
    position: int           # 1-indexed position in sequence
    residue: str            # Single-letter amino acid
    score: float            # Conservation score (lower = more conserved)
    qq_interval: Tuple[float, float]  # Bayesian confidence interval
    std_dev: float          # Standard deviation
    msa_data: int           # Number of sequences with data at this position
    msa_total: int          # Total sequences in MSA
    data_quality: str       # "sufficient" or "insufficient"


@dataclass
class Rate4SiteResult:
    """Complete Rate4Site analysis result."""
    scores: List[Rate4SiteScore]
    method: str             # "bayesian" or "ml"
    model: str              # Substitution model used (JTT, WAG, etc.)
    tree_method: str        # Tree construction method
    msa_sequences: int      # Number of sequences in MSA
    msa_length: int         # Alignment length
    raw_output: str         # Full Rate4Site output text


def run_rate4site(
    msa_fasta: str,
    tree_newick: Optional[str] = None,
    method: str = "bayesian",
    model: str = "JTT",
    query_name: Optional[str] = None,
    bayes_interval: int = 3,
    timeout: int = 300,
) -> Rate4SiteResult:
    """
    Run Rate4Site conservation analysis.

    Args:
        msa_fasta: Multiple sequence alignment in FASTA format
        tree_newick: Optional phylogenetic tree in Newick format
        method: "bayesian" (default, ConSurf standard) or "ml"
        model: Substitution model: JTT, WAG, mtREV, cpREV, Dayhoff
        query_name: Name of query sequence in MSA (first sequence if None)
        bayes_interval: Confidence interval width (default 3, ConSurf standard)
        timeout: Maximum execution time in seconds

    Returns:
        Rate4SiteResult with per-position scores

    Raises:
        RuntimeError: If Rate4Site binary not available or execution fails
        TimeoutError: If execution exceeds timeout
    """
    if not RATE4SITE_AVAILABLE:
        raise RuntimeError(
            "Rate4Site binary not found. Install Rate4Site and set RATE4SITE_PATH "
            "environment variable or add to PATH."
        )

    # Create temp directory for files
    work_dir = tempfile.mkdtemp(prefix="rate4site_")

    try:
        # Write MSA file
        msa_path = os.path.join(work_dir, "alignment.fasta")
        with open(msa_path, 'w') as f:
            f.write(msa_fasta)

        # Extract query name from first sequence if not provided
        if not query_name:
            for line in msa_fasta.split("\n"):
                if line.startswith(">"):
                    query_name = line[1:].split()[0]
                    break

        # Build command
        output_path = os.path.join(work_dir, "grades.txt")

        cmd = [
            RATE4SITE_PATH,
            "-s", msa_path,              # MSA file
            "-o", output_path,           # Output file
            "-a", query_name,            # Query sequence name
            "-Mw",                       # Weight matrix normalization
        ]

        # Method selection
        if method.lower() == "bayesian":
            cmd.append("-ib")  # Use empirical Bayesian method
        else:
            cmd.append("-im")  # Use maximum likelihood method

        # Substitution model
        model_flags = {
            "JTT": "-Mj",
            "WAG": "-Mw",
            "mtREV": "-Mr",
            "cpREV": "-Mc",
            "Dayhoff": "-Md",
        }
        if model.upper() in model_flags:
            cmd.append(model_flags[model.upper()])

        # Phylogenetic tree
        if tree_newick:
            tree_path = os.path.join(work_dir, "tree.nwk")
            with open(tree_path, 'w') as f:
                f.write(tree_newick)
            cmd.extend(["-t", tree_path])

        logger.info(f"Running Rate4Site: {' '.join(cmd)}")

        # Execute
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
                cwd=work_dir,
            )
        except subprocess.TimeoutExpired:
            raise TimeoutError(f"Rate4Site exceeded {timeout}s timeout")

        if result.returncode != 0:
            logger.error(f"Rate4Site stderr: {result.stderr}")
            raise RuntimeError(f"Rate4Site failed: {result.stderr}")

        # Parse output
        if not os.path.exists(output_path):
            # Try alternative output names
            for alt_name in ["alignment.fasta.grades", "grades.res", "r4s.res"]:
                alt_path = os.path.join(work_dir, alt_name)
                if os.path.exists(alt_path):
                    output_path = alt_path
                    break
            else:
                raise RuntimeError(f"Rate4Site output not found in {work_dir}")

        with open(output_path, 'r') as f:
            output_text = f.read()

        scores = parse_rate4site_output(output_text)

        # Count MSA stats
        msa_sequences = msa_fasta.count(">")
        msa_length = len(scores) if scores else 0

        return Rate4SiteResult(
            scores=scores,
            method=method.lower(),
            model=model.upper(),
            tree_method="neighbor-joining" if not tree_newick else "user-provided",
            msa_sequences=msa_sequences,
            msa_length=msa_length,
            raw_output=output_text,
        )

    finally:
        # Cleanup temp directory
        shutil.rmtree(work_dir, ignore_errors=True)


def parse_rate4site_output(output_text: str) -> List[Rate4SiteScore]:
    """
    Parse Rate4Site output file.

    Handles multiple output formats produced by different Rate4Site versions.

    Args:
        output_text: Contents of Rate4Site output file

    Returns:
        List of Rate4SiteScore objects
    """
    scores = []

    # Pattern for standard Rate4Site output:
    # POS  SEQ  SCORE  QQ-INTERVAL  STD  DATA
    # Or: pos  amino  score  [lower,upper]  std  X/Y
    pattern = re.compile(
        r'^\s*(\d+)\s+'       # Position
        r'([A-Z\-])\s+'       # Amino acid
        r'([-\d.]+)\s+'       # Score
        r'\[?([-\d.]+)\s*,\s*([-\d.]+)\]?\s+'  # Interval
        r'([-\d.]+)\s+'       # Std dev
        r'(\d+)/(\d+)',       # MSA data
        re.MULTILINE
    )

    for match in pattern.finditer(output_text):
        pos = int(match.group(1))
        residue = match.group(2)
        score = float(match.group(3))
        interval_low = float(match.group(4))
        interval_high = float(match.group(5))
        std = float(match.group(6))
        msa_data = int(match.group(7))
        msa_total = int(match.group(8))

        # Determine data quality (ConSurf uses 6 as threshold)
        data_quality = "sufficient" if msa_data >= 6 else "insufficient"

        scores.append(Rate4SiteScore(
            position=pos,
            residue=residue,
            score=score,
            qq_interval=(interval_low, interval_high),
            std_dev=std,
            msa_data=msa_data,
            msa_total=msa_total,
            data_quality=data_quality,
        ))

    # Alternative parsing for simpler output format
    if not scores:
        # Try simple format: pos aa score
        simple_pattern = re.compile(
            r'^\s*(\d+)\s+([A-Z\-])\s+([-\d.]+)',
            re.MULTILINE
        )
        for match in simple_pattern.finditer(output_text):
            pos = int(match.group(1))
            residue = match.group(2)
            score = float(match.group(3))

            scores.append(Rate4SiteScore(
                position=pos,
                residue=residue,
                score=score,
                qq_interval=(score - 0.5, score + 0.5),  # Estimate
                std_dev=0.0,
                msa_data=0,
                msa_total=0,
                data_quality="unknown",
            ))

    return scores


def assign_consurf_grades(
    scores: List[Rate4SiteScore],
    method: str = "layers",
) -> List[int]:
    """
    Assign ConSurf color grades (1-9) based on Rate4Site scores.

    Implements ConSurf's assign_colors_according_to_r4s_layers() algorithm.

    Grade interpretation:
    - 1-3: Highly conserved (maroon colors in ConSurf)
    - 4-5: Average conservation (light colors)
    - 6-9: Variable positions (cyan/turquoise colors)

    Args:
        scores: List of Rate4SiteScore objects
        method: Grade assignment method:
            - "layers": Equal division of score range into 9 layers (ConSurf default)
            - "percentile": Assign based on percentile rank

    Returns:
        List of grades (1-9) corresponding to input scores
    """
    if not scores:
        return []

    raw_scores = [s.score for s in scores]
    min_score = min(raw_scores)
    max_score = max(raw_scores)
    score_range = max_score - min_score

    grades = []

    if method == "percentile":
        # Percentile-based assignment
        import numpy as np
        percentiles = [(np.sum(np.array(raw_scores) <= s) / len(raw_scores)) * 100
                       for s in raw_scores]
        for p in percentiles:
            # Lower percentile = more conserved = lower grade
            grade = int(p / 11.11) + 1  # Maps 0-100 to 1-9
            grade = max(1, min(9, grade))
            grades.append(grade)

    else:  # "layers" method (ConSurf default)
        if score_range == 0:
            # All same score = average conservation
            return [5] * len(scores)

        for s in scores:
            # Normalize score to 0-1 range (lower score = more conserved)
            normalized = (s.score - min_score) / score_range

            # ConSurf maps lower scores to lower grades (more conserved)
            # Grade 1 = most conserved, Grade 9 = most variable
            grade = int(normalized * 8) + 1
            grade = max(1, min(9, grade))
            grades.append(grade)

    return grades


def scores_to_consurf_format(
    scores: List[Rate4SiteScore],
    grades: List[int],
    pdb_chain: str = "A",
) -> str:
    """
    Format scores in ConSurf GradesPE format.

    Args:
        scores: Rate4Site scores
        grades: Assigned grades (1-9)
        pdb_chain: Chain identifier for residue labels

    Returns:
        Tab-delimited string in ConSurf format
    """
    lines = [
        "# ConSurf Conservation Grades",
        "# Generated by Banta Lab Conservation Analyzer",
        "#",
        "# POS\tSEQ\t3LATOM\tSCORE\tCOLOR\tCI_LOW\tCI_HIGH\tMSA_DATA\tVARIETY",
    ]

    # 3-letter amino acid codes
    AA_3LETTER = {
        "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
        "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
        "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
        "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
        "-": "GAP", "X": "UNK",
    }

    for score, grade in zip(scores, grades):
        res_3 = AA_3LETTER.get(score.residue, "UNK")
        atom_label = f"{res_3}{score.position}:{pdb_chain}"

        # Mark insufficient data with asterisk (ConSurf convention)
        grade_str = str(grade)
        if score.data_quality == "insufficient":
            grade_str += "*"

        line = (
            f"{score.position}\t"
            f"{score.residue}\t"
            f"{atom_label}\t"
            f"{score.score:.3f}\t"
            f"{grade_str}\t"
            f"{score.qq_interval[0]:.3f}\t"
            f"{score.qq_interval[1]:.3f}\t"
            f"{score.msa_data}/{score.msa_total}\t"
            f"-"  # Variety field (not computed)
        )
        lines.append(line)

    return "\n".join(lines)


def create_conservation_pdb(
    pdb_content: str,
    scores: List[Rate4SiteScore],
    grades: List[int],
    chain: str = "A",
) -> str:
    """
    Create PDB file with conservation grades in B-factor column.

    This allows visualization in PyMOL/Chimera using B-factor coloring.

    Args:
        pdb_content: Original PDB content
        scores: Rate4Site scores
        grades: Assigned grades (1-9)
        chain: Chain to annotate

    Returns:
        Modified PDB content with conservation grades as B-factors
    """
    # Create position -> grade mapping
    pos_to_grade = {}
    for score, grade in zip(scores, grades):
        pos_to_grade[score.position] = grade

    modified_lines = []

    for line in pdb_content.split("\n"):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                line_chain = line[21]
                if line_chain == chain:
                    resnum = int(line[22:26].strip())
                    grade = pos_to_grade.get(resnum, 5.0)

                    # Replace B-factor (columns 61-66)
                    # Normalize grade to 10-90 range for better visualization
                    b_factor = grade * 10.0

                    # Reconstruct line with new B-factor
                    new_line = line[:60] + f"{b_factor:6.2f}" + line[66:]
                    modified_lines.append(new_line)
                    continue
            except (ValueError, IndexError):
                pass

        modified_lines.append(line)

    return "\n".join(modified_lines)


# Fallback conservation calculation without Rate4Site
def calculate_simple_conservation(msa_fasta: str) -> List[Rate4SiteScore]:
    """
    Calculate simple conservation scores from MSA using Shannon entropy.

    This is a fallback when Rate4Site binary is not available.
    Results are less accurate but provide approximate conservation.

    Args:
        msa_fasta: Multiple sequence alignment in FASTA format

    Returns:
        List of Rate4SiteScore objects
    """
    import math

    # Parse MSA
    sequences = []
    current_seq = []

    for line in msa_fasta.split("\n"):
        if line.startswith(">"):
            if current_seq:
                sequences.append("".join(current_seq))
            current_seq = []
        else:
            current_seq.append(line.strip().upper())

    if current_seq:
        sequences.append("".join(current_seq))

    if not sequences:
        return []

    # Ensure equal length (should be for valid MSA)
    max_len = max(len(s) for s in sequences)
    sequences = [s.ljust(max_len, "-") for s in sequences]

    query_seq = sequences[0]
    n_seqs = len(sequences)

    scores = []
    for pos in range(len(query_seq)):
        residue = query_seq[pos]
        if residue == "-":
            continue

        # Count amino acids at this position
        aa_counts = {}
        non_gap = 0
        for seq in sequences:
            aa = seq[pos] if pos < len(seq) else "-"
            if aa != "-":
                aa_counts[aa] = aa_counts.get(aa, 0) + 1
                non_gap += 1

        # Calculate Shannon entropy
        entropy = 0.0
        if non_gap > 0:
            for count in aa_counts.values():
                p = count / non_gap
                if p > 0:
                    entropy -= p * math.log2(p)

        # Normalize: max entropy for 20 amino acids is log2(20) ≈ 4.32
        max_entropy = math.log2(20)
        normalized_entropy = entropy / max_entropy if max_entropy > 0 else 0

        # Convert to conservation score (lower = more conserved)
        # ConSurf uses negative scores for conserved
        conservation_score = normalized_entropy * 2 - 1  # Range -1 to 1

        scores.append(Rate4SiteScore(
            position=pos + 1,  # 1-indexed
            residue=residue,
            score=conservation_score,
            qq_interval=(conservation_score - 0.3, conservation_score + 0.3),
            std_dev=0.3,
            msa_data=non_gap,
            msa_total=n_seqs,
            data_quality="sufficient" if non_gap >= 6 else "insufficient",
        ))

    return scores
