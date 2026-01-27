"""
MSA Preparation Utilities

Implements ConSurf prepareMSA.pm functionality in Python:
1. Filter BLAST hits by coverage and identity
2. Cluster with MMseqs2 (CD-HIT alternative, 95% identity)
3. Align with MUSCLE

ConSurf Parameters:
- FRAGMENT_REDUNDANCY_RATE = 95% (clustering identity)
- FRAGMENT_MINIMUM_LENGTH = 60% (minimum coverage)
- MINIMUM_FRAGMENTS_FOR_MSA = 5
- MAXIMUM_MODIFIED_PERCENT = 15%
"""

import subprocess
import tempfile
import shutil
import os
import logging
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple, Any
from pathlib import Path

logger = logging.getLogger(__name__)


# Check for external tools
MMSEQS_AVAILABLE = shutil.which("mmseqs") is not None
MUSCLE_AVAILABLE = shutil.which("muscle") is not None
CLUSTALW_AVAILABLE = shutil.which("clustalw2") is not None or shutil.which("clustalw") is not None


@dataclass
class FilteredHit:
    """BLAST hit that passed filtering criteria."""
    accession: str
    sequence: str
    identity: float
    coverage: float
    evalue: float


@dataclass
class MSAResult:
    """Multiple sequence alignment result."""
    alignment: str          # FASTA format alignment
    n_sequences: int
    alignment_length: int
    query_name: str
    method: str             # "muscle", "clustalw", or "simple"


def filter_blast_hits(
    hits: List[Any],  # List of BlastHit from blast_client
    query_length: int,
    min_coverage: float = 0.60,     # FRAGMENT_MINIMUM_LENGTH
    min_identity: float = 0.35,     # Minimum identity to include
    max_identity: float = 0.95,     # Maximum identity (exclude near-identical)
    max_modified_percent: float = 0.15,  # MAXIMUM_MODIFIED_PERCENT
) -> List[FilteredHit]:
    """
    Filter BLAST hits according to ConSurf criteria.

    ConSurf filtering removes:
    - Fragments too short (< 60% coverage)
    - Sequences too similar (> 95% identity, redundant)
    - Sequences too different (< 35% identity, noise)
    - Sequences with too many non-standard characters

    Args:
        hits: List of BlastHit objects from BLAST search
        query_length: Length of query sequence
        min_coverage: Minimum query coverage (default 0.60)
        min_identity: Minimum sequence identity (default 0.35)
        max_identity: Maximum sequence identity (default 0.95)
        max_modified_percent: Maximum non-standard amino acids (default 0.15)

    Returns:
        List of FilteredHit objects that passed all criteria
    """
    filtered = []
    seen_sequences = set()  # Avoid exact duplicates

    for hit in hits:
        # Skip exact duplicates
        seq_hash = hash(hit.sequence)
        if seq_hash in seen_sequences:
            continue
        seen_sequences.add(seq_hash)

        # Check coverage
        if hit.coverage < min_coverage * 100:
            logger.debug(f"Filtered {hit.accession}: coverage {hit.coverage:.1f}% < {min_coverage*100}%")
            continue

        # Check identity range
        if hit.identity < min_identity * 100:
            logger.debug(f"Filtered {hit.accession}: identity {hit.identity:.1f}% < {min_identity*100}%")
            continue

        if hit.identity > max_identity * 100:
            logger.debug(f"Filtered {hit.accession}: identity {hit.identity:.1f}% > {max_identity*100}%")
            continue

        # Check for non-standard amino acids
        standard_aa = set("ACDEFGHIKLMNPQRSTVWY")
        non_standard = sum(1 for c in hit.sequence.upper() if c not in standard_aa)
        non_standard_ratio = non_standard / len(hit.sequence) if hit.sequence else 1.0

        if non_standard_ratio > max_modified_percent:
            logger.debug(f"Filtered {hit.accession}: non-standard ratio {non_standard_ratio:.2f} > {max_modified_percent}")
            continue

        filtered.append(FilteredHit(
            accession=hit.accession,
            sequence=hit.sequence,
            identity=hit.identity,
            coverage=hit.coverage,
            evalue=hit.evalue,
        ))

    # Sort by E-value (best first)
    filtered.sort(key=lambda h: h.evalue)

    logger.info(f"Filtered {len(hits)} BLAST hits to {len(filtered)} after quality filtering")
    return filtered


def cluster_sequences_mmseqs(
    sequences: Dict[str, str],
    identity: float = 0.95,  # FRAGMENT_REDUNDANCY_RATE
    coverage: float = 0.80,
    timeout: int = 300,
) -> Dict[str, str]:
    """
    Cluster sequences using MMseqs2 to remove redundancy.

    MMseqs2 is faster than CD-HIT and produces equivalent results.
    This implements ConSurf's redundancy removal step.

    Args:
        sequences: Dict mapping ID to sequence
        identity: Sequence identity threshold (default 0.95 = 95%)
        coverage: Minimum coverage for clustering (default 0.80)
        timeout: Maximum execution time

    Returns:
        Dict of representative sequences (reduced redundancy)
    """
    if not MMSEQS_AVAILABLE:
        logger.warning("MMseqs2 not available, using simple identity filter")
        return _simple_redundancy_filter(sequences, identity)

    if len(sequences) < 2:
        return sequences

    work_dir = tempfile.mkdtemp(prefix="mmseqs_")

    try:
        # Write sequences to FASTA
        input_fasta = os.path.join(work_dir, "input.fasta")
        with open(input_fasta, 'w') as f:
            for seq_id, seq in sequences.items():
                f.write(f">{seq_id}\n{seq}\n")

        # Create MMseqs2 database
        db_path = os.path.join(work_dir, "db")
        subprocess.run(
            ["mmseqs", "createdb", input_fasta, db_path],
            capture_output=True,
            timeout=timeout,
            check=True,
        )

        # Cluster
        cluster_path = os.path.join(work_dir, "cluster")
        tmp_path = os.path.join(work_dir, "tmp")
        os.makedirs(tmp_path, exist_ok=True)

        subprocess.run(
            [
                "mmseqs", "cluster",
                db_path,
                cluster_path,
                tmp_path,
                "--min-seq-id", str(identity),
                "-c", str(coverage),
                "--cov-mode", "0",  # bidirectional coverage
            ],
            capture_output=True,
            timeout=timeout,
            check=True,
        )

        # Extract representative sequences
        rep_path = os.path.join(work_dir, "rep")
        subprocess.run(
            ["mmseqs", "createsubdb", cluster_path, db_path, rep_path],
            capture_output=True,
            timeout=timeout,
            check=True,
        )

        # Convert to FASTA
        output_fasta = os.path.join(work_dir, "representatives.fasta")
        subprocess.run(
            ["mmseqs", "convert2fasta", rep_path, output_fasta],
            capture_output=True,
            timeout=timeout,
            check=True,
        )

        # Parse representatives
        representatives = {}
        current_id = None
        current_seq = []

        with open(output_fasta, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    if current_id:
                        representatives[current_id] = "".join(current_seq)
                    current_id = line[1:].strip().split()[0]
                    current_seq = []
                else:
                    current_seq.append(line.strip())

            if current_id:
                representatives[current_id] = "".join(current_seq)

        logger.info(f"MMseqs2 clustering: {len(sequences)} -> {len(representatives)} sequences at {identity*100}% identity")
        return representatives

    except subprocess.CalledProcessError as e:
        logger.error(f"MMseqs2 failed: {e.stderr.decode() if e.stderr else str(e)}")
        return _simple_redundancy_filter(sequences, identity)

    except subprocess.TimeoutExpired:
        logger.error(f"MMseqs2 timed out after {timeout}s")
        return _simple_redundancy_filter(sequences, identity)

    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def _simple_redundancy_filter(
    sequences: Dict[str, str],
    identity_threshold: float,
) -> Dict[str, str]:
    """
    Simple redundancy filter without MMseqs2.

    Uses basic pairwise identity calculation (O(n²) - only for small sets).
    """
    if len(sequences) <= 1:
        return sequences

    items = list(sequences.items())
    representatives = {items[0][0]: items[0][1]}

    for seq_id, seq in items[1:]:
        is_redundant = False

        for rep_id, rep_seq in representatives.items():
            # Quick length check
            len_ratio = min(len(seq), len(rep_seq)) / max(len(seq), len(rep_seq))
            if len_ratio < 0.5:
                continue

            # Simple identity calculation (no alignment)
            min_len = min(len(seq), len(rep_seq))
            matches = sum(1 for i in range(min_len) if seq[i] == rep_seq[i])
            identity = matches / min_len

            if identity >= identity_threshold:
                is_redundant = True
                break

        if not is_redundant:
            representatives[seq_id] = seq

    logger.info(f"Simple filtering: {len(sequences)} -> {len(representatives)} sequences")
    return representatives


def run_muscle_alignment(
    sequences: Dict[str, str],
    query_name: str,
    timeout: int = 300,
) -> MSAResult:
    """
    Create multiple sequence alignment using MUSCLE.

    MUSCLE is the ConSurf default aligner for protein sequences.

    Args:
        sequences: Dict mapping ID to sequence (query should be included)
        query_name: Name of query sequence
        timeout: Maximum execution time

    Returns:
        MSAResult with aligned sequences in FASTA format
    """
    if not MUSCLE_AVAILABLE:
        logger.warning("MUSCLE not available, using simple alignment")
        return _simple_alignment(sequences, query_name)

    if len(sequences) < 2:
        # Single sequence - no alignment needed
        fasta = "\n".join(f">{k}\n{v}" for k, v in sequences.items())
        return MSAResult(
            alignment=fasta,
            n_sequences=len(sequences),
            alignment_length=max(len(s) for s in sequences.values()) if sequences else 0,
            query_name=query_name,
            method="none",
        )

    work_dir = tempfile.mkdtemp(prefix="muscle_")

    try:
        # Write input FASTA (ensure query is first)
        input_fasta = os.path.join(work_dir, "input.fasta")
        with open(input_fasta, 'w') as f:
            # Write query first
            if query_name in sequences:
                f.write(f">{query_name}\n{sequences[query_name]}\n")

            # Write others
            for seq_id, seq in sequences.items():
                if seq_id != query_name:
                    f.write(f">{seq_id}\n{seq}\n")

        # Run MUSCLE
        output_fasta = os.path.join(work_dir, "aligned.fasta")

        # Try MUSCLE 5 syntax first
        try:
            result = subprocess.run(
                ["muscle", "-align", input_fasta, "-output", output_fasta],
                capture_output=True,
                timeout=timeout,
            )
            if result.returncode != 0:
                raise subprocess.CalledProcessError(result.returncode, "muscle")
        except (subprocess.CalledProcessError, FileNotFoundError):
            # Try MUSCLE 3 syntax
            result = subprocess.run(
                ["muscle", "-in", input_fasta, "-out", output_fasta],
                capture_output=True,
                timeout=timeout,
            )
            if result.returncode != 0:
                raise subprocess.CalledProcessError(result.returncode, "muscle", result.stderr)

        # Read alignment
        with open(output_fasta, 'r') as f:
            alignment = f.read()

        # Count sequences and alignment length
        n_seqs = alignment.count(">")
        aligned_seqs = []
        for line in alignment.split("\n"):
            if not line.startswith(">"):
                aligned_seqs.append(line.strip())
        align_len = max(len(s) for s in aligned_seqs) if aligned_seqs else 0

        logger.info(f"MUSCLE alignment: {n_seqs} sequences, {align_len} columns")

        return MSAResult(
            alignment=alignment,
            n_sequences=n_seqs,
            alignment_length=align_len,
            query_name=query_name,
            method="muscle",
        )

    except subprocess.TimeoutExpired:
        logger.error(f"MUSCLE timed out after {timeout}s")
        return _simple_alignment(sequences, query_name)

    except subprocess.CalledProcessError as e:
        logger.error(f"MUSCLE failed: {e.stderr.decode() if e.stderr else str(e)}")
        return _simple_alignment(sequences, query_name)

    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def _simple_alignment(
    sequences: Dict[str, str],
    query_name: str,
) -> MSAResult:
    """
    Simple pairwise alignment fallback.

    Just pads sequences to equal length - not a real MSA.
    Only use when MUSCLE is unavailable.
    """
    if not sequences:
        return MSAResult(
            alignment="",
            n_sequences=0,
            alignment_length=0,
            query_name=query_name,
            method="none",
        )

    max_len = max(len(s) for s in sequences.values())

    lines = []
    # Query first
    if query_name in sequences:
        padded = sequences[query_name].ljust(max_len, "-")
        lines.append(f">{query_name}")
        lines.append(padded)

    for seq_id, seq in sequences.items():
        if seq_id != query_name:
            padded = seq.ljust(max_len, "-")
            lines.append(f">{seq_id}")
            lines.append(padded)

    return MSAResult(
        alignment="\n".join(lines),
        n_sequences=len(sequences),
        alignment_length=max_len,
        query_name=query_name,
        method="simple_padding",
    )


def prepare_msa_from_blast(
    query_sequence: str,
    query_name: str,
    blast_hits: List[Any],  # BlastHit objects
    max_sequences: int = 500,
    min_sequences: int = 5,     # MINIMUM_FRAGMENTS_FOR_MSA
    identity_cutoff: float = 0.95,  # FRAGMENT_REDUNDANCY_RATE
) -> MSAResult:
    """
    Complete MSA preparation pipeline from BLAST results.

    Implements full ConSurf prepareMSA workflow:
    1. Filter BLAST hits by quality
    2. Remove redundancy with clustering
    3. Create alignment with MUSCLE

    Args:
        query_sequence: Query protein sequence
        query_name: Identifier for query
        blast_hits: BLAST search results
        max_sequences: Maximum sequences for MSA
        min_sequences: Minimum sequences required
        identity_cutoff: Clustering identity threshold

    Returns:
        MSAResult with final alignment

    Raises:
        ValueError: If insufficient homologs found
    """
    # Step 1: Filter BLAST hits
    filtered = filter_blast_hits(
        blast_hits,
        query_length=len(query_sequence),
        min_coverage=0.60,
        min_identity=0.35,
        max_identity=0.95,
    )

    if len(filtered) < min_sequences:
        raise ValueError(
            f"Insufficient homologs: found {len(filtered)}, need at least {min_sequences}. "
            "Conservation analysis may be unreliable."
        )

    # Step 2: Prepare sequences dict
    sequences = {query_name: query_sequence}
    for hit in filtered[:max_sequences]:
        sequences[hit.accession] = hit.sequence

    # Step 3: Cluster to remove redundancy
    clustered = cluster_sequences_mmseqs(sequences, identity=identity_cutoff)

    # Ensure query is in clustered set
    if query_name not in clustered:
        clustered[query_name] = query_sequence

    if len(clustered) < min_sequences:
        logger.warning(
            f"After clustering: {len(clustered)} sequences (min {min_sequences}). "
            "Using filtered set instead."
        )
        clustered = sequences

    # Step 4: Align with MUSCLE
    msa = run_muscle_alignment(clustered, query_name)

    logger.info(
        f"MSA preparation complete: {msa.n_sequences} sequences, "
        f"{msa.alignment_length} columns, method={msa.method}"
    )

    return msa


def parse_fasta(fasta_content: str) -> Dict[str, str]:
    """Parse FASTA format into dict."""
    sequences = {}
    current_id = None
    current_seq = []

    for line in fasta_content.split("\n"):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id:
                sequences[current_id] = "".join(current_seq)
            current_id = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line)

    if current_id:
        sequences[current_id] = "".join(current_seq)

    return sequences


def write_fasta(sequences: Dict[str, str], line_width: int = 60) -> str:
    """Write sequences to FASTA format."""
    lines = []
    for seq_id, seq in sequences.items():
        lines.append(f">{seq_id}")
        for i in range(0, len(seq), line_width):
            lines.append(seq[i:i+line_width])
    return "\n".join(lines)
