"""
NCBI BLAST API Client

Async client for NCBI BLAST web service API.
Implements ConSurf-style homolog search with rate limiting and caching.

NCBI Guidelines:
- Maximum 3 requests per second
- Include email in requests for heavy usage
- Use RID (Request ID) for polling results

ConSurf Parameters Applied:
- BLAST_MAX_HOMOLOGS = 500
- E-value threshold = 0.0001
- Database = nr (non-redundant protein)
"""

import asyncio
import aiohttp
import xml.etree.ElementTree as ET
import hashlib
import time
import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any, Tuple
from datetime import datetime, timedelta
from pathlib import Path
import tempfile
import json
import os

logger = logging.getLogger(__name__)


@dataclass
class BlastHit:
    """Single BLAST hit result."""
    accession: str
    description: str
    sequence: str
    evalue: float
    identity: float      # Percentage identity (0-100)
    coverage: float      # Query coverage (0-100)
    alignment_length: int
    query_start: int
    query_end: int
    hit_start: int
    hit_end: int
    score: float


@dataclass
class BlastResult:
    """Complete BLAST search result."""
    query_sequence: str
    query_length: int
    hits: List[BlastHit]
    database: str
    program: str
    search_time: float   # Seconds
    rid: str             # NCBI Request ID
    cached: bool = False


# Simple file-based cache for BLAST results
_CACHE_DIR = Path(tempfile.gettempdir()) / "blast_cache"
_CACHE_TTL = timedelta(minutes=15)


def _get_cache_key(sequence: str, database: str, max_hits: int) -> str:
    """Generate cache key from parameters."""
    content = f"{sequence}:{database}:{max_hits}"
    return hashlib.md5(content.encode()).hexdigest()


def _get_cached_result(cache_key: str) -> Optional[BlastResult]:
    """Retrieve cached BLAST result if valid."""
    cache_file = _CACHE_DIR / f"{cache_key}.json"

    if not cache_file.exists():
        return None

    try:
        with open(cache_file, 'r') as f:
            data = json.load(f)

        # Check TTL
        cached_time = datetime.fromisoformat(data.get("_cached_at", "2000-01-01"))
        if datetime.now() - cached_time > _CACHE_TTL:
            cache_file.unlink()  # Delete expired cache
            return None

        # Reconstruct BlastResult
        hits = [BlastHit(**h) for h in data.get("hits", [])]
        return BlastResult(
            query_sequence=data["query_sequence"],
            query_length=data["query_length"],
            hits=hits,
            database=data["database"],
            program=data["program"],
            search_time=data["search_time"],
            rid=data["rid"],
            cached=True,
        )
    except Exception as e:
        logger.warning(f"Cache read failed: {e}")
        return None


def _save_to_cache(cache_key: str, result: BlastResult) -> None:
    """Save BLAST result to cache."""
    try:
        _CACHE_DIR.mkdir(parents=True, exist_ok=True)

        data = {
            "query_sequence": result.query_sequence,
            "query_length": result.query_length,
            "hits": [
                {
                    "accession": h.accession,
                    "description": h.description,
                    "sequence": h.sequence,
                    "evalue": h.evalue,
                    "identity": h.identity,
                    "coverage": h.coverage,
                    "alignment_length": h.alignment_length,
                    "query_start": h.query_start,
                    "query_end": h.query_end,
                    "hit_start": h.hit_start,
                    "hit_end": h.hit_end,
                    "score": h.score,
                }
                for h in result.hits
            ],
            "database": result.database,
            "program": result.program,
            "search_time": result.search_time,
            "rid": result.rid,
            "_cached_at": datetime.now().isoformat(),
        }

        cache_file = _CACHE_DIR / f"{cache_key}.json"
        with open(cache_file, 'w') as f:
            json.dump(data, f)

    except Exception as e:
        logger.warning(f"Cache write failed: {e}")


class NCBIBlastClient:
    """
    Async NCBI BLAST API client.

    Implements NCBI's REST API for sequence similarity search.
    Follows NCBI usage guidelines with rate limiting.

    ConSurf-equivalent parameters:
    - database="nr" (non-redundant protein)
    - program="blastp"
    - max_hits=500 (BLAST_MAX_HOMOLOGUES_TO_DISPLAY)
    - evalue=0.0001 (ESCORE default)

    Usage:
        client = NCBIBlastClient()
        result = await client.search(sequence, max_hits=500)
        for hit in result.hits:
            print(f"{hit.accession}: {hit.identity}% identity")
    """

    BASE_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"

    # Rate limiting: max 3 requests per second
    _last_request_time: float = 0
    _rate_limit_lock: asyncio.Lock = None
    MIN_REQUEST_INTERVAL = 0.34  # ~3 requests per second

    # Polling configuration
    POLL_INTERVAL = 5     # Seconds between status checks
    MAX_POLL_TIME = 600   # Maximum wait time (10 minutes)

    def __init__(self, email: Optional[str] = None, tool: str = "banta_lab_conservation"):
        """
        Initialize BLAST client.

        Args:
            email: Contact email (recommended for heavy usage)
            tool: Tool identifier for NCBI
        """
        self.email = email
        self.tool = tool
        self._rate_limit_lock = asyncio.Lock()

    async def _rate_limit(self) -> None:
        """Enforce rate limiting between requests."""
        async with self._rate_limit_lock:
            now = time.time()
            elapsed = now - NCBIBlastClient._last_request_time
            if elapsed < self.MIN_REQUEST_INTERVAL:
                await asyncio.sleep(self.MIN_REQUEST_INTERVAL - elapsed)
            NCBIBlastClient._last_request_time = time.time()

    async def search(
        self,
        sequence: str,
        database: str = "nr",
        program: str = "blastp",
        max_hits: int = 500,
        evalue: float = 0.0001,
        use_cache: bool = True,
        timeout: int = 600,
    ) -> BlastResult:
        """
        Run BLAST search against NCBI database.

        Args:
            sequence: Query protein sequence (one-letter codes)
            database: NCBI database (nr, swissprot, pdb, etc.)
            program: BLAST program (blastp, psiblast, etc.)
            max_hits: Maximum number of hits to return
            evalue: E-value threshold
            use_cache: Whether to use cached results
            timeout: Maximum wait time in seconds

        Returns:
            BlastResult with hits sorted by E-value

        Raises:
            TimeoutError: If search exceeds timeout
            RuntimeError: If BLAST search fails
        """
        # Clean sequence
        sequence = sequence.upper().replace(" ", "").replace("\n", "")

        # Check cache
        if use_cache:
            cache_key = _get_cache_key(sequence, database, max_hits)
            cached = _get_cached_result(cache_key)
            if cached:
                logger.info(f"Using cached BLAST result (RID: {cached.rid})")
                return cached

        start_time = time.time()

        # Submit job
        rid = await self._submit_search(sequence, database, program, max_hits, evalue)
        logger.info(f"BLAST job submitted: RID={rid}")

        # Poll for results
        xml_result = await self._poll_results(rid, timeout)

        # Parse results
        hits = self._parse_xml_results(xml_result, sequence)

        search_time = time.time() - start_time

        result = BlastResult(
            query_sequence=sequence,
            query_length=len(sequence),
            hits=hits[:max_hits],  # Ensure we don't exceed max_hits
            database=database,
            program=program,
            search_time=search_time,
            rid=rid,
        )

        # Cache result
        if use_cache:
            _save_to_cache(cache_key, result)

        logger.info(f"BLAST completed: {len(result.hits)} hits in {search_time:.1f}s")
        return result

    async def _submit_search(
        self,
        sequence: str,
        database: str,
        program: str,
        max_hits: int,
        evalue: float,
    ) -> str:
        """
        Submit BLAST search job.

        Returns:
            RID (Request ID) for polling results
        """
        await self._rate_limit()

        params = {
            "CMD": "Put",
            "PROGRAM": program,
            "DATABASE": database,
            "QUERY": sequence,
            "EXPECT": str(evalue),
            "HITLIST_SIZE": str(max_hits),
            "FORMAT_TYPE": "XML",
            "TOOL": self.tool,
        }

        if self.email:
            params["EMAIL"] = self.email

        async with aiohttp.ClientSession() as session:
            async with session.post(self.BASE_URL, data=params) as response:
                if response.status != 200:
                    raise RuntimeError(f"BLAST submit failed: HTTP {response.status}")

                text = await response.text()

                # Extract RID from response
                rid = None
                for line in text.split("\n"):
                    if line.strip().startswith("RID ="):
                        rid = line.split("=")[1].strip()
                        break

                if not rid:
                    raise RuntimeError(f"Failed to get BLAST RID: {text[:500]}")

                return rid

    async def _poll_results(self, rid: str, timeout: int) -> str:
        """
        Poll for BLAST results.

        Args:
            rid: Request ID from job submission
            timeout: Maximum wait time in seconds

        Returns:
            XML results string
        """
        start_time = time.time()

        while time.time() - start_time < timeout:
            await self._rate_limit()

            # Check status
            params = {
                "CMD": "Get",
                "RID": rid,
                "FORMAT_TYPE": "XML",
                "FORMAT_OBJECT": "SearchInfo",
            }

            async with aiohttp.ClientSession() as session:
                async with session.get(self.BASE_URL, params=params) as response:
                    text = await response.text()

                    # Check status
                    if "Status=WAITING" in text:
                        logger.debug(f"BLAST job {rid} still running...")
                        await asyncio.sleep(self.POLL_INTERVAL)
                        continue

                    elif "Status=FAILED" in text:
                        raise RuntimeError(f"BLAST job {rid} failed")

                    elif "Status=UNKNOWN" in text:
                        raise RuntimeError(f"BLAST job {rid} not found (expired?)")

                    elif "Status=READY" in text:
                        # Get actual results
                        break

                    else:
                        # May already be results
                        if "<BlastOutput>" in text:
                            return text
                        await asyncio.sleep(self.POLL_INTERVAL)

        # Fetch final results
        await self._rate_limit()

        params = {
            "CMD": "Get",
            "RID": rid,
            "FORMAT_TYPE": "XML",
        }

        async with aiohttp.ClientSession() as session:
            async with session.get(self.BASE_URL, params=params) as response:
                text = await response.text()

                if "<BlastOutput>" not in text:
                    raise TimeoutError(f"BLAST job {rid} timed out after {timeout}s")

                return text

    def _parse_xml_results(self, xml_content: str, query_sequence: str) -> List[BlastHit]:
        """
        Parse BLAST XML output.

        Args:
            xml_content: BLAST XML output
            query_sequence: Original query sequence

        Returns:
            List of BlastHit objects sorted by E-value
        """
        hits = []

        try:
            root = ET.fromstring(xml_content)
            query_len = len(query_sequence)

            # Navigate to hits
            iterations = root.findall(".//Iteration")

            for iteration in iterations:
                for hit in iteration.findall(".//Hit"):
                    try:
                        hit_id = hit.find("Hit_id").text or ""
                        hit_def = hit.find("Hit_def").text or ""
                        hit_accession = hit.find("Hit_accession").text or hit_id.split("|")[1] if "|" in hit_id else hit_id

                        # Get best HSP (High-scoring Segment Pair)
                        hsps = hit.findall(".//Hsp")
                        if not hsps:
                            continue

                        # Use first (best) HSP
                        hsp = hsps[0]

                        evalue = float(hsp.find("Hsp_evalue").text)
                        score = float(hsp.find("Hsp_bit-score").text)
                        identity_count = int(hsp.find("Hsp_identity").text)
                        align_len = int(hsp.find("Hsp_align-len").text)
                        query_from = int(hsp.find("Hsp_query-from").text)
                        query_to = int(hsp.find("Hsp_query-to").text)
                        hit_from = int(hsp.find("Hsp_hit-from").text)
                        hit_to = int(hsp.find("Hsp_hit-to").text)

                        # Get aligned sequence (remove gaps for storage)
                        hit_seq = hsp.find("Hsp_hseq").text or ""
                        hit_seq_clean = hit_seq.replace("-", "")

                        # Calculate metrics
                        identity = (identity_count / align_len * 100) if align_len > 0 else 0
                        coverage = ((query_to - query_from + 1) / query_len * 100) if query_len > 0 else 0

                        hits.append(BlastHit(
                            accession=hit_accession,
                            description=hit_def,
                            sequence=hit_seq_clean,
                            evalue=evalue,
                            identity=identity,
                            coverage=coverage,
                            alignment_length=align_len,
                            query_start=query_from,
                            query_end=query_to,
                            hit_start=hit_from,
                            hit_end=hit_to,
                            score=score,
                        ))

                    except (AttributeError, ValueError, TypeError) as e:
                        logger.warning(f"Failed to parse hit: {e}")
                        continue

            # Sort by E-value (lower is better)
            hits.sort(key=lambda h: h.evalue)

        except ET.ParseError as e:
            logger.error(f"Failed to parse BLAST XML: {e}")
            raise RuntimeError(f"Failed to parse BLAST results: {e}")

        return hits


async def fetch_sequences(
    accessions: List[str],
    batch_size: int = 100,
) -> Dict[str, str]:
    """
    Fetch full sequences for accessions from NCBI.

    Args:
        accessions: List of NCBI accession numbers
        batch_size: Number of accessions per request

    Returns:
        Dict mapping accession to sequence
    """
    sequences = {}

    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i + batch_size]

        # Use NCBI E-utilities
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = {
            "db": "protein",
            "id": ",".join(batch),
            "rettype": "fasta",
            "retmode": "text",
        }

        async with aiohttp.ClientSession() as session:
            async with session.get(url, params=params) as response:
                if response.status != 200:
                    logger.warning(f"Failed to fetch sequences: HTTP {response.status}")
                    continue

                text = await response.text()

                # Parse FASTA
                current_acc = None
                current_seq = []

                for line in text.split("\n"):
                    if line.startswith(">"):
                        if current_acc:
                            sequences[current_acc] = "".join(current_seq)

                        # Extract accession from header
                        header = line[1:].split()[0]
                        # Handle various formats: sp|P12345|NAME, P12345, etc.
                        parts = header.split("|")
                        if len(parts) >= 2:
                            current_acc = parts[1]
                        else:
                            current_acc = parts[0]
                        current_seq = []
                    else:
                        current_seq.append(line.strip())

                if current_acc:
                    sequences[current_acc] = "".join(current_seq)

        # Rate limiting
        await asyncio.sleep(0.34)

    return sequences


def extract_sequence_from_pdb(pdb_content: str, chain: str = "A") -> str:
    """
    Extract amino acid sequence from PDB content.

    Args:
        pdb_content: PDB file content
        chain: Chain ID to extract

    Returns:
        One-letter amino acid sequence
    """
    # Standard amino acid 3-letter to 1-letter mapping
    AA_MAP = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        # Non-standard (map to closest)
        "MSE": "M",  # Selenomethionine
        "HSD": "H", "HSE": "H", "HSP": "H",  # Histidine variants
        "CYX": "C",  # Oxidized cysteine
    }

    # Try SEQRES records first (complete sequence)
    seqres_sequence = []
    for line in pdb_content.split("\n"):
        if line.startswith("SEQRES") and len(line) > 11:
            seqres_chain = line[11]
            if seqres_chain == chain:
                residues = line[19:].split()
                for res in residues:
                    aa = AA_MAP.get(res.upper(), "X")
                    seqres_sequence.append(aa)

    if seqres_sequence:
        return "".join(seqres_sequence)

    # Fallback: extract from ATOM records
    residues = {}  # (resnum, icode) -> resname
    for line in pdb_content.split("\n"):
        if not line.startswith("ATOM"):
            continue

        line_chain = line[21] if len(line) > 21 else ""
        if line_chain != chain:
            continue

        try:
            resname = line[17:20].strip()
            resnum = int(line[22:26].strip())
            icode = line[26] if len(line) > 26 else " "

            key = (resnum, icode)
            if key not in residues:
                residues[key] = resname
        except (ValueError, IndexError):
            continue

    # Sort by residue number and convert
    sorted_keys = sorted(residues.keys())
    sequence = ""
    for key in sorted_keys:
        resname = residues[key]
        aa = AA_MAP.get(resname.upper(), "X")
        if aa != "X":  # Skip unknown residues
            sequence += aa

    return sequence
