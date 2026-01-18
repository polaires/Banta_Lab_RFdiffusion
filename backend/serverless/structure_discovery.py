# structure_discovery.py
"""
AI-Driven Structure Discovery Service

Unified interface for discovering protein structures with metal binding sites
from multiple databases. Provides natural language intent parsing and intelligent
aggregation of results from MetalPDB, RCSB PDB, UniProt, and AlphaFold.

This service is the primary interface for the AI-Driven Structure Discovery
Infrastructure in protein design workflows.

Example:
    >>> discovery = StructureDiscovery()
    >>> results = discovery.search_by_intent("terbium binding protein", limit=10)
    >>> template = discovery.get_best_template("TB", ligand="CIT")
    >>> coord_info = discovery.extract_coordination_info("1CA2", "ZN")
"""

import re
import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any, Tuple

# Import database adapters
from database_adapters import MetalPDBAdapter, UniProtAdapter, PubChemAdapter, AlphaFoldAdapter

# Import metal site fetcher for RCSB queries
from metal_site_fetcher import (
    query_metal_sites,
    query_metal_ligand_sites,
    generate_template_from_pdb,
    get_reference_template,
    REFERENCE_STRUCTURES,
)

# Import metal chemistry for coordination validation
from metal_chemistry import METAL_DATABASE, get_all_metals

# Import template system
from metal_ligand_templates import get_template_with_fallback

logger = logging.getLogger(__name__)

# =============================================================================
# CONSTANTS
# =============================================================================

# Metal name to code mapping for NL parsing
METAL_NAMES: Dict[str, str] = {
    # Transition metals
    "zinc": "ZN",
    "iron": "FE",
    "copper": "CU",
    "manganese": "MN",
    "cobalt": "CO",
    "nickel": "NI",
    "molybdenum": "MO",
    "tungsten": "W",
    # Alkaline earth metals
    "calcium": "CA",
    "magnesium": "MG",
    # Lanthanides (rare earth metals)
    "terbium": "TB",
    "europium": "EU",
    "gadolinium": "GD",
    "lanthanum": "LA",
    "cerium": "CE",
    "samarium": "SM",
    "ytterbium": "YB",
    "neodymium": "ND",
    "praseodymium": "PR",
    "dysprosium": "DY",
    "holmium": "HO",
    "erbium": "ER",
    "thulium": "TM",
    "lutetium": "LU",
    # Other
    "gold": "AU",
    "silver": "AG",
    "platinum": "PT",
    "cadmium": "CD",
    "mercury": "HG",
    "lead": "PB",
}

# Common ligand names to PDB codes
LIGAND_NAMES: Dict[str, str] = {
    # Cofactors
    "pqq": "PQQ",
    "pyrroloquinoline quinone": "PQQ",
    "citrate": "CIT",
    "citric acid": "CIT",
    "atp": "ATP",
    "adp": "ADP",
    "gtp": "GTP",
    "nad": "NAD",
    "nadh": "NAD",
    "fad": "FAD",
    "heme": "HEM",
    "haem": "HEM",
    # Amino acids as ligands
    "glutamate": "GLU",
    "aspartate": "ASP",
    "histidine": "HIS",
    "cysteine": "CYS",
}

# Protein type patterns for intent parsing
INTENT_PATTERNS: Dict[str, Dict[str, Any]] = {
    "dehydrogenase": {
        "keywords": ["dehydrogenase", "oxidoreductase", "reductase"],
        "metals": ["ZN", "FE", "MN", "CU"],
        "ligands": ["NAD", "FAD", "PQQ"],
    },
    "zinc_finger": {
        "keywords": ["zinc finger", "znf", "zf", "c2h2", "cys-his"],
        "metals": ["ZN"],
        "ligands": [],
    },
    "ef_hand": {
        "keywords": ["ef-hand", "ef hand", "calmodulin", "calcium binding"],
        "metals": ["CA"],
        "ligands": [],
    },
    "iron_sulfur": {
        "keywords": ["iron-sulfur", "fe-s", "ferredoxin", "iron sulfur"],
        "metals": ["FE"],
        "ligands": [],
    },
    "lanmodulin": {
        "keywords": ["lanmodulin", "lanthanide binding", "rare earth"],
        "metals": ["TB", "EU", "GD", "LA"],
        "ligands": ["CIT"],
    },
    "carbonic_anhydrase": {
        "keywords": ["carbonic anhydrase", "ca2"],
        "metals": ["ZN"],
        "ligands": [],
    },
    "superoxide_dismutase": {
        "keywords": ["superoxide dismutase", "sod", "superoxide"],
        "metals": ["ZN", "CU", "MN", "FE"],
        "ligands": [],
    },
    "metalloprotease": {
        "keywords": ["metalloprotease", "protease", "peptidase", "matrix metalloprotein"],
        "metals": ["ZN", "CO"],
        "ligands": [],
    },
}


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class DiscoveryResult:
    """
    Represents a structure discovery result.

    Attributes:
        pdb_id: PDB identifier
        metal: Metal element symbol
        resolution: Structure resolution in Angstroms (lower is better)
        coordination_number: Number of coordinating atoms
        geometry: Coordination geometry (e.g., "tetrahedral", "octahedral")
        source: Data source ("metalpdb", "rcsb", "uniprot", "alphafold")
        curated: Whether this is from a curated reference list
        ligands: List of non-protein ligand codes
        coordinating_residues: List of coordinating residue strings
        score: Ranking score (higher is better)
        uniprot_id: Optional UniProt accession
        protein_name: Optional protein name
        organism: Optional organism name
    """
    pdb_id: str
    metal: str
    resolution: float = 0.0
    coordination_number: int = 0
    geometry: str = ""
    source: str = ""
    curated: bool = False
    ligands: List[str] = field(default_factory=list)
    coordinating_residues: List[str] = field(default_factory=list)
    score: float = 0.0
    uniprot_id: str = ""
    protein_name: str = ""
    organism: str = ""

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        return {
            "pdb_id": self.pdb_id,
            "metal": self.metal,
            "resolution": self.resolution,
            "coordination_number": self.coordination_number,
            "geometry": self.geometry,
            "source": self.source,
            "curated": self.curated,
            "ligands": self.ligands,
            "coordinating_residues": self.coordinating_residues,
            "score": self.score,
            "uniprot_id": self.uniprot_id,
            "protein_name": self.protein_name,
            "organism": self.organism,
        }


# =============================================================================
# STRUCTURE DISCOVERY CLASS
# =============================================================================

class StructureDiscovery:
    """
    Unified structure discovery service for protein design.

    Aggregates results from multiple databases:
    - MetalPDB: Specialized metal binding site database
    - RCSB PDB: Primary experimental structure database
    - UniProt: Protein function and metal binding annotations
    - AlphaFold: Predicted structures

    Provides natural language intent parsing for intuitive queries
    and intelligent ranking of results.

    Example:
        >>> discovery = StructureDiscovery()
        >>> results = discovery.search_by_intent("terbium binding protein")
        >>> for r in results[:5]:
        ...     print(f"{r.pdb_id}: {r.metal} - {r.geometry}")
    """

    def __init__(
        self,
        timeout: int = 30,
        enable_fallback: bool = True,
    ):
        """
        Initialize the structure discovery service.

        Args:
            timeout: Request timeout in seconds for database queries
            enable_fallback: Whether to enable fallback behavior between databases
        """
        self.timeout = timeout
        self.enable_fallback = enable_fallback

        # Initialize adapters
        self.metalpdb = MetalPDBAdapter(timeout=timeout, enable_fallback=enable_fallback)
        self.uniprot = UniProtAdapter(timeout=timeout)
        self.pubchem = PubChemAdapter(timeout=timeout)
        self.alphafold = AlphaFoldAdapter(timeout=timeout)

    def search_by_intent(
        self,
        query: str,
        limit: int = 20,
    ) -> List[DiscoveryResult]:
        """
        Search for structures using natural language intent.

        Parses the query to extract:
        - Metal names (e.g., "zinc" -> ZN, "terbium" -> TB)
        - Ligand names (e.g., "PQQ", "citrate" -> CIT)
        - Protein types (e.g., "dehydrogenase", "zinc finger")

        Then queries appropriate databases based on the parsed intent.

        Args:
            query: Natural language query (e.g., "terbium binding protein",
                   "PQQ dehydrogenase", "zinc finger domain")
            limit: Maximum number of results to return

        Returns:
            List of DiscoveryResult objects, ranked by relevance and quality

        Example:
            >>> results = discovery.search_by_intent("terbium binding protein")
            >>> results = discovery.search_by_intent("PQQ dehydrogenase")
            >>> results = discovery.search_by_intent("zinc finger transcription factor")
        """
        # Parse intent
        metal_code, ligand_code, protein_type = self._parse_intent(query)

        logger.info(f"Parsed intent: metal={metal_code}, ligand={ligand_code}, type={protein_type}")

        results: List[DiscoveryResult] = []

        # Strategy 1: Query by metal + ligand if both present
        if metal_code and ligand_code:
            metal_ligand_results = self._search_metal_ligand(metal_code, ligand_code, limit)
            results.extend(metal_ligand_results)

        # Strategy 2: Query by metal only
        if metal_code:
            metal_results = self.search_by_metal(metal_code, limit=limit)
            results.extend(metal_results)

        # Strategy 3: Query by protein type (function search in UniProt)
        if protein_type:
            function_results = self._search_by_protein_type(protein_type, metal_code, limit)
            results.extend(function_results)

        # Strategy 4: If no specific metal/ligand parsed, do a keyword search
        if not metal_code and not ligand_code and not protein_type:
            keyword_results = self._search_by_keywords(query, limit)
            results.extend(keyword_results)

        # Deduplicate and rank
        results = self._deduplicate_results(results)
        results = self._rank_results(results)

        return results[:limit]

    def search_by_metal(
        self,
        metal: str,
        coordination_number: Optional[int] = None,
        geometry: Optional[str] = None,
        limit: int = 20,
    ) -> List[DiscoveryResult]:
        """
        Search for structures by metal type.

        Queries multiple databases in order of priority:
        1. Curated reference structures (high-quality templates)
        2. MetalPDB (specialized metal coordination database)
        3. RCSB PDB via metal_site_fetcher
        4. UniProt metal binding annotations

        Args:
            metal: Metal element symbol (e.g., "ZN", "FE", "TB")
            coordination_number: Optional filter by coordination number
            geometry: Optional filter by geometry (e.g., "tetrahedral")
            limit: Maximum number of results

        Returns:
            List of DiscoveryResult objects, ranked by resolution

        Example:
            >>> results = discovery.search_by_metal("ZN", coordination_number=4)
            >>> results = discovery.search_by_metal("TB", geometry="tricapped_trigonal_prismatic")
        """
        metal = metal.upper()
        results: List[DiscoveryResult] = []

        # Priority 1: Add curated reference structures
        curated_results = self._get_curated_references(metal)
        results.extend(curated_results)

        # Priority 2: Query MetalPDB
        try:
            metalpdb_results = self.metalpdb.search_by_metal(
                metal=metal,
                coordination_number=coordination_number,
                geometry=geometry,
                limit=limit,
            )
            for r in metalpdb_results:
                result = DiscoveryResult(
                    pdb_id=r.pdb_id if hasattr(r, 'pdb_id') else r.get('pdb_id', ''),
                    metal=metal,
                    resolution=r.resolution if hasattr(r, 'resolution') else r.get('resolution', 0.0),
                    coordination_number=r.coordination_number if hasattr(r, 'coordination_number') else r.get('coordination_number', 0),
                    geometry=r.geometry if hasattr(r, 'geometry') else r.get('geometry', ''),
                    source="metalpdb",
                    curated=False,
                    coordinating_residues=r.coordinating_residues if hasattr(r, 'coordinating_residues') else r.get('coordinating_residues', []),
                    ligands=r.ligands if hasattr(r, 'ligands') else r.get('ligands', []),
                )
                results.append(result)
        except Exception as e:
            logger.warning(f"MetalPDB query failed: {e}")

        # Priority 3: Query RCSB via metal_site_fetcher (if MetalPDB had limited results)
        if len(results) < limit and self.enable_fallback:
            try:
                rcsb_results = query_metal_sites(
                    metal=metal,
                    resolution_max=3.0,
                    limit=limit,
                )
                for r in rcsb_results:
                    result = DiscoveryResult(
                        pdb_id=r.get('pdb_id', ''),
                        metal=metal,
                        resolution=r.get('resolution', 0.0),
                        source="rcsb",
                        curated=False,
                    )
                    results.append(result)
            except Exception as e:
                logger.warning(f"RCSB query failed: {e}")

        # Priority 4: Query UniProt for metal binding proteins
        try:
            metal_name = self._get_metal_name(metal)
            uniprot_results = self.uniprot.search_by_metal_binding(
                metal=metal_name,
                limit=limit,
            )
            for r in uniprot_results:
                # Get PDB mappings for UniProt entries
                pdb_mappings = self.uniprot.get_pdb_mappings(r.get('uniprot_id', ''))
                for pdb in pdb_mappings[:2]:  # Limit to first 2 PDBs per protein
                    result = DiscoveryResult(
                        pdb_id=pdb.get('pdb_id', ''),
                        metal=metal,
                        source="uniprot",
                        curated=False,
                        uniprot_id=r.get('uniprot_id', ''),
                        protein_name=r.get('name', ''),
                        organism=r.get('organism', ''),
                    )
                    results.append(result)
        except Exception as e:
            logger.warning(f"UniProt query failed: {e}")

        # Deduplicate and rank
        results = self._deduplicate_results(results)
        results = self._rank_results(results)

        return results[:limit]

    def get_best_template(
        self,
        metal: str,
        ligand: Optional[str] = None,
        use_architector: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """
        Get the best template for a metal-ligand combination.

        Priority order:
        1. Curated templates (library templates with validated geometry)
        2. MetalPDB experimental structures
        3. RCSB PDB experimental structures
        4. Architector-generated coordinates (if enabled)

        Args:
            metal: Metal element symbol (e.g., "TB", "ZN")
            ligand: Optional ligand code (e.g., "CIT", "PQQ")
            use_architector: Whether to use Architector for calculated fallback

        Returns:
            Template dictionary with:
                - pdb_id: Source structure ID
                - pdb_content: PDB file content (if available)
                - metal: Metal code
                - ligand: Ligand code (if applicable)
                - coordination_number: Number of coordinating atoms
                - geometry: Coordination geometry
                - coordinating_residues: List of coordinating residue info
                - source: "curated", "metalpdb", "rcsb", or "architector"

        Example:
            >>> template = discovery.get_best_template("TB", ligand="CIT")
            >>> template = discovery.get_best_template("ZN")
        """
        metal = metal.upper()
        if ligand:
            ligand = ligand.upper()

        # Priority 1: Try curated templates
        if ligand:
            template_name = f"{ligand.lower()}_{metal.lower()}"
            template = get_template_with_fallback(
                template_name=template_name,
                metal=metal,
                ligand=ligand,
                fallback_enabled=use_architector,
            )
            if template:
                return template

        # Priority 2: Try MetalPDB/RCSB reference template
        reference_template = get_reference_template(metal, ligand)
        if reference_template:
            reference_template["source"] = "reference"
            return reference_template

        # Priority 3: Search for best structure and generate template
        if ligand:
            # Search for metal-ligand complexes
            ligand_results = query_metal_ligand_sites(
                metal=metal,
                ligand=ligand,
                limit=5,
            )
            if ligand_results:
                pdb_id = ligand_results[0].get('pdb_id', '')
                if pdb_id:
                    template = generate_template_from_pdb(pdb_id, metal)
                    if template:
                        template["ligand"] = ligand
                        return template

        # Priority 4: Get template from reference structures
        ref_pdb_ids = REFERENCE_STRUCTURES.get(metal, [])
        for pdb_id in ref_pdb_ids[:3]:
            template = generate_template_from_pdb(pdb_id, metal)
            if template:
                if ligand:
                    template["ligand"] = ligand
                return template

        # Priority 5: Architector fallback (calculated geometry)
        if use_architector:
            try:
                from metal_chemistry import get_coordination_number_range, get_preferred_donors

                ox_state = METAL_DATABASE.get(metal, {}).get("default_oxidation", 2)
                cn_range = get_coordination_number_range(metal, ox_state)
                donors = get_preferred_donors(metal, ox_state)

                # Generate calculated template
                template = {
                    "metal": metal,
                    "ligand": ligand,
                    "coordination_number": cn_range[1],  # Use max CN
                    "geometry": self._infer_geometry(cn_range[1]),
                    "source": "architector",
                    "calculated": True,
                    "warning": "Calculated template - validate geometry manually",
                }
                return template
            except Exception as e:
                logger.warning(f"Architector fallback failed: {e}")

        return None

    def extract_coordination_info(
        self,
        pdb_id: str,
        metal: str,
        site_index: int = 0,
    ) -> Optional[Dict[str, Any]]:
        """
        Extract detailed coordination information from a PDB structure.

        Queries MetalPDB first, then falls back to parsing the PDB file
        directly via metal_site_fetcher.

        Args:
            pdb_id: 4-character PDB identifier
            metal: Metal element symbol
            site_index: Which metal site to extract (for multi-metal structures)

        Returns:
            Dictionary with:
                - coordination_number: Number of coordinating atoms
                - geometry: Coordination geometry
                - coordinating_residues: List of residue info dicts
                - metal_coords: (x, y, z) coordinates
                - distances: Dict mapping residue to distance

        Example:
            >>> info = discovery.extract_coordination_info("1CA2", "ZN")
            >>> print(f"CN: {info['coordination_number']}, Geometry: {info['geometry']}")
        """
        pdb_id = pdb_id.upper()
        metal = metal.upper()

        # Try MetalPDB first (more detailed coordination info)
        try:
            site_details = self.metalpdb.get_site_details(pdb_id, metal, site_index)
            if site_details:
                return {
                    "pdb_id": pdb_id,
                    "metal": metal,
                    "coordination_number": site_details.coordination_number if hasattr(site_details, 'coordination_number') else site_details.get('coordination_number', 0),
                    "geometry": site_details.geometry if hasattr(site_details, 'geometry') else site_details.get('geometry', ''),
                    "coordinating_residues": site_details.coordinating_residues if hasattr(site_details, 'coordinating_residues') else site_details.get('coordinating_residues', []),
                    "metal_coords": site_details.metal_coords if hasattr(site_details, 'metal_coords') else site_details.get('metal_coords', (0.0, 0.0, 0.0)),
                    "distances": site_details.distances if hasattr(site_details, 'distances') else site_details.get('distances', {}),
                    "source": "metalpdb",
                }
        except Exception as e:
            logger.warning(f"MetalPDB site details failed: {e}")

        # Fallback to metal_site_fetcher (parses PDB directly)
        try:
            template = generate_template_from_pdb(
                pdb_id=pdb_id,
                metal=metal,
                metal_site_index=site_index,
            )
            if template:
                # Convert coordinating_residues format
                coord_residues = []
                distances = {}
                for res in template.get("coordinating_residues", []):
                    res_str = f"{res.get('residue', '')}{res.get('position', '')}"
                    coord_residues.append(res_str)
                    if res.get("distance"):
                        distances[res_str] = res.get("distance")

                return {
                    "pdb_id": pdb_id,
                    "metal": metal,
                    "coordination_number": template.get("coordination_number", 0),
                    "geometry": template.get("geometry", ""),
                    "coordinating_residues": coord_residues,
                    "metal_coords": template.get("metal_coords", (0.0, 0.0, 0.0)),
                    "distances": distances,
                    "source": "rcsb",
                }
        except Exception as e:
            logger.warning(f"metal_site_fetcher failed: {e}")

        return None

    # =========================================================================
    # INTERNAL PARSING METHODS
    # =========================================================================

    def _parse_intent(
        self,
        query: str,
    ) -> Tuple[Optional[str], Optional[str], Optional[str]]:
        """
        Parse natural language query to extract metal, ligand, and protein type.

        Args:
            query: Natural language query

        Returns:
            Tuple of (metal_code, ligand_code, protein_type)
            Any element may be None if not found in query
        """
        query_lower = query.lower()
        metal_code = None
        ligand_code = None
        protein_type = None

        # Extract metal name
        # Sort by length descending for longest-match
        sorted_metals = sorted(METAL_NAMES.items(), key=lambda x: len(x[0]), reverse=True)
        for name, code in sorted_metals:
            if name in query_lower:
                metal_code = code
                break

        # Also check for metal codes directly (e.g., "ZN", "TB")
        if not metal_code:
            for code in get_all_metals():
                # Match as whole word to avoid partial matches
                pattern = rf'\b{re.escape(code.lower())}\b'
                if re.search(pattern, query_lower):
                    metal_code = code
                    break

        # Extract ligand name
        # Sort by length descending for longest-match
        sorted_ligands = sorted(LIGAND_NAMES.items(), key=lambda x: len(x[0]), reverse=True)
        for name, code in sorted_ligands:
            if name in query_lower:
                ligand_code = code
                break

        # Check for 3-letter ligand codes directly
        if not ligand_code:
            # Look for uppercase 3-letter codes
            ligand_pattern = r'\b([A-Z]{3})\b'
            matches = re.findall(ligand_pattern, query)
            for match in matches:
                if match not in get_all_metals():  # Don't match metal codes
                    ligand_code = match
                    break

        # Extract protein type from patterns
        for ptype, info in INTENT_PATTERNS.items():
            for keyword in info["keywords"]:
                if keyword in query_lower:
                    protein_type = ptype
                    # If no metal found, suggest metals from this protein type
                    if not metal_code and info.get("metals"):
                        metal_code = info["metals"][0]
                    # If no ligand found, suggest ligands from this protein type
                    if not ligand_code and info.get("ligands"):
                        ligand_code = info["ligands"][0] if info["ligands"] else None
                    break
            if protein_type:
                break

        return metal_code, ligand_code, protein_type

    # =========================================================================
    # INTERNAL SEARCH METHODS
    # =========================================================================

    def _search_metal_ligand(
        self,
        metal: str,
        ligand: str,
        limit: int,
    ) -> List[DiscoveryResult]:
        """Search for metal-ligand complex structures."""
        results = []

        try:
            rcsb_results = query_metal_ligand_sites(
                metal=metal,
                ligand=ligand,
                limit=limit,
            )
            for r in rcsb_results:
                result = DiscoveryResult(
                    pdb_id=r.get('pdb_id', ''),
                    metal=metal,
                    resolution=r.get('resolution', 0.0),
                    source="rcsb",
                    curated=False,
                    ligands=[ligand],
                )
                results.append(result)
        except Exception as e:
            logger.warning(f"Metal-ligand search failed: {e}")

        return results

    def _search_by_protein_type(
        self,
        protein_type: str,
        metal: Optional[str],
        limit: int,
    ) -> List[DiscoveryResult]:
        """Search by protein type/function via UniProt."""
        results = []

        # Get keywords for this protein type
        pattern_info = INTENT_PATTERNS.get(protein_type, {})
        keywords = pattern_info.get("keywords", [protein_type])

        for keyword in keywords[:2]:  # Limit keyword searches
            try:
                uniprot_results = self.uniprot.search_by_function(
                    function_query=keyword,
                    limit=limit // 2,
                )
                for r in uniprot_results:
                    # Get PDB mappings
                    pdb_mappings = self.uniprot.get_pdb_mappings(r.get('uniprot_id', ''))
                    for pdb in pdb_mappings[:1]:  # First PDB only
                        result = DiscoveryResult(
                            pdb_id=pdb.get('pdb_id', ''),
                            metal=metal or "",
                            source="uniprot",
                            curated=False,
                            uniprot_id=r.get('uniprot_id', ''),
                            protein_name=r.get('name', ''),
                            organism=r.get('organism', ''),
                        )
                        results.append(result)
            except Exception as e:
                logger.warning(f"UniProt function search failed: {e}")

        return results

    def _search_by_keywords(
        self,
        query: str,
        limit: int,
    ) -> List[DiscoveryResult]:
        """General keyword search via UniProt."""
        results = []

        try:
            uniprot_results = self.uniprot.search_by_function(
                function_query=query,
                limit=limit,
            )
            for r in uniprot_results:
                pdb_mappings = self.uniprot.get_pdb_mappings(r.get('uniprot_id', ''))
                for pdb in pdb_mappings[:1]:
                    result = DiscoveryResult(
                        pdb_id=pdb.get('pdb_id', ''),
                        metal="",
                        source="uniprot",
                        curated=False,
                        uniprot_id=r.get('uniprot_id', ''),
                        protein_name=r.get('name', ''),
                        organism=r.get('organism', ''),
                    )
                    results.append(result)
        except Exception as e:
            logger.warning(f"Keyword search failed: {e}")

        return results

    def _get_curated_references(
        self,
        metal: str,
    ) -> List[DiscoveryResult]:
        """Get curated reference structures for a metal."""
        results = []

        ref_pdb_ids = REFERENCE_STRUCTURES.get(metal.upper(), [])
        for pdb_id in ref_pdb_ids:
            result = DiscoveryResult(
                pdb_id=pdb_id,
                metal=metal,
                source="curated",
                curated=True,
                score=100.0,  # High score for curated
            )
            results.append(result)

        return results

    # =========================================================================
    # INTERNAL RANKING METHODS
    # =========================================================================

    def _rank_results(
        self,
        results: List[DiscoveryResult],
    ) -> List[DiscoveryResult]:
        """
        Rank results by quality.

        Scoring factors:
        - Curated status: +100 points
        - Resolution: 10 - resolution (better resolution = higher score)
        - Source priority: curated > metalpdb > rcsb > uniprot
        """
        for result in results:
            score = 0.0

            # Curated bonus
            if result.curated:
                score += 100.0

            # Resolution score (lower is better, so invert)
            if result.resolution > 0:
                # Clamp resolution to reasonable range (0-10 for scoring)
                clamped_res = min(result.resolution, 10.0)
                score += max(0, 10 - clamped_res) * 5
            else:
                # Unknown resolution gets modest bonus (less than 5Å)
                score += 25

            # Source priority
            source_scores = {
                "curated": 50,
                "metalpdb": 30,
                "rcsb": 20,
                "uniprot": 10,
                "alphafold": 5,
            }
            score += source_scores.get(result.source, 0)

            # Coordination info bonus
            if result.coordination_number > 0:
                score += 5
            if result.geometry:
                score += 5

            result.score = score

        # Sort by score (descending)
        return sorted(results, key=lambda r: r.score, reverse=True)

    def _deduplicate_results(
        self,
        results: List[DiscoveryResult],
    ) -> List[DiscoveryResult]:
        """
        Remove duplicate PDB IDs, keeping the best entry for each.

        Preference: curated > higher resolution > first occurrence
        """
        seen: Dict[str, DiscoveryResult] = {}

        for result in results:
            pdb_id = result.pdb_id.upper()
            if not pdb_id:
                continue

            if pdb_id not in seen:
                seen[pdb_id] = result
            else:
                existing = seen[pdb_id]
                # Prefer curated
                if result.curated and not existing.curated:
                    seen[pdb_id] = result
                # Prefer better resolution (lower)
                elif result.resolution > 0 and (existing.resolution <= 0 or result.resolution < existing.resolution):
                    seen[pdb_id] = result
                # Prefer more info
                elif result.coordination_number > existing.coordination_number:
                    seen[pdb_id] = result

        return list(seen.values())

    # =========================================================================
    # UTILITY METHODS
    # =========================================================================

    def _get_metal_name(self, metal_code: str) -> str:
        """Convert metal code to name for UniProt queries."""
        code_to_name = {v: k for k, v in METAL_NAMES.items()}
        return code_to_name.get(metal_code.upper(), metal_code.lower())

    def _infer_geometry(self, coordination_number: int) -> str:
        """Infer geometry from coordination number."""
        geometry_map = {
            2: "linear",
            3: "trigonal_planar",
            4: "tetrahedral",
            5: "trigonal_bipyramidal",
            6: "octahedral",
            7: "pentagonal_bipyramidal",
            8: "square_antiprismatic",
            9: "tricapped_trigonal_prismatic",
        }
        return geometry_map.get(coordination_number, "unknown")
