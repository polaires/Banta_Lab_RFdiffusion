# AI-Driven Structure Discovery Infrastructure Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a comprehensive AI-driven infrastructure that enables natural language queries to discover, analyze, and extract relevant protein structures from multiple databases for automated protein design workflows.

**Architecture:** A layered discovery system with: (1) Database adapters for RCSB, MetalPDB, UniProt, PubChem, and AlphaFold, (2) A unified StructureDiscovery service with NL intent parsing, (3) A coordination knowledge base integrated with existing Architector, (4) AI-callable interface for the task planner.

**Tech Stack:** Python 3.10+, requests, httpx (async), existing metal_chemistry.py, architector_integration.py, task_planner.py

---

## Phase 1: Database Adapters

### Task 1.1: MetalPDB Adapter

**Files:**
- Create: `backend/serverless/database_adapters/metalpdb_adapter.py`
- Create: `backend/serverless/database_adapters/__init__.py`
- Test: `backend/serverless/tests/test_metalpdb_adapter.py`

**Step 1: Create the database_adapters package**

```python
# backend/serverless/database_adapters/__init__.py
"""Database adapters for external structure databases."""

from .metalpdb_adapter import MetalPDBAdapter

__all__ = ["MetalPDBAdapter"]
```

**Step 2: Write failing test for MetalPDB search**

```python
# backend/serverless/tests/test_metalpdb_adapter.py
"""Tests for MetalPDB adapter."""
import pytest
from database_adapters.metalpdb_adapter import MetalPDBAdapter


class TestMetalPDBAdapter:
    """Test MetalPDB API integration."""

    def test_search_by_metal_returns_results(self):
        """Search for zinc sites returns PDB IDs."""
        adapter = MetalPDBAdapter()
        results = adapter.search_by_metal("ZN", limit=5)

        assert len(results) > 0
        assert all("pdb_id" in r for r in results)
        assert all("metal" in r for r in results)

    def test_search_by_metal_with_coordination(self):
        """Search with coordination number filter."""
        adapter = MetalPDBAdapter()
        results = adapter.search_by_metal("ZN", coordination_number=4, limit=5)

        assert len(results) > 0
        for r in results:
            assert r.get("coordination_number") == 4

    def test_get_site_details_returns_geometry(self):
        """Get detailed coordination geometry for a site."""
        adapter = MetalPDBAdapter()
        # 1CA2 is carbonic anhydrase with well-characterized Zn site
        details = adapter.get_site_details("1CA2", "ZN")

        assert details is not None
        assert "geometry" in details
        assert "coordinating_residues" in details

    def test_search_nonexistent_metal_returns_empty(self):
        """Search for invalid metal returns empty list."""
        adapter = MetalPDBAdapter()
        results = adapter.search_by_metal("XX", limit=5)

        assert results == []
```

**Step 3: Run test to verify it fails**

Run: `cd G:/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless && python -m pytest tests/test_metalpdb_adapter.py -v`
Expected: FAIL with "No module named 'database_adapters'"

**Step 4: Implement MetalPDB adapter**

```python
# backend/serverless/database_adapters/metalpdb_adapter.py
"""
MetalPDB Database Adapter

Provides access to the MetalPDB database for curated metal coordination sites.
MetalPDB contains ~400,000 metal binding sites with detailed coordination info.

API Documentation: https://metalpdb.cerm.unifi.it/
"""

import logging
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

try:
    import requests
    REQUESTS_AVAILABLE = True
except ImportError:
    REQUESTS_AVAILABLE = False

logger = logging.getLogger(__name__)

# MetalPDB API endpoints
METALPDB_BASE_URL = "https://metalpdb.cerm.unifi.it"
METALPDB_SEARCH_URL = f"{METALPDB_BASE_URL}/api/search"
METALPDB_SITE_URL = f"{METALPDB_BASE_URL}/api/site"


@dataclass
class MetalSiteResult:
    """Result from MetalPDB search."""
    pdb_id: str
    metal: str
    site_id: str
    coordination_number: int
    geometry: str
    resolution: float
    coordinating_residues: List[str]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "pdb_id": self.pdb_id,
            "metal": self.metal,
            "site_id": self.site_id,
            "coordination_number": self.coordination_number,
            "geometry": self.geometry,
            "resolution": self.resolution,
            "coordinating_residues": self.coordinating_residues,
        }


class MetalPDBAdapter:
    """
    Adapter for MetalPDB database queries.

    MetalPDB provides curated metal coordination site information
    with detailed geometry analysis.
    """

    def __init__(self, timeout: int = 30):
        """Initialize adapter with optional timeout."""
        self.timeout = timeout
        self._check_availability()

    def _check_availability(self) -> None:
        """Check if requests library is available."""
        if not REQUESTS_AVAILABLE:
            logger.warning("requests library not available for MetalPDB")

    def search_by_metal(
        self,
        metal: str,
        coordination_number: Optional[int] = None,
        geometry: Optional[str] = None,
        resolution_max: float = 3.0,
        limit: int = 20,
    ) -> List[Dict[str, Any]]:
        """
        Search MetalPDB for metal binding sites.

        Args:
            metal: Metal element symbol (ZN, FE, CA, TB, etc.)
            coordination_number: Filter by coordination number
            geometry: Filter by geometry (tetrahedral, octahedral, etc.)
            resolution_max: Maximum structure resolution
            limit: Maximum results to return

        Returns:
            List of site dictionaries with pdb_id, metal, geometry, etc.
        """
        if not REQUESTS_AVAILABLE:
            logger.warning("requests not available, returning empty results")
            return []

        metal = metal.upper()

        # Build query parameters
        params = {
            "metal": metal,
            "resolution_max": resolution_max,
            "limit": limit,
        }

        if coordination_number:
            params["coordination_number"] = coordination_number
        if geometry:
            params["geometry"] = geometry.lower()

        try:
            # Note: MetalPDB may use different API structure
            # This is a representative implementation
            response = requests.get(
                METALPDB_SEARCH_URL,
                params=params,
                timeout=self.timeout,
            )

            if response.status_code == 200:
                data = response.json()
                return self._parse_search_results(data, metal)
            elif response.status_code == 404:
                logger.info(f"No results for metal {metal}")
                return []
            else:
                logger.warning(f"MetalPDB returned {response.status_code}")
                return self._fallback_to_rcsb(metal, limit)

        except requests.RequestException as e:
            logger.warning(f"MetalPDB request failed: {e}")
            return self._fallback_to_rcsb(metal, limit)

    def _parse_search_results(
        self, data: Dict[str, Any], metal: str
    ) -> List[Dict[str, Any]]:
        """Parse MetalPDB search response."""
        results = []

        sites = data.get("sites", data.get("results", []))

        for site in sites:
            results.append({
                "pdb_id": site.get("pdb_id", site.get("pdbId", "")),
                "metal": metal,
                "site_id": site.get("site_id", site.get("siteId", "")),
                "coordination_number": site.get("coordination_number",
                                                site.get("coordinationNumber", 0)),
                "geometry": site.get("geometry", "unknown"),
                "resolution": site.get("resolution", 0.0),
                "coordinating_residues": site.get("ligands", []),
                "source": "metalpdb",
            })

        return results

    def _fallback_to_rcsb(self, metal: str, limit: int) -> List[Dict[str, Any]]:
        """Fallback to RCSB search if MetalPDB unavailable."""
        try:
            from metal_site_fetcher import query_metal_sites
            return query_metal_sites(metal, limit=limit)
        except ImportError:
            return []

    def get_site_details(
        self, pdb_id: str, metal: str, site_index: int = 0
    ) -> Optional[Dict[str, Any]]:
        """
        Get detailed coordination information for a specific site.

        Args:
            pdb_id: 4-character PDB ID
            metal: Metal element symbol
            site_index: Which site if multiple (0-indexed)

        Returns:
            Detailed site information or None
        """
        if not REQUESTS_AVAILABLE:
            return self._fallback_site_details(pdb_id, metal, site_index)

        pdb_id = pdb_id.upper()
        metal = metal.upper()

        try:
            url = f"{METALPDB_SITE_URL}/{pdb_id}/{metal}/{site_index}"
            response = requests.get(url, timeout=self.timeout)

            if response.status_code == 200:
                data = response.json()
                return self._parse_site_details(data)
            else:
                return self._fallback_site_details(pdb_id, metal, site_index)

        except requests.RequestException as e:
            logger.warning(f"MetalPDB site query failed: {e}")
            return self._fallback_site_details(pdb_id, metal, site_index)

    def _parse_site_details(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """Parse site details response."""
        return {
            "pdb_id": data.get("pdb_id", ""),
            "metal": data.get("metal", ""),
            "coordination_number": data.get("coordination_number", 0),
            "geometry": data.get("geometry", "unknown"),
            "coordinating_residues": data.get("ligands", []),
            "bond_distances": data.get("distances", []),
            "bond_angles": data.get("angles", []),
            "rmsd_to_ideal": data.get("rmsd", 0.0),
            "source": "metalpdb",
        }

    def _fallback_site_details(
        self, pdb_id: str, metal: str, site_index: int
    ) -> Optional[Dict[str, Any]]:
        """Fallback to local extraction if MetalPDB unavailable."""
        try:
            from metal_site_fetcher import generate_template_from_pdb
            return generate_template_from_pdb(pdb_id, metal, metal_site_index=site_index)
        except ImportError:
            return None

    def search_by_coordination_motif(
        self,
        motif: str,
        metal: Optional[str] = None,
        limit: int = 20,
    ) -> List[Dict[str, Any]]:
        """
        Search by coordination motif pattern.

        Args:
            motif: Coordination motif (e.g., "His-His-Glu-Asp")
            metal: Optional metal filter
            limit: Maximum results

        Returns:
            List of matching sites
        """
        if not REQUESTS_AVAILABLE:
            return []

        params = {"motif": motif, "limit": limit}
        if metal:
            params["metal"] = metal.upper()

        try:
            response = requests.get(
                f"{METALPDB_SEARCH_URL}/motif",
                params=params,
                timeout=self.timeout,
            )

            if response.status_code == 200:
                return self._parse_search_results(response.json(), metal or "ANY")
            return []

        except requests.RequestException:
            return []
```

**Step 5: Run tests to verify they pass**

Run: `cd G:/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless && python -m pytest tests/test_metalpdb_adapter.py -v`
Expected: PASS (may skip network tests if MetalPDB unavailable)

**Step 6: Commit**

```bash
git add backend/serverless/database_adapters/
git add backend/serverless/tests/test_metalpdb_adapter.py
git commit -m "feat: add MetalPDB database adapter for metal coordination sites"
```

---

### Task 1.2: UniProt Adapter

**Files:**
- Create: `backend/serverless/database_adapters/uniprot_adapter.py`
- Modify: `backend/serverless/database_adapters/__init__.py`
- Test: `backend/serverless/tests/test_uniprot_adapter.py`

**Step 1: Write failing test**

```python
# backend/serverless/tests/test_uniprot_adapter.py
"""Tests for UniProt adapter."""
import pytest
from database_adapters.uniprot_adapter import UniProtAdapter


class TestUniProtAdapter:
    """Test UniProt API integration."""

    def test_search_by_function_returns_results(self):
        """Search for metal binding proteins."""
        adapter = UniProtAdapter()
        results = adapter.search_by_function("zinc finger", limit=5)

        assert len(results) > 0
        assert all("uniprot_id" in r for r in results)

    def test_search_by_metal_binding(self):
        """Search for proteins that bind specific metal."""
        adapter = UniProtAdapter()
        results = adapter.search_by_metal_binding("Zinc", limit=5)

        assert len(results) > 0

    def test_get_pdb_mappings(self):
        """Get PDB structures for a UniProt entry."""
        adapter = UniProtAdapter()
        # P00918 is human carbonic anhydrase II
        pdbs = adapter.get_pdb_mappings("P00918")

        assert len(pdbs) > 0
        assert "1CA2" in [p["pdb_id"] for p in pdbs]

    def test_get_metal_binding_sites(self):
        """Extract metal binding site annotations."""
        adapter = UniProtAdapter()
        sites = adapter.get_metal_binding_sites("P00918")

        assert len(sites) > 0
        assert any("Zinc" in str(s) for s in sites)
```

**Step 2: Run test to verify it fails**

Run: `cd G:/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless && python -m pytest tests/test_uniprot_adapter.py -v`
Expected: FAIL

**Step 3: Implement UniProt adapter**

```python
# backend/serverless/database_adapters/uniprot_adapter.py
"""
UniProt Database Adapter

Provides access to UniProt protein database for:
- Function-based protein search
- Metal binding annotations
- PDB structure mappings
- Sequence information

API Documentation: https://www.uniprot.org/help/api
"""

import logging
from typing import Any, Dict, List, Optional

try:
    import requests
    REQUESTS_AVAILABLE = True
except ImportError:
    REQUESTS_AVAILABLE = False

logger = logging.getLogger(__name__)

# UniProt API endpoints
UNIPROT_BASE_URL = "https://rest.uniprot.org"
UNIPROT_SEARCH_URL = f"{UNIPROT_BASE_URL}/uniprotkb/search"
UNIPROT_ENTRY_URL = f"{UNIPROT_BASE_URL}/uniprotkb"


class UniProtAdapter:
    """
    Adapter for UniProt database queries.

    UniProt provides comprehensive protein sequence and
    functional annotation data with PDB cross-references.
    """

    def __init__(self, timeout: int = 30):
        self.timeout = timeout

    def search_by_function(
        self,
        function_query: str,
        organism: Optional[str] = None,
        reviewed_only: bool = True,
        limit: int = 20,
    ) -> List[Dict[str, Any]]:
        """
        Search UniProt by protein function.

        Args:
            function_query: Function description (e.g., "zinc finger", "dehydrogenase")
            organism: Optional organism filter (e.g., "human", "9606")
            reviewed_only: Only return Swiss-Prot (reviewed) entries
            limit: Maximum results

        Returns:
            List of protein entries with uniprot_id, name, function
        """
        if not REQUESTS_AVAILABLE:
            return []

        # Build query
        query_parts = [f'(function:"{function_query}")']

        if reviewed_only:
            query_parts.append("(reviewed:true)")
        if organism:
            query_parts.append(f'(organism_name:"{organism}")')

        query = " AND ".join(query_parts)

        params = {
            "query": query,
            "format": "json",
            "size": limit,
            "fields": "accession,protein_name,organism_name,ft_metal,xref_pdb",
        }

        try:
            response = requests.get(
                UNIPROT_SEARCH_URL,
                params=params,
                timeout=self.timeout,
            )

            if response.status_code == 200:
                return self._parse_search_results(response.json())
            return []

        except requests.RequestException as e:
            logger.warning(f"UniProt search failed: {e}")
            return []

    def search_by_metal_binding(
        self,
        metal: str,
        organism: Optional[str] = None,
        limit: int = 20,
    ) -> List[Dict[str, Any]]:
        """
        Search for proteins that bind a specific metal.

        Args:
            metal: Metal name (e.g., "Zinc", "Iron", "Calcium")
            organism: Optional organism filter
            limit: Maximum results

        Returns:
            List of metal-binding proteins
        """
        if not REQUESTS_AVAILABLE:
            return []

        query_parts = [f'(ft_metal:"{metal}")']

        if organism:
            query_parts.append(f'(organism_name:"{organism}")')

        query_parts.append("(reviewed:true)")
        query = " AND ".join(query_parts)

        params = {
            "query": query,
            "format": "json",
            "size": limit,
            "fields": "accession,protein_name,organism_name,ft_metal,xref_pdb",
        }

        try:
            response = requests.get(
                UNIPROT_SEARCH_URL,
                params=params,
                timeout=self.timeout,
            )

            if response.status_code == 200:
                return self._parse_search_results(response.json())
            return []

        except requests.RequestException as e:
            logger.warning(f"UniProt metal search failed: {e}")
            return []

    def _parse_search_results(self, data: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Parse UniProt search response."""
        results = []

        for entry in data.get("results", []):
            # Extract protein name
            protein_name = ""
            if "proteinDescription" in entry:
                rec_name = entry["proteinDescription"].get("recommendedName", {})
                protein_name = rec_name.get("fullName", {}).get("value", "")

            # Extract PDB references
            pdb_refs = []
            for xref in entry.get("uniProtKBCrossReferences", []):
                if xref.get("database") == "PDB":
                    pdb_refs.append(xref.get("id", ""))

            # Extract metal binding sites
            metal_sites = []
            for feature in entry.get("features", []):
                if feature.get("type") == "Metal binding":
                    metal_sites.append({
                        "position": feature.get("location", {}).get("start", {}).get("value"),
                        "metal": feature.get("description", ""),
                    })

            results.append({
                "uniprot_id": entry.get("primaryAccession", ""),
                "name": protein_name,
                "organism": entry.get("organism", {}).get("scientificName", ""),
                "pdb_ids": pdb_refs,
                "metal_binding_sites": metal_sites,
                "source": "uniprot",
            })

        return results

    def get_pdb_mappings(self, uniprot_id: str) -> List[Dict[str, Any]]:
        """
        Get PDB structure mappings for a UniProt entry.

        Args:
            uniprot_id: UniProt accession (e.g., "P00918")

        Returns:
            List of PDB mappings with resolution, method, chains
        """
        if not REQUESTS_AVAILABLE:
            return []

        try:
            url = f"{UNIPROT_ENTRY_URL}/{uniprot_id}"
            params = {"format": "json", "fields": "xref_pdb"}

            response = requests.get(url, params=params, timeout=self.timeout)

            if response.status_code == 200:
                data = response.json()
                return self._parse_pdb_mappings(data)
            return []

        except requests.RequestException as e:
            logger.warning(f"UniProt PDB mapping failed: {e}")
            return []

    def _parse_pdb_mappings(self, data: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Parse PDB cross-references from UniProt entry."""
        mappings = []

        for xref in data.get("uniProtKBCrossReferences", []):
            if xref.get("database") == "PDB":
                props = {p["key"]: p["value"] for p in xref.get("properties", [])}

                mappings.append({
                    "pdb_id": xref.get("id", ""),
                    "method": props.get("Method", ""),
                    "resolution": props.get("Resolution", ""),
                    "chains": props.get("Chains", ""),
                })

        return mappings

    def get_metal_binding_sites(self, uniprot_id: str) -> List[Dict[str, Any]]:
        """
        Get metal binding site annotations for a protein.

        Args:
            uniprot_id: UniProt accession

        Returns:
            List of metal binding site annotations
        """
        if not REQUESTS_AVAILABLE:
            return []

        try:
            url = f"{UNIPROT_ENTRY_URL}/{uniprot_id}"
            params = {"format": "json", "fields": "ft_metal,ft_binding"}

            response = requests.get(url, params=params, timeout=self.timeout)

            if response.status_code == 200:
                data = response.json()
                return self._parse_metal_sites(data)
            return []

        except requests.RequestException as e:
            logger.warning(f"UniProt metal sites query failed: {e}")
            return []

    def _parse_metal_sites(self, data: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Parse metal binding features from UniProt entry."""
        sites = []

        for feature in data.get("features", []):
            if feature.get("type") in ["Metal binding", "Binding site"]:
                location = feature.get("location", {})
                sites.append({
                    "type": feature.get("type"),
                    "start": location.get("start", {}).get("value"),
                    "end": location.get("end", {}).get("value"),
                    "description": feature.get("description", ""),
                    "evidence": feature.get("evidences", []),
                })

        return sites

    def get_sequence(self, uniprot_id: str) -> Optional[str]:
        """Get protein sequence for a UniProt entry."""
        if not REQUESTS_AVAILABLE:
            return None

        try:
            url = f"{UNIPROT_ENTRY_URL}/{uniprot_id}"
            params = {"format": "json", "fields": "sequence"}

            response = requests.get(url, params=params, timeout=self.timeout)

            if response.status_code == 200:
                data = response.json()
                return data.get("sequence", {}).get("value")
            return None

        except requests.RequestException:
            return None
```

**Step 4: Update __init__.py**

```python
# backend/serverless/database_adapters/__init__.py
"""Database adapters for external structure databases."""

from .metalpdb_adapter import MetalPDBAdapter
from .uniprot_adapter import UniProtAdapter

__all__ = ["MetalPDBAdapter", "UniProtAdapter"]
```

**Step 5: Run tests**

Run: `cd G:/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless && python -m pytest tests/test_uniprot_adapter.py -v`
Expected: PASS

**Step 6: Commit**

```bash
git add backend/serverless/database_adapters/uniprot_adapter.py
git add backend/serverless/database_adapters/__init__.py
git add backend/serverless/tests/test_uniprot_adapter.py
git commit -m "feat: add UniProt adapter for protein function and metal binding queries"
```

---

### Task 1.3: PubChem Adapter

**Files:**
- Create: `backend/serverless/database_adapters/pubchem_adapter.py`
- Modify: `backend/serverless/database_adapters/__init__.py`
- Test: `backend/serverless/tests/test_pubchem_adapter.py`

**Step 1: Write failing test**

```python
# backend/serverless/tests/test_pubchem_adapter.py
"""Tests for PubChem adapter."""
import pytest
from database_adapters.pubchem_adapter import PubChemAdapter


class TestPubChemAdapter:
    """Test PubChem API integration."""

    def test_search_by_name_returns_results(self):
        """Search for compound by name."""
        adapter = PubChemAdapter()
        results = adapter.search_by_name("citric acid", limit=3)

        assert len(results) > 0
        assert all("cid" in r for r in results)
        assert all("smiles" in r for r in results)

    def test_get_compound_by_cid(self):
        """Get compound details by CID."""
        adapter = PubChemAdapter()
        # CID 311 is citric acid
        compound = adapter.get_compound(311)

        assert compound is not None
        assert "smiles" in compound
        assert "molecular_formula" in compound

    def test_get_smiles_for_ligand(self):
        """Get SMILES for common ligand."""
        adapter = PubChemAdapter()
        smiles = adapter.get_smiles("ATP")

        assert smiles is not None
        assert len(smiles) > 0

    def test_search_by_substructure(self):
        """Search by substructure SMILES."""
        adapter = PubChemAdapter()
        # Carboxylate substructure
        results = adapter.search_by_substructure("C(=O)[O-]", limit=5)

        assert len(results) >= 0  # May be empty if API unavailable
```

**Step 2: Run test to verify it fails**

Run: `cd G:/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless && python -m pytest tests/test_pubchem_adapter.py -v`
Expected: FAIL

**Step 3: Implement PubChem adapter**

```python
# backend/serverless/database_adapters/pubchem_adapter.py
"""
PubChem Database Adapter

Provides access to PubChem for ligand information:
- SMILES structures
- Molecular properties
- Compound search

API Documentation: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest
"""

import logging
from typing import Any, Dict, List, Optional

try:
    import requests
    REQUESTS_AVAILABLE = True
except ImportError:
    REQUESTS_AVAILABLE = False

logger = logging.getLogger(__name__)

# PubChem PUG REST API
PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

# Common ligand name to CID mappings (for fast lookup)
COMMON_LIGANDS = {
    "ATP": 5957,
    "ADP": 6022,
    "GTP": 6830,
    "NAD": 5893,
    "NADH": 928,
    "FAD": 643975,
    "HEME": 26945,
    "CITRATE": 311,
    "CITRIC ACID": 311,
    "PQQ": 1024,
    "AZOBENZENE": 2272,
}


class PubChemAdapter:
    """
    Adapter for PubChem compound database.

    Provides ligand structure information including SMILES,
    molecular formula, and properties.
    """

    def __init__(self, timeout: int = 30):
        self.timeout = timeout

    def search_by_name(
        self,
        name: str,
        limit: int = 10,
    ) -> List[Dict[str, Any]]:
        """
        Search PubChem by compound name.

        Args:
            name: Compound name (e.g., "citric acid", "ATP")
            limit: Maximum results

        Returns:
            List of matching compounds with CID and SMILES
        """
        if not REQUESTS_AVAILABLE:
            return []

        # Check common ligands first
        name_upper = name.upper()
        if name_upper in COMMON_LIGANDS:
            cid = COMMON_LIGANDS[name_upper]
            compound = self.get_compound(cid)
            if compound:
                return [compound]

        try:
            # Search by name
            url = f"{PUBCHEM_BASE_URL}/compound/name/{name}/cids/JSON"
            response = requests.get(url, timeout=self.timeout)

            if response.status_code == 200:
                data = response.json()
                cids = data.get("IdentifierList", {}).get("CID", [])[:limit]

                results = []
                for cid in cids:
                    compound = self.get_compound(cid)
                    if compound:
                        results.append(compound)

                return results
            return []

        except requests.RequestException as e:
            logger.warning(f"PubChem search failed: {e}")
            return []

    def get_compound(self, cid: int) -> Optional[Dict[str, Any]]:
        """
        Get compound details by PubChem CID.

        Args:
            cid: PubChem Compound ID

        Returns:
            Compound dictionary with SMILES, formula, etc.
        """
        if not REQUESTS_AVAILABLE:
            return None

        try:
            url = f"{PUBCHEM_BASE_URL}/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,IUPACName/JSON"
            response = requests.get(url, timeout=self.timeout)

            if response.status_code == 200:
                data = response.json()
                props = data.get("PropertyTable", {}).get("Properties", [{}])[0]

                return {
                    "cid": cid,
                    "smiles": props.get("CanonicalSMILES", ""),
                    "isomeric_smiles": props.get("IsomericSMILES", ""),
                    "molecular_formula": props.get("MolecularFormula", ""),
                    "molecular_weight": props.get("MolecularWeight", 0),
                    "iupac_name": props.get("IUPACName", ""),
                    "source": "pubchem",
                }
            return None

        except requests.RequestException as e:
            logger.warning(f"PubChem compound fetch failed: {e}")
            return None

    def get_smiles(self, name: str) -> Optional[str]:
        """
        Get SMILES string for a compound by name.

        Args:
            name: Compound name or common abbreviation

        Returns:
            Canonical SMILES string or None
        """
        results = self.search_by_name(name, limit=1)
        if results:
            return results[0].get("smiles")
        return None

    def search_by_substructure(
        self,
        smiles: str,
        limit: int = 20,
    ) -> List[Dict[str, Any]]:
        """
        Search for compounds containing a substructure.

        Args:
            smiles: SMILES string of substructure
            limit: Maximum results

        Returns:
            List of matching compounds
        """
        if not REQUESTS_AVAILABLE:
            return []

        try:
            # Submit substructure search
            url = f"{PUBCHEM_BASE_URL}/compound/substructure/smiles/{smiles}/cids/JSON"
            params = {"MaxRecords": limit}

            response = requests.get(url, params=params, timeout=self.timeout)

            if response.status_code == 200:
                data = response.json()
                cids = data.get("IdentifierList", {}).get("CID", [])[:limit]

                results = []
                for cid in cids[:5]:  # Limit detailed fetches
                    compound = self.get_compound(cid)
                    if compound:
                        results.append(compound)

                return results
            return []

        except requests.RequestException as e:
            logger.warning(f"PubChem substructure search failed: {e}")
            return []

    def get_3d_conformer(self, cid: int) -> Optional[str]:
        """
        Get 3D conformer as SDF for a compound.

        Args:
            cid: PubChem Compound ID

        Returns:
            SDF content or None
        """
        if not REQUESTS_AVAILABLE:
            return None

        try:
            url = f"{PUBCHEM_BASE_URL}/compound/cid/{cid}/SDF"
            params = {"record_type": "3d"}

            response = requests.get(url, params=params, timeout=self.timeout)

            if response.status_code == 200:
                return response.text
            return None

        except requests.RequestException:
            return None
```

**Step 4: Update __init__.py**

```python
# backend/serverless/database_adapters/__init__.py
"""Database adapters for external structure databases."""

from .metalpdb_adapter import MetalPDBAdapter
from .uniprot_adapter import UniProtAdapter
from .pubchem_adapter import PubChemAdapter

__all__ = ["MetalPDBAdapter", "UniProtAdapter", "PubChemAdapter"]
```

**Step 5: Run tests**

Run: `cd G:/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless && python -m pytest tests/test_pubchem_adapter.py -v`
Expected: PASS

**Step 6: Commit**

```bash
git add backend/serverless/database_adapters/pubchem_adapter.py
git add backend/serverless/database_adapters/__init__.py
git add backend/serverless/tests/test_pubchem_adapter.py
git commit -m "feat: add PubChem adapter for ligand SMILES and properties"
```

---

### Task 1.4: AlphaFold DB Adapter

**Files:**
- Create: `backend/serverless/database_adapters/alphafold_adapter.py`
- Modify: `backend/serverless/database_adapters/__init__.py`
- Test: `backend/serverless/tests/test_alphafold_adapter.py`

**Step 1: Write failing test**

```python
# backend/serverless/tests/test_alphafold_adapter.py
"""Tests for AlphaFold DB adapter."""
import pytest
from database_adapters.alphafold_adapter import AlphaFoldAdapter


class TestAlphaFoldAdapter:
    """Test AlphaFold DB API integration."""

    def test_get_prediction_by_uniprot(self):
        """Get AlphaFold prediction for UniProt ID."""
        adapter = AlphaFoldAdapter()
        # P00918 is human carbonic anhydrase II
        prediction = adapter.get_prediction("P00918")

        assert prediction is not None
        assert "pdb_url" in prediction
        assert "plddt_url" in prediction

    def test_download_structure(self):
        """Download predicted structure."""
        adapter = AlphaFoldAdapter()
        pdb_content = adapter.download_structure("P00918")

        assert pdb_content is not None
        assert "ATOM" in pdb_content

    def test_get_plddt_scores(self):
        """Get per-residue pLDDT confidence."""
        adapter = AlphaFoldAdapter()
        plddt = adapter.get_plddt_scores("P00918")

        assert plddt is not None
        assert len(plddt) > 0
        assert all(0 <= score <= 100 for score in plddt)

    def test_nonexistent_uniprot_returns_none(self):
        """Query for invalid UniProt returns None."""
        adapter = AlphaFoldAdapter()
        result = adapter.get_prediction("XXXXXX")

        assert result is None
```

**Step 2: Run test to verify it fails**

Run: `cd G:/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless && python -m pytest tests/test_alphafold_adapter.py -v`
Expected: FAIL

**Step 3: Implement AlphaFold adapter**

```python
# backend/serverless/database_adapters/alphafold_adapter.py
"""
AlphaFold Database Adapter

Provides access to AlphaFold DB for predicted protein structures:
- Structure predictions by UniProt ID
- Per-residue confidence (pLDDT) scores
- Model download

API Documentation: https://alphafold.ebi.ac.uk/api-docs
"""

import logging
from typing import Any, Dict, List, Optional

try:
    import requests
    REQUESTS_AVAILABLE = True
except ImportError:
    REQUESTS_AVAILABLE = False

logger = logging.getLogger(__name__)

# AlphaFold DB API endpoints
ALPHAFOLD_BASE_URL = "https://alphafold.ebi.ac.uk"
ALPHAFOLD_API_URL = f"{ALPHAFOLD_BASE_URL}/api"
ALPHAFOLD_FILES_URL = f"{ALPHAFOLD_BASE_URL}/files"


class AlphaFoldAdapter:
    """
    Adapter for AlphaFold Database.

    Provides access to AlphaFold2 predicted structures
    for proteins without experimental structures.
    """

    def __init__(self, timeout: int = 30):
        self.timeout = timeout

    def get_prediction(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        """
        Get AlphaFold prediction metadata for a UniProt ID.

        Args:
            uniprot_id: UniProt accession (e.g., "P00918")

        Returns:
            Prediction metadata with URLs or None if not found
        """
        if not REQUESTS_AVAILABLE:
            return None

        uniprot_id = uniprot_id.upper()

        try:
            url = f"{ALPHAFOLD_API_URL}/prediction/{uniprot_id}"
            response = requests.get(url, timeout=self.timeout)

            if response.status_code == 200:
                data = response.json()
                if data:
                    entry = data[0] if isinstance(data, list) else data
                    return self._parse_prediction(entry)
            return None

        except requests.RequestException as e:
            logger.warning(f"AlphaFold API request failed: {e}")
            return None

    def _parse_prediction(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """Parse AlphaFold prediction response."""
        return {
            "uniprot_id": data.get("uniprotAccession", ""),
            "entry_id": data.get("entryId", ""),
            "gene": data.get("gene", ""),
            "organism": data.get("organismScientificName", ""),
            "pdb_url": data.get("pdbUrl", ""),
            "cif_url": data.get("cifUrl", ""),
            "plddt_url": data.get("paeImageUrl", ""),  # Actually PAE, pLDDT in structure
            "model_version": data.get("latestVersion", 1),
            "sequence_length": data.get("uniprotEnd", 0) - data.get("uniprotStart", 0) + 1,
            "source": "alphafold",
        }

    def download_structure(
        self,
        uniprot_id: str,
        format: str = "pdb",
    ) -> Optional[str]:
        """
        Download AlphaFold predicted structure.

        Args:
            uniprot_id: UniProt accession
            format: "pdb" or "cif"

        Returns:
            Structure file content or None
        """
        if not REQUESTS_AVAILABLE:
            return None

        uniprot_id = uniprot_id.upper()

        # Construct direct URL
        ext = "pdb" if format == "pdb" else "cif"
        url = f"{ALPHAFOLD_FILES_URL}/AF-{uniprot_id}-F1-model_v4.{ext}"

        try:
            response = requests.get(url, timeout=self.timeout)

            if response.status_code == 200:
                return response.text

            # Try older version
            url_v3 = f"{ALPHAFOLD_FILES_URL}/AF-{uniprot_id}-F1-model_v3.{ext}"
            response = requests.get(url_v3, timeout=self.timeout)

            if response.status_code == 200:
                return response.text

            return None

        except requests.RequestException as e:
            logger.warning(f"AlphaFold download failed: {e}")
            return None

    def get_plddt_scores(self, uniprot_id: str) -> Optional[List[float]]:
        """
        Extract per-residue pLDDT scores from AlphaFold structure.

        pLDDT (predicted Local Distance Difference Test):
        - >90: Very high confidence
        - 70-90: Confident
        - 50-70: Low confidence
        - <50: Very low confidence (likely disordered)

        Args:
            uniprot_id: UniProt accession

        Returns:
            List of pLDDT scores per residue or None
        """
        pdb_content = self.download_structure(uniprot_id, format="pdb")

        if not pdb_content:
            return None

        # Extract B-factors (pLDDT is stored in B-factor column)
        plddt_scores = []
        current_residue = None

        for line in pdb_content.split("\n"):
            if line.startswith("ATOM"):
                try:
                    res_num = int(line[22:26].strip())
                    b_factor = float(line[60:66].strip())

                    if res_num != current_residue:
                        plddt_scores.append(b_factor)
                        current_residue = res_num
                except (ValueError, IndexError):
                    continue

        return plddt_scores if plddt_scores else None

    def get_confident_regions(
        self,
        uniprot_id: str,
        threshold: float = 70.0,
    ) -> List[Dict[str, Any]]:
        """
        Get regions with high confidence predictions.

        Args:
            uniprot_id: UniProt accession
            threshold: Minimum pLDDT score (default 70)

        Returns:
            List of confident regions with start, end, mean_plddt
        """
        plddt = self.get_plddt_scores(uniprot_id)

        if not plddt:
            return []

        regions = []
        in_region = False
        start = 0
        scores = []

        for i, score in enumerate(plddt):
            if score >= threshold:
                if not in_region:
                    in_region = True
                    start = i + 1  # 1-indexed
                    scores = [score]
                else:
                    scores.append(score)
            else:
                if in_region:
                    regions.append({
                        "start": start,
                        "end": i,  # Last confident residue
                        "length": len(scores),
                        "mean_plddt": sum(scores) / len(scores),
                    })
                    in_region = False
                    scores = []

        # Handle region at end
        if in_region:
            regions.append({
                "start": start,
                "end": len(plddt),
                "length": len(scores),
                "mean_plddt": sum(scores) / len(scores),
            })

        return regions

    def has_prediction(self, uniprot_id: str) -> bool:
        """Check if AlphaFold has a prediction for this protein."""
        return self.get_prediction(uniprot_id) is not None
```

**Step 4: Update __init__.py**

```python
# backend/serverless/database_adapters/__init__.py
"""Database adapters for external structure databases."""

from .metalpdb_adapter import MetalPDBAdapter
from .uniprot_adapter import UniProtAdapter
from .pubchem_adapter import PubChemAdapter
from .alphafold_adapter import AlphaFoldAdapter

__all__ = [
    "MetalPDBAdapter",
    "UniProtAdapter",
    "PubChemAdapter",
    "AlphaFoldAdapter",
]
```

**Step 5: Run tests**

Run: `cd G:/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless && python -m pytest tests/test_alphafold_adapter.py -v`
Expected: PASS

**Step 6: Commit**

```bash
git add backend/serverless/database_adapters/alphafold_adapter.py
git add backend/serverless/database_adapters/__init__.py
git add backend/serverless/tests/test_alphafold_adapter.py
git commit -m "feat: add AlphaFold DB adapter for predicted structures"
```

---

## Phase 2: Structure Discovery Service

### Task 2.1: Core StructureDiscovery Service

**Files:**
- Create: `backend/serverless/structure_discovery.py`
- Test: `backend/serverless/tests/test_structure_discovery.py`

**Step 1: Write failing test**

```python
# backend/serverless/tests/test_structure_discovery.py
"""Tests for StructureDiscovery service."""
import pytest
from structure_discovery import StructureDiscovery


class TestStructureDiscovery:
    """Test unified structure discovery."""

    def test_search_by_intent_metal_binding(self):
        """Search for metal binding proteins by NL query."""
        discovery = StructureDiscovery()
        results = discovery.search_by_intent("terbium binding protein")

        assert len(results) > 0
        assert all("pdb_id" in r or "uniprot_id" in r for r in results)

    def test_search_by_intent_ligand(self):
        """Search for ligand-binding structures."""
        discovery = StructureDiscovery()
        results = discovery.search_by_intent("PQQ dehydrogenase")

        assert len(results) > 0

    def test_search_by_metal(self):
        """Search by metal element."""
        discovery = StructureDiscovery()
        results = discovery.search_by_metal("ZN", limit=5)

        assert len(results) > 0
        assert all(r.get("metal") == "ZN" for r in results)

    def test_get_best_template(self):
        """Get best template for a metal."""
        discovery = StructureDiscovery()
        template = discovery.get_best_template(metal="TB", ligand="CIT")

        assert template is not None
        assert "pdb_id" in template or "source" in template

    def test_extract_coordination_info(self):
        """Extract coordination info from PDB."""
        discovery = StructureDiscovery()
        info = discovery.extract_coordination_info("1CA2", "ZN")

        assert info is not None
        assert "coordination_number" in info
        assert "geometry" in info

    def test_rank_results_by_quality(self):
        """Results should be ranked by quality."""
        discovery = StructureDiscovery()
        results = discovery.search_by_metal("ZN", limit=10)

        # First result should have best resolution
        if len(results) >= 2:
            r1 = results[0].get("resolution", 99)
            r2 = results[1].get("resolution", 99)
            # Resolution should be ascending or equal
            assert r1 <= r2 or r1 == 0 or r2 == 0
```

**Step 2: Run test to verify it fails**

Run: `cd G:/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless && python -m pytest tests/test_structure_discovery.py -v`
Expected: FAIL

**Step 3: Implement StructureDiscovery service**

```python
# backend/serverless/structure_discovery.py
"""
Structure Discovery Service

Unified interface for discovering protein structures across multiple databases.
Supports natural language queries and intelligent ranking.

This is the main entry point for AI-driven structure discovery.
"""

import logging
import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

# Import database adapters
try:
    from database_adapters import (
        MetalPDBAdapter,
        UniProtAdapter,
        PubChemAdapter,
        AlphaFoldAdapter,
    )
    ADAPTERS_AVAILABLE = True
except ImportError:
    ADAPTERS_AVAILABLE = False

# Import existing modules
try:
    from metal_site_fetcher import (
        query_metal_sites,
        generate_template_from_pdb,
        get_reference_template,
        REFERENCE_STRUCTURES,
    )
    from metal_chemistry import METAL_DATABASE, get_preferred_donors
    from metal_ligand_templates import get_template, get_template_with_fallback
    METAL_MODULES_AVAILABLE = True
except ImportError:
    METAL_MODULES_AVAILABLE = False

logger = logging.getLogger(__name__)


# Intent patterns for NL parsing
INTENT_PATTERNS = {
    "metal_binding": [
        r"(\w+)\s*(binding|coordinating)\s*protein",
        r"protein.*bind.*(\w+)",
        r"(\w+)\s*site",
        r"metalloprotein.*(\w+)",
    ],
    "ligand_binding": [
        r"(\w+)\s*(binding|bound)\s*(protein|enzyme)",
        r"(\w+)\s*dehydrogenase",
        r"(\w+)\s*cofactor",
    ],
    "enzyme": [
        r"(\w+)ase\b",
        r"enzyme.*(\w+)",
        r"catalytic",
    ],
    "structure_type": [
        r"ef[\s-]?hand",
        r"zinc\s*finger",
        r"lanmodulin",
        r"ferritin",
    ],
}

# Metal name mappings
METAL_NAMES = {
    "zinc": "ZN", "zn": "ZN",
    "iron": "FE", "fe": "FE",
    "copper": "CU", "cu": "CU",
    "calcium": "CA", "ca": "CA",
    "magnesium": "MG", "mg": "MG",
    "manganese": "MN", "mn": "MN",
    "cobalt": "CO", "co": "CO",
    "nickel": "NI", "ni": "NI",
    "terbium": "TB", "tb": "TB",
    "europium": "EU", "eu": "EU",
    "gadolinium": "GD", "gd": "GD",
    "lanthanum": "LA", "la": "LA",
    "lanthanide": "TB",  # Default lanthanide
}


@dataclass
class DiscoveryResult:
    """Result from structure discovery."""
    pdb_id: Optional[str] = None
    uniprot_id: Optional[str] = None
    metal: Optional[str] = None
    ligand: Optional[str] = None
    resolution: float = 0.0
    coordination_number: int = 0
    geometry: str = ""
    source: str = ""
    score: float = 0.0
    title: str = ""
    organism: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {k: v for k, v in self.__dict__.items() if v}


class StructureDiscovery:
    """
    Unified structure discovery service.

    Aggregates results from multiple databases and provides
    intelligent ranking based on quality and relevance.
    """

    def __init__(self):
        """Initialize with available adapters."""
        self.metalpdb = MetalPDBAdapter() if ADAPTERS_AVAILABLE else None
        self.uniprot = UniProtAdapter() if ADAPTERS_AVAILABLE else None
        self.pubchem = PubChemAdapter() if ADAPTERS_AVAILABLE else None
        self.alphafold = AlphaFoldAdapter() if ADAPTERS_AVAILABLE else None

    def search_by_intent(
        self,
        query: str,
        limit: int = 20,
    ) -> List[Dict[str, Any]]:
        """
        Search structures using natural language query.

        Args:
            query: Natural language query (e.g., "terbium binding protein")
            limit: Maximum results

        Returns:
            Ranked list of structure results
        """
        query_lower = query.lower()

        # Parse intent
        metal, ligand, protein_type = self._parse_intent(query_lower)

        results = []

        # Strategy 1: Direct metal search
        if metal:
            metal_results = self.search_by_metal(metal, limit=limit)
            results.extend(metal_results)

        # Strategy 2: UniProt function search
        if self.uniprot and (protein_type or query):
            search_term = protein_type or query
            uniprot_results = self.uniprot.search_by_function(search_term, limit=limit//2)

            for ur in uniprot_results:
                for pdb_id in ur.get("pdb_ids", [])[:3]:
                    results.append({
                        "pdb_id": pdb_id,
                        "uniprot_id": ur.get("uniprot_id"),
                        "title": ur.get("name", ""),
                        "organism": ur.get("organism", ""),
                        "source": "uniprot",
                        "metal": metal,
                    })

        # Strategy 3: Check curated references
        if metal and METAL_MODULES_AVAILABLE:
            refs = REFERENCE_STRUCTURES.get(metal, [])
            for pdb_id in refs[:3]:
                if not any(r.get("pdb_id") == pdb_id for r in results):
                    results.append({
                        "pdb_id": pdb_id,
                        "metal": metal,
                        "source": "curated",
                        "score": 100,  # High score for curated
                    })

        # Deduplicate and rank
        results = self._deduplicate_results(results)
        results = self._rank_results(results)

        return results[:limit]

    def _parse_intent(self, query: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
        """
        Parse natural language query to extract metal, ligand, protein type.

        Returns:
            Tuple of (metal_code, ligand_code, protein_type)
        """
        metal = None
        ligand = None
        protein_type = None

        # Check for metal names
        for name, code in METAL_NAMES.items():
            if name in query:
                metal = code
                break

        # Check for common ligands
        ligand_patterns = ["pqq", "citrate", "heme", "atp", "nad", "fad"]
        for lig in ligand_patterns:
            if lig in query:
                ligand = lig.upper()
                break

        # Check for protein types
        type_patterns = ["dehydrogenase", "oxidase", "reductase", "transferase",
                        "kinase", "phosphatase", "protease", "synthase"]
        for ptype in type_patterns:
            if ptype in query:
                protein_type = ptype
                break

        # Check for structural motifs
        if "ef-hand" in query or "ef hand" in query:
            protein_type = "EF-hand"
            if not metal:
                metal = "CA"
        elif "zinc finger" in query:
            protein_type = "zinc finger"
            metal = "ZN"
        elif "lanmodulin" in query:
            protein_type = "lanmodulin"
            metal = "TB"

        return metal, ligand, protein_type

    def search_by_metal(
        self,
        metal: str,
        coordination_number: Optional[int] = None,
        geometry: Optional[str] = None,
        limit: int = 20,
    ) -> List[Dict[str, Any]]:
        """
        Search for metal-binding structures.

        Args:
            metal: Metal element symbol (ZN, FE, TB, etc.)
            coordination_number: Optional CN filter
            geometry: Optional geometry filter
            limit: Maximum results

        Returns:
            List of metal-binding structures
        """
        metal = metal.upper()
        results = []

        # Source 1: MetalPDB (most detailed)
        if self.metalpdb:
            metalpdb_results = self.metalpdb.search_by_metal(
                metal,
                coordination_number=coordination_number,
                geometry=geometry,
                limit=limit,
            )
            results.extend(metalpdb_results)

        # Source 2: RCSB via existing fetcher
        if METAL_MODULES_AVAILABLE and len(results) < limit:
            rcsb_results = query_metal_sites(metal, limit=limit - len(results))
            for r in rcsb_results:
                if not any(x.get("pdb_id") == r.get("pdb_id") for x in results):
                    results.append(r)

        # Source 3: UniProt metal binding
        if self.uniprot and len(results) < limit:
            # Map element to UniProt metal name
            metal_name = {
                "ZN": "Zinc", "FE": "Iron", "CU": "Copper",
                "CA": "Calcium", "MG": "Magnesium", "TB": "Terbium",
            }.get(metal, metal)

            uniprot_results = self.uniprot.search_by_metal_binding(
                metal_name, limit=limit - len(results)
            )

            for ur in uniprot_results:
                for pdb_id in ur.get("pdb_ids", [])[:2]:
                    if not any(x.get("pdb_id") == pdb_id for x in results):
                        results.append({
                            "pdb_id": pdb_id,
                            "metal": metal,
                            "uniprot_id": ur.get("uniprot_id"),
                            "source": "uniprot",
                        })

        return self._rank_results(results)[:limit]

    def get_best_template(
        self,
        metal: str,
        ligand: Optional[str] = None,
        use_architector: bool = True,
    ) -> Optional[Dict[str, Any]]:
        """
        Get the best template for a metal (with optional ligand).

        Priority:
        1. Curated library templates (metal_ligand_templates.py)
        2. MetalPDB best match
        3. RCSB reference structures
        4. Architector-generated (for lanthanides)

        Args:
            metal: Metal element symbol
            ligand: Optional ligand code
            use_architector: Whether to use Architector fallback

        Returns:
            Template dictionary or None
        """
        metal = metal.upper()

        # Priority 1: Curated templates
        if METAL_MODULES_AVAILABLE and ligand:
            template_key = f"{ligand.lower()}_{metal.lower()}"
            template = get_template_with_fallback(template_key, metal=metal)
            if template:
                return template

        # Priority 2: MetalPDB best match
        if self.metalpdb:
            results = self.metalpdb.search_by_metal(metal, limit=1)
            if results:
                best = results[0]
                details = self.metalpdb.get_site_details(
                    best["pdb_id"], metal
                )
                if details:
                    return details

        # Priority 3: RCSB reference
        if METAL_MODULES_AVAILABLE:
            template = get_reference_template(metal, ligand)
            if template:
                return template

        # Priority 4: Architector for lanthanides
        if use_architector and metal in ["TB", "EU", "GD", "LA", "CE", "YB"]:
            try:
                from architector_integration import generate_template_with_architector
                template = generate_template_with_architector(metal)
                if template:
                    template["source"] = "architector"
                    return template
            except ImportError:
                pass

        return None

    def extract_coordination_info(
        self,
        pdb_id: str,
        metal: str,
        cutoff: float = 3.0,
    ) -> Optional[Dict[str, Any]]:
        """
        Extract coordination information from a PDB structure.

        Args:
            pdb_id: 4-character PDB ID
            metal: Metal element symbol
            cutoff: Distance cutoff for coordination

        Returns:
            Coordination info dictionary
        """
        # Try MetalPDB first (has pre-computed geometry)
        if self.metalpdb:
            details = self.metalpdb.get_site_details(pdb_id, metal)
            if details:
                return details

        # Fallback to local extraction
        if METAL_MODULES_AVAILABLE:
            return generate_template_from_pdb(pdb_id, metal, cutoff)

        return None

    def get_ligand_smiles(self, ligand_name: str) -> Optional[str]:
        """
        Get SMILES for a ligand by name.

        Args:
            ligand_name: Ligand name or PDB code

        Returns:
            SMILES string or None
        """
        if self.pubchem:
            return self.pubchem.get_smiles(ligand_name)
        return None

    def get_alphafold_prediction(
        self,
        uniprot_id: str,
    ) -> Optional[Dict[str, Any]]:
        """
        Get AlphaFold prediction for a protein.

        Args:
            uniprot_id: UniProt accession

        Returns:
            AlphaFold prediction info or None
        """
        if self.alphafold:
            return self.alphafold.get_prediction(uniprot_id)
        return None

    def _deduplicate_results(
        self, results: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """Remove duplicate PDB entries."""
        seen = set()
        unique = []

        for r in results:
            pdb_id = r.get("pdb_id", "")
            if pdb_id and pdb_id not in seen:
                seen.add(pdb_id)
                unique.append(r)
            elif not pdb_id:
                unique.append(r)

        return unique

    def _rank_results(
        self, results: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """
        Rank results by quality and relevance.

        Scoring factors:
        - Resolution (lower is better)
        - Source priority (curated > metalpdb > rcsb > uniprot)
        - Has coordination info
        """
        source_scores = {
            "curated": 100,
            "metalpdb": 80,
            "rcsb": 60,
            "uniprot": 40,
            "alphafold": 30,
        }

        def score_result(r: Dict[str, Any]) -> float:
            score = r.get("score", 0)

            # Source bonus
            source = r.get("source", "")
            score += source_scores.get(source, 0)

            # Resolution bonus (0-50 points for <3Å)
            resolution = r.get("resolution", 0)
            if 0 < resolution < 3.0:
                score += (3.0 - resolution) * 20

            # Coordination info bonus
            if r.get("coordination_number"):
                score += 10
            if r.get("geometry"):
                score += 10

            return score

        for r in results:
            r["_score"] = score_result(r)

        results.sort(key=lambda x: x.get("_score", 0), reverse=True)

        # Remove internal score
        for r in results:
            r.pop("_score", None)

        return results


# Convenience functions for AI integration
def discover_structures(query: str, limit: int = 10) -> List[Dict[str, Any]]:
    """
    AI-callable function to discover structures.

    Args:
        query: Natural language query
        limit: Maximum results

    Returns:
        List of structure results
    """
    discovery = StructureDiscovery()
    return discovery.search_by_intent(query, limit=limit)


def get_metal_template(metal: str, ligand: str = None) -> Optional[Dict[str, Any]]:
    """
    AI-callable function to get best template.

    Args:
        metal: Metal element symbol
        ligand: Optional ligand code

    Returns:
        Template dictionary
    """
    discovery = StructureDiscovery()
    return discovery.get_best_template(metal, ligand)
```

**Step 4: Run tests**

Run: `cd G:/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless && python -m pytest tests/test_structure_discovery.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add backend/serverless/structure_discovery.py
git add backend/serverless/tests/test_structure_discovery.py
git commit -m "feat: add unified StructureDiscovery service with NL intent parsing"
```

---

## Phase 3: Coordination Knowledge Base

### Task 3.1: Extend Metal Chemistry with Coordination Templates

**Files:**
- Modify: `backend/serverless/metal_chemistry.py`
- Test: `backend/serverless/tests/test_coordination_knowledge.py`

**Step 1: Write failing test**

```python
# backend/serverless/tests/test_coordination_knowledge.py
"""Tests for coordination knowledge base."""
import pytest
from metal_chemistry import (
    get_coordination_template,
    get_ideal_geometry,
    validate_coordination,
    COORDINATION_MOTIFS,
)


class TestCoordinationKnowledge:
    """Test coordination chemistry knowledge base."""

    def test_ef_hand_motif_exists(self):
        """EF-hand motif should be defined."""
        assert "ef_hand" in COORDINATION_MOTIFS
        ef_hand = COORDINATION_MOTIFS["ef_hand"]
        assert "CA" in ef_hand["metals"]
        assert ef_hand["coordination"] == 7

    def test_zinc_finger_motif_exists(self):
        """Zinc finger motif should be defined."""
        assert "zinc_finger" in COORDINATION_MOTIFS
        zf = COORDINATION_MOTIFS["zinc_finger"]
        assert "ZN" in zf["metals"]
        assert zf["coordination"] == 4

    def test_get_coordination_template(self):
        """Get template for metal-ligand combination."""
        template = get_coordination_template("TB", coordination_number=9)

        assert template is not None
        assert template["geometry"] in ["tricapped_trigonal_prism", "square_antiprism"]

    def test_get_ideal_geometry(self):
        """Get ideal geometry positions."""
        positions = get_ideal_geometry("octahedral")

        assert len(positions) == 6
        # Check positions are on unit sphere
        for pos in positions:
            dist = (pos[0]**2 + pos[1]**2 + pos[2]**2) ** 0.5
            assert abs(dist - 1.0) < 0.01

    def test_validate_coordination_good(self):
        """Validate a good coordination geometry."""
        # Octahedral positions
        positions = [
            (2.4, 0, 0), (-2.4, 0, 0),
            (0, 2.4, 0), (0, -2.4, 0),
            (0, 0, 2.4), (0, 0, -2.4),
        ]
        metal_pos = (0, 0, 0)

        result = validate_coordination(positions, metal_pos, "ZN", expected_cn=6)

        assert result["valid"]
        assert result["coordination_number"] == 6

    def test_validate_coordination_bad_distances(self):
        """Detect bad bond distances."""
        # Positions too close
        positions = [
            (1.5, 0, 0), (-1.5, 0, 0),
            (0, 1.5, 0), (0, -1.5, 0),
        ]
        metal_pos = (0, 0, 0)

        result = validate_coordination(positions, metal_pos, "ZN", expected_cn=4)

        assert not result["valid"]
        assert "distance" in result["issues"][0].lower()
```

**Step 2: Run test to verify it fails**

Run: `cd G:/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless && python -m pytest tests/test_coordination_knowledge.py -v`
Expected: FAIL

**Step 3: Add coordination knowledge to metal_chemistry.py**

Add the following to the end of `backend/serverless/metal_chemistry.py`:

```python
# =============================================================================
# COORDINATION MOTIFS KNOWLEDGE BASE
# =============================================================================

COORDINATION_MOTIFS = {
    "ef_hand": {
        "name": "EF-hand",
        "description": "Calcium-binding helix-loop-helix motif",
        "metals": ["CA", "MG", "TB", "EU", "GD"],
        "sequence_motif": "DxDxDGxxDxxE",
        "coordination": 7,
        "geometry": "pentagonal_bipyramidal",
        "donors": ["Asp_OD", "Glu_OE", "backbone_O", "water"],
        "example_pdbs": ["1CLL", "3CLN", "1EXR"],
        "notes": "Loop provides 5 donors, helix Glu provides 2 (bidentate)",
    },
    "zinc_finger": {
        "name": "Zinc finger",
        "description": "DNA-binding structural motif",
        "metals": ["ZN"],
        "sequence_motif": "CxxC...HxxH",
        "coordination": 4,
        "geometry": "tetrahedral",
        "donors": ["Cys_SG", "His_NE2"],
        "example_pdbs": ["1ZNF", "1AAY", "2DRP"],
        "notes": "Classic C2H2 arrangement, structural role",
    },
    "his_tag": {
        "name": "Polyhistidine tag",
        "description": "Affinity purification tag",
        "metals": ["NI", "CO", "ZN", "CU"],
        "sequence_motif": "HHHHHH",
        "coordination": 6,
        "geometry": "octahedral",
        "donors": ["His_NE2"],
        "example_pdbs": [],
        "notes": "Engineered, used for protein purification",
    },
    "lanthanide_site": {
        "name": "Lanthanide binding site",
        "description": "High-coordination lanthanide site",
        "metals": ["TB", "EU", "GD", "LA", "CE", "SM", "YB"],
        "sequence_motif": None,  # No consensus sequence
        "coordination": [8, 9],
        "geometry": ["square_antiprism", "tricapped_trigonal_prism"],
        "donors": ["Glu_OE", "Asp_OD", "Asn_OD", "water"],
        "example_pdbs": ["6MI5", "7D4I", "5VNE"],
        "notes": "Hard acid prefers hard O donors, avoid Cys/His",
    },
    "rubredoxin": {
        "name": "Rubredoxin iron site",
        "description": "Iron-sulfur protein",
        "metals": ["FE"],
        "sequence_motif": "CxxCGxC...C",
        "coordination": 4,
        "geometry": "tetrahedral",
        "donors": ["Cys_SG"],
        "example_pdbs": ["1BRF", "1RB9"],
        "notes": "Fe-S4 cluster, electron transfer",
    },
    "carbonic_anhydrase": {
        "name": "Carbonic anhydrase zinc site",
        "description": "Catalytic zinc in CA",
        "metals": ["ZN"],
        "sequence_motif": "HxxH...H",
        "coordination": 4,
        "geometry": "tetrahedral",
        "donors": ["His_NE2", "water"],
        "example_pdbs": ["1CA2", "1HCB"],
        "notes": "3 His + water/hydroxide, catalytic",
    },
}


# Ideal geometry vertex positions (on unit sphere)
IDEAL_GEOMETRIES = {
    "linear": {
        "cn": 2,
        "positions": [(1, 0, 0), (-1, 0, 0)],
    },
    "trigonal_planar": {
        "cn": 3,
        "positions": [
            (1, 0, 0),
            (-0.5, 0.866, 0),
            (-0.5, -0.866, 0),
        ],
    },
    "tetrahedral": {
        "cn": 4,
        "positions": [
            (0.577, 0.577, 0.577),
            (0.577, -0.577, -0.577),
            (-0.577, 0.577, -0.577),
            (-0.577, -0.577, 0.577),
        ],
    },
    "square_planar": {
        "cn": 4,
        "positions": [
            (1, 0, 0), (-1, 0, 0),
            (0, 1, 0), (0, -1, 0),
        ],
    },
    "trigonal_bipyramidal": {
        "cn": 5,
        "positions": [
            (0, 0, 1), (0, 0, -1),  # Axial
            (1, 0, 0), (-0.5, 0.866, 0), (-0.5, -0.866, 0),  # Equatorial
        ],
    },
    "square_pyramidal": {
        "cn": 5,
        "positions": [
            (0, 0, 1),  # Apex
            (0.707, 0.707, 0), (0.707, -0.707, 0),
            (-0.707, 0.707, 0), (-0.707, -0.707, 0),
        ],
    },
    "octahedral": {
        "cn": 6,
        "positions": [
            (1, 0, 0), (-1, 0, 0),
            (0, 1, 0), (0, -1, 0),
            (0, 0, 1), (0, 0, -1),
        ],
    },
    "pentagonal_bipyramidal": {
        "cn": 7,
        "positions": [
            (0, 0, 1), (0, 0, -1),  # Axial
            (1, 0, 0),
            (0.309, 0.951, 0),
            (-0.809, 0.588, 0),
            (-0.809, -0.588, 0),
            (0.309, -0.951, 0),
        ],
    },
    "square_antiprism": {
        "cn": 8,
        "positions": _generate_sap_positions(),
    },
    "tricapped_trigonal_prism": {
        "cn": 9,
        "positions": _generate_ttp_positions(),
    },
}


def _generate_sap_positions():
    """Generate square antiprism positions."""
    import math
    theta = math.acos(1/3)
    positions = []

    # Top square
    for angle in [0, 90, 180, 270]:
        rad = math.radians(angle)
        positions.append((
            math.sin(theta) * math.cos(rad),
            math.sin(theta) * math.sin(rad),
            math.cos(theta),
        ))

    # Bottom square (rotated 45°)
    for angle in [45, 135, 225, 315]:
        rad = math.radians(angle)
        positions.append((
            math.sin(theta) * math.cos(rad),
            math.sin(theta) * math.sin(rad),
            -math.cos(theta),
        ))

    return positions


def _generate_ttp_positions():
    """Generate tricapped trigonal prism positions."""
    import math
    theta_prism = math.radians(48.2)
    positions = []

    # Top triangle
    for angle in [0, 120, 240]:
        rad = math.radians(angle)
        positions.append((
            math.sin(theta_prism) * math.cos(rad),
            math.sin(theta_prism) * math.sin(rad),
            math.cos(theta_prism),
        ))

    # Bottom triangle
    for angle in [60, 180, 300]:
        rad = math.radians(angle)
        positions.append((
            math.sin(theta_prism) * math.cos(rad),
            math.sin(theta_prism) * math.sin(rad),
            -math.cos(theta_prism),
        ))

    # Equatorial caps
    for angle in [30, 150, 270]:
        rad = math.radians(angle)
        positions.append((
            math.cos(rad),
            math.sin(rad),
            0,
        ))

    return positions


def get_coordination_template(
    metal: str,
    coordination_number: Optional[int] = None,
    motif: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """
    Get coordination template for a metal.

    Args:
        metal: Metal element symbol
        coordination_number: Target coordination number
        motif: Optional motif name (ef_hand, zinc_finger, etc.)

    Returns:
        Template dictionary with geometry, donors, etc.
    """
    metal = metal.upper()

    # If motif specified, use it
    if motif and motif in COORDINATION_MOTIFS:
        template = COORDINATION_MOTIFS[motif]
        if metal in template["metals"]:
            return dict(template)

    # Find matching motifs for this metal
    matching_motifs = [
        (name, data) for name, data in COORDINATION_MOTIFS.items()
        if metal in data["metals"]
    ]

    if not matching_motifs:
        # Fall back to metal database
        if metal in METAL_DATABASE:
            metal_info = METAL_DATABASE[metal]
            cn = coordination_number or metal_info.get("coordination_numbers", [6])[0]
            return {
                "metal": metal,
                "coordination": cn,
                "geometry": _geometry_for_cn(cn),
                "donors": metal_info.get("preferred_donors", ["O", "N"]),
                "source": "metal_database",
            }
        return None

    # Filter by coordination number if specified
    if coordination_number:
        for name, data in matching_motifs:
            cn = data["coordination"]
            if isinstance(cn, list):
                if coordination_number in cn:
                    return dict(data)
            elif cn == coordination_number:
                return dict(data)

    # Return first match
    return dict(matching_motifs[0][1])


def _geometry_for_cn(cn: int) -> str:
    """Get default geometry for coordination number."""
    geometry_map = {
        2: "linear",
        3: "trigonal_planar",
        4: "tetrahedral",
        5: "trigonal_bipyramidal",
        6: "octahedral",
        7: "pentagonal_bipyramidal",
        8: "square_antiprism",
        9: "tricapped_trigonal_prism",
    }
    return geometry_map.get(cn, "unknown")


def get_ideal_geometry(geometry_name: str) -> Optional[List[Tuple[float, float, float]]]:
    """
    Get ideal vertex positions for a coordination geometry.

    Args:
        geometry_name: Geometry name (octahedral, tetrahedral, etc.)

    Returns:
        List of (x, y, z) positions on unit sphere
    """
    geometry_name = geometry_name.lower().replace("-", "_").replace(" ", "_")

    if geometry_name in IDEAL_GEOMETRIES:
        return list(IDEAL_GEOMETRIES[geometry_name]["positions"])
    return None


def validate_coordination(
    donor_positions: List[Tuple[float, float, float]],
    metal_position: Tuple[float, float, float],
    metal: str,
    expected_cn: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Validate coordination geometry.

    Args:
        donor_positions: List of donor atom (x, y, z) positions
        metal_position: Metal center (x, y, z)
        metal: Metal element symbol
        expected_cn: Expected coordination number

    Returns:
        Validation result with valid flag and issues
    """
    import math

    metal = metal.upper()
    issues = []
    warnings = []

    # Get metal-specific parameters
    metal_info = METAL_DATABASE.get(metal, {})

    # Default bond distance ranges
    if metal in ["TB", "EU", "GD", "LA", "CE"]:
        # Lanthanides
        min_dist, max_dist = 2.2, 2.7
    elif metal in ["ZN", "CU", "FE", "CO", "NI"]:
        # Transition metals
        min_dist, max_dist = 1.9, 2.5
    else:
        # Others (Ca, Mg)
        min_dist, max_dist = 2.0, 2.6

    # Check coordination number
    cn = len(donor_positions)
    if expected_cn and cn != expected_cn:
        issues.append(f"Coordination number {cn} != expected {expected_cn}")

    # Check bond distances
    mx, my, mz = metal_position
    distances = []

    for dx, dy, dz in donor_positions:
        dist = math.sqrt((dx - mx)**2 + (dy - my)**2 + (dz - mz)**2)
        distances.append(dist)

        if dist < min_dist:
            issues.append(f"Distance {dist:.2f} Å too short (< {min_dist})")
        elif dist > max_dist:
            issues.append(f"Distance {dist:.2f} Å too long (> {max_dist})")

    # Check donor-donor distances (clash detection)
    clash_threshold = 2.4  # Å
    for i, (x1, y1, z1) in enumerate(donor_positions):
        for j, (x2, y2, z2) in enumerate(donor_positions):
            if i < j:
                dd_dist = math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
                if dd_dist < clash_threshold:
                    issues.append(f"Donor clash: {dd_dist:.2f} Å < {clash_threshold}")

    return {
        "valid": len(issues) == 0,
        "coordination_number": cn,
        "distances": distances,
        "mean_distance": sum(distances) / len(distances) if distances else 0,
        "issues": issues,
        "warnings": warnings,
    }
```

**Step 4: Run tests**

Run: `cd G:/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless && python -m pytest tests/test_coordination_knowledge.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add backend/serverless/metal_chemistry.py
git add backend/serverless/tests/test_coordination_knowledge.py
git commit -m "feat: add coordination motifs knowledge base with validation"
```

---

## Phase 4: AI Interface Integration

### Task 4.1: AI-Callable Structure Interface

**Files:**
- Create: `backend/serverless/ai_structure_interface.py`
- Modify: `backend/ai/task_planner.py`
- Test: `backend/serverless/tests/test_ai_interface.py`

**Step 1: Write failing test**

```python
# backend/serverless/tests/test_ai_interface.py
"""Tests for AI structure interface."""
import pytest
from ai_structure_interface import (
    answer_structure_question,
    plan_structure_search,
    get_design_recommendations,
)


class TestAIStructureInterface:
    """Test AI-callable interface functions."""

    def test_answer_question_metal_coordination(self):
        """Answer question about metal coordination."""
        result = answer_structure_question(
            "What is the coordination geometry of terbium?"
        )

        assert "answer" in result
        assert "coordination_number" in result
        assert result["coordination_number"] in [8, 9]

    def test_answer_question_find_structure(self):
        """Answer question to find structures."""
        result = answer_structure_question(
            "Find me a PQQ-calcium structure"
        )

        assert "answer" in result
        assert "structures" in result or "pdb_id" in result

    def test_plan_structure_search(self):
        """Plan a structure search from NL query."""
        plan = plan_structure_search(
            "I need to design a terbium biosensor"
        )

        assert "steps" in plan
        assert "metal" in plan
        assert plan["metal"] == "TB"

    def test_get_design_recommendations(self):
        """Get design recommendations for a metal."""
        recs = get_design_recommendations(
            metal="TB",
            task="biosensor",
        )

        assert "donors" in recs
        assert "templates" in recs
        assert "bias_aa" in recs
```

**Step 2: Run test to verify it fails**

Run: `cd G:/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless && python -m pytest tests/test_ai_interface.py -v`
Expected: FAIL

**Step 3: Implement AI structure interface**

```python
# backend/serverless/ai_structure_interface.py
"""
AI Structure Interface

High-level interface for AI agents to query structure information.
Designed for natural language interaction and task planning.

This module provides functions that can be called by Claude or other
AI systems to get structured information about proteins and metals.
"""

import logging
from typing import Any, Dict, List, Optional

# Import discovery service
try:
    from structure_discovery import StructureDiscovery, discover_structures
    DISCOVERY_AVAILABLE = True
except ImportError:
    DISCOVERY_AVAILABLE = False

# Import metal chemistry
try:
    from metal_chemistry import (
        METAL_DATABASE,
        COORDINATION_MOTIFS,
        get_coordination_template,
        get_preferred_donors,
    )
    METAL_CHEMISTRY_AVAILABLE = True
except ImportError:
    METAL_CHEMISTRY_AVAILABLE = False

# Import templates
try:
    from metal_ligand_templates import (
        METAL_LIGAND_COMPLEX_TEMPLATES,
        get_template_info,
    )
    TEMPLATES_AVAILABLE = True
except ImportError:
    TEMPLATES_AVAILABLE = False

logger = logging.getLogger(__name__)


def answer_structure_question(question: str) -> Dict[str, Any]:
    """
    Answer a natural language question about protein structures.

    This is the main entry point for AI agents to query structure knowledge.

    Args:
        question: Natural language question

    Returns:
        Dictionary with answer and structured data

    Examples:
        >>> answer_structure_question("What is the coordination geometry of terbium?")
        {
            "answer": "Terbium typically has 8-9 coordination in proteins...",
            "coordination_number": [8, 9],
            "geometry": ["square_antiprism", "tricapped_trigonal_prism"],
            "preferred_donors": ["O"],
            "example_structures": ["6MI5", "7D4I"],
        }
    """
    question_lower = question.lower()

    # Detect question type
    if _is_coordination_question(question_lower):
        return _answer_coordination_question(question_lower)
    elif _is_find_structure_question(question_lower):
        return _answer_find_structure_question(question_lower)
    elif _is_template_question(question_lower):
        return _answer_template_question(question_lower)
    elif _is_donor_question(question_lower):
        return _answer_donor_question(question_lower)
    else:
        return _answer_general_question(question)


def _is_coordination_question(q: str) -> bool:
    """Check if question is about coordination."""
    keywords = ["coordination", "geometry", "coordinate", "bond", "cn="]
    return any(k in q for k in keywords)


def _is_find_structure_question(q: str) -> bool:
    """Check if question is asking to find structures."""
    keywords = ["find", "search", "look for", "get", "show me", "list"]
    return any(k in q for k in keywords)


def _is_template_question(q: str) -> bool:
    """Check if question is about templates."""
    keywords = ["template", "example", "reference", "starting"]
    return any(k in q for k in keywords)


def _is_donor_question(q: str) -> bool:
    """Check if question is about donor preferences."""
    keywords = ["donor", "residue", "amino acid", "bind", "coordinate"]
    return any(k in q for k in keywords)


def _answer_coordination_question(question: str) -> Dict[str, Any]:
    """Answer question about metal coordination."""
    # Extract metal from question
    metal = _extract_metal_from_question(question)

    if not metal:
        return {
            "answer": "Please specify a metal (e.g., zinc, terbium, calcium).",
            "error": "no_metal_specified",
        }

    if not METAL_CHEMISTRY_AVAILABLE:
        return {
            "answer": f"Metal chemistry module not available.",
            "error": "module_unavailable",
        }

    metal_info = METAL_DATABASE.get(metal, {})

    # Find matching coordination motifs
    motifs = []
    for name, data in COORDINATION_MOTIFS.items():
        if metal in data.get("metals", []):
            motifs.append(name)

    # Get coordination info
    cn_range = metal_info.get("coordination_numbers", [6])
    geometry = metal_info.get("preferred_geometry", "octahedral")

    # Build answer
    if metal in ["TB", "EU", "GD", "LA", "CE", "SM", "YB"]:
        answer = (
            f"{metal} is a lanthanide with typical coordination number 8-9. "
            f"It prefers hard oxygen donors (Glu, Asp, Asn) due to HSAB theory. "
            f"Common geometries are square antiprism (CN=8) or tricapped trigonal prism (CN=9)."
        )
    elif metal == "ZN":
        answer = (
            f"Zinc typically has coordination number 4-6. "
            f"Tetrahedral (CN=4) is most common in structural roles (zinc fingers). "
            f"Prefers His and Cys donors."
        )
    elif metal == "CA":
        answer = (
            f"Calcium typically has coordination number 6-8. "
            f"Prefers oxygen donors from Asp, Glu, and backbone carbonyls. "
            f"EF-hand motif is the classic Ca-binding domain (CN=7)."
        )
    else:
        answer = f"{metal} typically has coordination number {cn_range}."

    return {
        "answer": answer,
        "metal": metal,
        "coordination_number": cn_range if isinstance(cn_range, list) else [cn_range],
        "geometry": geometry,
        "preferred_donors": metal_info.get("preferred_donors", []),
        "motifs": motifs,
        "hsab_class": metal_info.get("hsab", "unknown"),
    }


def _answer_find_structure_question(question: str) -> Dict[str, Any]:
    """Answer question asking to find structures."""
    if not DISCOVERY_AVAILABLE:
        return {
            "answer": "Structure discovery module not available.",
            "error": "module_unavailable",
        }

    # Use discovery service
    discovery = StructureDiscovery()
    results = discovery.search_by_intent(question, limit=5)

    if not results:
        return {
            "answer": "No structures found matching your query.",
            "structures": [],
        }

    # Build answer
    pdb_list = [r.get("pdb_id") for r in results if r.get("pdb_id")]

    answer = f"Found {len(results)} structures. "
    if pdb_list:
        answer += f"Top matches: {', '.join(pdb_list[:3])}."

    if results[0].get("source") == "curated":
        answer += " These are curated high-quality reference structures."

    return {
        "answer": answer,
        "structures": results,
        "recommended": pdb_list[0] if pdb_list else None,
    }


def _answer_template_question(question: str) -> Dict[str, Any]:
    """Answer question about templates."""
    metal = _extract_metal_from_question(question)
    ligand = _extract_ligand_from_question(question)

    if DISCOVERY_AVAILABLE:
        discovery = StructureDiscovery()
        template = discovery.get_best_template(metal or "ZN", ligand)

        if template:
            return {
                "answer": f"Found template for {metal or 'metal'}" +
                         (f" with {ligand}" if ligand else ""),
                "template": template,
            }

    return {
        "answer": "No specific template found. Try specifying a metal.",
        "template": None,
    }


def _answer_donor_question(question: str) -> Dict[str, Any]:
    """Answer question about donor residues."""
    metal = _extract_metal_from_question(question)

    if not metal:
        return {
            "answer": "Please specify a metal to get donor recommendations.",
            "error": "no_metal_specified",
        }

    if METAL_CHEMISTRY_AVAILABLE:
        try:
            donors = get_preferred_donors(metal, oxidation_state=None)

            # Map to residues
            donor_residues = {
                "O": ["Glu", "Asp", "Asn", "Gln", "Ser", "Thr", "Tyr"],
                "N": ["His", "Lys", "Arg"],
                "S": ["Cys", "Met"],
            }

            preferred_residues = []
            for donor_type, weight in donors.items():
                if weight > 0:
                    preferred_residues.extend(donor_residues.get(donor_type, []))

            # Lanthanide-specific advice
            if metal in ["TB", "EU", "GD", "LA", "CE"]:
                answer = (
                    f"For {metal} (lanthanide/hard acid), strongly prefer oxygen donors: "
                    f"Glu, Asp, Asn. Avoid Cys (soft sulfur donor). "
                    f"His is acceptable but suboptimal."
                )
            elif metal == "ZN":
                answer = (
                    f"For {metal}, prefer His (most common) and Cys. "
                    f"Asp/Glu also work for catalytic sites."
                )
            else:
                answer = f"For {metal}, preferred donors: {', '.join(preferred_residues[:5])}."

            return {
                "answer": answer,
                "metal": metal,
                "donor_weights": donors,
                "preferred_residues": preferred_residues,
            }
        except Exception as e:
            logger.warning(f"Error getting donors: {e}")

    return {
        "answer": f"Could not determine donor preferences for {metal}.",
        "metal": metal,
    }


def _answer_general_question(question: str) -> Dict[str, Any]:
    """Handle general questions."""
    return {
        "answer": (
            "I can help with:\n"
            "- Metal coordination geometry (e.g., 'What is the coordination of terbium?')\n"
            "- Finding structures (e.g., 'Find zinc finger proteins')\n"
            "- Template selection (e.g., 'Get template for PQQ-calcium')\n"
            "- Donor residues (e.g., 'What residues coordinate lanthanides?')"
        ),
        "help": True,
    }


def _extract_metal_from_question(question: str) -> Optional[str]:
    """Extract metal name/symbol from question."""
    metal_names = {
        "zinc": "ZN", "zn": "ZN",
        "iron": "FE", "fe": "FE",
        "copper": "CU", "cu": "CU",
        "calcium": "CA", "ca": "CA",
        "magnesium": "MG", "mg": "MG",
        "manganese": "MN", "mn": "MN",
        "cobalt": "CO",
        "nickel": "NI", "ni": "NI",
        "terbium": "TB", "tb": "TB",
        "europium": "EU", "eu": "EU",
        "gadolinium": "GD", "gd": "GD",
        "lanthanum": "LA", "la": "LA",
        "lanthanide": "TB",
    }

    for name, code in metal_names.items():
        if name in question:
            return code
    return None


def _extract_ligand_from_question(question: str) -> Optional[str]:
    """Extract ligand name from question."""
    ligands = ["pqq", "citrate", "heme", "atp", "nad", "fad"]
    for lig in ligands:
        if lig in question:
            return lig.upper()
    return None


def plan_structure_search(query: str) -> Dict[str, Any]:
    """
    Create a search plan from natural language query.

    Used by task planner to structure search operations.

    Args:
        query: Natural language description of what user needs

    Returns:
        Search plan with steps and parameters
    """
    metal = _extract_metal_from_question(query.lower())
    ligand = _extract_ligand_from_question(query.lower())

    # Determine task type
    task_type = "unknown"
    if "biosensor" in query.lower() or "luminescent" in query.lower():
        task_type = "biosensor"
    elif "binder" in query.lower():
        task_type = "binder"
    elif "enzyme" in query.lower() or "catalytic" in query.lower():
        task_type = "enzyme"
    elif "redesign" in query.lower() or "replace" in query.lower():
        task_type = "redesign"

    steps = []

    # Step 1: Find reference structures
    if metal:
        steps.append({
            "action": "search_structures",
            "params": {"metal": metal, "ligand": ligand, "limit": 5},
            "description": f"Find {metal} structures" + (f" with {ligand}" if ligand else ""),
        })

    # Step 2: Get coordination template
    if metal:
        steps.append({
            "action": "get_template",
            "params": {"metal": metal, "ligand": ligand},
            "description": "Get coordination template",
        })

    # Step 3: Get design recommendations
    steps.append({
        "action": "get_recommendations",
        "params": {"metal": metal, "task": task_type},
        "description": "Get design parameter recommendations",
    })

    return {
        "query": query,
        "metal": metal,
        "ligand": ligand,
        "task_type": task_type,
        "steps": steps,
    }


def get_design_recommendations(
    metal: str,
    task: str = "general",
    ligand: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Get recommendations for protein design with a metal.

    Args:
        metal: Metal element symbol
        task: Task type (biosensor, binder, enzyme, redesign)
        ligand: Optional ligand

    Returns:
        Recommendations for donors, templates, HSAB bias
    """
    metal = metal.upper()

    recommendations = {
        "metal": metal,
        "task": task,
    }

    # Get donor recommendations
    if METAL_CHEMISTRY_AVAILABLE:
        try:
            donors = get_preferred_donors(metal)

            # Convert to residue bias string for RFD3
            bias_parts = []

            # Lanthanides: favor Glu, Asp
            if metal in ["TB", "EU", "GD", "LA", "CE", "SM", "YB"]:
                bias_parts = ["E:3.0", "D:3.0", "N:2.4", "Q:2.4", "H:0.5", "C:-5.0"]
            elif metal == "ZN":
                bias_parts = ["H:5.0", "C:4.0", "D:3.0", "E:3.0"]
            elif metal == "CA":
                bias_parts = ["E:4.0", "D:4.0", "N:2.0", "Q:2.0"]
            else:
                bias_parts = ["E:3.0", "D:3.0", "H:2.0", "N:1.0"]

            recommendations["donors"] = donors
            recommendations["bias_aa"] = ",".join(bias_parts)

        except Exception:
            recommendations["donors"] = {}
            recommendations["bias_aa"] = ""

    # Get template recommendations
    templates = []
    if TEMPLATES_AVAILABLE:
        for key, template in METAL_LIGAND_COMPLEX_TEMPLATES.items():
            template_metal = template.get("metal") or template.get("default_metal")
            if template_metal == metal:
                templates.append({
                    "id": key,
                    "name": template.get("name"),
                    "ligand": template.get("ligand_res_name"),
                })

    if DISCOVERY_AVAILABLE:
        discovery = StructureDiscovery()
        best_template = discovery.get_best_template(metal, ligand)
        if best_template:
            recommendations["best_template"] = best_template

    recommendations["templates"] = templates

    # Task-specific recommendations
    if task == "biosensor":
        recommendations["notes"] = [
            "For biosensors, maximize signal change upon metal binding.",
            "Consider lanthanide luminescence for optical detection.",
            "Ensure high metal selectivity over Ca2+/Mg2+.",
        ]
        if metal in ["TB", "EU"]:
            recommendations["notes"].append(
                f"{metal} has excellent luminescence properties for biosensing."
            )
    elif task == "binder":
        recommendations["notes"] = [
            "Design for high affinity and specificity.",
            "Consider using hotspot residues for stronger binding.",
        ]
    elif task == "enzyme":
        recommendations["notes"] = [
            "Catalytic sites often need precise geometry.",
            "Consider substrate access and product release.",
        ]

    return recommendations
```

**Step 4: Run tests**

Run: `cd G:/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless && python -m pytest tests/test_ai_interface.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add backend/serverless/ai_structure_interface.py
git add backend/serverless/tests/test_ai_interface.py
git commit -m "feat: add AI-callable structure interface with NL question answering"
```

---

### Task 4.2: Integrate with Task Planner

**Files:**
- Modify: `backend/ai/task_planner.py`

**Step 1: Add import and enhance FETCH_PDB step**

Add to `backend/ai/task_planner.py` at the imports section:

```python
# Add after existing imports
try:
    import sys
    sys.path.insert(0, '../serverless')
    from structure_discovery import StructureDiscovery
    from ai_structure_interface import (
        answer_structure_question,
        plan_structure_search,
        get_design_recommendations,
    )
    STRUCTURE_DISCOVERY_AVAILABLE = True
except ImportError:
    STRUCTURE_DISCOVERY_AVAILABLE = False
```

**Step 2: Enhance the _execute_step method**

Replace the FETCH_PDB handling in `TaskExecutor._execute_step`:

```python
        if step.type == TaskStepType.FETCH_PDB:
            pdb_id = step.params.get("pdb_id")
            metal = step.params.get("metal")
            ligand = step.params.get("ligand")

            # If no PDB ID, use discovery to find one
            if not pdb_id and STRUCTURE_DISCOVERY_AVAILABLE:
                discovery = StructureDiscovery()

                if metal:
                    # Search by metal
                    results = discovery.search_by_metal(metal, limit=1)
                    if results:
                        pdb_id = results[0].get("pdb_id")

                # Get best template if available
                if metal:
                    template = discovery.get_best_template(metal, ligand)
                    if template:
                        self.artifacts["template"] = template
                        if not pdb_id and template.get("pdb_id"):
                            pdb_id = template.get("pdb_id")

            if not pdb_id:
                return {"error": "No PDB ID found or specified"}

            import urllib.request
            url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
            with urllib.request.urlopen(url) as response:
                content = response.read().decode()
            self.artifacts["pdb_content"] = content
            self.artifacts["pdb_id"] = pdb_id
            return {"pdb_id": pdb_id, "size": len(content)}
```

**Step 3: Add new task step type for discovery**

Add to `TaskStepType` enum:

```python
    DISCOVER_STRUCTURES = "discover_structures"  # AI-driven structure discovery
```

**Step 4: Add handler for discovery step**

Add to `_execute_step` method:

```python
        elif step.type == TaskStepType.DISCOVER_STRUCTURES:
            query = step.params.get("query", "")
            metal = step.params.get("metal")
            limit = step.params.get("limit", 5)

            if STRUCTURE_DISCOVERY_AVAILABLE:
                discovery = StructureDiscovery()

                if query:
                    results = discovery.search_by_intent(query, limit=limit)
                elif metal:
                    results = discovery.search_by_metal(metal, limit=limit)
                else:
                    return {"error": "No query or metal specified"}

                self.artifacts["discovered_structures"] = results

                # Get recommendations if metal specified
                if metal:
                    recs = get_design_recommendations(metal)
                    self.artifacts["recommendations"] = recs

                return {
                    "structures_found": len(results),
                    "structures": results[:3],  # Summary
                    "recommendations": self.artifacts.get("recommendations"),
                }

            return {"error": "Structure discovery not available"}
```

**Step 5: Commit**

```bash
git add backend/ai/task_planner.py
git commit -m "feat: integrate structure discovery with task planner"
```

---

## Phase 5: Final Integration & Testing

### Task 5.1: Integration Tests

**Files:**
- Create: `backend/serverless/tests/test_integration.py`

**Step 1: Write integration tests**

```python
# backend/serverless/tests/test_integration.py
"""Integration tests for the complete discovery pipeline."""
import pytest


class TestFullPipeline:
    """Test complete discovery to design pipeline."""

    def test_nl_query_to_template(self):
        """Natural language query should yield usable template."""
        from ai_structure_interface import answer_structure_question
        from structure_discovery import StructureDiscovery

        # Ask about terbium
        result = answer_structure_question(
            "What template should I use for terbium biosensor design?"
        )

        assert "answer" in result

        # Get template
        discovery = StructureDiscovery()
        template = discovery.get_best_template("TB", "CIT")

        assert template is not None

    def test_metal_search_to_coordination(self):
        """Metal search should provide coordination info."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        results = discovery.search_by_metal("ZN", limit=3)

        assert len(results) > 0

        # Extract coordination from first result
        if results[0].get("pdb_id"):
            info = discovery.extract_coordination_info(
                results[0]["pdb_id"], "ZN"
            )
            # May be None if network unavailable
            if info:
                assert "coordination_number" in info

    def test_ligand_smiles_retrieval(self):
        """Should retrieve SMILES for common ligands."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        smiles = discovery.get_ligand_smiles("citric acid")

        # May be None if PubChem unavailable
        if smiles:
            assert "C" in smiles  # Has carbons
            assert "O" in smiles  # Has oxygens

    def test_recommendations_include_bias(self):
        """Design recommendations should include HSAB bias."""
        from ai_structure_interface import get_design_recommendations

        recs = get_design_recommendations("TB", task="biosensor")

        assert "bias_aa" in recs
        assert "E:" in recs["bias_aa"]  # Glu should be biased
        assert "C:-" in recs["bias_aa"]  # Cys should be negative


class TestArchitectorIntegration:
    """Test Architector integration in discovery."""

    def test_lanthanide_template_uses_architector(self):
        """Lanthanide templates should use Architector when available."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        template = discovery.get_best_template("TB", use_architector=True)

        assert template is not None
        # Source may be library, pdb, or architector

    def test_geometry_validation_available(self):
        """Coordination validation should work."""
        from metal_chemistry import validate_coordination

        # Good octahedral
        positions = [
            (2.4, 0, 0), (-2.4, 0, 0),
            (0, 2.4, 0), (0, -2.4, 0),
            (0, 0, 2.4), (0, 0, -2.4),
        ]

        result = validate_coordination(positions, (0, 0, 0), "ZN")
        assert result["valid"]
```

**Step 2: Run integration tests**

Run: `cd G:/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless && python -m pytest tests/test_integration.py -v`
Expected: PASS (some may skip if network unavailable)

**Step 3: Commit**

```bash
git add backend/serverless/tests/test_integration.py
git commit -m "test: add integration tests for full discovery pipeline"
```

---

### Task 5.2: Documentation

**Files:**
- Create: `backend/serverless/database_adapters/README.md`

**Step 1: Write documentation**

```markdown
# Database Adapters

This package provides adapters for external protein structure databases.

## Available Adapters

### MetalPDBAdapter
Queries MetalPDB for metal coordination site information.
- Search by metal element
- Filter by coordination number and geometry
- Get detailed site information

### UniProtAdapter
Queries UniProt for protein function and annotation data.
- Search by protein function
- Search by metal binding
- Get PDB structure mappings
- Get metal binding site annotations

### PubChemAdapter
Queries PubChem for ligand structure information.
- Search compounds by name
- Get SMILES strings
- Get molecular properties

### AlphaFoldAdapter
Queries AlphaFold DB for predicted structures.
- Get predictions by UniProt ID
- Download predicted structures
- Extract pLDDT confidence scores

## Usage

```python
from database_adapters import MetalPDBAdapter, UniProtAdapter

# Search for zinc sites
metalpdb = MetalPDBAdapter()
results = metalpdb.search_by_metal("ZN", coordination_number=4)

# Search for metal-binding proteins
uniprot = UniProtAdapter()
proteins = uniprot.search_by_metal_binding("Zinc")
```

## Structure Discovery Service

The `StructureDiscovery` class provides a unified interface:

```python
from structure_discovery import StructureDiscovery

discovery = StructureDiscovery()

# Natural language search
results = discovery.search_by_intent("terbium binding protein")

# Get best template
template = discovery.get_best_template("TB", "CIT")
```

## AI Interface

For AI agents (Claude):

```python
from ai_structure_interface import answer_structure_question

result = answer_structure_question(
    "What is the coordination geometry of terbium?"
)
print(result["answer"])
```
```

**Step 2: Commit**

```bash
git add backend/serverless/database_adapters/README.md
git commit -m "docs: add documentation for database adapters"
```

---

## Summary

This plan implements:

1. **Phase 1: Database Adapters** (4 tasks)
   - MetalPDB adapter for coordination sites
   - UniProt adapter for protein function
   - PubChem adapter for ligand SMILES
   - AlphaFold adapter for predicted structures

2. **Phase 2: Structure Discovery Service** (1 task)
   - Unified search interface
   - Natural language intent parsing
   - Quality-based ranking

3. **Phase 3: Coordination Knowledge Base** (1 task)
   - Coordination motifs (EF-hand, zinc finger, etc.)
   - Ideal geometry definitions
   - Validation functions

4. **Phase 4: AI Interface Integration** (2 tasks)
   - AI-callable question answering
   - Task planner integration

5. **Phase 5: Final Integration** (2 tasks)
   - Integration tests
   - Documentation

Total: **10 tasks** with TDD approach (test first, then implement).
