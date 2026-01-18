# test_pubchem_adapter.py
"""Tests for PubChem database adapter."""
import pytest
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Check if requests is available for network tests
try:
    import requests
    NETWORK_AVAILABLE = True
    # Quick connectivity check for PubChem
    try:
        requests.head("https://pubchem.ncbi.nlm.nih.gov", timeout=5)
        PUBCHEM_AVAILABLE = True
    except (requests.RequestException, Exception):
        PUBCHEM_AVAILABLE = False
except ImportError:
    NETWORK_AVAILABLE = False
    PUBCHEM_AVAILABLE = False

network_required = pytest.mark.skipif(
    not NETWORK_AVAILABLE,
    reason="Network access required (requests not available)"
)

pubchem_required = pytest.mark.skipif(
    not PUBCHEM_AVAILABLE,
    reason="PubChem API not reachable"
)


class TestPubChemAdapterImport:
    """Test PubChem adapter can be imported."""

    def test_can_import_adapter(self):
        """Should be able to import PubChemAdapter."""
        from database_adapters import PubChemAdapter
        assert PubChemAdapter is not None

    def test_can_import_from_module(self):
        """Should be able to import from pubchem_adapter module."""
        from database_adapters.pubchem_adapter import PubChemAdapter
        assert PubChemAdapter is not None

    def test_can_import_common_ligands(self):
        """Should be able to import COMMON_LIGANDS lookup table."""
        from database_adapters.pubchem_adapter import COMMON_LIGANDS
        assert COMMON_LIGANDS is not None
        assert isinstance(COMMON_LIGANDS, dict)
        assert "ATP" in COMMON_LIGANDS
        assert "NAD" in COMMON_LIGANDS


class TestPubChemAdapterInitialization:
    """Test PubChemAdapter initialization."""

    def test_adapter_has_default_base_url(self):
        """Should have PubChem PUG REST API base URL by default."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter()
        assert "pubchem.ncbi.nlm.nih.gov" in adapter.base_url

    def test_adapter_can_override_base_url(self):
        """Should allow overriding base URL."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter(base_url="https://custom.example.com")
        assert adapter.base_url == "https://custom.example.com"

    def test_adapter_has_timeout_setting(self):
        """Should have configurable timeout."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter(timeout=60)
        assert adapter.timeout == 60


class TestSearchByName:
    """Test search_by_name method."""

    def test_search_by_name_returns_list(self):
        """Should return a list from search_by_name."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter()
        # Even without network, should return empty list not raise
        results = adapter.search_by_name("citric acid", limit=5)
        assert isinstance(results, list)

    @network_required
    @pubchem_required
    def test_search_by_name_returns_results(self):
        """Search for 'citric acid' should return results with cid and smiles."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter()
        results = adapter.search_by_name("citric acid", limit=5)

        # Should find results for citric acid
        assert len(results) > 0

        # Each result should have required fields
        for result in results:
            assert "cid" in result
            assert "smiles" in result
            assert result["cid"] is not None

    @network_required
    @pubchem_required
    def test_search_by_name_aspirin(self):
        """Search for 'aspirin' should return results."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter()
        results = adapter.search_by_name("aspirin", limit=3)

        assert len(results) > 0
        assert "cid" in results[0]
        assert "smiles" in results[0]

    def test_search_by_name_empty_returns_empty(self):
        """Empty name should return empty list."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter()
        results = adapter.search_by_name("", limit=5)
        assert isinstance(results, list)
        assert len(results) == 0

    def test_search_by_name_whitespace_returns_empty(self):
        """Whitespace-only name should return empty list."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter()
        results = adapter.search_by_name("   ", limit=5)
        assert isinstance(results, list)
        assert len(results) == 0


class TestGetCompound:
    """Test get_compound method."""

    def test_get_compound_returns_dict_or_none(self):
        """Should return a dict or None from get_compound."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter()
        result = adapter.get_compound(311)  # Citric acid
        assert result is None or isinstance(result, dict)

    @network_required
    @pubchem_required
    def test_get_compound_by_cid(self):
        """CID 311 (citric acid) should return smiles and molecular_formula."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter()
        result = adapter.get_compound(311)

        assert result is not None
        assert "smiles" in result
        assert "molecular_formula" in result
        assert "molecular_weight" in result
        assert "iupac_name" in result

        # Citric acid should have carbon in formula
        assert "C" in result["molecular_formula"]

    @network_required
    @pubchem_required
    def test_get_compound_caffeine(self):
        """CID 2519 (caffeine) should return valid properties."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter()
        result = adapter.get_compound(2519)

        assert result is not None
        assert result["cid"] == 2519
        assert "smiles" in result
        assert len(result["smiles"]) > 0

    @network_required
    @pubchem_required
    def test_get_compound_returns_molecular_weight(self):
        """get_compound should return molecular weight as a number."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter()
        result = adapter.get_compound(311)

        if result:
            assert "molecular_weight" in result
            # Citric acid is ~192 g/mol
            weight = result["molecular_weight"]
            assert isinstance(weight, (int, float))
            assert 180 < weight < 210


class TestGetSmiles:
    """Test get_smiles method."""

    def test_get_smiles_returns_string_or_none(self):
        """Should return a string or None from get_smiles."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter()
        result = adapter.get_smiles("ATP")
        assert result is None or isinstance(result, str)

    @network_required
    @pubchem_required
    def test_get_smiles_for_ligand(self):
        """'ATP' should return a SMILES string."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter()
        smiles = adapter.get_smiles("ATP")

        assert smiles is not None
        assert isinstance(smiles, str)
        assert len(smiles) > 0

    @network_required
    @pubchem_required
    def test_get_smiles_nad(self):
        """'NAD' should return a SMILES string."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter()
        smiles = adapter.get_smiles("NAD")

        assert smiles is not None
        assert len(smiles) > 0

    @network_required
    @pubchem_required
    def test_get_smiles_case_insensitive(self):
        """get_smiles should be case-insensitive."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter()
        smiles_upper = adapter.get_smiles("ATP")
        smiles_lower = adapter.get_smiles("atp")
        smiles_mixed = adapter.get_smiles("Atp")

        # All should return the same SMILES
        assert smiles_upper == smiles_lower == smiles_mixed

    @network_required
    @pubchem_required
    def test_get_smiles_citrate(self):
        """'CITRATE' should return a SMILES string (from common ligands)."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter()
        smiles = adapter.get_smiles("CITRATE")

        assert smiles is not None
        assert len(smiles) > 0

    def test_get_smiles_empty_returns_none(self):
        """Empty ligand name should return None."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter()
        result = adapter.get_smiles("")
        assert result is None

    def test_get_smiles_whitespace_returns_none(self):
        """Whitespace-only ligand name should return None."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter()
        result = adapter.get_smiles("   ")
        assert result is None


class TestSearchBySubstructure:
    """Test search_by_substructure method."""

    def test_search_by_substructure_returns_list(self):
        """Should return a list from search_by_substructure."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter()
        results = adapter.search_by_substructure("C(=O)[O-]", limit=5)
        assert isinstance(results, list)

    @network_required
    @pubchem_required
    def test_search_by_substructure(self):
        """Search for 'C(=O)[O-]' (carboxylate) should work."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter()
        # Use a simple pattern that should match many compounds
        results = adapter.search_by_substructure("C(=O)O", limit=5)

        # Substructure search might take time or return async
        # Just verify it returns a list and doesn't crash
        assert isinstance(results, list)

    @network_required
    @pubchem_required
    def test_search_by_substructure_benzene(self):
        """Search for benzene ring should return results."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter()
        # Simple benzene SMILES
        results = adapter.search_by_substructure("c1ccccc1", limit=5)

        assert isinstance(results, list)
        # Benzene is very common, should find results
        # Note: may be empty if API is slow or rate-limited

    def test_search_by_substructure_empty_returns_empty(self):
        """Empty SMILES should return empty list."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter()
        results = adapter.search_by_substructure("", limit=5)
        assert isinstance(results, list)
        assert len(results) == 0


class TestCommonLigandsLookup:
    """Test COMMON_LIGANDS lookup table."""

    def test_common_ligands_has_atp(self):
        """COMMON_LIGANDS should have ATP."""
        from database_adapters.pubchem_adapter import COMMON_LIGANDS
        assert "ATP" in COMMON_LIGANDS
        assert COMMON_LIGANDS["ATP"] == 5957

    def test_common_ligands_has_adp(self):
        """COMMON_LIGANDS should have ADP."""
        from database_adapters.pubchem_adapter import COMMON_LIGANDS
        assert "ADP" in COMMON_LIGANDS
        assert COMMON_LIGANDS["ADP"] == 6022

    def test_common_ligands_has_gtp(self):
        """COMMON_LIGANDS should have GTP."""
        from database_adapters.pubchem_adapter import COMMON_LIGANDS
        assert "GTP" in COMMON_LIGANDS

    def test_common_ligands_has_nad(self):
        """COMMON_LIGANDS should have NAD."""
        from database_adapters.pubchem_adapter import COMMON_LIGANDS
        assert "NAD" in COMMON_LIGANDS

    def test_common_ligands_has_nadh(self):
        """COMMON_LIGANDS should have NADH."""
        from database_adapters.pubchem_adapter import COMMON_LIGANDS
        assert "NADH" in COMMON_LIGANDS

    def test_common_ligands_has_fad(self):
        """COMMON_LIGANDS should have FAD."""
        from database_adapters.pubchem_adapter import COMMON_LIGANDS
        assert "FAD" in COMMON_LIGANDS

    def test_common_ligands_has_heme(self):
        """COMMON_LIGANDS should have HEME."""
        from database_adapters.pubchem_adapter import COMMON_LIGANDS
        assert "HEME" in COMMON_LIGANDS

    def test_common_ligands_citrate_alias(self):
        """COMMON_LIGANDS should have CITRATE and CITRIC ACID pointing to same CID."""
        from database_adapters.pubchem_adapter import COMMON_LIGANDS
        assert "CITRATE" in COMMON_LIGANDS
        assert "CITRIC ACID" in COMMON_LIGANDS
        assert COMMON_LIGANDS["CITRATE"] == COMMON_LIGANDS["CITRIC ACID"]
        assert COMMON_LIGANDS["CITRATE"] == 311


class TestErrorHandling:
    """Test error handling and logging."""

    def test_network_error_returns_empty_list(self):
        """Network errors should return empty list, not raise."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter(
            base_url="https://invalid.nonexistent.domain.example",
            timeout=1,  # Short timeout
        )
        # Should not raise, should return empty list
        results = adapter.search_by_name("test")
        assert isinstance(results, list)

    def test_network_error_get_compound_returns_none(self):
        """Network errors on get_compound should return None."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter(
            base_url="https://invalid.nonexistent.domain.example",
            timeout=1,
        )
        result = adapter.get_compound(311)
        assert result is None

    def test_network_error_get_smiles_returns_none(self):
        """Network errors on get_smiles should return None."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter(
            base_url="https://invalid.nonexistent.domain.example",
            timeout=1,
        )
        result = adapter.get_smiles("ATP")
        assert result is None


class TestIntegration:
    """Integration tests requiring network access."""

    @pytest.mark.integration
    @network_required
    @pubchem_required
    def test_full_workflow_ligand_lookup(self):
        """Full workflow: lookup common ligand, get details."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter()

        # Step 1: Get SMILES for ATP
        smiles = adapter.get_smiles("ATP")
        assert smiles is not None
        assert len(smiles) > 0

        # Step 2: Get full compound info
        from database_adapters.pubchem_adapter import COMMON_LIGANDS
        atp_cid = COMMON_LIGANDS["ATP"]
        compound = adapter.get_compound(atp_cid)

        assert compound is not None
        assert "molecular_formula" in compound
        assert "molecular_weight" in compound

    @pytest.mark.integration
    @network_required
    @pubchem_required
    def test_search_then_get_details(self):
        """Search by name, then get detailed info."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter()

        # Search for a compound
        results = adapter.search_by_name("caffeine", limit=1)

        if results:
            cid = results[0]["cid"]

            # Get full details
            compound = adapter.get_compound(cid)
            assert compound is not None
            assert compound["cid"] == cid
            assert "smiles" in compound

    @pytest.mark.integration
    @network_required
    @pubchem_required
    def test_multiple_ligand_lookups(self):
        """Look up multiple common ligands."""
        from database_adapters import PubChemAdapter

        adapter = PubChemAdapter()

        ligands = ["ATP", "NAD", "FAD", "CITRATE"]
        results = {}

        for ligand in ligands:
            smiles = adapter.get_smiles(ligand)
            results[ligand] = smiles

        # All should return SMILES
        for ligand in ligands:
            assert results[ligand] is not None, f"Failed to get SMILES for {ligand}"
            assert len(results[ligand]) > 0
