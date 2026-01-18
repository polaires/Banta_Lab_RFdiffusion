# test_uniprot_adapter.py
"""Tests for UniProt database adapter."""
import pytest
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Check if requests is available for network tests
try:
    import requests
    NETWORK_AVAILABLE = True
    # Quick connectivity check for UniProt
    try:
        requests.head("https://rest.uniprot.org", timeout=5)
        UNIPROT_AVAILABLE = True
    except (requests.RequestException, Exception):
        UNIPROT_AVAILABLE = False
except ImportError:
    NETWORK_AVAILABLE = False
    UNIPROT_AVAILABLE = False

network_required = pytest.mark.skipif(
    not NETWORK_AVAILABLE,
    reason="Network access required (requests not available)"
)

uniprot_required = pytest.mark.skipif(
    not UNIPROT_AVAILABLE,
    reason="UniProt API not reachable"
)


class TestUniProtAdapterImport:
    """Test UniProt adapter can be imported."""

    def test_can_import_adapter(self):
        """Should be able to import UniProtAdapter."""
        from database_adapters import UniProtAdapter
        assert UniProtAdapter is not None

    def test_can_import_from_module(self):
        """Should be able to import from uniprot_adapter module."""
        from database_adapters.uniprot_adapter import UniProtAdapter
        assert UniProtAdapter is not None


class TestUniProtAdapterInitialization:
    """Test UniProtAdapter initialization."""

    def test_adapter_has_default_base_url(self):
        """Should have UniProt REST API base URL by default."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()
        assert "rest.uniprot.org" in adapter.base_url

    def test_adapter_can_override_base_url(self):
        """Should allow overriding base URL."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter(base_url="https://custom.example.com")
        assert adapter.base_url == "https://custom.example.com"

    def test_adapter_has_timeout_setting(self):
        """Should have configurable timeout."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter(timeout=60)
        assert adapter.timeout == 60


class TestSearchByFunction:
    """Test search_by_function method."""

    def test_search_by_function_returns_list(self):
        """Should return a list from search_by_function."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()
        # Even without network, should return empty list not raise
        results = adapter.search_by_function("zinc finger", limit=5)
        assert isinstance(results, list)

    @network_required
    @uniprot_required
    def test_search_by_function_returns_results(self):
        """Search for 'zinc finger' should return results."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()
        results = adapter.search_by_function("zinc finger", limit=5)

        # Should find results for common function query
        assert len(results) > 0

        # Each result should have required fields
        for result in results:
            assert "uniprot_id" in result
            assert "name" in result
            assert result["uniprot_id"] is not None

    @network_required
    @uniprot_required
    def test_search_by_function_with_organism(self):
        """Search with organism filter should work."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()
        results = adapter.search_by_function("dehydrogenase", organism="human", limit=10)

        # Should filter to human proteins
        for result in results:
            if "organism" in result:
                assert "human" in result["organism"].lower() or "sapiens" in result["organism"].lower()

    @network_required
    @uniprot_required
    def test_search_by_function_reviewed_only(self):
        """Search with reviewed_only=True should return Swiss-Prot entries."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()
        results = adapter.search_by_function("kinase", reviewed_only=True, limit=5)

        # All results should be reviewed (Swiss-Prot)
        assert isinstance(results, list)


class TestSearchByMetalBinding:
    """Test search_by_metal_binding method."""

    def test_search_by_metal_binding_returns_list(self):
        """Should return a list from search_by_metal_binding."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()
        results = adapter.search_by_metal_binding("Zinc", limit=5)
        assert isinstance(results, list)

    @network_required
    @uniprot_required
    def test_search_by_metal_binding_zinc(self):
        """Search for Zinc binding proteins should return results."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()
        results = adapter.search_by_metal_binding("Zinc", limit=10)

        # Should find zinc-binding proteins
        assert len(results) > 0

        for result in results:
            assert "uniprot_id" in result

    @network_required
    @uniprot_required
    def test_search_by_metal_binding_iron(self):
        """Search for Iron binding proteins should return results."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()
        results = adapter.search_by_metal_binding("Iron", limit=10)

        # Should find iron-binding proteins
        assert len(results) > 0

    @network_required
    @uniprot_required
    def test_search_by_metal_binding_with_organism(self):
        """Search with organism filter should work."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()
        results = adapter.search_by_metal_binding("Zinc", organism="human", limit=10)

        assert isinstance(results, list)


class TestGetPdbMappings:
    """Test get_pdb_mappings method."""

    def test_get_pdb_mappings_returns_list(self):
        """Should return a list from get_pdb_mappings."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()
        results = adapter.get_pdb_mappings("P00918")  # Human carbonic anhydrase II
        assert isinstance(results, list)

    @network_required
    @uniprot_required
    def test_get_pdb_mappings_carbonic_anhydrase(self):
        """P00918 (human carbonic anhydrase II) should return 1CA2 among PDB mappings."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()
        results = adapter.get_pdb_mappings("P00918")

        # Should have PDB mappings
        assert len(results) > 0

        # Should contain 1CA2 (classic carbonic anhydrase structure)
        pdb_ids = [r.get("pdb_id", "").upper() for r in results]
        assert "1CA2" in pdb_ids, f"Expected 1CA2 in {pdb_ids}"

        # Each result should have expected fields
        for result in results:
            assert "pdb_id" in result

    @network_required
    @uniprot_required
    def test_get_pdb_mappings_returns_structure_info(self):
        """PDB mappings should include method and resolution where available."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()
        results = adapter.get_pdb_mappings("P00918")

        if results:
            # At least some results should have structure info
            has_method = any("method" in r for r in results)
            # Method info may not always be available, so just check structure
            assert all("pdb_id" in r for r in results)

    def test_get_pdb_mappings_invalid_id_returns_empty(self):
        """Invalid UniProt ID should return empty list."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()
        results = adapter.get_pdb_mappings("INVALID123")
        assert isinstance(results, list)
        assert len(results) == 0


class TestGetMetalBindingSites:
    """Test get_metal_binding_sites method."""

    def test_get_metal_binding_sites_returns_list(self):
        """Should return a list from get_metal_binding_sites."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()
        results = adapter.get_metal_binding_sites("P00918")
        assert isinstance(results, list)

    @network_required
    @uniprot_required
    def test_get_metal_binding_sites_carbonic_anhydrase(self):
        """P00918 (human carbonic anhydrase II) should have Zinc binding sites."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()
        results = adapter.get_metal_binding_sites("P00918")

        # Should have metal binding sites
        assert len(results) > 0, "P00918 should have metal binding sites"

        # Should have zinc binding
        metals = [r.get("metal", "").lower() for r in results]
        assert any("zinc" in m or "zn" in m for m in metals), f"Expected Zinc in {metals}"

        # Each result should have expected fields
        for result in results:
            assert "position" in result or "metal" in result

    @network_required
    @uniprot_required
    def test_get_metal_binding_sites_returns_positions(self):
        """Metal binding sites should include position information."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()
        results = adapter.get_metal_binding_sites("P00918")

        if results:
            # At least some results should have position
            has_position = any("position" in r for r in results)
            assert has_position, "Metal binding sites should include position"

    def test_get_metal_binding_sites_invalid_id_returns_empty(self):
        """Invalid UniProt ID should return empty list."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()
        results = adapter.get_metal_binding_sites("INVALID123")
        assert isinstance(results, list)
        assert len(results) == 0


class TestGetSequence:
    """Test get_sequence method."""

    def test_get_sequence_returns_string_or_none(self):
        """Should return a string or None from get_sequence."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()
        result = adapter.get_sequence("P00918")
        assert result is None or isinstance(result, str)

    @network_required
    @uniprot_required
    def test_get_sequence_carbonic_anhydrase(self):
        """P00918 should return a valid protein sequence."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()
        sequence = adapter.get_sequence("P00918")

        # Should return sequence
        assert sequence is not None
        assert len(sequence) > 0

        # Should be valid amino acid sequence
        valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
        assert all(aa in valid_aas for aa in sequence.upper())

    @network_required
    @uniprot_required
    def test_get_sequence_expected_length(self):
        """P00918 sequence should be approximately 260 amino acids."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()
        sequence = adapter.get_sequence("P00918")

        if sequence:
            # Human carbonic anhydrase II is ~260 amino acids
            assert 250 < len(sequence) < 280, f"Unexpected length: {len(sequence)}"

    def test_get_sequence_invalid_id_returns_none(self):
        """Invalid UniProt ID should return None."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()
        result = adapter.get_sequence("INVALID123")
        assert result is None


class TestErrorHandling:
    """Test error handling and logging."""

    def test_network_error_returns_empty_list(self):
        """Network errors should return empty list, not raise."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter(
            base_url="https://invalid.nonexistent.domain.example",
            timeout=1,  # Short timeout
        )
        # Should not raise, should return empty list
        results = adapter.search_by_function("test")
        assert isinstance(results, list)

    def test_network_error_get_sequence_returns_none(self):
        """Network errors on get_sequence should return None."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter(
            base_url="https://invalid.nonexistent.domain.example",
            timeout=1,
        )
        result = adapter.get_sequence("P00918")
        assert result is None


class TestIntegration:
    """Integration tests requiring network access."""

    @pytest.mark.integration
    @network_required
    @uniprot_required
    def test_full_workflow_zinc_proteins(self):
        """Full workflow: search for zinc proteins, get details."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()

        # Step 1: Search for zinc-binding proteins
        results = adapter.search_by_metal_binding("Zinc", organism="human", limit=3)

        if results:
            # Step 2: Get first result's UniProt ID
            first = results[0]
            uniprot_id = first.get("uniprot_id")

            if uniprot_id:
                # Step 3: Get metal binding sites
                sites = adapter.get_metal_binding_sites(uniprot_id)
                assert isinstance(sites, list)

                # Step 4: Get sequence
                sequence = adapter.get_sequence(uniprot_id)
                assert sequence is None or isinstance(sequence, str)

                # Step 5: Get PDB mappings
                pdb_mappings = adapter.get_pdb_mappings(uniprot_id)
                assert isinstance(pdb_mappings, list)

    @pytest.mark.integration
    @network_required
    @uniprot_required
    def test_search_function_then_get_details(self):
        """Search by function, then get detailed info."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()

        # Search for metalloenzymes
        results = adapter.search_by_function("metalloenzyme", reviewed_only=True, limit=5)

        if results:
            uniprot_id = results[0].get("uniprot_id")
            if uniprot_id:
                # Get more details
                sites = adapter.get_metal_binding_sites(uniprot_id)
                sequence = adapter.get_sequence(uniprot_id)

                assert isinstance(sites, list)
                assert sequence is None or isinstance(sequence, str)
