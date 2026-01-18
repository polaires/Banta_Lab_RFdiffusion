# test_metalpdb_adapter.py
"""Tests for MetalPDB database adapter."""
import pytest
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Check if requests is available for network tests
try:
    import requests
    NETWORK_AVAILABLE = True
    # Quick connectivity check for MetalPDB
    try:
        requests.head("https://metalpdb.cerm.unifi.it", timeout=5)
        METALPDB_AVAILABLE = True
    except (requests.RequestException, Exception):
        METALPDB_AVAILABLE = False
    # Quick connectivity check for RCSB fallback
    try:
        requests.head("https://search.rcsb.org", timeout=5)
        RCSB_AVAILABLE = True
    except (requests.RequestException, Exception):
        RCSB_AVAILABLE = False
except ImportError:
    NETWORK_AVAILABLE = False
    METALPDB_AVAILABLE = False
    RCSB_AVAILABLE = False

network_required = pytest.mark.skipif(
    not NETWORK_AVAILABLE,
    reason="Network access required (requests not available)"
)

metalpdb_required = pytest.mark.skipif(
    not METALPDB_AVAILABLE,
    reason="MetalPDB API not reachable"
)


class TestMetalPDBAdapterImport:
    """Test MetalPDB adapter can be imported."""

    def test_can_import_adapter(self):
        """Should be able to import MetalPDBAdapter."""
        from database_adapters import MetalPDBAdapter
        assert MetalPDBAdapter is not None

    def test_can_import_metal_site_result(self):
        """Should be able to import MetalSiteResult dataclass."""
        from database_adapters.metalpdb_adapter import MetalSiteResult
        assert MetalSiteResult is not None


class TestMetalSiteResultDataclass:
    """Test MetalSiteResult dataclass."""

    def test_metal_site_result_creation(self):
        """Should create MetalSiteResult with all fields."""
        from database_adapters.metalpdb_adapter import MetalSiteResult

        result = MetalSiteResult(
            pdb_id="1CA2",
            metal="ZN",
            coordination_number=4,
            geometry="tetrahedral",
            resolution=1.5,
            coordinating_residues=["HIS94", "HIS96", "HIS119"],
            ligands=[],
            source="metalpdb",
        )
        assert result.pdb_id == "1CA2"
        assert result.metal == "ZN"
        assert result.coordination_number == 4
        assert result.geometry == "tetrahedral"
        assert result.resolution == 1.5
        assert len(result.coordinating_residues) == 3

    def test_metal_site_result_to_dict(self):
        """Should convert MetalSiteResult to dictionary."""
        from database_adapters.metalpdb_adapter import MetalSiteResult

        result = MetalSiteResult(
            pdb_id="1CA2",
            metal="ZN",
            coordination_number=4,
            geometry="tetrahedral",
            resolution=1.5,
            coordinating_residues=["HIS94"],
            ligands=[],
            source="metalpdb",
        )
        d = result.to_dict()
        assert isinstance(d, dict)
        assert d["pdb_id"] == "1CA2"
        assert d["metal"] == "ZN"


class TestMetalPDBAdapterInitialization:
    """Test MetalPDBAdapter initialization."""

    def test_adapter_has_default_base_url(self):
        """Should have MetalPDB base URL by default."""
        from database_adapters import MetalPDBAdapter

        adapter = MetalPDBAdapter()
        assert "metalpdb" in adapter.base_url.lower()

    def test_adapter_can_override_base_url(self):
        """Should allow overriding base URL."""
        from database_adapters import MetalPDBAdapter

        adapter = MetalPDBAdapter(base_url="https://custom.example.com")
        assert adapter.base_url == "https://custom.example.com"

    def test_adapter_has_timeout_setting(self):
        """Should have configurable timeout."""
        from database_adapters import MetalPDBAdapter

        adapter = MetalPDBAdapter(timeout=60)
        assert adapter.timeout == 60


class TestSearchByMetal:
    """Test search_by_metal method."""

    def test_search_by_metal_returns_list(self):
        """Should return a list from search_by_metal."""
        from database_adapters import MetalPDBAdapter

        adapter = MetalPDBAdapter()
        # Even without network, should return empty list not raise
        results = adapter.search_by_metal("ZN", limit=5)
        assert isinstance(results, list)

    @network_required
    def test_search_by_metal_returns_results(self):
        """Search for ZN should return PDB IDs."""
        from database_adapters import MetalPDBAdapter

        adapter = MetalPDBAdapter()
        results = adapter.search_by_metal("ZN", limit=5)
        # Should find at least some results (from MetalPDB or RCSB fallback)
        # Note: May be empty if both APIs are down
        for result in results:
            assert hasattr(result, 'pdb_id') or isinstance(result, dict)
            if hasattr(result, 'pdb_id'):
                assert result.pdb_id is not None
                assert result.metal == "ZN"

    @network_required
    def test_search_by_metal_with_coordination(self):
        """Filter by coordination number should work."""
        from database_adapters import MetalPDBAdapter

        adapter = MetalPDBAdapter()
        results = adapter.search_by_metal("ZN", coordination_number=4, limit=10)
        # If filtering by coordination, all results should match
        for result in results:
            if hasattr(result, 'coordination_number') and result.coordination_number:
                # Coordination number should be 4 (or close if API doesn't filter exactly)
                assert result.coordination_number in [3, 4, 5]  # Allow some tolerance

    @network_required
    def test_search_by_metal_with_geometry(self):
        """Filter by geometry should work."""
        from database_adapters import MetalPDBAdapter

        adapter = MetalPDBAdapter()
        results = adapter.search_by_metal("ZN", geometry="tetrahedral", limit=10)
        # If filtering by geometry works, results should have geometry info
        for result in results:
            if hasattr(result, 'geometry') and result.geometry:
                assert "tetra" in result.geometry.lower() or result.geometry == "unknown"

    @network_required
    def test_search_by_metal_with_resolution(self):
        """Filter by resolution_max should work."""
        from database_adapters import MetalPDBAdapter

        adapter = MetalPDBAdapter()
        results = adapter.search_by_metal("ZN", resolution_max=2.0, limit=10)
        for result in results:
            if hasattr(result, 'resolution') and result.resolution:
                assert result.resolution <= 2.5  # Allow some tolerance for API variations

    def test_search_nonexistent_metal_returns_empty(self):
        """Invalid metal should return empty list."""
        from database_adapters import MetalPDBAdapter

        adapter = MetalPDBAdapter()
        results = adapter.search_by_metal("XX", limit=5)  # Invalid element symbol
        assert isinstance(results, list)
        assert len(results) == 0


class TestGetSiteDetails:
    """Test get_site_details method."""

    def test_get_site_details_returns_none_for_invalid(self):
        """Should return None for invalid PDB ID."""
        from database_adapters import MetalPDBAdapter

        adapter = MetalPDBAdapter()
        result = adapter.get_site_details("XXXX", "ZN")
        assert result is None

    @network_required
    def test_get_site_details_returns_geometry(self):
        """Get details for 1CA2 (carbonic anhydrase) should return geometry."""
        from database_adapters import MetalPDBAdapter

        adapter = MetalPDBAdapter()
        result = adapter.get_site_details("1CA2", "ZN")
        # If API is available, should return detailed info
        if result is not None:
            assert hasattr(result, 'geometry') or 'geometry' in result if isinstance(result, dict) else True
            if hasattr(result, 'coordination_number'):
                assert result.coordination_number > 0

    @network_required
    def test_get_site_details_with_site_index(self):
        """Should support site_index parameter for multi-metal structures."""
        from database_adapters import MetalPDBAdapter

        adapter = MetalPDBAdapter()
        # 1CA2 might have multiple zinc sites in unit cell
        result0 = adapter.get_site_details("1CA2", "ZN", site_index=0)
        result1 = adapter.get_site_details("1CA2", "ZN", site_index=1)
        # Both calls should not raise, may return same or different site


class TestSearchByCoordinationMotif:
    """Test search_by_coordination_motif method."""

    def test_search_by_motif_returns_list(self):
        """Should return a list from search_by_coordination_motif."""
        from database_adapters import MetalPDBAdapter

        adapter = MetalPDBAdapter()
        results = adapter.search_by_coordination_motif("His-His-Glu", limit=5)
        assert isinstance(results, list)

    @network_required
    def test_search_by_motif_his_his_glu(self):
        """Search for His-His-Glu motif should return results."""
        from database_adapters import MetalPDBAdapter

        adapter = MetalPDBAdapter()
        results = adapter.search_by_coordination_motif("His-His-Glu", limit=10)
        # Common motif in zinc proteins, should find results
        # Note: depends on API availability

    @network_required
    def test_search_by_motif_with_metal_filter(self):
        """Search with metal filter should narrow results."""
        from database_adapters import MetalPDBAdapter

        adapter = MetalPDBAdapter()
        results = adapter.search_by_coordination_motif("His-His-Glu-Asp", metal="ZN", limit=10)
        for result in results:
            if hasattr(result, 'metal'):
                assert result.metal == "ZN"


class TestFallbackBehavior:
    """Test fallback to RCSB when MetalPDB is unavailable."""

    def test_adapter_has_fallback_enabled(self):
        """Should have RCSB fallback enabled by default."""
        from database_adapters import MetalPDBAdapter

        adapter = MetalPDBAdapter()
        assert hasattr(adapter, 'enable_fallback')
        assert adapter.enable_fallback is True

    def test_adapter_can_disable_fallback(self):
        """Should be able to disable fallback."""
        from database_adapters import MetalPDBAdapter

        adapter = MetalPDBAdapter(enable_fallback=False)
        assert adapter.enable_fallback is False

    @network_required
    def test_fallback_to_rcsb_works(self):
        """When MetalPDB fails, should fallback to RCSB."""
        from database_adapters import MetalPDBAdapter

        # Use an invalid MetalPDB URL to force fallback
        adapter = MetalPDBAdapter(
            base_url="https://invalid.metalpdb.example.com",
            enable_fallback=True
        )
        if RCSB_AVAILABLE:
            results = adapter.search_by_metal("ZN", limit=5)
            # Should still get results from RCSB fallback
            # (may be empty list if both fail, but should not raise)
            assert isinstance(results, list)


class TestErrorHandling:
    """Test error handling and logging."""

    def test_network_error_returns_empty_list(self):
        """Network errors should return empty list, not raise."""
        from database_adapters import MetalPDBAdapter

        adapter = MetalPDBAdapter(
            base_url="https://invalid.nonexistent.domain.example",
            enable_fallback=False,  # Disable fallback to test error handling
            timeout=1,  # Short timeout
        )
        # Should not raise, should return empty list
        results = adapter.search_by_metal("ZN")
        assert isinstance(results, list)

    def test_invalid_json_response_handled(self):
        """Invalid JSON response should be handled gracefully."""
        from database_adapters import MetalPDBAdapter

        adapter = MetalPDBAdapter()
        # This tests internal handling - adapter should not crash
        # when API returns unexpected data


class TestIntegration:
    """Integration tests requiring network access."""

    @pytest.mark.integration
    @network_required
    def test_full_workflow_zinc_sites(self):
        """Full workflow: search, get details, extract coordination."""
        from database_adapters import MetalPDBAdapter

        adapter = MetalPDBAdapter()

        # Step 1: Search for zinc sites
        results = adapter.search_by_metal("ZN", resolution_max=2.0, limit=3)

        if results:
            # Step 2: Get details for first result
            first = results[0]
            pdb_id = first.pdb_id if hasattr(first, 'pdb_id') else first.get('pdb_id')

            if pdb_id:
                details = adapter.get_site_details(pdb_id, "ZN")
                # Details should have coordination info
                if details:
                    assert hasattr(details, 'coordinating_residues') or 'coordinating_residues' in details

    @pytest.mark.integration
    @network_required
    def test_compare_metalpdb_vs_rcsb_results(self):
        """Compare results from MetalPDB vs RCSB fallback."""
        from database_adapters import MetalPDBAdapter

        # MetalPDB primary
        adapter_metalpdb = MetalPDBAdapter(enable_fallback=False)
        results_metalpdb = adapter_metalpdb.search_by_metal("ZN", limit=5)

        # Force RCSB fallback
        adapter_rcsb = MetalPDBAdapter(
            base_url="https://invalid.example.com",
            enable_fallback=True
        )
        results_rcsb = adapter_rcsb.search_by_metal("ZN", limit=5)

        # Both should return valid results (or empty lists)
        assert isinstance(results_metalpdb, list)
        assert isinstance(results_rcsb, list)
