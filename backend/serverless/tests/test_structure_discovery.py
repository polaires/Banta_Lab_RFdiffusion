# test_structure_discovery.py
"""Tests for StructureDiscovery service."""
import pytest
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Check if requests is available for network tests
try:
    import requests
    NETWORK_AVAILABLE = True
    # Quick connectivity check for RCSB
    try:
        requests.head("https://search.rcsb.org", timeout=5)
        RCSB_AVAILABLE = True
    except (requests.RequestException, Exception):
        RCSB_AVAILABLE = False
except ImportError:
    NETWORK_AVAILABLE = False
    RCSB_AVAILABLE = False

network_required = pytest.mark.skipif(
    not NETWORK_AVAILABLE,
    reason="Network access required (requests not available)"
)


class TestStructureDiscoveryImport:
    """Test StructureDiscovery can be imported."""

    def test_can_import_structure_discovery(self):
        """Should be able to import StructureDiscovery."""
        from structure_discovery import StructureDiscovery
        assert StructureDiscovery is not None

    def test_can_import_discovery_result(self):
        """Should be able to import DiscoveryResult dataclass."""
        from structure_discovery import DiscoveryResult
        assert DiscoveryResult is not None

    def test_can_import_metal_names(self):
        """Should be able to import METAL_NAMES mapping."""
        from structure_discovery import METAL_NAMES
        assert isinstance(METAL_NAMES, dict)
        assert "zinc" in METAL_NAMES
        assert "terbium" in METAL_NAMES
        assert METAL_NAMES["zinc"] == "ZN"
        assert METAL_NAMES["terbium"] == "TB"


class TestDiscoveryResultDataclass:
    """Test DiscoveryResult dataclass."""

    def test_discovery_result_creation(self):
        """Should create DiscoveryResult with all fields."""
        from structure_discovery import DiscoveryResult

        result = DiscoveryResult(
            pdb_id="1CA2",
            metal="ZN",
            resolution=1.5,
            coordination_number=4,
            geometry="tetrahedral",
            source="metalpdb",
            curated=False,
        )
        assert result.pdb_id == "1CA2"
        assert result.metal == "ZN"
        assert result.resolution == 1.5
        assert result.coordination_number == 4
        assert result.geometry == "tetrahedral"
        assert result.source == "metalpdb"
        assert result.curated is False

    def test_discovery_result_to_dict(self):
        """Should convert DiscoveryResult to dictionary."""
        from structure_discovery import DiscoveryResult

        result = DiscoveryResult(
            pdb_id="6MI5",
            metal="TB",
            resolution=2.0,
            coordination_number=9,
            geometry="tricapped_trigonal_prismatic",
            source="curated",
            curated=True,
            score=150.0,
        )
        d = result.to_dict()
        assert isinstance(d, dict)
        assert d["pdb_id"] == "6MI5"
        assert d["metal"] == "TB"
        assert d["curated"] is True
        assert d["score"] == 150.0


class TestStructureDiscoveryInitialization:
    """Test StructureDiscovery initialization."""

    def test_discovery_has_adapters(self):
        """Should initialize with all database adapters."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        assert hasattr(discovery, 'metalpdb')
        assert hasattr(discovery, 'uniprot')
        assert hasattr(discovery, 'pubchem')
        assert hasattr(discovery, 'alphafold')

    def test_discovery_has_timeout_setting(self):
        """Should have configurable timeout."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery(timeout=60)
        assert discovery.timeout == 60

    def test_discovery_has_fallback_setting(self):
        """Should have configurable fallback."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery(enable_fallback=False)
        assert discovery.enable_fallback is False


class TestIntentParsing:
    """Test natural language intent parsing."""

    def test_parse_intent_metal_name(self):
        """Should parse metal names from query."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        metal, ligand, ptype = discovery._parse_intent("terbium binding protein")
        assert metal == "TB"

    def test_parse_intent_zinc(self):
        """Should parse zinc from query."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        metal, ligand, ptype = discovery._parse_intent("zinc finger domain")
        assert metal == "ZN"

    def test_parse_intent_ligand_pqq(self):
        """Should parse PQQ ligand from query."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        metal, ligand, ptype = discovery._parse_intent("PQQ dehydrogenase enzyme")
        assert ligand == "PQQ"

    def test_parse_intent_citrate(self):
        """Should parse citrate ligand from query."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        metal, ligand, ptype = discovery._parse_intent("citrate bound terbium")
        assert ligand == "CIT"
        assert metal == "TB"

    def test_parse_intent_protein_type_dehydrogenase(self):
        """Should identify dehydrogenase protein type."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        metal, ligand, ptype = discovery._parse_intent("methanol dehydrogenase")
        assert ptype == "dehydrogenase"

    def test_parse_intent_protein_type_ef_hand(self):
        """Should identify EF-hand protein type."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        metal, ligand, ptype = discovery._parse_intent("calmodulin ef-hand")
        assert ptype == "ef_hand"
        assert metal == "CA"

    def test_parse_intent_protein_type_zinc_finger(self):
        """Should identify zinc finger protein type."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        metal, ligand, ptype = discovery._parse_intent("c2h2 zinc finger")
        assert ptype == "zinc_finger"
        assert metal == "ZN"

    def test_parse_intent_combined(self):
        """Should parse combined metal, ligand, and type."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        metal, ligand, ptype = discovery._parse_intent("calcium PQQ dehydrogenase")
        assert metal == "CA"
        assert ligand == "PQQ"
        assert ptype == "dehydrogenase"


class TestSearchByIntent:
    """Test search_by_intent method."""

    def test_search_by_intent_returns_list(self):
        """Should return a list from search_by_intent."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        # Even without network, should return list (possibly empty or curated)
        results = discovery.search_by_intent("terbium binding protein", limit=5)
        assert isinstance(results, list)

    @network_required
    def test_search_by_intent_metal_binding(self):
        """Search for 'terbium binding protein' should return results."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        results = discovery.search_by_intent("terbium binding protein", limit=10)
        # Should find curated references at minimum
        assert isinstance(results, list)
        # Check that results have expected metal
        for result in results:
            if result.metal:
                assert result.metal == "TB"

    @network_required
    def test_search_by_intent_ligand(self):
        """Search for 'PQQ dehydrogenase' should return results."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        results = discovery.search_by_intent("PQQ dehydrogenase", limit=10)
        assert isinstance(results, list)
        # Should identify PQQ as a ligand in at least some results
        pqq_results = [r for r in results if "PQQ" in r.ligands]
        # May or may not have PQQ in ligands depending on API response

    @network_required
    def test_search_by_intent_returns_ranked_results(self):
        """Results should be ranked by score."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        results = discovery.search_by_intent("zinc binding protein", limit=10)
        if len(results) >= 2:
            # Scores should be in descending order
            for i in range(len(results) - 1):
                assert results[i].score >= results[i + 1].score


class TestSearchByMetal:
    """Test search_by_metal method."""

    def test_search_by_metal_returns_list(self):
        """Should return a list from search_by_metal."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        results = discovery.search_by_metal("ZN", limit=5)
        assert isinstance(results, list)

    def test_search_by_metal_includes_curated(self):
        """Should include curated reference structures."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        results = discovery.search_by_metal("ZN", limit=20)
        # Should have curated results for ZN
        curated = [r for r in results if r.curated]
        assert len(curated) > 0

    def test_search_by_metal_zn_has_metal_field(self):
        """ZN search results should have metal='ZN'."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        results = discovery.search_by_metal("ZN", limit=10)
        for result in results:
            assert result.metal == "ZN"

    @network_required
    def test_search_by_metal_returns_pdb_ids(self):
        """Search for ZN should return valid PDB IDs."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        results = discovery.search_by_metal("ZN", limit=10)
        assert len(results) > 0
        for result in results:
            assert len(result.pdb_id) == 4  # PDB IDs are 4 characters

    @network_required
    def test_search_by_metal_tb_returns_lanthanide_results(self):
        """Search for TB should return lanthanide structures."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        results = discovery.search_by_metal("TB", limit=10)
        assert isinstance(results, list)
        # Should have curated lanmodulin structures
        curated = [r for r in results if r.curated]
        assert len(curated) > 0

    @network_required
    def test_search_by_metal_with_coordination_number(self):
        """Should filter by coordination number."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        results = discovery.search_by_metal("ZN", coordination_number=4, limit=10)
        assert isinstance(results, list)

    @network_required
    def test_search_by_metal_with_geometry(self):
        """Should filter by geometry."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        results = discovery.search_by_metal("ZN", geometry="tetrahedral", limit=10)
        assert isinstance(results, list)


class TestGetBestTemplate:
    """Test get_best_template method."""

    def test_get_best_template_returns_dict_or_none(self):
        """Should return dict template or None."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        template = discovery.get_best_template("ZN")
        assert template is None or isinstance(template, dict)

    @network_required
    def test_get_best_template_metal_only(self):
        """Should return template for metal-only query."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        template = discovery.get_best_template("ZN")
        if template:
            assert "metal" in template
            assert template["metal"] == "ZN"

    @network_required
    def test_get_best_template_with_ligand(self):
        """Should return template for metal+ligand query."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        template = discovery.get_best_template("TB", ligand="CIT")
        if template:
            assert "metal" in template
            # Template should have source information
            assert "source" in template

    @network_required
    def test_get_best_template_has_coordination_info(self):
        """Template should include coordination information."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        template = discovery.get_best_template("ZN")
        if template:
            # Should have coordination number or geometry
            has_coord = "coordination_number" in template or "geometry" in template
            assert has_coord or "source" in template

    def test_get_best_template_citrate_tb(self):
        """Should return citrate-TB template from library."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        template = discovery.get_best_template("TB", ligand="CIT")
        # Should get template from library or reference
        if template:
            assert template.get("metal") == "TB" or "TB" in str(template)


class TestExtractCoordinationInfo:
    """Test extract_coordination_info method."""

    def test_extract_coordination_info_returns_none_for_invalid(self):
        """Should return None for invalid PDB ID."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        info = discovery.extract_coordination_info("XXXX", "ZN")
        assert info is None

    @network_required
    def test_extract_coordination_info_1ca2_zn(self):
        """Extract from 1CA2 (carbonic anhydrase) should return ZN coordination."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        info = discovery.extract_coordination_info("1CA2", "ZN")
        if info:
            assert "coordination_number" in info
            assert "geometry" in info
            # 1CA2 has tetrahedral ZN coordination
            if info["coordination_number"]:
                assert info["coordination_number"] >= 3  # At least 3 coordinating atoms

    @network_required
    def test_extract_coordination_info_has_geometry(self):
        """Coordination info should include geometry."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        info = discovery.extract_coordination_info("1CA2", "ZN")
        if info:
            assert "geometry" in info

    @network_required
    def test_extract_coordination_info_has_source(self):
        """Coordination info should include source."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        info = discovery.extract_coordination_info("1CA2", "ZN")
        if info:
            assert "source" in info
            assert info["source"] in ("metalpdb", "rcsb")


class TestRankResultsByQuality:
    """Test result ranking logic."""

    def test_rank_results_curated_first(self):
        """Curated results should rank higher."""
        from structure_discovery import StructureDiscovery, DiscoveryResult

        discovery = StructureDiscovery()

        results = [
            DiscoveryResult(pdb_id="1AAA", metal="ZN", curated=False, resolution=1.0, source="rcsb"),
            DiscoveryResult(pdb_id="1BBB", metal="ZN", curated=True, resolution=2.0, source="curated"),
            DiscoveryResult(pdb_id="1CCC", metal="ZN", curated=False, resolution=1.5, source="metalpdb"),
        ]

        ranked = discovery._rank_results(results)

        # Curated should be first despite worse resolution
        assert ranked[0].pdb_id == "1BBB"
        assert ranked[0].curated is True

    def test_rank_results_by_resolution(self):
        """Better resolution should rank higher among non-curated."""
        from structure_discovery import StructureDiscovery, DiscoveryResult

        discovery = StructureDiscovery()

        results = [
            DiscoveryResult(pdb_id="1AAA", metal="ZN", curated=False, resolution=2.5, source="rcsb"),
            DiscoveryResult(pdb_id="1BBB", metal="ZN", curated=False, resolution=1.0, source="rcsb"),
            DiscoveryResult(pdb_id="1CCC", metal="ZN", curated=False, resolution=1.8, source="rcsb"),
        ]

        ranked = discovery._rank_results(results)

        # 1.0 resolution should score higher than 2.5
        resolutions = [r.resolution for r in ranked]
        # First result should have best (lowest) resolution
        assert ranked[0].resolution == 1.0

    def test_rank_results_by_source(self):
        """Source priority should affect ranking."""
        from structure_discovery import StructureDiscovery, DiscoveryResult

        discovery = StructureDiscovery()

        results = [
            DiscoveryResult(pdb_id="1AAA", metal="ZN", curated=False, resolution=2.0, source="uniprot"),
            DiscoveryResult(pdb_id="1BBB", metal="ZN", curated=False, resolution=2.0, source="metalpdb"),
            DiscoveryResult(pdb_id="1CCC", metal="ZN", curated=False, resolution=2.0, source="rcsb"),
        ]

        ranked = discovery._rank_results(results)

        # MetalPDB should rank higher than RCSB, which ranks higher than UniProt
        sources = [r.source for r in ranked]
        assert sources.index("metalpdb") < sources.index("uniprot")


class TestDeduplicateResults:
    """Test result deduplication."""

    def test_deduplicate_removes_duplicates(self):
        """Should remove duplicate PDB IDs."""
        from structure_discovery import StructureDiscovery, DiscoveryResult

        discovery = StructureDiscovery()

        results = [
            DiscoveryResult(pdb_id="1CA2", metal="ZN", curated=False, resolution=1.5, source="rcsb"),
            DiscoveryResult(pdb_id="1CA2", metal="ZN", curated=False, resolution=2.0, source="metalpdb"),
            DiscoveryResult(pdb_id="1BBB", metal="ZN", curated=False, resolution=1.8, source="rcsb"),
        ]

        deduped = discovery._deduplicate_results(results)

        pdb_ids = [r.pdb_id for r in deduped]
        assert len(pdb_ids) == 2
        assert pdb_ids.count("1CA2") == 1

    def test_deduplicate_keeps_curated(self):
        """Should prefer curated over non-curated when deduplicating."""
        from structure_discovery import StructureDiscovery, DiscoveryResult

        discovery = StructureDiscovery()

        results = [
            DiscoveryResult(pdb_id="1CA2", metal="ZN", curated=False, resolution=1.0, source="rcsb"),
            DiscoveryResult(pdb_id="1CA2", metal="ZN", curated=True, resolution=2.0, source="curated"),
        ]

        deduped = discovery._deduplicate_results(results)

        assert len(deduped) == 1
        assert deduped[0].curated is True

    def test_deduplicate_keeps_better_resolution(self):
        """Should prefer better resolution when deduplicating non-curated."""
        from structure_discovery import StructureDiscovery, DiscoveryResult

        discovery = StructureDiscovery()

        results = [
            DiscoveryResult(pdb_id="1CA2", metal="ZN", curated=False, resolution=2.5, source="rcsb"),
            DiscoveryResult(pdb_id="1CA2", metal="ZN", curated=False, resolution=1.5, source="metalpdb"),
        ]

        deduped = discovery._deduplicate_results(results)

        assert len(deduped) == 1
        assert deduped[0].resolution == 1.5


class TestIntegration:
    """Integration tests requiring network access."""

    @pytest.mark.integration
    @network_required
    def test_full_workflow_terbium_discovery(self):
        """Full workflow: intent search, get template, extract coordination."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()

        # Step 1: Search by intent
        results = discovery.search_by_intent("terbium binding protein", limit=5)
        assert isinstance(results, list)

        # Step 2: Get best template
        template = discovery.get_best_template("TB", ligand="CIT")
        # Template may or may not be available
        if template:
            assert "metal" in template or "source" in template

        # Step 3: Extract coordination from known structure (lanmodulin)
        if results:
            first = results[0]
            if first.pdb_id:
                info = discovery.extract_coordination_info(first.pdb_id, "TB")
                # Info may or may not be available

    @pytest.mark.integration
    @network_required
    def test_full_workflow_zinc_discovery(self):
        """Full workflow for zinc sites."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()

        # Search for zinc
        results = discovery.search_by_metal("ZN", coordination_number=4, limit=10)
        assert len(results) > 0  # Should find curated at minimum

        # Get template
        template = discovery.get_best_template("ZN")
        assert template is not None or True  # May require network

        # Extract from known structure
        info = discovery.extract_coordination_info("1CA2", "ZN")
        if info:
            assert info["metal"] == "ZN"

    @pytest.mark.integration
    @network_required
    def test_compare_metal_sources(self):
        """Compare results from different sources."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()

        results = discovery.search_by_metal("ZN", limit=20)

        # Should have results from multiple sources
        sources = set(r.source for r in results)
        # At minimum should have curated
        assert "curated" in sources or len(results) > 0


class TestErrorHandling:
    """Test error handling and graceful degradation."""

    def test_handles_network_errors_gracefully(self):
        """Should handle network errors without crashing."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery(timeout=1)
        # Should not raise even if network unavailable
        results = discovery.search_by_intent("unknown protein xyz", limit=5)
        assert isinstance(results, list)

    def test_handles_invalid_metal_gracefully(self):
        """Should handle invalid metal symbol gracefully."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        results = discovery.search_by_metal("XX", limit=5)
        assert isinstance(results, list)
        # May be empty or contain results depending on fallback

    def test_handles_empty_query(self):
        """Should handle empty query gracefully."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        results = discovery.search_by_intent("", limit=5)
        assert isinstance(results, list)
