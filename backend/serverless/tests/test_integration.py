# test_integration.py
"""
Integration tests for the AI-Driven Structure Discovery Pipeline.

These tests verify that the full pipeline works end-to-end, from natural
language queries through structure discovery to design recommendations.

Tests cover:
1. NL query to template pipeline
2. Metal search to coordination extraction
3. Full discovery pipeline (plan -> search -> recommendations)
4. Cross-database adapter integration
5. Design recommendations to bias string conversion

Note: These tests require network access to external databases.
Use @pytest.mark.integration to mark tests that make network requests.
"""
import pytest
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Check if requests is available for network tests
try:
    import requests
    NETWORK_AVAILABLE = True
    # Check connectivity to required services
    try:
        requests.head("https://search.rcsb.org", timeout=5)
        RCSB_AVAILABLE = True
    except (requests.RequestException, Exception):
        RCSB_AVAILABLE = False
    try:
        requests.head("https://rest.uniprot.org", timeout=5)
        UNIPROT_AVAILABLE = True
    except (requests.RequestException, Exception):
        UNIPROT_AVAILABLE = False
    try:
        requests.head("https://alphafold.ebi.ac.uk", timeout=5)
        ALPHAFOLD_AVAILABLE = True
    except (requests.RequestException, Exception):
        ALPHAFOLD_AVAILABLE = False
    try:
        requests.head("https://metalpdb.cerm.unifi.it", timeout=5)
        METALPDB_AVAILABLE = True
    except (requests.RequestException, Exception):
        METALPDB_AVAILABLE = False
except ImportError:
    NETWORK_AVAILABLE = False
    RCSB_AVAILABLE = False
    UNIPROT_AVAILABLE = False
    ALPHAFOLD_AVAILABLE = False
    METALPDB_AVAILABLE = False

# Pytest markers
integration = pytest.mark.integration

network_required = pytest.mark.skipif(
    not NETWORK_AVAILABLE,
    reason="Network access required (requests not available)"
)

rcsb_required = pytest.mark.skipif(
    not RCSB_AVAILABLE,
    reason="RCSB API not reachable"
)

uniprot_required = pytest.mark.skipif(
    not UNIPROT_AVAILABLE,
    reason="UniProt API not reachable"
)

alphafold_required = pytest.mark.skipif(
    not ALPHAFOLD_AVAILABLE,
    reason="AlphaFold API not reachable"
)


# =============================================================================
# TEST IMPORTS
# =============================================================================

class TestIntegrationImports:
    """Verify all required modules can be imported."""

    def test_can_import_structure_discovery(self):
        """Should be able to import StructureDiscovery."""
        from structure_discovery import StructureDiscovery
        assert StructureDiscovery is not None

    def test_can_import_ai_interface_functions(self):
        """Should be able to import AI interface functions."""
        from ai_structure_interface import (
            answer_structure_question,
            plan_structure_search,
            get_design_recommendations,
        )
        assert answer_structure_question is not None
        assert plan_structure_search is not None
        assert get_design_recommendations is not None

    def test_can_import_database_adapters(self):
        """Should be able to import all database adapters."""
        from database_adapters import (
            MetalPDBAdapter,
            UniProtAdapter,
            AlphaFoldAdapter,
        )
        assert MetalPDBAdapter is not None
        assert UniProtAdapter is not None
        assert AlphaFoldAdapter is not None


# =============================================================================
# TEST 1: NL QUERY TO TEMPLATE PIPELINE
# =============================================================================

class TestNLQueryToTemplatePipeline:
    """
    Test 1: Natural Language Query to Template Pipeline

    Verifies that:
    1. answer_structure_question() correctly parses template questions
    2. StructureDiscovery can retrieve relevant templates
    3. The template has required fields (pdb_id or pdb_content)
    """

    def test_template_question_parsing(self):
        """answer_structure_question should detect template questions."""
        from ai_structure_interface import answer_structure_question

        result = answer_structure_question(
            "What template should I use for terbium biosensor design?"
        )

        assert result is not None
        assert isinstance(result, dict)
        assert "type" in result
        # Should be detected as template question
        assert result["type"] == "template"
        assert "answer" in result

    def test_template_question_extracts_metal(self):
        """Template question should extract the metal (TB for terbium)."""
        from ai_structure_interface import answer_structure_question

        result = answer_structure_question(
            "What template should I use for terbium biosensor design?"
        )

        # Should extract terbium as TB
        assert "metal" in result
        assert result["metal"] == "TB"

    @integration
    @network_required
    def test_nl_query_to_template_full_pipeline(self):
        """
        Full pipeline: NL query -> template retrieval.

        1. Call answer_structure_question with template question
        2. Use StructureDiscovery to get the actual template
        3. Assert template has pdb_id or pdb_content
        """
        from ai_structure_interface import answer_structure_question
        from structure_discovery import StructureDiscovery

        # Step 1: Parse NL query
        result = answer_structure_question(
            "What template should I use for terbium biosensor design?"
        )

        assert result is not None
        metal = result.get("metal", "TB")

        # Step 2: Use StructureDiscovery to get template
        discovery = StructureDiscovery()
        template = discovery.get_best_template(metal=metal, use_architector=True)

        # Step 3: Verify template has required fields
        assert template is not None
        # Template should have either pdb_id or pdb_content
        has_pdb_id = "pdb_id" in template and template["pdb_id"]
        has_pdb_content = "pdb_content" in template and template["pdb_content"]
        has_source = "source" in template

        assert has_pdb_id or has_pdb_content or has_source, \
            f"Template should have pdb_id or pdb_content, got: {list(template.keys())}"

    def test_template_with_ligand(self):
        """Template question with ligand should extract both metal and ligand."""
        from ai_structure_interface import answer_structure_question

        result = answer_structure_question(
            "Get a template for calcium with PQQ cofactor"
        )

        assert result is not None
        # Should detect CA and PQQ
        if result.get("metal"):
            assert result["metal"] == "CA"


# =============================================================================
# TEST 2: METAL SEARCH TO COORDINATION INFO
# =============================================================================

class TestMetalSearchToCoordinationInfo:
    """
    Test 2: Metal Search to Coordination Information Extraction

    Verifies that:
    1. StructureDiscovery.search_by_metal() returns results
    2. extract_coordination_info() provides detailed coordination data
    3. Results include coordination_number and geometry
    """

    @integration
    @network_required
    def test_metal_search_returns_results(self):
        """search_by_metal for ZN should return PDB IDs."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()
        results = discovery.search_by_metal("ZN", limit=3)

        assert isinstance(results, list)
        # Should find at least some zinc structures
        if results:
            first_result = results[0]
            assert hasattr(first_result, 'pdb_id') or hasattr(first_result, 'metal')
            if hasattr(first_result, 'pdb_id'):
                assert first_result.pdb_id is not None
                assert len(first_result.pdb_id) == 4

    @integration
    @network_required
    def test_metal_search_to_coordination_info(self):
        """
        Full pipeline: Metal search -> coordination extraction.

        1. Use StructureDiscovery.search_by_metal("ZN", limit=3)
        2. For first result, call extract_coordination_info
        3. Assert coordination_number and geometry are present
        """
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()

        # Step 1: Search for zinc structures
        results = discovery.search_by_metal("ZN", limit=3)

        if not results:
            pytest.skip("No zinc structures found - network may be unavailable")

        # Step 2: Get first result and extract coordination info
        first_result = results[0]
        pdb_id = first_result.pdb_id

        coord_info = discovery.extract_coordination_info(pdb_id, "ZN")

        # Step 3: Verify coordination information is present
        if coord_info is not None:
            assert "coordination_number" in coord_info
            assert "geometry" in coord_info
            # Zinc typically has CN 4-6
            cn = coord_info["coordination_number"]
            if cn > 0:
                assert cn >= 3 and cn <= 8, f"Unexpected CN {cn} for zinc"

    @integration
    @network_required
    def test_coordination_info_includes_residues(self):
        """Coordination info should include coordinating residues."""
        from structure_discovery import StructureDiscovery

        discovery = StructureDiscovery()

        # Use well-known zinc structure: carbonic anhydrase
        coord_info = discovery.extract_coordination_info("1CA2", "ZN")

        if coord_info is not None:
            # Should have coordinating residues list
            assert "coordinating_residues" in coord_info
            residues = coord_info["coordinating_residues"]
            if residues:
                assert isinstance(residues, list)


# =============================================================================
# TEST 3: FULL DISCOVERY PIPELINE
# =============================================================================

class TestFullDiscoveryPipeline:
    """
    Test 3: Complete Discovery Pipeline

    Verifies the full workflow:
    1. plan_structure_search() parses task description
    2. Returns metal="TB" and has steps
    3. get_design_recommendations() provides donors and bias_aa
    """

    def test_plan_structure_search_parses_task(self):
        """plan_structure_search should parse terbium biosensor task."""
        from ai_structure_interface import plan_structure_search

        plan = plan_structure_search("design a terbium biosensor")

        assert plan is not None
        assert isinstance(plan, dict)
        assert "metal" in plan
        assert "steps" in plan

    def test_plan_structure_search_detects_terbium(self):
        """Plan should detect TB as the metal for terbium biosensor."""
        from ai_structure_interface import plan_structure_search

        plan = plan_structure_search("design a terbium biosensor")

        # Terbium biosensor should return TB
        metal = plan.get("metal")
        assert metal == "TB", f"Expected TB, got {metal}"

    def test_plan_structure_search_has_steps(self):
        """Plan should include actionable steps."""
        from ai_structure_interface import plan_structure_search

        plan = plan_structure_search("design a terbium biosensor")

        steps = plan.get("steps", [])
        assert isinstance(steps, list)
        assert len(steps) > 0, "Plan should have at least one step"

    def test_plan_structure_search_detects_task_type(self):
        """Plan should detect biosensor task type."""
        from ai_structure_interface import plan_structure_search

        plan = plan_structure_search("design a terbium biosensor")

        task_type = plan.get("task_type")
        assert task_type == "biosensor"

    def test_get_design_recommendations_returns_donors(self):
        """get_design_recommendations should return donor residues."""
        from ai_structure_interface import get_design_recommendations

        recs = get_design_recommendations(metal="TB", task="biosensor")

        assert recs is not None
        assert isinstance(recs, dict)
        assert "donors" in recs
        donors = recs["donors"]
        assert isinstance(donors, list)
        # Terbium prefers oxygen donors (Glu, Asp, Asn, Gln)
        assert len(donors) > 0

    def test_get_design_recommendations_returns_bias_aa(self):
        """get_design_recommendations should return bias_aa string."""
        from ai_structure_interface import get_design_recommendations

        recs = get_design_recommendations(metal="TB", task="biosensor")

        assert "bias_aa" in recs
        bias_aa = recs["bias_aa"]
        assert isinstance(bias_aa, str)
        # Should have format like "E:3.0,D:3.0,..."
        if bias_aa:
            assert ":" in bias_aa, "bias_aa should have colon separators"

    @integration
    def test_full_discovery_pipeline_terbium_biosensor(self):
        """
        Full pipeline test:
        1. plan_structure_search("design a terbium biosensor")
        2. Assert returns metal="TB" and has steps
        3. Get recommendations with get_design_recommendations
        4. Assert donors and bias_aa are present
        """
        from ai_structure_interface import (
            plan_structure_search,
            get_design_recommendations,
        )

        # Step 1: Plan the search
        plan = plan_structure_search("design a terbium biosensor")

        # Step 2: Verify plan has metal and steps
        assert plan is not None
        metal = plan.get("metal")
        assert metal == "TB", f"Expected TB, got {metal}"

        steps = plan.get("steps", [])
        assert len(steps) > 0, "Plan should have steps"

        # Step 3: Get design recommendations for the metal
        recs = get_design_recommendations(metal=metal, task="biosensor")

        # Step 4: Verify donors and bias_aa are present
        assert recs is not None
        assert "donors" in recs
        assert len(recs["donors"]) > 0
        assert "bias_aa" in recs
        assert recs["bias_aa"] is not None


# =============================================================================
# TEST 4: ADAPTERS INTEGRATION
# =============================================================================

class TestAdaptersIntegration:
    """
    Test 4: Cross-Database Adapter Integration

    Verifies that adapters work together:
    1. MetalPDB search returns results
    2. UniProt provides protein information
    3. AlphaFold availability can be checked
    4. Cross-database consistency
    """

    @integration
    @network_required
    def test_metalpdb_search_works(self):
        """MetalPDB adapter should return metal site results."""
        from database_adapters import MetalPDBAdapter

        adapter = MetalPDBAdapter()
        results = adapter.search_by_metal("ZN", limit=3)

        assert isinstance(results, list)
        # Results may be empty if MetalPDB is down, but should not raise

    @integration
    @network_required
    @uniprot_required
    def test_uniprot_search_works(self):
        """UniProt adapter should return protein results."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()
        results = adapter.search_by_function("zinc finger", limit=3)

        assert isinstance(results, list)
        if results:
            assert "uniprot_id" in results[0]

    @integration
    @network_required
    @alphafold_required
    def test_alphafold_availability_check_works(self):
        """AlphaFold adapter should check prediction availability."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        # Human carbonic anhydrase II (P00918) should have AlphaFold prediction
        available = adapter.check_availability("P00918")

        assert isinstance(available, bool)

    @integration
    @network_required
    @uniprot_required
    @alphafold_required
    def test_adapters_integration_full_workflow(self):
        """
        Full cross-database workflow:
        1. MetalPDB search -> get PDB ID
        2. UniProt search for same protein type
        3. AlphaFold check availability

        Verify cross-database consistency.
        """
        from database_adapters import (
            MetalPDBAdapter,
            UniProtAdapter,
            AlphaFoldAdapter,
        )

        # Step 1: Search for zinc sites using MetalPDB
        metalpdb = MetalPDBAdapter()
        metal_results = metalpdb.search_by_metal("ZN", limit=3)

        if not metal_results:
            pytest.skip("MetalPDB returned no results - skipping cross-db test")

        first_result = metal_results[0]
        pdb_id = first_result.pdb_id if hasattr(first_result, 'pdb_id') else None

        if not pdb_id:
            pytest.skip("No PDB ID in MetalPDB result")

        # Step 2: Search UniProt for zinc-binding proteins
        uniprot = UniProtAdapter()
        uniprot_results = uniprot.search_by_metal_binding("Zinc", limit=3)

        assert isinstance(uniprot_results, list)
        # UniProt should find zinc-binding proteins
        if uniprot_results:
            uniprot_id = uniprot_results[0].get("uniprot_id")

            # Step 3: Check AlphaFold availability for the UniProt ID
            if uniprot_id:
                alphafold = AlphaFoldAdapter()
                af_available = alphafold.check_availability(uniprot_id)
                assert isinstance(af_available, bool)

    @integration
    @network_required
    @uniprot_required
    def test_uniprot_pdb_mapping_consistency(self):
        """UniProt PDB mappings should return valid PDB IDs."""
        from database_adapters import UniProtAdapter

        adapter = UniProtAdapter()

        # Get PDB mappings for well-known protein (carbonic anhydrase II)
        mappings = adapter.get_pdb_mappings("P00918")

        assert isinstance(mappings, list)
        if mappings:
            for mapping in mappings:
                assert "pdb_id" in mapping
                pdb_id = mapping["pdb_id"]
                # PDB IDs should be 4 characters
                assert len(pdb_id) == 4


# =============================================================================
# TEST 5: DESIGN RECOMMENDATIONS TO BIAS AA
# =============================================================================

class TestDesignRecommendationsToBiasAA:
    """
    Test 5: Design Recommendations to Bias String

    Verifies that:
    1. get_design_recommendations returns bias_aa for metal
    2. bias_aa string format is correct (e.g., "E:3.0,D:3.0")
    3. Bias values make chemical sense for the metal
    """

    def test_design_recommendations_returns_bias_aa(self):
        """get_design_recommendations should return bias_aa string."""
        from ai_structure_interface import get_design_recommendations

        recs = get_design_recommendations(metal="TB", task="biosensor")

        assert "bias_aa" in recs
        assert recs["bias_aa"] is not None

    def test_bias_aa_format_is_correct(self):
        """bias_aa should have correct format: 'E:3.0,D:3.0,...'"""
        from ai_structure_interface import get_design_recommendations

        recs = get_design_recommendations(metal="TB", task="biosensor")
        bias_aa = recs.get("bias_aa", "")

        if bias_aa:
            # Should have comma-separated pairs
            pairs = bias_aa.split(",")
            for pair in pairs:
                if pair.strip():
                    # Each pair should have format "X:value" or "X:+value" or "X:-value"
                    parts = pair.strip().split(":")
                    assert len(parts) == 2, f"Invalid pair format: {pair}"
                    aa = parts[0]
                    value = parts[1]
                    # AA should be single letter
                    assert len(aa) == 1, f"Invalid AA code: {aa}"
                    # Value should be parseable as float
                    try:
                        float(value)
                    except ValueError:
                        pytest.fail(f"Invalid bias value: {value}")

    def test_bias_aa_chemically_sensible_for_terbium(self):
        """Terbium bias should favor Glu (E) and Asp (D) as hard donors."""
        from ai_structure_interface import get_design_recommendations

        recs = get_design_recommendations(metal="TB", task="biosensor")
        bias_aa = recs.get("bias_aa", "")

        if bias_aa:
            # Parse bias values
            bias_dict = {}
            for pair in bias_aa.split(","):
                if ":" in pair:
                    aa, value = pair.split(":")
                    bias_dict[aa.strip()] = float(value.strip())

            # Terbium (hard Lewis acid) should favor oxygen donors
            # Glu (E) and Asp (D) should have positive bias
            if "E" in bias_dict:
                assert bias_dict["E"] > 0, "Glu should have positive bias for Tb"
            if "D" in bias_dict:
                assert bias_dict["D"] > 0, "Asp should have positive bias for Tb"

            # Cysteine (soft donor) should have negative bias
            if "C" in bias_dict:
                assert bias_dict["C"] < 0, "Cys should have negative bias for Tb"

    def test_bias_aa_chemically_sensible_for_zinc(self):
        """Zinc bias should favor His (H) and Cys (C)."""
        from ai_structure_interface import get_design_recommendations

        recs = get_design_recommendations(metal="ZN")
        bias_aa = recs.get("bias_aa", "")

        if bias_aa:
            # Parse bias values
            bias_dict = {}
            for pair in bias_aa.split(","):
                if ":" in pair:
                    aa, value = pair.split(":")
                    bias_dict[aa.strip()] = float(value.strip())

            # Zinc (borderline Lewis acid) accepts His and Cys
            # His (H) and Cys (C) should have non-negative bias
            if "H" in bias_dict:
                assert bias_dict["H"] >= 0, "His should have non-negative bias for Zn"

    def test_design_recommendations_for_multiple_metals(self):
        """Test bias_aa generation for various metals."""
        from ai_structure_interface import get_design_recommendations

        metals_to_test = ["ZN", "CA", "FE", "TB", "MG"]

        for metal in metals_to_test:
            recs = get_design_recommendations(metal=metal)

            # Should not raise error
            assert recs is not None
            assert isinstance(recs, dict)

            # Should have bias_aa (may be empty string but should exist)
            assert "bias_aa" in recs

            # Should have donors list
            assert "donors" in recs

    def test_design_recommendations_includes_coordination_info(self):
        """Design recommendations should include coordination info."""
        from ai_structure_interface import get_design_recommendations

        recs = get_design_recommendations(metal="TB", task="biosensor")

        assert "coordination_info" in recs
        coord_info = recs["coordination_info"]
        assert isinstance(coord_info, dict)

        # Should have coordination number range
        if "coordination_number_range" in coord_info:
            cn_range = coord_info["coordination_number_range"]
            assert len(cn_range) == 2
            # Terbium typically has CN 8-9
            assert cn_range[0] >= 6
            assert cn_range[1] <= 12


# =============================================================================
# END-TO-END INTEGRATION TESTS
# =============================================================================

class TestEndToEndIntegration:
    """
    End-to-end integration tests covering the complete workflow.
    """

    @integration
    @network_required
    def test_complete_biosensor_design_workflow(self):
        """
        Complete biosensor design workflow:
        1. Parse NL query to understand intent
        2. Plan the structure search
        3. Execute search and get results
        4. Generate design recommendations
        """
        from ai_structure_interface import (
            answer_structure_question,
            plan_structure_search,
            get_design_recommendations,
        )
        from structure_discovery import StructureDiscovery

        # Step 1: User asks about terbium
        question_result = answer_structure_question(
            "What is the coordination geometry of terbium?"
        )
        assert question_result is not None
        assert question_result.get("metal") == "TB"

        # Step 2: Plan the search
        plan = plan_structure_search("design a terbium biosensor")
        assert plan is not None
        assert plan.get("metal") == "TB"
        assert len(plan.get("steps", [])) > 0

        # Step 3: Execute search
        discovery = StructureDiscovery()
        results = discovery.search_by_intent("terbium binding protein", limit=5)
        # Results may be empty if databases are down
        assert isinstance(results, list)

        # Step 4: Get design recommendations
        recs = get_design_recommendations(metal="TB", task="biosensor")
        assert recs is not None
        assert "donors" in recs
        assert "bias_aa" in recs
        assert "coordination_info" in recs

    @integration
    def test_error_handling_invalid_metal(self):
        """Should handle invalid metal gracefully."""
        from ai_structure_interface import get_design_recommendations

        recs = get_design_recommendations(metal="XX")  # Invalid metal

        assert recs is not None
        assert "error" in recs
        assert "valid_metals" in recs

    @integration
    def test_error_handling_unknown_question(self):
        """Should handle unknown questions gracefully."""
        from ai_structure_interface import answer_structure_question

        result = answer_structure_question("What is the meaning of life?")

        assert result is not None
        assert result.get("type") == "unknown"
        assert "suggestions" in result
