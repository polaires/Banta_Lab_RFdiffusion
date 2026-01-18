# test_ai_interface.py
"""Tests for AI-callable structure interface."""
import pytest
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Check if requests is available for network tests
try:
    import requests
    NETWORK_AVAILABLE = True
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


class TestAIInterfaceImport:
    """Test AI interface can be imported."""

    def test_can_import_answer_structure_question(self):
        """Should be able to import answer_structure_question."""
        from ai_structure_interface import answer_structure_question
        assert answer_structure_question is not None

    def test_can_import_plan_structure_search(self):
        """Should be able to import plan_structure_search."""
        from ai_structure_interface import plan_structure_search
        assert plan_structure_search is not None

    def test_can_import_get_design_recommendations(self):
        """Should be able to import get_design_recommendations."""
        from ai_structure_interface import get_design_recommendations
        assert get_design_recommendations is not None


class TestAnswerQuestionMetalCoordination:
    """Test answer_structure_question for coordination queries."""

    def test_answer_question_metal_coordination_terbium(self):
        """'What is the coordination geometry of terbium?' returns coordination_number in [8,9]."""
        from ai_structure_interface import answer_structure_question

        result = answer_structure_question("What is the coordination geometry of terbium?")

        assert result is not None
        assert isinstance(result, dict)
        assert "answer" in result
        assert "coordination_number" in result

        # Terbium should have CN 8-9
        cn = result["coordination_number"]
        assert isinstance(cn, list)
        assert 8 in cn or 9 in cn
        # Check that at least one value is in expected range
        assert any(n in [8, 9] for n in cn)

    def test_answer_question_metal_coordination_zinc(self):
        """Zinc coordination query should return CN 4-6."""
        from ai_structure_interface import answer_structure_question

        result = answer_structure_question("What is the coordination geometry of zinc?")

        assert result is not None
        assert "coordination_number" in result
        cn = result["coordination_number"]
        # Zinc typically has CN 4-6
        assert 4 in cn or 6 in cn

    def test_answer_question_metal_coordination_calcium(self):
        """Calcium coordination query should return CN 6-8."""
        from ai_structure_interface import answer_structure_question

        result = answer_structure_question("What coordination number does calcium have?")

        assert result is not None
        assert "metal" in result
        assert result["metal"] == "CA"
        assert "coordination_number" in result

    def test_answer_question_includes_hsab_class(self):
        """Coordination answer should include HSAB classification."""
        from ai_structure_interface import answer_structure_question

        result = answer_structure_question("What is the coordination geometry of terbium?")

        assert "hsab_class" in result
        assert result["hsab_class"] == "hard"  # Lanthanides are hard acids

    def test_answer_question_includes_common_residues(self):
        """Coordination answer should include common residues."""
        from ai_structure_interface import answer_structure_question

        result = answer_structure_question("coordination of zinc")

        assert "common_residues" in result
        assert isinstance(result["common_residues"], list)
        # Zinc should include His, Cys
        assert "HIS" in result["common_residues"] or "CYS" in result["common_residues"]

    def test_answer_question_type_is_coordination(self):
        """Coordination questions should return type='coordination'."""
        from ai_structure_interface import answer_structure_question

        result = answer_structure_question("What is the coordination number of copper?")

        assert result.get("type") == "coordination"


class TestAnswerQuestionFindStructure:
    """Test answer_structure_question for structure finding queries."""

    def test_answer_question_find_structure_pqq_calcium(self):
        """'Find me a PQQ-calcium structure' returns structures."""
        from ai_structure_interface import answer_structure_question

        result = answer_structure_question("Find me a PQQ-calcium structure")

        assert result is not None
        assert isinstance(result, dict)
        assert "type" in result
        assert result["type"] == "find_structure"
        assert "structures" in result

        # Should find at least the curated structures or return empty list
        assert isinstance(result["structures"], list)

    def test_answer_question_find_structure_terbium(self):
        """'Find terbium structures' should search for structures."""
        from ai_structure_interface import answer_structure_question

        result = answer_structure_question("Find terbium structures")

        assert result is not None
        assert result["type"] == "find_structure"
        assert "structures" in result

    def test_answer_question_find_structure_includes_count(self):
        """Find structure answer should include count."""
        from ai_structure_interface import answer_structure_question

        result = answer_structure_question("search for zinc structures")

        assert "count" in result
        assert isinstance(result["count"], int)

    def test_answer_question_find_structure_extracts_ligand(self):
        """Should extract ligand from query."""
        from ai_structure_interface import answer_structure_question

        result = answer_structure_question("Find structures with citrate and terbium")

        assert result.get("ligand") == "CIT"
        assert result.get("metal") == "TB"


class TestAnswerQuestionDonors:
    """Test answer_structure_question for donor queries."""

    def test_answer_question_donors_terbium(self):
        """'What residues coordinate terbium?' returns donor info."""
        from ai_structure_interface import answer_structure_question

        result = answer_structure_question("What residues coordinate terbium?")

        assert result is not None
        assert result.get("type") == "donors"
        assert "donor_weights" in result
        assert "common_residues" in result

        # Terbium prefers O donors
        donors = result["donor_weights"]
        assert "O" in donors
        assert donors["O"] > donors.get("S", -10)  # O should be preferred over S

    def test_answer_question_donors_zinc(self):
        """Zinc donor query should include His, Cys."""
        from ai_structure_interface import answer_structure_question

        result = answer_structure_question("Which amino acids coordinate zinc?")

        assert result is not None
        assert "common_residues" in result
        residues = result["common_residues"]
        # Zinc commonly uses His and Cys
        assert "HIS" in residues or "CYS" in residues

    def test_answer_question_donors_includes_bias_string(self):
        """Donor answer should include LigandMPNN bias string."""
        from ai_structure_interface import answer_structure_question

        result = answer_structure_question("Which amino acids coordinate calcium?")

        assert "bias_aa" in result
        assert isinstance(result["bias_aa"], str)
        # Calcium should bias towards E, D (Glu, Asp)
        assert "E:" in result["bias_aa"] or "D:" in result["bias_aa"]


class TestAnswerQuestionTemplates:
    """Test answer_structure_question for template queries."""

    def test_answer_question_template_terbium(self):
        """'Get a template for terbium' should return template info."""
        from ai_structure_interface import answer_structure_question

        result = answer_structure_question("Get a starting structure template for terbium")

        assert result is not None
        assert result.get("type") == "template"

    def test_answer_question_template_citrate_tb(self):
        """Template query for citrate-terbium should find citrate_tb template."""
        from ai_structure_interface import answer_structure_question

        result = answer_structure_question("Find a template for citrate terbium complex")

        assert result is not None
        if result.get("template_name"):
            assert "citrate" in result["template_name"].lower()

    def test_answer_question_template_lists_available(self):
        """Template query should list available templates."""
        from ai_structure_interface import answer_structure_question

        result = answer_structure_question("What templates are available?")

        assert result is not None
        assert "available_templates" in result or "template_name" in result


class TestAnswerQuestionUnknown:
    """Test answer_structure_question for unknown queries."""

    def test_answer_question_unknown_returns_suggestions(self):
        """Unknown questions should return suggestions."""
        from ai_structure_interface import answer_structure_question

        result = answer_structure_question("This is a random question")

        assert result is not None
        assert result.get("type") == "unknown" or "suggestions" in result or "metal" in result


class TestPlanStructureSearch:
    """Test plan_structure_search function."""

    def test_plan_structure_search_terbium_biosensor(self):
        """'I need to design a terbium biosensor' returns metal='TB'."""
        from ai_structure_interface import plan_structure_search

        result = plan_structure_search("I need to design a terbium biosensor")

        assert result is not None
        assert isinstance(result, dict)
        assert result.get("metal") == "TB"
        assert "steps" in result
        assert isinstance(result["steps"], list)
        assert len(result["steps"]) > 0

    def test_plan_structure_search_detects_task_type(self):
        """Should detect biosensor task type."""
        from ai_structure_interface import plan_structure_search

        result = plan_structure_search("I need to design a terbium biosensor")

        assert result.get("task_type") == "biosensor"

    def test_plan_structure_search_detects_ligand(self):
        """Should detect ligand in task description."""
        from ai_structure_interface import plan_structure_search

        result = plan_structure_search("Design a citrate-binding terbium sensor")

        assert result.get("ligand") == "CIT"
        assert result.get("metal") == "TB"

    def test_plan_structure_search_enzyme_task(self):
        """Should detect enzyme task type."""
        from ai_structure_interface import plan_structure_search

        result = plan_structure_search("I want to engineer a zinc enzyme")

        assert result.get("task_type") == "enzyme"
        assert result.get("metal") == "ZN"

    def test_plan_structure_search_includes_steps(self):
        """Plan should include actionable steps."""
        from ai_structure_interface import plan_structure_search

        result = plan_structure_search("Create a calcium binding site")

        steps = result.get("steps", [])
        assert len(steps) >= 2
        # Should include a search step
        assert any("search" in step.lower() or "find" in step.lower() for step in steps)

    def test_plan_structure_search_includes_databases(self):
        """Plan should include database recommendations."""
        from ai_structure_interface import plan_structure_search

        result = plan_structure_search("Design a terbium biosensor")

        assert "databases" in result
        assert isinstance(result["databases"], list)
        assert len(result["databases"]) >= 1

    def test_plan_structure_search_lanthanide_strategy(self):
        """Lanthanide search should use specialized strategy."""
        from ai_structure_interface import plan_structure_search

        result = plan_structure_search("europium luminescent sensor design")

        assert result.get("metal") == "EU"
        assert result.get("search_strategy") == "lanthanide_specialized"


class TestGetDesignRecommendations:
    """Test get_design_recommendations function."""

    def test_get_design_recommendations_tb_biosensor(self):
        """metal='TB', task='biosensor' returns donors, templates, bias_aa."""
        from ai_structure_interface import get_design_recommendations

        result = get_design_recommendations(metal="TB", task="biosensor")

        assert result is not None
        assert isinstance(result, dict)

        # Check required fields
        assert "donors" in result
        assert "templates" in result
        assert "bias_aa" in result

        # Check donors for terbium (hard acid, should prefer Glu, Asp)
        donors = result["donors"]
        assert isinstance(donors, list)
        assert "Glu" in donors or "Asp" in donors

        # Check bias string is non-empty
        assert isinstance(result["bias_aa"], str)
        assert len(result["bias_aa"]) > 0

    def test_get_design_recommendations_includes_avoid_residues(self):
        """Should include residues to avoid."""
        from ai_structure_interface import get_design_recommendations

        result = get_design_recommendations(metal="TB")

        assert "avoid_residues" in result
        # Terbium (hard acid) should avoid Cys (soft donor)
        assert "Cys" in result["avoid_residues"]

    def test_get_design_recommendations_zinc(self):
        """Zinc recommendations should include His, Cys donors."""
        from ai_structure_interface import get_design_recommendations

        result = get_design_recommendations(metal="ZN")

        donors = result.get("donors", [])
        # Zinc is borderline, can use N and S donors
        assert "His" in donors or "Cys" in donors

    def test_get_design_recommendations_includes_hsab_class(self):
        """Should include HSAB classification."""
        from ai_structure_interface import get_design_recommendations

        result = get_design_recommendations(metal="CA")

        assert "hsab_class" in result
        assert result["hsab_class"] == "hard"  # Calcium is hard acid

    def test_get_design_recommendations_includes_coordination_info(self):
        """Should include coordination information."""
        from ai_structure_interface import get_design_recommendations

        result = get_design_recommendations(metal="TB")

        assert "coordination_info" in result
        coord_info = result["coordination_info"]
        assert "coordination_number_range" in coord_info
        # Terbium has CN 8-9
        cn_range = coord_info["coordination_number_range"]
        assert 8 in cn_range or cn_range[0] == 8 or cn_range[1] == 9

    def test_get_design_recommendations_with_ligand(self):
        """Should find templates when ligand is specified."""
        from ai_structure_interface import get_design_recommendations

        result = get_design_recommendations(metal="TB", ligand="CIT")

        templates = result.get("templates", [])
        # Should find citrate_tb template
        template_names = [t["name"] for t in templates]
        assert any("citrate" in name.lower() for name in template_names)

    def test_get_design_recommendations_task_recommendations(self):
        """Should include task-specific recommendations."""
        from ai_structure_interface import get_design_recommendations

        result = get_design_recommendations(metal="EU", task="biosensor")

        task_recs = result.get("task_recommendations", [])
        # Should mention luminescence for Eu biosensor
        assert len(task_recs) > 0
        assert any("luminescent" in r.lower() or "biosensor" in r.lower() for r in task_recs)

    def test_get_design_recommendations_unknown_metal(self):
        """Should return error for unknown metal."""
        from ai_structure_interface import get_design_recommendations

        result = get_design_recommendations(metal="XX")

        assert "error" in result
        assert "valid_metals" in result


class TestMetalExtraction:
    """Test internal metal extraction helpers."""

    def test_extract_metal_from_question_names(self):
        """Should extract metal from natural language names."""
        from ai_structure_interface import _extract_metal_from_question

        assert _extract_metal_from_question("terbium binding") == "TB"
        assert _extract_metal_from_question("zinc finger") == "ZN"
        assert _extract_metal_from_question("calcium ions") == "CA"
        assert _extract_metal_from_question("copper enzyme") == "CU"

    def test_extract_metal_from_question_codes(self):
        """Should extract metal from element codes."""
        from ai_structure_interface import _extract_metal_from_question

        assert _extract_metal_from_question("TB binding site") == "TB"
        assert _extract_metal_from_question("ZN coordination") == "ZN"

    def test_extract_metal_case_insensitive(self):
        """Metal extraction should be case insensitive."""
        from ai_structure_interface import _extract_metal_from_question

        assert _extract_metal_from_question("TERBIUM") == "TB"
        assert _extract_metal_from_question("Zinc") == "ZN"
        assert _extract_metal_from_question("EUROPIUM") == "EU"


class TestLigandExtraction:
    """Test internal ligand extraction helpers."""

    def test_extract_ligand_from_question_names(self):
        """Should extract ligand from natural language names."""
        from ai_structure_interface import _extract_ligand_from_question

        assert _extract_ligand_from_question("citrate complex") == "CIT"
        assert _extract_ligand_from_question("pqq dehydrogenase") == "PQQ"
        assert _extract_ligand_from_question("atp binding") == "ATP"

    def test_extract_ligand_from_question_codes(self):
        """Should extract ligand from PDB codes."""
        from ai_structure_interface import _extract_ligand_from_question

        assert _extract_ligand_from_question("CIT ligand") == "CIT"
        assert _extract_ligand_from_question("find PQQ") == "PQQ"


class TestQuestionTypeDetection:
    """Test internal question type detection."""

    def test_is_coordination_question(self):
        """Should detect coordination questions."""
        from ai_structure_interface import _is_coordination_question

        assert _is_coordination_question("what is the coordination geometry")
        assert _is_coordination_question("coordination number of zinc")
        assert _is_coordination_question("how many ligands bind")
        assert not _is_coordination_question("find structures")

    def test_is_find_structure_question(self):
        """Should detect find structure questions."""
        from ai_structure_interface import _is_find_structure_question

        assert _is_find_structure_question("find me a structure")
        assert _is_find_structure_question("search for proteins")
        assert _is_find_structure_question("get pdb structures")
        assert not _is_find_structure_question("coordination number")

    def test_is_donor_question(self):
        """Should detect donor questions."""
        from ai_structure_interface import _is_donor_question

        assert _is_donor_question("what residues coordinate zinc")
        assert _is_donor_question("which amino acids coordinate")
        assert _is_donor_question("preferred donors for terbium")
        assert not _is_donor_question("find structures")

    def test_is_template_question(self):
        """Should detect template questions."""
        from ai_structure_interface import _is_template_question

        assert _is_template_question("get a template for")
        assert _is_template_question("starting structure for design")
        assert _is_template_question("scaffold for biosensor")
        assert not _is_template_question("coordination number")


class TestIntegration:
    """Integration tests for the AI interface."""

    def test_full_workflow_terbium_biosensor_design(self):
        """Full workflow: plan -> recommendations -> answer questions."""
        from ai_structure_interface import (
            answer_structure_question,
            plan_structure_search,
            get_design_recommendations,
        )

        # Step 1: Plan the search
        plan = plan_structure_search("I want to design a terbium-citrate biosensor")
        assert plan["metal"] == "TB"
        assert plan["ligand"] == "CIT"
        assert plan["task_type"] == "biosensor"

        # Step 2: Get design recommendations
        recs = get_design_recommendations(metal="TB", task="biosensor", ligand="CIT")
        assert "donors" in recs
        assert "bias_aa" in recs
        assert len(recs["donors"]) > 0

        # Step 3: Answer coordination question
        coord_answer = answer_structure_question(
            "What is the coordination geometry of terbium?"
        )
        assert coord_answer["coordination_number"] == [8, 9]
        assert coord_answer["hsab_class"] == "hard"

        # Step 4: Answer donor question
        donor_answer = answer_structure_question(
            "What residues should coordinate terbium?"
        )
        assert "common_residues" in donor_answer

    def test_workflow_zinc_enzyme_design(self):
        """Workflow for zinc enzyme design."""
        from ai_structure_interface import (
            answer_structure_question,
            plan_structure_search,
            get_design_recommendations,
        )

        # Plan
        plan = plan_structure_search("Design a catalytic zinc site")
        assert plan["metal"] == "ZN"

        # Get recommendations
        recs = get_design_recommendations(metal="ZN", task="enzyme")
        assert "His" in recs["donors"] or "Cys" in recs["donors"]

        # Answer questions
        answer = answer_structure_question("coordination of zinc")
        assert 4 in answer["coordination_number"] or 6 in answer["coordination_number"]

    @network_required
    def test_find_structures_with_network(self):
        """Find structures should work with network access."""
        from ai_structure_interface import answer_structure_question

        result = answer_structure_question("Find zinc structures in PDB")

        assert result["type"] == "find_structure"
        # With network, should find structures
        if RCSB_AVAILABLE:
            assert result["count"] >= 0  # May find structures


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_question(self):
        """Should handle empty question gracefully."""
        from ai_structure_interface import answer_structure_question

        result = answer_structure_question("")

        assert result is not None
        assert "answer" in result

    def test_no_metal_specified_coordination(self):
        """Should handle coordination question without metal."""
        from ai_structure_interface import answer_structure_question

        result = answer_structure_question("What is the coordination geometry?")

        assert result is not None
        assert "error" in result or "suggestions" in result or "answer" in result

    def test_invalid_metal_recommendations(self):
        """Should handle invalid metal gracefully."""
        from ai_structure_interface import get_design_recommendations

        result = get_design_recommendations(metal="INVALID")

        assert "error" in result

    def test_case_sensitivity(self):
        """Should handle different cases."""
        from ai_structure_interface import answer_structure_question

        result1 = answer_structure_question("TERBIUM coordination")
        result2 = answer_structure_question("terbium coordination")
        result3 = answer_structure_question("Terbium coordination")

        assert result1["metal"] == result2["metal"] == result3["metal"] == "TB"


class TestResponseFormat:
    """Test response format consistency."""

    def test_answer_always_has_answer_field(self):
        """All answers should have 'answer' field."""
        from ai_structure_interface import answer_structure_question

        questions = [
            "What is the coordination geometry of zinc?",
            "Find terbium structures",
            "What donors does calcium prefer?",
            "Get a template for biosensor design",
            "Random question",
        ]

        for q in questions:
            result = answer_structure_question(q)
            assert "answer" in result, f"Missing 'answer' for: {q}"

    def test_answer_always_has_type_field(self):
        """All answers should have 'type' field."""
        from ai_structure_interface import answer_structure_question

        questions = [
            "What is the coordination geometry of zinc?",
            "Find terbium structures",
            "What donors does calcium prefer?",
        ]

        for q in questions:
            result = answer_structure_question(q)
            assert "type" in result, f"Missing 'type' for: {q}"

    def test_plan_always_has_steps(self):
        """All plans should have 'steps' field."""
        from ai_structure_interface import plan_structure_search

        tasks = [
            "Design a terbium biosensor",
            "Create zinc enzyme",
            "Make calcium binding site",
        ]

        for t in tasks:
            result = plan_structure_search(t)
            assert "steps" in result, f"Missing 'steps' for: {t}"
            assert isinstance(result["steps"], list)

    def test_recommendations_always_has_donors(self):
        """All recommendations should have 'donors' field."""
        from ai_structure_interface import get_design_recommendations

        metals = ["TB", "ZN", "CA", "FE", "CU"]

        for m in metals:
            result = get_design_recommendations(metal=m)
            assert "donors" in result, f"Missing 'donors' for: {m}"
