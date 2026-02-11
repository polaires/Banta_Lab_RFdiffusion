"""
Tests for scaffold_search module.

Tests the search, validation, scoring, and ranking pipeline with mocked
network calls to metal_site_fetcher.
"""

import sys
import os
import unittest
from unittest.mock import patch, MagicMock

# Ensure serverless directory is on path
sys.path.insert(0, os.path.join(
    "G:", os.sep, "Github_local_repo", "Banta_Lab_RFdiffusion", "backend", "serverless"
))

from scaffold_search import (
    ScaffoldCandidate,
    ScaffoldSearchResult,
    search_scaffold_candidates,
    rank_candidates,
    _score_candidate,
    _get_hsab_compatible_metals,
    _resolve_ligand_code,
    _search_pdb_hits,
    _validate_candidate,
    SCAFFOLD_SEARCH_THRESHOLD,
)


class TestResolveLeadCode(unittest.TestCase):
    """Test ligand code resolution."""

    def test_known_ligand_name(self):
        self.assertEqual(_resolve_ligand_code("citrate"), "CIT")

    def test_already_code(self):
        self.assertEqual(_resolve_ligand_code("CIT"), "CIT")

    def test_unknown_fallback(self):
        # Unknown names get uppercased
        result = _resolve_ligand_code("xyz")
        self.assertEqual(result, "XYZ")


class TestHSABCompatibleMetals(unittest.TestCase):
    """Test compatible metal discovery."""

    def test_tb_compatible(self):
        compatible = _get_hsab_compatible_metals("TB")
        # TB is hard, O-preferring — CA, MG, EU etc. should be compatible
        # Exact list depends on METAL_DATABASE content
        self.assertIsInstance(compatible, list)
        # TB should not be in its own compatible list
        self.assertNotIn("TB", compatible)

    def test_unknown_metal(self):
        compatible = _get_hsab_compatible_metals("XX")
        # Unknown metal — no compatible metals
        self.assertEqual(compatible, [])


class TestSearchPDBHits(unittest.TestCase):
    """Test two-phase PDB search with mocked API."""

    @patch("scaffold_search.query_metal_ligand_sites")
    def test_exact_match_found(self, mock_query):
        mock_query.return_value = [
            {"pdb_id": "3C9H", "metal": "TB", "ligand": "CIT", "resolution": 1.8},
        ]
        hits = _search_pdb_hits("CIT", "TB", limit=10)
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["pdb_id"], "3C9H")

    @patch("scaffold_search.query_metal_ligand_sites")
    @patch("scaffold_search._get_hsab_compatible_metals")
    def test_compatible_metal_fallback(self, mock_compat, mock_query):
        # Phase 1: no exact hits
        # Phase 2: compatible metal CA has hits
        mock_query.side_effect = [
            [],  # exact TB search
            [{"pdb_id": "1XYZ", "metal": "CA", "ligand": "CIT", "resolution": 2.0}],
            [],  # MG search — no hits
        ]
        mock_compat.return_value = ["CA", "MG"]
        hits = _search_pdb_hits("CIT", "TB", limit=10)
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["metal"], "CA")

    @patch("scaffold_search.query_metal_ligand_sites")
    def test_deduplication(self, mock_query):
        mock_query.return_value = [
            {"pdb_id": "3C9H", "metal": "TB", "ligand": "CIT", "resolution": 1.8},
            {"pdb_id": "3C9H", "metal": "TB", "ligand": "CIT", "resolution": 1.8},
        ]
        hits = _search_pdb_hits("CIT", "TB", limit=10)
        self.assertEqual(len(hits), 1)


class TestValidateCandidate(unittest.TestCase):
    """Test spatial proximity validation."""

    @patch("scaffold_search.find_metal_ligand_active_site")
    def test_valid_site(self, mock_find):
        mock_find.return_value = {
            "pdb_id": "3C9H",
            "metal": "TB",
            "metal_chain": "A",
            "metal_resnum": "201",
            "metal_coords": (10.0, 20.0, 30.0),
            "ligand": "CIT",
            "ligand_distance": 2.4,
            "coordination_number": 9,
            "coordinating_atoms": [
                {"atom_name": "OE1", "res_name": "GLU", "distance": 2.3},
                {"atom_name": "OD1", "res_name": "ASP", "distance": 2.4},
            ],
            "protein_donors": ["A10:GLU OE1@2.30", "A42:ASP OD1@2.40"],
            "ligand_donors": ["O1@2.20", "O5@2.50", "O6@2.60"],
        }
        candidate = _validate_candidate("3C9H", "TB", "TB", "CIT")
        self.assertIsNotNone(candidate)
        self.assertEqual(candidate.pdb_id, "3C9H")
        self.assertFalse(candidate.needs_substitution)
        self.assertEqual(candidate.coordination_number, 9)
        self.assertEqual(len(candidate.ligand_donors), 3)

    @patch("scaffold_search.find_metal_ligand_active_site")
    def test_no_site_found(self, mock_find):
        mock_find.return_value = None
        candidate = _validate_candidate("XXXX", "TB", "TB", "CIT")
        self.assertIsNone(candidate)

    @patch("scaffold_search.find_metal_ligand_active_site")
    def test_substitution_flagged(self, mock_find):
        mock_find.return_value = {
            "coordination_number": 7,
            "protein_donors": ["A10:GLU OE1@2.30"],
            "ligand_donors": ["O1@2.20"],
            "ligand_distance": 2.4,
            "coordinating_atoms": [],
        }
        candidate = _validate_candidate("1XYZ", "CA", "TB", "CIT")
        self.assertIsNotNone(candidate)
        self.assertTrue(candidate.needs_substitution)


class TestScoring(unittest.TestCase):
    """Test the composite scoring function."""

    def _make_candidate(self, **kwargs):
        defaults = {
            "pdb_id": "TEST",
            "source_metal": "TB",
            "target_metal": "TB",
            "ligand_code": "CIT",
            "coordination_number": 9,
            "protein_donors": ["A10:GLU OE1@2.30", "A42:ASP OD1@2.40",
                             "A77:ASN ND2@2.50", "A99:GLU OE2@2.35"],
            "ligand_donors": ["O1@2.20", "O5@2.50", "O6@2.60"],
            "ligand_distance": 2.4,
            "coordinating_atoms": [
                {"distance": 2.3}, {"distance": 2.4},
                {"distance": 2.5}, {"distance": 2.35},
            ],
            "resolution": 1.5,
        }
        defaults.update(kwargs)
        return ScaffoldCandidate(**defaults)

    def test_high_quality_candidate_scores_well(self):
        candidate = self._make_candidate()
        score = _score_candidate(candidate, "TB")
        # Good CN, good donors, good resolution — should score high
        self.assertGreater(score, 60.0)
        self.assertEqual(score, candidate.total_score)

    def test_low_cn_reduces_score(self):
        good = self._make_candidate(coordination_number=9)
        bad = self._make_candidate(coordination_number=3)
        _score_candidate(good, "TB")
        _score_candidate(bad, "TB")
        self.assertGreater(good.score_cn, bad.score_cn)

    def test_no_ligand_donors_reduces_score(self):
        good = self._make_candidate(ligand_donors=["O1@2.2", "O5@2.5", "O6@2.6"])
        bad = self._make_candidate(ligand_donors=[])
        _score_candidate(good, "TB")
        _score_candidate(bad, "TB")
        self.assertGreater(good.score_lig_donors, bad.score_lig_donors)

    def test_poor_resolution_reduces_score(self):
        good = self._make_candidate(resolution=1.0)
        bad = self._make_candidate(resolution=4.0)
        _score_candidate(good, "TB")
        _score_candidate(bad, "TB")
        self.assertGreater(good.score_resolution, bad.score_resolution)

    def test_score_range(self):
        candidate = self._make_candidate()
        score = _score_candidate(candidate, "TB")
        self.assertGreaterEqual(score, 0.0)
        self.assertLessEqual(score, 100.0)


class TestRankCandidates(unittest.TestCase):
    """Test ranking sorts by total_score descending."""

    def test_ranking_order(self):
        c1 = ScaffoldCandidate(
            pdb_id="A", source_metal="TB", target_metal="TB", ligand_code="CIT",
            coordination_number=9, resolution=1.5,
            protein_donors=["A10:GLU OE1@2.3", "A42:ASP OD1@2.4",
                          "A77:ASN ND2@2.5", "A99:GLU OE2@2.35"],
            ligand_donors=["O1@2.2", "O5@2.5", "O6@2.6"],
            coordinating_atoms=[{"distance": 2.3}, {"distance": 2.4}],
        )
        c2 = ScaffoldCandidate(
            pdb_id="B", source_metal="TB", target_metal="TB", ligand_code="CIT",
            coordination_number=3, resolution=4.0,
            protein_donors=[], ligand_donors=[],
            coordinating_atoms=[],
        )
        ranked = rank_candidates([c2, c1], "TB")
        self.assertEqual(ranked[0].pdb_id, "A")
        self.assertGreater(ranked[0].total_score, ranked[1].total_score)


class TestSearchScaffoldCandidates(unittest.TestCase):
    """Integration test for the main entry point."""

    @patch("scaffold_search.find_metal_ligand_active_site")
    @patch("scaffold_search.query_metal_ligand_sites")
    def test_full_pipeline_recommends_scaffold(self, mock_query, mock_find):
        mock_query.return_value = [
            {"pdb_id": "3C9H", "metal": "TB", "ligand": "CIT", "resolution": 1.8},
        ]
        mock_find.return_value = {
            "coordination_number": 9,
            "protein_donors": ["A10:GLU OE1@2.30", "A42:ASP OD1@2.40",
                             "A77:ASN ND2@2.50", "A99:GLU OE2@2.35"],
            "ligand_donors": ["O1@2.20", "O5@2.50", "O6@2.60"],
            "ligand_distance": 2.4,
            "coordinating_atoms": [
                {"distance": 2.3}, {"distance": 2.4},
                {"distance": 2.5}, {"distance": 2.35},
            ],
        }

        result = search_scaffold_candidates(
            metal="TB", ligand_name="citrate"
        )

        self.assertTrue(result.searched)
        self.assertEqual(result.query_metal, "TB")
        self.assertEqual(result.ligand_code, "CIT")
        self.assertEqual(result.num_pdb_hits, 1)
        self.assertEqual(result.num_validated, 1)
        self.assertIsNotNone(result.best_candidate)
        self.assertEqual(result.best_candidate.pdb_id, "3C9H")
        # With good data, score should be above threshold
        self.assertGreaterEqual(result.best_candidate.total_score, SCAFFOLD_SEARCH_THRESHOLD)
        self.assertEqual(result.recommended_action, "scaffold")

    @patch("scaffold_search.query_metal_ligand_sites")
    def test_no_hits_returns_de_novo(self, mock_query):
        mock_query.return_value = []

        result = search_scaffold_candidates(
            metal="TB", ligand_name="obscureligand"
        )

        self.assertTrue(result.searched)
        self.assertEqual(result.num_pdb_hits, 0)
        self.assertEqual(result.recommended_action, "de_novo")
        self.assertIn("No PDB structures found", result.reason)

    @patch("scaffold_search.find_metal_ligand_active_site")
    @patch("scaffold_search.query_metal_ligand_sites")
    def test_all_validation_fails_returns_de_novo(self, mock_query, mock_find):
        mock_query.return_value = [
            {"pdb_id": "XXXX", "metal": "TB", "ligand": "CIT", "resolution": 2.0},
        ]
        mock_find.return_value = None  # No site found

        result = search_scaffold_candidates(
            metal="TB", ligand_name="citrate"
        )

        self.assertEqual(result.num_pdb_hits, 1)
        self.assertEqual(result.num_validated, 0)
        self.assertEqual(result.recommended_action, "de_novo")

    def test_to_dict_serialization(self):
        result = ScaffoldSearchResult(
            searched=True,
            query_metal="TB",
            query_ligand="citrate",
            ligand_code="CIT",
        )
        d = result.to_dict()
        self.assertTrue(d["searched"])
        self.assertEqual(d["query_metal"], "TB")
        self.assertIsNone(d["best_candidate"])

    def test_candidate_to_dict(self):
        c = ScaffoldCandidate(
            pdb_id="3C9H",
            source_metal="CA",
            target_metal="TB",
            ligand_code="CIT",
            needs_substitution=True,
            total_score=65.3,
        )
        d = c.to_dict()
        self.assertEqual(d["pdb_id"], "3C9H")
        self.assertTrue(d["needs_substitution"])
        self.assertEqual(d["total_score"], 65.3)


if __name__ == "__main__":
    unittest.main()
