# test_coordination_knowledge.py
"""Tests for coordination motifs knowledge base and geometry validation."""
import pytest
import sys
import os
import math

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class TestCoordinationMotifsImport:
    """Test coordination motifs can be imported."""

    def test_can_import_coordination_motifs(self):
        """Should be able to import COORDINATION_MOTIFS."""
        from metal_chemistry import COORDINATION_MOTIFS
        assert COORDINATION_MOTIFS is not None
        assert isinstance(COORDINATION_MOTIFS, dict)

    def test_can_import_ideal_geometries(self):
        """Should be able to import IDEAL_GEOMETRIES."""
        from metal_chemistry import IDEAL_GEOMETRIES
        assert IDEAL_GEOMETRIES is not None
        assert isinstance(IDEAL_GEOMETRIES, dict)

    def test_can_import_get_coordination_template(self):
        """Should be able to import get_coordination_template."""
        from metal_chemistry import get_coordination_template
        assert get_coordination_template is not None
        assert callable(get_coordination_template)

    def test_can_import_get_ideal_geometry(self):
        """Should be able to import get_ideal_geometry."""
        from metal_chemistry import get_ideal_geometry
        assert get_ideal_geometry is not None
        assert callable(get_ideal_geometry)

    def test_can_import_validate_coordination(self):
        """Should be able to import validate_coordination."""
        from metal_chemistry import validate_coordination
        assert validate_coordination is not None
        assert callable(validate_coordination)


class TestEFHandMotif:
    """Test EF-hand coordination motif."""

    def test_ef_hand_motif_exists(self):
        """ef_hand should exist in COORDINATION_MOTIFS."""
        from metal_chemistry import COORDINATION_MOTIFS

        assert "ef_hand" in COORDINATION_MOTIFS
        ef_hand = COORDINATION_MOTIFS["ef_hand"]
        assert "CA" in ef_hand["metals"]
        assert ef_hand["coordination_number"] == 7
        assert ef_hand["geometry"] == "pentagonal_bipyramidal"

    def test_ef_hand_has_example_pdbs(self):
        """ef_hand should have example PDB IDs."""
        from metal_chemistry import COORDINATION_MOTIFS

        ef_hand = COORDINATION_MOTIFS["ef_hand"]
        assert "example_pdbs" in ef_hand
        assert "1CLL" in ef_hand["example_pdbs"]
        assert "3CLN" in ef_hand["example_pdbs"]


class TestZincFingerMotif:
    """Test zinc finger coordination motif."""

    def test_zinc_finger_motif_exists(self):
        """zinc_finger should exist in COORDINATION_MOTIFS."""
        from metal_chemistry import COORDINATION_MOTIFS

        assert "zinc_finger" in COORDINATION_MOTIFS
        zinc_finger = COORDINATION_MOTIFS["zinc_finger"]
        assert "ZN" in zinc_finger["metals"]
        assert zinc_finger["coordination_number"] == 4
        assert zinc_finger["geometry"] == "tetrahedral"

    def test_zinc_finger_has_example_pdbs(self):
        """zinc_finger should have example PDB IDs."""
        from metal_chemistry import COORDINATION_MOTIFS

        zinc_finger = COORDINATION_MOTIFS["zinc_finger"]
        assert "1ZNF" in zinc_finger["example_pdbs"]
        assert "1AAY" in zinc_finger["example_pdbs"]


class TestLanthanideSite:
    """Test lanthanide site coordination motif."""

    def test_lanthanide_site_motif_exists(self):
        """lanthanide_site should exist in COORDINATION_MOTIFS."""
        from metal_chemistry import COORDINATION_MOTIFS

        assert "lanthanide_site" in COORDINATION_MOTIFS
        lnt_site = COORDINATION_MOTIFS["lanthanide_site"]
        assert "TB" in lnt_site["metals"]
        assert "EU" in lnt_site["metals"]
        assert "GD" in lnt_site["metals"]
        assert "LA" in lnt_site["metals"]

    def test_lanthanide_site_coordination_numbers(self):
        """lanthanide_site should support CN 8 and 9."""
        from metal_chemistry import COORDINATION_MOTIFS

        lnt_site = COORDINATION_MOTIFS["lanthanide_site"]
        assert isinstance(lnt_site["coordination_number"], list)
        assert 8 in lnt_site["coordination_number"]
        assert 9 in lnt_site["coordination_number"]

    def test_lanthanide_site_geometries(self):
        """lanthanide_site should have appropriate geometries."""
        from metal_chemistry import COORDINATION_MOTIFS

        lnt_site = COORDINATION_MOTIFS["lanthanide_site"]
        assert isinstance(lnt_site["geometry"], list)
        assert "square_antiprism" in lnt_site["geometry"]
        assert "tricapped_trigonal_prism" in lnt_site["geometry"]


class TestOtherMotifs:
    """Test other coordination motifs."""

    def test_his_tag_motif_exists(self):
        """his_tag should exist with multiple supported metals."""
        from metal_chemistry import COORDINATION_MOTIFS

        assert "his_tag" in COORDINATION_MOTIFS
        his_tag = COORDINATION_MOTIFS["his_tag"]
        assert "NI" in his_tag["metals"]
        assert "CO" in his_tag["metals"]
        assert "ZN" in his_tag["metals"]
        assert "CU" in his_tag["metals"]
        assert his_tag["coordination_number"] == 6
        assert his_tag["geometry"] == "octahedral"

    def test_rubredoxin_motif_exists(self):
        """rubredoxin should exist with FE."""
        from metal_chemistry import COORDINATION_MOTIFS

        assert "rubredoxin" in COORDINATION_MOTIFS
        rubredoxin = COORDINATION_MOTIFS["rubredoxin"]
        assert "FE" in rubredoxin["metals"]
        assert rubredoxin["coordination_number"] == 4
        assert rubredoxin["geometry"] == "tetrahedral"

    def test_carbonic_anhydrase_motif_exists(self):
        """carbonic_anhydrase should exist with ZN."""
        from metal_chemistry import COORDINATION_MOTIFS

        assert "carbonic_anhydrase" in COORDINATION_MOTIFS
        ca = COORDINATION_MOTIFS["carbonic_anhydrase"]
        assert "ZN" in ca["metals"]
        assert ca["coordination_number"] == 4
        assert ca["geometry"] == "tetrahedral"
        assert "1CA2" in ca["example_pdbs"]


class TestIdealGeometries:
    """Test ideal geometry positions."""

    def test_all_geometries_present(self):
        """All required geometries should be present."""
        from metal_chemistry import IDEAL_GEOMETRIES

        required = [
            "linear", "trigonal_planar", "tetrahedral", "square_planar",
            "trigonal_bipyramidal", "square_pyramidal", "octahedral",
            "pentagonal_bipyramidal", "square_antiprism", "tricapped_trigonal_prism"
        ]
        for geom in required:
            assert geom in IDEAL_GEOMETRIES, f"Missing geometry: {geom}"

    def test_linear_geometry(self):
        """Linear geometry should have 2 positions."""
        from metal_chemistry import IDEAL_GEOMETRIES

        linear = IDEAL_GEOMETRIES["linear"]
        assert linear["coordination_number"] == 2
        assert len(linear["positions"]) == 2

    def test_tetrahedral_geometry(self):
        """Tetrahedral geometry should have 4 positions."""
        from metal_chemistry import IDEAL_GEOMETRIES

        tetra = IDEAL_GEOMETRIES["tetrahedral"]
        assert tetra["coordination_number"] == 4
        assert len(tetra["positions"]) == 4

    def test_octahedral_geometry(self):
        """Octahedral geometry should have 6 positions on unit sphere."""
        from metal_chemistry import IDEAL_GEOMETRIES

        octa = IDEAL_GEOMETRIES["octahedral"]
        assert octa["coordination_number"] == 6
        assert len(octa["positions"]) == 6

        # Check positions are on unit sphere (distance from origin = 1)
        for x, y, z in octa["positions"]:
            dist = math.sqrt(x*x + y*y + z*z)
            assert abs(dist - 1.0) < 0.01, f"Position ({x}, {y}, {z}) not on unit sphere"

    def test_square_antiprism_geometry(self):
        """Square antiprism should have 8 positions."""
        from metal_chemistry import IDEAL_GEOMETRIES

        sap = IDEAL_GEOMETRIES["square_antiprism"]
        assert sap["coordination_number"] == 8
        assert len(sap["positions"]) == 8

    def test_tricapped_trigonal_prism_geometry(self):
        """Tricapped trigonal prism should have 9 positions."""
        from metal_chemistry import IDEAL_GEOMETRIES

        ttp = IDEAL_GEOMETRIES["tricapped_trigonal_prism"]
        assert ttp["coordination_number"] == 9
        assert len(ttp["positions"]) == 9


class TestGetCoordinationTemplate:
    """Test get_coordination_template function."""

    def test_get_coordination_template_basic(self):
        """Should return template dict for valid metal."""
        from metal_chemistry import get_coordination_template

        template = get_coordination_template("ZN")
        assert template is not None
        assert template["metal"] == "ZN"
        assert "coordination_number" in template
        assert "geometries" in template
        assert "typical_residues" in template

    def test_get_coordination_template_with_cn(self):
        """Should return template for specific coordination number."""
        from metal_chemistry import get_coordination_template

        template = get_coordination_template("TB", coordination_number=9)
        assert template["metal"] == "TB"
        assert template["coordination_number"] == 9
        # Should match lanthanide_site motif
        assert "lanthanide_site" in template["motifs"]
        # Should have appropriate geometries
        geoms = template["geometries"]
        assert "square_antiprism" in geoms or "tricapped_trigonal_prism" in geoms

    def test_get_coordination_template_zinc_finger(self):
        """Should return zinc_finger motif for ZN with CN=4."""
        from metal_chemistry import get_coordination_template

        template = get_coordination_template("ZN", coordination_number=4)
        assert template["coordination_number"] == 4
        # Should include motifs that match
        assert any(m in template["motifs"] for m in ["zinc_finger", "carbonic_anhydrase"])
        assert "tetrahedral" in template["geometries"]

    def test_get_coordination_template_invalid_metal(self):
        """Should raise ValueError for invalid metal."""
        from metal_chemistry import get_coordination_template

        with pytest.raises(ValueError, match="Unknown metal"):
            get_coordination_template("XX")


class TestGetIdealGeometry:
    """Test get_ideal_geometry function."""

    def test_get_ideal_geometry_octahedral(self):
        """Should return 6 positions for octahedral."""
        from metal_chemistry import get_ideal_geometry

        positions = get_ideal_geometry("octahedral")
        assert len(positions) == 6

        # All positions should be on unit sphere
        for x, y, z in positions:
            dist = math.sqrt(x*x + y*y + z*z)
            assert abs(dist - 1.0) < 0.01

    def test_get_ideal_geometry_tetrahedral(self):
        """Should return 4 positions for tetrahedral."""
        from metal_chemistry import get_ideal_geometry

        positions = get_ideal_geometry("tetrahedral")
        assert len(positions) == 4

        # All positions should be on unit sphere
        for x, y, z in positions:
            dist = math.sqrt(x*x + y*y + z*z)
            assert abs(dist - 1.0) < 0.01

    def test_get_ideal_geometry_case_insensitive(self):
        """Should handle case-insensitive geometry names."""
        from metal_chemistry import get_ideal_geometry

        positions1 = get_ideal_geometry("octahedral")
        positions2 = get_ideal_geometry("OCTAHEDRAL")
        positions3 = get_ideal_geometry("Octahedral")
        assert positions1 == positions2 == positions3

    def test_get_ideal_geometry_invalid(self):
        """Should raise ValueError for invalid geometry."""
        from metal_chemistry import get_ideal_geometry

        with pytest.raises(ValueError, match="Unknown geometry"):
            get_ideal_geometry("invalid_geometry")


class TestValidateCoordination:
    """Test validate_coordination function."""

    def test_validate_coordination_good(self):
        """Valid octahedral ZN coordination should pass."""
        from metal_chemistry import validate_coordination

        # Metal at origin, ligands at 2.0 A in octahedral arrangement
        metal_pos = (0.0, 0.0, 0.0)
        dist = 2.1  # Within ZN-N range (1.95-2.20)
        positions = [
            (dist, 0.0, 0.0),
            (-dist, 0.0, 0.0),
            (0.0, dist, 0.0),
            (0.0, -dist, 0.0),
            (0.0, 0.0, dist),
            (0.0, 0.0, -dist),
        ]

        result = validate_coordination(positions, metal_pos, "ZN")
        assert result["valid"] is True
        assert result["coordination_number"] == 6
        assert result["geometry_match"] == "octahedral"
        assert len(result["issues"]) == 0

    def test_validate_coordination_bad_distances(self):
        """Positions too close should fail with distance issue."""
        from metal_chemistry import validate_coordination

        metal_pos = (0.0, 0.0, 0.0)
        # Distances way too short (0.5 A)
        dist = 0.5
        positions = [
            (dist, 0.0, 0.0),
            (-dist, 0.0, 0.0),
            (0.0, dist, 0.0),
            (0.0, -dist, 0.0),
        ]

        result = validate_coordination(positions, metal_pos, "ZN")
        assert result["valid"] is False
        assert result["coordination_number"] == 4
        # Should have distance issues
        assert any("distance" in issue for issue in result["issues"])

    def test_validate_coordination_far_distances(self):
        """Positions too far should fail with distance issue."""
        from metal_chemistry import validate_coordination

        metal_pos = (0.0, 0.0, 0.0)
        # Distances too long (5.0 A)
        dist = 5.0
        positions = [
            (dist, 0.0, 0.0),
            (-dist, 0.0, 0.0),
            (0.0, dist, 0.0),
            (0.0, -dist, 0.0),
        ]

        result = validate_coordination(positions, metal_pos, "ZN")
        assert result["valid"] is False
        assert any("distance" in issue for issue in result["issues"])

    def test_validate_coordination_returns_distances(self):
        """Should return calculated distances."""
        from metal_chemistry import validate_coordination

        metal_pos = (0.0, 0.0, 0.0)
        positions = [
            (2.0, 0.0, 0.0),
            (0.0, 2.5, 0.0),
        ]

        result = validate_coordination(positions, metal_pos, "ZN")
        assert "distances" in result
        assert len(result["distances"]) == 2
        assert abs(result["distances"][0] - 2.0) < 0.01
        assert abs(result["distances"][1] - 2.5) < 0.01

    def test_validate_coordination_expected_cn(self):
        """Should flag mismatched expected CN."""
        from metal_chemistry import validate_coordination

        metal_pos = (0.0, 0.0, 0.0)
        dist = 2.1
        positions = [
            (dist, 0.0, 0.0),
            (-dist, 0.0, 0.0),
            (0.0, dist, 0.0),
            (0.0, -dist, 0.0),
        ]

        # Expect 6 but provide 4
        result = validate_coordination(positions, metal_pos, "ZN", expected_cn=6)
        assert result["valid"] is False
        assert any("expected 6" in issue for issue in result["issues"])

    def test_validate_coordination_unknown_metal(self):
        """Should handle unknown metal gracefully."""
        from metal_chemistry import validate_coordination

        metal_pos = (0.0, 0.0, 0.0)
        positions = [(2.0, 0.0, 0.0)]

        result = validate_coordination(positions, metal_pos, "XX")
        assert result["valid"] is False
        assert any("Unknown metal" in issue for issue in result["issues"])


class TestGeometryPositionsOnUnitSphere:
    """Test that geometry positions are properly normalized to unit sphere."""

    def test_all_geometries_on_unit_sphere(self):
        """All geometry positions should be on unit sphere."""
        from metal_chemistry import IDEAL_GEOMETRIES

        for geom_name, geom_data in IDEAL_GEOMETRIES.items():
            for i, (x, y, z) in enumerate(geom_data["positions"]):
                dist = math.sqrt(x*x + y*y + z*z)
                assert abs(dist - 1.0) < 0.05, (
                    f"{geom_name} position {i} ({x:.3f}, {y:.3f}, {z:.3f}) "
                    f"has distance {dist:.3f}, expected 1.0"
                )

    def test_geometry_cn_matches_positions(self):
        """Coordination number should match position count."""
        from metal_chemistry import IDEAL_GEOMETRIES

        for geom_name, geom_data in IDEAL_GEOMETRIES.items():
            expected_cn = geom_data["coordination_number"]
            actual_cn = len(geom_data["positions"])
            assert expected_cn == actual_cn, (
                f"{geom_name}: CN is {expected_cn} but has {actual_cn} positions"
            )
