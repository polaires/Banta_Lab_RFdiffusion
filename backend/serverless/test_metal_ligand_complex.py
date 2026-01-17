#!/usr/bin/env python3
"""
Test cases for metal-ligand complex dimer design.

Test complexes:
1. PQQ-Ca (pyrroloquinoline quinone with calcium)
2. Citrate-Tb (citrate with terbium)

Run with: pytest test_metal_ligand_complex.py -v
"""

import os
import sys
import pytest
from typing import Dict, Any

# Setup path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# PQQ-Ca test complex (minimal PDB)
# PQQ coordinating Ca via O5, N6, O7
PQQ_CA_TEST_PDB = """HETATM    1  C1  PQQ L   1      48.000  50.000  48.000  1.00  0.00           C
HETATM    2  C2  PQQ L   1      49.200  50.000  48.500  1.00  0.00           C
HETATM    3  O5  PQQ L   1      50.500  51.000  50.000  1.00  0.00           O
HETATM    4  N6  PQQ L   1      49.500  50.500  50.500  1.00  0.00           N
HETATM    5  O7  PQQ L   1      50.000  49.000  50.000  1.00  0.00           O
HETATM    6  C3  PQQ L   1      52.000  50.000  48.000  1.00  0.00           C
HETATM    7  O8  PQQ L   1      53.000  50.500  48.000  1.00  0.00           O
HETATM   10 CA   CA  M   2      50.000  50.000  52.400  1.00  0.00          CA
END
"""

# Citrate-Tb test complex
# Citrate coordinating Tb via 4 carboxylate oxygens
CITRATE_TB_TEST_PDB = """HETATM    1  C1  CIT L   1      50.000  50.000  48.000  1.00  0.00           C
HETATM    2  O1  CIT L   1      50.500  51.000  49.500  1.00  0.00           O
HETATM    3  O2  CIT L   1      49.500  49.000  49.500  1.00  0.00           O
HETATM    4  C2  CIT L   1      50.000  50.000  46.500  1.00  0.00           C
HETATM    5  O3  CIT L   1      51.200  50.500  49.200  1.00  0.00           O
HETATM    6  O4  CIT L   1      48.800  49.500  49.200  1.00  0.00           O
HETATM    7  C3  CIT L   1      50.000  51.500  46.000  1.00  0.00           C
HETATM    8  O5  CIT L   1      50.500  52.000  46.500  1.00  0.00           O
HETATM    9  O6  CIT L   1      49.500  51.500  45.000  1.00  0.00           O
HETATM   10 TB   TB  M   2      50.000  50.000  50.350  1.00  0.00          TB
END
"""

# Simple test complex without coordinating atoms for negative testing
SIMPLE_TEST_PDB = """HETATM    1  C1  LIG L   1      50.000  50.000  48.000  1.00  0.00           C
HETATM    2  C2  LIG L   1      51.000  50.000  48.000  1.00  0.00           C
HETATM   10 CA   CA  M   2      50.000  50.000  55.000  1.00  0.00          CA
END
"""


class TestComplexParser:
    """Test metal-ligand complex parsing from inference_utils."""

    def test_parse_pqq_ca_complex(self):
        """Should correctly parse PQQ-Ca complex."""
        from inference_utils import parse_metal_ligand_complex

        complexes = parse_metal_ligand_complex(PQQ_CA_TEST_PDB)

        assert len(complexes) == 1
        c = complexes[0]
        assert c.metal_code == "CA"
        assert c.ligand_res_name == "PQQ"
        assert c.metal_coords is not None
        # Should detect O5, N6, O7 coordinating Ca
        assert len(c.coordination_bonds) >= 2, f"Expected at least 2 bonds, got {len(c.coordination_bonds)}"

    def test_parse_citrate_tb_complex(self):
        """Should correctly parse Citrate-Tb complex."""
        from inference_utils import parse_metal_ligand_complex

        complexes = parse_metal_ligand_complex(CITRATE_TB_TEST_PDB)

        assert len(complexes) == 1
        c = complexes[0]
        assert c.metal_code == "TB"
        assert c.ligand_res_name == "CIT"
        # Citrate should have multiple O donors
        assert len(c.coordination_bonds) >= 2

    def test_parse_no_metal(self):
        """Should return empty list for PDB without metal."""
        from inference_utils import parse_metal_ligand_complex

        pdb_no_metal = """HETATM    1  C1  LIG L   1      50.000  50.000  48.000  1.00  0.00           C
END
"""
        complexes = parse_metal_ligand_complex(pdb_no_metal)
        assert len(complexes) == 0

    def test_parse_metal_too_far(self):
        """Should not detect coordination when metal is too far from ligand."""
        from inference_utils import parse_metal_ligand_complex

        complexes = parse_metal_ligand_complex(SIMPLE_TEST_PDB, coordination_cutoff=3.0)

        # Metal at 55Å z vs ligand at 48Å z = 7Å apart, beyond cutoff
        if len(complexes) > 0:
            # If complex is found, it should have no coordination bonds
            assert len(complexes[0].coordination_bonds) == 0

    def test_detect_available_sites_octahedral(self):
        """Should detect available coordination sites for Ca complex."""
        from inference_utils import parse_metal_ligand_complex, detect_available_coordination_sites

        complexes = parse_metal_ligand_complex(PQQ_CA_TEST_PDB)
        if not complexes:
            pytest.skip("No complex found in test PDB")

        sites = detect_available_coordination_sites(complexes[0], target_coordination=6)

        # Should have some available sites for protein coordination
        assert len(sites) >= 1, "Expected at least 1 available coordination site"
        assert all("coords" in s for s in sites)
        assert all("preferred_donors" in s for s in sites)


class TestTemplateLibrary:
    """Test metal-ligand complex templates."""

    def test_pqq_ca_template_exists(self):
        """PQQ-Ca template should exist."""
        from metal_ligand_templates import METAL_LIGAND_COMPLEX_TEMPLATES

        assert "pqq_ca" in METAL_LIGAND_COMPLEX_TEMPLATES
        template = METAL_LIGAND_COMPLEX_TEMPLATES["pqq_ca"]
        assert template["metal"] == "CA"
        assert "coordination" in template
        assert template["ligand_res_name"] == "PQQ"

    def test_citrate_templates_exist(self):
        """Citrate-lanthanide templates should exist."""
        from metal_ligand_templates import METAL_LIGAND_COMPLEX_TEMPLATES

        assert "citrate_tb" in METAL_LIGAND_COMPLEX_TEMPLATES
        assert "citrate_eu" in METAL_LIGAND_COMPLEX_TEMPLATES
        assert "citrate_gd" in METAL_LIGAND_COMPLEX_TEMPLATES

    def test_list_templates(self):
        """Should list all available templates."""
        from metal_ligand_templates import list_templates

        templates = list_templates()
        assert isinstance(templates, list)
        assert len(templates) >= 4  # pqq_ca, citrate_tb, citrate_eu, citrate_gd
        assert "pqq_ca" in templates

    def test_get_template(self):
        """Should retrieve template by name."""
        from metal_ligand_templates import get_template

        template = get_template("pqq_ca")
        assert template is not None
        assert template["metal"] == "CA"

        # Case insensitivity
        template2 = get_template("PQQ_CA")
        assert template2 is not None

        # Unknown template
        assert get_template("unknown_template") is None

    def test_get_template_info(self):
        """Should get summary info for template."""
        from metal_ligand_templates import get_template_info

        info = get_template_info("pqq_ca")
        assert info is not None
        assert "metal" in info
        assert "ligand" in info
        assert "coordination_number" in info
        assert "preferred_donors" in info

    def test_generate_pqq_ca_pdb(self):
        """Should generate valid PQQ-Ca PDB from template."""
        from metal_ligand_templates import generate_complex_pdb

        pdb = generate_complex_pdb("pqq_ca")
        assert pdb is not None
        assert "PQQ" in pdb or "CA" in pdb  # Either ligand or metal should be present
        assert "HETATM" in pdb
        assert "END" in pdb

    def test_generate_citrate_tb_pdb(self):
        """Should generate valid Citrate-Tb PDB from template."""
        from metal_ligand_templates import generate_complex_pdb

        pdb = generate_complex_pdb("citrate_tb")
        assert pdb is not None
        assert "HETATM" in pdb
        assert "TB" in pdb

    def test_get_coordination_info(self):
        """Should return coordination info for RFD3 config."""
        from metal_ligand_templates import get_template_coordination_info

        info = get_template_coordination_info("pqq_ca")
        assert "ligand_name" in info
        assert "metal_code" in info
        assert "target_coordination" in info
        assert "protein_sites" in info


class TestMetalValidation:
    """Test metal-ligand complex validation."""

    def test_validate_pqq_ca_complex(self):
        """Should validate PQQ-Ca complex structure."""
        from metal_validation import validate_metal_ligand_complex_site

        result = validate_metal_ligand_complex_site(
            PQQ_CA_TEST_PDB,
            metal="CA",
            ligand_name="PQQ",
            expected_ligand_donors=3,
            expected_protein_donors=0,  # No protein in test PDB
        )

        assert "coordination_number" in result
        assert "ligand_coordination" in result
        assert result["ligand_coordination"] >= 2  # O5, N6, O7 within cutoff

    def test_validate_citrate_tb_complex(self):
        """Should validate Citrate-Tb complex structure."""
        from metal_validation import validate_metal_ligand_complex_site

        result = validate_metal_ligand_complex_site(
            CITRATE_TB_TEST_PDB,
            metal="TB",
            ligand_name="CIT",
            expected_ligand_donors=4,
            expected_protein_donors=0,
        )

        assert "coordination_number" in result
        assert "ligand_coordination" in result
        assert result["ligand_coordination"] >= 2  # Carboxylate oxygens

    def test_validate_missing_metal(self):
        """Should report error for missing metal."""
        from metal_validation import validate_metal_ligand_complex_site

        pdb_no_metal = """HETATM    1  O1  PQQ L   1      50.000  50.000  50.000  1.00  0.00           O
END
"""
        result = validate_metal_ligand_complex_site(
            pdb_no_metal,
            metal="CA",
            ligand_name="PQQ",
        )

        assert result["success"] is False
        assert "error" in result or result["coordination_number"] == 0


class TestHandlerIntegration:
    """Test handler function integration."""

    def test_handler_task_registered(self):
        """Handler should accept interface_metal_ligand_design task."""
        # Import handler to check task dispatch
        import handler as h

        # Check that the function exists
        assert hasattr(h, "handle_interface_metal_ligand_design")

    def test_handler_requires_input(self):
        """Handler should require complex_pdb or template_name."""
        from handler import handle_interface_metal_ligand_design

        result = handle_interface_metal_ligand_design({})

        assert result["status"] == "failed"
        assert "error" in result

    def test_handler_unknown_template(self):
        """Handler should reject unknown template."""
        from handler import handle_interface_metal_ligand_design

        result = handle_interface_metal_ligand_design({
            "template_name": "nonexistent_template",
        })

        assert result["status"] == "failed"
        assert "Unknown template" in result["error"]

    def test_pqq_ca_design_mock(self):
        """Should generate mock design for PQQ-Ca."""
        from handler import handle_interface_metal_ligand_design

        result = handle_interface_metal_ligand_design({
            "template_name": "pqq_ca",
            "approach": "joint",
            "chain_length": "60-80",
            "num_designs": 1,
            "use_mock": True,
        })

        assert result["status"] == "completed", f"Failed: {result.get('error')}"
        assert "designs" in result["result"]
        assert len(result["result"]["designs"]) >= 1
        assert result["result"]["metal"] == "CA"

    def test_citrate_tb_design_mock(self):
        """Should generate mock design for Citrate-Tb."""
        from handler import handle_interface_metal_ligand_design

        result = handle_interface_metal_ligand_design({
            "template_name": "citrate_tb",
            "approach": "joint",
            "chain_length": "60-80",
            "num_designs": 1,
            "chain_a_donors": ["Glu", "Asp"],
            "chain_b_donors": ["Glu", "Asp", "Asn"],
            "use_mock": True,
        })

        assert result["status"] == "completed", f"Failed: {result.get('error')}"
        assert "designs" in result["result"]
        assert result["result"]["metal"] == "TB"

    def test_custom_complex_pdb(self):
        """Should accept custom complex PDB."""
        from handler import handle_interface_metal_ligand_design

        result = handle_interface_metal_ligand_design({
            "complex_pdb": PQQ_CA_TEST_PDB,
            "approach": "joint",
            "chain_length": "60-80",
            "num_designs": 1,
            "use_mock": True,
        })

        assert result["status"] == "completed", f"Failed: {result.get('error')}"
        assert "designs" in result["result"]


class TestMetalLigandComplexDataclass:
    """Test MetalLigandComplex dataclass methods."""

    def test_get_ligand_coordinating_atoms(self):
        """Should return list of coordinating atom names."""
        from inference_utils import parse_metal_ligand_complex

        complexes = parse_metal_ligand_complex(PQQ_CA_TEST_PDB)
        if not complexes:
            pytest.skip("No complex found")

        atoms = complexes[0].get_ligand_coordinating_atoms()
        assert isinstance(atoms, list)

    def test_get_complex_pdb(self):
        """Should generate PDB content from complex."""
        from inference_utils import parse_metal_ligand_complex

        complexes = parse_metal_ligand_complex(PQQ_CA_TEST_PDB)
        if not complexes:
            pytest.skip("No complex found")

        pdb = complexes[0].get_complex_pdb()
        assert "HETATM" in pdb
        assert "CA" in pdb
        assert "PQQ" in pdb


class TestDatabaseDrivenTemplates:
    """Test database-first template retrieval."""

    def test_get_template_prefers_database(self):
        """Template retrieval should try database first."""
        from metal_ligand_templates import get_template_with_fallback

        # This should try PDB first, then fall back to library
        template = get_template_with_fallback("pqq_ca")

        assert template is not None
        # Should have either source=pdb or source=library or source=calculated
        assert template.get("source") in ["pdb", "calculated", "library"]

    def test_fallback_to_calculated(self):
        """Should fall back to calculated if PDB unavailable."""
        from metal_ligand_templates import get_template_with_fallback

        # Non-existent complex - must fall back
        template = get_template_with_fallback(
            "unknown_complex",
            metal="TB",
            fallback_enabled=True,
        )

        # Should return None or a calculated template
        if template:
            assert template["source"] == "calculated"


# Standalone test runner
if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
