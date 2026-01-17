# test_metal_dimer_pipeline.py
"""
Integration Tests for Metal-Ligand Dimer Design Pipeline

Tests the complete pipeline from template retrieval through amino acid bias
generation and HSAB validation for metal-mediated dimer design.

Test Classes:
- TestTemplateRetrievalPipeline: Template retrieval with fallback
- TestBiasGenerationPipeline: Amino acid bias generation for different metals
- TestLigandAnalysisPipeline: Ligand donor analysis (requires RDKit)
- TestValidationPipeline: HSAB validation and coordination number checks
- TestEndToEndDesign: Complete workflow tests for lanthanide and zinc
"""

import pytest
from typing import Dict, Any

# Import modules under test
from metal_chemistry import (
    get_amino_acid_bias,
    get_coordination_number_range,
    validate_coordination_chemistry,
    get_hsab_class,
    is_lanthanide,
    METAL_DATABASE,
)
from design_types import METAL_PRESETS, get_recommended_config, DesignType
from metal_ligand_templates import (
    get_template,
    get_template_with_fallback,
    list_templates,
    get_template_info,
)
from metal_site_fetcher import get_reference_template, REFERENCE_STRUCTURES

# Check RDKit availability for conditional tests
try:
    from ligand_donors import (
        identify_donors_from_smiles,
        score_ligand_metal_compatibility,
        HAS_RDKIT,
    )
    RDKIT_AVAILABLE = HAS_RDKIT
except ImportError:
    RDKIT_AVAILABLE = False


# =============================================================================
# TEST: TEMPLATE RETRIEVAL PIPELINE
# =============================================================================

class TestTemplateRetrievalPipeline:
    """Test template retrieval with fallback mechanisms."""

    def test_library_template_available(self):
        """Built-in templates should be available without network."""
        templates = list_templates()
        assert len(templates) > 0
        assert "pqq_ca" in templates
        assert "citrate_tb" in templates

    def test_get_template_returns_dict(self):
        """get_template should return a dictionary with required keys."""
        template = get_template("pqq_ca")
        assert template is not None
        assert isinstance(template, dict)
        assert "name" in template
        assert "coordination" in template

    def test_citrate_lanthanide_templates_exist(self):
        """Citrate templates for lanthanides should exist."""
        for metal in ["tb", "eu", "gd"]:
            template = get_template(f"citrate_{metal}")
            assert template is not None, f"Missing template: citrate_{metal}"

    def test_template_with_fallback_library_source(self):
        """get_template_with_fallback should use library when available."""
        template = get_template_with_fallback("pqq_ca")
        assert template is not None
        # Should have source metadata
        assert "source" in template
        assert template["source"] in ["library", "pdb", "calculated"]

    def test_template_fallback_for_unknown_template(self):
        """get_template_with_fallback should fallback for unknown templates."""
        # Request an unknown template with a known metal
        template = get_template_with_fallback(
            "unknown_template",
            metal="ZN",
            fallback_enabled=True,
        )
        # Should return calculated fallback
        if template:
            assert template["source"] == "calculated"
            assert "warning" in template

    def test_template_info_returns_coordination_info(self):
        """get_template_info should return coordination details."""
        info = get_template_info("citrate_tb")
        assert info is not None
        assert "coordination_number" in info
        assert "protein_sites" in info
        assert "geometry" in info

    def test_reference_structures_defined(self):
        """Reference structures should be defined for common metals."""
        assert "ZN" in REFERENCE_STRUCTURES
        assert "TB" in REFERENCE_STRUCTURES
        assert "FE" in REFERENCE_STRUCTURES
        assert "CA" in REFERENCE_STRUCTURES


# =============================================================================
# TEST: BIAS GENERATION PIPELINE
# =============================================================================

class TestBiasGenerationPipeline:
    """Test amino acid bias generation respects HSAB chemistry."""

    def test_lanthanide_excludes_cysteine(self):
        """Lanthanide bias should exclude cysteine (Cys bias <= -3.0)."""
        for metal in ["TB", "EU", "GD", "LA"]:
            bias = get_amino_acid_bias(metal, 3)

            # Parse the bias string to get Cys weight
            cys_weight = None
            for part in bias.split(","):
                if part.startswith("C:"):
                    cys_weight = float(part.split(":")[1])
                    break

            assert cys_weight is not None, f"Cys not found in bias for {metal}"
            assert cys_weight <= -3.0, (
                f"Lanthanide {metal} should exclude Cys (weight <= -3.0), "
                f"got {cys_weight}"
            )

    def test_zinc_prefers_cysteine(self):
        """Zinc bias should prefer cysteine (Cys bias > 0)."""
        bias = get_amino_acid_bias("ZN", 2)

        cys_weight = None
        for part in bias.split(","):
            if part.startswith("C:"):
                cys_weight = float(part.split(":")[1])
                break

        assert cys_weight is not None, "Cys not found in zinc bias"
        assert cys_weight > 0, f"Zinc should prefer Cys (weight > 0), got {cys_weight}"

    def test_zinc_prefers_histidine(self):
        """Zinc bias should prefer histidine (His bias > 0)."""
        bias = get_amino_acid_bias("ZN", 2)

        his_weight = None
        for part in bias.split(","):
            if part.startswith("H:"):
                his_weight = float(part.split(":")[1])
                break

        assert his_weight is not None, "His not found in zinc bias"
        assert his_weight > 0, f"Zinc should prefer His (weight > 0), got {his_weight}"

    def test_lanthanide_favors_glu_asp(self):
        """Lanthanide bias should favor Glu and Asp (O donors)."""
        bias = get_amino_acid_bias("TB", 3)

        glu_weight = None
        asp_weight = None
        for part in bias.split(","):
            if part.startswith("E:"):
                glu_weight = float(part.split(":")[1])
            elif part.startswith("D:"):
                asp_weight = float(part.split(":")[1])

        assert glu_weight is not None, "Glu not found in lanthanide bias"
        assert asp_weight is not None, "Asp not found in lanthanide bias"
        assert glu_weight > 0, f"Lanthanide should favor Glu, got {glu_weight}"
        assert asp_weight > 0, f"Lanthanide should favor Asp, got {asp_weight}"

    def test_metal_presets_have_correct_bias(self):
        """METAL_PRESETS should have bias_AA matching HSAB chemistry."""
        # Lanthanide preset should exclude Cys
        lanthanide = METAL_PRESETS["lanthanide"]
        assert "excluded_residues" in lanthanide
        assert "C" in lanthanide["excluded_residues"]

        # Zinc preset should prefer Cys
        zinc = METAL_PRESETS["zinc"]
        assert "C" in zinc["preferred_residues"]

    def test_bias_format_is_valid_ligandmpnn(self):
        """Bias strings should follow LigandMPNN format: AA:weight,AA:weight."""
        for metal in ["ZN", "TB", "FE", "CA"]:
            ox_state = METAL_DATABASE[metal]["default_oxidation"]
            bias = get_amino_acid_bias(metal, ox_state)

            for part in bias.split(","):
                assert ":" in part, f"Invalid format in {metal} bias: {part}"
                aa, weight = part.split(":")
                assert len(aa) == 1, f"AA code should be single letter: {aa}"
                assert aa.isupper(), f"AA code should be uppercase: {aa}"
                # Should parse as float without error
                float(weight)


# =============================================================================
# TEST: LIGAND ANALYSIS PIPELINE (Requires RDKit)
# =============================================================================

class TestLigandAnalysisPipeline:
    """Test ligand donor analysis with HSAB scoring (requires RDKit)."""

    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not available")
    def test_identify_donors_from_citrate(self):
        """Citrate SMILES should identify carboxylate O donors."""
        # Citrate trianion SMILES
        citrate_smiles = "OC(CC(=O)[O-])(CC(=O)[O-])C(=O)[O-]"
        donors = identify_donors_from_smiles(citrate_smiles)

        assert len(donors) > 0, "No donors identified from citrate"

        # Should have O donors
        o_donors = [d for d in donors if d["element"] == "O"]
        assert len(o_donors) > 0, "No O donors found in citrate"

        # Should have carboxylate donors
        carboxylate_donors = [d for d in donors if d["type"] == "carboxylate"]
        assert len(carboxylate_donors) > 0, "No carboxylate donors found"

        # All citrate donors should be hard HSAB
        for donor in o_donors:
            assert donor["hsab"] == "hard", f"O donor should be hard: {donor}"

    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not available")
    def test_identify_donors_from_pqq(self):
        """PQQ SMILES should identify O and N donors."""
        # PQQ SMILES (pyrroloquinoline quinone)
        pqq_smiles = "OC(=O)c1[nH]c2c(c1)C(=O)C(=O)c3nc(cc(C(O)=O)c23)C(O)=O"
        donors = identify_donors_from_smiles(pqq_smiles)

        assert len(donors) > 0, "No donors identified from PQQ"

        # Should have both O and N donors
        elements = set(d["element"] for d in donors)
        assert "O" in elements, "PQQ should have O donors"
        assert "N" in elements, "PQQ should have N donors"

    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not available")
    def test_score_citrate_terbium_compatibility(self):
        """Citrate should be highly compatible with terbium (hard-hard)."""
        citrate_smiles = "OC(CC(=O)[O-])(CC(=O)[O-])C(=O)[O-]"
        score = score_ligand_metal_compatibility(citrate_smiles, "TB")

        # Hard ligand with hard metal should score high
        assert score >= 0.7, f"Citrate-Tb should be highly compatible, got {score}"

    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not available")
    def test_score_thiol_lanthanide_incompatibility(self):
        """Thiol ligand should be incompatible with lanthanide (soft-hard)."""
        thiol_smiles = "CCS"  # Ethyl thiol
        score = score_ligand_metal_compatibility(thiol_smiles, "TB")

        # Soft ligand with hard metal should score low
        assert score <= 0.3, f"Thiol-Tb should be incompatible, got {score}"

    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not available")
    def test_score_thiol_zinc_compatibility(self):
        """Thiol ligand should be compatible with zinc (soft-borderline)."""
        thiol_smiles = "CCS"  # Ethyl thiol
        score = score_ligand_metal_compatibility(thiol_smiles, "ZN")

        # Soft ligand with borderline metal should be reasonably compatible
        assert score >= 0.5, f"Thiol-Zn should be compatible, got {score}"


# =============================================================================
# TEST: VALIDATION PIPELINE
# =============================================================================

class TestValidationPipeline:
    """Test HSAB validation and coordination number validation."""

    def test_hsab_validation_lanthanide_with_cys(self):
        """Lanthanide with Cys coordination should fail HSAB validation."""
        residues = ["GLU", "ASP", "CYS", "GLU", "ASP", "ASN", "GLN", "GLU"]
        result = validate_coordination_chemistry("TB", residues, 3)

        assert result["hsab_compatible"] == False, (
            "Lanthanide with Cys should fail HSAB validation"
        )
        assert "CYS" in result.get("incompatible_residues", [])

    def test_hsab_validation_zinc_with_his_cys(self):
        """Zinc with His/Cys coordination should pass HSAB validation."""
        residues = ["HIS", "HIS", "CYS", "CYS"]
        result = validate_coordination_chemistry("ZN", residues, 2)

        assert result["valid"] == True
        assert result["hsab_compatible"] == True

    def test_coordination_number_lanthanide(self):
        """Lanthanide coordination number should be 8-9."""
        min_cn, max_cn = get_coordination_number_range("TB", 3)

        assert min_cn >= 7, f"Lanthanide min CN should be >= 7, got {min_cn}"
        assert 8 in range(min_cn, max_cn + 1), "CN 8 should be valid for lanthanide"
        assert 9 in range(min_cn, max_cn + 1), "CN 9 should be valid for lanthanide"

    def test_coordination_number_zinc(self):
        """Zinc coordination number should include 4."""
        min_cn, max_cn = get_coordination_number_range("ZN", 2)

        assert 4 in range(min_cn, max_cn + 1), "CN 4 should be valid for zinc"

    def test_validate_too_few_residues(self):
        """Too few coordinating residues should produce warnings."""
        residues = ["HIS"]  # Only 1 residue for zinc
        result = validate_coordination_chemistry("ZN", residues, 2)

        # Should either be invalid or have warnings
        assert result["valid"] == False or len(result.get("warnings", [])) > 0

    def test_validate_lanthanide_with_only_oxygen(self):
        """Lanthanide with only O donors should pass validation."""
        residues = ["GLU", "GLU", "ASP", "ASP", "ASN", "GLN", "GLU", "ASP"]
        result = validate_coordination_chemistry("TB", residues, 3)

        assert result["valid"] == True
        assert result["hsab_compatible"] == True
        assert result["coordination_number"] == 8

    def test_is_lanthanide_function(self):
        """is_lanthanide should correctly identify lanthanides."""
        lanthanides = ["TB", "EU", "GD", "LA", "CE", "SM", "YB"]
        non_lanthanides = ["ZN", "FE", "CA", "MG", "CU"]

        for metal in lanthanides:
            assert is_lanthanide(metal) == True, f"{metal} should be lanthanide"

        for metal in non_lanthanides:
            assert is_lanthanide(metal) == False, f"{metal} should not be lanthanide"


# =============================================================================
# TEST: END-TO-END DESIGN WORKFLOWS
# =============================================================================

class TestEndToEndDesign:
    """End-to-end workflow tests for complete design pipeline."""

    def test_lanthanide_dimer_design_workflow(self):
        """Test complete workflow for lanthanide-mediated dimer design."""
        # Step 1: Get template
        template = get_template("citrate_tb")
        assert template is not None

        # Step 2: Get coordination info
        coord_info = template.get("coordination", {})
        expected_cn = coord_info.get("metal_coordination_number", 9)
        assert expected_cn >= 8

        # Step 3: Get amino acid bias
        bias = get_amino_acid_bias("TB", 3)

        # Verify Cys is excluded (weight <= -3.0)
        cys_weight = None
        for part in bias.split(","):
            if part.startswith("C:"):
                cys_weight = float(part.split(":")[1])
                break
        assert cys_weight is not None
        assert cys_weight <= -3.0

        # Step 4: Validate example coordination
        residues = ["GLU", "GLU", "ASP", "ASP", "ASN", "GLN", "GLU", "ASP"]
        result = validate_coordination_chemistry("TB", residues, 3)
        assert result["hsab_compatible"] == True

        # Step 5: Get design config
        config = get_recommended_config(
            DesignType.METAL_MEDIATED_DIMER,
            metal_type="lanthanide"
        )
        assert config.sequence_tool == "ligand_mpnn"
        assert config.use_atom_context == True

    def test_zinc_dimer_design_workflow(self):
        """Test complete workflow for zinc-mediated dimer design."""
        # Step 1: Get coordination number range
        min_cn, max_cn = get_coordination_number_range("ZN", 2)
        assert 4 in range(min_cn, max_cn + 1)

        # Step 2: Get amino acid bias
        bias = get_amino_acid_bias("ZN", 2)

        # Verify Cys and His are preferred (weight > 0)
        cys_weight = None
        his_weight = None
        for part in bias.split(","):
            if part.startswith("C:"):
                cys_weight = float(part.split(":")[1])
            elif part.startswith("H:"):
                his_weight = float(part.split(":")[1])

        assert cys_weight is not None and cys_weight > 0
        assert his_weight is not None and his_weight > 0

        # Step 3: Validate example coordination
        residues = ["HIS", "HIS", "CYS", "CYS"]
        result = validate_coordination_chemistry("ZN", residues, 2)
        assert result["valid"] == True
        assert result["hsab_compatible"] == True

        # Step 4: Get design config
        config = get_recommended_config(
            DesignType.METAL_MEDIATED_DIMER,
            metal_type="zinc"
        )
        assert config.sequence_tool == "ligand_mpnn"
        assert config.symmetry_type == "C2"

    def test_hsab_class_consistency(self):
        """HSAB class should be consistent across pipeline."""
        # Lanthanides should be hard
        for metal in ["TB", "EU", "GD", "LA"]:
            hsab = get_hsab_class(metal, 3)
            assert hsab == "hard", f"{metal} should be hard acid"

        # Zinc should be borderline
        assert get_hsab_class("ZN", 2) == "borderline"

        # Calcium should be hard
        assert get_hsab_class("CA", 2) == "hard"

    def test_template_fallback_chain(self):
        """Template retrieval should work through fallback chain."""
        # Known template should work directly
        template1 = get_template_with_fallback("citrate_tb")
        assert template1 is not None
        assert template1["source"] in ["library", "pdb"]

        # Unknown template with known metal should fallback
        template2 = get_template_with_fallback(
            "unknown_complex",
            metal="TB",
            fallback_enabled=True,
        )
        if template2:
            assert template2["source"] == "calculated"

    def test_metal_preset_matches_chemistry(self):
        """METAL_PRESETS should match underlying chemistry."""
        # Lanthanide preset
        lanthanide = METAL_PRESETS["lanthanide"]
        tb_bias = get_amino_acid_bias("TB", 3)

        # Both should exclude Cys
        assert "C" in lanthanide.get("excluded_residues", [])
        assert "C:-" in tb_bias

        # Zinc preset
        zinc = METAL_PRESETS["zinc"]
        zn_bias = get_amino_acid_bias("ZN", 2)

        # Both should prefer Cys
        assert "C" in zinc.get("preferred_residues", [])
        cys_in_bias = any(
            part.startswith("C:") and float(part.split(":")[1]) > 0
            for part in zn_bias.split(",")
        )
        assert cys_in_bias

    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not available")
    def test_ligand_metal_compatibility_integration(self):
        """Test ligand-metal compatibility across pipeline."""
        # Citrate should be compatible with lanthanides
        citrate = "OC(CC(=O)[O-])(CC(=O)[O-])C(=O)[O-]"
        tb_score = score_ligand_metal_compatibility(citrate, "TB")

        # Should be highly compatible (hard-hard)
        assert tb_score >= 0.7

        # Verify this matches with validation
        result = validate_coordination_chemistry(
            "TB",
            ["GLU", "GLU", "ASP", "ASP", "ASN", "GLN", "GLU", "ASP"],
            3
        )
        assert result["hsab_compatible"] == True


# =============================================================================
# TEST: DESIGN CONFIG INTEGRATION
# =============================================================================

class TestDesignConfigIntegration:
    """Test design configuration matches metal chemistry."""

    def test_metal_dimer_config_uses_ligandmpnn(self):
        """Metal dimer design should use LigandMPNN for atom context."""
        config = get_recommended_config(DesignType.METAL_MEDIATED_DIMER)
        assert config.sequence_tool == "ligand_mpnn"
        assert config.use_atom_context == True

    def test_metal_dimer_config_is_symmetric(self):
        """Metal dimer design should be configured for symmetry."""
        config = get_recommended_config(DesignType.METAL_MEDIATED_DIMER)
        assert config.use_symmetry == True
        assert config.symmetry_type == "C2"

    def test_metal_dimer_config_with_lanthanide_preset(self):
        """Lanthanide preset should override bias_AA appropriately."""
        config = get_recommended_config(
            DesignType.METAL_MEDIATED_DIMER,
            metal_type="lanthanide"
        )

        # Should have a bias_AA string
        assert config.bias_AA is not None

        # Should have Cys penalized
        if "C:" in config.bias_AA:
            for part in config.bias_AA.split(","):
                if part.startswith("C:"):
                    weight = float(part.split(":")[1])
                    assert weight < 0, "Lanthanide config should penalize Cys"

    def test_metal_dimer_config_with_zinc_preset(self):
        """Zinc preset should favor His/Cys coordination."""
        config = get_recommended_config(
            DesignType.METAL_MEDIATED_DIMER,
            metal_type="zinc"
        )

        # Should have a bias_AA string
        assert config.bias_AA is not None

        # Should favor Cys
        cys_positive = False
        for part in config.bias_AA.split(","):
            if part.startswith("C:"):
                weight = float(part.split(":")[1])
                if weight > 0:
                    cys_positive = True

        assert cys_positive, "Zinc config should favor Cys"
