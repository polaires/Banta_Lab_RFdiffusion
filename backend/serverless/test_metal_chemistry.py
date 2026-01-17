# test_metal_chemistry.py
"""
Tests for metal chemistry database module.

Tests HSAB (Hard-Soft Acid-Base) theory compliance and metal coordination chemistry.
"""

import pytest
from metal_chemistry import (
    get_hsab_class,
    get_preferred_donors,
    get_bond_distance_range,
    get_amino_acid_bias,
    get_coordination_number_range,
    validate_coordination_chemistry,
    METAL_DATABASE,
)


# =============================================================================
# HSAB CLASSIFICATION TESTS
# =============================================================================

class TestHSABClassification:
    """Tests for HSAB (Hard-Soft Acid-Base) classification."""

    def test_zinc_is_borderline(self):
        """Zn2+ is a borderline acid."""
        assert get_hsab_class("ZN", 2) == "borderline"

    def test_iron2_is_borderline(self):
        """Fe2+ is a borderline acid."""
        assert get_hsab_class("FE", 2) == "borderline"

    def test_iron3_is_hard(self):
        """Fe3+ is a hard acid."""
        assert get_hsab_class("FE", 3) == "hard"

    def test_copper1_is_soft(self):
        """Cu1+ is a soft acid."""
        assert get_hsab_class("CU", 1) == "soft"

    def test_copper2_is_borderline(self):
        """Cu2+ is a borderline acid."""
        assert get_hsab_class("CU", 2) == "borderline"

    def test_calcium_is_hard(self):
        """Ca2+ is a hard acid."""
        assert get_hsab_class("CA", 2) == "hard"

    def test_magnesium_is_hard(self):
        """Mg2+ is a hard acid."""
        assert get_hsab_class("MG", 2) == "hard"

    def test_manganese_is_borderline(self):
        """Mn2+ is a borderline acid."""
        assert get_hsab_class("MN", 2) == "borderline"

    def test_cobalt_is_borderline(self):
        """Co2+ is a borderline acid."""
        assert get_hsab_class("CO", 2) == "borderline"

    def test_nickel_is_borderline(self):
        """Ni2+ is a borderline acid."""
        assert get_hsab_class("NI", 2) == "borderline"

    # Lanthanides - all hard acids
    def test_terbium_is_hard(self):
        """Tb3+ is a hard acid."""
        assert get_hsab_class("TB", 3) == "hard"

    def test_europium_is_hard(self):
        """Eu3+ is a hard acid."""
        assert get_hsab_class("EU", 3) == "hard"

    def test_gadolinium_is_hard(self):
        """Gd3+ is a hard acid."""
        assert get_hsab_class("GD", 3) == "hard"

    def test_lanthanum_is_hard(self):
        """La3+ is a hard acid."""
        assert get_hsab_class("LA", 3) == "hard"

    def test_cerium_is_hard(self):
        """Ce3+ is a hard acid."""
        assert get_hsab_class("CE", 3) == "hard"

    def test_samarium_is_hard(self):
        """Sm3+ is a hard acid."""
        assert get_hsab_class("SM", 3) == "hard"

    def test_ytterbium_is_hard(self):
        """Yb3+ is a hard acid."""
        assert get_hsab_class("YB", 3) == "hard"

    def test_unknown_metal_raises_error(self):
        """Unknown metal should raise ValueError."""
        with pytest.raises(ValueError, match="Unknown metal"):
            get_hsab_class("XX", 2)

    def test_invalid_oxidation_state_raises_error(self):
        """Invalid oxidation state should raise ValueError."""
        with pytest.raises(ValueError, match="oxidation state"):
            get_hsab_class("ZN", 5)

    def test_case_insensitive(self):
        """Metal symbols should be case insensitive."""
        assert get_hsab_class("zn", 2) == "borderline"
        assert get_hsab_class("Zn", 2) == "borderline"
        assert get_hsab_class("ZN", 2) == "borderline"


# =============================================================================
# PREFERRED DONORS TESTS
# =============================================================================

class TestPreferredDonors:
    """Tests for preferred donor atoms based on HSAB theory."""

    def test_zinc_prefers_sulfur_and_nitrogen(self):
        """Zn2+ (borderline) prefers S (Cys) and N (His)."""
        donors = get_preferred_donors("ZN", 2, "catalytic")
        assert "S" in donors
        assert "N" in donors
        assert donors["S"] > 0  # Positive weight
        assert donors["N"] > 0  # Positive weight

    def test_lanthanide_only_oxygen_donors(self):
        """Lanthanides (hard acids) ONLY accept O donors."""
        for metal in ["TB", "EU", "GD", "LA", "CE", "SM", "YB"]:
            donors = get_preferred_donors(metal, 3, "structural")
            assert "O" in donors
            assert donors["O"] > 0
            # S should be excluded or heavily penalized
            if "S" in donors:
                assert donors["S"] < 0  # Negative weight (penalized)

    def test_copper1_prefers_soft_donors(self):
        """Cu1+ (soft) prefers S and N donors."""
        donors = get_preferred_donors("CU", 1, "catalytic")
        assert "S" in donors
        assert donors["S"] > 0

    def test_calcium_prefers_oxygen(self):
        """Ca2+ (hard) prefers O donors."""
        donors = get_preferred_donors("CA", 2, "structural")
        assert "O" in donors
        assert donors["O"] > 0

    def test_iron3_prefers_hard_donors(self):
        """Fe3+ (hard) prefers O donors."""
        donors = get_preferred_donors("FE", 3, "catalytic")
        assert "O" in donors
        assert donors["O"] > 0


# =============================================================================
# AMINO ACID BIAS TESTS
# =============================================================================

class TestAminoAcidBias:
    """Tests for amino acid bias string generation."""

    def test_zinc_bias_includes_his_cys(self):
        """Zinc bias should favor H (His) and C (Cys)."""
        bias = get_amino_acid_bias("ZN", 2, "catalytic")
        assert "H:" in bias
        assert "C:" in bias
        # Should have positive weights
        assert "H:2" in bias or "H:1" in bias or "H:3" in bias

    def test_lanthanide_excludes_cys(self):
        """Lanthanide bias should heavily penalize C (Cys)."""
        for metal in ["TB", "EU", "GD"]:
            bias = get_amino_acid_bias(metal, 3, "structural")
            # Cys should be penalized with negative weight
            assert "C:-" in bias
            # Should favor E (Glu) and D (Asp)
            assert "E:" in bias
            assert "D:" in bias

    def test_lanthanide_cys_penalty_is_strong(self):
        """Lanthanide Cys penalty should be -5.0 or stronger."""
        bias = get_amino_acid_bias("TB", 3, "structural")
        # Parse the bias to check Cys weight
        parts = bias.split(",")
        for part in parts:
            if part.startswith("C:"):
                weight = float(part.split(":")[1])
                assert weight <= -5.0, f"Cys penalty {weight} should be <= -5.0"

    def test_calcium_omits_cys(self):
        """Calcium bias should penalize cysteine."""
        bias = get_amino_acid_bias("CA", 2, "structural")
        assert "C:-" in bias

    def test_iron2_allows_cys(self):
        """Fe2+ (borderline) allows Cys coordination."""
        bias = get_amino_acid_bias("FE", 2, "catalytic")
        # Should not heavily penalize Cys
        if "C:" in bias:
            parts = bias.split(",")
            for part in parts:
                if part.startswith("C:"):
                    weight = float(part.split(":")[1])
                    assert weight >= 0, f"Fe2+ should not penalize Cys"

    def test_bias_format_is_valid(self):
        """Bias string should match LigandMPNN format: AA:weight,AA:weight."""
        bias = get_amino_acid_bias("ZN", 2, "catalytic")
        parts = bias.split(",")
        for part in parts:
            assert ":" in part, f"Invalid format: {part}"
            aa, weight = part.split(":")
            assert len(aa) == 1, f"AA code should be single letter: {aa}"
            assert aa.isupper(), f"AA code should be uppercase: {aa}"
            float(weight)  # Should not raise


# =============================================================================
# BOND DISTANCE TESTS
# =============================================================================

class TestBondDistanceRange:
    """Tests for metal-donor bond distance ranges."""

    def test_zinc_nitrogen_distance(self):
        """Zn-N bond should be around 2.0-2.2 Angstroms."""
        min_dist, max_dist = get_bond_distance_range("ZN", "N", 2)
        assert 1.9 <= min_dist <= 2.1
        assert 2.1 <= max_dist <= 2.4

    def test_zinc_sulfur_distance(self):
        """Zn-S bond should be around 2.2-2.4 Angstroms."""
        min_dist, max_dist = get_bond_distance_range("ZN", "S", 2)
        assert 2.1 <= min_dist <= 2.3
        assert 2.3 <= max_dist <= 2.6

    def test_lanthanide_oxygen_distance(self):
        """Lanthanide-O bond should be around 2.3-2.6 Angstroms."""
        min_dist, max_dist = get_bond_distance_range("TB", "O", 3)
        assert 2.2 <= min_dist <= 2.5
        assert 2.4 <= max_dist <= 2.8

    def test_calcium_oxygen_distance(self):
        """Ca-O bond should be around 2.3-2.6 Angstroms."""
        min_dist, max_dist = get_bond_distance_range("CA", "O", 2)
        assert 2.2 <= min_dist <= 2.5
        assert 2.4 <= max_dist <= 2.8

    def test_copper_sulfur_distance(self):
        """Cu-S bond should be around 2.1-2.4 Angstroms."""
        min_dist, max_dist = get_bond_distance_range("CU", "S", 1)
        assert 2.0 <= min_dist <= 2.3
        assert 2.2 <= max_dist <= 2.5


# =============================================================================
# COORDINATION NUMBER TESTS
# =============================================================================

class TestCoordinationNumberRange:
    """Tests for coordination number ranges."""

    def test_zinc_coordination_number(self):
        """Zn2+ typically has CN 4-6."""
        min_cn, max_cn = get_coordination_number_range("ZN", 2)
        assert min_cn >= 3
        assert max_cn <= 7
        assert 4 in range(min_cn, max_cn + 1)

    def test_lanthanide_coordination_number(self):
        """Lanthanides typically have CN 8-9."""
        min_cn, max_cn = get_coordination_number_range("TB", 3)
        assert min_cn >= 6
        assert max_cn >= 8
        assert 8 in range(min_cn, max_cn + 1)
        assert 9 in range(min_cn, max_cn + 1)

    def test_calcium_coordination_number(self):
        """Ca2+ typically has CN 6-8."""
        min_cn, max_cn = get_coordination_number_range("CA", 2)
        assert 6 in range(min_cn, max_cn + 1)
        assert 7 in range(min_cn, max_cn + 1)

    def test_copper_coordination_number(self):
        """Cu2+ typically has CN 4-6."""
        min_cn, max_cn = get_coordination_number_range("CU", 2)
        assert 4 in range(min_cn, max_cn + 1)


# =============================================================================
# COORDINATION VALIDATION TESTS
# =============================================================================

class TestCoordinationValidation:
    """Tests for coordination chemistry validation."""

    def test_valid_zinc_coordination(self):
        """Valid zinc coordination with His and Cys."""
        residues = ["HIS", "HIS", "CYS", "CYS"]
        result = validate_coordination_chemistry("ZN", residues, 2)
        assert result["valid"] == True
        assert result["hsab_compatible"] == True

    def test_valid_lanthanide_coordination(self):
        """Valid lanthanide coordination with only O donors."""
        residues = ["GLU", "GLU", "ASP", "ASP", "GLU", "ASP", "ASN", "GLN"]
        result = validate_coordination_chemistry("TB", residues, 3)
        assert result["valid"] == True
        assert result["hsab_compatible"] == True

    def test_invalid_lanthanide_with_cys(self):
        """Lanthanide with Cys should fail HSAB validation."""
        residues = ["GLU", "GLU", "ASP", "CYS", "GLU", "ASP", "ASN", "GLN"]
        result = validate_coordination_chemistry("TB", residues, 3)
        assert result["hsab_compatible"] == False
        assert "cysteine" in result.get("warnings", [""])[0].lower() or \
               "soft" in result.get("warnings", [""])[0].lower() or \
               "incompatible" in result.get("warnings", [""])[0].lower()

    def test_invalid_lanthanide_with_his(self):
        """Lanthanide with His (N donor) may trigger warning."""
        residues = ["GLU", "GLU", "ASP", "HIS", "GLU", "ASP", "ASN", "GLN"]
        result = validate_coordination_chemistry("TB", residues, 3)
        # His is not ideal for lanthanides but may be tolerated
        # At minimum, it should be flagged in warnings
        if result["hsab_compatible"]:
            assert len(result.get("warnings", [])) > 0

    def test_valid_calcium_coordination(self):
        """Valid calcium coordination with O donors."""
        residues = ["GLU", "GLU", "ASP", "ASP", "ASN", "GLN"]
        result = validate_coordination_chemistry("CA", residues, 2)
        assert result["valid"] == True

    def test_calcium_with_cys_warning(self):
        """Calcium with Cys should trigger HSAB warning."""
        residues = ["GLU", "GLU", "ASP", "CYS"]
        result = validate_coordination_chemistry("CA", residues, 2)
        # Should flag HSAB incompatibility
        assert result["hsab_compatible"] == False or len(result.get("warnings", [])) > 0

    def test_coordination_number_too_low(self):
        """Too few coordinating residues should be flagged."""
        residues = ["HIS"]  # Only 1 residue
        result = validate_coordination_chemistry("ZN", residues, 2)
        assert result["valid"] == False or len(result.get("warnings", [])) > 0

    def test_coordination_number_too_high(self):
        """Too many coordinating residues should be flagged."""
        residues = ["HIS"] * 10  # 10 residues for zinc (way too many)
        result = validate_coordination_chemistry("ZN", residues, 2)
        assert len(result.get("warnings", [])) > 0


# =============================================================================
# METAL DATABASE STRUCTURE TESTS
# =============================================================================

class TestMetalDatabaseStructure:
    """Tests for the METAL_DATABASE structure."""

    def test_database_contains_required_metals(self):
        """Database should contain all required metals."""
        required = ["ZN", "FE", "CU", "MN", "CO", "NI", "CA", "MG",
                    "TB", "EU", "GD", "LA", "CE", "SM", "YB"]
        for metal in required:
            assert metal in METAL_DATABASE, f"Missing metal: {metal}"

    def test_database_entry_has_required_fields(self):
        """Each metal entry should have required fields."""
        required_fields = ["hsab_class", "preferred_donors", "coordination_numbers"]
        for metal, data in METAL_DATABASE.items():
            for field in required_fields:
                assert field in data, f"Metal {metal} missing field: {field}"

    def test_lanthanides_all_marked_hard(self):
        """All lanthanides should be marked as hard acids."""
        lanthanides = ["TB", "EU", "GD", "LA", "CE", "SM", "YB"]
        for metal in lanthanides:
            data = METAL_DATABASE[metal]
            # Check if there's a 3+ oxidation state marked as hard
            if isinstance(data["hsab_class"], dict):
                assert data["hsab_class"].get(3) == "hard", f"{metal} 3+ should be hard"
            else:
                assert data["hsab_class"] == "hard", f"{metal} should be hard"


# =============================================================================
# EDGE CASES AND ERROR HANDLING
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases and error handling."""

    def test_empty_residue_list(self):
        """Empty residue list should return valid=False."""
        result = validate_coordination_chemistry("ZN", [], 2)
        assert result["valid"] == False

    def test_non_coordinating_residue(self):
        """Non-coordinating residues like ALA should be flagged."""
        residues = ["ALA", "ALA", "ALA", "ALA"]
        result = validate_coordination_chemistry("ZN", residues, 2)
        assert result["valid"] == False or result["hsab_compatible"] == False

    def test_mixed_case_residue_names(self):
        """Residue names should be case insensitive."""
        residues = ["his", "His", "HIS", "cys"]
        result = validate_coordination_chemistry("ZN", residues, 2)
        # Should process without error
        assert "valid" in result

    def test_three_letter_and_one_letter_codes(self):
        """Should handle both 3-letter and 1-letter codes."""
        # Test that the bias string uses 1-letter codes
        bias = get_amino_acid_bias("ZN", 2, "catalytic")
        # All parts should be single letter codes
        for part in bias.split(","):
            aa = part.split(":")[0]
            assert len(aa) == 1
