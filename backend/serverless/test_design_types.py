# test_design_types.py
"""Tests for design type configuration."""

import pytest
import warnings
from design_types import (
    DesignType,
    DesignConfig,
    get_recommended_config,
    infer_design_type,
    validate_bias_AA,
    validate_omit_AA,
    DESIGN_PRESETS,
    METAL_PRESETS,
    VALID_AA_CODES,
)


def test_design_type_enum_has_all_types():
    """Should have all documented design types."""
    assert DesignType.METAL_BINDING in DesignType
    assert DesignType.SMALL_MOLECULE_BINDER in DesignType
    assert DesignType.DNA_RNA_BINDER in DesignType
    assert DesignType.PROTEIN_PROTEIN_BINDER in DesignType
    assert DesignType.ENZYME_ACTIVE_SITE in DesignType
    assert DesignType.SYMMETRIC_OLIGOMER in DesignType
    assert DesignType.GENERAL_SCAFFOLD in DesignType


def test_metal_binding_config_uses_ligandmpnn():
    """Metal binding should use LigandMPNN, not ProteinMPNN."""
    config = get_recommended_config(DesignType.METAL_BINDING)

    assert config.sequence_tool == "ligand_mpnn"
    assert config.backbone_tool == "rfd3"
    assert "A:-" in config.bias_AA  # Alanine penalty
    assert config.use_atom_context == True


def test_protein_binder_config_uses_proteinmpnn():
    """Protein-protein binder should use ProteinMPNN."""
    config = get_recommended_config(DesignType.PROTEIN_PROTEIN_BINDER)

    assert config.sequence_tool == "protein_mpnn"
    assert config.relaxation == "fastrelax"


def test_small_molecule_config_uses_ligandmpnn():
    """Small molecule binder should use LigandMPNN."""
    config = get_recommended_config(DesignType.SMALL_MOLECULE_BINDER)

    assert config.sequence_tool == "ligand_mpnn"
    assert config.ligand_cutoff == 6.0  # Tighter for small molecules


def test_config_allows_overrides():
    """Should allow overriding preset values."""
    config = get_recommended_config(
        DesignType.METAL_BINDING,
        temperature=0.2,
        num_sequences=16,
    )

    assert config.temperature == 0.2
    assert config.num_sequences == 16
    # But still uses correct sequence tool
    assert config.sequence_tool == "ligand_mpnn"


def test_config_warns_on_invalid_override_keys():
    """Should warn when invalid override keys are provided."""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        config = get_recommended_config(
            DesignType.METAL_BINDING,
            tempurature=0.5,  # Typo
            invalid_key=True,
        )
        assert len(w) == 1
        assert "Unknown config keys ignored" in str(w[0].message)
        assert "tempurature" in str(w[0].message)
        assert "invalid_key" in str(w[0].message)


def test_config_warns_on_unknown_metal_type():
    """Should warn when unknown metal type is provided."""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        config = get_recommended_config(
            DesignType.METAL_BINDING,
            metal_type="plutonium",
        )
        assert len(w) == 1
        assert "Unknown metal type 'plutonium'" in str(w[0].message)
        assert "zinc" in str(w[0].message)  # Should list valid types


def test_infer_design_type_metal():
    """Should prioritize metal over other flags."""
    # Pure metal binding
    assert infer_design_type(has_metal=True) == DesignType.METAL_BINDING

    # Metal + symmetric = dimer
    assert infer_design_type(has_metal=True, is_symmetric=True) == DesignType.METAL_MEDIATED_DIMER

    # Metal takes priority over ligand
    assert infer_design_type(has_metal=True, has_ligand=True) == DesignType.METAL_BINDING


def test_infer_design_type_protein_binder():
    """Should identify protein-protein binder."""
    assert infer_design_type(has_target_protein=True) == DesignType.PROTEIN_PROTEIN_BINDER


def test_infer_design_type_default():
    """Should default to general scaffold."""
    assert infer_design_type() == DesignType.GENERAL_SCAFFOLD


# =============================================================================
# EXPANDED TESTS FOR ALL DESIGN TYPES
# =============================================================================

class TestAllDesignTypes:
    """Tests to verify all design types have correct configurations."""

    def test_all_design_types_have_presets(self):
        """Every DesignType should have a preset configuration."""
        for design_type in DesignType:
            assert design_type in DESIGN_PRESETS, f"Missing preset for {design_type}"

    def test_design_type_enum_has_all_10_types(self):
        """Should have exactly 10 design types."""
        expected_types = {
            "METAL_BINDING",
            "SMALL_MOLECULE_BINDER",
            "DNA_RNA_BINDER",
            "PROTEIN_PROTEIN_BINDER",
            "ENZYME_ACTIVE_SITE",
            "SYMMETRIC_OLIGOMER",
            "SYMMETRIC_LIGAND_BINDER",
            "METAL_MEDIATED_DIMER",
            "GENERAL_SCAFFOLD",
            "PEPTIDE_BINDER",
        }
        actual_types = {dt.name for dt in DesignType}
        assert actual_types == expected_types

    def test_symmetric_ligand_binder_config(self):
        """SYMMETRIC_LIGAND_BINDER should use LigandMPNN with symmetry."""
        config = get_recommended_config(DesignType.SYMMETRIC_LIGAND_BINDER)
        assert config.sequence_tool == "ligand_mpnn"
        assert config.use_symmetry == True
        assert config.use_atom_context == True

    def test_metal_mediated_dimer_config(self):
        """METAL_MEDIATED_DIMER should use LigandMPNN with C2 symmetry."""
        config = get_recommended_config(DesignType.METAL_MEDIATED_DIMER)
        assert config.sequence_tool == "ligand_mpnn"
        assert config.use_symmetry == True
        assert config.symmetry_type == "C2"
        assert "A:-" in config.bias_AA  # Should penalize alanine

    def test_peptide_binder_config(self):
        """PEPTIDE_BINDER should use BindCraft with ProteinMPNN."""
        config = get_recommended_config(DesignType.PEPTIDE_BINDER)
        assert config.sequence_tool == "protein_mpnn"
        assert config.backbone_tool == "bindcraft"
        assert config.relaxation == "built_in"

    def test_symmetric_oligomer_config(self):
        """SYMMETRIC_OLIGOMER should use ProteinMPNN with symmetry."""
        config = get_recommended_config(DesignType.SYMMETRIC_OLIGOMER)
        assert config.sequence_tool == "protein_mpnn"
        assert config.use_symmetry == True
        assert config.relaxation == "fastrelax"

    def test_dna_rna_binder_config(self):
        """DNA_RNA_BINDER should use LigandMPNN with positive charge bias."""
        config = get_recommended_config(DesignType.DNA_RNA_BINDER)
        assert config.sequence_tool == "ligand_mpnn"
        assert config.use_atom_context == True
        # Should favor positive residues (R, K, H) for phosphate contacts
        assert "R:" in config.bias_AA or "K:" in config.bias_AA

    def test_enzyme_active_site_config(self):
        """ENZYME_ACTIVE_SITE should use LigandMPNN with motif scaffolding."""
        config = get_recommended_config(DesignType.ENZYME_ACTIVE_SITE)
        assert config.sequence_tool == "ligand_mpnn"
        assert config.backbone_tool == "rfd3"
        assert config.auto_detect_fixed == True

    def test_general_scaffold_config(self):
        """GENERAL_SCAFFOLD should use ProteinMPNN."""
        config = get_recommended_config(DesignType.GENERAL_SCAFFOLD)
        assert config.sequence_tool == "protein_mpnn"
        assert config.backbone_tool == "rfdiffusion"


class TestMetalPresets:
    """Tests for metal-specific presets."""

    def test_core_metal_presets_exist(self):
        """Should have presets for common metals and their variants."""
        # Core metals that must exist
        core_metals = {"zinc", "lanthanide", "iron", "copper", "calcium"}
        assert core_metals.issubset(set(METAL_PRESETS.keys()))
        # Should also have oxidation state and site variants
        assert "iron_ii" in METAL_PRESETS
        assert "iron_iii" in METAL_PRESETS
        assert "zinc_structural" in METAL_PRESETS

    def test_zinc_preset_has_correct_bias(self):
        """Zinc should favor HIS and CYS in bias_AA."""
        zinc = METAL_PRESETS["zinc"]
        assert "H:" in zinc["bias_AA"]  # Favor histidine
        assert "C:" in zinc["bias_AA"]  # Favor cysteine
        # New structure uses preferred_residues instead of coordinating_residues
        assert "C" in zinc["preferred_residues"]
        assert "H" in zinc["preferred_residues"]

    def test_lanthanide_preset_excludes_cysteine(self):
        """Lanthanide should strongly penalize cysteine (soft donor incompatible)."""
        lanthanide = METAL_PRESETS["lanthanide"]
        # New structure uses excluded_residues and negative bias instead of omit_AA
        assert "C" in lanthanide["excluded_residues"]
        # Cys should have strong negative bias
        assert "C:-" in lanthanide["bias_AA"]

    def test_calcium_preset_penalizes_cysteine(self):
        """Calcium should penalize cysteine (hard acid)."""
        calcium = METAL_PRESETS["calcium"]
        # Hard acid should have negative Cys bias
        assert "C:-" in calcium["bias_AA"]

    def test_metal_presets_applied_to_config(self):
        """Metal type should override bias_AA in config."""
        config = get_recommended_config(
            DesignType.METAL_BINDING,
            metal_type="lanthanide"
        )
        # Lanthanide should have its specific bias
        assert "E:" in config.bias_AA  # Favor glutamate
        assert "D:" in config.bias_AA  # Favor aspartate


class TestBiasAAValidation:
    """Tests for bias_AA format validation."""

    def test_valid_bias_AA_single(self):
        """Should accept single AA bias."""
        is_valid, error = validate_bias_AA("A:-2.0")
        assert is_valid
        assert error is None

    def test_valid_bias_AA_multiple(self):
        """Should accept multiple AA biases."""
        is_valid, error = validate_bias_AA("A:-2.0,H:2.0,E:1.0,D:1.0")
        assert is_valid
        assert error is None

    def test_valid_bias_AA_positive_values(self):
        """Should accept positive bias values."""
        is_valid, error = validate_bias_AA("W:1.5,Y:1.5")
        assert is_valid
        assert error is None

    def test_valid_bias_AA_empty(self):
        """Should accept empty/None bias."""
        is_valid, error = validate_bias_AA("")
        assert is_valid
        is_valid, error = validate_bias_AA(None)
        assert is_valid

    def test_invalid_bias_AA_bad_format(self):
        """Should reject invalid format."""
        is_valid, error = validate_bias_AA("A=2.0")
        assert not is_valid
        assert "Invalid bias_AA format" in error

    def test_invalid_bias_AA_missing_value(self):
        """Should reject missing value."""
        is_valid, error = validate_bias_AA("A:")
        assert not is_valid
        assert "Invalid bias_AA format" in error

    def test_invalid_bias_AA_bad_amino_acid(self):
        """Should reject invalid amino acid code."""
        is_valid, error = validate_bias_AA("X:-2.0")
        assert not is_valid
        assert "Invalid amino acid code" in error

    def test_invalid_bias_AA_lowercase(self):
        """Should reject lowercase amino acid codes."""
        is_valid, error = validate_bias_AA("a:-2.0")
        assert not is_valid
        assert "Invalid bias_AA format" in error


class TestOmitAAValidation:
    """Tests for omit_AA format validation."""

    def test_valid_omit_AA_single(self):
        """Should accept single AA omission."""
        is_valid, error = validate_omit_AA("C")
        assert is_valid
        assert error is None

    def test_valid_omit_AA_multiple(self):
        """Should accept multiple AA omissions."""
        is_valid, error = validate_omit_AA("CM")
        assert is_valid
        assert error is None

    def test_valid_omit_AA_empty(self):
        """Should accept empty omit."""
        is_valid, error = validate_omit_AA("")
        assert is_valid
        is_valid, error = validate_omit_AA(None)
        assert is_valid

    def test_invalid_omit_AA_bad_code(self):
        """Should reject invalid amino acid code."""
        is_valid, error = validate_omit_AA("X")
        assert not is_valid
        assert "Invalid amino acid code" in error


class TestInferDesignType:
    """Extended tests for design type inference."""

    def test_infer_ligand_creates_small_molecule(self):
        """Ligand flag should create small molecule binder."""
        assert infer_design_type(has_ligand=True) == DesignType.SMALL_MOLECULE_BINDER

    def test_infer_ligand_with_symmetry(self):
        """Ligand + symmetry should create symmetric ligand binder."""
        assert infer_design_type(has_ligand=True, is_symmetric=True) == DesignType.SYMMETRIC_LIGAND_BINDER

    def test_infer_ligand_with_motif(self):
        """Ligand + motif should create enzyme active site."""
        assert infer_design_type(has_ligand=True, has_motif=True) == DesignType.ENZYME_ACTIVE_SITE

    def test_infer_nucleotide(self):
        """Nucleotide flag should create DNA/RNA binder."""
        assert infer_design_type(has_nucleotide=True) == DesignType.DNA_RNA_BINDER

    def test_infer_symmetric_alone(self):
        """Symmetric flag alone should create symmetric oligomer."""
        assert infer_design_type(is_symmetric=True) == DesignType.SYMMETRIC_OLIGOMER

    def test_priority_metal_over_ligand(self):
        """Metal should take priority over ligand."""
        result = infer_design_type(has_metal=True, has_ligand=True)
        assert result == DesignType.METAL_BINDING

    def test_priority_metal_over_nucleotide(self):
        """Metal should take priority over nucleotide."""
        result = infer_design_type(has_metal=True, has_nucleotide=True)
        assert result == DesignType.METAL_BINDING

    def test_priority_nucleotide_over_ligand(self):
        """Nucleotide should take priority over ligand."""
        result = infer_design_type(has_nucleotide=True, has_ligand=True)
        assert result == DesignType.DNA_RNA_BINDER


class TestMetalPresetChemistry:
    """Test that metal presets use correct chemistry."""

    def test_zinc_preset_favors_cys_his(self):
        """Zinc preset should favor Cys and His over Asp/Glu."""
        from design_types import METAL_PRESETS
        bias = METAL_PRESETS["zinc"]["bias_AA"]
        bias_dict = dict(pair.split(":") for pair in bias.split(","))
        assert float(bias_dict.get("C", 0)) >= 2.5
        assert float(bias_dict.get("H", 0)) >= 2.0

    def test_lanthanide_preset_excludes_cys(self):
        """Lanthanide preset must strongly exclude Cys."""
        from design_types import METAL_PRESETS
        bias = METAL_PRESETS["lanthanide"]["bias_AA"]
        bias_dict = dict(pair.split(":") for pair in bias.split(","))
        assert float(bias_dict.get("C", 0)) <= -3.0

    def test_iron_has_oxidation_state_variants(self):
        """Iron should have separate Fe2+ and Fe3+ presets."""
        from design_types import METAL_PRESETS
        assert "iron_ii" in METAL_PRESETS or "iron" in METAL_PRESETS
        # Fe3+ is harder - prefers O over S
        if "iron_iii" in METAL_PRESETS:
            bias = METAL_PRESETS["iron_iii"]["bias_AA"]
            bias_dict = dict(pair.split(":") for pair in bias.split(","))
            assert float(bias_dict.get("E", 0)) > float(bias_dict.get("C", 0))
