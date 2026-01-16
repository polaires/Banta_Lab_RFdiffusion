# test_design_types.py
"""Tests for design type configuration."""

import pytest
import warnings
from design_types import (
    DesignType,
    DesignConfig,
    get_recommended_config,
    infer_design_type,
    DESIGN_PRESETS,
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
