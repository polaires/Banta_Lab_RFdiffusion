# test_design_types.py
"""Tests for design type configuration."""

import pytest
from design_types import (
    DesignType,
    DesignConfig,
    get_recommended_config,
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
