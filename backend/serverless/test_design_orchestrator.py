# test_design_orchestrator.py
"""Tests for DesignOrchestrator."""

import pytest
from design_orchestrator import DesignOrchestrator
from design_types import DesignType


def test_orchestrator_initializes_with_design_type():
    """Should initialize with inferred config."""
    orchestrator = DesignOrchestrator(design_type=DesignType.METAL_BINDING)

    assert orchestrator.config.sequence_tool == "ligand_mpnn"
    assert orchestrator.config.use_atom_context == True


def test_orchestrator_builds_mpnn_request():
    """Should build valid MPNN request."""
    orchestrator = DesignOrchestrator(
        design_type=DesignType.METAL_BINDING,
        metal_type="zinc"
    )

    pdb = "ATOM...HETATM...ZN..."
    request = orchestrator.build_mpnn_request(pdb)

    assert request["model_type"] == "ligand_mpnn"
    assert "A:-" in request["bias_AA"]
    assert request["temperature"] == 0.1


def test_orchestrator_builds_rfd3_config():
    """Should build valid RFD3 config for backbone generation."""
    orchestrator = DesignOrchestrator(
        design_type=DesignType.METAL_MEDIATED_DIMER,
        metal_type="zinc"
    )

    rfd3_config = orchestrator.build_rfd3_config(
        ligand_code="ZN",
        chain_length="60-80",
    )

    assert "symmetry" in rfd3_config
    assert rfd3_config["symmetry"]["type"] == "C2"
    assert rfd3_config["ligand"] == "ZN"


def test_orchestrator_recommends_validation():
    """Should recommend appropriate validation method."""
    # Metal binding
    orchestrator = DesignOrchestrator(design_type=DesignType.METAL_BINDING)
    assert "af2" in orchestrator.get_validation_method()

    # Protein-protein
    orchestrator = DesignOrchestrator(design_type=DesignType.PROTEIN_PROTEIN_BINDER)
    assert "multimer" in orchestrator.get_validation_method()
