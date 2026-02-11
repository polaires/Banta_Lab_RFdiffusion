# test_unified_design.py
"""Tests for unified design API endpoint."""

import pytest
import sys
from unittest.mock import MagicMock

# Mock external dependencies that may not be available in test environment
sys.modules['runpod'] = MagicMock()


def test_design_request_infers_type():
    """Should infer design type from request parameters."""
    from handler import infer_design_type_from_request
    from design_types import DesignType

    # Metal binding
    request = {"metal": "ZN", "pdb_content": "..."}
    design_type = infer_design_type_from_request(request)
    assert design_type == DesignType.METAL_BINDING

    # Protein binder
    request = {"target_pdb": "...", "pdb_content": "..."}
    design_type = infer_design_type_from_request(request)
    assert design_type == DesignType.PROTEIN_PROTEIN_BINDER

    # Metal dimer
    request = {"metal": "ZN", "symmetry": "C2", "pdb_content": "..."}
    design_type = infer_design_type_from_request(request)
    assert design_type == DesignType.METAL_MEDIATED_DIMER


def test_design_request_infers_small_molecule():
    """Should infer small molecule binder from ligand parameter."""
    from handler import infer_design_type_from_request
    from design_types import DesignType

    request = {"ligand": "AZO", "pdb_content": "..."}
    design_type = infer_design_type_from_request(request)
    assert design_type == DesignType.SMALL_MOLECULE_BINDER


def test_design_request_defaults_to_general():
    """Should default to general scaffold when no specific parameters."""
    from handler import infer_design_type_from_request
    from design_types import DesignType

    request = {"pdb_content": "..."}
    design_type = infer_design_type_from_request(request)
    assert design_type == DesignType.GENERAL_SCAFFOLD


def test_design_request_validates_temperature():
    """Should reject invalid temperature."""
    from handler import handle_unified_design

    result = handle_unified_design({
        "pdb_content": "ATOM...",
        "temperature": "not_a_number"
    })

    assert result["status"] == "failed"
    assert "temperature" in result["error"]


def test_design_request_validates_bias_AA():
    """Should reject invalid bias_AA format."""
    from handler import handle_unified_design

    result = handle_unified_design({
        "pdb_content": "ATOM...",
        "bias_AA": "invalid_format"
    })

    assert result["status"] == "failed"
    assert "bias_aa" in result["error"].lower()


def test_design_request_validates_bias_AA_bad_amino_acid():
    """Should reject invalid amino acid code in bias_AA."""
    from handler import handle_unified_design

    result = handle_unified_design({
        "pdb_content": "ATOM...",
        "bias_AA": "X:-2.0"  # X is not a valid AA code
    })

    assert result["status"] == "failed"
    assert "amino acid" in result["error"].lower()


def test_design_request_validates_omit_AA():
    """Should reject invalid omit_AA codes."""
    from handler import handle_unified_design

    result = handle_unified_design({
        "pdb_content": "ATOM...",
        "omit_AA": "XZ"  # Invalid AA codes
    })

    assert result["status"] == "failed"
    assert "amino acid" in result["error"].lower()


def test_design_request_accepts_valid_bias_AA():
    """Should accept valid bias_AA format."""
    from handler import handle_unified_design

    # This will fail on MPNN execution but should pass validation
    result = handle_unified_design({
        "pdb_content": "ATOM...",
        "bias_AA": "A:-2.0,H:2.0,E:1.0"
    })

    # If it fails, should NOT be due to bias_AA validation
    if result.get("status") == "failed":
        assert "bias_AA" not in result.get("error", "").lower()


def test_design_request_accepts_valid_omit_AA():
    """Should accept valid omit_AA codes."""
    from handler import handle_unified_design

    result = handle_unified_design({
        "pdb_content": "ATOM...",
        "omit_AA": "CM"  # Valid: Cysteine and Methionine
    })

    # If it fails, should NOT be due to omit_AA validation
    if result.get("status") == "failed":
        assert "omit_AA" not in result.get("error", "").lower()
        assert "amino acid" not in result.get("error", "").lower()
