"""Tests for analysis type definitions."""
import pytest
from analysis_types import (
    AnalysisStatus,
    AnalysisResult,
    DesignType,
    detect_design_type,
)


def test_analysis_status_success():
    result = AnalysisResult.success({"score": 0.85})
    assert result.status == AnalysisStatus.SUCCESS
    assert result.metrics == {"score": 0.85}
    assert result.reason is None


def test_analysis_status_not_applicable():
    result = AnalysisResult.not_applicable("no ligand present")
    assert result.status == AnalysisStatus.NOT_APPLICABLE
    assert result.metrics is None
    assert result.reason == "no ligand present"


def test_analysis_status_skipped():
    result = AnalysisResult.skipped("GNINA binary not found")
    assert result.status == AnalysisStatus.SKIPPED
    assert result.reason == "GNINA binary not found"


def test_analysis_result_to_dict():
    result = AnalysisResult.success({"plddt": 0.9})
    d = result.to_dict()
    assert d["status"] == "success"
    assert d["metrics"]["plddt"] == 0.9


def test_detect_design_type_ligand_dimer():
    dtype = detect_design_type(
        has_ligand=True,
        has_metal=False,
        chain_count=2
    )
    assert dtype == DesignType.LIGAND_INTERFACE_DIMER


def test_detect_design_type_metal_dimer():
    dtype = detect_design_type(
        has_ligand=False,
        has_metal=True,
        chain_count=2
    )
    assert dtype == DesignType.METAL_INTERFACE_DIMER


def test_detect_design_type_metal_ligand_dimer():
    dtype = detect_design_type(
        has_ligand=True,
        has_metal=True,
        chain_count=2
    )
    assert dtype == DesignType.METAL_LIGAND_INTERFACE_DIMER


def test_detect_design_type_monomer():
    dtype = detect_design_type(
        has_ligand=False,
        has_metal=False,
        chain_count=1
    )
    assert dtype == DesignType.MONOMER
