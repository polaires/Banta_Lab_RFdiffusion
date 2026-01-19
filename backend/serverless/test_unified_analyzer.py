"""Tests for UnifiedDesignAnalyzer."""
import pytest
import json
from unittest.mock import patch, MagicMock
from unified_analyzer import UnifiedDesignAnalyzer
from analysis_types import AnalysisStatus, DesignType


# Sample PDB content for testing
SAMPLE_DIMER_PDB = """ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00  0.00           O
ATOM      5  N   ALA B   1      10.000   0.000   0.000  1.00  0.00           N
ATOM      6  CA  ALA B   1      11.458   0.000   0.000  1.00  0.00           C
END
"""

SAMPLE_MONOMER_PDB = """ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
END
"""


class TestUnifiedDesignAnalyzer:
    """Test suite for UnifiedDesignAnalyzer."""

    def test_init(self):
        analyzer = UnifiedDesignAnalyzer()
        assert analyzer is not None
        assert hasattr(analyzer, "analyze")

    def test_detect_design_type_dimer(self):
        analyzer = UnifiedDesignAnalyzer()
        dtype = analyzer._detect_design_type(SAMPLE_DIMER_PDB)
        assert dtype == DesignType.PROTEIN_DIMER

    def test_detect_design_type_monomer(self):
        analyzer = UnifiedDesignAnalyzer()
        dtype = analyzer._detect_design_type(SAMPLE_MONOMER_PDB)
        assert dtype == DesignType.MONOMER

    def test_analyze_returns_structured_json(self):
        analyzer = UnifiedDesignAnalyzer()
        result = analyzer.analyze(
            pdb_content=SAMPLE_DIMER_PDB,
            design_params={"name": "test"}
        )
        assert "design_id" in result
        assert "design_type" in result
        assert "timestamp" in result
        assert "analyses" in result
        assert isinstance(result["analyses"], dict)

    def test_analyze_includes_applicable_analyses(self):
        analyzer = UnifiedDesignAnalyzer()
        result = analyzer.analyze(
            pdb_content=SAMPLE_DIMER_PDB,
            design_params={}
        )
        # For a simple dimer, interface analysis should be applicable
        assert "interface_quality" in result["analyses"]

    def test_analyze_marks_inapplicable_analyses(self):
        analyzer = UnifiedDesignAnalyzer()
        result = analyzer.analyze(
            pdb_content=SAMPLE_DIMER_PDB,
            design_params={}
        )
        # No metal present, so metal analysis should be not_applicable
        metal_result = result["analyses"].get("metal_coordination", {})
        assert metal_result.get("status") in ["not_applicable", "skipped"]

    def test_analyze_with_metal_type(self):
        analyzer = UnifiedDesignAnalyzer()
        result = analyzer.analyze(
            pdb_content=SAMPLE_DIMER_PDB,
            design_params={},
            metal_type="TB"
        )
        assert result is not None
        assert "metal_coordination" in result["analyses"]
