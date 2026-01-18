# Unified Design Analyzer Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a unified local testing and analysis infrastructure for ligand/metal/metal-ligand dimer designs with persistent learning.

**Architecture:** UnifiedDesignAnalyzer orchestrates existing analysis tools (GNINA, metal_validation, binding_analysis) into a single entry point producing structured JSON. DesignHistoryManager handles persistent storage with session tracking, filter presets, and automatic lesson synthesis.

**Tech Stack:** Python 3.12, existing backend modules (binding_analysis.py, metal_chemistry.py, metal_validation.py, tebl_analysis.py, hotspot_detection.py), JSON for storage, pytest for testing.

---

## Phase 1: Foundation

### Task 1: Create Directory Structure

**Files:**
- Create: `experiments/design_history/`
- Create: `experiments/design_history/runs/`
- Create: `experiments/design_history/sessions/`
- Create: `experiments/design_history/lessons/`
- Create: `experiments/design_history/exports/`
- Create: `experiments/design_history/filter_presets/`

**Step 1: Create directory structure**

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion
mkdir -p experiments/design_history/runs
mkdir -p experiments/design_history/sessions
mkdir -p experiments/design_history/lessons
mkdir -p experiments/design_history/exports
mkdir -p experiments/design_history/filter_presets
```

**Step 2: Create .gitkeep files to preserve empty directories**

```bash
touch experiments/design_history/runs/.gitkeep
touch experiments/design_history/sessions/.gitkeep
touch experiments/design_history/lessons/.gitkeep
touch experiments/design_history/exports/.gitkeep
```

**Step 3: Create initial index.json**

Create file `experiments/design_history/index.json`:
```json
{
  "version": "1.0.0",
  "created": "2026-01-18T00:00:00Z",
  "designs": []
}
```

**Step 4: Commit**

```bash
git add experiments/design_history/
git commit -m "feat: add design_history directory structure"
```

---

### Task 2: Create Filter Presets

**Files:**
- Create: `experiments/design_history/filter_presets/default.json`
- Create: `experiments/design_history/filter_presets/relaxed.json`
- Create: `experiments/design_history/filter_presets/stringent.json`
- Create: `experiments/design_history/filter_presets/metal_specific.json`
- Create: `experiments/design_history/filter_presets/ligand_specific.json`

**Step 1: Create default.json**

Create file `experiments/design_history/filter_presets/default.json`:
```json
{
  "name": "default",
  "description": "Standard thresholds for general design filtering",
  "filters": {
    "structure_confidence": {
      "plddt": {"min": 0.8, "direction": "higher_better"},
      "pae_mean": {"max": 10.0, "direction": "lower_better"}
    },
    "interface_quality": {
      "shape_complementarity": {"min": 0.6, "direction": "higher_better"},
      "dG": {"max": 0, "direction": "lower_better"},
      "dSASA": {"min": 500, "direction": "higher_better"},
      "interface_hbonds": {"min": 3, "direction": "higher_better"},
      "unsatisfied_hbonds": {"max": 4, "direction": "lower_better"},
      "surface_hydrophobicity": {"max": 0.35, "direction": "lower_better"}
    },
    "metal_coordination": {
      "coordination_distance": {"max": 2.5, "direction": "lower_better"},
      "hsab_compatible": {"equals": true}
    },
    "ligand_binding": {
      "gnina_affinity": {"max": -6.0, "direction": "lower_better"}
    },
    "symmetry": {
      "c2_score": {"min": 0.7, "direction": "higher_better"}
    }
  }
}
```

**Step 2: Create relaxed.json**

Create file `experiments/design_history/filter_presets/relaxed.json`:
```json
{
  "name": "relaxed",
  "description": "Permissive thresholds for early exploration",
  "filters": {
    "structure_confidence": {
      "plddt": {"min": 0.7, "direction": "higher_better"},
      "pae_mean": {"max": 15.0, "direction": "lower_better"}
    },
    "interface_quality": {
      "shape_complementarity": {"min": 0.5, "direction": "higher_better"},
      "dG": {"max": 5, "direction": "lower_better"},
      "dSASA": {"min": 300, "direction": "higher_better"},
      "interface_hbonds": {"min": 1, "direction": "higher_better"},
      "unsatisfied_hbonds": {"max": 6, "direction": "lower_better"},
      "surface_hydrophobicity": {"max": 0.45, "direction": "lower_better"}
    },
    "metal_coordination": {
      "coordination_distance": {"max": 2.8, "direction": "lower_better"},
      "hsab_compatible": {"equals": true}
    },
    "ligand_binding": {
      "gnina_affinity": {"max": -4.0, "direction": "lower_better"}
    },
    "symmetry": {
      "c2_score": {"min": 0.5, "direction": "higher_better"}
    }
  }
}
```

**Step 3: Create stringent.json**

Create file `experiments/design_history/filter_presets/stringent.json`:
```json
{
  "name": "stringent",
  "description": "High bar for final candidates",
  "filters": {
    "structure_confidence": {
      "plddt": {"min": 0.85, "direction": "higher_better"},
      "pae_mean": {"max": 6.0, "direction": "lower_better"}
    },
    "interface_quality": {
      "shape_complementarity": {"min": 0.7, "direction": "higher_better"},
      "dG": {"max": -10, "direction": "lower_better"},
      "dSASA": {"min": 700, "direction": "higher_better"},
      "interface_hbonds": {"min": 5, "direction": "higher_better"},
      "unsatisfied_hbonds": {"max": 2, "direction": "lower_better"},
      "surface_hydrophobicity": {"max": 0.30, "direction": "lower_better"}
    },
    "metal_coordination": {
      "coordination_distance": {"max": 2.3, "direction": "lower_better"},
      "hsab_compatible": {"equals": true}
    },
    "ligand_binding": {
      "gnina_affinity": {"max": -8.0, "direction": "lower_better"}
    },
    "symmetry": {
      "c2_score": {"min": 0.85, "direction": "higher_better"}
    }
  }
}
```

**Step 4: Create metal_specific.json**

Create file `experiments/design_history/filter_presets/metal_specific.json`:
```json
{
  "name": "metal_specific",
  "description": "Optimized for metal coordination designs",
  "filters": {
    "structure_confidence": {
      "plddt": {"min": 0.8, "direction": "higher_better"}
    },
    "metal_coordination": {
      "coordination_distance": {"max": 2.5, "direction": "lower_better"},
      "coordination_number": {"min": 4, "direction": "higher_better"},
      "geometry_rmsd": {"max": 15.0, "direction": "lower_better"},
      "hsab_compatible": {"equals": true},
      "carboxylate_donors": {"min": 4, "direction": "higher_better"}
    },
    "tebl": {
      "has_antenna": {"equals": true},
      "predicted_efficiency": {"min": 0.3, "direction": "higher_better"}
    }
  }
}
```

**Step 5: Create ligand_specific.json**

Create file `experiments/design_history/filter_presets/ligand_specific.json`:
```json
{
  "name": "ligand_specific",
  "description": "Optimized for small molecule binding designs",
  "filters": {
    "structure_confidence": {
      "plddt": {"min": 0.8, "direction": "higher_better"}
    },
    "ligand_binding": {
      "gnina_affinity": {"max": -6.0, "direction": "lower_better"},
      "gnina_cnn_score": {"min": 0.5, "direction": "higher_better"}
    },
    "interface_quality": {
      "shape_complementarity": {"min": 0.6, "direction": "higher_better"},
      "contacts": {"min": 10, "direction": "higher_better"}
    },
    "symmetry": {
      "c2_score": {"min": 0.7, "direction": "higher_better"}
    }
  }
}
```

**Step 6: Commit**

```bash
git add experiments/design_history/filter_presets/
git commit -m "feat: add BindCraft-inspired filter presets"
```

---

### Task 3: Create Analysis Status Types

**Files:**
- Create: `backend/serverless/analysis_types.py`
- Test: `backend/serverless/test_analysis_types.py`

**Step 1: Write the failing test**

Create file `backend/serverless/test_analysis_types.py`:
```python
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
```

**Step 2: Run test to verify it fails**

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless
python -m pytest test_analysis_types.py -v
```

Expected: FAIL with "ModuleNotFoundError: No module named 'analysis_types'"

**Step 3: Write minimal implementation**

Create file `backend/serverless/analysis_types.py`:
```python
"""
Analysis Type Definitions

Provides structured types for design analysis results with clear status handling.
"""
from enum import Enum
from typing import Dict, Any, Optional
from dataclasses import dataclass


class AnalysisStatus(Enum):
    """Status of an analysis run."""
    SUCCESS = "success"
    NOT_APPLICABLE = "not_applicable"
    SKIPPED = "skipped"


class DesignType(Enum):
    """Types of protein designs."""
    MONOMER = "monomer"
    PROTEIN_DIMER = "protein_dimer"
    LIGAND_INTERFACE_DIMER = "ligand_interface_dimer"
    METAL_INTERFACE_DIMER = "metal_interface_dimer"
    METAL_LIGAND_INTERFACE_DIMER = "metal_ligand_interface_dimer"


@dataclass
class AnalysisResult:
    """Result of a single analysis module."""
    status: AnalysisStatus
    metrics: Optional[Dict[str, Any]] = None
    reason: Optional[str] = None

    @classmethod
    def success(cls, metrics: Dict[str, Any]) -> "AnalysisResult":
        """Create a successful analysis result."""
        return cls(status=AnalysisStatus.SUCCESS, metrics=metrics)

    @classmethod
    def not_applicable(cls, reason: str) -> "AnalysisResult":
        """Create a not-applicable result (analysis doesn't apply to this design)."""
        return cls(status=AnalysisStatus.NOT_APPLICABLE, reason=reason)

    @classmethod
    def skipped(cls, reason: str) -> "AnalysisResult":
        """Create a skipped result (analysis could apply but couldn't run)."""
        return cls(status=AnalysisStatus.SKIPPED, reason=reason)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        result = {"status": self.status.value}
        if self.metrics is not None:
            result["metrics"] = self.metrics
        if self.reason is not None:
            result["reason"] = self.reason
        return result


def detect_design_type(
    has_ligand: bool,
    has_metal: bool,
    chain_count: int,
) -> DesignType:
    """
    Detect design type from structural properties.

    Args:
        has_ligand: Whether a small molecule ligand is present
        has_metal: Whether a metal ion is present
        chain_count: Number of protein chains

    Returns:
        DesignType enum value
    """
    is_dimer = chain_count >= 2

    if has_ligand and has_metal and is_dimer:
        return DesignType.METAL_LIGAND_INTERFACE_DIMER
    elif has_ligand and is_dimer:
        return DesignType.LIGAND_INTERFACE_DIMER
    elif has_metal and is_dimer:
        return DesignType.METAL_INTERFACE_DIMER
    elif is_dimer:
        return DesignType.PROTEIN_DIMER
    else:
        return DesignType.MONOMER
```

**Step 4: Run test to verify it passes**

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless
python -m pytest test_analysis_types.py -v
```

Expected: All tests PASS

**Step 5: Commit**

```bash
git add backend/serverless/analysis_types.py backend/serverless/test_analysis_types.py
git commit -m "feat: add analysis type definitions with status handling"
```

---

### Task 4: Create UnifiedDesignAnalyzer Core

**Files:**
- Create: `backend/serverless/unified_analyzer.py`
- Test: `backend/serverless/test_unified_analyzer.py`

**Step 1: Write the failing test**

Create file `backend/serverless/test_unified_analyzer.py`:
```python
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
```

**Step 2: Run test to verify it fails**

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless
python -m pytest test_unified_analyzer.py -v
```

Expected: FAIL with "ModuleNotFoundError: No module named 'unified_analyzer'"

**Step 3: Write minimal implementation**

Create file `backend/serverless/unified_analyzer.py`:
```python
"""
Unified Design Analyzer

Orchestrates all analysis tools into a single entry point producing structured JSON.
Auto-detects design type and runs appropriate analyses.
"""
import os
import json
from datetime import datetime
from typing import Dict, Any, Optional, List
from dataclasses import dataclass, field

from analysis_types import AnalysisResult, AnalysisStatus, DesignType, detect_design_type


class UnifiedDesignAnalyzer:
    """
    Unified analyzer that orchestrates all design analysis tools.

    Usage:
        analyzer = UnifiedDesignAnalyzer()
        result = analyzer.analyze(pdb_content, design_params, metal_type="TB")
    """

    def __init__(self):
        """Initialize analyzer with available analysis modules."""
        self._analysis_modules = self._discover_modules()

    def _discover_modules(self) -> Dict[str, bool]:
        """Discover which analysis modules are available."""
        modules = {
            "binding_analysis": False,
            "metal_validation": False,
            "metal_chemistry": False,
            "tebl_analysis": False,
            "hotspot_detection": False,
            "topology_validation": False,
            "gnina": False,
        }

        # Check for Python modules
        try:
            import binding_analysis
            modules["binding_analysis"] = True
        except ImportError:
            pass

        try:
            import metal_validation
            modules["metal_validation"] = True
        except ImportError:
            pass

        try:
            import metal_chemistry
            modules["metal_chemistry"] = True
        except ImportError:
            pass

        try:
            import tebl_analysis
            modules["tebl_analysis"] = True
        except ImportError:
            pass

        try:
            import hotspot_detection
            modules["hotspot_detection"] = True
        except ImportError:
            pass

        try:
            import topology_validation
            modules["topology_validation"] = True
        except ImportError:
            pass

        # Check for GNINA binary
        gnina_paths = ["/usr/local/bin/gnina", "/app/gnina", "gnina"]
        for path in gnina_paths:
            if os.path.exists(path) or os.system(f"which {path} > /dev/null 2>&1") == 0:
                modules["gnina"] = True
                break

        return modules

    def _detect_design_type(
        self,
        pdb_content: str,
        metal_type: Optional[str] = None,
        ligand_sdf: Optional[str] = None,
    ) -> DesignType:
        """Detect design type from PDB content and provided metadata."""
        # Count chains
        chains = set()
        has_hetatm = False

        for line in pdb_content.split("\n"):
            if line.startswith("ATOM"):
                chain_id = line[21] if len(line) > 21 else ""
                if chain_id.strip():
                    chains.add(chain_id)
            elif line.startswith("HETATM"):
                has_hetatm = True

        chain_count = len(chains)
        has_ligand = ligand_sdf is not None or has_hetatm
        has_metal = metal_type is not None

        return detect_design_type(
            has_ligand=has_ligand,
            has_metal=has_metal,
            chain_count=chain_count,
        )

    def _generate_design_id(self, design_type: DesignType) -> str:
        """Generate unique design ID."""
        timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
        type_short = design_type.value.replace("_", "-")
        return f"{timestamp}_{type_short}"

    def analyze(
        self,
        pdb_content: str,
        design_params: Dict[str, Any],
        pdb_path: Optional[str] = None,
        ligand_sdf: Optional[str] = None,
        metal_type: Optional[str] = None,
        metal_chain: str = "L",
        metal_resnum: int = 1,
    ) -> Dict[str, Any]:
        """
        Run comprehensive analysis on a design.

        Args:
            pdb_content: PDB file content as string
            design_params: Parameters used to generate the design
            pdb_path: Optional path to PDB file (for tools that need file path)
            ligand_sdf: Optional SDF content for ligand
            metal_type: Optional metal type code (e.g., "TB", "ZN")
            metal_chain: Chain ID of metal (default "L")
            metal_resnum: Residue number of metal (default 1)

        Returns:
            Structured metrics JSON with all analysis results
        """
        # Detect design type
        design_type = self._detect_design_type(pdb_content, metal_type, ligand_sdf)
        design_id = self._generate_design_id(design_type)

        # Initialize result structure
        result = {
            "design_id": design_id,
            "design_type": design_type.value,
            "timestamp": datetime.now().isoformat(),
            "design_params": design_params,
            "analyses": {},
        }

        # Run structure confidence analysis
        result["analyses"]["structure_confidence"] = self._analyze_structure_confidence(
            pdb_content
        ).to_dict()

        # Run interface analysis (for dimers)
        if design_type != DesignType.MONOMER:
            result["analyses"]["interface_quality"] = self._analyze_interface(
                pdb_content
            ).to_dict()
        else:
            result["analyses"]["interface_quality"] = AnalysisResult.not_applicable(
                "monomer design - no interface"
            ).to_dict()

        # Run metal coordination analysis
        if metal_type:
            result["analyses"]["metal_coordination"] = self._analyze_metal_coordination(
                pdb_content, metal_type, metal_chain, metal_resnum
            ).to_dict()
        else:
            result["analyses"]["metal_coordination"] = AnalysisResult.not_applicable(
                "no metal type specified"
            ).to_dict()

        # Run ligand binding analysis
        if ligand_sdf:
            result["analyses"]["ligand_binding"] = self._analyze_ligand_binding(
                pdb_content, pdb_path, ligand_sdf
            ).to_dict()
        else:
            result["analyses"]["ligand_binding"] = AnalysisResult.not_applicable(
                "no ligand provided"
            ).to_dict()

        # Run topology validation
        result["analyses"]["topology"] = self._analyze_topology(pdb_content).to_dict()

        # Run symmetry analysis (for dimers)
        if design_type != DesignType.MONOMER:
            result["analyses"]["symmetry"] = self._analyze_symmetry(pdb_content).to_dict()
        else:
            result["analyses"]["symmetry"] = AnalysisResult.not_applicable(
                "monomer design - no symmetry analysis"
            ).to_dict()

        return result

    def _analyze_structure_confidence(self, pdb_content: str) -> AnalysisResult:
        """Extract structure confidence metrics (pLDDT, pAE)."""
        # Try to extract pLDDT from B-factor column (common convention)
        plddt_values = []

        for line in pdb_content.split("\n"):
            if line.startswith("ATOM"):
                try:
                    bfactor = float(line[60:66].strip())
                    # pLDDT is typically stored as 0-100 in B-factor
                    if 0 <= bfactor <= 100:
                        plddt_values.append(bfactor / 100.0)
                except (ValueError, IndexError):
                    pass

        if plddt_values:
            return AnalysisResult.success({
                "plddt": sum(plddt_values) / len(plddt_values),
                "plddt_min": min(plddt_values),
                "plddt_max": max(plddt_values),
            })
        else:
            return AnalysisResult.skipped("could not extract pLDDT from B-factors")

    def _analyze_interface(self, pdb_content: str) -> AnalysisResult:
        """Analyze protein-protein interface."""
        if not self._analysis_modules.get("binding_analysis"):
            return AnalysisResult.skipped("binding_analysis module not available")

        try:
            from binding_analysis import analyze_interface
            result = analyze_interface(pdb_content)

            if result.get("status") == "error":
                return AnalysisResult.skipped(result.get("error", "unknown error"))

            return AnalysisResult.success({
                "dSASA": result.get("dSASA_int", 0),
                "contacts": result.get("contacts", 0),
                "interface_residues": result.get("nres_int", 0),
                "interface_hbonds": result.get("hbonds_int", 0),
                "packstat": result.get("packstat", 0),
            })
        except Exception as e:
            return AnalysisResult.skipped(f"interface analysis failed: {str(e)}")

    def _analyze_metal_coordination(
        self,
        pdb_content: str,
        metal_type: str,
        metal_chain: str,
        metal_resnum: int,
    ) -> AnalysisResult:
        """Analyze metal coordination geometry."""
        if not self._analysis_modules.get("metal_validation"):
            return AnalysisResult.skipped("metal_validation module not available")

        try:
            from metal_validation import validate_lanthanide_site

            result = validate_lanthanide_site(
                pdb_content=pdb_content,
                metal=metal_type,
                metal_chain=metal_chain,
                metal_resnum=metal_resnum,
            )

            return AnalysisResult.success({
                "coordination_number": result.get("coordination_number", 0),
                "geometry_type": result.get("geometry_type", "unknown"),
                "geometry_rmsd": result.get("geometry_rmsd", 0),
                "coordination_distance": result.get("mean_distance", 0),
                "hsab_compatible": result.get("hsab_compatible", False),
                "coordinating_residues": result.get("coordinating_residues", []),
            })
        except Exception as e:
            return AnalysisResult.skipped(f"metal coordination analysis failed: {str(e)}")

    def _analyze_ligand_binding(
        self,
        pdb_content: str,
        pdb_path: Optional[str],
        ligand_sdf: str,
    ) -> AnalysisResult:
        """Analyze protein-ligand binding with GNINA."""
        if not self._analysis_modules.get("gnina"):
            return AnalysisResult.skipped("GNINA binary not available")

        if not self._analysis_modules.get("binding_analysis"):
            return AnalysisResult.skipped("binding_analysis module not available")

        try:
            from binding_analysis import run_gnina_scoring

            result = run_gnina_scoring(
                receptor_pdb=pdb_content,
                ligand_sdf=ligand_sdf,
            )

            if result.get("status") == "error":
                return AnalysisResult.skipped(result.get("error", "GNINA failed"))

            return AnalysisResult.success({
                "gnina_affinity": result.get("affinity", 0),
                "gnina_cnn_score": result.get("cnn_score", 0),
                "gnina_cnn_affinity": result.get("cnn_affinity", 0),
            })
        except Exception as e:
            return AnalysisResult.skipped(f"ligand binding analysis failed: {str(e)}")

    def _analyze_topology(self, pdb_content: str) -> AnalysisResult:
        """Validate topology (ring closure, etc.)."""
        if not self._analysis_modules.get("topology_validation"):
            return AnalysisResult.skipped("topology_validation module not available")

        try:
            from topology_validation import validate_topology

            result = validate_topology(pdb_content)

            return AnalysisResult.success({
                "valid": result.get("valid", False),
                "issues": result.get("issues", []),
            })
        except Exception as e:
            return AnalysisResult.skipped(f"topology validation failed: {str(e)}")

    def _analyze_symmetry(self, pdb_content: str) -> AnalysisResult:
        """Analyze C2 symmetry for dimers."""
        # Basic symmetry analysis - calculate RMSD between chains
        try:
            chain_coords = {}
            for line in pdb_content.split("\n"):
                if line.startswith("ATOM") and " CA " in line:
                    chain_id = line[21]
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    if chain_id not in chain_coords:
                        chain_coords[chain_id] = []
                    chain_coords[chain_id].append((x, y, z))

            chains = list(chain_coords.keys())
            if len(chains) < 2:
                return AnalysisResult.not_applicable("need at least 2 chains for symmetry")

            # Simple symmetry score based on chain length similarity
            len_a = len(chain_coords[chains[0]])
            len_b = len(chain_coords[chains[1]])
            length_ratio = min(len_a, len_b) / max(len_a, len_b) if max(len_a, len_b) > 0 else 0

            return AnalysisResult.success({
                "c2_score": length_ratio,
                "chain_lengths": {chains[0]: len_a, chains[1]: len_b},
            })
        except Exception as e:
            return AnalysisResult.skipped(f"symmetry analysis failed: {str(e)}")
```

**Step 4: Run test to verify it passes**

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless
python -m pytest test_unified_analyzer.py -v
```

Expected: All tests PASS

**Step 5: Commit**

```bash
git add backend/serverless/unified_analyzer.py backend/serverless/test_unified_analyzer.py
git commit -m "feat: add UnifiedDesignAnalyzer core with analysis orchestration"
```

---

### Task 5: Create DesignHistoryManager

**Files:**
- Create: `backend/serverless/design_history.py`
- Test: `backend/serverless/test_design_history.py`

**Step 1: Write the failing test**

Create file `backend/serverless/test_design_history.py`:
```python
"""Tests for DesignHistoryManager."""
import pytest
import json
import os
import tempfile
import shutil
from design_history import DesignHistoryManager


class TestDesignHistoryManager:
    """Test suite for DesignHistoryManager."""

    @pytest.fixture
    def temp_history_dir(self):
        """Create temporary history directory."""
        temp_dir = tempfile.mkdtemp()
        # Create required subdirectories
        os.makedirs(os.path.join(temp_dir, "runs"))
        os.makedirs(os.path.join(temp_dir, "sessions"))
        os.makedirs(os.path.join(temp_dir, "lessons"))
        os.makedirs(os.path.join(temp_dir, "exports"))
        os.makedirs(os.path.join(temp_dir, "filter_presets"))
        # Create index.json
        with open(os.path.join(temp_dir, "index.json"), "w") as f:
            json.dump({"version": "1.0.0", "designs": []}, f)
        yield temp_dir
        shutil.rmtree(temp_dir)

    def test_init(self, temp_history_dir):
        manager = DesignHistoryManager(temp_history_dir)
        assert manager is not None
        assert manager.history_dir == temp_history_dir

    def test_start_session(self, temp_history_dir):
        manager = DesignHistoryManager(temp_history_dir)
        session = manager.start_session("test_exploration")
        assert session is not None
        assert "test_exploration" in session.session_id

    def test_save_run(self, temp_history_dir):
        manager = DesignHistoryManager(temp_history_dir)
        session = manager.start_session("test")

        run_id = manager.save_run(
            session=session,
            params={"temperature": 0.1},
            outputs={"pdb": "ATOM..."},
            metrics={"design_id": "test_001", "analyses": {}}
        )

        assert run_id is not None
        # Verify files were created
        run_dir = os.path.join(temp_history_dir, "runs", run_id)
        assert os.path.exists(run_dir)
        assert os.path.exists(os.path.join(run_dir, "input", "params.json"))
        assert os.path.exists(os.path.join(run_dir, "analysis", "metrics.json"))

    def test_load_index(self, temp_history_dir):
        manager = DesignHistoryManager(temp_history_dir)
        index = manager.load_index()
        assert "designs" in index
        assert isinstance(index["designs"], list)

    def test_get_session_stats(self, temp_history_dir):
        manager = DesignHistoryManager(temp_history_dir)
        session = manager.start_session("test")

        # Save a few runs
        for i in range(3):
            manager.save_run(
                session=session,
                params={},
                outputs={},
                metrics={
                    "design_id": f"test_{i}",
                    "analyses": {},
                    "filter_results": {"default": {"pass": i % 2 == 0}}
                }
            )

        stats = manager.get_session_stats(session)
        assert stats["total_designs"] == 3
        assert "acceptance_rate" in stats

    def test_export_metrics_csv(self, temp_history_dir):
        manager = DesignHistoryManager(temp_history_dir)
        session = manager.start_session("test")

        manager.save_run(
            session=session,
            params={},
            outputs={},
            metrics={
                "design_id": "test_001",
                "design_type": "metal_dimer",
                "analyses": {
                    "metal_coordination": {
                        "status": "success",
                        "metrics": {"coordination_distance": 2.3}
                    }
                }
            }
        )

        csv_path = manager.export_metrics_csv()
        assert os.path.exists(csv_path)
        with open(csv_path, "r") as f:
            content = f.read()
            assert "design_id" in content
            assert "test_001" in content
```

**Step 2: Run test to verify it fails**

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless
python -m pytest test_design_history.py -v
```

Expected: FAIL with "ModuleNotFoundError: No module named 'design_history'"

**Step 3: Write minimal implementation**

Create file `backend/serverless/design_history.py`:
```python
"""
Design History Manager

Handles persistent storage of design runs with session tracking,
filter evaluation, and CSV export.
"""
import os
import json
import csv
from datetime import datetime
from typing import Dict, Any, Optional, List
from dataclasses import dataclass, field


@dataclass
class DesignSession:
    """Represents a design exploration session."""
    session_id: str
    name: str
    started: str
    run_ids: List[str] = field(default_factory=list)


class DesignHistoryManager:
    """
    Manages persistent storage of design history.

    Usage:
        manager = DesignHistoryManager("experiments/design_history")
        session = manager.start_session("tb_exploration")
        run_id = manager.save_run(session, params, outputs, metrics)
    """

    def __init__(self, history_dir: str):
        """
        Initialize history manager.

        Args:
            history_dir: Path to design_history directory
        """
        self.history_dir = history_dir
        self._ensure_structure()

    def _ensure_structure(self):
        """Ensure directory structure exists."""
        subdirs = ["runs", "sessions", "lessons", "exports", "filter_presets"]
        for subdir in subdirs:
            os.makedirs(os.path.join(self.history_dir, subdir), exist_ok=True)

        # Ensure index.json exists
        index_path = os.path.join(self.history_dir, "index.json")
        if not os.path.exists(index_path):
            with open(index_path, "w") as f:
                json.dump({
                    "version": "1.0.0",
                    "created": datetime.now().isoformat(),
                    "designs": []
                }, f, indent=2)

    def load_index(self) -> Dict[str, Any]:
        """Load the design index."""
        index_path = os.path.join(self.history_dir, "index.json")
        with open(index_path, "r") as f:
            return json.load(f)

    def _save_index(self, index: Dict[str, Any]):
        """Save the design index."""
        index_path = os.path.join(self.history_dir, "index.json")
        with open(index_path, "w") as f:
            json.dump(index, f, indent=2)

    def start_session(self, name: str) -> DesignSession:
        """
        Start a new design exploration session.

        Args:
            name: Descriptive name for the session

        Returns:
            DesignSession object
        """
        timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
        session_id = f"{timestamp}_{name}"

        session = DesignSession(
            session_id=session_id,
            name=name,
            started=datetime.now().isoformat(),
        )

        # Save session file
        session_path = os.path.join(
            self.history_dir, "sessions", f"{session_id}.json"
        )
        with open(session_path, "w") as f:
            json.dump({
                "session_id": session.session_id,
                "name": session.name,
                "started": session.started,
                "run_ids": session.run_ids,
            }, f, indent=2)

        return session

    def save_run(
        self,
        session: DesignSession,
        params: Dict[str, Any],
        outputs: Dict[str, Any],
        metrics: Dict[str, Any],
    ) -> str:
        """
        Save a design run to history.

        Args:
            session: Current design session
            params: Design parameters
            outputs: Design outputs (PDB content, sequences, etc.)
            metrics: Analysis metrics from UnifiedDesignAnalyzer

        Returns:
            Run ID
        """
        # Generate run ID
        run_id = metrics.get("design_id", datetime.now().strftime("%Y-%m-%d_%H%M%S"))

        # Create run directory
        run_dir = os.path.join(self.history_dir, "runs", run_id)
        os.makedirs(run_dir, exist_ok=True)
        os.makedirs(os.path.join(run_dir, "input"), exist_ok=True)
        os.makedirs(os.path.join(run_dir, "output"), exist_ok=True)
        os.makedirs(os.path.join(run_dir, "analysis"), exist_ok=True)

        # Save params
        with open(os.path.join(run_dir, "input", "params.json"), "w") as f:
            json.dump(params, f, indent=2)

        # Save outputs
        for name, content in outputs.items():
            ext = "pdb" if "pdb" in name.lower() else "json"
            output_path = os.path.join(run_dir, "output", f"{name}.{ext}")
            if isinstance(content, str):
                with open(output_path, "w") as f:
                    f.write(content)
            else:
                with open(output_path, "w") as f:
                    json.dump(content, f, indent=2)

        # Save metrics
        with open(os.path.join(run_dir, "analysis", "metrics.json"), "w") as f:
            json.dump(metrics, f, indent=2)

        # Save metadata
        meta = {
            "run_id": run_id,
            "session_id": session.session_id,
            "timestamp": datetime.now().isoformat(),
            "design_type": metrics.get("design_type", "unknown"),
        }
        with open(os.path.join(run_dir, "meta.json"), "w") as f:
            json.dump(meta, f, indent=2)

        # Update session
        session.run_ids.append(run_id)
        self._update_session(session)

        # Update index
        index = self.load_index()
        index["designs"].append({
            "run_id": run_id,
            "session_id": session.session_id,
            "design_type": metrics.get("design_type", "unknown"),
            "timestamp": meta["timestamp"],
            "filter_pass": metrics.get("filter_results", {}).get("default", {}).get("pass"),
        })
        self._save_index(index)

        return run_id

    def _update_session(self, session: DesignSession):
        """Update session file."""
        session_path = os.path.join(
            self.history_dir, "sessions", f"{session.session_id}.json"
        )
        with open(session_path, "w") as f:
            json.dump({
                "session_id": session.session_id,
                "name": session.name,
                "started": session.started,
                "run_ids": session.run_ids,
            }, f, indent=2)

    def get_session_stats(self, session: DesignSession) -> Dict[str, Any]:
        """
        Get statistics for a session.

        Args:
            session: Design session

        Returns:
            Session statistics including acceptance rate
        """
        total = len(session.run_ids)
        passing = 0

        for run_id in session.run_ids:
            metrics_path = os.path.join(
                self.history_dir, "runs", run_id, "analysis", "metrics.json"
            )
            if os.path.exists(metrics_path):
                with open(metrics_path, "r") as f:
                    metrics = json.load(f)
                    if metrics.get("filter_results", {}).get("default", {}).get("pass"):
                        passing += 1

        return {
            "session_id": session.session_id,
            "total_designs": total,
            "passing_designs": passing,
            "acceptance_rate": passing / total if total > 0 else 0,
        }

    def export_metrics_csv(self, output_path: Optional[str] = None) -> str:
        """
        Export all metrics to CSV.

        Args:
            output_path: Optional custom output path

        Returns:
            Path to exported CSV
        """
        if output_path is None:
            output_path = os.path.join(
                self.history_dir, "exports", "all_metrics.csv"
            )

        index = self.load_index()
        rows = []

        for design in index["designs"]:
            run_id = design["run_id"]
            metrics_path = os.path.join(
                self.history_dir, "runs", run_id, "analysis", "metrics.json"
            )

            if os.path.exists(metrics_path):
                with open(metrics_path, "r") as f:
                    metrics = json.load(f)

                row = {
                    "design_id": run_id,
                    "design_type": metrics.get("design_type", ""),
                    "timestamp": metrics.get("timestamp", ""),
                }

                # Flatten analyses
                for analysis_name, analysis_data in metrics.get("analyses", {}).items():
                    if analysis_data.get("status") == "success":
                        for metric_name, value in analysis_data.get("metrics", {}).items():
                            if not isinstance(value, (list, dict)):
                                row[f"{analysis_name}_{metric_name}"] = value

                rows.append(row)

        if rows:
            # Get all column names
            all_columns = set()
            for row in rows:
                all_columns.update(row.keys())

            # Write CSV
            with open(output_path, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=sorted(all_columns))
                writer.writeheader()
                writer.writerows(rows)

        return output_path
```

**Step 4: Run test to verify it passes**

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless
python -m pytest test_design_history.py -v
```

Expected: All tests PASS

**Step 5: Commit**

```bash
git add backend/serverless/design_history.py backend/serverless/test_design_history.py
git commit -m "feat: add DesignHistoryManager with session tracking and CSV export"
```

---

## Phase 2: Complete Analysis Coverage

### Task 6: Add Filter Evaluation

**Files:**
- Create: `backend/serverless/filter_evaluator.py`
- Test: `backend/serverless/test_filter_evaluator.py`

**Step 1: Write the failing test**

Create file `backend/serverless/test_filter_evaluator.py`:
```python
"""Tests for FilterEvaluator."""
import pytest
import json
import os
import tempfile
from filter_evaluator import FilterEvaluator


class TestFilterEvaluator:
    """Test suite for FilterEvaluator."""

    @pytest.fixture
    def sample_preset(self):
        return {
            "name": "test",
            "filters": {
                "structure_confidence": {
                    "plddt": {"min": 0.8, "direction": "higher_better"}
                },
                "metal_coordination": {
                    "coordination_distance": {"max": 2.5, "direction": "lower_better"}
                }
            }
        }

    @pytest.fixture
    def temp_preset_dir(self, sample_preset):
        temp_dir = tempfile.mkdtemp()
        with open(os.path.join(temp_dir, "test.json"), "w") as f:
            json.dump(sample_preset, f)
        yield temp_dir
        import shutil
        shutil.rmtree(temp_dir)

    def test_init(self, temp_preset_dir):
        evaluator = FilterEvaluator(temp_preset_dir)
        assert evaluator is not None

    def test_load_preset(self, temp_preset_dir):
        evaluator = FilterEvaluator(temp_preset_dir)
        preset = evaluator.load_preset("test")
        assert preset["name"] == "test"
        assert "filters" in preset

    def test_evaluate_passing(self, temp_preset_dir):
        evaluator = FilterEvaluator(temp_preset_dir)
        metrics = {
            "analyses": {
                "structure_confidence": {
                    "status": "success",
                    "metrics": {"plddt": 0.9}
                },
                "metal_coordination": {
                    "status": "success",
                    "metrics": {"coordination_distance": 2.3}
                }
            }
        }
        result = evaluator.evaluate(metrics, "test")
        assert result["pass"] is True
        assert result["failed_filters"] == []

    def test_evaluate_failing(self, temp_preset_dir):
        evaluator = FilterEvaluator(temp_preset_dir)
        metrics = {
            "analyses": {
                "structure_confidence": {
                    "status": "success",
                    "metrics": {"plddt": 0.7}  # Below threshold
                },
                "metal_coordination": {
                    "status": "success",
                    "metrics": {"coordination_distance": 2.3}
                }
            }
        }
        result = evaluator.evaluate(metrics, "test")
        assert result["pass"] is False
        assert len(result["failed_filters"]) > 0

    def test_evaluate_skipped_analysis(self, temp_preset_dir):
        evaluator = FilterEvaluator(temp_preset_dir)
        metrics = {
            "analyses": {
                "structure_confidence": {
                    "status": "success",
                    "metrics": {"plddt": 0.9}
                },
                "metal_coordination": {
                    "status": "not_applicable",
                    "reason": "no metal"
                }
            }
        }
        result = evaluator.evaluate(metrics, "test")
        # Should pass since metal analysis is not applicable
        assert result["pass"] is True

    def test_evaluate_all_presets(self, temp_preset_dir):
        evaluator = FilterEvaluator(temp_preset_dir)
        metrics = {
            "analyses": {
                "structure_confidence": {
                    "status": "success",
                    "metrics": {"plddt": 0.9}
                }
            }
        }
        results = evaluator.evaluate_all_presets(metrics)
        assert "test" in results
```

**Step 2: Run test to verify it fails**

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless
python -m pytest test_filter_evaluator.py -v
```

Expected: FAIL with "ModuleNotFoundError: No module named 'filter_evaluator'"

**Step 3: Write minimal implementation**

Create file `backend/serverless/filter_evaluator.py`:
```python
"""
Filter Evaluator

Evaluates design metrics against filter presets (BindCraft-style).
"""
import os
import json
from typing import Dict, Any, List, Optional


class FilterEvaluator:
    """
    Evaluates metrics against filter presets.

    Usage:
        evaluator = FilterEvaluator("experiments/design_history/filter_presets")
        result = evaluator.evaluate(metrics, "default")
    """

    def __init__(self, presets_dir: str):
        """
        Initialize evaluator.

        Args:
            presets_dir: Path to filter_presets directory
        """
        self.presets_dir = presets_dir
        self._presets_cache: Dict[str, Dict] = {}

    def load_preset(self, name: str) -> Dict[str, Any]:
        """
        Load a filter preset by name.

        Args:
            name: Preset name (without .json extension)

        Returns:
            Preset dictionary
        """
        if name in self._presets_cache:
            return self._presets_cache[name]

        preset_path = os.path.join(self.presets_dir, f"{name}.json")
        with open(preset_path, "r") as f:
            preset = json.load(f)

        self._presets_cache[name] = preset
        return preset

    def list_presets(self) -> List[str]:
        """List available preset names."""
        presets = []
        for filename in os.listdir(self.presets_dir):
            if filename.endswith(".json"):
                presets.append(filename[:-5])
        return presets

    def evaluate(
        self,
        metrics: Dict[str, Any],
        preset_name: str,
    ) -> Dict[str, Any]:
        """
        Evaluate metrics against a filter preset.

        Args:
            metrics: Metrics from UnifiedDesignAnalyzer
            preset_name: Name of filter preset to use

        Returns:
            Evaluation result with pass/fail and failed filters
        """
        preset = self.load_preset(preset_name)
        failed_filters = []

        for analysis_name, filters in preset.get("filters", {}).items():
            analysis_data = metrics.get("analyses", {}).get(analysis_name, {})

            # Skip if analysis was not applicable or skipped
            if analysis_data.get("status") in ["not_applicable", "skipped"]:
                continue

            # Check each filter
            analysis_metrics = analysis_data.get("metrics", {})
            for metric_name, threshold in filters.items():
                value = analysis_metrics.get(metric_name)

                if value is None:
                    continue

                # Check threshold
                passed = self._check_threshold(value, threshold)
                if not passed:
                    failed_filters.append({
                        "analysis": analysis_name,
                        "metric": metric_name,
                        "value": value,
                        "threshold": threshold,
                    })

        return {
            "preset": preset_name,
            "pass": len(failed_filters) == 0,
            "failed_filters": failed_filters,
        }

    def _check_threshold(
        self,
        value: Any,
        threshold: Dict[str, Any],
    ) -> bool:
        """
        Check if value passes threshold.

        Args:
            value: Metric value
            threshold: Threshold definition

        Returns:
            True if passes, False if fails
        """
        if "min" in threshold:
            if value < threshold["min"]:
                return False

        if "max" in threshold:
            if value > threshold["max"]:
                return False

        if "equals" in threshold:
            if value != threshold["equals"]:
                return False

        return True

    def evaluate_all_presets(
        self,
        metrics: Dict[str, Any],
    ) -> Dict[str, Dict[str, Any]]:
        """
        Evaluate metrics against all available presets.

        Args:
            metrics: Metrics from UnifiedDesignAnalyzer

        Returns:
            Dictionary of preset_name -> evaluation result
        """
        results = {}
        for preset_name in self.list_presets():
            try:
                results[preset_name] = self.evaluate(metrics, preset_name)
            except Exception as e:
                results[preset_name] = {
                    "preset": preset_name,
                    "pass": False,
                    "error": str(e),
                }
        return results
```

**Step 4: Run test to verify it passes**

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless
python -m pytest test_filter_evaluator.py -v
```

Expected: All tests PASS

**Step 5: Commit**

```bash
git add backend/serverless/filter_evaluator.py backend/serverless/test_filter_evaluator.py
git commit -m "feat: add FilterEvaluator for BindCraft-style threshold evaluation"
```

---

### Task 7: Add Lesson Trigger Detection

**Files:**
- Create: `backend/serverless/lesson_detector.py`
- Test: `backend/serverless/test_lesson_detector.py`

**Step 1: Write the failing test**

Create file `backend/serverless/test_lesson_detector.py`:
```python
"""Tests for LessonDetector."""
import pytest
from lesson_detector import LessonDetector, LessonTrigger


class TestLessonDetector:
    """Test suite for LessonDetector."""

    def test_init(self):
        detector = LessonDetector()
        assert detector is not None

    def test_detect_failure_pattern(self):
        detector = LessonDetector()
        history = [
            {"outcome": "failure", "params": {"loop_length": 6}, "design_type": "metal_dimer"},
            {"outcome": "failure", "params": {"loop_length": 6}, "design_type": "metal_dimer"},
            {"outcome": "failure", "params": {"loop_length": 6}, "design_type": "metal_dimer"},
        ]
        trigger = detector.check_triggers(
            new_result=history[-1],
            history=history
        )
        assert trigger is not None
        assert trigger.trigger_type == "failure_pattern"

    def test_detect_breakthrough(self):
        detector = LessonDetector()
        history = [
            {"outcome": "failure", "metrics": {"gnina_affinity": -4}},
            {"outcome": "failure", "metrics": {"gnina_affinity": -5}},
        ]
        new_result = {"outcome": "success", "metrics": {"gnina_affinity": -8}}

        trigger = detector.check_triggers(new_result, history)
        assert trigger is not None
        assert trigger.trigger_type == "breakthrough"

    def test_detect_improvement(self):
        detector = LessonDetector()
        history = [
            {"metrics": {"coordination_distance": 3.0}},
            {"metrics": {"coordination_distance": 2.8}},
        ]
        new_result = {"metrics": {"coordination_distance": 2.3}}

        trigger = detector.check_triggers(new_result, history)
        assert trigger is not None
        assert trigger.trigger_type == "improvement"

    def test_no_trigger_for_normal_result(self):
        detector = LessonDetector()
        history = [
            {"outcome": "success", "metrics": {"gnina_affinity": -6}},
            {"outcome": "failure", "metrics": {"gnina_affinity": -4}},
        ]
        new_result = {"outcome": "success", "metrics": {"gnina_affinity": -6.5}}

        trigger = detector.check_triggers(new_result, history)
        # Small improvement, not significant enough
        assert trigger is None or trigger.trigger_type != "breakthrough"
```

**Step 2: Run test to verify it fails**

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless
python -m pytest test_lesson_detector.py -v
```

Expected: FAIL with "ModuleNotFoundError: No module named 'lesson_detector'"

**Step 3: Write minimal implementation**

Create file `backend/serverless/lesson_detector.py`:
```python
"""
Lesson Detector

Detects significant events (failure patterns, breakthroughs, improvements)
that warrant synthesizing lessons learned.
"""
from dataclasses import dataclass
from typing import Dict, Any, List, Optional
from collections import Counter


@dataclass
class LessonTrigger:
    """Represents a detected lesson trigger."""
    trigger_type: str  # "failure_pattern", "breakthrough", "improvement"
    description: str
    relevant_designs: List[str]
    metrics_involved: List[str]


class LessonDetector:
    """
    Detects significant events that should trigger lesson synthesis.

    Triggers:
    - failure_pattern: 3+ consecutive failures with similar characteristics
    - breakthrough: New best achieved on key metrics
    - improvement: Significant improvement (>15%) from previous best
    """

    def __init__(
        self,
        failure_threshold: int = 3,
        improvement_threshold: float = 0.15,
    ):
        """
        Initialize detector.

        Args:
            failure_threshold: Number of similar failures to trigger pattern detection
            improvement_threshold: Relative improvement to trigger (0.15 = 15%)
        """
        self.failure_threshold = failure_threshold
        self.improvement_threshold = improvement_threshold

        # Key metrics to track for breakthroughs/improvements
        self.tracked_metrics = [
            ("gnina_affinity", "lower_better"),
            ("coordination_distance", "lower_better"),
            ("plddt", "higher_better"),
            ("shape_complementarity", "higher_better"),
            ("c2_score", "higher_better"),
        ]

    def check_triggers(
        self,
        new_result: Dict[str, Any],
        history: List[Dict[str, Any]],
    ) -> Optional[LessonTrigger]:
        """
        Check if new result triggers lesson synthesis.

        Args:
            new_result: Latest design result
            history: List of previous results

        Returns:
            LessonTrigger if significant event detected, None otherwise
        """
        # Check failure pattern
        failure_trigger = self._check_failure_pattern(new_result, history)
        if failure_trigger:
            return failure_trigger

        # Check breakthrough
        breakthrough_trigger = self._check_breakthrough(new_result, history)
        if breakthrough_trigger:
            return breakthrough_trigger

        # Check improvement
        improvement_trigger = self._check_improvement(new_result, history)
        if improvement_trigger:
            return improvement_trigger

        return None

    def _check_failure_pattern(
        self,
        new_result: Dict[str, Any],
        history: List[Dict[str, Any]],
    ) -> Optional[LessonTrigger]:
        """Check for repeated failure patterns."""
        if new_result.get("outcome") != "failure":
            return None

        # Get recent failures
        recent_failures = [
            r for r in history[-10:]
            if r.get("outcome") == "failure"
        ]

        if len(recent_failures) < self.failure_threshold - 1:
            return None

        # Check for common parameters in failures
        param_patterns = self._find_common_params(recent_failures + [new_result])

        if param_patterns:
            return LessonTrigger(
                trigger_type="failure_pattern",
                description=f"Repeated failures with common parameters: {param_patterns}",
                relevant_designs=[r.get("design_id", "unknown") for r in recent_failures],
                metrics_involved=list(param_patterns.keys()),
            )

        # Check for common design type
        design_types = [r.get("design_type") for r in recent_failures + [new_result]]
        type_counts = Counter(design_types)
        most_common = type_counts.most_common(1)

        if most_common and most_common[0][1] >= self.failure_threshold:
            return LessonTrigger(
                trigger_type="failure_pattern",
                description=f"Repeated failures for design type: {most_common[0][0]}",
                relevant_designs=[r.get("design_id", "unknown") for r in recent_failures],
                metrics_involved=["design_type"],
            )

        return None

    def _find_common_params(
        self,
        results: List[Dict[str, Any]],
    ) -> Dict[str, Any]:
        """Find common parameters across results."""
        if not results:
            return {}

        common = {}
        first_params = results[0].get("params", {})

        for key, value in first_params.items():
            if all(r.get("params", {}).get(key) == value for r in results[1:]):
                common[key] = value

        return common

    def _check_breakthrough(
        self,
        new_result: Dict[str, Any],
        history: List[Dict[str, Any]],
    ) -> Optional[LessonTrigger]:
        """Check if result achieves new best on key metrics."""
        new_metrics = new_result.get("metrics", {})
        breakthrough_metrics = []

        for metric_name, direction in self.tracked_metrics:
            new_value = self._get_nested_metric(new_metrics, metric_name)
            if new_value is None:
                continue

            # Get best previous value
            best_prev = None
            for result in history:
                prev_value = self._get_nested_metric(
                    result.get("metrics", {}), metric_name
                )
                if prev_value is not None:
                    if best_prev is None:
                        best_prev = prev_value
                    elif direction == "lower_better":
                        best_prev = min(best_prev, prev_value)
                    else:
                        best_prev = max(best_prev, prev_value)

            if best_prev is None:
                continue

            # Check if new value is significantly better
            if direction == "lower_better":
                improvement = (best_prev - new_value) / abs(best_prev) if best_prev != 0 else 0
                is_better = new_value < best_prev
            else:
                improvement = (new_value - best_prev) / abs(best_prev) if best_prev != 0 else 0
                is_better = new_value > best_prev

            if is_better and improvement > self.improvement_threshold * 2:
                breakthrough_metrics.append(metric_name)

        if breakthrough_metrics:
            return LessonTrigger(
                trigger_type="breakthrough",
                description=f"New best achieved: {', '.join(breakthrough_metrics)}",
                relevant_designs=[new_result.get("design_id", "unknown")],
                metrics_involved=breakthrough_metrics,
            )

        return None

    def _check_improvement(
        self,
        new_result: Dict[str, Any],
        history: List[Dict[str, Any]],
    ) -> Optional[LessonTrigger]:
        """Check for significant improvement from recent results."""
        if len(history) < 2:
            return None

        new_metrics = new_result.get("metrics", {})
        improved_metrics = []

        for metric_name, direction in self.tracked_metrics:
            new_value = self._get_nested_metric(new_metrics, metric_name)
            if new_value is None:
                continue

            # Get recent average
            recent_values = []
            for result in history[-5:]:
                prev_value = self._get_nested_metric(
                    result.get("metrics", {}), metric_name
                )
                if prev_value is not None:
                    recent_values.append(prev_value)

            if not recent_values:
                continue

            avg_recent = sum(recent_values) / len(recent_values)

            # Check improvement
            if direction == "lower_better":
                improvement = (avg_recent - new_value) / abs(avg_recent) if avg_recent != 0 else 0
            else:
                improvement = (new_value - avg_recent) / abs(avg_recent) if avg_recent != 0 else 0

            if improvement > self.improvement_threshold:
                improved_metrics.append(metric_name)

        if improved_metrics:
            return LessonTrigger(
                trigger_type="improvement",
                description=f"Significant improvement in: {', '.join(improved_metrics)}",
                relevant_designs=[new_result.get("design_id", "unknown")],
                metrics_involved=improved_metrics,
            )

        return None

    def _get_nested_metric(
        self,
        metrics: Dict[str, Any],
        metric_name: str,
    ) -> Optional[float]:
        """Get metric value, searching in nested analyses."""
        # Direct lookup
        if metric_name in metrics:
            return metrics[metric_name]

        # Search in analyses
        for analysis_name, analysis_data in metrics.items():
            if isinstance(analysis_data, dict):
                nested_metrics = analysis_data.get("metrics", {})
                if metric_name in nested_metrics:
                    return nested_metrics[metric_name]

        return None
```

**Step 4: Run test to verify it passes**

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless
python -m pytest test_lesson_detector.py -v
```

Expected: All tests PASS

**Step 5: Commit**

```bash
git add backend/serverless/lesson_detector.py backend/serverless/test_lesson_detector.py
git commit -m "feat: add LessonDetector for failure patterns and breakthroughs"
```

---

### Task 8: Create CLI Entry Point

**Files:**
- Create: `backend/serverless/analyze_design_cli.py`

**Step 1: Create CLI script**

Create file `backend/serverless/analyze_design_cli.py`:
```python
#!/usr/bin/env python3
"""
CLI for Unified Design Analyzer

Usage:
    python -m analyze_design_cli output.pdb --params params.json --metal TB
    python -m analyze_design_cli output.pdb --ligand mol.sdf --session my_session
"""
import argparse
import json
import os
import sys
from pathlib import Path

from unified_analyzer import UnifiedDesignAnalyzer
from design_history import DesignHistoryManager
from filter_evaluator import FilterEvaluator
from lesson_detector import LessonDetector


def main():
    parser = argparse.ArgumentParser(
        description="Analyze protein design and save to history"
    )
    parser.add_argument("pdb_file", help="Path to PDB file to analyze")
    parser.add_argument("--params", help="Path to design parameters JSON")
    parser.add_argument("--ligand", help="Path to ligand SDF file")
    parser.add_argument("--metal", help="Metal type code (e.g., TB, ZN)")
    parser.add_argument("--metal-chain", default="L", help="Metal chain ID")
    parser.add_argument("--metal-resnum", type=int, default=1, help="Metal residue number")
    parser.add_argument("--session", help="Session name (creates new if not exists)")
    parser.add_argument(
        "--history-dir",
        default="experiments/design_history",
        help="Path to design history directory"
    )
    parser.add_argument("--no-save", action="store_true", help="Don't save to history")
    parser.add_argument("--json", action="store_true", help="Output JSON only")

    args = parser.parse_args()

    # Read PDB file
    with open(args.pdb_file, "r") as f:
        pdb_content = f.read()

    # Read params if provided
    params = {}
    if args.params:
        with open(args.params, "r") as f:
            params = json.load(f)

    # Read ligand if provided
    ligand_sdf = None
    if args.ligand:
        with open(args.ligand, "r") as f:
            ligand_sdf = f.read()

    # Initialize analyzer
    analyzer = UnifiedDesignAnalyzer()

    # Run analysis
    metrics = analyzer.analyze(
        pdb_content=pdb_content,
        design_params=params,
        pdb_path=args.pdb_file,
        ligand_sdf=ligand_sdf,
        metal_type=args.metal,
        metal_chain=args.metal_chain,
        metal_resnum=args.metal_resnum,
    )

    # Evaluate against filter presets
    presets_dir = os.path.join(args.history_dir, "filter_presets")
    if os.path.exists(presets_dir):
        evaluator = FilterEvaluator(presets_dir)
        metrics["filter_results"] = evaluator.evaluate_all_presets(metrics)

    # Save to history if requested
    if not args.no_save:
        history = DesignHistoryManager(args.history_dir)

        # Get or create session
        session_name = args.session or "cli_session"
        session = history.start_session(session_name)

        # Save run
        run_id = history.save_run(
            session=session,
            params=params,
            outputs={"pdb": pdb_content},
            metrics=metrics,
        )

        # Check for lesson triggers
        index = history.load_index()
        detector = LessonDetector()

        # Build history for detector
        design_history = []
        for design in index["designs"][-20:]:
            metrics_path = os.path.join(
                args.history_dir, "runs", design["run_id"],
                "analysis", "metrics.json"
            )
            if os.path.exists(metrics_path):
                with open(metrics_path) as f:
                    design_history.append(json.load(f))

        trigger = detector.check_triggers(metrics, design_history)

        if not args.json:
            print(f"\nSaved to: {run_id}")
            if trigger:
                print(f"\n*** LESSON TRIGGER: {trigger.trigger_type} ***")
                print(f"    {trigger.description}")

    # Output
    if args.json:
        print(json.dumps(metrics, indent=2))
    else:
        print(f"\nDesign Analysis: {metrics['design_id']}")
        print(f"Type: {metrics['design_type']}")
        print("\nAnalyses:")
        for name, data in metrics["analyses"].items():
            status = data["status"]
            if status == "success":
                print(f"  {name}: {status}")
                for k, v in data.get("metrics", {}).items():
                    if not isinstance(v, (list, dict)):
                        print(f"    - {k}: {v}")
            else:
                print(f"  {name}: {status}")
                if data.get("reason"):
                    print(f"    ({data['reason']})")

        if "filter_results" in metrics:
            print("\nFilter Results:")
            for preset, result in metrics["filter_results"].items():
                status = "PASS" if result.get("pass") else "FAIL"
                print(f"  {preset}: {status}")


if __name__ == "__main__":
    main()
```

**Step 2: Test CLI manually**

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless
python analyze_design_cli.py --help
```

Expected: Help message displays

**Step 3: Commit**

```bash
git add backend/serverless/analyze_design_cli.py
git commit -m "feat: add CLI entry point for design analysis"
```

---

### Task 9: Update CLAUDE.md with Design History Reference

**Files:**
- Modify: `G:\Github_local_repo\Banta_Lab_RFdiffusion\CLAUDE.md`

**Step 1: Add design history section to CLAUDE.md**

Add the following section to the project's CLAUDE.md:

```markdown
## Design History & Lessons

Local design experiments are tracked in `experiments/design_history/`.

**Current lessons learned:** See `experiments/design_history/lessons/current_summary.md`

When running local design tests:
1. Use `UnifiedDesignAnalyzer` for comprehensive metrics
2. Results auto-save to design_history with full provenance
3. Lessons auto-update when significant patterns detected

**Quick reference:**
- Recent runs: `experiments/design_history/index.json`
- Filter presets: `experiments/design_history/filter_presets/`
- Lesson triggers: failure patterns (3+ similar), breakthroughs, significant improvements

**CLI usage:**
```bash
cd backend/serverless
python analyze_design_cli.py output.pdb --metal TB --session my_exploration
```
```

**Step 2: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: add design history reference to CLAUDE.md"
```

---

### Task 10: Create Initial Lessons Summary

**Files:**
- Create: `experiments/design_history/lessons/current_summary.md`

**Step 1: Create initial summary file**

Create file `experiments/design_history/lessons/current_summary.md`:
```markdown
# Design Lessons (Auto-updated)

> This file is automatically refreshed when significant patterns are detected.
> Last updated: 2026-01-18

## Metal Coordination

*No lessons recorded yet. Run local designs to build knowledge base.*

## Ligand Interface

*No lessons recorded yet.*

## Failure Patterns to Avoid

*No failure patterns detected yet.*

---

## How Lessons Are Captured

Lessons are automatically synthesized when:
1. **Failure pattern detected**: 3+ consecutive failures with similar parameters
2. **Breakthrough success**: New best achieved on key metrics
3. **Meaningful improvement**: Significant improvement (>15%) from previous best

To contribute to lessons, run local designs using:
```bash
python analyze_design_cli.py your_design.pdb --metal TB --session exploration_name
```
```

**Step 2: Commit**

```bash
git add experiments/design_history/lessons/current_summary.md
git commit -m "docs: add initial lessons summary template"
```

---

## Phase 3: Integration Tests

### Task 11: Create Integration Test

**Files:**
- Create: `backend/serverless/test_unified_integration.py`

**Step 1: Write integration test**

Create file `backend/serverless/test_unified_integration.py`:
```python
"""Integration tests for the unified design analysis system."""
import pytest
import os
import json
import tempfile
import shutil

from unified_analyzer import UnifiedDesignAnalyzer
from design_history import DesignHistoryManager
from filter_evaluator import FilterEvaluator
from lesson_detector import LessonDetector


# More realistic test PDB
REALISTIC_DIMER_PDB = """HEADER    TEST DIMER
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 85.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 85.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 85.00           C
ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00 85.00           O
ATOM      5  N   GLU A   2       3.300   1.500   0.000  1.00 90.00           N
ATOM      6  CA  GLU A   2       3.900   2.800   0.000  1.00 90.00           C
ATOM      7  C   GLU A   2       5.400   2.800   0.000  1.00 90.00           C
ATOM      8  O   GLU A   2       6.000   1.700   0.000  1.00 90.00           O
ATOM      9  N   ALA B   1      20.000   0.000   0.000  1.00 85.00           N
ATOM     10  CA  ALA B   1      21.458   0.000   0.000  1.00 85.00           C
ATOM     11  C   ALA B   1      22.009   1.420   0.000  1.00 85.00           C
ATOM     12  O   ALA B   1      21.251   2.390   0.000  1.00 85.00           O
ATOM     13  N   ASP B   2      23.300   1.500   0.000  1.00 88.00           N
ATOM     14  CA  ASP B   2      23.900   2.800   0.000  1.00 88.00           C
ATOM     15  C   ASP B   2      25.400   2.800   0.000  1.00 88.00           C
ATOM     16  O   ASP B   2      26.000   1.700   0.000  1.00 88.00           O
END
"""


class TestUnifiedIntegration:
    """Integration tests for the complete analysis pipeline."""

    @pytest.fixture
    def temp_history_dir(self):
        """Create complete temporary history setup."""
        temp_dir = tempfile.mkdtemp()

        # Create directory structure
        os.makedirs(os.path.join(temp_dir, "runs"))
        os.makedirs(os.path.join(temp_dir, "sessions"))
        os.makedirs(os.path.join(temp_dir, "lessons"))
        os.makedirs(os.path.join(temp_dir, "exports"))
        os.makedirs(os.path.join(temp_dir, "filter_presets"))

        # Create index
        with open(os.path.join(temp_dir, "index.json"), "w") as f:
            json.dump({"version": "1.0.0", "designs": []}, f)

        # Create default filter preset
        default_preset = {
            "name": "default",
            "filters": {
                "structure_confidence": {
                    "plddt": {"min": 0.8, "direction": "higher_better"}
                }
            }
        }
        with open(os.path.join(temp_dir, "filter_presets", "default.json"), "w") as f:
            json.dump(default_preset, f)

        yield temp_dir
        shutil.rmtree(temp_dir)

    def test_full_pipeline(self, temp_history_dir):
        """Test complete analysis pipeline."""
        # Initialize components
        analyzer = UnifiedDesignAnalyzer()
        history = DesignHistoryManager(temp_history_dir)
        evaluator = FilterEvaluator(os.path.join(temp_history_dir, "filter_presets"))

        # Start session
        session = history.start_session("integration_test")

        # Analyze design
        metrics = analyzer.analyze(
            pdb_content=REALISTIC_DIMER_PDB,
            design_params={"test": True},
        )

        # Evaluate filters
        metrics["filter_results"] = evaluator.evaluate_all_presets(metrics)

        # Save to history
        run_id = history.save_run(
            session=session,
            params={"test": True},
            outputs={"pdb": REALISTIC_DIMER_PDB},
            metrics=metrics,
        )

        # Verify
        assert run_id is not None
        assert os.path.exists(os.path.join(temp_history_dir, "runs", run_id))

        # Check index was updated
        index = history.load_index()
        assert len(index["designs"]) == 1
        assert index["designs"][0]["run_id"] == run_id

        # Check session stats
        stats = history.get_session_stats(session)
        assert stats["total_designs"] == 1

    def test_multiple_designs_session(self, temp_history_dir):
        """Test multiple designs in a session."""
        analyzer = UnifiedDesignAnalyzer()
        history = DesignHistoryManager(temp_history_dir)

        session = history.start_session("multi_test")

        # Run multiple designs
        for i in range(5):
            metrics = analyzer.analyze(
                pdb_content=REALISTIC_DIMER_PDB,
                design_params={"iteration": i},
            )
            history.save_run(
                session=session,
                params={"iteration": i},
                outputs={},
                metrics=metrics,
            )

        # Check stats
        stats = history.get_session_stats(session)
        assert stats["total_designs"] == 5

        # Check CSV export
        csv_path = history.export_metrics_csv()
        assert os.path.exists(csv_path)

    def test_lesson_detection_integration(self, temp_history_dir):
        """Test lesson detection with real history."""
        history = DesignHistoryManager(temp_history_dir)
        detector = LessonDetector()

        session = history.start_session("lesson_test")

        # Simulate failures
        history_list = []
        for i in range(4):
            metrics = {
                "design_id": f"fail_{i}",
                "design_type": "metal_dimer",
                "outcome": "failure",
                "params": {"loop_length": 6},
                "metrics": {"coordination_distance": 3.5},
            }
            history.save_run(session, {"loop_length": 6}, {}, metrics)
            history_list.append(metrics)

        # Check for failure pattern
        trigger = detector.check_triggers(history_list[-1], history_list[:-1])
        assert trigger is not None
        assert trigger.trigger_type == "failure_pattern"
```

**Step 2: Run integration tests**

```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless
python -m pytest test_unified_integration.py -v
```

Expected: All tests PASS

**Step 3: Commit**

```bash
git add backend/serverless/test_unified_integration.py
git commit -m "test: add integration tests for unified analysis pipeline"
```

---

## Summary

### Files Created
1. `experiments/design_history/` - Directory structure
2. `experiments/design_history/filter_presets/*.json` - Filter presets
3. `backend/serverless/analysis_types.py` - Type definitions
4. `backend/serverless/unified_analyzer.py` - Core analyzer
5. `backend/serverless/design_history.py` - History manager
6. `backend/serverless/filter_evaluator.py` - Filter evaluation
7. `backend/serverless/lesson_detector.py` - Lesson detection
8. `backend/serverless/analyze_design_cli.py` - CLI entry point
9. `experiments/design_history/lessons/current_summary.md` - Lessons template

### Files Modified
1. `CLAUDE.md` - Added design history reference

### Test Files
1. `backend/serverless/test_analysis_types.py`
2. `backend/serverless/test_unified_analyzer.py`
3. `backend/serverless/test_design_history.py`
4. `backend/serverless/test_filter_evaluator.py`
5. `backend/serverless/test_lesson_detector.py`
6. `backend/serverless/test_unified_integration.py`

### Run All Tests
```bash
cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless
python -m pytest test_analysis_types.py test_unified_analyzer.py test_design_history.py test_filter_evaluator.py test_lesson_detector.py test_unified_integration.py -v
```

---

*Plan created: 2026-01-18*
*Ready for execution*
