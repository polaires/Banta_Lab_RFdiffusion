"""
Tests for the modular pipeline architecture.

Tests pipeline_types, inference_backend, pipeline_modules,
workflow_runner, and iteration_strategies WITHOUT requiring
GPU or Foundry imports.
"""

import sys
import os
import json
import pytest
from unittest.mock import MagicMock, patch
from pathlib import Path

# Add parent dir to path so we can import backend modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


# ============================================================
# Test pipeline_types.py
# ============================================================

class TestPipelineTypes:
    """Test StepContext, result dataclasses, and PipelineStep protocol."""

    def test_step_context_defaults(self):
        from pipeline_types import StepContext
        ctx = StepContext()
        assert ctx.backbone_pdbs == []
        assert ctx.sequences == []
        assert ctx.predictions == []
        assert ctx.params == {}
        assert ctx.metadata == {}
        assert ctx.backend is None

    def test_step_context_get_param(self):
        from pipeline_types import StepContext
        ctx = StepContext(params={"num_designs": 10, "metal": "TB"})
        assert ctx.get_param("num_designs") == 10
        assert ctx.get_param("metal") == "TB"
        assert ctx.get_param("missing") is None
        assert ctx.get_param("missing", 42) == 42

    def test_step_context_to_dict(self):
        from pipeline_types import StepContext
        ctx = StepContext(
            params={"query": "test"},
            backbone_pdbs=["PDB1", "PDB2"],
        )
        d = ctx.to_dict()
        assert d["params"] == {"query": "test"}
        assert d["num_backbones"] == 2
        assert d["num_sequences"] == 0
        assert len(d["backbone_pdbs"]) == 2

    def test_step_context_to_dict_truncates_large_pdbs(self):
        from pipeline_types import StepContext
        ctx = StepContext(backbone_pdbs=[f"PDB_{i}" for i in range(20)])
        d = ctx.to_dict()
        assert len(d["backbone_pdbs"]) == 10
        assert d["backbone_pdbs_truncated"] is True

    def test_sequence_result(self):
        from pipeline_types import SequenceResult
        sr = SequenceResult(sequence="MEVTI", score=-1.5, backbone_index=0)
        assert sr.sequence == "MEVTI"
        assert sr.score == -1.5
        assert sr.backbone_index == 0

    def test_prediction_result(self):
        from pipeline_types import PredictionResult
        pr = PredictionResult(
            sequence="MEVTI",
            plddt=0.85,
            ptm=0.9,
            passed=True,
        )
        assert pr.plddt == 0.85
        assert pr.ptm == 0.9
        assert pr.passed is True

    def test_design_candidate(self):
        from pipeline_types import DesignCandidate
        dc = DesignCandidate(backbone_pdb="PDB", sequence="MEVTI")
        assert dc.backbone_pdb == "PDB"
        assert dc.rank == 0

    def test_pipeline_step_protocol(self):
        """Verify that a class with the right attributes satisfies the protocol."""
        from pipeline_types import PipelineStep, StepContext

        class FakeStep:
            name = "fake"
            description = "test"
            input_keys = []
            output_keys = []
            optional_keys = []
            def run(self, context: StepContext) -> StepContext:
                return context

        step = FakeStep()
        assert isinstance(step, PipelineStep)


# ============================================================
# Test inference_backend.py
# ============================================================

class TestInferenceBackend:
    """Test InProcessBackend, HTTPBackend, and create_backend factory."""

    def test_filter_params_strips_task(self):
        from inference_backend import _filter_params
        def dummy(a, b, c=None):
            pass
        result = _filter_params({"task": "rfd3", "a": 1, "b": 2, "c": 3}, dummy)
        assert "task" not in result
        assert result == {"a": 1, "b": 2, "c": 3}

    def test_filter_params_drops_unknown_keys(self):
        from inference_backend import _filter_params
        def dummy(a, b):
            pass
        result = _filter_params({"a": 1, "b": 2, "unknown": 99}, dummy)
        assert result == {"a": 1, "b": 2}

    def test_filter_params_allows_all_with_kwargs(self):
        from inference_backend import _filter_params
        def dummy(a, **kwargs):
            pass
        result = _filter_params({"a": 1, "b": 2, "c": 3}, dummy)
        assert result == {"a": 1, "b": 2, "c": 3}

    def test_create_backend_unknown_mode(self):
        from inference_backend import create_backend
        with pytest.raises(ValueError, match="Unknown backend mode"):
            create_backend("invalid_mode")

    def test_http_backend_init(self):
        from inference_backend import HTTPBackend
        backend = HTTPBackend(api_url="http://test:8000", timeout=30)
        assert backend.api_url == "http://test:8000"
        assert backend.timeout == 30

    def test_http_backend_raises_on_failure(self):
        from inference_backend import HTTPBackend
        backend = HTTPBackend()
        mock_call = MagicMock(return_value={"status": "failed", "error": "boom"})
        with patch.dict("sys.modules", {"api_client": MagicMock(call_api=mock_call, API_URL="http://test:8000", DEFAULT_TIMEOUT=600)}):
            # Re-import to pick up mock
            backend2 = HTTPBackend.__new__(HTTPBackend)
            backend2.api_url = "http://test:8000"
            backend2.timeout = 600
            # Patch the lazy import inside run_rfd3
            with patch("api_client.call_api", mock_call):
                try:
                    backend2.run_rfd3({"contig": "80"})
                    assert False, "Should have raised"
                except RuntimeError as e:
                    assert "boom" in str(e)


# ============================================================
# Test filter_evaluator.py FILTER_PRESETS
# ============================================================

class TestFilterPresets:
    """Test centralized FILTER_PRESETS format."""

    def test_filter_presets_exist(self):
        from filter_evaluator import FILTER_PRESETS
        assert "scout_relaxed" in FILTER_PRESETS
        assert "metal_binding" in FILTER_PRESETS
        assert "production_strict" in FILTER_PRESETS
        assert len(FILTER_PRESETS) == 8

    def test_filter_presets_use_min_max(self):
        from filter_evaluator import FILTER_PRESETS
        for preset_name, thresholds in FILTER_PRESETS.items():
            for metric, threshold in thresholds.items():
                assert isinstance(threshold, dict), (
                    f"Preset {preset_name!r}, metric {metric!r}: "
                    f"expected dict with min/max, got {type(threshold).__name__}"
                )
                assert "min" in threshold or "max" in threshold, (
                    f"Preset {preset_name!r}, metric {metric!r}: "
                    f"dict must have 'min' or 'max' key"
                )

    def test_scout_relaxed_thresholds(self):
        from filter_evaluator import FILTER_PRESETS
        scout = FILTER_PRESETS["scout_relaxed"]
        assert scout["ptm"]["min"] == 0.6
        assert scout["plddt"]["min"] == 0.65

    def test_metal_binding_has_max_for_rmsd(self):
        from filter_evaluator import FILTER_PRESETS
        metal = FILTER_PRESETS["metal_binding"]
        assert "max" in metal["geometry_rmsd"]
        assert "min" in metal["plddt"]


# ============================================================
# Test pipeline_modules registry
# ============================================================

class TestModuleRegistry:
    """Test MODULE_REGISTRY and module imports."""

    def test_registry_has_all_modules(self):
        from pipeline_modules import MODULE_REGISTRY
        expected = {
            "intent_parser", "design_configurator", "scaffolder",
            "backbone_generator", "sequence_designer", "structure_predictor",
            "analyzer", "reporter",
        }
        assert set(MODULE_REGISTRY.keys()) == expected

    def test_modules_have_protocol_attributes(self):
        from pipeline_modules import MODULE_REGISTRY
        for name, cls in MODULE_REGISTRY.items():
            instance = cls()
            assert hasattr(instance, "name"), f"{name} missing 'name'"
            assert hasattr(instance, "description"), f"{name} missing 'description'"
            assert hasattr(instance, "input_keys"), f"{name} missing 'input_keys'"
            assert hasattr(instance, "output_keys"), f"{name} missing 'output_keys'"
            assert hasattr(instance, "run"), f"{name} missing 'run'"
            assert instance.name == name, f"{name}: name={instance.name!r}"


# ============================================================
# Test BackboneGenerator module
# ============================================================

class TestBackboneGenerator:
    """Test BackboneGenerator with mock backend."""

    def _make_mock_backend(self, designs):
        backend = MagicMock()
        backend.run_rfd3.return_value = {
            "designs": [{"pdb_content": d} for d in designs]
        }
        return backend

    def test_basic_generation(self):
        from pipeline_types import StepContext
        from pipeline_modules.backbone_generator import BackboneGenerator

        backend = self._make_mock_backend(["PDB_A", "PDB_B"])
        ctx = StepContext(
            backend=backend,
            rfd3_config={"contig": "80-120", "num_designs": 2},
        )

        gen = BackboneGenerator()
        result = gen.run(ctx)

        assert len(result.backbone_pdbs) == 2
        assert result.backbone_pdbs[0] == "PDB_A"
        backend.run_rfd3.assert_called_once()

    def test_hetatm_merge(self):
        from pipeline_modules.backbone_generator import _merge_hetatm

        design = "ATOM      1  N   ALA A   1\nEND"
        input_pdb = "HETATM    1 ZN    ZN X   1      10.0  10.0  10.0\nATOM      1  N   ALA A   1"

        merged = _merge_hetatm(design, input_pdb)
        assert "HETATM" in merged
        assert "ZN" in merged
        # HETATM should be before END
        lines = merged.split("\n")
        hetatm_idx = next(i for i, l in enumerate(lines) if l.startswith("HETATM"))
        end_idx = next(i for i, l in enumerate(lines) if l.startswith("END"))
        assert hetatm_idx < end_idx

    def test_no_backend_raises(self):
        from pipeline_types import StepContext
        from pipeline_modules.backbone_generator import BackboneGenerator

        ctx = StepContext()
        gen = BackboneGenerator()
        with pytest.raises(RuntimeError, match="No inference backend"):
            gen.run(ctx)

    def test_sweep_aware(self):
        from pipeline_types import StepContext
        from pipeline_modules.backbone_generator import BackboneGenerator

        backend = MagicMock()
        backend.run_rfd3.return_value = {"designs": [{"pdb_content": "PDB"}]}

        ctx = StepContext(
            backend=backend,
            rfd3_config={"contig": "80-120"},
            metadata={
                "sweep_configs": [
                    {"cfg_scale": 1.5},
                    {"cfg_scale": 2.5},
                ],
                "sweep_designs_per_config": 1,
            },
        )

        gen = BackboneGenerator()
        result = gen.run(ctx)

        assert backend.run_rfd3.call_count == 2
        assert len(result.backbone_pdbs) == 2
        assert "sweep_results" in result.metadata["backbone_generator"]


# ============================================================
# Test SequenceDesigner module
# ============================================================

class TestSequenceDesigner:
    """Test SequenceDesigner with mock backend."""

    def test_basic_design(self):
        from pipeline_types import StepContext
        from pipeline_modules.sequence_designer import SequenceDesigner

        backend = MagicMock()
        backend.run_mpnn.return_value = {
            "sequences": [
                {"sequence": "MEVTI", "score": -1.2},
                {"sequence": "MKLIP", "score": -1.5},
            ]
        }
        ctx = StepContext(
            backend=backend,
            backbone_pdbs=["PDB_A"],
        )

        designer = SequenceDesigner(num_sequences=2)
        result = designer.run(ctx)

        assert len(result.sequences) == 2
        assert result.sequences[0].sequence == "MEVTI"
        assert result.sequences[0].score == -1.2

    def test_fasta_parsing(self):
        from pipeline_modules.sequence_designer import _parse_fasta

        fasta = ">seq_0001, score=-1.234\nMEVTI\n>seq_0002, score=-0.5\nMKLIP\n"
        results = _parse_fasta(fasta, backbone_index=0)

        assert len(results) == 2
        assert results[0].sequence == "MEVTI"
        assert abs(results[0].score - (-1.234)) < 0.001
        assert results[1].sequence == "MKLIP"

    def test_no_backbones_raises(self):
        from pipeline_types import StepContext
        from pipeline_modules.sequence_designer import SequenceDesigner

        ctx = StepContext(backend=MagicMock())
        with pytest.raises(ValueError, match="No backbone PDBs"):
            SequenceDesigner().run(ctx)


# ============================================================
# Test StructurePredictor module
# ============================================================

class TestStructurePredictor:
    """Test StructurePredictor with mock backend."""

    def test_basic_prediction(self):
        from pipeline_types import StepContext, SequenceResult
        from pipeline_modules.structure_predictor import StructurePredictor

        backend = MagicMock()
        backend.run_rf3.return_value = {
            "predictions": [{
                "mean_plddt": 0.85,
                "ptm": 0.9,
                "pdb_content": "ATOM...",
            }]
        }
        ctx = StepContext(
            backend=backend,
            sequences=[SequenceResult(sequence="MEVTI", backbone_index=0)],
        )

        predictor = StructurePredictor()
        result = predictor.run(ctx)

        assert len(result.predictions) == 1
        assert result.predictions[0].plddt == 0.85
        assert result.predictions[0].ptm == 0.9
        assert result.predictions[0].passed is True  # 0.9 >= 0.6 and 0.85 >= 0.65

    def test_failing_prediction(self):
        from pipeline_types import StepContext, SequenceResult
        from pipeline_modules.structure_predictor import StructurePredictor

        backend = MagicMock()
        backend.run_rf3.return_value = {
            "predictions": [{"mean_plddt": 0.3, "ptm": 0.2}]
        }
        ctx = StepContext(
            backend=backend,
            sequences=[SequenceResult(sequence="MEVTI", backbone_index=0)],
        )

        result = StructurePredictor().run(ctx)
        assert result.predictions[0].passed is False


# ============================================================
# Test Analyzer module
# ============================================================

class TestAnalyzer:
    """Test Analyzer with filter presets."""

    def test_passing_prediction(self):
        from pipeline_types import StepContext, PredictionResult
        from pipeline_modules.analyzer import Analyzer

        ctx = StepContext(predictions=[
            PredictionResult(sequence="MEVTI", plddt=0.8, ptm=0.85, passed=True),
        ])
        result = Analyzer(filter_preset="scout_relaxed").run(ctx)

        assert result.analysis["passed"] == 1
        assert result.analysis["pass_rate"] == 1.0
        assert result.analysis["candidates"][0]["passed"] is True

    def test_failing_prediction_min_threshold(self):
        from pipeline_types import StepContext, PredictionResult
        from pipeline_modules.analyzer import Analyzer

        ctx = StepContext(predictions=[
            PredictionResult(sequence="MEVTI", plddt=0.5, ptm=0.4, passed=False),
        ])
        result = Analyzer(filter_preset="scout_relaxed").run(ctx)

        assert result.analysis["passed"] == 0
        assert result.analysis["candidates"][0]["passed"] is False
        # Should have failed on both ptm and plddt
        failed = result.analysis["candidates"][0]["failed_filters"]
        assert len(failed) == 2

    def test_max_threshold_rmsd(self):
        """Verify that max thresholds correctly reject high RMSD values.

        Note: production_strict uses 'backbone_rmsd' but PredictionResult has 'rmsd'.
        The analyzer checks getattr(pred, metric), so we need the metric key to match
        a field on PredictionResult. We test with a custom preset that uses 'rmsd'.
        """
        from pipeline_types import StepContext, PredictionResult
        from pipeline_modules.analyzer import Analyzer

        # Temporarily inject a test preset that uses 'rmsd' (matches PredictionResult field)
        from filter_evaluator import FILTER_PRESETS
        FILTER_PRESETS["_test_rmsd"] = {
            "ptm": {"min": 0.6},
            "rmsd": {"max": 1.5},
        }
        try:
            ctx = StepContext(predictions=[
                PredictionResult(
                    sequence="MEVTI", plddt=0.9, ptm=0.9,
                    rmsd=5.0,  # Too high for max of 1.5
                    passed=True,
                ),
            ])
            result = Analyzer(filter_preset="_test_rmsd").run(ctx)

            candidate = result.analysis["candidates"][0]
            assert candidate["passed"] is False
            failed_strs = candidate["failed_filters"]
            assert any("rmsd" in f for f in failed_strs)
        finally:
            del FILTER_PRESETS["_test_rmsd"]

    def test_no_predictions(self):
        from pipeline_types import StepContext
        from pipeline_modules.analyzer import Analyzer

        ctx = StepContext()
        result = Analyzer().run(ctx)
        assert result.analysis["status"] == "no_predictions"

    def test_quality_tiers(self):
        from pipeline_types import StepContext, PredictionResult
        from pipeline_modules.analyzer import Analyzer

        ctx = StepContext(predictions=[
            PredictionResult(sequence="A", plddt=0.9, ptm=0.9, passed=True),   # tier 3
            PredictionResult(sequence="B", plddt=0.75, ptm=0.75, passed=True), # tier 2
            PredictionResult(sequence="C", plddt=0.62, ptm=0.62, passed=True), # tier 1
            PredictionResult(sequence="D", plddt=0.3, ptm=0.3, passed=False),  # tier 0
        ])
        result = Analyzer(filter_preset="scout_relaxed").run(ctx)

        tiers = [c["quality_tier"] for c in result.analysis["candidates"]]
        assert 3 in tiers
        assert 0 in tiers

    def test_imports_from_filter_evaluator(self):
        """Verify analyzer uses centralized FILTER_PRESETS, not a local copy."""
        from pipeline_modules.analyzer import FILTER_PRESETS as analyzer_presets
        from filter_evaluator import FILTER_PRESETS as evaluator_presets
        assert analyzer_presets is evaluator_presets


# ============================================================
# Test Reporter module
# ============================================================

class TestReporter:
    """Test Reporter module."""

    def test_basic_report(self):
        from pipeline_types import StepContext, SequenceResult, PredictionResult
        from pipeline_modules.reporter import Reporter

        ctx = StepContext(
            backbone_pdbs=["PDB1", "PDB2"],
            sequences=[SequenceResult(sequence="MEVTI", score=-1.0, backbone_index=0)],
            predictions=[PredictionResult(sequence="MEVTI", plddt=0.8, ptm=0.85, passed=True)],
            analysis={
                "filter_preset": "scout_relaxed",
                "passed": 1,
                "total_predictions": 1,
                "pass_rate": 1.0,
                "candidates": [{"rank": 1, "ptm": 0.85, "plddt": 0.8, "quality_tier": 3, "passed": True}],
            },
        )

        result = Reporter().run(ctx)
        assert result.report is not None
        assert "PIPELINE REPORT" in result.report
        assert "Backbones generated: 2" in result.report
        assert "Total sequences: 1" in result.report

    def test_empty_context_report(self):
        from pipeline_types import StepContext
        from pipeline_modules.reporter import Reporter

        ctx = StepContext()
        result = Reporter().run(ctx)
        assert result.report is not None
        assert "PIPELINE REPORT" in result.report


# ============================================================
# Test WorkflowRunner
# ============================================================

class TestWorkflowRunner:
    """Test WorkflowRunner with mock modules."""

    def test_sequential_run(self):
        from pipeline_types import StepContext
        from workflow_runner import WorkflowRunner

        step1 = MagicMock()
        step1.name = "step_a"
        step1.run.side_effect = lambda ctx: ctx

        step2 = MagicMock()
        step2.name = "step_b"
        step2.run.side_effect = lambda ctx: ctx

        ctx = StepContext()
        runner = WorkflowRunner()
        result = runner.run([step1, step2], ctx)

        step1.run.assert_called_once()
        step2.run.assert_called_once()
        assert result.metadata["step_status"]["step_a"] == "completed"
        assert result.metadata["step_status"]["step_b"] == "completed"

    def test_step_failure_propagates(self):
        from pipeline_types import StepContext
        from workflow_runner import WorkflowRunner

        step1 = MagicMock()
        step1.name = "failing_step"
        step1.run.side_effect = RuntimeError("boom")

        ctx = StepContext()
        runner = WorkflowRunner()
        with pytest.raises(RuntimeError, match="boom"):
            runner.run([step1], ctx)

        assert ctx.metadata["step_status"]["failing_step"] == "failed"

    def test_on_step_complete_callback(self):
        from pipeline_types import StepContext
        from workflow_runner import WorkflowRunner

        step1 = MagicMock()
        step1.name = "step_a"
        step1.run.side_effect = lambda ctx: ctx

        callbacks = []
        def on_complete(name, ctx):
            callbacks.append(name)

        ctx = StepContext()
        WorkflowRunner().run([step1], ctx, on_step_complete=on_complete)
        assert callbacks == ["step_a"]

    def test_from_spec(self):
        from workflow_runner import WorkflowRunner

        spec = {
            "name": "test_workflow",
            "backend": "http",
            "params": {"num_designs": 5},
            "steps": [
                {"module": "backbone_generator", "params": {"seed": 42}},
                {"module": "sequence_designer"},
            ],
        }

        with patch("inference_backend.HTTPBackend"):
            runner, steps, ctx = WorkflowRunner.from_spec(spec)

        assert len(steps) == 2
        assert steps[0].name == "backbone_generator"
        assert steps[1].name == "sequence_designer"
        assert ctx.params["num_designs"] == 5
        assert ctx.metadata["workflow_name"] == "test_workflow"

    def test_from_spec_unknown_module(self):
        from workflow_runner import WorkflowRunner

        spec = {
            "name": "bad",
            "steps": [{"module": "nonexistent"}],
        }

        with patch("inference_backend.create_backend", return_value=MagicMock()):
            with pytest.raises(ValueError, match="Unknown module"):
                WorkflowRunner.from_spec(spec)

    def test_get_progress(self):
        from pipeline_types import StepContext
        from workflow_runner import WorkflowRunner

        ctx = StepContext(
            backbone_pdbs=["PDB1"],
            metadata={
                "workflow_name": "test",
                "workflow_steps": ["backbone_generator", "sequence_designer"],
                "step_status": {"backbone_generator": "completed"},
                "step_timings": {"backbone_generator": 2.5},
                "current_step": "sequence_designer",
                "current_step_index": 1,
            },
        )

        runner = WorkflowRunner()
        progress = runner.get_progress(ctx)

        assert progress["status"] == "running"
        assert progress["current_step"] == "sequence_designer"
        assert progress["total_steps"] == 2
        assert "backbone_generator" in progress["step_results"]


# ============================================================
# Test IterationStrategies
# ============================================================

class TestIterationStrategies:
    """Test ScoutStrategy and SweepStrategy."""

    def test_sweep_grid_generation(self):
        from iteration_strategies import SweepStrategy

        strategy = SweepStrategy(grid={
            "cfg_scale": [1.5, 2.0],
            "length": ["100-120", "120-150"],
        })

        configs = strategy.generate_grid()
        assert len(configs) == 4  # 2 x 2
        # Check all combinations exist
        cfg_scales = {c["cfg_scale"] for c in configs}
        assert cfg_scales == {1.5, 2.0}

    def test_sweep_static_grid(self):
        from iteration_strategies import SweepStrategy

        configs = SweepStrategy.generate_sweep_grid(
            sizes=["100-120"],
            cfg_scales=[1.5, 2.0, 2.5],
            base_params={"step_scale": 1.5},
        )

        assert len(configs) == 3
        assert all(c["step_scale"] == 1.5 for c in configs)
        assert configs[0]["cfg_scale"] == 1.5

    def test_sweep_run_stores_configs(self):
        from pipeline_types import StepContext
        from iteration_strategies import SweepStrategy

        ctx = StepContext(
            rfd3_config={"contig": "80-120"},
        )
        strategy = SweepStrategy(
            grid={"cfg_scale": [1.5, 2.5]},
            designs_per_config=5,
        )

        result = strategy.run(ctx)
        assert "sweep_configs" in result.metadata
        assert len(result.metadata["sweep_configs"]) == 2
        assert result.metadata["sweep_designs_per_config"] == 5

    def test_scout_filters_backbones(self):
        from pipeline_types import StepContext
        from iteration_strategies import ScoutStrategy

        backend = MagicMock()
        # MPNN returns a sequence
        backend.run_mpnn.return_value = {
            "sequences": [{"sequence": "MEVTI"}]
        }
        # RF3: first backbone passes, second fails
        backend.run_rf3.side_effect = [
            {"predictions": [{"ptm": 0.8, "mean_plddt": 0.85}]},
            {"predictions": [{"ptm": 0.3, "mean_plddt": 0.4}]},
        ]

        ctx = StepContext(
            backend=backend,
            backbone_pdbs=["GOOD_PDB", "BAD_PDB"],
        )

        scout = ScoutStrategy(ptm_threshold=0.6, plddt_threshold=0.65)
        result = scout.run(ctx)

        assert len(result.backbone_pdbs) == 1
        assert result.backbone_pdbs[0] == "GOOD_PDB"
        assert result.metadata["scout_filter"]["original_backbones"] == 2
        assert result.metadata["scout_filter"]["passing_backbones"] == 1


# ============================================================
# Test rfd3_config_generator scaffold_to_rfd3_params
# ============================================================

class TestScaffoldToRfd3Params:
    """Test the scaffold_to_rfd3_params bridge function."""

    def test_dict_input(self):
        from rfd3_config_generator import scaffold_to_rfd3_params

        scaffold = {
            "pdb_content": "ATOM...",
            "metal_type": "TB",
            "ligand_code": "CIT",
            "residues": ["A35", "A40"],
        }

        params = scaffold_to_rfd3_params(scaffold)
        assert params["pdb_content"] == "ATOM..."
        assert params["select_fixed_atoms"]["X1"] == "all"
        assert params["select_fixed_atoms"]["L1"] == "all"
        assert params["hotspots"] == ["L1"]

    def test_no_metal_no_ligand(self):
        from rfd3_config_generator import scaffold_to_rfd3_params

        scaffold = {"pdb_content": "ATOM..."}
        params = scaffold_to_rfd3_params(scaffold)
        assert "pdb_content" in params
        assert "select_fixed_atoms" not in params


# ============================================================
# Run
# ============================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
