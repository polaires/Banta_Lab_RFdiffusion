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
