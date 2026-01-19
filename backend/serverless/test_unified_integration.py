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
            json.dump({"version": "1.0.0", "created": "", "designs": []}, f)

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
