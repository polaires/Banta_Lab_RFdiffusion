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
