"""Tests for PipelineHandler and PipelineSession."""
import pytest
import json
import os
import tempfile
import shutil
from datetime import datetime

from pipeline_session import (
    PipelineSession,
    PipelineSessionManager,
    PipelineMode,
    PipelineDesign,
    FilterThresholds,
    SweepConfig,
    DesignStatus,
    PipelineProgress,
    generate_sweep_configs,
)
from pipeline_handler import (
    PipelineHandler,
    handle_pipeline_design,
    handle_pipeline_status,
    handle_pipeline_cancel,
    handle_pipeline_export,
)


class TestFilterThresholds:
    """Test suite for FilterThresholds."""

    def test_default_values(self):
        filters = FilterThresholds()
        assert filters.plddt == 0.80
        assert filters.ptm == 0.80
        assert filters.pae == 5.0
        assert filters.plddt_discard == 0.60

    def test_strict_preset(self):
        filters = FilterThresholds.strict()
        assert filters.plddt == 0.80
        assert filters.ptm == 0.80
        assert filters.pae == 5.0

    def test_relaxed_preset(self):
        filters = FilterThresholds.relaxed()
        assert filters.plddt == 0.70
        assert filters.ptm == 0.65
        assert filters.pae == 10.0

    def test_to_dict(self):
        filters = FilterThresholds(plddt=0.85, ptm=0.90, pae=4.0)
        d = filters.to_dict()
        assert d["plddt"] == 0.85
        assert d["ptm"] == 0.90
        assert d["pae"] == 4.0


class TestSweepConfig:
    """Test suite for SweepConfig."""

    def test_init(self):
        config = SweepConfig(
            name="small_low_cfg",
            contig_size="small",
            contig_range="50-70",
            cfg_scale=1.5,
            num_designs=10,
        )
        assert config.name == "small_low_cfg"
        assert config.contig_size == "small"
        assert config.contig_range == "50-70"
        assert config.cfg_scale == 1.5
        assert config.num_designs == 10

    def test_to_dict(self):
        config = SweepConfig(
            name="medium_mid_cfg",
            contig_size="medium",
            contig_range="70-90",
            cfg_scale=2.0,
            num_designs=20,
        )
        d = config.to_dict()
        assert d["name"] == "medium_mid_cfg"
        assert d["contig_size"] == "medium"


class TestPipelineDesign:
    """Test suite for PipelineDesign."""

    def test_init(self):
        design = PipelineDesign(
            name="test_design_001",
            config_name="small_low_cfg",
            sequence="MVKLSTG",
            pdb_content="ATOM...",
            plddt=0.85,
            ptm=0.90,
            pae=4.5,
            status=DesignStatus.PASS,
        )
        assert design.name == "test_design_001"
        assert design.status == DesignStatus.PASS

    def test_to_dict(self):
        design = PipelineDesign(
            name="test_design",
            config_name="config",
            sequence="MVKLSTG",
            pdb_content="ATOM...",
            plddt=0.85,
            ptm=0.90,
            pae=4.5,
            status=DesignStatus.REVIEW,
        )
        d = design.to_dict()
        assert d["status"] == "review"
        assert d["name"] == "test_design"

    def test_from_dict(self):
        data = {
            "name": "test_design",
            "config_name": "config",
            "sequence": "MVKLSTG",
            "pdb_content": "ATOM...",
            "plddt": 0.85,
            "ptm": 0.90,
            "pae": 4.5,
            "status": "pass",
            "timestamp": "2026-01-23T10:00:00",
        }
        design = PipelineDesign.from_dict(data)
        assert design.status == DesignStatus.PASS
        assert design.plddt == 0.85


class TestPipelineProgress:
    """Test suite for PipelineProgress."""

    def test_default_values(self):
        progress = PipelineProgress()
        assert progress.current_config == 0
        assert progress.total_generated == 0
        assert progress.pass_rate == 0.0

    def test_pass_rate_calculation(self):
        progress = PipelineProgress(
            total_generated=100,
            total_passing=30,
        )
        assert progress.pass_rate == 0.30

    def test_to_dict(self):
        progress = PipelineProgress(
            current_config=3,
            total_configs=9,
            current_design=7,
            total_generated=27,
            total_passing=12,
        )
        d = progress.to_dict()
        assert d["current_config"] == 3
        assert d["total_configs"] == 9
        assert "pass_rate" in d


class TestPipelineSession:
    """Test suite for PipelineSession."""

    @pytest.fixture
    def temp_session_dir(self):
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir)

    @pytest.fixture
    def sample_session(self, temp_session_dir):
        return PipelineSession(
            session_dir=temp_session_dir,
            session_id="test_session_001",
            mode=PipelineMode.SWEEP,
            metal="TB",
            ligand="CIT",
            filters=FilterThresholds.strict(),
        )

    def test_session_creates_structure(self, sample_session, temp_session_dir):
        assert os.path.exists(os.path.join(temp_session_dir, "configs"))
        assert os.path.exists(os.path.join(temp_session_dir, "results"))
        assert os.path.exists(os.path.join(temp_session_dir, "passing_designs"))
        assert os.path.exists(os.path.join(temp_session_dir, "review_designs"))

    def test_session_saves_metadata(self, sample_session, temp_session_dir):
        meta_path = os.path.join(temp_session_dir, "meta.json")
        assert os.path.exists(meta_path)

        with open(meta_path, "r") as f:
            meta = json.load(f)

        assert meta["session_id"] == "test_session_001"
        assert meta["mode"] == "sweep"
        assert meta["metal"] == "TB"
        assert meta["ligand"] == "CIT"

    def test_evaluate_design_pass(self, sample_session):
        status = sample_session.evaluate_design(
            plddt=0.90,
            ptm=0.85,
            pae=3.0,
        )
        assert status == DesignStatus.PASS

    def test_evaluate_design_review(self, sample_session):
        # Below strict but above relaxed
        status = sample_session.evaluate_design(
            plddt=0.75,  # Below 0.80 (strict), above 0.70 (relaxed)
            ptm=0.70,    # Below 0.80, above 0.65
            pae=8.0,     # Above 5.0, below 10.0
        )
        assert status == DesignStatus.REVIEW

    def test_evaluate_design_fail(self, sample_session):
        status = sample_session.evaluate_design(
            plddt=0.55,  # Below discard threshold
            ptm=0.50,
            pae=15.0,
        )
        assert status == DesignStatus.FAIL

    def test_add_passing_design(self, sample_session, temp_session_dir):
        design = PipelineDesign(
            name="passing_001",
            config_name="small_low_cfg",
            sequence="MVKLSTGAEI",
            pdb_content="ATOM 1...",
            plddt=0.90,
            ptm=0.85,
            pae=3.0,
            status=DesignStatus.PASS,
        )

        stored = sample_session.add_design(design)

        assert stored is True
        assert sample_session.progress.total_generated == 1
        assert sample_session.progress.total_passing == 1

        # Check PDB was saved
        pdb_path = os.path.join(temp_session_dir, "passing_designs", "passing_001.pdb")
        assert os.path.exists(pdb_path)

    def test_add_failed_design_not_stored(self, sample_session, temp_session_dir):
        design = PipelineDesign(
            name="failed_001",
            config_name="small_low_cfg",
            sequence="MVKLSTGAEI",
            pdb_content="ATOM 1...",
            plddt=0.50,
            ptm=0.40,
            pae=20.0,
            status=DesignStatus.FAIL,
        )

        stored = sample_session.add_design(design)

        assert stored is False
        assert sample_session.progress.total_generated == 1
        assert sample_session.progress.total_failed == 1
        assert sample_session.progress.total_passing == 0

    def test_export_fasta(self, sample_session):
        # Add some designs
        sample_session.add_design(PipelineDesign(
            name="pass_001",
            config_name="config1",
            sequence="MVKLSTG",
            pdb_content="ATOM...",
            plddt=0.90,
            ptm=0.85,
            pae=3.0,
            status=DesignStatus.PASS,
        ))
        sample_session.add_design(PipelineDesign(
            name="review_001",
            config_name="config1",
            sequence="AEITKLM",
            pdb_content="ATOM...",
            plddt=0.75,
            ptm=0.70,
            pae=8.0,
            status=DesignStatus.REVIEW,
        ))

        # Export only passing
        fasta = sample_session.export_fasta(include_review=False)
        assert ">pass_001" in fasta
        assert "MVKLSTG" in fasta
        assert "review_001" not in fasta

        # Export including review
        fasta_all = sample_session.export_fasta(include_review=True)
        assert ">pass_001" in fasta_all
        assert ">review_001" in fasta_all

    def test_finalize_session(self, sample_session, temp_session_dir):
        sample_session.add_design(PipelineDesign(
            name="test_001",
            config_name="config1",
            sequence="MVKLSTG",
            pdb_content="ATOM...",
            plddt=0.90,
            ptm=0.85,
            pae=3.0,
            status=DesignStatus.PASS,
        ))

        summary = sample_session.finalize(status="completed")

        assert summary["status"] == "completed"
        assert summary["statistics"]["total_passing"] == 1

        # Check summary file exists
        summary_path = os.path.join(temp_session_dir, "summary.json")
        assert os.path.exists(summary_path)

    def test_save_and_load_progress(self, sample_session):
        sample_session.progress.current_config = 5
        sample_session.progress.total_generated = 50
        sample_session.save_progress()

        # Create new session that loads progress
        sample_session.progress = PipelineProgress()
        loaded = sample_session.load_progress()

        assert loaded.current_config == 5
        assert loaded.total_generated == 50


class TestPipelineSessionManager:
    """Test suite for PipelineSessionManager."""

    @pytest.fixture
    def temp_base_dir(self):
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir)

    @pytest.fixture
    def manager(self, temp_base_dir):
        return PipelineSessionManager(temp_base_dir)

    def test_create_session(self, manager):
        session = manager.create_session(
            mode=PipelineMode.SWEEP,
            metal="TB",
            ligand="CIT",
        )

        assert session is not None
        assert session.mode == PipelineMode.SWEEP
        assert session.metal == "TB"
        assert "TB_CIT" in session.session_id

    def test_get_session(self, manager):
        # Create a session
        session = manager.create_session(
            mode=PipelineMode.PRODUCTION,
            metal="EU",
            ligand="PQQ",
        )
        session_id = session.session_id

        # Retrieve it
        retrieved = manager.get_session(session_id)

        assert retrieved is not None
        assert retrieved.session_id == session_id
        assert retrieved.mode == PipelineMode.PRODUCTION

    def test_get_nonexistent_session(self, manager):
        retrieved = manager.get_session("nonexistent_session")
        assert retrieved is None

    def test_list_sessions(self, manager):
        # Create multiple sessions
        manager.create_session(PipelineMode.SWEEP, "TB", "CIT")
        manager.create_session(PipelineMode.PRODUCTION, "EU", "PQQ")

        sessions = manager.list_sessions()

        assert len(sessions) == 2
        assert all("session_id" in s for s in sessions)


class TestGenerateSweepConfigs:
    """Test suite for generate_sweep_configs helper."""

    def test_default_configs(self):
        configs = generate_sweep_configs()

        # 3 sizes x 3 CFG values = 9 configs
        assert len(configs) == 9

    def test_config_names(self):
        configs = generate_sweep_configs()
        names = [c.name for c in configs]

        assert "small_low_cfg" in names
        assert "medium_mid_cfg" in names
        assert "large_high_cfg" in names

    def test_custom_designs_per_config(self):
        configs = generate_sweep_configs(designs_per_config=20)

        for config in configs:
            assert config.num_designs == 20


class TestPipelineHandler:
    """Test suite for PipelineHandler."""

    @pytest.fixture
    def temp_base_dir(self):
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir)

    @pytest.fixture
    def handler(self, temp_base_dir):
        return PipelineHandler(temp_base_dir)

    def test_start_sweep(self, handler):
        result = handler.start_sweep(
            metal="TB",
            ligand="CIT",
            motif_pdb="ATOM 1...",
            designs_per_config=5,
        )

        assert result["status"] == "running"
        assert result["mode"] == "sweep"
        assert result["total_configs"] == 9  # Default 3x3
        assert "session_id" in result

    def test_start_sweep_with_custom_configs(self, handler):
        custom_configs = [
            {
                "name": "custom_1",
                "contig_size": "small",
                "contig_range": "40-50",
                "cfg_scale": 1.5,
                "num_designs": 5,
            }
        ]

        result = handler.start_sweep(
            metal="TB",
            ligand="CIT",
            motif_pdb="ATOM 1...",
            sweep_configs=custom_configs,
        )

        assert result["total_configs"] == 1

    def test_start_production(self, handler):
        result = handler.start_production(
            metal="TB",
            ligand="CIT",
            motif_pdb="ATOM 1...",
            production_config={
                "contig_size": "medium",
                "contig_range": "70-90",
                "cfg_scale": 2.0,
            },
            num_designs=100,
        )

        assert result["status"] == "running"
        assert result["mode"] == "production"
        assert result["total_configs"] == 1
        assert result["num_designs"] == 100

    def test_get_status_not_found(self, handler):
        result = handler.get_status("nonexistent_session")

        assert result["status"] == "not_found"

    def test_get_status_running(self, handler):
        # Start a session first
        start_result = handler.start_sweep(
            metal="TB",
            ligand="CIT",
            motif_pdb="ATOM 1...",
        )
        session_id = start_result["session_id"]

        # Get status
        status = handler.get_status(session_id)

        assert status["session_id"] == session_id
        assert status["status"] == "running"
        assert "current_config" in status

    def test_cancel_session(self, handler):
        # Start a session
        start_result = handler.start_sweep(
            metal="TB",
            ligand="CIT",
            motif_pdb="ATOM 1...",
        )
        session_id = start_result["session_id"]

        # Cancel it
        cancel_result = handler.cancel(session_id)

        assert cancel_result["status"] == "cancelled"
        assert cancel_result["session_id"] == session_id

    def test_export_fasta_empty_session(self, handler):
        # Start a session
        start_result = handler.start_sweep(
            metal="TB",
            ligand="CIT",
            motif_pdb="ATOM 1...",
        )
        session_id = start_result["session_id"]

        # Export (should be empty)
        export_result = handler.export_fasta(session_id)

        assert export_result["status"] == "success"
        assert export_result["num_sequences"] == 0


class TestHandlerFunctions:
    """Test suite for handler entry point functions."""

    @pytest.fixture
    def temp_base_dir(self):
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir)

    def test_handle_pipeline_design_sweep(self, temp_base_dir):
        result = handle_pipeline_design({
            "mode": "sweep",
            "metal": "TB",
            "ligand": "CIT",
            "motif_pdb": "ATOM 1...",
            "base_dir": temp_base_dir,
            "designs_per_config": 5,
        })

        assert result["status"] == "completed"
        assert result["result"]["mode"] == "sweep"

    def test_handle_pipeline_design_production(self, temp_base_dir):
        result = handle_pipeline_design({
            "mode": "production",
            "metal": "TB",
            "ligand": "CIT",
            "motif_pdb": "ATOM 1...",
            "base_dir": temp_base_dir,
            "production_config": {
                "contig_size": "medium",
                "contig_range": "70-90",
                "cfg_scale": 2.0,
            },
            "num_designs": 50,
        })

        assert result["status"] == "completed"
        assert result["result"]["mode"] == "production"

    def test_handle_pipeline_design_invalid_mode(self, temp_base_dir):
        result = handle_pipeline_design({
            "mode": "invalid",
            "metal": "TB",
            "ligand": "CIT",
            "motif_pdb": "ATOM 1...",
            "base_dir": temp_base_dir,
        })

        assert result["status"] == "failed"
        assert "Unknown mode" in result["error"]

    def test_handle_pipeline_status(self, temp_base_dir):
        # Create a session first
        design_result = handle_pipeline_design({
            "mode": "sweep",
            "metal": "TB",
            "ligand": "CIT",
            "motif_pdb": "ATOM 1...",
            "base_dir": temp_base_dir,
        })
        session_id = design_result["result"]["session_id"]

        # Get status
        status_result = handle_pipeline_status({
            "session_id": session_id,
            "base_dir": temp_base_dir,
        })

        assert status_result["status"] == "completed"
        assert status_result["result"]["session_id"] == session_id

    def test_handle_pipeline_cancel(self, temp_base_dir):
        # Create a session
        design_result = handle_pipeline_design({
            "mode": "sweep",
            "metal": "TB",
            "ligand": "CIT",
            "motif_pdb": "ATOM 1...",
            "base_dir": temp_base_dir,
        })
        session_id = design_result["result"]["session_id"]

        # Cancel
        cancel_result = handle_pipeline_cancel({
            "session_id": session_id,
            "base_dir": temp_base_dir,
        })

        assert cancel_result["status"] == "completed"
        assert cancel_result["result"]["status"] == "cancelled"

    def test_handle_pipeline_export(self, temp_base_dir):
        # Create a session
        design_result = handle_pipeline_design({
            "mode": "sweep",
            "metal": "TB",
            "ligand": "CIT",
            "motif_pdb": "ATOM 1...",
            "base_dir": temp_base_dir,
        })
        session_id = design_result["result"]["session_id"]

        # Export
        export_result = handle_pipeline_export({
            "session_id": session_id,
            "include_review": True,
            "base_dir": temp_base_dir,
        })

        assert export_result["status"] == "completed"
        assert "filename" in export_result["result"]


class TestEdgeCases:
    """Test edge cases and error handling."""

    @pytest.fixture
    def temp_base_dir(self):
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir)

    def test_session_with_custom_filters(self, temp_base_dir):
        manager = PipelineSessionManager(temp_base_dir)

        custom_filters = FilterThresholds(
            plddt=0.85,
            ptm=0.90,
            pae=3.0,
        )

        session = manager.create_session(
            mode=PipelineMode.SWEEP,
            metal="TB",
            ligand="CIT",
            filters=custom_filters,
        )

        # Test that custom filters are used
        # High confidence should pass
        assert session.evaluate_design(0.90, 0.92, 2.5) == DesignStatus.PASS
        # Below custom threshold should not pass strict
        assert session.evaluate_design(0.82, 0.85, 4.0) != DesignStatus.PASS

    def test_best_design_tracking(self, temp_base_dir):
        session = PipelineSession(
            session_dir=temp_base_dir,
            session_id="test",
            mode=PipelineMode.SWEEP,
            metal="TB",
            ligand="CIT",
            filters=FilterThresholds.strict(),
        )

        # Add design with plddt 0.85
        session.add_design(PipelineDesign(
            name="design_1",
            config_name="config",
            sequence="MVKLSTG",
            pdb_content="ATOM...",
            plddt=0.85,
            ptm=0.85,
            pae=4.0,
            status=DesignStatus.PASS,
        ))

        assert session.progress.best_design["name"] == "design_1"
        assert session.progress.best_design["plddt"] == 0.85

        # Add better design
        session.add_design(PipelineDesign(
            name="design_2",
            config_name="config",
            sequence="MVKLSTG",
            pdb_content="ATOM...",
            plddt=0.92,
            ptm=0.90,
            pae=3.0,
            status=DesignStatus.PASS,
        ))

        assert session.progress.best_design["name"] == "design_2"
        assert session.progress.best_design["plddt"] == 0.92

    def test_review_zone_boundary(self, temp_base_dir):
        session = PipelineSession(
            session_dir=temp_base_dir,
            session_id="test",
            mode=PipelineMode.SWEEP,
            metal="TB",
            ligand="CIT",
            filters=FilterThresholds.strict(),
        )

        # Exactly at relaxed threshold - should be REVIEW
        status = session.evaluate_design(0.70, 0.65, 10.0)
        assert status == DesignStatus.REVIEW

        # Just below relaxed threshold - should be FAIL
        status = session.evaluate_design(0.69, 0.65, 10.0)
        assert status == DesignStatus.FAIL
