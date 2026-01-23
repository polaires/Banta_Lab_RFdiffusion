"""
Pipeline Handler

Orchestrates parameter sweeps and production runs for metal-ligand design.
Provides real-time progress updates and intelligent design filtering.
"""
import os
import json
import tempfile
from datetime import datetime
from typing import Dict, Any, Optional, List, Generator, Callable

from pipeline_session import (
    PipelineSession,
    PipelineSessionManager,
    PipelineMode,
    PipelineDesign,
    FilterThresholds,
    SweepConfig,
    DesignStatus,
    generate_sweep_configs,
)

# Import inference utilities (available in serverless environment)
try:
    from inference_utils import (
        run_rfd3_inference,
        run_rf3_inference,
        run_mpnn_inference,
    )
    INFERENCE_AVAILABLE = True
except ImportError:
    INFERENCE_AVAILABLE = False

# Import unified analyzer for design evaluation
try:
    from unified_analyzer import UnifiedDesignAnalyzer
    ANALYZER_AVAILABLE = True
except ImportError:
    ANALYZER_AVAILABLE = False

# Import filter evaluator
try:
    from filter_evaluator import FilterEvaluator
    FILTER_EVALUATOR_AVAILABLE = True
except ImportError:
    FILTER_EVALUATOR_AVAILABLE = False


# Default batch size for progress updates
DEFAULT_BATCH_SIZE = 10


class PipelineHandler:
    """
    Handles pipeline execution for metal-ligand design.

    Supports:
    - Parameter sweeps across configurations
    - Production runs with single configuration
    - Real-time progress tracking
    - Intelligent design filtering
    """

    def __init__(self, base_dir: str):
        """
        Initialize pipeline handler.

        Args:
            base_dir: Base directory for pipeline sessions
        """
        self.base_dir = base_dir
        self.session_manager = PipelineSessionManager(base_dir)

        if ANALYZER_AVAILABLE:
            self.analyzer = UnifiedDesignAnalyzer()
        else:
            self.analyzer = None

    def start_sweep(
        self,
        metal: str,
        ligand: str,
        motif_pdb: str,
        sweep_configs: Optional[List[Dict[str, Any]]] = None,
        filters: Optional[Dict[str, float]] = None,
        designs_per_config: int = 10,
        batch_size: int = DEFAULT_BATCH_SIZE,
    ) -> Dict[str, Any]:
        """
        Start a parameter sweep pipeline.

        Args:
            metal: Metal code (e.g., "TB")
            ligand: Ligand code (e.g., "CIT")
            motif_pdb: PDB content for the metal-ligand motif
            sweep_configs: Optional custom configurations
            filters: Filter thresholds {plddt, ptm, pae}
            designs_per_config: Number of designs per configuration
            batch_size: Batch size for progress updates

        Returns:
            Initial status with session_id
        """
        # Parse filters
        if filters:
            filter_thresholds = FilterThresholds(**filters)
        else:
            filter_thresholds = FilterThresholds.strict()

        # Create session
        session = self.session_manager.create_session(
            mode=PipelineMode.SWEEP,
            metal=metal,
            ligand=ligand,
            filters=filter_thresholds,
        )

        # Generate or parse sweep configurations
        if sweep_configs:
            configs = [
                SweepConfig(**cfg) for cfg in sweep_configs
            ]
        else:
            configs = generate_sweep_configs(
                metal=metal,
                ligand=ligand,
                designs_per_config=designs_per_config,
            )

        # Save configurations
        for i, config in enumerate(configs):
            session.save_config(config, i)

        # Update progress
        session.progress.total_configs = len(configs)
        session.progress.designs_per_config = designs_per_config
        session.save_progress()

        # Store motif for later use
        motif_path = os.path.join(session.session_dir, "motif.pdb")
        with open(motif_path, "w") as f:
            f.write(motif_pdb)

        return {
            "session_id": session.session_id,
            "status": "running",
            "mode": "sweep",
            "total_configs": len(configs),
            "designs_per_config": designs_per_config,
        }

    def start_production(
        self,
        metal: str,
        ligand: str,
        motif_pdb: str,
        production_config: Dict[str, Any],
        num_designs: int = 1000,
        filters: Optional[Dict[str, float]] = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
    ) -> Dict[str, Any]:
        """
        Start a production pipeline run.

        Args:
            metal: Metal code (e.g., "TB")
            ligand: Ligand code (e.g., "CIT")
            motif_pdb: PDB content for the metal-ligand motif
            production_config: Single configuration to run
            num_designs: Total number of designs to generate
            filters: Filter thresholds {plddt, ptm, pae}
            batch_size: Batch size for progress updates

        Returns:
            Initial status with session_id
        """
        # Parse filters
        if filters:
            filter_thresholds = FilterThresholds(**filters)
        else:
            filter_thresholds = FilterThresholds.strict()

        # Create session
        session = self.session_manager.create_session(
            mode=PipelineMode.PRODUCTION,
            metal=metal,
            ligand=ligand,
            filters=filter_thresholds,
        )

        # Create single config
        config = SweepConfig(
            name="production",
            contig_size=production_config.get("contig_size", "medium"),
            contig_range=production_config.get("contig_range", "70-90"),
            cfg_scale=production_config.get("cfg_scale", 2.0),
            num_designs=num_designs,
        )
        session.save_config(config, 0)

        # Update progress
        session.progress.total_configs = 1
        session.progress.designs_per_config = num_designs
        session.save_progress()

        # Store motif
        motif_path = os.path.join(session.session_dir, "motif.pdb")
        with open(motif_path, "w") as f:
            f.write(motif_pdb)

        return {
            "session_id": session.session_id,
            "status": "running",
            "mode": "production",
            "total_configs": 1,
            "num_designs": num_designs,
        }

    def process_batch(
        self,
        session: PipelineSession,
        config: SweepConfig,
        motif_pdb: str,
        start_idx: int,
        batch_size: int,
    ) -> List[PipelineDesign]:
        """
        Process a batch of designs.

        Args:
            session: Active pipeline session
            config: Current sweep configuration
            motif_pdb: PDB content for the motif
            start_idx: Starting design index
            batch_size: Number of designs in batch

        Returns:
            List of generated designs
        """
        if not INFERENCE_AVAILABLE:
            raise RuntimeError("Inference modules not available")

        designs = []

        for i in range(batch_size):
            if session.cancelled:
                break

            design_idx = start_idx + i
            design_name = f"{config.name}_{design_idx:03d}"

            # Generate backbone with RFD3
            rfd3_params = {
                "contig": config.contig_range,
                "step_scale": config.cfg_scale,
                "num_designs": 1,
                "pdb_content": motif_pdb,
                "seed": design_idx,  # Use index as seed for reproducibility
            }

            try:
                rfd3_result = run_rfd3_inference(**rfd3_params)
                if rfd3_result.get("status") != "completed":
                    continue

                backbone_pdb = rfd3_result.get("result", {}).get(
                    "designs", [{}]
                )[0].get("content", "")

                if not backbone_pdb:
                    continue

                # Generate sequence with LigandMPNN
                mpnn_params = {
                    "pdb_content": backbone_pdb,
                    "num_sequences": 3,  # Generate 3 sequences per backbone
                    "temperature": 0.1,
                    "model_type": "ligand_mpnn",
                }

                mpnn_result = run_mpnn_inference(**mpnn_params)
                if mpnn_result.get("status") != "completed":
                    continue

                sequences = mpnn_result.get("result", {}).get("sequences", [])
                if not sequences:
                    continue

                # Predict structure with RF3 for each sequence
                for seq_idx, seq_data in enumerate(sequences):
                    sequence = seq_data.get("content", "").strip()
                    if not sequence:
                        continue

                    rf3_params = {"sequence": sequence}
                    rf3_result = run_rf3_inference(**rf3_params)

                    if rf3_result.get("status") != "completed":
                        continue

                    predictions = rf3_result.get("result", {}).get("predictions", [])
                    if not predictions:
                        continue

                    prediction = predictions[0]
                    pred_pdb = prediction.get("content", "")

                    # Extract confidence metrics
                    confidences = rf3_result.get("result", {}).get("confidences", {})
                    summary = confidences.get("summary_confidences", {})

                    plddt = summary.get("overall_plddt", 0)
                    ptm = summary.get("ptm", 0)
                    pae = summary.get("overall_pae", 100)

                    # Normalize pLDDT if needed (0-100 to 0-1)
                    if plddt > 1:
                        plddt = plddt / 100.0

                    # Evaluate against filters
                    status = session.evaluate_design(plddt, ptm, pae)

                    # Create design object
                    full_name = f"{design_name}_seq{seq_idx}"
                    design = PipelineDesign(
                        name=full_name,
                        config_name=config.name,
                        sequence=sequence,
                        pdb_content=pred_pdb,
                        plddt=plddt,
                        ptm=ptm,
                        pae=pae,
                        status=status,
                    )

                    # Add to session
                    session.add_design(design)
                    designs.append(design)

                    # Update current design counter
                    session.progress.current_design = design_idx + 1

            except Exception as e:
                print(f"Error processing design {design_name}: {e}")
                continue

        return designs

    def run_sweep(
        self,
        session_id: str,
        progress_callback: Optional[Callable[[Dict[str, Any]], None]] = None,
    ) -> Dict[str, Any]:
        """
        Execute a sweep pipeline.

        Args:
            session_id: Session identifier
            progress_callback: Optional callback for progress updates

        Returns:
            Final results
        """
        session = self.session_manager.get_session(session_id)
        if not session:
            return {"status": "failed", "error": "Session not found"}

        # Load motif
        motif_path = os.path.join(session.session_dir, "motif.pdb")
        with open(motif_path, "r") as f:
            motif_pdb = f.read()

        # Load configurations
        configs_dir = os.path.join(session.session_dir, "configs")
        configs = []
        for filename in sorted(os.listdir(configs_dir)):
            if filename.endswith(".json"):
                filepath = os.path.join(configs_dir, filename)
                with open(filepath, "r") as f:
                    config_data = json.load(f)
                configs.append(SweepConfig(**config_data))

        # Process each configuration
        for config_idx, config in enumerate(configs):
            if session.cancelled:
                break

            session.progress.current_config = config_idx + 1
            session.progress.current_design = 0
            session.save_progress()

            # Process in batches
            for batch_start in range(0, config.num_designs, DEFAULT_BATCH_SIZE):
                if session.cancelled:
                    break

                batch_size = min(
                    DEFAULT_BATCH_SIZE,
                    config.num_designs - batch_start
                )

                self.process_batch(
                    session=session,
                    config=config,
                    motif_pdb=motif_pdb,
                    start_idx=batch_start,
                    batch_size=batch_size,
                )

                # Report progress
                if progress_callback:
                    progress_callback(session.progress.to_dict())

        # Finalize session
        status = "cancelled" if session.cancelled else "completed"
        summary = session.finalize(status=status)

        return {
            "status": status,
            "results": session.get_results(),
            "summary": summary,
        }

    def run_production(
        self,
        session_id: str,
        progress_callback: Optional[Callable[[Dict[str, Any]], None]] = None,
    ) -> Dict[str, Any]:
        """
        Execute a production pipeline.

        Args:
            session_id: Session identifier
            progress_callback: Optional callback for progress updates

        Returns:
            Final results
        """
        # Production is just a sweep with single config
        return self.run_sweep(session_id, progress_callback)

    def get_status(self, session_id: str) -> Dict[str, Any]:
        """
        Get current pipeline status.

        Args:
            session_id: Session identifier

        Returns:
            Current status and progress
        """
        session = self.session_manager.get_session(session_id)
        if not session:
            return {"status": "not_found", "error": "Session not found"}

        # Load metadata for status
        meta_path = os.path.join(session.session_dir, "meta.json")
        with open(meta_path, "r") as f:
            meta = json.load(f)

        pipeline_status = meta.get("status", "unknown")

        response = {
            "session_id": session_id,
            "status": pipeline_status,
            "mode": session.mode.value,
            **session.progress.to_dict(),
        }

        # Include results if completed
        if pipeline_status == "completed":
            response["results"] = session.get_results()

            # Load summary if available
            summary_path = os.path.join(session.session_dir, "summary.json")
            if os.path.exists(summary_path):
                with open(summary_path, "r") as f:
                    response["summary"] = json.load(f)

        return response

    def cancel(self, session_id: str) -> Dict[str, Any]:
        """
        Cancel a running pipeline.

        Args:
            session_id: Session identifier

        Returns:
            Cancellation status
        """
        session = self.session_manager.get_session(session_id)
        if not session:
            return {"status": "failed", "error": "Session not found"}

        session.cancel()

        return {
            "status": "cancelled",
            "session_id": session_id,
        }

    def export_fasta(
        self,
        session_id: str,
        include_review: bool = False,
    ) -> Dict[str, Any]:
        """
        Export designs as FASTA.

        Args:
            session_id: Session identifier
            include_review: Include gray zone designs

        Returns:
            FASTA content and filename
        """
        session = self.session_manager.get_session(session_id)
        if not session:
            return {"status": "failed", "error": "Session not found"}

        fasta_content = session.export_fasta(include_review=include_review)

        return {
            "status": "success",
            "filename": f"{session_id}_sequences.fasta",
            "content": fasta_content,
            "num_sequences": len([
                d for d in session.designs
                if d.status == DesignStatus.PASS or
                (include_review and d.status == DesignStatus.REVIEW)
            ]),
        }


def handle_pipeline_design(input_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Handler function for pipeline_design task.

    Args:
        input_data: Task input containing:
            - mode: "sweep" or "production"
            - metal: Metal code
            - ligand: Ligand code
            - motif_pdb: PDB content
            - sweep_configs: (optional) Custom sweep configurations
            - production_config: (for production mode) Single configuration
            - num_designs: (for production mode) Total designs
            - filters: Filter thresholds

    Returns:
        Task result with session_id and initial status
    """
    # Determine base directory
    base_dir = input_data.get(
        "base_dir",
        os.environ.get("DESIGN_HISTORY_DIR", "/tmp/design_history")
    )

    handler = PipelineHandler(base_dir)
    mode = input_data.get("mode", "sweep")

    if mode == "sweep":
        result = handler.start_sweep(
            metal=input_data.get("metal", "TB"),
            ligand=input_data.get("ligand", "CIT"),
            motif_pdb=input_data.get("motif_pdb", ""),
            sweep_configs=input_data.get("sweep_configs"),
            filters=input_data.get("filters"),
            designs_per_config=input_data.get("designs_per_config", 10),
        )
    elif mode == "production":
        result = handler.start_production(
            metal=input_data.get("metal", "TB"),
            ligand=input_data.get("ligand", "CIT"),
            motif_pdb=input_data.get("motif_pdb", ""),
            production_config=input_data.get("production_config", {}),
            num_designs=input_data.get("num_designs", 100),
            filters=input_data.get("filters"),
        )
    else:
        return {"status": "failed", "error": f"Unknown mode: {mode}"}

    return {"status": "completed", "result": result}


def handle_pipeline_status(input_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Handler function for pipeline_status task.

    Args:
        input_data: Task input containing:
            - session_id: Pipeline session identifier

    Returns:
        Current pipeline status and progress
    """
    base_dir = input_data.get(
        "base_dir",
        os.environ.get("DESIGN_HISTORY_DIR", "/tmp/design_history")
    )

    handler = PipelineHandler(base_dir)
    result = handler.get_status(input_data.get("session_id", ""))

    return {"status": "completed", "result": result}


def handle_pipeline_cancel(input_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Handler function for pipeline_cancel task.

    Args:
        input_data: Task input containing:
            - session_id: Pipeline session identifier

    Returns:
        Cancellation result
    """
    base_dir = input_data.get(
        "base_dir",
        os.environ.get("DESIGN_HISTORY_DIR", "/tmp/design_history")
    )

    handler = PipelineHandler(base_dir)
    result = handler.cancel(input_data.get("session_id", ""))

    return {"status": "completed", "result": result}


def handle_pipeline_export(input_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Handler function for pipeline_export task.

    Args:
        input_data: Task input containing:
            - session_id: Pipeline session identifier
            - include_review: Include gray zone designs

    Returns:
        FASTA export content
    """
    base_dir = input_data.get(
        "base_dir",
        os.environ.get("DESIGN_HISTORY_DIR", "/tmp/design_history")
    )

    handler = PipelineHandler(base_dir)
    result = handler.export_fasta(
        session_id=input_data.get("session_id", ""),
        include_review=input_data.get("include_review", False),
    )

    return {"status": "completed", "result": result}
