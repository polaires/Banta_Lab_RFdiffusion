"""
Pipeline Session Manager

Extends DesignHistoryManager for pipeline-specific session management.
Handles parameter sweeps, production runs, and intelligent design filtering.
"""
import os
import json
import shutil
from datetime import datetime, timedelta
from typing import Dict, Any, Optional, List
from dataclasses import dataclass, field, asdict
from enum import Enum

from design_history import DesignHistoryManager, DesignSession


class PipelineMode(str, Enum):
    """Pipeline execution modes."""
    SINGLE = "single"
    SWEEP = "sweep"
    PRODUCTION = "production"


class DesignStatus(str, Enum):
    """Design filtering status."""
    PASS = "pass"
    REVIEW = "review"  # Gray zone - needs review
    FAIL = "fail"


@dataclass
class FilterThresholds:
    """Filter thresholds for design evaluation."""
    plddt: float = 0.80
    ptm: float = 0.80
    pae: float = 5.0

    # Discard threshold - don't even store designs below this
    plddt_discard: float = 0.60

    def to_dict(self) -> Dict[str, float]:
        return asdict(self)

    @classmethod
    def strict(cls) -> "FilterThresholds":
        """Strict filtering preset."""
        return cls(plddt=0.80, ptm=0.80, pae=5.0)

    @classmethod
    def relaxed(cls) -> "FilterThresholds":
        """Relaxed filtering preset."""
        return cls(plddt=0.70, ptm=0.65, pae=10.0)


@dataclass
class SweepConfig:
    """Configuration for a single sweep iteration."""
    name: str
    contig_size: str  # "small", "medium", "large"
    contig_range: str  # e.g., "50-70"
    cfg_scale: float  # 1.5, 2.0, 2.5
    num_designs: int = 10

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class PipelineDesign:
    """A single design result from the pipeline."""
    name: str
    config_name: str
    sequence: str
    pdb_content: str
    plddt: float
    ptm: float
    pae: float
    status: DesignStatus
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())

    # Optional metrics
    coordination_number: Optional[int] = None
    geometry_rmsd: Optional[float] = None
    interface_area: Optional[float] = None

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d["status"] = self.status.value
        return d

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "PipelineDesign":
        data["status"] = DesignStatus(data["status"])
        return cls(**data)


@dataclass
class PipelineProgress:
    """Real-time progress tracking for pipeline execution."""
    current_config: int = 0
    total_configs: int = 0
    current_design: int = 0
    designs_per_config: int = 10
    total_generated: int = 0
    total_passing: int = 0
    total_review: int = 0
    total_failed: int = 0
    best_design: Optional[Dict[str, Any]] = None

    @property
    def pass_rate(self) -> float:
        if self.total_generated == 0:
            return 0.0
        return self.total_passing / self.total_generated

    def to_dict(self) -> Dict[str, Any]:
        return {
            "current_config": self.current_config,
            "total_configs": self.total_configs,
            "current_design": self.current_design,
            "designs_per_config": self.designs_per_config,
            "total_generated": self.total_generated,
            "total_passing": self.total_passing,
            "total_review": self.total_review,
            "total_failed": self.total_failed,
            "pass_rate": self.pass_rate,
            "best_design": self.best_design,
        }


class PipelineSession:
    """
    Manages a single pipeline execution session.

    Directory structure:
        pipeline_sessions/{session_id}/
            meta.json           - Session metadata
            progress.json       - Current progress state
            configs/            - Sweep configurations
            results/            - All design results
            passing_designs/    - Designs that pass filters
            review_designs/     - Designs in gray zone
            summary.json        - Final summary
    """

    def __init__(
        self,
        session_dir: str,
        session_id: str,
        mode: PipelineMode,
        metal: str,
        ligand: str,
        filters: FilterThresholds,
    ):
        self.session_dir = session_dir
        self.session_id = session_id
        self.mode = mode
        self.metal = metal
        self.ligand = ligand
        self.filters = filters
        self.progress = PipelineProgress()
        self.designs: List[PipelineDesign] = []
        self.cancelled = False

        self._ensure_structure()
        self._save_metadata()

    def _ensure_structure(self):
        """Create session directory structure."""
        subdirs = ["configs", "results", "passing_designs", "review_designs"]
        for subdir in subdirs:
            os.makedirs(os.path.join(self.session_dir, subdir), exist_ok=True)

    def _save_metadata(self):
        """Save session metadata."""
        meta = {
            "session_id": self.session_id,
            "mode": self.mode.value,
            "metal": self.metal,
            "ligand": self.ligand,
            "filters": self.filters.to_dict(),
            "started": datetime.now().isoformat(),
            "status": "running",
        }
        meta_path = os.path.join(self.session_dir, "meta.json")
        with open(meta_path, "w") as f:
            json.dump(meta, f, indent=2)

    def save_progress(self):
        """Save current progress state."""
        progress_path = os.path.join(self.session_dir, "progress.json")
        with open(progress_path, "w") as f:
            json.dump(self.progress.to_dict(), f, indent=2)

    def load_progress(self) -> PipelineProgress:
        """Load progress state (for resuming)."""
        progress_path = os.path.join(self.session_dir, "progress.json")
        if os.path.exists(progress_path):
            with open(progress_path, "r") as f:
                data = json.load(f)
            self.progress = PipelineProgress(**{
                k: v for k, v in data.items()
                if k in PipelineProgress.__dataclass_fields__
            })
        return self.progress

    def save_config(self, config: SweepConfig, index: int):
        """Save a sweep configuration."""
        config_path = os.path.join(
            self.session_dir, "configs", f"{index:02d}_{config.name}.json"
        )
        with open(config_path, "w") as f:
            json.dump(config.to_dict(), f, indent=2)

    def evaluate_design(
        self,
        plddt: float,
        ptm: float,
        pae: float,
    ) -> DesignStatus:
        """
        Evaluate a design against filter thresholds.

        Returns:
            DesignStatus indicating pass/review/fail
        """
        # Discard immediately if below discard threshold
        if plddt < self.filters.plddt_discard:
            return DesignStatus.FAIL

        # Check strict thresholds
        passes_strict = (
            plddt >= self.filters.plddt and
            ptm >= self.filters.ptm and
            pae <= self.filters.pae
        )

        if passes_strict:
            return DesignStatus.PASS

        # Check relaxed thresholds for gray zone
        relaxed = FilterThresholds.relaxed()
        passes_relaxed = (
            plddt >= relaxed.plddt and
            ptm >= relaxed.ptm and
            pae <= relaxed.pae
        )

        if passes_relaxed:
            return DesignStatus.REVIEW

        return DesignStatus.FAIL

    def add_design(self, design: PipelineDesign) -> bool:
        """
        Add a design to the session.

        Returns:
            True if design was stored (pass or review), False if discarded
        """
        # Update progress
        self.progress.total_generated += 1

        if design.status == DesignStatus.FAIL:
            self.progress.total_failed += 1
            # Don't store failed designs, just count them
            return False

        # Store the design
        self.designs.append(design)

        if design.status == DesignStatus.PASS:
            self.progress.total_passing += 1
            subdir = "passing_designs"
        else:  # REVIEW
            self.progress.total_review += 1
            subdir = "review_designs"

        # Save PDB
        pdb_path = os.path.join(
            self.session_dir, subdir, f"{design.name}.pdb"
        )
        with open(pdb_path, "w") as f:
            f.write(design.pdb_content)

        # Save metrics
        metrics_path = os.path.join(
            self.session_dir, "results", f"{design.name}.json"
        )
        with open(metrics_path, "w") as f:
            # Don't store full PDB in JSON
            design_dict = design.to_dict()
            design_dict.pop("pdb_content", None)
            json.dump(design_dict, f, indent=2)

        # Update best design
        if (
            self.progress.best_design is None or
            design.plddt > self.progress.best_design.get("plddt", 0)
        ):
            self.progress.best_design = {
                "name": design.name,
                "plddt": design.plddt,
                "ptm": design.ptm,
                "pae": design.pae,
            }

        # Save progress
        self.save_progress()

        return True

    def get_results(self) -> List[Dict[str, Any]]:
        """Get all design results (without full PDB content)."""
        results = []
        results_dir = os.path.join(self.session_dir, "results")

        for filename in sorted(os.listdir(results_dir)):
            if filename.endswith(".json"):
                filepath = os.path.join(results_dir, filename)
                with open(filepath, "r") as f:
                    results.append(json.load(f))

        return results

    def get_passing_designs(self) -> List[Dict[str, Any]]:
        """Get only passing designs."""
        return [
            r for r in self.get_results()
            if r.get("status") == DesignStatus.PASS.value
        ]

    def export_fasta(self, include_review: bool = False) -> str:
        """
        Export passing sequences as FASTA.

        Args:
            include_review: Include gray zone designs

        Returns:
            FASTA formatted string
        """
        fasta_lines = []

        for design in self.designs:
            if design.status == DesignStatus.PASS:
                include = True
            elif design.status == DesignStatus.REVIEW and include_review:
                include = True
            else:
                include = False

            if include:
                header = (
                    f">{design.name} "
                    f"pLDDT={design.plddt:.2f} "
                    f"pTM={design.ptm:.2f} "
                    f"PAE={design.pae:.1f} "
                    f"config={design.config_name}"
                )
                fasta_lines.append(header)
                # Wrap sequence at 80 characters
                seq = design.sequence
                for i in range(0, len(seq), 80):
                    fasta_lines.append(seq[i:i+80])

        return "\n".join(fasta_lines)

    def finalize(self, status: str = "completed"):
        """
        Finalize the session and generate summary.

        Args:
            status: Final status (completed, cancelled, failed)
        """
        summary = {
            "session_id": self.session_id,
            "mode": self.mode.value,
            "status": status,
            "completed": datetime.now().isoformat(),
            "filters": self.filters.to_dict(),
            "statistics": {
                "total_generated": self.progress.total_generated,
                "total_passing": self.progress.total_passing,
                "total_review": self.progress.total_review,
                "total_failed": self.progress.total_failed,
                "pass_rate": self.progress.pass_rate,
            },
            "best_design": self.progress.best_design,
        }

        # Update metadata
        meta_path = os.path.join(self.session_dir, "meta.json")
        with open(meta_path, "r") as f:
            meta = json.load(f)
        meta["status"] = status
        meta["completed"] = summary["completed"]
        with open(meta_path, "w") as f:
            json.dump(meta, f, indent=2)

        # Save summary
        summary_path = os.path.join(self.session_dir, "summary.json")
        with open(summary_path, "w") as f:
            json.dump(summary, f, indent=2)

        # Save FASTA export
        fasta_content = self.export_fasta(include_review=True)
        if fasta_content:
            fasta_path = os.path.join(self.session_dir, "passing_sequences.fasta")
            with open(fasta_path, "w") as f:
                f.write(fasta_content)

        return summary

    def cancel(self):
        """Cancel the pipeline execution."""
        self.cancelled = True
        self.finalize(status="cancelled")


class PipelineSessionManager:
    """
    Manages multiple pipeline sessions.

    Directory structure:
        {base_dir}/
            pipeline_sessions/
                {session_id}/
                    ...
    """

    def __init__(self, base_dir: str):
        self.base_dir = base_dir
        self.sessions_dir = os.path.join(base_dir, "pipeline_sessions")
        os.makedirs(self.sessions_dir, exist_ok=True)

    def create_session(
        self,
        mode: PipelineMode,
        metal: str,
        ligand: str,
        filters: Optional[FilterThresholds] = None,
    ) -> PipelineSession:
        """
        Create a new pipeline session.

        Args:
            mode: Pipeline execution mode
            metal: Metal code (e.g., "TB")
            ligand: Ligand code (e.g., "CIT")
            filters: Filter thresholds (defaults to strict)

        Returns:
            New PipelineSession
        """
        if filters is None:
            filters = FilterThresholds.strict()

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        session_id = f"pipeline_{timestamp}_{metal}_{ligand}"
        session_dir = os.path.join(self.sessions_dir, session_id)

        return PipelineSession(
            session_dir=session_dir,
            session_id=session_id,
            mode=mode,
            metal=metal,
            ligand=ligand,
            filters=filters,
        )

    def get_session(self, session_id: str) -> Optional[PipelineSession]:
        """
        Get an existing session by ID.

        Args:
            session_id: Session identifier

        Returns:
            PipelineSession if exists, None otherwise
        """
        session_dir = os.path.join(self.sessions_dir, session_id)
        meta_path = os.path.join(session_dir, "meta.json")

        if not os.path.exists(meta_path):
            return None

        with open(meta_path, "r") as f:
            meta = json.load(f)

        filters = FilterThresholds(**meta.get("filters", {}))

        session = PipelineSession(
            session_dir=session_dir,
            session_id=session_id,
            mode=PipelineMode(meta["mode"]),
            metal=meta["metal"],
            ligand=meta["ligand"],
            filters=filters,
        )

        # Load existing progress
        session.load_progress()

        return session

    def list_sessions(self) -> List[Dict[str, Any]]:
        """List all pipeline sessions with basic info."""
        sessions = []

        for session_id in os.listdir(self.sessions_dir):
            meta_path = os.path.join(
                self.sessions_dir, session_id, "meta.json"
            )
            if os.path.exists(meta_path):
                with open(meta_path, "r") as f:
                    meta = json.load(f)
                sessions.append({
                    "session_id": session_id,
                    "mode": meta.get("mode"),
                    "status": meta.get("status"),
                    "started": meta.get("started"),
                    "metal": meta.get("metal"),
                    "ligand": meta.get("ligand"),
                })

        # Sort by started time descending
        sessions.sort(key=lambda x: x.get("started", ""), reverse=True)

        return sessions

    def cleanup_old_sessions(self, max_age_days: int = 7):
        """
        Clean up old sessions based on retention policy.

        - Failed designs: Already not stored
        - Review designs: Delete after max_age_days unless bookmarked
        - Passing designs: Keep indefinitely

        Args:
            max_age_days: Maximum age for review designs
        """
        cutoff = datetime.now() - timedelta(days=max_age_days)

        for session_id in os.listdir(self.sessions_dir):
            session_dir = os.path.join(self.sessions_dir, session_id)
            meta_path = os.path.join(session_dir, "meta.json")

            if not os.path.exists(meta_path):
                continue

            with open(meta_path, "r") as f:
                meta = json.load(f)

            # Check if session is old enough
            started = datetime.fromisoformat(meta.get("started", ""))
            if started > cutoff:
                continue

            # Delete review designs
            review_dir = os.path.join(session_dir, "review_designs")
            if os.path.exists(review_dir):
                shutil.rmtree(review_dir)
                os.makedirs(review_dir)  # Recreate empty

            # Update results to remove review designs
            results_dir = os.path.join(session_dir, "results")
            if os.path.exists(results_dir):
                for filename in os.listdir(results_dir):
                    filepath = os.path.join(results_dir, filename)
                    with open(filepath, "r") as f:
                        result = json.load(f)
                    if result.get("status") == DesignStatus.REVIEW.value:
                        os.remove(filepath)


def generate_sweep_configs(
    metal: str = "TB",
    ligand: str = "CIT",
    designs_per_config: int = 10,
) -> List[SweepConfig]:
    """
    Generate standard sweep configurations (3 sizes x 3 CFG values = 9 configs).

    Args:
        metal: Metal code
        ligand: Ligand code
        designs_per_config: Number of designs per configuration

    Returns:
        List of SweepConfig objects
    """
    sizes = {
        "small": "50-70",
        "medium": "70-90",
        "large": "90-120",
    }

    cfg_scales = [1.5, 2.0, 2.5]

    configs = []
    for size_name, contig_range in sizes.items():
        for cfg in cfg_scales:
            cfg_name = f"low_cfg" if cfg == 1.5 else f"mid_cfg" if cfg == 2.0 else "high_cfg"
            name = f"{size_name}_{cfg_name}"
            configs.append(SweepConfig(
                name=name,
                contig_size=size_name,
                contig_range=contig_range,
                cfg_scale=cfg,
                num_designs=designs_per_config,
            ))

    return configs
