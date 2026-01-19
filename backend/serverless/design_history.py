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

# Default filter name used for filter evaluation
DEFAULT_FILTER = "default"


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
        try:
            with open(index_path, "r") as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError) as e:
            # Return empty index if corrupted/missing
            return {"version": "1.0.0", "created": "", "designs": []}

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
            "filter_pass": metrics.get("filter_results", {}).get(DEFAULT_FILTER, {}).get("pass"),
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
                try:
                    with open(metrics_path, "r") as f:
                        metrics = json.load(f)
                        if metrics.get("filter_results", {}).get(DEFAULT_FILTER, {}).get("pass"):
                            passing += 1
                except (json.JSONDecodeError, IOError) as e:
                    # Skip corrupted metrics files
                    pass

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
                try:
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
                except (json.JSONDecodeError, IOError) as e:
                    # Skip corrupted metrics files
                    pass

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
