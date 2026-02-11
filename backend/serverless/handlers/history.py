"""Design history and lesson detection handlers.

Tasks: save_design_history, check_lessons
"""

import os
import json
import traceback
from typing import Dict, Any, List


def handle_save_design_history(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """Save a completed pipeline run to persistent design history."""
    try:
        from design_history import DesignHistoryManager

        history_dir = job_input.get("history_dir", "/app/design_history")
        session_name = job_input.get("session_name", "nl_pipeline")

        manager = DesignHistoryManager(history_dir)
        session = manager.start_session(session_name)

        run_id = manager.save_run(
            session=session,
            params=job_input.get("design_params", {}),
            outputs=job_input.get("design_outputs", {}),
            metrics=job_input.get("design_metrics", {}),
        )

        return {
            "status": "completed",
            "result": {
                "run_id": run_id,
                "session_id": session.session_id,
                "history_dir": history_dir,
            },
        }
    except Exception as e:
        return {
            "status": "failed",
            "error": f"Failed to save design history: {str(e)}",
        }


def handle_check_lessons(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """Check if this design run triggers any lesson patterns."""
    try:
        from lesson_detector import LessonDetector
        from design_history import DesignHistoryManager

        history_dir = job_input.get("history_dir", "/app/design_history")
        new_result = job_input.get("result", {})

        # Load recent history
        manager = DesignHistoryManager(history_dir)
        index = manager.load_index()
        history: List[Dict[str, Any]] = []

        # Load metrics from recent runs
        recent_designs = index.get("designs", [])[-20:]
        for design_entry in recent_designs:
            run_id = design_entry.get("run_id", "")
            metrics_path = os.path.join(
                history_dir, "runs", run_id, "analysis", "metrics.json"
            )
            if os.path.exists(metrics_path):
                try:
                    with open(metrics_path, "r") as f:
                        metrics = json.load(f)
                    history.append(metrics)
                except (json.JSONDecodeError, IOError):
                    pass

        detector = LessonDetector(
            failure_threshold=job_input.get("failure_threshold", 3),
            improvement_threshold=job_input.get("improvement_threshold", 0.15),
        )

        trigger = detector.check_triggers(new_result, history)

        return {
            "status": "completed",
            "result": {
                "trigger_detected": trigger is not None,
                "trigger": {
                    "type": trigger.trigger_type,
                    "description": trigger.description,
                    "relevant_designs": trigger.relevant_designs,
                    "metrics_involved": trigger.metrics_involved,
                } if trigger else None,
                "history_count": len(history),
            },
        }
    except Exception as e:
        return {
            "status": "failed",
            "error": f"Failed to check lessons: {str(e)}",
        }
