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
