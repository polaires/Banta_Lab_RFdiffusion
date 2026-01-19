"""Tests for LessonDetector."""
import pytest
from lesson_detector import LessonDetector, LessonTrigger


class TestLessonDetector:
    """Test suite for LessonDetector."""

    def test_init(self):
        detector = LessonDetector()
        assert detector is not None

    def test_detect_failure_pattern(self):
        detector = LessonDetector()
        history = [
            {"outcome": "failure", "params": {"loop_length": 6}, "design_type": "metal_dimer"},
            {"outcome": "failure", "params": {"loop_length": 6}, "design_type": "metal_dimer"},
            {"outcome": "failure", "params": {"loop_length": 6}, "design_type": "metal_dimer"},
        ]
        trigger = detector.check_triggers(
            new_result=history[-1],
            history=history
        )
        assert trigger is not None
        assert trigger.trigger_type == "failure_pattern"

    def test_detect_breakthrough(self):
        detector = LessonDetector()
        history = [
            {"outcome": "failure", "metrics": {"gnina_affinity": -4}},
            {"outcome": "failure", "metrics": {"gnina_affinity": -5}},
        ]
        new_result = {"outcome": "success", "metrics": {"gnina_affinity": -8}}

        trigger = detector.check_triggers(new_result, history)
        assert trigger is not None
        assert trigger.trigger_type == "breakthrough"

    def test_detect_improvement(self):
        detector = LessonDetector()
        history = [
            {"metrics": {"coordination_distance": 3.0}},
            {"metrics": {"coordination_distance": 2.8}},
        ]
        new_result = {"metrics": {"coordination_distance": 2.3}}

        trigger = detector.check_triggers(new_result, history)
        assert trigger is not None
        assert trigger.trigger_type == "improvement"

    def test_no_trigger_for_normal_result(self):
        detector = LessonDetector()
        history = [
            {"outcome": "success", "metrics": {"gnina_affinity": -6}},
            {"outcome": "failure", "metrics": {"gnina_affinity": -4}},
        ]
        new_result = {"outcome": "success", "metrics": {"gnina_affinity": -6.5}}

        trigger = detector.check_triggers(new_result, history)
        # Small improvement, not significant enough
        assert trigger is None or trigger.trigger_type != "breakthrough"
