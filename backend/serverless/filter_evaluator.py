"""
Filter Evaluator

Evaluates design metrics against filter presets (BindCraft-style).
"""
import os
import json
from typing import Dict, Any, List, Optional


# Centralized filter presets — threshold dicts for common design types.
# These can be used directly by the Analyzer module (M7) without
# requiring JSON files on disk.
#
# Each metric uses {"min": X} (value must be >= X) or {"max": X} (value must be <= X).
# "min" is for quality scores (ptm, plddt, iptm) — higher is better.
# "max" is for distance/error metrics (rmsd, pae) — lower is better.
FILTER_PRESETS = {
    "ppi_binder": {
        "inter_chain_pae": {"max": 1.5},
        "ptm": {"min": 0.8},
        "ca_rmsd": {"max": 2.5},
    },
    "small_molecule": {
        "backbone_rmsd": {"max": 1.5},
        "ligand_rmsd": {"max": 5.0},
        "iptm": {"min": 0.8},
    },
    "metal_binding": {
        "coordination_number": {"min": 6},
        "geometry_rmsd": {"max": 2.0},
        "plddt": {"min": 0.70},
    },
    "metal_binding_strict": {
        "coordination_number": {"min": 8},
        "geometry_rmsd": {"max": 0.8},
        "plddt": {"min": 0.85},
    },
    "dna_binder": {
        "dna_aligned_rmsd": {"max": 5.0},
    },
    "enzyme_scaffold": {
        "motif_all_atom_rmsd": {"max": 1.5},
    },
    "scout_relaxed": {
        "ptm": {"min": 0.6},
        "plddt": {"min": 0.65},
    },
    "production_strict": {
        "ptm": {"min": 0.8},
        "plddt": {"min": 0.75},
        "backbone_rmsd": {"max": 1.5},
    },
}


class FilterEvaluator:
    """
    Evaluates metrics against filter presets.

    Usage:
        evaluator = FilterEvaluator("experiments/design_history/filter_presets")
        result = evaluator.evaluate(metrics, "default")
    """

    def __init__(self, presets_dir: str):
        """
        Initialize evaluator.

        Args:
            presets_dir: Path to filter_presets directory
        """
        self.presets_dir = presets_dir
        self._presets_cache: Dict[str, Dict] = {}

    def load_preset(self, name: str) -> Dict[str, Any]:
        """
        Load a filter preset by name.

        Args:
            name: Preset name (without .json extension)

        Returns:
            Preset dictionary
        """
        if name in self._presets_cache:
            return self._presets_cache[name]

        preset_path = os.path.join(self.presets_dir, f"{name}.json")
        with open(preset_path, "r") as f:
            preset = json.load(f)

        self._presets_cache[name] = preset
        return preset

    def list_presets(self) -> List[str]:
        """List available preset names."""
        presets = []
        for filename in os.listdir(self.presets_dir):
            if filename.endswith(".json"):
                presets.append(filename[:-5])
        return presets

    def evaluate(
        self,
        metrics: Dict[str, Any],
        preset_name: str,
    ) -> Dict[str, Any]:
        """
        Evaluate metrics against a filter preset.

        Args:
            metrics: Metrics from UnifiedDesignAnalyzer
            preset_name: Name of filter preset to use

        Returns:
            Evaluation result with pass/fail and failed filters
        """
        preset = self.load_preset(preset_name)
        failed_filters = []

        for analysis_name, filters in preset.get("filters", {}).items():
            analysis_data = metrics.get("analyses", {}).get(analysis_name, {})

            # Skip if analysis was not applicable or skipped
            if analysis_data.get("status") in ["not_applicable", "skipped"]:
                continue

            # Check each filter
            analysis_metrics = analysis_data.get("metrics", {})
            for metric_name, threshold in filters.items():
                value = analysis_metrics.get(metric_name)

                if value is None:
                    continue

                # Check threshold
                passed = self._check_threshold(value, threshold)
                if not passed:
                    failed_filters.append({
                        "analysis": analysis_name,
                        "metric": metric_name,
                        "value": value,
                        "threshold": threshold,
                    })

        return {
            "preset": preset_name,
            "pass": len(failed_filters) == 0,
            "failed_filters": failed_filters,
        }

    def _check_threshold(
        self,
        value: Any,
        threshold: Dict[str, Any],
    ) -> bool:
        """
        Check if value passes threshold.

        Args:
            value: Metric value
            threshold: Threshold definition

        Returns:
            True if passes, False if fails
        """
        if "min" in threshold:
            if value < threshold["min"]:
                return False

        if "max" in threshold:
            if value > threshold["max"]:
                return False

        if "equals" in threshold:
            if value != threshold["equals"]:
                return False

        return True

    def evaluate_all_presets(
        self,
        metrics: Dict[str, Any],
    ) -> Dict[str, Dict[str, Any]]:
        """
        Evaluate metrics against all available presets.

        Args:
            metrics: Metrics from UnifiedDesignAnalyzer

        Returns:
            Dictionary of preset_name -> evaluation result
        """
        results = {}
        for preset_name in self.list_presets():
            try:
                results[preset_name] = self.evaluate(metrics, preset_name)
            except Exception as e:
                results[preset_name] = {
                    "preset": preset_name,
                    "pass": False,
                    "error": str(e),
                }
        return results
