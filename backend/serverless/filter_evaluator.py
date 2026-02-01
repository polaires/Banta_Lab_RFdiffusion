"""
Filter Evaluator

Evaluates design metrics against filter presets (BindCraft-style).
Includes chemistry-aware preset generation from metal_chemistry.py data.
"""
import os
import json
from typing import Dict, Any, List, Optional, Tuple


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


# Tier multipliers for chemistry-aware preset generation.
# cn_offset: added to the empirical CN base (formal_cn_min - 2)
# rmsd_base: geometry RMSD threshold for CN 4-6 metals, scaled up for CN 7+ metals
_TIER_PARAMS = {
    "relaxed": {"cn_offset": -1, "rmsd_base": 2.5, "rmsd_high_cn": 3.0, "plddt": 0.60, "ptm": 0.50, "lc": 2},
    "standard": {"cn_offset": 0, "rmsd_base": 1.5, "rmsd_high_cn": 2.0, "plddt": 0.70, "ptm": 0.60, "lc": 3},
    "strict": {"cn_offset": 1, "rmsd_base": 0.8, "rmsd_high_cn": 1.2, "plddt": 0.80, "ptm": 0.75, "lc": 5},
}


def generate_chemistry_aware_preset(
    metal: Optional[str] = None,
    ligand_name: Optional[str] = None,
    tier: str = "standard",
    design_type: str = "metal_binding",
) -> Tuple[str, Dict[str, Any]]:
    """
    Generate filter thresholds from metal_chemistry.py data.

    For metals in METAL_DATABASE, derives CN thresholds from formal coordination
    range, adjusted by empirical observation that monomer designs typically achieve
    formal_cn_min - 2 total coordination (protein + ligand combined).

    Args:
        metal: Metal element code (TB, ZN, CA, etc.)
        ligand_name: Ligand residue name (CIT, PQQ, etc.) — adds ligand_contacts filter
        tier: Quality tier — "relaxed", "standard", or "strict"
        design_type: Fallback design type if metal not in database

    Returns:
        Tuple of (preset_name, preset_dict) where preset_dict uses
        {"min": X} / {"max": X} format compatible with FILTER_PRESETS.
    """
    if tier not in _TIER_PARAMS:
        tier = "standard"
    tp = _TIER_PARAMS[tier]

    # If no metal or metal not known, fall back to existing hardcoded presets
    if not metal:
        fallback = _fallback_preset_name(design_type)
        return fallback, dict(FILTER_PRESETS.get(fallback, {}))

    try:
        from metal_chemistry import (
            METAL_DATABASE,
            get_coordination_number_range,
            get_hsab_class_simple,
        )
    except ImportError:
        fallback = _fallback_preset_name(design_type)
        return fallback, dict(FILTER_PRESETS.get(fallback, {}))

    metal_upper = metal.upper()
    if metal_upper not in METAL_DATABASE:
        fallback = _fallback_preset_name(design_type)
        return fallback, dict(FILTER_PRESETS.get(fallback, {}))

    metal_data = METAL_DATABASE[metal_upper]
    default_ox = metal_data.get("default_oxidation", 2)

    try:
        formal_cn_min, formal_cn_max = get_coordination_number_range(metal_upper, default_ox)
    except (ValueError, KeyError):
        formal_cn_min, formal_cn_max = 4, 6

    try:
        hsab_class = get_hsab_class_simple(metal_upper)
    except (ValueError, KeyError):
        hsab_class = "borderline"

    # Empirical CN base: formal_cn_min - 2 for monomers
    # (protein provides ~3-5 donors, ligand provides ~3-4, total rarely reaches formal max)
    cn_base = max(3, formal_cn_min - 2)
    cn_threshold = max(3, cn_base + tp["cn_offset"])

    # Geometry RMSD: higher CN metals have more flexible coordination spheres
    if formal_cn_min >= 7:
        rmsd_threshold = tp["rmsd_high_cn"]
    else:
        rmsd_threshold = tp["rmsd_base"]

    # Build preset
    preset: Dict[str, Any] = {
        "coordination_number": {"min": cn_threshold},
        "plddt": {"min": tp["plddt"]},
        "ptm": {"min": tp["ptm"]},
    }

    # Only add geometry_rmsd if it's commonly computed
    if rmsd_threshold is not None:
        preset["geometry_rmsd"] = {"max": rmsd_threshold}

    # Add ligand_contacts filter when ligand is specified
    if ligand_name:
        preset["ligand_contacts"] = {"min": tp["lc"]}

    # Build descriptive name
    name = f"metal_{metal_upper.lower()}_{tier}"

    return name, preset


def _fallback_preset_name(design_type: str) -> str:
    """Map design_type to a hardcoded FILTER_PRESETS key."""
    if design_type in ("metal", "metal_ligand", "metal_monomer", "metal_ligand_monomer", "metal_binding"):
        return "metal_binding"
    if design_type in ("enzyme", "scaffold", "enzyme_scaffold"):
        return "enzyme_scaffold"
    return "scout_relaxed"


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
