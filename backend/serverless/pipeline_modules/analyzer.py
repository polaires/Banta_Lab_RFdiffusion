"""
Analyzer Module (M7)

Runs comprehensive analysis and filtering on design predictions.
Wraps UnifiedDesignAnalyzer + FilterEvaluator — does NOT rewrite them.

Reads: predictions (and backbone_pdbs) from context
Writes: analysis to context
"""

import logging
from typing import Any, Dict, List, Optional

from pipeline_types import StepContext
from filter_evaluator import FILTER_PRESETS

logger = logging.getLogger(__name__)


class Analyzer:
    """Run comprehensive analysis and filtering on design predictions."""

    name = "analyzer"
    description = "Analyze predicted structures with metrics, RMSD, and quality filtering"
    input_keys = ["predictions"]
    output_keys = ["analysis"]
    optional_keys = ["backbone_pdbs", "design_intent", "resolved_ligand"]

    def __init__(self, filter_preset: str = "auto", **kwargs):
        """
        Args:
            filter_preset: Name of filter preset or "auto" to detect from design type.
        """
        self.filter_preset = filter_preset
        self.overrides = kwargs

    @classmethod
    def with_preset(cls, preset_name: str) -> "Analyzer":
        """Factory: create Analyzer with a specific filter preset."""
        return cls(filter_preset=preset_name)

    def run(self, context: StepContext) -> StepContext:
        """Run analysis on all predictions."""
        predictions = context.predictions
        if not predictions:
            logger.warning("No predictions to analyze")
            context.analysis = {"status": "no_predictions", "candidates": []}
            return context

        logger.info(f"Analyzing {len(predictions)} prediction(s)")

        # Determine filter preset
        preset_name = self._resolve_preset(context)
        thresholds = FILTER_PRESETS.get(preset_name, FILTER_PRESETS.get("scout_relaxed", {}))

        # Run analysis per prediction
        candidates: List[Dict[str, Any]] = []
        passed_count = 0

        for i, pred in enumerate(predictions):
            candidate = self._analyze_one(context, pred, thresholds, rank=i + 1)
            candidates.append(candidate)
            if candidate.get("passed"):
                passed_count += 1

        # Sort by quality tier then ptm
        candidates.sort(key=lambda c: (-c.get("quality_tier", 0), -c.get("ptm", 0)))

        # Re-rank after sorting
        for i, c in enumerate(candidates):
            c["rank"] = i + 1

        pass_rate = passed_count / len(predictions) if predictions else 0.0

        context.analysis = {
            "status": "completed",
            "filter_preset": preset_name,
            "thresholds": thresholds,
            "total_predictions": len(predictions),
            "passed": passed_count,
            "pass_rate": pass_rate,
            "candidates": candidates,
            "best_candidate": candidates[0] if candidates else None,
        }

        logger.info(f"Analysis complete: {passed_count}/{len(predictions)} passed ({pass_rate:.0%})")

        context.metadata["analyzer"] = {
            "filter_preset": preset_name,
            "passed": passed_count,
            "total": len(predictions),
            "pass_rate": pass_rate,
        }
        return context

    def _resolve_preset(self, context: StepContext) -> str:
        """Determine which filter preset to use."""
        if self.filter_preset != "auto":
            return self.filter_preset

        preset = context.get_param("filter_preset")
        if preset:
            return preset

        # Auto-detect from design intent
        intent = context.design_intent
        if intent:
            metal = getattr(intent, "metal_type", None)
            goal = getattr(intent, "design_goal", None)
            is_scaffold = getattr(intent, "is_scaffolding", False)

            if is_scaffold:
                return "enzyme_scaffold"
            if metal:
                return "metal_binding"
            if goal == "ppi":
                return "ppi_binder"

        return "scout_relaxed"

    def _analyze_one(
        self,
        context: StepContext,
        pred: Any,
        thresholds: Dict[str, Any],
        rank: int,
    ) -> Dict[str, Any]:
        """Analyze a single prediction."""
        result: Dict[str, Any] = {
            "rank": rank,
            "sequence": pred.sequence[:50] + "..." if len(pred.sequence) > 50 else pred.sequence,
            "full_sequence": pred.sequence,
            "plddt": pred.plddt,
            "ptm": pred.ptm,
            "iptm": pred.iptm,
            "rmsd": pred.rmsd,
            "predictor": pred.predictor,
            "backbone_index": pred.backbone_index,
        }

        # Determine pass/fail against thresholds
        # Thresholds use {"min": X} (value >= X) or {"max": X} (value <= X)
        passed = True
        failed_filters: List[str] = []

        for metric, threshold in thresholds.items():
            value = getattr(pred, metric, None)
            if value is None:
                value = result.get(metric)
            if value is None:
                continue

            if isinstance(threshold, dict):
                if "min" in threshold and value < threshold["min"]:
                    passed = False
                    failed_filters.append(f"{metric}={value:.3f} < min {threshold['min']}")
                if "max" in threshold and value > threshold["max"]:
                    passed = False
                    failed_filters.append(f"{metric}={value:.3f} > max {threshold['max']}")
            elif isinstance(threshold, (int, float)):
                # Legacy flat value — treat as min threshold
                if value < threshold:
                    passed = False
                    failed_filters.append(f"{metric}={value:.3f} < {threshold}")

        result["passed"] = passed
        result["failed_filters"] = failed_filters

        # Quality tier (0-3)
        result["quality_tier"] = self._quality_tier(pred)

        # Run unified analyzer if available
        if pred.predicted_pdb:
            result["has_structure"] = True
            try:
                detailed = self._run_unified_analysis(context, pred)
                if detailed:
                    result["detailed_analysis"] = detailed
            except Exception as e:
                logger.debug(f"Unified analysis failed: {e}")
        else:
            result["has_structure"] = False

        return result

    def _quality_tier(self, pred: Any) -> int:
        """Assign quality tier 0-3 based on confidence metrics."""
        ptm = pred.ptm
        plddt = pred.plddt
        if ptm >= 0.85 and plddt >= 0.80:
            return 3  # Excellent
        elif ptm >= 0.7 and plddt >= 0.70:
            return 2  # Good
        elif ptm >= 0.6 and plddt >= 0.60:
            return 1  # Marginal
        return 0  # Poor

    def _run_unified_analysis(self, context: StepContext, pred: Any) -> Optional[Dict[str, Any]]:
        """Run UnifiedDesignAnalyzer if available."""
        try:
            from unified_analyzer import UnifiedDesignAnalyzer
            analyzer = UnifiedDesignAnalyzer()

            # Get reference PDB for RMSD
            ref_pdb = None
            if context.backbone_pdbs and pred.backbone_index < len(context.backbone_pdbs):
                ref_pdb = context.backbone_pdbs[pred.backbone_index]

            return analyzer.analyze(
                predicted_pdb=pred.predicted_pdb,
                reference_pdb=ref_pdb,
                sequence=pred.sequence,
            )
        except ImportError:
            return None
        except Exception as e:
            logger.debug(f"Unified analysis error: {e}")
            return None
