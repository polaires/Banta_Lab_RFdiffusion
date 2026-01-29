"""
Reporter Module (M8)

Generates human-readable reports and recommendations from pipeline results.

Reads: analysis, predictions, sequences, metadata from context
Writes: report to context
"""

import logging
from typing import Any, Dict, List, Optional

from pipeline_types import StepContext

logger = logging.getLogger(__name__)


class Reporter:
    """Generate summary reports and recommendations."""

    name = "reporter"
    description = "Generate human-readable summary report with recommendations"
    input_keys = []
    output_keys = ["report"]
    optional_keys = ["analysis", "predictions", "sequences", "design_intent",
                     "backbone_pdbs", "scaffold_result"]

    def __init__(self, **kwargs):
        self.overrides = kwargs

    def run(self, context: StepContext) -> StepContext:
        """Generate report from context data."""
        sections: List[str] = []

        # Header
        sections.append("=" * 60)
        sections.append("PROTEIN DESIGN PIPELINE REPORT")
        sections.append("=" * 60)

        # Design intent summary
        self._add_intent_section(sections, context)

        # Scaffold info
        self._add_scaffold_section(sections, context)

        # Backbone generation
        self._add_backbone_section(sections, context)

        # Sequence design
        self._add_sequence_section(sections, context)

        # Predictions and analysis
        self._add_analysis_section(sections, context)

        # Recommendations
        self._add_recommendations(sections, context)

        # Footer
        sections.append("")
        sections.append("=" * 60)

        context.report = "\n".join(sections)
        context.metadata["reporter"] = {"sections": len(sections)}
        return context

    def _add_intent_section(self, sections: List[str], context: StepContext) -> None:
        intent = context.design_intent
        if not intent:
            return
        sections.append("")
        sections.append("## Design Intent")
        metal = getattr(intent, "metal_type", None)
        ligand = getattr(intent, "ligand_name", None)
        goal = getattr(intent, "design_goal", None)
        mode = "scaffolding" if getattr(intent, "is_scaffolding", False) else "de novo"
        if metal:
            sections.append(f"  Metal: {metal}")
        if ligand:
            sections.append(f"  Ligand: {ligand}")
        if goal:
            sections.append(f"  Goal: {goal}")
        sections.append(f"  Mode: {mode}")

    def _add_scaffold_section(self, sections: List[str], context: StepContext) -> None:
        scaffold = context.scaffold_result
        if not scaffold:
            return
        sections.append("")
        sections.append("## Scaffold")
        if isinstance(scaffold, dict):
            desc = scaffold.get("pocket_description", "N/A")
            pdb_id = context.get_param("pdb_id", "N/A")
            sections.append(f"  Source PDB: {pdb_id}")
            sections.append(f"  Pocket: {desc}")

    def _add_backbone_section(self, sections: List[str], context: StepContext) -> None:
        if not context.backbone_pdbs:
            return
        sections.append("")
        sections.append("## Backbone Generation")
        sections.append(f"  Backbones generated: {len(context.backbone_pdbs)}")

    def _add_sequence_section(self, sections: List[str], context: StepContext) -> None:
        if not context.sequences:
            return
        sections.append("")
        sections.append("## Sequence Design")
        sections.append(f"  Total sequences: {len(context.sequences)}")
        if context.sequences:
            scores = [s.score for s in context.sequences if s.score != 0.0]
            if scores:
                sections.append(f"  Score range: {min(scores):.3f} to {max(scores):.3f}")

    def _add_analysis_section(self, sections: List[str], context: StepContext) -> None:
        analysis = context.analysis
        if not analysis:
            if context.predictions:
                sections.append("")
                sections.append("## Predictions")
                sections.append(f"  Total: {len(context.predictions)}")
                passed = sum(1 for p in context.predictions if p.passed)
                sections.append(f"  Passed: {passed}/{len(context.predictions)}")
            return

        sections.append("")
        sections.append("## Analysis Results")
        sections.append(f"  Filter preset: {analysis.get('filter_preset', 'N/A')}")
        sections.append(f"  Passed: {analysis.get('passed', 0)}/{analysis.get('total_predictions', 0)}")
        sections.append(f"  Pass rate: {analysis.get('pass_rate', 0):.0%}")

        # Top candidates
        candidates = analysis.get("candidates", [])
        top = [c for c in candidates if c.get("passed")][:5]
        if top:
            sections.append("")
            sections.append("  Top Candidates:")
            for c in top:
                sections.append(
                    f"    #{c['rank']}: pTM={c.get('ptm', 0):.3f}, "
                    f"pLDDT={c.get('plddt', 0):.3f}, "
                    f"tier={c.get('quality_tier', 0)}"
                )

    def _add_recommendations(self, sections: List[str], context: StepContext) -> None:
        sections.append("")
        sections.append("## Recommendations")
        recommendations: List[str] = []

        analysis = context.analysis or {}
        pass_rate = analysis.get("pass_rate", 0)
        total = analysis.get("total_predictions", len(context.predictions))

        if total == 0:
            recommendations.append("- No predictions were generated. Check backbone generation.")
        elif pass_rate == 0:
            recommendations.append("- No designs passed filters. Consider:")
            recommendations.append("  - Relaxing filter thresholds")
            recommendations.append("  - Increasing num_designs for more backbone diversity")
            recommendations.append("  - Adjusting cfg_scale or step_scale")
        elif pass_rate < 0.2:
            recommendations.append("- Low pass rate. Consider increasing num_designs.")
        elif pass_rate > 0.5:
            recommendations.append("- Good pass rate. Consider tightening filters for higher quality.")

        # Check if scout would help
        if len(context.backbone_pdbs) > 10 and pass_rate < 0.3:
            recommendations.append("- Use scout_mode to pre-filter backbones before full sequence design.")

        if not recommendations:
            recommendations.append("- Results look reasonable. Proceed with experimental validation.")

        sections.extend(recommendations)
