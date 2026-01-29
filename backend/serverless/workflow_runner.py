"""
Workflow Runner (M9)

JSON workflow spec interpreter and step executor.
The AI assistant can compose workflows either:
  1. As Python: runner.run([step1, step2, ...], ctx)
  2. As JSON: WorkflowRunner.from_spec(json_spec)

Supports optional iteration strategies (scout, sweep).
"""

import json
import logging
import time
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

from pipeline_types import StepContext, PipelineStep
from inference_backend import create_backend

logger = logging.getLogger(__name__)


class WorkflowRunner:
    """Executes a sequence of pipeline modules with optional strategies."""

    def __init__(
        self,
        strategy: Optional[str] = None,
        strategy_params: Optional[Dict[str, Any]] = None,
        checkpoint_dir: Optional[Path] = None,
    ):
        """
        Args:
            strategy: "scout", "sweep", or None (plain sequential)
            strategy_params: Parameters for the strategy (e.g., ptm_threshold)
            checkpoint_dir: If set, save context after each step
        """
        self.strategy = strategy
        self.strategy_params = strategy_params or {}
        self.checkpoint_dir = checkpoint_dir

    def run(
        self,
        steps: List[Any],
        context: StepContext,
        on_step_complete: Optional[Callable] = None,
    ) -> StepContext:
        """
        Execute steps sequentially. Returns final context.

        Args:
            steps: List of PipelineStep instances
            context: Initial StepContext
            on_step_complete: Optional callback(step_name, context) after each step

        Returns:
            Final StepContext with all results
        """
        total_start = time.time()
        context.metadata["workflow_start"] = total_start
        context.metadata["workflow_steps"] = [s.name for s in steps]
        context.metadata["step_timings"] = {}
        context.metadata["step_status"] = {}

        # Apply strategy wrapping if configured
        if self.strategy:
            steps = self._apply_strategy(steps, context)

        for i, step in enumerate(steps):
            step_name = step.name
            logger.info(f"[{i + 1}/{len(steps)}] Running: {step_name}")
            context.metadata["current_step"] = step_name
            context.metadata["current_step_index"] = i

            step_start = time.time()
            try:
                context = step.run(context)
                elapsed = time.time() - step_start
                context.metadata["step_timings"][step_name] = elapsed
                context.metadata["step_status"][step_name] = "completed"
                logger.info(f"[{i + 1}/{len(steps)}] {step_name} completed ({elapsed:.1f}s)")
            except Exception as e:
                elapsed = time.time() - step_start
                context.metadata["step_timings"][step_name] = elapsed
                context.metadata["step_status"][step_name] = "failed"
                context.metadata["step_error"] = str(e)
                logger.error(f"[{i + 1}/{len(steps)}] {step_name} failed: {e}")
                raise

            # Checkpoint
            if self.checkpoint_dir:
                self._save_checkpoint(step_name, context)

            # Callback
            if on_step_complete:
                on_step_complete(step_name, context)

        context.metadata["workflow_elapsed"] = time.time() - total_start
        return context

    @classmethod
    def from_spec(
        cls,
        spec: Dict[str, Any],
    ) -> Tuple["WorkflowRunner", List[Any], StepContext]:
        """
        Parse a JSON workflow spec into runner + steps + context.

        This is the AI assistant's primary interface:
            spec = AI generates JSON
            runner, steps, ctx = WorkflowRunner.from_spec(spec)
            result = runner.run(steps, ctx)

        Args:
            spec: JSON workflow specification dict

        Returns:
            Tuple of (WorkflowRunner, steps list, StepContext)
        """
        from pipeline_modules import MODULE_REGISTRY

        # Create backend
        backend_mode = spec.get("backend", "auto")
        backend = create_backend(backend_mode)

        # Create context
        ctx = StepContext(
            backend=backend,
            params=spec.get("params", {}),
            metadata={"workflow_name": spec.get("name", "unnamed")},
        )

        # Working directory
        working_dir = spec.get("working_dir")
        if working_dir:
            ctx.working_dir = Path(working_dir)

        # Pre-populated data (e.g., backbone_pdbs passed directly)
        if "backbone_pdbs" in spec.get("params", {}):
            ctx.backbone_pdbs = spec["params"].pop("backbone_pdbs")

        # Build step instances from spec
        steps = []
        for step_spec in spec.get("steps", []):
            module_name = step_spec["module"]
            if module_name not in MODULE_REGISTRY:
                raise ValueError(
                    f"Unknown module: {module_name!r}. "
                    f"Available: {list(MODULE_REGISTRY.keys())}"
                )
            module_cls = MODULE_REGISTRY[module_name]
            step_params = step_spec.get("params", {})
            step = module_cls(**step_params)
            steps.append(step)

        # Create runner with strategy
        strategy_spec = spec.get("strategy", {})
        strategy = strategy_spec.get("type") if strategy_spec else None
        runner = cls(
            strategy=strategy,
            strategy_params=strategy_spec,
            checkpoint_dir=Path(spec["checkpoint_dir"]) if spec.get("checkpoint_dir") else None,
        )

        return runner, steps, ctx

    def _apply_strategy(self, steps: List[Any], context: StepContext) -> List[Any]:
        """Wrap steps with iteration strategy if configured."""
        if self.strategy == "scout":
            return self._apply_scout_strategy(steps, context)
        elif self.strategy == "sweep":
            return self._apply_sweep_strategy(steps, context)
        return steps

    def _apply_scout_strategy(self, steps: List[Any], context: StepContext) -> List[Any]:
        """
        Scout strategy: insert a scout round between backbone_generator and
        sequence_designer to filter bad backbones early.
        """
        try:
            from iteration_strategies import ScoutStrategy
            ptm_threshold = self.strategy_params.get("ptm_threshold", 0.6)
            scout = ScoutStrategy(ptm_threshold=ptm_threshold)

            # Find the insertion point
            new_steps = []
            for step in steps:
                new_steps.append(step)
                if step.name == "backbone_generator":
                    new_steps.append(scout)
            return new_steps
        except ImportError:
            logger.warning("Scout strategy requested but iteration_strategies not available")
            return steps

    def _apply_sweep_strategy(self, steps: List[Any], context: StepContext) -> List[Any]:
        """
        Sweep strategy: wrap the core pipeline in a parameter sweep loop.
        """
        try:
            from iteration_strategies import SweepStrategy
            grid = self.strategy_params.get("grid", {})
            designs_per_config = self.strategy_params.get("designs_per_config", 10)
            sweep = SweepStrategy(grid=grid, designs_per_config=designs_per_config)

            # Insert sweep wrapper before backbone_generator
            new_steps = []
            for step in steps:
                if step.name == "backbone_generator":
                    new_steps.append(sweep)
                new_steps.append(step)
            return new_steps
        except ImportError:
            logger.warning("Sweep strategy requested but iteration_strategies not available")
            return steps

    def _save_checkpoint(self, step_name: str, context: StepContext) -> None:
        """Save context checkpoint after a step."""
        if not self.checkpoint_dir:
            return
        try:
            self.checkpoint_dir.mkdir(parents=True, exist_ok=True)
            checkpoint_path = self.checkpoint_dir / f"checkpoint_{step_name}.json"
            with open(checkpoint_path, "w") as f:
                json.dump(context.to_dict(), f, indent=2, default=str)
            logger.debug(f"Checkpoint saved: {checkpoint_path}")
        except Exception as e:
            logger.warning(f"Failed to save checkpoint for {step_name}: {e}")

    def get_progress(self, context: StepContext) -> Dict[str, Any]:
        """
        Get workflow progress information (for frontend polling).

        Returns a WorkflowProgress-compatible dict.
        """
        steps = context.metadata.get("workflow_steps", [])
        step_status = context.metadata.get("step_status", {})
        current_step = context.metadata.get("current_step", "")
        current_idx = context.metadata.get("current_step_index", 0)

        step_results = {}
        for step_name, status in step_status.items():
            step_results[step_name] = {
                "status": status,
                "summary": self._step_summary(step_name, context),
                "timing": context.metadata.get("step_timings", {}).get(step_name),
            }

        all_done = all(step_status.get(s) == "completed" for s in steps)
        any_failed = any(step_status.get(s) == "failed" for s in steps)

        return {
            "status": "failed" if any_failed else ("completed" if all_done else "running"),
            "workflow_name": context.metadata.get("workflow_name", ""),
            "current_step": current_step,
            "current_step_index": current_idx,
            "total_steps": len(steps),
            "step_results": step_results,
        }

    def _step_summary(self, step_name: str, context: StepContext) -> str:
        """Generate a short summary string for a completed step."""
        if step_name == "backbone_generator":
            return f"Generated {len(context.backbone_pdbs)} backbone(s)"
        elif step_name == "sequence_designer":
            return f"Designed {len(context.sequences)} sequence(s)"
        elif step_name == "structure_predictor":
            passed = sum(1 for p in context.predictions if p.passed)
            return f"Predicted {len(context.predictions)} structures, {passed} passed"
        elif step_name == "analyzer":
            analysis = context.analysis or {}
            return f"Pass rate: {analysis.get('pass_rate', 0):.0%}"
        elif step_name == "reporter":
            return "Report generated"
        elif step_name == "intent_parser":
            intent = context.design_intent
            if intent:
                return f"Parsed: {getattr(intent, 'metal_type', '')} + {getattr(intent, 'ligand_name', '')}"
            return "Intent parsed"
        elif step_name == "design_configurator":
            return "Configs generated"
        elif step_name == "scaffolder":
            return f"Scaffold from {context.get_param('pdb_id', 'PDB')}"
        return "Completed"


def run_workflow_spec(spec: Dict[str, Any], on_step_complete: Optional[Callable] = None) -> Dict[str, Any]:
    """
    Top-level convenience function for handler.py integration.

    Args:
        spec: JSON workflow specification
        on_step_complete: Optional callback per step

    Returns:
        Context as dict (for API response)
    """
    runner, steps, ctx = WorkflowRunner.from_spec(spec)
    result_ctx = runner.run(steps, ctx, on_step_complete=on_step_complete)
    return result_ctx.to_dict()
