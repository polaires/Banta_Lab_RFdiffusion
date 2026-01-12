"""
AI Task Planner and Executor

Infrastructure for AI-driven multi-step protein engineering tasks.
Designed for text-based AI feedback (Claude API integration).

Architecture:
┌─────────────────────────────────────────────────────────────────┐
│                         AI Agent (Claude)                        │
│   - Receives natural language request                            │
│   - Analyzes and breaks into TaskPlan                           │
│   - Reviews results as text, suggests iterations                 │
└───────────────────────────┬─────────────────────────────────────┘
                            │ TaskPlan (JSON)
                            ▼
┌─────────────────────────────────────────────────────────────────┐
│                      Task Executor                               │
│   - Executes each TaskStep sequentially                         │
│   - Handles RFD3, RF3, MPNN, Analysis tasks                     │
│   - Collects results, formats for AI review                     │
└───────────────────────────┬─────────────────────────────────────┘
                            │ TaskResult (JSON/Text)
                            ▼
┌─────────────────────────────────────────────────────────────────┐
│                      Result Reporter                             │
│   - Formats results as human/AI readable text                   │
│   - Extracts key metrics for iteration decisions                │
│   - Stores artifacts (PDB, sequences) for next steps            │
└─────────────────────────────────────────────────────────────────┘

This enables Claude to:
1. Plan: "Design azobenzene binder" → [step1, step2, step3...]
2. Execute: Backend runs each step
3. Review: Claude gets text summary of results
4. Iterate: Claude suggests refinements
"""

import json
import uuid
from typing import Dict, List, Any, Optional, Literal
from dataclasses import dataclass, field, asdict
from enum import Enum
from datetime import datetime


class TaskStepType(str, Enum):
    """Types of executable task steps."""
    FETCH_PDB = "fetch_pdb"           # Fetch structure from RCSB
    ANALYZE = "analyze"                # Analyze structure
    ANALYZE_LIGAND = "analyze_ligand" # Analyze ligand binding site
    RFD3_DESIGN = "rfd3_design"       # Run RFD3 diffusion
    RF3_PREDICT = "rf3_predict"       # Run RF3 folding
    MPNN_DESIGN = "mpnn_design"       # Run MPNN sequence design
    MODIFY_PDB = "modify_pdb"         # Modify PDB (add/remove/replace atoms)
    PLACE_LIGAND = "place_ligand"     # Place ligand at interface/position
    EVALUATE = "evaluate"              # Evaluate design metrics
    COMPARE = "compare"                # Compare structures (RMSD, etc.)
    EXTRACT_LIGAND = "extract_ligand" # Extract ligand from structure
    CUSTOM = "custom"                  # Custom operation


class TaskStatus(str, Enum):
    """Status of a task or step."""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    WAITING_FOR_AI = "waiting_for_ai"  # Needs AI review


@dataclass
class TaskStep:
    """A single step in a task plan."""
    id: str
    type: TaskStepType
    name: str
    description: str
    params: Dict[str, Any]
    depends_on: List[str] = field(default_factory=list)  # IDs of prerequisite steps
    status: TaskStatus = TaskStatus.PENDING
    result: Optional[Dict[str, Any]] = None
    error: Optional[str] = None
    started_at: Optional[str] = None
    completed_at: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class TaskPlan:
    """A complete task plan with multiple steps."""
    id: str
    name: str
    description: str
    goal: str  # Natural language goal from user
    steps: List[TaskStep]
    status: TaskStatus = TaskStatus.PENDING
    created_at: str = field(default_factory=lambda: datetime.now().isoformat())
    ai_context: Optional[str] = None  # Context for AI to understand the task
    iteration: int = 0  # For multi-round refinement

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d['steps'] = [s.to_dict() if isinstance(s, TaskStep) else s for s in self.steps]
        return d

    def to_ai_summary(self) -> str:
        """Generate a text summary for AI review."""
        lines = [
            f"# Task Plan: {self.name}",
            f"**Goal:** {self.goal}",
            f"**Status:** {self.status.value}",
            f"**Iteration:** {self.iteration}",
            "",
            "## Steps:",
        ]

        for i, step in enumerate(self.steps, 1):
            status_icon = {
                "pending": "⏳",
                "running": "🔄",
                "completed": "✅",
                "failed": "❌",
                "waiting_for_ai": "🤖"
            }.get(step.status.value, "❓")

            lines.append(f"{i}. {status_icon} **{step.name}** ({step.type.value})")
            lines.append(f"   {step.description}")

            if step.result:
                # Summarize result (limit size for AI context)
                result_summary = self._summarize_result(step.result)
                lines.append(f"   Result: {result_summary}")

            if step.error:
                lines.append(f"   ⚠️ Error: {step.error}")

            lines.append("")

        return "\n".join(lines)

    def _summarize_result(self, result: Dict[str, Any], max_len: int = 200) -> str:
        """Create a brief summary of a step result."""
        parts = []

        if "designs" in result:
            count = len(result["designs"])
            parts.append(f"{count} designs generated")
        if "sequences" in result:
            count = len(result["sequences"])
            parts.append(f"{count} sequences designed")
        if "rmsd" in result:
            parts.append(f"RMSD: {result['rmsd']:.2f} Å")
        if "plddt" in result:
            parts.append(f"pLDDT: {result['plddt']:.2f}")
        if "binding_sites" in result:
            count = len(result["binding_sites"])
            parts.append(f"{count} binding sites found")
        if "ligand_code" in result:
            parts.append(f"Ligand: {result['ligand_code']}")
        if "center" in result:
            c = result["center"]
            parts.append(f"Center: ({c['x']}, {c['y']}, {c['z']})")
        if "nearby_residues" in result:
            count = len(result["nearby_residues"])
            parts.append(f"{count} nearby residues")
        if "placed" in result and result["placed"]:
            parts.append("Ligand placed successfully")
        if "extracted" in result and result["extracted"]:
            parts.append(f"Ligand extracted ({result.get('atom_count', 0)} atoms)")
        if "atom_count" in result:
            parts.append(f"{result['atom_count']} atoms")
        if "residue_count" in result:
            parts.append(f"{result['residue_count']} residues")
        if "chain_balance" in result:
            parts.append(f"Chain balance: {result['chain_balance']:.2%}")

        if parts:
            return "; ".join(parts)

        # Fallback: stringify and truncate
        s = json.dumps(result)
        if len(s) > max_len:
            s = s[:max_len] + "..."
        return s


@dataclass
class TaskResult:
    """Result of executing a task plan."""
    plan_id: str
    status: TaskStatus
    steps_completed: int
    steps_total: int
    final_artifacts: Dict[str, Any]  # PDB content, sequences, etc.
    metrics: Dict[str, Any]  # Key metrics for evaluation
    ai_summary: str  # Human/AI readable summary
    suggestions: List[str]  # Suggested next steps
    errors: List[str]

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


class TaskPlanner:
    """
    Plans complex protein engineering tasks.

    Used by AI (Claude) to break down natural language requests
    into executable TaskPlans.
    """

    # Templates for common task patterns
    TASK_TEMPLATES = {
        "metal_replacement": [
            {"type": "fetch_pdb", "name": "Fetch structure", "params": ["pdb_id"]},
            {"type": "analyze", "name": "Analyze metal site", "params": ["target_metal"]},
            {"type": "modify_pdb", "name": "Replace metal", "params": ["old_metal", "new_metal"]},
            {"type": "rfd3_design", "name": "Redesign pocket", "params": ["partial_t", "contig"]},
            {"type": "mpnn_design", "name": "Design sequences", "params": ["num_sequences"]},
            {"type": "rf3_predict", "name": "Validate folding", "params": ["sequence"]},
            {"type": "evaluate", "name": "Check metrics", "params": ["criteria"]},
        ],
        "symmetric_binder": [
            {"type": "fetch_pdb", "name": "Fetch ligand structure", "params": ["ligand_file"]},
            {"type": "rfd3_design", "name": "Generate symmetric scaffold", "params": ["symmetry", "length"]},
            {"type": "modify_pdb", "name": "Place ligand at interface", "params": ["ligand", "position"]},
            {"type": "rfd3_design", "name": "Refine binding pocket", "params": ["partial_t"]},
            {"type": "mpnn_design", "name": "Design sequences", "params": ["ligand_aware"]},
            {"type": "rf3_predict", "name": "Validate structure", "params": []},
            {"type": "evaluate", "name": "Evaluate binding", "params": ["binding_metrics"]},
        ],
        "dimerization_binder": [
            # Step 1: Analyze ligand from reference structure
            {"type": "fetch_pdb", "name": "Fetch reference with ligand", "params": ["pdb_id"]},
            {"type": "analyze_ligand", "name": "Analyze ligand binding mode", "params": ["ligand_code"]},
            {"type": "extract_ligand", "name": "Extract ligand coordinates", "params": ["ligand_code"]},
            # Step 2: Design C2 symmetric dimer scaffold
            {"type": "rfd3_design", "name": "Generate C2 symmetric dimer", "params": ["symmetry", "length"]},
            # Step 3: Place ligand at dimer interface (symmetry axis)
            {"type": "place_ligand", "name": "Place ligand at interface", "params": ["ligand_coords", "interface_center"]},
            # Step 4: Refine binding pocket around ligand
            {"type": "rfd3_design", "name": "Refine interface pocket", "params": ["partial_t", "select_buried"]},
            # Step 5: Design sequences
            {"type": "mpnn_design", "name": "Design sequences (ligand-aware)", "params": ["model_type"]},
            # Step 6: Validate
            {"type": "rf3_predict", "name": "Predict structure", "params": []},
            {"type": "evaluate", "name": "Evaluate dimerization", "params": ["binding_metrics", "symmetry_metrics"]},
        ],
        "photoswitch_binder": [
            # Specialized template for photoswitchable binders (like azobenzene)
            {"type": "fetch_pdb", "name": "Fetch structure with photoswitch", "params": ["pdb_id"]},
            {"type": "analyze_ligand", "name": "Analyze trans conformation", "params": ["ligand_code", "conformation"]},
            # Generate dimer for trans conformation
            {"type": "rfd3_design", "name": "Generate C2 dimer (trans)", "params": ["symmetry", "length"]},
            {"type": "place_ligand", "name": "Place trans conformer", "params": ["ligand_coords"]},
            {"type": "rfd3_design", "name": "Refine for trans", "params": ["partial_t"]},
            {"type": "mpnn_design", "name": "Design sequences", "params": []},
            # Evaluate both conformations
            {"type": "evaluate", "name": "Compare trans vs cis binding", "params": ["compare_conformations"]},
        ],
        "de_novo": [
            {"type": "rfd3_design", "name": "Generate structure", "params": ["length"]},
            {"type": "mpnn_design", "name": "Design sequences", "params": []},
            {"type": "rf3_predict", "name": "Predict folding", "params": []},
            {"type": "evaluate", "name": "Check designability", "params": []},
        ],
    }

    @classmethod
    def create_plan(
        cls,
        goal: str,
        template: str,
        params: Dict[str, Any],
        name: Optional[str] = None,
    ) -> TaskPlan:
        """
        Create a TaskPlan from a template.

        Args:
            goal: Natural language goal
            template: Template name (metal_replacement, symmetric_binder, etc.)
            params: Parameters for the template
            name: Optional custom name

        Returns:
            TaskPlan ready for execution
        """
        template_steps = cls.TASK_TEMPLATES.get(template, [])

        if not template_steps:
            raise ValueError(f"Unknown template: {template}")

        plan_id = f"plan-{uuid.uuid4().hex[:8]}"
        steps = []

        for i, tpl in enumerate(template_steps):
            step_id = f"step-{i+1}"
            step_params = {}

            # Fill in parameters from the provided params dict
            for param_name in tpl.get("params", []):
                if param_name in params:
                    step_params[param_name] = params[param_name]

            step = TaskStep(
                id=step_id,
                type=TaskStepType(tpl["type"]),
                name=tpl["name"],
                description=f"Step {i+1}: {tpl['name']}",
                params=step_params,
                depends_on=[f"step-{i}"] if i > 0 else [],
            )
            steps.append(step)

        return TaskPlan(
            id=plan_id,
            name=name or f"{template} task",
            description=f"Task plan for: {goal}",
            goal=goal,
            steps=steps,
        )

    @classmethod
    def create_custom_plan(
        cls,
        goal: str,
        steps: List[Dict[str, Any]],
        name: Optional[str] = None,
    ) -> TaskPlan:
        """
        Create a custom TaskPlan from a list of step definitions.

        This is used when the AI generates a custom plan for complex tasks.

        Args:
            goal: Natural language goal
            steps: List of step definitions
            name: Optional custom name

        Returns:
            TaskPlan ready for execution
        """
        plan_id = f"plan-{uuid.uuid4().hex[:8]}"
        task_steps = []

        for i, step_def in enumerate(steps):
            step_id = step_def.get("id", f"step-{i+1}")

            step = TaskStep(
                id=step_id,
                type=TaskStepType(step_def.get("type", "custom")),
                name=step_def.get("name", f"Step {i+1}"),
                description=step_def.get("description", ""),
                params=step_def.get("params", {}),
                depends_on=step_def.get("depends_on", [f"step-{i}"] if i > 0 else []),
            )
            task_steps.append(step)

        return TaskPlan(
            id=plan_id,
            name=name or "Custom task",
            description=f"Custom task plan for: {goal}",
            goal=goal,
            steps=task_steps,
        )


class TaskExecutor:
    """
    Executes TaskPlans step by step.

    This class interfaces with the actual inference functions
    and collects results for AI review.
    """

    def __init__(self):
        self.artifacts = {}  # Store intermediate results

    async def execute_plan(
        self,
        plan: TaskPlan,
        progress_callback: Optional[callable] = None
    ) -> TaskResult:
        """
        Execute a task plan.

        Args:
            plan: The TaskPlan to execute
            progress_callback: Optional callback for progress updates

        Returns:
            TaskResult with execution summary
        """
        plan.status = TaskStatus.RUNNING
        errors = []
        metrics = {}

        for step in plan.steps:
            # Check dependencies
            deps_satisfied = all(
                self._get_step_by_id(plan, dep_id).status == TaskStatus.COMPLETED
                for dep_id in step.depends_on
            )

            if not deps_satisfied:
                step.status = TaskStatus.PENDING
                continue

            step.status = TaskStatus.RUNNING
            step.started_at = datetime.now().isoformat()

            if progress_callback:
                progress_callback(step.id, step.name, "running")

            try:
                result = await self._execute_step(step, plan)
                step.result = result
                step.status = TaskStatus.COMPLETED
                step.completed_at = datetime.now().isoformat()

                # Extract metrics
                if "rmsd" in result:
                    metrics["rmsd"] = result["rmsd"]
                if "plddt" in result:
                    metrics["plddt"] = result["plddt"]

                if progress_callback:
                    progress_callback(step.id, step.name, "completed")

            except Exception as e:
                step.status = TaskStatus.FAILED
                step.error = str(e)
                step.completed_at = datetime.now().isoformat()
                errors.append(f"{step.name}: {str(e)}")

                if progress_callback:
                    progress_callback(step.id, step.name, "failed")

        # Determine overall status
        if all(s.status == TaskStatus.COMPLETED for s in plan.steps):
            plan.status = TaskStatus.COMPLETED
        elif any(s.status == TaskStatus.FAILED for s in plan.steps):
            plan.status = TaskStatus.FAILED
        else:
            plan.status = TaskStatus.PENDING

        # Generate suggestions based on results
        suggestions = self._generate_suggestions(plan, metrics)

        return TaskResult(
            plan_id=plan.id,
            status=plan.status,
            steps_completed=sum(1 for s in plan.steps if s.status == TaskStatus.COMPLETED),
            steps_total=len(plan.steps),
            final_artifacts=self.artifacts,
            metrics=metrics,
            ai_summary=plan.to_ai_summary(),
            suggestions=suggestions,
            errors=errors,
        )

    async def _execute_step(self, step: TaskStep, plan: TaskPlan) -> Dict[str, Any]:
        """Execute a single step."""
        # Import inference functions
        from inference_utils import (
            run_rfd3_inference,
            run_rf3_inference,
            run_mpnn_inference,
            analyze_structure,
        )

        if step.type == TaskStepType.FETCH_PDB:
            pdb_id = step.params.get("pdb_id")
            import urllib.request
            url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
            with urllib.request.urlopen(url) as response:
                content = response.read().decode()
            self.artifacts["pdb_content"] = content
            return {"pdb_id": pdb_id, "size": len(content)}

        elif step.type == TaskStepType.ANALYZE:
            pdb_content = self.artifacts.get("pdb_content") or step.params.get("pdb_content")
            result = analyze_structure(pdb_content)
            return result.get("result", result)

        elif step.type == TaskStepType.RFD3_DESIGN:
            result = run_rfd3_inference(**step.params)
            if result.get("status") == "completed" and result.get("result", {}).get("designs"):
                self.artifacts["designs"] = result["result"]["designs"]
                self.artifacts["pdb_content"] = result["result"]["designs"][0]["content"]
            return result.get("result", result)

        elif step.type == TaskStepType.RF3_PREDICT:
            sequence = step.params.get("sequence") or self.artifacts.get("sequence")
            result = run_rf3_inference(sequence=sequence)
            return result.get("result", result)

        elif step.type == TaskStepType.MPNN_DESIGN:
            pdb_content = self.artifacts.get("pdb_content") or step.params.get("pdb_content")
            result = run_mpnn_inference(
                pdb_content=pdb_content,
                **{k: v for k, v in step.params.items() if k != "pdb_content"}
            )
            if result.get("result", {}).get("sequences"):
                self.artifacts["sequences"] = result["result"]["sequences"]
                self.artifacts["sequence"] = result["result"]["sequences"][0]["sequence"]
            return result.get("result", result)

        elif step.type == TaskStepType.MODIFY_PDB:
            # Handle PDB modifications (metal replacement, ligand placement, etc.)
            return await self._modify_pdb(step.params)

        elif step.type == TaskStepType.ANALYZE_LIGAND:
            return await self._analyze_ligand(step.params)

        elif step.type == TaskStepType.EXTRACT_LIGAND:
            return await self._extract_ligand(step.params)

        elif step.type == TaskStepType.PLACE_LIGAND:
            return await self._place_ligand(step.params)

        elif step.type == TaskStepType.EVALUATE:
            return await self._evaluate(step.params)

        elif step.type == TaskStepType.COMPARE:
            return await self._compare_structures(step.params)

        else:
            raise ValueError(f"Unknown step type: {step.type}")

    async def _modify_pdb(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Modify PDB content (metal replacement, ligand addition, etc.)."""
        pdb_content = self.artifacts.get("pdb_content", "")

        if "old_metal" in params and "new_metal" in params:
            old = params["old_metal"]
            new = params["new_metal"]
            # Simple text replacement for metal
            modified = pdb_content.replace(f" {old} ", f" {new} ")
            modified = modified.replace(f" {old}\n", f" {new}\n")
            self.artifacts["pdb_content"] = modified
            return {"modified": True, "old_metal": old, "new_metal": new}

        return {"modified": False, "reason": "Unknown modification type"}

    async def _analyze_ligand(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Analyze ligand binding site and properties."""
        pdb_content = self.artifacts.get("pdb_content", "")
        ligand_code = params.get("ligand_code", "")

        lines = pdb_content.split("\n")
        hetatm_lines = [l for l in lines if l.startswith("HETATM")]
        ligand_atoms = [l for l in hetatm_lines if ligand_code in l[17:20]]

        if not ligand_atoms:
            return {"found": False, "ligand_code": ligand_code, "error": "Ligand not found"}

        # Parse ligand coordinates
        coords = []
        atom_info = []
        for line in ligand_atoms:
            try:
                atom_name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append((x, y, z))
                atom_info.append({"name": atom_name, "x": x, "y": y, "z": z})
            except (ValueError, IndexError):
                continue

        if not coords:
            return {"found": False, "error": "Could not parse ligand coordinates"}

        # Calculate center of mass
        center_x = sum(c[0] for c in coords) / len(coords)
        center_y = sum(c[1] for c in coords) / len(coords)
        center_z = sum(c[2] for c in coords) / len(coords)

        # Find nearby protein residues (within 6 Angstroms)
        atom_lines = [l for l in lines if l.startswith("ATOM")]
        nearby_residues = set()
        for line in atom_lines:
            try:
                ax = float(line[30:38])
                ay = float(line[38:46])
                az = float(line[46:54])
                dist = ((ax - center_x)**2 + (ay - center_y)**2 + (az - center_z)**2) ** 0.5
                if dist < 6.0:
                    res_name = line[17:20].strip()
                    res_num = line[22:26].strip()
                    chain = line[21]
                    nearby_residues.add(f"{chain}:{res_name}{res_num}")
            except (ValueError, IndexError):
                continue

        result = {
            "found": True,
            "ligand_code": ligand_code,
            "atom_count": len(coords),
            "center": {"x": round(center_x, 2), "y": round(center_y, 2), "z": round(center_z, 2)},
            "atoms": atom_info[:10],  # First 10 atoms
            "nearby_residues": list(nearby_residues)[:15],  # First 15 residues
        }

        # Store for later use
        self.artifacts["ligand_info"] = result
        self.artifacts["ligand_center"] = (center_x, center_y, center_z)

        return result

    async def _extract_ligand(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Extract ligand HETATM records from structure."""
        pdb_content = self.artifacts.get("pdb_content", "")
        ligand_code = params.get("ligand_code", "")

        lines = pdb_content.split("\n")
        ligand_lines = [l for l in lines if l.startswith("HETATM") and ligand_code in l[17:20]]

        if not ligand_lines:
            return {"extracted": False, "error": f"Ligand {ligand_code} not found"}

        ligand_pdb = "\n".join(ligand_lines)
        self.artifacts["ligand_pdb"] = ligand_pdb

        return {
            "extracted": True,
            "ligand_code": ligand_code,
            "atom_count": len(ligand_lines),
            "pdb_lines": len(ligand_lines),
        }

    async def _place_ligand(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Place ligand at specified position (e.g., dimer interface center)."""
        ligand_pdb = self.artifacts.get("ligand_pdb", "")
        pdb_content = self.artifacts.get("pdb_content", "")

        if not ligand_pdb:
            return {"placed": False, "error": "No ligand extracted"}

        # Get target position (default: center of current structure)
        target_pos = params.get("position")
        if not target_pos:
            # Calculate center of the protein structure
            lines = pdb_content.split("\n")
            atom_lines = [l for l in lines if l.startswith("ATOM")]
            coords = []
            for line in atom_lines:
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
                except (ValueError, IndexError):
                    continue
            if coords:
                target_pos = {
                    "x": sum(c[0] for c in coords) / len(coords),
                    "y": sum(c[1] for c in coords) / len(coords),
                    "z": sum(c[2] for c in coords) / len(coords),
                }
            else:
                target_pos = {"x": 0, "y": 0, "z": 0}

        # Get current ligand center
        ligand_center = self.artifacts.get("ligand_center")
        if not ligand_center:
            # Calculate from ligand PDB
            ligand_lines = ligand_pdb.split("\n")
            coords = []
            for line in ligand_lines:
                if line.startswith("HETATM"):
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append((x, y, z))
                    except (ValueError, IndexError):
                        continue
            if coords:
                ligand_center = (
                    sum(c[0] for c in coords) / len(coords),
                    sum(c[1] for c in coords) / len(coords),
                    sum(c[2] for c in coords) / len(coords),
                )
            else:
                return {"placed": False, "error": "Could not calculate ligand center"}

        # Calculate translation
        dx = target_pos["x"] - ligand_center[0]
        dy = target_pos["y"] - ligand_center[1]
        dz = target_pos["z"] - ligand_center[2]

        # Translate ligand atoms
        translated_lines = []
        for line in ligand_pdb.split("\n"):
            if line.startswith("HETATM"):
                try:
                    x = float(line[30:38]) + dx
                    y = float(line[38:46]) + dy
                    z = float(line[46:54]) + dz
                    new_line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
                    translated_lines.append(new_line)
                except (ValueError, IndexError):
                    translated_lines.append(line)
            else:
                translated_lines.append(line)

        translated_ligand = "\n".join(translated_lines)

        # Combine protein with translated ligand
        combined = pdb_content.rstrip()
        if not combined.endswith("\n"):
            combined += "\n"
        combined += translated_ligand

        self.artifacts["pdb_content"] = combined
        self.artifacts["ligand_placed"] = True

        return {
            "placed": True,
            "target_position": target_pos,
            "translation": {"dx": round(dx, 2), "dy": round(dy, 2), "dz": round(dz, 2)},
        }

    async def _compare_structures(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Compare two structures (RMSD, etc.)."""
        pdb1 = params.get("pdb1") or self.artifacts.get("pdb_content")
        pdb2 = params.get("pdb2")

        if not pdb1 or not pdb2:
            return {"compared": False, "error": "Need two structures to compare"}

        # Simple atom count comparison
        lines1 = [l for l in pdb1.split("\n") if l.startswith("ATOM")]
        lines2 = [l for l in pdb2.split("\n") if l.startswith("ATOM")]

        return {
            "compared": True,
            "structure1_atoms": len(lines1),
            "structure2_atoms": len(lines2),
            "atom_difference": abs(len(lines1) - len(lines2)),
        }

    async def _evaluate(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Evaluate design metrics."""
        metrics = {}

        if self.artifacts.get("pdb_content"):
            pdb_content = self.artifacts["pdb_content"]
            lines = pdb_content.split("\n")
            atoms = [l for l in lines if l.startswith("ATOM")]
            hetatms = [l for l in lines if l.startswith("HETATM")]
            metrics["atom_count"] = len(atoms)
            metrics["hetatm_count"] = len(hetatms)

            # Count residues
            residues = set()
            for line in atoms:
                res_num = line[22:26].strip()
                chain = line[21]
                residues.add(f"{chain}_{res_num}")
            metrics["residue_count"] = len(residues)

            # Check for ligand
            if self.artifacts.get("ligand_placed"):
                metrics["ligand_placed"] = True

        if self.artifacts.get("sequences"):
            metrics["sequence_count"] = len(self.artifacts["sequences"])
            # Sequence length
            if self.artifacts["sequences"]:
                first_seq = self.artifacts["sequences"][0]
                if isinstance(first_seq, dict):
                    metrics["sequence_length"] = len(first_seq.get("sequence", ""))
                else:
                    metrics["sequence_length"] = len(first_seq)

        # Binding metrics evaluation
        if "binding_metrics" in params:
            # Check if ligand is properly positioned at interface
            if self.artifacts.get("ligand_info"):
                ligand_info = self.artifacts["ligand_info"]
                metrics["ligand_center"] = ligand_info.get("center")
                metrics["nearby_residue_count"] = len(ligand_info.get("nearby_residues", []))

        # Symmetry metrics evaluation
        if "symmetry_metrics" in params:
            pdb_content = self.artifacts.get("pdb_content", "")
            lines = pdb_content.split("\n")
            atom_lines = [l for l in lines if l.startswith("ATOM")]
            # Check chain distribution (C2 should have equal atoms per chain)
            chain_counts = {}
            for line in atom_lines:
                chain = line[21]
                chain_counts[chain] = chain_counts.get(chain, 0) + 1
            metrics["chains"] = chain_counts
            if len(chain_counts) == 2:
                counts = list(chain_counts.values())
                metrics["chain_balance"] = min(counts) / max(counts) if max(counts) > 0 else 0

        return metrics

    def _get_step_by_id(self, plan: TaskPlan, step_id: str) -> Optional[TaskStep]:
        """Get a step by ID."""
        for step in plan.steps:
            if step.id == step_id:
                return step
        return None

    def _generate_suggestions(self, plan: TaskPlan, metrics: Dict[str, Any]) -> List[str]:
        """Generate suggestions based on results."""
        suggestions = []

        if metrics.get("rmsd", float("inf")) > 3.0:
            suggestions.append("High RMSD detected. Consider reducing partial_t or increasing num_timesteps.")

        if metrics.get("plddt", 0) < 0.7:
            suggestions.append("Low pLDDT score. The design may need refinement.")

        if not self.artifacts.get("sequences"):
            suggestions.append("No sequences generated yet. Run MPNN to design sequences.")

        if plan.status == TaskStatus.COMPLETED:
            suggestions.append("Task completed. Review the results and iterate if needed.")

        return suggestions


# Convenience function for API endpoint
def plan_task(
    goal: str,
    template: Optional[str] = None,
    params: Optional[Dict[str, Any]] = None,
    custom_steps: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    """
    Plan a protein engineering task.

    This is the main entry point for AI-driven task planning.

    Args:
        goal: Natural language goal
        template: Optional template name
        params: Parameters for the template
        custom_steps: Custom step definitions (if not using template)

    Returns:
        TaskPlan as dictionary
    """
    if custom_steps:
        plan = TaskPlanner.create_custom_plan(goal, custom_steps)
    elif template:
        plan = TaskPlanner.create_plan(goal, template, params or {})
    else:
        # AI should call with either template or custom_steps
        raise ValueError("Either template or custom_steps must be provided")

    return plan.to_dict()


def plan_dimerization_binder(
    ligand_pdb_id: str,
    ligand_code: str,
    monomer_length: str = "60-80",
    symmetry: str = "C2",
    partial_t: int = 12,
) -> Dict[str, Any]:
    """
    Create a plan for designing a dimerization binder.

    This is specialized for designing proteins that:
    1. Form symmetric dimers (C2, D2, etc.)
    2. Bind a small molecule at the dimer interface
    3. Use the ligand to mediate/stabilize dimerization

    Args:
        ligand_pdb_id: PDB ID containing the ligand (e.g., "3R1V" for azobenzene)
        ligand_code: Three-letter ligand code (e.g., "AZB")
        monomer_length: Length of each monomer (e.g., "60-80")
        symmetry: Symmetry type (C2, C3, D2, etc.)
        partial_t: Noise level for pocket refinement

    Returns:
        TaskPlan as dictionary
    """
    params = {
        "pdb_id": ligand_pdb_id,
        "ligand_code": ligand_code,
        "length": monomer_length,
        "symmetry": {"id": symmetry},
        "partial_t": partial_t,
        "model_type": "ligand_mpnn",
        "binding_metrics": True,
        "symmetry_metrics": True,
    }

    return plan_task(
        goal=f"Design {symmetry} symmetric dimer that binds {ligand_code} at the interface",
        template="dimerization_binder",
        params=params,
    )


def plan_photoswitch_binder(
    photoswitch_pdb_id: str,
    photoswitch_code: str,
    monomer_length: str = "60-80",
) -> Dict[str, Any]:
    """
    Create a plan for designing a photoswitch-responsive binder.

    Specialized for photoswitchable ligands like azobenzene that undergo
    trans↔cis isomerization upon light exposure.

    Args:
        photoswitch_pdb_id: PDB ID containing the photoswitch (e.g., "3R1V")
        photoswitch_code: Three-letter code (e.g., "AZB" for azobenzene)
        monomer_length: Length of each monomer

    Returns:
        TaskPlan as dictionary
    """
    params = {
        "pdb_id": photoswitch_pdb_id,
        "ligand_code": photoswitch_code,
        "conformation": "trans",  # Design for trans, then test cis
        "length": monomer_length,
        "symmetry": {"id": "C2"},
        "partial_t": 12,
        "compare_conformations": True,
    }

    return plan_task(
        goal=f"Design C2 dimer responsive to {photoswitch_code} photoswitching",
        template="photoswitch_binder",
        params=params,
    )
