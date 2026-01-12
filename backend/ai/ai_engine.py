"""
AI-Driven Protein Engineering Engine

This module implements a hybrid AI system that:
1. Uses a comprehensive handbook as context for Claude
2. Detects task type from natural language
3. Generates RFD3 parameters dynamically
4. Provides reasoning and form configuration

Works with ANY RFdiffusion task, not just metal binding.
"""

import json
import os
import re
from typing import Dict, List, Any, Optional, Literal
from dataclasses import dataclass, asdict, field
from enum import Enum
from pathlib import Path

# Import existing chemistry rules
from .recommender import MetalChemistryRules, DesignStrategyRules, MetalProperties

# Try to import anthropic
ANTHROPIC_AVAILABLE = False
try:
    import anthropic
    ANTHROPIC_AVAILABLE = True
except ImportError:
    pass


class TaskType(str, Enum):
    """Supported RFdiffusion design tasks."""
    DE_NOVO = "de_novo"
    BINDER = "binder"
    METAL_REDESIGN = "metal_redesign"
    ENZYME = "enzyme"
    REFINEMENT = "refinement"
    SYMMETRIC = "symmetric"
    SCAFFOLD = "scaffold"
    UNKNOWN = "unknown"


@dataclass
class FormFieldConfig:
    """Configuration for a dynamic form field."""
    id: str
    label: str
    type: Literal["text", "number", "select", "range", "structure", "residue_picker"]
    required: bool = False
    default_value: Any = None
    suggested_value: Any = None
    options: Optional[List[Dict[str, str]]] = None  # For select type
    range_config: Optional[Dict[str, float]] = None  # min, max, step
    help_text: str = ""
    ai_reasoning: Optional[str] = None


@dataclass
class AIResponse:
    """Response from AI analysis."""
    success: bool
    task_type: TaskType
    params: Dict[str, Any]
    reasoning: str
    confidence: float
    clarifying_questions: List[str] = field(default_factory=list)
    form_config: List[FormFieldConfig] = field(default_factory=list)
    error: Optional[str] = None
    raw_response: Optional[str] = None


class HandbookLoader:
    """Loads and caches the AI handbook."""

    _cache: Optional[str] = None
    _cache_path: Optional[str] = None

    @classmethod
    def load(cls, handbook_path: Optional[str] = None) -> str:
        """Load the handbook from file."""
        if handbook_path is None:
            # Default path relative to this file
            current_dir = Path(__file__).parent
            handbook_path = current_dir.parent.parent / "docs" / "ai-handbook.md"

        # Check cache
        if cls._cache is not None and cls._cache_path == str(handbook_path):
            return cls._cache

        # Load from file
        try:
            with open(handbook_path, 'r', encoding='utf-8') as f:
                cls._cache = f.read()
                cls._cache_path = str(handbook_path)
                return cls._cache
        except FileNotFoundError:
            # Return minimal handbook if file not found
            return cls._get_minimal_handbook()

    @classmethod
    def _get_minimal_handbook(cls) -> str:
        """Return minimal handbook content if file not found."""
        return """
# RFdiffusion AI Handbook (Minimal)

## Supported Tasks
- de_novo: Design new proteins from scratch (use length parameter)
- metal_redesign: Redesign metal binding sites (use partial_t, ligand)
- refinement: Improve existing structures (use low partial_t)
- binder: Design protein binders (use hotspots, high step_scale)
- symmetric: Design oligomers (use symmetry config)

## Key Parameters
- partial_t: 0-30, controls noise level (higher = more change)
- num_timesteps: 50-500, quality (higher = better)
- step_scale: 0.5-3.0, diversity (higher = more diverse)
- gamma_0: 0.1-1.0, stability (higher = more stable)
"""


class AIEngine:
    """
    Hybrid AI engine using Handbook + Claude for parameter decisions.

    This engine can handle ANY RFdiffusion design task by:
    1. Parsing user intent from natural language
    2. Consulting the handbook for parameter knowledge
    3. Generating task-appropriate parameters
    4. Creating dynamic form configuration
    """

    # Task detection keywords
    TASK_KEYWORDS = {
        TaskType.DE_NOVO: ["create", "design", "generate", "new protein", "from scratch", "de novo"],
        TaskType.BINDER: ["binder", "bind to", "target", "interface", "ppi", "interaction"],
        TaskType.METAL_REDESIGN: ["metal", "lanthanide", "coordination", "tb", "gd", "eu", "zn", "fe", "ca", "terbium", "gadolinium", "europium", "zinc", "iron", "calcium"],
        TaskType.ENZYME: ["enzyme", "catalytic", "active site", "catalysis", "reaction", "substrate"],
        TaskType.REFINEMENT: ["refine", "optimize", "improve", "partial", "adjust", "tweak"],
        TaskType.SYMMETRIC: ["symmetric", "dimer", "trimer", "oligomer", "c2", "c3", "c4", "d2", "homo"],
        TaskType.SCAFFOLD: ["scaffold", "motif", "graft", "functional site", "transplant"],
    }

    def __init__(self, anthropic_api_key: Optional[str] = None):
        """
        Initialize the AI engine.

        Args:
            anthropic_api_key: API key for Claude. If not provided,
                              falls back to rule-based recommendations.
        """
        self.handbook = HandbookLoader.load()
        self.has_llm = ANTHROPIC_AVAILABLE and anthropic_api_key is not None

        if self.has_llm:
            self.client = anthropic.Anthropic(api_key=anthropic_api_key)
        else:
            self.client = None

    def detect_task_type(self, user_input: str) -> TaskType:
        """
        Detect the design task type from user input.

        Args:
            user_input: Natural language description of the task

        Returns:
            Detected TaskType
        """
        input_lower = user_input.lower()

        # Score each task type by keyword matches
        scores = {}
        for task_type, keywords in self.TASK_KEYWORDS.items():
            score = sum(1 for kw in keywords if kw in input_lower)
            if score > 0:
                scores[task_type] = score

        if not scores:
            return TaskType.UNKNOWN

        # Return highest scoring task
        return max(scores, key=scores.get)

    def analyze_request(
        self,
        user_input: str,
        structure_analysis: Optional[Dict[str, Any]] = None,
        conversation_history: Optional[List[Dict[str, str]]] = None,
        pdb_content: Optional[str] = None,
    ) -> AIResponse:
        """
        Main entry point for AI-driven parameter generation.

        Args:
            user_input: Natural language description of what user wants
            structure_analysis: Optional analysis of uploaded structure
            conversation_history: Previous messages in conversation
            pdb_content: Optional PDB file content

        Returns:
            AIResponse with parameters, reasoning, and form config
        """
        # Detect task type
        task_type = self.detect_task_type(user_input)

        # If we have LLM, use it for full analysis
        if self.has_llm:
            return self._analyze_with_llm(
                user_input, task_type, structure_analysis,
                conversation_history, pdb_content
            )

        # Fall back to rule-based analysis
        return self._analyze_with_rules(
            user_input, task_type, structure_analysis, pdb_content
        )

    def _analyze_with_llm(
        self,
        user_input: str,
        detected_task: TaskType,
        structure_analysis: Optional[Dict[str, Any]],
        conversation_history: Optional[List[Dict[str, str]]],
        pdb_content: Optional[str],
    ) -> AIResponse:
        """Use Claude API for intelligent parameter generation."""

        system_prompt = f"""You are an expert protein engineer assistant for RFdiffusion3 (RFD3).

You have access to a comprehensive handbook with all parameter knowledge:

<handbook>
{self.handbook}
</handbook>

Your role:
1. Understand the user's design intent from their natural language request
2. Detect the appropriate task type (de_novo, binder, metal_redesign, enzyme, refinement, symmetric, scaffold)
3. Generate optimal RFD3 parameters with clear reasoning
4. Ask clarifying questions if the intent is ambiguous
5. Explain your decisions in plain English

ALWAYS output a valid JSON response with this exact structure:
{{
  "task_type": "string (de_novo|binder|metal_redesign|enzyme|refinement|symmetric|scaffold)",
  "params": {{
    "contig": "optional string",
    "length": "optional string like '80-100'",
    "partial_t": "optional float 0-30",
    "num_timesteps": "optional int 50-500",
    "step_scale": "optional float 0.5-3.0",
    "gamma_0": "optional float 0.1-1.0",
    "num_designs": "optional int 1-20",
    "ligand": "optional string for metal/ligand",
    "symmetry": "optional object with id",
    "select_hotspots": "optional dict",
    "unindex": "optional string of residue IDs"
  }},
  "reasoning": "string explaining why you chose these parameters",
  "confidence": "float 0-1 indicating confidence in recommendation",
  "questions": ["array of clarifying questions if needed, empty if clear"],
  "form_fields": [
    {{
      "id": "field_id",
      "label": "Display Label",
      "type": "number|text|select|range",
      "required": true|false,
      "suggested_value": "value",
      "help_text": "explanation",
      "ai_reasoning": "why this value"
    }}
  ]
}}

Important guidelines:
- Only include parameters relevant to the detected task
- For metal_redesign, always include 'ligand' parameter
- For de_novo without structure, use 'length' not 'contig'
- For partial diffusion, 'contig' defines what to redesign
- Set realistic num_designs (5-10 for most tasks)
- Provide helpful form_fields for user adjustment
"""

        # Build user message
        user_message = f"User request: {user_input}\n\n"

        if detected_task != TaskType.UNKNOWN:
            user_message += f"Detected task type: {detected_task.value}\n\n"

        if structure_analysis:
            user_message += f"Structure analysis:\n{json.dumps(structure_analysis, indent=2)}\n\n"

        if pdb_content:
            # Extract basic info from PDB
            lines = pdb_content.split('\n')[:50]  # First 50 lines
            user_message += f"PDB excerpt (first 50 lines):\n{''.join(lines[:10])}\n...\n\n"

        if conversation_history:
            user_message += "Previous conversation:\n"
            for msg in conversation_history[-5:]:  # Last 5 messages
                user_message += f"{msg['role']}: {msg['content']}\n"

        try:
            response = self.client.messages.create(
                model="claude-sonnet-4-20250514",
                max_tokens=2048,
                system=system_prompt,
                messages=[{"role": "user", "content": user_message}]
            )

            response_text = response.content[0].text

            # Parse JSON from response
            json_match = re.search(r'\{[\s\S]*\}', response_text)
            if json_match:
                parsed = json.loads(json_match.group())

                # Convert form_fields to FormFieldConfig objects
                form_configs = []
                for field in parsed.get("form_fields", []):
                    form_configs.append(FormFieldConfig(
                        id=field.get("id", ""),
                        label=field.get("label", ""),
                        type=field.get("type", "text"),
                        required=field.get("required", False),
                        suggested_value=field.get("suggested_value"),
                        help_text=field.get("help_text", ""),
                        ai_reasoning=field.get("ai_reasoning"),
                    ))

                return AIResponse(
                    success=True,
                    task_type=TaskType(parsed.get("task_type", "unknown")),
                    params=parsed.get("params", {}),
                    reasoning=parsed.get("reasoning", ""),
                    confidence=parsed.get("confidence", 0.7),
                    clarifying_questions=parsed.get("questions", []),
                    form_config=form_configs,
                    raw_response=response_text,
                )
            else:
                return AIResponse(
                    success=False,
                    task_type=detected_task,
                    params={},
                    reasoning="Could not parse AI response",
                    confidence=0.0,
                    error="Failed to parse JSON from response",
                    raw_response=response_text,
                )

        except Exception as e:
            # Fall back to rules on LLM error
            result = self._analyze_with_rules(
                user_input, detected_task, structure_analysis, pdb_content
            )
            result.error = f"LLM failed ({str(e)}), using rule-based fallback"
            return result

    def _analyze_with_rules(
        self,
        user_input: str,
        detected_task: TaskType,
        structure_analysis: Optional[Dict[str, Any]],
        pdb_content: Optional[str],
    ) -> AIResponse:
        """Rule-based parameter generation without LLM."""

        params = {}
        reasoning_parts = []
        form_fields = []

        # Get basic structure info
        has_structure = pdb_content is not None or structure_analysis is not None

        if detected_task == TaskType.DE_NOVO:
            # De novo design
            params = {
                "length": "80-100",
                "num_designs": 10,
                "num_timesteps": 200,
                "step_scale": 1.5,
                "gamma_0": 0.6,
            }
            reasoning_parts.append("Detected de novo design task.")
            reasoning_parts.append("Using balanced parameters for new protein generation.")

            form_fields = [
                FormFieldConfig(
                    id="length", label="Protein Length", type="text",
                    required=True, suggested_value="80-100",
                    help_text="Range like '80-100' or single value",
                    ai_reasoning="Standard length for de novo design"
                ),
                FormFieldConfig(
                    id="num_designs", label="Number of Designs", type="range",
                    required=False, suggested_value=10,
                    range_config={"min": 1, "max": 20, "step": 1},
                    help_text="More designs = better sampling"
                ),
            ]

        elif detected_task == TaskType.METAL_REDESIGN:
            # Metal binding redesign
            target_metal = self._extract_metal_from_input(user_input)

            params = {
                "ligand": target_metal or "TB",
                "partial_t": 15.0,
                "num_timesteps": 200,
                "step_scale": 1.5,
                "gamma_0": 0.6,
                "num_designs": 10,
            }

            # Add contig if we have structure info
            if structure_analysis:
                num_res = structure_analysis.get("num_residues", 100)
                chain = structure_analysis.get("chains", ["A"])[0]
                params["contig"] = f"{chain}1-{num_res}"

            reasoning_parts.append(f"Detected metal binding redesign task targeting {target_metal or 'lanthanide'}.")
            reasoning_parts.append("Using moderate partial diffusion (15.0) for coordination change.")

            form_fields = [
                FormFieldConfig(
                    id="ligand", label="Target Metal", type="select",
                    required=True, suggested_value=target_metal or "TB",
                    options=[
                        {"value": "TB", "label": "Terbium (TB) - Luminescent"},
                        {"value": "GD", "label": "Gadolinium (GD) - MRI"},
                        {"value": "EU", "label": "Europium (EU) - Red emission"},
                        {"value": "CA", "label": "Calcium (CA) - Signaling"},
                        {"value": "ZN", "label": "Zinc (ZN) - Catalytic"},
                        {"value": "FE", "label": "Iron (FE) - Redox"},
                    ],
                    help_text="Metal ion to bind"
                ),
                FormFieldConfig(
                    id="partial_t", label="Redesign Level", type="range",
                    required=False, suggested_value=15.0,
                    range_config={"min": 5, "max": 25, "step": 1},
                    help_text="Higher = more aggressive redesign",
                    ai_reasoning="15 is good for moderate coordination changes"
                ),
            ]

        elif detected_task == TaskType.REFINEMENT:
            # Structure refinement
            params = {
                "partial_t": 6.0,
                "num_timesteps": 150,
                "step_scale": 1.2,
                "gamma_0": 0.8,
                "num_designs": 5,
            }

            if structure_analysis:
                num_res = structure_analysis.get("num_residues", 100)
                chain = structure_analysis.get("chains", ["A"])[0]
                params["contig"] = f"{chain}1-{num_res}"

            reasoning_parts.append("Detected refinement task.")
            reasoning_parts.append("Using low partial_t (6.0) for conservative improvements.")

            form_fields = [
                FormFieldConfig(
                    id="partial_t", label="Refinement Level", type="range",
                    required=False, suggested_value=6.0,
                    range_config={"min": 3, "max": 12, "step": 1},
                    help_text="Lower = more conservative",
                    ai_reasoning="Low noise for subtle improvements"
                ),
            ]

        elif detected_task == TaskType.BINDER:
            # Binder design
            params = {
                "step_scale": 3.0,
                "gamma_0": 0.2,
                "num_timesteps": 250,
                "num_designs": 20,
            }

            reasoning_parts.append("Detected protein binder design task.")
            reasoning_parts.append("Using high step_scale (3.0) and low gamma_0 (0.2) for diverse sampling.")

            form_fields = [
                FormFieldConfig(
                    id="binder_length", label="Binder Length", type="text",
                    required=True, suggested_value="60-80",
                    help_text="Length of binder protein to design"
                ),
                FormFieldConfig(
                    id="hotspots", label="Target Hotspots", type="text",
                    required=False, suggested_value="",
                    help_text="Residues to focus on (e.g., A45,A50,A55)"
                ),
            ]

        elif detected_task == TaskType.SYMMETRIC:
            # Symmetric design
            symmetry_type = self._extract_symmetry_from_input(user_input)

            params = {
                "length": "60-80",
                "symmetry": {"id": symmetry_type, "is_symmetric_motif": True},
                "num_timesteps": 250,
                "step_scale": 1.5,
                "num_designs": 10,
            }

            reasoning_parts.append(f"Detected symmetric oligomer design ({symmetry_type}).")

            form_fields = [
                FormFieldConfig(
                    id="symmetry", label="Symmetry Type", type="select",
                    required=True, suggested_value=symmetry_type,
                    options=[
                        {"value": "C2", "label": "C2 - Dimer"},
                        {"value": "C3", "label": "C3 - Trimer"},
                        {"value": "C4", "label": "C4 - Tetramer"},
                        {"value": "D2", "label": "D2 - 4 subunits"},
                        {"value": "D3", "label": "D3 - 6 subunits"},
                    ],
                    help_text="Symmetry group for oligomer"
                ),
                FormFieldConfig(
                    id="length", label="Subunit Length", type="text",
                    required=True, suggested_value="60-80",
                    help_text="Length of each symmetric subunit"
                ),
            ]

        else:
            # Unknown task - ask for clarification
            params = {}
            reasoning_parts.append("Could not determine the specific design task.")
            reasoning_parts.append("Please provide more details about your design goal.")

            return AIResponse(
                success=True,
                task_type=TaskType.UNKNOWN,
                params=params,
                reasoning=" ".join(reasoning_parts),
                confidence=0.3,
                clarifying_questions=[
                    "What type of protein are you trying to design?",
                    "Do you have an input structure to modify?",
                    "What is the main goal (binding, catalysis, stability)?",
                ],
                form_config=[],
            )

        return AIResponse(
            success=True,
            task_type=detected_task,
            params=params,
            reasoning=" ".join(reasoning_parts),
            confidence=0.8,
            clarifying_questions=[],
            form_config=form_fields,
        )

    def _extract_metal_from_input(self, user_input: str) -> Optional[str]:
        """Extract metal ion from user input."""
        input_lower = user_input.lower()

        metal_mapping = {
            "terbium": "TB", "tb": "TB",
            "gadolinium": "GD", "gd": "GD",
            "europium": "EU", "eu": "EU",
            "calcium": "CA", "ca": "CA",
            "zinc": "ZN", "zn": "ZN",
            "iron": "FE", "fe": "FE",
            "copper": "CU", "cu": "CU",
            "manganese": "MN", "mn": "MN",
            "magnesium": "MG", "mg": "MG",
            "cobalt": "CO", "co": "CO",
            "nickel": "NI", "ni": "NI",
            "lanthanum": "LA", "la": "LA",
        }

        for keyword, metal in metal_mapping.items():
            if keyword in input_lower:
                return metal

        return None

    def _extract_symmetry_from_input(self, user_input: str) -> str:
        """Extract symmetry type from user input."""
        input_lower = user_input.lower()

        if "dimer" in input_lower or "c2" in input_lower:
            return "C2"
        elif "trimer" in input_lower or "c3" in input_lower:
            return "C3"
        elif "tetramer" in input_lower or "c4" in input_lower:
            return "C4"
        elif "d2" in input_lower:
            return "D2"
        elif "d3" in input_lower:
            return "D3"
        elif "c5" in input_lower:
            return "C5"
        elif "c6" in input_lower:
            return "C6"

        return "C3"  # Default to trimer

    def refine_parameters(
        self,
        current_params: Dict[str, Any],
        user_feedback: str,
        task_type: TaskType,
    ) -> AIResponse:
        """
        Refine parameters based on user feedback.

        Args:
            current_params: Current parameter values
            user_feedback: User's adjustment request
            task_type: Current task type

        Returns:
            Updated AIResponse with refined parameters
        """
        feedback_lower = user_feedback.lower()
        new_params = current_params.copy()
        reasoning_parts = []

        # Adjustment patterns
        if "more aggressive" in feedback_lower or "higher" in feedback_lower:
            if "partial_t" in new_params:
                new_params["partial_t"] = min(new_params["partial_t"] + 5, 25)
                reasoning_parts.append("Increased partial_t for more aggressive redesign.")
            if "step_scale" in new_params:
                new_params["step_scale"] = min(new_params["step_scale"] + 0.5, 3.0)
                reasoning_parts.append("Increased step_scale for more diversity.")

        elif "more conservative" in feedback_lower or "lower" in feedback_lower:
            if "partial_t" in new_params:
                new_params["partial_t"] = max(new_params["partial_t"] - 5, 3)
                reasoning_parts.append("Decreased partial_t for more conservative redesign.")
            if "step_scale" in new_params:
                new_params["step_scale"] = max(new_params["step_scale"] - 0.3, 0.5)

        elif "more stable" in feedback_lower:
            new_params["gamma_0"] = min(new_params.get("gamma_0", 0.6) + 0.15, 1.0)
            new_params["num_timesteps"] = min(new_params.get("num_timesteps", 200) + 100, 500)
            reasoning_parts.append("Increased gamma_0 and timesteps for stability.")

        elif "more diverse" in feedback_lower:
            new_params["step_scale"] = min(new_params.get("step_scale", 1.5) + 0.5, 3.0)
            new_params["gamma_0"] = max(new_params.get("gamma_0", 0.6) - 0.2, 0.1)
            reasoning_parts.append("Adjusted for more diverse sampling.")

        elif "faster" in feedback_lower or "quick" in feedback_lower:
            new_params["num_timesteps"] = max(new_params.get("num_timesteps", 200) - 100, 50)
            new_params["num_designs"] = max(new_params.get("num_designs", 10) // 2, 3)
            reasoning_parts.append("Reduced timesteps and designs for faster generation.")

        elif "higher quality" in feedback_lower:
            new_params["num_timesteps"] = min(new_params.get("num_timesteps", 200) + 150, 500)
            new_params["num_designs"] = min(new_params.get("num_designs", 10) + 5, 20)
            reasoning_parts.append("Increased timesteps and designs for higher quality.")

        if not reasoning_parts:
            reasoning_parts.append("Kept parameters unchanged - please specify what to adjust.")

        return AIResponse(
            success=True,
            task_type=task_type,
            params=new_params,
            reasoning=" ".join(reasoning_parts),
            confidence=0.85,
            clarifying_questions=[],
            form_config=[],
        )


# Convenience function for API endpoint
def analyze_user_request(
    user_input: str,
    structure_analysis: Optional[Dict[str, Any]] = None,
    pdb_content: Optional[str] = None,
    conversation_history: Optional[List[Dict[str, str]]] = None,
    anthropic_api_key: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Main entry point for AI-driven analysis.

    Returns dict-serializable response for API.
    """
    engine = AIEngine(anthropic_api_key)
    response = engine.analyze_request(
        user_input=user_input,
        structure_analysis=structure_analysis,
        conversation_history=conversation_history,
        pdb_content=pdb_content,
    )

    # Convert to dict for JSON serialization
    return {
        "success": response.success,
        "task_type": response.task_type.value,
        "params": response.params,
        "reasoning": response.reasoning,
        "confidence": response.confidence,
        "clarifying_questions": response.clarifying_questions,
        "form_config": [asdict(f) for f in response.form_config],
        "error": response.error,
    }
