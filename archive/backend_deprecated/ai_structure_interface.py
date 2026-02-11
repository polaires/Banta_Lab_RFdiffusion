# ai_structure_interface.py
"""
AI-Callable Structure Interface

Provides natural language interface for AI/LLM agents (like Claude) to query
metal coordination chemistry, find protein structures, and get design recommendations.

This module bridges natural language queries to the underlying structure discovery
and metal chemistry databases.

Example:
    >>> answer = answer_structure_question("What is the coordination geometry of terbium?")
    >>> plan = plan_structure_search("I need to design a terbium biosensor")
    >>> recs = get_design_recommendations("TB", task="biosensor")
"""

import re
import logging
from typing import Dict, List, Optional, Any, Tuple

# Import structure discovery service
from structure_discovery import StructureDiscovery, METAL_NAMES, LIGAND_NAMES

# Import metal chemistry database
from metal_chemistry import (
    METAL_DATABASE,
    COORDINATION_MOTIFS,
    get_coordination_template,
    get_preferred_donors,
    get_coordination_number_range,
    get_amino_acid_bias,
    get_common_residues,
    get_hsab_class,
    is_lanthanide,
)

# Import metal-ligand templates
from metal_ligand_templates import (
    METAL_LIGAND_COMPLEX_TEMPLATES,
    get_template,
    get_template_info,
    recommend_template,
    list_templates,
)

logger = logging.getLogger(__name__)


# =============================================================================
# CONSTANTS FOR NL PARSING
# =============================================================================

# Question type patterns - ORDER MATTERS for detection priority!
# Find structure patterns (checked FIRST)
FIND_STRUCTURE_KEYWORDS = [
    "find", "search", "get me", "show me", "retrieve", "look for",
    "structures", "pdbs", "examples", "list",
]

# Template patterns (checked SECOND)
TEMPLATE_KEYWORDS = [
    "template", "starting point", "scaffold", "starting structure",
    "reference structure", "basis structure",
]

# Donor patterns (checked THIRD)
DONOR_KEYWORDS = [
    "what residues", "which amino acids", "preferred donors",
    "what amino acid", "which residue", "coordinating residues",
    "amino acids coordinate", "residues coordinate",
]

# Coordination patterns (checked LAST - broadest)
COORDINATION_KEYWORDS = [
    "coordination", "geometry", "coordinate", "coordination number",
    "cn of", "how many ligands", "coordination sphere",
]

DESIGN_KEYWORDS = [
    "design", "engineer", "create", "build", "make", "develop",
    "biosensor", "sensor", "binding site", "active site",
]


# =============================================================================
# MAIN PUBLIC FUNCTIONS
# =============================================================================

def answer_structure_question(question: str) -> Dict[str, Any]:
    """
    Main entry point for natural language questions about metal coordination.

    Detects the question type and routes to appropriate handler.

    Args:
        question: Natural language question about metal coordination chemistry.

    Returns:
        Dict with:
            - answer: Human-readable answer string
            - type: Question type detected ("coordination", "find_structure",
                    "template", "donors", "unknown")
            - Additional structured data based on question type

    Example:
        >>> result = answer_structure_question("What is the coordination geometry of terbium?")
        >>> print(result["answer"])
        "Terbium typically has coordination number 8-9 with geometry..."
        >>> print(result["coordination_number"])
        [8, 9]
    """
    question_lower = question.lower()

    # Detect question type - ORDER MATTERS!
    # Check most specific patterns first, broadest (coordination) last
    if _is_find_structure_question(question_lower):
        return _answer_find_structure_question(question)
    elif _is_template_question(question_lower):
        return _answer_template_question(question)
    elif _is_donor_question(question_lower):
        return _answer_donor_question(question)
    elif _is_coordination_question(question_lower):
        return _answer_coordination_question(question)
    else:
        # Try to detect metal and give general info
        metal = _extract_metal_from_question(question)
        if metal:
            return _answer_general_metal_question(question, metal)

        return {
            "answer": "I couldn't determine what you're asking about. "
                      "Try asking about coordination geometry, finding structures, "
                      "templates, or preferred donors for a specific metal.",
            "type": "unknown",
            "suggestions": [
                "What is the coordination geometry of terbium?",
                "Find me structures with calcium and PQQ",
                "What are the preferred donors for zinc?",
                "Get a template for terbium biosensor design",
            ],
        }


def plan_structure_search(task_description: str) -> Dict[str, Any]:
    """
    Plan a structure search based on a task description.

    Parses the task to extract metal, ligand, and application context,
    then returns a structured plan for finding relevant structures.

    Args:
        task_description: Natural language description of design task.

    Returns:
        Dict with:
            - steps: List of recommended steps
            - metal: Detected metal code (or None)
            - ligand: Detected ligand code (or None)
            - task_type: Detected task type (biosensor, enzyme, etc.)
            - search_strategy: Recommended search approach
            - databases: Databases to query

    Example:
        >>> plan = plan_structure_search("I need to design a terbium biosensor")
        >>> print(plan["metal"])
        "TB"
        >>> print(plan["steps"])
        ["Search for terbium binding proteins", ...]
    """
    task_lower = task_description.lower()

    # Extract metal and ligand
    metal = _extract_metal_from_question(task_description)
    ligand = _extract_ligand_from_question(task_description)

    # Detect task type
    task_type = _detect_task_type(task_lower)

    # Build plan
    steps = []
    search_strategy = "general"
    databases = ["metalpdb", "rcsb"]

    if metal:
        # Add metal-specific search step
        steps.append(f"Search for {_get_metal_full_name(metal)} binding proteins")

        # Get coordination info for planning
        try:
            ox_state = METAL_DATABASE.get(metal, {}).get("default_oxidation", 2)
            cn_range = get_coordination_number_range(metal, ox_state)
            donors = get_preferred_donors(metal, ox_state)

            # Add coordination planning step
            donor_list = [d for d, w in sorted(donors.items(), key=lambda x: -x[1]) if w > 0]
            steps.append(
                f"Target coordination number: {cn_range[0]}-{cn_range[1]}, "
                f"preferred donors: {', '.join(donor_list)}"
            )
        except (ValueError, KeyError):
            pass

        # Add metal-specific databases
        if is_lanthanide(metal):
            databases.append("uniprot")  # For lanmodulin-like proteins
            search_strategy = "lanthanide_specialized"

    if ligand:
        steps.append(f"Search for structures containing {ligand} ligand")
        search_strategy = "metal_ligand_complex"

    if task_type == "biosensor":
        steps.append("Look for high-affinity binding sites")
        steps.append("Consider luminescent metals (Tb, Eu) for optical readout")
    elif task_type == "enzyme":
        steps.append("Search for catalytic metal sites")
        steps.append("Review existing enzyme mechanisms")
    elif task_type == "structural":
        steps.append("Search for structural metal sites (zinc fingers, etc.)")

    # Add template recommendation step
    if metal:
        template_name = recommend_template(metal, ligand)
        if template_name:
            steps.append(f"Use '{template_name}' as starting template")

    # Final steps
    steps.append("Extract coordination information from best structures")
    steps.append("Generate design recommendations using get_design_recommendations()")

    return {
        "steps": steps,
        "metal": metal,
        "ligand": ligand,
        "task_type": task_type or "general",
        "search_strategy": search_strategy,
        "databases": databases,
        "raw_task": task_description,
    }


def get_design_recommendations(
    metal: str,
    task: Optional[str] = None,
    ligand: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Get design recommendations for a metal binding site.

    Returns information about preferred donors, available templates,
    amino acid biasing for sequence design, and residues to avoid.

    Args:
        metal: Metal element code (e.g., "TB", "ZN", "CA")
        task: Optional task description for context (e.g., "biosensor", "enzyme")
        ligand: Optional ligand code (e.g., "CIT", "PQQ")

    Returns:
        Dict with:
            - donors: List of preferred donor residues
            - templates: List of available templates
            - bias_aa: Amino acid bias string for LigandMPNN
            - avoid_residues: Residues to avoid
            - coordination_info: Coordination chemistry details
            - hsab_class: HSAB classification
            - warnings: Any warnings or caveats

    Example:
        >>> recs = get_design_recommendations(metal="TB", task="biosensor")
        >>> print(recs["donors"])
        ["Glu", "Asp", "Asn"]
        >>> print(recs["bias_aa"])
        "E:3.0,D:3.0,N:2.4,Q:2.4,H:0.5,C:-5.0,A:-2.0"
    """
    metal = metal.upper()
    warnings = []

    # Validate metal
    if metal not in METAL_DATABASE:
        return {
            "error": f"Unknown metal: {metal}",
            "valid_metals": list(METAL_DATABASE.keys()),
        }

    metal_info = METAL_DATABASE[metal]
    ox_state = metal_info.get("default_oxidation", 2)

    # Get HSAB class
    hsab_class = get_hsab_class(metal, ox_state)

    # Get preferred donors
    try:
        donor_weights = get_preferred_donors(metal, ox_state)
        # Sort by weight descending
        sorted_donors = sorted(donor_weights.items(), key=lambda x: -x[1])
        preferred_elements = [d for d, w in sorted_donors if w > 0]
    except ValueError:
        donor_weights = {}
        preferred_elements = []
        warnings.append(f"Could not get donor preferences for {metal}")

    # Map donor elements to residues
    donor_residues = []
    if "O" in preferred_elements:
        donor_residues.extend(["Glu", "Asp", "Asn", "Gln"])
    if "N" in preferred_elements:
        donor_residues.append("His")
    if "S" in preferred_elements:
        donor_residues.extend(["Cys", "Met"])

    # Determine residues to avoid
    avoid_residues = []
    if hsab_class == "hard":
        avoid_residues.extend(["Cys", "Met"])  # Soft donors bad for hard acids
    elif hsab_class == "soft":
        avoid_residues.extend(["Glu", "Asp", "Asn", "Gln"])  # Hard donors less optimal

    # Get coordination info
    try:
        cn_range = get_coordination_number_range(metal, ox_state)
        coordination_info = {
            "coordination_number_range": cn_range,
            "common_residues": get_common_residues(metal),
            "is_lanthanide": is_lanthanide(metal),
        }
    except ValueError:
        cn_range = (4, 6)  # Default
        coordination_info = {"coordination_number_range": cn_range}
        warnings.append(f"Could not get coordination range for {metal}")

    # Get amino acid bias string for LigandMPNN
    try:
        bias_aa = get_amino_acid_bias(metal, ox_state)
    except ValueError:
        bias_aa = ""
        warnings.append(f"Could not generate bias string for {metal}")

    # Find available templates
    templates = []

    # Check library templates
    if ligand:
        template_name = recommend_template(metal, ligand)
        if template_name:
            templates.append({
                "name": template_name,
                "source": "library",
                "info": get_template_info(template_name),
            })

    # Check for metal-specific templates
    for tname in list_templates():
        tinfo = get_template_info(tname)
        if tinfo and tinfo.get("metal") == metal:
            if tname not in [t["name"] for t in templates]:
                templates.append({
                    "name": tname,
                    "source": "library",
                    "info": tinfo,
                })

    # Get coordination template info
    try:
        coord_template = get_coordination_template(metal)
        if coord_template.get("example_pdbs"):
            for pdb in coord_template["example_pdbs"]:
                templates.append({
                    "name": f"pdb_{pdb}",
                    "source": "pdb",
                    "pdb_id": pdb,
                    "info": {
                        "geometries": coord_template.get("geometries", []),
                        "coordination_number": coord_template.get("coordination_number"),
                    },
                })
    except ValueError:
        pass

    # Task-specific recommendations
    task_recommendations = []
    if task:
        task_lower = task.lower()
        if "biosensor" in task_lower or "sensor" in task_lower:
            if is_lanthanide(metal):
                task_recommendations.append(
                    "Lanthanides (Tb, Eu) are excellent for luminescent biosensors"
                )
            task_recommendations.append(
                "Design for high binding affinity and selectivity"
            )
        elif "enzyme" in task_lower or "catalytic" in task_lower:
            task_recommendations.append(
                "Consider substrate access and reaction mechanism"
            )
        elif "structural" in task_lower:
            task_recommendations.append(
                "Focus on stable, tetrahedral or octahedral geometries"
            )

    return {
        "metal": metal,
        "donors": donor_residues,
        "donor_weights": donor_weights,
        "templates": templates,
        "bias_aa": bias_aa,
        "avoid_residues": avoid_residues,
        "coordination_info": coordination_info,
        "hsab_class": hsab_class,
        "oxidation_state": ox_state,
        "task_recommendations": task_recommendations,
        "warnings": warnings,
    }


# =============================================================================
# INTERNAL QUESTION TYPE DETECTION
# =============================================================================

def _is_coordination_question(question: str) -> bool:
    """Check if question is about coordination chemistry."""
    question_lower = question.lower()
    return any(kw in question_lower for kw in COORDINATION_KEYWORDS)


def _is_find_structure_question(question: str) -> bool:
    """Check if question is about finding structures."""
    question_lower = question.lower()
    return any(kw in question_lower for kw in FIND_STRUCTURE_KEYWORDS)


def _is_template_question(question: str) -> bool:
    """Check if question is about templates."""
    question_lower = question.lower()
    return any(kw in question_lower for kw in TEMPLATE_KEYWORDS)


def _is_donor_question(question: str) -> bool:
    """Check if question is about donor atoms/residues."""
    question_lower = question.lower()
    return any(kw in question_lower for kw in DONOR_KEYWORDS)


# =============================================================================
# INTERNAL QUESTION ANSWERING
# =============================================================================

def _answer_coordination_question(question: str) -> Dict[str, Any]:
    """Answer questions about coordination geometry."""
    metal = _extract_metal_from_question(question)

    if not metal:
        return {
            "answer": "Please specify a metal (e.g., zinc, terbium, calcium) "
                      "to get coordination information.",
            "type": "coordination",
            "error": "no_metal_specified",
        }

    # Get metal info
    if metal not in METAL_DATABASE:
        return {
            "answer": f"Unknown metal: {metal}. Supported metals: {', '.join(METAL_DATABASE.keys())}",
            "type": "coordination",
            "error": "unknown_metal",
        }

    metal_info = METAL_DATABASE[metal]
    ox_state = metal_info.get("default_oxidation", 2)

    # Get coordination data
    cn_range = metal_info["coordination_numbers"].get(ox_state, (4, 6))
    hsab_class = get_hsab_class(metal, ox_state)
    common_residues = metal_info.get("common_residues", [])
    description = metal_info.get("description", "")

    # Get geometry info
    try:
        coord_template = get_coordination_template(metal, cn_range[1])
        geometries = coord_template.get("geometries", [])
    except ValueError:
        geometries = []

    # Build answer
    metal_name = metal_info.get("name", metal)
    geometry_str = ", ".join(geometries) if geometries else "varies"
    residue_str = ", ".join(common_residues)

    answer = (
        f"{metal_name} ({metal}{'{}+'.format(ox_state) if ox_state else ''}) typically has "
        f"coordination number {cn_range[0]}-{cn_range[1]} with {geometry_str} geometry. "
        f"It is classified as a {hsab_class} Lewis acid (HSAB theory). "
        f"Common coordinating residues: {residue_str}. "
        f"{description}"
    )

    return {
        "answer": answer,
        "type": "coordination",
        "metal": metal,
        "metal_name": metal_name,
        "coordination_number": list(cn_range),
        "geometries": geometries,
        "hsab_class": hsab_class,
        "common_residues": common_residues,
        "oxidation_state": ox_state,
    }


def _answer_find_structure_question(question: str) -> Dict[str, Any]:
    """Answer questions about finding structures."""
    metal = _extract_metal_from_question(question)
    ligand = _extract_ligand_from_question(question)

    # Initialize discovery service
    discovery = StructureDiscovery()

    # Search for structures
    if metal and ligand:
        query = f"{_get_metal_full_name(metal)} {ligand}"
        results = discovery.search_by_intent(query, limit=10)
    elif metal:
        results = discovery.search_by_metal(metal, limit=10)
    elif ligand:
        results = discovery.search_by_intent(ligand, limit=10)
    else:
        return {
            "answer": "Please specify a metal and/or ligand to search for structures.",
            "type": "find_structure",
            "error": "no_search_terms",
        }

    # Build answer
    if not results:
        answer = f"No structures found for {'metal=' + metal if metal else ''} {'ligand=' + ligand if ligand else ''}."
        return {
            "answer": answer,
            "type": "find_structure",
            "structures": [],
            "count": 0,
        }

    structures = []
    for r in results[:10]:
        structures.append({
            "pdb_id": r.pdb_id,
            "metal": r.metal,
            "resolution": r.resolution,
            "geometry": r.geometry,
            "coordination_number": r.coordination_number,
            "source": r.source,
            "curated": r.curated,
        })

    # Build readable answer
    pdb_list = ", ".join([s["pdb_id"] for s in structures[:5]])
    answer = (
        f"Found {len(structures)} structures"
        f"{' for ' + _get_metal_full_name(metal) if metal else ''}"
        f"{' with ' + ligand if ligand else ''}. "
        f"Top PDB IDs: {pdb_list}"
    )

    if structures and structures[0].get("curated"):
        answer += " (includes curated reference structures)"

    return {
        "answer": answer,
        "type": "find_structure",
        "structures": structures,
        "count": len(structures),
        "metal": metal,
        "ligand": ligand,
    }


def _answer_template_question(question: str) -> Dict[str, Any]:
    """Answer questions about templates."""
    metal = _extract_metal_from_question(question)
    ligand = _extract_ligand_from_question(question)

    # Find recommended template
    if metal:
        template_name = recommend_template(metal, ligand)
    else:
        template_name = None

    if template_name:
        info = get_template_info(template_name)
        template = get_template(template_name)

        answer = (
            f"Recommended template: '{template_name}' - {info.get('description', '')}. "
            f"Metal: {info.get('metal')}, Ligand: {info.get('ligand')}, "
            f"Coordination: {info.get('coordination_number')}, "
            f"Geometry: {info.get('geometry')}."
        )

        return {
            "answer": answer,
            "type": "template",
            "template_name": template_name,
            "template_info": info,
            "metal": metal,
            "ligand": ligand,
        }
    else:
        # List available templates
        available = list_templates()
        answer = (
            f"No specific template found for {'metal=' + metal if metal else 'unknown metal'}. "
            f"Available templates: {', '.join(available)}"
        )

        return {
            "answer": answer,
            "type": "template",
            "available_templates": available,
            "metal": metal,
            "ligand": ligand,
        }


def _answer_donor_question(question: str) -> Dict[str, Any]:
    """Answer questions about preferred donor residues."""
    metal = _extract_metal_from_question(question)

    if not metal:
        return {
            "answer": "Please specify a metal to get donor preferences.",
            "type": "donors",
            "error": "no_metal_specified",
        }

    if metal not in METAL_DATABASE:
        return {
            "answer": f"Unknown metal: {metal}",
            "type": "donors",
            "error": "unknown_metal",
        }

    metal_info = METAL_DATABASE[metal]
    ox_state = metal_info.get("default_oxidation", 2)

    # Get donors
    try:
        donors = get_preferred_donors(metal, ox_state)
        hsab_class = get_hsab_class(metal, ox_state)
        common_residues = get_common_residues(metal)
        bias_aa = get_amino_acid_bias(metal, ox_state)
    except ValueError as e:
        return {
            "answer": f"Error getting donor info: {e}",
            "type": "donors",
            "error": str(e),
        }

    # Sort donors by weight
    sorted_donors = sorted(donors.items(), key=lambda x: -x[1])
    donor_str = ", ".join([f"{d}({w:+.1f})" for d, w in sorted_donors])
    residue_str = ", ".join(common_residues)

    metal_name = metal_info.get("name", metal)

    answer = (
        f"{metal_name} ({hsab_class} acid) prefers donors in order: {donor_str}. "
        f"Best coordinating residues: {residue_str}. "
        f"LigandMPNN bias string: {bias_aa}"
    )

    return {
        "answer": answer,
        "type": "donors",
        "metal": metal,
        "donor_weights": donors,
        "common_residues": common_residues,
        "hsab_class": hsab_class,
        "bias_aa": bias_aa,
    }


def _answer_general_metal_question(question: str, metal: str) -> Dict[str, Any]:
    """Answer general questions about a metal."""
    if metal not in METAL_DATABASE:
        return {
            "answer": f"Unknown metal: {metal}",
            "type": "general",
            "error": "unknown_metal",
        }

    metal_info = METAL_DATABASE[metal]
    ox_state = metal_info.get("default_oxidation", 2)

    # Gather all info
    cn_range = metal_info["coordination_numbers"].get(ox_state, (4, 6))
    hsab_class = get_hsab_class(metal, ox_state)
    common_residues = metal_info.get("common_residues", [])
    description = metal_info.get("description", "")

    metal_name = metal_info.get("name", metal)
    residue_str = ", ".join(common_residues)

    answer = (
        f"{metal_name} ({metal}): {description} "
        f"CN: {cn_range[0]}-{cn_range[1]}, HSAB: {hsab_class}, "
        f"Residues: {residue_str}."
    )

    return {
        "answer": answer,
        "type": "general",
        "metal": metal,
        "metal_name": metal_name,
        "coordination_number_range": cn_range,
        "hsab_class": hsab_class,
        "common_residues": common_residues,
        "description": description,
    }


# =============================================================================
# INTERNAL EXTRACTION HELPERS
# =============================================================================

def _extract_metal_from_question(question: str) -> Optional[str]:
    """Extract metal code from natural language question."""
    question_lower = question.lower()

    # Check for metal names first (longer match preferred)
    sorted_metals = sorted(METAL_NAMES.items(), key=lambda x: len(x[0]), reverse=True)
    for name, code in sorted_metals:
        if name in question_lower:
            return code

    # Check for metal codes (2 letters, uppercase)
    for code in METAL_DATABASE.keys():
        # Match as whole word
        pattern = rf'\b{re.escape(code)}\b'
        if re.search(pattern, question, re.IGNORECASE):
            return code

    return None


def _extract_ligand_from_question(question: str) -> Optional[str]:
    """Extract ligand code from natural language question."""
    question_lower = question.lower()

    # Check for ligand names first
    sorted_ligands = sorted(LIGAND_NAMES.items(), key=lambda x: len(x[0]), reverse=True)
    for name, code in sorted_ligands:
        if name in question_lower:
            return code

    # Check for 3-letter codes in uppercase
    ligand_pattern = r'\b([A-Z]{3})\b'
    matches = re.findall(ligand_pattern, question)
    for match in matches:
        # Exclude metal codes
        if match not in METAL_DATABASE:
            return match

    return None


def _get_metal_full_name(metal_code: str) -> str:
    """Get full name for a metal code."""
    code_to_name = {v: k for k, v in METAL_NAMES.items()}
    return code_to_name.get(metal_code.upper(), metal_code)


def _detect_task_type(task_lower: str) -> Optional[str]:
    """Detect the type of design task from description."""
    if any(kw in task_lower for kw in ["biosensor", "sensor", "detect", "luminescent"]):
        return "biosensor"
    elif any(kw in task_lower for kw in ["enzyme", "catalytic", "catalyst", "reaction"]):
        return "enzyme"
    elif any(kw in task_lower for kw in ["structural", "stabilize", "fold", "zinc finger"]):
        return "structural"
    elif any(kw in task_lower for kw in ["binding", "affinity", "selectivity"]):
        return "binding"
    return None
