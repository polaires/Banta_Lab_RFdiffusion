"""
Natural Language Design Parser

Parses natural language protein design requests into structured DesignIntent
using Claude API.

Example inputs:
- "Design a protein to bind citrate with terbium"
- "Create a zinc finger protein"
- "I want a PQQ-binding dehydrogenase with calcium"

Uses Claude API for semantic understanding and metal_chemistry.py for validation.
"""

import logging
import json
import re
import requests
from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional, Tuple

# Import metal chemistry for validation
from metal_chemistry import (
    METAL_DATABASE,
    get_hsab_class,
    get_amino_acid_bias,
    is_lanthanide,
)

logger = logging.getLogger(__name__)


# =============================================================================
# Data Classes
# =============================================================================

# =============================================================================
# PDB ID and Scaffolding Detection
# =============================================================================

# PDB ID pattern: 4 characters starting with a digit (e.g., 4CVB, 1ABC, 7XYZ)
PDB_ID_PATTERN = re.compile(r'\b([0-9][A-Za-z0-9]{3})\b')

# Keywords that indicate scaffolding intent (rather than de novo design)
SCAFFOLD_KEYWORDS = [
    "scaffold", "scaffolding", "motif", "pocket of", "site of",
    "from structure", "from pdb", "active site of", "binding site of",
    "extract from", "based on", "template from"
]

# Keywords that indicate user wants ALL interacting residues (not just metal coordination)
FULL_POCKET_KEYWORDS = [
    "all interacting", "all residues", "keep all", "preserve all",
    "full pocket", "complete pocket", "entire pocket", "whole pocket",
    "all contacts", "all interactions", "binding pocket"
]

# Keywords that indicate stability optimization
STABILITY_KEYWORDS = [
    "more stable", "stable", "stability", "stabilize", "stabilization",
    "robust", "thermostable", "thermal stability", "improve stability"
]


@dataclass
class DesignIntent:
    """
    Structured representation of a protein design request.

    Extracted from natural language input by Claude API.
    """
    # Core design parameters
    metal_type: Optional[str] = None       # Metal symbol (TB, ZN, CA, etc.)
    ligand_name: Optional[str] = None      # Ligand name (citrate, PQQ, etc.)
    design_goal: str = "binding"           # binding | catalysis | sensing | structural
    target_topology: str = "monomer"       # monomer | dimer | symmetric | custom

    # Size parameters
    chain_length_min: int = 80             # Minimum residue count
    chain_length_max: int = 120            # Maximum residue count

    # Scaffolding parameters
    design_mode: str = "de_novo"           # "de_novo" | "scaffold"
    source_pdb_id: Optional[str] = None    # e.g., "4CVB" for scaffold mode
    pocket_description: Optional[str] = None  # e.g., "PQQ-Ca pocket"

    # Extended pocket and stability options
    include_all_contacts: bool = False     # If True, include all ligand-contacting residues
    stability_focus: bool = False          # If True, optimize for protein stability

    # Confidence and metadata
    confidence: float = 0.0                # Parser confidence (0-1)
    raw_query: str = ""                    # Original user query

    # Warnings and suggestions
    warnings: List[str] = field(default_factory=list)
    suggestions: List[str] = field(default_factory=list)

    # Parsed entities for debugging
    parsed_entities: Dict[str, Any] = field(default_factory=dict)

    @property
    def is_scaffolding(self) -> bool:
        """Check if this is a scaffolding request."""
        return self.design_mode == "scaffold" and self.source_pdb_id is not None

    @property
    def chain_length_range(self) -> str:
        """Return chain length as RFD3 contig range string."""
        return f"{self.chain_length_min}-{self.chain_length_max}"

    @property
    def has_metal(self) -> bool:
        """Check if design includes a metal."""
        return self.metal_type is not None

    @property
    def has_ligand(self) -> bool:
        """Check if design includes a ligand."""
        return self.ligand_name is not None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return asdict(self)


# =============================================================================
# Metal/Ligand Name Normalization
# =============================================================================

# Common metal name variations
METAL_ALIASES: Dict[str, str] = {
    # Lanthanides
    "terbium": "TB",
    "tb": "TB",
    "tb3+": "TB",
    "europium": "EU",
    "eu": "EU",
    "eu3+": "EU",
    "gadolinium": "GD",
    "gd": "GD",
    "gd3+": "GD",
    "lanthanum": "LA",
    "la": "LA",
    "la3+": "LA",
    "cerium": "CE",
    "ce": "CE",
    "ce3+": "CE",
    "samarium": "SM",
    "sm": "SM",
    "ytterbium": "YB",
    "yb": "YB",

    # Transition metals
    "zinc": "ZN",
    "zn": "ZN",
    "zn2+": "ZN",
    "iron": "FE",
    "fe": "FE",
    "fe2+": "FE",
    "fe3+": "FE",
    "copper": "CU",
    "cu": "CU",
    "cu2+": "CU",
    "manganese": "MN",
    "mn": "MN",
    "cobalt": "CO",
    "co": "CO",
    "nickel": "NI",
    "ni": "NI",

    # Alkaline earth
    "calcium": "CA",
    "ca": "CA",
    "ca2+": "CA",
    "magnesium": "MG",
    "mg": "MG",
    "mg2+": "MG",

    # Shorthand
    "ln": "TB",  # Default lanthanide
    "lanthanide": "TB",
}

# Common ligand name variations
LIGAND_ALIASES: Dict[str, str] = {
    # Citrate variations
    "citrate": "citrate",
    "citric acid": "citrate",
    "citric": "citrate",
    "cit": "citrate",

    # PQQ variations
    "pqq": "pqq",
    "pyrroloquinoline quinone": "pqq",
    "pyrroloquinoline": "pqq",

    # Common cofactors
    "heme": "heme",
    "hem": "heme",
    "protoporphyrin": "heme",

    # Other ligands
    "atp": "atp",
    "adp": "adp",
    "nad": "nad",
    "nadp": "nadp",
    "fad": "fad",
    "fmn": "fmn",
}


def normalize_metal(name: str) -> Optional[str]:
    """Normalize metal name to standard symbol."""
    if not name:
        return None

    name_lower = name.lower().strip()

    # Check aliases
    if name_lower in METAL_ALIASES:
        return METAL_ALIASES[name_lower]

    # Check if already a valid symbol
    name_upper = name.upper().strip()
    if name_upper in METAL_DATABASE:
        return name_upper

    return None


def normalize_ligand(name: str) -> Optional[str]:
    """Normalize ligand name to standard form."""
    if not name:
        return None

    name_lower = name.lower().strip()

    # Check aliases
    if name_lower in LIGAND_ALIASES:
        return LIGAND_ALIASES[name_lower]

    # Return as-is if not found (might be valid)
    return name_lower


# =============================================================================
# Claude API Parser
# =============================================================================

# System prompt for Claude to parse design requests
PARSER_SYSTEM_PROMPT = """You are a protein design assistant that extracts structured information from natural language design requests.

Your task is to parse user queries about protein design and extract key parameters.

For each query, extract:
1. metal_type: The metal ion involved (use standard symbols: ZN, FE, CA, TB, EU, etc.)
2. ligand_name: The ligand/cofactor name (citrate, PQQ, heme, etc.)
3. design_goal: One of: binding, catalysis, sensing, structural
4. target_topology: One of: monomer, dimer, symmetric, custom
5. chain_length: Suggested protein size (small: 60-80, medium: 80-120, large: 120-160)
6. confidence: Your confidence in the parsing (0.0 to 1.0)

Common metal-ligand combinations:
- Citrate with lanthanides (Tb, Eu, Gd) - luminescent biosensors
- PQQ with calcium - quinoprotein dehydrogenases
- Heme with iron - electron transfer, catalysis
- Zinc fingers - DNA binding

Return ONLY valid JSON with these fields. No explanation text."""


PARSER_USER_TEMPLATE = """Parse this protein design request:

"{query}"

Return JSON with fields: metal_type, ligand_name, design_goal, target_topology, chain_length_min, chain_length_max, confidence, reasoning
"""


class NLDesignParser:
    """
    Natural Language Design Parser using Claude API.

    Parses natural language protein design requests into structured DesignIntent.
    """

    def __init__(
        self,
        api_key: str,
        api_url: str = "https://yinli.one/v1/messages",
        model: str = "claude-sonnet-4-5-20250929",
        timeout: int = 30,
    ):
        """
        Initialize the parser.

        Args:
            api_key: Claude API key
            api_url: API endpoint URL
            model: Model to use for parsing
            timeout: Request timeout in seconds
        """
        self.api_key = api_key
        self.api_url = api_url
        self.model = model
        self.timeout = timeout

    def parse(self, query: str) -> DesignIntent:
        """
        Parse a natural language design query into structured intent.

        Args:
            query: Natural language design request

        Returns:
            DesignIntent with extracted parameters
        """
        intent = DesignIntent(raw_query=query)

        try:
            # Call Claude API
            parsed = self._call_claude_api(query)

            if parsed:
                # Extract and normalize metal
                if parsed.get("metal_type"):
                    metal = normalize_metal(parsed["metal_type"])
                    if metal:
                        intent.metal_type = metal
                        intent.parsed_entities["metal"] = parsed["metal_type"]
                    else:
                        intent.warnings.append(
                            f"Unknown metal: {parsed['metal_type']}"
                        )

                # Extract and normalize ligand
                if parsed.get("ligand_name"):
                    ligand = normalize_ligand(parsed["ligand_name"])
                    if ligand:
                        intent.ligand_name = ligand
                        intent.parsed_entities["ligand"] = parsed["ligand_name"]

                # Extract design goal
                if parsed.get("design_goal"):
                    goal = parsed["design_goal"].lower()
                    if goal in ["binding", "catalysis", "sensing", "structural"]:
                        intent.design_goal = goal

                # Extract topology
                if parsed.get("target_topology"):
                    topo = parsed["target_topology"].lower()
                    if topo in ["monomer", "dimer", "symmetric", "custom"]:
                        intent.target_topology = topo

                # Extract chain length
                if parsed.get("chain_length_min"):
                    try:
                        intent.chain_length_min = int(parsed["chain_length_min"])
                    except (ValueError, TypeError):
                        pass

                if parsed.get("chain_length_max"):
                    try:
                        intent.chain_length_max = int(parsed["chain_length_max"])
                    except (ValueError, TypeError):
                        pass

                # Ensure valid range
                if intent.chain_length_min > intent.chain_length_max:
                    intent.chain_length_min, intent.chain_length_max = (
                        intent.chain_length_max, intent.chain_length_min
                    )

                # Extract confidence
                if parsed.get("confidence"):
                    try:
                        intent.confidence = float(parsed["confidence"])
                    except (ValueError, TypeError):
                        intent.confidence = 0.5

                # Store reasoning
                if parsed.get("reasoning"):
                    intent.parsed_entities["reasoning"] = parsed["reasoning"]

        except Exception as e:
            logger.error(f"Error parsing query: {e}")
            intent.warnings.append(f"Parsing error: {str(e)}")
            intent.confidence = 0.0

        # Validate and add warnings
        self._validate_intent(intent)

        return intent

    def _call_claude_api(self, query: str) -> Optional[Dict[str, Any]]:
        """
        Call Claude API to parse the query.

        Returns parsed JSON or None on error.
        """
        headers = {
            "Content-Type": "application/json",
            "x-api-key": self.api_key,
            "anthropic-version": "2023-06-01",
        }

        payload = {
            "model": self.model,
            "max_tokens": 1024,
            "system": PARSER_SYSTEM_PROMPT,
            "messages": [
                {
                    "role": "user",
                    "content": PARSER_USER_TEMPLATE.format(query=query)
                }
            ]
        }

        try:
            response = requests.post(
                self.api_url,
                headers=headers,
                json=payload,
                timeout=self.timeout,
            )
            response.raise_for_status()

            result = response.json()

            # Extract text from response
            content = result.get("content", [])
            if content and len(content) > 0:
                text = content[0].get("text", "")

                # Parse JSON from response
                # Handle potential markdown code blocks
                if "```json" in text:
                    text = text.split("```json")[1].split("```")[0]
                elif "```" in text:
                    text = text.split("```")[1].split("```")[0]

                return json.loads(text.strip())

        except requests.RequestException as e:
            logger.error(f"API request failed: {e}")
        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse JSON response: {e}")
        except Exception as e:
            logger.error(f"Unexpected error in API call: {e}")

        return None

    def _validate_intent(self, intent: DesignIntent) -> None:
        """
        Validate the parsed intent and add warnings/suggestions.

        Checks HSAB compatibility, coordination requirements, etc.
        """
        # Check if we have enough information
        if not intent.has_metal and not intent.has_ligand:
            intent.warnings.append(
                "No metal or ligand detected. Design may need more specifics."
            )
            intent.suggestions.append(
                "Try specifying a metal (e.g., 'with terbium') or ligand (e.g., 'bind citrate')"
            )
            return

        # Validate metal if present
        if intent.has_metal:
            metal = intent.metal_type

            # Check if metal is in database
            if metal not in METAL_DATABASE:
                intent.warnings.append(f"Metal {metal} not in chemistry database")
            else:
                # Get HSAB class for suggestions
                try:
                    metal_info = METAL_DATABASE[metal]
                    default_ox = metal_info.get("default_oxidation", 2)
                    hsab_class = get_hsab_class(metal, default_ox)

                    intent.parsed_entities["hsab_class"] = hsab_class
                    intent.parsed_entities["default_oxidation"] = default_ox

                    # Add HSAB-based suggestions
                    if hsab_class == "hard":
                        intent.suggestions.append(
                            f"{metal} is a hard acid - prefer Glu/Asp/Asn donors"
                        )
                        if is_lanthanide(metal):
                            intent.suggestions.append(
                                f"Lanthanide {metal} prefers O donors, exclude Cys"
                            )
                    elif hsab_class == "soft":
                        intent.suggestions.append(
                            f"{metal} is a soft acid - consider Cys/Met donors"
                        )
                    else:
                        intent.suggestions.append(
                            f"{metal} is borderline - His/Cys/Glu all viable"
                        )

                    # Get amino acid bias for later use
                    bias = get_amino_acid_bias(metal, default_ox)
                    intent.parsed_entities["suggested_bias_aa"] = bias

                except ValueError as e:
                    intent.warnings.append(f"Metal validation error: {e}")

        # Suggest templates based on metal-ligand combination
        if intent.has_metal and intent.has_ligand:
            metal = intent.metal_type
            ligand = intent.ligand_name

            # Known good combinations
            if ligand == "citrate" and is_lanthanide(metal):
                intent.suggestions.append(
                    f"Citrate-{metal} is a known template - using optimized config"
                )
                intent.parsed_entities["template_match"] = f"citrate_{metal.lower()}"

            elif ligand == "pqq" and metal == "CA":
                intent.suggestions.append(
                    "PQQ-Ca is a known template - using optimized config"
                )
                intent.parsed_entities["template_match"] = "pqq_ca"

        # Adjust chain length based on topology
        if intent.target_topology == "dimer":
            if intent.chain_length_max < 100:
                intent.suggestions.append(
                    "Dimers typically need longer chains (100-150) for interface"
                )

        # Adjust for lanthanide coordination requirements
        if intent.has_metal and is_lanthanide(intent.metal_type):
            if intent.chain_length_max < 100:
                intent.suggestions.append(
                    f"Lanthanide {intent.metal_type} needs 8-9 coordination sites - consider larger protein"
                )


# =============================================================================
# Fallback Parser (No API)
# =============================================================================

class SimpleFallbackParser:
    """
    Simple rule-based fallback parser when API is unavailable.

    Uses keyword matching for basic extraction.
    """

    def parse(self, query: str) -> DesignIntent:
        """Parse using simple keyword matching."""
        intent = DesignIntent(raw_query=query)
        query_lower = query.lower()

        # Detect PDB ID (4 characters starting with digit, e.g., 4CVB)
        pdb_match = PDB_ID_PATTERN.search(query.upper())
        if pdb_match:
            intent.source_pdb_id = pdb_match.group(1)
            intent.parsed_entities["pdb_id"] = intent.source_pdb_id

        # Detect scaffolding intent
        is_scaffold = any(kw in query_lower for kw in SCAFFOLD_KEYWORDS)
        if is_scaffold and intent.source_pdb_id:
            intent.design_mode = "scaffold"
            intent.parsed_entities["design_mode"] = "scaffold"
            # Try to extract pocket description (text between ligand/metal and PDB ID)
            intent.pocket_description = self._extract_pocket_description(query)
        else:
            intent.design_mode = "de_novo"

        # Extract metal
        for alias, metal in METAL_ALIASES.items():
            if alias in query_lower:
                intent.metal_type = metal
                intent.parsed_entities["metal_match"] = alias
                break

        # Extract ligand
        for alias, ligand in LIGAND_ALIASES.items():
            if alias in query_lower:
                intent.ligand_name = ligand
                intent.parsed_entities["ligand_match"] = alias
                break

        # Detect design goal
        if any(word in query_lower for word in ["sensor", "sensing", "detect"]):
            intent.design_goal = "sensing"
        elif any(word in query_lower for word in ["catalysis", "enzyme", "react"]):
            intent.design_goal = "catalysis"
        elif any(word in query_lower for word in ["structural", "scaffold", "fold"]):
            intent.design_goal = "structural"
        else:
            intent.design_goal = "binding"

        # Detect topology
        if any(word in query_lower for word in ["dimer", "two chain", "interface"]):
            intent.target_topology = "dimer"
        elif any(word in query_lower for word in ["symmetric", "symmetry", "c3", "c4"]):
            intent.target_topology = "symmetric"

        # Detect full pocket request (keep all interacting residues)
        if any(kw in query_lower for kw in FULL_POCKET_KEYWORDS):
            intent.include_all_contacts = True
            intent.parsed_entities["include_all_contacts"] = True

        # Detect stability optimization request
        if any(kw in query_lower for kw in STABILITY_KEYWORDS):
            intent.stability_focus = True
            intent.parsed_entities["stability_focus"] = True

        # Set confidence based on matches
        matches = sum([
            intent.has_metal,
            intent.has_ligand,
            intent.design_goal != "binding",
            intent.target_topology != "monomer",
            intent.is_scaffolding,  # Scaffolding detection
            intent.include_all_contacts,  # Full pocket detection
            intent.stability_focus,  # Stability detection
        ])
        intent.confidence = 0.3 + (matches * 0.15)

        intent.warnings.append("Used fallback parser (API unavailable)")

        return intent

    def _extract_pocket_description(self, query: str) -> Optional[str]:
        """
        Extract pocket description from query.

        Examples:
        - "scaffold the PQQ-Ca pocket of 4CVB" -> "PQQ-Ca pocket"
        - "scaffold the active site of 1ABC" -> "active site"
        """
        query_lower = query.lower()

        # Try common patterns
        patterns = [
            r'the\s+(\S+-\S+\s+pocket)',     # "the PQQ-Ca pocket"
            r'the\s+(\S+\s+pocket)',          # "the zinc pocket"
            r'(active\s+site)',               # "active site"
            r'(binding\s+site)',              # "binding site"
            r'(coordination\s+site)',         # "coordination site"
        ]

        for pattern in patterns:
            match = re.search(pattern, query_lower)
            if match:
                return match.group(1)

        return None


# =============================================================================
# Factory Function
# =============================================================================

def create_parser(
    api_key: Optional[str] = None,
    api_url: str = "https://yinli.one/v1/messages",
    use_fallback: bool = True,
) -> NLDesignParser:
    """
    Create an appropriate parser based on configuration.

    Args:
        api_key: Claude API key (None for fallback)
        api_url: API endpoint URL
        use_fallback: Use fallback parser if API unavailable

    Returns:
        Parser instance
    """
    if api_key:
        return NLDesignParser(api_key=api_key, api_url=api_url)
    elif use_fallback:
        return SimpleFallbackParser()
    else:
        raise ValueError("API key required when fallback disabled")


# =============================================================================
# Test Function
# =============================================================================

def test_parser():
    """Test the parser with example queries."""
    test_queries = [
        # De novo design queries
        "Design a protein to bind citrate with terbium",
        "Create a zinc finger DNA-binding protein",
        "I want a PQQ-binding dehydrogenase with calcium",
        "Make a luminescent sensor for europium",
        "Design a dimer interface with lanthanide binding",
        # Scaffolding queries (NEW)
        "Scaffold the PQQ-Ca pocket of 4CVB",
        "Scaffold the active site of 1ABC with zinc",
        "Extract the citrate binding site from 3XYZ",
    ]

    # Use fallback parser for testing without API
    parser = SimpleFallbackParser()

    for query in test_queries:
        print(f"\nQuery: {query}")
        intent = parser.parse(query)
        print(f"  Design Mode: {intent.design_mode}")
        if intent.source_pdb_id:
            print(f"  Source PDB: {intent.source_pdb_id}")
        if intent.pocket_description:
            print(f"  Pocket: {intent.pocket_description}")
        print(f"  Metal: {intent.metal_type}")
        print(f"  Ligand: {intent.ligand_name}")
        print(f"  Goal: {intent.design_goal}")
        print(f"  Topology: {intent.target_topology}")
        print(f"  Confidence: {intent.confidence:.2f}")
        if intent.warnings:
            print(f"  Warnings: {intent.warnings}")


if __name__ == "__main__":
    test_parser()
