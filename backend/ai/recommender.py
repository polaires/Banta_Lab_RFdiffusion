"""
Enhanced Hybrid AI Recommendation Engine

This module implements a two-layer recommendation system:
1. Chemistry Rule Base: Hard-coded rules from research/literature
2. LLM Refinement: Claude API for natural language interpretation and parameter tuning

The rule base ensures scientifically valid recommendations while
the LLM layer provides flexibility for edge cases and user interaction.
"""

import json
from typing import Dict, List, Any, Optional, Tuple
from dataclasses import dataclass, asdict
from enum import Enum

# Try to import anthropic for LLM integration
ANTHROPIC_AVAILABLE = False
try:
    import anthropic
    ANTHROPIC_AVAILABLE = True
except ImportError:
    pass


# ============== Chemistry Rule Base ==============

class MetalClass(Enum):
    """Classification of metal ions by binding behavior."""
    HARD_ACID = "hard_acid"         # Prefers O, N donors (lanthanides, Ca, Mg)
    BORDERLINE = "borderline"        # Flexible (Fe, Zn, Cu, Ni)
    SOFT_ACID = "soft_acid"          # Prefers S donors (Ag, Au, Hg)


@dataclass
class MetalProperties:
    """Properties of a metal ion for binding prediction."""
    name: str
    symbol: str
    metal_class: MetalClass
    typical_coordination: List[int]
    preferred_geometries: List[str]
    preferred_donors: List[str]
    avoid_donors: List[str]
    ionic_radius: float  # Angstroms for 3+ or most common oxidation state
    bond_distance_range: Tuple[float, float]
    special_notes: Optional[str] = None


class MetalChemistryRules:
    """
    Chemistry rule base for metal binding preferences.

    Based on hard-soft acid-base (HSAB) theory and crystallographic data.
    """

    # Metal properties database
    METALS = {
        # Lanthanides (hard acids, high coordination, O preference)
        "LA": MetalProperties(
            name="Lanthanum", symbol="LA", metal_class=MetalClass.HARD_ACID,
            typical_coordination=[9, 10, 12], preferred_geometries=["tricapped_trigonal_prism"],
            preferred_donors=["O_carboxylate", "O_carbonyl", "O_water", "O_hydroxyl"],
            avoid_donors=["S_thiolate", "S_thioether"],
            ionic_radius=1.16, bond_distance_range=(2.4, 2.8),
            special_notes="Largest lanthanide, highest coordination numbers"
        ),
        "TB": MetalProperties(
            name="Terbium", symbol="TB", metal_class=MetalClass.HARD_ACID,
            typical_coordination=[8, 9], preferred_geometries=["square_antiprism", "tricapped_trigonal_prism"],
            preferred_donors=["O_carboxylate", "O_carbonyl", "O_water", "O_amide"],
            avoid_donors=["S_thiolate", "S_thioether"],
            ionic_radius=1.04, bond_distance_range=(2.3, 2.6),
            special_notes="Mid-series lanthanide, luminescent properties"
        ),
        "GD": MetalProperties(
            name="Gadolinium", symbol="GD", metal_class=MetalClass.HARD_ACID,
            typical_coordination=[8, 9], preferred_geometries=["square_antiprism", "tricapped_trigonal_prism"],
            preferred_donors=["O_carboxylate", "O_carbonyl", "O_water"],
            avoid_donors=["S_thiolate", "S_thioether"],
            ionic_radius=1.05, bond_distance_range=(2.3, 2.6),
            special_notes="MRI contrast agent, paramagnetic"
        ),
        "EU": MetalProperties(
            name="Europium", symbol="EU", metal_class=MetalClass.HARD_ACID,
            typical_coordination=[8, 9], preferred_geometries=["square_antiprism", "bicapped_trigonal_prism"],
            preferred_donors=["O_carboxylate", "O_carbonyl", "O_water"],
            avoid_donors=["S_thiolate", "S_thioether"],
            ionic_radius=1.07, bond_distance_range=(2.3, 2.6),
            special_notes="Red luminescence, biosensing applications"
        ),

        # Alkaline earth metals
        "CA": MetalProperties(
            name="Calcium", symbol="CA", metal_class=MetalClass.HARD_ACID,
            typical_coordination=[6, 7, 8], preferred_geometries=["pentagonal_bipyramidal", "square_antiprism"],
            preferred_donors=["O_carboxylate", "O_carbonyl", "O_water", "O_hydroxyl"],
            avoid_donors=["S_thiolate"],
            ionic_radius=1.00, bond_distance_range=(2.3, 2.6),
            special_notes="EF-hand proteins, signaling"
        ),
        "MG": MetalProperties(
            name="Magnesium", symbol="MG", metal_class=MetalClass.HARD_ACID,
            typical_coordination=[6], preferred_geometries=["octahedral"],
            preferred_donors=["O_carboxylate", "O_water", "O_phosphate"],
            avoid_donors=["S_thiolate"],
            ionic_radius=0.72, bond_distance_range=(2.0, 2.3),
            special_notes="Catalysis, ATP binding"
        ),

        # Transition metals (borderline)
        "FE": MetalProperties(
            name="Iron", symbol="FE", metal_class=MetalClass.BORDERLINE,
            typical_coordination=[4, 6], preferred_geometries=["tetrahedral", "octahedral"],
            preferred_donors=["S_thiolate", "N_imidazole", "O_carboxylate", "N_porphyrin"],
            avoid_donors=[],
            ionic_radius=0.64, bond_distance_range=(1.9, 2.5),
            special_notes="Iron-sulfur clusters, heme proteins"
        ),
        "ZN": MetalProperties(
            name="Zinc", symbol="ZN", metal_class=MetalClass.BORDERLINE,
            typical_coordination=[4, 5, 6], preferred_geometries=["tetrahedral", "trigonal_bipyramidal"],
            preferred_donors=["S_thiolate", "N_imidazole", "O_carboxylate"],
            avoid_donors=[],
            ionic_radius=0.74, bond_distance_range=(1.9, 2.4),
            special_notes="Zinc fingers, catalytic sites"
        ),
        "CU": MetalProperties(
            name="Copper", symbol="CU", metal_class=MetalClass.BORDERLINE,
            typical_coordination=[4, 5, 6], preferred_geometries=["square_planar", "square_pyramidal"],
            preferred_donors=["N_imidazole", "S_thiolate", "S_thioether"],
            avoid_donors=[],
            ionic_radius=0.73, bond_distance_range=(1.9, 2.3),
            special_notes="Electron transfer, Type 1/2/3 copper"
        ),
        "MN": MetalProperties(
            name="Manganese", symbol="MN", metal_class=MetalClass.BORDERLINE,
            typical_coordination=[4, 6], preferred_geometries=["tetrahedral", "octahedral"],
            preferred_donors=["O_carboxylate", "N_imidazole", "O_water"],
            avoid_donors=[],
            ionic_radius=0.83, bond_distance_range=(2.0, 2.4),
            special_notes="OEC in photosynthesis, SOD"
        ),
        "CO": MetalProperties(
            name="Cobalt", symbol="CO", metal_class=MetalClass.BORDERLINE,
            typical_coordination=[4, 6], preferred_geometries=["tetrahedral", "octahedral"],
            preferred_donors=["N_imidazole", "O_carboxylate", "S_thiolate"],
            avoid_donors=[],
            ionic_radius=0.61, bond_distance_range=(1.9, 2.3),
            special_notes="B12 cofactor, catalysis"
        ),
        "NI": MetalProperties(
            name="Nickel", symbol="NI", metal_class=MetalClass.BORDERLINE,
            typical_coordination=[4, 6], preferred_geometries=["square_planar", "octahedral"],
            preferred_donors=["S_thiolate", "N_imidazole", "O_carboxylate"],
            avoid_donors=[],
            ionic_radius=0.69, bond_distance_range=(1.9, 2.2),
            special_notes="Hydrogenases, urease"
        ),
    }

    @classmethod
    def get_metal_properties(cls, metal: str) -> Optional[MetalProperties]:
        """Get properties for a metal ion."""
        return cls.METALS.get(metal.upper())

    @classmethod
    def is_compatible_donor(cls, metal: str, donor_type: str) -> Tuple[bool, str]:
        """
        Check if a donor type is compatible with a metal.

        Returns:
            Tuple of (is_compatible, explanation)
        """
        props = cls.get_metal_properties(metal)
        if props is None:
            return True, "Unknown metal - assuming compatible"

        if donor_type in props.avoid_donors:
            return False, f"{donor_type} is not recommended for {props.name}"

        if donor_type in props.preferred_donors:
            return True, f"{donor_type} is preferred for {props.name}"

        return True, f"{donor_type} is acceptable for {props.name}"

    @classmethod
    def calculate_coordination_delta(
        cls,
        current_metal: str,
        target_metal: str,
        current_coordination: int
    ) -> Dict[str, Any]:
        """
        Calculate the coordination change needed for metal swap.

        Returns:
            Dict with coordination analysis and recommendations
        """
        current_props = cls.get_metal_properties(current_metal)
        target_props = cls.get_metal_properties(target_metal)

        if target_props is None:
            return {
                "error": f"Unknown target metal: {target_metal}",
                "delta": 0,
                "target_coordination": current_coordination,
            }

        target_min = min(target_props.typical_coordination)
        target_max = max(target_props.typical_coordination)

        delta = target_min - current_coordination

        return {
            "current_metal": current_metal,
            "target_metal": target_metal,
            "current_coordination": current_coordination,
            "target_coordination_range": target_props.typical_coordination,
            "delta": delta,
            "needs_increase": delta > 0,
            "needs_decrease": delta < 0,
            "recommendation": (
                f"Increase coordination by {delta}" if delta > 0 else
                f"Decrease coordination by {-delta}" if delta < 0 else
                "Coordination is appropriate"
            ),
        }


class DesignStrategyRules:
    """
    Rules for choosing RFD3 design strategies based on the design goal.
    """

    STRATEGIES = {
        "increase_coordination": {
            "description": "Add more coordinating residues around metal site",
            "partial_t_range": (12.0, 18.0),
            "recommended_partial_t": 15.0,
            "fix_backbone": True,
            "allow_length_change": True,
            "length_expansion": 1.2,  # Allow 20% more residues
            "num_timesteps": 200,
            "step_scale": 1.5,
            "gamma_0": 0.6,
        },
        "decrease_coordination": {
            "description": "Remove coordinating residues or increase site volume",
            "partial_t_range": (8.0, 12.0),
            "recommended_partial_t": 10.0,
            "fix_backbone": True,
            "allow_length_change": False,
            "num_timesteps": 200,
            "step_scale": 1.5,
            "gamma_0": 0.6,
        },
        "change_donor_type": {
            "description": "Replace donor residues (e.g., Cys to Asp/Glu)",
            "partial_t_range": (10.0, 15.0),
            "recommended_partial_t": 12.0,
            "fix_backbone": True,
            "allow_length_change": False,
            "num_timesteps": 200,
            "step_scale": 1.5,
            "gamma_0": 0.6,
        },
        "expand_pocket": {
            "description": "Make binding pocket larger",
            "partial_t_range": (15.0, 20.0),
            "recommended_partial_t": 18.0,
            "fix_backbone": False,
            "allow_length_change": True,
            "num_timesteps": 250,
            "step_scale": 1.3,
            "gamma_0": 0.7,
        },
        "refine_geometry": {
            "description": "Fine-tune coordination geometry without major changes",
            "partial_t_range": (5.0, 10.0),
            "recommended_partial_t": 8.0,
            "fix_backbone": True,
            "allow_length_change": False,
            "num_timesteps": 150,
            "step_scale": 1.8,
            "gamma_0": 0.5,
        },
    }

    @classmethod
    def select_strategy(
        cls,
        coordination_delta: int,
        donor_change_needed: bool,
        pocket_expansion_needed: bool
    ) -> str:
        """
        Select the appropriate design strategy based on analysis.

        Returns:
            Strategy name from STRATEGIES
        """
        if pocket_expansion_needed:
            return "expand_pocket"

        if donor_change_needed:
            return "change_donor_type"

        if coordination_delta > 2:
            return "increase_coordination"
        elif coordination_delta < -2:
            return "decrease_coordination"
        elif abs(coordination_delta) <= 2 and coordination_delta != 0:
            return "change_donor_type"  # Minor coordination change with donor swap
        else:
            return "refine_geometry"

    @classmethod
    def get_strategy_params(cls, strategy_name: str) -> Dict[str, Any]:
        """Get parameters for a design strategy."""
        return cls.STRATEGIES.get(strategy_name, cls.STRATEGIES["refine_geometry"])


class AIRecommender:
    """
    Enhanced hybrid AI recommendation engine.

    Combines rule-based chemistry knowledge with optional LLM refinement.
    """

    def __init__(self, anthropic_api_key: Optional[str] = None):
        """
        Initialize the recommender.

        Args:
            anthropic_api_key: Optional API key for Claude integration
        """
        self.has_llm = ANTHROPIC_AVAILABLE and anthropic_api_key is not None
        if self.has_llm:
            self.client = anthropic.Anthropic(api_key=anthropic_api_key)
        else:
            self.client = None

    def analyze_structure(
        self,
        pdb_content: str,
        current_analysis: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Analyze structure and extract relevant features for recommendation.

        Args:
            pdb_content: PDB file content
            current_analysis: Output from coordination geometry analyzer

        Returns:
            Dict with analysis features
        """
        if not current_analysis.get("success"):
            return current_analysis

        metal = current_analysis["metal"]["element"]
        coordination = current_analysis["coordination"]
        donor_analysis = current_analysis["donor_analysis"]

        # Get metal properties
        metal_props = MetalChemistryRules.get_metal_properties(metal)

        # Analyze donor compatibility
        incompatible_donors = []
        for donor_type in donor_analysis.get("donor_list", []):
            compatible, reason = MetalChemistryRules.is_compatible_donor(metal, donor_type)
            if not compatible:
                incompatible_donors.append({
                    "donor_type": donor_type,
                    "reason": reason,
                })

        return {
            "metal": metal,
            "metal_properties": asdict(metal_props) if metal_props else None,
            "coordination_number": coordination["number"],
            "geometry": coordination["geometry"],
            "donor_types": donor_analysis.get("types", {}),
            "incompatible_donors": incompatible_donors,
            "suggestions_from_analysis": current_analysis.get("suggestions", []),
        }

    def generate_baseline_parameters(
        self,
        current_analysis: Dict[str, Any],
        target_metal: str,
        user_goal: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Generate baseline RFD3 parameters using chemistry rules.

        Args:
            current_analysis: Output from coordination analysis
            target_metal: Target metal ion (e.g., "TB")
            user_goal: Optional natural language goal description

        Returns:
            Dict with RFD3 parameters and reasoning
        """
        if not current_analysis.get("success"):
            return {
                "success": False,
                "error": current_analysis.get("error", "Analysis failed"),
            }

        current_metal = current_analysis["metal"]["element"]
        current_coord = current_analysis["coordination"]["number"]

        # Get target metal properties
        target_props = MetalChemistryRules.get_metal_properties(target_metal)
        if target_props is None:
            return {
                "success": False,
                "error": f"Unknown target metal: {target_metal}",
            }

        # Calculate coordination change needed
        coord_analysis = MetalChemistryRules.calculate_coordination_delta(
            current_metal, target_metal, current_coord
        )

        # Analyze donor changes needed
        current_donors = current_analysis["donor_analysis"].get("donor_list", [])
        needs_donor_change = any(
            donor in target_props.avoid_donors for donor in current_donors
        )

        # Check if pocket expansion needed (lanthanides are larger)
        current_props = MetalChemistryRules.get_metal_properties(current_metal)
        needs_expansion = (
            target_props.ionic_radius > (current_props.ionic_radius if current_props else 0.8) + 0.3
        )

        # Select design strategy
        strategy = DesignStrategyRules.select_strategy(
            coord_analysis["delta"],
            needs_donor_change,
            needs_expansion
        )
        strategy_params = DesignStrategyRules.get_strategy_params(strategy)

        # Build coordinating residue list
        coord_residues = current_analysis["coordination"]["coordinating_atoms"]
        residue_ids = [f"A{atom['residue_number']}" for atom in coord_residues]

        # Generate RFD3 parameters
        rfd3_params = {
            "ligand": target_metal.upper(),
            "partial_t": strategy_params["recommended_partial_t"],
            "unindex": ",".join(residue_ids),
            "select_fixed_atoms": {
                res_id: "BKBN" for res_id in residue_ids
            } if strategy_params["fix_backbone"] else {},
            "num_timesteps": strategy_params["num_timesteps"],
            "step_scale": strategy_params["step_scale"],
            "gamma_0": strategy_params["gamma_0"],
            "num_designs": 10,
        }

        # Generate reasoning
        reasoning = []
        reasoning.append(f"Converting from {current_metal} to {target_metal}")
        reasoning.append(f"Strategy selected: {strategy} - {strategy_params['description']}")
        reasoning.append(coord_analysis["recommendation"])

        if needs_donor_change:
            reasoning.append(
                f"Donor change required: {target_metal} prefers {target_props.preferred_donors}"
            )

        if needs_expansion:
            reasoning.append(
                f"Pocket expansion needed: {target_metal} is larger "
                f"({target_props.ionic_radius}Å vs {current_props.ionic_radius if current_props else '?'}Å)"
            )

        return {
            "success": True,
            "parameters": rfd3_params,
            "strategy": strategy,
            "reasoning": reasoning,
            "coordination_analysis": coord_analysis,
            "target_properties": asdict(target_props),
            "evaluation_criteria": {
                "target_coordination": target_props.typical_coordination,
                "target_distance_range": target_props.bond_distance_range,
                "preferred_geometry": target_props.preferred_geometries,
                "preferred_donors": target_props.preferred_donors,
            },
        }

    async def refine_with_llm(
        self,
        baseline_params: Dict[str, Any],
        user_description: str,
        structure_context: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Refine baseline parameters using Claude API.

        Args:
            baseline_params: Output from generate_baseline_parameters
            user_description: User's natural language description
            structure_context: Additional structure context

        Returns:
            Dict with refined parameters
        """
        if not self.has_llm:
            return {
                **baseline_params,
                "llm_refined": False,
                "note": "LLM refinement not available - using baseline parameters",
            }

        # Build LLM prompt
        prompt = f"""You are an expert protein engineer specializing in metal binding site design.

## Task
Refine RFdiffusion3 (RFD3) parameters for a metal binding site redesign.

## User Goal
{user_description}

## Current Analysis
- Current metal: {baseline_params.get('coordination_analysis', {}).get('current_metal', 'Unknown')}
- Target metal: {baseline_params.get('parameters', {}).get('ligand', 'Unknown')}
- Coordination delta: {baseline_params.get('coordination_analysis', {}).get('delta', 0)}

## Structure Context
{json.dumps(structure_context, indent=2)}

## Baseline Parameters (from chemistry rules)
{json.dumps(baseline_params.get('parameters', {}), indent=2)}

## Reasoning
{chr(10).join(baseline_params.get('reasoning', []))}

## Your Task
1. Evaluate if the baseline parameters are appropriate for the user's goal
2. Suggest refinements if needed (adjust partial_t, num_timesteps, etc.)
3. Provide additional insights for achieving the user's goal

Respond with a JSON object containing:
- "refined_parameters": {{...}}  (adjusted RFD3 parameters)
- "refinements_made": ["list of changes and why"]
- "additional_suggestions": ["list of tips for better results"]
- "confidence": 0.0-1.0 (how confident you are in these parameters)

Important: Keep parameters within valid ranges:
- partial_t: 5.0-20.0
- num_timesteps: 50-500
- step_scale: 0.5-3.0
- gamma_0: 0.1-1.0
"""

        try:
            message = self.client.messages.create(
                model="claude-sonnet-4-20250514",
                max_tokens=1024,
                messages=[
                    {"role": "user", "content": prompt}
                ]
            )

            # Parse response
            response_text = message.content[0].text

            # Try to extract JSON from response
            import re
            json_match = re.search(r'\{[\s\S]*\}', response_text)
            if json_match:
                refined = json.loads(json_match.group())

                return {
                    **baseline_params,
                    "llm_refined": True,
                    "parameters": refined.get("refined_parameters", baseline_params["parameters"]),
                    "llm_refinements": refined.get("refinements_made", []),
                    "llm_suggestions": refined.get("additional_suggestions", []),
                    "llm_confidence": refined.get("confidence", 0.8),
                }
            else:
                return {
                    **baseline_params,
                    "llm_refined": False,
                    "note": "Could not parse LLM response - using baseline parameters",
                }

        except Exception as e:
            return {
                **baseline_params,
                "llm_refined": False,
                "error": f"LLM refinement failed: {str(e)}",
            }


# Convenience function for synchronous use
def generate_rfd3_parameters(
    current_analysis: Dict[str, Any],
    target_metal: str,
    user_description: Optional[str] = None,
    anthropic_api_key: Optional[str] = None
) -> Dict[str, Any]:
    """
    Generate RFD3 parameters for metal binding site redesign.

    This is the main entry point for the AI recommendation system.

    Args:
        current_analysis: Output from coordination geometry analyzer
        target_metal: Target metal ion (e.g., "TB", "GD", "EU")
        user_description: Optional natural language description
        anthropic_api_key: Optional API key for LLM refinement

    Returns:
        Dict with RFD3 parameters and reasoning
    """
    recommender = AIRecommender(anthropic_api_key)

    return recommender.generate_baseline_parameters(
        current_analysis,
        target_metal,
        user_description
    )
