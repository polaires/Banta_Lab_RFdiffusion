"""
Case Study: Converting Rubredoxin (1BRF) to Terbium Binding

This script demonstrates the AI-guided protein engineering workflow:
1. Fetch 1BRF from RCSB
2. Analyze the iron binding site
3. Generate RFD3 parameters for terbium conversion
4. Show evaluation criteria for output assessment

Run this script to test the AI protein engineering system:
    python test_1brf_case_study.py
"""

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from utils.pdb_fetch import fetch_pdb, parse_pdb_content
from utils.coordination import (
    analyze_coordination_geometry,
    suggest_lanthanide_conversion,
    MetalChemistryRules,
)
from ai.recommender import generate_rfd3_parameters, MetalChemistryRules as AIMetalRules


def print_section(title):
    """Print a section header."""
    print("\n" + "=" * 60)
    print(f"  {title}")
    print("=" * 60)


def main():
    print_section("Case Study: Rubredoxin (1BRF) → Terbium Binding")

    # Step 1: Fetch 1BRF
    print_section("Step 1: Fetching 1BRF from RCSB")

    result = fetch_pdb("1BRF")
    if not result["success"]:
        print(f"ERROR: Could not fetch 1BRF - {result['error']}")
        return

    pdb_content = result["content"]
    print(f"✓ Fetched 1BRF from {result['source']}")

    # Parse basic info
    info = parse_pdb_content(pdb_content)
    print(f"  Title: {info['title']}")
    print(f"  Chains: {info['chains']}")
    print(f"  Atoms: {info['num_atoms']}")
    print(f"  Residues: {info['num_residues']}")
    print(f"  Metals found: {info['metals']}")

    # Step 2: Analyze metal binding site
    print_section("Step 2: Analyzing Iron Binding Site")

    # In 1BRF, the iron is typically:
    # - Chain A (or only chain)
    # - Residue name: FE or FE2
    # - Coordinated by Cys5, Cys8, Cys38, Cys41

    metal_info = info['metals'][0] if info['metals'] else None
    if metal_info:
        metal_chain = metal_info['chain']
        metal_residue = metal_info['residue']
        metal_resnum = int(metal_info['res_num'])
    else:
        # Default for 1BRF
        metal_chain = "A"
        metal_residue = "FE"
        metal_resnum = 1
        print("  (Using default metal location: A:FE1)")

    analysis = analyze_coordination_geometry(
        pdb_content,
        metal_chain,
        metal_residue,
        metal_resnum,
        distance_cutoff=3.0
    )

    if not analysis.get("success"):
        print(f"ERROR: Could not analyze metal site - {analysis.get('error')}")
        # Try with FE2
        metal_residue = "FE2"
        analysis = analyze_coordination_geometry(
            pdb_content, metal_chain, metal_residue, metal_resnum
        )

    if analysis.get("success"):
        print(f"\n✓ Metal Analysis Results:")
        print(f"  Metal: {analysis['metal']['element']} (Chain {analysis['metal']['chain']})")
        print(f"  Coordination Number: {analysis['coordination']['number']}")
        print(f"  Geometry: {analysis['coordination']['geometry']}")
        print(f"  Geometry RMSD: {analysis['coordination']['geometry_rmsd']}Å")

        print(f"\n  Coordinating Atoms:")
        for atom in analysis['coordination']['coordinating_atoms']:
            print(f"    - {atom['chain']}:{atom['residue']}{atom['residue_number']}:{atom['atom']} "
                  f"({atom['donor_type']}) - {atom['distance']}Å")

        print(f"\n  Donor Type Summary:")
        for dtype, count in analysis['donor_analysis']['types'].items():
            print(f"    - {dtype}: {count}")

        print(f"\n  Bond Distance Statistics:")
        print(f"    Average: {analysis['bond_analysis']['average_distance']}Å")
        print(f"    Range: {analysis['bond_analysis']['min_distance']} - {analysis['bond_analysis']['max_distance']}Å")

        if analysis.get('suggestions'):
            print(f"\n  Suggestions:")
            for sug in analysis['suggestions']:
                print(f"    ⚠ {sug}")
    else:
        print(f"ERROR: {analysis.get('error')}")
        return

    # Step 3: Generate RFD3 Parameters for Terbium
    print_section("Step 3: AI Parameter Recommendation for Terbium Conversion")

    target_metal = "TB"

    # Get target metal properties
    tb_props = AIMetalRules.get_metal_properties(target_metal)
    print(f"\n  Target Metal Properties ({tb_props.name}):")
    print(f"    Typical Coordination: {tb_props.typical_coordination}")
    print(f"    Preferred Geometries: {tb_props.preferred_geometries}")
    print(f"    Preferred Donors: {tb_props.preferred_donors}")
    print(f"    Avoid Donors: {tb_props.avoid_donors}")
    print(f"    Ionic Radius: {tb_props.ionic_radius}Å")
    print(f"    Bond Distance Range: {tb_props.bond_distance_range}Å")

    # Generate parameters
    recommendation = generate_rfd3_parameters(
        current_analysis=analysis,
        target_metal=target_metal,
        user_description="Convert rubredoxin iron binding site to bind terbium with high coordination"
    )

    if recommendation.get("success"):
        print(f"\n✓ AI Recommendation:")
        print(f"  Strategy: {recommendation['strategy']}")

        print(f"\n  Reasoning:")
        for reason in recommendation['reasoning']:
            print(f"    • {reason}")

        print(f"\n  Coordination Analysis:")
        coord_analysis = recommendation['coordination_analysis']
        print(f"    Current: {coord_analysis['current_coordination']}")
        print(f"    Target Range: {coord_analysis['target_coordination_range']}")
        print(f"    Delta: {coord_analysis['delta']}")
        print(f"    Recommendation: {coord_analysis['recommendation']}")

        print(f"\n  Recommended RFD3 Parameters:")
        params = recommendation['parameters']
        for key, value in params.items():
            print(f"    {key}: {value}")

        print(f"\n  Evaluation Criteria:")
        eval_criteria = recommendation['evaluation_criteria']
        print(f"    Target Coordination: {eval_criteria['target_coordination']}")
        print(f"    Target Distance Range: {eval_criteria['target_distance_range']}Å")
        print(f"    Preferred Geometry: {eval_criteria['preferred_geometry']}")
        print(f"    Preferred Donors: {eval_criteria['preferred_donors']}")
    else:
        print(f"ERROR: {recommendation.get('error')}")

    # Step 4: Show API Usage
    print_section("Step 4: API Usage Example")

    print("""
To use via the FastAPI backend:

1. Start the server:
   uvicorn main:app --host 0.0.0.0 --port 8000

2. Fetch and analyze 1BRF:
   POST /api/analyze/metal-binding
   {
       "pdb_id": "1BRF",
       "metal_chain": "A",
       "metal_residue": "FE",
       "metal_resnum": 1,
       "distance_cutoff": 3.0
   }

3. Get AI recommendations:
   POST /api/ai/recommend-parameters
   {
       "pdb_id": "1BRF",
       "metal_chain": "A",
       "metal_residue": "FE",
       "metal_resnum": 1,
       "target_metal": "TB",
       "user_description": "Convert to terbium with 8-9 coordination"
   }

4. Run RFD3 design with recommended parameters:
   POST /api/rfd3/design
   {
       "pdb_content": "<1BRF pdb content>",
       "ligand": "TB",
       "partial_t": 12.0,
       "unindex": "A5,A8,A38,A41",
       "select_fixed_atoms": {"A5": "BKBN", "A8": "BKBN", "A38": "BKBN", "A41": "BKBN"},
       "num_timesteps": 200,
       "step_scale": 1.5,
       "gamma_0": 0.6,
       "num_designs": 10
   }

5. Evaluate design output:
   POST /api/evaluate/design
   {
       "pdb_content": "<design output>",
       "target_metal": "TB",
       "metal_chain": "A",
       "metal_resnum": 1
   }
""")

    print_section("Case Study Complete")
    print("\nThe AI protein engineering system has:")
    print("  ✓ Analyzed the iron binding site in rubredoxin")
    print("  ✓ Identified the need to change from 4-coord Fe to 8-9 coord Tb")
    print("  ✓ Recognized that Cys (S-donor) should be replaced with Asp/Glu (O-donor)")
    print("  ✓ Generated appropriate RFD3 parameters for the conversion")
    print("  ✓ Provided evaluation criteria for assessing the design output")


if __name__ == "__main__":
    main()
