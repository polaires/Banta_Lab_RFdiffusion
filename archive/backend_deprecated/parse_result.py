#!/usr/bin/env python3
import json
import sys

filename = sys.argv[1] if len(sys.argv) > 1 else "/tmp/lanthanide_result3.json"
with open(filename) as f:
    data = json.load(f)

result = data.get("output", {}).get("result", {})

print("=== Lanthanide Design Test Results ===")
print(f"Approach: {result.get('approach')}")
print(f"Metal: {result.get('metal')}")
print(f"Template type: {result.get('template_type')}")
print(f"Target coordination: {result.get('target_coordination')}")
print()

# Check best design metrics
if result.get("designs"):
    design = result["designs"][0]
    metrics = design.get("metrics", {})
    print("=== Design Metrics ===")
    print(f"Coordination number: {metrics.get('coordination_number')}")
    print(f"Chain A donors: {metrics.get('chain_a_donors')}")
    print(f"Chain B donors: {metrics.get('chain_b_donors')}")
    seq_id = metrics.get("sequence_identity", 0)
    print(f"Sequence identity: {seq_id:.1f}%")
    print(f"Is heterodimer: {metrics.get('is_heterodimer')}")
    print()

# Check validation results
validation = result.get("validation", {})
if validation:
    print("=== Lanthanide Validation ===")
    print(f"Validation success: {validation.get('success')}")
    print(f"Coordination number: {validation.get('coordination_number')}")
    print(f"Quality score: {validation.get('quality_score')}")
    print(f"Quality rating: {validation.get('quality_rating')}")
    print(f"Carboxylate count: {validation.get('carboxylate_count')}")
    print(f"Chain contribution: {validation.get('chain_contribution')}")
    if validation.get("suggestions"):
        print("Suggestions:")
        for s in validation["suggestions"][:3]:
            print(f"  - {s[:80]}...")
    print()

    # TEBL details
    tebl = validation.get("tebl_details", {})
    if tebl:
        print("=== TEBL Analysis (from validation) ===")
        print(f"TEBL ready: {tebl.get('has_antenna')}")
        print(f"Trp chain: {tebl.get('trp_chain')}")
        print(f"Trp residue: {tebl.get('trp_resnum')}")
        print(f"Trp-metal distance: {tebl.get('trp_metal_distance')} A")
        print(f"Energy transfer efficiency: {tebl.get('energy_transfer_efficiency')}")

# Check TEBL analysis
tebl_analysis = result.get("tebl_analysis", {})
if tebl_analysis:
    print()
    print("=== TEBL Signal Prediction ===")
    print(f"Has antenna: {tebl_analysis.get('has_antenna')}")
    print(f"Signal strength: {tebl_analysis.get('signal_strength')}")
    print(f"Predicted efficiency: {tebl_analysis.get('predicted_efficiency')}")
    if tebl_analysis.get("best_antenna"):
        best = tebl_analysis["best_antenna"]
        print(f"Best antenna: {best.get('chain')}{best.get('residue_number')} at {best.get('distance')} A")
