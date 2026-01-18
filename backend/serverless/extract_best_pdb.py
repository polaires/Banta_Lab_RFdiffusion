#!/usr/bin/env python3
"""Extract best PDB from design result JSON."""
import json
import sys

def main():
    input_file = sys.argv[1] if len(sys.argv) > 1 else "/tmp/design_result.json"
    output_file = sys.argv[2] if len(sys.argv) > 2 else "citrate_tb_best_v6.pdb"

    with open(input_file) as f:
        data = json.load(f)

    if data.get("status") != "COMPLETED":
        print(f"Error: {data.get('error', 'Unknown')}")
        sys.exit(1)

    result = data["output"]["result"]
    designs = result["designs"]
    best_idx = result["best_design"]
    best = designs[best_idx]

    print(f"Best design: #{best_idx + 1}")
    print(f"Score: {best['score']}")

    if "pipeline_validation" in best:
        pv = best["pipeline_validation"]
        print(f"Pipeline passed: {best.get('pipeline_passed', 'N/A')}")
        print(f"Pipeline score: {best.get('pipeline_score', 'N/A')}")
        if pv.get("clashes"):
            clashes = pv["clashes"]
            print(f"Clashes: {clashes['total_clash_count']}")
            print(f"Worst overlap: {clashes['worst_overlap']:.2f}Å")
        if pv.get("geometry"):
            geom = pv["geometry"]
            print(f"Min C-M distance: {geom['min_carbon_metal_distance']:.2f}Å")

    with open(output_file, "w") as pf:
        pf.write(best["pdb_content"])
    print(f"Saved to {output_file}")

if __name__ == "__main__":
    main()
