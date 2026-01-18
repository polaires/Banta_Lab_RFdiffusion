#!/usr/bin/env python3
"""Analyze design results and save best one."""
import json
import sys

def analyze_designs(json_path, output_pdb_path):
    with open(json_path) as f:
        d = json.load(f)

    if d.get("status") != "COMPLETED":
        print(f"Error: {d}")
        return None

    designs = d["output"]["result"]["designs"]
    print(f"Total designs: {len(designs)}\n")

    best_idx = -1
    best_score = -999

    for i, design in enumerate(designs):
        v = design.get("validation", {})
        score = v.get("quality_score", 0)
        coord = v.get("coordination_number", 0)
        ligand_coord = v.get("ligand_coordination", 0)
        protein_coord = v.get("protein_coordination", 0)
        hbonds = v.get("ligand_protein_hbond_count", 0)

        donors = v.get("donor_residues", [])
        donor_str = ", ".join([f"{r['resname']} {r['chain']}:{r['resnum']}" for r in donors[:4]])

        ligand_d = v.get("ligand_donors", [])
        ligand_dists = [f"{l['atom']}:{l['distance']:.2f}A" for l in ligand_d[:4]]

        print(f"Design {i+1}: Score {score}, Coord {coord}/9 (Lig:{ligand_coord} Prot:{protein_coord}), H-bonds:{hbonds}")
        print(f"  Protein donors: {donor_str}")
        print(f"  Ligand O dists: {', '.join(ligand_dists)}")
        print()

        if score > best_score:
            best_score = score
            best_idx = i

    print(f"\nBest design: #{best_idx + 1} with score {best_score}")

    # Save best design
    best_pdb = designs[best_idx]["pdb_content"]
    with open(output_pdb_path, "w") as f:
        f.write(best_pdb)
    print(f"Saved best design to: {output_pdb_path}")

    return best_idx

if __name__ == "__main__":
    json_path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/citrate_tb_designs.json"
    output_path = sys.argv[2] if len(sys.argv) > 2 else "/tmp/citrate_tb_best.pdb"
    analyze_designs(json_path, output_path)
