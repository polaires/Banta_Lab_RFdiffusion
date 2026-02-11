"""Test motif scaffolding with a coordinating Asp residue."""

import requests
import json
from pathlib import Path
import math

API_URL = "http://localhost:8000/runsync"
INPUT_DIR = Path("G:/Github_local_repo/Banta_Lab_RFdiffusion/experiments/Dy_TriNOx_scaffold/inputs")
OUTPUT_DIR = Path("G:/Github_local_repo/Banta_Lab_RFdiffusion/experiments/Dy_TriNOx_scaffold/outputs_v2")
OUTPUT_DIR.mkdir(exist_ok=True)

# Load input PDB with Asp motif
input_pdb = INPUT_DIR / "Dy_TriNOx_with_coord_Asp.pdb"
with open(input_pdb, 'r') as f:
    pdb_content = f.read()

# Contig: Metal (X1), Ligand (L1), Asp motif (M1), and new protein
# Format: include fixed chains, then scaffold
# M1 = the Asp at the coordination site

payload = {
    "input": {
        "task": "rfd3",
        "pdb_content": pdb_content,
        # Include Dy (X1), TriNOx (L1), Asp motif (M1), then add protein scaffold
        # Contig format for scaffolding around motif
        "contig": "X1,L1,5-15,M1,5-15",
        "ligand": "DY,UNL",
        "num_designs": 3,
        "select_fixed_atoms": {
            "X1": "all",
            "L1": "all",
            "M1": "all",  # Fix the Asp motif
        },
        "select_buried": {"X1": "all"},
        "use_classifier_free_guidance": True,
        "cfg_scale": 2.5,
    }
}

print("Testing motif scaffolding with fixed Asp at coordination site...")
print(f"Contig: {payload['input']['contig']}")

response = requests.post(API_URL, json=payload, timeout=300)
result = response.json()

print(f"\nResponse status: {result.get('status')}")

if result.get("status") == "FAILED":
    print(f"Error: {result.get('error')}")
    print(f"Full: {json.dumps(result, indent=2)[:2000]}")
elif result.get("status") == "COMPLETED":
    output = result.get("output", {})
    res = output.get("result", {})
    designs = res.get("designs", [])
    print(f"Generated {len(designs)} designs")

    for i, design_data in enumerate(designs):
        design = design_data["content"]
        output_file = OUTPUT_DIR / f"test_motif_scaffold_{i}.pdb"
        with open(output_file, 'w') as f:
            f.write(design)
        print(f"Saved: {output_file}")
else:
    print(f"Full: {json.dumps(result, indent=2)[:2000]}")
