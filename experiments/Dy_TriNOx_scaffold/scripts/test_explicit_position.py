"""Test explicit origin positioning to place protein near Dy apical site."""

import requests
import json
from pathlib import Path
import math

API_URL = "http://localhost:8000/runsync"
INPUT_DIR = Path("G:/Github_local_repo/Banta_Lab_RFdiffusion/experiments/Dy_TriNOx_scaffold/inputs")
OUTPUT_DIR = Path("G:/Github_local_repo/Banta_Lab_RFdiffusion/experiments/Dy_TriNOx_scaffold/outputs_v2")
OUTPUT_DIR.mkdir(exist_ok=True)

# Load input PDB
input_pdb = INPUT_DIR / "Dy_TriNOx_split.pdb"
with open(input_pdb, 'r') as f:
    pdb_content = f.read()

# Dy is at (-0.001, 7.832, 3.094)
# Apical site is in +Z direction
# Position protein origin above Dy at Z ~ 6-7

# Test with explicit ori_token
payload = {
    "input": {
        "task": "rfd3",
        "pdb_content": pdb_content,
        "contig": "X1,L1,/0,40-60",
        "ligand": "DY,UNL",
        "num_designs": 3,
        "select_fixed_atoms": {
            "X1": "all",
            "L1": "all",
        },
        # Position protein above Dy apical site
        "ori_token": [0.0, 8.0, 6.0],  # Above Dy in Z direction
        "select_buried": {"X1": "all"},
        "use_classifier_free_guidance": True,
        "cfg_scale": 2.5,
    }
}

print("Testing explicit origin positioning...")
print(f"Original Dy position: (-0.001, 7.832, 3.094)")
print(f"ori_token: [0.0, 8.0, 6.0] - placing protein above Dy")

response = requests.post(API_URL, json=payload, timeout=300)
result = response.json()

print(f"\nResponse status: {result.get('status')}")

if result.get("status") == "COMPLETED":
    output = result.get("output", {})
    res = output.get("result", {})
    designs = res.get("designs", [])
    print(f"Generated {len(designs)} designs")

    for i, design_data in enumerate(designs):
        design = design_data["content"]
        output_file = OUTPUT_DIR / f"test_explicit_pos_{i}.pdb"
        with open(output_file, 'w') as f:
            f.write(design)
        print(f"Saved: {output_file}")

        # Analyze coordination
        print(f"\n  Analyzing design {i}:")
        lines = design.split('\n')

        # Find Dy position
        dy_pos = None
        protein_atoms = []

        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                if res_name == "DY" or atom_name == "DY":
                    dy_pos = (x, y, z)
                elif res_name not in ["DY", "UNL"]:
                    protein_atoms.append({
                        "name": atom_name,
                        "res_name": res_name,
                        "x": x, "y": y, "z": z
                    })

        if dy_pos:
            print(f"    Dy position: ({dy_pos[0]:.2f}, {dy_pos[1]:.2f}, {dy_pos[2]:.2f})")

            # Find closest protein atom
            min_dist = float('inf')
            closest = None
            for a in protein_atoms:
                d = math.sqrt((a["x"]-dy_pos[0])**2 + (a["y"]-dy_pos[1])**2 + (a["z"]-dy_pos[2])**2)
                if d < min_dist:
                    min_dist = d
                    closest = a

            print(f"    Closest protein atom: {closest['res_name']} {closest['name']} at {min_dist:.2f} A")

            # Count atoms within coordination distance
            within_3 = sum(1 for a in protein_atoms
                         if math.sqrt((a["x"]-dy_pos[0])**2 + (a["y"]-dy_pos[1])**2 + (a["z"]-dy_pos[2])**2) <= 3.0)
            within_5 = sum(1 for a in protein_atoms
                         if math.sqrt((a["x"]-dy_pos[0])**2 + (a["y"]-dy_pos[1])**2 + (a["z"]-dy_pos[2])**2) <= 5.0)

            print(f"    Atoms within 3A: {within_3}")
            print(f"    Atoms within 5A: {within_5}")

else:
    print(f"Error: {result.get('error', 'Unknown')}")
    print(f"Full: {json.dumps(result, indent=2)[:2000]}")
