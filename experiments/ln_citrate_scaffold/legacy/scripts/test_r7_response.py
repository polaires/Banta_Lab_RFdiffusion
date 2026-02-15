"""Test R7 response structure"""

import requests
import json
from pathlib import Path

API_URL = "http://localhost:8000/runsync"
BASE_DIR = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\experiments\ln_citrate_scaffold")
INPUT_PDB = BASE_DIR / "inputs" / "citrate_ln_only.pdb"
OUTPUT_DIR = BASE_DIR / "outputs" / "round_07_scaffold"

def load_pdb_content(pdb_path: Path) -> str:
    with open(pdb_path, 'r') as f:
        return f.read()

pdb_content = load_pdb_content(INPUT_PDB)

payload = {
    "input": {
        "task": "rfd3",
        "pdb_content": pdb_content,
        "contig": "110-130",
        "ligand": "CIT,TB",
        "num_designs": 2,  # Just 2 for testing
        "step_scale": 1.5,
        "num_timesteps": 200,
        "select_fixed_atoms": {"X1": "all", "L1": "all"},
        "select_buried": {"X1": "all"},
        "select_hbond_acceptor": {"L1": "O1,O2,O3,O4,O5,O6"},
        "infer_ori_strategy": "com",
        "use_classifier_free_guidance": True,
        "cfg_scale": 2.0,
    }
}

print("Sending test request...")
response = requests.post(API_URL, json=payload, timeout=300)
result = response.json()

print("\n=== Response Structure ===")
print(f"Top-level keys: {result.keys()}")
print(f"Status: {result.get('status')}")

output = result.get("output", {})
print(f"\nOutput keys: {output.keys()}")
print(f"Output status: {output.get('status')}")

res = output.get("result", {})
print(f"\nResult keys: {res.keys()}")

if "designs" in res:
    designs = res["designs"]
    print(f"\nNumber of designs: {len(designs)}")
    if designs:
        print(f"First design keys: {designs[0].keys()}")
        content = designs[0].get("content", "")
        print(f"Content preview: {content[:200]}...")

        # Save to file
        for i, design in enumerate(designs):
            pdb_path = OUTPUT_DIR / f"r7_test_{i:03d}.pdb"
            with open(pdb_path, 'w') as f:
                f.write(design["content"])
            print(f"Saved: {pdb_path}")
else:
    print("\nNo 'designs' in result")
    print(f"Full result: {json.dumps(res, indent=2)[:500]}")
