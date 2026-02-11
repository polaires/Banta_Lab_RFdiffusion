"""Test hotspot conditioning on Dy-TriNOx."""

import requests
import json
from pathlib import Path

API_URL = "http://localhost:8000/runsync"
INPUT_DIR = Path("G:/Github_local_repo/Banta_Lab_RFdiffusion/experiments/Dy_TriNOx_scaffold/inputs")

# Load input PDB
input_pdb = INPUT_DIR / "Dy_TriNOx_split.pdb"
with open(input_pdb, 'r') as f:
    pdb_content = f.read()

# Test with hotspot on Dy
payload = {
    "input": {
        "task": "rfd3",
        "pdb_content": pdb_content,
        "contig": "X1,L1,/0,40-60",
        "ligand": "DY,UNL",
        "num_designs": 1,
        "select_fixed_atoms": {
            "X1": "all",
            "L1": "all",
        },
        "select_hotspots": {"X1": "1"},
        "select_buried": {"X1": "all"},
        "use_classifier_free_guidance": True,
        "cfg_scale": 2.5,
        "infer_ori_strategy": "hotspots",
    }
}

print("Testing hotspot on Dy...")
print(f"Payload: {json.dumps(payload, indent=2)}")

response = requests.post(API_URL, json=payload, timeout=300)
result = response.json()

print(f"\nResponse status: {result.get('status')}")
print(f"Full response: {json.dumps(result, indent=2)[:2000]}")
