"""
Test R7 scaffolding v3: Generate protein around ligand using conditioning only
"""

import requests
import json
import base64
from pathlib import Path

API_URL = "http://localhost:8000/runsync"
BASE_DIR = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\experiments\ln_citrate_scaffold")
INPUT_PDB = BASE_DIR / "inputs" / "citrate_ln_only.pdb"

def load_pdb_base64(pdb_path: Path) -> str:
    with open(pdb_path, 'rb') as f:
        return base64.b64encode(f.read()).decode('utf-8')

pdb_b64 = load_pdb_base64(INPUT_PDB)

# Test configurations
tests = [
    {
        "name": "test1_simple_with_ligand",
        "payload": {
            "input": {
                "task": "rfd3",
                "pdb_base64": pdb_b64,
                "contig": "100-120",
                "num_designs": 1,
                "step_scale": 1.5,
            }
        }
    },
    {
        "name": "test2_buried_conditioning",
        "payload": {
            "input": {
                "task": "rfd3",
                "pdb_base64": pdb_b64,
                "contig": "100-120",
                "num_designs": 1,
                "step_scale": 1.5,
                "select_buried": {"L": "all"},
            }
        }
    },
    {
        "name": "test3_hbond_conditioning",
        "payload": {
            "input": {
                "task": "rfd3",
                "pdb_base64": pdb_b64,
                "contig": "100-120",
                "num_designs": 1,
                "step_scale": 1.5,
                "select_hbond_acceptor": {"L": "O1,O2,O3,O4,O5,O6"},
            }
        }
    },
    {
        "name": "test4_all_conditioning_with_cfg",
        "payload": {
            "input": {
                "task": "rfd3",
                "pdb_base64": pdb_b64,
                "contig": "100-120",
                "num_designs": 1,
                "step_scale": 1.5,
                "select_buried": {"L": "all"},
                "select_hbond_acceptor": {"L": "O1,O2,O3,O4,O5,O6"},
                "select_hbond_donor": {"L": "O7"},
                "use_classifier_free_guidance": True,
                "cfg_scale": 2.0,
            }
        }
    },
]

for test in tests:
    print(f"\n{'='*60}")
    print(f"Testing: {test['name']}")
    print("="*60)

    try:
        response = requests.post(API_URL, json=test['payload'], timeout=300)
        result = response.json()

        print(f"Status: {result.get('status')}")

        if result.get('status') == 'FAILED':
            output = result.get('output', {})
            error = output.get('error', 'No error message')
            print(f"Error: {error[:200] if error else 'None'}")
        elif result.get('status') == 'COMPLETED':
            print("SUCCESS!")
            # Check if output has PDB
            output = result.get('output', {})
            if 'pdb_base64' in output:
                print("PDB generated successfully")

    except Exception as e:
        print(f"Exception: {e}")
