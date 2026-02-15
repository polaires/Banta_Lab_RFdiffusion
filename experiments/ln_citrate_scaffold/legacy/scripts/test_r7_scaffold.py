"""
Test R7 scaffolding with detailed error reporting
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

# Read the input PDB to understand structure
print("Input PDB content:")
with open(INPUT_PDB) as f:
    print(f.read())

# Test simple contig first
pdb_b64 = load_pdb_base64(INPUT_PDB)

# Try different contig formats
contigs_to_test = [
    # Format 1: Chain letter, residue range
    "L1,X1,/0,100-120",
    # Format 2: Just chains
    "L,X,/0,100-120",
    # Format 3: Skip the HETATM chains entirely, just generate protein
    "100-120",
]

for contig in contigs_to_test:
    print(f"\n{'='*60}")
    print(f"Testing contig: {contig}")
    print("="*60)

    payload = {
        "input": {
            "task": "rfd3",
            "pdb_base64": pdb_b64,
            "contig": contig,
            "num_designs": 1,
            "step_scale": 1.5,
        }
    }

    # Only add select_fixed_atoms if we reference those chains
    if "L" in contig or "X" in contig:
        payload["input"]["select_fixed_atoms"] = {"X": "all", "L": "all"}

    try:
        response = requests.post(API_URL, json=payload, timeout=120)
        result = response.json()

        print(f"Status: {result.get('status')}")

        if result.get('status') == 'FAILED':
            output = result.get('output', {})
            error = output.get('error', 'No error message')
            print(f"Error: {error}")
        elif result.get('status') == 'COMPLETED':
            print("SUCCESS!")
            break

    except Exception as e:
        print(f"Exception: {e}")
