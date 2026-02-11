#!/usr/bin/env python3
"""Test CIF to PDB conversion for RF3 outputs."""
import json
import requests

# Read the actual RF3 CIF output
with open("../../experiments/ln_citrate_scaffold/scaffolding_3c9h/outputs/batch_01_pdbs/seq_0019.pdb") as f:
    cif_content = f.read()

print(f"Loaded CIF: {len(cif_content)} chars")
print(f"First 100 chars: {cif_content[:100]}")

# Test via API
params = {
    "input": {
        "task": "analyze_design",
        "pdb_content": cif_content,
        "metal_type": "TB"
    }
}

resp = requests.post("http://localhost:8000/runsync", json=params, timeout=120)
result = resp.json()

output = result.get("output", result)
print(f"\nStatus: {output.get('status')}")
print(f"CN: {output.get('result', {}).get('metrics', {}).get('coordination_number')}")
print(f"Auto-detected: {output.get('result', {}).get('auto_detected', {})}")
