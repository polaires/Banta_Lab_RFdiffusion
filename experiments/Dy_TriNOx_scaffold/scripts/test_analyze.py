#!/usr/bin/env python3
"""Test analyze API on a single design."""

import json
import subprocess
from pathlib import Path

PDB_FILE = Path("/mnt/g/Github_local_repo/Banta_Lab_RFdiffusion/experiments/Dy_TriNOx_scaffold/outputs_v4/v4_small_000.pdb")
API_URL = "http://localhost:8000/runsync"

# Read PDB
with open(PDB_FILE) as f:
    pdb_content = f.read()

print(f"PDB file size: {len(pdb_content)} chars")
print(f"First 200 chars:\n{pdb_content[:200]}")

# Analyze
params = {
    "task": "analyze",
    "pdb_content": pdb_content,
    "metal_type": "DY",
    "ligand_name": "UNL"
}

payload = json.dumps({"input": params})
print(f"\nPayload size: {len(payload)} chars")

result = subprocess.run(
    ["curl", "-s", "-X", "POST", API_URL, "-H", "Content-Type: application/json", "-d", payload],
    capture_output=True, text=True, timeout=120
)

print(f"\nStdout length: {len(result.stdout)}")
print(f"Stderr: {result.stderr[:500] if result.stderr else 'None'}")

try:
    response = json.loads(result.stdout)
    print(f"\nResponse status: {response.get('status')}")

    if response.get("status") == "COMPLETED":
        output = response.get("output", {})
        result_data = output.get("result", {})
        print(f"\nDesign type: {result_data.get('design_type')}")
        print(f"Detected metal: {result_data.get('detected_metal')}")
        print(f"Detected ligand: {result_data.get('detected_ligand')}")

        analyses = result_data.get("analyses", {})
        print(f"\nAnalyses available: {list(analyses.keys())}")

        # Metal coordination
        metal = analyses.get("metal_coordination", {})
        print(f"\nMetal coordination:")
        print(f"  CN: {metal.get('coordination_number')}")
        print(f"  Distances: {metal.get('coord_distances')}")

        # Structure confidence
        conf = analyses.get("structure_confidence", {})
        print(f"\nStructure confidence:")
        print(f"  pLDDT mean: {conf.get('plddt_mean')}")
    else:
        print(f"\nError: {response.get('error', response.get('output', {}).get('error'))}")

except json.JSONDecodeError as e:
    print(f"\nJSON decode error: {e}")
    print(f"Raw response: {result.stdout[:1000]}")
