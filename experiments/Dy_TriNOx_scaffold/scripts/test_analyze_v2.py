#!/usr/bin/env python3
"""Test analyze API and examine output format."""

import json
import subprocess
from pathlib import Path

PDB_FILE = Path("/mnt/g/Github_local_repo/Banta_Lab_RFdiffusion/experiments/Dy_TriNOx_scaffold/outputs_v4/v4_small_000.pdb")
API_URL = "http://localhost:8000/runsync"

# Read PDB
with open(PDB_FILE) as f:
    pdb_content = f.read()

# Analyze with target_ligands specified
params = {
    "task": "analyze",
    "pdb_content": pdb_content,
    "target_ligands": ["DY", "UNL"]
}

payload = json.dumps({"input": params})

result = subprocess.run(
    ["curl", "-s", "-X", "POST", API_URL, "-H", "Content-Type: application/json", "-d", payload],
    capture_output=True, text=True, timeout=120
)

try:
    response = json.loads(result.stdout)
    print(f"Response status: {response.get('status')}")
    print(f"\nFull response:")
    print(json.dumps(response, indent=2))

except json.JSONDecodeError as e:
    print(f"JSON decode error: {e}")
    print(f"Raw response: {result.stdout[:2000]}")
