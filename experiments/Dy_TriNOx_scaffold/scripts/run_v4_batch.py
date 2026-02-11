#!/usr/bin/env python3
"""Run full v4 RFD3 batch for Dy-TriNOx scaffold design."""

import json
import subprocess
import os
from pathlib import Path

# Configuration
BASE_DIR = Path("/mnt/g/Github_local_repo/Banta_Lab_RFdiffusion/experiments/Dy_TriNOx_scaffold")
OUTPUT_DIR = BASE_DIR / "outputs_v4"
INPUT_PDB = BASE_DIR / "inputs" / "Dy_TriNOx_UNL.pdb"
API_URL = "http://localhost:8000/runsync"

# V4 configs with Round 7b-style settings
CONFIGS = [
    {"name": "v4_small", "contig": "X1,L1,/0,80-100", "num_designs": 25},
    {"name": "v4_medium", "contig": "X1,L1,/0,100-120", "num_designs": 25},
    {"name": "v4_large", "contig": "X1,L1,/0,120-140", "num_designs": 25},
    {"name": "v4_extended", "contig": "X1,L1,/0,140-160", "num_designs": 25},
]

BASE_PARAMS = {
    "task": "rfd3",
    "use_classifier_free_guidance": True,
    "cfg_scale": 2.0,
    "step_scale": 1.5,
    "num_timesteps": 200,
    "select_fixed_atoms": {"X1": "all", "L1": "all"},
    "select_buried": {"X1": "all"},
    "select_hbond_acceptor": {"L1": "O1,O2,O3"},
    "infer_ori_strategy": "com"
}


def run_rfd3_batch(name: str, contig: str, num_designs: int, pdb_content: str) -> list:
    """Run a single RFD3 batch and return designs."""
    params = {**BASE_PARAMS, "contig": contig, "num_designs": num_designs, "pdb_content": pdb_content}

    payload = json.dumps({"input": params})

    result = subprocess.run(
        ["curl", "-s", "-X", "POST", API_URL, "-H", "Content-Type: application/json", "-d", payload],
        capture_output=True, text=True
    )

    try:
        response = json.loads(result.stdout)
        if response.get("status") == "COMPLETED":
            return response.get("output", {}).get("result", {}).get("designs", [])
        else:
            print(f"  ERROR: {response.get('output', {}).get('error', 'Unknown error')[:100]}")
            return []
    except json.JSONDecodeError as e:
        print(f"  ERROR: Failed to parse response: {e}")
        return []


def main():
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Read input PDB
    with open(INPUT_PDB) as f:
        pdb_content = f.read()

    total_designs = 0

    for config in CONFIGS:
        name = config["name"]
        contig = config["contig"]
        num_designs = config["num_designs"]

        print(f"\n{'='*60}")
        print(f"Running {name}: {contig}, {num_designs} designs")
        print(f"{'='*60}")

        designs = run_rfd3_batch(name, contig, num_designs, pdb_content)

        # Save designs
        for i, design in enumerate(designs):
            filename = OUTPUT_DIR / f"{name}_{i:03d}.pdb"
            with open(filename, "w") as f:
                f.write(design["content"])

        print(f"  Saved {len(designs)} designs to {OUTPUT_DIR}")
        total_designs += len(designs)

    print(f"\n{'='*60}")
    print(f"TOTAL: {total_designs} designs saved to {OUTPUT_DIR}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
