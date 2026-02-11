#!/usr/bin/env python3
"""Prefilter v4 designs using static analysis before full validation."""

import json
import subprocess
from pathlib import Path

# Configuration
BASE_DIR = Path("/mnt/g/Github_local_repo/Banta_Lab_RFdiffusion/experiments/Dy_TriNOx_scaffold")
OUTPUT_DIR = BASE_DIR / "outputs_v4"
API_URL = "http://localhost:8000/runsync"

# Prefiltering thresholds
MIN_METAL_CN = 7  # Dy needs 8-9, TriNOx provides ~7, want at least 1 protein donor
MIN_PLDDT = 0.6   # Relaxed threshold for RFD3 output (before ESMFold)


def analyze_design(pdb_content: str, design_name: str) -> dict:
    """Run static analysis on a design."""
    params = {
        "task": "analyze",
        "pdb_content": pdb_content,
        "metal_type": "DY",
        "ligand_name": "UNL"
    }

    payload = json.dumps({"input": params})

    result = subprocess.run(
        ["curl", "-s", "-X", "POST", API_URL, "-H", "Content-Type: application/json", "-d", payload],
        capture_output=True, text=True, timeout=120
    )

    try:
        response = json.loads(result.stdout)
        if response.get("status") == "COMPLETED":
            return response.get("output", {}).get("result", {})
        else:
            print(f"  {design_name}: ERROR - {response.get('output', {}).get('error', 'Unknown')[:80]}")
            return {}
    except json.JSONDecodeError as e:
        print(f"  {design_name}: Parse error - {e}")
        return {}


def passes_prefilter(result: dict, design_name: str) -> bool:
    """Check if design passes prefilter criteria."""
    if not result:
        return False

    analyses = result.get("analyses", {})

    # Get metal coordination
    metal = analyses.get("metal_coordination", {})
    cn = metal.get("coordination_number", 0)

    # Get structure confidence (pLDDT)
    confidence = analyses.get("structure_confidence", {})
    plddt = confidence.get("plddt_mean", 0)

    # Get ligand binding info
    ligand = analyses.get("ligand_binding", {})
    protein_contacts = ligand.get("protein_contacts", 0) if ligand else 0

    # Check criteria
    passes = cn >= MIN_METAL_CN and plddt >= MIN_PLDDT

    status = "PASS" if passes else "FAIL"
    print(f"  {design_name}: {status} | CN={cn}, pLDDT={plddt:.2f}, protein_contacts={protein_contacts}")

    return passes


def main():
    # Find all design files
    design_files = sorted(OUTPUT_DIR.glob("v4_*.pdb"))

    if not design_files:
        print("No design files found in outputs_v4/")
        return

    print(f"Found {len(design_files)} designs to analyze")
    print(f"Prefilter criteria: CN >= {MIN_METAL_CN}, pLDDT >= {MIN_PLDDT}")
    print("=" * 60)

    passing_designs = []
    failing_designs = []

    for pdb_file in design_files:
        design_name = pdb_file.stem

        with open(pdb_file) as f:
            pdb_content = f.read()

        result = analyze_design(pdb_content, design_name)

        if passes_prefilter(result, design_name):
            passing_designs.append({
                "file": str(pdb_file),
                "name": design_name,
                "result": result
            })
        else:
            failing_designs.append(design_name)

    print("=" * 60)
    print(f"PREFILTER RESULTS:")
    print(f"  Passing: {len(passing_designs)}/{len(design_files)}")
    print(f"  Failing: {len(failing_designs)}/{len(design_files)}")

    if passing_designs:
        print(f"\nDesigns passing prefilter:")
        for d in passing_designs:
            print(f"  - {d['name']}")

        # Save passing designs list
        passing_file = OUTPUT_DIR / "prefilter_passing.json"
        with open(passing_file, "w") as f:
            json.dump(passing_designs, f, indent=2, default=str)
        print(f"\nSaved passing designs to {passing_file}")


if __name__ == "__main__":
    main()
