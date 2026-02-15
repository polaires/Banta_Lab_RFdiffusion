#!/usr/bin/env python3
"""
Direct Round 2 execution script for Ln-citrate scaffold design.
Runs RFD3 via local Docker API.

Key changes from R1:
- Using simpler input (just TB + CIT)
- Focus on burial conditioning without explicit motif residues
- Let RFD3 figure out coordination
"""

import json
import os
import sys
import time
from datetime import datetime
from pathlib import Path

import requests

API_URL = "http://localhost:8000/runsync"
EXPERIMENT_DIR = Path(__file__).parent.parent
OUTPUT_DIR = EXPERIMENT_DIR / "outputs" / "round_02"
INPUT_DIR = EXPERIMENT_DIR / "inputs"

# Ensure output dir exists
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_pdb(pdb_path: Path) -> str:
    """Load PDB file content."""
    with open(pdb_path) as f:
        return f.read()


def run_rfd3(config: dict, pdb_content: str, name: str, num_designs: int = 8) -> dict:
    """Run RFD3 design via API."""

    # Build the RFD3 request
    # Remove description fields
    clean_config = {k: v for k, v in config.items() if not k.startswith("_")}

    # The API expects the config in a specific format
    # NOTE: Use 'num_designs' not 'n_designs'
    payload = {
        "input": {
            "task": "rfd3",
            "pdb_content": pdb_content,
            "num_designs": num_designs,
            "contig": clean_config.get("contig"),
            "ligand": clean_config.get("ligand"),
            "unindex": clean_config.get("unindex"),
            "select_fixed_atoms": clean_config.get("select_fixed_atoms"),
            "select_buried": clean_config.get("select_buried"),
            "select_partially_buried": clean_config.get("select_partially_buried"),
            "select_hbond_acceptor": clean_config.get("select_hbond_acceptor"),
            "ori_token": clean_config.get("ori_token"),
            "infer_ori_strategy": clean_config.get("infer_ori_strategy"),
            "is_non_loopy": clean_config.get("is_non_loopy", True),
            "plddt_enhanced": clean_config.get("plddt_enhanced", True),
        }
    }

    # Remove None values
    payload["input"] = {k: v for k, v in payload["input"].items() if v is not None}

    print(f"\n  Sending request for {name}...")
    print(f"  Config: contig={clean_config.get('contig')}, num_designs={num_designs}")
    print(f"  Fixed atoms: {clean_config.get('select_fixed_atoms')}")

    try:
        response = requests.post(API_URL, json=payload, timeout=600)
        response.raise_for_status()
        result = response.json()
        return result
    except requests.exceptions.Timeout:
        return {"error": "Timeout after 600s"}
    except requests.exceptions.RequestException as e:
        return {"error": str(e)}


def save_designs(result: dict, config_name: str) -> list:
    """Save generated designs to PDB files."""
    saved = []

    if "error" in result:
        print(f"  ERROR: {result['error']}")
        return saved

    # Extract designs from response
    # Format: output.result.designs[].content
    output = result.get("output", {})
    result_data = output.get("result", {})
    designs = result_data.get("designs", [])

    if not designs:
        # Try alternative paths
        designs = output.get("designs", [])

    if not designs:
        designs = result.get("designs", [])

    if not designs:
        print(f"  No designs in response.")
        print(f"  Response structure: {list(result.keys())}")
        if "output" in result:
            print(f"  output keys: {list(output.keys())}")
            if "result" in output:
                print(f"  result keys: {list(result_data.keys())}")
        return saved

    for i, design in enumerate(designs):
        # API returns 'content' field for PDB
        pdb_content = design.get("content") or design.get("pdb_content") or design.get("pdb")
        if not pdb_content:
            print(f"  Design {i} has no content. Keys: {design.keys()}")
            continue

        filename = f"{config_name}_{i:03d}.pdb"
        filepath = OUTPUT_DIR / filename

        with open(filepath, "w") as f:
            f.write(pdb_content)

        saved.append(filepath)
        print(f"  Saved: {filename}")

    return saved


def main():
    print("="*60)
    print("Ln-Citrate Scaffold Design - Round 2")
    print("="*60)
    print(f"\nTimestamp: {datetime.now().isoformat()}")
    print(f"Output dir: {OUTPUT_DIR}")
    print("\nKey changes from R1:")
    print("  - Using simpler input (just TB + CIT, no motif residues)")
    print("  - Focus on burial conditioning")
    print("  - Let RFD3 figure out coordination naturally")

    # Load configs
    config_path = EXPERIMENT_DIR / "configs" / "round_02" / "r2_all_configs.json"
    with open(config_path) as f:
        all_configs = json.load(f)

    # Remove meta
    configs = {k: v for k, v in all_configs.items() if not k.startswith("_")}

    print(f"\nConfigs to run: {list(configs.keys())}")

    # Load input PDB (simple: just citrate + TB)
    input_pdb = INPUT_DIR / "citrate_ln_only.pdb"
    pdb_content = load_pdb(input_pdb)
    print(f"Loaded input: {input_pdb.name} ({len(pdb_content)} bytes)")

    # Test API health first
    print("\nChecking API health...")
    try:
        health = requests.post(API_URL, json={"input": {"task": "health"}}, timeout=30)
        health_data = health.json()
        if health_data.get("output", {}).get("result", {}).get("healthy"):
            print("  API healthy, GPU available")
        else:
            print("  WARNING: API may not be healthy")
    except Exception as e:
        print(f"  ERROR: Cannot connect to API: {e}")
        print("  Make sure Docker is running!")
        sys.exit(1)

    # Run each config
    all_results = {}
    total_saved = 0
    num_designs = 8  # More designs this round for better statistics

    for config_name, config in configs.items():
        print(f"\n{'='*40}")
        print(f"Running: {config_name}")
        print(f"Description: {config.get('_description', 'N/A')}")
        print(f"Hypothesis: {config.get('_hypothesis', 'N/A')}")

        start_time = time.time()
        result = run_rfd3(config, pdb_content, config_name, num_designs=num_designs)
        elapsed = time.time() - start_time

        saved = save_designs(result, config_name)
        total_saved += len(saved)

        all_results[config_name] = {
            "n_generated": len(saved),
            "elapsed_seconds": round(elapsed, 1),
            "files": [str(f.name) for f in saved]
        }

        print(f"  Generated: {len(saved)} designs in {elapsed:.1f}s")

    # Save summary
    summary = {
        "round": 2,
        "timestamp": datetime.now().isoformat(),
        "num_designs_per_config": num_designs,
        "total_generated": total_saved,
        "configs": all_results
    }

    summary_path = OUTPUT_DIR / "run_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n{'='*60}")
    print("Round 2 Complete")
    print(f"{'='*60}")
    print(f"Total designs: {total_saved}")
    print(f"Summary: {summary_path}")
    print(f"\nNext: Analyze with simple_analyze.py --round 2")


if __name__ == "__main__":
    main()
