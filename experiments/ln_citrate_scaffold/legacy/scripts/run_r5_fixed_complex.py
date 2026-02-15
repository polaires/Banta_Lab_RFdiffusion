#!/usr/bin/env python3
"""
Round 5: RFD3 backbone generation with FIXED TB-citrate complex.

CRITICAL FIX: In Round 2b, citrate was NOT fixed, causing it to drift 17A away
from TB during diffusion. This round fixes BOTH TB and citrate to maintain
the original coordination geometry (TB-O distances 2.35-2.5A).
"""

import json
import sys
import time
from datetime import datetime
from pathlib import Path

import requests

API_URL = "http://localhost:8000/runsync"
EXPERIMENT_DIR = Path(__file__).parent.parent
CONFIG_FILE = EXPERIMENT_DIR / "configs" / "round_05" / "r5_configs.json"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs" / "round_05"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

DESIGNS_PER_CONFIG = 16


def load_pdb(pdb_path: Path) -> str:
    """Load PDB file content."""
    with open(pdb_path) as f:
        return f.read()


def run_rfd3(config: dict, config_name: str, pdb_content: str, design_idx: int) -> dict:
    """
    Run RFD3 via API.

    CRITICAL: select_fixed_atoms must include BOTH TB and CIT
    to maintain coordination geometry.
    """
    # Build payload from config
    payload = {
        "input": {
            "task": "rfd3",
            "pdb_content": pdb_content,
            "contig": config.get("contig"),
            "ligand": config.get("ligand"),
        }
    }

    # Add all conditioning parameters
    for key in ["select_fixed_atoms", "select_buried", "select_partially_buried",
                "select_hbond_acceptor", "select_hbond_donor",
                "infer_ori_strategy", "use_classifier_free_guidance", "cfg_scale",
                "is_non_loopy", "plddt_enhanced"]:
        if key in config and config[key]:
            payload["input"][key] = config[key]

    print(f"\n  Running RFD3 design {design_idx+1}/{DESIGNS_PER_CONFIG}...")
    print(f"  Fixed atoms: {config.get('select_fixed_atoms', {})}")

    try:
        response = requests.post(API_URL, json=payload, timeout=600)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.Timeout:
        return {"error": "Timeout after 600s"}
    except requests.exceptions.RequestException as e:
        return {"error": str(e)}


def save_result(result: dict, config_name: str, design_idx: int) -> dict:
    """Save RFD3 result to PDB file."""
    saved = {"config": config_name, "design_idx": design_idx}

    if "error" in result:
        print(f"  ERROR: {result['error']}")
        saved["error"] = result["error"]
        return saved

    # Extract PDB content
    output = result.get("output", {})
    result_data = output.get("result", {})

    pdb_content = None
    for key in ["pdb_content", "pdb", "structure", "content"]:
        if key in result_data:
            pdb_content = result_data[key]
            break

    if not pdb_content:
        predictions = result_data.get("predictions", [])
        if predictions and isinstance(predictions[0], dict):
            pdb_content = predictions[0].get("pdb_content") or predictions[0].get("content")

    if not pdb_content:
        print(f"  No structure in response")
        saved["error"] = "No structure returned"
        return saved

    # Save PDB
    name = f"{config_name}_{design_idx:03d}"
    pdb_file = OUTPUT_DIR / f"{name}.pdb"
    with open(pdb_file, "w") as f:
        f.write(pdb_content)

    saved["pdb"] = pdb_file.name
    print(f"  Saved: {pdb_file.name}")

    return saved


def main():
    print("="*60)
    print("Ln-Citrate Scaffold - Round 5: FIXED TB-Citrate Complex")
    print("="*60)
    print(f"\nTimestamp: {datetime.now().isoformat()}")
    print(f"Output dir: {OUTPUT_DIR}")
    print("\nCRITICAL FIX:")
    print("  - Round 2b had CIT NOT fixed -> drifted 17A from TB")
    print("  - Round 5 fixes BOTH TB AND CIT atoms")
    print("  - This maintains coordination geometry (TB-O: 2.35-2.5A)")

    # Check API health
    print("\nChecking API health...")
    try:
        health = requests.post(API_URL, json={"input": {"task": "health"}}, timeout=30)
        health_data = health.json()
        if health_data.get("output", {}).get("result", {}).get("healthy"):
            print("  API healthy")
        else:
            print("  WARNING: API may not be healthy")
    except Exception as e:
        print(f"  ERROR: Cannot connect to API: {e}")
        sys.exit(1)

    # Load configs
    with open(CONFIG_FILE) as f:
        all_configs = json.load(f)

    # Filter out metadata
    configs = {k: v for k, v in all_configs.items() if not k.startswith("_")}

    print(f"\nConfigs to run: {list(configs.keys())}")
    print(f"Designs per config: {DESIGNS_PER_CONFIG}")

    # Process each config
    all_results = {}
    total_designs = 0

    for config_name, config in configs.items():
        print(f"\n{'='*50}")
        print(f"Config: {config_name}")
        print(f"Description: {config.get('_description', 'N/A')}")

        # Load input PDB
        input_path = EXPERIMENT_DIR / "configs" / "round_05" / config["input"]
        if not input_path.exists():
            # Try relative to experiment dir
            input_path = EXPERIMENT_DIR / config["input"].lstrip("../")

        if not input_path.exists():
            print(f"  ERROR: Input file not found: {config['input']}")
            continue

        pdb_content = load_pdb(input_path)
        print(f"  Input: {input_path.name}")

        config_results = []

        for i in range(DESIGNS_PER_CONFIG):
            start_time = time.time()
            result = run_rfd3(config, config_name, pdb_content, i)
            elapsed = time.time() - start_time

            saved = save_result(result, config_name, i)
            saved["elapsed_seconds"] = round(elapsed, 1)
            config_results.append(saved)

            if "pdb" in saved:
                total_designs += 1

        all_results[config_name] = {
            "description": config.get("_description"),
            "hypothesis": config.get("_hypothesis"),
            "n_designs": len([r for r in config_results if "pdb" in r]),
            "results": config_results
        }

    # Save summary
    summary = {
        "round": 5,
        "task": "RFD3 with FIXED TB-citrate complex",
        "timestamp": datetime.now().isoformat(),
        "critical_fix": "Both TB and CIT now fixed to maintain coordination",
        "total_designs": total_designs,
        "configs": all_results
    }

    summary_path = OUTPUT_DIR / "run_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n{'='*60}")
    print("Round 5 Complete")
    print(f"{'='*60}")
    print(f"Total designs: {total_designs}")
    print(f"Summary: {summary_path}")
    print(f"\nNext: Analyze backbones for TB-citrate distance and coordination")


if __name__ == "__main__":
    main()
