#!/usr/bin/env python3
"""
Round 6 execution script - Two approaches to achieve CN=8-9.

Approaches:
- r6_a: Partial diffusion from R5 best (partial_t=0.5)
- r6_b: Aggressive partial diffusion (partial_t=0.7)
- r6_c: Homodimer with C2 symmetry
- r6_d: Extended de novo with stronger H-bond conditioning

Key parameters:
- step_scale: 1.5 (η from paper)
- num_timesteps: 200
- cfg_scale: 2.0 or 2.5
- 32 designs per config
"""

import json
import sys
import time
from datetime import datetime
from pathlib import Path

import requests

API_URL = "http://localhost:8000/runsync"
EXPERIMENT_DIR = Path(__file__).parent.parent
OUTPUT_DIR = EXPERIMENT_DIR / "outputs" / "round_06"
INPUT_DIR = EXPERIMENT_DIR / "inputs"

# Ensure output dir exists
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_pdb(pdb_path: Path) -> str:
    """Load PDB file content."""
    with open(pdb_path) as f:
        return f.read()


def run_rfd3(config: dict, pdb_content: str, name: str, num_designs: int = 32) -> dict:
    """Run RFD3 design via API with full parameter support."""

    # Build the RFD3 request - remove description fields
    clean_config = {k: v for k, v in config.items() if not k.startswith("_")}

    payload = {
        "input": {
            "task": "rfd3",
            "pdb_content": pdb_content,
            "num_designs": num_designs,
            "contig": clean_config.get("contig"),
            "ligand": clean_config.get("ligand"),
            # Fixed atoms
            "select_fixed_atoms": clean_config.get("select_fixed_atoms"),
            # Burial conditioning
            "select_buried": clean_config.get("select_buried"),
            "select_partially_buried": clean_config.get("select_partially_buried"),
            # H-bond conditioning
            "select_hbond_acceptor": clean_config.get("select_hbond_acceptor"),
            "select_hbond_donor": clean_config.get("select_hbond_donor"),
            # Hotspots
            "hotspots": clean_config.get("hotspots"),
            # CFG parameters (CRITICAL)
            "use_classifier_free_guidance": clean_config.get("use_classifier_free_guidance", True),
            "cfg_scale": clean_config.get("cfg_scale", 2.0),
            # Partial diffusion
            "partial_t": clean_config.get("partial_t"),
            # Symmetry
            "sym": clean_config.get("sym"),
            "diffusion_batch_size": clean_config.get("diffusion_batch_size"),
            # Diffusion parameters
            "step_scale": clean_config.get("step_scale", 1.5),
            "num_timesteps": clean_config.get("num_timesteps", 200),
            # Origin
            "infer_ori_strategy": clean_config.get("infer_ori_strategy"),
            # Other params
            "is_non_loopy": clean_config.get("is_non_loopy", False),
            "plddt_enhanced": clean_config.get("plddt_enhanced", True),
        }
    }

    # Remove None values
    payload["input"] = {k: v for k, v in payload["input"].items() if v is not None}

    print(f"\n  Sending request for {name}...")
    print(f"  Config: contig={clean_config.get('contig')}")
    print(f"  Partial_t: {clean_config.get('partial_t', 'N/A')}")
    print(f"  Symmetry: {clean_config.get('sym', 'None')}")
    print(f"  CFG: scale={clean_config.get('cfg_scale', 2.0)}, step_scale={clean_config.get('step_scale', 1.5)}")
    print(f"  num_designs={num_designs}")

    try:
        # Extended timeout for complex designs
        response = requests.post(API_URL, json=payload, timeout=1800)  # 30 min timeout
        response.raise_for_status()
        result = response.json()
        return result
    except requests.exceptions.Timeout:
        return {"error": "Timeout after 1800s"}
    except requests.exceptions.RequestException as e:
        return {"error": str(e)}


def save_designs(result: dict, config_name: str) -> list:
    """Save generated designs to PDB files."""
    saved = []

    if "error" in result:
        print(f"  ERROR: {result['error']}")
        return saved

    # Extract designs from response
    output = result.get("output", {})
    result_data = output.get("result", {})
    designs = result_data.get("designs", [])

    if not designs:
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
    print("="*70)
    print("Ln-Citrate Scaffold Design - Round 6 (Partial Diffusion + Dimer)")
    print("="*70)
    print(f"\nTimestamp: {datetime.now().isoformat()}")
    print(f"Output dir: {OUTPUT_DIR}")
    print("\nApproaches:")
    print("  - r6_a: Partial diffusion from R5 best (partial_t=0.5)")
    print("  - r6_b: Aggressive partial diffusion (partial_t=0.7)")
    print("  - r6_c: Homodimer with C2 symmetry")
    print("  - r6_d: Extended de novo with stronger H-bond conditioning")

    # Load configs
    config_path = EXPERIMENT_DIR / "configs" / "round_06" / "r6_configs.json"
    with open(config_path) as f:
        all_configs = json.load(f)

    # Remove meta
    configs = {k: v for k, v in all_configs.items() if not k.startswith("_")}
    n_designs = all_configs.get("_meta", {}).get("n_designs_per_config", 32)

    print(f"\nConfigs to run: {list(configs.keys())}")
    print(f"Designs per config: {n_designs}")
    print(f"Total designs: {len(configs) * n_designs}")

    # Track results
    results_summary = {
        "round": "6",
        "timestamp": datetime.now().isoformat(),
        "configs": {}
    }

    for config_name, config in configs.items():
        print(f"\n{'='*60}")
        print(f"Running: {config_name}")
        print(f"{'='*60}")

        # Determine input PDB based on config
        input_path_str = config.get("input", "../inputs/citrate_ln_only.pdb")
        # Resolve relative to configs/round_06 directory
        config_dir = EXPERIMENT_DIR / "configs" / "round_06"
        input_path = (config_dir / input_path_str).resolve()

        if not input_path.exists():
            print(f"  ERROR: Input file not found: {input_path}")
            results_summary["configs"][config_name] = {"error": f"Input not found: {input_path}"}
            continue

        pdb_content = load_pdb(input_path)
        print(f"  Input: {input_path.name} ({len(pdb_content)} bytes)")

        start_time = time.time()
        result = run_rfd3(config, pdb_content, config_name, num_designs=n_designs)
        elapsed = time.time() - start_time

        saved = save_designs(result, config_name)

        results_summary["configs"][config_name] = {
            "designs_saved": len(saved),
            "elapsed_seconds": round(elapsed, 1),
            "input": input_path.name,
            "partial_t": config.get("partial_t"),
            "sym": config.get("sym"),
        }

        print(f"\n  Completed: {len(saved)} designs saved in {elapsed:.1f}s")

    # Save summary
    summary_path = OUTPUT_DIR / "run_summary.json"
    with open(summary_path, "w") as f:
        json.dump(results_summary, f, indent=2)

    print(f"\n{'='*70}")
    print("ROUND 6 COMPLETE")
    print(f"{'='*70}")
    print(f"Summary saved: {summary_path}")

    total_designs = sum(c.get("designs_saved", 0) for c in results_summary["configs"].values())
    print(f"Total designs generated: {total_designs}")


if __name__ == "__main__":
    main()
