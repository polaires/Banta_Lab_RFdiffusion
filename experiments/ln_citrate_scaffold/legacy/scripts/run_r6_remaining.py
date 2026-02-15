#!/usr/bin/env python3
"""
Run remaining Round 6 configs: r6_c (dimer) and r6_d (extended).
"""

import json
import time
from datetime import datetime
from pathlib import Path

import requests

API_URL = "http://localhost:8000/runsync"
EXPERIMENT_DIR = Path(__file__).parent.parent
OUTPUT_DIR = EXPERIMENT_DIR / "outputs" / "round_06"
INPUT_DIR = EXPERIMENT_DIR / "inputs"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_pdb(pdb_path: Path) -> str:
    with open(pdb_path) as f:
        return f.read()


def run_rfd3(config: dict, pdb_content: str, name: str, num_designs: int = 32) -> dict:
    clean_config = {k: v for k, v in config.items() if not k.startswith("_")}

    payload = {
        "input": {
            "task": "rfd3",
            "pdb_content": pdb_content,
            "num_designs": num_designs,
            "contig": clean_config.get("contig"),
            "ligand": clean_config.get("ligand"),
            "select_fixed_atoms": clean_config.get("select_fixed_atoms"),
            "select_buried": clean_config.get("select_buried"),
            "select_partially_buried": clean_config.get("select_partially_buried"),
            "select_hbond_acceptor": clean_config.get("select_hbond_acceptor"),
            "select_hbond_donor": clean_config.get("select_hbond_donor"),
            "hotspots": clean_config.get("hotspots"),
            "use_classifier_free_guidance": clean_config.get("use_classifier_free_guidance", True),
            "cfg_scale": clean_config.get("cfg_scale", 2.0),
            "sym": clean_config.get("sym"),
            "diffusion_batch_size": clean_config.get("diffusion_batch_size"),
            "step_scale": clean_config.get("step_scale", 1.5),
            "num_timesteps": clean_config.get("num_timesteps", 200),
            "infer_ori_strategy": clean_config.get("infer_ori_strategy"),
            "is_non_loopy": clean_config.get("is_non_loopy", False),
            "plddt_enhanced": clean_config.get("plddt_enhanced", True),
        }
    }

    payload["input"] = {k: v for k, v in payload["input"].items() if v is not None}

    print(f"\n  Sending request for {name}...")
    print(f"  Config: contig={clean_config.get('contig')}")
    print(f"  Symmetry: {clean_config.get('sym', 'None')}")
    print(f"  num_designs={num_designs}")

    try:
        response = requests.post(API_URL, json=payload, timeout=1800)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.Timeout:
        return {"error": "Timeout after 1800s"}
    except requests.exceptions.RequestException as e:
        return {"error": str(e)}


def save_designs(result: dict, config_name: str) -> list:
    saved = []

    if "error" in result:
        print(f"  ERROR: {result['error']}")
        return saved

    output = result.get("output", {})
    result_data = output.get("result", {})
    designs = result_data.get("designs", [])

    if not designs:
        designs = output.get("designs", [])
    if not designs:
        designs = result.get("designs", [])

    if not designs:
        print(f"  No designs in response.")
        return saved

    for i, design in enumerate(designs):
        pdb_content = design.get("content") or design.get("pdb_content") or design.get("pdb")
        if not pdb_content:
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
    print("Round 6 - Remaining Configs (r6_c dimer, r6_d extended)")
    print("="*60)
    print(f"\nTimestamp: {datetime.now().isoformat()}")

    # Load input PDB
    input_path = INPUT_DIR / "citrate_ln_only.pdb"
    pdb_content = load_pdb(input_path)
    print(f"Input: {input_path.name} ({len(pdb_content)} bytes)")

    # Configs to run
    configs = {
        "r6_c_dimer_symmetric": {
            "contig": "50-70,/0,50-70",
            "ligand": "CIT,TB",
            "sym": "C2",
            "diffusion_batch_size": 1,
            "select_fixed_atoms": {"X1": "all", "L1": "all"},
            "select_buried": {"X1": "all", "L1": "O2,O5,O7"},
            "select_hbond_acceptor": {"L1": "O1,O3,O4,O6"},
            "hotspots": {"X1": "all"},
            "infer_ori_strategy": "com",
            "use_classifier_free_guidance": True,
            "cfg_scale": 2.0,
            "step_scale": 1.5,
            "num_timesteps": 200,
            "plddt_enhanced": True
        },
        "r6_d_denovo_extended": {
            "contig": "110-130",
            "ligand": "CIT,TB",
            "select_fixed_atoms": {"X1": "all", "L1": "all"},
            "select_buried": {"X1": "all"},
            "select_partially_buried": {"L1": "O2,O5,O7"},
            "select_hbond_acceptor": {"L1": "O1,O2,O3,O4,O5,O6,O7"},
            "select_hbond_donor": {"L1": "O7"},
            "infer_ori_strategy": "com",
            "use_classifier_free_guidance": True,
            "cfg_scale": 2.5,
            "step_scale": 1.5,
            "num_timesteps": 200,
            "is_non_loopy": True,
            "plddt_enhanced": True
        }
    }

    results = {}
    for config_name, config in configs.items():
        print(f"\n{'='*60}")
        print(f"Running: {config_name}")
        print(f"{'='*60}")

        start_time = time.time()
        result = run_rfd3(config, pdb_content, config_name, num_designs=32)
        elapsed = time.time() - start_time

        saved = save_designs(result, config_name)
        results[config_name] = {
            "designs_saved": len(saved),
            "elapsed_seconds": round(elapsed, 1)
        }

        print(f"\n  Completed: {len(saved)} designs in {elapsed:.1f}s")

    print(f"\n{'='*60}")
    print("REMAINING CONFIGS COMPLETE")
    print(f"{'='*60}")
    for name, res in results.items():
        print(f"  {name}: {res['designs_saved']} designs")


if __name__ == "__main__":
    main()
