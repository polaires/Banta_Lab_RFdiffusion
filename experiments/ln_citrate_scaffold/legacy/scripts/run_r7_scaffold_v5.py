"""
Round 7v5: Corrected scaffolding - use pdb_content not pdb_base64

The API handler expects 'pdb_content' (raw PDB string) not 'pdb_base64'.
When pdb_content is provided, it's written to a temp file and set as spec["input"].
"""

import requests
import json
import base64
from pathlib import Path
import time

API_URL = "http://localhost:8000/runsync"
BASE_DIR = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\experiments\ln_citrate_scaffold")
INPUT_PDB = BASE_DIR / "inputs" / "citrate_ln_only.pdb"
OUTPUT_DIR = BASE_DIR / "outputs" / "round_07_scaffold"

def load_pdb_content(pdb_path: Path) -> str:
    """Load PDB as string content."""
    with open(pdb_path, 'r') as f:
        return f.read()

def run_rfd3(config: dict, name: str) -> dict:
    """Submit RFD3 job."""
    # Load PDB as string content (not base64!)
    pdb_content = load_pdb_content(INPUT_PDB)

    payload = {
        "input": {
            "task": "rfd3",
            "pdb_content": pdb_content,  # Raw string, not base64!
            **config
        }
    }

    print(f"\n{'='*60}")
    print(f"Running: {name}")
    print(f"Contig: {config.get('contig')}")
    print(f"CFG: {config.get('use_classifier_free_guidance', False)}, scale: {config.get('cfg_scale', 'N/A')}")
    print("="*60)

    try:
        response = requests.post(API_URL, json=payload, timeout=600)
        result = response.json()

        print(f"Status: {result.get('status')}")

        if result.get('status') == 'failed':
            output = result.get('output', {})
            error = output.get('error', str(result))
            print(f"Error: {error[:300]}")
        return result

    except Exception as e:
        print(f"Exception: {e}")
        return {"error": str(e)}

def save_results(result: dict, name: str):
    """Save RFD3 output PDBs."""
    output = result.get("output", result)  # Handle both formats
    res = output.get("result", {})

    if "designs" in res:
        for i, design in enumerate(res["designs"]):
            pdb_path = OUTPUT_DIR / f"{name}_{i:03d}.pdb"
            with open(pdb_path, 'w') as f:
                f.write(design["content"])
            print(f"Saved: {pdb_path.name}")

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Configs using correct format
    # Ligand chains: X for TB (residue 1), L for citrate (residue 1)
    configs = [
        {
            "name": "r7v5a_hbond_buried",
            "config": {
                "contig": "110-130",
                "ligand": "CIT,TB",
                "num_designs": 8,
                "step_scale": 1.5,
                "num_timesteps": 200,
                "select_fixed_atoms": {
                    "X1": "all",
                    "L1": "all"
                },
                "select_buried": {
                    "X1": "all"
                },
                "select_partially_buried": {
                    "L1": "O2,O5,O7"
                },
                "select_hbond_acceptor": {
                    "L1": "O1,O2,O3,O4,O5,O6,O7"
                },
                "select_hbond_donor": {
                    "L1": "O7"
                },
                "infer_ori_strategy": "com",
                "use_classifier_free_guidance": True,
                "cfg_scale": 2.0,
                "plddt_enhanced": True
            }
        },
        {
            "name": "r7v5b_fully_buried",
            "config": {
                "contig": "110-130",
                "ligand": "CIT,TB",
                "num_designs": 8,
                "step_scale": 1.5,
                "num_timesteps": 200,
                "select_fixed_atoms": {
                    "X1": "all",
                    "L1": "all"
                },
                "select_buried": {
                    "X1": "all",
                    "L1": "all"
                },
                "select_hbond_acceptor": {
                    "L1": "O1,O2,O3,O4,O5,O6"
                },
                "select_hbond_donor": {
                    "L1": "O7"
                },
                "infer_ori_strategy": "com",
                "use_classifier_free_guidance": True,
                "cfg_scale": 2.0,
                "plddt_enhanced": True
            }
        },
        {
            "name": "r7v5c_larger_scaffold",
            "config": {
                "contig": "130-160",
                "ligand": "CIT,TB",
                "num_designs": 8,
                "step_scale": 1.5,
                "num_timesteps": 200,
                "select_fixed_atoms": {
                    "X1": "all",
                    "L1": "all"
                },
                "select_buried": {
                    "X1": "all"
                },
                "select_partially_buried": {
                    "L1": "O2,O5,O7"
                },
                "select_hbond_acceptor": {
                    "L1": "O1,O2,O3,O4,O5,O6,O7"
                },
                "select_hbond_donor": {
                    "L1": "O7"
                },
                "infer_ori_strategy": "com",
                "use_classifier_free_guidance": True,
                "cfg_scale": 2.0,
                "plddt_enhanced": True
            }
        },
        {
            "name": "r7v5d_high_cfg",
            "config": {
                "contig": "110-130",
                "ligand": "CIT,TB",
                "num_designs": 8,
                "step_scale": 1.5,
                "num_timesteps": 200,
                "select_fixed_atoms": {
                    "X1": "all",
                    "L1": "all"
                },
                "select_buried": {
                    "X1": "all"
                },
                "select_hbond_acceptor": {
                    "L1": "O1,O2,O3,O4,O5,O6,O7"
                },
                "select_hbond_donor": {
                    "L1": "O7"
                },
                "infer_ori_strategy": "com",
                "use_classifier_free_guidance": True,
                "cfg_scale": 3.0,
                "plddt_enhanced": True
            }
        }
    ]

    # Test API health
    print("Testing API health...")
    try:
        health_response = requests.post(
            API_URL,
            json={"input": {"task": "health"}},
            timeout=30
        )
        print(f"Health: {health_response.json().get('output', {}).get('result', {}).get('healthy')}")
    except Exception as e:
        print(f"API not available: {e}")
        return

    # Run each config
    results = {}

    for item in configs:
        name = item["name"]
        result = run_rfd3(item["config"], name)
        results[name] = result

        status = result.get("status") or result.get("output", {}).get("status")
        if status == "completed":
            save_results(result, name)

        time.sleep(2)

    # Summary
    print("\n" + "="*60)
    print("ROUND 7v5 SUMMARY")
    print("="*60)

    for name, result in results.items():
        status = result.get("status") or result.get("output", {}).get("status", "UNKNOWN")
        print(f"{name}: {status}")

    print(f"\nOutputs saved to: {OUTPUT_DIR}")

    # Save configs
    config_path = OUTPUT_DIR / "r7v5_configs_used.json"
    with open(config_path, 'w') as f:
        json.dump(configs, f, indent=2)

if __name__ == "__main__":
    main()
