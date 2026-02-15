"""
Round 7v4: Correct ligand-centric scaffolding using proper chain+residue format

Key fixes from R6 configs:
1. Use "L1" and "X1" format (chain + residue number) not just "L"
2. Include "ligand": "CIT,TB" parameter
3. Use infer_ori_strategy for positioning
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

def load_pdb_base64(pdb_path: Path) -> str:
    with open(pdb_path, 'rb') as f:
        return base64.b64encode(f.read()).decode('utf-8')

def run_rfd3(config: dict, name: str) -> dict:
    """Submit RFD3 job."""
    pdb_b64 = load_pdb_base64(INPUT_PDB)

    payload = {
        "input": {
            "task": "rfd3",
            "pdb_base64": pdb_b64,
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

        if result.get('status') == 'FAILED':
            output = result.get('output', {})
            error = output.get('error', str(output))
            print(f"Error: {error[:300]}")
        return result

    except Exception as e:
        print(f"Exception: {e}")
        return {"error": str(e)}

def save_results(result: dict, name: str):
    """Save RFD3 output PDBs."""
    output = result.get("output", {})

    if "pdb_files" in output:
        for i, pdb_b64 in enumerate(output["pdb_files"]):
            pdb_path = OUTPUT_DIR / f"{name}_{i:03d}.pdb"
            pdb_content = base64.b64decode(pdb_b64)
            with open(pdb_path, 'wb') as f:
                f.write(pdb_content)
            print(f"Saved: {pdb_path.name}")
    elif "pdb_base64" in output:
        pdb_path = OUTPUT_DIR / f"{name}_000.pdb"
        pdb_content = base64.b64decode(output["pdb_base64"])
        with open(pdb_path, 'wb') as f:
            f.write(pdb_content)
        print(f"Saved: {pdb_path.name}")

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Configs using correct L1/X1 format (ligand at residue 1 in input PDB)
    configs = [
        {
            "name": "r7v4a_hbond_buried",
            "config": {
                "contig": "110-130",
                "ligand": "CIT,TB",
                "num_designs": 8,
                "step_scale": 1.5,
                "num_timesteps": 200,
                "select_fixed_atoms": {
                    "X1": "all",   # TB at chain X, residue 1
                    "L1": "all"    # Citrate at chain L, residue 1
                },
                "select_buried": {
                    "X1": "all"    # Bury TB
                },
                "select_partially_buried": {
                    "L1": "O2,O5,O7"  # Partial burial for citrate arms
                },
                "select_hbond_acceptor": {
                    "L1": "O1,O2,O3,O4,O5,O6,O7"  # All citrate oxygens
                },
                "select_hbond_donor": {
                    "L1": "O7"    # Citrate hydroxyl
                },
                "infer_ori_strategy": "com",
                "use_classifier_free_guidance": True,
                "cfg_scale": 2.0,
                "plddt_enhanced": True
            }
        },
        {
            "name": "r7v4b_fully_buried",
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
                    "L1": "all"    # Fully bury citrate
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
            "name": "r7v4c_larger_scaffold",
            "config": {
                "contig": "130-160",   # Larger protein for bigger pocket
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
            "name": "r7v4d_high_cfg",
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
                "cfg_scale": 3.0,   # Higher CFG for stronger conditioning
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

        if result.get("status") == "COMPLETED":
            save_results(result, name)

        time.sleep(2)

    # Summary
    print("\n" + "="*60)
    print("ROUND 7v4 SUMMARY")
    print("="*60)

    for name, result in results.items():
        status = result.get("status", "UNKNOWN")
        if "error" in result:
            status = f"FAILED"
        print(f"{name}: {status}")

    print(f"\nOutputs saved to: {OUTPUT_DIR}")

    # Save configs for reference
    config_path = OUTPUT_DIR / "r7v4_configs_used.json"
    with open(config_path, 'w') as f:
        json.dump(configs, f, indent=2)

if __name__ == "__main__":
    main()
