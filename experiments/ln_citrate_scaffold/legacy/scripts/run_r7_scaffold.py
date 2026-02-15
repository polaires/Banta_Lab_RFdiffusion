"""
Round 7: Motif Scaffolding for Tb-Citrate Binding Protein

Strategy: Fix TB-citrate + H-bond coordinating residues, scaffold around them.
This preserves the citrate binding pocket while exploring protein structure.
"""

import requests
import json
import base64
from pathlib import Path
import time

# Configuration
API_URL = "http://localhost:8000/runsync"
BASE_DIR = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\experiments\ln_citrate_scaffold")
INPUT_PDB = BASE_DIR / "inputs" / "tb_citrate_motif_scaffold.pdb"
OUTPUT_DIR = BASE_DIR / "outputs" / "round_07_scaffold"

def load_pdb_base64(pdb_path: Path) -> str:
    """Load PDB and encode as base64."""
    with open(pdb_path, 'rb') as f:
        return base64.b64encode(f.read()).decode('utf-8')

def run_rfd3_scaffold(config: dict, name: str) -> dict:
    """Submit scaffolding job to RFD3 API."""

    # Load input PDB
    pdb_b64 = load_pdb_base64(INPUT_PDB)

    payload = {
        "input": {
            "task": "rfd3",
            "pdb_base64": pdb_b64,
            "contig": config["contig"],
            "num_designs": 8,
            "step_scale": config.get("step_scale", 1.5),
            "gamma_0": 0.6,
            "num_timesteps": 200,
        }
    }

    # Add conditioning parameters
    if "select_fixed_atoms" in config:
        payload["input"]["select_fixed_atoms"] = config["select_fixed_atoms"]

    if "select_hbond_acceptor" in config:
        payload["input"]["select_hbond_acceptor"] = config["select_hbond_acceptor"]

    if "select_hbond_donor" in config:
        payload["input"]["select_hbond_donor"] = config["select_hbond_donor"]

    if "select_buried" in config:
        payload["input"]["select_buried"] = config["select_buried"]

    if "unindex" in config:
        payload["input"]["unindex"] = config["unindex"]

    if config.get("use_classifier_free_guidance"):
        payload["input"]["use_classifier_free_guidance"] = True
        payload["input"]["cfg_scale"] = config.get("cfg_scale", 2.0)

    print(f"\n{'='*60}")
    print(f"Running: {name}")
    print(f"Contig: {config['contig']}")
    print(f"CFG: {config.get('use_classifier_free_guidance', False)}, scale: {config.get('cfg_scale', 'N/A')}")
    print(f"{'='*60}")

    try:
        response = requests.post(
            API_URL,
            json=payload,
            timeout=600  # 10 min timeout for scaffolding
        )

        if response.status_code == 200:
            result = response.json()
            print(f"Status: {result.get('status', 'unknown')}")
            return result
        else:
            print(f"Error: {response.status_code}")
            print(response.text)
            return {"error": response.text}

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
    # Ensure output directory exists
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load config
    config_path = BASE_DIR / "configs" / "round_07" / "r7_scaffold_config.json"
    with open(config_path) as f:
        full_config = json.load(f)

    configs = full_config["round_07_scaffold"]["configs"]

    # Test API health first
    print("Testing API health...")
    try:
        health_response = requests.post(
            API_URL,
            json={"input": {"task": "health"}},
            timeout=30
        )
        print(f"Health check: {health_response.json()}")
    except Exception as e:
        print(f"API not available: {e}")
        print("Please start Docker container first!")
        return

    # Run each scaffolding config
    results = {}

    for config in configs:
        name = config["name"]
        result = run_rfd3_scaffold(config, name)
        results[name] = result

        if "error" not in result:
            save_results(result, name)

        # Brief pause between jobs
        time.sleep(2)

    # Summary
    print("\n" + "="*60)
    print("ROUND 7 SCAFFOLDING SUMMARY")
    print("="*60)

    for name, result in results.items():
        status = "SUCCESS" if "error" not in result else "FAILED"
        print(f"{name}: {status}")

    print(f"\nOutputs saved to: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
