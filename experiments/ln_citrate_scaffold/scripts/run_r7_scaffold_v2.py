"""
Round 7v2: Ligand-centric Scaffolding for Tb-Citrate Binding Protein

Strategy: Fix TB-citrate complex, let RFD3 build protein around it with:
1. RASA conditioning to bury the citrate
2. H-bond conditioning for citrate oxygens
3. CFG scale 2.0 for strong conditioning adherence

This is the "fixed ligand" workflow from RFD3 skill.
"""

import requests
import json
import base64
from pathlib import Path
import time

# Configuration
API_URL = "http://localhost:8000/runsync"
BASE_DIR = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\experiments\ln_citrate_scaffold")
INPUT_PDB = BASE_DIR / "inputs" / "citrate_ln_only.pdb"
OUTPUT_DIR = BASE_DIR / "outputs" / "round_07_scaffold"

def load_pdb_base64(pdb_path: Path) -> str:
    """Load PDB and encode as base64."""
    with open(pdb_path, 'rb') as f:
        return base64.b64encode(f.read()).decode('utf-8')

def run_rfd3_ligand_scaffold(config: dict, name: str) -> dict:
    """Submit ligand-centric scaffolding job to RFD3 API."""

    pdb_b64 = load_pdb_base64(INPUT_PDB)

    payload = {
        "input": {
            "task": "rfd3",
            "pdb_base64": pdb_b64,
            "contig": config["contig"],
            "num_designs": config.get("num_designs", 8),
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

    if config.get("use_classifier_free_guidance"):
        payload["input"]["use_classifier_free_guidance"] = True
        payload["input"]["cfg_scale"] = config.get("cfg_scale", 2.0)

    print(f"\n{'='*60}")
    print(f"Running: {name}")
    print(f"Contig: {config['contig']}")
    print(f"CFG: {config.get('use_classifier_free_guidance', False)}, scale: {config.get('cfg_scale', 'N/A')}")
    print(f"Conditioning: fixed={config.get('select_fixed_atoms')}")
    print(f"             buried={config.get('select_buried')}")
    print(f"             hbond_acceptor={config.get('select_hbond_acceptor')}")
    print(f"{'='*60}")

    try:
        response = requests.post(
            API_URL,
            json=payload,
            timeout=600
        )

        if response.status_code == 200:
            result = response.json()
            status = result.get('status', 'unknown')
            print(f"Status: {status}")

            # Check for errors in output
            output = result.get('output', {})
            if 'error' in output:
                print(f"Error in output: {output['error']}")

            return result
        else:
            print(f"HTTP Error: {response.status_code}")
            print(response.text[:500])
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
    else:
        print(f"No PDB output found in result")
        print(f"Output keys: {output.keys()}")

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Define scaffolding configurations
    # Input PDB has: TB on chain X, citrate on chain L
    configs = [
        {
            "name": "r7v2a_fixed_buried",
            "description": "Fix TB+citrate, bury citrate, H-bond conditioning",
            "contig": "X,L,/0,100-130",  # Fix X and L, generate 100-130 residue protein
            "select_fixed_atoms": {
                "X": "all",  # Fix TB
                "L": "all"   # Fix citrate
            },
            "select_buried": {
                "L": "all"   # Bury citrate
            },
            "select_hbond_acceptor": {
                "L": "O1,O2,O3,O4,O5,O6"  # Citrate carboxylate oxygens
            },
            "select_hbond_donor": {
                "L": "O7"    # Citrate hydroxyl
            },
            "use_classifier_free_guidance": True,
            "cfg_scale": 2.0,
            "num_designs": 8,
            "step_scale": 1.5
        },
        {
            "name": "r7v2b_fixed_partial_burial",
            "description": "Fix TB+citrate, partial burial for accessibility",
            "contig": "X,L,/0,100-130",
            "select_fixed_atoms": {
                "X": "all",
                "L": "all"
            },
            "select_partially_buried": {
                "L": "all"   # Partial burial for better access
            },
            "select_hbond_acceptor": {
                "L": "O1,O2,O3,O4,O5,O6"
            },
            "select_hbond_donor": {
                "L": "O7"
            },
            "use_classifier_free_guidance": True,
            "cfg_scale": 2.0,
            "num_designs": 8,
            "step_scale": 1.5
        },
        {
            "name": "r7v2c_larger_protein",
            "description": "Larger protein for bigger pocket",
            "contig": "X,L,/0,130-160",  # Larger protein
            "select_fixed_atoms": {
                "X": "all",
                "L": "all"
            },
            "select_buried": {
                "L": "all"
            },
            "select_hbond_acceptor": {
                "L": "O1,O2,O3,O4,O5,O6"
            },
            "select_hbond_donor": {
                "L": "O7"
            },
            "use_classifier_free_guidance": True,
            "cfg_scale": 2.0,
            "num_designs": 8,
            "step_scale": 1.5
        },
        {
            "name": "r7v2d_high_cfg",
            "description": "Higher CFG scale for stronger conditioning",
            "contig": "X,L,/0,100-130",
            "select_fixed_atoms": {
                "X": "all",
                "L": "all"
            },
            "select_buried": {
                "L": "all"
            },
            "select_hbond_acceptor": {
                "L": "O1,O2,O3,O4,O5,O6"
            },
            "select_hbond_donor": {
                "L": "O7"
            },
            "use_classifier_free_guidance": True,
            "cfg_scale": 3.0,  # Higher CFG
            "num_designs": 8,
            "step_scale": 1.5
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
        health = health_response.json()
        print(f"Health: {health.get('output', {}).get('result', {}).get('healthy', False)}")
    except Exception as e:
        print(f"API not available: {e}")
        return

    # Run each config
    results = {}

    for config in configs:
        name = config["name"]
        result = run_rfd3_ligand_scaffold(config, name)
        results[name] = result

        if "error" not in result and result.get("status") == "COMPLETED":
            save_results(result, name)

        time.sleep(2)

    # Summary
    print("\n" + "="*60)
    print("ROUND 7v2 SCAFFOLDING SUMMARY")
    print("="*60)

    for name, result in results.items():
        status = result.get("status", "UNKNOWN")
        if "error" in result:
            status = f"FAILED: {result['error'][:50]}"
        print(f"{name}: {status}")

    print(f"\nOutputs saved to: {OUTPUT_DIR}")

    # Save config for reference
    config_path = OUTPUT_DIR / "r7v2_configs_used.json"
    with open(config_path, 'w') as f:
        json.dump(configs, f, indent=2)
    print(f"Config saved: {config_path.name}")

if __name__ == "__main__":
    main()
