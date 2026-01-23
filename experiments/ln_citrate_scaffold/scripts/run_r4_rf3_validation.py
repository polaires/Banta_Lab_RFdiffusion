#!/usr/bin/env python3
"""
Round 4: RoseTTAFold3 structure prediction for LigandMPNN sequences.

Goal: Validate that designed sequences fold correctly and maintain metal coordination.
Input: Best sequences from R3c (LigandMPNN with no bias, metal context)
"""

import json
import sys
import time
from datetime import datetime
from pathlib import Path

import requests

API_URL = "http://localhost:8000/runsync"
EXPERIMENT_DIR = Path(__file__).parent.parent
MPNN_OUTPUT_DIR = EXPERIMENT_DIR / "outputs" / "round_03c_mpnn_natural"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs" / "round_04_rf3"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Select best sequences from each backbone (first 2 per backbone for speed)
SEQUENCES_PER_BACKBONE = 2


def load_fasta(fasta_path: Path) -> tuple:
    """Load sequence from FASTA file."""
    with open(fasta_path) as f:
        lines = f.readlines()
    header = lines[0].strip()[1:]  # Remove >
    sequence = "".join(line.strip() for line in lines[1:])
    return header, sequence


def run_rf3(sequence: str, name: str, ligand_smiles: str = None) -> dict:
    """
    Run RoseTTAFold3 structure prediction.

    Args:
        sequence: Protein sequence
        name: Design name
        ligand_smiles: Optional ligand SMILES (citrate)
    """
    # Citrate SMILES
    citrate_smiles = "C(C(=O)[O-])(CC(=O)[O-])(CC(=O)[O-])O"

    payload = {
        "input": {
            "task": "rf3",
            "sequence": sequence,
            # Include citrate as ligand context
            "ligand_smiles": citrate_smiles if ligand_smiles is None else ligand_smiles,
            # Include terbium as metal
            "metal": "TB",
        }
    }

    print(f"\n  Running RF3 on {name}...")
    print(f"  Sequence length: {len(sequence)}")

    try:
        response = requests.post(API_URL, json=payload, timeout=600)
        response.raise_for_status()
        result = response.json()
        return result
    except requests.exceptions.Timeout:
        return {"error": "Timeout after 600s"}
    except requests.exceptions.RequestException as e:
        return {"error": str(e)}


def save_prediction(result: dict, name: str) -> dict:
    """Save RF3 prediction results."""
    saved = {"name": name}

    if "error" in result:
        print(f"  ERROR: {result['error']}")
        saved["error"] = result["error"]
        return saved

    # Extract prediction from response
    output = result.get("output", {})
    result_data = output.get("result", {})

    # Try different paths for PDB content
    pdb_content = None
    for key in ["pdb_content", "pdb", "structure", "content"]:
        if key in result_data:
            pdb_content = result_data[key]
            break

    if not pdb_content:
        # Check if it's in a predictions list
        predictions = result_data.get("predictions", [])
        if predictions and isinstance(predictions[0], dict):
            pdb_content = predictions[0].get("pdb_content") or predictions[0].get("content")

    if not pdb_content:
        print(f"  No structure in response")
        print(f"  Result keys: {list(result_data.keys()) if result_data else 'None'}")
        saved["error"] = "No structure returned"
        return saved

    # Save PDB
    pdb_file = OUTPUT_DIR / f"{name}_rf3.pdb"
    with open(pdb_file, "w") as f:
        f.write(pdb_content)

    saved["pdb"] = pdb_file.name

    # Extract confidence metrics if available
    for metric in ["plddt", "ptm", "iptm", "confidence"]:
        if metric in result_data:
            saved[metric] = result_data[metric]

    print(f"  Saved: {pdb_file.name}")
    if "plddt" in saved:
        print(f"  pLDDT: {saved['plddt']:.2f}")

    return saved


def main():
    print("="*60)
    print("Ln-Citrate Scaffold - Round 4: RF3 Structure Prediction")
    print("="*60)
    print(f"\nTimestamp: {datetime.now().isoformat()}")
    print(f"Input dir: {MPNN_OUTPUT_DIR}")
    print(f"Output dir: {OUTPUT_DIR}")
    print(f"\nStrategy:")
    print(f"  - Predict structures for top {SEQUENCES_PER_BACKBONE} sequences per backbone")
    print(f"  - Include citrate ligand context")
    print(f"  - Validate metal coordination geometry")

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

    # Find all FASTA files from R3c
    fasta_files = sorted(MPNN_OUTPUT_DIR.glob("*_seq*.fasta"))
    if not fasta_files:
        print(f"ERROR: No FASTA files found in {MPNN_OUTPUT_DIR}")
        sys.exit(1)

    print(f"\nFound {len(fasta_files)} sequences")

    # Group by backbone
    backbones = {}
    for fasta in fasta_files:
        # Extract backbone name (everything before _seq)
        name = fasta.stem
        backbone = name.rsplit("_seq", 1)[0]
        if backbone not in backbones:
            backbones[backbone] = []
        backbones[backbone].append(fasta)

    print(f"Backbones: {list(backbones.keys())}")

    # Process top sequences from each backbone
    all_results = {}
    total_predictions = 0

    for backbone, fastas in backbones.items():
        print(f"\n{'='*40}")
        print(f"Backbone: {backbone}")

        # Take first N sequences
        for fasta in fastas[:SEQUENCES_PER_BACKBONE]:
            header, sequence = load_fasta(fasta)
            name = fasta.stem

            start_time = time.time()
            result = run_rf3(sequence, name)
            elapsed = time.time() - start_time

            saved = save_prediction(result, name)
            saved["elapsed_seconds"] = round(elapsed, 1)
            saved["sequence"] = sequence

            all_results[name] = saved

            if "pdb" in saved:
                total_predictions += 1

            print(f"  Completed in {elapsed:.1f}s")

    # Save summary
    summary = {
        "round": 4,
        "task": "RF3 structure prediction",
        "timestamp": datetime.now().isoformat(),
        "total_predictions": total_predictions,
        "sequences_per_backbone": SEQUENCES_PER_BACKBONE,
        "results": all_results
    }

    summary_path = OUTPUT_DIR / "run_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n{'='*60}")
    print("Round 4 Complete")
    print(f"{'='*60}")
    print(f"Total predictions: {total_predictions}")
    print(f"Summary: {summary_path}")
    print(f"\nNext: Analyze predicted structures for metal coordination")


if __name__ == "__main__":
    main()
