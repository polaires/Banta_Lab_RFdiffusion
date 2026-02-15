#!/usr/bin/env python3
"""
Run LigandMPNN on best Round 6 designs with strong Glu/Asp bias.

Top designs to process:
1. r6_c_dimer_symmetric_025 (CN=6, 3 carboxylates) - BEST
2. r6_c_dimer_symmetric_012 (CN=5, 4 carboxylates)
3. r6_b_partial_aggressive_031 (CN=6, 1 carboxylate)
4. r6_c_dimer_symmetric_010 (CN=4, 2 carboxylates)

Parameters:
- model_type: ligand_mpnn
- temperature: 0.1 (conservative)
- bias_AA: E:3.0,D:3.0,N:2.0,Q:2.0,C:-5.0 (prefer carboxylates, avoid Cys)
- num_seqs: 8 per design
"""

import json
import time
from datetime import datetime
from pathlib import Path

import requests

API_URL = "http://localhost:8000/runsync"
EXPERIMENT_DIR = Path(__file__).parent.parent
INPUT_DIR = EXPERIMENT_DIR / "outputs" / "round_06"
OUTPUT_DIR = EXPERIMENT_DIR / "outputs" / "round_06_mpnn"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Top designs to process
TOP_DESIGNS = [
    "r6_c_dimer_symmetric_025",  # Best: CN=6, 3 carboxylates
    "r6_c_dimer_symmetric_012",  # CN=5, 4 carboxylates
    "r6_b_partial_aggressive_031",  # CN=6, 1 carboxylate
    "r6_c_dimer_symmetric_010",  # CN=4, 2 carboxylates
]


def load_pdb(pdb_path: Path) -> str:
    with open(pdb_path) as f:
        return f.read()


def run_ligandmpnn(pdb_content: str, name: str, num_seqs: int = 8) -> dict:
    """Run LigandMPNN with carboxylate bias."""

    payload = {
        "input": {
            "task": "mpnn",
            "pdb_content": pdb_content,
            "model_type": "ligand_mpnn",
            "temperature": 0.1,
            "num_seqs": num_seqs,
            # Strong bias for carboxylates (hard acid donors for Ln3+)
            # Penalize Cys (soft sulfur donor)
            "bias_AA": "E:3.0,D:3.0,N:2.0,Q:2.0,C:-5.0",
            # Use atom and side chain context for metal coordination
            "ligand_mpnn_use_atom_context": 1,
            "ligand_mpnn_use_side_chain_context": 1,
        }
    }

    print(f"\n  Running LigandMPNN for {name}...")
    print(f"  Params: temperature=0.1, num_seqs={num_seqs}")
    print(f"  Bias: E:3.0, D:3.0, N:2.0, Q:2.0, C:-5.0")

    try:
        response = requests.post(API_URL, json=payload, timeout=600)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.Timeout:
        return {"error": "Timeout after 600s"}
    except requests.exceptions.RequestException as e:
        return {"error": str(e)}


def parse_fasta_content(fasta_content: str) -> list:
    """Parse FASTA content string into list of (header, sequence) tuples."""
    sequences = []
    current_header = None
    current_seq = []

    for line in fasta_content.strip().split('\n'):
        line = line.strip()
        if line.startswith('>'):
            if current_header is not None:
                sequences.append((current_header, ''.join(current_seq)))
            current_header = line[1:]  # Remove '>'
            current_seq = []
        else:
            current_seq.append(line)

    if current_header is not None:
        sequences.append((current_header, ''.join(current_seq)))

    return sequences


def save_sequences(result: dict, design_name: str) -> list:
    """Save sequences to FASTA files."""
    saved = []

    if "error" in result:
        print(f"  ERROR: {result['error']}")
        return saved

    output = result.get("output", {})
    result_data = output.get("result", {})
    sequences_data = result_data.get("sequences", [])

    if not sequences_data:
        sequences_data = output.get("sequences", [])
    if not sequences_data:
        sequences_data = result.get("sequences", [])

    if not sequences_data:
        print(f"  No sequences in response.")
        return saved

    # API returns sequences as FASTA file objects with 'content' field
    all_seqs = []
    for seq_obj in sequences_data:
        if isinstance(seq_obj, dict) and "content" in seq_obj:
            # Parse FASTA content
            parsed = parse_fasta_content(seq_obj["content"])
            all_seqs.extend(parsed)
        elif isinstance(seq_obj, dict):
            # Try other formats
            seq = seq_obj.get("sequence") or seq_obj.get("seq")
            if seq:
                all_seqs.append((f"seq", seq))
        elif isinstance(seq_obj, str):
            all_seqs.append((f"seq", seq_obj))

    if not all_seqs:
        print(f"  No sequences parsed from response.")
        return saved

    # Save all sequences to a single FASTA with renamed headers
    all_fasta_path = OUTPUT_DIR / f"{design_name}_all_seqs.fasta"
    with open(all_fasta_path, "w") as f:
        for i, (header, seq) in enumerate(all_seqs):
            f.write(f">{design_name}_seq{i:02d}\n")
            f.write(f"{seq}\n")
            saved.append(f"{design_name}_seq{i:02d}")

    print(f"  Saved {len(saved)} sequences to {all_fasta_path.name}")

    # Also save individual FASTA files
    for i, (header, seq) in enumerate(all_seqs):
        fasta_path = OUTPUT_DIR / f"{design_name}_seq{i:02d}.fasta"
        with open(fasta_path, "w") as f:
            f.write(f">{design_name}_seq{i:02d}\n")
            f.write(f"{seq}\n")

    return saved


def main():
    print("="*60)
    print("Round 6 - LigandMPNN Sequence Design")
    print("="*60)
    print(f"\nTimestamp: {datetime.now().isoformat()}")
    print(f"Output dir: {OUTPUT_DIR}")
    print(f"\nDesigns to process: {len(TOP_DESIGNS)}")
    print("Bias: Strong Glu/Asp preference, Cys penalty")

    results = {}
    for design_name in TOP_DESIGNS:
        pdb_path = INPUT_DIR / f"{design_name}.pdb"

        if not pdb_path.exists():
            print(f"\n  ERROR: {pdb_path} not found")
            continue

        print(f"\n{'='*60}")
        print(f"Processing: {design_name}")
        print(f"{'='*60}")

        pdb_content = load_pdb(pdb_path)
        print(f"  Input: {pdb_path.name} ({len(pdb_content)} bytes)")

        start_time = time.time()
        result = run_ligandmpnn(pdb_content, design_name, num_seqs=8)
        elapsed = time.time() - start_time

        saved = save_sequences(result, design_name)
        results[design_name] = {
            "sequences_saved": len(saved),
            "elapsed_seconds": round(elapsed, 1)
        }

        print(f"\n  Completed: {len(saved)} sequences in {elapsed:.1f}s")

    # Save summary
    summary_path = OUTPUT_DIR / "mpnn_summary.json"
    with open(summary_path, "w") as f:
        json.dump({
            "timestamp": datetime.now().isoformat(),
            "designs": results,
            "parameters": {
                "model_type": "ligand_mpnn",
                "temperature": 0.1,
                "bias_AA": "E:3.0,D:3.0,N:2.0,Q:2.0,C:-5.0",
                "num_seqs": 8
            }
        }, f, indent=2)

    print(f"\n{'='*60}")
    print("LIGANDMPNN COMPLETE")
    print(f"{'='*60}")
    for name, res in results.items():
        print(f"  {name}: {res['sequences_saved']} sequences")
    print(f"\nSummary saved: {summary_path}")


if __name__ == "__main__":
    main()
