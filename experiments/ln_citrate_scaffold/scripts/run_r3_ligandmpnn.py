#!/usr/bin/env python3
"""
Round 3: LigandMPNN sequence design on best RFD3 backbones.

Goal: Design sequences with strong Glu/Asp bias for lanthanide coordination.
Input: Best backbones from R2/R2b with CN=3-5

Key parameters:
- model_type: ligand_mpnn (metal-aware)
- bias_AA: "D:5.0,E:4.0,N:2.0,Q:2.0" (favor hard O donors)
- omit_AA: "C" (exclude soft S donor - HSAB principle)
- temperature: 0.1 (default, low diversity for high confidence)
- ligand_mpnn_use_atom_context: 1 (use metal/ligand context)
"""

import json
import sys
import time
from datetime import datetime
from pathlib import Path

import requests

API_URL = "http://localhost:8000/runsync"
EXPERIMENT_DIR = Path(__file__).parent.parent
OUTPUT_DIR = EXPERIMENT_DIR / "outputs" / "round_03c_mpnn_natural"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Best backbones from R2/R2b analysis
BEST_BACKBONES = [
    # R2b best - pure carboxylate coordination
    ("round_02b", "r2b_b_cfg_hbond_007"),  # CN=3, GLU+GLU+ASP - best quality
    ("round_02b", "r2b_b_cfg_hbond_005"),  # CN=5, 2 ASP + HIS
    ("round_02b", "r2b_b_cfg_hbond_012"),  # CN=3, 2 GLU + ASN
    # R2 best - for comparison
    ("round_02", "r2a_simple_buried_metal_007"),  # CN=5, 3 GLU - excellent
    ("round_02", "r2b_hbond_citrate_005"),  # CN=5, mixed O/N
]


def load_pdb(pdb_path: Path) -> str:
    """Load PDB file content."""
    with open(pdb_path) as f:
        return f.read()


def run_ligandmpnn(pdb_content: str, name: str, num_seqs: int = 8) -> dict:
    """
    Run LigandMPNN sequence design via API.

    Uses lanthanide-optimized biasing:
    - Strong bias for Asp (D) and Glu (E) - hard oxygen donors
    - Moderate bias for Asn (N) and Gln (Q) - additional O donors
    - Exclude Cys (C) - soft sulfur donor, bad for Ln³⁺
    """

    payload = {
        "input": {
            "task": "mpnn",  # API uses "mpnn" not "ligandmpnn"
            "pdb_content": pdb_content,
            "num_seqs": num_seqs,
            "model_type": "ligand_mpnn",
            # NO global bias - let LigandMPNN use metal context naturally
            # Just exclude cysteine (soft donor, bad for Ln)
            "omit_AA": "C",
            # Use ligand/metal context - this should naturally favor O-donors near metal
            "ligand_mpnn_use_atom_context": 1,
            "ligand_mpnn_use_side_chain_context": 1,
            # Sampling
            "temperature": 0.1,
            # Side chain packing
            "pack_side_chains": 1,
        }
    }

    print(f"\n  Running LigandMPNN on {name}...")
    print(f"  No global bias - using metal context | Omit: C")
    print(f"  Sequences requested: {num_seqs}")

    try:
        response = requests.post(API_URL, json=payload, timeout=300)
        response.raise_for_status()
        result = response.json()
        return result
    except requests.exceptions.Timeout:
        return {"error": "Timeout after 300s"}
    except requests.exceptions.RequestException as e:
        return {"error": str(e)}


def parse_fasta_content(fasta_content: str) -> list:
    """Parse FASTA format content into list of (name, sequence) tuples."""
    sequences = []
    current_name = None
    current_seq = []

    for line in fasta_content.strip().split('\n'):
        if line.startswith('>'):
            if current_name is not None:
                sequences.append((current_name, ''.join(current_seq)))
            current_name = line[1:].strip()
            current_seq = []
        else:
            current_seq.append(line.strip())

    if current_name is not None:
        sequences.append((current_name, ''.join(current_seq)))

    return sequences


def save_sequences(result: dict, backbone_name: str) -> list:
    """Save designed sequences to FASTA files."""
    saved = []

    if "error" in result:
        print(f"  ERROR: {result['error']}")
        return saved

    # Extract sequences from response - API returns FASTA in result.sequences[0].content
    output = result.get("output", {})
    result_data = output.get("result", {})
    sequences_list = result_data.get("sequences", [])

    if not sequences_list:
        print(f"  No sequences in response.")
        print(f"  Response keys: {list(result.keys())}")
        if "output" in result:
            print(f"  output keys: {list(output.keys())}")
            if "result" in output:
                print(f"  result keys: {list(result_data.keys())}")
        return saved

    # API returns [{'filename': 'sequences.fasta', 'content': '>design_1\nSEQ\n...'}]
    fasta_content = None
    for item in sequences_list:
        if isinstance(item, dict) and 'content' in item:
            fasta_content = item['content']
            break

    if not fasta_content:
        print(f"  Could not find FASTA content in response")
        return saved

    # Parse FASTA content
    parsed_seqs = parse_fasta_content(fasta_content)

    # Save combined FASTA file
    combined_fasta = OUTPUT_DIR / f"{backbone_name}_all_seqs.fasta"
    with open(combined_fasta, "w") as f:
        for name, seq in parsed_seqs:
            f.write(f">{backbone_name}_{name}\n{seq}\n")

    # Save individual sequences and analyze
    for i, (name, sequence) in enumerate(parsed_seqs):
        # Count coordinating residues
        counts = count_coordinating_residues(sequence)

        # Save individual FASTA
        fasta_file = OUTPUT_DIR / f"{backbone_name}_seq{i:02d}.fasta"
        with open(fasta_file, "w") as f:
            f.write(f">{backbone_name}_seq{i:02d} D={counts['D']} E={counts['E']} total_O={counts['total_O_donors']}\n")
            f.write(f"{sequence}\n")

        saved.append({
            "fasta": fasta_file.name,
            "sequence": sequence,
            "D_count": counts["D"],
            "E_count": counts["E"],
            "total_O_donors": counts["total_O_donors"],
            "C_count": counts["C"]
        })

        print(f"  Saved: {fasta_file.name} (D={counts['D']}, E={counts['E']}, C={counts['C']})")

    return saved


def count_coordinating_residues(sequence: str) -> dict:
    """Count potential metal-coordinating residues in sequence."""
    counts = {
        "D": sequence.count("D"),  # Asp - carboxylate O
        "E": sequence.count("E"),  # Glu - carboxylate O
        "N": sequence.count("N"),  # Asn - amide O
        "Q": sequence.count("Q"),  # Gln - amide O
        "H": sequence.count("H"),  # His - imidazole N
        "S": sequence.count("S"),  # Ser - hydroxyl O
        "T": sequence.count("T"),  # Thr - hydroxyl O
        "C": sequence.count("C"),  # Cys - thiol S (bad for Ln)
    }
    counts["total_O_donors"] = counts["D"] + counts["E"] + counts["N"] + counts["Q"]
    counts["total_hard"] = counts["total_O_donors"] + counts["S"] + counts["T"]
    return counts


def main():
    print("="*60)
    print("Ln-Citrate Scaffold - Round 3: LigandMPNN Sequence Design")
    print("="*60)
    print(f"\nTimestamp: {datetime.now().isoformat()}")
    print(f"Output dir: {OUTPUT_DIR}")
    print("\nStrategy:")
    print("  - Use metal-aware LigandMPNN model")
    print("  - Bias toward Asp/Glu (hard O donors for Ln3+)")
    print("  - Exclude Cys (soft S donor - HSAB violation)")
    print("  - 8 sequences per backbone")

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

    # Process each backbone
    all_results = {}
    total_seqs = 0
    num_seqs = 8

    for round_dir, backbone_name in BEST_BACKBONES:
        pdb_path = EXPERIMENT_DIR / "outputs" / round_dir / f"{backbone_name}.pdb"

        if not pdb_path.exists():
            print(f"\n  WARNING: {pdb_path} not found, skipping")
            continue

        print(f"\n{'='*40}")
        print(f"Processing: {backbone_name}")
        print(f"Source: {round_dir}")

        pdb_content = load_pdb(pdb_path)

        start_time = time.time()
        result = run_ligandmpnn(pdb_content, backbone_name, num_seqs=num_seqs)
        elapsed = time.time() - start_time

        saved = save_sequences(result, backbone_name)
        total_seqs += len(saved)

        all_results[backbone_name] = {
            "source": round_dir,
            "n_sequences": len(saved),
            "elapsed_seconds": round(elapsed, 1),
            "sequences": saved
        }

        print(f"  Generated: {len(saved)} sequences in {elapsed:.1f}s")

    # Save summary
    summary = {
        "round": 3,
        "task": "LigandMPNN sequence design",
        "timestamp": datetime.now().isoformat(),
        "bias_AA": "D:5.0,E:4.0,N:2.0,Q:2.0",
        "omit_AA": "C",
        "temperature": 0.1,
        "total_sequences": total_seqs,
        "backbones": all_results
    }

    summary_path = OUTPUT_DIR / "run_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n{'='*60}")
    print("Round 3 Complete")
    print(f"{'='*60}")
    print(f"Total sequences: {total_seqs}")
    print(f"Summary: {summary_path}")
    print(f"\nNext: Analyze sequences and run AF2/AF3 validation")


if __name__ == "__main__":
    main()
