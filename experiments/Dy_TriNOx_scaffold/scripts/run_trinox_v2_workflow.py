"""
Run LigandMPNN + RF3 validation on v2 Dy-TriNOx scaffolds.
"""

import requests
import json
from pathlib import Path
from datetime import datetime
import math

API_URL = "http://localhost:8000/runsync"
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR.parent / "outputs_v2"


def run_ligandmpnn(pdb_content: str, num_seqs: int = 8, temperature: float = 0.2) -> list:
    """Run LigandMPNN on a scaffold."""
    payload = {
        "input": {
            "task": "mpnn",
            "pdb_content": pdb_content,
            "num_seqs": num_seqs,
            "ligand_mpnn_use_atom_context": 1,
            "temperature": temperature,
        }
    }

    try:
        response = requests.post(API_URL, json=payload, timeout=300)
        result = response.json()

        output = result.get("output", {})
        res = output.get("result", {})
        sequences = res.get("sequences", [])

        extracted = []
        for seq_data in sequences:
            if isinstance(seq_data, dict) and "content" in seq_data:
                content = seq_data["content"]
                lines = content.strip().split("\n")
                for j in range(0, len(lines), 2):
                    if j + 1 < len(lines) and lines[j].startswith(">"):
                        seq = lines[j + 1]
                        extracted.append(seq)
        return extracted
    except Exception as e:
        print(f"    MPNN Error: {e}")
        return []


def run_rf3_validation(sequence: str, name: str) -> dict:
    """Run RF3 structure prediction."""
    payload = {
        "input": {
            "task": "rf3",
            "sequence": sequence,
            "name": name
        }
    }

    try:
        response = requests.post(API_URL, json=payload, timeout=600)
        result = response.json()

        output = result.get("output", {})
        res = output.get("result", {})
        confidences = res.get("confidences", {})

        return {
            "plddt": confidences.get("mean_plddt", 0),
            "ptm": confidences.get("ptm", 0),
            "pae": confidences.get("overall_pae", 999),
        }
    except Exception as e:
        print(f"    RF3 Error: {e}")
        return {"plddt": 0, "ptm": 0, "pae": 999}


def main():
    print("=" * 70)
    print("Dy-TriNOx v2 WORKFLOW")
    print("LigandMPNN + RF3 Validation")
    print("=" * 70)

    # Find all v2 design PDBs
    pdb_files = sorted([f for f in OUTPUT_DIR.glob("v2_*.pdb") if not f.name.startswith("test_")])
    print(f"\nFound {len(pdb_files)} v2 scaffolds")

    if not pdb_files:
        print("ERROR: No v2 scaffold PDBs found!")
        return

    all_sequences = []

    # Step 1: LigandMPNN
    print("\n" + "=" * 70)
    print("STEP 1: LIGANDMPNN SEQUENCE DESIGN (8 seqs per backbone)")
    print("=" * 70)

    for pdb_file in pdb_files:
        print(f"\n{pdb_file.stem}...")

        with open(pdb_file, 'r') as f:
            pdb_content = f.read()

        sequences = run_ligandmpnn(pdb_content, num_seqs=8, temperature=0.2)

        if sequences:
            for i, seq in enumerate(sequences):
                name = f"{pdb_file.stem}_seq{i+1:02d}"
                all_sequences.append({
                    "name": name,
                    "seq": seq,
                    "backbone": pdb_file.stem,
                    "length": len(seq)
                })
                print(f"  {name}: {len(seq)} aa")
        else:
            print(f"  No sequences generated")

    # Save all sequences
    if all_sequences:
        fasta_out = OUTPUT_DIR / "trinox_v2_sequences.fasta"
        with open(fasta_out, 'w') as f:
            for s in all_sequences:
                f.write(f">{s['name']}\n{s['seq']}\n")
        print(f"\nSaved {len(all_sequences)} sequences to {fasta_out}")

    # Step 2: RF3 Validation
    print("\n" + "=" * 70)
    print("STEP 2: RF3 STRUCTURE VALIDATION")
    print("=" * 70)

    results = []

    for i, s in enumerate(all_sequences):
        print(f"\n[{i+1}/{len(all_sequences)}] Validating {s['name']}...")

        conf = run_rf3_validation(s["seq"], s["name"])
        s.update(conf)
        results.append(s)

        print(f"  pLDDT: {conf['plddt']:.4f}, pTM: {conf['ptm']:.4f}, PAE: {conf['pae']:.2f}")

    # Results Summary
    print("\n" + "=" * 70)
    print("RESULTS SUMMARY - Sorted by pLDDT")
    print("=" * 70)

    sorted_results = sorted(results, key=lambda x: x['plddt'], reverse=True)

    print(f"\n{'Design':<30} {'Backbone':<18} {'Len':>5} {'pLDDT':>8} {'pTM':>8} {'PAE':>8}")
    print("-" * 85)
    for r in sorted_results[:20]:
        print(f"{r['name']:<30} {r['backbone']:<18} {r['length']:>5} {r['plddt']:>8.4f} {r['ptm']:>8.4f} {r['pae']:>8.2f}")

    # Quality thresholds
    high_quality = [r for r in results if r['plddt'] >= 0.85 and r['ptm'] >= 0.80]
    medium_quality = [r for r in results if r['plddt'] >= 0.75 and r['ptm'] >= 0.70]

    print(f"\n" + "=" * 70)
    print("QUALITY SUMMARY")
    print("=" * 70)
    print(f"Total sequences: {len(results)}")
    print(f"High quality (pLDDT>=0.85, pTM>=0.80): {len(high_quality)}")
    print(f"Medium quality (pLDDT>=0.75, pTM>=0.70): {len(medium_quality)}")

    if sorted_results:
        best = sorted_results[0]
        print(f"\nBEST DESIGN: {best['name']}")
        print(f"  Backbone: {best['backbone']}")
        print(f"  Length: {best['length']} aa")
        print(f"  pLDDT: {best['plddt']:.4f}")
        print(f"  pTM: {best['ptm']:.4f}")
        print(f"  PAE: {best['pae']:.2f}")

    # Save results
    results_file = OUTPUT_DIR / "trinox_v2_workflow_results.json"
    with open(results_file, 'w') as f:
        json.dump({
            "timestamp": datetime.now().isoformat(),
            "total_backbones": len(pdb_files),
            "total_sequences": len(results),
            "high_quality_count": len(high_quality),
            "medium_quality_count": len(medium_quality),
            "results": sorted_results,
        }, f, indent=2)
    print(f"\nFull results saved to {results_file}")


if __name__ == "__main__":
    main()
