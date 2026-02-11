"""
Resume RF3 validation from sequence 47 onward.
"""

import requests
import json
from pathlib import Path

API_URL = "http://localhost:8000/runsync"
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR.parent / "outputs"

# Results from sequences 1-46 (already completed)
COMPLETED_RESULTS = [
    {"name": "r1a_small_000_seq01", "backbone": "r1a_small_000", "length": 63, "plddt": 0.8776, "ptm": 0.8595, "pae": 2.46},
    {"name": "r1a_small_000_seq02", "backbone": "r1a_small_000", "length": 63, "plddt": 0.8751, "ptm": 0.8593, "pae": 2.51},
    {"name": "r1a_small_000_seq03", "backbone": "r1a_small_000", "length": 63, "plddt": 0.8778, "ptm": 0.8728, "pae": 2.02},
    {"name": "r1a_small_000_seq04", "backbone": "r1a_small_000", "length": 63, "plddt": 0.8897, "ptm": 0.8733, "pae": 1.96},
    {"name": "r1a_small_000_seq05", "backbone": "r1a_small_000", "length": 63, "plddt": 0.8454, "ptm": 0.8318, "pae": 2.84},
    {"name": "r1a_small_000_seq06", "backbone": "r1a_small_000", "length": 63, "plddt": 0.8540, "ptm": 0.8362, "pae": 2.99},
    {"name": "r1a_small_000_seq07", "backbone": "r1a_small_000", "length": 63, "plddt": 0.8850, "ptm": 0.8379, "pae": 2.69},
    {"name": "r1a_small_000_seq08", "backbone": "r1a_small_000", "length": 63, "plddt": 0.8707, "ptm": 0.8422, "pae": 2.47},
    {"name": "r1a_small_001_seq01", "backbone": "r1a_small_001", "length": 63, "plddt": 0.8860, "ptm": 0.7958, "pae": 3.10},
    {"name": "r1a_small_001_seq02", "backbone": "r1a_small_001", "length": 63, "plddt": 0.8710, "ptm": 0.7698, "pae": 3.33},
    {"name": "r1a_small_001_seq03", "backbone": "r1a_small_001", "length": 63, "plddt": 0.9184, "ptm": 0.8742, "pae": 1.84},
    {"name": "r1a_small_001_seq04", "backbone": "r1a_small_001", "length": 63, "plddt": 0.8968, "ptm": 0.7733, "pae": 3.12},
    {"name": "r1a_small_001_seq05", "backbone": "r1a_small_001", "length": 63, "plddt": 0.9072, "ptm": 0.8522, "pae": 2.17},
    {"name": "r1a_small_001_seq06", "backbone": "r1a_small_001", "length": 63, "plddt": 0.9108, "ptm": 0.8553, "pae": 2.18},
    {"name": "r1a_small_001_seq07", "backbone": "r1a_small_001", "length": 63, "plddt": 0.9050, "ptm": 0.8248, "pae": 2.38},
    {"name": "r1a_small_001_seq08", "backbone": "r1a_small_001", "length": 63, "plddt": 0.8976, "ptm": 0.8202, "pae": 2.62},
    {"name": "r1a_small_002_seq01", "backbone": "r1a_small_002", "length": 63, "plddt": 0.8364, "ptm": 0.6627, "pae": 6.36},
    {"name": "r1a_small_002_seq02", "backbone": "r1a_small_002", "length": 63, "plddt": 0.8487, "ptm": 0.7016, "pae": 5.43},
    {"name": "r1a_small_002_seq03", "backbone": "r1a_small_002", "length": 63, "plddt": 0.8358, "ptm": 0.6669, "pae": 6.02},
    {"name": "r1a_small_002_seq04", "backbone": "r1a_small_002", "length": 63, "plddt": 0.8708, "ptm": 0.7325, "pae": 4.79},
    {"name": "r1a_small_002_seq05", "backbone": "r1a_small_002", "length": 63, "plddt": 0.8649, "ptm": 0.7392, "pae": 4.40},
    {"name": "r1a_small_002_seq06", "backbone": "r1a_small_002", "length": 63, "plddt": 0.8532, "ptm": 0.6213, "pae": 6.74},
    {"name": "r1a_small_002_seq07", "backbone": "r1a_small_002", "length": 63, "plddt": 0.8639, "ptm": 0.7238, "pae": 4.87},
    {"name": "r1a_small_002_seq08", "backbone": "r1a_small_002", "length": 63, "plddt": 0.8826, "ptm": 0.8110, "pae": 3.08},
    {"name": "r1a_small_003_seq01", "backbone": "r1a_small_003", "length": 63, "plddt": 0.8113, "ptm": 0.5274, "pae": 8.24},
    {"name": "r1a_small_003_seq02", "backbone": "r1a_small_003", "length": 63, "plddt": 0.8208, "ptm": 0.5354, "pae": 7.92},
    {"name": "r1a_small_003_seq03", "backbone": "r1a_small_003", "length": 63, "plddt": 0.8180, "ptm": 0.5012, "pae": 9.27},
    {"name": "r1a_small_003_seq04", "backbone": "r1a_small_003", "length": 63, "plddt": 0.8092, "ptm": 0.5309, "pae": 7.91},
    {"name": "r1a_small_003_seq05", "backbone": "r1a_small_003", "length": 63, "plddt": 0.8099, "ptm": 0.5236, "pae": 8.03},
    {"name": "r1a_small_003_seq06", "backbone": "r1a_small_003", "length": 63, "plddt": 0.8039, "ptm": 0.5177, "pae": 8.54},
    {"name": "r1a_small_003_seq07", "backbone": "r1a_small_003", "length": 63, "plddt": 0.8150, "ptm": 0.5124, "pae": 8.87},
    {"name": "r1a_small_003_seq08", "backbone": "r1a_small_003", "length": 63, "plddt": 0.8006, "ptm": 0.4874, "pae": 9.52},
    {"name": "r1a_small_004_seq01", "backbone": "r1a_small_004", "length": 63, "plddt": 0.8087, "ptm": 0.6016, "pae": 6.56},
    {"name": "r1a_small_004_seq02", "backbone": "r1a_small_004", "length": 63, "plddt": 0.8221, "ptm": 0.6402, "pae": 5.94},
    {"name": "r1a_small_004_seq03", "backbone": "r1a_small_004", "length": 63, "plddt": 0.7929, "ptm": 0.5829, "pae": 6.85},
    {"name": "r1a_small_004_seq04", "backbone": "r1a_small_004", "length": 63, "plddt": 0.8226, "ptm": 0.6217, "pae": 5.96},
    {"name": "r1a_small_004_seq05", "backbone": "r1a_small_004", "length": 63, "plddt": 0.8256, "ptm": 0.6593, "pae": 5.39},
    {"name": "r1a_small_004_seq06", "backbone": "r1a_small_004", "length": 63, "plddt": 0.8264, "ptm": 0.6904, "pae": 4.66},
    {"name": "r1a_small_004_seq07", "backbone": "r1a_small_004", "length": 63, "plddt": 0.8034, "ptm": 0.6798, "pae": 5.46},
    {"name": "r1a_small_004_seq08", "backbone": "r1a_small_004", "length": 63, "plddt": 0.8045, "ptm": 0.6233, "pae": 5.70},
    {"name": "r1b_medium_000_seq01", "backbone": "r1b_medium_000", "length": 108, "plddt": 0.9117, "ptm": 0.9164, "pae": 2.33},
    {"name": "r1b_medium_000_seq02", "backbone": "r1b_medium_000", "length": 108, "plddt": 0.8995, "ptm": 0.9146, "pae": 2.53},
    {"name": "r1b_medium_000_seq03", "backbone": "r1b_medium_000", "length": 108, "plddt": 0.8922, "ptm": 0.9121, "pae": 2.47},
    {"name": "r1b_medium_000_seq04", "backbone": "r1b_medium_000", "length": 108, "plddt": 0.9223, "ptm": 0.9303, "pae": 2.04},
    {"name": "r1b_medium_000_seq05", "backbone": "r1b_medium_000", "length": 108, "plddt": 0.8797, "ptm": 0.9153, "pae": 2.43},
    {"name": "r1b_medium_000_seq06", "backbone": "r1b_medium_000", "length": 108, "plddt": 0.8848, "ptm": 0.9185, "pae": 2.48},
]


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
    from datetime import datetime

    print("=" * 70)
    print("RESUMING RF3 VALIDATION FROM SEQUENCE 47")
    print("=" * 70)

    # Read all sequences from FASTA
    fasta_file = OUTPUT_DIR / "trinox_all_sequences.fasta"
    all_sequences = []

    with open(fasta_file, 'r') as f:
        lines = f.readlines()
        for i in range(0, len(lines), 2):
            if i + 1 < len(lines):
                name = lines[i].strip()[1:]  # Remove '>'
                seq = lines[i + 1].strip()
                backbone = "_".join(name.split("_")[:-1])
                all_sequences.append({
                    "name": name,
                    "seq": seq,
                    "backbone": backbone,
                    "length": len(seq)
                })

    print(f"Loaded {len(all_sequences)} sequences from FASTA")

    # Start from sequence 47 (index 46)
    start_idx = 46
    remaining_sequences = all_sequences[start_idx:]

    print(f"Already completed: {start_idx} sequences")
    print(f"Remaining to validate: {len(remaining_sequences)} sequences")

    # Use previously completed results
    results = list(COMPLETED_RESULTS)

    # Validate remaining sequences
    for i, s in enumerate(remaining_sequences):
        total_idx = start_idx + i + 1
        print(f"\n[{total_idx}/104] Validating {s['name']}...")

        conf = run_rf3_validation(s["seq"], s["name"])

        result_entry = {
            "name": s["name"],
            "backbone": s["backbone"],
            "length": s["length"],
            "plddt": conf["plddt"],
            "ptm": conf["ptm"],
            "pae": conf["pae"]
        }
        results.append(result_entry)

        print(f"  pLDDT: {conf['plddt']:.4f}, pTM: {conf['ptm']:.4f}, PAE: {conf['pae']:.2f}")

    # ================================================================
    # RESULTS SUMMARY
    # ================================================================
    print("\n" + "=" * 70)
    print("RESULTS SUMMARY - Sorted by pLDDT")
    print("=" * 70)

    sorted_results = sorted(results, key=lambda x: x['plddt'], reverse=True)

    print(f"\n{'Design':<30} {'Backbone':<18} {'Len':>5} {'pLDDT':>8} {'pTM':>8} {'PAE':>8}")
    print("-" * 85)
    for r in sorted_results[:20]:  # Show top 20
        print(f"{r['name']:<30} {r['backbone']:<18} {r['length']:>5} {r['plddt']:>8.4f} {r['ptm']:>8.4f} {r['pae']:>8.2f}")

    # Quality thresholds
    high_quality = [r for r in results if r['plddt'] >= 0.85 and r['ptm'] >= 0.80]
    medium_quality = [r for r in results if r['plddt'] >= 0.75 and r['ptm'] >= 0.70]

    print(f"\n" + "=" * 70)
    print("QUALITY SUMMARY")
    print("=" * 70)
    print(f"Total sequences: {len(results)}")
    print(f"High quality (pLDDT≥0.85, pTM≥0.80): {len(high_quality)}")
    print(f"Medium quality (pLDDT≥0.75, pTM≥0.70): {len(medium_quality)}")

    if sorted_results:
        best = sorted_results[0]
        print(f"\nBEST DESIGN: {best['name']}")
        print(f"  Backbone: {best['backbone']}")
        print(f"  Length: {best['length']} aa")
        print(f"  pLDDT: {best['plddt']:.4f}")
        print(f"  pTM: {best['ptm']:.4f}")
        print(f"  PAE: {best['pae']:.2f}")

    # Save full results
    results_file = OUTPUT_DIR / "trinox_workflow_results.json"
    with open(results_file, 'w') as f:
        json.dump({
            "timestamp": datetime.now().isoformat(),
            "total_backbones": 13,
            "total_sequences": len(results),
            "high_quality_count": len(high_quality),
            "medium_quality_count": len(medium_quality),
            "results": sorted_results,
        }, f, indent=2)
    print(f"\nFull results saved to {results_file}")


if __name__ == "__main__":
    main()
