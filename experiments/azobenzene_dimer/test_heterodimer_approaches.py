#!/usr/bin/env python3
"""
Test script for Heterodimer Design Approaches

Tests all three approaches for cis-azobenzene heterodimer design:
- Approach 5: Joint Multi-Chain Diffusion
- Approach 7: Asymmetric RASA
- Approach 8: Induced Dimerization

Generates 10 designs per approach (30 total) and compares results.
"""

import json
import os
import sys
import time
from datetime import datetime
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed

# API endpoint
API_URL = "http://localhost:8000/runsync"

# Cis-azobenzene SMILES (Z-isomer, bent conformer)
CIS_AZOBENZENE_SMILES = "c1ccc(/N=N\\c2ccccc2)cc1"

# Test configuration
NUM_DESIGNS_PER_APPROACH = 10
CHAIN_LENGTH = "50-70"
BASE_SEED = 100


def call_api(payload):
    """Make API call and return result."""
    try:
        response = requests.post(
            API_URL,
            json=payload,
            headers={"Content-Type": "application/json"},
            timeout=600  # 10 minute timeout per design
        )
        return response.json()
    except Exception as e:
        return {"status": "error", "error": str(e)}


def run_approach(approach_name, approach_type, design_idx):
    """Run a single design for an approach."""
    print(f"[{approach_name}] Starting design {design_idx + 1}/{NUM_DESIGNS_PER_APPROACH}...")

    payload = {
        "input": {
            "task": "interface_ligand_design",
            "ligand_smiles": CIS_AZOBENZENE_SMILES,
            "approach": approach_type,
            "chain_length": CHAIN_LENGTH,
            "num_designs": 1,
            "seed": BASE_SEED + (design_idx * 100),
        }
    }

    start_time = time.time()
    result = call_api(payload)
    elapsed = time.time() - start_time

    return {
        "approach": approach_name,
        "design_index": design_idx,
        "elapsed_seconds": elapsed,
        "result": result,
    }


def extract_metrics(result):
    """Extract key metrics from design result."""
    # API returns {"id": ..., "status": "COMPLETED", "output": {"status": "completed", ...}}
    # Need to extract from the output field
    output = result.get("output", result)  # Fall back to result if no output wrapper

    if output.get("status") != "completed":
        return {
            "success": False,
            "error": output.get("error", "Unknown error"),
        }

    dimer = output.get("result", {}).get("dimer", {})
    metrics = dimer.get("metrics", {})
    anti_homo = dimer.get("anti_homodimerization", {})

    return {
        "success": True,
        "affinity": metrics.get("affinity"),
        "contacts_a": metrics.get("contacts_a"),
        "contacts_b": metrics.get("contacts_b"),
        "has_clashes": metrics.get("has_clashes"),
        "separable": metrics.get("separable"),
        "sequence_identity": metrics.get("sequence_identity"),
        "is_heterodimer": metrics.get("is_heterodimer"),
        "passes_anti_homo": anti_homo.get("passes_anti_homodimerization"),
        "ab_affinity": anti_homo.get("affinities", {}).get("ab_heterodimer"),
        "aa_affinity": anti_homo.get("affinities", {}).get("aa_homodimer"),
        "bb_affinity": anti_homo.get("affinities", {}).get("bb_homodimer"),
        "selectivity": anti_homo.get("selectivity_kcal"),
        "anti_homo_score": anti_homo.get("anti_homo_score"),
        # H-bond metrics for azobenzene N5/N6
        "n5_hbonds": metrics.get("n5_hbonds", 0),
        "n6_hbonds": metrics.get("n6_hbonds", 0),
        "total_ligand_hbonds": metrics.get("total_ligand_hbonds", 0),
    }


def save_design_pdb(approach_name, design_idx, result, output_dir):
    """Save design PDB to file."""
    output = result.get("output", result)  # Fall back to result if no output wrapper
    if output.get("status") != "completed":
        return None

    pdb_content = output.get("result", {}).get("dimer", {}).get("pdb_content")
    if not pdb_content:
        return None

    filename = f"{approach_name}_design_{design_idx:02d}.pdb"
    filepath = os.path.join(output_dir, filename)

    with open(filepath, 'w') as f:
        f.write(pdb_content)

    return filepath


def run_all_approaches():
    """Run all heterodimer approaches and compare results."""
    print("=" * 70)
    print("HETERODIMER DESIGN TEST - CIS-AZOBENZENE")
    print("=" * 70)
    print(f"Ligand: cis-azobenzene (bent conformer)")
    print(f"SMILES: {CIS_AZOBENZENE_SMILES}")
    print(f"Chain length: {CHAIN_LENGTH}")
    print(f"Designs per approach: {NUM_DESIGNS_PER_APPROACH}")
    print("=" * 70)

    # Create output directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = f"outputs/heterodimer_test_{timestamp}"
    os.makedirs(output_dir, exist_ok=True)

    # Define approaches to test
    approaches = [
        ("joint", "joint"),
        ("asymmetric_rasa", "asymmetric_rasa"),
        ("induced", "induced"),
    ]

    all_results = []

    # Test each approach
    for approach_name, approach_type in approaches:
        print(f"\n{'='*50}")
        print(f"Testing approach: {approach_name.upper()}")
        print(f"{'='*50}")

        approach_results = []
        approach_start = time.time()

        for design_idx in range(NUM_DESIGNS_PER_APPROACH):
            result = run_approach(approach_name, approach_type, design_idx)
            approach_results.append(result)

            # Extract and print metrics
            metrics = extract_metrics(result["result"])
            if metrics["success"]:
                print(f"  Design {design_idx + 1}: affinity={metrics['affinity']}, "
                      f"identity={metrics['sequence_identity']:.1f}%, "
                      f"anti_homo={metrics['passes_anti_homo']}, "
                      f"n5_hb={metrics['n5_hbonds']}, n6_hb={metrics['n6_hbonds']}")

                # Save PDB
                save_design_pdb(approach_name, design_idx, result["result"], output_dir)
            else:
                print(f"  Design {design_idx + 1}: FAILED - {metrics.get('error', 'Unknown')}")

        approach_elapsed = time.time() - approach_start
        print(f"\n{approach_name} completed in {approach_elapsed:.1f}s")

        all_results.extend(approach_results)

    # Analyze results
    print("\n" + "=" * 70)
    print("RESULTS SUMMARY")
    print("=" * 70)

    summary = {}
    for approach_name, _ in approaches:
        approach_data = [r for r in all_results if r["approach"] == approach_name]
        metrics_list = [extract_metrics(r["result"]) for r in approach_data]

        successful = [m for m in metrics_list if m["success"]]
        anti_homo_passed = [m for m in successful if m.get("passes_anti_homo")]
        heterodimers = [m for m in successful if m.get("is_heterodimer")]

        affinities = [m["affinity"] for m in successful if m.get("affinity") is not None]
        identities = [m["sequence_identity"] for m in successful if m.get("sequence_identity") is not None]

        # H-bond metrics
        n5_hbonds_list = [m["n5_hbonds"] for m in successful]
        n6_hbonds_list = [m["n6_hbonds"] for m in successful]
        total_hbonds_list = [m["total_ligand_hbonds"] for m in successful]
        designs_with_n5 = len([h for h in n5_hbonds_list if h >= 1])
        designs_with_n6 = len([h for h in n6_hbonds_list if h >= 1])
        designs_with_both = len([i for i, (n5, n6) in enumerate(zip(n5_hbonds_list, n6_hbonds_list)) if n5 >= 1 and n6 >= 1])

        summary[approach_name] = {
            "total": len(approach_data),
            "successful": len(successful),
            "success_rate": len(successful) / len(approach_data) * 100 if approach_data else 0,
            "anti_homo_passed": len(anti_homo_passed),
            "anti_homo_rate": len(anti_homo_passed) / len(successful) * 100 if successful else 0,
            "heterodimers": len(heterodimers),
            "heterodimer_rate": len(heterodimers) / len(successful) * 100 if successful else 0,
            "avg_affinity": sum(affinities) / len(affinities) if affinities else None,
            "best_affinity": min(affinities) if affinities else None,
            "avg_identity": sum(identities) / len(identities) if identities else None,
            # H-bond statistics
            "avg_n5_hbonds": sum(n5_hbonds_list) / len(n5_hbonds_list) if n5_hbonds_list else 0,
            "avg_n6_hbonds": sum(n6_hbonds_list) / len(n6_hbonds_list) if n6_hbonds_list else 0,
            "avg_total_hbonds": sum(total_hbonds_list) / len(total_hbonds_list) if total_hbonds_list else 0,
            "designs_with_n5_hbond": designs_with_n5,
            "designs_with_n6_hbond": designs_with_n6,
            "designs_with_both_hbonds": designs_with_both,
            "n5_hbond_rate": designs_with_n5 / len(successful) * 100 if successful else 0,
            "n6_hbond_rate": designs_with_n6 / len(successful) * 100 if successful else 0,
            "both_hbond_rate": designs_with_both / len(successful) * 100 if successful else 0,
        }

        print(f"\n{approach_name.upper()}:")
        print(f"  Success rate: {summary[approach_name]['success_rate']:.1f}% ({summary[approach_name]['successful']}/{summary[approach_name]['total']})")
        print(f"  Anti-homodimerization passed: {summary[approach_name]['anti_homo_passed']}/{summary[approach_name]['successful']} ({summary[approach_name]['anti_homo_rate']:.1f}%)")
        print(f"  True heterodimers (identity < 70%): {summary[approach_name]['heterodimers']}/{summary[approach_name]['successful']} ({summary[approach_name]['heterodimer_rate']:.1f}%)")
        if summary[approach_name]['avg_affinity'] is not None:
            print(f"  Average affinity: {summary[approach_name]['avg_affinity']:.2f} kcal/mol")
            print(f"  Best affinity: {summary[approach_name]['best_affinity']:.2f} kcal/mol")
        if summary[approach_name]['avg_identity'] is not None:
            print(f"  Average sequence identity: {summary[approach_name]['avg_identity']:.1f}%")
        # H-bond summary
        print(f"  --- H-bond Metrics (Azobenzene N5/N6) ---")
        print(f"  Designs with N5 H-bond: {designs_with_n5}/{len(successful)} ({summary[approach_name]['n5_hbond_rate']:.1f}%)")
        print(f"  Designs with N6 H-bond: {designs_with_n6}/{len(successful)} ({summary[approach_name]['n6_hbond_rate']:.1f}%)")
        print(f"  Designs with BOTH H-bonds: {designs_with_both}/{len(successful)} ({summary[approach_name]['both_hbond_rate']:.1f}%)")
        print(f"  Average total ligand H-bonds: {summary[approach_name]['avg_total_hbonds']:.2f}")

    # Save full results
    results_file = os.path.join(output_dir, "results.json")
    with open(results_file, 'w') as f:
        json.dump({
            "timestamp": timestamp,
            "config": {
                "ligand_smiles": CIS_AZOBENZENE_SMILES,
                "chain_length": CHAIN_LENGTH,
                "num_designs_per_approach": NUM_DESIGNS_PER_APPROACH,
            },
            "summary": summary,
            "all_results": [
                {
                    "approach": r["approach"],
                    "design_index": r["design_index"],
                    "elapsed_seconds": r["elapsed_seconds"],
                    "metrics": extract_metrics(r["result"]),
                }
                for r in all_results
            ]
        }, f, indent=2)

    print(f"\n{'='*70}")
    print(f"Results saved to: {output_dir}/")
    print(f"  - results.json (full metrics)")
    print(f"  - *.pdb (design structures)")
    print(f"{'='*70}")

    # Determine best approach (considering anti-homo rate, H-bond rate, and affinity)
    best_approach = max(summary.items(),
                        key=lambda x: (
                            x[1]["anti_homo_rate"],
                            x[1]["both_hbond_rate"],  # Prioritize designs with H-bonds to both N5 and N6
                            -x[1].get("avg_affinity", 0) if x[1].get("avg_affinity") else 0
                        ))

    print(f"\nBEST APPROACH: {best_approach[0].upper()}")
    print(f"  Anti-homodimerization rate: {best_approach[1]['anti_homo_rate']:.1f}%")
    print(f"  Both N5+N6 H-bond rate: {best_approach[1]['both_hbond_rate']:.1f}%")
    if best_approach[1]['avg_affinity']:
        print(f"  Average affinity: {best_approach[1]['avg_affinity']:.2f} kcal/mol")

    return summary


def test_single_approach(approach_type="joint"):
    """Test a single approach for debugging."""
    print(f"Testing single design with approach: {approach_type}")

    result = run_approach(approach_type, approach_type, 0)
    metrics = extract_metrics(result["result"])

    print(f"\nResult: {json.dumps(metrics, indent=2)}")

    if result["result"].get("status") == "completed":
        pdb = result["result"].get("result", {}).get("dimer", {}).get("pdb_content", "")
        if pdb:
            print(f"\nPDB preview (first 500 chars):\n{pdb[:500]}...")

    return result


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Test heterodimer design approaches")
    parser.add_argument("--single", type=str, help="Test single approach (joint, asymmetric_rasa, induced)")
    parser.add_argument("--designs", type=int, default=10, help="Number of designs per approach")

    args = parser.parse_args()

    if args.single:
        test_single_approach(args.single)
    else:
        NUM_DESIGNS_PER_APPROACH = args.designs
        run_all_approaches()
