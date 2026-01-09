#!/usr/bin/env python3
"""
Local Testing Script for RFdiffusion Serverless

Usage:
    python test_local.py [test_name]

Tests:
    health   - Check if server is running and GPU available
    rfd3     - Test RFD3 protein design
    rf3      - Test RF3 structure prediction
    mpnn     - Test MPNN sequence design
    all      - Run all tests

Examples:
    python test_local.py health
    python test_local.py rfd3
    python test_local.py all
"""

import requests
import json
import sys
import time

BASE_URL = "http://localhost:8000"

def call_api(task: str, params: dict = None, timeout: int = 300):
    """Call the local RunPod API"""
    payload = {
        "input": {
            "task": task,
            **(params or {})
        }
    }

    print(f"\n{'='*60}")
    print(f"Testing: {task.upper()}")
    print(f"{'='*60}")
    print(f"Request: {json.dumps(payload, indent=2)[:500]}...")

    try:
        start = time.time()
        response = requests.post(
            f"{BASE_URL}/runsync",
            json=payload,
            headers={"Content-Type": "application/json"},
            timeout=timeout
        )
        elapsed = time.time() - start

        print(f"\nStatus Code: {response.status_code}")
        print(f"Time: {elapsed:.2f}s")

        result = response.json()

        # Pretty print result (truncate large outputs)
        result_str = json.dumps(result, indent=2)
        if len(result_str) > 2000:
            print(f"Response (truncated):\n{result_str[:2000]}...")
        else:
            print(f"Response:\n{result_str}")

        # Check for errors
        if result.get("status") == "failed":
            print(f"\n❌ FAILED: {result.get('error')}")
            if result.get("traceback"):
                print(f"\nTraceback:\n{result.get('traceback')}")
            return False
        else:
            print(f"\n✅ SUCCESS")
            return True

    except requests.exceptions.ConnectionError:
        print(f"\n❌ Connection Error: Is the container running?")
        print(f"   Start with: docker-compose -f docker-compose.local.yml up")
        return False
    except requests.exceptions.Timeout:
        print(f"\n❌ Timeout after {timeout}s")
        return False
    except Exception as e:
        print(f"\n❌ Error: {e}")
        return False


def test_health():
    """Test health endpoint"""
    return call_api("health")


def test_rfd3():
    """Test RFD3 protein design"""
    return call_api("rfd3", {
        "contig": "100",  # Simple 100-residue de novo design
        "num_designs": 1,
        "seed": 42
    }, timeout=600)


def test_rfd3_binder():
    """Test RFD3 binder design with input PDB"""
    # Example: design a binder to a target
    # You would need to provide actual PDB content here
    sample_pdb = """ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.246   2.390   0.000  1.00  0.00           O
END
"""
    return call_api("rfd3", {
        "contig": "A1-4/0 50-50",  # Fixed A1-4, then design 50 residues
        "num_designs": 1,
        "pdb_content": sample_pdb,
        "seed": 42
    }, timeout=600)


def test_rf3():
    """Test RF3 structure prediction"""
    return call_api("rf3", {
        "sequence": "MKFLILLFNILCLFPVLAADNHGVGPQGASGVDPITFDINSNQTGVQLTLGN",
        "name": "test_protein"
    }, timeout=600)


def test_mpnn():
    """Test MPNN sequence design"""
    # Simple test PDB
    sample_pdb = """ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.246   2.390   0.000  1.00  0.00           O
ATOM      5  N   GLY A   2       3.332   1.554   0.000  1.00  0.00           N
ATOM      6  CA  GLY A   2       3.988   2.847   0.000  1.00  0.00           C
ATOM      7  C   GLY A   2       5.504   2.720   0.000  1.00  0.00           C
ATOM      8  O   GLY A   2       6.088   1.628   0.000  1.00  0.00           O
END
"""
    return call_api("mpnn", {
        "pdb_content": sample_pdb,
        "num_sequences": 4,
        "temperature": 0.1
    }, timeout=300)


def run_all_tests():
    """Run all tests"""
    results = {}

    tests = [
        ("health", test_health),
        ("rfd3", test_rfd3),
        ("rf3", test_rf3),
        ("mpnn", test_mpnn),
    ]

    for name, test_fn in tests:
        results[name] = test_fn()
        time.sleep(1)  # Brief pause between tests

    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    for name, passed in results.items():
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"  {name}: {status}")

    all_passed = all(results.values())
    print(f"\nOverall: {'✅ ALL TESTS PASSED' if all_passed else '❌ SOME TESTS FAILED'}")
    return all_passed


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    test_name = sys.argv[1].lower()

    test_map = {
        "health": test_health,
        "rfd3": test_rfd3,
        "rfd3_binder": test_rfd3_binder,
        "rf3": test_rf3,
        "mpnn": test_mpnn,
        "all": run_all_tests,
    }

    if test_name not in test_map:
        print(f"Unknown test: {test_name}")
        print(f"Available tests: {', '.join(test_map.keys())}")
        sys.exit(1)

    success = test_map[test_name]()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
