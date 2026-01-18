"""
Test FastRelax integration for resolving ligand-protein clashes.

Run after rebuilding the Docker container with PyRosetta:
    python test_fastrelax.py

This script tests the fastrelax API task on designs with known clashes.
"""

import json
import requests
from pathlib import Path

API_URL = "http://localhost:8000/runsync"
AZOBENZENE_SMILES = "c1ccc(cc1)N=Nc2ccccc2"


def test_fastrelax(pdb_path: str) -> dict:
    """Run FastRelax via API on a design with known clashes."""
    with open(pdb_path, 'r') as f:
        pdb_content = f.read()

    request = {
        "input": {
            "task": "fastrelax",
            "pdb_content": pdb_content,
            "ligand_smiles": AZOBENZENE_SMILES,
            "max_iter": 200,
            "interface_only": True,
            "run_gnina_before": True,
            "run_gnina_after": True,
        }
    }

    try:
        print("Sending FastRelax request...")
        response = requests.post(API_URL, json=request, timeout=600)
        result = response.json()

        # Handle RunPod wrapper
        if "output" in result:
            result = result["output"]

        return result
    except Exception as e:
        return {"error": str(e)}


def main():
    # Designs with known clashes (positive GNINA affinity)
    designs_to_test = [
        "45_65_s62_dimer.pdb",  # User said: "azobenzene penetrate the dimer"
        "symmetric_40_60_dimer.pdb",  # User said: "azobenzene too close to backbone"
        "35_55_s50_dimer.pdb",  # User said: "also clashing"
    ]

    symmetric_dir = Path(__file__).parent / "outputs" / "approach_symmetric"
    results = []

    for design_name in designs_to_test:
        pdb_path = symmetric_dir / design_name
        if not pdb_path.exists():
            print(f"File not found: {pdb_path}")
            continue

        print(f"\n{'='*60}")
        print(f"Testing FastRelax: {design_name}")
        print("="*60)

        result = test_fastrelax(str(pdb_path))

        if result.get("status") == "error":
            print(f"Error: {result.get('error')}")
            results.append({
                "design": design_name,
                "status": "error",
                "error": result.get("error")
            })
            continue

        # Extract results
        improvement = result.get("improvement", {})

        print(f"\nRosetta Energy:")
        print(f"  Before: {result.get('energy_before', 'N/A')}")
        print(f"  After:  {result.get('energy_after', 'N/A')}")
        print(f"  Change: {result.get('energy_change', 'N/A')}")

        gnina_before = result.get("gnina_before", {})
        gnina_after = result.get("gnina_after", {})

        print(f"\nGNINA Affinity:")
        print(f"  Before: {gnina_before.get('best_affinity', 'N/A')} kcal/mol")
        print(f"  After:  {gnina_after.get('best_affinity', 'N/A')} kcal/mol")
        print(f"  Improved: {improvement.get('gnina_improved', 'N/A')}")

        print(f"\nCNN Score:")
        print(f"  Before: {gnina_before.get('best_cnn_score', 'N/A')}")
        print(f"  After:  {gnina_after.get('best_cnn_score', 'N/A')}")
        print(f"  Improved: {improvement.get('cnn_improved', 'N/A')}")

        # Save relaxed PDB
        if result.get("status") == "completed" and result.get("relaxed_pdb"):
            relaxed_path = pdb_path.with_stem(pdb_path.stem + "_relaxed")
            with open(relaxed_path, 'w') as f:
                f.write(result["relaxed_pdb"])
            print(f"\nRelaxed structure saved: {relaxed_path.name}")

        results.append({
            "design": design_name,
            "status": result.get("status"),
            "rosetta_before": result.get("energy_before"),
            "rosetta_after": result.get("energy_after"),
            "rosetta_change": result.get("energy_change"),
            "gnina_before": gnina_before.get("best_affinity"),
            "gnina_after": gnina_after.get("best_affinity"),
            "gnina_improved": improvement.get("gnina_improved"),
            "cnn_before": gnina_before.get("best_cnn_score"),
            "cnn_after": gnina_after.get("best_cnn_score"),
            "full_result": result
        })

    # Save results
    output_file = Path(__file__).parent / "outputs" / "fastrelax_test_results.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)

    print(f"\n\n{'='*60}")
    print("SUMMARY")
    print("="*60)

    for r in results:
        if r.get("status") == "completed":
            gnina_delta = (r.get("gnina_after") or 0) - (r.get("gnina_before") or 0)
            improved = "YES" if r.get("gnina_improved") else "NO"
            print(f"\n{r['design']}:")
            print(f"  GNINA: {r.get('gnina_before'):.2f} -> {r.get('gnina_after'):.2f} ({gnina_delta:+.2f})")
            print(f"  Rosetta: {r.get('rosetta_before'):.1f} -> {r.get('rosetta_after'):.1f} ({r.get('rosetta_change'):+.1f})")
            print(f"  Improved: {improved}")

    print(f"\n\nResults saved to: {output_file}")


if __name__ == "__main__":
    main()
