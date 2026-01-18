"""
Run GNINA scoring on valid symmetric designs via API.
"""

import json
import requests
from pathlib import Path

API_URL = "http://localhost:8000/runsync"
AZOBENZENE_SMILES = "c1ccc(cc1)N=Nc2ccccc2"

def run_gnina_scoring(pdb_path: str) -> dict:
    """Run GNINA scoring via API."""
    with open(pdb_path, 'r') as f:
        pdb_content = f.read()

    request = {
        "input": {
            "task": "binding_eval",
            "pdb_content": pdb_content,
            "ligand_smiles": AZOBENZENE_SMILES,
            "chain_a": "A",
            "chain_b": "B",
            "run_gnina": True,
            "check_clashes": True,
            "whole_protein_search": True,
            "docking_box_size": 30.0
        }
    }

    try:
        response = requests.post(API_URL, json=request, timeout=300)
        result = response.json()

        # Handle RunPod wrapper
        if "output" in result:
            result = result["output"]

        return result
    except Exception as e:
        return {"error": str(e)}


def main():
    # Top valid designs from validation
    valid_designs = [
        "45_65_s62_dimer.pdb",
        "symmetric_40_60_dimer.pdb",
        "35_55_s50_dimer.pdb",
        "45_65_s61_dimer.pdb",
        "sym_40_60_seed49_dimer.pdb",
    ]

    symmetric_dir = Path(__file__).parent / "outputs" / "approach_symmetric"

    results = []

    for design in valid_designs:
        pdb_path = symmetric_dir / design
        if not pdb_path.exists():
            print(f"File not found: {pdb_path}")
            continue

        print(f"\n{'='*60}")
        print(f"Scoring: {design}")
        print("="*60)

        result = run_gnina_scoring(str(pdb_path))

        # Extract key metrics
        gnina = result.get("gnina_scoring", {}).get("result", {})
        affinity = gnina.get("best_affinity", "N/A")
        cnn_score = gnina.get("best_cnn_score", "N/A")

        interface = result.get("interface_analysis", {}).get("metrics", {})
        contacts = interface.get("contacts", "N/A")
        hbonds = interface.get("hbonds_int", "N/A")

        clash = result.get("clash_check", {})
        has_clashes = clash.get("has_clashes", "N/A")

        composite = result.get("summary", {}).get("composite_score", {})
        total_score = composite.get("total", "N/A")
        quality = composite.get("quality", "N/A")

        print(f"GNINA Affinity: {affinity} kcal/mol")
        print(f"CNN Score: {cnn_score}")
        print(f"Interface Contacts: {contacts}")
        print(f"H-bonds: {hbonds}")
        print(f"Has Clashes: {has_clashes}")
        print(f"Composite Score: {total_score}/100 ({quality})")

        results.append({
            "design": design,
            "affinity": affinity,
            "cnn_score": cnn_score,
            "contacts": contacts,
            "hbonds": hbonds,
            "has_clashes": has_clashes,
            "composite_score": total_score,
            "quality": quality,
            "full_result": result
        })

    # Save results
    output_file = Path(__file__).parent / "outputs" / "gnina_scoring_results.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)

    print(f"\n\n{'='*60}")
    print("SUMMARY")
    print("="*60)

    # Sort by affinity (most negative first)
    sorted_results = sorted(
        [r for r in results if isinstance(r.get("affinity"), (int, float))],
        key=lambda x: x["affinity"]
    )

    print("\nRanked by GNINA affinity:")
    for i, r in enumerate(sorted_results, 1):
        print(f"  {i}. {r['design']}: {r['affinity']:.2f} kcal/mol (CNN: {r.get('cnn_score', 'N/A')})")

    print(f"\nResults saved to: {output_file}")


if __name__ == "__main__":
    main()
