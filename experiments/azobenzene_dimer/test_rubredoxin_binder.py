"""
Test the protein_binder_design pipeline with Rubredoxin (1RB9).

Rubredoxin is a small iron-sulfur protein (~52 residues) that makes
a good test target for binder design.
"""

import json
import requests
from pathlib import Path

API_URL = "http://localhost:8000/runsync"


def get_rubredoxin_pdb():
    """Fetch and clean rubredoxin PDB from RCSB."""
    import urllib.request

    url = "https://files.rcsb.org/download/1RB9.pdb"
    response = urllib.request.urlopen(url)
    pdb_content = response.read().decode('utf-8')

    # Clean the PDB: keep only chain A heavy atoms, handle alt conformations
    seen = set()
    clean_lines = []

    for line in pdb_content.split('\n'):
        if line.startswith('ATOM'):
            chain = line[21]
            if chain != 'A':
                continue
            # Skip hydrogens
            atom_name = line[12:16].strip()
            if atom_name.startswith('H') or (len(atom_name) > 1 and atom_name[0].isdigit() and atom_name[1] == 'H'):
                continue
            # Handle alternate conformations - take first one only
            alt_loc = line[16]
            res_num = line[22:26].strip()
            key = (res_num, atom_name)
            if key in seen:
                continue
            seen.add(key)
            # Remove alternate location indicator
            if alt_loc != ' ':
                line = line[:16] + ' ' + line[17:]
            clean_lines.append(line)

    clean_lines.append('END')
    return '\n'.join(clean_lines)


def test_rubredoxin_binder(
    num_designs: int = 2,
    binder_length: str = "40-60",
    quality_threshold: str = "relaxed",
    hotspots: list = None,
    run_esm_scoring: bool = True,
):
    """Run the protein binder design pipeline on rubredoxin."""

    print("Fetching rubredoxin PDB (1RB9)...")
    target_pdb = get_rubredoxin_pdb()
    print(f"Target PDB: {len(target_pdb)} characters, chain A")

    # Default hotspots: surface residues near the iron center
    # Residues 6, 9, 39, 42 are cysteines coordinating iron - nearby residues are good targets
    if hotspots is None:
        hotspots = ["A10", "A11", "A12", "A40", "A41"]  # Surface residues

    request = {
        "input": {
            "task": "protein_binder_design",
            "target_pdb": target_pdb,
            "hotspots": hotspots,
            "binder_length": binder_length,
            "num_designs": num_designs,
            "quality_threshold": quality_threshold,
            "run_mpnn": True,
            "run_esm_scoring": run_esm_scoring,
            "run_fastrelax": False,  # Skip for faster testing
        }
    }

    try:
        print(f"\nSending protein_binder_design request...")
        print(f"  Target: Rubredoxin (1RB9, 52 residues)")
        print(f"  Hotspots: {hotspots}")
        print(f"  Binder length: {binder_length}")
        print(f"  Designs requested: {num_designs}")
        print(f"  Quality threshold: {quality_threshold}")
        print()

        response = requests.post(API_URL, json=request, timeout=1800)
        result = response.json()

        # Handle RunPod wrapper
        if "output" in result:
            result = result["output"]

        return result
    except requests.exceptions.Timeout:
        return {"status": "error", "error": "Request timed out after 30 minutes"}
    except Exception as e:
        return {"status": "error", "error": str(e)}


def main():
    print("=" * 60)
    print("Testing protein_binder_design with Rubredoxin (1RB9)")
    print("=" * 60)
    print()

    result = test_rubredoxin_binder(
        num_designs=2,
        binder_length="40-60",
        quality_threshold="relaxed",
        run_esm_scoring=True,  # ESM now has HF_TOKEN set
    )

    if result.get("status") == "error":
        print(f"\nError: {result.get('error')}")
        return

    # Print statistics
    stats = result.get("statistics", {})
    print("\n" + "=" * 60)
    print("PIPELINE STATISTICS")
    print("=" * 60)
    for key, value in stats.items():
        print(f"  {key}: {value}")

    # Print design results
    designs = result.get("designs", [])
    print("\n" + "=" * 60)
    print(f"DESIGNS RETURNED: {len(designs)}")
    print("=" * 60)

    output_dir = Path(__file__).parent / "outputs" / "rubredoxin_binder_test"
    output_dir.mkdir(parents=True, exist_ok=True)

    for i, design in enumerate(designs):
        rank = design.get("rank", i + 1)
        print(f"\nDesign {rank}:")

        for key in ["interface_contacts", "interface_hbonds", "buried_sasa",
                    "esm_perplexity", "esm_confidence", "rosetta_energy"]:
            value = design.get(key)
            if value is not None:
                if isinstance(value, float):
                    print(f"  {key}: {value:.2f}")
                else:
                    print(f"  {key}: {value}")

        # Print binder sequence
        binder_seq = design.get("binder_sequence")
        if binder_seq:
            print(f"  Binder sequence ({len(binder_seq)} aa): {binder_seq[:50]}...")

        # Save PDB
        pdb_content = design.get("pdb_content")
        if pdb_content:
            pdb_path = output_dir / f"rubredoxin_binder_rank_{rank}.pdb"
            with open(pdb_path, 'w') as f:
                f.write(pdb_content)
            print(f"  Saved: {pdb_path.name}")

    # Save full results
    results_path = output_dir / "rubredoxin_binder_results.json"
    results_for_json = {
        "status": result.get("status"),
        "statistics": stats,
        "designs": [
            {k: v for k, v in d.items() if k != "pdb_content"}
            for d in designs
        ]
    }
    with open(results_path, 'w') as f:
        json.dump(results_for_json, f, indent=2, default=str)
    print(f"\nResults saved to: {results_path}")

    print("\n" + "=" * 60)
    print("TEST COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
