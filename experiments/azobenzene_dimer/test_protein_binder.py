"""
Test the protein_binder_design API task for PROTEIN-PROTEIN binding.

This tests the full protein-protein binder design pipeline:
1. Structure Generation (RFD3 with hotspots)
2. Sequence Design (ProteinMPNN)
3. Sequence Scoring (ESM-3)
4. Structure Refinement (FastRelax)
5. Interface Analysis
6. Quality Filtering
7. Ranking

Run after rebuilding the Docker container:
    python test_protein_binder.py

Note: Requires a target PDB file. A sample mini-protein is included below.
"""

import json
import requests
from pathlib import Path

API_URL = "http://localhost:8000/runsync"

# Sample mini-protein target (a simple 3-helix bundle, ~60 residues)
# This is a de novo designed protein for testing purposes
SAMPLE_TARGET_PDB = """HEADER    DE NOVO PROTEIN
TITLE     SAMPLE TARGET FOR BINDER DESIGN TESTING
ATOM      1  N   MET A   1      -0.525   1.362   0.000  1.00  0.00           N
ATOM      2  CA  MET A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   MET A   1       1.520   0.000   0.000  1.00  0.00           C
ATOM      4  O   MET A   1       2.153   1.055   0.000  1.00  0.00           O
ATOM      5  N   GLU A   2       2.061  -1.218   0.000  1.00  0.00           N
ATOM      6  CA  GLU A   2       3.502  -1.414   0.000  1.00  0.00           C
ATOM      7  C   GLU A   2       4.200  -0.065   0.000  1.00  0.00           C
ATOM      8  O   GLU A   2       5.430  -0.008   0.000  1.00  0.00           O
ATOM      9  N   ALA A   3       3.445   1.031   0.000  1.00  0.00           N
ATOM     10  CA  ALA A   3       3.970   2.390   0.000  1.00  0.00           C
ATOM     11  C   ALA A   3       5.490   2.496   0.000  1.00  0.00           C
ATOM     12  O   ALA A   3       6.135   3.547   0.000  1.00  0.00           O
ATOM     13  N   LEU A   4       6.031   1.280   0.000  1.00  0.00           N
ATOM     14  CA  LEU A   4       7.472   1.095   0.000  1.00  0.00           C
ATOM     15  C   LEU A   4       8.170   2.443   0.000  1.00  0.00           C
ATOM     16  O   LEU A   4       9.400   2.500   0.000  1.00  0.00           O
ATOM     17  N   LYS A   5       7.415   3.537   0.000  1.00  0.00           N
ATOM     18  CA  LYS A   5       7.940   4.896   0.000  1.00  0.00           C
ATOM     19  C   LYS A   5       9.460   5.002   0.000  1.00  0.00           C
ATOM     20  O   LYS A   5      10.105   6.053   0.000  1.00  0.00           O
ATOM     21  N   ARG A   6      10.001   3.786   0.000  1.00  0.00           N
ATOM     22  CA  ARG A   6      11.442   3.600   0.000  1.00  0.00           C
ATOM     23  C   ARG A   6      12.140   4.949   0.000  1.00  0.00           C
ATOM     24  O   ARG A   6      13.370   5.006   0.000  1.00  0.00           O
ATOM     25  N   ILE A   7      11.385   6.043   0.000  1.00  0.00           N
ATOM     26  CA  ILE A   7      11.910   7.402   0.000  1.00  0.00           C
ATOM     27  C   ILE A   7      13.430   7.508   0.000  1.00  0.00           C
ATOM     28  O   ILE A   7      14.075   8.559   0.000  1.00  0.00           O
ATOM     29  N   GLU A   8      13.971   6.292   0.000  1.00  0.00           N
ATOM     30  CA  GLU A   8      15.412   6.106   0.000  1.00  0.00           C
ATOM     31  C   GLU A   8      16.110   7.455   0.000  1.00  0.00           C
ATOM     32  O   GLU A   8      17.340   7.512   0.000  1.00  0.00           O
ATOM     33  N   ALA A   9      15.355   8.549   0.000  1.00  0.00           N
ATOM     34  CA  ALA A   9      15.880   9.908   0.000  1.00  0.00           C
ATOM     35  C   ALA A   9      17.400  10.014   0.000  1.00  0.00           C
ATOM     36  O   ALA A   9      18.045  11.065   0.000  1.00  0.00           O
ATOM     37  N   LEU A  10      17.941   8.798   0.000  1.00  0.00           N
ATOM     38  CA  LEU A  10      19.382   8.612   0.000  1.00  0.00           C
ATOM     39  C   LEU A  10      20.080   9.961   0.000  1.00  0.00           C
ATOM     40  O   LEU A  10      21.310  10.018   0.000  1.00  0.00           O
END
"""


def load_target_pdb(pdb_path: str = None) -> str:
    """Load target PDB content from file or use sample."""
    if pdb_path:
        path = Path(pdb_path)
        if path.exists():
            return path.read_text()
        else:
            print(f"Warning: PDB file not found: {pdb_path}")
            print("Using sample mini-protein target...")
    return SAMPLE_TARGET_PDB


def test_protein_binder_design(
    target_pdb: str,
    hotspots: list = None,
    binder_length: str = "60-80",
    num_designs: int = 3,
    quality_threshold: str = "relaxed",
    run_mpnn: bool = True,
    run_esm_scoring: bool = True,
    run_fastrelax: bool = True,
) -> dict:
    """Run the full protein-protein binder design pipeline."""

    request = {
        "input": {
            "task": "protein_binder_design",
            "target_pdb": target_pdb,
            "binder_length": binder_length,
            "num_designs": num_designs,
            "quality_threshold": quality_threshold,
            "run_mpnn": run_mpnn,
            "run_esm_scoring": run_esm_scoring,
            "run_fastrelax": run_fastrelax,
        }
    }

    # Add hotspots if specified
    if hotspots:
        request["input"]["hotspots"] = hotspots

    try:
        print(f"Sending protein_binder_design request (timeout=1800s)...")
        print(f"  Binder length: {binder_length}")
        print(f"  Hotspots: {hotspots or 'auto (all surface residues)'}")
        print(f"  Designs requested: {num_designs}")
        print(f"  Quality threshold: {quality_threshold}")
        print(f"  Run MPNN: {run_mpnn}")
        print(f"  Run ESM scoring: {run_esm_scoring}")
        print(f"  Run FastRelax: {run_fastrelax}")
        print()

        response = requests.post(API_URL, json=request, timeout=1800)
        result = response.json()

        # Handle RunPod wrapper
        if "output" in result:
            result = result["output"]

        return result
    except requests.exceptions.Timeout:
        return {"error": "Request timed out after 30 minutes"}
    except Exception as e:
        return {"error": str(e)}


def main():
    print("=" * 60)
    print("Testing protein_binder_design pipeline")
    print("(Protein-Protein Binding)")
    print("=" * 60)
    print()

    # Load target PDB (use sample or provide your own path)
    target_pdb = load_target_pdb()
    print(f"Target PDB loaded ({len(target_pdb)} characters)")

    # Test with relaxed threshold (more designs should pass)
    result = test_protein_binder_design(
        target_pdb=target_pdb,
        hotspots=["A5", "A6", "A7"],  # Target residues 5-7 as hotspots
        binder_length="40-60",
        num_designs=3,
        quality_threshold="relaxed",
        run_mpnn=True,
        run_esm_scoring=True,
        run_fastrelax=True,
    )

    if result.get("status") == "error":
        print(f"Error: {result.get('error')}")
        return

    # Print statistics
    stats = result.get("statistics", {})
    print("\n" + "=" * 60)
    print("PIPELINE STATISTICS")
    print("=" * 60)
    print(f"  Generated:        {stats.get('generated', 0)}")
    print(f"  MPNN designed:    {stats.get('mpnn_designed', 0)}")
    print(f"  ESM scored:       {stats.get('esm_scored', 0)}")
    print(f"  Relaxed:          {stats.get('relaxed', 0)}")
    print(f"  Interface scored: {stats.get('interface_scored', 0)}")
    print(f"  Passed filters:   {stats.get('passed_filters', 0)}")
    print(f"  Returned:         {stats.get('returned', 0)}")

    # Print design results
    designs = result.get("designs", [])
    print("\n" + "=" * 60)
    print(f"TOP {len(designs)} DESIGNS")
    print("=" * 60)

    output_dir = Path(__file__).parent / "outputs" / "protein_binder_test"
    output_dir.mkdir(parents=True, exist_ok=True)

    for design in designs:
        rank = design.get("rank", "?")
        contacts = design.get("interface_contacts")
        hbonds = design.get("interface_hbonds")
        buried_sasa = design.get("buried_sasa")
        esm_perp = design.get("esm_perplexity")
        esm_conf = design.get("esm_confidence")
        rosetta_e = design.get("rosetta_energy")

        print(f"\nRank {rank}:")
        print(f"  Interface Contacts: {contacts}" if contacts is not None else "  Interface Contacts: N/A")
        print(f"  H-bonds:            {hbonds}" if hbonds is not None else "  H-bonds:            N/A")
        print(f"  Buried SASA:        {buried_sasa:.1f} Å²" if buried_sasa else "  Buried SASA:        N/A")
        print(f"  ESM Perplexity:     {esm_perp:.2f}" if esm_perp else "  ESM Perplexity:     N/A")
        print(f"  ESM Confidence:     {esm_conf:.3f}" if esm_conf else "  ESM Confidence:     N/A")
        print(f"  Rosetta Energy:     {rosetta_e:.1f}" if rosetta_e else "  Rosetta Energy:     N/A")

        # Print binder sequence
        binder_seq = design.get("binder_sequence")
        if binder_seq:
            print(f"  Binder Sequence:    {binder_seq[:40]}..." if len(binder_seq) > 40 else f"  Binder Sequence:    {binder_seq}")

        # Save PDB
        pdb_content = design.get("pdb_content")
        if pdb_content:
            pdb_path = output_dir / f"rank_{rank}_complex.pdb"
            with open(pdb_path, 'w') as f:
                f.write(pdb_content)
            print(f"  Saved: {pdb_path.name}")

    # Save full results
    results_path = output_dir / "protein_binder_results.json"
    # Remove PDB content for cleaner JSON
    results_for_json = {
        "status": result.get("status"),
        "statistics": stats,
        "thresholds_used": result.get("thresholds_used"),
        "designs": [
            {k: v for k, v in d.items() if k != "pdb_content"}
            for d in designs
        ]
    }
    with open(results_path, 'w') as f:
        json.dump(results_for_json, f, indent=2, default=str)
    print(f"\nFull results saved to: {results_path}")

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    if designs:
        best = designs[0]
        print(f"Best design (Rank 1):")
        print(f"  Interface Contacts: {best.get('interface_contacts', 'N/A')}")
        print(f"  H-bonds:            {best.get('interface_hbonds', 'N/A')}")
        print(f"  ESM Confidence:     {best.get('esm_confidence', 'N/A')}")

        # Calculate pass rates
        generated = stats.get("generated", 1)
        passed = stats.get("passed_filters", 0)
        pass_rate = (passed / generated * 100) if generated > 0 else 0
        print(f"\nPass rate: {passed}/{generated} ({pass_rate:.1f}%)")
    else:
        print("No designs passed filters. Try using 'relaxed' threshold.")


if __name__ == "__main__":
    main()
