"""
Dy-TriNOx scaffold design v2 - with explicit coordination constraints.

The v1 designs had no metal coordination because the protein was positioned
too far from the Dy center (~10A away instead of ~2.5A).

This version uses:
1. select_hotspots on the Dy to force close contact
2. H-bond conditioning for the apical oxygen atoms
3. Smaller proteins to ensure tight binding pocket
"""

import requests
import json
from pathlib import Path
from datetime import datetime

API_URL = "http://localhost:8000/runsync"
SCRIPT_DIR = Path(__file__).parent
INPUT_DIR = SCRIPT_DIR.parent / "inputs"
OUTPUT_DIR = SCRIPT_DIR.parent / "outputs_v2"

# Create output dir
OUTPUT_DIR.mkdir(exist_ok=True)


def run_rfd3_design(pdb_content: str, config: dict, timeout: int = 1200) -> list:
    """Run RFD3 design with coordination constraints."""

    # Build contig: include metal (X1) and ligand (L1), then add protein chain
    protein_length = config.get("contig", "50-70")
    contig = f"X1,L1,/0,{protein_length}"

    payload = {
        "input": {
            "task": "rfd3",
            "pdb_content": pdb_content,
            "contig": contig,
            "ligand": "DY,UNL",
            "num_designs": config.get("num_designs", 5),

            # Fix both metal and ligand in place
            "select_fixed_atoms": {
                "X1": "all",
                "L1": "all",
            },

            # CRITICAL: Use hotspots on Dy to force close contact
            "select_hotspots": config.get("select_hotspots", {"X1": "1"}),

            # Burial conditioning
            "select_buried": config.get("select_buried", {"X1": "all"}),

            # H-bond conditioning if specified
            # This can help pull protein donors toward the metal
            **({"select_hbond_acceptor": config["select_hbond_acceptor"]}
               if "select_hbond_acceptor" in config else {}),

            # CFG for better conditioning
            "use_classifier_free_guidance": True,
            "cfg_scale": config.get("cfg_scale", 2.5),

            # Position protein near the metal
            "infer_ori_strategy": "hotspots",
        }
    }

    print(f"  Contig: {contig}")
    print(f"  Hotspots: {config.get('select_hotspots', {'X1': '1'})}")

    try:
        response = requests.post(API_URL, json=payload, timeout=timeout)
        result = response.json()

        if result.get("status") != "COMPLETED":
            error = result.get("output", {}).get("error", "Unknown error")
            print(f"  ERROR: {error}")
            return []

        output = result.get("output", {})
        res = output.get("result", {})
        predictions = res.get("predictions", [])

        designs = []
        for pred in predictions:
            if "content" in pred:
                designs.append(pred["content"])

        return designs

    except Exception as e:
        print(f"  Error: {e}")
        return []


def main():
    print("=" * 70)
    print("Dy-TriNOx SCAFFOLD DESIGN v2")
    print("With explicit coordination constraints")
    print("=" * 70)

    # Load input PDB
    input_pdb = INPUT_DIR / "Dy_TriNOx_split.pdb"
    with open(input_pdb, 'r') as f:
        pdb_content = f.read()

    # Design configs with different hotspot/conditioning strategies
    configs = [
        # Strategy 1: Hotspot on Dy only, small protein
        {
            "name": "v2a_hotspot_dy",
            "contig": "40-60",
            "num_designs": 5,
            "select_hotspots": {"X1": "1"},  # Hotspot on Dy
            "select_buried": {"X1": "all"},
            "cfg_scale": 2.5,
        },
        # Strategy 2: Hotspot on Dy + burial of ligand oxygens
        {
            "name": "v2b_hotspot_burial",
            "contig": "50-70",
            "num_designs": 5,
            "select_hotspots": {"X1": "1"},
            "select_buried": {"X1": "all", "L1": "O1,O2,O3"},  # Bury oxygens
            "cfg_scale": 3.0,
        },
        # Strategy 3: Very small protein for tight binding
        {
            "name": "v2c_minimal",
            "contig": "30-45",
            "num_designs": 5,
            "select_hotspots": {"X1": "1"},
            "select_buried": {"X1": "all"},
            "cfg_scale": 2.5,
        },
    ]

    all_designs = []
    summary = {
        "timestamp": datetime.now().isoformat(),
        "input_pdb": str(input_pdb),
        "total_designs": 0,
        "designs": [],
        "configs": configs,
    }

    for config in configs:
        print(f"\n{'='*70}")
        print(f"Config: {config['name']}")
        print(f"Length: {config['contig']} residues")
        print(f"={'='*70}")

        designs = run_rfd3_design(pdb_content, config)

        print(f"Generated {len(designs)} designs")

        for i, design in enumerate(designs):
            design_name = f"{config['name']}_{i:03d}"
            design_file = OUTPUT_DIR / f"{design_name}.pdb"

            with open(design_file, 'w') as f:
                f.write(design)

            all_designs.append(design_name)
            summary["designs"].append(design_name)
            print(f"  Saved: {design_name}")

    summary["total_designs"] = len(all_designs)

    # Save summary
    summary_file = OUTPUT_DIR / "run_summary_v2.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"\n{'='*70}")
    print(f"COMPLETE: Generated {len(all_designs)} designs")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
