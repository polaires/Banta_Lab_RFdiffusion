"""
Run Dy-TriNOx protein design using RFdiffusion3 via local Docker API.

Key differences from citrate-Ln:
- Only 1-2 coordination sites needed (TriNOx provides 7 of 8-9)
- No H-bond conditioning (no free carboxylates on TriNOx)
- Larger ligand requires larger pocket

Usage:
    python run_trinox_design.py --quick   # Quick test with 2 designs
    python run_trinox_design.py           # Full campaign
"""

import argparse
import json
import requests
from pathlib import Path
from datetime import datetime

API_URL = "http://localhost:8000/runsync"
SCRIPT_DIR = Path(__file__).parent
INPUT_DIR = SCRIPT_DIR.parent / "inputs"
OUTPUT_DIR = SCRIPT_DIR.parent / "outputs"


def run_rfd3_design(pdb_content: str, config: dict, timeout: int = 1200) -> list:
    """
    Run RFD3 design via local API.

    Returns list of PDB content strings.

    Strategy (per RFD3 skill reference):
    - pdb_content: Contains Dy (chain X) and TriNOx (chain L)
    - contig: "X,L,/0,70-90" - include ligand chains, then design protein
    - ligand: "DY,UNL" - tell RFD3 which residues are ligands

    Per the RFD3 docs: "L,/0,80-120" means ligand chain L + new protein chain.
    """
    # Build contig: include metal (X1) and ligand (L1), then add protein chain
    # Format: "X1,L1,/0,70-90" = metal + ligand + chain break + new 70-90 res protein
    protein_length = config.get("contig", "80-100")
    contig = f"X1,L1,/0,{protein_length}"

    payload = {
        "input": {
            "task": "rfd3",
            "pdb_content": pdb_content,
            "contig": contig,  # Include ligand chains in contig
            "ligand": "DY,UNL",  # Residue names of metal and ligand
            "num_designs": config.get("num_designs", 5),

            # Fix both metal and ligand in place
            "select_fixed_atoms": {
                "X1": "all",  # Dy metal on chain X
                "L1": "all",  # TriNOx ligand on chain L
            },

            # Bury the metal for pocket formation
            "select_buried": config.get("select_buried", {"X1": "all"}),

            # CFG guidance for conditioning
            "use_classifier_free_guidance": config.get("use_cfg", True),
            "cfg_scale": config.get("cfg_scale", 2.0),

            # Origin strategy - center on ligand
            "infer_ori_strategy": "com",
        }
    }

    print(f"  Submitting to API: {API_URL}")
    print(f"  Contig: {config.get('contig')}, Designs: {config.get('num_designs')}")

    try:
        response = requests.post(API_URL, json=payload, timeout=timeout)
        result = response.json()

        status = result.get("status", result.get("output", {}).get("status", "unknown"))
        print(f"  Status: {status}")

        if status.upper() == "COMPLETED":
            output = result.get("output", {})
            res = output.get("result", {})
            pdbs = res.get("pdbs", res.get("designs", []))

            # Extract PDB content
            pdb_contents = []
            for pdb_data in pdbs:
                if isinstance(pdb_data, dict) and "content" in pdb_data:
                    pdb_contents.append(pdb_data["content"])
                elif isinstance(pdb_data, str):
                    pdb_contents.append(pdb_data)

            return pdb_contents
        else:
            print(f"  Failed: {result}")
            return []

    except requests.exceptions.ConnectionError:
        print(f"  ERROR: Cannot connect to API at {API_URL}")
        print(f"  Make sure the Docker container is running!")
        print(f"  Try: docker start rfd3-api")
        return []
    except Exception as e:
        print(f"  Exception: {e}")
        return []


def run_quick_test():
    """Run a quick test with 2 designs."""
    print("=" * 70)
    print("QUICK TEST - Dy-TriNOx Protein Design")
    print("=" * 70)

    # Read PDB with Dy (chain X) and TriNOx (chain L, UNL residue)
    pdb_path = INPUT_DIR / "Dy_TriNOx_split.pdb"
    if not pdb_path.exists():
        print(f"ERROR: PDB not found at {pdb_path}")
        return []

    with open(pdb_path, 'r') as f:
        pdb_content = f.read()
    print(f"Loaded PDB: {pdb_path}")

    print(f"Output: {OUTPUT_DIR}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    config = {
        "name": "quick_test",
        "contig": "70-90",
        "num_designs": 2,
        "cfg_scale": 2.0,
        "select_buried": {"X1": "all"},  # Bury Dy metal
        "use_cfg": True,
    }

    print(f"\nRunning: {config['name']}")
    pdbs = run_rfd3_design(pdb_content, config)

    if pdbs:
        print(f"\nGenerated {len(pdbs)} designs")
        for i, pdb in enumerate(pdbs):
            out_path = OUTPUT_DIR / f"quick_test_{i:03d}.pdb"
            with open(out_path, 'w') as f:
                f.write(pdb)
            print(f"  Saved: {out_path.name}")
        return [OUTPUT_DIR / f"quick_test_{i:03d}.pdb" for i in range(len(pdbs))]
    else:
        print("\nNo designs generated.")
        return []


def run_full_campaign():
    """Run the full design campaign with multiple configurations."""
    print("=" * 70)
    print("FULL CAMPAIGN - Dy-TriNOx Protein Design")
    print("=" * 70)

    # Read PDB with Dy (chain X) and TriNOx (chain L, UNL residue)
    pdb_path = INPUT_DIR / "Dy_TriNOx_split.pdb"
    if not pdb_path.exists():
        print(f"ERROR: PDB not found at {pdb_path}")
        return

    with open(pdb_path, 'r') as f:
        pdb_content = f.read()

    print(f"\nInput PDB: {pdb_path}")
    print(f"Output: {OUTPUT_DIR}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Configuration matrix
    # NOTE: Contig format "X,L,/0,LENGTH" includes metal+ligand chains, then protein
    # Dy on chain X (X1), TriNOx on chain L (L1)
    configs = [
        # Small proteins - direct coordination
        {
            "name": "r1a_small",
            "contig": "60-80",
            "num_designs": 5,
            "cfg_scale": 2.0,
            "select_buried": {"X1": "all"},  # Bury Dy metal
        },
        # Medium proteins - better burial
        {
            "name": "r1b_medium",
            "contig": "90-120",
            "num_designs": 5,
            "cfg_scale": 2.0,
            "select_buried": {"X1": "all"},  # Bury Dy metal
        },
        # Large proteins - full encapsulation
        {
            "name": "r1c_large",
            "contig": "130-160",
            "num_designs": 3,
            "cfg_scale": 2.5,
            "select_buried": {"X1": "all", "L1": "O1,O2,O3"},  # Bury Dy + O atoms
        },
    ]

    all_designs = []

    for config in configs:
        print(f"\n{'='*60}")
        print(f"Running: {config['name']}")
        print(f"  Protein length: {config['contig']}")
        print(f"  Designs: {config['num_designs']}")
        print(f"  CFG: {config['cfg_scale']}")
        print("=" * 60)

        pdbs = run_rfd3_design(pdb_content, config)

        if pdbs:
            print(f"  Generated {len(pdbs)} designs")
            for i, pdb in enumerate(pdbs):
                name = f"{config['name']}_{i:03d}"
                out_path = OUTPUT_DIR / f"{name}.pdb"
                with open(out_path, 'w') as f:
                    f.write(pdb)
                print(f"    Saved: {name}.pdb")
                all_designs.append(name)

    # Save summary
    summary = {
        "timestamp": datetime.now().isoformat(),
        "input_pdb": str(pdb_path),
        "total_designs": len(all_designs),
        "designs": all_designs,
        "configs": configs,
    }

    summary_path = OUTPUT_DIR / "run_summary.json"
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"\n{'='*70}")
    print(f"CAMPAIGN COMPLETE: {len(all_designs)} designs")
    print(f"Summary: {summary_path}")
    print("=" * 70)


def main():
    parser = argparse.ArgumentParser(
        description="Run Dy-TriNOx protein design campaign"
    )
    parser.add_argument(
        "--quick", "-q",
        action="store_true",
        help="Run quick test with 2 designs"
    )

    args = parser.parse_args()

    if args.quick:
        run_quick_test()
    else:
        run_full_campaign()


if __name__ == "__main__":
    main()
