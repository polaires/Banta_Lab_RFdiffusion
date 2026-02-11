"""
Dy-TriNOx scaffold design v2 - Motif scaffolding approach.

Uses a fixed Asp residue at the coordination site, then scaffolds around it.
This guarantees proper metal coordination (~2.4A Asp-Dy distance).
"""

import requests
import json
from pathlib import Path
from datetime import datetime
import math

API_URL = "http://localhost:8000/runsync"
SCRIPT_DIR = Path(__file__).parent
INPUT_DIR = SCRIPT_DIR.parent / "inputs"
OUTPUT_DIR = SCRIPT_DIR.parent / "outputs_v2"

OUTPUT_DIR.mkdir(exist_ok=True)


def run_rfd3_design(pdb_content: str, config: dict, timeout: int = 600) -> list:
    """Run RFD3 motif scaffolding design."""

    contig = config.get("contig")

    payload = {
        "input": {
            "task": "rfd3",
            "pdb_content": pdb_content,
            "contig": contig,
            "ligand": "DY,UNL",
            "num_designs": config.get("num_designs", 5),
            "select_fixed_atoms": {
                "X1": "all",  # Dy
                "L1": "all",  # TriNOx
                "M1": "all",  # Coordinating Asp
            },
            "select_buried": config.get("select_buried", {"X1": "all"}),
            "use_classifier_free_guidance": True,
            "cfg_scale": config.get("cfg_scale", 2.5),
        }
    }

    print(f"  Contig: {contig}")

    try:
        response = requests.post(API_URL, json=payload, timeout=timeout)
        result = response.json()

        if result.get("status") != "COMPLETED":
            error = result.get("error", "Unknown error")
            print(f"  ERROR: {error}")
            return []

        output = result.get("output", {})
        res = output.get("result", {})
        designs = res.get("designs", [])

        return [d["content"] for d in designs]

    except Exception as e:
        print(f"  Error: {e}")
        return []


def analyze_design(pdb_content: str) -> dict:
    """Analyze coordination in a design."""
    atoms = []
    for line in pdb_content.split('\n'):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms.append({
                    "name": atom_name,
                    "res_name": res_name,
                    "x": x, "y": y, "z": z
                })
            except:
                pass

    # Find Dy
    dy = None
    for a in atoms:
        if a["res_name"] == "DY" or a["name"] == "DY":
            dy = a
            break

    if not dy:
        return {"error": "No Dy"}

    # Find coordinating atoms
    protein_atoms = [a for a in atoms if a["res_name"] not in ["DY", "UNL"]]
    coordinating = []
    for a in protein_atoms:
        d = math.sqrt((a["x"]-dy["x"])**2 + (a["y"]-dy["y"])**2 + (a["z"]-dy["z"])**2)
        if d <= 3.0:
            coordinating.append({"res": a["res_name"], "atom": a["name"], "dist": round(d, 2)})

    return {
        "num_coordinating": len(coordinating),
        "coordinating": coordinating,
        "num_protein_atoms": len(protein_atoms)
    }


def main():
    print("=" * 70)
    print("Dy-TriNOx SCAFFOLD DESIGN v2")
    print("Motif scaffolding with fixed Asp coordination")
    print("=" * 70)

    # Load input PDB with Asp motif
    input_pdb = INPUT_DIR / "Dy_TriNOx_with_coord_Asp.pdb"
    with open(input_pdb, 'r') as f:
        pdb_content = f.read()

    # Design configs with different scaffold lengths
    configs = [
        # Small scaffolds - tight binding pocket
        {
            "name": "v2_small",
            "contig": "X1,L1,3-8,M1,3-8",
            "num_designs": 5,
            "cfg_scale": 2.5,
        },
        # Medium scaffolds
        {
            "name": "v2_medium",
            "contig": "X1,L1,8-15,M1,8-15",
            "num_designs": 5,
            "cfg_scale": 2.5,
        },
        # Longer N-terminal linker
        {
            "name": "v2_nterm",
            "contig": "X1,L1,15-25,M1,3-8",
            "num_designs": 5,
            "cfg_scale": 2.5,
        },
        # Longer C-terminal linker
        {
            "name": "v2_cterm",
            "contig": "X1,L1,3-8,M1,15-25",
            "num_designs": 5,
            "cfg_scale": 2.5,
        },
    ]

    all_designs = []
    results = []

    for config in configs:
        print(f"\n{'='*70}")
        print(f"Config: {config['name']}")
        print(f"{'='*70}")

        designs = run_rfd3_design(pdb_content, config)

        print(f"Generated {len(designs)} designs")

        for i, design in enumerate(designs):
            design_name = f"{config['name']}_{i:03d}"
            design_file = OUTPUT_DIR / f"{design_name}.pdb"

            with open(design_file, 'w') as f:
                f.write(design)

            # Analyze
            analysis = analyze_design(design)

            all_designs.append(design_name)
            results.append({
                "name": design_name,
                "config": config["name"],
                **analysis
            })

            coord_info = ", ".join([f"{c['res']} {c['atom']}:{c['dist']}A" for c in analysis.get("coordinating", [])])
            print(f"  {design_name}: {analysis.get('num_coordinating', 0)} coord atoms - {coord_info}")

    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"Total designs: {len(all_designs)}")

    designs_with_coord = sum(1 for r in results if r.get("num_coordinating", 0) > 0)
    print(f"Designs with coordination (<3A): {designs_with_coord}")

    # Save summary
    summary = {
        "timestamp": datetime.now().isoformat(),
        "input_pdb": str(input_pdb),
        "total_designs": len(all_designs),
        "designs_with_coordination": designs_with_coord,
        "configs": configs,
        "results": results,
    }

    summary_file = OUTPUT_DIR / "run_summary_v2.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"\nSummary saved to {summary_file}")


if __name__ == "__main__":
    main()
