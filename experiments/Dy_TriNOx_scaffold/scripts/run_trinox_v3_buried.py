"""
Dy-TriNOx scaffold design v3 - With partial burial of TriNOx-Dy complex.

Key changes from v2:
- Use select_partially_buried on TriNOx to create protein pocket
- Larger scaffolds (60-100 aa) to wrap around the complex
- Still use fixed Asp motif for coordination
"""

import requests
import json
from pathlib import Path
from datetime import datetime
import math

API_URL = "http://localhost:8000/runsync"
SCRIPT_DIR = Path(__file__).parent
INPUT_DIR = SCRIPT_DIR.parent / "inputs"
OUTPUT_DIR = SCRIPT_DIR.parent / "outputs_v3"

OUTPUT_DIR.mkdir(exist_ok=True)


def run_rfd3_design(pdb_content: str, config: dict, timeout: int = 600) -> list:
    """Run RFD3 with burial conditioning."""

    contig = config.get("contig")

    payload = {
        "input": {
            "task": "rfd3",
            "pdb_content": pdb_content,
            "contig": contig,
            "ligand": "DY,UNL",
            "num_designs": config.get("num_designs", 5),
            "select_fixed_atoms": {
                "X1": "all",  # Dy - fixed
                "L1": "all",  # TriNOx - fixed coordinates
                "M1": "all",  # Coordinating Asp - fixed
            },
            # CRITICAL: Partially bury the TriNOx to create protein pocket
            "select_partially_buried": {
                "L1": "all",  # TriNOx should be in a pocket
            },
            # Also bury the Dy
            "select_buried": {
                "X1": "all",  # Dy should be buried
            },
            "use_classifier_free_guidance": True,
            "cfg_scale": config.get("cfg_scale", 2.5),
        }
    }

    print(f"  Contig: {contig}")
    print(f"  Burial: TriNOx=partially_buried, Dy=buried")

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
    """Analyze coordination and burial in a design."""
    atoms = []
    for line in pdb_content.split('\n'):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms.append({
                    "name": atom_name,
                    "res_name": res_name,
                    "chain": chain,
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

    # Count protein atoms
    protein_atoms = [a for a in atoms if a["res_name"] in
                    ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                     "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
                     "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]]

    # Count protein residues
    protein_residues = set()
    for a in atoms:
        if a["res_name"] in ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                             "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
                             "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]:
            protein_residues.add(a["res_name"])

    # Find TriNOx atoms
    trinox_atoms = [a for a in atoms if a["res_name"] == "UNL"]

    # Check protein atoms near TriNOx (within 5A) - indicates burial
    trinox_nearby = 0
    for ta in trinox_atoms:
        for pa in protein_atoms:
            d = math.sqrt((ta["x"]-pa["x"])**2 + (ta["y"]-pa["y"])**2 + (ta["z"]-pa["z"])**2)
            if d <= 5.0:
                trinox_nearby += 1
                break  # Count each trinox atom only once

    # Find coordinating atoms to Dy
    coordinating = []
    for a in atoms:
        if a["res_name"] in ["ASP", "GLU", "HIS", "CYS"]:
            if a["name"] in ["OD1", "OD2", "OE1", "OE2", "ND1", "NE2", "SG"]:
                d = math.sqrt((a["x"]-dy["x"])**2 + (a["y"]-dy["y"])**2 + (a["z"]-dy["z"])**2)
                if d <= 3.0:
                    coordinating.append({"res": a["res_name"], "atom": a["name"], "dist": round(d, 2)})

    return {
        "num_residues": len(protein_residues),
        "num_protein_atoms": len(protein_atoms),
        "trinox_atoms_near_protein": trinox_nearby,
        "total_trinox_atoms": len(trinox_atoms),
        "burial_ratio": round(trinox_nearby / max(len(trinox_atoms), 1), 2),
        "num_coordinating": len(coordinating),
        "coordinating": coordinating,
    }


def main():
    print("=" * 70)
    print("Dy-TriNOx SCAFFOLD DESIGN v3")
    print("With PARTIAL BURIAL of TriNOx-Dy complex")
    print("=" * 70)

    # Load input PDB with Asp motif
    input_pdb = INPUT_DIR / "Dy_TriNOx_with_coord_Asp.pdb"
    with open(input_pdb, 'r') as f:
        pdb_content = f.read()

    # Design configs with LARGER scaffolds and burial
    configs = [
        # Medium scaffolds with burial
        {
            "name": "v3_medium",
            "contig": "X1,L1,20-30,M1,20-30",  # 40-60 new residues
            "num_designs": 5,
            "cfg_scale": 2.5,
        },
        # Larger scaffolds
        {
            "name": "v3_large",
            "contig": "X1,L1,30-50,M1,30-50",  # 60-100 new residues
            "num_designs": 5,
            "cfg_scale": 2.5,
        },
        # Asymmetric - more N-terminal
        {
            "name": "v3_nterm",
            "contig": "X1,L1,40-60,M1,15-25",  # 55-85 new residues
            "num_designs": 5,
            "cfg_scale": 2.5,
        },
        # Asymmetric - more C-terminal
        {
            "name": "v3_cterm",
            "contig": "X1,L1,15-25,M1,40-60",  # 55-85 new residues
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
            print(f"  {design_name}: {analysis.get('num_protein_atoms', 0)} atoms, "
                  f"burial={analysis.get('burial_ratio', 0)}, "
                  f"coord={analysis.get('num_coordinating', 0)} - {coord_info}")

    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"Total designs: {len(all_designs)}")

    designs_with_coord = sum(1 for r in results if r.get("num_coordinating", 0) > 0)
    designs_with_burial = sum(1 for r in results if r.get("burial_ratio", 0) > 0.3)
    print(f"Designs with coordination (<3A): {designs_with_coord}")
    print(f"Designs with good burial (>30%): {designs_with_burial}")

    # Save summary
    summary = {
        "timestamp": datetime.now().isoformat(),
        "input_pdb": str(input_pdb),
        "total_designs": len(all_designs),
        "designs_with_coordination": designs_with_coord,
        "designs_with_burial": designs_with_burial,
        "configs": configs,
        "results": results,
    }

    summary_file = OUTPUT_DIR / "run_summary_v3.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"\nSummary saved to {summary_file}")


if __name__ == "__main__":
    main()
