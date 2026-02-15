"""
Run a single R7 scaffolding design for quick testing.
"""

import requests
import json
from pathlib import Path

API_URL = "http://localhost:8000/runsync"
BASE_DIR = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\experiments\ln_citrate_scaffold")
INPUT_PDB = BASE_DIR / "inputs" / "citrate_ln_only.pdb"
OUTPUT_DIR = BASE_DIR / "outputs" / "round_07_scaffold"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def load_pdb_content(pdb_path: Path) -> str:
    with open(pdb_path, 'r') as f:
        return f.read()

pdb_content = load_pdb_content(INPUT_PDB)

# Single design config
config = {
    "task": "rfd3",
    "pdb_content": pdb_content,
    "contig": "110-130",
    "ligand": "CIT,TB",
    "num_designs": 1,  # Single design
    "step_scale": 1.5,
    "num_timesteps": 200,
    "select_fixed_atoms": {"X1": "all", "L1": "all"},
    "select_buried": {"X1": "all"},
    "select_hbond_acceptor": {"L1": "O1,O2,O3,O4,O5,O6"},
    "select_hbond_donor": {"L1": "O7"},
    "infer_ori_strategy": "com",
    "use_classifier_free_guidance": True,
    "cfg_scale": 2.0,
}

print("Sending single design request...")
print("This may take 2-3 minutes...")

try:
    response = requests.post(
        API_URL,
        json={"input": config},
        timeout=600
    )
    result = response.json()

    print(f"\nStatus: {result.get('status')}")

    # Handle both structures
    output = result.get("output", {})
    inner_status = output.get("status")
    print(f"Inner status: {inner_status}")

    res = output.get("result", {})
    designs = res.get("designs", [])

    print(f"Number of designs: {len(designs)}")

    for i, design in enumerate(designs):
        content = design.get("content", "")
        filename = design.get("filename", f"design_{i}.pdb")

        # Save PDB
        pdb_path = OUTPUT_DIR / f"r7_single_{i:03d}.pdb"
        with open(pdb_path, 'w') as f:
            f.write(content)
        print(f"Saved: {pdb_path}")

        # Count ligand atoms
        cit_count = content.count("CIT L")
        tb_count = content.count("TB  X") or content.count("TB X")
        protein_atoms = len([l for l in content.split('\n') if l.startswith('ATOM')])
        print(f"  Protein atoms: {protein_atoms}")
        print(f"  Citrate atoms: {cit_count}")
        print(f"  TB atoms: {tb_count}")

except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
