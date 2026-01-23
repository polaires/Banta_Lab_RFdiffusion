"""
Round 7b - Improved scaffolding with more designs and higher CFG.

Improvements over R7:
1. Generate 5 scaffold backbones (not just 1)
2. Try 2 scaffold size ranges (110-130, 130-150)
3. Higher CFG scale (2.5) for stronger conditioning
4. More LigandMPNN samples (16 per backbone)
"""

import requests
import json
from pathlib import Path

API_URL = "http://localhost:8000/runsync"
INPUT_DIR = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\experiments\ln_citrate_scaffold\inputs")
OUTPUT_DIR = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\experiments\ln_citrate_scaffold\outputs\round_07b_improved")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Read input PDB
pdb_path = INPUT_DIR / "tb_citrate_motif_scaffold.pdb"
with open(pdb_path, 'r') as f:
    pdb_content = f.read()

print("=" * 70)
print("ROUND 7b - IMPROVED SCAFFOLDING")
print("=" * 70)

# Configurations to test
configs = [
    {"name": "r7b_small", "contig": "110-130", "num_designs": 3, "cfg_scale": 2.5},
    {"name": "r7b_medium", "contig": "130-150", "num_designs": 2, "cfg_scale": 2.5},
]

all_designs = []

for cfg in configs:
    print(f"\n{'='*60}")
    print(f"Running: {cfg['name']}")
    print(f"Contig: {cfg['contig']}, Designs: {cfg['num_designs']}, CFG: {cfg['cfg_scale']}")
    print("=" * 60)

    payload = {
        "input": {
            "task": "rfd3",
            "pdb_content": pdb_content,
            "contig": cfg["contig"],
            "ligand": "CIT,TB",
            "num_designs": cfg["num_designs"],
            # Fix TB and citrate in place
            "select_fixed_atoms": {
                "X1": "all",  # TB
                "L1": "all",  # Citrate
            },
            # Bury TB for pocket formation
            "select_buried": {"X1": "all"},
            # H-bond conditioning for citrate oxygens
            "select_hbond_acceptor": {"L1": "O1,O2,O3,O4,O5,O6"},
            "select_hbond_donor": {"L1": "O7"},
            # Stronger CFG for better conditioning
            "use_classifier_free_guidance": True,
            "cfg_scale": cfg["cfg_scale"],
            "infer_ori_strategy": "com",
        }
    }

    try:
        response = requests.post(API_URL, json=payload, timeout=1200)
        result = response.json()
        status = result.get("status", result.get("output", {}).get("status", "unknown"))
        print(f"Status: {status}")

        if status.upper() == "COMPLETED":
            output = result.get("output", {})
            res = output.get("result", {})
            # API can return 'pdbs' or 'designs'
            pdbs = res.get("pdbs", res.get("designs", []))

            print(f"Generated {len(pdbs)} designs")

            for i, pdb_data in enumerate(pdbs):
                design_name = f"{cfg['name']}_{i:03d}"
                pdb_out = OUTPUT_DIR / f"{design_name}.pdb"
                with open(pdb_out, 'w') as f:
                    f.write(pdb_data.get("content", pdb_data) if isinstance(pdb_data, dict) else pdb_data)
                print(f"  Saved: {design_name}.pdb")
                all_designs.append(design_name)
        else:
            print(f"Failed: {result}")
    except Exception as e:
        print(f"Exception: {e}")

print(f"\n{'='*70}")
print(f"SCAFFOLDING COMPLETE: {len(all_designs)} designs")
print("=" * 70)

# Now run LigandMPNN on all designs
if all_designs:
    print("\n" + "=" * 70)
    print("RUNNING LIGANDMPNN ON ALL DESIGNS")
    print("=" * 70)

    all_sequences = []

    for design_name in all_designs:
        pdb_file = OUTPUT_DIR / f"{design_name}.pdb"
        with open(pdb_file, 'r') as f:
            design_pdb = f.read()

        print(f"\nRunning LigandMPNN on {design_name}...")

        mpnn_payload = {
            "input": {
                "task": "ligandmpnn",
                "pdb_content": design_pdb,
                "num_seqs": 4,  # 4 sequences per backbone
                "ligand_mpnn_use_atom_context": 1,
                "temperature": 0.2,
            }
        }

        try:
            response = requests.post(API_URL, json=mpnn_payload, timeout=600)
            result = response.json()

            output = result.get("output", {})
            res = output.get("result", {})
            sequences = res.get("sequences", [])

            if sequences:
                for seq_data in sequences:
                    if isinstance(seq_data, dict) and "content" in seq_data:
                        # Parse FASTA content
                        content = seq_data["content"]
                        lines = content.strip().split("\n")
                        for j in range(0, len(lines), 2):
                            if j + 1 < len(lines) and lines[j].startswith(">"):
                                header = lines[j]
                                seq = lines[j + 1]
                                new_name = f"{design_name}_seq{len(all_sequences) + 1}"
                                all_sequences.append({"name": new_name, "seq": seq})
                                print(f"  {new_name}: {len(seq)} aa")
                    elif isinstance(seq_data, str):
                        all_sequences.append({"name": f"{design_name}_seq", "seq": seq_data})
        except Exception as e:
            print(f"  LigandMPNN error: {e}")

    # Save all sequences
    if all_sequences:
        fasta_out = OUTPUT_DIR / "r7b_all_sequences.fasta"
        with open(fasta_out, 'w') as f:
            for s in all_sequences:
                f.write(f">{s['name']}\n{s['seq']}\n")
        print(f"\nSaved {len(all_sequences)} sequences to {fasta_out}")

print("\n" + "=" * 70)
print("ROUND 7b COMPLETE")
print("=" * 70)
