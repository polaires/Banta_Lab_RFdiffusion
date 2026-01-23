"""
Run LigandMPNN and RF3 validation on R7b scaffolds.
"""

import requests
from pathlib import Path

API_URL = "http://localhost:8000/runsync"
OUTPUT_DIR = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\experiments\ln_citrate_scaffold\outputs\round_07b_improved")

# Find all R7b PDBs
pdb_files = sorted(OUTPUT_DIR.glob("r7b_*.pdb"))
print(f"Found {len(pdb_files)} R7b scaffolds")

all_sequences = []

# Run LigandMPNN on each scaffold
print("\n" + "=" * 60)
print("LIGANDMPNN SEQUENCE DESIGN")
print("=" * 60)

for pdb_file in pdb_files:
    print(f"\n{pdb_file.stem}...")

    with open(pdb_file, 'r') as f:
        pdb_content = f.read()

    payload = {
        "input": {
            "task": "mpnn",
            "pdb_content": pdb_content,
            "num_seqs": 2,  # 2 sequences per backbone
            "ligand_mpnn_use_atom_context": 1,
            "temperature": 0.2,
        }
    }

    try:
        response = requests.post(API_URL, json=payload, timeout=300)
        result = response.json()

        output = result.get("output", {})
        res = output.get("result", {})
        sequences = res.get("sequences", [])

        if sequences:
            for seq_data in sequences:
                if isinstance(seq_data, dict) and "content" in seq_data:
                    content = seq_data["content"]
                    lines = content.strip().split("\n")
                    for j in range(0, len(lines), 2):
                        if j + 1 < len(lines) and lines[j].startswith(">"):
                            seq = lines[j + 1]
                            name = f"{pdb_file.stem}_seq{len(all_sequences) + 1}"
                            all_sequences.append({"name": name, "seq": seq, "backbone": pdb_file.stem})
                            print(f"  {name}: {len(seq)} aa")
    except Exception as e:
        print(f"  Error: {e}")

# Save all sequences
if all_sequences:
    fasta_out = OUTPUT_DIR / "r7b_all_sequences.fasta"
    with open(fasta_out, 'w') as f:
        for s in all_sequences:
            f.write(f">{s['name']}\n{s['seq']}\n")
    print(f"\nSaved {len(all_sequences)} sequences to {fasta_out}")

# Run RF3 validation
print("\n" + "=" * 60)
print("RF3 VALIDATION")
print("=" * 60)

results = []

for s in all_sequences[:8]:  # Validate first 8
    print(f"\nValidating {s['name']}...")

    payload = {
        "input": {
            "task": "rf3",
            "sequence": s["seq"],
            "name": s["name"]
        }
    }

    try:
        response = requests.post(API_URL, json=payload, timeout=300)
        result = response.json()

        output = result.get("output", {})
        res = output.get("result", {})
        confidences = res.get("confidences", {})

        plddt = confidences.get("mean_plddt", 0)
        ptm = confidences.get("ptm", 0)
        pae = confidences.get("overall_pae", 999)

        print(f"  pLDDT: {plddt:.4f}, pTM: {ptm:.4f}, PAE: {pae:.2f}")

        results.append({
            "name": s["name"],
            "seq": s["seq"],
            "backbone": s["backbone"],
            "plddt": plddt,
            "ptm": ptm,
            "pae": pae,
        })
    except Exception as e:
        print(f"  Error: {e}")

# Summary
print("\n" + "=" * 60)
print("SUMMARY - Sorted by pLDDT")
print("=" * 60)

sorted_results = sorted(results, key=lambda x: x['plddt'], reverse=True)

print(f"{'Design':<25} {'Backbone':<15} {'pLDDT':>8} {'pTM':>8} {'PAE':>8}")
print("-" * 70)
for r in sorted_results:
    print(f"{r['name']:<25} {r['backbone']:<15} {r['plddt']:>8.4f} {r['ptm']:>8.4f} {r['pae']:>8.2f}")

# Compare with R7 best
print("\n" + "=" * 60)
print("COMPARISON WITH R7 BEST (r7_design_1)")
print("=" * 60)
print(f"R7 best: pLDDT=0.865, pTM=0.891, PAE=2.34")

if sorted_results:
    best = sorted_results[0]
    print(f"R7b best: pLDDT={best['plddt']:.3f}, pTM={best['ptm']:.3f}, PAE={best['pae']:.2f}")

    if best['plddt'] > 0.865:
        print(f"\n*** R7b IMPROVEMENT: {best['name']} beats R7! ***")
    elif best['plddt'] > 0.80:
        print(f"\n*** R7b {best['name']} is viable (pLDDT > 0.80) ***")
