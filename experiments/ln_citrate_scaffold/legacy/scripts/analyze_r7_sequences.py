"""
Analyze R7 LigandMPNN sequences and run RF3 validation.
"""

import requests
from pathlib import Path

API_URL = "http://localhost:8000/runsync"
BASE_DIR = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\experiments\ln_citrate_scaffold")
OUTPUT_DIR = BASE_DIR / "outputs" / "round_07_scaffold"

# Read sequences
fasta_path = OUTPUT_DIR / "r7_mpnn_sequences.fasta"
sequences = {}
current_name = None
with open(fasta_path) as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            current_name = line[1:]
            sequences[current_name] = ""
        elif current_name:
            sequences[current_name] += line

print("=" * 60)
print("R7 SEQUENCE ANALYSIS")
print("=" * 60)

for name, seq in sequences.items():
    total = len(seq)
    ala = seq.count('A')
    glu = seq.count('E')
    asp = seq.count('D')
    arg = seq.count('R')
    lys = seq.count('K')
    ser = seq.count('S')
    thr = seq.count('T')

    print(f"\n{name}:")
    print(f"  Length: {total}")
    print(f"  Ala: {ala} ({100*ala/total:.1f}%)")
    print(f"  E+D (carboxylates): {glu+asp} ({100*(glu+asp)/total:.1f}%)")
    print(f"  R+K (positively charged): {arg+lys} ({100*(arg+lys)/total:.1f}%)")
    print(f"  S+T (H-bond donors): {ser+thr} ({100*(ser+thr)/total:.1f}%)")

# Run RF3 on best candidates (first few)
print("\n" + "=" * 60)
print("RF3 VALIDATION")
print("=" * 60)

rf3_results = []

for name, seq in list(sequences.items())[:4]:  # First 4
    print(f"\nValidating {name}...")

    payload = {
        "input": {
            "task": "rf3",
            "sequence": seq,
            "name": name
        }
    }

    try:
        response = requests.post(API_URL, json=payload, timeout=180)
        result = response.json()

        output = result.get("output", {})
        res = output.get("result", {})

        plddt = res.get("plddt", 0)
        rmsd = res.get("rmsd", "N/A")

        print(f"  pLDDT: {plddt:.2f}" if isinstance(plddt, float) else f"  pLDDT: {plddt}")
        print(f"  RMSD: {rmsd}")

        rf3_results.append({
            "name": name,
            "seq": seq,
            "plddt": plddt,
            "rmsd": rmsd,
        })

        # Save PDB if available
        pdb_content = res.get("pdb_content")
        if pdb_content:
            pdb_path = OUTPUT_DIR / f"{name}_rf3.pdb"
            with open(pdb_path, 'w') as f:
                f.write(pdb_content)
            print(f"  Saved: {pdb_path.name}")

    except Exception as e:
        print(f"  Error: {e}")

# Summary
print("\n" + "=" * 60)
print("RF3 SUMMARY")
print("=" * 60)

for r in rf3_results:
    print(f"{r['name']}: pLDDT={r['plddt']}, RMSD={r['rmsd']}")

# Find best
best = min(rf3_results, key=lambda x: float(x['rmsd']) if isinstance(x['rmsd'], (int, float)) else 999)
print(f"\nBest: {best['name']} with RMSD={best['rmsd']}")
