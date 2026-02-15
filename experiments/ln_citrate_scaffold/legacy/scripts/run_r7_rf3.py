"""
Run RF3 validation on R7 sequences.
"""

import requests
from pathlib import Path

API_URL = "http://localhost:8000/runsync"
OUTPUT_DIR = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\experiments\ln_citrate_scaffold\outputs\round_07_scaffold")

# Best sequences from analysis (lower Ala, good E+D/S+T)
sequences = {
    "r7_design_1": "TLRLTITFAPGDLKVTFFDAETGEKLGTYVGRDAIIAANNELRDAGVWHTEVVARADAPEKAVRDSVPHRTLSIEQIAPNTIVAVVETADPAAFVAEETAEVAALGGSLTYEVL",
    "r7_design_4": "RLRVTITYAPGEKLVRFFDAETELLGTYVGMDAILAANADAAAGKWITEEVALADRRPEKAVRDSVPHRILSLERIAPNTRAVVETDDPAAFAAEETAEVAAMGGSMTWEVL",
    "r7_design_7": "TLRVTITWAPGEKLVTFFDAETKEKLGTYVGRDAILAAKNELTAAGKWHTEEVALADRKPEKAVRESVPHKILSLEQVAPDTVRAVIETDDPDALIAAETAAVAAMGGSMTAERL",
    "r7_design_8": "RLRLTITYAPGDLLVRFFNAETGELLGTFVGRDAILAANEELAAGVWHTEEVARADRPEDAVAETTPHILSLEQVAPNTVVAEVDTADPDALVAEETAAVAALGGSLTYERL",
}

print("RF3 VALIDATION - R7 Sequences")
print("=" * 60)

results = []

for name, seq in sequences.items():
    print(f"\nValidating {name} ({len(seq)} aa)...")

    payload = {
        "input": {
            "task": "rf3",
            "sequence": seq,
            "name": name
        }
    }

    try:
        response = requests.post(API_URL, json=payload, timeout=300)
        result = response.json()

        output = result.get("output", {})
        res = output.get("result", {})

        plddt = res.get("plddt", res.get("mean_plddt", 0))
        rmsd = res.get("rmsd", res.get("backbone_rmsd", "N/A"))

        print(f"  pLDDT: {plddt:.2f}" if isinstance(plddt, (int, float)) else f"  pLDDT: {plddt}")
        print(f"  RMSD to backbone: {rmsd}")

        results.append({
            "name": name,
            "seq": seq,
            "plddt": plddt,
            "rmsd": rmsd
        })

        # Save PDB
        pdb_content = res.get("pdb_content")
        if pdb_content:
            pdb_path = OUTPUT_DIR / f"{name}_rf3.pdb"
            with open(pdb_path, 'w') as f:
                f.write(pdb_content)
            print(f"  Saved: {pdb_path.name}")

    except Exception as e:
        print(f"  Error: {e}")
        results.append({"name": name, "seq": seq, "plddt": 0, "rmsd": 999})

# Summary
print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)

# Sort by RMSD
sorted_results = sorted(results, key=lambda x: float(x['rmsd']) if isinstance(x['rmsd'], (int, float)) else 999)

for r in sorted_results:
    rmsd_str = f"{r['rmsd']:.2f}Å" if isinstance(r['rmsd'], (int, float)) else str(r['rmsd'])
    print(f"{r['name']}: pLDDT={r['plddt']}, RMSD={rmsd_str}")

# Best candidate
best = sorted_results[0]
if isinstance(best['rmsd'], (int, float)) and best['rmsd'] < 2.0:
    print(f"\n*** EXCELLENT: {best['name']} has RMSD {best['rmsd']:.2f}Å ***")
    print(f"Sequence: {best['seq']}")

    # Save best as FASTA for AF3 submission
    best_fasta = OUTPUT_DIR / "BEST_R7.fasta"
    with open(best_fasta, 'w') as f:
        f.write(f">{best['name']}_RMSD_{best['rmsd']:.2f}A\n")
        f.write(f"{best['seq']}\n")
    print(f"Saved: {best_fasta}")
