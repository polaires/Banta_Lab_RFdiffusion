"""
Run RF3 validation on R7 sequences - corrected response parsing.
"""

import requests
from pathlib import Path

API_URL = "http://localhost:8000/runsync"
OUTPUT_DIR = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\experiments\ln_citrate_scaffold\outputs\round_07_scaffold")

# Best sequences from LigandMPNN (from r7_mpnn_sequences.fasta)
sequences = {
    "r7_design_1": "TLRLTITFAPGDLKVTFFDAETGEKLGTYVGRDAIIAANNELRDAGVWHTEVVARADAPEKAVRDSVPHRTLSIEQIAPNTIVAVVETADPAAFVAEETAEVAALGGSLTYEVL",
    "r7_design_4": "RLRVTITYAPGEKLVRFFDAETELLGTYVGMDAILAANADAAAGKWITEEVALADRRPEKAVRDSVPHRILSLERIAPNTRAVVETDDPAAFAAEETAEVAAMGGSMTWEVL",
    "r7_design_7": "TLRVTITWAPGEKLVTFFDAETKEKLGTYVGRDAILAAKNELTAAGKWHTEEVALADRKPEKAVRESVPHKILSLEQVAPDTVRAVIETDDPDALIAAETAAVAAMGGSMTAERL",
    "r7_design_8": "RLRLTITYAPGDLLVRFFNAETGELLGTFVGRDAILAANEELAAGVWHTEEVARADRPEDAVAETTPHILSLEQVAPNTVVAEVDTADPDALVAEETAAVAALGGSLTYERL",
}

print("RF3 VALIDATION - R7 Sequences (v2)")
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

        # Debug: print full response structure
        print(f"  Response status: {result.get('status')}")

        output = result.get("output", {})
        res = output.get("result", {})

        # RF3 response can have different structures
        # Check for predictions array first (newer format)
        predictions = res.get("predictions", [])
        if predictions and len(predictions) > 0:
            pred = predictions[0]
            plddt = pred.get("mean_plddt", pred.get("overall_plddt", 0))
            ptm = pred.get("ptm", 0)
            pdb_content = pred.get("pdb_string", pred.get("pdb_content", ""))
        else:
            # Try direct fields (older format)
            plddt = res.get("mean_plddt", res.get("plddt", res.get("overall_plddt", 0)))
            ptm = res.get("ptm", 0)
            pdb_content = res.get("pdb_string", res.get("pdb_content", ""))

        # Format pLDDT
        if isinstance(plddt, (int, float)):
            print(f"  pLDDT: {plddt:.3f}")
        else:
            print(f"  pLDDT: {plddt}")

        if isinstance(ptm, (int, float)) and ptm > 0:
            print(f"  pTM: {ptm:.3f}")

        results.append({
            "name": name,
            "seq": seq,
            "plddt": plddt if isinstance(plddt, (int, float)) else 0,
            "ptm": ptm if isinstance(ptm, (int, float)) else 0,
        })

        # Save PDB
        if pdb_content:
            pdb_path = OUTPUT_DIR / f"{name}_rf3.pdb"
            with open(pdb_path, 'w') as f:
                f.write(pdb_content)
            print(f"  Saved: {pdb_path.name}")
        else:
            print("  No PDB content in response")

    except Exception as e:
        print(f"  Error: {e}")
        import traceback
        traceback.print_exc()
        results.append({"name": name, "seq": seq, "plddt": 0, "ptm": 0})

# Summary
print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)

# Sort by pLDDT
sorted_results = sorted(results, key=lambda x: x['plddt'], reverse=True)

for r in sorted_results:
    print(f"{r['name']}: pLDDT={r['plddt']:.3f}, pTM={r['ptm']:.3f}")

# Best candidate
best = sorted_results[0]
if best['plddt'] > 0.7:
    print(f"\n*** BEST: {best['name']} with pLDDT {best['plddt']:.3f} ***")
    print(f"Sequence ({len(best['seq'])} aa): {best['seq']}")

    # Save best as FASTA for AF3 submission
    best_fasta = OUTPUT_DIR / "BEST_R7_for_AF3.fasta"
    with open(best_fasta, 'w') as f:
        f.write(f">{best['name']}_pLDDT_{best['plddt']:.2f}\n")
        f.write(f"{best['seq']}\n")
    print(f"Saved: {best_fasta}")
