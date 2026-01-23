"""
Run RF3 validation on R7 sequences - correct response parsing.
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

print("RF3 VALIDATION - R7 Sequences (v3)")
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

        # RF3 returns predictions array with metrics and CIF content
        predictions = res.get("predictions", [])

        plddt = 0
        ptm = 0
        pae = 999
        cif_content = ""

        if predictions and len(predictions) > 0:
            pred = predictions[0]
            # Metrics are directly in predictions[0]
            plddt = pred.get("mean_plddt", pred.get("overall_plddt", 0))
            ptm = pred.get("ptm", 0)
            pae = pred.get("overall_pae", 999)
            cif_content = pred.get("content", "")
            filename = pred.get("filename", "")

        print(f"  mean_pLDDT: {plddt:.4f}" if isinstance(plddt, (int, float)) else f"  pLDDT: {plddt}")
        print(f"  pTM: {ptm:.4f}" if isinstance(ptm, (int, float)) else f"  pTM: {ptm}")
        print(f"  PAE: {pae:.2f}" if isinstance(pae, (int, float)) else f"  PAE: {pae}")

        results.append({
            "name": name,
            "seq": seq,
            "plddt": plddt if isinstance(plddt, (int, float)) else 0,
            "ptm": ptm if isinstance(ptm, (int, float)) else 0,
            "pae": pae if isinstance(pae, (int, float)) else 999,
        })

        # Save CIF (RF3 outputs mmCIF format)
        if cif_content:
            cif_path = OUTPUT_DIR / f"{name}_rf3.cif"
            with open(cif_path, 'w') as f:
                f.write(cif_content)
            print(f"  Saved: {cif_path.name}")

    except Exception as e:
        print(f"  Error: {e}")
        import traceback
        traceback.print_exc()
        results.append({"name": name, "seq": seq, "plddt": 0, "ptm": 0, "pae": 999})

# Summary
print("\n" + "=" * 60)
print("SUMMARY - Sorted by pLDDT")
print("=" * 60)

# Sort by pLDDT (descending)
sorted_results = sorted(results, key=lambda x: x['plddt'], reverse=True)

print(f"{'Design':<15} {'pLDDT':>8} {'pTM':>8} {'PAE':>8}")
print("-" * 45)
for r in sorted_results:
    print(f"{r['name']:<15} {r['plddt']:>8.4f} {r['ptm']:>8.4f} {r['pae']:>8.2f}")

# Best candidate
best = sorted_results[0]
if best['plddt'] > 0.7:
    print(f"\n*** BEST: {best['name']} ***")
    print(f"    pLDDT={best['plddt']:.4f}, pTM={best['ptm']:.4f}, PAE={best['pae']:.2f}")
    print(f"    Sequence ({len(best['seq'])} aa):")
    print(f"    {best['seq']}")

    # Save best as FASTA for AF3 submission
    best_fasta = OUTPUT_DIR / "BEST_R7_for_AF3.fasta"
    with open(best_fasta, 'w') as f:
        f.write(f">{best['name']}_pLDDT{best['plddt']:.2f}_pTM{best['ptm']:.2f}\n")
        f.write(f"{best['seq']}\n")
    print(f"    Saved: {best_fasta}")

# Selection criteria check
print("\n" + "=" * 60)
print("SELECTION CRITERIA")
print("=" * 60)
print("Passing thresholds: pLDDT > 0.80, pTM > 0.80, PAE < 3.0")
print("-" * 60)

passing = [r for r in results if r['plddt'] > 0.80 and r['ptm'] > 0.80 and r['pae'] < 3.0]
if passing:
    print(f"✓ {len(passing)} designs pass all thresholds:")
    for r in passing:
        print(f"  - {r['name']}: pLDDT={r['plddt']:.3f}, pTM={r['ptm']:.3f}, PAE={r['pae']:.2f}")
else:
    print("No designs pass all strict thresholds.")
    # Relaxed criteria
    relaxed = [r for r in results if r['plddt'] > 0.75]
    if relaxed:
        print(f"Relaxed (pLDDT > 0.75): {len(relaxed)} designs")
        for r in relaxed:
            print(f"  - {r['name']}: pLDDT={r['plddt']:.3f}, pTM={r['ptm']:.3f}")
