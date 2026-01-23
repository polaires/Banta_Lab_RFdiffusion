"""
Run RF3 validation on R7 sequences - FINAL with correct response parsing.
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

print("RF3 VALIDATION - R7 Sequences (FINAL)")
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

        # RF3 response structure:
        # - predictions: list of {filename, content} for structure files
        # - confidences: dict with all metrics
        predictions = res.get("predictions", [])
        confidences = res.get("confidences", {})

        # Extract metrics from confidences dict
        plddt = confidences.get("mean_plddt", confidences.get("overall_plddt", 0))
        ptm = confidences.get("ptm", 0)
        pae = confidences.get("overall_pae", 999)

        print(f"  mean_pLDDT: {plddt:.4f}")
        print(f"  pTM: {ptm:.4f}")
        print(f"  PAE: {pae:.2f}")

        results.append({
            "name": name,
            "seq": seq,
            "plddt": plddt,
            "ptm": ptm,
            "pae": pae,
        })

        # Save CIF file(s)
        for pred in predictions:
            filename = pred.get("filename", "")
            content = pred.get("content", "")
            if content and filename.endswith(".cif"):
                cif_path = OUTPUT_DIR / f"{name}_rf3.cif"
                with open(cif_path, 'w') as f:
                    f.write(content)
                print(f"  Saved: {cif_path.name}")
                break  # Just save first CIF

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
    print(f"[PASS] {len(passing)} designs pass all thresholds:")
    for r in passing:
        print(f"  - {r['name']}: pLDDT={r['plddt']:.3f}, pTM={r['ptm']:.3f}, PAE={r['pae']:.2f}")
else:
    print("No designs pass all strict thresholds.")
    relaxed = [r for r in results if r['plddt'] > 0.75 and r['ptm'] > 0.75]
    if relaxed:
        print(f"Relaxed (pLDDT > 0.75, pTM > 0.75): {len(relaxed)} designs")
        for r in relaxed:
            print(f"  - {r['name']}: pLDDT={r['plddt']:.3f}, pTM={r['ptm']:.3f}")

# AF3 submission prep
print("\n" + "=" * 60)
print("AF3 SUBMISSION PREPARATION")
print("=" * 60)

# All designs with pLDDT > 0.80
af3_candidates = [r for r in sorted_results if r['plddt'] > 0.80]
print(f"Candidates for AF3 (pLDDT > 0.80): {len(af3_candidates)}")

af3_fasta = OUTPUT_DIR / "R7_AF3_candidates.fasta"
with open(af3_fasta, 'w') as f:
    for r in af3_candidates:
        f.write(f">{r['name']}_pLDDT{r['plddt']:.2f}\n")
        f.write(f"{r['seq']}\n")
print(f"Saved all candidates: {af3_fasta}")
print("\nFor AF3 submission, include:")
print("  - Protein sequence (chain A)")
print("  - Tb³⁺ ion (chain B)")
print("  - Citrate ligand (chain C)")
