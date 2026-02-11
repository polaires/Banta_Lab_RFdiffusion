"""
Save final results from completed RF3 validation.
"""

import json
from pathlib import Path
from datetime import datetime

OUTPUT_DIR = Path("G:/Github_local_repo/Banta_Lab_RFdiffusion/experiments/Dy_TriNOx_scaffold/outputs")

# All 104 results from the completed RF3 validation
ALL_RESULTS = [
    # Small designs (r1a_small) - 63 aa
    {"name": "r1a_small_000_seq01", "backbone": "r1a_small_000", "length": 63, "plddt": 0.8776, "ptm": 0.8595, "pae": 2.46},
    {"name": "r1a_small_000_seq02", "backbone": "r1a_small_000", "length": 63, "plddt": 0.8751, "ptm": 0.8593, "pae": 2.51},
    {"name": "r1a_small_000_seq03", "backbone": "r1a_small_000", "length": 63, "plddt": 0.8778, "ptm": 0.8728, "pae": 2.02},
    {"name": "r1a_small_000_seq04", "backbone": "r1a_small_000", "length": 63, "plddt": 0.8897, "ptm": 0.8733, "pae": 1.96},
    {"name": "r1a_small_000_seq05", "backbone": "r1a_small_000", "length": 63, "plddt": 0.8454, "ptm": 0.8318, "pae": 2.84},
    {"name": "r1a_small_000_seq06", "backbone": "r1a_small_000", "length": 63, "plddt": 0.8540, "ptm": 0.8362, "pae": 2.99},
    {"name": "r1a_small_000_seq07", "backbone": "r1a_small_000", "length": 63, "plddt": 0.8850, "ptm": 0.8379, "pae": 2.69},
    {"name": "r1a_small_000_seq08", "backbone": "r1a_small_000", "length": 63, "plddt": 0.8707, "ptm": 0.8422, "pae": 2.47},
    {"name": "r1a_small_001_seq01", "backbone": "r1a_small_001", "length": 63, "plddt": 0.8860, "ptm": 0.7958, "pae": 3.10},
    {"name": "r1a_small_001_seq02", "backbone": "r1a_small_001", "length": 63, "plddt": 0.8710, "ptm": 0.7698, "pae": 3.33},
    {"name": "r1a_small_001_seq03", "backbone": "r1a_small_001", "length": 63, "plddt": 0.9184, "ptm": 0.8742, "pae": 1.84},
    {"name": "r1a_small_001_seq04", "backbone": "r1a_small_001", "length": 63, "plddt": 0.8968, "ptm": 0.7733, "pae": 3.12},
    {"name": "r1a_small_001_seq05", "backbone": "r1a_small_001", "length": 63, "plddt": 0.9072, "ptm": 0.8522, "pae": 2.17},
    {"name": "r1a_small_001_seq06", "backbone": "r1a_small_001", "length": 63, "plddt": 0.9108, "ptm": 0.8553, "pae": 2.18},
    {"name": "r1a_small_001_seq07", "backbone": "r1a_small_001", "length": 63, "plddt": 0.9050, "ptm": 0.8248, "pae": 2.38},
    {"name": "r1a_small_001_seq08", "backbone": "r1a_small_001", "length": 63, "plddt": 0.8976, "ptm": 0.8202, "pae": 2.62},
    {"name": "r1a_small_002_seq01", "backbone": "r1a_small_002", "length": 63, "plddt": 0.8364, "ptm": 0.6627, "pae": 6.36},
    {"name": "r1a_small_002_seq02", "backbone": "r1a_small_002", "length": 63, "plddt": 0.8487, "ptm": 0.7016, "pae": 5.43},
    {"name": "r1a_small_002_seq03", "backbone": "r1a_small_002", "length": 63, "plddt": 0.8358, "ptm": 0.6669, "pae": 6.02},
    {"name": "r1a_small_002_seq04", "backbone": "r1a_small_002", "length": 63, "plddt": 0.8708, "ptm": 0.7325, "pae": 4.79},
    {"name": "r1a_small_002_seq05", "backbone": "r1a_small_002", "length": 63, "plddt": 0.8649, "ptm": 0.7392, "pae": 4.40},
    {"name": "r1a_small_002_seq06", "backbone": "r1a_small_002", "length": 63, "plddt": 0.8532, "ptm": 0.6213, "pae": 6.74},
    {"name": "r1a_small_002_seq07", "backbone": "r1a_small_002", "length": 63, "plddt": 0.8639, "ptm": 0.7238, "pae": 4.87},
    {"name": "r1a_small_002_seq08", "backbone": "r1a_small_002", "length": 63, "plddt": 0.8826, "ptm": 0.8110, "pae": 3.08},
    {"name": "r1a_small_003_seq01", "backbone": "r1a_small_003", "length": 63, "plddt": 0.8113, "ptm": 0.5274, "pae": 8.24},
    {"name": "r1a_small_003_seq02", "backbone": "r1a_small_003", "length": 63, "plddt": 0.8208, "ptm": 0.5354, "pae": 7.92},
    {"name": "r1a_small_003_seq03", "backbone": "r1a_small_003", "length": 63, "plddt": 0.8180, "ptm": 0.5012, "pae": 9.27},
    {"name": "r1a_small_003_seq04", "backbone": "r1a_small_003", "length": 63, "plddt": 0.8092, "ptm": 0.5309, "pae": 7.91},
    {"name": "r1a_small_003_seq05", "backbone": "r1a_small_003", "length": 63, "plddt": 0.8099, "ptm": 0.5236, "pae": 8.03},
    {"name": "r1a_small_003_seq06", "backbone": "r1a_small_003", "length": 63, "plddt": 0.8039, "ptm": 0.5177, "pae": 8.54},
    {"name": "r1a_small_003_seq07", "backbone": "r1a_small_003", "length": 63, "plddt": 0.8150, "ptm": 0.5124, "pae": 8.87},
    {"name": "r1a_small_003_seq08", "backbone": "r1a_small_003", "length": 63, "plddt": 0.8006, "ptm": 0.4874, "pae": 9.52},
    {"name": "r1a_small_004_seq01", "backbone": "r1a_small_004", "length": 63, "plddt": 0.8087, "ptm": 0.6016, "pae": 6.56},
    {"name": "r1a_small_004_seq02", "backbone": "r1a_small_004", "length": 63, "plddt": 0.8221, "ptm": 0.6402, "pae": 5.94},
    {"name": "r1a_small_004_seq03", "backbone": "r1a_small_004", "length": 63, "plddt": 0.7929, "ptm": 0.5829, "pae": 6.85},
    {"name": "r1a_small_004_seq04", "backbone": "r1a_small_004", "length": 63, "plddt": 0.8226, "ptm": 0.6217, "pae": 5.96},
    {"name": "r1a_small_004_seq05", "backbone": "r1a_small_004", "length": 63, "plddt": 0.8256, "ptm": 0.6593, "pae": 5.39},
    {"name": "r1a_small_004_seq06", "backbone": "r1a_small_004", "length": 63, "plddt": 0.8264, "ptm": 0.6904, "pae": 4.66},
    {"name": "r1a_small_004_seq07", "backbone": "r1a_small_004", "length": 63, "plddt": 0.8034, "ptm": 0.6798, "pae": 5.46},
    {"name": "r1a_small_004_seq08", "backbone": "r1a_small_004", "length": 63, "plddt": 0.8045, "ptm": 0.6233, "pae": 5.70},

    # Medium designs (r1b_medium) - 108 aa
    {"name": "r1b_medium_000_seq01", "backbone": "r1b_medium_000", "length": 108, "plddt": 0.9117, "ptm": 0.9164, "pae": 2.33},
    {"name": "r1b_medium_000_seq02", "backbone": "r1b_medium_000", "length": 108, "plddt": 0.8995, "ptm": 0.9146, "pae": 2.53},
    {"name": "r1b_medium_000_seq03", "backbone": "r1b_medium_000", "length": 108, "plddt": 0.8922, "ptm": 0.9121, "pae": 2.47},
    {"name": "r1b_medium_000_seq04", "backbone": "r1b_medium_000", "length": 108, "plddt": 0.9223, "ptm": 0.9303, "pae": 2.04},
    {"name": "r1b_medium_000_seq05", "backbone": "r1b_medium_000", "length": 108, "plddt": 0.8797, "ptm": 0.9153, "pae": 2.43},
    {"name": "r1b_medium_000_seq06", "backbone": "r1b_medium_000", "length": 108, "plddt": 0.8848, "ptm": 0.9185, "pae": 2.48},
    {"name": "r1b_medium_000_seq07", "backbone": "r1b_medium_000", "length": 108, "plddt": 0.9021, "ptm": 0.9237, "pae": 2.03},
    {"name": "r1b_medium_000_seq08", "backbone": "r1b_medium_000", "length": 108, "plddt": 0.9297, "ptm": 0.9396, "pae": 1.86},
    {"name": "r1b_medium_001_seq01", "backbone": "r1b_medium_001", "length": 108, "plddt": 0.8698, "ptm": 0.8343, "pae": 4.17},
    {"name": "r1b_medium_001_seq02", "backbone": "r1b_medium_001", "length": 108, "plddt": 0.8672, "ptm": 0.8552, "pae": 3.48},
    {"name": "r1b_medium_001_seq03", "backbone": "r1b_medium_001", "length": 108, "plddt": 0.8670, "ptm": 0.8586, "pae": 3.41},
    {"name": "r1b_medium_001_seq04", "backbone": "r1b_medium_001", "length": 108, "plddt": 0.8481, "ptm": 0.8462, "pae": 3.74},
    {"name": "r1b_medium_001_seq05", "backbone": "r1b_medium_001", "length": 108, "plddt": 0.8415, "ptm": 0.8501, "pae": 4.05},
    {"name": "r1b_medium_001_seq06", "backbone": "r1b_medium_001", "length": 108, "plddt": 0.8849, "ptm": 0.8791, "pae": 3.08},
    {"name": "r1b_medium_001_seq07", "backbone": "r1b_medium_001", "length": 108, "plddt": 0.8616, "ptm": 0.8493, "pae": 3.58},
    {"name": "r1b_medium_001_seq08", "backbone": "r1b_medium_001", "length": 108, "plddt": 0.7734, "ptm": 0.6572, "pae": 8.56},
    {"name": "r1b_medium_002_seq01", "backbone": "r1b_medium_002", "length": 108, "plddt": 0.7918, "ptm": 0.6351, "pae": 9.16},
    {"name": "r1b_medium_002_seq02", "backbone": "r1b_medium_002", "length": 108, "plddt": 0.8494, "ptm": 0.8193, "pae": 4.17},
    {"name": "r1b_medium_002_seq03", "backbone": "r1b_medium_002", "length": 108, "plddt": 0.8284, "ptm": 0.7017, "pae": 7.36},
    {"name": "r1b_medium_002_seq04", "backbone": "r1b_medium_002", "length": 108, "plddt": 0.8489, "ptm": 0.7887, "pae": 5.02},
    {"name": "r1b_medium_002_seq05", "backbone": "r1b_medium_002", "length": 108, "plddt": 0.8400, "ptm": 0.8065, "pae": 4.46},
    {"name": "r1b_medium_002_seq06", "backbone": "r1b_medium_002", "length": 108, "plddt": 0.8488, "ptm": 0.8165, "pae": 4.24},
    {"name": "r1b_medium_002_seq07", "backbone": "r1b_medium_002", "length": 108, "plddt": 0.8522, "ptm": 0.8410, "pae": 3.99},
    {"name": "r1b_medium_002_seq08", "backbone": "r1b_medium_002", "length": 108, "plddt": 0.8405, "ptm": 0.7703, "pae": 5.13},
    {"name": "r1b_medium_003_seq01", "backbone": "r1b_medium_003", "length": 108, "plddt": 0.8468, "ptm": 0.6838, "pae": 6.66},
    {"name": "r1b_medium_003_seq02", "backbone": "r1b_medium_003", "length": 108, "plddt": 0.8579, "ptm": 0.8338, "pae": 4.31},
    {"name": "r1b_medium_003_seq03", "backbone": "r1b_medium_003", "length": 108, "plddt": 0.8893, "ptm": 0.8204, "pae": 4.23},
    {"name": "r1b_medium_003_seq04", "backbone": "r1b_medium_003", "length": 108, "plddt": 0.8601, "ptm": 0.8064, "pae": 3.90},
    {"name": "r1b_medium_003_seq05", "backbone": "r1b_medium_003", "length": 108, "plddt": 0.8786, "ptm": 0.8381, "pae": 3.42},
    {"name": "r1b_medium_003_seq06", "backbone": "r1b_medium_003", "length": 108, "plddt": 0.8545, "ptm": 0.7917, "pae": 4.04},
    {"name": "r1b_medium_003_seq07", "backbone": "r1b_medium_003", "length": 108, "plddt": 0.8692, "ptm": 0.8307, "pae": 3.91},
    {"name": "r1b_medium_003_seq08", "backbone": "r1b_medium_003", "length": 108, "plddt": 0.8331, "ptm": 0.6046, "pae": 8.71},
    {"name": "r1b_medium_004_seq01", "backbone": "r1b_medium_004", "length": 108, "plddt": 0.8162, "ptm": 0.7213, "pae": 5.31},
    {"name": "r1b_medium_004_seq02", "backbone": "r1b_medium_004", "length": 108, "plddt": 0.8617, "ptm": 0.7643, "pae": 4.87},
    {"name": "r1b_medium_004_seq03", "backbone": "r1b_medium_004", "length": 108, "plddt": 0.8518, "ptm": 0.7622, "pae": 4.76},
    {"name": "r1b_medium_004_seq04", "backbone": "r1b_medium_004", "length": 108, "plddt": 0.8527, "ptm": 0.7503, "pae": 4.56},
    {"name": "r1b_medium_004_seq05", "backbone": "r1b_medium_004", "length": 108, "plddt": 0.8356, "ptm": 0.7410, "pae": 4.34},
    {"name": "r1b_medium_004_seq06", "backbone": "r1b_medium_004", "length": 108, "plddt": 0.8559, "ptm": 0.7741, "pae": 3.97},
    {"name": "r1b_medium_004_seq07", "backbone": "r1b_medium_004", "length": 108, "plddt": 0.8438, "ptm": 0.6963, "pae": 6.33},
    {"name": "r1b_medium_004_seq08", "backbone": "r1b_medium_004", "length": 108, "plddt": 0.8702, "ptm": 0.8146, "pae": 3.65},

    # Large designs (r1c_large) - 131 aa
    {"name": "r1c_large_000_seq01", "backbone": "r1c_large_000", "length": 131, "plddt": 0.8668, "ptm": 0.7023, "pae": 7.59},
    {"name": "r1c_large_000_seq02", "backbone": "r1c_large_000", "length": 131, "plddt": 0.8583, "ptm": 0.7648, "pae": 6.78},
    {"name": "r1c_large_000_seq03", "backbone": "r1c_large_000", "length": 131, "plddt": 0.8606, "ptm": 0.7142, "pae": 8.28},
    {"name": "r1c_large_000_seq04", "backbone": "r1c_large_000", "length": 131, "plddt": 0.8485, "ptm": 0.6901, "pae": 8.91},
    {"name": "r1c_large_000_seq05", "backbone": "r1c_large_000", "length": 131, "plddt": 0.8846, "ptm": 0.7582, "pae": 6.17},
    {"name": "r1c_large_000_seq06", "backbone": "r1c_large_000", "length": 131, "plddt": 0.8799, "ptm": 0.7359, "pae": 6.47},
    {"name": "r1c_large_000_seq07", "backbone": "r1c_large_000", "length": 131, "plddt": 0.8807, "ptm": 0.7376, "pae": 6.54},
    {"name": "r1c_large_000_seq08", "backbone": "r1c_large_000", "length": 131, "plddt": 0.8964, "ptm": 0.8311, "pae": 4.29},
    {"name": "r1c_large_001_seq01", "backbone": "r1c_large_001", "length": 131, "plddt": 0.8404, "ptm": 0.5578, "pae": 9.89},
    {"name": "r1c_large_001_seq02", "backbone": "r1c_large_001", "length": 131, "plddt": 0.8170, "ptm": 0.5813, "pae": 10.62},
    {"name": "r1c_large_001_seq03", "backbone": "r1c_large_001", "length": 131, "plddt": 0.8125, "ptm": 0.6232, "pae": 10.51},
    {"name": "r1c_large_001_seq04", "backbone": "r1c_large_001", "length": 131, "plddt": 0.7922, "ptm": 0.5678, "pae": 9.88},
    {"name": "r1c_large_001_seq05", "backbone": "r1c_large_001", "length": 131, "plddt": 0.8014, "ptm": 0.4857, "pae": 12.71},
    {"name": "r1c_large_001_seq06", "backbone": "r1c_large_001", "length": 131, "plddt": 0.8162, "ptm": 0.6484, "pae": 8.53},
    {"name": "r1c_large_001_seq07", "backbone": "r1c_large_001", "length": 131, "plddt": 0.8195, "ptm": 0.5839, "pae": 10.36},
    {"name": "r1c_large_001_seq08", "backbone": "r1c_large_001", "length": 131, "plddt": 0.8379, "ptm": 0.6032, "pae": 9.43},
    {"name": "r1c_large_002_seq01", "backbone": "r1c_large_002", "length": 131, "plddt": 0.7740, "ptm": 0.4183, "pae": 11.78},
    {"name": "r1c_large_002_seq02", "backbone": "r1c_large_002", "length": 131, "plddt": 0.7490, "ptm": 0.4935, "pae": 10.69},
    {"name": "r1c_large_002_seq03", "backbone": "r1c_large_002", "length": 131, "plddt": 0.7715, "ptm": 0.4900, "pae": 10.83},
    {"name": "r1c_large_002_seq04", "backbone": "r1c_large_002", "length": 131, "plddt": 0.7999, "ptm": 0.4972, "pae": 10.86},
    {"name": "r1c_large_002_seq05", "backbone": "r1c_large_002", "length": 131, "plddt": 0.7662, "ptm": 0.5144, "pae": 11.12},
    {"name": "r1c_large_002_seq06", "backbone": "r1c_large_002", "length": 131, "plddt": 0.7842, "ptm": 0.4849, "pae": 11.16},
    {"name": "r1c_large_002_seq07", "backbone": "r1c_large_002", "length": 131, "plddt": 0.7678, "ptm": 0.4291, "pae": 11.22},
    {"name": "r1c_large_002_seq08", "backbone": "r1c_large_002", "length": 131, "plddt": 0.7805, "ptm": 0.5624, "pae": 11.46},
]

def main():
    results = ALL_RESULTS
    sorted_results = sorted(results, key=lambda x: x['plddt'], reverse=True)

    # Quality thresholds
    high_quality = [r for r in results if r['plddt'] >= 0.85 and r['ptm'] >= 0.80]
    medium_quality = [r for r in results if r['plddt'] >= 0.75 and r['ptm'] >= 0.70]

    print("=" * 70)
    print("QUALITY SUMMARY")
    print("=" * 70)
    print(f"Total sequences: {len(results)}")
    print(f"High quality (pLDDT>=0.85, pTM>=0.80): {len(high_quality)}")
    print(f"Medium quality (pLDDT>=0.75, pTM>=0.70): {len(medium_quality)}")

    # Best design
    best = sorted_results[0]
    print(f"\nBEST DESIGN: {best['name']}")
    print(f"  Backbone: {best['backbone']}")
    print(f"  Length: {best['length']} aa")
    print(f"  pLDDT: {best['plddt']:.4f}")
    print(f"  pTM: {best['ptm']:.4f}")
    print(f"  PAE: {best['pae']:.2f}")

    # Save full results
    results_file = OUTPUT_DIR / "trinox_workflow_results.json"
    with open(results_file, 'w') as f:
        json.dump({
            "timestamp": datetime.now().isoformat(),
            "total_backbones": 13,
            "total_sequences": len(results),
            "high_quality_count": len(high_quality),
            "medium_quality_count": len(medium_quality),
            "results": sorted_results,
        }, f, indent=2)
    print(f"\nFull results saved to {results_file}")

    # Print top 20
    print("\n" + "=" * 70)
    print("TOP 20 DESIGNS")
    print("=" * 70)
    print(f"{'Design':<30} {'Backbone':<18} {'Len':>5} {'pLDDT':>8} {'pTM':>8} {'PAE':>8}")
    print("-" * 85)
    for r in sorted_results[:20]:
        print(f"{r['name']:<30} {r['backbone']:<18} {r['length']:>5} {r['plddt']:>8.4f} {r['ptm']:>8.4f} {r['pae']:>8.2f}")

    # Breakdown by backbone size
    print("\n" + "=" * 70)
    print("BREAKDOWN BY SIZE")
    print("=" * 70)

    small = [r for r in results if r['length'] == 63]
    medium = [r for r in results if r['length'] == 108]
    large = [r for r in results if r['length'] == 131]

    small_hq = [r for r in small if r['plddt'] >= 0.85 and r['ptm'] >= 0.80]
    medium_hq = [r for r in medium if r['plddt'] >= 0.85 and r['ptm'] >= 0.80]
    large_hq = [r for r in large if r['plddt'] >= 0.85 and r['ptm'] >= 0.80]

    print(f"Small (63 aa): {len(small)} total, {len(small_hq)} high quality")
    print(f"Medium (108 aa): {len(medium)} total, {len(medium_hq)} high quality")
    print(f"Large (131 aa): {len(large)} total, {len(large_hq)} high quality")

if __name__ == "__main__":
    main()
