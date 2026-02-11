"""
Test UnifiedDesignAnalyzer on the best v3 design.

Enhanced version that tests:
- Metal-ligand complex validation (separates ligand vs protein donors)
- PLIP-based ligand interaction analysis (H-bonds, hydrophobic, pi-stacking)
- Fixed pLDDT detection for RFD3 PDBs
"""

import sys
import json
from pathlib import Path

# Add serverless path for imports
serverless_path = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless")
sys.path.insert(0, str(serverless_path))

from unified_analyzer import UnifiedDesignAnalyzer

# Path to best design
PDB_PATH = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\experiments\Dy_TriNOx_scaffold\outputs_v3\v3_nterm_004.pdb")

def main():
    print("=" * 70)
    print("UNIFIED DESIGN ANALYZER TEST")
    print("=" * 70)

    # Read PDB content
    with open(PDB_PATH, 'r') as f:
        pdb_content = f.read()

    print(f"\nAnalyzing: {PDB_PATH.name}")
    print(f"PDB size: {len(pdb_content)} bytes")

    # Initialize analyzer
    analyzer = UnifiedDesignAnalyzer()
    print(f"\nAvailable modules: {analyzer._analysis_modules}")

    # Design params (what was used to generate this design)
    design_params = {
        "contig": "X1,L1,30-50,M1,30-50",
        "ligand": "DY,UNL",
        "select_fixed_atoms": {"X1": "all", "L1": "all", "M1": "all"},
        "select_partially_buried": {"L1": "all"},
        "select_buried": {"X1": "all"},
        "cfg_scale": 2.5,
    }

    # Run analysis
    # Note: metal_type="DY" for Dysprosium, ligand is UNL (TriNOx)
    print("\nRunning unified analysis...")
    result = analyzer.analyze(
        pdb_content=pdb_content,
        design_params=design_params,
        pdb_path=str(PDB_PATH),
        metal_type="DY",  # Dysprosium
    )

    # Print results
    print("\n" + "=" * 70)
    print("ANALYSIS RESULTS")
    print("=" * 70)

    print(f"\nDesign ID: {result.get('design_id')}")
    print(f"Design Type: {result.get('design_type')}")
    print(f"Timestamp: {result.get('timestamp')}")

    # Auto-detected info
    if "auto_detected" in result:
        print(f"\nAuto-detected:")
        for key, val in result["auto_detected"].items():
            print(f"  {key}: {val}")

    # PDB sanitization
    if "pdb_sanitization" in result:
        sanitization = result["pdb_sanitization"]
        print(f"\nPDB Sanitization ({sanitization['count']} issues fixed):")
        for issue in sanitization["issues_fixed"][:5]:
            print(f"  - {issue}")
        if sanitization['count'] > 5:
            print(f"  ... and {sanitization['count'] - 5} more")

    # Analysis results
    print("\nAnalyses:")
    for analysis_name, analysis_result in result.get("analyses", {}).items():
        status = analysis_result.get("status", "unknown")
        print(f"\n  [{analysis_name}] - {status}")

        if status == "success":
            metrics = analysis_result.get("metrics", {})
            for metric_name, metric_value in metrics.items():
                if isinstance(metric_value, float):
                    print(f"    {metric_name}: {metric_value:.3f}")
                elif isinstance(metric_value, list) and len(metric_value) > 3:
                    print(f"    {metric_name}: [{len(metric_value)} items]")
                else:
                    print(f"    {metric_name}: {metric_value}")
        elif status == "skipped":
            print(f"    Reason: {analysis_result.get('reason', 'unknown')}")
        elif status == "not_applicable":
            print(f"    Reason: {analysis_result.get('reason', 'unknown')}")

    # Save full results to JSON
    output_file = PDB_PATH.parent / "validation" / f"{PDB_PATH.stem}_unified_analysis.json"
    output_file.parent.mkdir(exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"\n\nFull results saved to: {output_file}")


if __name__ == "__main__":
    main()
