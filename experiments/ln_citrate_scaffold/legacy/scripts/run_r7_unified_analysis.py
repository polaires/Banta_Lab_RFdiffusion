"""
Run comprehensive unified analysis on R7 designs.

Uses:
1. Metal-ligand complex validation (TB-citrate coordination accounting)
2. RF3 structure confidence
3. Sequence composition analysis
"""

import sys
import os
from pathlib import Path

# Add serverless dir to path for metal_validation module
serverless_path = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless")
sys.path.insert(0, str(serverless_path))

# Also add utils path
utils_path = serverless_path / "utils"
if utils_path.exists():
    sys.path.insert(0, str(utils_path))

from metal_validation import validate_metal_ligand_complex_site

OUTPUT_DIR = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\experiments\ln_citrate_scaffold\outputs\round_07_scaffold")

# Best R7 sequences and their RF3 validation results
r7_designs = {
    "r7_design_1": {
        "sequence": "TLRLTITFAPGDLKVTFFDAETGEKLGTYVGRDAIIAANNELRDAGVWHTEVVARADAPEKAVRDSVPHRTLSIEQIAPNTIVAVVETADPAAFVAEETAEVAALGGSLTYEVL",
        "rf3_plddt": 0.8646,
        "rf3_ptm": 0.8913,
        "rf3_pae": 2.34,
    },
    "r7_design_7": {
        "sequence": "TLRVTITWAPGEKLVTFFDAETKEKLGTYVGRDAILAAKNELTAAGKWHTEEVALADRKPEKAVRESVPHKILSLEQVAPDTVRAVIETDDPDALIAAETAAVAAMGGSMTAERL",
        "rf3_plddt": 0.8276,
        "rf3_ptm": 0.7157,
        "rf3_pae": 5.58,
    },
    "r7_design_4": {
        "sequence": "RLRVTITYAPGEKLVRFFDAETELLGTYVGMDAILAANADAAAGKWITEEVALADRRPEKAVRDSVPHRILSLERIAPNTRAVVETDDPAAFAAEETAEVAAMGGSMTWEVL",
        "rf3_plddt": 0.7977,
        "rf3_ptm": 0.5400,
        "rf3_pae": 12.36,
    },
    "r7_design_8": {
        "sequence": "RLRLTITYAPGDLLVRFFNAETGELLGTFVGRDAILAANEELAAGVWHTEEVARADRPEDAVAETTPHILSLEQVAPNTVVAEVDTADPDALVAEETAAVAALGGSLTYERL",
        "rf3_plddt": 0.8139,
        "rf3_ptm": 0.5932,
        "rf3_pae": 10.16,
    },
}


def analyze_sequence_composition(seq):
    """Analyze sequence composition for metal-binding residues."""
    length = len(seq)
    composition = {}

    # Key residues for lanthanide binding
    binding_residues = {
        "E": "Glu (carboxylate)",
        "D": "Asp (carboxylate)",
        "H": "His (imidazole)",
        "C": "Cys (thiol)",
        "S": "Ser (hydroxyl)",
        "T": "Thr (hydroxyl)",
    }

    for aa, name in binding_residues.items():
        count = seq.count(aa)
        composition[aa] = {
            "count": count,
            "percent": round(100 * count / length, 1),
            "name": name,
        }

    # Total carboxylate (E+D)
    carbox_count = seq.count("E") + seq.count("D")
    composition["carboxylate_total"] = {
        "count": carbox_count,
        "percent": round(100 * carbox_count / length, 1),
    }

    # Alanine content (check for low diversity)
    ala_count = seq.count("A")
    composition["A"] = {
        "count": ala_count,
        "percent": round(100 * ala_count / length, 1),
        "flag": "high" if ala_count / length > 0.25 else "ok",
    }

    return composition


def run_metal_ligand_validation(pdb_content):
    """
    Run metal-ligand complex validation.

    For TB-citrate:
    - Citrate provides 4 donors to TB (O1, O3, O5, O7)
    - TB has coordination capacity of 9
    - Protein should provide 5 additional donors
    - Citrate has 3 free carboxylates for H-bonds (O2, O4, O6)
    """
    result = validate_metal_ligand_complex_site(
        pdb_content=pdb_content,
        metal="TB",
        ligand_name="CIT",
        expected_ligand_donors=4,      # Citrate -> TB
        expected_protein_donors=5,     # Protein -> TB (to reach CN=9)
        metal_coordination_target=9,   # Tb3+ coordination capacity
        distance_cutoff=3.0,           # Angstroms
        ligand_hbond_cutoff=3.5,       # For citrate-protein H-bonds
        ligand_hbond_acceptors=["O2", "O4", "O6"],  # Free carboxylates
    )
    return result


print("=" * 70)
print("R7 UNIFIED ANALYSIS")
print("=" * 70)

# Check if we have the RFD3 scaffold PDB
scaffold_pdb = OUTPUT_DIR / "r7_single_000.pdb"
if scaffold_pdb.exists():
    with open(scaffold_pdb, 'r') as f:
        scaffold_content = f.read()

    print("\n[1] METAL-LIGAND COMPLEX VALIDATION")
    print("-" * 60)
    print("Analyzing TB-citrate coordination in RFD3 scaffold...")

    validation = run_metal_ligand_validation(scaffold_content)

    if validation.get("success"):
        print(f"  Quality Score: {validation['quality_score']}/150")
        print(f"  Metal Coordination: {validation['coordination_number']}/{validation['metal_coordination_target']}")
        print(f"  - Ligand donors: {validation['ligand_coordination']}/{validation['ligand_budget']['expected']}")
        print(f"  - Protein donors: {validation['protein_coordination']}/{validation['protein_budget']['expected']}")
        print(f"  - Chain A: {validation['chain_a_donors']} donors")
        print(f"  - Chain B: {validation['chain_b_donors']} donors")
        print(f"  Ligand-Protein H-bonds: {validation['ligand_protein_hbond_count']}")

        if validation.get("ligand_donors"):
            print("\n  Citrate -> TB coordination:")
            for donor in validation["ligand_donors"]:
                print(f"    {donor['atom']}: {donor['distance']} A")

        if validation.get("donor_residues"):
            print("\n  Protein -> TB coordination:")
            for donor in validation["donor_residues"][:5]:  # First 5
                print(f"    {donor['chain']}{donor['resnum']} {donor['resname']} {donor['atom']}: {donor['distance']} A")
    else:
        print(f"  Error: {validation.get('error', 'Unknown')}")
        if validation.get("issues"):
            for issue in validation["issues"]:
                print(f"  - {issue}")
else:
    print("Scaffold PDB not found. Skipping metal-ligand validation.")

print("\n[2] SEQUENCE COMPOSITION ANALYSIS")
print("-" * 60)

for name, data in r7_designs.items():
    seq = data["sequence"]
    comp = analyze_sequence_composition(seq)

    print(f"\n{name} ({len(seq)} aa):")
    print(f"  RF3: pLDDT={data['rf3_plddt']:.3f}, pTM={data['rf3_ptm']:.3f}, PAE={data['rf3_pae']:.1f}")
    print(f"  Carboxylate (E+D): {comp['carboxylate_total']['count']} ({comp['carboxylate_total']['percent']}%)")
    print(f"  Glu (E): {comp['E']['count']}, Asp (D): {comp['D']['count']}")
    print(f"  His (H): {comp['H']['count']}, Ser (S): {comp['S']['count']}, Thr (T): {comp['T']['count']}")
    print(f"  Alanine: {comp['A']['count']} ({comp['A']['percent']}%) [{comp['A']['flag']}]")

print("\n[3] DESIGN RANKING")
print("-" * 60)

# Rank by combined score: RF3 pLDDT * pTM
rankings = []
for name, data in r7_designs.items():
    combined_score = data["rf3_plddt"] * data["rf3_ptm"]
    rankings.append({
        "name": name,
        "plddt": data["rf3_plddt"],
        "ptm": data["rf3_ptm"],
        "pae": data["rf3_pae"],
        "combined": combined_score,
    })

rankings.sort(key=lambda x: x["combined"], reverse=True)

print(f"{'Rank':<5} {'Design':<15} {'pLDDT':>8} {'pTM':>8} {'PAE':>8} {'Score':>8}")
print("-" * 55)
for i, r in enumerate(rankings, 1):
    print(f"{i:<5} {r['name']:<15} {r['plddt']:>8.3f} {r['ptm']:>8.3f} {r['pae']:>8.1f} {r['combined']:>8.3f}")

print("\n[4] RECOMMENDATIONS FOR AF3")
print("-" * 60)

best = rankings[0]
print(f"Best design for AF3: {best['name']}")
print(f"  Combined score: {best['combined']:.4f}")
print(f"  Rationale: Highest pLDDT * pTM indicates best fold confidence")
print("\nExpected AF3 behavior:")
print("  - Citrate should remain coordinated to TB")
print("  - Protein should form H-bonds with citrate carboxylates")
print("  - Look for additional Glu/Asp coordination to TB")

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
