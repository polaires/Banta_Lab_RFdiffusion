"""
Compare UnifiedDesignAnalyzer vs direct metal_validation call.
"""

import sys
import json
from pathlib import Path

# Add serverless path for imports
serverless_path = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless")
sys.path.insert(0, str(serverless_path))

from unified_analyzer import UnifiedDesignAnalyzer
from metal_validation import validate_metal_ligand_complex_site

# Path to best design
PDB_PATH = Path(r"G:\Github_local_repo\Banta_Lab_RFdiffusion\experiments\Dy_TriNOx_scaffold\outputs_v3\v3_nterm_004.pdb")

# Dy-TriNOx specific parameters
DY_TRINOX_PARAMS = {
    "metal": "DY",
    "ligand_name": "UNL",
    "expected_ligand_donors": 3,  # O1, O2, O3 phenolates
    "expected_protein_donors": 2,  # At least 2 from protein
    "metal_coordination_target": 8,
    "distance_cutoff": 3.0,
    "ligand_hbond_cutoff": 3.5,
    "ligand_hbond_acceptors": ["N1", "N2", "N3", "N4"],
}


def main():
    print("=" * 70)
    print("ANALYZER COMPARISON")
    print("=" * 70)

    # Read PDB content
    with open(PDB_PATH, 'r') as f:
        pdb_content = f.read()

    print(f"\nAnalyzing: {PDB_PATH.name}")

    # 1. UnifiedDesignAnalyzer
    print("\n" + "-" * 70)
    print("1. UNIFIED DESIGN ANALYZER")
    print("-" * 70)

    analyzer = UnifiedDesignAnalyzer()
    result = analyzer.analyze(
        pdb_content=pdb_content,
        design_params={"source": "v3_nterm_004"},
        metal_type="DY",
    )

    metal_coord = result["analyses"]["metal_coordination"]
    print(f"  Status: {metal_coord['status']}")
    if metal_coord['status'] == 'success':
        metrics = metal_coord['metrics']
        print(f"  Coordination number: {metrics['coordination_number']}")
        print(f"  Geometry type: {metrics['geometry_type']}")
        print(f"  Coordinating residues: {metrics['coordinating_residues']}")
    else:
        print(f"  Error: {metal_coord.get('reason', 'unknown')}")

    # 2. Direct validate_metal_ligand_complex_site call
    print("\n" + "-" * 70)
    print("2. DIRECT METAL-LIGAND COMPLEX VALIDATION")
    print("-" * 70)

    result2 = validate_metal_ligand_complex_site(
        pdb_content=pdb_content,
        **DY_TRINOX_PARAMS
    )

    print(f"  Success: {result2.get('success')}")
    print(f"  Total coordination: {result2.get('coordination_number')}")
    print(f"  Ligand donors: {result2.get('ligand_coordination')}/{DY_TRINOX_PARAMS['expected_ligand_donors']}")
    print(f"  Protein donors: {result2.get('protein_coordination')}/{DY_TRINOX_PARAMS['expected_protein_donors']}")

    if result2.get('ligand_donors'):
        print(f"\n  Ligand donor details:")
        for d in result2['ligand_donors']:
            print(f"    {d['atom']}: {d['distance']}A")

    if result2.get('donor_residues'):
        print(f"\n  Protein donor details:")
        for d in result2['donor_residues']:
            print(f"    {d['chain']}{d['resnum']} {d['resname']} {d['atom']}: {d['distance']}A")

    if result2.get('issues'):
        print(f"\n  Issues:")
        for issue in result2['issues']:
            print(f"    - {issue}")

    # 3. Manual contact analysis (our script)
    print("\n" + "-" * 70)
    print("3. MANUAL TRINOX-PROTEIN CONTACT ANALYSIS")
    print("-" * 70)

    # Quick manual analysis
    import math

    atoms = []
    for line in pdb_content.split('\n'):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21]
                res_num = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms.append({
                    "name": atom_name, "res_name": res_name,
                    "chain": chain, "res_num": res_num,
                    "x": x, "y": y, "z": z
                })
            except:
                pass

    # TriNOx atoms
    trinox_o = [a for a in atoms if a["res_name"] == "UNL" and a["name"] in ["O1", "O2", "O3"]]
    trinox_n = [a for a in atoms if a["res_name"] == "UNL" and a["name"] in ["N1", "N2", "N3", "N4"]]
    trinox_c = [a for a in atoms if a["res_name"] == "UNL" and a["name"].startswith("C")]

    # Protein atoms
    protein_aa = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
                  "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    protein_atoms = [a for a in atoms if a["res_name"] in protein_aa]

    # H-bond donors from protein
    hbond_donors = ["N", "NE", "NH1", "NH2", "ND2", "NE2", "NZ", "OG", "OG1", "OH", "NE1"]

    def dist(a1, a2):
        return math.sqrt((a1["x"]-a2["x"])**2 + (a1["y"]-a2["y"])**2 + (a1["z"]-a2["z"])**2)

    # Count H-bonds to TriNOx (O and N)
    hbond_acceptors = trinox_o + trinox_n
    hbond_count = 0
    for ta in hbond_acceptors:
        for pa in protein_atoms:
            if pa["name"] in hbond_donors:
                if dist(ta, pa) <= 3.5:
                    hbond_count += 1
                    break

    # Count hydrophobic contacts to TriNOx carbons
    hydrophobic_res = ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "TYR", "PRO"]
    hydrophobic_atoms = [a for a in protein_atoms if a["res_name"] in hydrophobic_res]

    contacted_carbons = set()
    for tc in trinox_c:
        for pa in hydrophobic_atoms:
            if dist(tc, pa) <= 4.5:
                contacted_carbons.add(tc["name"])
                break

    carbon_coverage = len(contacted_carbons) / len(trinox_c) if trinox_c else 0

    print(f"  H-bonds to TriNOx: {hbond_count}")
    print(f"  Hydrophobic coverage: {carbon_coverage:.0%} ({len(contacted_carbons)}/{len(trinox_c)} carbons)")

    # Aromatic stacking
    aromatic_residues = [a for a in protein_atoms if a["res_name"] in ["PHE", "TYR", "TRP"]]
    aromatic_contacts = 0
    for tc in trinox_c:
        for pa in aromatic_residues:
            if dist(tc, pa) <= 5.0:
                aromatic_contacts += 1
                break

    print(f"  Aromatic stacking contacts: {aromatic_contacts}")

    # Summary
    print("\n" + "=" * 70)
    print("COMPARISON SUMMARY")
    print("=" * 70)

    print(f"""
  UnifiedDesignAnalyzer:
    - Uses validate_lanthanide_site (simpler, for pure metal sites)
    - Found CN=11 but no residue details
    - Does NOT analyze TriNOx-protein interactions
    - Does NOT separate ligand vs protein donors

  Direct metal_validation:
    - Uses validate_metal_ligand_complex_site
    - Separates ligand donors ({result2.get('ligand_coordination')}) from protein donors ({result2.get('protein_coordination')})
    - More appropriate for metal-ligand complexes

  Manual analysis:
    - Focuses on TriNOx-protein interactions (per user's priority)
    - H-bonds: {hbond_count}
    - Hydrophobic coverage: {carbon_coverage:.0%}
    - Aromatic stacking: {aromatic_contacts}

  RECOMMENDATION:
    UnifiedDesignAnalyzer should be enhanced to:
    1. Use validate_metal_ligand_complex_site for metal-ligand complexes
    2. Add TriNOx-protein contact analysis
    3. Fix pLDDT extraction for RFD3 PDBs
""")


if __name__ == "__main__":
    main()
