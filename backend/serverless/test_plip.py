#!/usr/bin/env python
"""Test PLIP integration in the serverless container."""

from shared.interaction_analysis import analyze_all_interactions

# Simple test PDB with protein and ligand
test_pdb = """ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00  0.00           O
HETATM    5  C1  UNL A   2       5.000   0.000   0.000  1.00  0.00           C
HETATM    6  O1  UNL A   2       5.000   1.500   0.000  1.00  0.00           O
END
"""

try:
    print("Testing PLIP integration...")
    result = analyze_all_interactions(test_pdb, "UNL")
    print(f"Analysis method: {result.analysis_method}")
    print(f"Total H-bonds: {len(result.hydrogen_bonds)}")
    print(f"Total hydrophobic: {len(result.hydrophobic_contacts)}")
    print(f"Total pi-stacking: {len(result.pi_stacking)}")
    print(f"Total salt bridges: {len(result.salt_bridges)}")
    print(f"Total halogen bonds: {len(result.halogen_bonds)}")
    print(f"Key residues: {result.key_residues}")
    print("Success!")
except Exception as e:
    import traceback
    print(f"Error: {e}")
    traceback.print_exc()
