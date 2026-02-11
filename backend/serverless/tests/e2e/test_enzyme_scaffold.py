"""Test Enzyme Scaffold Design approach.

This uses the CORRECT approach that works in the UI:
- Use `length` parameter (NOT contig with fixed residue ranges)
- Use `unindex` to mark catalytic residues as flexible
- Use `ligand` with comma-separated codes
- Use `select_fixed_atoms` to fix ligand/metal position

This is different from the blocker-based approach I was using before.
"""
import asyncio
import sys
import math
sys.path.insert(0, '.')

from scaffolding_workflow import ScaffoldingWorkflow
from api_client import run_rfd3_api


def check_ligand_clashes(pdb_content: str, ligand_name: str, threshold: float = 2.0):
    """Check for clashes between protein and ligand."""
    ligand_atoms = []
    for line in pdb_content.split('\n'):
        if line.startswith('HETATM') and line[17:20].strip() == ligand_name:
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                ligand_atoms.append((x, y, z))
            except:
                pass

    if not ligand_atoms:
        return {"has_clashes": False, "clash_count": 0, "min_distance": float('inf')}

    protein_atoms = []
    for line in pdb_content.split('\n'):
        if line.startswith('ATOM'):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                protein_atoms.append((x, y, z))
            except:
                pass

    if not protein_atoms:
        return {"has_clashes": False, "clash_count": 0, "min_distance": float('inf')}

    clashes = 0
    min_dist = float('inf')
    for px, py, pz in protein_atoms:
        for lx, ly, lz in ligand_atoms:
            dist = math.sqrt((px-lx)**2 + (py-ly)**2 + (pz-lz)**2)
            min_dist = min(min_dist, dist)
            if dist < threshold:
                clashes += 1

    return {
        "has_clashes": clashes > 0,
        "clash_count": clashes,
        "min_distance": min_dist
    }


async def test_enzyme_scaffold_design():
    """Test the Enzyme Scaffold Design approach (same as UI)."""
    print("=" * 60)
    print("ENZYME SCAFFOLD DESIGN (Correct Approach)")
    print("=" * 60)
    print("Using: length + unindex + ligand (NOT contig with blockers)")
    print()

    # Extract theozyme from crystal structure
    workflow = ScaffoldingWorkflow()
    scaffold = await workflow.run(
        pdb_id='4CVB',
        ligand_code='pqq',
        metal='CA',
        include_all_ligand_contacts=True,
        fixed_atom_type='BKBN',  # Fix backbone of catalytic residues
    )

    if not scaffold.success:
        print(f"[FAIL] Scaffold extraction failed: {scaffold.error_message}")
        return None

    print(f"[OK] Theozyme extracted from 4CVB")
    print(f"     Length: {scaffold.length}")
    print(f"     Unindex: {scaffold.unindex}")
    print(f"     Ligand codes: {scaffold.ligand_codes}")
    print(f"     Fixed atoms: {scaffold.fixed_atoms}")
    print(f"     Catalytic residues: {len(scaffold.coordinating_residues)}")

    # Run RFD3 with Enzyme Scaffold Design parameters
    num_designs = 3
    print(f"\nRunning RFD3 ({num_designs} designs)...")

    rfd3_result = run_rfd3_api(
        # Enzyme Scaffold Design parameters
        length=scaffold.length,
        unindex=scaffold.unindex,
        ligand=scaffold.ligand_codes,
        pdb_content=scaffold.motif_pdb,
        select_fixed_atoms=scaffold.fixed_atoms,
        # Quality settings
        num_designs=num_designs,
        is_non_loopy=True,
        cfg_scale=2.5,
        # Optional: bury the metal
        select_buried=scaffold.rasa_targets if scaffold.rasa_targets else None,
        timeout=600
    )

    if rfd3_result["status"] != "completed":
        print(f"[FAIL] RFD3 failed: {rfd3_result.get('error')}")
        return None

    designs = rfd3_result["result"].get("designs", [])
    print(f"[OK] Generated {len(designs)} designs")

    # Analyze designs for clashes
    import os
    os.makedirs("design_results", exist_ok=True)

    passed = 0
    for i, design in enumerate(designs):
        pdb = design.get("content") or design.get("pdb_content", "")
        if not pdb:
            print(f"  Design {i+1}: No PDB content")
            continue

        # Check for PQQ clashes
        clashes = check_ligand_clashes(pdb, "PQQ")
        # Also check for CA clashes
        ca_clashes = check_ligand_clashes(pdb, "CA")

        status = "CLASH" if clashes["has_clashes"] else "PASS"
        print(f"  Design {i+1}: {status}")
        print(f"     PQQ: min_dist={clashes['min_distance']:.2f} A, clashes={clashes['clash_count']}")
        print(f"     CA:  min_dist={ca_clashes['min_distance']:.2f} A, clashes={ca_clashes['clash_count']}")

        if not clashes["has_clashes"]:
            passed += 1
            # Save passing design
            with open(f"design_results/enzyme_scaffold_{i+1}.pdb", "w") as f:
                f.write(pdb)
            print(f"     Saved to enzyme_scaffold_{i+1}.pdb")

    print(f"\nResult: {passed}/{len(designs)} passed ({100*passed/len(designs) if designs else 0:.0f}%)")
    return rfd3_result


if __name__ == "__main__":
    result = asyncio.run(test_enzyme_scaffold_design())

    print("\n" + "=" * 60)
    print("TEST COMPLETE")
    print("=" * 60)
    print()
    print("KEY INSIGHT: The Enzyme Scaffold Design approach uses:")
    print("  - length: Total scaffold length (e.g., '80-120')")
    print("  - unindex: Catalytic residues marked as FLEXIBLE")
    print("  - ligand: Comma-separated 3-letter codes")
    print("  - select_fixed_atoms: Fix ligand/metal position")
    print()
    print("This is different from the blocker-based approach which tried to")
    print("design linkers between fixed crystal structure residues.")
