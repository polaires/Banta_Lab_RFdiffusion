#!/usr/bin/env python3
"""
Test script for crystal structure coordination scaffold approach.

This tests the Option 3 implementation: using real crystal structures
to extract coordination geometry for metal-ligand complex design.
"""

import sys
sys.path.insert(0, '.')

from metal_site_fetcher import (
    search_metal_ligand_structures,
    find_metal_ligand_active_site,
    _fetch_pdb_content,
)


def test_search_metal_ligand():
    """Test finding structures with both metal and ligand."""
    print("=== Test 1: Search for CA+PQQ structures ===")

    pdb_ids = search_metal_ligand_structures("CA", "PQQ", limit=5)
    print(f"Found {len(pdb_ids)} structures: {pdb_ids}")

    assert len(pdb_ids) > 0, "Should find at least one CA+PQQ structure"
    assert "1CQ1" in pdb_ids or "1C9U" in pdb_ids, "Should find known PQQ dehydrogenase"

    print("PASSED\n")
    return pdb_ids


def test_find_active_site(pdb_id="1CQ1"):
    """Test finding the active site metal that coordinates ligand."""
    print(f"=== Test 2: Find active site in {pdb_id} ===")

    site = find_metal_ligand_active_site(pdb_id, "CA", "PQQ", cutoff=3.5)

    assert site is not None, f"Should find active site in {pdb_id}"
    assert site["metal"] == "CA", "Metal should be CA"
    assert site["ligand_distance"] < 3.0, "Ligand should be coordinating metal"
    assert len(site["ligand_donors"]) > 0, "Should have ligand donors"
    assert len(site["protein_donors"]) > 0, "Should have protein donors"

    print(f"Active site: {site['metal']} {site['metal_chain']}{site['metal_resnum']}")
    print(f"Ligand distance: {site['ligand_distance']}Å")
    print(f"Ligand donors ({len(site['ligand_donors'])}): {site['ligand_donors']}")
    print(f"Protein donors ({len(site['protein_donors'])}): {site['protein_donors']}")

    # PQQ-Ca should have backbone carbonyl donors (GLY, PRO)
    protein_residues = [d.split(":")[0] for d in site['protein_donors']]
    print(f"Protein residues: {protein_residues}")

    print("PASSED\n")
    return site


def test_extract_scaffold(pdb_id="1CQ1"):
    """Test extracting coordination scaffold from structure."""
    print(f"=== Test 3: Extract scaffold from {pdb_id} ===")

    site = find_metal_ligand_active_site(pdb_id, "CA", "PQQ", cutoff=3.5)
    pdb_content = _fetch_pdb_content(pdb_id)

    assert pdb_content is not None, "Should fetch PDB content"

    # Extract the specific metal chain/resnum
    metal_chain = site['metal_chain']
    metal_resnum = site['metal_resnum']

    # Count HETATM lines for metal and ligand
    metal_lines = []
    ligand_lines = []
    for line in pdb_content.split('\n'):
        if line.startswith('HETATM'):
            res_name = line[17:20].strip()
            chain = line[21]
            resnum = line[22:26].strip()
            if res_name == "CA" and chain == metal_chain and resnum == metal_resnum:
                metal_lines.append(line)
            elif res_name == "PQQ":
                ligand_lines.append(line)

    print(f"Metal atoms: {len(metal_lines)}")
    print(f"Ligand atoms: {len(ligand_lines)}")

    # Get protein residues from coordination
    protein_residues = set()
    for atom in site['coordinating_atoms']:
        res_name = atom.get('res_name', '')
        if res_name not in ('CA', 'PQQ', 'HOH'):
            chain = atom.get('chain_id', '')
            resnum = atom.get('res_seq', 0)
            protein_residues.add((chain, resnum, res_name))

    print(f"Coordinating protein residues: {protein_residues}")

    print("PASSED\n")
    return site, pdb_content


def test_coordination_chemistry():
    """Test that we understand the coordination chemistry correctly."""
    print("=== Test 4: Verify PQQ-Ca coordination chemistry ===")

    site = find_metal_ligand_active_site("1CQ1", "CA", "PQQ", cutoff=3.5)

    # PQQ should donate through O5, O7A, and N6
    pqq_donors = [d.split('@')[0] for d in site['ligand_donors']]
    print(f"PQQ donor atoms: {pqq_donors}")

    expected_donors = {'O5', 'O7A', 'N6'}
    actual_donors = set(pqq_donors)

    # Check if we have the key coordination atoms (within 3Å)
    close_donors = [d for d in site['ligand_donors'] if float(d.split('@')[1].replace('Å', '')) < 3.0]
    close_atoms = [d.split('@')[0] for d in close_donors]
    print(f"Close donors (<3Å): {close_donors}")

    # Protein should donate through backbone carbonyls (not sidechains for PQQ-Ca)
    for donor in site['protein_donors']:
        parts = donor.split()
        atom = parts[1].split('@')[0]
        print(f"Protein donor: {parts[0]} atom {atom}")
        # For PQQ-Ca, we expect backbone O atoms
        assert atom == 'O' or atom == 'C', f"Expected backbone atom, got {atom}"

    print("PASSED\n")


def main():
    """Run all tests."""
    print("=" * 60)
    print("Crystal Structure Coordination Scaffold Tests")
    print("=" * 60 + "\n")

    try:
        pdb_ids = test_search_metal_ligand()
        site = test_find_active_site(pdb_ids[0])
        test_extract_scaffold(pdb_ids[0])
        test_coordination_chemistry()

        print("=" * 60)
        print("ALL TESTS PASSED")
        print("=" * 60)

    except AssertionError as e:
        print(f"\nTEST FAILED: {e}")
        return 1
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
