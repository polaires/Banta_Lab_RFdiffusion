# test_metal_coordination.py
"""Tests for metal coordination site detection."""

import pytest


ZINC_PDB = """
ATOM      1  N   HIS A  63      10.000  10.000  10.000  1.00 50.00           N
ATOM      2  CA  HIS A  63      11.000  10.000  10.000  1.00 50.00           C
ATOM      3  NE2 HIS A  63      12.500  10.000   7.230  1.00 50.00           N
ATOM     50  N   CYS B  39      10.000  10.000   0.000  1.00 50.00           N
ATOM     51  CA  CYS B  39      11.000  10.000   0.000  1.00 50.00           C
ATOM     52  SG  CYS B  39      12.300  10.000   3.440  1.00 50.00           S
HETATM  100  ZN  ZN  X   1      12.000  10.000   5.000  1.00 50.00          ZN
END
"""


def test_detect_coordinating_residues():
    """Should detect His and Cys coordinating Zn."""
    from inference_utils import detect_coordinating_residues

    sites = detect_coordinating_residues(ZINC_PDB, cutoff=4.0)

    assert len(sites) == 1
    site = sites[0]
    assert site.metal_code == "ZN"

    # Should find both coordinating residues
    chains_residues = [(r.chain, r.resnum) for r in site.coordinating_residues]
    assert ("A", 63) in chains_residues
    assert ("B", 39) in chains_residues


def test_get_fixed_positions_string():
    """Should generate fixed_positions for LigandMPNN."""
    from inference_utils import detect_coordinating_residues

    sites = detect_coordinating_residues(ZINC_PDB, cutoff=4.0)
    fixed = sites[0].get_fixed_positions_list()

    assert "A63" in fixed
    assert "B39" in fixed


def test_detect_coordinating_residues_empty_pdb():
    """Should return empty list for empty PDB."""
    from inference_utils import detect_coordinating_residues

    sites = detect_coordinating_residues("", cutoff=4.0)
    assert sites == []


def test_detect_coordinating_residues_no_metal():
    """Should return empty list when no metals present."""
    from inference_utils import detect_coordinating_residues

    pdb_no_metal = """
ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 50.00           N
ATOM      2  CA  ALA A   1      11.000  10.000  10.000  1.00 50.00           C
END
"""
    sites = detect_coordinating_residues(pdb_no_metal, cutoff=4.0)
    assert sites == []


def test_detect_coordinating_residues_invalid_cutoff():
    """Should raise ValueError for invalid cutoff."""
    from inference_utils import detect_coordinating_residues

    with pytest.raises(ValueError, match="cutoff must be positive"):
        detect_coordinating_residues(ZINC_PDB, cutoff=0)

    with pytest.raises(ValueError, match="cutoff must be positive"):
        detect_coordinating_residues(ZINC_PDB, cutoff=-1.0)
