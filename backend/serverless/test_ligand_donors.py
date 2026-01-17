# test_ligand_donors.py
"""Tests for ligand donor atom identification."""
import pytest

# Check if RDKit is available
try:
    from rdkit import Chem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

# Skip marker for tests requiring RDKit
requires_rdkit = pytest.mark.skipif(
    not HAS_RDKIT,
    reason="RDKit not available - ligand analysis requires RDKit"
)


@requires_rdkit
class TestDonorIdentification:
    """Test identification of coordinating atoms from ligand."""

    def test_citrate_donors(self):
        """Citrate should have 6 potential O donors (3 carboxylates)."""
        from ligand_donors import identify_donors_from_smiles

        citrate_smiles = "C(C(=O)O)(CC(=O)O)(C(=O)O)O"
        donors = identify_donors_from_smiles(citrate_smiles)

        # 3 carboxylates = 6 oxygen atoms potentially coordinating
        oxygen_donors = [d for d in donors if d["element"] == "O"]
        assert len(oxygen_donors) >= 4  # At least 4 accessible

    def test_pqq_donors(self):
        """PQQ should have multiple O and N donors."""
        from ligand_donors import identify_donors_from_smiles

        # Simplified PQQ core
        pqq_smiles = "O=C1C(=O)c2cc3C(=O)NC(=O)c3c(c2N1)C(=O)O"
        donors = identify_donors_from_smiles(pqq_smiles)

        # Should have O and N donors
        elements = set(d["element"] for d in donors)
        assert "O" in elements
        assert "N" in elements

    def test_donor_type_classification(self):
        """Should classify donor types correctly."""
        from ligand_donors import identify_donors_from_smiles

        # Carboxylic acid
        acetic = "CC(=O)O"
        donors = identify_donors_from_smiles(acetic)

        carboxylate_donors = [d for d in donors if d["type"] == "carboxylate"]
        assert len(carboxylate_donors) >= 1


@requires_rdkit
class TestMetalCompatibility:
    """Test metal-donor compatibility scoring."""

    def test_citrate_lanthanide_compatibility(self):
        """Citrate should be highly compatible with lanthanides."""
        from ligand_donors import score_ligand_metal_compatibility

        citrate_smiles = "C(C(=O)O)(CC(=O)O)(C(=O)O)O"
        score = score_ligand_metal_compatibility(citrate_smiles, "TB")

        assert score >= 0.8  # High compatibility

    def test_thiol_lanthanide_incompatibility(self):
        """Thiol ligands should be incompatible with lanthanides."""
        from ligand_donors import score_ligand_metal_compatibility

        thiol_smiles = "CCCS"  # Propanethiol
        score = score_ligand_metal_compatibility(thiol_smiles, "TB")

        assert score < 0.3  # Low compatibility


@requires_rdkit
class TestDenticity:
    """Test denticity determination."""

    def test_citrate_bidentate_carboxylates(self):
        """Citrate carboxylates can be bidentate."""
        from ligand_donors import analyze_denticity

        citrate_smiles = "C(C(=O)O)(CC(=O)O)(C(=O)O)O"
        denticity = analyze_denticity(citrate_smiles)

        # Should have potential for bidentate
        bidentate_groups = [d for d in denticity if d["max_denticity"] >= 2]
        assert len(bidentate_groups) >= 1


class TestGracefulDegradation:
    """Test that functions degrade gracefully without RDKit."""

    def test_identify_donors_no_rdkit(self):
        """identify_donors_from_smiles should return empty list if RDKit unavailable."""
        from ligand_donors import identify_donors_from_smiles, HAS_RDKIT as MODULE_HAS_RDKIT

        if MODULE_HAS_RDKIT:
            pytest.skip("RDKit is available, skipping graceful degradation test")

        donors = identify_donors_from_smiles("CC(=O)O")
        assert donors == []

    def test_score_compatibility_no_rdkit(self):
        """score_ligand_metal_compatibility should return 0.0 if RDKit unavailable."""
        from ligand_donors import score_ligand_metal_compatibility, HAS_RDKIT as MODULE_HAS_RDKIT

        if MODULE_HAS_RDKIT:
            pytest.skip("RDKit is available, skipping graceful degradation test")

        score = score_ligand_metal_compatibility("CC(=O)O", "TB")
        assert score == 0.0

    def test_analyze_denticity_no_rdkit(self):
        """analyze_denticity should return empty list if RDKit unavailable."""
        from ligand_donors import analyze_denticity, HAS_RDKIT as MODULE_HAS_RDKIT

        if MODULE_HAS_RDKIT:
            pytest.skip("RDKit is available, skipping graceful degradation test")

        denticity = analyze_denticity("CC(=O)O")
        assert denticity == []
