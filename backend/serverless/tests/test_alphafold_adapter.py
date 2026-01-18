# test_alphafold_adapter.py
"""Tests for AlphaFold database adapter."""
import pytest
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Check if requests is available for network tests
try:
    import requests
    NETWORK_AVAILABLE = True
    # Quick connectivity check for AlphaFold
    try:
        requests.head("https://alphafold.ebi.ac.uk", timeout=5)
        ALPHAFOLD_AVAILABLE = True
    except (requests.RequestException, Exception):
        ALPHAFOLD_AVAILABLE = False
except ImportError:
    NETWORK_AVAILABLE = False
    ALPHAFOLD_AVAILABLE = False

network_required = pytest.mark.skipif(
    not NETWORK_AVAILABLE,
    reason="Network access required (requests not available)"
)

alphafold_required = pytest.mark.skipif(
    not ALPHAFOLD_AVAILABLE,
    reason="AlphaFold API not reachable"
)


class TestAlphaFoldAdapterImport:
    """Test AlphaFold adapter can be imported."""

    def test_can_import_adapter(self):
        """Should be able to import AlphaFoldAdapter."""
        from database_adapters import AlphaFoldAdapter
        assert AlphaFoldAdapter is not None

    def test_can_import_from_module(self):
        """Should be able to import from alphafold_adapter module."""
        from database_adapters.alphafold_adapter import AlphaFoldAdapter
        assert AlphaFoldAdapter is not None

    def test_can_import_constants(self):
        """Should be able to import constants."""
        from database_adapters.alphafold_adapter import (
            ALPHAFOLD_BASE_URL,
            ALPHAFOLD_FILES_URL,
            DEFAULT_TIMEOUT,
        )
        assert ALPHAFOLD_BASE_URL is not None
        assert ALPHAFOLD_FILES_URL is not None
        assert DEFAULT_TIMEOUT == 30


class TestAlphaFoldAdapterInitialization:
    """Test AlphaFoldAdapter initialization."""

    def test_adapter_has_default_base_url(self):
        """Should have AlphaFold API base URL by default."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        assert "alphafold.ebi.ac.uk" in adapter.base_url

    def test_adapter_has_default_files_url(self):
        """Should have AlphaFold files URL by default."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        assert "alphafold.ebi.ac.uk" in adapter.files_url

    def test_adapter_can_override_base_url(self):
        """Should allow overriding base URL."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter(base_url="https://custom.example.com")
        assert adapter.base_url == "https://custom.example.com"

    def test_adapter_has_timeout_setting(self):
        """Should have configurable timeout."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter(timeout=60)
        assert adapter.timeout == 60


class TestGetPrediction:
    """Test get_prediction method."""

    def test_get_prediction_returns_dict_or_none(self):
        """Should return a dict or None from get_prediction."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        result = adapter.get_prediction("P00918")
        assert result is None or isinstance(result, dict)

    @network_required
    @alphafold_required
    def test_get_prediction_returns_metadata(self):
        """P00918 (human carbonic anhydrase) should return prediction metadata."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        result = adapter.get_prediction("P00918")

        # Should have prediction for this well-known protein
        assert result is not None, "P00918 should have AlphaFold prediction"

        # Should have key metadata fields
        assert "pdb_url" in result
        assert "model_version" in result
        assert "uniprot_id" in result

        # Model version should be a positive integer
        assert isinstance(result["model_version"], int)
        assert result["model_version"] > 0

    @network_required
    @alphafold_required
    def test_get_prediction_has_structure_urls(self):
        """Prediction should include PDB and/or CIF URLs."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        result = adapter.get_prediction("P00918")

        if result:
            # Should have at least one structure URL
            has_pdb = result.get("pdb_url") and len(result["pdb_url"]) > 0
            has_cif = result.get("cif_url") and len(result["cif_url"]) > 0
            assert has_pdb or has_cif, "Should have PDB or CIF URL"

    def test_get_prediction_empty_id_returns_none(self):
        """Empty UniProt ID should return None."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        result = adapter.get_prediction("")
        assert result is None

    def test_get_prediction_whitespace_id_returns_none(self):
        """Whitespace-only UniProt ID should return None."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        result = adapter.get_prediction("   ")
        assert result is None

    @network_required
    @alphafold_required
    def test_nonexistent_returns_none(self):
        """Invalid/nonexistent UniProt ID should return None."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        result = adapter.get_prediction("INVALID12345XYZ")
        assert result is None


class TestDownloadStructure:
    """Test download_structure method."""

    def test_download_structure_returns_string_or_none(self):
        """Should return a string or None from download_structure."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        result = adapter.download_structure("P00918")
        assert result is None or isinstance(result, str)

    @network_required
    @alphafold_required
    def test_download_structure_pdb(self):
        """Download PDB format should return content starting with PDB header."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        result = adapter.download_structure("P00918", format="pdb")

        assert result is not None, "P00918 should have AlphaFold structure"
        assert len(result) > 0

        # PDB file should start with HEADER or contain ATOM records
        lines = result.strip().split("\n")
        first_line = lines[0] if lines else ""
        has_header = first_line.startswith("HEADER")
        has_atom = any(line.startswith("ATOM") for line in lines[:100])

        assert has_header or has_atom, "PDB content should have HEADER or ATOM records"

    @network_required
    @alphafold_required
    def test_download_structure_cif(self):
        """Download CIF format should return valid mmCIF content."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        result = adapter.download_structure("P00918", format="cif")

        assert result is not None, "P00918 should have AlphaFold structure in CIF format"
        assert len(result) > 0

        # CIF file should contain data_ block
        assert "data_" in result or "_atom_site" in result, "CIF should have data_ block or _atom_site"

    @network_required
    @alphafold_required
    def test_download_structure_default_is_pdb(self):
        """Default format should be PDB."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        result = adapter.download_structure("P00918")

        if result:
            # Default should be PDB format
            lines = result.strip().split("\n")
            has_pdb_content = any(
                line.startswith(("HEADER", "ATOM", "HETATM", "TER", "END"))
                for line in lines[:100]
            )
            assert has_pdb_content, "Default format should be PDB"

    def test_download_structure_invalid_format(self):
        """Invalid format should return None."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        result = adapter.download_structure("P00918", format="xyz")
        assert result is None

    def test_download_structure_empty_id_returns_none(self):
        """Empty UniProt ID should return None."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        result = adapter.download_structure("")
        assert result is None

    @network_required
    @alphafold_required
    def test_download_structure_nonexistent_returns_none(self):
        """Nonexistent UniProt ID should return None."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        result = adapter.download_structure("INVALID12345XYZ")
        assert result is None


class TestGetPlddtScores:
    """Test get_plddt_scores method."""

    def test_get_plddt_scores_returns_list_or_none(self):
        """Should return a list or None from get_plddt_scores."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        result = adapter.get_plddt_scores("P00918")
        assert result is None or isinstance(result, list)

    @network_required
    @alphafold_required
    def test_get_plddt_scores(self):
        """P00918 should return pLDDT scores as list of floats."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        result = adapter.get_plddt_scores("P00918")

        assert result is not None, "P00918 should have pLDDT scores"
        assert len(result) > 0, "Should have at least one score"

        # All scores should be floats in valid range (0-100)
        for score in result:
            assert isinstance(score, (int, float))
            assert 0 <= score <= 100, f"pLDDT score {score} out of range"

    @network_required
    @alphafold_required
    def test_get_plddt_scores_length_matches_sequence(self):
        """pLDDT scores count should approximately match protein length."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        scores = adapter.get_plddt_scores("P00918")

        if scores:
            # Human carbonic anhydrase II is ~260 amino acids
            assert 200 < len(scores) < 300, f"Unexpected score count: {len(scores)}"

    def test_get_plddt_scores_empty_id_returns_none(self):
        """Empty UniProt ID should return None."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        result = adapter.get_plddt_scores("")
        assert result is None


class TestCheckAvailability:
    """Test check_availability method."""

    def test_check_availability_returns_bool(self):
        """Should return a boolean from check_availability."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        result = adapter.check_availability("P00918")
        assert isinstance(result, bool)

    @network_required
    @alphafold_required
    def test_check_availability(self):
        """P00918 should return True (prediction exists)."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        result = adapter.check_availability("P00918")

        assert result is True, "P00918 should have AlphaFold prediction"

    @network_required
    @alphafold_required
    def test_check_availability_nonexistent(self):
        """Invalid UniProt ID should return False."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        result = adapter.check_availability("INVALID12345XYZ")

        assert result is False, "Invalid ID should return False"

    def test_check_availability_empty_id_returns_false(self):
        """Empty UniProt ID should return False."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        result = adapter.check_availability("")
        assert result is False

    def test_check_availability_whitespace_id_returns_false(self):
        """Whitespace-only UniProt ID should return False."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()
        result = adapter.check_availability("   ")
        assert result is False


class TestErrorHandling:
    """Test error handling and logging."""

    def test_network_error_get_prediction_returns_none(self):
        """Network errors on get_prediction should return None."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter(
            base_url="https://invalid.nonexistent.domain.example",
            timeout=1,
        )
        result = adapter.get_prediction("P00918")
        assert result is None

    def test_network_error_download_structure_returns_none(self):
        """Network errors on download_structure should return None."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter(
            files_url="https://invalid.nonexistent.domain.example",
            timeout=1,
        )
        result = adapter.download_structure("P00918")
        assert result is None

    def test_network_error_check_availability_returns_false(self):
        """Network errors on check_availability should return False."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter(
            base_url="https://invalid.nonexistent.domain.example",
            timeout=1,
        )
        result = adapter.check_availability("P00918")
        assert result is False


class TestIntegration:
    """Integration tests requiring network access."""

    @pytest.mark.integration
    @network_required
    @alphafold_required
    def test_full_workflow_get_prediction_and_structure(self):
        """Full workflow: get prediction metadata, then download structure."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()

        # Step 1: Check availability
        available = adapter.check_availability("P00918")
        assert available, "P00918 should be available"

        # Step 2: Get prediction metadata
        prediction = adapter.get_prediction("P00918")
        assert prediction is not None
        assert "pdb_url" in prediction

        # Step 3: Download structure
        structure = adapter.download_structure("P00918", format="pdb")
        assert structure is not None
        assert len(structure) > 0

        # Step 4: Get pLDDT scores
        scores = adapter.get_plddt_scores("P00918")
        assert scores is not None
        assert len(scores) > 0

    @pytest.mark.integration
    @network_required
    @alphafold_required
    def test_multiple_proteins(self):
        """Test multiple well-known proteins."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()

        # Test several well-characterized proteins
        test_proteins = [
            "P00918",  # Human carbonic anhydrase II
            "P04637",  # Human p53
            "P68871",  # Hemoglobin beta
        ]

        for uniprot_id in test_proteins:
            available = adapter.check_availability(uniprot_id)
            # Most SwissProt proteins should have predictions
            if available:
                prediction = adapter.get_prediction(uniprot_id)
                assert prediction is not None, f"Failed to get prediction for {uniprot_id}"

    @pytest.mark.integration
    @network_required
    @alphafold_required
    def test_both_formats(self):
        """Test downloading both PDB and CIF formats."""
        from database_adapters import AlphaFoldAdapter

        adapter = AlphaFoldAdapter()

        # Download PDB
        pdb_content = adapter.download_structure("P00918", format="pdb")
        assert pdb_content is not None

        # Download CIF
        cif_content = adapter.download_structure("P00918", format="cif")
        assert cif_content is not None

        # Both should have content
        assert len(pdb_content) > 0
        assert len(cif_content) > 0

        # They should be different formats
        assert pdb_content != cif_content
