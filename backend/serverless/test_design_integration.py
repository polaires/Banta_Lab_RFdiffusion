# test_design_integration.py
"""
Integration tests for the unified design system.

Run with: pytest test_design_integration.py -v -m integration

Unit tests can run locally without runpod:
    pytest test_design_integration.py -v -k "not test_unified_design and not test_zinc_design"
"""

import os
import pytest
import requests
from pathlib import Path
from collections import Counter

API_URL = os.getenv("DESIGN_API_URL", "http://localhost:8000/runsync")
INTEGRATION_TEST_TIMEOUT = 300  # seconds

pytestmark = pytest.mark.integration

# Check if handler can be imported (requires runpod)
try:
    from handler import infer_design_type_from_request
    HANDLER_AVAILABLE = True
except ImportError:
    HANDLER_AVAILABLE = False
    infer_design_type_from_request = None


@pytest.fixture
def zinc_pdb():
    """Load zinc dimer test PDB."""
    paths = [
        Path(__file__).parent.parent.parent / "design_zn_dimer.pdb",
        Path(__file__).parent.parent.parent / "azobenzene_dimer_design.pdb",
        Path("/mnt/g/Github_local_repo/Banta_Lab_RFdiffusion/design_zn_dimer.pdb"),
        Path("/mnt/g/Github_local_repo/Banta_Lab_RFdiffusion/azobenzene_dimer_design.pdb"),
    ]
    for p in paths:
        if p.exists():
            return p.read_text()
    pytest.skip("Test PDB not found")


def parse_sequences(result):
    """Extract sequence strings from result."""
    sequences = []
    for seq in result.get("sequences", []):
        if isinstance(seq, dict):
            for line in seq.get("content", "").split("\n"):
                if not line.startswith(">"):
                    sequences.append(line.strip())
    return sequences


class TestMetalBindingDesign:
    """Tests for metal-binding design."""

    def test_unified_design_with_metal_uses_ligandmpnn(self, zinc_pdb):
        """Unified design with metal should use LigandMPNN."""
        payload = {
            "input": {
                "task": "design",
                "metal": "ZN",
                "pdb_content": zinc_pdb,
                "num_sequences": 2,
            }
        }

        response = requests.post(API_URL, json=payload, timeout=INTEGRATION_TEST_TIMEOUT)
        result = response.json()

        assert "error" not in result, f"API error: {result.get('error')}. Status: {response.status_code}"
        output = result.get("output", {})

        # Should have selected LigandMPNN
        assert output.get("workflow", {}).get("sequence_tool") == "ligand_mpnn"
        assert output.get("design_type") == "METAL_BINDING"

    def test_zinc_design_reduces_alanine(self, zinc_pdb):
        """Zinc design should have low alanine content."""
        payload = {
            "input": {
                "task": "design",
                "metal": "ZN",
                "pdb_content": zinc_pdb,
                "num_sequences": 4,
            }
        }

        response = requests.post(API_URL, json=payload, timeout=INTEGRATION_TEST_TIMEOUT)
        result = response.json()

        assert "error" not in result, f"API error: {result.get('error')}. Status: {response.status_code}"
        output = result.get("output", {})
        sequences = parse_sequences(output)

        if not sequences:
            pytest.skip("No sequences returned")

        combined = "".join(sequences)
        total = len(combined)
        ala_pct = (combined.count("A") / total * 100) if total > 0 else 0

        print(f"Alanine: {ala_pct:.1f}%")
        assert ala_pct < 25, f"Alanine too high: {ala_pct:.1f}%"

    def test_zinc_design_has_coordinating_residues(self, zinc_pdb):
        """Zinc design should have H/C/E/D residues."""
        payload = {
            "input": {
                "task": "design",
                "metal": "ZN",
                "pdb_content": zinc_pdb,
                "num_sequences": 4,
            }
        }

        response = requests.post(API_URL, json=payload, timeout=INTEGRATION_TEST_TIMEOUT)
        result = response.json()

        assert "error" not in result, f"API error: {result.get('error')}. Status: {response.status_code}"
        output = result.get("output", {})
        sequences = parse_sequences(output)

        if not sequences:
            pytest.skip("No sequences returned")

        combined = "".join(sequences)
        aa_counts = Counter(combined)

        coordinating = sum(aa_counts.get(aa, 0) for aa in "HCED")
        assert coordinating > 0, "Should have coordinating residues"


@pytest.mark.skipif(not HANDLER_AVAILABLE, reason="handler requires runpod module")
class TestDesignTypeInference:
    """Tests for design type inference logic."""

    def test_protein_binder_uses_proteinmpnn(self):
        """Protein binder should use ProteinMPNN."""
        from design_types import DesignType

        request = {"target_pdb": "..."}
        design_type = infer_design_type_from_request(request)

        assert design_type == DesignType.PROTEIN_PROTEIN_BINDER

    def test_metal_with_symmetry_creates_dimer(self):
        """Metal + symmetry should create metal-mediated dimer."""
        from design_types import DesignType

        request = {"metal": "ZN", "symmetry": "C2"}
        design_type = infer_design_type_from_request(request)

        assert design_type == DesignType.METAL_MEDIATED_DIMER

    def test_ligand_creates_small_molecule_binder(self):
        """Ligand parameter should create small molecule binder."""
        from design_types import DesignType

        request = {"ligand": "AZO"}
        design_type = infer_design_type_from_request(request)

        assert design_type == DesignType.SMALL_MOLECULE_BINDER

    def test_nucleotide_creates_dna_rna_binder(self):
        """DNA/RNA parameter should create nucleotide binder."""
        from design_types import DesignType

        request = {"dna": True}
        design_type = infer_design_type_from_request(request)

        assert design_type == DesignType.DNA_RNA_BINDER


class TestConfigOverrides:
    """Tests for config override handling."""

    def test_temperature_override_applied(self):
        """Temperature override should be applied."""
        from design_orchestrator import DesignOrchestrator
        from design_types import DesignType

        orchestrator = DesignOrchestrator(
            design_type=DesignType.METAL_BINDING,
            temperature=0.5
        )

        assert orchestrator.config.temperature == 0.5

    def test_metal_preset_applied(self):
        """Metal-specific preset should be applied."""
        from design_orchestrator import DesignOrchestrator
        from design_types import DesignType

        orchestrator = DesignOrchestrator(
            design_type=DesignType.METAL_BINDING,
            metal_type="lanthanide"
        )

        # Lanthanide preset should penalize cysteine (soft donor incompatible)
        bias_aa = orchestrator.config.bias_AA or ""
        assert "C:-" in bias_aa, f"Lanthanide should penalize cysteine, got bias_AA: {bias_aa}"


class TestMetalCoordinationDetection:
    """Tests for metal coordination site detection."""

    def test_detect_zinc_coordination(self):
        """Should detect zinc coordination sites."""
        from inference_utils import detect_coordinating_residues

        pdb_content = """
ATOM      1  N   HIS A  63      10.000  10.000  10.000  1.00 50.00           N
ATOM      2  CA  HIS A  63      11.000  10.000  10.000  1.00 50.00           C
ATOM      3  NE2 HIS A  63      12.500  10.000   7.230  1.00 50.00           N
ATOM     50  N   CYS B  39      10.000  10.000   0.000  1.00 50.00           N
ATOM     51  CA  CYS B  39      11.000  10.000   0.000  1.00 50.00           C
ATOM     52  SG  CYS B  39      12.300  10.000   3.440  1.00 50.00           S
HETATM  100  ZN  ZN  X   1      12.000  10.000   5.000  1.00 50.00          ZN
END
"""
        sites = detect_coordinating_residues(pdb_content, cutoff=4.0)

        assert len(sites) == 1
        assert sites[0].metal_code == "ZN"

        fixed = sites[0].get_fixed_positions_list()
        assert "A63" in fixed
        assert "B39" in fixed


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
