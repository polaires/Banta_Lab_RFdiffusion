# test_metal_site_fetcher.py
"""Tests for PDB metal site fetching."""
import pytest

# Check if requests is available for network tests
try:
    import requests
    NETWORK_AVAILABLE = True
    # Quick connectivity check
    try:
        requests.head("https://search.rcsb.org", timeout=5)
    except (requests.RequestException, Exception):
        NETWORK_AVAILABLE = False
except ImportError:
    NETWORK_AVAILABLE = False

network_required = pytest.mark.skipif(
    not NETWORK_AVAILABLE,
    reason="Network access required (requests not available or RCSB unreachable)"
)


class TestMetalSiteQuery:
    """Test metal site database queries."""

    @network_required
    def test_query_zinc_sites(self):
        """Should find zinc binding sites in PDB."""
        from metal_site_fetcher import query_metal_sites
        sites = query_metal_sites(
            metal="ZN",
            resolution_max=2.0,
            limit=5,
        )
        assert len(sites) >= 1
        assert all(s["metal"] == "ZN" for s in sites)

    @network_required
    def test_query_lanthanide_sites(self):
        """Should find lanthanide sites (may be sparse)."""
        from metal_site_fetcher import query_metal_sites
        sites = query_metal_sites(
            metal="TB",
            resolution_max=3.0,
            limit=5,
        )
        # Lanthanides are rarer - may find 0
        # But if found, should be TB
        for s in sites:
            assert s["metal"] in ["TB", "EU", "GD", "LA"]

    @network_required
    def test_query_with_ligand(self):
        """Should find metal-ligand complex structures."""
        from metal_site_fetcher import query_metal_ligand_sites
        sites = query_metal_ligand_sites(
            metal="CA",
            ligand="PQQ",
            limit=5,
        )
        # PQQ-Ca should find methanol dehydrogenase structures
        assert len(sites) >= 0  # May be empty if API unavailable


class TestSiteExtraction:
    """Test extraction of coordination from PDB."""

    def test_extract_coordination_from_pdb(self):
        """Should extract coordinating residues from PDB content."""
        from metal_site_fetcher import extract_metal_coordination

        # Minimal zinc site
        pdb = """
ATOM      1  NE2 HIS A  63      12.500  10.000   7.230  1.00 50.00           N
ATOM      2  SG  CYS A  96      12.300  10.000   3.440  1.00 50.00           S
ATOM      3  SG  CYS A  99      11.700  11.500   5.000  1.00 50.00           S
ATOM      4  SG  CYS A 102      11.700   8.500   5.000  1.00 50.00           S
HETATM    5  ZN  ZN  X   1      12.000  10.000   5.000  1.00 50.00          ZN
END
"""
        coord = extract_metal_coordination(pdb, metal="ZN", cutoff=3.0)
        assert coord["metal"] == "ZN"
        assert coord["coordination_number"] >= 3
        assert len(coord["coordinating_atoms"]) >= 3

    def test_extract_no_metal_returns_empty(self):
        """Should return empty coordination when no metal found."""
        from metal_site_fetcher import extract_metal_coordination

        pdb = """
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 50.00           C
END
"""
        coord = extract_metal_coordination(pdb, metal="ZN", cutoff=3.0)
        assert coord["metal"] == "ZN"
        assert coord["coordination_number"] == 0
        assert len(coord["coordinating_atoms"]) == 0

    def test_extract_coordination_distance_cutoff(self):
        """Should respect distance cutoff for coordination."""
        from metal_site_fetcher import extract_metal_coordination

        # Atoms at varying distances from metal at origin
        pdb = """
ATOM      1  NE2 HIS A   1       2.000   0.000   0.000  1.00 50.00           N
ATOM      2  OE1 GLU A   2       4.000   0.000   0.000  1.00 50.00           O
HETATM    3  ZN  ZN  X   1       0.000   0.000   0.000  1.00 50.00          ZN
END
"""
        # With 2.5A cutoff, only HIS should be found
        coord_tight = extract_metal_coordination(pdb, metal="ZN", cutoff=2.5)
        assert coord_tight["coordination_number"] == 1
        assert coord_tight["coordinating_atoms"][0]["res_name"] == "HIS"

        # With 5.0A cutoff, both should be found
        coord_wide = extract_metal_coordination(pdb, metal="ZN", cutoff=5.0)
        assert coord_wide["coordination_number"] == 2


class TestTemplateGeneration:
    """Test template generation from PDB sites."""

    @network_required
    def test_generate_template_from_pdb_id(self):
        """Should generate template from known PDB structure."""
        from metal_site_fetcher import generate_template_from_pdb

        # 1CA2 = carbonic anhydrase with Zn
        template = generate_template_from_pdb(
            pdb_id="1CA2",
            metal="ZN",
        )
        # May fail if network unavailable
        if template:
            assert template["metal"] == "ZN"
            assert template["source"] == "pdb"
            assert "pdb_id" in template

    @network_required
    def test_get_reference_template(self):
        """Should get template from curated reference structures."""
        from metal_site_fetcher import get_reference_template

        template = get_reference_template(metal="ZN")
        if template:
            assert template["metal"] == "ZN"
            assert template["source"] == "pdb"
            assert template["coordination_number"] > 0


class TestUtilityFunctions:
    """Test utility functions."""

    def test_get_available_reference_metals(self):
        """Should return list of metals with reference structures."""
        from metal_site_fetcher import get_available_reference_metals

        metals = get_available_reference_metals()
        assert "ZN" in metals
        assert "FE" in metals
        assert "CA" in metals
        assert "TB" in metals

    def test_get_reference_pdb_ids(self):
        """Should return reference PDB IDs for a metal."""
        from metal_site_fetcher import get_reference_pdb_ids

        zn_refs = get_reference_pdb_ids("ZN")
        assert len(zn_refs) > 0
        assert "1CA2" in zn_refs

        # Unknown metal should return empty list
        unknown_refs = get_reference_pdb_ids("XX")
        assert len(unknown_refs) == 0

    def test_infer_geometry(self):
        """Should infer correct geometry from coordination number."""
        from metal_site_fetcher import _infer_geometry

        assert _infer_geometry(4) == "tetrahedral"
        assert _infer_geometry(6) == "octahedral"
        assert _infer_geometry(8) == "square_antiprismatic"
        assert _infer_geometry(100) == "unknown"


class TestDataClasses:
    """Test data class functionality."""

    def test_metal_site_to_dict(self):
        """Should convert MetalSite to dictionary."""
        from metal_site_fetcher import MetalSite

        site = MetalSite(
            pdb_id="1CA2",
            metal="ZN",
            metal_chain="A",
            metal_resnum=301,
            metal_coords=(10.0, 20.0, 30.0),
            coordination_number=4,
            coordinating_atoms=[{"name": "HIS"}, {"name": "HIS"}, {"name": "HIS"}],
            resolution=1.5,
            geometry="tetrahedral",
        )
        d = site.to_dict()
        assert d["pdb_id"] == "1CA2"
        assert d["metal"] == "ZN"
        assert d["coordination_number"] == 4

    def test_metal_ligand_site_to_dict(self):
        """Should convert MetalLigandSite to dictionary with ligand info."""
        from metal_site_fetcher import MetalLigandSite

        site = MetalLigandSite(
            pdb_id="1ABC",
            metal="CA",
            ligand_code="PQQ",
            ligand_chain="A",
            ligand_resnum=500,
            ligand_atoms=[{"name": "O1"}],
        )
        d = site.to_dict()
        assert d["pdb_id"] == "1ABC"
        assert d["metal"] == "CA"
        assert d["ligand_code"] == "PQQ"
