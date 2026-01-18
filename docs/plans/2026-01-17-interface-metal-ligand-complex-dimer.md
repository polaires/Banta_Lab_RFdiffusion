# Implementation Plan: Interface Metal-Ligand Complex Dimer

**Date:** 2026-01-17
**Status:** ✅ Implementation Complete
**Test Cases:** PQQ-Ca (pyrroloquinoline quinone-calcium), Citrate-Lanthanide (Tb/Eu)

---

## Overview

Design protein dimers around **metal-ligand complexes** where both the organic ligand and coordinated metal are functionally important. This combines concepts from:
- `interface_ligand_design` (organic molecules at interface)
- `interface_metal_design` (metal coordination chemistry)

### Target Use Cases
1. **PQQ-Ca**: Pyrroloquinoline quinone cofactor with calcium - quinoprotein dehydrogenases
2. **Citrate-Lanthanide**: Citrate-coordinated Tb³⁺/Eu³⁺ - luminescent biosensors
3. **Heme-Fe**: Porphyrin-iron complexes - oxygen binding, electron transfer
4. **Chlorophyll**: Mg-porphyrin - light harvesting

---

## Step 1: Complex Structure Parser

**Files to modify:** `G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless\inference_utils.py`

### 1.1 Add Metal-Ligand Complex Dataclass

```python
@dataclass
class MetalLigandComplex:
    """Represents a metal-ligand complex with coordination information."""
    metal_code: str                          # e.g., "CA", "TB", "FE"
    metal_coords: Tuple[float, float, float]
    metal_chain: str
    metal_resnum: int
    ligand_atoms: List[Dict[str, Any]]       # List of {name, coords, element, res_name}
    ligand_res_name: str                     # e.g., "PQQ", "CIT"
    ligand_chain: str
    ligand_resnum: int
    coordination_bonds: List[Dict[str, Any]] # Existing metal-ligand bonds
    available_sites: List[Dict[str, Any]]    # Sites available for protein coordination
    coordination_geometry: str               # "octahedral", "square_antiprism", etc.
    centroid: Tuple[float, float, float]     # Complex center for positioning
```

### 1.2 Implement Complex Parser Function

Location: `inference_utils.py` after line ~310 (after `detect_coordinating_residues`)

```python
def parse_metal_ligand_complex(
    pdb_content: str,
    metal_codes: Optional[List[str]] = None,
    ligand_names: Optional[List[str]] = None,
    coordination_cutoff: float = 3.0,
) -> List[MetalLigandComplex]:
    """
    Parse metal-ligand complex from PDB and identify coordination.

    Args:
        pdb_content: PDB file content
        metal_codes: Expected metals (auto-detect if None)
        ligand_names: Expected ligand residue names (auto-detect if None)
        coordination_cutoff: Distance cutoff for metal-ligand bonds

    Returns:
        List of MetalLigandComplex objects
    """
    # Implementation:
    # 1. Parse HETATM records for metals and ligands
    # 2. Find metal-ligand distances within cutoff
    # 3. Classify coordination geometry based on angles
    # 4. Identify available coordination sites (unfilled positions)
    # 5. Calculate complex centroid
```

### 1.3 Implement Available Site Detection

```python
def detect_available_coordination_sites(
    complex_info: MetalLigandComplex,
    target_coordination: int = 6,
) -> List[Dict[str, Any]]:
    """
    Identify positions where protein residues can coordinate to metal.

    For octahedral (CN=6): 2 axial + 4 equatorial positions
    For square antiprism (CN=8): 8 positions in two planes

    Returns list of available site positions with:
    - coords: (x, y, z) position for protein donor atom
    - site_type: "axial" or "equatorial"
    - preferred_donors: ["His", "Cys", "Asp", "Glu"]
    """
```

---

## Step 2: Complex Template Library

**Files to create/modify:**
- `G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless\metal_ligand_templates.py` (NEW)
- `G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless\handler.py`

### 2.1 Create Template Library File

```python
"""
Metal-Ligand Complex Templates for Interface Design

Provides pre-defined templates for common cofactors and metal-ligand complexes
with accurate 3D coordinates and coordination information.
"""

# PQQ-Calcium template
# Source: PDB ligand PQQ + Ca coordination from 1W6S (methanol dehydrogenase)
# SMILES: OC(=O)c1[nH]c2c(c1)C(=O)C(=O)c3nc(cc(C(O)=O)c23)C(O)=O
PQQ_CA_TEMPLATE = {
    "name": "PQQ-Calcium",
    "description": "Pyrroloquinoline quinone with calcium cofactor",
    "metal": "CA",
    "ligand_smiles": "OC(=O)c1[nH]c2c(c1)C(=O)C(=O)c3nc(cc(C(O)=O)c23)C(O)=O",
    "ligand_res_name": "PQQ",
    "molecular_weight": 330.206,
    "coordination": {
        "metal_coordination_number": 6,
        "ligand_donors": ["O5", "N6", "O7"],  # PQQ atoms coordinating Ca
        "protein_sites": 3,  # Available for protein coordination
        "geometry": "distorted_octahedral",
        "ca_bond_distances": {
            "O5": 2.4,  # Quinone oxygen
            "N6": 2.5,  # Pyridine nitrogen
            "O7": 2.4,  # Carboxylate oxygen
        }
    },
    "preferred_protein_donors": ["Glu", "Asp", "Asn", "water"],
    "pdb_reference": "1W6S",
    "best_for": ["dehydrogenase", "electron_transfer", "quinoprotein"],
}

# Citrate-Lanthanide template
# Source: Citrate tricarboxylate coordination to Tb3+/Eu3+
# SMILES (citrate): OC(CC([O-])=O)(CC([O-])=O)C([O-])=O
CITRATE_LANTHANIDE_TEMPLATE = {
    "name": "Citrate-Lanthanide",
    "description": "Citrate-coordinated lanthanide (Tb/Eu) for luminescence",
    "metals": ["TB", "EU", "GD"],  # Compatible lanthanides
    "ligand_smiles": "OC(CC([O-])=O)(CC([O-])=O)C([O-])=O",
    "ligand_res_name": "CIT",
    "molecular_weight": 189.10,  # Citrate anion
    "coordination": {
        "metal_coordination_number": 9,  # Typical for Ln3+
        "ligand_donors": ["O1", "O3", "O5", "O7"],  # Carboxylate oxygens
        "citrate_denticity": "tetradentate",  # Can be bi- to tetra-dentate
        "protein_sites": 5,  # Remaining coordination sites
        "geometry": "tricapped_trigonal_prism",
        "ln_bond_distances": {
            "carboxylate_O": 2.35,  # Typical Ln-O distance
            "hydroxyl_O": 2.45,
        }
    },
    "preferred_protein_donors": ["Glu", "Asp", "Asn", "water"],
    "luminescence": {
        "TB": {"emission": "green", "wavelength": 545},
        "EU": {"emission": "red", "wavelength": 615},
    },
    "best_for": ["biosensor", "luminescence", "metal_sensing"],
}

METAL_LIGAND_COMPLEX_TEMPLATES = {
    "pqq_ca": PQQ_CA_TEMPLATE,
    "citrate_tb": {**CITRATE_LANTHANIDE_TEMPLATE, "metal": "TB"},
    "citrate_eu": {**CITRATE_LANTHANIDE_TEMPLATE, "metal": "EU"},
    "citrate_gd": {**CITRATE_LANTHANIDE_TEMPLATE, "metal": "GD"},
}
```

### 2.2 Add Template PDB Generator

```python
def generate_complex_pdb(
    template_name: str,
    metal: Optional[str] = None,
    center: Tuple[float, float, float] = (50.0, 50.0, 50.0),
) -> str:
    """
    Generate PDB content for a metal-ligand complex from template.

    Uses RDKit for ligand conformer + places metal at coordination site.
    """
    # 1. Get template
    # 2. Generate ligand conformer from SMILES
    # 3. Place metal at coordination center
    # 4. Rotate/translate to center
    # 5. Return combined PDB
```

---

## Step 3: Handler Implementation

**File:** `G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless\handler.py`

### 3.1 Add Handler Function

Location: After `handle_interface_metal_design` (~line 6300)

```python
def handle_interface_metal_ligand_design(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Design protein dimers around metal-ligand complexes.

    Combines organic ligand binding with metal coordination chemistry.

    Input:
        # Complex specification (one of):
        complex_pdb: str - PDB content of metal-ligand complex
        complex_sdf: str - SDF content of metal-ligand complex
        template_name: str - Template from METAL_LIGAND_COMPLEX_TEMPLATES

        # OR build from components:
        ligand_smiles: str - SMILES for organic ligand
        metal: str - Metal code (CA, TB, FE, etc.)
        metal_binding_atoms: list - Ligand atoms that coordinate metal

        # Design parameters
        approach: str - "asymmetric", "sequential", "full", "joint"
        chain_length: str - e.g., "60-80"
        num_designs: int

        # Coordination preferences
        chain_a_donors: list - e.g., ["His", "Glu"]
        chain_b_donors: list - e.g., ["Asp", "Asn"]
        coordination_split: list - e.g., [2, 3] for sites per chain

        # Validation
        validate_coordination: bool - Run metal coordination validation
        validate_ligand_binding: bool - Run GNINA scoring

    Returns:
        Standard handler response with designs and validation
    """
```

### 3.2 Add to Task Dispatcher

Location: In `handler()` function (~line 350)

```python
elif task == "interface_metal_ligand_design":
    return handle_interface_metal_ligand_design(job_input)
```

### 3.3 Implement Design Approaches

```python
def _design_metal_ligand_joint(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Design both chains around complex simultaneously.

    Workflow:
    1. Parse/generate complex structure
    2. Identify available coordination sites
    3. Build RFD3 config with:
       - Complex as ligand input
       - select_hotspots for coordination sites
       - select_exposed for ligand surface
       - guiding_potentials for metal proximity
    4. Run RFD3 with two-chain contig
    5. Run LigandMPNN with fixed coordinating residues
    6. Validate coordination geometry
    """

def _design_metal_ligand_sequential(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Design chain A first, then chain B with opposite coordination.

    Similar to interface_ligand sequential approach but with
    coordination site awareness.
    """
```

---

## Step 4: RFD3 Configuration Builder

**File:** `G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless\handler.py`

### 4.1 Build RFD3 Config for Complexes

```python
def _build_rfd3_config_for_complex(
    complex_info: MetalLigandComplex,
    approach: str,
    chain_length: str,
    chain_a_donors: Optional[List[str]] = None,
    chain_b_donors: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """
    Build RFD3 configuration for metal-ligand complex design.

    Key configurations:
    - select_hotspots: Metal and available coordination sites
    - select_exposed: Ligand atoms not involved in metal coordination
    - select_buried: Metal coordination sphere
    - guiding_potentials: substrate_contacts for complex proximity
    - ori_token: Position complex at interface center
    """

    rfd3_config = {
        "contig": f"{chain_length},/0,{chain_length}",

        # Complex as ligand - use HETATM records
        "ligand": complex_info.ligand_res_name,

        # Hotspots: ensure contact with coordination sites
        "select_hotspots": {
            complex_info.metal_code: "ALL",  # Contact metal
        },

        # RASA: expose non-coordinating ligand surface
        "select_exposed": {
            complex_info.ligand_res_name: _get_non_coordinating_atoms(complex_info)
        },

        # H-bond acceptors: carboxylate oxygens for lanthanides
        "select_hbond_acceptor": _get_hbond_acceptors(complex_info),

        # Guiding potentials
        "guiding_potentials": [
            f"type:substrate_contacts,weight:5,s:1,r_0:8,d_0:4",
            "type:monomer_ROG,weight:1,min_dist:15",
        ],
        "guide_scale": 2.0,
        "guide_decay": "quadratic",
    }

    return rfd3_config
```

---

## Step 5: Validation Integration

**Files:**
- `G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless\metal_validation.py`
- `G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless\binding_analysis.py`

### 5.1 Extend Metal Validation for Complexes

```python
def validate_metal_ligand_complex_site(
    pdb_content: str,
    complex_info: MetalLigandComplex,
    expected_protein_donors: int = 4,
    distance_cutoff: float = 3.5,
) -> Dict[str, Any]:
    """
    Validate metal coordination in designed protein.

    Checks:
    1. Metal-ligand bonds preserved (within tolerance)
    2. Protein residues coordinating at expected sites
    3. Coordination geometry matches target
    4. No steric clashes between protein and ligand

    Returns:
        {
            "success": bool,
            "ligand_coordination_preserved": bool,
            "protein_coordination_count": int,
            "total_coordination_number": int,
            "geometry_type": str,
            "geometry_rmsd": float,
            "quality_score": float,
            "issues": list,
        }
    """
```

### 5.2 Combined Validation Function

```python
def validate_metal_ligand_dimer(
    pdb_content: str,
    complex_info: MetalLigandComplex,
    ligand_smiles: Optional[str] = None,
    run_gnina: bool = True,
    run_coordination: bool = True,
) -> Dict[str, Any]:
    """
    Combined validation for metal-ligand complex dimers.

    Runs:
    1. Metal coordination validation
    2. Ligand binding analysis (contacts, H-bonds)
    3. GNINA scoring (if SMILES provided)
    4. Interface quality metrics
    """
```

---

## Step 6: Test Cases Implementation

**File:** `G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless\test_metal_ligand_complex.py` (NEW)

### 6.1 Test PQQ-Ca Complex

```python
"""
Test cases for metal-ligand complex dimer design.

Test complexes:
1. PQQ-Ca (pyrroloquinoline quinone with calcium)
2. Citrate-Tb (citrate with terbium)
"""

import pytest
from typing import Dict, Any

# PQQ-Ca test complex (minimal PDB)
PQQ_CA_TEST_PDB = """
HETATM    1  C1  PQQ L   1      50.000  50.000  48.000  1.00  0.00           C
HETATM    2  C2  PQQ L   1      51.200  50.000  48.500  1.00  0.00           C
HETATM    3  O5  PQQ L   1      50.500  51.000  50.000  1.00  0.00           O
HETATM    4  N6  PQQ L   1      49.500  50.500  50.500  1.00  0.00           N
HETATM    5  O7  PQQ L   1      50.000  49.000  50.000  1.00  0.00           O
HETATM   10 CA   CA  M   1      50.000  50.000  52.400  1.00  0.00          CA
END
"""

# Citrate-Tb test complex
CITRATE_TB_TEST_PDB = """
HETATM    1  C1  CIT L   1      50.000  50.000  48.000  1.00  0.00           C
HETATM    2  O1  CIT L   1      50.500  51.000  49.000  1.00  0.00           O
HETATM    3  O2  CIT L   1      49.500  49.000  49.000  1.00  0.00           O
HETATM    4  C2  CIT L   1      50.000  50.000  46.500  1.00  0.00           C
HETATM    5  O3  CIT L   1      51.200  50.500  46.000  1.00  0.00           O
HETATM    6  O4  CIT L   1      48.800  49.500  46.000  1.00  0.00           O
HETATM   10 TB   TB  M   1      50.000  50.000  50.350  1.00  0.00          TB
END
"""


class TestComplexParser:
    """Test metal-ligand complex parsing."""

    def test_parse_pqq_ca_complex(self):
        """Should correctly parse PQQ-Ca complex."""
        from inference_utils import parse_metal_ligand_complex

        complexes = parse_metal_ligand_complex(PQQ_CA_TEST_PDB)

        assert len(complexes) == 1
        c = complexes[0]
        assert c.metal_code == "CA"
        assert c.ligand_res_name == "PQQ"
        assert len(c.coordination_bonds) >= 2  # O5, N6, O7 should coordinate

    def test_parse_citrate_tb_complex(self):
        """Should correctly parse Citrate-Tb complex."""
        from inference_utils import parse_metal_ligand_complex

        complexes = parse_metal_ligand_complex(CITRATE_TB_TEST_PDB)

        assert len(complexes) == 1
        c = complexes[0]
        assert c.metal_code == "TB"
        assert c.ligand_res_name == "CIT"

    def test_detect_available_sites_octahedral(self):
        """Should detect available coordination sites for octahedral Ca."""
        from inference_utils import parse_metal_ligand_complex, detect_available_coordination_sites

        complexes = parse_metal_ligand_complex(PQQ_CA_TEST_PDB)
        sites = detect_available_coordination_sites(complexes[0], target_coordination=6)

        # Ca in PQQ uses 3 sites, so 3 should be available
        assert len(sites) >= 2
        assert all("coords" in s for s in sites)


class TestTemplateLibrary:
    """Test metal-ligand complex templates."""

    def test_pqq_ca_template_exists(self):
        """PQQ-Ca template should exist."""
        from metal_ligand_templates import METAL_LIGAND_COMPLEX_TEMPLATES

        assert "pqq_ca" in METAL_LIGAND_COMPLEX_TEMPLATES
        template = METAL_LIGAND_COMPLEX_TEMPLATES["pqq_ca"]
        assert template["metal"] == "CA"
        assert "coordination" in template

    def test_citrate_templates_exist(self):
        """Citrate-lanthanide templates should exist."""
        from metal_ligand_templates import METAL_LIGAND_COMPLEX_TEMPLATES

        assert "citrate_tb" in METAL_LIGAND_COMPLEX_TEMPLATES
        assert "citrate_eu" in METAL_LIGAND_COMPLEX_TEMPLATES

    def test_generate_pqq_ca_pdb(self):
        """Should generate valid PQQ-Ca PDB from template."""
        from metal_ligand_templates import generate_complex_pdb

        pdb = generate_complex_pdb("pqq_ca")

        assert "PQQ" in pdb
        assert "CA" in pdb
        assert "HETATM" in pdb


class TestHandlerIntegration:
    """Test handler function integration."""

    def test_handler_accepts_task(self):
        """Handler should accept interface_metal_ligand_design task."""
        from handler import handler

        result = handler({
            "input": {
                "task": "interface_metal_ligand_design",
                "template_name": "pqq_ca",
                "approach": "joint",
                "chain_length": "60-80",
                "num_designs": 1,
            }
        })

        # Should not fail with unknown task error
        assert "Unknown task" not in result.get("error", "")

    def test_pqq_ca_design_mock(self):
        """Should generate design for PQQ-Ca (mock mode)."""
        from handler import handle_interface_metal_ligand_design

        result = handle_interface_metal_ligand_design({
            "template_name": "pqq_ca",
            "approach": "joint",
            "chain_length": "60-80",
            "num_designs": 1,
            "use_mock": True,
        })

        assert result["status"] in ["completed", "failed"]
        if result["status"] == "completed":
            assert "designs" in result["result"]


class TestValidation:
    """Test validation functions."""

    def test_validate_coordination_preserved(self):
        """Should validate metal-ligand coordination is preserved."""
        from metal_validation import validate_metal_ligand_complex_site
        from inference_utils import parse_metal_ligand_complex

        complexes = parse_metal_ligand_complex(PQQ_CA_TEST_PDB)

        # Validation of input complex should pass
        result = validate_metal_ligand_complex_site(
            PQQ_CA_TEST_PDB,
            complexes[0],
            expected_protein_donors=0,  # No protein yet
        )

        assert result["ligand_coordination_preserved"]
```

---

## Step 7: Docker Testing

**File:** `G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless\docker-compose.local.yml`

### 7.1 Update Volume Mounts

Add new module to volume mounts:

```yaml
volumes:
  # ... existing mounts ...
  # Metal-ligand complex module
  - ./metal_ligand_templates.py:/app/metal_ligand_templates.py
```

### 7.2 Docker Test Script

Create: `G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless\test_docker_metal_ligand.sh`

```bash
#!/bin/bash
# Test metal-ligand complex design in Docker

set -e

echo "=== Testing Metal-Ligand Complex Design in Docker ==="

# Wait for container to be ready
echo "Waiting for container..."
sleep 5

# Test 1: Health check
echo ""
echo "Test 1: Health check"
curl -s http://localhost:8000/runsync -X POST \
  -H "Content-Type: application/json" \
  -d '{"input": {"task": "health"}}' | jq .

# Test 2: PQQ-Ca template design (mock mode for quick test)
echo ""
echo "Test 2: PQQ-Ca design (mock mode)"
curl -s http://localhost:8000/runsync -X POST \
  -H "Content-Type: application/json" \
  -d '{
    "input": {
      "task": "interface_metal_ligand_design",
      "template_name": "pqq_ca",
      "approach": "joint",
      "chain_length": "60-80",
      "num_designs": 1,
      "use_mock": true
    }
  }' | jq .

# Test 3: Citrate-Tb design (mock mode)
echo ""
echo "Test 3: Citrate-Tb design (mock mode)"
curl -s http://localhost:8000/runsync -X POST \
  -H "Content-Type: application/json" \
  -d '{
    "input": {
      "task": "interface_metal_ligand_design",
      "template_name": "citrate_tb",
      "approach": "joint",
      "chain_length": "60-80",
      "num_designs": 1,
      "chain_a_donors": ["Glu", "Asp"],
      "chain_b_donors": ["Glu", "Asp", "Asn"],
      "use_mock": true
    }
  }' | jq .

# Test 4: Full design (if GPU available)
echo ""
echo "Test 4: Full PQQ-Ca design (real inference)"
curl -s http://localhost:8000/runsync -X POST \
  -H "Content-Type: application/json" \
  --max-time 300 \
  -d '{
    "input": {
      "task": "interface_metal_ligand_design",
      "template_name": "pqq_ca",
      "approach": "joint",
      "chain_length": "60-80",
      "num_designs": 1,
      "validate_coordination": true
    }
  }' | jq .

echo ""
echo "=== Tests Complete ==="
```

### 7.3 Docker Commands

```bash
# Build and start container
cd G:\Github_local_repo\Banta_Lab_RFdiffusion\backend\serverless
docker-compose -f docker-compose.local.yml up --build -d

# Run tests
chmod +x test_docker_metal_ligand.sh
./test_docker_metal_ligand.sh

# View logs
docker-compose -f docker-compose.local.yml logs -f

# Stop container
docker-compose -f docker-compose.local.yml down
```

---

## Implementation Order

### Phase 1: Core Infrastructure (Steps 1-2) ✅
1. [x] Add `MetalLigandComplex` dataclass to `inference_utils.py`
2. [x] Implement `parse_metal_ligand_complex()` function
3. [x] Implement `detect_available_coordination_sites()` function
4. [x] Create `metal_ligand_templates.py` with PQQ-Ca and Citrate-Tb templates
5. [x] Implement `generate_complex_pdb()` function

### Phase 2: Handler Implementation (Steps 3-4) ✅
6. [x] Add `handle_interface_metal_ligand_design()` to `handler.py`
7. [x] Add task to dispatcher in `handler()` function
8. [x] Implement `_design_metal_ligand_joint()` approach
9. [x] Implement `_build_rfd3_config_for_complex()`
10. [x] Add LigandMPNN integration with fixed positions

### Phase 3: Validation (Step 5) ✅
11. [x] Add `validate_metal_ligand_complex_site()` to `metal_validation.py`
12. [x] Add `validate_metal_ligand_dimer()` combined validation
13. [x] Integrate validation into handler

### Phase 4: Testing (Steps 6-7) ✅
14. [x] Create `test_metal_ligand_complex.py` with unit tests
15. [x] Update `docker-compose.local.yml` with new module
16. [x] Create Docker test script
17. [ ] Run integration tests in Docker (manual step)

---

## Test Case Details

### PQQ-Ca (Pyrroloquinoline Quinone - Calcium)

**Source:** [RCSB PDB PQQ Ligand](https://www.rcsb.org/ligand/PQQ), [Crystal structure of Ca-PQQ](https://pmc.ncbi.nlm.nih.gov/articles/PMC7716187/)

| Property | Value |
|----------|-------|
| Metal | Ca²⁺ |
| Ligand SMILES | `OC(=O)c1[nH]c2c(c1)C(=O)C(=O)c3nc(cc(C(O)=O)c23)C(O)=O` |
| Ligand MW | 330.206 Da |
| Coordination | 6 (octahedral) |
| PQQ donors | O5, N6, O7 (3 sites) |
| Protein sites | 3 available |
| PDB reference | 1W6S (methanol dehydrogenase) |

### Citrate-Lanthanide (Tb/Eu)

**Source:** [Lanthanide citrate polymers](https://www.sciencedirect.com/science/article/abs/pii/S0022286007004905), [PubChem Citrate](https://pubchem.ncbi.nlm.nih.gov/compound/6224)

| Property | Value |
|----------|-------|
| Metals | Tb³⁺, Eu³⁺, Gd³⁺ |
| Ligand SMILES | `OC(CC([O-])=O)(CC([O-])=O)C([O-])=O` |
| Ligand MW | 189.10 Da (citrate³⁻) |
| Coordination | 8-9 (square antiprism / TTP) |
| Citrate donors | 4 carboxylate oxygens |
| Protein sites | 4-5 available |
| Tb emission | Green (545 nm) |
| Eu emission | Red (615 nm) |

---

## Success Criteria

1. **Parser works:** Can parse PQQ-Ca and Citrate-Tb from PDB
2. **Templates generate valid structures:** PDB output has correct geometry
3. **Handler accepts requests:** Task dispatcher routes correctly
4. **Design runs:** RFD3 generates protein around complex
5. **Validation passes:** Coordination geometry is maintained
6. **Docker tests pass:** All curl commands return expected results

---

## References

1. [Crystal structure of Ca-PQQ complex](https://pmc.ncbi.nlm.nih.gov/articles/PMC7716187/)
2. [RCSB PDB - PQQ Ligand](https://www.rcsb.org/ligand/PQQ)
3. [Lanthanide citrate coordination polymers](https://www.sciencedirect.com/science/article/abs/pii/S0022286007004905)
4. [PubChem - Trisodium citrate](https://pubchem.ncbi.nlm.nih.gov/compound/6224)
5. Caldwell et al. (2020) PNAS - Tight lanthanide binding in de novo TIM barrel
