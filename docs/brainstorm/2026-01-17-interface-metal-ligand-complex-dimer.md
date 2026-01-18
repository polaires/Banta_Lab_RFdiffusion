# Brainstorming: Interface Metal-Ligand Complex Dimer Design

## User Request
Design protein dimers with **metal-ligand complexes** at the interface, combining the concepts from:
1. **Interface Ligand Dimer** - organic molecules (azobenzene) at dimer interface
2. **Interface Metal Dimer** - metal ions (Zn, Tb, Fe) coordinated by protein residues

The goal is to support coordination compounds like **heme-iron**, **zinc-porphyrin**, **metal-cofactor complexes** where both the organic ligand AND the metal are critical for function.

---

## Problem Analysis

### What is a Metal-Ligand Complex?
A metal-ligand complex consists of:
1. **Central metal ion** (Fe, Zn, Cu, Mn, Co, etc.)
2. **Organic ligand(s)** that coordinate to the metal (porphyrin, bipyridine, salen, etc.)
3. **Additional coordination sites** that can be occupied by protein residues

**Examples:**
- **Heme (Fe-porphyrin)**: Iron coordinated by 4 porphyrin nitrogens + 1-2 axial ligands from protein (His, Cys)
- **Zinc-porphyrin**: Zinc in porphyrin ring + 1 axial ligand
- **Metal-bipyridine**: Metal + 2-3 bidentate bipyridine ligands + protein coordination
- **Chlorophyll**: Mg-porphyrin derivative
- **Cobalamin (B12)**: Cobalt in corrin ring + axial ligands

### Current System Gaps

| Feature | Interface Ligand | Interface Metal | Metal-Ligand Complex |
|---------|-----------------|-----------------|---------------------|
| Input format | SMILES string | Metal code (ZN) | **BOTH needed** |
| 3D structure | RDKit conformer | Single atom + template | **Pre-computed complex** |
| Binding mode | H-bonds, pi-stacking | Coordination bonds | **Mixed: coord + secondary** |
| Fixed positions | None | Template residues | **Metal coord sites** |
| RASA control | Yes (select_exposed) | Limited | **Yes + coord geometry** |

---

## Proposed Solutions

### Option 1: Pre-computed Complex PDB Input (Recommended)

**Approach:** User provides the metal-ligand complex as a PDB/SDF file with accurate 3D coordinates.

**Workflow:**
1. User provides complex structure from:
   - Cambridge Structural Database (CSD)
   - Protein Data Bank (PDB) - extract from existing proteins
   - Computational chemistry (DFT-optimized, Architector)
   - RDKit + metal placement

2. System processes the complex:
   - Identify metal center(s) and their coordination geometry
   - Identify organic ligand atoms
   - Determine available coordination sites (for protein binding)
   - Calculate hotspot positions for RFD3

3. Design proteins using existing approaches:
   - Use `interface_ligand` workflow with complex as "ligand"
   - Add `select_hotspots` for metal coordination sites
   - Apply `select_hbond_acceptor/donor` for ligand interactions

**Pros:**
- Accurate complex geometry from validated sources
- Flexible - works with ANY metal-ligand complex
- Minimal changes to existing codebase
- User has full control over complex structure

**Cons:**
- Requires user to source/prepare complex structure
- Learning curve for structure preparation

### Option 2: Template Library Expansion

**Approach:** Build a library of common metal-ligand complex templates (like `lanthanide_templates.py`).

**Workflow:**
1. Pre-define common complexes:
   ```python
   METAL_LIGAND_TEMPLATES = {
       "heme_b": {
           "metal": "FE",
           "ligand_smiles": "...", # Protoporphyrin IX
           "coordination": 6,
           "protein_sites": 2,  # Axial His positions
           "axial_geometry": "octahedral",
       },
       "zinc_porphyrin": {...},
       "chlorophyll_a": {...},
       "cobalamin": {...},
   }
   ```

2. Generate 3D coordinates from template + metal placement

3. Use existing workflows with template-derived structure

**Pros:**
- Easy for common cases (heme, chlorophyll)
- Consistent, validated geometries
- Lower user barrier

**Cons:**
- Limited to pre-defined complexes
- Maintenance burden for template library
- Complex biochemistry to capture correctly

### Option 3: SMILES + Metal Annotation

**Approach:** Extend SMILES input to annotate metal binding sites.

**Workflow:**
1. User provides:
   - Ligand SMILES (e.g., porphyrin without metal)
   - Metal type (FE, ZN, etc.)
   - Binding atom indices or SMARTS pattern

2. System:
   - Generates ligand conformer
   - Places metal at binding site with correct geometry
   - Adds available coordination sites for protein

**Pros:**
- Flexible - any ligand + metal combination
- Familiar SMILES workflow

**Cons:**
- Accurate metal placement is complex
- Requires sophisticated chemistry (bond lengths, angles)
- May not capture subtle geometry requirements

### Option 4: Architector Integration (Advanced)

**Approach:** Use [Architector](https://github.com/lanl/Architector) for metal complex generation.

Architector is a tool specifically designed for generating 3D structures of metal complexes with accurate coordination geometry.

**Workflow:**
1. User specifies:
   - Metal type and oxidation state
   - Ligand SMILES or name
   - Coordination number/geometry preference

2. Architector generates optimized 3D structure

3. System uses generated complex in design workflow

**Pros:**
- Chemically accurate metal placement
- Handles complex geometries (octahedral, tetrahedral, etc.)
- Academic tool with active development

**Cons:**
- External dependency
- May require installation/licensing
- Integration complexity

---

## Recommended Implementation Strategy

### Phase 1: Pre-computed Complex Input (MVP)

Add support for user-provided metal-ligand complex structures:

```python
def handle_interface_metal_ligand_design(job_input: Dict[str, Any]) -> Dict[str, Any]:
    """
    Design protein dimers around metal-ligand complexes.

    Input:
        complex_pdb: str (required) - PDB content of metal-ligand complex
        OR
        complex_sdf: str (required) - SDF content of metal-ligand complex

        # Auto-detection or manual specification
        metal_code: str (optional) - Override auto-detected metal (FE, ZN, etc.)
        coordination_atoms: list (optional) - Atoms available for protein coordination

        # Design parameters (inherited from interface_ligand_design)
        approach: str - "asymmetric", "sequential", "full"
        chain_length: str
        num_designs: int

        # From interface_metal_design
        chain_a_donors: list - Preferred coordinating residues for chain A
        chain_b_donors: list - Preferred coordinating residues for chain B
    """
```

**Key Implementation Tasks:**

1. **Complex Structure Parser**
   ```python
   def parse_metal_ligand_complex(pdb_content: str) -> Dict[str, Any]:
       """
       Parse metal-ligand complex and extract:
       - Metal center(s): element, coordinates, charge
       - Ligand atoms: coordinates, atom types
       - Coordination geometry: existing bonds, available sites
       - Hotspot positions: where protein can coordinate
       """
   ```

2. **Available Site Detection**
   ```python
   def detect_coordination_sites(complex_info: Dict) -> List[Dict]:
       """
       Identify positions where protein residues can coordinate.

       For heme: 2 axial positions above/below porphyrin plane
       For zinc-porphyrin: 1 axial position
       For incomplete complexes: unfilled coordination sites
       """
   ```

3. **RFD3 Configuration Builder**
   ```python
   def build_rfd3_config_for_complex(
       complex_info: Dict,
       approach: str,
       chain_length: str,
   ) -> Dict:
       """
       Build RFD3 configuration that:
       - Places complex at interface (ori_token)
       - Sets hotspots for coordination sites
       - Configures RASA for ligand binding surfaces
       - Sets up guiding potentials for metal coordination
       """
   ```

### Phase 2: Common Template Library

Add templates for frequently used complexes:

```python
METAL_LIGAND_COMPLEX_TEMPLATES = {
    # Heme variants
    "heme_b": {
        "description": "Heme B (protoporphyrin IX + Fe)",
        "pdb_template": "...",  # Pre-computed coordinates
        "metal": "FE",
        "oxidation_states": [2, 3],  # Fe2+ or Fe3+
        "protein_coordination": {
            "num_sites": 2,
            "geometry": "axial",
            "preferred_donors": ["His", "Cys", "Met"],
            "bond_distance": 2.0,  # Angstroms
        },
        "best_for": ["oxygen_binding", "electron_transfer", "catalysis"],
    },

    "heme_c": {
        "description": "Heme C (covalently attached)",
        "pdb_template": "...",
        "covalent_attachment": ["Cys-vinyl"],
        ...
    },

    # Zinc complexes
    "zinc_porphyrin": {
        "description": "Zinc tetraphenylporphyrin",
        "metal": "ZN",
        "protein_coordination": {
            "num_sites": 1,
            "geometry": "axial",
            "preferred_donors": ["His", "Glu", "Asp"],
        },
    },

    # Other cofactors
    "chlorophyll_a": {...},
    "cobalamin_b12": {...},
    "fad_flavin": {...},
    "nad_nicotinamide": {...},
}
```

### Phase 3: Structure Sources Integration

Integrate external structure databases:

1. **PDB Ligand Extraction**
   - Fetch ligand + coordinated metal from PDB entries
   - Example: Extract heme from 1MBO (myoglobin)

2. **CSD Integration** (if licensed)
   - Query Cambridge Structural Database for complexes
   - High-quality experimental geometries

3. **Architector Integration**
   - Generate complexes computationally
   - Support custom ligands

---

## Technical Considerations

### Accurate Structure Sources

| Source | Quality | Availability | Best For |
|--------|---------|--------------|----------|
| **PDB** | Good (but may need optimization) | Free | Heme, common cofactors |
| **CSD** | Excellent | Licensed | Novel complexes |
| **Architector** | Good | Free | Custom complexes |
| **DFT (Gaussian, ORCA)** | Excellent | Licensed | Precise geometries |
| **RDKit + Manual** | Fair | Free | Simple cases |

### Recommended Structure Preparation Workflow

1. **For common cofactors (heme, chlorophyll):**
   - Extract from high-resolution PDB structure
   - Optional: DFT optimization of extracted geometry

2. **For novel complexes:**
   - Use Architector to generate initial geometry
   - Validate coordination geometry
   - Optional: QM/MM refinement

3. **Quality checks:**
   - Bond lengths within expected ranges
   - Coordination geometry matches expected (octahedral, tetrahedral, etc.)
   - No steric clashes between ligand and expected protein positions

### Coordination Geometry Validation

```python
COORDINATION_GEOMETRIES = {
    "octahedral": {
        "coordination_number": 6,
        "angles": {90, 180},
        "example_metals": ["Fe", "Co", "Cr", "Mn"],
    },
    "tetrahedral": {
        "coordination_number": 4,
        "angles": {109.5},
        "example_metals": ["Zn", "Cu(I)"],
    },
    "square_planar": {
        "coordination_number": 4,
        "angles": {90, 180},
        "example_metals": ["Pt", "Pd", "Ni", "Cu(II)"],
    },
    "trigonal_bipyramidal": {
        "coordination_number": 5,
        "angles": {90, 120, 180},
        "example_metals": ["Fe", "Cu"],
    },
}
```

---

## RFD3 Parameter Mapping

### From Interface Ligand Design
- `ligand_smiles` → N/A (use complex_pdb instead)
- `select_exposed` → Yes, for ligand surface accessibility
- `select_buried` → Yes, for metal-coordination pocket
- `ori_token` → Position complex at interface

### From Interface Metal Design
- `template_name` → Extended to complex templates
- `select_fixed_atoms` → Fix metal coordination geometry
- `select_hotspots` → Coordination sites for protein
- `guiding_potentials` → `substrate_contacts` for complex binding
- `chain_a_donors` / `chain_b_donors` → Coordinating residue preferences

### New Parameters Needed
- `complex_pdb` / `complex_sdf` → Input structure
- `coordination_sites` → Available positions for protein
- `axial_preference` → Which axial site to fill first
- `coordination_geometry_constraint` → Enforce geometry

---

## Example Use Cases

### 1. Heme-binding Heterodimer

**Goal:** Design two proteins that sandwich a heme, each providing an axial ligand.

```json
{
  "task": "interface_metal_ligand_design",
  "complex_pdb": "<heme PDB content>",
  "approach": "joint",
  "chain_length": "60-80",
  "chain_a_donors": ["His"],  // Proximal His
  "chain_b_donors": ["His"],  // Distal His (or open for O2)
  "coordination_split": [1, 1],  // One axial site per chain
  "num_designs": 5
}
```

### 2. Light-harvesting Zinc Porphyrin Dimer

**Goal:** Design proteins around zinc porphyrin for energy transfer studies.

```json
{
  "task": "interface_metal_ligand_design",
  "template_name": "zinc_porphyrin",
  "approach": "asymmetric",
  "chain_a_donors": ["His"],  // Single axial ligand
  "add_antenna_residues": true,  // Trp for energy transfer
  "num_designs": 5
}
```

### 3. Enzyme Mimic with Metal-Salen Complex

**Goal:** Design artificial metalloenzyme using metal-salen catalyst.

```json
{
  "task": "interface_metal_ligand_design",
  "complex_sdf": "<Mn-salen SDF content>",
  "metal_code": "MN",
  "approach": "full",
  "chain_length": "80-100",
  "include_substrate_channel": true,
  "num_designs": 10
}
```

---

## Questions to Resolve

1. **Structure Source Priority:**
   - Should we prioritize pre-computed templates or user-provided structures?
   - What's the minimum viable set of templates for Phase 2?

2. **Coordination Validation:**
   - How strictly should we enforce coordination geometry?
   - Should we auto-fix sub-optimal geometries?

3. **LigandMPNN Handling:**
   - How does LigandMPNN handle metal-ligand complexes?
   - Do we need custom training data for metalloporphyrins?

4. **Axial Site Assignment:**
   - For heme with 2 axial sites, how do we assign chain A vs B?
   - Should this be automatic or user-specified?

5. **Covalent Attachments:**
   - Some cofactors (heme c, FAD) are covalently attached
   - Do we need to support `covalent_bonds` parameter?

---

## Recommended Next Steps

1. **Immediate (Phase 1):**
   - [ ] Implement `parse_metal_ligand_complex()` function
   - [ ] Add `handle_interface_metal_ligand_design()` handler
   - [ ] Create test cases with heme from PDB 1MBO

2. **Short-term (Phase 2):**
   - [ ] Build template library for heme_b, zinc_porphyrin
   - [ ] Add coordination geometry validation
   - [ ] Integration test with full design workflow

3. **Medium-term (Phase 3):**
   - [ ] Architector integration for custom complexes
   - [ ] PDB ligand extraction utility
   - [ ] Frontend UI for complex design

---

## References

1. Caldwell et al. (2020) PNAS - Tight lanthanide binding in de novo TIM barrel
2. Cambridge Structural Database - https://www.ccdc.cam.ac.uk/
3. Architector - https://github.com/lanl/Architector
4. RosettaCommons RFD3 documentation
5. [Transition metal porphyrin complexes](https://en.wikipedia.org/wiki/Transition_metal_porphyrin_complexes)

---

*Brainstorming document created: 2026-01-17*
*Status: Ready for discussion and refinement*
