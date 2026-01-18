# Interface Metal-Ligand Dimer Design Workflow

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Design protein homodimers that bind metal-ligand complexes (e.g., citrate-terbium) at the interface, with proper coordination chemistry and clash-free structures.

**Architecture:** 7-stage pipeline: (1) Parse complex → (2) RFD3 backbone generation with RASA/hotspot conditioning → (3) Restore template HETATM → (4) LigandMPNN sequence design with HSAB bias → (5) Apply sequence + ligand-aware FastRelax → (6) Restore HETATM → (7) Validation pipeline.

**Tech Stack:** RFdiffusion3 (Foundry API), LigandMPNN, PyRosetta, RDKit, rdkit-to-params, biotite

---

## Pipeline Overview

```
INPUT: Metal-Ligand Complex (PDB or template name)
         │
         ▼
┌─────────────────────────────────────────────────────────┐
│ Stage 1: Parse & Prepare Complex                        │
│ • Parse metal/ligand from PDB                           │
│ • Get HSAB chemistry (bias, coordination range)         │
│ • Set coordination split per chain                      │
└─────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────┐
│ Stage 2: RFdiffusion3 Backbone Generation               │
│ • contig: "60-80,/0,60-80" (two new chains)            │
│ • select_buried: bury metal + coordinating atoms        │
│ • select_fixed_atoms: lock metal position               │
│ • hotspots: ensure contact with metal+ligand            │
│ • select_hbond_acceptor: H-bond to carboxylates         │
└─────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────┐
│ Stage 3: Restore Original Metal-Ligand Geometry         │
│ • Remove RFD3's distorted HETATM                        │
│ • Replace with template HETATM (correct geometry)       │
└─────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────┐
│ Stage 4: LigandMPNN Sequence Design                     │
│ • fix_pdb_for_mpnn() - fix NaN coords                   │
│ • HSAB bias: "E:3.0,D:3.0,N:2.4,Q:2.4,H:0.5,C:-5.0"    │
│ • pack_side_chains + pack_with_ligand_context           │
└─────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────┐
│ Stage 5: Apply Sequence & Ligand-Aware FastRelax        │
│ • apply_designed_sequence_to_backbone()                 │
│ • fastrelax_with_ligand_in_pose():                      │
│   - Load ligand from PDB (preserves coords)             │
│   - Params.from_mol() for correct atom names            │
│   - fa_rep naturally avoids ligand clashes              │
└─────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────┐
│ Stage 6: Restore HETATM (if lost)                       │
│ • Re-append template HETATM if missing                  │
│ • Renumber atom serials                                 │
└─────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────┐
│ Stage 7: Validation Pipeline                            │
│ • Coordination: donors, geometry, HSAB compliance       │
│ • Clashes: ligand-protein, metal-protein, internal      │
│ • Scoring: +10/coord, +20/chain, -15/clash              │
└─────────────────────────────────────────────────────────┘
         │
         ▼
OUTPUT: Validated designs with scores and metrics
```

---

## Task 1: Parse Metal-Ligand Complex

**Files:**
- Reference: `handler.py:6375-6572` (`handle_interface_metal_ligand_design`)
- Reference: `metal_ligand_templates.py` (template definitions)

**Input Options:**
1. `complex_pdb` - User-provided PDB with HETATM
2. `template_name` - e.g., "citrate_tb", "pqq_ca"
3. `ligand_smiles` + `metal` - Build from components (not yet implemented)

**Key Functions:**
```python
# Parse complex from PDB
complexes = parse_metal_ligand_complex(complex_pdb)
complex_info = complexes[0]  # MetalLigandComplex object
metal = complex_info.metal_code  # "TB"
ligand_name = complex_info.ligand_res_name  # "CIT"

# Or use template with database fallback
template = get_template_with_fallback(template_name, metal, ligand)
complex_pdb = generate_complex_pdb(template_name, metal, center=(0.0, 0.0, 0.0))
```

**HSAB Chemistry Integration:**
```python
from metal_chemistry import get_amino_acid_bias, get_hsab_class, get_coordination_number_range

hsab_bias = get_amino_acid_bias(metal, oxidation_state)  # "E:3.0,D:3.0,..."
hsab_class = get_hsab_class(metal, oxidation_state)  # "hard"
cn_min, cn_max = get_coordination_number_range(metal, oxidation_state)  # (8, 9) for Tb³⁺
```

---

## Task 2: RFdiffusion3 Backbone Generation

**Files:**
- Reference: `handler.py:6730-6870` (`_design_metal_ligand_joint`)
- Reference: `inference_utils.py` (`run_rfd3_inference`)

**Contig Construction:**
```python
# Two new chains around the ligand
contig = f"{min_len}-{max_len},/0,{min_len}-{max_len}"  # e.g., "60-80,/0,60-80"
```

**Critical Conditioning Parameters:**

| Parameter | Purpose | Example |
|-----------|---------|---------|
| `select_buried` | Force protein contact | `{M1: "ALL", L1: "O2,O5,O7"}` |
| `select_fixed_atoms` | Lock metal position | `{M1: "ALL"}` |
| `hotspots` | Ensure contact | `["M1", "L1"]` |
| `select_hbond_acceptor` | H-bond to carboxylates | `{L1: "O1,O3,O4,O6"}` |

**RFD3 Call:**
```python
result = run_rfd3_inference(
    contig=contig,
    pdb_content=complex_pdb,
    num_designs=1,
    seed=design_seed,
    ligand=ligand_name,  # "CIT"
    select_buried={"M1": "ALL", "L1": "O2,O5,O7"},
    select_fixed_atoms={"M1": "ALL"},
    hotspots=["M1", "L1"],
    select_hbond_acceptor={"L1": "O1,O3,O4,O6"},
)
```

---

## Task 3: Restore Original Metal-Ligand Geometry

**Files:**
- Reference: `handler.py:6892-6930`

**Problem:** RFD3's internal ligand processing can distort metal-ligand geometry.

**Solution:** Always replace RFD3's HETATM with original template:
```python
# Remove ALL HETATM from RFD3 output
atom_lines = [line for line in pdb_content.split('\n')
              if line.startswith('ATOM') or line.startswith('TER')]

# Get template HETATM (correct geometry)
hetatm_lines = [line for line in complex_pdb.split('\n')
                if line.startswith('HETATM')]

# Renumber and append
max_atom_num = max(int(line[6:11].strip()) for line in atom_lines if line.startswith('ATOM'))
renumbered_hetatm = []
for i, line in enumerate(hetatm_lines):
    new_atom_num = max_atom_num + i + 1
    renumbered_hetatm.append(f"{line[:6]}{new_atom_num:5d}{line[11:]}")

pdb_content = '\n'.join(atom_lines) + '\n' + '\n'.join(renumbered_hetatm) + '\nEND\n'
```

---

## Task 4: LigandMPNN Sequence Design

**Files:**
- Reference: `handler.py:6933-7136`
- Reference: `inference_utils.py:1295-1400` (`fix_pdb_for_mpnn`)

**PDB Fixes (Critical for NaN Prevention):**
```python
from inference_utils import fix_incomplete_backbone, fix_pdb_for_mpnn

pdb_for_mpnn = fix_incomplete_backbone(pdb_content)
pdb_for_mpnn = fix_pdb_for_mpnn(pdb_for_mpnn)  # Fixes:
# 1. Renumbers ALL atom serial numbers sequentially
# 2. Renames duplicate atom names (O -> O1, O2)
# 3. Skips atoms with NaN coordinates
# 4. Skips malformed lines
```

**LigandMPNN Configuration:**
```python
mpnn_input = {
    "pdb_content": pdb_for_mpnn,
    "num_sequences": 1,
    "temperature": 0.1,
    "model_type": "ligand_mpnn",
    "remove_waters": True,
    "bias_AA": "E:3.0,D:3.0,N:2.4,Q:2.4,H:0.5,C:-5.0",  # HSAB bias
    "omit_AA": "C",  # Exclude Cys for lanthanides
    "ligand_cutoff_for_score": 4.0,
    "model_noise_level": "010",
    "pack_side_chains": True,
    "pack_with_ligand_context": True,
    "number_of_packs_per_design": 4,
}
mpnn_result = handle_mpnn(mpnn_input)
```

---

## Task 5: Apply Sequence & Ligand-Aware FastRelax

**Files:**
- Reference: `handler.py:7017-7085`
- Reference: `rosetta_utils.py:521-800` (`fastrelax_with_ligand_in_pose`)
- Reference: `rosetta_utils.py:148-230` (`generate_ligand_params_from_pdb`)

**Step 5.1: Apply Designed Sequence**
```python
from inference_utils import apply_designed_sequence_to_backbone

# Extract sequence from LigandMPNN FASTA output
fasta_content = sequences[0].get("content", "")
lines = [l for l in fasta_content.strip().split('\n') if not l.startswith('>')]
designed_seq = ''.join(lines).upper()

mutated_pdb = apply_designed_sequence_to_backbone(pdb_content, designed_seq)
```

**Step 5.2: Ligand-Aware FastRelax (CRITICAL FIX)**

The key insight: PyRosetta's `fa_rep` only avoids clashes with atoms IN the pose. We must load the ligand into the pose.

```python
from rosetta_utils import fastrelax_with_ligand_in_pose

CITRATE_SMILES = "OC(=O)CC(O)(C(=O)O)CC(=O)O"

relax_result = fastrelax_with_ligand_in_pose(
    mutated_pdb,
    ligand_smiles=CITRATE_SMILES,
    ligand_residue_name="CIT",
    max_iter=200,
    repack_only=True,
)
```

**How `generate_ligand_params_from_pdb` Works (PROPER FIX):**
```python
from rdkit import Chem
from rdkit_to_params import Params

# OLD (WRONG): from_smiles_w_pdbfile generates NEW RDKit coordinates
# p = Params.from_smiles_w_pdbfile(pdb_path, smiles, False, name)

# NEW (CORRECT): Load mol directly from PDB - preserves atom names AND coordinates
mol = Chem.MolFromPDBFile(pdb_path, removeHs=False)
p = Params.from_mol(mol, name=name)
# Result: Ligand at original PDB coordinates, fa_rep can detect real clashes
```

---

## Task 6: Restore HETATM (if lost)

**Files:**
- Reference: `handler.py:7097-7122`

**Problem:** FastRelax may not output HETATM records.

**Solution:**
```python
hetatm_count = sum(1 for line in pdb_content.split('\n') if line.startswith('HETATM'))
if hetatm_count == 0 and complex_pdb:
    # Extract HETATM from template
    hetatm_lines = [line for line in complex_pdb.split('\n') if line.startswith('HETATM')]

    # Find max atom serial
    max_atom_num = max(int(line[6:11].strip()) for line in pdb_content.split('\n')
                       if line.startswith('ATOM'))

    # Renumber and append
    renumbered = [f"{line[:6]}{max_atom_num + i + 1:5d}{line[11:]}"
                  for i, line in enumerate(hetatm_lines)]

    pdb_content = pdb_content.rstrip() + '\n' + '\n'.join(renumbered) + '\nEND\n'
```

---

## Task 7: Validation Pipeline

**Files:**
- Reference: `handler.py:7142-7206`
- Reference: `validation_pipeline/` directory

**Validation Components:**

### 7.1 Coordination Validation
```python
validation_result = _validate_metal_ligand_design(
    pdb_content, complex_info, target_coord,
    template_coord_info=coord_info
)
# Returns: coordination_number, chain_a_donors, chain_b_donors, nearby_residues
```

### 7.2 Pipeline Validation
```python
from validation_pipeline import ValidationPipeline

pipeline = ValidationPipeline()
pipeline_report = pipeline.validate(
    pdb_content,
    metal=complex_info.metal_code,
    ligand=complex_info.ligand_res_name,
    target_coordination=target_coord,
    expected_ligand_donors=ligand_donors,
)
# Returns: clashes, coordination, geometry, summary
```

### 7.3 Scoring
```python
score = 0
if coordination_number >= 4:
    score += coordination_number * 10
if chain_a_donors >= 1:
    score += 20
if chain_b_donors >= 1:
    score += 20

# Subtract for clashes
if pipeline_report.clashes:
    score -= pipeline_report.clashes.total_clash_count * 15
```

---

## API Usage Example

**Request:**
```json
{
  "input": {
    "task": "interface_metal_ligand_design",
    "template_name": "citrate_tb",
    "contig_str": "60-80,/0,60-80",
    "num_designs": 3,
    "approach": "joint",
    "validate_coordination": true
  }
}
```

**Or with custom PDB:**
```json
{
  "input": {
    "task": "interface_metal_ligand_design",
    "complex_pdb": "HETATM    1  C1  CIT...",
    "contig_str": "60-80,/0,60-80",
    "num_designs": 3
  }
}
```

**Response:**
```json
{
  "status": "COMPLETED",
  "output": {
    "result": {
      "designs": [
        {
          "pdb_content": "ATOM...",
          "score": 45,
          "pipeline_validation": {
            "clashes": {"total_clash_count": 4, "worst_overlap": 0.40},
            "coordination": {"total": 5, "ligand_donors": 3, "protein_donors": 2}
          },
          "pipeline_passed": false
        }
      ],
      "best_design": 0
    }
  }
}
```

---

## Key Bug Fixes Documented

### Fix 1: LigandMPNN NaN Coordinates
**Problem:** RFD3 output has duplicate atom names/serials causing NaN in LigandMPNN.
**Solution:** `fix_pdb_for_mpnn()` renumbers atoms and renames duplicates.

### Fix 2: Ligand Coordinate Mismatch in FastRelax
**Problem:** `Params.from_smiles_w_pdbfile()` generates RDKit coordinates, not PDB coordinates.
**Solution:** Use `Chem.MolFromPDBFile()` + `Params.from_mol()` to preserve exact PDB coordinates.

**Before:** 24-38 clashes (ligand at wrong position)
**After:** 1-6 clashes (ligand at correct position, fa_rep works)

---

## Testing

**Run test:**
```bash
python3 -c "
import json, requests
with open('test_data/citrate_tb_complex.pdb') as f:
    pdb = f.read()
r = requests.post('http://localhost:8000/runsync', json={
    'input': {
        'task': 'interface_metal_ligand_design',
        'complex_pdb': pdb,
        'contig_str': '60-80,/0,60-80',
        'num_designs': 1
    }
})
print(r.json()['status'])
"
```

**Expected:** `COMPLETED` with designs showing low clash counts.

---

## Future Improvements

1. **Clash tolerance tuning:** 0.04Å overlap is within VDW tolerance, consider adjusting threshold
2. **Coordination geometry detection:** Currently shows "unknown", needs improvement
3. **Sequential approach:** Design chains one at a time for better control
4. **Multi-metal support:** Extend to bimetallic complexes
