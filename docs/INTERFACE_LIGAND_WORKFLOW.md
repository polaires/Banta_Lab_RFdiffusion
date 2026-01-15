# Interface Ligand Design Workflow

## Problem Statement
Design a C2 symmetric protein dimer with azobenzene ligand at the interface between the two chains.

**Target Metrics:**
- Both chains contact the ligand (≥3 contacts each)
- No interchain clashes
- GNINA affinity: -5 to -8 kcal/mol (useful binding)

---

## Approaches Tried (Chronological)

### Approach 1: Native RFD3 Symmetry with Seed Atom
**Hypothesis:** Add a glycine seed residue to provide protein entities for RFD3's native symmetry.

**Implementation:**
- Generate glycine seed at (+12, 0, +3) offset
- Use contig "A1,50-70" to extend from seed
- Enable native C2 symmetry

**Result:** FAILED
- Error: "Input has no possible symmetry. Multiplicity: 1"
- RFD3 requires multiplicity matching symmetry order (2 for C2)
- Single seed + ligand = multiplicity 1

**Lesson:** Native RFD3 symmetry requires either pre-existing symmetric input OR pure de novo (no input structure).

---

### Approach 2: RASA Conditioning with Half-Binding
**Hypothesis:** Use `select_exposed` on one side of ligand to prevent 360° wrapping.

**Implementation:**
- `select_exposed: {UNL: "C9,C10,C11,C12,C13,C14"}` (second phenyl ring)
- `select_partially_buried: {UNL: "N7,N8"}` (azo bridge)
- Post-processing C2 symmetry

**Result:** PARTIAL SUCCESS (geometry) / FAILED (binding)
- Geometry: Valid interface, no clashes
- Contacts: Asymmetric (8 vs 2)
- Affinity: -1.44 kcal/mol (very weak)

**Lesson:** RASA `select_exposed` prevents BOTH chains from binding those atoms after symmetry, creating asymmetric contacts.

---

### Approach 3: Natural Design + Post-Processing Symmetry
**Hypothesis:** Remove RASA conditioning; let natural design + post-processing symmetry create interface.

**Implementation:**
- Ligand oriented along Z-axis (C2 axis)
- No RASA conditioning
- Post-processing C2 symmetry with translation optimization
- Translation optimization: maximize (contacts_A + contacts_B + symmetry_bonus) - clash_penalty

**Result:** PARTIAL SUCCESS
- Geometry: Valid interface (≥3 contacts each chain)
- Affinity: -0.69 to -1.25 kcal/mol (very weak)

**Lesson:** Post-processing symmetry creates geometric validity but not proper binding pockets.

---

### Approach 4: Multiple Designs with Filtering
**Hypothesis:** Generate multiple designs, filter by GNINA affinity.

**Implementation:**
- Generate 5-10 designs per batch
- Evaluate each with GNINA
- Rank by composite score: affinity + contacts + symmetry

**Results by Chain Length:**

| Chain Length | Best Affinity | Valid Designs | Notes |
|--------------|---------------|---------------|-------|
| 50-70 | -0.94 kcal/mol | 1/10 | Most designs have positive (clashing) affinities |
| 80-100 | -3.78 kcal/mol | 7/10 | Better pocket flexibility |

**Key Finding:** High contact counts correlate with WORSE affinities (steric clashes).

| Contacts | Typical Affinity | Interpretation |
|----------|------------------|----------------|
| 0 | -1 to -2 | No binding, no clashes |
| 8-12 | -2 to -4 | Sweet spot: some binding, minimal clashes |
| 16-20 | +5 to +30 | Severe steric clashes |

**Lesson:** Longer chains (80-100) provide flexibility for better pocket formation. Moderate contacts are better than maximum contacts.

---

## Decision Tree for Future Sessions

```
User wants: Symmetric dimer with ligand at interface
│
├─ Can use native RFD3 symmetry?
│   └─ NO: Ligand-only input has no protein entities
│
├─ Try RASA conditioning?
│   └─ CAUTION: select_exposed affects BOTH chains after symmetry
│
├─ Use post-processing symmetry?
│   └─ YES, but with caveats:
│       ├─ Monomer wraps 360° around ligand
│       ├─ Symmetry creates overlapping half-pockets
│       └─ Result: geometric validity, weak binding
│
└─ Optimization strategy:
    ├─ Use longer chains (80-100+ residues)
    ├─ Generate multiple designs (10+)
    ├─ Filter by GNINA affinity, NOT contact count
    └─ Target moderate contacts (8-12), not maximum
```

---

## Critical Insights

### 1. Design Problem vs Geometry Problem
The user correctly identified: **This is a DESIGN problem, not a GEOMETRY problem.**

Post-processing symmetry can create geometric validity (ligand at interface, both chains contact it) but cannot fix fundamental design issues (no proper binding pocket).

### 2. Contact Count Paradox
More contacts ≠ better binding. High contact counts often indicate:
- Protein too close to ligand
- Steric clashes (positive GNINA affinity)
- Poor pocket shape

### 3. GNINA Affinity Thresholds
| Affinity | Interpretation |
|----------|----------------|
| < -8 | Strong (drug-like) |
| < -5 | Moderate (useful) |
| < -2 | Weak (marginal) |
| -2 to 0 | Very weak (thermal noise) |
| > 0 | Repulsive (clashes) |

### 4. Chain Length Matters
Longer chains (80-100+ residues) provide:
- More flexibility for pocket formation
- Better accommodation of ligand without clashing
- Higher success rate for valid designs

---

## Files Modified

1. **`backend/serverless/inference_utils.py`**
   - `_apply_post_symmetry()`: Translation optimization with clash detection
   - Interface ligand handling: ligand orientation, post-processing symmetry
   - Removed RASA conditioning for de novo designs (ineffective)

2. **`backend/serverless/handler.py`**
   - Added `to_python_types()` conversion for JSON serialization

3. **`backend/serverless/binding_analysis.py`**
   - Enhanced numpy type conversion (float32, bool_)

4. **`test_multi_design.js`**
   - Multi-design generation and GNINA filtering

---

## Best Results Achieved

| Metric | Value |
|--------|-------|
| Best Affinity | -3.78 kcal/mol |
| Contacts | 10 (6 + 4) |
| Symmetry | 80% |
| Chain Length | 80-100 residues |
| Valid Interface | YES |

**Gap to Target:** Need -5 to -8 kcal/mol; achieved -3.78 kcal/mol (1.2-4.2 kcal/mol gap)

---

## Approach 5: Binder-to-Target Mode (Priority 3)

**Hypothesis:** Use H-bond and burial conditioning to create explicit binding interactions.

**Implementation:**
- `select_hbond_acceptor: {UNL: "N7,N8"}` (azo nitrogens accept H-bonds)
- `select_buried: {UNL: "N7,N8"}` (bury azo bridge)
- `select_partially_buried: {UNL: "C4,C9"}` (partial burial of adjacent carbons)

**Result:** PARTIAL SUCCESS
- Best affinity: -3.30 kcal/mol (slightly worse than natural mode's -3.78)
- Best design: 100% symmetry (7 + 7 contacts)
- One design achieved an actual H-bond

**Comparison:**
| Mode | Best Affinity | Best Symmetry | H-bonds |
|------|---------------|---------------|---------|
| Natural | -3.78 kcal/mol | 80% | 0 |
| Binder | -3.30 kcal/mol | 100% | 1 |

**Lesson:** H-bond conditioning for de novo small molecule design has limited effect. The conditioning may work better for binder design against existing protein structures.

---

## Final Summary: All Approaches

| Approach | Best Affinity | Valid Designs | Key Issue |
|----------|---------------|---------------|-----------|
| Native symmetry + seed | FAILED | N/A | Multiplicity mismatch |
| RASA half-binding | -1.44 kcal/mol | YES | Asymmetric contacts |
| Natural + post-symmetry (50-70) | -0.94 kcal/mol | 1/10 | Short chains, clashes |
| Natural + post-symmetry (80-100) | **-3.78 kcal/mol** | 7/10 | Best approach |
| Binder mode (80-100) | -3.30 kcal/mol | 4/10 | Similar to natural |

**Recommended Approach:** Natural mode + 80-100 residue chains + 10 designs + filter by GNINA affinity

---

## Recommendations for AI Workflow

1. **Start with documentation review** - Check existing troubleshooting docs before attempting solutions
2. **Test incrementally** - One change at a time, verify each step
3. **Use quantitative metrics** - GNINA affinity, not just contact counts
4. **Generate multiple designs** - Statistical filtering is essential
5. **Record failed approaches** - Prevents re-attempting failed strategies
6. **Identify root cause** - "Design problem vs geometry problem" insight was key

---

## Running the Demo

A unified demo script is available for testing the interface ligand workflow.

### Quick Start

```bash
# Navigate to project root
cd G:\Github_local_repo\Banta_Lab_RFdiffusion

# Ensure the RFD3 server is running on localhost:8000

# Run demo with natural mode (best affinity results)
node scripts/demos/demo_azobenzene_interface.js

# Run demo with binder mode (H-bond conditioning)
node scripts/demos/demo_azobenzene_interface.js binder

# Customize number of designs
node scripts/demos/demo_azobenzene_interface.js natural 5
node scripts/demos/demo_azobenzene_interface.js binder 10
```

### Expected Output

```
============================================================
  DEMO: Azobenzene Dimerization Binder
============================================================

Mode: NATURAL
Designs: 10
Seed: <random>
Chain length: 80-100 residues
Ligand: Azobenzene (c1ccc(cc1)N=Nc2ccccc2)

Step 1: Generating designs...
Generated 10 designs in ~50s

Step 2: Analyzing geometry and binding...
  Design 1/10... -2.45 kcal/mol, 8 contacts, VALID
  ...

============================================================
  RESULTS SUMMARY
============================================================

Rank | Design | Affinity   | Contacts | Symmetry | Valid
-----|--------|------------|----------|----------|------
   1 |    3   |  -3.78    |    10    |    80%   |  YES
   ...

============================================================
  BEST DESIGN
============================================================
  Design #3
  GNINA Affinity: -3.78 kcal/mol
  Binding Quality: WEAK (marginal)
  ...
```

### Mode Comparison

| Mode | Best Affinity | Symmetry | H-bonds | Use Case |
|------|---------------|----------|---------|----------|
| natural | -3.78 kcal/mol | 80% | 0 | General interface design |
| binder | -3.30 kcal/mol | 100% | 1 | When H-bond interactions are important |

### Output Files

- `azobenzene_demo_natural_best.pdb` - Best design from natural mode
- `azobenzene_demo_binder_best.pdb` - Best design from binder mode

---

## Approach 6: Separable Dimer via Two-Step Asymmetric Design (RECOMMENDED)

**Date:** January 2026

### Problem with Post-Processing Symmetry

The previous approaches used post-processing symmetry which:
1. Creates a single monomer that wraps 360° around the ligand
2. Applies C2 symmetry to create a dimer
3. **Problem:** The chains are NOT separable - they're entangled around the ligand

### Solution: Independent Asymmetric Design + Combine

**Key Insight:** Contig-based binder design (`A1-76/0 80`) doesn't work for ligands - RFD3's binder mode is designed for protein-protein interactions, not adding chains around ligands with existing protein context.

**Working Approach:** Design both chains independently around the same ligand position, then combine:

1. **Step 1 (Asymmetric):** Design Chain A that binds only ONE side of the ligand
   - Use `select_exposed` to keep opposite side atoms accessible
   - Chain A leaves space for Chain B

2. **Step 2 (Independent):** Extract ligand from Chain A, design Chain B around ONLY the ligand
   - Extract the ligand PDB from Chain A's output (preserves exact 3D position)
   - Design Chain B using `length` (NOT contig) with the extracted ligand as input
   - Use OPPOSITE `select_exposed` atoms so Chain B binds the other side

3. **Step 3 (Combine):** Merge Chain A + Chain B + Ligand into final dimer
   - Extract protein-only from each design
   - Relabel chains (A, B, L)
   - Combine into single PDB

### Key Technical Details

#### Azobenzene Atom Naming
```
C1-C6:  First phenyl ring (left side)
N7-N8:  Azo bridge (-N=N-)
C9-C14: Second phenyl ring (right side)
```

#### RASA Conditioning for Asymmetric Binding
```python
# Step 1: Chain A binds C9-C14 (right), keeps C1-C6 exposed
rfd3_input["select_exposed"] = {"UNL": "C1,C2,C3,C4,C5,C6"}

# Step 2: Chain B binds C1-C6 (left), keeps C9-C14 exposed
# Uses extracted ligand PDB as input, NOT Chain A PDB
rfd3_input["select_exposed"] = {"UNL": "C9,C10,C11,C12,C13,C14"}
```

#### Why Contig-Based Sequential Design Failed

The original approach tried:
```python
contig = f"A1-{chain_a_length}/0 {chain_b_length}"  # e.g., "A1-76/0 80"
pdb_content = chain_a_pdb  # Chain A + ligand
```

**Problem:** RFD3 only output Chain A (76 residues) without generating Chain B. The contig-based binder design mode works for protein-protein interactions, not for adding chains around ligands.

#### Working Implementation

```python
# Step 2: Design Chain B around ONLY the ligand
ligand_pdb = _extract_ligand_pdb(chain_a_pdb, "UNL")  # Just the ligand

rfd3_input = {
    "task": "rfd3",
    "length": chain_length,      # Use length, NOT contig
    "pdb_content": ligand_pdb,   # Just the ligand, NOT Chain A
    "ligand": "UNL",
    "select_exposed": {"UNL": "C9,C10,C11,C12,C13,C14"},
}

# After RFD3 generates Chain B:
chain_b_protein = _extract_protein_pdb(chain_b_pdb)
chain_b_protein = _relabel_chain(chain_b_protein, "B")
dimer_pdb = _combine_chains(chain_a_protein, chain_b_protein, ligand_pdb)
```

#### Helper Functions Added

- `_extract_ligand_pdb()` - Extract HETATM records for ligand
- `_extract_protein_pdb()` - Extract ATOM records for protein
- `_relabel_chain()` - Change chain ID in PDB
- `_combine_chains()` - Merge Chain A + Chain B + Ligand

### API Usage

#### Full Dimer Design (Recommended)
```json
{
  "task": "interface_ligand_design",
  "approach": "full",
  "ligand_smiles": "c1ccc(cc1)N=Nc2ccccc2",
  "chain_length": "60-80",
  "num_designs": 5,
  "side": "left",
  "seed": 42
}
```

#### Response Structure
```json
{
  "status": "completed",
  "result": {
    "approach": "full",
    "chain_a": {
      "pdb_content": "...",
      "design_index": 0,
      "metrics": {...}
    },
    "dimer": {
      "pdb_content": "...",
      "design_index": 2,
      "affinity": -4.2,
      "contacts_a": 8,
      "contacts_b": 6,
      "has_clashes": false,
      "separable": true
    }
  }
}
```

### Workflow Diagram

```
Step 1: Asymmetric Chain A
┌─────────────────────────────────────────────┐
│                                             │
│   ┌─────┐                                   │
│   │Chain│    ┌───────────┐                  │
│   │  A  │◄───│ Azobenzene│   (C1-C6 exposed)│
│   └─────┘    └───────────┘                  │
│                                             │
└─────────────────────────────────────────────┘
         ▼
Step 2: Extract ligand, design Chain B independently
┌─────────────────────────────────────────────┐
│                                             │
│               ┌───────────┐   ┌─────┐       │
│               │ Azobenzene│───►│Chain│      │
│               └───────────┘    │  B  │      │
│               (from Step 1)    └─────┘      │
│               (C9-C14 exposed)              │
└─────────────────────────────────────────────┘
         ▼
Step 3: Combine Chain A + Chain B + Ligand
┌─────────────────────────────────────────────┐
│                                             │
│   ┌─────┐    ┌───────────┐   ┌─────┐        │
│   │Chain│◄───│ Azobenzene│───►│Chain│       │
│   │  A  │    └───────────┘    │  B  │       │
│   └─────┘                     └─────┘       │
│                                             │
└─────────────────────────────────────────────┘
         ▼
Result: Separable Dimer
- Both chains contact ligand
- Chains can physically separate
- No entanglement around ligand
```

### Backend Implementation

**File:** `backend/serverless/handler.py`

```python
# Key functions:
# - handle_interface_ligand_design(): Main entry point
# - _design_asymmetric_binder(): Step 1 - one-sided binder with select_exposed
# - _design_sequential_binder(): Step 2 - Extract ligand, design Chain B, combine
# - _design_full_dimer(): Orchestrates both steps

# Helper functions for combining chains:
# - _extract_ligand_pdb(): Extract ligand HETATM records
# - _extract_protein_pdb(): Extract protein ATOM records
# - _relabel_chain(): Change chain ID in PDB
# - _combine_chains(): Merge Chain A + Chain B + Ligand
```

### Frontend Demo (AI Assistant)

Trigger the demo by entering "AZOB" in the AI Design Assistant:

1. Displays ligand analysis (molecular weight, rings, symmetry)
2. Interview questions for design preferences
3. Two-step execution with progress indicators
4. Evaluation showing contacts_a, contacts_b, affinity, separability
5. 3D viewer with cartoon representation

### Achieved Results (January 2026)

| Metric | Target | Achieved | Notes |
|--------|--------|----------|-------|
| contacts_a | ≥3 | **26** | Excellent contact count |
| contacts_b | ≥3 | **20** | Strong binding |
| has_clashes | false | **false** | Clean geometry |
| separable | true | **true** | Topology verified |
| affinity | < -3 kcal/mol | **-9.78 kcal/mol** | Drug-like binding |

**Comparison with Previous Approaches:**

| Approach | Best Affinity | Contacts A | Contacts B | Separable |
|----------|---------------|------------|------------|-----------|
| Post-symmetry (80-100) | -3.78 kcal/mol | 6 | 4 | NO |
| Cleavable Monomer | -7.05 kcal/mol | 21 | 0 | NO |
| **Independent Design** | **-9.78 kcal/mol** | **26** | **20** | **YES** |

### Testing Locally (Requires GPU + Checkpoints)

```bash
# Build Docker image
cd backend/serverless
docker build -t foundry-rfd3-serverless .

# Run with GPU and checkpoint volume
docker run --gpus all \
  -v /path/to/checkpoints:/runpod-volume/checkpoints \
  foundry-rfd3-serverless \
  python -c "
from handler import handle_interface_ligand_design
result = handle_interface_ligand_design({
    'ligand_smiles': 'c1ccc(cc1)N=Nc2ccccc2',
    'approach': 'full',
    'chain_length': '60-80',
    'num_designs': 1,
    'seed': 42,
})
print(result)
"
```

**Note:** Local testing requires:
- NVIDIA GPU with CUDA
- RFD3 model checkpoints (~10GB)
- Foundry installation

For testing without GPU, deploy to RunPod Serverless.
