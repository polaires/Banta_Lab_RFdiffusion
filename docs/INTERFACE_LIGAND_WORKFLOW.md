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
node demo_azobenzene_interface.js

# Run demo with binder mode (H-bond conditioning)
node demo_azobenzene_interface.js binder

# Customize number of designs
node demo_azobenzene_interface.js natural 5
node demo_azobenzene_interface.js binder 10
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
