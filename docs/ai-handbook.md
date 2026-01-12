# RFdiffusion AI Handbook

**Version:** 1.0.0
**Purpose:** Comprehensive knowledge base for AI-driven protein engineering

This handbook provides the context and rules for an AI system to:
1. Identify what type of protein design task the user wants
2. Generate appropriate RFD3 parameters
3. Explain decisions in plain English
4. Evaluate results and suggest improvements

---

## 1. TASK IDENTIFICATION

### 1.1 Supported Design Tasks

| Task Type | Keywords | Requires Input Structure | Primary Parameters |
|-----------|----------|--------------------------|-------------------|
| `de_novo` | "create", "design", "generate", "new protein" | No | length, symmetry |
| `binder` | "binder", "bind to", "target", "PPI", "interface" | Yes (target) | hotspots, step_scale |
| `metal_redesign` | "metal", "lanthanide", "coordination", metal names | Yes | partial_t, ligand |
| `enzyme` | "enzyme", "catalytic", "active site", "substrate" | Yes (template) | partial_t, hotspots |
| `refinement` | "refine", "optimize", "improve", "partial" | Yes | partial_t (low) |
| `symmetric` | "symmetric", "dimer", "trimer", "oligomer", "Cn", "Dn" | Optional | symmetry config |
| `scaffold` | "scaffold", "motif", "functional site", "graft" | Yes (motif) | contig, partial_t |

### 1.2 Intent Recognition Patterns

```
"Design a 100 residue protein" → de_novo, length=100
"Create a binder for HER2" → binder, need target structure
"Convert iron site to terbium" → metal_redesign, FE→TB
"Make this bind calcium" → metal_redesign, target=CA
"Improve stability" → refinement, low partial_t
"Design a C3 symmetric trimer" → symmetric, symmetry=C3
"Refine this structure" → refinement
"Graft EF-hand motif" → scaffold
```

### 1.3 Ambiguity Handling

When intent is unclear, ask clarifying questions about:
- **What**: Specific design goal (binding? stability? function?)
- **How much change**: Conservative vs aggressive redesign
- **Constraints**: What must be preserved?

---

## 2. PARAMETER KNOWLEDGE BASE

### 2.1 Core RFD3 Parameters

| Parameter | Type | Range | Default | Effect | Guidance |
|-----------|------|-------|---------|--------|----------|
| `contig` | string | - | required | Defines chains/residues | Format: "A1-50,0 B1-30" |
| `length` | string | "50-200" | - | De novo protein length | Use for no input structure |
| `partial_t` | float | 0-30 | 10 | Noise level (Angstroms) | Higher = more change |
| `num_timesteps` | int | 50-500 | 200 | Denoising steps | Higher = better quality |
| `step_scale` | float | 0.5-3.0 | 1.5 | Step size multiplier | Higher = more diversity |
| `gamma_0` | float | 0.1-1.0 | 0.6 | Self-conditioning strength | Higher = more stable |
| `noise_scale` | float | 0.9-1.1 | 1.0 | Noise amplitude | Usually leave at 1.0 |
| `num_designs` | int | 1-20 | 5 | Output count | More = better sampling |
| `seed` | int | 0-999999 | random | Random seed | Set for reproducibility |
| `ligand` | string | - | - | Metal/ligand ID | "TB", "ZN", "ATP", etc. |

### 2.2 Hotspot Parameters

| Parameter | Type | Format | Effect |
|-----------|------|--------|--------|
| `select_hotspots` | dict | {"A15-20": "ALL"} | Fix residues in binding |
| `unindex` | string | "A15,A20,A25" | Residues to redesign |
| `select_fixed_atoms` | dict | {"A15": "BKBN"} | Fix backbone only |

**Atom Selection Values:**
- `ALL`: Fix entire residue
- `BKBN`: Fix backbone (N, CA, C, O)
- `TIP`: Fix sidechain tip only
- `SC`: Fix sidechain

### 2.3 Symmetry Parameters

| Symmetry | Subunits | Use Case | Memory |
|----------|----------|----------|--------|
| C2 | 2 | Homodimer | Low |
| C3 | 3 | Trimer | Medium |
| C4-C6 | 4-6 | Higher cyclic | Medium |
| D2 | 4 | Dihedral | Medium |
| D3 | 6 | Dihedral | High |
| T | 12 | Tetrahedral | Very High |
| O | 24 | Octahedral | Very High |
| I | 60 | Icosahedral | Extreme |

**For T, O, I:** Enable `low_memory_mode: true`

### 2.4 Classifier-Free Guidance (CFG)

```json
{
  "cfg": {
    "enabled": true,
    "scale": 1.5,      // 0.5-3.0, higher = more guided
    "t_max": 0.8,      // 0-1, fraction of steps with CFG
    "features": ["active_donor", "active_acceptor"]
  }
}
```

---

## 3. METAL CHEMISTRY DATABASE

### 3.1 Metal Classification (HSAB Theory)

**Hard Acids** (prefer O donors):
- Lanthanides: La, Ce, Pr, Nd, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu
- Alkaline earth: Ca, Mg, Sr, Ba
- Group 3: Sc, Y

**Borderline** (flexible):
- Fe, Zn, Cu, Co, Ni, Mn

**Soft Acids** (prefer S donors):
- Ag, Au, Hg, Cd, Pt, Pd

### 3.2 Metal Properties Reference

| Metal | Coordination | Geometry | Preferred Donors | Bond Distance | Notes |
|-------|--------------|----------|------------------|---------------|-------|
| TB | 8-9 | Square antiprism | Asp, Glu, Asn, Gln | 2.3-2.5 Å | Luminescent |
| GD | 8-9 | Square antiprism | Asp, Glu, Asn | 2.3-2.5 Å | MRI contrast |
| EU | 8-9 | Square antiprism | Asp, Glu | 2.3-2.6 Å | Red luminescence |
| CA | 6-8 | Octahedral/PBP | Asp, Glu, Asn | 2.3-2.6 Å | EF-hand signaling |
| MG | 6 | Octahedral | Asp, Glu | 2.0-2.3 Å | ATP binding |
| ZN | 4-6 | Tetrahedral | His, Cys, Asp | 1.9-2.4 Å | Zinc fingers |
| FE | 4-6 | Tetrahedral/Oct | Cys, His, Asp | 1.9-2.5 Å | Iron-sulfur |
| CU | 4-6 | Square planar | His, Cys, Met | 1.9-2.3 Å | Electron transfer |
| MN | 4-6 | Tetrahedral/Oct | Asp, His | 2.0-2.4 Å | OEC, SOD |

### 3.3 Metal Conversion Strategies

**Transition Metal → Lanthanide (e.g., Fe→Tb):**
- Increase coordination: +4 ligands typically needed
- Replace S donors with O donors
- Expand binding pocket
- **Parameters:** partial_t=15-18, step_scale=1.8

**Same-class swap (e.g., Zn→Fe):**
- Maintain coordination number
- Adjust donor types if needed
- **Parameters:** partial_t=8-12

**Lanthanide → Lanthanide (e.g., Tb→Gd):**
- Minimal changes needed
- Same coordination preferences
- **Parameters:** partial_t=5-8

---

## 4. TASK-SPECIFIC RECIPES

### 4.1 De Novo Design

**When:** No input structure, designing from scratch

```json
{
  "length": "80-100",
  "num_designs": 10,
  "num_timesteps": 200,
  "step_scale": 1.5,
  "gamma_0": 0.6
}
```

**Quality preset:**
- Fast: num_timesteps=50, num_designs=3
- Balanced: num_timesteps=200, num_designs=10
- High: num_timesteps=500, num_designs=20

### 4.2 Protein Binder Design

**When:** Designing a protein to bind a target

```json
{
  "contig": "A1-100,0 60-80",  // Target chain + binder length
  "select_hotspots": {"A45-55": "ALL"},  // Interface residues
  "step_scale": 3.0,  // High for diversity
  "gamma_0": 0.2,  // Low for exploration
  "num_designs": 20
}
```

**Key considerations:**
- Identify target hotspot residues first
- Use high step_scale for diverse sampling
- Low gamma_0 allows more novel solutions

### 4.3 Metal Binding Redesign

**When:** Converting metal binding site to new metal

1. Analyze current coordination (coordination_analysis endpoint)
2. Calculate coordination delta needed
3. Select strategy based on delta:

| Coordination Change | Strategy | partial_t |
|---------------------|----------|-----------|
| +4 or more | expand_pocket | 18-20 |
| +2 to +3 | increase_coordination | 15-18 |
| -1 to +1 | refine_geometry | 5-10 |
| -2 to -3 | decrease_coordination | 8-12 |

```json
{
  "ligand": "TB",
  "partial_t": 15.0,
  "contig": "A1-100",
  "unindex": "A15,A20,A25,A30",  // Coordinating residues
  "select_fixed_atoms": {"A15": "BKBN", "A20": "BKBN"},
  "num_timesteps": 200,
  "step_scale": 1.5
}
```

### 4.4 Enzyme Active Site Design

**When:** Designing catalytic function

```json
{
  "contig": "A1-150",
  "ligand": "substrate_id",
  "select_hotspots": {"A45": "ALL", "A90": "ALL"},  // Catalytic residues
  "partial_t": 12.0,
  "step_scale": 1.3,
  "gamma_0": 0.7
}
```

### 4.5 Structure Refinement

**When:** Improving an existing design

```json
{
  "contig": "A1-100",
  "partial_t": 5.0,  // Low noise
  "num_timesteps": 150,
  "step_scale": 1.2,
  "gamma_0": 0.8  // High for stability
}
```

### 4.6 Symmetric Oligomer

**When:** Designing homo-oligomers

```json
{
  "length": "60-80",
  "symmetry": {
    "id": "C3",
    "is_symmetric_motif": true
  },
  "num_timesteps": 250,
  "step_scale": 1.5
}
```

---

## 5. PARAMETER DECISION LOGIC

### 5.1 partial_t Selection

```
IF task == "refinement":
    partial_t = 5-8
ELIF task == "metal_redesign":
    IF coordination_delta > 3:
        partial_t = 15-20
    ELIF coordination_delta > 0:
        partial_t = 10-15
    ELSE:
        partial_t = 8-12
ELIF task == "binder":
    partial_t = 10-15 (if partial diffusion mode)
ELIF task == "de_novo":
    partial_t = None (not used)
```

### 5.2 Quality vs Speed Trade-offs

| Priority | num_timesteps | step_scale | gamma_0 | num_designs |
|----------|---------------|------------|---------|-------------|
| Speed | 50-100 | 1.8 | 0.5 | 3-5 |
| Balanced | 200 | 1.5 | 0.6 | 10 |
| Quality | 300-500 | 1.2 | 0.7 | 15-20 |
| Exploration | 200 | 2.5-3.0 | 0.2-0.4 | 20 |

### 5.3 User Intent Mapping

| User Says | Parameter Adjustments |
|-----------|----------------------|
| "conservative" | partial_t-=3, gamma_0+=0.1 |
| "aggressive" | partial_t+=5, step_scale+=0.5 |
| "stable" | gamma_0=0.8, num_timesteps+=100 |
| "diverse" | step_scale=2.5, gamma_0=0.3 |
| "quick" | num_timesteps=100, num_designs=5 |
| "high quality" | num_timesteps=400, num_designs=15 |

---

## 6. EVALUATION CRITERIA

### 6.1 Per-Task Success Metrics

**Metal Binding:**
- Coordination number matches target range
- Geometry matches preferred type
- Bond distances within range
- Correct donor types (no avoided donors)

**Binder Design:**
- Interface SASA > 800 Å²
- Shape complementarity > 0.65
- No clashes at interface

**De Novo:**
- pLDDT > 0.7
- No significant clashes
- Reasonable Ramachandran

**Refinement:**
- RMSD from input < 2.0 Å
- pLDDT improved or maintained
- No new clashes

### 6.2 Designability Assessment

Use RF3 (RosettaFold3) to refold designed sequences:

| RMSD | Verdict | Meaning |
|------|---------|---------|
| < 1.0 Å | Excellent | Will fold as designed |
| 1.0-2.0 Å | Good | Likely to fold correctly |
| 2.0-3.0 Å | Moderate | Some deviation expected |
| > 3.0 Å | Poor | May not fold as intended |

### 6.3 Failure Modes & Remediation

| Failure | Cause | Fix |
|---------|-------|-----|
| Low coordination | partial_t too low | Increase by 3-5 |
| Wrong geometry | Not enough sampling | More designs, higher step_scale |
| Clashes | gamma_0 too low | Increase to 0.7+ |
| Too different | partial_t too high | Decrease by 3-5 |
| No metal | Pocket too small | Increase partial_t, expand |

---

## 7. CONVERSATION PATTERNS

### 7.1 Task Detection Response

When task type is detected, respond with:
1. Detected task type and confidence
2. Required inputs (structure, parameters)
3. Suggested approach with reasoning
4. Any clarifying questions needed

### 7.2 Parameter Explanation Template

"I'm setting **{parameter}** to **{value}** because {reason}.
This will {effect}. If you want {alternative_outcome}, I can adjust to {alternative_value}."

### 7.3 Results Interpretation

After design completes, explain:
1. Key metrics and what they mean
2. Whether targets were met
3. Suggestions for improvement if needed
4. Next steps (validation, MPNN, etc.)

---

## 8. WORKFLOW INTEGRATION

### 8.1 Complete Design Pipeline

1. **Input Phase**: Get structure (PDB ID or upload) + user goal
2. **Analysis Phase**: Analyze structure, detect metals, binding sites
3. **Planning Phase**: Determine task type, generate parameters
4. **Confirmation Phase**: Show plan to user, allow adjustments
5. **Execution Phase**: Run RFD3 job
6. **Evaluation Phase**: Analyze results against criteria
7. **Iteration Phase**: Suggest refinements if needed

### 8.2 Multi-Step Workflows

**Metal Redesign + Validation:**
1. RFD3 with metal binding → Structure
2. MPNN on structure → Sequences
3. RF3 refold top sequences → Validation
4. RMSD check → Designability

**Binder Design + Optimization:**
1. RFD3 binder design → Initial binder
2. Evaluate interface metrics
3. RFD3 refinement if needed → Optimized binder
4. MPNN sequence design → Final sequences

---

## 9. ERROR HANDLING

### 9.1 Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| "Missing contig or length" | No structure spec | Add contig or length param |
| "Symmetry batch_size must be 1" | Memory | Set batch_size=1 |
| "Metal not found" | Wrong residue number | Check PDB for metal position |
| "GPU OOM" | Large structure | Enable low_memory_mode |

### 9.2 Graceful Degradation

If optimal params fail:
1. Try with lower num_timesteps
2. Reduce num_designs
3. Enable low_memory_mode
4. Switch to conservative strategy

---

## 10. APPENDIX

### 10.1 Three-Letter Metal Codes

```
TB = Terbium      GD = Gadolinium   EU = Europium
CA = Calcium      MG = Magnesium    ZN = Zinc
FE = Iron         CU = Copper       MN = Manganese
CO = Cobalt       NI = Nickel       LA = Lanthanum
```

### 10.2 Common Ligand Codes

```
ATP = Adenosine triphosphate
NAD = Nicotinamide adenine dinucleotide
FAD = Flavin adenine dinucleotide
HEM = Heme
```

### 10.3 Residue Donor Classifications

**Oxygen donors (hard):** Asp (OD1/OD2), Glu (OE1/OE2), Asn (OD1), Gln (OE1), Ser (OG), Thr (OG1), Tyr (OH)

**Nitrogen donors:** His (ND1/NE2), Lys (NZ), Arg (NH1/NH2)

**Sulfur donors (soft):** Cys (SG), Met (SD)

---

*End of Handbook*
