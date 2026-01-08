# RFD3 Advanced Features Research & Implementation Roadmap

## Executive Summary

This document analyzes the advanced features available in RosettaCommons Foundry RFD3 (RFdiffusion3) and identifies opportunities to make these capabilities accessible to users with general protein engineering knowledge.

**Current State:** Our tool exposes basic contig specification, number of designs, seed control, and PDB upload.

**Opportunity:** RFD3 offers extensive advanced features including symmetry design, hotspot binding, partial diffusion refinement, classifier-free guidance, and sophisticated selection mechanisms that could significantly expand the tool's capabilities.

---

## Feature Gap Analysis

### Currently Implemented (Our Tool)

| Feature | Implementation | UI Element |
|---------|---------------|------------|
| De novo design | ✅ Full | Length input (e.g., "100") |
| Binder design | ✅ Basic | Contig string (e.g., "A1-100/0 50-100") |
| Scaffold design | ✅ Basic | Contig string with fixed regions |
| Number of designs | ✅ Full | 1-10 slider |
| Seed control | ✅ Full | Checkbox + number input |
| PDB input | ✅ Full | Drag & drop file upload |
| Visual contig builder | ✅ Full | Interactive segment editor |
| Preset examples | ✅ Basic | Quick example buttons |

### Not Yet Implemented (High-Value Features)

| Feature | RFD3 Support | User Value | Implementation Complexity |
|---------|-------------|------------|--------------------------|
| **Symmetry Design** | Full | Very High | Medium |
| **Hotspot Residues** | Full | Very High | Low |
| **Partial Diffusion** | Full | High | Low |
| **Diffusion Parameters** | Full | Medium | Low |
| **Classifier-Free Guidance** | Full | Medium | Medium |
| **Ligand Binding** | Full | Very High | Medium |
| **Selection Mini-Language** | Full | High | Medium |
| **RASA-based Selection** | Full | Medium | Medium |
| **H-bond Specifications** | Full | Medium | Medium |

---

## Detailed Feature Specifications

### 1. Symmetry Design (Priority: HIGH)

Enables design of symmetric oligomeric proteins - critical for enzyme engineering, viral capsids, and protein assemblies.

#### Supported Symmetry Types

| Type | Subunits | Use Cases |
|------|----------|-----------|
| **Cyclic (Cn)** | n | Homo-oligomers (C2 dimers, C3 trimers, C4 tetramers, etc.) |
| **Dihedral (Dn)** | 2n | Higher-order assemblies with two rotation axes |
| **Tetrahedral (T)** | 12 | Cage-like structures |
| **Octahedral (O)** | 24 | Cubic symmetry, larger cages |
| **Icosahedral (I)** | 60 | Viral capsid-like structures |

#### Configuration Schema
```json
{
  "symmetry": {
    "id": "C3",
    "is_unsym_motif": "Y1-11,Z16-25",
    "is_symmetric_motif": true
  }
}
```

#### Technical Notes
- Must set `diffusion_batch_size=1` (memory constraint)
- Use `low_memory_mode=True` for T, O, I symmetries
- Sampler kind must be set to "symmetry"

#### Proposed UI
```
[Symmetry Design Panel]
┌─────────────────────────────────────────────────────────────┐
│ Symmetry Type: [Cyclic ▾]                                   │
│                                                             │
│ Number of Subunits: [3] ← (auto-calculates based on type)   │
│                                                             │
│ ○ C2 (Dimer)   ○ C3 (Trimer)   ○ C4 (Tetramer)            │
│ ○ C6 (Hexamer) ○ D2 (4 chains) ○ D3 (6 chains)            │
│                                                             │
│ □ Advanced: Asymmetric motifs (optional)                    │
│   Chains to exclude from symmetry: [____]                   │
│                                                             │
│ ℹ️ Memory note: Higher symmetry requires more GPU memory    │
└─────────────────────────────────────────────────────────────┘
```

---

### 2. Hotspot Residue Selection (Priority: HIGH)

Define specific residues that must be within binding distance (~4.5Å) of the designed structure - essential for functional binder design.

#### Configuration
```json
{
  "select_hotspots": {
    "A15-20": "ALL",
    "A25": "TIP",
    "A30-35": "BKBN"
  }
}
```

#### Selection Shorthands
- `"ALL"` - All atoms in residue
- `"TIP"` - Tip/functional atom only (e.g., Lys NH3+)
- `"BKBN"` - Backbone atoms only (N, CA, C, O)
- Custom: `"N,CA,C,O,CB"` - Specific atoms

#### Proposed UI
```
[Hotspot Residue Panel]
┌─────────────────────────────────────────────────────────────┐
│ Hotspot Residues (design must contact these)                │
│                                                             │
│ Chain: [A ▾]  Residues: [15-20] Atoms: [All ▾]  [+ Add]    │
│                                                             │
│ Current Hotspots:                                          │
│ ┌─────────────────────────────────────────────────────────┐│
│ │ A:15-20 (All atoms)                              [×]    ││
│ │ A:25 (Tip atom)                                  [×]    ││
│ │ A:30-35 (Backbone)                               [×]    ││
│ └─────────────────────────────────────────────────────────┘│
│                                                             │
│ ℹ️ Selected residues will be within ~4.5Å of designed      │
│    protein for functional binding interactions              │
└─────────────────────────────────────────────────────────────┘
```

---

### 3. Partial Diffusion / Structure Refinement (Priority: HIGH)

Start from existing structure and refine it - perfect for optimizing existing designs or experimental structures.

#### Configuration
```json
{
  "partial_t": 10.0
}
```

#### Parameter Range
- **5-10 Å**: Minor refinement, preserve most of original structure
- **10-15 Å**: Moderate refinement, more flexibility
- **15-20 Å**: Major redesign while maintaining topology

#### Proposed UI
```
[Partial Diffusion Panel]
┌─────────────────────────────────────────────────────────────┐
│ □ Enable Partial Diffusion (refine existing structure)      │
│                                                             │
│ Noise Level: [10.0 Å] ─────●────────────────                │
│              5Å          10Å           15Å          20Å     │
│              Minor       Moderate      Major                │
│              Refinement  Refinement    Redesign             │
│                                                             │
│ ℹ️ Lower values preserve more of the input structure        │
│    Higher values allow more structural changes              │
│                                                             │
│ CA-RMSD to input will be logged for assessment              │
└─────────────────────────────────────────────────────────────┘
```

---

### 4. Diffusion Sampling Parameters (Priority: MEDIUM)

Fine-tune the diffusion process for specific design goals.

#### Key Parameters

| Parameter | Default | Range | Effect |
|-----------|---------|-------|--------|
| `num_timesteps` | 200 | 50-500 | More steps = higher quality, slower |
| `step_scale` | 1.5 | 0.5-3.0 | Higher = less diverse, more designable |
| `noise_scale` | 1.003 | 0.9-1.1 | Diversity during inference |
| `gamma_0` | 0.6 | 0.3-1.0 | Lower = more designable, less diverse |

#### Proposed UI
```
[Advanced Sampling Options]
┌─────────────────────────────────────────────────────────────┐
│ ▼ Sampling Parameters                                       │
│                                                             │
│ Diversity vs Designability: ──────────●───────────          │
│                             Diverse        Designable       │
│                             (creative)     (realistic)      │
│                                                             │
│ Quality vs Speed: ──────●─────────────────────              │
│                   Fast    Balanced      High Quality        │
│                   (50)    (200)         (500 steps)         │
│                                                             │
│ □ Show detailed parameters                                  │
│   • num_timesteps: [200]                                    │
│   • step_scale: [1.5]                                       │
│   • noise_scale: [1.003]                                    │
│   • gamma_0: [0.6]                                          │
└─────────────────────────────────────────────────────────────┘
```

---

### 5. Classifier-Free Guidance (Priority: MEDIUM)

Use guidance to steer designs toward specific properties without explicit classifiers.

#### Configuration
```json
{
  "inference_sampler": {
    "use_classifier_free_guidance": true,
    "cfg_scale": 1.5,
    "cfg_t_max": 0.8,
    "cfg_features": ["active_donor", "active_acceptor", "ref_atomwise_rasa"]
  }
}
```

#### Proposed UI
```
[Guidance Options]
┌─────────────────────────────────────────────────────────────┐
│ □ Enable Classifier-Free Guidance                           │
│                                                             │
│ Guidance Strength: [1.5] ─────●────────────                 │
│                    Weak      Medium     Strong              │
│                                                             │
│ Guide toward:                                               │
│ □ H-bond donors (polar surfaces)                            │
│ □ H-bond acceptors (binding interfaces)                     │
│ □ Surface accessibility patterns                            │
│                                                             │
│ ℹ️ Guidance helps bias designs toward desired properties    │
└─────────────────────────────────────────────────────────────┘
```

---

### 6. Ligand Binding Design (Priority: HIGH)

Design proteins that bind specific small molecules or cofactors.

#### Configuration
```json
{
  "ligand": "ATP",
  "select_hotspots": {
    "ligand": "ALL"
  }
}
```

#### Supported Ligand Types
- Chemical component IDs (e.g., "ATP", "NAD", "HEM")
- Custom ligand from input structure
- Metal ions (e.g., "ZN", "MG", "CA")

#### Proposed UI
```
[Ligand Binding Design]
┌─────────────────────────────────────────────────────────────┐
│ □ Include Ligand in Design                                  │
│                                                             │
│ Ligand Source:                                              │
│ ○ From input structure (auto-detected)                      │
│ ○ Standard ligand: [ATP ▾]                                  │
│     Common: ATP, NAD, FAD, HEM, ZN, MG, CA                 │
│                                                             │
│ Binding Requirements:                                       │
│ □ Coordinate all ligand atoms                               │
│ □ Focus on functional groups only                           │
│                                                             │
│ ℹ️ Designed protein will have binding site for ligand       │
└─────────────────────────────────────────────────────────────┘
```

---

### 7. RASA-Based Selection (Priority: MEDIUM)

Select residues based on solvent accessibility - useful for designing interfaces vs core regions.

#### Configuration
```json
{
  "select_buried": true,
  "select_partially_buried": "A1-50",
  "select_exposed": false
}
```

#### Categories
- **Buried**: Core residues, low solvent accessibility
- **Partially Buried**: Interface residues
- **Exposed**: Surface residues

#### Proposed UI
```
[Surface Accessibility Selection]
┌─────────────────────────────────────────────────────────────┐
│ Design based on solvent accessibility:                      │
│                                                             │
│ ○ All residues (no filter)                                  │
│ ○ Core only (buried residues)                               │
│ ○ Interface (partially buried)                              │
│ ○ Surface only (exposed)                                    │
│ ○ Custom selection per region                               │
│                                                             │
│ Preview: [Visual showing selected regions]                  │
└─────────────────────────────────────────────────────────────┘
```

---

### 8. H-Bond Specifications (Priority: MEDIUM)

Define hydrogen bonding requirements for designed interfaces.

#### Configuration
```json
{
  "select_hbond_donor": "A15-20",
  "select_hbond_acceptor": "A25-30"
}
```

#### Proposed UI
```
[Hydrogen Bonding Requirements]
┌─────────────────────────────────────────────────────────────┐
│ □ Specify H-bond requirements                               │
│                                                             │
│ H-bond Donors (must donate):                               │
│ Chain: [A ▾]  Residues: [15-20]                  [+ Add]   │
│                                                             │
│ H-bond Acceptors (must accept):                            │
│ Chain: [A ▾]  Residues: [25-30]                  [+ Add]   │
│                                                             │
│ ℹ️ Ensures proper hydrogen bonding at interface             │
└─────────────────────────────────────────────────────────────┘
```

---

## Implementation Recommendations

### Phase 1: Quick Wins (Low Complexity, High Value)

1. **Hotspot Residue Selection**
   - Add simple residue range input
   - Map to `select_hotspots` config
   - Estimated: 2-3 hours

2. **Partial Diffusion**
   - Add slider for `partial_t` parameter
   - Show/hide based on PDB upload
   - Estimated: 1-2 hours

3. **Diffusion Presets**
   - "Quick & Diverse", "Balanced", "High Quality" presets
   - Map to num_timesteps, step_scale, gamma_0
   - Estimated: 2 hours

### Phase 2: Major Features (Medium Complexity)

4. **Symmetry Design Panel**
   - Add symmetry type selector
   - Handle batch_size=1 constraint
   - Add visual preview of symmetry
   - Estimated: 4-6 hours

5. **Ligand Binding Support**
   - Auto-detect ligands from input PDB
   - Add standard ligand dropdown
   - Configure ligand hotspots
   - Estimated: 4-6 hours

6. **Classifier-Free Guidance**
   - Add guidance toggle and scale slider
   - Feature selection checkboxes
   - Estimated: 3-4 hours

### Phase 3: Advanced Features (Higher Complexity)

7. **Full Selection Mini-Language**
   - Visual selection builder
   - Atom-level specificity
   - Preview of selections
   - Estimated: 6-8 hours

8. **RASA & H-Bond Selection**
   - Surface accessibility calculator/viewer
   - H-bond network visualization
   - Estimated: 6-8 hours

---

## User Experience Guidelines

### Making Advanced Features Accessible

1. **Progressive Disclosure**
   - Basic mode: Just contig + number of designs
   - Intermediate: Show preset buttons for common workflows
   - Advanced: Collapsible panels for detailed parameters

2. **Context-Aware Options**
   - Only show relevant options based on design type
   - Hide symmetry for binder design (usually asymmetric)
   - Show partial diffusion only when PDB is uploaded

3. **Helpful Defaults**
   - Pre-populate with sensible defaults
   - Show recommended ranges inline
   - Provide tooltips with scientific context

4. **Visual Feedback**
   - Preview symmetry with 3D visualization
   - Highlight hotspot residues in structure viewer
   - Show affected regions in contig visualization

5. **Workflow Templates**
   - "Design symmetric enzyme" preset
   - "Optimize existing binder" preset
   - "Add ligand binding site" preset

---

## Backend Changes Required

### RFD3Request Model Updates

```python
class RFD3Request(BaseModel):
    """RFdiffusion3 design request - extended"""
    contig: Optional[str] = None
    num_designs: int = Field(default=1, ge=1, le=10)
    pdb_content: Optional[str] = None
    seed: Optional[int] = None

    # NEW: Symmetry
    symmetry: Optional[Dict[str, Any]] = None  # {id, is_unsym_motif, is_symmetric_motif}

    # NEW: Partial diffusion
    partial_t: Optional[float] = Field(default=None, ge=0.0, le=30.0)

    # NEW: Hotspots
    select_hotspots: Optional[Dict[str, str]] = None  # {"A15-20": "ALL", ...}

    # NEW: Diffusion parameters
    num_timesteps: Optional[int] = Field(default=200, ge=50, le=500)
    step_scale: Optional[float] = Field(default=1.5, ge=0.5, le=3.0)
    noise_scale: Optional[float] = Field(default=1.003, ge=0.9, le=1.1)
    gamma_0: Optional[float] = Field(default=0.6, ge=0.3, le=1.0)

    # NEW: Classifier-free guidance
    use_cfg: Optional[bool] = False
    cfg_scale: Optional[float] = Field(default=1.5, ge=0.5, le=3.0)
    cfg_features: Optional[List[str]] = None

    # NEW: Ligand
    ligand: Optional[str] = None

    # NEW: Selection options
    select_buried: Optional[bool] = None
    select_exposed: Optional[bool] = None
    select_hbond_donor: Optional[str] = None
    select_hbond_acceptor: Optional[str] = None
```

### Engine Configuration Mapping

```python
def build_rfd3_config(request: RFD3Request) -> Dict[str, Any]:
    """Build RFD3 inference configuration from request"""
    config = {}

    # Base specification
    spec = {}
    if request.contig:
        if request.contig.replace("-", "").isdigit():
            spec["length"] = request.contig
        else:
            spec["contig"] = request.contig

    # Hotspots
    if request.select_hotspots:
        spec["select_hotspots"] = request.select_hotspots

    # Partial diffusion
    if request.partial_t is not None:
        spec["partial_t"] = request.partial_t

    # Ligand
    if request.ligand:
        spec["ligand"] = request.ligand

    # Symmetry
    if request.symmetry:
        spec["symmetry"] = request.symmetry
        config["diffusion_batch_size"] = 1  # Required for symmetry
        if request.symmetry.get("id", "").startswith(("T", "O", "I")):
            config["low_memory_mode"] = True
    else:
        config["diffusion_batch_size"] = request.num_designs

    # Sampling parameters
    sampler = {}
    if request.num_timesteps != 200:
        sampler["num_timesteps"] = request.num_timesteps
    if request.step_scale != 1.5:
        sampler["step_scale"] = request.step_scale
    if request.noise_scale != 1.003:
        sampler["noise_scale"] = request.noise_scale
    if request.gamma_0 != 0.6:
        sampler["gamma_0"] = request.gamma_0

    # CFG
    if request.use_cfg:
        sampler["use_classifier_free_guidance"] = True
        sampler["cfg_scale"] = request.cfg_scale
        if request.cfg_features:
            sampler["cfg_features"] = request.cfg_features

    if sampler:
        config["inference_sampler"] = sampler

    config["specification"] = spec
    return config
```

---

## API Changes Summary

### New API Types (api.ts)

```typescript
export interface RFD3Request {
  contig: string;
  num_designs?: number;
  pdb_content?: string;
  seed?: number;

  // Symmetry
  symmetry?: {
    id: string;  // "C2", "C3", "D2", "T", "O", "I"
    is_unsym_motif?: string;
    is_symmetric_motif?: boolean;
  };

  // Partial diffusion
  partial_t?: number;

  // Hotspots
  select_hotspots?: Record<string, string>;

  // Sampling
  num_timesteps?: number;
  step_scale?: number;
  noise_scale?: number;
  gamma_0?: number;

  // CFG
  use_cfg?: boolean;
  cfg_scale?: number;
  cfg_features?: string[];

  // Ligand
  ligand?: string;

  // Selections
  select_buried?: boolean;
  select_exposed?: boolean;
  select_hbond_donor?: string;
  select_hbond_acceptor?: string;
}
```

---

## References

- [RFD3 Input Documentation](https://github.com/RosettaCommons/foundry/blob/production/models/rfd3/docs/input.md)
- [RFD3 Symmetry Documentation](https://github.com/RosettaCommons/foundry/blob/production/models/rfd3/docs/symmetry.md)
- IPD Official Pipeline Examples

---

## Next Steps

1. Review this document with stakeholders
2. Prioritize features based on user needs
3. Start with Phase 1 quick wins
4. Implement backend configuration mapping
5. Design and implement UI components
6. Add comprehensive tooltips and documentation
7. Test with real-world design scenarios
