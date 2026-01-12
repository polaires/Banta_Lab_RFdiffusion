# RFD3 Design System: Comprehensive Instruction Guide

## Purpose of This Document

This document helps you (Claude Code) understand RFD3's **full capability space** before implementing any protein design workflow. The goal is to prevent a common failure mode: implementing features in isolation when they should be combined.

---

## RFD3 Core Capabilities (From Paper)

### What RFD3 Does
RFD3 (RoseTTAFold Diffusion 3) is a generative model that creates protein structures **conditioned on complex molecular contexts**:

- **Joint protein-ligand sampling**: Protein structure and ligand conformation are sampled TOGETHER (not rigid docking)
- **Multi-component design**: Can design around ligands, nucleic acids, metal ions, other proteins simultaneously
- **Atom-level constraints**: RASA (burial), H-bonds, fixed coordinates at individual atom resolution
- **Symmetry at inference**: C2, C3, C4, C5, C7, D2, D3 symmetries applied without retraining

### Key Performance Numbers
| Task | Success Rate | Notes |
|------|--------------|-------|
| Monomeric ligand binders | ~15-25% | AF3 validation |
| Dimeric DNA binders | 6.67% | Symmetry + nucleic acid |
| Enzyme motif scaffolding | Variable | Depends on motif complexity |
| Custom (non-CCD) ligands | Lower | Less training data |

### What Makes RFD3 Different from RFdiffusion (v1/v2)
| Feature | RFdiffusion v1/v2 | RFD3 |
|---------|-------------------|------|
| Ligand handling | Rigid (fixed position) | Co-diffused (flexible) |
| Atom-level control | Limited | Full (RASA, H-bonds) |
| Input format | Contig only | JSON/YAML with rich options |
| Symmetry + ligand | Separate | Combined in single run |

### Critical Insight: Features COMPOSE
The paper demonstrates that RFD3 features are **orthogonal and composable**:

```
symmetry + ligand           ✓ Works (dimeric binder around ligand)
symmetry + DNA              ✓ Works (6.67% success demonstrated)
ligand + RASA               ✓ Works (control burial depth)
symmetry + motif            ✓ Works (symmetric enzyme)
symmetry + ligand + RASA    ✓ Works (all three together)
```

**This means:** When a user asks for "dimer with ligand at interface", you MUST use symmetry AND ligand AND RASA together, not choose one.

---

## CRITICAL: Read This First

### The #1 Mistake to Avoid

```
❌ WRONG THINKING:
   User wants: "ligand binder" → Use ligand conditioning
   User wants: "dimer" → Use symmetry
   User wants: "ligand at interface of dimer" → ???
   
   Claude Code often does: Implement ligand binding (monomer) because 
   that's what the examples show.

✅ CORRECT THINKING:
   User wants: "ligand at interface of dimer"
   → This requires COMBINING symmetry + ligand conditioning
   → Search/verify: "Can RFD3 do symmetry + ligand together?"
   → Answer: YES, these are orthogonal features that compose
```

### Before Writing Any Code, Answer These Questions

1. **What is the user's actual end goal?** (Not "run RFD3" but "protein that does X")
2. **What RFD3 features are needed to achieve that goal?**
3. **Do these features need to be combined in a single run?**
4. **How will I validate the output matches the goal?**

---

## RFD3 Input Format (OFFICIAL)

### Basic Structure
RFD3 accepts JSON or YAML with **design names as top-level keys**:

```json
{
    "design-1": {
        "input": "path/to/structure.pdb",
        "contig": "50-80,/0,A1-100",
        "ligand": "HAX,OAA"
    },
    "design-2": {
        // ... independent configuration
    }
}
```

### CLI Command
```bash
rfd3 design out_dir=<path/to/outdir> inputs=<path/to/inputs.json>
```

### Key CLI Arguments
| Argument | Default | Description |
|----------|---------|-------------|
| `n_batches` | 1 | Number of batches per input key |
| `diffusion_batch_size` | 8 | Designs per batch |
| `inference_sampler.num_timesteps` | 200 | Diffusion steps |
| `inference_sampler.step_scale` | 1.5 | Step size scale (higher = less diverse, more designable) |
| `inference_sampler.kind` | "default" | Set to "symmetry" for symmetric designs |
| `inference_sampler.center_option` | "all" | Options: "all", "motif", "diffuse" |

---

## RFD3 Capability Matrix

### Core Features (Can Be Combined)

| Feature | Description | Key Parameters |
|---------|-------------|----------------|
| **Unconditional** | Generate novel folds | `contig: "50-100"` (no chain label = design) |
| **Binder design** | Bind to target protein | `contig: "50-70,\0,A1-200"` (A=target from input) |
| **Symmetric oligomers** | C2, C3, D2, etc. | `symmetry: {...}` (see symmetry docs) |
| **Ligand conditioning** | Design around small molecule | `ligand: "CCD_NAME"` or `ligand: "L1"` (by index) |
| **RASA conditioning** | Control atom burial | `select_buried`, `select_partially_buried`, `select_exposed` |
| **H-bond specification** | Force specific H-bonds | `select_hbond_donor`, `select_hbond_acceptor` |
| **Hotspots** | Target specific residues | `select_hotspots: "A42,A67"` |
| **Origin control** | Position protein COM | `ori_token: [x,y,z]` or `infer_ori_strategy: "hotspots"` |
| **Partial diffusion** | Refine existing structure | `partial_t: 15.0` (noise in Å, 5-15 recommended) |

### Feature Composition Rules

```
✅ VALID COMBINATIONS:
   - Symmetry + Ligand → Symmetric binder around ligand
   - Symmetry + Motif → Symmetric enzyme
   - Ligand + RASA → Control how buried ligand is
   - Binder + Hotspots → Targeted interface design
   - DNA + Symmetry → Dimeric DNA binder (6.67% success rate reported)

❌ CONFLICTING COMBINATIONS:
   - Unconditional + Binder (contradictory)
   - Multiple symmetry types simultaneously
   - Overlapping RASA selections (same atom can't be buried AND exposed)
```

---

## InputSpecification Fields (OFFICIAL)

| Field | Type | Description |
|-------|------|-------------|
| `input` | `str` | Path to PDB/CIF file |
| `contig` | `InputSelection` | Indexed motif spec, e.g., `"A1-80,10,\0,B5-12"` |
| `unindex` | `InputSelection` | Unindexed motif components (position unknown) |
| `length` | `str` | Total design length, e.g., `"50-100"` or `100` |
| `ligand` | `str` | Ligand(s) by CCD name or index, e.g., `"HAX,OAA"` or `"L1"` |
| `select_fixed_atoms` | `InputSelection` | Atoms with fixed coordinates |
| `select_unfixed_sequence` | `InputSelection` | Where sequence can change |
| `select_buried` | `InputSelection` | RASA: buried atoms |
| `select_partially_buried` | `InputSelection` | RASA: partially buried atoms |
| `select_exposed` | `InputSelection` | RASA: exposed atoms |
| `select_hbond_donor` | `InputSelection` | H-bond donor atoms |
| `select_hbond_acceptor` | `InputSelection` | H-bond acceptor atoms |
| `select_hotspots` | `InputSelection` | Target residues to contact |
| `symmetry` | `SymmetryConfig` | Symmetry specification |
| `ori_token` | `list[float]` | `[x,y,z]` origin override for COM placement |
| `infer_ori_strategy` | `str` | `"com"` or `"hotspots"` |
| `partial_t` | `float` | Noise (Å) for partial diffusion (5-15 recommended) |
| `is_non_loopy` | `bool` | Default True. Fewer loops, more helices |

---

## The InputSelection Mini-Language

Fields marked `InputSelection` accept **boolean**, **contig string**, or **dictionary**:

```yaml
# Boolean
select_fixed_atoms: true  # Fix all atoms from input

# Contig string
select_fixed_atoms: "A1-10,B5-8"  # Fix specific residues

# Dictionary (most expressive)
select_fixed_atoms:
  A1-2: BKBN      # Backbone only (N,CA,C,O)
  A3: N,CA,C,O,CB # Specific atoms
  B5-7: ALL       # All atoms
  B10: TIP        # Tip atom only
  LIG: ''         # UNFIX ligand atoms (empty string!)
```

**Shortcuts:**
- `ALL` = all atoms
- `BKBN` = backbone (N, CA, C, O)
- `TIP` = tip/functional atoms
- `''` (empty) = select nothing / unfix

---

## Contig String Format

```
contig: "A40-60,70,A120-170,\0,B3-45,60-80"
         │       │  │          │  │        │
         │       │  │          │  │        └─ Design 60-80 residues
         │       │  │          │  └─ Input residues B3-45
         │       │  │          └─ Chain break (no peptide bond)
         │       │  └─ Input residues A120-170
         │       └─ Design exactly 70 residues
         └─ Input residues A40-60
```

**Rules:**
- **With chain label** (e.g., `A40-60`): Residues from input structure
- **Without chain label** (e.g., `70` or `50-80`): Design new residues
- **`\0`**: Chain break (separate chains)
- **Comma**: Continue same chain
- **Range** (e.g., `50-80`): Random length within range

---

## Mapping User Goals to RFD3 Features

### Decision Tree

```
User Request
    │
    ├─► "Bind to [protein/DNA/ligand]"
    │       │
    │       ├─► Single chain? → Binder design mode
    │       │
    │       └─► Symmetric/multivalent? → Symmetry + Binder/Ligand
    │
    ├─► "Enzyme that catalyzes X"
    │       │
    │       └─► Motif scaffolding (unindex) + (optional) Symmetry
    │
    ├─► "Symmetric assembly"
    │       │
    │       ├─► With function? → Symmetry + [Ligand/Motif/Binder]
    │       │
    │       └─► Just structure? → Symmetry only
    │
    └─► "Ligand binder"
            │
            ├─► Monomer? → Ligand conditioning
            │
            ├─► Dimer with ligand at interface? 
            │       → Symmetry + Ligand + ori_token at ligand
            │
            └─► Ligand-induced dimerization?
                    → Two-state design (advanced)
```

### Common User Requests → RFD3 Configuration

#### 1. "Design a binder for azobenzene"
```json
{
  "azo_binder": {
    "input": "azobenzene.pdb",
    "contig": "60-80",
    "ligand": "AZO",
    "select_partially_buried": {"AZO": "ALL"}
  }
}
```

#### 2. "Design a DIMER that binds azobenzene at the INTERFACE"
```json
{
  "azo_c2_dimer": {
    "input": "azobenzene_centered.pdb",
    "contig": "50-70",
    "ligand": "AZO",
    "symmetry": {"type": "C2"},
    "select_partially_buried": {"AZO": "N1,N2"},
    "ori_token": [0.0, 0.0, 0.0],
    "infer_ori_strategy": "com"
  }
}
```
**Key differences:** 
- `symmetry: {"type": "C2"}` creates the dimer
- `ori_token` + ligand at origin centers dimer on ligand
- `select_partially_buried` keeps ligand at interface (not fully buried)

#### 3. "Design an enzyme with catalytic triad (unindexed motif)"
```json
{
  "enzyme_design": {
    "input": "active_site.pdb",
    "contig": "80-120",
    "unindex": "A42,A87,A123",
    "select_fixed_atoms": {
      "A42": "TIP",
      "A87": "TIP", 
      "A123": "TIP"
    }
  }
}
```

#### 4. "Symmetric enzyme (trimeric) with active site"
```json
{
  "symmetric_enzyme": {
    "input": "motif.pdb",
    "contig": "60-80",
    "symmetry": {"type": "C3"},
    "unindex": "A42,A87",
    "select_fixed_atoms": {"A42": "ALL", "A87": "ALL"}
  }
}
```

#### 5. "Binder for protein target with hotspot residues"
```json
{
  "protein_binder": {
    "input": "target.pdb",
    "contig": "60-80,\\0,A1-200",
    "select_hotspots": "A42,A67,A103",
    "infer_ori_strategy": "hotspots"
  }
}
```

---

## Validation: Does Output Match User Goal?

### Critical Validation Checks

**Don't just check "did RFD3 run successfully"—check "does output satisfy user intent":**

```python
def validate_design_meets_goal(design, user_goal):
    """
    ALWAYS implement goal-specific validation, not just structural metrics.
    """
    
    if user_goal == "ligand_at_interface":
        # Check ligand contacts BOTH chains
        contacts_A = get_contacts(design, ligand_chain="L", protein_chain="A")
        contacts_B = get_contacts(design, ligand_chain="L", protein_chain="B")
        
        if len(contacts_A) < 3:
            return False, "Ligand doesn't contact chain A sufficiently"
        if len(contacts_B) < 3:
            return False, "Ligand doesn't contact chain B sufficiently"
        
        # Check ligand is at interface, not buried in one chain
        ligand_sasa = calculate_sasa(design, selection="ligand")
        if ligand_sasa < 5:
            return False, "Ligand completely buried (not at interface)"
            
        return True, "Ligand properly at interface"
    
    elif user_goal == "symmetric_dimer":
        # Check both chains exist
        chains = get_chains(design)
        if len(chains) < 2:
            return False, "Not a dimer - only one chain generated"
        
        # Check symmetry
        rmsd = calculate_symmetry_rmsd(design, symmetry="C2")
        if rmsd > 1.5:
            return False, f"C2 symmetry broken (RMSD={rmsd})"
            
        return True, "Valid symmetric dimer"
```

### Red Flags That Indicate Wrong Implementation

| User Said | Output Shows | Problem |
|-----------|--------------|---------|
| "dimer" | Single chain | Forgot `symmetry` field |
| "ligand at interface" | Ligand buried inside one chain | Missing `select_partially_buried` on ligand |
| "symmetric enzyme" | Asymmetric structure | `inference_sampler.kind` not set to "symmetry" |
| "bind DNA" | No DNA in output | DNA not included in input |

---

## Ligand Handling Deep Dive

### Specifying Ligands

Ligands are specified by **CCD name** (Chemical Component Dictionary) or **index**:

```json
{
  "my_design": {
    "input": "structure_with_ligand.pdb",
    "ligand": "HAX"           // By CCD 3-letter code
    // OR
    "ligand": "HAX,OAA"       // Multiple ligands
    // OR  
    "ligand": "L1"            // By residue index in input
  }
}
```

### Custom Ligands (Not in CCD)

For ligands not in the Chemical Component Dictionary:

1. **Use unique residue name** with `L:` prefix to avoid CCD clashes
2. **Provide complete atom definitions** in the input PDB/CIF

```python
# When preparing input structure with custom ligand:
# Rename residue to avoid CCD lookup
ligand_resname = "L:AZO"  # L: prefix = custom, no CCD lookup

# Or provide full CIF with atom definitions
```

### RASA Control for Ligands

Control how buried/exposed ligand atoms are:

```json
{
  "ligand_design": {
    "input": "with_ligand.pdb",
    "contig": "60-80",
    "ligand": "AZO",
    "select_partially_buried": {"AZO": "N1,N2"},  // These atoms at interface
    "select_buried": {"AZO": "C3,C4,C5,C6"}       // These atoms buried
  }
}
```

**RASA options:**
- `select_buried`: Atom SASA ≈ 0 (inside protein)
- `select_partially_buried`: Atom at interface (some exposure)
- `select_exposed`: Atom on surface

### Ligand Orientation for Symmetric Designs

For C2 symmetric designs with ligand at interface:

```python
def orient_ligand_for_c2(ligand_mol, axis_atoms=["N1", "N2"]):
    """
    Orient ligand so the C2 axis passes through specified atoms.
    
    For azobenzene: N=N bond should be along the C2 axis.
    Each phenyl ring then contacts one monomer of the dimer.
    
           Chain A
              \
        Ph────N═N────Ph    ← N=N along C2 axis (z-axis)
              /
           Chain B
    """
    # Get axis atom positions
    n1_pos = ligand_mol.GetAtomPosition("N1")
    n2_pos = ligand_mol.GetAtomPosition("N2")
    
    # Align N1-N2 vector with z-axis
    current_axis = n2_pos - n1_pos
    target_axis = np.array([0, 0, 1])
    
    rotation = calculate_rotation_matrix(current_axis, target_axis)
    ligand_mol.Transform(rotation)
    
    # Center on midpoint of N=N bond (this becomes ori_token origin)
    center = (n1_pos + n2_pos) / 2
    ligand_mol.Translate(-center)
    
    return ligand_mol
```

---

## Common Pitfalls and Solutions

### Pitfall 1: Implementing Features in Isolation

```python
# ❌ BAD: Separate functions that can't compose
def design_symmetric_protein(symmetry):
    return rfd3_run(symmetry=symmetry)

def design_ligand_binder(ligand):
    return rfd3_run(ligand=ligand)

# User asks for "symmetric ligand binder" → which function to call??
```

```python
# ✅ GOOD: Single function with composable parameters
def design_protein(
    input_pdb: str,
    contig: str,
    symmetry: Optional[dict] = None,
    ligand: Optional[str] = None,
    select_partially_buried: Optional[dict] = None,
    select_hotspots: Optional[str] = None,
    # ... all features as optional parameters
):
    config = {
        "input": input_pdb,
        "contig": contig
    }
    if symmetry:
        config["symmetry"] = symmetry
    if ligand:
        config["ligand"] = ligand
    if select_partially_buried:
        config["select_partially_buried"] = select_partially_buried
    # ... combine all specified features
    return rfd3_run(config)
```

### Pitfall 2: Not Validating Against User Intent

```python
# ❌ BAD: Only check if RFD3 succeeded
if rfd3_output.exists():
    print("Success!")  # But is it actually what user wanted?

# ✅ GOOD: Validate against user's goal
design = load_structure(rfd3_output)
goal_met, reason = validate_design_meets_goal(design, user_goal="ligand_at_interface")
if not goal_met:
    print(f"Design failed user requirements: {reason}")
    # Try different parameters or approach
```

### Pitfall 3: Wrong Contig Syntax

```python
# ❌ WRONG: Using chain labels for NEW residues
contig = "A50-70"  # This means: take residues 50-70 FROM INPUT chain A

# ✅ CORRECT: No chain label for designed residues
contig = "50-70"   # This means: DESIGN 50-70 new residues

# ❌ WRONG: Wrong chain break syntax
contig = "50-70/0/A1-100"  # Wrong separator

# ✅ CORRECT: Use \0 for chain breaks
contig = "50-70,\\0,A1-100"  # Design 50-70, then chain break, then input A1-100
```

### Pitfall 4: Forgetting to Unfix Ligand Atoms for Diffusion

```json
// ❌ WRONG: Ligand atoms stay fixed (rigid docking)
{
  "design": {
    "input": "with_ligand.pdb",
    "ligand": "AZO"
    // ligand atoms fixed by default!
  }
}

// ✅ CORRECT: Explicitly unfix ligand for co-diffusion
{
  "design": {
    "input": "with_ligand.pdb",
    "ligand": "AZO",
    "select_fixed_atoms": {"AZO": ""}  // Empty string = unfix all ligand atoms
  }
}
```

---

## Pre-Implementation Checklist

Before writing code for any RFD3 workflow:

- [ ] **Parsed user goal**: What is the biological/functional outcome they want?
- [ ] **Identified all required features**: Symmetry? Ligand? Motif? RASA?
- [ ] **Checked feature compatibility**: Can these features combine?
- [ ] **Designed validation criteria**: How will I verify output matches goal?
- [ ] **Correct contig syntax**: Chain labels for input, no labels for design
- [ ] **Ligand handling**: CCD name correct? Need to unfix atoms for diffusion?
- [ ] **Symmetry setup**: `inference_sampler.kind: symmetry` in CLI args?
- [ ] **Origin placement**: Using `ori_token` or `infer_ori_strategy` correctly?

---

## Quick Reference: RFD3 JSON Schema (OFFICIAL)

```json
{
  "design_name": {
    // INPUT STRUCTURE
    "input": "path/to/structure.pdb",     // PDB or CIF
    
    // WHAT TO DESIGN
    "contig": "50-70,\\0,A1-100",          // Design spec + input residues
    "length": "80-120",                    // Total length constraint
    
    // LIGAND
    "ligand": "HAX,OAA",                   // By CCD name or index
    
    // MOTIF SCAFFOLDING
    "unindex": "A42,A87,A123",             // Unindexed motif residues
    
    // ATOM CONTROL
    "select_fixed_atoms": {                // What coordinates to fix
      "A1-10": "BKBN",
      "AZO": ""                            // Empty = unfix for diffusion
    },
    "select_unfixed_sequence": "A20-35",   // Where sequence can change
    
    // RASA (burial) CONTROL
    "select_buried": {"AZO": "C1,C2"},
    "select_partially_buried": {"AZO": "N1,N2"},
    "select_exposed": "A42",
    
    // H-BONDS
    "select_hbond_donor": {"A42": "N"},
    "select_hbond_acceptor": {"AZO": "N1"},
    
    // INTERFACE DESIGN
    "select_hotspots": "A42,A67",
    
    // SYMMETRY
    "symmetry": {"type": "C2"},            // See symmetry docs
    
    // POSITIONING
    "ori_token": [0.0, 0.0, 0.0],          // Custom origin
    "infer_ori_strategy": "hotspots",      // "com" or "hotspots"
    
    // PARTIAL DIFFUSION
    "partial_t": 15.0,                     // Noise in Å (5-15 recommended)
    
    // OPTIONS
    "is_non_loopy": true,                  // Fewer loops
    "plddt_enhanced": true                 // Better confidence
  }
}
```

CLI for symmetric designs:
```bash
rfd3 design \
  out_dir=outputs \
  inputs=my_design.json \
  inference_sampler.kind=symmetry \
  n_batches=10 \
  diffusion_batch_size=8
```

---

## Summary: The Mental Model

```
┌────────────────────────────────────────────────────────────────────┐
│                     RFD3 IS A FEATURE COMPOSER                     │
│                                                                    │
│  Think of RFD3 like Photoshop layers, not separate tools:         │
│                                                                    │
│    Layer 1: Base structure (contig)          ← Always needed       │
│    Layer 2: Input structure (input)          ← If using existing   │
│    Layer 3: Symmetry (symmetry)              ← Adds oligomer       │
│    Layer 4: Ligand (ligand)                  ← Adds binding site   │
│    Layer 5: Motif (unindex)                  ← Adds active site    │
│    Layer 6: RASA (select_*_buried)           ← Controls burial     │
│    Layer 7: H-bonds (select_hbond_*)         ← Fine control        │
│    Layer 8: Origin (ori_token)               ← Positioning         │
│                                                                    │
│  Layers STACK. Don't implement them as separate workflows.        │
│                                                                    │
└────────────────────────────────────────────────────────────────────┘
```

**When user asks for complex design:**
1. Identify ALL layers needed
2. Compose them in ONE RFD3 config JSON
3. Use correct CLI args (especially `inference_sampler.kind=symmetry`)
4. Validate output has ALL requested features
5. If validation fails, adjust parameters (don't remove features)

---

## Appendix: InputSelection Atom Shortcuts

| Shortcut | Atoms Selected | Use Case |
|----------|----------------|----------|
| `ALL` | All atoms in residue | Fix everything |
| `BKBN` | N, CA, C, O | Backbone only |
| `TIP` | Functional tip atoms | Catalytic residues |
| `""` (empty) | Nothing | UNFIX atoms (for diffusion) |
| `"N,CA,C,O,CB"` | Listed atoms | Custom selection |

---

## Appendix: Symmetry Types

| Type | Description | Example Use |
|------|-------------|-------------|
| `C2` | 2-fold cyclic | Dimers |
| `C3` | 3-fold cyclic | Trimers |
| `C4` | 4-fold cyclic | Tetramers |
| `D2` | 2-fold dihedral | Tetramers (222) |
| `D3` | 3-fold dihedral | Hexamers (32) |

Specify in config:
```json
"symmetry": {"type": "C2"}
```

And use CLI flag:
```bash
inference_sampler.kind=symmetry
```

---

## Appendix: FAQ from Official Docs

**Q: Can I guide on secondary structure?**
A: Not directly, but use `is_non_loopy: true` for fewer loops (more helices).

**Q: Do I need `select_fixed_atoms` & `select_unfixed_sequence` every time?**
A: No, defaults apply when input is present.

**Q: Why "Input provided but unused" warning?**
A: You gave an input PDB/CIF but no `contig`, `unindex`, `ligand`, or `partial_t` to use it.

**Q: What do logged b-factors mean?**
A: Sequence head confidence. Use `spectrum b` in PyMOL to visualize. High entropy = uncertain sequence assignment.
