# RFdiffusion3 Comprehensive Skill Document

> **Purpose**: Complete reference for RFdiffusion3 protein design via Foundry API based on official documentation. Use this skill when working on any RFD3-related task.

---

## 1. Core Concepts

### What is RFdiffusion3?
RFdiffusion3 (RFD3) is a diffusion-based generative model for protein structure design, accessed through the Rosetta Foundry API. It generates novel protein backbones conditioned on various constraints.

### Key Terminology
| Term | Definition |
|------|------------|
| **Contig** | Specification of chain lengths and fixed residues (e.g., `"60-80,/0,60-80"`) |
| **Diffusion** | Iterative denoising process that generates structures |
| **Hotspot** | Key interface residues to preserve during design |
| **RASA** | Relative Accessible Surface Area - controls burial/exposure |
| **Conditioning** | Constraints applied during diffusion (fixed atoms, H-bonds, etc.) |

### API vs CLI Differences
**CRITICAL**: RFdiffusion3 via Foundry API differs from RFdiffusion v1/v2 CLI:

| Feature | RFD3/Foundry | RFdiffusion v1/v2 CLI |
|---------|--------------|----------------------|
| Guiding potentials | ❌ NOT SUPPORTED | ✅ substrate_contacts, olig_contacts, monomer_ROG |
| Symmetry | ✅ Built-in | ✅ Built-in |
| Partial diffusion | ✅ Yes | ✅ Yes |
| H-bond conditioning | ✅ Yes | ❌ No |
| RASA conditioning | ✅ Yes | ❌ No |

---

## 2. Input Specification Parameters

### Contig Syntax
```
Format: "length_range,/chain_break,length_range" or "ChainResidue,length,ChainResidue"

Examples:
- "60-80"                    # Single chain, 60-80 residues
- "60-80,/0,60-80"          # Two chains, 60-80 residues each
- "A10-50,/0,40-60"         # Keep A10-50, add new 40-60 residue chain
- "5-12,A10,3-8,A15,15-35"  # Scaffold around fixed residues A10, A15
```

**Syntax Rules:**
- Comma (`,`) = element separator
- Forward slash `/0` = chain break with zero gap
- Hyphen (`-`) = length range (e.g., `60-80`)
- `ChainResnum` = keep specific residue (e.g., `A10`)

### Fixed Atoms (select_fixed_atoms)
```python
# Fix specific atoms during diffusion
select_fixed_atoms = {
    "A10": "CA,C,N,O",      # Backbone atoms
    "A15": "CA,C,N,O,CB",   # Backbone + CB
    "L1": "NI"              # Metal atom (e.g., Nickel)
}
```

### Hotspots (select_hotspots)
```python
# Preserve key interface residues
select_hotspots = {
    "B": "45,48,52,89,93"   # Residues on chain B to maintain
}
```

### RASA Conditioning (Burial/Exposure)
```python
# Control surface accessibility
select_buried = {"A": "10,15,20"}    # Force burial
select_exposed = {"A": "1,50,100"}   # Force exposure

# Ranges supported
select_buried = {"A": "10-20"}       # Bury residues 10-20
```

### H-Bond Conditioning
```python
# Specify H-bond donors/acceptors
select_hbond_donor = {"A10": "N"}           # Backbone NH as donor
select_hbond_acceptor = {"UNL": "O1,O2"}    # Ligand oxygens as acceptors
```

**H-Bond Success Rate**: ~37% when properly specified (per official docs)

### Symmetry Options
```python
# Point group symmetry
symmetry = "c2"   # C2 dimer
symmetry = "c3"   # C3 trimer
symmetry = "d2"   # D2 tetramer

# Symmetry constraints:
# - diffusion_batch_size MUST be 1
# - For T, O, I: use low_memory_mode=True
```

**Available Symmetry Types:**
| Type | Description | Subunits |
|------|-------------|----------|
| C2-Cn | Cyclic | n copies around axis |
| D2-Dn | Dihedral | 2n copies |
| T | Tetrahedral | 12 copies |
| O | Octahedral | 24 copies |
| I | Icosahedral | 60 copies |

**Important Symmetry Notes** (from official docs):
- Contigs must specify something precisely symmetric
- Motifs should be pre-symmetrized to RFdiffusion's canonical axes
- Single monomer is generated, then n-1 copies propagated automatically

### Partial Diffusion
```python
# Start from existing structure, add noise, then denoise
partial_T = 20              # Noise timesteps (0-50 typical)
partial_diffusion = True    # Enable partial diffusion mode
```

### Inference Parameters
```python
num_designs = 10            # Number of designs to generate
diffusion_steps = 50        # Denoising steps (default: 50)
diffusion_batch_size = 1    # Batch size (must be 1 for symmetry)
seed = 42                   # Random seed for reproducibility
```

---

## 3. Design Types & Official Examples

### 3.1 Small Molecule Binder Design
**Goal**: Design protein to bind small molecule ligand.

From official `sm_binder_design.json`:
```json
{
  "contig": "A1-74/0 70-100",
  "select_fixed_atoms": {
    "A44": "all",
    "A56": "all",
    "A60": "all"
  },
  "select_hotspots": {
    "A": "44,56,60"
  }
}
```

**Key Points:**
- Fix atoms on residues that interact with ligand
- Use hotspots to preserve binding interface
- Contig specifies existing binding domain + new scaffold

### 3.2 Protein Binder Design
**Goal**: Design protein that binds target protein.

From official `protein_binder_design.json`:
```json
{
  "contig": "B1-100/0 70-100",
  "select_hotspots": {
    "B": "25,29,33"
  }
}
```

**Key Points:**
- Target protein on chain B
- Hotspots identify key interface residues on target
- New binder designed as 70-100 residue chain

### 3.3 Enzyme Design (Active Site Scaffolding)
**Goal**: Design scaffold around catalytic residues.

From official `enzyme_design.json`:
```json
{
  "contig": "5-20/A10/3-10/A25/3-10/A40/10-30",
  "select_fixed_atoms": {
    "A10": "all",
    "A25": "all",
    "A40": "all"
  }
}
```

**Key Points:**
- Theozyme residues (A10, A25, A40) are fixed
- Variable-length loops connect fixed residues
- "Unindexed" residues in contig allow flexibility

### 3.4 Symmetric Oligomer Design
**Goal**: Design symmetric assembly.

From official `symmetry.md`:
```python
# C3 trimer example
{
    "contig": "100",
    "symmetry": "c3",
    "diffusion_batch_size": 1
}
```

**Key Points:**
- Single monomer contig (symmetry handles copies)
- `diffusion_batch_size` must be 1 for symmetry
- For T/O/I symmetry, use `low_memory_mode=True`

### 3.5 Nucleic Acid Binder Design
From official docs - similar to protein binder:
```json
{
  "contig": "A1-20/0 60-80",
  "select_hotspots": {
    "A": "5,8,12,15"
  }
}
```

---

## 4. Guiding Potentials - NOT SUPPORTED

**CRITICAL**: The following potentials from RFdiffusion v1/v2 CLI are **NOT available** in RFD3/Foundry API:

| Potential | Purpose | Alternative in RFD3 |
|-----------|---------|---------------------|
| `substrate_contacts` | Metal/ligand proximity | Use `select_fixed_atoms` |
| `olig_contacts` | Interface optimization | Use `symmetry` constraints |
| `monomer_ROG` | Compactness | Rely on contig length constraints |

If your workflow requires these potentials, you must use the RFdiffusion v1/v2 CLI directly, not the Foundry API.

---

## 5. PPI Design Tutorial (Official)

From `ppi_design_tutorial.md`:

### Step 1: Prepare Target
- Obtain target PDB structure
- Identify interface residues (hotspots)
- Clean PDB (remove waters, select relevant chains)

### Step 2: Design Binder
```python
{
    "pdb_content": target_pdb,
    "contig": "B1-200/0 60-100",  # Target chain B, design 60-100 binder
    "select_hotspots": {
        "B": "45,48,52,89,93"
    },
    "num_designs": 100  # Generate many candidates
}
```

### Step 3: Sequence Design
- Use ProteinMPNN or LigandMPNN for sequence optimization
- Generate multiple sequences per backbone

### Step 4: Validation
- Structure prediction (ESMFold/AlphaFold2)
- Binding affinity prediction
- Experimental validation

---

## 6. Inference Configuration

From official `rfdiffusion3.yaml`:

```yaml
rfdiffusion3:
  engine: rfdiffusion3
  version: 1.0.0

  default_params:
    diffusion_steps: 50
    diffusion_batch_size: 1

  symmetry_params:
    low_memory_mode: false  # Set true for T/O/I

  conditioning:
    hbond_success_rate: 0.37  # Expected H-bond satisfaction
```

---

## 7. Best Practices

### Contig Design
1. **Be specific about lengths** - Use ranges that make biological sense
2. **Chain breaks** - Use `/0` for separate chains
3. **Fixed residues** - Include in contig with chain+resnum format

### Hotspot Selection
1. **Interface residues** - Choose residues that make key contacts
2. **Energetically important** - Prioritize buried, H-bonding residues
3. **Typically 3-10 hotspots** - More isn't always better

### Symmetry Usage
1. **Always set `diffusion_batch_size=1`** when using symmetry
2. **Pre-symmetrize motifs** to canonical axes
3. **Use `low_memory_mode=True`** for T/O/I symmetry

### H-Bond Conditioning
1. **Expect ~37% success rate** - not all H-bonds will form
2. **Specify both donor and acceptor** when possible
3. **Use for critical interactions** (catalytic, binding)

---

## 8. Quick Reference Card

```
CONTIG SYNTAX:
  "60-80"              → Single chain, 60-80 residues
  "60-80/0 60-80"      → Two chains (official notation)
  "60-80,/0,60-80"     → Two chains (alternative notation)
  "A10-50/0 40-60"     → Keep A10-50, add new chain
  "5-12/A10/3-8/A15"   → Scaffold around fixed residues

SYMMETRY:
  c2, c3, c4, ...      → Cyclic
  d2, d3, d4, ...      → Dihedral
  t, o, i              → Platonic (use low_memory_mode)

KEY PARAMETERS:
  select_fixed_atoms   → Fix backbone/sidechain atoms
  select_hotspots      → Preserve interface residues
  select_buried        → Force burial (low RASA)
  select_exposed       → Force exposure (high RASA)
  select_hbond_*       → H-bond conditioning

REQUIRED FOR SYMMETRY:
  diffusion_batch_size = 1

NOT SUPPORTED IN FOUNDRY API:
  guiding_potentials (substrate_contacts, olig_contacts, monomer_ROG)
```

---

## 9. Troubleshooting

### Common Issues

| Issue | Cause | Solution |
|-------|-------|----------|
| "Extra inputs not permitted" | Unsupported parameter | Remove guiding_potentials |
| Asymmetric output | Missing symmetry flag | Add `symmetry="c2"` etc. |
| Memory error with T/O/I | Large assembly | Use `low_memory_mode=True` |
| H-bonds not forming | Expected behavior | ~37% success rate is normal |
| Motif not fixed | Wrong atom selection | Verify `select_fixed_atoms` syntax |

### Debugging Tips
1. **Check contig parsing** - Verify chain breaks interpreted correctly
2. **Validate PDB input** - Ensure clean PDB with proper chain IDs
3. **Start simple** - Test with minimal constraints first
4. **Generate many designs** - Use `num_designs=100` for diversity

---

## 10. Official Documentation Sources

| Resource | URL |
|----------|-----|
| Foundry API | https://github.com/RosettaCommons/foundry |
| RFdiffusion3 Docs | foundry/docs/rfdiffusion3/ |
| Input Reference | foundry/docs/rfdiffusion3/input.md |
| Symmetry Guide | foundry/docs/rfdiffusion3/symmetry.md |
| Design Examples | foundry/docs/rfdiffusion3/examples/ |

---

*Document generated: 2026-01-15*
*Source: Official RFdiffusion3/Foundry documentation from RosettaCommons*
