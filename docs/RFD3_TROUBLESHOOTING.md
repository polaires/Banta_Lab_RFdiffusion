# RFD3 Pre-Implementation Checklist & Troubleshooting

## STOP! Before Writing Code, Complete This Checklist

### Step 1: Parse the User's TRUE Goal

Don't implement the literal words—understand the biological intent.

| User Says | Literal Interpretation | TRUE Goal |
|-----------|----------------------|-----------|
| "ligand binder" | Protein that touches ligand | Protein with functional binding pocket |
| "dimerization binder" | Two proteins | TWO-CHAIN protein with specific interface |
| "azobenzene at interface" | Azobenzene somewhere | Azobenzene BETWEEN two chains, contacting BOTH |
| "photoswitchable" | Has azobenzene | Binding changes with cis/trans isomerization |

**Write down the user's goal in one sentence:**
```
"User wants a [oligomeric state] protein where [ligand/target] is [position] 
and the function is [what it does]."
```

Example for our case:
```
"User wants a C2 DIMER protein where AZOBENZENE is AT THE INTERFACE 
BETWEEN THE TWO CHAINS and the function is DIMERIZATION CONTROLLED BY LIGAND."
```

---

### Step 2: Identify Required RFD3 Features

Check ALL that apply to user's goal:

```
□ Unconditional generation (just make a fold)
□ Binder design (bind to existing protein target)
■ Symmetric oligomer (C2, C3, D2, etc.)  ← NEEDED for "dimer"
■ Ligand conditioning                      ← NEEDED for "azobenzene"
□ Nucleic acid binding
□ Enzyme/motif scaffolding (unindex)
■ RASA control (burial depth)              ← NEEDED for "at interface"
□ H-bond specification
■ Origin positioning                        ← NEEDED to center on ligand
```

**If multiple boxes checked → Features must be COMBINED in single config!**

---

### Step 3: Write the RFD3 Config BEFORE Implementation

```json
// Fill this out BEFORE writing Python code:
{
  "design_name": {
    "input": "____",              // Path to input structure
    "contig": "____",             // Design spec
    "symmetry": {},               // From Step 2 if needed
    "ligand": "____",             // CCD name if needed
    "select_partially_buried": {},// RASA if needed
    "ori_token": []               // Origin if needed
  }
}
```

For our azobenzene dimer case:
```json
{
  "azobenzene_c2_dimer": {
    "input": "azobenzene_at_origin.pdb",
    "contig": "50-70",
    "symmetry": {"type": "C2"},
    "ligand": "AZO",
    "select_partially_buried": {
      "AZO": "N1,N2"
    },
    "ori_token": [0.0, 0.0, 0.0]
  }
}
```

---

### Step 4: Define Validation Criteria BEFORE Running

**What MUST be true for the output to be "correct"?**

```python
# Write these checks BEFORE implementation:

def output_is_correct(design):
    checks = {
        "has_two_chains": len(get_chains(design)) >= 2,
        "ligand_contacts_chain_A": count_contacts("AZO", "A") >= 3,
        "ligand_contacts_chain_B": count_contacts("AZO", "B") >= 3,
        "ligand_at_interface": ligand_sasa > 5,  # Not completely buried
        "c2_symmetry_maintained": symmetry_rmsd < 1.5,
    }
    
    failed = [k for k, v in checks.items() if not v]
    return len(failed) == 0, failed
```

**If you can't write validation criteria, you don't understand the goal yet.**

---

## Troubleshooting: When Output Doesn't Match Goal

### Problem: "Got a monomer instead of dimer"

**Symptom:** Output has only chain A, no chain B

**Cause:** Forgot `symmetry` field OR forgot CLI flag

**Fix:**
```json
// WRONG - no symmetry
{"contig": "50-70", "ligand": "AZO"}

// RIGHT - add symmetry field
{"contig": "50-70", "ligand": "AZO", "symmetry": {"type": "C2"}}
```

AND use correct CLI:
```bash
rfd3 design out_dir=out inputs=config.json inference_sampler.kind=symmetry
```

---

### Problem: "Ligand is buried inside one chain, not at interface"

**Symptom:** Ligand contacts only chain A, SASA ≈ 0

**Cause:** Didn't use `select_partially_buried` and/or wrong origin

**Fix:**
```json
// WRONG - ligand gets fully buried
{
  "input": "with_ligand.pdb",
  "contig": "50-70",
  "ligand": "AZO",
  "symmetry": {"type": "C2"}
}

// RIGHT - control burial and position
{
  "input": "with_ligand_at_origin.pdb",
  "contig": "50-70",
  "ligand": "AZO",
  "symmetry": {"type": "C2"},
  "select_partially_buried": {"AZO": "N1,N2"},  // Keep at interface
  "ori_token": [0.0, 0.0, 0.0]                  // Center protein on ligand
}
```

---

### Problem: "Input provided but unused" warning

**Symptom:** RFD3 logs warning about unused input

**Cause:** You provided `input` but no `contig`, `unindex`, `ligand`, or `partial_t`

**Fix:** You must specify what to DO with the input:
```json
// WRONG - input given but nothing references it
{
  "input": "structure.pdb"
}

// RIGHT - specify what to use from input
{
  "input": "structure.pdb",
  "contig": "A1-100",           // Use residues from input
  // OR
  "ligand": "HAX",              // Use ligand from input
  // OR
  "partial_t": 10.0             // Partial diffusion on input
}
```

---

### Problem: "Atom X not found" or CCD lookup error

**Symptom:** RFD3 crashes looking for atoms or CCD entry

**Cause:** Ligand residue name conflicts with CCD, or custom ligand not properly defined

**Fix:**
```python
# WRONG - "LIG" or "AZO" might exist in CCD with different atoms
residue_name = "LIG"

# RIGHT - prefix with L: to skip CCD lookup
residue_name = "L:AZO"

# OR - ensure your input PDB/CIF has complete atom definitions
# that match what you reference in config
```

---

### Problem: "Symmetry is broken in output"

**Symptom:** Chain A and B aren't symmetric, high RMSD

**Cause:** Ligand not aligned with symmetry axis, or wrong CLI args

**Fix:**
```python
def orient_for_c2(ligand):
    """
    For C2 symmetry, the ligand's C2 axis must align with 
    the protein's C2 axis.
    
    For azobenzene: N=N bond IS the C2 axis.
    """
    # Align N=N bond with z-axis
    n1, n2 = get_atom_positions(ligand, ["N1", "N2"])
    align_vector_to_axis(ligand, n2 - n1, target=[0, 0, 1])
    
    # Center on N=N midpoint (place at origin)
    center_on_point(ligand, (n1 + n2) / 2)
    
    return ligand
```

And ensure CLI uses:
```bash
inference_sampler.kind=symmetry
```

---

### Problem: "Ligand atoms are fixed/rigid, not co-diffusing"

**Symptom:** Ligand stays exactly in input conformation

**Cause:** Ligand atoms fixed by default; need to explicitly unfix

**Fix:**
```json
// WRONG - ligand atoms stay fixed
{
  "input": "with_ligand.pdb",
  "ligand": "AZO"
}

// RIGHT - unfix ligand atoms for co-diffusion
{
  "input": "with_ligand.pdb",
  "ligand": "AZO",
  "select_fixed_atoms": {"AZO": ""}  // Empty string = unfix all
}
```

---

### Problem: "Design looks good structurally but won't work biologically"

**Symptom:** High pLDDT, good RMSD, but nonsensical design

**Cause:** Validated structure, not function

**Fix:** Add biological validation:
```python
def validate_biology(design):
    # 1. Binding site is accessible
    pocket_volume = calculate_pocket_volume(design)
    if pocket_volume < ligand_volume:
        return False, "Pocket too small for ligand"
    
    # 2. Interface has reasonable chemistry
    interface_residues = get_interface_residues(design)
    hydrophobic_fraction = count_hydrophobic(interface_residues) / len(interface_residues)
    if hydrophobic_fraction < 0.3:
        return False, "Interface too polar for hydrophobic ligand"
    
    # 3. No clashes
    clashes = detect_clashes(design)
    if clashes > 0:
        return False, f"{clashes} atomic clashes detected"
    
    return True, "Biologically reasonable"
```

---

## Quick Debugging Commands

```bash
# Check what chains exist in output
grep "^ATOM" design.pdb | cut -c22 | sort -u

# Count ligand contacts per chain (atoms within 4Å)
pymol -c -d "load design.pdb; select contacts_A, chain A within 4 of resn AZO; print(len(contacts_A))"

# Check symmetry RMSD
pymol -c -d "load design.pdb; align chain A, chain B; print(cmd.get_rmsd())"

# Visualize ligand position
pymol design.pdb -d "show surface; color white, chain A; color gray, chain B; color red, resn AZO"

# Check if ligand is at interface (partial SASA)
pymol -c -d "load design.pdb; get_area resn AZO"
```

---

## RFD3 CLI Quick Reference

```bash
# Basic run
rfd3 design out_dir=outputs inputs=config.json

# With symmetry
rfd3 design out_dir=outputs inputs=config.json inference_sampler.kind=symmetry

# More designs
rfd3 design out_dir=outputs inputs=config.json n_batches=10 diffusion_batch_size=8

# Debugging (keep guideposts, trajectories)
rfd3 design out_dir=outputs inputs=config.json \
  cleanup_guideposts=False \
  dump_trajectories=True

# Quick unconditional test
rfd3 design out_dir=test inputs=null specification.length=100
```

---

## The Golden Rule

```
╔══════════════════════════════════════════════════════════════════╗
║                                                                  ║
║   If the user asks for X AND Y, don't implement X OR Y.          ║
║   Implement X AND Y together.                                    ║
║                                                                  ║
║   RFD3 features COMPOSE. Use them together.                      ║
║                                                                  ║
╚══════════════════════════════════════════════════════════════════╝
```
