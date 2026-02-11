# Dy-TriNOx Protein Design Campaign

## Overview

Design proteins to coordinate the Dy-TriNOx complex at its single open apical coordination site.

## The Dy-TriNOx Complex

```
            Open Site (Z+)
                 ↑
                 |
    N2 (apical) ← Dy3+ → (needs 1-2 donors here)
           /    |    \
         O1    O2    O3  (phenolate oxygens)
          |     |     |
         N1    N3    N4  (amine nitrogens)
          |     |     |
        Arm1  Arm2  Arm3 (aromatic + tert-butyl)
```

**TriNOx provides 7 donors:**
- 3× phenolate O (O1, O2, O3)
- 3× amine N (N1, N3, N4)
- 1× apical N (N2)

**Protein must provide:** 1-2 O/N donors at position (-0.007, 7.835, 5.494)

## Why Conformers Are NOT Critical Here

Unlike flexible ligands (e.g., azobenzene photoswitches), TriNOx is a **rigid tripodal scaffold**:

1. **Aromatic framework locked**: The three phenyl arms are geometrically constrained
2. **Metal position fixed**: Dy sits in a pre-defined coordination cage
3. **Open site direction defined**: The +Z direction is the only accessible coordination vector

**Where we DO need diversity:**
- Protein backbone conformations (via multiple RFD3 runs)
- Side chain rotamers at metal site (via LigandMPNN sampling)

## Workflow

### Step 1: Prepare Input (DONE)
- `inputs/Dy_TriNOx_clean.pdb` - Complex without dummy X atom
- Chain X = Dy (metal)
- Chain L = TriNOx ligand

### Step 2: Run RFD3 Designs

```bash
# Via Foundry API
python run_trinox_design.py --config configs/trinox_scaffold_config.json

# Or direct RFD3 call:
rfdiffusion3 \
  --input_pdb inputs/Dy_TriNOx_clean.pdb \
  --contig "80-100" \
  --select_fixed_atoms "X:all,L:all" \
  --select_buried "X:all" \
  --use_classifier_free_guidance True \
  --cfg_scale 2.0 \
  --num_designs 20
```

### Step 3: Sequence Design with LigandMPNN

```bash
ligandmpnn \
  --pdb_path outputs/r1a_small_000.pdb \
  --ligand_file inputs/Dy_TriNOx_clean.pdb \
  --temperature 0.1 \
  --num_seq 8 \
  --bias_AA "E:2.0,D:1.5,H:1.5" \
  --bias_by_distance_from_ligand
```

### Step 4: Validate with AlphaFold3

Include in AF3 job:
- Protein sequence
- Dy3+ ion
- TriNOx ligand (as custom SDF or CCD)

## Key Differences from Citrate Campaign

| Aspect | Citrate-Ln | Dy-TriNOx |
|--------|------------|-----------|
| Coordination needed | 3-4 donors | **1-2 donors** |
| Ligand conformers | Might help | **Not needed** |
| Protein conformers | Via seeds | **Via seeds** |
| Pocket size | Small | **Large (for TriNOx)** |
| Steric challenge | Low | **High (tert-butyl)** |

## Expected Outcomes

**Success criteria:**
- [ ] Dy-O/N distance: 2.2-2.7 Å
- [ ] Metal burial: >70% SASA reduction
- [ ] Protein pLDDT: >0.80
- [ ] Complete 8-9 coordination sphere

## Files

```
Dy_TriNOx_scaffold/
├── README.md                          # This file
├── inputs/
│   └── Dy_TriNOx_clean.pdb           # Cleaned complex (no dummy atom)
├── configs/
│   └── trinox_scaffold_config.json   # RFD3 configurations
├── outputs/                          # RFD3 outputs (to be created)
└── analysis/                         # Analysis results (to be created)
```
