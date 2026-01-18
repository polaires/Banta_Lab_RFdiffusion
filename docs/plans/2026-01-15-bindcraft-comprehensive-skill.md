# BindCraft Comprehensive Skill Document

> **Purpose**: Complete reference for BindCraft protein binder design pipeline based on official documentation from martinpacesa/BindCraft. Use this skill when designing de novo protein binders or understanding the AF2 backpropagation + MPNN + PyRosetta workflow.

---

## 1. Core Concepts

### What is BindCraft?
BindCraft is a computational pipeline for de novo protein binder design that combines:
- **AlphaFold2 backpropagation**: Optimizes sequences by gradient descent through AF2
- **ProteinMPNN**: Generates diverse sequences compatible with designed backbones
- **PyRosetta**: Structural relaxation and interface analysis

### Pipeline Overview
```
Target PDB → AF2 Hallucination → MPNN Sequence Design → AF2 Validation → PyRosetta Relaxation → Filtering
                    ↓                    ↓                   ↓                    ↓
              Design trajectory    Multiple seqs      Predict complex       Score interface
```

### Key Innovation
Unlike RFdiffusion (backbone-first), BindCraft uses **AF2 backpropagation** to simultaneously optimize sequence AND structure by computing gradients through AlphaFold2's neural network.

---

## 2. Design Algorithms

### Algorithm Options
| Algorithm | Stages | Method | Speed | Memory |
|-----------|--------|--------|-------|--------|
| **2stage** | Logits→PSSM | Semi-greedy sampling | Fastest | Low |
| **3stage** | Logits→Softmax→One-hot | Standard progression | Medium | Medium |
| **4stage** | Logits→Softmax→One-hot→PSSM | Most comprehensive (default) | Slower | Higher |
| **greedy** | Random mutations | Loss-reducing mutations | Slowest | Lowest |
| **mcmc** | Probabilistic | Temperature-based sampling | Variable | Low |

### 4-Stage Design (Default)
1. **Soft iterations** (75): All amino acids considered at all positions
2. **Temporary iterations** (45): Softmax-weighted amino acid probabilities
3. **Hard iterations** (5): Single amino acid per position (one-hot)
4. **Greedy iterations** (15): Random PSSM mutations reducing loss

---

## 3. Configuration Files

### Three Required JSON Files

#### 1. Target Settings (`settings_target/*.json`)
```json
{
    "design_path": "/path/to/output/",
    "binder_name": "MyBinder",
    "starting_pdb": "/path/to/target.pdb",
    "chains": "A",
    "target_hotspot_residues": "56,58,62",
    "lengths": [65, 150],
    "number_of_final_designs": 100
}
```

| Parameter | Description |
|-----------|-------------|
| `design_path` | Output directory for designs |
| `binder_name` | Prefix for output files |
| `starting_pdb` | Target protein PDB file |
| `chains` | Which chains to target |
| `target_hotspot_residues` | Specific residues to bind (null = auto-select) |
| `lengths` | [min, max] binder length range |
| `number_of_final_designs` | Target count of passing designs |

#### 2. Advanced Settings (`settings_advanced/*.json`)
See Section 4 for complete parameter reference.

#### 3. Filter Settings (`settings_filters/*.json`)
See Section 5 for complete filter reference.

---

## 4. Advanced Settings Reference

### Design Algorithm
```json
{
    "design_algorithm": "4stage",
    "soft_iterations": 75,
    "temporary_iterations": 45,
    "hard_iterations": 5,
    "greedy_iterations": 15,
    "greedy_percentage": 1
}
```

### AlphaFold2 Configuration
```json
{
    "use_multimer_design": true,
    "num_recycles_design": 1,
    "num_recycles_validation": 3,
    "sample_models": true,
    "rm_template_seq_design": false,
    "rm_template_sc_design": false,
    "predict_initial_guess": false,
    "predict_bigbang": false
}
```

| Parameter | Description |
|-----------|-------------|
| `use_multimer_design` | AF2-multimer vs AF2-ptm |
| `num_recycles_design` | AF2 recycles during design |
| `num_recycles_validation` | AF2 recycles for validation |
| `sample_models` | Randomize AF2 params (prevents overfitting) |
| `rm_template_seq_design` | Remove target template sequences |
| `rm_template_sc_design` | Remove target sidechains |
| `predict_initial_guess` | Bias with binder atom positions |
| `predict_bigbang` | Position bias for large complexes (>600 AA) |

### Loss Function Weights
```json
{
    "weights_plddt": 0.1,
    "weights_pae_intra": 0.4,
    "weights_pae_inter": 0.1,
    "weights_con_intra": 1.0,
    "weights_con_inter": 1.0,
    "weights_helicity": -0.3,
    "weights_iptm": 0.05,
    "weights_rg": 0.3,
    "weights_termini_loss": 0.1
}
```

| Weight | Description | Effect |
|--------|-------------|--------|
| `weights_plddt` | Confidence score | Higher = more confident |
| `weights_pae_intra` | Intra-chain alignment error | Higher = better internal packing |
| `weights_pae_inter` | Inter-chain alignment error | Higher = better interface |
| `weights_con_intra` | Internal contacts | Higher = more compact |
| `weights_con_inter` | Interface contacts | Higher = more interface contacts |
| `weights_helicity` | Helix propensity | Negative = favor β-sheets |
| `weights_iptm` | Interface pTM | Higher = better interface prediction |
| `weights_rg` | Radius of gyration | Higher = more compact |
| `weights_termini_loss` | N/C terminus distance | Higher = closer termini |

### Contact Parameters
```json
{
    "intra_contact_distance": 14.0,
    "inter_contact_distance": 20.0,
    "intra_contact_number": 2,
    "inter_contact_number": 2
}
```

### MPNN Settings
```json
{
    "enable_mpnn": true,
    "mpnn_fix_interface": true,
    "num_seqs": 20,
    "max_mpnn_sequences": 2,
    "sampling_temp": 0.1,
    "backbone_noise": 0.00,
    "mpnn_weights": "soluble",
    "model_path": "v_48_020",
    "save_mpnn_fasta": false
}
```

### Beta Sheet Optimization
```json
{
    "optimise_beta": true,
    "optimise_beta_extra_soft": 0,
    "optimise_beta_extra_temp": 0,
    "optimise_beta_recycles_design": 3,
    "optimise_beta_recycles_valid": 3
}
```

### Amino Acid Control
```json
{
    "omit_AAs": "C",
    "force_reject_AA": false
}
```

### Output Management
```json
{
    "save_design_animations": true,
    "save_design_trajectory_plots": true,
    "remove_unrelaxed_trajectory": true,
    "remove_unrelaxed_complex": true,
    "remove_binder_monomer": true,
    "zip_animations": true,
    "zip_plots": true,
    "save_trajectory_pickle": false
}
```

### Acceptance Control
```json
{
    "max_trajectories": null,
    "acceptance_rate": 0.01,
    "start_monitoring": 600,
    "enable_rejection_check": true
}
```

---

## 5. Filter Settings Reference

### Structural Confidence Filters
| Filter | Threshold | Direction | Description |
|--------|-----------|-----------|-------------|
| `pLDDT` | ≥0.8 | Higher better | Overall confidence |
| `pTM` | ≥0.55 | Higher better | Structure accuracy |
| `i_pTM` | ≥0.5 | Higher better | Interface confidence |
| `i_pAE` | ≤0.35 | Lower better | Interface alignment error |
| `Binder_pLDDT` | ≥0.8 | Higher better | Binder monomer confidence |

### Interface Quality Filters
| Filter | Threshold | Direction | Description |
|--------|-----------|-----------|-------------|
| `ShapeComplementarity` | ≥0.6 | Higher better | Surface fit |
| `dG` | ≤0 | Lower better | Binding energy (negative = favorable) |
| `dSASA` | ≥1 | Higher better | Buried surface area |
| `n_InterfaceResidues` | ≥7 | Higher better | Contact residue count |
| `n_InterfaceHbonds` | ≥3 | Higher better | Hydrogen bonds |
| `n_InterfaceUnsatHbonds` | ≤4 | Lower better | Unsatisfied H-bonds |
| `Surface_Hydrophobicity` | ≤0.35 | Lower better | Surface hydrophobicity |

### Structural Composition Filters
| Filter | Threshold | Direction | Description |
|--------|-----------|-----------|-------------|
| `Binder_Loop%` | ≤90% | Lower better | Loop content |
| `Binder_RMSD` | ≤3.5 | Lower better | Structural deviation |
| `HotspotRMSD` | ≤6.0 | Lower better | Binding site deviation |

### Amino Acid Composition Filters
| Filter | Threshold | Description |
|--------|-----------|-------------|
| `InterfaceAAs_K` | ≤3 | Max lysines at interface |
| `InterfaceAAs_M` | ≤3 | Max methionines at interface |

### Filter Prefix Convention
- **No prefix**: Average across all models
- **N_prefix** (e.g., `1_pLDDT`): Per-model statistics (models 1-5)

---

## 6. Available Presets

### Advanced Settings Presets
| Preset | Use Case |
|--------|----------|
| `default_4stage_multimer.json` | Standard protein binders (default) |
| `default_4stage_multimer_hardtarget.json` | Difficult targets |
| `default_4stage_multimer_flexible.json` | Flexible interface |
| `betasheet_4stage_multimer.json` | β-sheet rich binders |
| `peptide_3stage_multimer.json` | Short peptide binders |
| `*_mpnn.json` variants | With MPNN interface fixing |

### Filter Presets
| Preset | Use Case |
|--------|----------|
| `default_filters.json` | Standard stringency |
| `relaxed_filters.json` | Less stringent (more designs pass) |
| `peptide_filters.json` | Peptide-specific thresholds |
| `no_filters.json` | No filtering (debugging) |

---

## 7. Loss Functions

### Standard Loss Functions
1. **pLDDT Loss**: Maximize predicted confidence
2. **PAE Loss**: Minimize alignment error (intra/inter)
3. **Contact Loss**: Maximize residue contacts (intra/inter)

### Advanced Loss Functions
```python
# Radius of Gyration Loss
# Penalizes extended conformations
rg_loss = (actual_rg - theoretical_rg) / theoretical_rg

# Interface pTM Loss
# Optimizes binding interface confidence
iptm_loss = 1 - interface_ptm

# Helix Loss
# Encourages helical structure (or β-sheets if negative)
helix_loss = distance_geometry_penalty

# Termini Distance Loss
# Promotes compact N/C terminus arrangement
termini_loss = distance(N_terminus, C_terminus)
```

---

## 8. PyRosetta Metrics

### Interface Analysis (`score_interface()`)
| Metric | Description | Good Value |
|--------|-------------|------------|
| `ShapeComplementarity` | Surface geometric fit | >0.6 |
| `dG` | Interface binding energy | <0 (negative) |
| `dSASA` | Buried surface area change | >500 Å² |
| `n_InterfaceHbonds` | Hydrogen bond count | >5 |
| `InterfaceUnsatHbonds` | Buried polar without H-bonds | <5 |
| `PackStat` | Interface packing quality | >0.6 |
| `dG/dSASA` | Energy per interface area | <-0.01 |
| `Interface_Hydrophobicity` | Hydrophobic content | 0.3-0.5 |
| `Binder_Energy_Score` | Rosetta monomer energy | Lower better |
| `Surface_Hydrophobicity` | Binder surface polarity | <0.35 |

### Structural Validation
| Function | Purpose |
|----------|---------|
| `pr_relax()` | FastRelax structure optimization |
| `align_pdbs()` | Superimpose structures |
| `unaligned_rmsd()` | Calculate RMSD without alignment |

---

## 9. Common Workflows

### Workflow 1: Standard Binder Design
```bash
# Using SLURM
sbatch ./bindcraft.slurm \
    --settings './settings_target/MyTarget.json' \
    --filters './settings_filters/default_filters.json' \
    --advanced './settings_advanced/default_4stage_multimer.json'

# Direct Python
python -u ./bindcraft.py \
    --settings './settings_target/MyTarget.json' \
    --filters './settings_filters/default_filters.json' \
    --advanced './settings_advanced/default_4stage_multimer.json'
```

### Workflow 2: Difficult Target
```bash
python -u ./bindcraft.py \
    --settings './settings_target/HardTarget.json' \
    --filters './settings_filters/relaxed_filters.json' \
    --advanced './settings_advanced/default_4stage_multimer_hardtarget.json'
```

### Workflow 3: Beta-Sheet Binder
```bash
python -u ./bindcraft.py \
    --settings './settings_target/MyTarget.json' \
    --filters './settings_filters/default_filters.json' \
    --advanced './settings_advanced/betasheet_4stage_multimer.json'
```

### Workflow 4: Peptide Design
```bash
python -u ./bindcraft.py \
    --settings './settings_target/PeptideTarget.json' \
    --filters './settings_filters/peptide_filters.json' \
    --advanced './settings_advanced/peptide_3stage_multimer.json'
```

---

## 10. Output Structure

### Directory Layout
```
design_path/
├── trajectories/           # Design trajectory animations
├── mpnn_sequences/         # MPNN-generated sequences
├── predictions/            # AF2 validation predictions
├── relaxed/               # PyRosetta relaxed structures
├── plots/                 # Design metric plots
├── final_designs/         # Designs passing all filters
│   ├── *.pdb              # Final relaxed structures
│   └── metrics.csv        # All metrics for final designs
└── logs/                  # Design logs and statistics
```

### Output Metrics CSV
Contains all computed metrics for each design:
- AF2 metrics: pLDDT, pTM, i_pTM, pAE, i_pAE
- MPNN metrics: score, sequence recovery
- Interface metrics: dG, dSASA, H-bonds, shape complementarity
- Structural metrics: RMSD, secondary structure percentages

---

## 11. Best Practices

### Target Preparation
1. **Trim PDB**: Remove unnecessary regions for speed
2. **Clean structure**: Remove waters, alternate conformations
3. **Identify hotspots**: Key binding residues (or let AF2 auto-select)
4. **Check target size**: Large targets (>600 AA) need `predict_bigbang=true`

### Design Strategy
1. **Start with defaults**: `default_4stage_multimer.json`
2. **Generate many**: Aim for 100+ passing designs minimum
3. **Screen top 5-20**: Experimentally validate best candidates
4. **Difficult targets**: Use `hardtarget` presets, expect 1000+ trajectories

### Filter Adjustment
1. **Too few designs**: Use `relaxed_filters.json`
2. **Too many designs**: Tighten thresholds manually
3. **Specific requirements**: Modify individual filter values

### Memory Management
- **32GB GPU recommended** for local installation
- **Use `greedy` algorithm** if memory limited
- **Reduce `num_recycles`** for large complexes
- **Enable `predict_bigbang`** for >600 AA targets

---

## 12. Known Limitations

1. **AF2 struggles with**:
   - Hydrophilic interfaces
   - Flexible/disordered regions
   - Very small binding sites

2. **Design artifacts**:
   - "Squashed" trajectories (automatically filtered)
   - Overfitting to AF2 confidence metrics

3. **ipTM not predictive of affinity**:
   - Good for binary binding prediction
   - Cannot rank by binding strength
   - High ipTM ≠ high affinity

4. **PyRosetta license**:
   - Commercial license required for non-academic use

---

## 13. Comparison with RFdiffusion

| Aspect | BindCraft | RFdiffusion |
|--------|-----------|-------------|
| **Approach** | AF2 backpropagation | Diffusion model |
| **Sequence** | Simultaneous | MPNN post-design |
| **Validation** | Built-in AF2 | Separate step |
| **Speed** | Slower (gradient through AF2) | Faster (forward pass) |
| **Memory** | Higher (backprop) | Lower |
| **Interface quality** | Often better metrics | Varies |
| **Structural diversity** | Limited by AF2 | Higher diversity |

### When to Use BindCraft
- Need high-confidence AF2 predictions
- Want integrated validation pipeline
- Targeting well-defined binding sites
- Quality over quantity

### When to Use RFdiffusion
- Need structural diversity
- Large-scale screening
- Novel fold exploration
- Memory constrained

---

## 14. Quick Reference Card

```
EXECUTION:
  python bindcraft.py --settings X.json --filters Y.json --advanced Z.json

TARGET SETTINGS:
  chains: "A"                      → Target chain(s)
  target_hotspot_residues: "56"    → Binding site (null = auto)
  lengths: [65, 150]               → Binder length range
  number_of_final_designs: 100     → Target design count

ALGORITHM OPTIONS:
  2stage, 3stage, 4stage (default), greedy, mcmc

KEY WEIGHTS:
  weights_con_inter: 1.0           → Interface contacts
  weights_plddt: 0.1               → Confidence
  weights_helicity: -0.3           → Negative = β-sheets

KEY FILTERS:
  pLDDT ≥ 0.8                      → Confidence
  i_pTM ≥ 0.5                      → Interface confidence
  dG ≤ 0                           → Binding energy
  ShapeComplementarity ≥ 0.6       → Surface fit

PRESETS:
  default_4stage_multimer.json     → Standard
  betasheet_4stage_multimer.json   → β-sheet binders
  peptide_3stage_multimer.json     → Peptides
  *_hardtarget.json                → Difficult targets
  *_flexible.json                  → Flexible interface
```

---

## 15. Concepts to Apply to Banta Lab Development

### Applicable Concepts
1. **Multi-stage optimization**: Soft → temp → hard → greedy progression
2. **Combined loss functions**: pLDDT + pAE + contacts + custom losses
3. **Integrated validation**: AF2 prediction after MPNN design
4. **PyRosetta interface scoring**: dG, dSASA, H-bonds, shape complementarity
5. **Filter-based selection**: Multi-metric thresholds for design quality
6. **Acceptance rate monitoring**: Auto-terminate unsuccessful campaigns

### Potential Improvements for Banta Lab
1. **Add AF2 validation step** after RFD3 + LigandMPNN
2. **Implement PyRosetta interface metrics** for protein binders
3. **Add multi-stage design** to metal site optimization
4. **Use loss function weights** for balancing coordination vs stability
5. **Add acceptance rate monitoring** for long design runs

---

*Document generated: 2026-01-15*
*Source: Official BindCraft repository (martinpacesa/BindCraft)*
