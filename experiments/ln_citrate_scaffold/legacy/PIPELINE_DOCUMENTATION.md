# TB-Citrate Binding Protein Design Pipeline

## Overview

This pipeline automates the design of terbium-citrate binding proteins using RFdiffusion3 (RFD3) scaffolding, LigandMPNN sequence design, and RF3 structure validation.

## Pipeline Stages

```
┌─────────────────────────────────────────────────────────────────────────┐
│  PARAMETER SWEEP (10 trials × 9 configs = 90 designs)                   │
│  ┌─────────────┐    ┌─────────────┐    ┌─────────────┐                 │
│  │ Small       │    │ Medium      │    │ Large       │                 │
│  │ 100-120 aa  │    │ 110-130 aa  │    │ 130-150 aa  │                 │
│  ├─────────────┤    ├─────────────┤    ├─────────────┤                 │
│  │ CFG 1.5     │    │ CFG 1.5     │    │ CFG 1.5     │                 │
│  │ CFG 2.0     │    │ CFG 2.0     │    │ CFG 2.0     │                 │
│  │ CFG 2.5     │    │ CFG 2.5     │    │ CFG 2.5     │                 │
│  └─────────────┘    └─────────────┘    └─────────────┘                 │
│                              │                                          │
│                              ▼                                          │
│                    ┌─────────────────┐                                  │
│                    │ RANK BY:        │                                  │
│                    │ - Pass rate     │                                  │
│                    │ - Avg pLDDT     │                                  │
│                    │ - Avg pTM       │                                  │
│                    └────────┬────────┘                                  │
│                              │                                          │
│                              ▼                                          │
│  PRODUCTION RUN (1000 designs with best config)                        │
│  ┌──────────────────────────────────────────────────────────────────┐  │
│  │ Best Config → RFD3 → MPNN → RF3 Filter → Save Passing Designs   │  │
│  └──────────────────────────────────────────────────────────────────┘  │
└─────────────────────────────────────────────────────────────────────────┘
```

## Filter Criteria

### Stage 3: RF3 Validation (Primary Filter)

| Metric | Strict | Relaxed | Description |
|--------|--------|---------|-------------|
| pLDDT | > 0.80 | > 0.75 | Fold confidence (per-residue) |
| pTM | > 0.80 | > 0.70 | Template matching score |
| PAE | < 5.0 | < 10.0 | Predicted aligned error (Å) |

### Stage 4: Sequence Analysis (Secondary Filter)

| Metric | Target | Flag | Description |
|--------|--------|------|-------------|
| Carboxylate (E+D) | > 15% | < 10% | Metal coordination donors |
| Alanine (A) | < 20% | > 25% | Low diversity indicator |
| Binding Score | > 60 | < 40 | Combined binding potential |

## Key Parameters

### RFD3 Scaffolding

```python
{
    "task": "rfd3",
    "contig": "130-150",              # Protein length range
    "ligand": "CIT,TB",               # Ligands to include
    "select_fixed_atoms": {
        "X1": "all",                  # Fix TB position
        "L1": "all"                   # Fix citrate position
    },
    "select_buried": {"X1": "all"},   # Bury TB for pocket
    "select_hbond_acceptor": {
        "L1": "O1,O2,O3,O4,O5,O6"    # Citrate H-bond acceptors
    },
    "select_hbond_donor": {
        "L1": "O7"                    # Citrate hydroxyl donor
    },
    "use_classifier_free_guidance": True,
    "cfg_scale": 2.5,                 # Conditioning strength
    "infer_ori_strategy": "com",      # Center of mass orientation
}
```

### LigandMPNN Sequence Design

```python
{
    "task": "mpnn",
    "ligand_mpnn_use_atom_context": 1,  # Ligand-aware design
    "temperature": 0.2,                  # Lower = more deterministic
    "num_seqs": 8,                       # Sequences per backbone
}
```

## Usage

### 1. Parameter Sweep (Find Best Config)

```bash
python pipeline_tb_citrate.py --mode sweep --trials 10
```

This tests 9 configurations with 10 trials each:
- 3 size ranges: small (100-120), medium (110-130), large (130-150)
- 3 CFG scales: 1.5, 2.0, 2.5

Output: `sweep_results.json` with rankings

### 2. Production Run (Generate Designs)

```bash
python pipeline_tb_citrate.py --mode production --num-designs 1000 --config large_high_cfg
```

Output:
- `all_passing.fasta` - All designs passing strict filters
- `production_summary.json` - Statistics and best designs
- Individual FASTA files for each passing design

## Results from Manual Optimization

| Config | Contig | CFG | Best pLDDT | Best pTM | Pass Rate |
|--------|--------|-----|------------|----------|-----------|
| R7 | 110-130 | 2.0 | 0.865 | 0.891 | ~25% |
| **R7b** | **130-150** | **2.5** | **0.874** | **0.938** | **~50%** |

**Best configuration: `large_high_cfg` (130-150 aa, CFG 2.5)**

## Expected Pass Rates

Based on our experiments:
- Strict filters (pLDDT>0.80, pTM>0.80): ~30-50%
- Relaxed filters (pLDDT>0.75, pTM>0.70): ~60-80%

For 1000 generated designs:
- Expect 300-500 designs passing strict filters
- Top 10-20 designs for AF3 validation

## Files Structure

```
experiments/ln_citrate_scaffold/
├── inputs/
│   └── tb_citrate_motif_scaffold.pdb    # Input complex
├── scripts/
│   └── pipeline_tb_citrate.py           # Automated pipeline
├── outputs/
│   ├── round_07_scaffold/               # R7 results
│   ├── round_07b_improved/              # R7b results
│   ├── sweep_YYYYMMDD_HHMMSS/          # Parameter sweep results
│   └── production_YYYYMMDD_HHMMSS/     # Production run results
└── PIPELINE_DOCUMENTATION.md            # This file
```

## Next Steps After Pipeline

1. **AF3 Validation**: Submit top 10-20 designs to AlphaFold3
2. **Structural Analysis**: Check citrate-TB coordination preservation
3. **Experimental Validation**: Express and test binding affinity

## Troubleshooting

| Issue | Solution |
|-------|----------|
| 0 designs generated | Check Docker is running: `docker ps` |
| Low pass rate (<10%) | Increase CFG scale or adjust contig range |
| High Ala content (>25%) | Lower MPNN temperature or increase num_seqs |
| Connection errors | Restart Docker: `docker compose up -d` |
