# Tb-Citrate Scaffold Production Report

**Date**: 2026-02-02
**Pipeline**: RFD3 → Scout → LigandMPNN → RF3 → Stability → Analysis
**Filter Preset**: metal_tb_standard (CN≥6, pLDDT≥0.70, geoRMSD≤2.0)

## Production Parameters

```json
{
  "cfg_scale": 1.5,
  "contig": "100-130",
  "bury_ligand": true,
  "num_timesteps": 200,
  "step_scale": 1.0,
  "mpnn_temperature": 0.05,
  "mpnn_bias_AA": "A:-2.0",
  "mpnn_omit_AA": "C"
}
```

## Pipeline Funnel

| Stage | Count | Rate |
|-------|-------|------|
| Backbones generated | 500 | — |
| After scout pre-filter | 500 | 100% |
| Sequences designed | 1,800 | 4 seq/bb |
| After stability filter | 899 | 50.0% |
| Pass metal_binding | 173 | 19.2% |
| Pass strict | 0 | 0% |

## Per-Batch Results

| Batch | BB | Seqs | Stability | Pass | Rate | Notes |
|-------|----|------|-----------|------|------|-------|
| B1 | 50 | 200 | 98 | 18 | 18.4% | Baseline |
| B2 | 50 | 200 | 105 | 22 | 21.0% | — |
| B3 | 50 | 100 | 50 | 10 | 20.0% | 2 seq/bb test |
| B4 | 50 | 100 | 49 | 15 | 30.6% | Best rate (B1-B4) |
| B5 | 50 | 200 | 90 | 11 | 12.2% | Partial co-diffusion |
| B6 | 50 | 200 | 103 | 10 | 9.7% | Co-diffusion confirmed harmful |
| B7 | 50 | 200 | 108 | 21 | 19.4% | Reverted to all-fixed |
| B8 | 50 | 200 | 91 | 16 | 17.6% | — |
| B9 | 50 | 200 | 101 | 20 | 20.0% | — |
| B10 | 50 | 200 | 104 | 30 | 28.8% | Best batch overall |

**Excluding B5-B6 (co-diffusion experiment)**: 152/808 = 18.8% avg pass rate

## Coordination Number Distribution (173 passing)

| CN | Count | % |
|----|-------|---|
| 6 | 92 | 53.2% |
| 7 | 62 | 35.8% |
| 8 | 16 | 9.2% |
| 9 | 2 | 1.2% |
| 10 | 1 | 0.6% |

## Top 10 Candidates

| Rank | Design | pTM | pLDDT | CN | LC | Ala% | Score |
|------|--------|-----|-------|----|----|------|-------|
| 1 | b01_seq_0180 | 0.908 | 0.864 | 7 | 14 | 3.9% | 0.622 |
| 2 | b09_seq_0138 | 0.901 | 0.864 | 7 | 13 | 3.1% | 0.620 |
| 3 | b01_seq_0094 | 0.895 | 0.867 | 6 | 8 | 0.9% | 0.618 |
| 4 | b09_seq_0016 | 0.893 | 0.868 | 6 | 9 | 0.0% | 0.617 |
| 5 | b02_seq_0177 | 0.904 | 0.849 | 7 | 13 | 3.1% | 0.617 |
| 6 | b01_seq_0176 | 0.884 | 0.871 | 7 | 12 | 0.9% | 0.615 |
| 7 | b03_seq_0079 | 0.894 | 0.855 | 7 | 17 | 0.9% | 0.614 |
| 8 | b01_seq_0096 | 0.885 | 0.861 | 6 | 9 | 1.9% | 0.613 |
| 9 | b05_seq_0114 | 0.878 | 0.869 | 7 | 4 | 3.5% | 0.612 |
| 10 | b02_seq_0130 | 0.888 | 0.853 | 7 | 16 | 0.9% | 0.611 |

## Key Findings

### Partial Ligand Co-Diffusion (B5-B6)
- Setting `ligand_fix_atoms="O2,O5,O7"` to only fix coordinating oxygens while co-diffusing the rest
- Result: Pass rate dropped from ~21% to ~10%
- Cause: RFD3 displaces non-fixed ligand atoms, distorting citrate geometry
- Lesson: Keep all ligand atoms fixed for small rigid ligands like citrate

### H-Bond Acceptor Conflict
- Coordinating oxygens (O2,O5,O7) should NOT be H-bond acceptors
- Their lone pairs are occupied by metal coordination bonds
- Fixed `get_ligand_hbond_atoms()` to exclude coordinating atoms when `ligand_fix_atoms` is set
- Impact unclear since B7-B10 ran without `ligand_fix_atoms` (all fixed, all as H-bond acceptors)

### Strict Filter Gap
- 0/173 designs passed strict filter across 10 batches
- Strict requires: geoRMSD≤1.0, pLDDT≥0.80, CN≥8, protein_coordination≥5, ligand_contacts≥5
- Most designs have CN=6-7 (TB empirically achieves 6-9 in monomers vs formal 8-9)
- geoRMSD is null (no reference geometry available for comparison)
- Recommendation: Relax strict CN threshold to ≥6 for TB, or accept standard-passing as sufficient

## Output Files

```
final/
├── top_candidates/           # Top 10 organized
│   ├── 01_b01_seq_0180/
│   │   ├── rf3_prediction.pdb
│   │   ├── sequence.fasta
│   │   └── metrics.json
│   ├── 02_b09_seq_0138/
│   │   └── ...
│   └── ... (10 directories)
├── top_10_summary.json       # Summary with all metrics
├── promising/                # All 173 passing PDBs
└── PRODUCTION_REPORT.md      # This file

analysis/
├── production_batch_01-10.json  # Per-batch detailed results
├── best_params.json             # Locked production parameters
└── cumulative_lessons.md        # All lessons learned
```

## Next Steps

1. **AF3 Cross-Validation**: Run top 10 through AlphaFold3 with Ca-citrate (primary) and Tb-citrate (secondary)
   - AF3 inputs prepared at `~/af_input/tb_citrate_validation/`
   - Requires stopping Docker container to free GPU
2. **Experimental Testing**: Order gene synthesis for top candidates
3. **Iteration**: If AF3 validation identifies issues, adjust parameters and run additional batches
