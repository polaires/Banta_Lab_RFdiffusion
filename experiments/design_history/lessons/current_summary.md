# Design Lessons (Auto-updated)

> This file is automatically refreshed when significant patterns are detected.
> Last updated: 2026-02-02

## Metal Coordination

### Tb-Citrate (Rounds 1-8)

**Best approach**: NL pipeline automated defaults via `metal_binding_design` single mode.

| Metric | Expert-Tuned (R7b) | NL Pipeline (R8) | Delta |
|--------|---------------------|-------------------|-------|
| Avg pTM | 0.771 | 0.902 | +17% |
| Avg pLDDT | 0.801 | 0.867 | +8% |
| Avg PAE | 6.22 | 3.61 | -42% |
| Pass rate | 85% | 100% | +15pp |
| Within-backbone std | 0.126 | 0.008 | 15.7x better |
| Alanine % | 17.1% | 6.7% | 2.5x less |

**Key lessons**:
1. **MPNN bias is critical**: `bias_AA: "A:-2.0"` + `omit_AA: "C"` reduced alanine from 17% to 6.7%. Without bias, MPNN defaults to alanine-heavy sequences that fold poorly.
2. **Ligand-aware RF3 is transformative**: Providing `ligand_smiles` + `metal` to RF3 yields uniformly high confidence (worst pTM 0.86). Sequence-only RF3 produces 3/20 designs below pTM 0.5.
3. **Auto-templates outperform hand-crafted motifs**: Fewer constraints = more designable space. Hand-crafted 5-residue motif scaffolds over-constrain the problem.
4. **Lower CFG is sufficient with good MPNN**: CFG 2.0 + quality sequences beats CFG 2.5 + poor sequences.
5. **Side chain packing matters**: `pack_side_chains: true` in MPNN improves coordination geometry.
6. **Bury both metal and ligand**: `bury_ligand: true` didn't hurt H-bond access despite initial concerns.
7. **MPNN temperature 0.1 > 0.2 for metals**: Lower temperature produces more consistent sequences when combined with bias.

## Production Statistics (Final — 10 Batches, 2026-02-02)

| Metric | Cumulative |
|--------|-----------|
| Backbones | 500 |
| Sequences | 1,800 |
| After stability | 899 (50%) |
| Pass metal_binding | 173 (19.2%) |
| Pass strict | 0 |
| CN=8+ designs | 21 |
| Best composite score | 0.622 (b01_seq_0180) |

**Per-batch pass rates**: B1: 18.4%, B2: 21.0%, B3: 20.0%, B4: 30.6%, B5: 12.2%*, B6: 9.7%*, B7: 19.4%, B8: 17.6%, B9: 20.0%, B10: 28.8%
(*B5-B6 used partial co-diffusion — harmful, reverted)

**Excluding co-diffusion experiment (B5-B6)**: 152/808 = 18.8% avg pass rate

**CN distribution**: CN=6 (53%), CN=7 (36%), CN=8 (9%), CN=9+ (2%)

## AF3 Cross-Validation (2026-02-02)

All 173 metal_binding-passing designs screened through AlphaFold3 v3.0.1 (noMSA, RTX 5090).

### Tiered Screening Results

| Phase | Designs | Seeds | Ca Pass | Tb Pass | Both Pass |
|-------|---------|-------|---------|---------|-----------|
| 1. Screen | 172 | 1 | 42 (24.4%) | 132 (76.7%) | 41 (23.8%) |
| 2. Validate | 65 | 3 | 53 (81.5%) | 64 (98.5%) | 52 (80.0%) |
| 3. Final | 10 | 5 | 10 (100%) | 10 (100%) | 10 (100%) |

**Total: 571 AF3 jobs in 145 min** (batch mode via `--input_dir`, ~8s/job)

### Top 10 Candidates for Gene Synthesis

| Rank | Design | RF3 pTM | Ca iPTM | Tb iPTM | Holo-Apo Delta |
|------|--------|---------|---------|---------|----------------|
| 1 | b09_seq_0016 | 0.893 | 0.91 | 0.91 | +0.09 |
| 2 | b01_seq_0094 | 0.895 | 0.89 | 0.89 | +0.13 |
| 3 | b07_seq_0015 | 0.867 | 0.89 | 0.83 | +0.29 |
| 4 | b01_seq_0180 | 0.908 | 0.89 | 0.82 | +0.06 |
| 5 | b02_seq_0105 | 0.818 | 0.92 | 0.90 | +0.15 |

### Key AF3 Findings

1. **Ca and Tb produce identical AF3 predictions**: Mean iPTM 0.703 vs 0.705. Ca is a reliable Tb surrogate.
2. **RF3 moderately predicts AF3**: r=0.40 correlation. RF3 useful for triage but doesn't fully predict AF3.
3. **CN is NOT predictive**: CN=6 designs pass Ca at 27% vs CN=7 at 19%. High CN reflects RF3 scoring biases.
4. **Multi-seed averaging transforms results**: 1-seed 24.4% → 3-seed 81.5% → 5-seed 100% pass rate.
5. **Strong holo-apo delta confirms binding**: Avg +0.201 (Phase 2), +0.151 (Phase 3). AF3 recognizes metal-ligand binding.

## Ligand Interface

### Citrate Coordination

- **CONFIRMED (B5-B6)**: Partial ligand co-diffusion (`ligand_fix_atoms: "O2,O5,O7"`) HURTS pass rate (10.9% vs 18.8%). All-fixed is correct.
- Coordinating oxygens (O2, O5, O7) should NOT be H-bond acceptors — their lone pairs are occupied by metal coordination bonds. Only non-coordinating oxygens (O1, O3, O4, O6) should be set as `select_hbond_acceptor`.
- However, when `ligand_fix_atoms` is NOT set (production default), all O atoms are used as H-bond acceptors. This produces ~19% pass rate and is acceptable.
- Ligand SMILES for RF3: `OC(=O)CC(O)(CC(O)=O)C(O)=O`

## Failure Patterns to Avoid

1. **No MPNN amino acid bias** -> 17%+ alanine content -> poor folding (pTM variance 0.13+)
2. **Sequence-only RF3 for metal designs** -> 15% of designs below pTM 0.5 threshold
3. **Over-constrained motif scaffolds** -> Reduces designable space, doesn't improve quality
4. **RF3 CIF res_name assumption** -> CIF output uses `L:0`/`L:1` not standard codes; use atom_name-based detection
5. **H-bond on coordinating atoms** -> Conditioning H-bonds on metal-coordinating O atoms conflicts with metal coordination geometry. Exclude them from `select_hbond_acceptor`.
6. **50 backbone API calls timeout** -> Sub-batch at 10/call with 900s timeout

---

## How Lessons Are Captured

Lessons are automatically synthesized when:
1. **Failure pattern detected**: 3+ consecutive failures with similar parameters
2. **Breakthrough success**: New best achieved on key metrics
3. **Meaningful improvement**: Significant improvement (>15%) from previous best

To contribute to lessons, run local designs using:
```bash
python analyze_design_cli.py your_design.pdb --metal TB --session exploration_name
```
