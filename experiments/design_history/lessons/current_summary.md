# Design Lessons (Auto-updated)

> This file is automatically refreshed when significant patterns are detected.
> Last updated: 2026-01-31

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

## Ligand Interface

### Citrate Coordination

- Auto-assignment of all O atoms as H-bond acceptors is effective (O1-O7)
- Missing explicit H-bond donor conditioning (O7) didn't hurt when combined with burial
- Ligand SMILES for RF3: `OC(=O)CC(O)(CC(O)=O)C(O)=O`

## Failure Patterns to Avoid

1. **No MPNN amino acid bias** -> 17%+ alanine content -> poor folding (pTM variance 0.13+)
2. **Sequence-only RF3 for metal designs** -> 15% of designs below pTM 0.5 threshold
3. **Over-constrained motif scaffolds** -> Reduces designable space, doesn't improve quality
4. **RF3 CIF res_name assumption** -> CIF output uses `L:0`/`L:1` not standard codes; use atom_name-based detection

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
