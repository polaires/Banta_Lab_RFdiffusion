# 3C9H Scaffolding Production Plan (REVISED)

**Date**: 2026-02-03 (Updated after structural analysis)
**Goal**: Produce 400 Tb-citrate designs by scaffolding the 3C9H Mg-citrate coordination motif
**Hypothesis**: Scaffolding a crystallographically-validated coordination geometry will improve metal binding success

---

## Structural Analysis Summary

### 3C9H Binding Site (Mg-Citrate, 1.9 A resolution)

**Direct Mg Coordination (CN=6):**
| Donor | Distance | Source |
|-------|----------|--------|
| Ser41.O (backbone!) | 2.27 A | Protein |
| Asp43.OD1 | 2.32 A | Protein |
| Citrate O2 | 2.48 A | Ligand |
| Citrate O5 | 2.61 A | Ligand |
| Citrate O7 | 2.36 A | Ligand |
| Water (HOH414) | 2.28 A | Solvent |

**Key insight**: The protein contribution is only Ser41 backbone O and Asp43 sidechain - a **3-residue span** (41-43). This is compact and extractable!

**Citrate H-bond network** (provides selectivity):
- Thr70, Ser94, Ser214, Thr216, Tyr238, Arg263
- Spatially close despite sequence separation (Ser41-Thr70 = 5.8 A in space)

### Tb vs Mg Considerations

| Property | Mg2+ | Tb3+ |
|----------|------|------|
| Ionic radius | 0.72 A | 1.04 A |
| Preferred CN | 6 (octahedral) | 8-9 |
| HSAB | Hard | Hard |
| Donor preference | O > N | O > N |

**Challenge**: The water position (HOH414) in 3C9H could be replaced by additional protein donors for Tb, which RFD3 can design.

---

## Motif Files Created

```
experiments/ln_citrate_scaffold/scaffolding_3c9h/inputs/
├── 3c9h.pdb                      # Full structure (334 aa × 2 chains)
├── 3c9h_ultra_minimal_motif.pdb  # Ser41-Leu42-Asp43 + TB + CIT (original numbering)
├── 3c9h_minimal_motif.pdb        # + Thr70 (8 residues, 2 segments)
└── 3c9h_rfd3_motif.pdb           # ★ RFD3-formatted: Chain M (Ser1-Leu2-Asp3), X (TB), L (CIT)
```

**Use `3c9h_rfd3_motif.pdb` for all RFD3 runs** — properly formatted with:
- Chain M: Protein motif (renumbered 1-3)
- Chain X: Metal (TB)
- Chain L: Ligand (CIT)

---

## RFD3 Paper Key Findings (Applied)

### Sources
- [RFD3 Paper: bioRxiv 2025.09.18.676967](https://www.biorxiv.org/content/10.1101/2025.09.18.676967v2)
- [IPD Announcement](https://www.ipd.uw.edu/2025/12/rfdiffusion3-now-available/)
- [Ailurus Deep Dive](https://www.ailurus.bio/post/the-atomic-era-of-protein-design-a-deep-dive-into-rfdiffusion3)

### Enzyme/Motif Scaffolding Performance
- **Success rate**: RFD3 scaffolds catalytic motifs in **90% of cases** (cysteine hydrolase benchmark)
- **Multi-island performance**: Outperforms RFD2 on **37/41** AME (Atomic Motif Enzyme) benchmarks
- **Pass rate improvement**: 15% vs 4% for multi-island motifs compared to RFD2
- **Experimental validation**: 35/190 cysteine hydrolase designs showed catalytic activity (Cys-His-Asp triad)

### Unindexed Atomic Motifs
The `unindex` parameter allows RFD3 to place the motif at any position in the generated protein:
- During training, subsets of residues are duplicated with index features removed
- Network learns that unindexed residues superimpose on indexed residues in unnoised structures
- By providing coordinates for unindexed residues, RFD3 models proteins conditioned on known substructure at unknown indices
- **Advantage**: No need to pre-specify where motif should appear in sequence

### Key Parameters for Enzyme Scaffolding
| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `n_batches` | ~100 | Paper recommendation for motif scaffolding |
| `sequences/backbone` | 8 | Higher than PPI (4) due to constraint complexity |
| `cfg_scale` | 2.0 | Required with H-bond + RASA conditioning |
| `step_scale` | 1.5 | Production mode (1.0 for exploration diversity) |
| `diffusion_batch_size` | 1 | For motif scaffolding |

### Success Thresholds (from Paper)
| Metric | Threshold | Note |
|--------|-----------|------|
| Motif all-atom RMSD | < 1.5 Å | Primary success criterion |
| Backbone RMSD | ≤ 1.5 Å | For ligand binding |
| Ligand RMSD | < 5 Å | Backbone-aligned |
| Interface min PAE | ≤ 1.5 | For binding |
| iPTM | ≥ 0.8 | Interface confidence |

---

## Background & Motivation

### De Novo Production Results (B1-B10)
| Metric | Value |
|--------|-------|
| Total backbones | 500 |
| Total sequences | 1,800 |
| Pass metal_tb_standard | 173 (19.2%) |
| Top pTM | 0.91 |
| Top pLDDT | 0.87 |
| Typical CN | 6-7 |
| AF3 validation | 10/10 passed both Ca and Tb |

**Limitation identified**: High global confidence (pTM/pLDDT) but potentially low *local* confidence at the binding pocket. De novo pockets are "invented" — the model is confident about the fold but uncertain about the novel coordination geometry.

### Why 3C9H Scaffolding?

**3C9H** is an ABC transporter substrate-binding protein (355 aa, 1.9 Å resolution) that natively binds **Mg²⁺-citrate**. This is an ideal scaffold because:

1. **Proven citrate binding** — The pocket geometry is crystallographically validated
2. **HSAB compatible** — Mg (hard) → Tb (hard), same O-donor preference
3. **Natural fold** — High baseline confidence for overall structure
4. **Established coordination** — Known coordinating residues (Asp/Glu) preserved

### Key Question

Will scaffolding 3C9H produce designs with:
- Higher local pLDDT at binding site residues?
- Better AF3 iptm for metal-citrate interface?
- Similar or better overall pass rate?

---

## Lessons from De Novo Production (Applied)

### Critical Parameters (Locked)
```json
{
  "mpnn_temperature": 0.05,
  "mpnn_bias_AA": "A:-2.0",
  "mpnn_omit_AA": "C",
  "rf3_with_ligand_context": true
}
```

### Lessons Applied to Scaffolding

| Lesson | De Novo Finding | Scaffolding Adaptation |
|--------|-----------------|------------------------|
| **MPNN bias critical** | A:-2.0 reduces Ala from 17%→6.7% | Same — applies to scaffold designs |
| **Omit Cys** | Tb is hard acid, Cys is soft | Same — omit C in all MPNN |
| **RF3 with ligand** | +15% pass rate with metal context | Same — always use ligand-aware RF3 |
| **All-fixed ligand** | Partial co-diffusion hurt rate (21%→10%) | Same — fix all citrate atoms |
| **CFG 1.5** | Sufficient with good MPNN | Start at 1.5, may need 2.0 for motif |
| **4 seq/bb** | Better coverage than 2 seq/bb | Same |

### New Considerations for Scaffolding

| Factor | De Novo | Scaffolding | Impact |
|--------|---------|-------------|--------|
| Contig | `100-130` (full protein) | `50-80` scaffold around motif | Smaller, focused design space |
| Fixed residues | None (metal+ligand only) | Coordinating residues from 3C9H | Constrained backbone |
| Pocket geometry | Invented by RFD3 | Crystallographic (1.9 Å) | Higher local confidence expected |
| Metal substitution | N/A | MG → TB | HSAB compatible |

---

## Phase 1: Exploration (100 backbones)

### 1.1 Motif Scaffolding Approach (from RFD3 Paper + Skill Reference)

**Key RFD3 capabilities for enzyme/motif scaffolding:**
- **Unindexed atomic motifs**: Motif can be placed anywhere in sequence (no index features)
- **Multi-island motifs**: 15% pass rate vs 4% for RFD2 (outperforms on 37/41 AME benchmarks)
- **Native sidechain diffusion**: 14-atom representation overlaps with fixed motif atoms
- **H-bond conditioning**: Improves active site geometry (33% → 37% with CFG)
- **RASA conditioning**: `select_buried` for ligand pocket creation

**Paper validation**: "RFD3 successfully scaffolded the catalytic motif in 90% of cases" for cysteine hydrolase design.

Using the ultra-minimal motif (Ser41-Leu42-Asp43 + TB + CIT):

```
RFD3 Motif Scaffolding with Unindexed Placement:

    The motif (M chain) can be placed anywhere optimal in the new protein.
    RFD3 will find the best structural context via learned superposition.

    [NEW SCAFFOLD ~100-140 residues]
              ↓
         ┌─────────┐
         │  Ser1   │ ← Fixed: backbone O coordinates TB (2.27 Å)
         │  Leu2   │ ← Fixed: provides hydrophobic context
         │  Asp3   │ ← Fixed: OD1 coordinates TB (2.32 Å)
         └────┬────┘
              │
           TB─CIT (fixed, buried)
```

**Two scaffolding modes:**

1. **Indexed scaffolding** - Motif at specified position in contig:
```
contig: "30-50,M1-3,30-50"

Where:
  30-50    = Design 30-50 new residues (N-terminal segment)
  M1-3     = Fixed motif residues (Ser-Leu-Asp on chain M)
  30-50    = Design 30-50 new residues (C-terminal segment)

Total: 63-103 residues (including 3 fixed)
```

2. **Unindexed scaffolding** - RFD3 decides optimal placement:
```python
{
    "input": "3c9h_rfd3_motif.pdb",
    "contig": "100-140",  # Just specify total length
    "unindex": {"M": "all"},  # Motif position not constrained
    "select_fixed_atoms": {"M": "all", "X": "all", "L": "all"}
}
```

**How unindexed works** (from paper): During training, subsets of residues are duplicated with index features removed. The network learns that unindexed residues are always superimposed on indexed residues in unnoised structures. By providing coordinates for unindexed residues, RFD3 models the distribution of proteins conditioned on inclusion of a known substructure at unknown indices.

### 1.2 Exploration Configs (Based on RFD3 Best Practices)

**Paper recommendation**: ~100 backbones per motif, 8 sequences/backbone
**Success threshold**: Motif all-atom RMSD < 1.5 Å (validated with Chai-1 or AF3)

| Config | Contig | Unindex | CFG | step_scale | Description |
|--------|--------|---------|-----|------------|-------------|
| E1 | `35-45,M1-3,35-45` | No | 2.0 | 1.0 | Indexed, small (73-93 res), diverse |
| E2 | `45-55,M1-3,45-55` | No | 2.0 | 1.5 | Indexed, medium (93-113 res), designable |
| E3 | `100-130` | Yes | 2.0 | 1.0 | **Unindexed**, diverse |
| E4 | `100-130` | Yes | 2.0 | 1.5 | **Unindexed**, designable |

**Step scale trade-off**:
- η=1.0: Higher diversity, lower designability (exploration)
- η=1.5: Higher designability, lower diversity (production)

### 1.3 Fixed Atoms & Conditioning Configuration

**Input PDB**: `3c9h_rfd3_motif.pdb` (properly formatted with chain M/X/L)

```json
{
  "input": "3c9h_rfd3_motif.pdb",

  "select_fixed_atoms": {
    "M": "all",     // Protein motif (Ser1-Leu2-Asp3 on chain M)
    "X": "all",     // Metal (TB)
    "L": "all"      // Ligand (CIT)
  },

  "unindex": {
    "M": "all"      // Let RFD3 find optimal sequence position
  },

  "select_buried": {
    "X": "all",     // Bury metal in pocket
    "L": "all"      // Bury ligand in pocket
  },

  "select_hbond_acceptor": {
    "L": "O1,O3,O4,O6"   // Non-coordinating citrate oxygens for H-bond network
  },

  "use_classifier_free_guidance": true,
  "cfg_scale": 2.0,

  "infer_ori_strategy": "com",  // Center on fixed atoms CoM

  "n_batches": 100,
  "diffusion_batch_size": 1,  // For motif scaffolding
  "num_timesteps": 200,
  "step_scale": 1.5,
  "gamma_0": 0.6
}
```

**Citrate oxygen assignments** (from crystal structure):
- **Coordinating TB** (exclude from H-bond acceptors): O2 (2.48 Å), O5 (2.61 Å), O7 (2.36 Å)
- **Available for H-bonds**: O1, O3, O4, O6

### 1.4 Motif PDB Preparation (COMPLETED)

The motif has been properly relabeled in `3c9h_rfd3_motif.pdb`:
- **Chain M**: Protein motif (Ser1-Leu2-Asp3, renumbered from 41-43)
- **Chain X**: Metal (TB at original MG position)
- **Chain L**: Ligand (CIT with all atoms)

File location: `experiments/ln_citrate_scaffold/scaffolding_3c9h/inputs/3c9h_rfd3_motif.pdb`

### 1.5 Exploration Metrics

**Key metrics for motif scaffolding**:

1. **Motif all-atom RMSD** (primary): Must be < 1.5 Å for success
2. **Binding site local pLDDT**: Compare to global pLDDT
3. **Metal coordination number**: Target CN ≥ 6 for Tb

```python
def compute_binding_site_plddt(rf3_cif, coordinating_residues):
    """Extract per-residue pLDDT for coordinating residues only."""
    # Parse B_iso_or_equiv column from CIF
    # Filter to coordinating residue atoms
    # Return mean pLDDT for binding site
    pass

# Compare:
# - Global pLDDT (whole protein average)
# - Binding site pLDDT (coordinating residues only)
# - Binding site delta = binding_site_plddt - global_plddt
```

**Success criteria for exploration** (from RFD3 paper):
- **Motif RMSD < 1.5 Å** (must achieve — primary criterion)
- Pass rate ≥ 15% (RFD3 paper baseline for multi-island)
- Binding site pLDDT ≥ global pLDDT (the key hypothesis)
- At least one config produces CN ≥ 6

---

## Phase 2: Production (400 backbones)

### 2.1 Production Parameters (RFD3 Foundry API Format)

Based on exploration, lock best parameters. Initial estimate using **unindexed scaffolding** (E3/E4 approach):

**RFD3 Configuration:**
```json
{
  "input": "3c9h_rfd3_motif.pdb",
  "contig": "100-130",

  "select_fixed_atoms": {
    "M": "all",
    "X": "all",
    "L": "all"
  },

  "unindex": {
    "M": "all"
  },

  "select_buried": {
    "X": "all",
    "L": "all"
  },

  "select_hbond_acceptor": {
    "L": "O1,O3,O4,O6"
  },

  "use_classifier_free_guidance": true,
  "cfg_scale": 2.0,
  "infer_ori_strategy": "com",

  "n_batches": 400,
  "diffusion_batch_size": 1,
  "num_timesteps": 200,
  "step_scale": 1.5,
  "gamma_0": 0.6
}
```

**LigandMPNN Configuration:**
```json
{
  "temperature": 0.05,
  "bias_AA": "A:-2.0",
  "omit_AA": "C",
  "fixed_residues_multi": "M1,M2,M3",
  "sequences_per_backbone": 8
}
```

**Key difference from de novo**: The Ser-Leu-Asp motif (chain M) provides a crystallographically-validated starting geometry for the metal coordination sphere. RFD3's unindexed scaffolding finds the optimal position for this motif in the generated protein.

### 2.2 Batch Structure

| Batch | Backbones | Sequences | Notes |
|-------|-----------|-----------|-------|
| S1-S4 | 50 each | 200 each | Baseline scaffolding |
| S5 | 50 | 200 | Adjust based on S1-S4 lessons |
| S6-S8 | 50 each | 200 each | Continue with best params |
| **Total** | **400** | **1,600** | |

### 2.3 Pipeline Per Batch

```
1. RFD3 scaffold generation (10/sub-batch × 5 = 50 backbones)
   ↓
2. Scout filter (1 seq/bb, pTM ≥ 0.6, pLDDT ≥ 0.65)
   ↓ ~60% pass
3. LigandMPNN (4 seq/bb on passing backbones)
   ↓
4. Stability filter (Ala < 25%, pTM ≥ 0.6)
   ↓ ~50% pass
5. RF3 validation (ligand + metal context)
   ↓
6. Metal binding filter (CN ≥ 6, pLDDT ≥ 0.70, geoRMSD ≤ 2.0)
   ↓ ~15-20% pass
7. Analysis with binding site pLDDT extraction
```

### 2.4 Expected Funnel

| Stage | Count | Rate | Notes |
|-------|-------|------|-------|
| Backbones | 400 | — | |
| After scout | ~240 | 60% | Scaffolding may have lower scout pass |
| Sequences | ~960 | 4 seq/bb | |
| After stability | ~480 | 50% | |
| Pass metal_binding | **60-80** | **12-17%** | Lower than de novo due to constraints |

**Note**: Expect lower pass rate than de novo (19%) because:
- Motif constraints reduce designable space
- Fixed coordinating residues may conflict with MPNN optimization
- Scaffold topology may not fit all TB coordination geometries

---

## Phase 3: Analysis & Comparison

### 3.1 De Novo vs Scaffolding Comparison

| Metric | De Novo (B1-B10) | Scaffolding (S1-S8) | Interpretation |
|--------|------------------|---------------------|----------------|
| Pass rate | 19.2% | TBD | Lower expected |
| Mean pTM | 0.88 | TBD | Similar expected |
| Mean pLDDT | 0.86 | TBD | Similar expected |
| **Binding site pLDDT** | TBD (extract) | TBD | **Key metric** |
| CN distribution | 53% CN=6, 36% CN=7 | TBD | May shift |
| AF3 Ca-citrate iptm | 0.84-0.92 | TBD | **Key metric** |
| AF3 Tb-citrate iptm | 0.82-0.91 | TBD | **Key metric** |

### 3.2 New Analysis: Binding Site Local Confidence

For all passing designs (de novo + scaffolding):

```python
analysis = {
    "design_id": "b01_seq_0180",
    "approach": "de_novo",  # or "scaffold_3c9h"
    "global_plddt": 0.864,
    "binding_site_residues": ["E31", "D53", "E93"],
    "binding_site_plddt": 0.72,  # Example — may be lower
    "binding_site_delta": -0.144,  # global - site
    "af3_ca_iptm": 0.89,
    "af3_tb_iptm": 0.90,
}
```

**Hypothesis test**:
- H0: binding_site_plddt(scaffold) = binding_site_plddt(de_novo)
- H1: binding_site_plddt(scaffold) > binding_site_plddt(de_novo)

### 3.3 AF3 Cross-Validation

Select top 10-20 scaffolding designs for AF3 validation:
- Primary: Ca-citrate (strict thresholds)
- Secondary: Tb-citrate (relaxed thresholds)

Compare AF3 interface confidence vs de novo designs.

---

## Implementation Checklist

### Setup
- [ ] Fetch 3C9H from RCSB
- [ ] Extract motif with `extract_metal_motif()`
- [ ] Verify MG → TB substitution works
- [ ] Inspect coordinating residues

### Exploration (Phase 1)
- [ ] Run E1-E4 configs (50 designs total)
- [ ] Analyze pass rates per config
- [ ] Extract binding site pLDDT for comparison
- [ ] Select best config for production

### Production (Phase 2)
- [ ] Lock best parameters
- [ ] Run S1-S4 batches (200 backbones)
- [ ] Mid-production checkpoint after S4
- [ ] Run S5-S8 batches (200 backbones)
- [ ] Full funnel analysis

### Analysis (Phase 3)
- [ ] Extract binding site pLDDT from de novo designs
- [ ] Compare de novo vs scaffolding metrics
- [ ] AF3 validation of top scaffolding designs
- [ ] Write comparison report

---

## Risk Mitigation

| Risk | Mitigation |
|------|------------|
| Low pass rate (<10%) | Relax cfg_scale, increase contig range, reduce fixed residues |
| MG → TB substitution fails | Manually edit motif PDB to place TB at MG position |
| Pocket too constrained | Try partial motif (only 2-3 coordinating residues) |
| No CN improvement | Accept that de novo may be superior for this target |

---

## Timeline Estimate

| Phase | Tasks | Duration |
|-------|-------|----------|
| Setup | Fetch 3C9H, extract motif | 1 hour |
| Exploration | E1-E4 (50 designs) | 2-3 hours |
| Production S1-S4 | 200 backbones, 800 sequences | 4-6 hours |
| Mid-checkpoint | Analysis, parameter adjustment | 1 hour |
| Production S5-S8 | 200 backbones, 800 sequences | 4-6 hours |
| Analysis | Comparison, AF3 validation | 2-3 hours |
| **Total** | | **~16-20 hours** |

---

## Success Criteria

### Primary (Must achieve)
1. ≥50 designs pass metal_tb_standard filter
2. Binding site pLDDT measurably extracted for all designs

### Secondary (Hypothesis test)
3. Scaffolding binding_site_pLDDT > de novo binding_site_pLDDT (p < 0.05)
4. ≥5 scaffolding designs pass AF3 Ca+Tb validation

### Exploratory
5. Any scaffolding design achieves CN ≥ 8
6. Any scaffolding design achieves pTM ≥ 0.92 (higher than de novo best)

---

## Files & Outputs

```
experiments/ln_citrate_scaffold/
├── scaffolding_3c9h/
│   ├── SCAFFOLDING_3C9H_PLAN.md      # This file
│   ├── inputs/
│   │   ├── 3c9h.pdb                  # Full PDB from RCSB
│   │   └── 3c9h_motif.pdb            # Extracted motif
│   ├── exploration/
│   │   ├── e1_configs.json
│   │   ├── e1_results.json
│   │   └── ...
│   ├── production/
│   │   ├── s01_results.json
│   │   ├── s01_promising/
│   │   └── ...
│   ├── analysis/
│   │   ├── binding_site_plddt_comparison.json
│   │   ├── de_novo_vs_scaffolding.json
│   │   └── af3_validation/
│   └── final/
│       ├── top_candidates/
│       └── PRODUCTION_REPORT.md
```

---

## Commands to Run

### Fetch 3C9H and extract motif
```bash
cd backend/serverless
python -c "
from scaffolding_workflow import ScaffoldingWorkflow
from metal_site_fetcher import _fetch_pdb_content

pdb = _fetch_pdb_content('3C9H')
with open('../../experiments/ln_citrate_scaffold/scaffolding_3c9h/inputs/3c9h.pdb', 'w') as f:
    f.write(pdb)

workflow = ScaffoldingWorkflow()
result = workflow.extract_metal_motif(pdb, metal='TB', ligand='CIT', contact_distance=5.0)
print(f'Coordinating residues: {result.coordinating_residues}')
print(f'Fixed atoms: {result.fixed_atoms}')
with open('../../experiments/ln_citrate_scaffold/scaffolding_3c9h/inputs/3c9h_motif.pdb', 'w') as f:
    f.write(result.motif_pdb)
"
```

### Run exploration batch
```bash
python production_runner.py explore-scaffold --config e1 --num-designs 15
```

### Run production batch
```bash
python production_runner.py produce-scaffold --batch 1 --num-backbones 50
```

---

## Appendix: 3C9H Structure Details

- **PDB ID**: 3C9H
- **Organism**: Agrobacterium tumefaciens (C58)
- **Function**: ABC transporter substrate-binding protein (citrate transporter)
- **Resolution**: 1.90 Å
- **Size**: 355 amino acids
- **Metal**: Mg²⁺ (will substitute → Tb³⁺)
- **Ligand**: Citrate (CIT)
- **Coordination**: Mg coordinated by citrate + protein Asp/Glu residues

### HSAB Compatibility Check
| Metal | HSAB Class | CN Range | Donor Preference |
|-------|------------|----------|------------------|
| Mg²⁺ | Hard | 4-6 | O > N |
| Tb³⁺ | Hard | 8-9 | O > N |

✓ Compatible for substitution (both hard, both O-preferring)

### Pocket Size Concern
- Mg ionic radius: 0.72 Å
- Tb ionic radius: 1.04 Å (CN=8)
- **44% larger** — may need pocket expansion
- Mitigation: RFD3 scaffold generation can expand pocket while preserving coordination geometry
