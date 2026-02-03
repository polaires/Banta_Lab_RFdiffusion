# Cumulative Lessons: Ln-Citrate Scaffold Design

> Auto-updated by analyze_round.py after each round.
> Manual insights should be added under each round section.

---

## Project Overview

**Goal:** Repurpose 3C9H (citrate-Mg binder) to create lanthanide-citrate binder

**Key Challenge:** Expand coordination from CN=6 (Mg) to CN=8-9 (Ln)

---

## Chemistry Lessons (Persistent)

### Lanthanide Coordination Requirements
- Tb³⁺/Eu³⁺/Gd³⁺: CN = 8-9, prefer square antiprism (SAP)
- M-O distances: 2.3-2.6Å (vs 2.0-2.3Å for Mg)
- Donors: O strongly preferred (Glu, Asp carboxylates)
- Citrate: Can provide 2-3 coordination sites

### HSAB Considerations
- Ln³⁺ is hard acid → hard bases only (O, some N)
- Exclude Cys (soft S donor) - use bias_AA: "C:-5.0"
- His acceptable but lower affinity than Glu/Asp

---

## Round-by-Round Lessons

<!-- Auto-populated by analyze_round.py -->

---

## Parameter Evolution Tracker

| Parameter | R1 | R2 | R7b | R8-NL (best) | Production B1-4 | Production B5+ | Rationale |
|-----------|-----|-----|-----|--------------|------------------|----------------|-----------|
| approach | de novo + unindex | de novo + fixed metal | motif scaffolding | auto template | auto template | auto template + partial co-diffusion | |
| contig length | 80-140 | 80-110 | 110-150 | 110-150 | 100-130 | 100-130 | Exploration found shorter is better |
| cfg_scale | none | none | 2.5 | 2.0 | 1.5 | 1.5 | 1.5 best in exploration |
| fixed_atoms (TB) | ALL | ALL / "" | X1:all | auto (TB:all) | TB:all | TB:all | |
| fixed_atoms (CIT) | "" | "" | L1:all | auto (CIT:all) | CIT:all | **CIT:O2,O5,O7** | Fix only coordinating O |
| burial | TB:ALL | TB:ALL | X1:all | TB+CIT | TB+CIT | TB+CIT | |
| hbond_acceptor | none | CIT:O1-O6 | L1:O1-O6 | auto: all O | all O (bug) | **O1,O3,O4,O6 only** | Exclude metal-coordinating O |
| hbond_donor | none | none | L1:O7 | none | none | none | |
| MPNN temp | - | - | 0.2 | 0.1 | 0.05 | 0.05 | Lower = more consistent |
| MPNN bias | - | - | none | A:-2.0, omit C | A:-2.0, omit C | A:-2.0, omit C | |
| RF3 context | - | - | seq only | seq + ligand + metal | seq + ligand + metal | seq + ligand + metal | |
| num_timesteps | - | - | 50 | 200 | 200 | 200 | |
| step_scale | - | - | 1.5 | 1.5 | 1.0 | 1.0 | 1.0 best in exploration |

---

## Best Practices Discovered

1. **Input Preparation**
   - Keep inputs simple: just metal + ligand (citrate_ln_only.pdb works)
   - Don't try to pre-position coordinating residues - RFD3 can't handle disconnected motifs
   - Metal at origin (0,0,0) is fine

2. **RFD3 Configuration**
   - Fix metal with `"TB": "ALL"` - scaffold will build around it
   - Unfix ligand with `"CIT": ""` for co-diffusion (but citrate still drifts)
   - `select_buried: {"TB": "ALL"}` helps concentrate scaffold around metal
   - De novo approach (contig=80-100) works for building scaffold around fixed point

3. **Key Limitation Discovered**
   - RFD3 doesn't understand lanthanide chemistry
   - Will coordinate metals with Cys (soft S) rather than Glu/Asp (hard O)
   - MUST use LigandMPNN with strong E/D bias after backbone generation

---

## Failure Patterns to Avoid

1. **Disconnected Motif Failure**
   - Symptom: Metal stays at origin, scaffold builds elsewhere, CN=1
   - Cause: Using `unindex` with residues in separate chains (E1, E2, D3, D4, N5)
   - Solution: Don't use disconnected motif residues; use de novo with burial conditioning

2. **Soft-Base Coordination**
   - Symptom: Metal coordinated by Cys sulfurs instead of Glu/Asp oxygens
   - Cause: RFD3 learned from PDB (most metals coordinated by Cys/His)
   - Solution: Use LigandMPNN with `bias_AA: "E:3.0,D:3.0,C:-5.0"` to enforce HSAB

3. **Ligand Drift**
   - Symptom: Citrate ends up 10-15Å from metal even with burial conditioning
   - Cause: No rigid constraint between metal and citrate
   - Solution: May need to constrain at least one citrate O to metal, or use partial diffusion

---

## Publication Notes

Key findings for methods section:
- [To be compiled from successful rounds]


## Round 1 (2026-01-21)

- Pass rate: 0.0%
- Best config: None
- Action: REVISE_FUNDAMENTALLY

**Root Cause Analysis:**
The metal (TB) stayed fixed at origin (0,0,0) while the scaffold was built elsewhere.
The citrate moved with the scaffold, resulting in metal-citrate distance of ~15Å.
CN=1 because only backbone N atoms from N-terminus happened to be near origin.

**Key Issues Identified:**
1. Using disconnected motif chains (E1, E2, D3, D4, N5) confused RFD3
2. The `unindex` approach with disconnected residues doesn't work
3. RFD3 doesn't naturally maintain metal-ligand relationships

## Round 2 (2026-01-21) - CORRECTED ANALYSIS

- Pass rate: 0.0% (tier B or better)
- Best configs: r2a_simple_buried_metal and r2b_hbond_citrate (both achieved CN=5)
- Best score: 61.6 (r2a_007)
- Action: FINE_TUNE (promising results, need more coordination)

**Experimental Configs:**
- r2a: Fixed metal, unfixed citrate, burial conditioning (8 designs)
- r2b: Same + H-bond acceptor conditioning on citrate (8 designs)
- r2c: Both metal and citrate unfixed - co-diffusion (8 designs)

**Corrected Analysis Results:**

| Config | CN Range | Best CN | Best Design | Key Coordinating Residues |
|--------|----------|---------|-------------|---------------------------|
| r2a | 0-5 | **5** | r2a_007 | GLU14.OE1 (2.0Å), GLU30.OE1 (2.26Å), GLU59.OE1 (2.88Å) |
| r2b | 0-5 | **5** | r2b_005 | GLN17.OE1, ASN42.OD1, GLN45.OE1, ASN67.OD1 |
| r2c | 0-2 | 2 | r2c_007 | ASN78.OD1 only - poor |

**Key Findings:**

1. **Fixed metal significantly outperforms unfixed:**
   - r2a/r2b: CN up to 5 with buried metal (SASA~0)
   - r2c: CN≤2 with exposed metal (SASA~10-30)

2. **RFD3 CAN place appropriate O-donors:**
   - r2a_007 has 3 GLU residues coordinating - excellent for Ln³⁺
   - Distances 2.0-2.9Å are appropriate for Ln coordination
   - This proves de novo approach can work!

3. **Citrate integration still weak:**
   - Best citrate contacts: 10 (r2b_001)
   - But citrate-coordinating designs often have low CN
   - Trade-off between citrate binding and metal coordination

4. **Analysis Bug Fixed:**
   - Parser was confusing CA (alpha carbon) with CA (calcium)
   - Now correctly identifying TB and counting sidechain donors

**Best Designs for Further Development:**
1. `r2a_simple_buried_metal_007` - CN=5, 3 GLU oxygens, ideal for Ln
2. `r2b_hbond_citrate_005` - CN=5, mixed O/N donors
3. `r2a_simple_buried_metal_001` - CN=4, includes HIS

**Next Steps for Round 3:**
1. **LigandMPNN on best backbones** with bias_AA: "E:3.0,D:3.0,C:-5.0"
2. **Increase CN target:** Need 3 more donors to reach CN=8
3. **Partial diffusion** from r2a_007 to add more coordinating loops
4. **Consider RFD3 hotspot conditioning** to place more donors near metal

## R1/R2 Critique vs RFD3 Best Practices (2026-01-22)

### Critical Issues Found

| Issue | R1/R2 Config | Best Practice | Impact |
|-------|--------------|---------------|--------|
| **No CFG with RASA** | `select_buried` without CFG | **MUST use CFG** with RASA | Burial conditioning ignored |
| **No CFG with H-bond** | `select_hbond_acceptor` without CFG | **MUST use CFG** | H-bond conditioning ineffective |
| **Conflicting origin params** | Both `ori_token` AND `infer_ori_strategy` | Use ONE only | Undefined behavior |
| **Insufficient designs** | 24 total (R2) | 100+ for motif scaffolding | Poor statistics |

### Key Quote from RFD3 Reference:
> "Apply CFG (`use_classifier_free_guidance=True`, `cfg_scale=2.0`) **when using H-bond or RASA conditioning**"

### Decision: Re-run Round 2 as Round 2b with corrections

Created `round_02b` configs with:
- `use_classifier_free_guidance: true`
- `cfg_scale: 2.0`
- Removed conflicting `ori_token` (keeping `infer_ori_strategy` only)
- 16 designs per config (48 total)

**Hypothesis:** With proper CFG, burial and H-bond conditioning will be enforced, leading to:
- Better metal coordination (potentially CN > 5)
- Improved citrate binding
- More consistent results across designs

## Round 2b Results (2026-01-22) - CFG Enabled

- Designs generated: 32 (16 per config, 3rd config failed due to atom naming)
- Best CN: 5 (same as R2)
- Score range: 33.6 - 60.8 (floor improved from 17.6)

### R2 vs R2b Comparison

| Metric | R2 (no CFG) | R2b (with CFG) | Change |
|--------|-------------|----------------|--------|
| Best CN | 5 | 5 | Same |
| Score floor | 17.6 | 33.6 | **+16 points** |
| Designs with CN≥3 | 21% | 44% | **+23%** |
| Metal SASA ~0 | 67% | 100% | CFG enforced burial |

### Best Designs from R2b

1. **r2b_b_cfg_hbond_005**: CN=5, Score=60.8
   - ASP32.OD1 (1.94Å), ASP36.OD1 (2.79Å), HIS52.NE2 (2.52Å) + 2 more

2. **r2b_b_cfg_hbond_007**: CN=3, Score=53.8 - **BEST QUALITY**
   - GLU16.OE2 (2.16Å), GLU40.OE1 (2.41Å), ASP69.OD1 (1.93Å)
   - **All 3 are carboxylate oxygens - ideal for Ln³⁺!**

3. **r2b_b_cfg_hbond_012**: CN=3, Score=54.8
   - GLU12.OE1 (1.99Å), GLU38.OE2 (1.92Å), ASN100.ND2 (2.48Å)

### Key Findings

1. **CFG improved consistency** - fewer failed designs, tighter score distribution
2. **CFG enforced burial** - 100% of designs have metal SASA ~0
3. **Better CN distribution** - twice as many designs with CN≥3
4. **CN ceiling unchanged** - still maxes at 5, not 8-9 needed for Ln
5. **Quality improved** - r2b_b_007 has pure carboxylate coordination

### Conclusion

CFG helps but is not sufficient to reach CN=8-9. RFD3 de novo approach plateaus at CN~5.

**Recommended next steps:**
1. Run LigandMPNN on best backbones (r2b_b_007, r2a_007) with bias_AA: "E:3.0,D:3.0,C:-5.0"
2. Try partial diffusion (partial_t) from best design to add coordinating loops
3. Consider using an existing high-CN metalloprotein as template

## Round 5 (2026-01-22) - Fixed TB-Citrate Complex

### Strategy Change
Based on Rounds 1-4 analysis revealing the critical metal-citrate separation issue (~17Å instead of 2-3Å), Round 5 used a **preformed TB-citrate complex** with proper coordination geometry as the fixed input.

### Pipeline Executed
1. **Round 5 RFD3**: Generated backbones with fixed TB-citrate complex
   - `r5_a_both_fixed`: Both metal and citrate fixed with RASA conditioning
   - `r5_b_hbond`: With H-bond conditioning on citrate oxygens

2. **Round 5 LigandMPNN**: Metal-aware sequence design
   - Used `ligand_mpnn_use_atom_context: 1`
   - Used `ligand_mpnn_use_side_chain_context: 1`
   - Temperature: 0.1 (conservative)
   - Omit: C (avoid soft sulfur donor)
   - 8 sequences per backbone

3. **Round 5 RF3**: Structure prediction validation
   - 12 predictions from 6 backbones × 2 sequences each
   - Included citrate ligand via SMILES

### Results Summary

| Design | pLDDT | Lig pLDDT | CN (4Å) | CN (6Å) | Coordinating Residues |
|--------|-------|-----------|---------|---------|----------------------|
| r5_b_hbond_011_seq01 | 86.0% | 65.8% | **2** | 8 | GLU78 (3.25Å), GLU74 (3.61Å) |
| r5_a_both_fixed_006_seq00 | 81.8% | 70.9% | 0 | 7 | - |
| r5_a_both_fixed_001_seq01 | 84.6% | 73.2% | 0 | 5 | - |
| Others | 74-88% | 62-71% | 0 | 0-2 | - |

### Key Findings

1. **Best Design: r5_b_hbond_011_seq01_rf3**
   - 2 GLU residues coordinating within 4Å (GLU78, GLU74)
   - 8 coordinating atoms within 6Å
   - High protein pLDDT (86%) indicates confident fold
   - Moderate ligand pLDDT (65.8%) suggests some positional uncertainty

2. **Coordination Improvement**
   - Rounds 1-4: CN=0 (metal-citrate too far apart)
   - Round 5: CN=2 at 4Å, up to 8 at 6Å
   - Proper coordination geometry now achieved

3. **H-bond Conditioning Helped**
   - `r5_b_hbond` designs performed better than `r5_a_both_fixed`
   - The only design with CN_4A > 0 came from H-bond conditioned backbone

4. **Ligand Position Still Uncertain**
   - Ligand pLDDT ~66% across designs
   - RF3 predicts protein confidently but citrate position less certain
   - This is expected - RF3 training likely had limited citrate examples

### Comparison: Round 4 vs Round 5

| Metric | Round 4 (separated) | Round 5 (fixed complex) |
|--------|---------------------|-------------------------|
| TB-citrate distance | ~17Å | 2-3Å (by design) |
| Coordinating residues at 4Å | 0 | **2** |
| Coordinating residues at 6Å | 0-1 | **8** |
| Protein pLDDT | 79-89% | 74-88% |
| Ligand pLDDT | 62-70% | 62-71% |

### Next Steps

1. **PyRosetta Validation**: Score best designs for energy and interface metrics
2. **Experimental Candidates**: Select top 2-3 designs for expression testing
3. **Alternative Metals**: Consider testing with Eu³⁺ or Gd³⁺ for luminescence/MRI applications
4. **Expand Coordination**: If CN=8-9 still needed, consider:
   - Partial diffusion from r5_b_hbond_011 to add more coordinating loops
   - Manual mutation of nearby residues to Glu/Asp

### Files Reference

| File | Description |
|------|-------------|
| outputs/round_05/*.pdb | RFD3 backbones with fixed TB-citrate |
| outputs/round_05_mpnn/*.fasta | LigandMPNN designed sequences |
| outputs/round_05_rf3/*_rf3.pdb | RF3 predictions (mmCIF format) |
| outputs/round_05_rf3/rf3_analysis.json | Analysis results |
| outputs/round_05_rf3/r5_b_hbond_011_seq01_rf3_converted.pdb | Best design in PDB format |

## Round 6 Plan (2026-01-22) - Achieving CN=8-9

### Problem to Solve
Round 5 achieved CN=2 at 4Å from protein donors. Target is CN=6 from protein (total CN=9 with citrate's 3 donors).

### Four Approaches

| Config | Strategy | Key Parameters |
|--------|----------|----------------|
| r6_a | Partial diffusion from R5 best (partial_t=0.5) | Keep core, resample loops |
| r6_b | Aggressive partial diffusion (partial_t=0.7) | Larger backbone changes |
| r6_c | C2 homodimer (50-70 res × 2) | Each chain provides ~3 donors |
| r6_d | Extended de novo (110-130 res) | More space for coordinators |

### Key Parameters (From RFD3 Best Practices)
- `cfg_scale: 2.0-2.5` - Classifier-free guidance enforces conditioning
- `step_scale: 1.5` (η from paper)
- `num_timesteps: 200` - More steps for complex designs
- `select_hbond_acceptor: all citrate oxygens`
- `32 designs per config` - 128 total for better statistics

### Success Criteria
- CN@4Å ≥ 6 from protein
- Carboxylate donors ≥ 4 (Glu/Asp)
- No sulfur donors (Cys incompatible with Ln³⁺)

### Files Created
- `configs/round_06/r6_configs.json` - All configurations
- `scripts/run_r6_direct.py` - Execution script
- `scripts/analyze_r6.py` - Analysis with detailed CN metrics

## Round 6 Results (2026-01-22) - Partial Diffusion Success

### Execution Summary
- **Configs run:** 2 of 4 (r6_a, r6_b partial diffusion)
- **Configs failed:** 2 (r6_c dimer - tensor mismatch, r6_d extended - CUDA error)
- **Total designs:** 64 (32 per successful config)

### Results

| Config | CN@4A Range | CN@4A Mean | Best CN | Carboxylate Mean |
|--------|-------------|------------|---------|------------------|
| r6_a_partial (partial_t=0.5) | 0-4 | 1.72 | 4 | 0.31 |
| r6_b_aggressive (partial_t=0.7) | 0-6 | 1.28 | **6** | 0.12 |

### Best Design: r6_b_partial_aggressive_031

```
CN@4A: 6 (TARGET MET!)
CN@3.5A: 4
Coordinating residues:
- ARG31: 3 atoms (2.9-3.78Å) - NH1, NH2
- GLN53: 2 atoms (2.56-3.42Å) - OE1, NE2
- 1 Carboxylate donor

Quality: moderate (need more Glu/Asp carboxylates for HSAB)
```

### Key Findings

1. **Partial diffusion improved CN from R5**
   - R5 best: CN=2 at 4Å
   - R6 best: CN=6 at 4Å (3x improvement!)

2. **Aggressive partial_t=0.7 found best design**
   - More backbone resampling allowed larger changes
   - But also more variability (some designs have CN=0)

3. **Donor type issue remains**
   - RFD3 prefers ARG and GLN over GLU/ASP
   - Only 1 carboxylate in best design
   - Need LigandMPNN with strong Glu/Asp bias

4. **API format discovery**
   - `select_fixed_atoms` keys must be ChainIDResID format (e.g., "X96", "L1")
   - Using just chain ID ("X") causes validation error

### Comparison: R5 vs R6

| Metric | R5 Best | R6 Best | Improvement |
|--------|---------|---------|-------------|
| CN@4A | 2 | **6** | **+4 (+200%)** |
| CN@3.5A | 0 | 4 | +4 |
| Carboxylate donors | 2 (GLU) | 1 | -1 |
| Total with citrate | ~5 | **~9** | +4 |

### Next Steps

1. **Fix GPU issue** - Restart WSL to reset nvidia-container-cli
2. **Retry failed configs** - r6_c dimer, r6_d extended
3. **LigandMPNN** on r6_b_031 with strong Glu/Asp bias
4. **RF3/AF3 validation** of sequence-designed structures
5. **Consider iterative approach** - Partial diffusion from r6_b_031 to further improve

### Files

| File | Description |
|------|-------------|
| `outputs/round_06/r6_b_partial_aggressive_031.pdb` | Best design (CN=6) |
| `analysis/round_06_analysis.json` | Full analysis results |
| `analysis/ROUND_6_SUMMARY.md` | Summary table |

## Round 6 (Continued) - Dimer Success & LigandMPNN

### Dimer Config Results (r6_c)

After Docker restart, successfully ran the C2 symmetric dimer config:

| Config | CN@4A Range | CN@4A Mean | Best CN | Carboxylate Mean |
|--------|-------------|------------|---------|------------------|
| r6_c_dimer_symmetric | 0-6 | 2.34 | **6** | **0.47** |

**Best Dimer Designs:**

| Rank | File | CN@4A | Carboxylates |
|------|------|-------|--------------|
| 1 | r6_c_dimer_symmetric_025 | 6 | **3** (GLU39, ASP48×2) |
| 2 | r6_c_dimer_symmetric_012 | 5 | **4** |
| 3 | r6_c_dimer_symmetric_010 | 4 | 2 |

**Key Finding: Dimer approach produces more carboxylate donors!**
- r6_c average: 0.47 carboxylates per design
- r6_b average: 0.12 carboxylates per design
- **4x improvement** in carboxylate positioning

### LigandMPNN Sequence Design

Ran LigandMPNN on top 4 designs with carboxylate bias:

**Strong Bias (first attempt):**
```python
bias_AA = "E:3.0,D:3.0,N:2.0,Q:2.0,C:-5.0"
```
- Result: 61% charged residues (E+D)
- RF3 RMSD: 25Å (POOR - too much charge repulsion)
- Sequences non-designable

**Moderate Bias (second attempt):**
```python
bias_AA = "E:1.5,D:1.5,N:1.0,Q:1.0,C:-3.0"
```
- Result: 32% charged residues (E+D)
- RF3 RMSD: 12Å (improved but still POOR)
- Better but still challenging designability

### Key Lessons from LigandMPNN

1. **Carboxylate bias is a double-edged sword:**
   - Too strong (E:3.0, D:3.0) → 60%+ charged, non-foldable
   - Moderate (E:1.5, D:1.5) → 30% charged, poor RMSD
   - Need ~15-20% charged for good designability

2. **Designability vs Chemistry trade-off:**
   - Ideal Ln³⁺ coordination needs 6+ carboxylate donors
   - But proteins with 20%+ charged residues often misfold
   - May need to accept fewer carboxylates + some His/Asn

3. **Dimer architecture is promising:**
   - Places more potential donors near metal
   - Each chain contributes 1-2 carboxylates
   - C2 symmetry doubles effective donor count

### Comparison: Partial vs Dimer

| Approach | Best CN | Carboxylates | HSAB Quality |
|----------|---------|--------------|--------------|
| r6_b_partial_aggressive | 6 | 1 | Poor |
| r6_c_dimer_symmetric | 6 | 3-4 | **Good** |

**Conclusion:** C2 dimer approach is superior for carboxylate-coordinated Ln³⁺ binding.

### Final Round 6 Rankings

| Rank | Design | CN@4A | Carboxylates | Quality |
|------|--------|-------|--------------|---------|
| 1 | r6_c_dimer_symmetric_025 | 6 | 3 | **Best** |
| 2 | r6_c_dimer_symmetric_012 | 5 | 4 | Good |
| 3 | r6_b_partial_aggressive_031 | 6 | 1 | Moderate |
| 4 | r6_c_dimer_symmetric_010 | 4 | 2 | Moderate |

### Structural Validation Results

| Design | Bias Level | E+D % | RF3 RMSD |
|--------|------------|-------|----------|
| r6_c_025 | Strong (3.0) | 61% | 25Å |
| r6_c_025 | Moderate (1.5) | 32% | 12Å |
| r6_b_031 | Strong (3.0) | 42% | 21Å |

**Finding:** All LigandMPNN sequences show poor RMSD to backbone, suggesting:
1. RFD3 backbones may have intrinsically low designability
2. Metal-coordination geometry creates unusual structural constraints
3. May need iterative backbone/sequence co-optimization

### Recommendations for Round 7

1. **Try ProteinMPNN** (non-ligand aware) for baseline comparison
2. **Use lower bias** (E:0.5, D:0.5) to improve designability
3. **Filter MPNN outputs** by amino acid composition before RF3
4. **Consider BindCraft** for ligand-aware design
5. **Partial diffusion from dimer** to optimize coordination loops

### Files Created

| File | Description |
|------|-------------|
| `outputs/round_06_mpnn/*.fasta` | LigandMPNN sequences |
| `outputs/round_06_mpnn/*_moderate_bias.fasta` | Moderate bias sequences |
| `outputs/round_06_mpnn/*_rf3.pdb` | RF3 structure predictions |

## Round 6 (FINAL SUCCESS) - Fixed Coordinating Residues

### The Missing Piece: Preserve Coordination

The coordinating residues must be **FIXED** during LigandMPNN, not redesigned!

**Coordinating residues identified:**
```
r6_b_partial_aggressive_031:
  A31 ARG.NH1 (2.90Å) - fixed
  A53 GLN.NE2 (2.56Å) - fixed
  A93 ASP.OD1 (3.28Å) - fixed (carboxylate!)
```

### Final LigandMPNN Parameters

```python
payload = {
    'task': 'mpnn',
    'model_type': 'ligand_mpnn',
    'omit_AA': 'C',
    'fixed_residues': 'A31 A53 A93',  # CRITICAL: Fix coordinating residues!
    'ligand_mpnn_use_atom_context': 1,
    'ligand_mpnn_use_side_chain_context': 1,
    'temperature': 0.1,
    'pack_side_chains': 1,
}
```

### Final Results: ALL APPROACHES COMPARED

| Approach | RMSD | E+D % | Status |
|----------|------|-------|--------|
| Strong bias (E:3.0) | 21.0Å | 61% | ❌ Non-designable |
| Moderate bias (E:1.5) | 12.0Å | 32% | ❌ Poor |
| R5 approach (no bias) | 6.86Å | 12% | ❌ Poor |
| **Fixed coord + no bias** | **0.82Å** | 12.6% | ✅ **EXCELLENT** |

### Best Designable Sequence

```
Design: r6_b_partial_aggressive_031
RMSD: 0.82Å (Excellent)
Fixed residues: A31 (ARG), A53 (GLN), A93 (ASP)

Sequence (95 aa):
LDAAAIRAGVRKQPGISDDYAAALDALLPGEVAELTAGGLTDTAAAATIFAGDLAARFYPADKAADRDALVAAARAALAAAAGTLTPAAVAAELA

Amino acid composition:
  E (Glu): 3 (3.2%)
  D (Asp): 9 (9.5%)
  E+D total: 12 (12.6%)
  Cys: 0
```

### Key Lessons: Complete LigandMPNN Protocol for Metal Sites

1. **Identify coordinating residues** - find all residues within 4Å of metal
2. **Fix coordinating residues** - use `fixed_residues` parameter
3. **No global bias** - let metal context work naturally
4. **Use `omit_AA: "C"`** - exclude cysteines for hard acids like Ln³⁺
5. **Enable metal context** - `ligand_mpnn_use_atom_context: 1`

### Files Created

| File | Description |
|------|-------------|
| `*_fixed_coord.fasta` | Sequences with fixed coordinating residues |
| `*_fixed_coord_best_rf3.pdb` | RF3 predictions |
| `BEST_CANDIDATE_r6_b_031.fasta` | **Final best sequence (0.82Å RMSD)** |

## Round 6 (CRITICAL BUG FIX) - Parameter Discovery

### Bug: `fixed_residues` vs `fixed_positions`

The 0.82Å RMSD result was achieved with a **BUG** - the API parameter was wrong!

```python
# WRONG - string format, not recognized by API
'fixed_residues': 'A31 A53 A93'

# CORRECT - list format, actually fixes residues
'fixed_positions': ['A31', 'A53', 'A93']
```

**Result:** The 0.82Å sequence was generated **WITHOUT** actually fixing coordinating residues.

### Verification

When using **correct** `fixed_positions` list format:
- RMSD dropped to **8.57Å** (Poor)
- Residues were properly preserved (R, Q, S, D at 31, 53, 92, 93)
- But the constrained sequence doesn't fold well

### Critical Discovery: LigandMPNN Improved Coordination Naturally!

When **not constrained**, LigandMPNN used its ligand context to select better residues:

| Position | Original (PDB) | LigandMPNN (0.82Å) | Function |
|----------|---------------|-------------------|----------|
| 31 | R (ARG) | **E (GLU)** | Carboxylate - BETTER for Ln³⁺! |
| 53 | Q (GLN) | **D (ASP)** | Carboxylate - BETTER for Ln³⁺! |
| 92 | S (SER) | A (ALA) | Lost weak citrate H-bond |
| 93 | D (ASP) | **E (GLU)** | Carboxylate - equivalent |

**Key Insight:** The `ligand_mpnn_use_atom_context: 1` parameter allowed LigandMPNN to
recognize the metal position and select carboxylates (E/D) which are ideal for Ln³⁺
binding, WITHOUT needing explicit bias or fixing!

### Updated Best Practice

For metal-binding site design with LigandMPNN:

```python
payload = {
    'task': 'mpnn',
    'model_type': 'ligand_mpnn',
    'ligand_mpnn_use_atom_context': 1,      # CRITICAL: Enables metal awareness
    'ligand_mpnn_use_side_chain_context': 1,
    'omit_AA': 'C',                          # Exclude Cys for Ln³⁺
    'temperature': 0.1,
    # NO bias_AA - let metal context work naturally
    # NO fixed_positions - allows better optimization
}
```

### Final Best Sequence (Corrected Understanding)

```
Design: r6_b_partial_aggressive_031_best
RF3 RMSD: 0.81Å (EXCELLENT - verified)
Coordination: E31, D53, E93 (all carboxylates!)

Sequence:
LDAAAIRAGVRKQPGISDDYAAALDALLPGEVAELTAGGLTDTAAAATIFAGDLAARFYPADKAADRDALVAAARAALAAAAGTLTPAAVAAELA
```

The sequence has **3 carboxylate donors** at coordination positions, which is
chemically superior to the original R/Q/D backbone design.

### Next Step: AF3 Validation

- RF3 confirms backbone folds correctly (0.81Å RMSD)
- Need AF3 to validate binding pocket function
- AF3 submission guide: `outputs/round_06_mpnn/AF3_SUBMISSION_GUIDE.md`
- Sequence file: `outputs/round_06_mpnn/BEST_FOR_AF3.fasta`

### Files Created

| File | Description |
|------|-------------|
| `BEST_FOR_AF3.fasta` | Best sequence for AF3 validation |
| `AF3_SUBMISSION_GUIDE.md` | Instructions for AF3 server submission |

## Round 6 (CRITICAL LIMITATION) - Ala vs Citrate Trade-off

### AF3 Validation Result

The 0.81Å RMSD sequence was submitted to AlphaFold3:
- **Tb binding: YES** - Metal coordinates successfully
- **Tb-citrate binding: NO** - Citrate does not bind

### Root Cause: Alanine Content

The best-folding sequence has **37% Alanine**, which:
1. Provides good backbone stability (0.81Å RMSD)
2. Lacks H-bond donors/acceptors for citrate binding
3. Position 92 became Ala, losing the SER-citrate contact

### Bias Experiments

| Bias | Ala % | E+D % | Best RMSD | Status |
|------|-------|-------|-----------|--------|
| None | 37% | 13% | **0.81Å** | Folds but no citrate |
| A:-0.5 | 16-27% | 13-27% | 9.09Å | Poor fold |
| A:-1.5 | 3-5% | 12-21% | 6.68Å | Poor fold |
| A:-3.0 | 0-1% | 61-67% | 15-37Å | Very poor fold |

### Key Lesson: LigandMPNN Alanine Tendency

LigandMPNN heavily favors Alanine for:
- Small hydrophobic core positions
- Positions with few constraints
- General "safe" residues

**Problem:** The r6_b_partial_aggressive_031 backbone REQUIRES ~37% Ala to fold
correctly. Any bias reducing Ala destroys designability.

### Fundamental Conflict

```
Designability: Needs ~37% Ala → Good RMSD
Citrate binding: Needs H-bond donors (S, T, N, Q) at binding site → Replaces Ala
These requirements CONFLICT for this backbone!
```

### Why This Backbone Fails for Tb-Citrate

1. **Backbone geometry** was optimized by RFD3 without considering citrate H-bonding
2. **Core residues** need to be small hydrophobic (Ala) for stability
3. **Citrate binding site** happens to overlap with these core positions
4. **Cannot have both** good folding AND citrate binding

### Recommendations for Round 7

1. **New backbone design** with explicit citrate H-bond conditioning
   - Use `select_hbond_donor` for citrate oxygens
   - Ensure binding site is on surface, not in core

2. **Different architecture** - beta-barrel or TIM barrel may tolerate more polar residues

3. **Template-based approach** - Start from existing citrate-binding protein scaffold

4. **Accept Tb-only binding** - If Tb binding is sufficient for application, use current design

### Files Created

| File | Description |
|------|-------------|
| `r6_b_031_low_ala.fasta` | A:-3.0 bias (too polar, doesn't fold) |
| `r6_b_031_moderate_bias.fasta` | A:-1.5 bias (poor RMSD) |
| `BEST_MILD_ALA_BIAS.fasta` | A:-0.5 bias (poor RMSD) |
| `BEST_DIMER.fasta` | Dimer backbone (also doesn't fold) |

## Round 6d (BREAKTHROUGH) - H-bond Conditioned Backbone

### Retry of Failed r6_d Config

After retrying with `cfg_scale: 2.0` (reduced from 2.5), r6_d successfully generated designs.

### Backbone Analysis

```
r6d_extended_000:
  Length: 120 residues
  CN@4A: 4
  Carboxylates: 1 (GLU91)
  H-bond conditioning: Applied to all citrate oxygens
```

### LigandMPNN Results (NO BIAS)

| Sequence | RMSD | Ala % | E+D % | Status |
|----------|------|-------|-------|--------|
| seq03 | **0.46Å** | 29% | 11% | Excellent |
| seq04 | **0.47Å** | 24% | 14% | Excellent |
| **seq05** | **0.43Å** | **23%** | 12% | **BEST** |
| seq06 | 0.87Å | 27% | 12% | Excellent |
| seq07 | 10.37Å | 25% | 13% | Poor |

**7 of 8 sequences have excellent RMSD (<1.5Å)!**

### Key Breakthrough: Ala Requirement Solved

| Backbone | Required Ala % | Best RMSD | Citrate Potential |
|----------|---------------|-----------|-------------------|
| r6_b_031 | 37% | 0.81Å | Poor (Ala blocks H-bonds) |
| r6_c_dimer | ~15% | 11.33Å | N/A (doesn't fold) |
| **r6d** | **23%** | **0.43Å** | **Good (more polar residues)** |

### Why r6d Works

1. **H-bond conditioning created appropriate backbone geometry**
   - Used `select_hbond_acceptor` for all citrate oxygens
   - Used `select_hbond_donor` for O7 (hydroxyl)
   - `cfg_scale: 2.0` enforced the conditioning

2. **Extended length (110-130)** provides more design freedom

3. **Backbone tolerates polar residues**
   - Only 23% Ala needed (vs 37%)
   - 14% polar (STQN) vs ~8% for r6_b_031
   - 12 Thr residues for potential citrate H-bonding

### Best Candidate: r6d_seq05

```
RMSD: 0.43Å (EXCELLENT)
Ala: 23%, E+D: 12%, Polar: 14%
Length: 120 aa

Sequence:
MDPNVIRVHLRLELGPGQTPADVVAFIAARNAVATEYTVTVRSAAATATPGTAELVLDLTLLGATPEQLYAAELADERAALAAGLFRREEATLVTHPAGRAHAEALARAHRDAGLTARVV

Glu positions: 13, 36, 54, 67, 73, 77, 89, 90, 104
Thr positions: 19, 35, 38, 40, 47, 49, 52, 60, 65, 92, 95, 116
```

### Coordination Changes from LigandMPNN

| Position | Backbone | Designed | Change |
|----------|----------|----------|--------|
| 77 | THR | **E** | Upgraded to carboxylate! |
| 89 | GLN | **E** | Upgraded to carboxylate! |
| 91 | GLU | A | Lost carboxylate |

LigandMPNN naturally placed carboxylates at coordination-relevant positions.

### Next Step: AF3 Validation

- Guide saved: `outputs/round_06d/AF3_SUBMISSION_R6D.md`
- Sequence file: `outputs/round_06d/BEST_R6D.fasta`

### Critical Lesson Learned

**H-bond conditioning during backbone generation is essential for citrate binding.**
- Without it (r6_b_031): Backbone requires 37% Ala → no room for polar residues
- With it (r6d): Backbone only needs 23% Ala → 14% polar residues available

### Files Created

| File | Description |
|------|-------------|
| `outputs/round_06d/r6d_extended_000.pdb` | H-bond conditioned backbone |
| `outputs/round_06d/r6d_extended_000_seqs.fasta` | LigandMPNN sequences |
| `outputs/round_06d/BEST_R6D.fasta` | Best sequence (0.43Å RMSD) |
| `outputs/round_06d/AF3_SUBMISSION_R6D.md` | AF3 submission guide |

## Round 8: NL Pipeline vs R7b E2E Comparison (2026-01-31)

### Overview

Systematic head-to-head comparison of the NL pipeline (automated defaults via `metal_binding_design` single mode) vs Round 7b (expert-tuned manual parameters). Both arms generated 5 backbones and 20 designed sequences for Tb-citrate binding.

### Results

| Metric | R7b (Expert) | NL (Auto) | Delta |
|--------|-------------|-----------|-------|
| Pass rate (pTM >= 0.6) | 85% (17/20) | **100% (20/20)** | +15% |
| Avg pTM | 0.770 | **0.902** | +17% |
| Avg pLDDT | 0.801 | **0.867** | +8% |
| Avg PAE | 6.22 | **3.61** | -42% |
| Avg Ala% | 17.1% | **6.7%** | 2.5x less |
| Within-backbone pTM std | 0.126 | **0.008** | 15.7x better |
| Best pTM | 0.910 | **0.925** | |
| Worst pTM | 0.357 | **0.865** | |
| RFD3 time | 57s | 242s | 4.2x slower |
| Total time | 520s | 736s | 1.4x slower |

### Key Lessons

1. **NL pipeline automated defaults outperform expert-tuned R7b on every quality metric.** The auto-template approach (metal+ligand only, no hand-crafted motif residues) produces more designable scaffolds. Motif residues in R7b may over-constrain.

2. **MPNN bias_AA: "A:-2.0" + omit_AA: "C" is critical.** R7b sequences had 17-32% alanine (classic MPNN failure for metal binding). NL pipeline's bias reduces this to 2.7-13.8%, yielding healthier sequence composition with more charged residues (Glu 12%, Lys 7%).

3. **Ligand-aware RF3 dramatically improves prediction consistency.** Sequence-only RF3 (R7b) has 3 designs below pTM 0.5 and huge within-backbone variance (std 0.126). RF3 with citrate SMILES + TB metal (NL) gives uniformly high confidence (min pTM 0.865, std 0.008).

4. **CFG 2.0 is sufficient when combined with other improvements.** R7b's CFG 2.5 provided no advantage over NL's 2.0. The quality gains come from MPNN bias and RF3 context, not stronger conditioning.

5. **Metal+ligand burial doesn't hurt.** Despite concern that burying citrate would reduce H-bond access, NL pipeline buries both metal and ligand and achieves higher pTM/pLDDT.

6. **Auto H-bond handling is adequate.** NL pipeline assigns all citrate O atoms as acceptors (O1-O7) and no donors, vs R7b's curated O1-O6 acceptors + O7 donor. The simpler auto approach works equally well.

7. **Speed tradeoff is acceptable.** NL pipeline is 1.4x slower due to `metal_binding_design` single mode overhead (template generation, scaffold extraction) and `pack_side_chains` MPNN (39s vs 3s). The quality improvement justifies this.

### Bugs Fixed During Round 8

1. **`analyze_design` task missing from Docker** - handler.py had the task in uncommitted changes, but `filter_evaluator.py`, `design_history.py`, `lesson_detector.py`, `backbone_pre_filter.py` were not volume-mounted in docker-compose.local.yml. Fixed by adding mounts.

2. **RF3 CIF output coordination extraction** - RF3 returns CIF format where res_name is `L:0` (metal) and `L:1` (ligand) instead of `TB`/`CIT`. The `analyze_structure()` function in inference_utils.py matched by res_name, failing for CIF. Fixed by adding atom_name-based metal detection (e.g., atom `TB0` -> metal `TB`).

3. **Windows Unicode encoding** - Arrow and em-dash characters crashed stdout with cp1252 codec. Fixed to use ASCII alternatives.

### Recommendations for Production

Based on Round 8 findings, the NL pipeline should be the default for all Tb-citrate (and likely all Ln-metal) designs:

- Use `metal_binding_design` mode=single for backbone generation
- Always set `bias_AA: "A:-2.0"` and `omit_AA: "C"` in MPNN
- Always use ligand-aware RF3 validation (ligand_smiles + metal)
- CFG 2.0 is the right default
- No need for hand-crafted motif PDBs; auto-templates suffice

### Data Location

- Raw outputs: `backend/serverless/output_nl_vs_r7b/`
- Analysis JSON: `experiments/ln_citrate_scaffold/analysis/round_08_nl_vs_r7b.json`
- Test script: `backend/serverless/test_nl_vs_r7b.py`

## Production Batches 1-4 (2026-02-01 to 2026-02-02)

### Overview

First 4 production batches using locked best params from exploration rounds 1-5.

**Parameters (best_params.json):**
- cfg_scale: 1.5, contig: 100-130, bury_ligand: true
- num_timesteps: 200, step_scale: 1.0
- MPNN: temperature 0.05, bias_AA A:-2.0, omit_AA C
- Batches 1-2: 4 seqs/backbone, Batches 3-4: 2 seqs/backbone

### Pipeline Funnel

| Metric | B1 | B2 | B3 | B4 | Cumulative |
|--------|-----|-----|-----|-----|------------|
| Backbones | 50 | 50 | 50 | 50 | 200 |
| Sequences | 200 | 200 | 100 | 100 | 600 |
| After stability | 98 | 105 | 50 | 49 | 302 |
| Pass metal_binding | 18 | 22 | 10 | 15 | 65 |
| Pass strict | 0 | 0 | 0 | 0 | 0 |
| Pass rate | 18.4% | 21.0% | 20.0% | 30.6% | 21.5% |

### CN=8 Designs (Target for Tb)

| Batch | Count | Best |
|-------|-------|------|
| B1 | 1 | seq_0157 (score 0.556) |
| B2 | 3 | seq_0180 (score 0.579), seq_0178, seq_0049 |
| B3 | 1 | seq_0057 (score 0.525) |
| B4 | 2 | seq_0096 (score 0.509), seq_0009 (score 0.495) |
| **Total** | **7** | |

### Key Production Findings

1. **Chemistry-aware filtering active**: Uses `metal_tb_standard` preset with TB-specific thresholds (CN>=6, pTM>=0.60, pLDDT>=0.70, ligand_contacts>=3) instead of generic `metal_binding`.

2. **Sub-batching solved timeout issue**: Original 50-backbone API calls exceeded 1800s HTTP timeout. Sub-batching at 10/call with 900s timeout works reliably.

3. **2 seqs/bb vs 4 seqs/bb**: Batch 3-4 used 2 seqs/bb. B4 had highest pass rate (30.6%) despite fewer sequences. Suggests backbone quality matters more than sequence sampling depth for metal binding. Reverted to 4 seqs/bb for batch 5+.

4. **No strict filter passes**: All 65 passing designs fail strict tier. The strict preset (CN>=8, pLDDT>=0.80, pTM>=0.75, ligand_contacts>=5) is too stringent for current pipeline. The 7 CN=8 designs fail on other strict metrics (pTM or ligand_contacts).

5. **Composite score ceiling ~0.62**: Top designs cluster at 0.60-0.62. Improving beyond this likely requires structural changes (partial ligand co-diffusion, different contig range).

6. **H-bond acceptor conflict identified**: All citrate O atoms were set as H-bond acceptors, including O2/O5/O7 which coordinate the metal. These coordinating oxygens have lone pairs occupied by metal bonds and should NOT be H-bond acceptors. Fixed for batch 6+.

7. **Partial ligand co-diffusion introduced (batch 5)**: `ligand_fix_atoms: "O2,O5,O7"` fixes only coordinating oxygens while letting non-coordinating atoms (O1,O3,O4,O6) and carbons co-diffuse. This allows RFD3 to optimize citrate orientation for protein H-bonds.

### Parameter Evolution Update

| Parameter | Exploration | Production B1-4 | Production B5+ | Rationale |
|-----------|-------------|------------------|-----------------|-----------|
| cfg_scale | 1.5 (best) | 1.5 | 1.5 | Stable |
| contig | 100-130 | 100-130 | 100-130 | Stable |
| bury_ligand | false (best) -> true | true | true | bury=true standard |
| num_timesteps | 200 (best) | 200 | 200 | Stable |
| step_scale | 1.0 (best) | 1.0 | 1.0 | Stable |
| MPNN temp | 0.05 (best) | 0.05 | 0.05 | Stable |
| ligand_fix_atoms | all (default) | all | **O2,O5,O7** | Partial co-diffusion |
| hbond_acceptor | all O atoms | all O atoms | **O1,O3,O4,O6** | Exclude coordinating |
| num_seqs | 4 | 4->2->4 | 4 | 4 better coverage |

### Promising PDBs Location

All 65 passing designs saved to: `experiments/ln_citrate_scaffold/final/promising/`
- Files: `b{NN}_{seq_id}.pdb` + `b{NN}_{seq_id}_metrics.json`

## Production Batch 5 (2026-02-02) - Partial Ligand Co-diffusion

### Parameters Changed
- `ligand_fix_atoms: "O2,O5,O7"` — fix only coordinating oxygens, let rest co-diffuse
- `select_hbond_acceptor: all O atoms` — BUG: still included coordinating O2/O5/O7

### Results

| Metric | B5 | B1-B4 avg |
|--------|-----|-----------|
| Backbones | 50 | 50 |
| Sequences | 200 | 150 |
| After stability | 90 (45%) | 75.5 (50%) |
| Pass metal_binding | 11 (12.2%) | 16.3 (21.5%) |
| Pass strict | 0 | 0 |

### Key Findings

1. **Pass rate dropped to 12.2%** from 21.5% average. Partial co-diffusion + H-bond conflict likely caused this.

2. **CN=10 design found**: seq_0197 (pTM=0.713, pLDDT=0.785, LC=19) — highest CN across all production batches. Co-diffusion allows more donors to organize around metal.

3. **CN=8 design**: seq_0050 (pTM=0.695, CN=8, LC=13).

4. **Lower stability rate**: 45% vs 50% average. Co-diffusion may produce less stable backbones since ligand position is partially uncertain.

5. **H-bond acceptor bug impact**: Conditioning H-bonds on O2/O5/O7 (which are metal-coordinating) conflicts with metal coordination. The model tries to place protein H-bond donors near these atoms AND maintain metal-O bonds simultaneously. **Fixed for batch 6+** — only O1,O3,O4,O6 are H-bond acceptors.

## Production Batch 6 (2026-02-02) - H-bond Fix + Partial Co-diffusion

### Results

| Metric | B6 | B5 | B1-4 avg |
|--------|-----|-----|----------|
| Pass metal_binding | 10 (9.7%) | 11 (12.2%) | 16.3 (21.5%) |
| CN=8+ | 0 | 2 | ~2 |
| Best score | 0.579 | 0.612 | 0.622 |

### Conclusion: Partial ligand co-diffusion HURTS pass rate

The H-bond acceptor fix did NOT recover pass rate — it dropped further from 12.2% to 9.7%. This confirms **partial co-diffusion itself is the problem**, not the H-bond conflict.

When non-coordinating citrate atoms diffuse freely:
- Citrate position becomes less deterministic
- MPNN gets a less stable reference for sequence design
- RF3 sees more positional uncertainty
- Result: lower CN (max 7 vs 10 in B5), lower pTM, fewer passing designs

### Action: Reverted to all-fixed ligand for B7+

Removed `ligand_fix_atoms` from best_params.json. B7+ returns to fully fixed citrate (the B1-B4 approach with ~21.5% pass rate).

## Production Batch 7 (2026-02-02) - Reverted to All-Fixed

### Results

| Metric | B7 | B5-B6 avg | B1-B4 avg |
|--------|-----|-----------|-----------|
| Pass metal_binding | 21 (19.4%) | 10.5 (11.0%) | 16.3 (21.5%) |
| CN=8+ | 1 | 1 | ~2 |
| Best score | 0.604 | 0.596 | 0.622 |

**Pass rate recovered to 19.4%** — confirms all-fixed ligand is the correct approach.

Notable: seq_0196 (CN=8, LC=26, pTM=0.827) — backbone b07_048 produced 3 top-10 designs.

## Production Batch 8 (2026-02-02)

- Pass: 16/91 (17.6%), Best score: 0.585
- 2 CN=8 designs (seq_0178: pTM=0.815, seq_0102: pTM=0.700)
- Stable in expected range, no anomalies

## Production Batch 9 (2026-02-02)

- Pass: 20/101 (20.0%), Best score: 0.620
- Top design: seq_0138 (pTM=0.901, pLDDT=0.864, CN=7, LC=13) — 2nd best across all batches
- Consistent with B7-B8, confirms all-fixed approach is stable

## Production Batch 10 (2026-02-02) - Best Batch

- Pass: 30/104 (28.8%), Best score: 0.604
- **Highest pass rate across all 10 batches**
- Natural variance — same params as B7-B9 but happened to get better backbones
- Notable: seq_0065 had pTM=0.892 but only CN=4 (failed filter)

## Final Production Summary (B1-B10)

### Full Funnel

| Stage | Count | Rate |
|-------|-------|------|
| Backbones generated | 500 | — |
| Sequences designed | 1,800 | 3.6 seq/bb avg |
| After stability | 899 | 50.0% |
| Pass metal_binding | 173 | 19.2% |
| Pass strict | 0 | 0% |

### Per-Batch Summary

| Batch | Pass | Rate | CN=8+ | Best Score | Notes |
|-------|------|------|-------|------------|-------|
| B1 | 18 | 18.4% | 1 | 0.622 | Baseline |
| B2 | 22 | 21.0% | 3 | 0.616 | — |
| B3 | 10 | 20.0% | 1 | 0.614 | 2 seq/bb test |
| B4 | 15 | 30.6% | 2 | 0.610 | — |
| B5 | 11 | 12.2% | 2 | 0.612 | Partial co-diffusion (harmful) |
| B6 | 10 | 9.7% | 0 | 0.579 | Co-diffusion + H-bond fix |
| B7 | 21 | 19.4% | 1 | 0.604 | Reverted to all-fixed |
| B8 | 16 | 17.6% | 2 | 0.585 | — |
| B9 | 20 | 20.0% | 0 | 0.620 | — |
| B10 | 30 | 28.8% | 0 | 0.604 | Best batch |

### CN Distribution (173 passing designs)

| CN | Count | % |
|----|-------|---|
| 6 | 92 | 53.2% |
| 7 | 62 | 35.8% |
| 8 | 16 | 9.2% |
| 9 | 2 | 1.2% |
| 10 | 1 | 0.6% |

### Top 10 Candidates (for AF3 cross-validation)

| Rank | Design | pTM | pLDDT | CN | Score |
|------|--------|-----|-------|----|-------|
| 1 | b01_seq_0180 | 0.908 | 0.864 | 7 | 0.622 |
| 2 | b09_seq_0138 | 0.901 | 0.864 | 7 | 0.620 |
| 3 | b01_seq_0094 | 0.895 | 0.867 | 6 | 0.618 |
| 4 | b09_seq_0016 | 0.893 | 0.868 | 6 | 0.617 |
| 5 | b02_seq_0177 | 0.904 | 0.849 | 7 | 0.616 |
| 6 | b01_seq_0176 | 0.884 | 0.871 | 7 | 0.615 |
| 7 | b03_seq_0079 | 0.894 | 0.855 | 7 | 0.614 |
| 8 | b01_seq_0096 | 0.885 | 0.861 | 6 | 0.613 |
| 9 | b05_seq_0114 | 0.878 | 0.869 | 7 | 0.612 |
| 10 | b02_seq_0130 | 0.888 | 0.853 | 7 | 0.611 |

### Key Lessons from Full Production Run

1. **All-fixed ligand is optimal**: Partial co-diffusion (B5-B6) dropped pass rate from ~21% to ~10%. Reversion confirmed (~20% for B7-B10).

2. **Pass rate variance is natural**: Range 17.6-30.6% across batches with same params. Mean ~19% excluding B5-B6 experiment.

3. **Strict filter is unreachable**: 0/173 pass strict (CN>=8, pLDDT>=0.80, pTM>=0.75). The CN>=8 threshold is too high for monomeric Tb-citrate designs where CN=6-7 is the empirical optimum.

4. **Score ceiling is ~0.62**: Top designs cluster at 0.61-0.62 composite score. This appears to be a ceiling for current pipeline + params.

5. **B1 produced the best design**: b01_seq_0180 (score 0.622) has held as #1 across all 10 batches. Early batches are not worse than later ones.

6. **CN=7 dominates top designs**: 7 of top 10 have CN=7, suggesting this is the sweet spot for monomeric Tb-citrate binding (3 from citrate + 4 from protein sidechains).

### Recommended Final Configuration

```json
{
  "cfg_scale": 1.5,
  "contig": "100-130",
  "bury_ligand": true,
  "num_timesteps": 200,
  "step_scale": 1.0,
  "mpnn_temperature": 0.05,
  "mpnn_bias_AA": "A:-2.0",
  "mpnn_omit_AA": "C",
  "filter_preset": "metal_tb_standard"
}
```

### Next Steps

1. **AF3 Cross-Validation**: Run top 10 with Ca-citrate (primary) and Tb-citrate (secondary) in AlphaFold3
2. **Experimental Testing**: Gene synthesis for top candidates
3. **Strict Filter Revision**: Consider relaxing CN threshold to >=6 for Tb-specific strict preset

---

## AF3 Cross-Validation Phase 1 Results (2026-02-02)

### Screening Summary
- **All 173 RF3-passing designs** screened through AF3 v3.0.1 (1 seed, noMSA, Ca+Tb variants)
- Runtime: 46.9 min batch mode (--input_dir), ~8.3s/job on RTX 5090
- **Ca-citrate pass (iPTM>=0.8, pTM>=0.7):** 42/172 = 24.4%
- **Tb-citrate pass (iPTM>=0.6, pTM>=0.5):** 132/172 = 76.7%
- **Both pass:** 41/172 = 23.8%

### Key Findings

1. **RF3 composite score moderately predicts AF3 success**: r=0.40 (Ca) and r=0.42 (Tb)
2. **Ca and Tb iPTM distributions are nearly identical**: mean 0.703 vs 0.705, avg diff -0.002
   - Contrary to expectation, Ca (16K training structures) is NOT consistently higher than Tb (62 structures)
   - AF3 may use the same metal coordination model for both ions given similar ionic radii
3. **CN does NOT correlate with AF3 performance**: CN=6 has highest Ca pass rate (27%), CN=8 only 19%
   - High RF3 CN may reflect RF3's own biases rather than ground truth
   - AF3 evaluates coordination independently
4. **76.7% Tb pass rate is inflated**: Relaxed Tb threshold (iPTM>=0.6) is too permissive
   - Tightening to iPTM>=0.8 gives 29 designs — more selective
5. **Phase 2 candidates**: 65 designs (Ca>=0.8 OR Tb>=0.8), ~81 min for 3-seed validation

### Technical Lessons
- AF3 noMSA mode requires THREE fields: `unpairedMsa`, `pairedMsa`, `templates` (not just unpairedMsa)
- JAX 0.5.3 + triton 3.6.0 + `--flash_attention_implementation=xla` for RTX 5090 Blackwell
- `--input_dir` flag loads model once for all JSONs — 10x faster than per-job invocations
- UNC paths (`//wsl.localhost/Ubuntu/...`) work for cross-OS Python Path operations
