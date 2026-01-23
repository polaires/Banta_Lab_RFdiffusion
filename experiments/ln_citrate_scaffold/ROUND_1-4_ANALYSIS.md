# Ln-Citrate Scaffold Design: Rounds 1-4 Analysis

## Summary

Rounds 1-4 of the lanthanide-citrate scaffold design revealed a **critical flaw in the input structure**: the terbium metal and citrate ligand were positioned ~17 Angstroms apart in the RFD3 backbone designs, when they should be 2-3 Angstroms for proper coordination.

## Round-by-Round Results

### Rounds 1-2: RFD3 Backbone Generation
- Generated backbones with various conditioning (RASA burial, H-bond)
- Achieved max CN=5 (from protein residues)
- **Issue identified**: Metal and citrate not properly co-located

### Round 2b: CFG-Enhanced Generation
- Added classifier-free guidance (cfg_scale=2.0)
- Improved consistency: floor score 17.6 -> 33.6
- Designs with CN>=3: 21% -> 44%
- **Same underlying issue**: Metal-citrate separation

### Round 3: LigandMPNN Sequence Design
- Used metal-aware LigandMPNN with atom context
- No global bias (natural sequences), omit C only
- Generated sequences with D=1-4, E=15-20 per design
- Reasonable D/E content for lanthanide coordination

### Round 4: RF3 Structure Prediction
- 6 predictions completed successfully
- High protein pLDDT: 79-89% (good folding)
- Ligand pLDDT: 62-70% (uncertain position)
- **No coordinating residues near citrate** (0-1 within 6A)

## Root Cause Analysis

### The Problem
In the original input, the terbium was placed at origin (0, 0, 0) while citrate was ~17 Angstroms away. This means:
1. RFD3 wasn't designing a coordinated metal-citrate complex
2. The metal and citrate were in separate regions of the scaffold
3. Coordinating residues were designed around the METAL, not the METAL-CITRATE complex

### Distance Measurements (r2b_b_cfg_hbond_007.pdb)
| Atom Pair | Distance |
|-----------|----------|
| TB to Citrate center | 16.69 A |
| Expected for coordination | 2.0-3.0 A |

## Proposed Fix Strategy

### Option 1: Use Preformed Metal-Citrate Complex (Recommended)

1. **Start with proper coordination geometry**
   - Model citrate-Tb complex with Tb coordinated by citrate carboxylates
   - Use expected Tb-O distances: 2.35-2.45 A
   - Position citrate to provide ~3-4 coordination sites

2. **Input as fixed complex**
   ```json
   {
     "input": "tb_citrate_complex.pdb",
     "contig": "L,M,/0,80-120",
     "select_fixed_atoms": {"L": "all", "M": "all"},
     "select_buried": {"L": "all", "M": "all"},
     "use_classifier_free_guidance": true,
     "cfg_scale": 2.0
   }
   ```

3. **Design scaffold to complete coordination**
   - Target: Tb needs 8-9 coordination
   - Citrate provides ~3-4, protein provides 4-5
   - Use H-bond conditioning for protein carboxylates

### Option 2: Template-Based Approach

1. Use 3C9H (original Mg-citrate binder) as template
2. Replace Mg with Tb (same coordination geometry concept)
3. Run RFD3 partial diffusion to adapt for Tb

### Option 3: Two-Stage Design

1. **Stage 1**: Design citrate-only binder
   - Use citrate alone, buried, with H-bond conditioning

2. **Stage 2**: Introduce metal
   - Dock Tb into the citrate-protein interface
   - Validate with RF3/AF3

## Technical Lessons Learned

1. **CFG is required** for RASA/H-bond conditioning (floor score +16)
2. **LigandMPNN bias is very sensitive** - no global bias with atom context works best
3. **RF3 doesn't predict metals** when not explicitly in training set
4. **Input structure geometry is critical** - garbage in, garbage out

## Files Reference

| File | Description |
|------|-------------|
| outputs/round_02b/*.pdb | RFD3 backbones with CFG |
| outputs/round_03c_mpnn_natural/*.fasta | Designed sequences |
| outputs/round_04_rf3/*_rf3.pdb | RF3 predictions (mmCIF) |
| outputs/round_04_rf3/rf3_analysis.json | Analysis results |

## Next Steps

1. [ ] Model proper Tb-citrate coordination complex
2. [ ] Re-run RFD3 with fixed metal-citrate complex
3. [ ] Validate with RF3 (protein only) + docking (metal-citrate)
4. [ ] Or use AF3 for full prediction with metal/ligand

## Conclusion

The design workflow is sound, but the input structure was flawed. By properly co-locating the metal and citrate before scaffold generation, we should achieve the target coordination geometry of CN=8-9 for the lanthanide.
