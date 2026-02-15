# Round 6 Plan: Achieving CN=8-9 for Tb-Citrate Scaffold

**Date:** 2026-01-22
**Status:** Ready to run

---

## Problem Statement

Round 5 achieved **CN=2** at 4Å from protein donors (target: **CN=6** to reach total CN=9 with citrate's 3 donors). The scaffold successfully hosts the TB-citrate complex, but needs more coordinating residues.

---

## Round 6 Strategy

### Four Approaches

| Config | Approach | Hypothesis |
|--------|----------|------------|
| **r6_a_partial_resample** | Partial diffusion from R5 best (partial_t=0.5) | Keep working scaffold core, resample loops to add coordinators |
| **r6_b_partial_aggressive** | More aggressive partial diffusion (partial_t=0.7) | Allow larger changes to enable loop insertions |
| **r6_c_dimer_symmetric** | C2 homodimer (50-70 res × 2) | Each chain provides ~3 coordinators = 6 total |
| **r6_d_denovo_extended** | Extended de novo (110-130 res) | Longer scaffold provides space for more coordinators |

### Key Parameters (From RFD3 Best Practices)

- `cfg_scale: 2.0-2.5` - Classifier-free guidance for conditioning enforcement
- `step_scale: 1.5` (η from paper) - Improves sampling
- `num_timesteps: 200` - More diffusion steps for complex designs
- `select_hbond_acceptor: {L: O1,O2,O3,O4,O5,O6,O7}` - All citrate oxygens
- `select_fixed_atoms: {X: ALL, L: ALL}` - Lock TB and citrate positions
- `32 designs per config` - 128 total designs for better statistics

---

## Input Files

| Config | Input | Status |
|--------|-------|--------|
| r6_a, r6_b | `outputs/round_05/r5_b_hbond_011.pdb` | ✓ Exists |
| r6_c, r6_d | `inputs/citrate_ln_only.pdb` | ✓ Exists |

---

## How to Run

### Step 1: Start Docker API

```bash
# In WSL or Docker environment
cd G:/Github_local_repo/Banta_Lab_RFdiffusion
docker-compose up -d
# Wait for API to be healthy
curl http://localhost:8000/health
```

### Step 2: Run Round 6 Designs

```bash
python G:/Github_local_repo/Banta_Lab_RFdiffusion/experiments/ln_citrate_scaffold/scripts/run_r6_direct.py
```

### Step 3: Analyze Results

```bash
python G:/Github_local_repo/Banta_Lab_RFdiffusion/experiments/ln_citrate_scaffold/scripts/analyze_r6.py
```

---

## Success Criteria

| Metric | Target | Rationale |
|--------|--------|-----------|
| CN@4Å (protein) | ≥6 | With citrate (3), total = 9 |
| Carboxylate donors | ≥4 | Glu/Asp for HSAB compatibility |
| Sulfur donors | 0 | Cys is soft, Ln is hard acid |
| Quality | "excellent" | CN≥6 AND carboxylate≥4 |

---

## Expected Outcomes

1. **Partial diffusion (r6_a, r6_b)**: Should maintain the good features of R5 best while adding coordinating loops near the metal.

2. **Dimer (r6_c)**: C2 symmetry naturally places residues from two chains around the metal. Each chain only needs to contribute 3 donors.

3. **Extended de novo (r6_d)**: Longer scaffold provides more opportunities for coordinating residues.

---

## Next Steps After Round 6

1. **LigandMPNN** on best backbones with bias `E:3.0,D:3.0,C:-5.0`
2. **RF3/AF3 validation** of sequence-designed structures
3. **PyRosetta FastRelax** for energy minimization
4. **Experimental candidates** for expression testing

---

## File Locations

```
configs/round_06/r6_configs.json     # RFD3 configurations
scripts/run_r6_direct.py             # Execution script
scripts/analyze_r6.py                # Analysis script
outputs/round_06/                    # Output PDBs (after run)
analysis/round_06_analysis.json      # Analysis results (after analysis)
analysis/ROUND_6_SUMMARY.md          # Human-readable summary (after analysis)
```
