# Ln-Citrate Binder Design: Iterative RFD3 Exploration Plan

## Project Goal

Repurpose 3C9H (citrate-Mg binder) to create a lanthanide-citrate binder through systematic RFD3 backbone exploration with iterative parameter optimization.

**Target Publication:** High-impact journal (Nature Methods, PNAS, etc.)
**Current Focus:** Phase 1 - RFD3 backbone generation and optimization

---

## Chemistry Context

| Property | Mg²⁺ (3C9H) | Ln³⁺ (Target) | Design Implication |
|----------|-------------|---------------|-------------------|
| CN | 4-6 | 8-9 | Expand coordination shell |
| M-O Distance | 2.0-2.3Å | 2.3-2.6Å | Adjust expected distances |
| Geometry | Octahedral | SAP/TTP | Different spatial arrangement |
| Citrate Role | 2-3 donors | 2-3 donors | Keep citrate, add protein donors |

**Key Challenge:** Expand coordination from 6 → 8-9 while preserving citrate-binding ability.

---

## Directory Structure

```
experiments/ln_citrate_scaffold/
├── DESIGN_PLAN.md              # This file
├── inputs/                     # Input PDB structures
│   ├── 3c9h_original.pdb       # Original citrate-Mg structure
│   ├── 3c9h_ln_prepped.pdb     # Ln-replaced structure
│   ├── citrate_ln_motif.pdb    # Minimal motif for unindexed design
│   └── coordinating_motif.pdb  # 5-residue coordination motif
├── configs/                    # RFD3 JSON configurations
│   ├── round_01/               # First exploration round
│   ├── round_02/               # Refined based on R1 lessons
│   └── ...
├── outputs/                    # Generated PDBs (gitignored)
│   ├── round_01/
│   ├── round_02/
│   └── ...
├── analysis/                   # Analysis results
│   ├── round_01_analysis.json
│   ├── round_02_analysis.json
│   └── cumulative_lessons.md
├── scripts/                    # Automation scripts
│   ├── run_round.py
│   ├── analyze_round.py
│   └── generate_next_config.py
└── final/                      # Production-ready configs + best designs
    ├── production_config.json
    └── top_candidates/
```

---

## Phase 1: RFD3 Backbone Exploration

### Round Structure

Each round follows this cycle:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  ROUND N CYCLE                                                               │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  1. CONFIGURE                                                               │
│     └─ Generate configs based on Round N-1 lessons                          │
│                                                                             │
│  2. EXECUTE                                                                 │
│     └─ Run RFD3 via local Docker                                           │
│     └─ Generate 8-16 designs per config variant                            │
│                                                                             │
│  3. ANALYZE                                                                 │
│     └─ Run UnifiedDesignAnalyzer on all outputs                            │
│     └─ Score against metal_specific filter preset                          │
│     └─ Compute aggregate statistics                                        │
│                                                                             │
│  4. JUDGE (Decision Tree)                                                   │
│     ├─ IF pass_rate < 10% → Revise fundamentally                           │
│     ├─ IF pass_rate 10-30% → Adjust key parameters                         │
│     ├─ IF pass_rate 30-50% → Fine-tune                                     │
│     └─ IF pass_rate > 50% → Move to LigandMPNN                             │
│                                                                             │
│  5. DOCUMENT                                                                │
│     └─ Record lessons learned                                               │
│     └─ Update cumulative_lessons.md                                        │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Round 1: Baseline Exploration

### Objective
Establish baseline performance across three fundamental strategies.

### Configs to Test

#### 1A: Conservative Partial Diffusion (Preserve 3C9H fold)
```json
{
  "r1_partial_conservative": {
    "input": "inputs/3c9h_ln_prepped.pdb",
    "contig": "A1-180",
    "ligand": "CIT,TB",
    "partial_t": 8.0,
    "select_fixed_atoms": {
      "A15": "TIP",
      "A45": "TIP",
      "A78": "TIP",
      "TB": "ALL",
      "CIT": "C1,C2,C3,C4,C5,C6"
    },
    "select_buried": {"TB": "ALL"},
    "select_partially_buried": {"CIT": "O1,O2,O3,O4,O5,O6"},
    "is_non_loopy": true,
    "plddt_enhanced": true
  }
}
```

#### 1B: Moderate Partial Diffusion (More flexibility)
```json
{
  "r1_partial_moderate": {
    "input": "inputs/3c9h_ln_prepped.pdb",
    "contig": "A1-180",
    "ligand": "CIT,TB",
    "partial_t": 12.0,
    "select_fixed_atoms": {
      "A15": "TIP",
      "A45": "TIP",
      "A78": "TIP",
      "TB": "ALL",
      "CIT": ""
    },
    "select_unfixed_sequence": "A20-40,A85-105",
    "select_buried": {"TB": "ALL"},
    "select_partially_buried": {"CIT": "ALL"},
    "is_non_loopy": true,
    "plddt_enhanced": true
  }
}
```

#### 1C: Unindexed Motif Scaffolding (De novo around motif)
```json
{
  "r1_unindex_motif": {
    "input": "inputs/coordinating_motif.pdb",
    "contig": "100-120",
    "length": "105-115",
    "ligand": "CIT,TB",
    "unindex": "E1,E2,D3,D4,N5",
    "select_fixed_atoms": {
      "E1": "TIP",
      "E2": "TIP",
      "D3": "TIP",
      "D4": "TIP",
      "N5": "TIP",
      "TB": "ALL",
      "CIT": ""
    },
    "select_buried": {"TB": "ALL"},
    "select_partially_buried": {"CIT": "O1,O2,O3"},
    "ori_token": [0.0, 0.0, 0.0],
    "is_non_loopy": true,
    "plddt_enhanced": true
  }
}
```

### Execution

```bash
# Start Docker (use docker-wsl skill)
# Then run:
cd backend/serverless
python -c "
import requests
import json

configs = json.load(open('../../experiments/ln_citrate_scaffold/configs/round_01/all_configs.json'))
for name, config in configs.items():
    response = requests.post('http://localhost:8000/runsync', json={
        'input': {
            'task': 'rfd3',
            'config': config,
            'n_designs': 8
        }
    })
    # Save outputs...
"
```

### Analysis Criteria (Round 1)

| Metric | Threshold | Weight | Notes |
|--------|-----------|--------|-------|
| CN (donors within 3.0Å) | ≥7 | 25% | Critical for Ln |
| Ln-O distances | 2.3-2.6Å | 20% | Must be in range |
| Metal SASA | <10Ų | 15% | Should be buried |
| Citrate contacts | ≥2 to Ln | 15% | Citrate must bridge |
| Backbone clashes | None <2.0Å | 15% | Basic quality |
| SS content | ≥50% | 10% | Designability |

### Decision Criteria (Round 1)

```python
def judge_round_1(results):
    """
    Decide next steps based on Round 1 results.
    """
    pass_rates = {
        'partial_conservative': calc_pass_rate(results['1A']),
        'partial_moderate': calc_pass_rate(results['1B']),
        'unindex_motif': calc_pass_rate(results['1C'])
    }

    best_strategy = max(pass_rates, key=pass_rates.get)
    best_rate = pass_rates[best_strategy]

    if best_rate < 0.10:
        return {
            'action': 'REVISE_FUNDAMENTALLY',
            'notes': [
                'Re-examine input structure preparation',
                'Check if Ln-citrate geometry is physically reasonable',
                'Consider alternative coordinating residue arrangements'
            ]
        }
    elif best_rate < 0.30:
        return {
            'action': 'ADJUST_KEY_PARAMETERS',
            'focus_strategy': best_strategy,
            'adjustments': [
                'Vary partial_t (±3Å)',
                'Adjust select_fixed_atoms (BKBN vs TIP)',
                'Try different contig lengths'
            ]
        }
    elif best_rate < 0.50:
        return {
            'action': 'FINE_TUNE',
            'focus_strategy': best_strategy,
            'adjustments': [
                'Add H-bond conditioning',
                'Refine RASA selections',
                'Increase n_designs for statistics'
            ]
        }
    else:
        return {
            'action': 'PROCEED_TO_MPNN',
            'top_backbones': get_top_n(results[best_strategy], n=5)
        }
```

---

## Round 2+: Iterative Refinement

Based on Round 1 lessons, subsequent rounds will:

### If Partial Diffusion Wins:
- Test partial_t range: [6, 8, 10, 12, 15]
- Test fixed atom strategies:
  - TIP only (functional atoms)
  - BKBN (backbone constrained)
  - ALL (full residue fixed)
- Test which residues to fix vs. design

### If Unindexed Motif Wins:
- Test scaffold lengths: [80-100, 100-120, 120-140]
- Test motif compositions:
  - 4 Glu (bidentate) → CN=8
  - 3 Glu + 2 Asp → mixed bidentate/monodentate
  - 2 Glu + 2 Asp + 1 Asn → include amide
- Test with/without H-bond conditioning

### Parameter Sweep Template

```json
{
  "r2_sweep": {
    "base_config": "...",
    "sweep_parameters": {
      "partial_t": [6, 8, 10, 12, 15],
      "contig_length": ["80-100", "100-120", "120-140"],
      "fixed_atom_mode": ["TIP", "BKBN", "ALL"]
    },
    "n_designs_per_combo": 4,
    "total_expected": 180
  }
}
```

---

## Quality Judgment Framework

### Tier System

| Tier | CN | Geometry RMSD | Metal SASA | Citrate Contacts | Action |
|------|-----|---------------|------------|------------------|--------|
| **S** | ≥8 | <0.8Å | <2Ų | ≥3 | → LigandMPNN immediately |
| **A** | ≥8 | 0.8-1.2Å | 2-5Ų | ≥3 | → LigandMPNN (priority) |
| **B** | 7-8 | 1.2-1.5Å | 5-10Ų | 2-3 | → More RFD3 refinement |
| **C** | 6-7 | 1.5-2.0Å | 10-20Ų | 1-2 | → Parameter adjustment needed |
| **F** | <6 | >2.0Å | >20Ų | <1 | → Fundamental revision |

### Automated Scoring Script

```python
# experiments/ln_citrate_scaffold/scripts/score_design.py

import sys
sys.path.insert(0, '../../backend/serverless')
from unified_analyzer import UnifiedDesignAnalyzer

def score_ln_citrate_design(pdb_path):
    """
    Score a single design against Ln-citrate criteria.

    Returns:
        {
            'tier': 'S'|'A'|'B'|'C'|'F',
            'score': float,
            'metrics': {...},
            'pass': bool,
            'issues': [...]
        }
    """
    analyzer = UnifiedDesignAnalyzer(metal='TB')
    result = analyzer.analyze(pdb_path)

    # Extract key metrics
    cn = result.get('coordination', {}).get('coordination_number', 0)
    geom_rmsd = result.get('coordination', {}).get('geometry_rmsd', 999)
    metal_sasa = result.get('burial', {}).get('metal_sasa', 999)
    cit_contacts = result.get('ligand_interface', {}).get('contacts_to_metal', 0)

    # Determine tier
    if cn >= 8 and geom_rmsd < 0.8 and metal_sasa < 2 and cit_contacts >= 3:
        tier = 'S'
    elif cn >= 8 and geom_rmsd < 1.2 and metal_sasa < 5 and cit_contacts >= 3:
        tier = 'A'
    elif cn >= 7 and geom_rmsd < 1.5 and metal_sasa < 10 and cit_contacts >= 2:
        tier = 'B'
    elif cn >= 6 and geom_rmsd < 2.0 and metal_sasa < 20 and cit_contacts >= 1:
        tier = 'C'
    else:
        tier = 'F'

    # Compute composite score (0-100)
    score = (
        min(30, cn * 3.75) +                    # CN: max 30
        max(0, 25 - geom_rmsd * 16.67) +        # Geometry: max 25
        max(0, 20 - metal_sasa * 1.0) +         # Burial: max 20
        min(15, cit_contacts * 5) +             # Citrate: max 15
        min(10, result.get('ss_content', 0) * 10)  # SS: max 10
    )

    return {
        'tier': tier,
        'score': score,
        'metrics': {
            'cn': cn,
            'geometry_rmsd': geom_rmsd,
            'metal_sasa': metal_sasa,
            'citrate_contacts': cit_contacts
        },
        'pass': tier in ['S', 'A', 'B'],
        'issues': identify_issues(cn, geom_rmsd, metal_sasa, cit_contacts)
    }
```

---

## Lesson Documentation Template

After each round, document findings in `analysis/round_XX_lessons.md`:

```markdown
# Round XX Lessons

## Summary
- Designs tested: N
- Pass rate: X%
- Best tier achieved: X
- Top score: XX/100

## What Worked
1. [Specific parameter that improved results]
2. [Configuration that produced best designs]

## What Failed
1. [Parameter/approach that consistently failed]
2. [Unexpected issues encountered]

## Key Insights
- [Insight about coordination geometry]
- [Insight about scaffold flexibility]
- [Insight about citrate positioning]

## Next Round Recommendations
1. [Specific parameter to try]
2. [New strategy to explore]
3. [What to avoid]

## Parameter Evolution
| Parameter | R1 Value | RXX Value | Change Reason |
|-----------|----------|-----------|---------------|
| partial_t | 8.0 | 10.0 | Better geometry at higher noise |
| ... | ... | ... | ... |
```

---

## Decision Tree: When to Proceed to LigandMPNN

```
                    Round Complete
                          │
                          ▼
              ┌───── Pass Rate? ─────┐
              │                      │
         < 30%                    ≥ 30%
              │                      │
              ▼                      ▼
    More RFD3 rounds          ┌── Best Tier? ──┐
    (max 5 rounds)            │                │
              │            S or A           B only
              │               │                │
              ▼               ▼                ▼
    If R5 still < 30%:   Proceed to      1 more round
    → Fundamental         LigandMPNN     to try for A
      revision
```

### LigandMPNN Trigger Conditions

Proceed to LigandMPNN when ANY of these are met:

1. **≥3 Tier A or S designs** from a single round
2. **Pass rate ≥50%** (Tier B or better)
3. **5 rounds completed** with at least one Tier B design

---

## CLI Workflow Commands

### Round Execution

```bash
# 1. Start Docker
# (use docker-wsl skill)

# 2. Generate configs
cd experiments/ln_citrate_scaffold
python scripts/generate_configs.py --round 1

# 3. Run RFD3
python scripts/run_round.py --round 1 --n-designs 8

# 4. Analyze results
python scripts/analyze_round.py --round 1 --save

# 5. Generate judgment report
python scripts/judge_round.py --round 1

# 6. Document lessons
# (manually update analysis/round_01_lessons.md)

# 7. Generate next round configs
python scripts/generate_next_config.py --from-round 1
```

### Quick Analysis

```bash
# Single design
cd backend/serverless
python analyze_design_cli.py \
  ../../experiments/ln_citrate_scaffold/outputs/round_01/design_001.pdb \
  --metal TB \
  --session ln_citrate_r1

# Full round
python -c "
from pathlib import Path
from unified_analyzer import UnifiedDesignAnalyzer
import json

outputs = list(Path('../../experiments/ln_citrate_scaffold/outputs/round_01').glob('*.pdb'))
analyzer = UnifiedDesignAnalyzer(metal='TB')

results = []
for pdb in outputs:
    result = analyzer.analyze(str(pdb))
    results.append({'file': pdb.name, **result})

with open('../../experiments/ln_citrate_scaffold/analysis/round_01_analysis.json', 'w') as f:
    json.dump(results, f, indent=2)
"
```

---

## Final Production Parameters

After iterative refinement, the final production config will be documented in `final/production_config.json`:

```json
{
  "_meta": {
    "derived_from_rounds": [1, 2, 3, 4],
    "total_designs_tested": 320,
    "final_pass_rate": "67%",
    "best_tier_achieved": "S",
    "key_lessons": [
      "partial_t=10 optimal for this system",
      "TIP fixing better than BKBN for coordinators",
      "105-115 residue scaffold length ideal"
    ]
  },
  "production_config": {
    "input": "inputs/optimized_motif.pdb",
    "contig": "105-115",
    "ligand": "CIT,TB",
    "partial_t": 10.0,
    "unindex": "E1,E2,D3,D4,N5",
    "select_fixed_atoms": {
      "E1": "TIP",
      "E2": "TIP",
      "D3": "TIP",
      "D4": "TIP",
      "N5": "TIP",
      "TB": "ALL",
      "CIT": ""
    },
    "select_buried": {"TB": "ALL"},
    "select_partially_buried": {
      "CIT": "O1,O2,O3",
      "E1": "OE1,OE2",
      "E2": "OE1,OE2"
    },
    "select_hbond_acceptor": {"CIT": "O1,O3,O5"},
    "ori_token": [0.0, 0.0, 0.0],
    "is_non_loopy": true,
    "plddt_enhanced": true
  },
  "cli_args": {
    "n_batches": 50,
    "diffusion_batch_size": 8,
    "inference_sampler.num_timesteps": 300,
    "inference_sampler.step_scale": 1.8
  }
}
```

---

## Cumulative Lessons Log

`analysis/cumulative_lessons.md` will track all lessons across rounds:

```markdown
# Cumulative Lessons: Ln-Citrate Scaffold Design

## Round 1 (Baseline)
- [Lessons from R1]

## Round 2 (First Refinement)
- [Lessons from R2]

## Round 3+ (Continued Refinement)
- [Lessons from subsequent rounds]

## Final Synthesis
- [Overall insights for publication methods section]
```

---

## Timeline Estimate

| Phase | Rounds | Designs/Round | Total Designs |
|-------|--------|---------------|---------------|
| Baseline | 1 | 24 (3 configs × 8) | 24 |
| Refinement | 2-3 | 40-60 | 100-180 |
| Fine-tuning | 4-5 | 80-100 | 80-100 |
| Production | Final | 200-400 | 200-400 |
| **Total** | 5-6 | - | **400-700** |

---

## Next Steps

1. **Prepare Input Structures**
   - Download 3C9H from PDB
   - Replace Mg with Tb
   - Position coordinating residues

2. **Create Round 1 Configs**
   - Implement configs 1A, 1B, 1C

3. **Set Up Analysis Pipeline**
   - Test UnifiedDesignAnalyzer with metal=TB
   - Verify filter preset works

4. **Execute Round 1**
   - Start local Docker
   - Run all three config variants
   - Analyze and judge

5. **Document and Iterate**
   - Record lessons
   - Generate Round 2 configs based on findings
