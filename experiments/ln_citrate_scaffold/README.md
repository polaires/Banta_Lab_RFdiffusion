# Ln-Citrate Scaffold Design Experiment

Systematic RFD3 exploration for repurposing 3C9H (citrate-Mg binder) to lanthanide-citrate binder.

## Quick Start

```bash
# 1. Start Docker (from Windows)
# Use docker-wsl skill or:
wsl -d Ubuntu -u root -- service docker start

# 2. Run Round 1
cd experiments/ln_citrate_scaffold
python scripts/run_round.py --round 1 --n-designs 8

# 3. Analyze results
python scripts/analyze_round.py --round 1 --save

# 4. Score individual design
python scripts/score_design.py outputs/round_01/design.pdb
```

## Directory Structure

```
ln_citrate_scaffold/
├── DESIGN_PLAN.md      # Full methodology (READ THIS FIRST)
├── README.md           # This file
├── inputs/             # Input PDB structures
├── configs/            # RFD3 JSON configurations by round
├── outputs/            # Generated PDBs (gitignored)
├── analysis/           # Analysis results + lessons
├── scripts/            # Automation scripts
└── final/              # Production configs + top candidates
```

## Workflow

```
Configure → Execute → Analyze → Judge → Document → Iterate
```

See `DESIGN_PLAN.md` for detailed methodology.

## Key Files

| File | Purpose |
|------|---------|
| `DESIGN_PLAN.md` | Full design methodology |
| `configs/round_01/r1_all_configs.json` | Round 1 configurations |
| `scripts/run_round.py` | Execute RFD3 designs |
| `scripts/analyze_round.py` | Analyze + judge results |
| `scripts/score_design.py` | Score single design |
| `analysis/cumulative_lessons.md` | Learning tracker |

## Tier System

| Tier | CN | Geometry | Metal SASA | Citrate | Action |
|------|-----|----------|------------|---------|--------|
| S | ≥8 | <0.8Å | <2Ų | ≥3 | → LigandMPNN |
| A | ≥8 | <1.2Å | <5Ų | ≥3 | → LigandMPNN |
| B | ≥7 | <1.5Å | <10Ų | ≥2 | More RFD3 |
| C | ≥6 | <2.0Å | <20Ų | ≥1 | Adjust params |
| F | <6 | >2.0Å | >20Ų | <1 | Fundamental revision |

## Current Status

- [ ] Round 1: Baseline exploration
- [ ] Round 2: Parameter refinement
- [ ] Round 3+: Fine-tuning
- [ ] LigandMPNN sequence design
- [ ] Validation pipeline
