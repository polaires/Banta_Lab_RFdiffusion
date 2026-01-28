# Banta Lab RFdiffusion - Development Guide

## Project Overview

AI-powered protein design platform using RFdiffusion3 (Foundry API), LigandMPNN, ESMFold, RoseTTAFold3, and PyRosetta. Generates protein structures from natural language input via a React frontend and Python serverless backend. Supports de novo design, motif scaffolding, metal-binding proteins, enzyme engineering, and multi-method structure validation.

## Skills Reference

Detailed instructions are in Claude skills - invoke these instead of reading inline docs:

| Skill | Purpose |
|-------|---------|
| `docker-wsl` | Start/stop Docker in WSL, passwords, commands |
| `rfd3-reference` | RFdiffusion3 contig syntax, parameters |
| `ligandmpnn-reference` | Sequence design parameters |
| `bindcraft-reference` | Binder design pipeline |
| `analyze-design` | Analyze design results, save to history, check quality |

## Infrastructure

### Production: RunPod Serverless
- **Endpoint ID:** `k26hmhdfqllx57`
- **API:** `https://api.runpod.ai/v2/k26hmhdfqllx57/run`
- **Deployment:** Push to `main` branch → RunPod auto-builds from GitHub
- **GitHub:** `polaires/Banta_Lab_RFdiffusion` (main branch)
- **GPU:** 48GB
- **Network Volume:** `foundry-checkpoints`
- **Dockerfile:** `backend/serverless/Dockerfile`

### Local Development
- **Backend:** WSL Ubuntu + Docker (see `docker-wsl` skill)
- **Frontend:** Next.js at `frontend/` (port 3000)
- **API Port:** 8000
- **GPU:** NVIDIA RTX 3090 (24GB VRAM)
- **Compose:** `backend/serverless/docker-compose.local.yml`

### Checkpoints
| Environment | Path |
|-------------|------|
| RunPod | `/runpod-volume/checkpoints` |
| Local (WSL) | `/home/polaire/.foundry/checkpoints` |

## API Endpoints

All endpoints via `POST http://localhost:8000/runsync` (local) or RunPod API:

| Task | Purpose |
|------|---------|
| `health` | Health check |
| `rfd3` | Backbone diffusion |
| `rf3` | RoseTTAFold3 structure prediction |
| `mpnn` | Sequence design |
| `validate_design` | Full validation pipeline (MPNN + ESMFold + RF3) |
| `ai_design` | NL-driven end-to-end design pipeline |
| `interface_ligand_design` | Ligand-binding dimer |
| `interface_metal_design` | Metal-coordinated dimer |
| `binding_eval` | GNINA scoring |
| `fastrelax` | PyRosetta refinement |

## Directory Structure

```
backend/serverless/               # Main API handler (RunPod serverless)
  handler.py                      #   API router (all tasks)
  inference_utils.py              #   RFD3, RF3, MPNN inference engines
  esmfold_utils.py                #   ESMFold + RF3 structure validation (Kabsch RMSD)
  ai_design_pipeline.py           #   NL → backbone → sequence → validation pipeline
  design_validation_pipeline.py   #   Post-design validation (MPNN + ESMFold + RF3)
  scaffolding_workflow.py         #   Motif-based scaffold design from PDB
  nl_design_parser.py             #   Natural language → DesignIntent
  ligand_resolver.py              #   Ligand name → SMILES/SDF resolution
  rfd3_config_generator.py        #   Intent → RFD3 config generation
  design_rules.py                 #   Decision engine (stability profiles, bias)
  enzyme_chemistry.py             #   Enzyme class detection + activity preservation
  minimal_motif_selector.py       #   Evidence-based motif residue selection
  unified_analyzer.py             #   Comprehensive design analysis (geometry, contacts, etc.)
  metal_validation.py             #   Metal coordination validation
  rosetta_utils.py                #   PyRosetta refinement + scoring
  shared/                         #   Shared utilities (interaction_analysis, etc.)
  run_50_designs.py               #   CLI: batch design + RF3 validation script
  check_rmsd_rf3.py               #   CLI: RF3 RMSD validation script
frontend/src/                     # React components (Next.js)
  components/ai/                  #   AI design panel + pipeline workflow
  hooks/                          #   useAIDesign hook
docs/                             # All documentation
  getting-started/                #   Setup guides
  archive/                        #   Old/superseded docs
experiments/                      # Research projects (local experiment data)
  azobenzene_dimer/               #   Completed dimer research (archived)
  lanthanide_metal/               #   Active metal binding work
  Dy_TriNOx_scaffold/             #   Dy-TriNOx metalloenzyme scaffolding
  ln_citrate_scaffold/            #   Lanthanide-citrate scaffolding
  design_history/                 #   Design run history + lessons learned
scripts/                          # Utility scripts
  demos/                          #   Demo implementations
  tests/                          #   Test scripts
```

## Current Focus

- Motif scaffolding for metalloenzymes (PQQ-Ca, Dy-TriNOx)
- Dual-method structure validation (ESMFold + RF3 best-of)
- AI design pipeline with NL input → validated sequences
- Enzyme activity preservation during redesign

## Documentation

- **Design guides:** `docs/` (RFD3_DESIGN_INSTRUCTION.md, INTERFACE_LIGAND_WORKFLOW.md, etc.)
- **API reference:** `backend/serverless/README.md`
- **Research:** `experiments/*/docs/`

## Git Workflow

1. Make changes locally
2. Test with local Docker (`docker-wsl` skill)
3. Push to `main` → RunPod auto-rebuilds
4. Monitor build at RunPod console

## File Organization Rules

**New files go here:**
| File Type | Location |
|-----------|----------|
| Backend modules (deployed to container) | `backend/serverless/` |
| Reusable CLI tools (batch runners, analysis) | `backend/serverless/` (committed) |
| Backend Python tests | `backend/serverless/test_*.py` |
| One-off experiment scripts (local only) | `experiments/<project>/` or gitignored |
| JS integration tests | `scripts/tests/` |
| Demo scripts | `scripts/demos/` |
| New experiments | `experiments/<project>/` |
| Generated PDB outputs | `experiments/<project>/outputs/` (gitignored) |
| Documentation | `docs/` |

**Deployed vs local distinction:**
Scripts in `backend/serverless/` fall into two categories:
1. **Deployed modules** — imported by `handler.py` or other modules, shipped in Docker container
2. **CLI/experiment scripts** — run locally via `python script.py`, not imported by container code

Local-only experiment runners (`analyze_batch.py`, `production_run.py`, `sweep_analysis.py`, etc.) are gitignored. Reusable CLI tools (`run_50_designs.py`, `check_rmsd_rf3.py`, `analyze_design_cli.py`) are committed.

**Never put in root:**
- Test files (use `scripts/tests/` or `backend/serverless/`)
- PDB files (gitignored, use `experiments/`)
- Temp files (gitignored)

**Gitignored (don't track):**
- `*.pdb` - Generated structures
- `experiments/*/archive/` - Archived outputs
- `experiments/*/outputs/` - Generated experiment outputs
- `backend/serverless/design_results/` - Design output data
- `backend/serverless/output_*/` - Batch run output directories
- `backend/serverless/outputs/` - General output directory
- `backend/serverless/analyze_batch.py`, `analyze_results.py`, `api_client.py`, `compare_scaffolding.py`, `debug_clash.py`, `production_run.py`, `sweep_analysis.py`, `validate_designs.py` - Local experiment runner scripts
- `backend/serverless/sweep_config_*.json` - Experiment parameter configs
- `tmpclaude-*` - Temp working files
- `.env*.local` - Environment secrets

**Skills:** Use Claude skills for detailed references (docker-wsl, rfd3-reference, etc.)

## Structure Validation Architecture

Both `validate_design` and `ai_design` pipelines run dual-method validation:

1. **ESMFold** — Fast single-sequence prediction via HuggingFace Transformers
2. **RF3 (RoseTTAFold3)** — Ligand-aware prediction via in-container inference

**Best-of logic:** For each candidate, compute `best_plddt = max(ESMFold, RF3)` and `best_rmsd = min(ESMFold, RF3)`. Pass/fail uses the best metrics — a design rescued by either method is viable.

**Graceful degradation:** If RF3 is unavailable, pipelines fall back to ESMFold-only without errors. Controlled by `use_rf3=True` (default) parameter.

**Key files:**
- `esmfold_utils.py` — All validation functions (`validate_structure_esmfold()`, `validate_structure_rf3()`, Kabsch RMSD, CIF/PDB parsers)
- `design_validation_pipeline.py` — Post-design pipeline with `SequenceCandidate` dataclass (RF3 + best-of fields)
- `ai_design_pipeline.py` — End-to-end NL pipeline with `ValidationResult` dataclass (RF3 + best-of fields)

**Thresholds** (in `ESMFOLD_THRESHOLDS`): strict, standard, relaxed, exploratory. The AI pipeline defaults to "exploratory" (lenient RMSD) for early-stage designs.

## Design History & Lessons

Local design experiments are tracked in `experiments/design_history/`.

**Current lessons learned:** See `experiments/design_history/lessons/current_summary.md`

### Auto-Analysis Workflow

**IMPORTANT:** After receiving design results from Docker API, always run analysis:

1. Extract PDB content from API response
2. Run `UnifiedDesignAnalyzer` on each design
3. Evaluate against filter presets
4. Save to design history
5. Check for lesson triggers

Use the `analyze-design` skill for detailed instructions.

### Manual Analysis

When running local design tests:
1. Use `UnifiedDesignAnalyzer` for comprehensive metrics
2. Results auto-save to design_history with full provenance
3. Lessons auto-update when significant patterns detected

**Quick reference:**
- Recent runs: `experiments/design_history/index.json`
- Filter presets: `experiments/design_history/filter_presets/`
- Lesson triggers: failure patterns (3+ similar), breakthroughs, significant improvements

**CLI usage:**
```bash
cd backend/serverless
python analyze_design_cli.py output.pdb --metal TB --session my_exploration
```
