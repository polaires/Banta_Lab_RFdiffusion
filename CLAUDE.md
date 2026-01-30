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
| `ai_parse` | NL intent parsing (Claude API + keyword fallback) |
| `scaffold_search` | RCSB PDB scaffold discovery (metal+ligand) |
| `scout_filter` | Scout: 1 seq/backbone → RF3 → filter weak backbones |
| `save_design_history` | Persist pipeline run to design history |
| `check_lessons` | Detect failure patterns, breakthroughs, improvements |
| `workflow_run` | Modular pipeline via JSON workflow spec |
| `metal_binding_design` | Metal binding sweep/production (round 7b) |
| `interface_ligand_design` | Ligand-binding dimer |
| `interface_metal_design` | Metal-coordinated dimer |
| `binding_eval` | GNINA scoring |
| `fastrelax` | PyRosetta refinement |

## Directory Structure

```
backend/serverless/               # Main API handler (RunPod serverless)
  handler.py                      #   API router (all tasks incl. scout_filter, save_design_history, check_lessons)
  inference_utils.py              #   RFD3, RF3, MPNN inference engines
  esmfold_utils.py                #   ESMFold + RF3 structure validation (Kabsch RMSD)
  ai_design_pipeline.py           #   DEPRECATED — monolithic NL pipeline (kept as reference)
  design_validation_pipeline.py   #   Post-design validation (MPNN + ESMFold + RF3)
  scaffolding_workflow.py         #   Motif-based scaffold design from PDB
  scaffold_search.py              #   RCSB PDB scaffold discovery (metal+ligand)
  nl_design_parser.py             #   Natural language → DesignIntent
  ligand_resolver.py              #   Ligand name → SMILES/SDF resolution
  rfd3_config_generator.py        #   Intent → RFD3 config + scaffold_to_rfd3_params()
  design_rules.py                 #   Decision engine (stability profiles, bias)
  design_history.py               #   DesignHistoryManager — persistent run storage
  lesson_detector.py              #   LessonDetector — failure/breakthrough/improvement detection
  filter_evaluator.py             #   FILTER_PRESETS + FilterEvaluator class
  unified_analyzer.py             #   Comprehensive design analysis (geometry, contacts, etc.)
  pipeline_types.py               #   StepContext, PipelineStep protocol, result dataclasses
  inference_backend.py            #   InProcessBackend, HTTPBackend, create_backend()
  workflow_runner.py              #   WorkflowRunner + run_workflow_spec()
  iteration_strategies.py         #   ScoutStrategy, SweepStrategy
  pipeline_modules/               #   Composable pipeline step modules
    __init__.py                   #     MODULE_REGISTRY (all 9 modules)
    intent_parser.py              #     M1: NL parse + ligand resolve
    design_configurator.py        #     M2: Intent → RFD3/MPNN configs
    scaffolder.py                 #     M3: PDB motif extraction
    scaffold_searcher.py          #     M4: RCSB PDB scaffold search
    backbone_generator.py         #     M5: RFD3 backbone generation
    sequence_designer.py          #     M6: MPNN sequence design
    structure_predictor.py        #     M7: RF3/ESMFold prediction
    analyzer.py                   #     M8: Analysis + filtering
    reporter.py                   #     M9: Report generation
  shared/                         #   Shared utilities (interaction_analysis, etc.)
  run_50_designs.py               #   CLI: batch design + RF3 validation script
  check_rmsd_rf3.py               #   CLI: RF3 RMSD validation script
frontend/src/                     # React components (Next.js)
  components/ai/                  #   AI assistant panel + pipeline workflow
  components/pipeline/            #   PipelineRunner, StepCard, preview components
  hooks/                          #   useDesignJob hook
  lib/api.ts                      #   FoundryAPI client (all backend tasks)
  lib/pipeline-types.ts           #   PipelineStepDefinition, StepResult TS interfaces
  lib/pipelines/                  #   Pipeline definitions
    natural-language.ts           #     NL pipeline (10 steps, see below)
    shared-steps.ts               #     Reusable step factories
    metal-mpnn-bias.ts            #     Metal-aware MPNN bias config
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

## Modular Pipeline Architecture

Composable "Lego" modules that the AI assistant can chain via Python or JSON workflow specs.

### Modules (in `backend/serverless/pipeline_modules/`)

| Module | Key | Wraps |
|--------|-----|-------|
| IntentParser | `intent_parser` | NLDesignParser + LigandResolver |
| DesignConfigurator | `design_configurator` | design_rules + RFD3ConfigGenerator |
| Scaffolder | `scaffolder` | ScaffoldingWorkflow |
| BackboneGenerator | `backbone_generator` | inference_utils.run_rfd3_inference |
| SequenceDesigner | `sequence_designer` | inference_utils.run_mpnn_inference |
| StructurePredictor | `structure_predictor` | inference_utils.run_rf3_inference |
| Analyzer | `analyzer` | UnifiedDesignAnalyzer + FILTER_PRESETS |
| Reporter | `reporter` | Summary generation |

### Key Files

| File | Purpose |
|------|---------|
| `pipeline_types.py` | `StepContext` (shared state), `PipelineStep` protocol, result dataclasses |
| `inference_backend.py` | `InProcessBackend` (GPU) / `HTTPBackend` (Docker API) / `create_backend("auto")` |
| `workflow_runner.py` | `WorkflowRunner.from_spec(json)` — JSON workflow interpreter, checkpointing, progress |
| `iteration_strategies.py` | `ScoutStrategy` (pre-filter backbones), `SweepStrategy` (parameter grid) |
| `filter_evaluator.py` | `FILTER_PRESETS` dict — centralized min/max thresholds for all design types |

### Composing Workflows

**JSON spec** (AI generates, handler.py executes via `workflow_run` task):
```json
{
  "name": "metal_binding_design",
  "params": {"pdb_id": "4CVB", "ligand": "PQQ", "metal": "CA", "num_designs": 50},
  "steps": [
    {"module": "scaffolder"},
    {"module": "design_configurator"},
    {"module": "backbone_generator"},
    {"module": "sequence_designer"},
    {"module": "structure_predictor"},
    {"module": "analyzer", "params": {"filter_preset": "metal_binding"}},
    {"module": "reporter"}
  ],
  "strategy": {"type": "scout", "ptm_threshold": 0.6}
}
```

**Python** (direct composition):
```python
from pipeline_modules import BackboneGenerator, SequenceDesigner, StructurePredictor, Analyzer
from pipeline_types import StepContext
from inference_backend import create_backend

ctx = StepContext(backend=create_backend("auto"), params={"num_designs": 10})
ctx.rfd3_config = {"contig": "80-120", "cfg_scale": 2.0}
for step in [BackboneGenerator(), SequenceDesigner(), StructurePredictor(), Analyzer()]:
    ctx = step.run(ctx)
```

### Design Rules

- Modules **wrap** existing code — never rewrite `inference_utils.py`, `design_rules.py`, etc.
- `StepContext` holds PDB/FASTA as **strings** (matches actual inference patterns)
- `backend` field on StepContext abstracts in-process vs HTTP execution
- `FILTER_PRESETS` uses `{"min": X}` for quality scores, `{"max": X}` for error metrics
- Frontend uses `'workflow'` job type with `WorkflowProgressCard` for step-by-step display

## Frontend Natural Language Pipeline

The primary user-facing pipeline is defined in `frontend/src/lib/pipelines/natural-language.ts`. It uses the `PipelineRunner` component which provides step-by-step execution with review gates, abort support, and design selection between steps.

### Pipeline Flow (10 steps)

```
parse_intent → resolve_structure → scaffold_search → configure
  → rfd3_nl (backbone generation)
  → scout_filter (optional, general RFD3 path only)
  → mpnn_nl (sequence design)
  → rf3_nl (validation)
  → analysis
  → save_history (automatic)
  → check_lessons (automatic)
```

### Step Details

| Step | Factory | Key | Behavior |
|------|---------|-----|----------|
| Parse Intent | inline | `parse_intent` | AI parser (Claude API) or keyword fallback |
| Resolve Structure | inline | `resolve_structure` | Fetch PDB from RCSB or accept upload |
| Scaffold Search | `createScaffoldSearchStep` | `scaffold_search_nl` | Auto-discover PDB templates with metal+ligand |
| Configure | inline | `configure` | Merge intent + scaffold + params → RFD3 config |
| Structure Gen | inline | `rfd3_nl` | Metal path: sweep (RFD3+MPNN+RF3 in one). General: raw RFD3 |
| Scout Filter | `createScoutFilterStep` | `scout_filter_nl` | 1 seq/backbone → RF3 → filter. Skips for metal sweep or ≤1 backbone |
| MPNN | `createMpnnStep` | `mpnn_nl` | LigandMPNN sequence design with metal bias |
| RF3 Validation | `createRf3Step` | `rf3_nl` | RoseTTAFold3 structure prediction per sequence |
| Analysis | inline | `analysis` | Score, rank, optional metal evaluation |
| Save History | `createSaveHistoryStep` | `save_history_nl` | Persist to backend design_history (non-fatal) |
| Check Lessons | `createCheckLessonsStep` | `check_lessons_nl` | Detect failure patterns/breakthroughs (non-fatal) |

### Preview Components

| Component | Used By | Shows |
|-----------|---------|-------|
| `IntentResultPreview` | parse_intent | Parsed design type, metal, ligand, confidence |
| `ScaffoldSearchResultPreview` | scaffold_search | Candidates table with scores, PDB links |
| `ScoutResultPreview` | scout_filter | Pass/fail bar, per-backbone pTM/pLDDT table |
| `LessonResultPreview` | check_lessons | Colored trigger banners (red/green/blue) |

### Key Design Decisions

- **Scout only on general RFD3 path**: Metal binding sweep already validates internally. Scout is auto-skipped when `session_id` is present (indicating sweep was used).
- **History/lessons are automatic and non-fatal**: Save history and check lessons run silently. Errors don't block the pipeline.
- **All new steps are optional**: Users can skip scout (proceed with all backbones), history, and lessons.
- **Shared step factories** in `shared-steps.ts` are reused across pipelines via `createXxxStep({ id: 'custom_id' })`.

## Current Focus

- **Frontend modular pipeline** as primary design interface (NL → PipelineRunner)
- Scout mode for backbone pre-filtering before full MPNN
- Design history persistence and lesson detection for iterative improvement
- Scaffold search for auto-discovering PDB templates with metal+ligand
- Motif scaffolding for metalloenzymes (PQQ-Ca, Dy-TriNOx)
- Dual-method structure validation (ESMFold + RF3 best-of)

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

Design history is now integrated into the frontend pipeline as automatic steps (`save_history_nl` → `check_lessons_nl`). The backend tasks (`save_design_history`, `check_lessons`) use `DesignHistoryManager` and `LessonDetector` classes.

Local design experiments are also tracked in `experiments/design_history/`.

**Current lessons learned:** See `experiments/design_history/lessons/current_summary.md`

### Pipeline Integration (Automatic)

The NL pipeline automatically saves results and checks for lesson triggers after analysis:
- `save_history_nl`: Persists params, outputs, metrics via `DesignHistoryManager`
- `check_lessons_nl`: Detects failure patterns (3+ similar), breakthroughs, improvements via `LessonDetector`

Both steps are non-fatal — errors don't block the pipeline.

### Manual Analysis

When running local design tests:
1. Use `UnifiedDesignAnalyzer` for comprehensive metrics
2. Results auto-save to design_history with full provenance
3. Lessons auto-update when significant patterns detected

Use the `analyze-design` skill for detailed instructions.

**Quick reference:**
- Recent runs: `experiments/design_history/index.json`
- Filter presets: `experiments/design_history/filter_presets/`
- Lesson triggers: failure patterns (3+ similar), breakthroughs, significant improvements
- Backend classes: `design_history.py` (DesignHistoryManager), `lesson_detector.py` (LessonDetector)

**CLI usage:**
```bash
cd backend/serverless
python analyze_design_cli.py output.pdb --metal TB --session my_exploration
```
