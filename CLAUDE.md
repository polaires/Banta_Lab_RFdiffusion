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

**Core Inference:**
| Task | Purpose |
|------|---------|
| `rfd3` | RFD3 backbone diffusion |
| `rf3` | RoseTTAFold3 structure prediction |
| `mpnn` | LigandMPNN sequence design |
| `rmsd` | RMSD calculation |

**Design Orchestration:**
| Task | Purpose |
|------|---------|
| `design` | Unified design endpoint (DesignOrchestrator) |
| `ai_parse` | NL intent parsing (Claude API + keyword fallback) |
| `scout_filter` | 1 seq/backbone → RF3 → filter weak backbones |
| `workflow_run` | Modular pipeline via JSON workflow spec |

**Specialized Pipelines:**
| Task | Purpose |
|------|---------|
| `validate_design` | Full validation pipeline (MPNN + ESMFold + RF3) |
| `metal_binding_design` | Metal binding sweep/production (round 7b) |
| `interface_ligand_design` | Ligand-binding dimer |
| `interface_metal_design` | Metal-coordinated dimer |
| `interface_metal_ligand_design` | Tri-component (metal+ligand+dimer) |
| `cleavable_monomer` | Cleavable monomer algorithm |
| `protein_binder_design` | Binder design mode |

**Analysis & Scoring:**
| Task | Purpose |
|------|---------|
| `analyze_design` | Unified metrics + FILTER_PRESETS |
| `scaffold_search` | RCSB PDB scaffold discovery (metal+ligand) |
| `analyze_conservation` | ConSurf-style evolutionary scoring |
| `interaction_analysis` | Comprehensive interaction analysis |
| `detect_hotspots` | SASA hotspot detection |
| `binding_eval` | GNINA scoring |
| `resolve_ligand` | Ligand name → SMILES/SDF resolution |
| `analyze_ligand_features` | Ligand chemistry feature extraction |

**Refinement:**
| Task | Purpose |
|------|---------|
| `fastrelax` | PyRosetta refinement |

**Pipeline Management:**
| Task | Purpose |
|------|---------|
| `pipeline_design` | Parameter sweep orchestration |
| `pipeline_status` | Sweep status check |
| `pipeline_cancel` | Sweep cancellation |
| `pipeline_export` | Result export |

**History & Learning:**
| Task | Purpose |
|------|---------|
| `save_design_history` | Persist pipeline run to design history |
| `check_lessons` | Detect failure patterns/breakthroughs |

**ESM3:**
| Task | Purpose |
|------|---------|
| `esm3_score` | ESM3 sequence scoring |
| `esm3_generate` | ESM3 sequence generation |
| `esm3_embed` | ESM3 embeddings |

**Utility:**
| Task | Purpose |
|------|---------|
| `health` | Health check |
| `download_checkpoints` | Model checkpoint download |

## Directory Structure

### Backend Serverless (~116 .py files, ~78K lines)

The `backend/serverless/` directory contains both **deployed modules** (shipped in Docker) and **local-only scripts** (CLI tools, tests, experiment runners). The Dockerfile COPY commands define the exact deployment boundary.

#### Deployed Modules (42 files → Docker container)

```
backend/serverless/
  # === CORE (always imported by handler) ===
  handler.py                      # API router + 6 monolithic design workflows (~10.5K lines)
  inference_utils.py              # RFD3, RF3, MPNN inference engines (164K)
  design_types.py                 # DesignType enum, parameter validation
  design_orchestrator.py          # Unified design endpoint orchestration
  binding_analysis.py             # GNINA scoring, interface analysis (77K)
  rosetta_utils.py                # PyRosetta FastRelax, clash scoring (70K)
  cleavage_utils.py               # Cleavable monomer algorithm, homodimer
  hotspot_detection.py            # SASA-based hotspot detection
  esm_utils.py                    # ESM3 scoring/generation
  esmfold_utils.py                # ESMFold + RF3 validation, Kabsch RMSD
  conformer_utils.py              # xTB/CREST conformer generation
  topology_validation.py          # Topology checks

  # === EXTRACTED HANDLERS (focused task modules, imported by handler.py) ===
  handlers/
    __init__.py                   #   Re-exports all extracted handle_* functions
    history.py                    #   save_design_history, check_lessons
    utility.py                    #   download_checkpoints, delete_file
    esm3.py                       #   esm3_score, esm3_generate, esm3_embed
    ligand.py                     #   resolve_ligand, analyze_ligand_features
    analysis.py                   #   validate_design, binding_eval, fastrelax, detect_hotspots, conservation, interaction_analysis
    orchestration.py              #   ai_parse, scaffold_search, scout_filter
    design_analysis.py            #   analyze_design (UnifiedDesignAnalyzer + filter evaluation)

  # === ANALYSIS ===
  unified_analyzer.py             # Comprehensive design scoring (61K)
  analysis_types.py               # StructureType enum (MONOMER, METAL_MONOMER, etc.), detection utils
  design_validation_pipeline.py   # MPNN + ESMFold + RF3 post-design validation
  filter_evaluator.py             # FILTER_PRESETS — centralized quality thresholds
  geometry_validation.py          # 3D geometry validation utils

  # === METAL & LANTHANIDE ===
  metal_chemistry.py              # HSAB theory, coordination chemistry DB (43K)
  metal_site_fetcher.py           # PDB metal reference template fetching (32K)
  metal_ligand_templates.py       # PQQ-Ca, Citrate-Tb complex templates (41K)
  metal_validation.py             # Lanthanide site validation
  metal_binding_pipeline.py       # Metal binding sweep/production (Round 7b)
  lanthanide_templates.py         # EF-hand, C4 symmetric templates (98K)
  tebl_analysis.py                # TEBL signal prediction, Trp antenna design
  architector_integration.py      # Architector docking/scoring
  hbplus_utils.py                 # HBPlus hydrogen bond analysis

  # === LIGAND ===
  ligand_donors.py                # RDKit donor atom identification (circular import w/ ligand_features)
  ligand_features.py              # Self-growing ligand knowledge base
  ligand_resolver.py              # Name → SMILES/SDF resolution (PubChem, ChemSpider)

  # === NL PIPELINE ===
  nl_design_parser.py             # Natural language → DesignIntent
  rfd3_config_generator.py        # Intent → RFD3 config
  design_rules.py                 # Decision engine (stability profiles, bias)
  enzyme_chemistry.py             # Enzyme-specific chemistry models
  backbone_pre_filter.py          # Pre-filtering for RFD3 backbones
  conservation_analyzer.py        # ConSurf-style evolutionary conservation

  # === SCAFFOLDING ===
  scaffolding_workflow.py         # Motif-based scaffold design from PDB (90K)
  scaffold_search.py              # RCSB PDB scaffold discovery (metal+ligand)

  # === MODULAR PIPELINE SYSTEM ===
  pipeline_types.py               # StepContext, PipelineStep protocol
  inference_backend.py            # InProcessBackend / HTTPBackend abstraction
  workflow_runner.py              # JSON workflow spec interpreter
  iteration_strategies.py         # ScoutStrategy, SweepStrategy
  design_history.py               # DesignHistoryManager — persistent run storage
  lesson_detector.py              # LessonDetector — failure/breakthrough detection

  # === SWEEP ORCHESTRATION ===
  pipeline_handler.py             # Parameter sweep/production run orchestration (4 handler tasks)
  pipeline_session.py             # Session state, design filtering, sweep configs

  # === MOTIF SELECTION ===
  minimal_motif_selector.py       # Multi-source evidence motif extraction (M-CSA+PLIP+ConSurf)
  mcsa_client_lite.py             # M-CSA catalytic residue queries

  pipeline_modules/               # Composable step modules (9 modules)
    __init__.py                   #   MODULE_REGISTRY
    intent_parser.py              #   M1: NL parse + ligand resolve
    design_configurator.py        #   M2: Intent → RFD3/MPNN configs
    scaffolder.py                 #   M3: PDB motif extraction
    scaffold_searcher.py          #   M4: RCSB PDB scaffold search
    backbone_generator.py         #   M5: RFD3 backbone generation
    sequence_designer.py          #   M6: MPNN sequence design
    structure_predictor.py        #   M7: RF3/ESMFold prediction
    analyzer.py                   #   M8: Analysis + filtering
    reporter.py                   #   M9: Report generation

  database_adapters/              # External database integrations
    pubchem_adapter.py            #   PubChem SMILES lookup (used by ligand_resolver)
    metalpdb_adapter.py           #   MetalPDB metal site queries
    alphafold_adapter.py          #   AlphaFold structure download
    uniprot_adapter.py            #   UniProt sequence/domain data

  shared/                         # Shared utilities (interaction_analysis)
  utils/                          # SASA calculator, DSSP
  validation_pipeline/            # Geometry, clash detection, coordination analysis
```

#### Archived (deprecated/dead code, moved to `archive/backend_deprecated/`)

```
  ai_design_pipeline.py           # DEPRECATED — monolithic NL pipeline, replaced by pipeline_modules/
  ai_structure_interface.py       # Wrapper of deployed modules, never used by handler.py
  structure_discovery.py          # Multi-DB abstraction, replaced by scaffold_search.py
  parse_result.py                 # Trivial JSON utility, not imported by anything
```

#### CLI Tools (committed in `cli/`, local-only, never deployed)

```
  cli/
    README.md                     # Usage instructions
    run_50_designs.py             # Batch RFD3 + RF3 validation
    check_rmsd_rf3.py             # RF3 RMSD validation
    analyze_design_cli.py         # Single design analyzer
    analyze_coordination.py       # Coordination geometry analyzer
    analyze_designs.py            # Minimal design analyzer
    analyze_metal_dimer.py        # Metal dimer analysis
    analyze_tb_coordination.py    # Terbium coordination specific
    extract_best_pdb.py           # PDB extraction utility
    validate_foldability.py       # Foldability validation
    validate_tb_foldability.py    # Tb-specific foldability
```

#### Test Files (organized in `tests/` subdirectory)

```
  tests/
    conftest.py              # Adds serverless dir to sys.path for all tests
    unit/                    # Fast, no Docker needed (26 files)
      test_design_types.py, test_filter_evaluator.py, test_metal_chemistry.py, ...
    integration/             # Needs running Docker at localhost:8000 (8 files)
      test_local.py, test_ligandmpnn_features.py, test_design_integration.py, ...
    e2e/                     # Full pipelines, 10-50 min, GPU (10+ files)
      test_e2e_pqq_ca.py, test_production_pipeline.py, test_nl_vs_r7b.py, ...
    fixtures/                # Test data (JSON, PDB inputs)
      test_input.json, test_cif_metal.json, ...
    outputs/                 # Gitignored — test-generated files go here
```

**Running tests:**
```bash
pytest backend/serverless/tests/unit/              # Fast unit tests only
pytest backend/serverless/tests/integration/       # Needs Docker running
pytest backend/serverless/tests/e2e/               # Full pipelines (slow)
pytest backend/serverless/tests/                   # Everything
```

#### Gitignored Scripts (local experiment runners, not tracked)

```
  analyze_batch.py, analyze_results.py, api_client.py,
  compare_scaffolding.py, debug_clash.py, sweep_analysis.py,
  validate_designs.py, production_run.py,
  production_runner.py, scaffolding_production_runner.py,
  run_batch_interactive.py, ml_donor_prediction.py
```

### Frontend, Docs, Experiments, Scripts

```
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
  design_outputs/                 #   Organized batch outputs (gitignored)
    production/                   #     Production batches (b01-b10)
    exploration/                  #     Exploration rounds (r1-r5)
    other/                        #     Misc outputs (nl_pqq_ca, etc.)
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

### Deployment Boundary (Critical)

The **Dockerfile** (`backend/serverless/Dockerfile`) defines exactly what ships to RunPod. If you create a new backend module that handler.py imports, you **MUST** also add a `COPY` line to the Dockerfile or it will silently fail in production (handler uses `try/except` for most imports).

**Checklist for new backend modules:**
1. Create the `.py` file in `backend/serverless/`
2. Add `COPY filename.py /app/` to the Dockerfile (or `COPY dir/ /app/dir/` for directories)
3. If it's a new API task handler, add it to `handlers/` and re-export from `handlers/__init__.py`
4. Route the task in `handler.py`'s dispatcher (the `if task == ...` chain)
5. Add to the "Deployed Modules" section of this file

### File Placement

| File Type | Location |
|-----------|----------|
| Deployed backend modules | `backend/serverless/` + add to Dockerfile |
| New API task handlers | `backend/serverless/handlers/` + re-export in `__init__.py` |
| Reusable CLI tools (committed, not deployed) | `backend/serverless/cli/` |
| One-off experiment scripts (local only) | `backend/serverless/` (gitignored) or `experiments/<project>/` |
| Unit tests (fast, mocked) | `backend/serverless/tests/unit/` |
| Integration tests (need Docker) | `backend/serverless/tests/integration/` |
| E2E/experiment tests (full pipeline) | `backend/serverless/tests/e2e/` |
| Test data (JSON, PDB fixtures) | `backend/serverless/tests/fixtures/` |
| Test output files | `backend/serverless/tests/outputs/` (gitignored) |
| JS integration tests | `scripts/tests/` |
| Demo scripts | `scripts/demos/` |
| New experiments | `experiments/<project>/` |
| Generated PDB outputs | `experiments/<project>/outputs/` (gitignored) |
| Batch run outputs | `experiments/design_outputs/` (gitignored) |
| Documentation | `docs/` |

### Test Organization Rules

**NEVER place `test_*.py` files directly in `backend/serverless/`.** All tests go in the `tests/` subdirectory:

- **`tests/unit/`** — Pure pytest, mock all external dependencies, no Docker/network, < 5 sec
- **`tests/integration/`** — Calls local Docker API (localhost:8000), needs running container, 30-120 sec
- **`tests/e2e/`** — Full pipeline runs, produces output files, needs GPU, 10-50 min

**Test file conventions:**
- Import modules with bare imports (`from handler import ...`) — `conftest.py` handles `sys.path`
- Never use `sys.path.insert()` in individual test files
- Never hardcode absolute Windows paths — use `Path(__file__).parent` or `conftest.FIXTURES_DIR`
- Write all test outputs to `conftest.TEST_OUTPUT_DIR` (gitignored), never to the source directory
- Test data fixtures go in `tests/fixtures/`, not alongside production code

### Three Categories in `backend/serverless/`

1. **Deployed modules** — COPY'd in Dockerfile, imported by handler.py or other deployed modules
2. **Committed CLI tools** — tracked in git, run locally via `python script.py`, never imported by container
3. **Gitignored scripts** — local experiment runners, not tracked, not deployed

### Never Put in Root

- Test files (use `scripts/tests/` or `backend/serverless/`)
- PDB files (gitignored, use `experiments/`)
- Temp files (gitignored)

### Gitignored (don't track)

- `*.pdb`, `*.fasta` - Generated structures and sequences
- `experiments/design_outputs/` - All batch run outputs (organized by type)
- `experiments/*/archive/` - Archived experiment outputs
- `experiments/*/outputs/` - Generated experiment outputs
- `backend/serverless/test_data/` - Test data
- `backend/serverless/sweep_config_*.json` - Experiment parameter configs
- `tmpclaude-*` - Temp working files
- `.env*.local` - Environment secrets
- Local experiment runner scripts in `backend/serverless/`:
  - `analyze_batch.py`, `analyze_results.py`, `api_client.py`
  - `compare_scaffolding.py`, `debug_clash.py`, `sweep_analysis.py`
  - `validate_designs.py`, `production_run.py`
  - `production_runner.py`, `scaffolding_production_runner.py`
  - `run_batch_interactive.py`, `ml_donor_prediction.py`

### Known Architecture Issues

1. **Circular import**: `ligand_donors.py` ↔ `ligand_features.py` (low risk — runtime function calls only, not module-level).
2. **handler.py still has 6 monolithic design workflows (~8.5K lines)**: `handle_protein_binder_design` (741 lines), `handle_cleavable_monomer` (287 lines), `handle_interface_ligand_design` (2363 lines), `handle_interface_metal_design` (1188 lines), `handle_interface_metal_ligand_design` (3701 lines), `handle_metal_binding_design` (733 lines). Each has 10-60 inline helper functions. Extracting these requires moving helpers as cohesive units.
3. **Sweep orchestration overlap**: `pipeline_handler.py` and `metal_binding_pipeline.py` both implement sweep/production patterns with different state management. Could be consolidated.

### Recently Resolved Issues

- **DesignType duplication** (Feb 2026): `analysis_types.py` renamed to `StructureType` (structural composition: MONOMER, METAL_MONOMER, etc.). `design_types.py` keeps `DesignType` (design intent: METAL_BINDING, ENZYME_ACTIVE_SITE, etc.). No collision.
- **handler.py partial decomposition** (Feb 2026): 20 task handlers (~1650 lines) extracted to `handlers/` modules. Remaining 6 monolithic design workflows still inline (see above).
- **VALIDATION_PIPELINE_AVAILABLE shadowing** (Feb 2026): Two different modules overwrote the same flag. Renamed to `GEOMETRY_VALIDATION_AVAILABLE` for the geometry/clash module. Removed dead `design_validation_pipeline` imports (now handled lazily in `handlers/analysis.py`).
- **Output directory cleanup** (Feb 2026): 18 `output_*` dirs + loose PDBs moved from `backend/serverless/` to `experiments/design_outputs/` (organized by production/exploration/other). Tracked junk files (design_result_raw.json, .fasta, STRUCTURE_DISCOVERY_README.md) removed. Stale `experiments/` and `docs/` subdirs inside serverless removed.

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
- `unified_analyzer.py` — Comprehensive design scoring with filter presets

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
