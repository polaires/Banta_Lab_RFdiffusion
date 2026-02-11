# Pipeline Orchestration Analysis & Modularization Plan

## 1. Current State: Inventory of All Pipelines/Workflows

### 1.1 Python Pipeline Classes (8 files)

| # | File | Class | Role | Lines |
|---|------|-------|------|-------|
| 1 | `ai_design_pipeline.py` | `AIDesignPipeline` | Full NL-to-report orchestrator (8 stages) | ~800 |
| 2 | `design_orchestrator.py` | `DesignOrchestrator` | Tool selection & config builder | ~400 |
| 3 | `design_validation_pipeline.py` | `DesignValidationPipeline` | Post-design validation (5-stage) | ~600 |
| 4 | `metal_binding_pipeline.py` | `MetalBindingPipeline` | Metal site optimization with quality tiers | ~500 |
| 5 | `production_run.py` | `ProductionRun` (CLI) | 9-stage production pipeline with checkpointing | ~900 |
| 6 | `pipeline_handler.py` | `PipelineHandler` | Sweep/production run orchestrator via API | ~650 |
| 7 | `pipeline_session.py` | `PipelineSessionManager` | Session/progress infrastructure | ~620 |
| 8 | `scaffolding_workflow.py` | `ScaffoldingWorkflow` | PDB motif extraction & theozyme building | ~1850 |

### 1.2 Supporting Decision Modules (5 files)

| File | Purpose |
|------|---------|
| `design_types.py` | 10 `DesignType` enums, presets, `infer_design_type()` |
| `design_rules.py` | 8 decision points → `DesignDecisions` dataclass |
| `rfd3_config_generator.py` | `RFD3Config` (30+ fields), `LigandMPNNConfig` |
| `nl_design_parser.py` | NL → `DesignIntent` (20+ fields) |
| `inference_utils.py` | `run_rfd3_inference()`, `run_mpnn_inference()` |

---

## 2. Overlap Analysis: What Duplicates What

```
                    NL Parse  Scaffold  RFD3  MPNN  Validate  Filter  Report
AIDesignPipeline      [X]       [X]     [X]   [X]    [X]      [X]     [X]
DesignOrchestrator    [ ]       [ ]     [X]   [X]    [X]      [ ]     [ ]
MetalBindingPipeline  [ ]       [X]     [X]   [X]    [X]      [X]     [ ]
ProductionRun (CLI)   [ ]       [X]     [X]   [X]    [X]      [X]     [X]
PipelineHandler       [ ]       [ ]     [X]   [X]    [X]      [X]     [ ]
ValidationPipeline    [ ]       [ ]     [ ]   [X]    [X]      [X]     [ ]
```

**Key overlaps:**
- RFD3 → MPNN → Validate loop is duplicated in **5 places**
- Filtering logic exists in **4 places** with different thresholds
- Scaffold construction in **3 places** (AIDesignPipeline, MetalBindingPipeline, ProductionRun)
- Report generation in **2 places** (AIDesignPipeline, ProductionRun)

---

## 3. Decision Tree: 5-Layer Cascade

```
Layer 0: USER INPUT
    │
    ├── Natural Language ("Design a protein to bind citrate with terbium")
    ├── Structured Config (JSON/CLI params)
    └── Partial Override (NL + specific param overrides)
    │
    ▼
Layer 1: INTENT PARSING  (nl_design_parser.py)
    │
    │  NLDesignParser (Claude API) or SimpleFallbackParser (regex)
    │
    │  Decisions at this layer:
    │  ├── target_molecule → ligand identification
    │  ├── metal_ion → metal type detection
    │  ├── pdb_id → reference structure
    │  ├── enzyme_class → catalytic requirements
    │  ├── design_goal → binding / scaffolding / engineering
    │  ├── symmetry → symmetry group
    │  ├── chain_count → monomer vs oligomer
    │  └── constraints → user-specified overrides
    │
    │  Output: DesignIntent (20+ fields)
    │
    ▼
Layer 2: DESIGN CLASSIFICATION  (design_types.py + design_rules.py)
    │
    │  infer_design_type(intent) → DesignType enum
    │
    │  Priority cascade:
    │  1. Has symmetry?           → SYMMETRIC_*
    │  2. Has enzyme class?       → ENZYME_SCAFFOLD
    │  3. Has metal + ligand?     → METAL_LIGAND_COMPLEX
    │  4. Has metal only?         → METAL_BINDING
    │  5. Has ligand (large)?     → COFACTOR_BINDING
    │  6. Has ligand (small)?     → SMALL_MOLECULE_BINDING
    │  7. Has target protein?     → PROTEIN_BINDER
    │  8. Has DNA/RNA?            → NUCLEIC_ACID_BINDER
    │  9. Else                    → UNCONDITIONAL
    │
    │  make_design_decisions(intent) → DesignDecisions
    │
    │  8 decision points:
    │  ┌─────────────────────────────────────────────────────────┐
    │  │ D1: Chain Length (LigandSize → base/burial/symmetric)   │
    │  │ D2: Hotspot Strategy (hydrophobic/polar/coordinating)   │
    │  │ D3: Burial Strategy (full/partial/exposed/none)         │
    │  │ D4: Orientation (com/hotspots/explicit/none)            │
    │  │ D5: CFG Scale (0.0 to 3.0 based on conditioning)       │
    │  │ D6: MPNN Bias (metal→HSAB, ligand→aromatic/polar/etc)  │
    │  │ D7: Stability Profile (balanced/focused/ultra_stable)   │
    │  │ D8: Enzyme Preservation (catalytic residues/channels)   │
    │  └─────────────────────────────────────────────────────────┘
    │
    │  Output: DesignType + DesignDecisions
    │
    ▼
Layer 3: CONFIG GENERATION  (rfd3_config_generator.py)
    │
    │  RFD3ConfigGenerator.generate(intent, decisions) → RFD3Config + LigandMPNNConfig
    │
    │  RFD3Config (30+ fields):
    │  ├── contig, input_pdb, n_batches
    │  ├── step_scale, gamma_0, num_timesteps
    │  ├── cfg_scale, use_classifier_free_guidance
    │  ├── select_hotspots, select_fixed_atoms
    │  ├── select_buried, select_exposed, select_partially_buried
    │  ├── select_hbond_donor, select_hbond_acceptor
    │  ├── ori_token / infer_ori_strategy
    │  └── symmetry, diffusion_batch_size
    │
    │  LigandMPNNConfig:
    │  ├── model_type, temperature, num_sequences
    │  ├── bias_AA, omit_AA, fixed_residues
    │  └── ligand_mpnn_use_atom_context, pack_side_chains
    │
    │  Output: RFD3Config + LigandMPNNConfig
    │
    ▼
Layer 4: EXECUTION & VALIDATION
    │
    │  ┌─ Scaffolding (optional) ─────────────────────────┐
    │  │  ScaffoldingWorkflow.run()                        │
    │  │  ├── Normalize ligand code                        │
    │  │  ├── Fetch PDB                                    │
    │  │  ├── Find active site (with fallback strategies)  │
    │  │  ├── Extract coordinating residues                │
    │  │  ├── Build theozyme PDB                           │
    │  │  └── Determine conditioning params                │
    │  │  Output: ScaffoldResult (motif_pdb, length,       │
    │  │          unindex, fixed_atoms, rasa, hbonds)      │
    │  └──────────────────────────────────────────────────-┘
    │           │
    │           ▼
    │  ┌─ RFD3 Backbone Generation ───────────────────────┐
    │  │  run_rfd3_inference(rfd3_config)                  │
    │  │  Output: N backbone PDB files                     │
    │  └──────────────────────────────────────────────────-┘
    │           │
    │           ▼ (for each backbone)
    │  ┌─ MPNN Sequence Design ───────────────────────────┐
    │  │  run_mpnn_inference(mpnn_config, backbone_pdb)    │
    │  │  Output: M sequences per backbone                 │
    │  └──────────────────────────────────────────────────-┘
    │           │
    │           ▼ (for each sequence)
    │  ┌─ Structure Prediction ───────────────────────────┐
    │  │  ESMFold (fast, single-chain)                     │
    │  │  or RF3/AF3 (accurate, multi-chain/ligand)        │
    │  │  Output: predicted PDB + confidence scores        │
    │  └──────────────────────────────────────────────────-┘
    │           │
    │           ▼
    │  ┌─ Analysis ───────────────────────────────────────┐
    │  │  UnifiedDesignAnalyzer                            │
    │  │  ├── Backbone RMSD                                │
    │  │  ├── pLDDT, pTM, iPTM                             │
    │  │  ├── Metal coordination (CN, geometry, SASA)      │
    │  │  ├── Ligand contacts                              │
    │  │  └── Interface metrics (pAE, contacts)            │
    │  └──────────────────────────────────────────────────-┘
    │           │
    │           ▼
    │  ┌─ Filtering ──────────────────────────────────────┐
    │  │  Per application type:                            │
    │  │  ├── PPI:    pAE≤1.5, pTM≥0.8, RMSD<2.5Å        │
    │  │  ├── Ligand: RMSD≤1.5Å, lig RMSD≤5Å, iPTM≥0.8   │
    │  │  ├── Metal:  CN≥6, geom RMSD<2.0, pLDDT≥0.70     │
    │  │  └── DNA:    DNA-aligned RMSD<5Å                  │
    │  │  Quality tier: S/A/B/C/F                          │
    │  └──────────────────────────────────────────────────-┘
    │           │
    │           ▼
    │  ┌─ Report/Export ──────────────────────────────────┐
    │  │  ├── Summary statistics                           │
    │  │  ├── Top candidates ranked                        │
    │  │  ├── PDB files saved                              │
    │  │  └── Design history updated                       │
    │  └──────────────────────────────────────────────────-┘
    │
    │  Output: Ranked designs with full metrics
    │
    ▼
Layer 5: ITERATION (optional)
    │
    │  Scout Mode (ProductionRun):
    │  ├── Generate 1 sequence per backbone first
    │  ├── Skip backbones below pTM threshold
    │  └── Full sequence design only for passing backbones
    │
    │  Parameter Sweep (PipelineHandler / MetalBindingPipeline):
    │  ├── 3 sizes × 3 CFG scales = 9 configs
    │  ├── Evaluate each config on small batch
    │  └── Select best config for production run
    │
    │  Iterative Optimization (MetalBindingPipeline):
    │  ├── Round 1: Scout sweep
    │  ├── Round 2: Production with best params
    │  └── Round 3: Fine-tune (adjust bias, temperature)
```

---

## 4. Modular "Lego Block" Architecture

### 4.1 Proposed Modules (13 blocks)

Each block has a clear **Input → Output** contract and can be assembled independently.

```
┌─────────────────────────────────────────────────────────────┐
│                    MODULE CATALOG                             │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  PARSING LAYER                                               │
│  ┌──────────────────┐  ┌──────────────────┐                 │
│  │ M1: NL Parser    │  │ M2: Config       │                 │
│  │                  │  │    Validator      │                 │
│  │ IN:  text str    │  │ IN:  raw dict     │                 │
│  │ OUT: DesignIntent│  │ OUT: validated    │                 │
│  │                  │  │      config dict  │                 │
│  └──────────────────┘  └──────────────────┘                 │
│                                                              │
│  DECISION LAYER                                              │
│  ┌──────────────────┐  ┌──────────────────┐                 │
│  │ M3: Type         │  │ M4: Decision     │                 │
│  │    Classifier    │  │    Engine        │                 │
│  │                  │  │                  │                 │
│  │ IN:  DesignIntent│  │ IN:  DesignIntent│                 │
│  │ OUT: DesignType  │  │      DesignType  │                 │
│  │                  │  │ OUT: Design-     │                 │
│  │                  │  │      Decisions   │                 │
│  └──────────────────┘  └──────────────────┘                 │
│                                                              │
│  CONFIG LAYER                                                │
│  ┌──────────────────┐  ┌──────────────────┐                 │
│  │ M5: RFD3 Config  │  │ M6: MPNN Config  │                 │
│  │    Generator     │  │    Generator     │                 │
│  │                  │  │                  │                 │
│  │ IN:  DesignIntent│  │ IN:  DesignIntent│                 │
│  │      Decisions   │  │      Decisions   │                 │
│  │ OUT: RFD3Config  │  │ OUT: MPNNConfig  │                 │
│  └──────────────────┘  └──────────────────┘                 │
│                                                              │
│  PREPARATION LAYER                                           │
│  ┌──────────────────┐  ┌──────────────────┐                 │
│  │ M7: Scaffolding  │  │ M8: PDB Utils    │                 │
│  │                  │  │                  │                 │
│  │ IN:  PDB + metal │  │ IN:  PDB content │                 │
│  │      + ligand    │  │ OUT: parsed/     │                 │
│  │ OUT: Scaffold-   │  │      sanitized/  │                 │
│  │      Result      │  │      extracted   │                 │
│  └──────────────────┘  └──────────────────┘                 │
│                                                              │
│  EXECUTION LAYER                                             │
│  ┌──────────────────┐  ┌──────────────────┐                 │
│  │ M9: RFD3 Runner  │  │ M10: MPNN Runner │                 │
│  │                  │  │                  │                 │
│  │ IN:  RFD3Config  │  │ IN:  MPNNConfig  │                 │
│  │      + PDB       │  │      + backbone  │                 │
│  │ OUT: backbone    │  │ OUT: sequences   │                 │
│  │      PDB files   │  │      + scores    │                 │
│  └──────────────────┘  └──────────────────┘                 │
│                                                              │
│  VALIDATION LAYER                                            │
│  ┌──────────────────┐  ┌──────────────────┐                 │
│  │ M11: Structure   │  │ M12: Analyzer    │                 │
│  │    Predictor     │  │    & Filter      │                 │
│  │                  │  │                  │                 │
│  │ IN:  sequence    │  │ IN:  design PDB  │                 │
│  │      + context   │  │      + ref PDB   │                 │
│  │ OUT: predicted   │  │      + thresholds│                 │
│  │      PDB + conf  │  │ OUT: metrics +   │                 │
│  │                  │  │      tier + pass/ │                 │
│  │                  │  │      fail        │                 │
│  └──────────────────┘  └──────────────────┘                 │
│                                                              │
│  OUTPUT LAYER                                                │
│  ┌──────────────────┐                                       │
│  │ M13: Reporter    │                                       │
│  │                  │                                       │
│  │ IN:  all results │                                       │
│  │      + metadata  │                                       │
│  │ OUT: report +    │                                       │
│  │      exports     │                                       │
│  └──────────────────┘                                       │
│                                                              │
│  META / ORCHESTRATION                                        │
│  ┌──────────────────┐  ┌──────────────────┐                 │
│  │ M14: Session     │  │ M15: Pipeline    │                 │
│  │    Manager       │  │    Composer      │                 │
│  │                  │  │                  │                 │
│  │ IN:  pipeline    │  │ IN:  module list │                 │
│  │      config      │  │      + config    │                 │
│  │ OUT: checkpoints │  │ OUT: executable  │                 │
│  │      progress    │  │      pipeline    │                 │
│  │      resume      │  │                  │                 │
│  └──────────────────┘  └──────────────────┘                 │
└─────────────────────────────────────────────────────────────┘
```

### 4.2 Module Interface Contracts

```python
# Each module implements this protocol:
class PipelineModule(Protocol):
    name: str                          # Unique module identifier
    version: str                       # Semantic version
    input_schema: dict                 # JSON schema for inputs
    output_schema: dict                # JSON schema for outputs

    async def run(self, inputs: dict, context: PipelineContext) -> dict:
        """Execute module with validated inputs, return outputs."""
        ...

    def validate_inputs(self, inputs: dict) -> ValidationResult:
        """Check inputs before execution."""
        ...

# Pipeline context carries shared state:
@dataclass
class PipelineContext:
    session_id: str
    working_dir: Path
    checkpoint_dir: Path
    design_type: DesignType
    metadata: dict                     # Arbitrary key-value store
    log: logging.Logger
```

### 4.3 Module Dependency Graph

```
M1 (NL Parser)
  └──→ M3 (Type Classifier)
         └──→ M4 (Decision Engine)
                ├──→ M5 (RFD3 Config Gen) ──→ M9 (RFD3 Runner)
                └──→ M6 (MPNN Config Gen) ──→ M10 (MPNN Runner)

M7 (Scaffolding) ──→ M9 (RFD3 Runner)
  ↑
  M8 (PDB Utils)

M9 (RFD3 Runner) ──→ M10 (MPNN Runner) ──→ M11 (Structure Predictor)
                                              └──→ M12 (Analyzer & Filter)
                                                     └──→ M13 (Reporter)

M14 (Session Manager) wraps any pipeline
M15 (Pipeline Composer) assembles any subset of M1-M13
```

---

## 5. Pre-Built Pipeline Recipes (Assembly Examples)

### 5.1 Full NL-to-Report Pipeline

```python
pipeline = PipelineComposer.build([
    "nl_parser",           # M1
    "type_classifier",     # M3
    "decision_engine",     # M4
    "rfd3_config_gen",     # M5
    "mpnn_config_gen",     # M6
    "scaffolding",         # M7 (conditional: only if motif-based)
    "rfd3_runner",         # M9
    "mpnn_runner",         # M10
    "structure_predictor", # M11
    "analyzer_filter",     # M12
    "reporter",            # M13
], session_manager=True, checkpoint=True)

result = await pipeline.run({"text": "Design a protein to bind citrate with terbium"})
```

### 5.2 Quick Validation Only (already have backbone)

```python
pipeline = PipelineComposer.build([
    "mpnn_runner",         # M10
    "structure_predictor", # M11
    "analyzer_filter",     # M12
])

result = await pipeline.run({
    "backbone_pdb": pdb_content,
    "mpnn_config": {"model_type": "ligand_mpnn", "temperature": 0.1}
})
```

### 5.3 Parameter Sweep

```python
pipeline = PipelineComposer.build([
    "rfd3_runner",         # M9
    "mpnn_runner",         # M10 (scout: 1 seq)
    "structure_predictor", # M11
    "analyzer_filter",     # M12
], sweep_mode=True, sweep_params={
    "rfd3_config.chain_length": ["80-100", "120-150", "150-200"],
    "rfd3_config.cfg_scale": [1.5, 2.0, 2.5],
})
```

### 5.4 Scaffolding + RFD3 Only (no sequence design)

```python
pipeline = PipelineComposer.build([
    "scaffolding",         # M7
    "rfd3_runner",         # M9
])

result = await pipeline.run({
    "pdb_id": "4CVB",
    "ligand": "PQQ",
    "metal": "CA",
    "rfd3_config": {"step_scale": 1.5, "n_batches": 50}
})
```

### 5.5 Metal Binding with Quality Tiers

```python
pipeline = PipelineComposer.build([
    "scaffolding",         # M7
    "rfd3_runner",         # M9
    "mpnn_runner",         # M10
    "structure_predictor", # M11
    "analyzer_filter",     # M12 (with metal quality tier thresholds)
    "reporter",            # M13
], filter_preset="metal_binding_quality_tiers")
```

---

## 6. Consolidation Plan: From 8 Files → Modular System

### Phase 1: Extract Core Modules

| Current Code | Becomes Module | What Changes |
|---|---|---|
| `nl_design_parser.py` | **M1: NL Parser** | No change, already clean |
| `design_types.py` | **M3: Type Classifier** | Extract `infer_design_type()` |
| `design_rules.py` | **M4: Decision Engine** | Extract `make_design_decisions()` |
| `rfd3_config_generator.py` | **M5 + M6: Config Generators** | Split RFD3 and MPNN generation |
| `scaffolding_workflow.py` | **M7: Scaffolding** + **M8: PDB Utils** | Split PDB utils (fetch, parse, sanitize) from scaffolding logic |
| `inference_utils.py` | **M9 + M10: Runners** | Already clean, wrap with module interface |
| `design_validation_pipeline.py` | **M11 + M12** | Split prediction from analysis/filtering |
| Report generation (in multiple files) | **M13: Reporter** | Extract and unify |
| `pipeline_session.py` | **M14: Session Manager** | Keep as-is, add checkpoint protocol |

### Phase 2: Build Pipeline Composer (M15)

```python
class PipelineComposer:
    """Assembles modules into executable pipelines."""

    @staticmethod
    def build(
        modules: List[str],
        session_manager: bool = False,
        checkpoint: bool = False,
        sweep_mode: bool = False,
        sweep_params: dict = None,
        filter_preset: str = "default",
    ) -> ExecutablePipeline:
        """Build a pipeline from module names."""
        ...

    @staticmethod
    def from_recipe(recipe_name: str) -> ExecutablePipeline:
        """Build from a named recipe (full_nl, validation_only, sweep, etc.)."""
        ...

    @staticmethod
    def from_design_type(design_type: DesignType) -> ExecutablePipeline:
        """Auto-select modules based on design type."""
        ...
```

### Phase 3: Deprecate Overlapping Pipelines

| Deprecated File | Replaced By |
|---|---|
| `ai_design_pipeline.py` | `PipelineComposer.from_recipe("full_nl")` |
| `metal_binding_pipeline.py` | `PipelineComposer.from_recipe("metal_sweep")` |
| `production_run.py` | `PipelineComposer.from_recipe("production")` with CLI wrapper |
| `pipeline_handler.py` | `PipelineComposer` + API handler (thin wrapper) |
| `design_orchestrator.py` | Absorbed into M3+M4+M5+M6 modules |

### Phase 4: AI Assistant Integration

The AI assistant can compose pipelines dynamically:

```python
# AI assistant receives: "Run just the sequence design on this backbone"
# AI composes:
pipeline = PipelineComposer.build(["mpnn_runner"])
result = await pipeline.run({"backbone_pdb": user_pdb, "mpnn_config": auto_config})

# AI assistant receives: "Sweep 3 scaffold sizes for this metal binding design"
# AI composes:
pipeline = PipelineComposer.build(
    ["scaffolding", "rfd3_runner", "mpnn_runner", "structure_predictor", "analyzer_filter"],
    sweep_mode=True,
    sweep_params={"rfd3_config.chain_length": ["80-100", "120-150", "150-200"]}
)
```

---

## 7. Scaffolding Module Deep-Dive (Needs Most Refactoring)

`scaffolding_workflow.py` (1858 lines) currently mixes multiple concerns:

### Current Responsibilities (to split):

```
scaffolding_workflow.py
├── PDB fetching (→ M8: PDB Utils)
│   ├── fetch_pdb_async()
│   └── PDB format parsing
│
├── Active site detection (→ M7: Scaffolding)
│   ├── find_active_site_residues()
│   ├── _extract_residues_from_pdb()
│   └── 4-level fallback strategy
│
├── Theozyme construction (→ M7: Scaffolding)
│   ├── _extract_theozyme_pdb()
│   └── _build_minimal_theozyme()
│
├── Metal substitution (→ M8: PDB Utils)
│   ├── normalize_ligand_code()
│   ├── is_metal_substitution_compatible()
│   └── HSAB compatibility tables
│
├── Conditioning computation (→ M7: Scaffolding)
│   ├── _determine_rasa_conditioning()
│   ├── _determine_hbond_conditioning()
│   └── _calculate_scaffold_length()
│
└── Backbone continuity check (→ M8: PDB Utils)
    └── check_backbone_continuity()
```

### Proposed Split:

**M8: PDB Utils** (~600 lines)
- PDB fetch, parse, sanitize
- Metal substitution / HSAB tables
- Backbone continuity check
- Coordinate extraction utilities

**M7: Scaffolding** (~1200 lines)
- Active site detection (with fallback)
- Theozyme construction
- Conditioning parameter computation
- ScaffoldResult assembly

---

## 8. Filter Thresholds: Unified Registry

Currently scattered across 4 files. Consolidate into `M12`:

```python
FILTER_PRESETS = {
    "ppi_binder": {
        "inter_chain_pae": {"max": 1.5},
        "ptm": {"min": 0.8},
        "ca_rmsd": {"max": 2.5},
    },
    "small_molecule": {
        "backbone_rmsd": {"max": 1.5},
        "ligand_rmsd": {"max": 5.0},
        "iptm": {"min": 0.8},
        "interface_pae": {"max": 1.5},
    },
    "metal_binding": {
        "coordination_number": {"min": 6},
        "geometry_rmsd": {"max": 2.0},
        "plddt": {"min": 0.70},
        "quality_tiers": True,  # Enable S/A/B/C/F grading
    },
    "metal_binding_strict": {
        "coordination_number": {"min": 8},
        "geometry_rmsd": {"max": 0.8},
        "metal_sasa": {"max": 2.0},
        "plddt": {"min": 0.85},
    },
    "dna_binder": {
        "dna_aligned_rmsd": {"max": 5.0},
    },
    "enzyme_scaffold": {
        "motif_all_atom_rmsd": {"max": 1.5},
    },
    "pipeline_relaxed": {  # For sweep/scout phases
        "ptm": {"min": 0.6},
        "plddt": {"min": 0.65},
    },
    "pipeline_strict": {   # For final filtering
        "ptm": {"min": 0.8},
        "plddt": {"min": 0.75},
        "backbone_rmsd": {"max": 1.5},
    },
}
```

---

## 9. Implementation Priority

| Priority | Task | Rationale |
|---|---|---|
| **P0** | Define `PipelineModule` protocol + `PipelineContext` | Foundation for all modules |
| **P0** | Split `scaffolding_workflow.py` → M7 + M8 | Largest file, most mixed concerns |
| **P1** | Build `PipelineComposer` (M15) | Enables lego assembly |
| **P1** | Unify filter presets into M12 | Removes 4-way duplication |
| **P1** | Wrap existing runners (M9, M10) with module interface | Minimal code change |
| **P2** | Build pre-built recipes | Replace deprecated pipelines |
| **P2** | Add checkpoint protocol to M14 | Enable resume for any pipeline |
| **P3** | Deprecate `ai_design_pipeline.py`, `metal_binding_pipeline.py`, `production_run.py` | After recipes prove equivalent |
| **P3** | AI assistant dynamic composition | After composer is stable |

---

## 10. Summary

**Current state**: 8 pipeline files with significant overlap, especially in the RFD3→MPNN→Validate core loop (duplicated 5x).

**Target state**: 15 modular blocks with clear input/output contracts, assembled via `PipelineComposer` into pre-built recipes or dynamically by the AI assistant.

**Key insight**: The decision cascade (NL → Intent → Type → Decisions → Config → Execution → Validation → Report) is already well-structured in the supporting modules. The problem is that 4 different pipeline classes each re-implement the execution and validation layers. Consolidation means keeping the decision cascade as-is and unifying execution into composable modules.
