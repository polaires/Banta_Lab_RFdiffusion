# Unified Design Analyzer & Local Testing Infrastructure

> **Purpose**: Create a comprehensive local testing and analysis infrastructure for ligand/metal/metal-ligand dimer designs with persistent learning capabilities.

---

## 1. Problem Statement

### Current Pain Points
- **Fragmented analysis**: Tools exist (GNINA, PLIP, metal_validation, binding_analysis) but require manual invocation
- **No feedback loop**: Can't compare designs across runs or identify patterns
- **Missing validation**: Critical checks aren't automated
- **No learning persistence**: Insights lost between sessions

### Goals
1. Unified analysis pipeline for all design types (ligand, metal, metal-ligand dimers)
2. Structured JSON output for future AI agent decision-making
3. Local storage with persistent history
4. Automatic lesson synthesis on significant events
5. CLAUDE.md integration with skill-style references

---

## 2. Architecture Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                     Local Design Workflow                        │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│   Design Run ──► UnifiedDesignAnalyzer ──► DesignHistoryManager │
│       │                   │                        │             │
│       │         ┌─────────┴─────────┐              │             │
│       │         ▼                   ▼              ▼             │
│       │    ┌─────────┐      ┌──────────┐    ┌───────────┐       │
│       │    │Existing │      │   New    │    │  Storage  │       │
│       │    │ Tools   │      │ Integr.  │    │           │       │
│       │    ├─────────┤      ├──────────┤    ├───────────┤       │
│       │    │ GNINA   │      │ ESM      │    │ index.json│       │
│       │    │ metal_* │      │ pLDDT    │    │ runs/     │       │
│       │    │ binding │      │ symmetry │    │ lessons/  │       │
│       │    │ PLIP    │      │ TEBL     │    │ sessions/ │       │
│       │    │ hotspot │      │ PyRosetta│    │ exports/  │       │
│       │    └─────────┘      └──────────┘    └───────────┘       │
│       │                                           │             │
│       │                              ┌────────────┴───────────┐ │
│       │                              │    Lesson Triggers     │ │
│       │                              │  - failure patterns    │ │
│       │                              │  - breakthroughs       │ │
│       │                              │  - improvements        │ │
│       │                              └────────────┬───────────┘ │
│       │                                           │             │
│       │                              ┌────────────▼───────────┐ │
│       │                              │  current_summary.md    │ │
│       │                              │   (CLAUDE.md refs)     │ │
│       │                              └────────────────────────┘ │
└─────────────────────────────────────────────────────────────────┘
```

---

## 3. Core Components

### 3.1 UnifiedDesignAnalyzer

**Entry point:**
```python
analyze_design(pdb_path, design_type, metadata) → StructuredMetricsJSON
```

**Auto-detection:** Infers design type from input (ligand present? metal present? dimer?) and runs appropriate analyses.

**Analysis modules:**
| Analysis Type | Source | Status |
|--------------|--------|--------|
| Binding affinity (GNINA) | `binding_analysis.py` | Existing |
| Interface metrics | `binding_analysis.py` | Existing |
| Shape complementarity | `binding_analysis.py` | Existing |
| Metal coordination | `metal_chemistry.py` + `metal_validation.py` | Existing |
| HSAB validation | `metal_chemistry.py` | Existing |
| Clash detection | `binding_analysis.py` | Existing |
| Topology validation | `topology_validation.py` | Existing |
| PLIP interactions | `test_plip.py` patterns | Needs formalization |
| Hotspot detection | `hotspot_detection.py` | Existing (unused) |
| Conformer quality | `conformer_utils.py` | Existing |
| ESM perplexity | Imported | Needs formalization |
| pLDDT extraction | Manual | Needs implementation |
| Symmetry scoring | `inference_utils.py` | Needs extraction |
| TEBL metrics | `tebl_analysis.py` | Existing |
| PyRosetta interface | `rosetta_utils.py` | Needs enhancement |

**Applicability handling:**
```json
{
  "analyses": {
    "gnina_affinity": {
      "status": "not_applicable",
      "reason": "no small molecule ligand present"
    },
    "metal_coordination": {
      "status": "success",
      "metrics": { "distance": 2.3, "geometry": "octahedral" }
    }
  }
}
```

Three statuses:
- `success` - analysis ran, metrics included
- `not_applicable` - analysis doesn't apply to this design type
- `skipped` - analysis could apply but couldn't run (missing dependency)

---

### 3.2 Filter Presets (BindCraft-inspired)

**Location:** `experiments/design_history/filter_presets/`

**Thresholds:**
| Metric | Default | Relaxed | Stringent |
|--------|---------|---------|-----------|
| pLDDT | ≥0.8 | ≥0.7 | ≥0.85 |
| Shape Complementarity | ≥0.6 | ≥0.5 | ≥0.7 |
| dG (binding energy) | ≤0 | ≤5 | ≤-10 |
| dSASA (buried surface) | ≥500 | ≥300 | ≥700 |
| Interface H-bonds | ≥3 | ≥1 | ≥5 |
| Unsatisfied H-bonds | ≤4 | ≤6 | ≤2 |
| Surface Hydrophobicity | ≤0.35 | ≤0.45 | ≤0.30 |
| GNINA affinity (ligand) | ≤-6 | ≤-4 | ≤-8 |
| Coordination distance (metal) | ≤2.5Å | ≤2.8Å | ≤2.3Å |

---

### 3.3 Local Storage Structure

```
experiments/design_history/
├── index.json                      # Quick lookup with pass/fail status
├── filter_presets/
│   ├── default.json
│   ├── relaxed.json
│   ├── stringent.json
│   ├── metal_specific.json
│   └── ligand_specific.json
├── runs/
│   └── <run_id>/
│       ├── input/
│       │   ├── params.json         # Design parameters
│       │   └── ligand.sdf          # Input ligand if applicable
│       ├── output/
│       │   ├── design.pdb          # Generated structure
│       │   └── mpnn_seqs.fa        # Designed sequences
│       ├── analysis/
│       │   ├── metrics.json        # Full UnifiedDesignAnalyzer output
│       │   └── filter_result.json  # Pass/fail per preset
│       └── meta.json               # Timestamp, design_type, outcome_tag
├── exports/
│   └── all_metrics.csv             # BindCraft-style CSV export
├── sessions/
│   └── <session_id>/
│       ├── session_stats.json      # Acceptance rate, alerts
│       └── run_ids.json            # Designs in this campaign
└── lessons/
    ├── current_summary.md          # Latest synthesized lessons
    └── archive/
        └── <date>_summary.md       # Previous summaries
```

---

### 3.4 Lesson Detection & Synthesis

**Trigger events (significance-driven):**

1. **Failure pattern** - 3+ failures with similar parameters or error modes
2. **Breakthrough success** - Design crosses key thresholds for first time
3. **Meaningful improvement** - Metric improves significantly from previous best

**Detection mechanism:**
```python
def check_lesson_trigger(new_result, history_index):
    # Check failure patterns
    recent_failures = get_recent_by_outcome(history_index, "failure", limit=5)
    if detect_common_pattern(recent_failures):
        return LessonTrigger("failure_pattern", pattern_description)

    # Check breakthrough
    if is_new_best(new_result, history_index):
        return LessonTrigger("breakthrough", improvement_description)

    # Check significant improvement
    if significant_improvement(new_result, history_index, threshold=0.15):
        return LessonTrigger("improvement", delta_description)

    return None
```

**Lesson file format (`current_summary.md`):**
```markdown
# Design Lessons (Auto-updated: 2026-01-18)

## Metal Coordination
- Tb prefers 8-coordinate with Asp/Glu carboxylates (learned: 2026-01-15)
- 6-residue loops insufficient for lanthanide coordination (learned: 2026-01-17)

## Ligand Interface
- Chain length 80-100 residues optimal for interface flexibility (learned: 2026-01-10)
- GNINA threshold < -6 correlates with stable binding (learned: 2026-01-12)

## Failure Patterns to Avoid
- High RASA + tight coordination = steric clashes (observed: 3 failures)
```

---

### 3.5 CLAUDE.md Integration

**Addition to project CLAUDE.md:**
```markdown
## Design History & Lessons

Local design experiments are tracked in `experiments/design_history/`.

**Current lessons learned:** See `experiments/design_history/lessons/current_summary.md`

When running local design tests:
1. Use `UnifiedDesignAnalyzer` for comprehensive metrics
2. Results auto-save to design_history with full provenance
3. Lessons auto-update when significant patterns detected

**Quick reference:**
- Recent runs: `experiments/design_history/index.json`
- Lesson triggers: failure patterns (3+ similar), breakthroughs, significant improvements
```

---

## 4. Structured Output Format

**Full metrics JSON:**
```json
{
  "design_id": "2026-01-18_metal_dimer_001",
  "design_type": "metal_interface_dimer",
  "timestamp": "2026-01-18T14:30:00Z",
  "analyses": {
    "structure_confidence": {
      "status": "success",
      "metrics": { "plddt": 0.85, "pae_mean": 4.2 }
    },
    "metal_coordination": {
      "status": "success",
      "metrics": {
        "distance": 2.3,
        "geometry": "octahedral",
        "hsab_compatible": true,
        "coordinating_residues": ["D45", "E48", "D52", "E55"]
      }
    },
    "interface_quality": {
      "status": "success",
      "metrics": {
        "dG": -12.5,
        "dSASA": 650,
        "shape_complementarity": 0.72,
        "interface_hbonds": 5,
        "unsatisfied_hbonds": 2
      }
    },
    "gnina_affinity": {
      "status": "not_applicable",
      "reason": "no small molecule ligand present"
    },
    "symmetry": {
      "status": "success",
      "metrics": { "c2_score": 0.87 }
    }
  },
  "filter_results": {
    "default": { "pass": true, "failed_filters": [] },
    "stringent": { "pass": false, "failed_filters": ["dG <= -15"] }
  },
  "session": "2026-01-18_tb_dimer",
  "session_stats": {
    "designs_in_session": 12,
    "passing_default": 4,
    "acceptance_rate": 0.33
  }
}
```

---

## 5. Existing Tools Integration Map

| Analysis | Existing Code | Integration Action |
|----------|---------------|-------------------|
| GNINA scoring | `binding_analysis.py:run_gnina_scoring()` | Direct call |
| Interface metrics | `binding_analysis.py:analyze_interface()` | Direct call |
| Shape complementarity | `binding_analysis.py:calculate_shape_complementarity()` | Direct call |
| Metal coordination | `metal_chemistry.py` + `metal_validation.py` | Combine |
| HSAB validation | `metal_chemistry.py:validate_coordination()` | Direct call |
| Clash detection | `binding_analysis.py:check_steric_clashes()` | Direct call |
| Topology | `topology_validation.py` | Direct call |
| PLIP interactions | `test_plip.py` patterns | Formalize module |
| Hotspot detection | `hotspot_detection.py` | Wire in (unused) |
| Conformer quality | `conformer_utils.py` | Wire in |
| ESM perplexity | Imported | Formalize function |
| pLDDT extraction | Manual | Implement extractor |
| Symmetry scoring | `inference_utils.py` | Extract standalone |
| TEBL metrics | `tebl_analysis.py` | Wire in |
| PyRosetta interface | `rosetta_utils.py` | Add `score_interface()` |

---

## 6. Implementation Phases

### Phase 1: Foundation
- Create `UnifiedDesignAnalyzer` class with plugin architecture
- Implement applicability detection
- Wire existing tools: GNINA, metal_validation, binding_analysis, topology
- Define JSON output schema
- Create `DesignHistoryManager` with save/load/index

### Phase 2: Complete Analysis Coverage
- Formalize PLIP as module
- Wire hotspot_detection
- Add ESM perplexity scoring
- Extract symmetry scoring
- Wire TEBL metrics
- Add pLDDT/pAE extraction
- Implement PyRosetta `score_interface()`

### Phase 3: Learning System
- Implement lesson trigger detection
- Create lesson synthesis logic
- Set up `lessons/current_summary.md` auto-update
- Add CLAUDE.md reference integration

### Phase 4: Developer Experience
- CLI wrapper (`python -m analyze_design`)
- Filter preset management
- Session tracking with acceptance rate
- CSV export functionality
- Integration with existing test scripts

---

## 7. Usage Example

```python
from unified_analyzer import UnifiedDesignAnalyzer
from design_history import DesignHistoryManager

analyzer = UnifiedDesignAnalyzer()
history = DesignHistoryManager("experiments/design_history")

# Start a design session
session = history.start_session("tb_dimer_exploration")

# Run your design (existing workflow)
result_pdb = run_local_design(params)

# Analyze comprehensively
metrics = analyzer.analyze(
    pdb_path=result_pdb,
    design_params=params,
    metal_type="Tb"
)

# Auto-save with provenance
run_id = history.save_run(
    session=session,
    params=params,
    outputs={"pdb": result_pdb},
    metrics=metrics
)

# Check for lesson triggers
trigger = history.check_lesson_trigger(metrics)
if trigger:
    print(f"Lesson detected: {trigger.type}")
    # Auto-updates lessons/current_summary.md

# Export for analysis
history.export_metrics_csv()
```

**CLI shortcut:**
```bash
python -m analyze_design output.pdb --params params.json --metal Tb --session tb_exploration
```

---

## 8. Success Criteria

1. **Single entry point**: One function call analyzes any design type
2. **Comprehensive metrics**: All existing tools unified in structured output
3. **Queryable history**: Can ask "show all Tb designs with coordination < 2.5Å"
4. **Automatic learning**: Significant events trigger lesson synthesis
5. **BindCraft parity**: Filter presets, acceptance rate monitoring, CSV export
6. **Minimal friction**: Integrates into existing local dev workflow seamlessly

---

*Design validated: 2026-01-18*
*Ready for implementation planning*
