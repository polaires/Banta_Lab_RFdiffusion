# Metal-Ligand Complex Validation Pipeline

## Problem Statement

Current approach is reactive - we fix issues one by one as they're discovered visually.
Need a systematic validation and correction pipeline that:
1. Automatically detects geometry issues
2. Accepts feedback (from AI API or human)
3. Applies correction algorithms
4. Iterates until quality passes

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                    DESIGN GENERATION                             │
│  (RFDiffusion → ProteinMPNN → LigandMPNN → FastRelax)           │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                 VALIDATION PIPELINE                              │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐          │
│  │ Geometry     │  │ Clash        │  │ Coordination │          │
│  │ Validator    │  │ Detector     │  │ Checker      │          │
│  └──────────────┘  └──────────────┘  └──────────────┘          │
│         │                 │                 │                   │
│         └─────────────────┼─────────────────┘                   │
│                           ▼                                     │
│              ┌──────────────────────┐                           │
│              │ Validation Report    │                           │
│              │ (JSON structure)     │                           │
│              └──────────────────────┘                           │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                 FEEDBACK INTERFACE                               │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐          │
│  │ AI API       │  │ Human Text   │  │ Automated    │          │
│  │ (Claude/GPT) │  │ Feedback     │  │ Rules        │          │
│  └──────────────┘  └──────────────┘  └──────────────┘          │
│         │                 │                 │                   │
│         └─────────────────┼─────────────────┘                   │
│                           ▼                                     │
│              ┌──────────────────────┐                           │
│              │ Action Items         │                           │
│              │ (structured JSON)    │                           │
│              └──────────────────────┘                           │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                 CORRECTION ENGINE                                │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐          │
│  │ Position     │  │ Clash        │  │ Template     │          │
│  │ Optimizer    │  │ Resolver     │  │ Regenerator  │          │
│  └──────────────┘  └──────────────┘  └──────────────┘          │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
                    [Loop back to Validation]
```

## Module Specifications

### 1. Geometry Validator (`geometry_validator.py`)

```python
class GeometryValidator:
    """Validate metal-ligand-protein geometry."""

    def validate(self, pdb_content: str, metal: str, ligand: str) -> ValidationReport:
        """Run all geometry checks."""
        return ValidationReport(
            metal_ligand_distances=self._check_metal_ligand_distances(),
            carbon_metal_distances=self._check_carbons_far_from_metal(),
            ligand_bond_lengths=self._check_internal_bond_lengths(),
            coordination_geometry=self._check_coordination_geometry(),
        )

    def _check_metal_ligand_distances(self) -> List[DistanceCheck]:
        """Verify coordinating atoms are at correct distance (2.3-2.6Å for Ln-O)."""

    def _check_carbons_far_from_metal(self) -> List[DistanceCheck]:
        """Verify ALL carbons are >3Å from metal."""

    def _check_internal_bond_lengths(self) -> List[BondCheck]:
        """Verify ligand bond lengths are chemically reasonable."""
```

### 2. Clash Detector (`clash_detector.py`)

```python
class ClashDetector:
    """Detect steric clashes between ligand, metal, and protein."""

    VDW_RADII = {
        'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80,
        'H': 1.20, 'TB': 1.75, 'CA': 1.74, ...
    }

    def detect_clashes(self, pdb_content: str) -> List[Clash]:
        """Find all atom pairs with overlapping van der Waals radii."""

    def get_ligand_protein_clashes(self) -> List[Clash]:
        """Specifically check ligand vs protein clashes."""

    def get_metal_protein_clashes(self) -> List[Clash]:
        """Check metal vs protein backbone clashes."""
```

### 3. Coordination Checker (`coordination_checker.py`)

```python
class CoordinationChecker:
    """Verify metal coordination sphere is correct."""

    def check_coordination(self, pdb_content: str, metal: str) -> CoordinationReport:
        """Analyze metal coordination."""
        return CoordinationReport(
            total_coordination=self._count_coordinating_atoms(),
            ligand_donors=self._identify_ligand_donors(),
            protein_donors=self._identify_protein_donors(),
            geometry=self._classify_geometry(),  # octahedral, etc.
            hsab_compliance=self._check_hsab_compatibility(),
        )
```

### 4. Feedback Processor (`feedback_processor.py`)

```python
class FeedbackProcessor:
    """Process text/JSON feedback into actionable corrections."""

    def process_text_feedback(self, feedback: str) -> List[CorrectionAction]:
        """Parse natural language feedback into structured actions.

        Example input: "The citrate is clashing with the protein backbone"
        Example output: [CorrectionAction(type='resolve_clash', target='ligand')]
        """

    def process_ai_analysis(self, ai_response: dict) -> List[CorrectionAction]:
        """Process structured AI API response."""

    def process_validation_report(self, report: ValidationReport) -> List[CorrectionAction]:
        """Auto-generate corrections from validation failures."""
```

### 5. Correction Engine (`correction_engine.py`)

```python
class CorrectionEngine:
    """Apply corrections to fix identified issues."""

    def apply_corrections(
        self,
        pdb_content: str,
        actions: List[CorrectionAction]
    ) -> str:
        """Apply all corrections and return fixed PDB."""

    def resolve_clash(self, pdb: str, clash: Clash) -> str:
        """Move atoms to resolve steric clash."""
        # Options:
        # 1. Translate ligand away from clash
        # 2. Rotate ligand to alternative orientation
        # 3. Flag for protein redesign

    def optimize_position(self, pdb: str, target: str) -> str:
        """Optimize ligand/metal position using energy minimization."""

    def regenerate_with_constraints(self, constraints: dict) -> str:
        """Re-run design with additional constraints."""
```

## Validation Report Schema

```json
{
  "overall_pass": false,
  "quality_score": 65,
  "issues": [
    {
      "type": "clash",
      "severity": "critical",
      "description": "Ligand C3 clashes with protein VAL A:45 CB",
      "atoms": ["CIT:C3", "VAL:A:45:CB"],
      "distance": 1.8,
      "required_distance": 3.2,
      "suggested_action": "translate_ligand"
    },
    {
      "type": "coordination",
      "severity": "warning",
      "description": "Metal coordination incomplete: 4/9",
      "suggested_action": "redesign_binding_site"
    }
  ],
  "geometry": {
    "metal_ligand_distances": {"O2": 2.48, "O5": 2.61, "O7": 2.36},
    "carbon_distances_ok": true,
    "bond_lengths_ok": true
  },
  "coordination": {
    "total": 4,
    "ligand_donors": 3,
    "protein_donors": 1,
    "target": 9
  }
}
```

## Usage Example

```python
from validation_pipeline import ValidationPipeline

pipeline = ValidationPipeline()

# Generate design
pdb = run_design_workflow(template="citrate_tb")

# Validate
report = pipeline.validate(pdb)

if not report.overall_pass:
    # Option 1: Auto-correct
    corrected_pdb = pipeline.auto_correct(pdb, report)

    # Option 2: Get AI feedback
    ai_feedback = call_claude_api(
        prompt=f"Analyze this structure: {report.to_summary()}",
        model="claude-sonnet"
    )
    actions = pipeline.process_feedback(ai_feedback)
    corrected_pdb = pipeline.apply_corrections(pdb, actions)

    # Option 3: Human feedback
    human_feedback = "Move the citrate 2Å away from the protein"
    actions = pipeline.process_text_feedback(human_feedback)
    corrected_pdb = pipeline.apply_corrections(pdb, actions)

# Iterate until pass
while not pipeline.validate(corrected_pdb).overall_pass:
    corrected_pdb = pipeline.auto_correct(corrected_pdb)
```

## Implementation Priority

1. **Phase 1: Core Validation** (immediate)
   - Geometry validator
   - Clash detector
   - Structured report output

2. **Phase 2: Auto-Correction** (next)
   - Position optimization
   - Clash resolution algorithms
   - Constraint-based regeneration

3. **Phase 3: Feedback Interface** (future)
   - Text feedback parser
   - AI API integration
   - Iterative refinement loop

## Files to Create

```
backend/serverless/
├── validation_pipeline/
│   ├── __init__.py
│   ├── geometry_validator.py
│   ├── clash_detector.py
│   ├── coordination_checker.py
│   ├── feedback_processor.py
│   ├── correction_engine.py
│   └── report.py
└── tests/
    └── test_validation_pipeline.py
```
