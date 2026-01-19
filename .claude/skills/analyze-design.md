# Analyze Design Skill

Use this skill after receiving design results from Docker API or when you have PDB files to analyze.

## When to Use

- After Docker API returns design results (`interface_metal_design`, `interface_ligand_design`, `rfd3`, etc.)
- When user asks to analyze a design
- After any local design run completes
- When user wants to check design quality metrics

## How to Analyze

### From Docker API Response

When you receive a design response with PDB content:

```python
import sys
sys.path.insert(0, 'G:/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless')

from unified_analyzer import UnifiedDesignAnalyzer
from design_history import DesignHistoryManager
from filter_evaluator import FilterEvaluator
from lesson_detector import LessonDetector

# Initialize
analyzer = UnifiedDesignAnalyzer()
history = DesignHistoryManager('G:/Github_local_repo/Banta_Lab_RFdiffusion/experiments/design_history')
evaluator = FilterEvaluator('G:/Github_local_repo/Banta_Lab_RFdiffusion/experiments/design_history/filter_presets')

# Start or continue session
session = history.start_session("docker_design_session")

# For each design in the response
for design in response['result']['designs']:
    pdb_content = design['content']

    # Run analysis
    metrics = analyzer.analyze(
        pdb_content=pdb_content,
        design_params=request_params,  # The original request params
        metal_type="TB" if 'metal' in request_params else None,
    )

    # Evaluate against filter presets
    metrics["filter_results"] = evaluator.evaluate_all_presets(metrics)

    # Save to history
    run_id = history.save_run(session, request_params, {"pdb": pdb_content}, metrics)

    # Check for lesson triggers
    detector = LessonDetector()
    trigger = detector.check_triggers(metrics, [...previous_designs...])

    if trigger:
        print(f"LESSON TRIGGER: {trigger.trigger_type} - {trigger.description}")
```

### From PDB File

```bash
cd G:/Github_local_repo/Banta_Lab_RFdiffusion/backend/serverless
python analyze_design_cli.py /path/to/design.pdb --metal TB --session my_session
```

### CLI Options

| Option | Description |
|--------|-------------|
| `pdb_file` | Path to PDB file (required) |
| `--params FILE` | Design parameters JSON |
| `--ligand FILE` | Ligand SDF file |
| `--metal TYPE` | Metal type code (TB, ZN, CA, etc.) |
| `--metal-chain ID` | Metal chain ID (default: L) |
| `--session NAME` | Session name for grouping designs |
| `--no-save` | Skip saving to history |
| `--json` | Output JSON only |

## Output

The analyzer provides:

1. **Design Type Detection**: monomer, protein_dimer, ligand_interface_dimer, metal_interface_dimer, metal_ligand_interface_dimer

2. **Analyses**:
   - `structure_confidence`: pLDDT scores
   - `interface_quality`: dSASA, contacts, H-bonds
   - `metal_coordination`: coordination number, geometry, HSAB compatibility
   - `ligand_binding`: GNINA affinity scores
   - `symmetry`: C2 symmetry score
   - `topology`: structural validation

3. **Filter Results**: Pass/fail against default, relaxed, stringent, metal_specific, ligand_specific presets

4. **Lesson Triggers**: Alerts for failure patterns, breakthroughs, or significant improvements

## Example Workflow

After running Docker design:

```
User: Run interface_metal_design with TB and contig 70-100
Claude: [Runs Docker API call, receives designs]
Claude: Let me analyze these designs...
[Invokes /analyze-design]
Claude: Analysis complete:
  - Design 1: pLDDT 0.87, coordination distance 2.3Å - PASS (default)
  - Design 2: pLDDT 0.72, coordination distance 2.8Å - FAIL (pLDDT below 0.8)

  Saved to session: docker_design_session
  No lesson triggers detected.
```

## Related

- `experiments/design_history/` - Storage location
- `experiments/design_history/lessons/current_summary.md` - Learned lessons
- `docker-wsl` skill - Start/stop Docker for local testing
