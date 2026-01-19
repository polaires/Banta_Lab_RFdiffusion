---
name: analyze-design-reminder
enabled: true
event: stop
pattern: .*
conditions:
  - field: stop_reason
    operator: regex_match
    pattern: .*
---

## Before Completing Design Work

If you received design results from Docker API during this session:

**Check:** Did you run analysis on the designs using `/analyze-design`?

**Quick analysis checklist:**
- [ ] Ran `UnifiedDesignAnalyzer` on PDB content
- [ ] Evaluated against filter presets (default, stringent, etc.)
- [ ] Saved to design history with session tracking
- [ ] Checked for lesson triggers

**If designs were not analyzed**, consider:
```python
# Quick analysis
from unified_analyzer import UnifiedDesignAnalyzer
analyzer = UnifiedDesignAnalyzer()
metrics = analyzer.analyze(pdb_content, design_params, metal_type="TB")
```

Or use CLI:
```bash
cd backend/serverless
python analyze_design_cli.py design.pdb --metal TB --session docker_session
```

**Skip this** if no designs were generated in this session.
