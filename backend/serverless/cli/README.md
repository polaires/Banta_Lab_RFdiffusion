# CLI Tools

Local-only command-line tools for design analysis and batch operations.
These are **never deployed** to the RunPod container.

## Usage

All scripts add the parent `serverless/` directory to `sys.path` automatically,
so they can import backend modules directly.

```bash
cd backend/serverless
python cli/analyze_design_cli.py output.pdb --metal TB
python cli/run_50_designs.py --num-designs 50
python cli/check_rmsd_rf3.py results/
```
