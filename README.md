# Banta Lab RFdiffusion

AI-powered protein design platform using RFdiffusion3, LigandMPNN, and PyRosetta.

## Features

- **Backbone Design:** Generate protein backbones with RFdiffusion3
- **Sequence Optimization:** Design sequences with LigandMPNN
- **Ligand Binding:** Design proteins that bind small molecules
- **Metal Coordination:** Design lanthanide-binding proteins
- **Structure Refinement:** PyRosetta FastRelax for clash resolution

## Quick Start

### Production API
```bash
curl -X POST https://api.runpod.ai/v2/k26hmhdfqllx57/run \
  -H "Content-Type: application/json" \
  -H "Authorization: Bearer YOUR_API_KEY" \
  -d '{"input": {"task": "health"}}'
```

### Local Development
See [CLAUDE.md](CLAUDE.md) for setup instructions.

## Documentation

- [Development Guide](CLAUDE.md)
- [RFD3 Design Instructions](docs/RFD3_DESIGN_INSTRUCTION.md)
- [Interface Ligand Workflow](docs/INTERFACE_LIGAND_WORKFLOW.md)
- [Backend API](backend/serverless/README.md)

## Project Structure

```
backend/serverless/   # Python serverless API
frontend/             # Next.js React frontend
docs/                 # Documentation
experiments/          # Research case studies
scripts/              # Utility scripts
```

## License

Proprietary - Banta Lab
