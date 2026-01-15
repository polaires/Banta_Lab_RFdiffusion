# Banta Lab RFdiffusion - Development Guide

## Project Overview

AI-powered protein design platform using RFdiffusion3 (Foundry API), LigandMPNN, and PyRosetta. Generates protein structures from natural language input via a React frontend and Python serverless backend.

## Skills Reference

Detailed instructions are in Claude skills - invoke these instead of reading inline docs:

| Skill | Purpose |
|-------|---------|
| `docker-wsl` | Start/stop Docker in WSL, passwords, commands |
| `rfd3-reference` | RFdiffusion3 contig syntax, parameters |
| `ligandmpnn-reference` | Sequence design parameters |
| `bindcraft-reference` | Binder design pipeline |

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
| `mpnn` | Sequence design |
| `interface_ligand_design` | Ligand-binding dimer |
| `interface_metal_design` | Metal-coordinated dimer |
| `binding_eval` | GNINA scoring |
| `fastrelax` | PyRosetta refinement |

## Directory Structure

```
backend/serverless/   # Main API handler (RunPod serverless)
frontend/src/         # React components (Next.js)
docs/                 # All documentation
  getting-started/    # Setup guides
  archive/            # Old/superseded docs
experiments/          # Research projects
  azobenzene_dimer/   # Completed dimer research (archived)
  lanthanide_metal/   # Active metal binding work
scripts/              # Utility scripts
  demos/              # Demo implementations
  tests/              # Test scripts
```

## Current Focus

- Lanthanide metal interface design (`experiments/lanthanide_metal/`)
- Frontend UI improvements

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
| Backend Python code | `backend/serverless/` |
| Backend Python tests | `backend/serverless/test_*.py` |
| JS integration tests | `scripts/tests/` |
| Demo scripts | `scripts/demos/` |
| New experiments | `experiments/<project>/` |
| Generated PDB outputs | `experiments/<project>/outputs/` (gitignored) |
| Documentation | `docs/` |

**Never put in root:**
- Test files (use `scripts/tests/` or `backend/serverless/`)
- PDB files (gitignored, use `experiments/`)
- Temp files (gitignored)

**Gitignored (don't track):**
- `*.pdb` - Generated structures
- `experiments/*/archive/` - Archived outputs
- `tmpclaude-*` - Temp working files
- `.env*.local` - Environment secrets

**Skills:** Use Claude skills for detailed references (docker-wsl, rfd3-reference, etc.)
