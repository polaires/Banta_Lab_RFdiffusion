# Repository Reorganization Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Clean up repository clutter, reorganize files into proper directories, update CLAUDE.md to slim project overview, and create proper README.md.

**Architecture:** Delete 60+ temp files, move completed research to archive, relocate demo/test scripts to scripts/ directory, update documentation references, and rewrite CLAUDE.md as minimal project overview pointing to skills.

**Tech Stack:** Git, Bash, Markdown

---

## Task 1: Delete Temporary Files

**Files:**
- Delete: All `tmpclaude-*-cwd` files (60 files across repo)
- Delete: All `nul` files (Windows artifacts)
- Delete: `docker.cmd` in root (duplicate)

**Step 1: Delete tmpclaude files**

```bash
cd "G:\Github_local_repo\Banta_Lab_RFdiffusion"
find . -name "tmpclaude-*-cwd" -type f -delete 2>/dev/null || Get-ChildItem -Recurse -Filter "tmpclaude-*-cwd" | Remove-Item -Force
```

**Step 2: Delete nul files**

```bash
rm -f "./nul" "./frontend/nul" "./backend/serverless/nul" 2>/dev/null
```

**Step 3: Delete duplicate docker.cmd**

```bash
rm -f "./docker.cmd"
```

**Step 4: Verify cleanup**

```bash
find . -name "tmpclaude-*" 2>/dev/null | wc -l
# Expected: 0
```

---

## Task 2: Create Directory Structure

**Files:**
- Create: `experiments/azobenzene_dimer/archive/`
- Create: `experiments/azobenzene_dimer/archive/outputs/`
- Create: `scripts/`
- Create: `scripts/demos/`
- Create: `scripts/tests/`
- Create: `docs/getting-started/`
- Create: `docs/archive/`

**Step 1: Create all directories**

```bash
cd "G:\Github_local_repo\Banta_Lab_RFdiffusion"
mkdir -p experiments/azobenzene_dimer/archive/outputs
mkdir -p scripts/demos
mkdir -p scripts/tests
mkdir -p docs/getting-started
mkdir -p docs/archive
```

**Step 2: Verify directories exist**

```bash
ls -d experiments/azobenzene_dimer/archive scripts/demos scripts/tests docs/getting-started docs/archive
# Expected: All directories listed
```

---

## Task 3: Move Untracked PDB and JSON Files to Archive

**Files:**
- Move: 23 `*.pdb` files from root → `experiments/azobenzene_dimer/archive/outputs/`
- Move: `test_result*.json`, `test_full*.json`, `test_sequential.json` → `experiments/azobenzene_dimer/archive/outputs/`

**Step 1: Move PDB files**

```bash
cd "G:\Github_local_repo\Banta_Lab_RFdiffusion"
mv *.pdb experiments/azobenzene_dimer/archive/outputs/ 2>/dev/null || true
```

**Step 2: Move JSON test result files**

```bash
mv test_result.json test_result2.json test_full.json test_full_approach.json test_sequential.json experiments/azobenzene_dimer/archive/outputs/ 2>/dev/null || true
```

**Step 3: Verify files moved**

```bash
ls experiments/azobenzene_dimer/archive/outputs/*.pdb | wc -l
# Expected: ~23
ls *.pdb 2>/dev/null | wc -l
# Expected: 0
```

---

## Task 4: Move Untracked Script Files to Archive

**Files:**
- Move: `analyze_azobenzene.js` → `experiments/azobenzene_dimer/archive/`
- Move: `analyze_result.py` → `experiments/azobenzene_dimer/archive/`
- Move: `compare_structures.py` → `experiments/azobenzene_dimer/archive/`
- Move: `demo_interface_ligand.js` → `experiments/azobenzene_dimer/archive/`

**Step 1: Move untracked analysis scripts**

```bash
cd "G:\Github_local_repo\Banta_Lab_RFdiffusion"
mv analyze_azobenzene.js experiments/azobenzene_dimer/archive/ 2>/dev/null || true
mv analyze_result.py experiments/azobenzene_dimer/archive/ 2>/dev/null || true
mv compare_structures.py experiments/azobenzene_dimer/archive/ 2>/dev/null || true
mv demo_interface_ligand.js experiments/azobenzene_dimer/archive/ 2>/dev/null || true
```

**Step 2: Verify files moved**

```bash
ls experiments/azobenzene_dimer/archive/
# Expected: analyze_azobenzene.js, analyze_result.py, compare_structures.py, demo_interface_ligand.js, outputs/
```

---

## Task 5: Move Docs to Archive

**Files:**
- Move: `RF3_official_README.md` → `docs/archive/`
- Move: `ipd_design_pipeline_collab.ipynb` → `docs/archive/`

**Step 1: Move outdated docs**

```bash
cd "G:\Github_local_repo\Banta_Lab_RFdiffusion"
mv RF3_official_README.md docs/archive/ 2>/dev/null || true
mv ipd_design_pipeline_collab.ipynb docs/archive/ 2>/dev/null || true
```

**Step 2: Verify files moved**

```bash
ls docs/archive/
# Expected: RF3_official_README.md, ipd_design_pipeline_collab.ipynb
```

---

## Task 6: Git Move Tracked Demo/Test Files

**Files:**
- Git mv: `demo_azobenzene_interface.js` → `scripts/demos/`
- Git mv: `demo_cleavable_monomer.js` → `scripts/demos/`
- Git mv: `test_azobenzene_ligand_first.js` → `scripts/tests/`

**Step 1: Move tracked files with git**

```bash
cd "G:\Github_local_repo\Banta_Lab_RFdiffusion"
git mv demo_azobenzene_interface.js scripts/demos/
git mv demo_cleavable_monomer.js scripts/demos/
git mv test_azobenzene_ligand_first.js scripts/tests/
```

**Step 2: Verify git status shows renames**

```bash
git status
# Expected: renamed: demo_azobenzene_interface.js -> scripts/demos/demo_azobenzene_interface.js
# Expected: renamed: demo_cleavable_monomer.js -> scripts/demos/demo_cleavable_monomer.js
# Expected: renamed: test_azobenzene_ligand_first.js -> scripts/tests/test_azobenzene_ligand_first.js
```

---

## Task 7: Update Documentation References

**Files:**
- Modify: `docs/INTERFACE_LIGAND_WORKFLOW.md` (update script paths)

**Step 1: Update demo script references**

In `docs/INTERFACE_LIGAND_WORKFLOW.md`, replace all occurrences of:
- `node demo_azobenzene_interface.js` → `node scripts/demos/demo_azobenzene_interface.js`

**Step 2: Verify changes**

```bash
grep -n "demo_azobenzene_interface" docs/INTERFACE_LIGAND_WORKFLOW.md
# Expected: All references should show scripts/demos/ path
```

---

## Task 8: Update .gitignore

**Files:**
- Modify: `.gitignore`

**Step 1: Add new patterns to .gitignore**

Append to `.gitignore`:
```
# Temporary Claude working files
tmpclaude-*

# Experiment archives (large outputs)
experiments/*/archive/

# Environment files
.env*.local

# Build artifacts
frontend/.next/
frontend/tsconfig.tsbuildinfo
```

**Step 2: Verify .gitignore updated**

```bash
tail -15 .gitignore
# Expected: New patterns visible
```

---

## Task 9: Write New CLAUDE.md

**Files:**
- Overwrite: `CLAUDE.md`

**Step 1: Write slimmed CLAUDE.md (~150 lines)**

```markdown
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
  design-guides/      # RFD3 workflow docs (move existing here)
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

- **Setup:** `docs/getting-started/`
- **Design guides:** `docs/` (RFD3_DESIGN_INSTRUCTION.md, INTERFACE_LIGAND_WORKFLOW.md, etc.)
- **API reference:** `backend/serverless/README.md`
- **Research:** `experiments/*/docs/`

## Git Workflow

1. Make changes locally
2. Test with local Docker (`docker-wsl` skill)
3. Push to `main` → RunPod auto-rebuilds
4. Monitor build at RunPod console

## File Conventions

- **Tracked:** Source code, documentation, configuration
- **Ignored:** PDB outputs, test results, temp files, archives
- **Skills:** Use Claude skills for detailed references (docker-wsl, rfd3-reference, etc.)
```

**Step 2: Verify CLAUDE.md is ~100-150 lines**

```bash
wc -l CLAUDE.md
# Expected: ~100-150 lines
```

---

## Task 10: Write New README.md

**Files:**
- Overwrite: `README.md`

**Step 1: Write proper README.md**

```markdown
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
```

**Step 2: Verify README.md exists**

```bash
cat README.md | head -20
# Expected: New README content
```

---

## Task 11: Commit All Changes

**Files:**
- All modified and moved files

**Step 1: Stage all changes**

```bash
cd "G:\Github_local_repo\Banta_Lab_RFdiffusion"
git add -A
```

**Step 2: Review staged changes**

```bash
git status
# Expected: Many files staged (moves, deletes, new files)
```

**Step 3: Commit**

```bash
git commit -m "chore: reorganize repository structure

- Delete 60+ temp files (tmpclaude-*, nul)
- Move demo/test scripts to scripts/
- Archive completed dimer research
- Slim down CLAUDE.md to project overview
- Add proper README.md
- Update .gitignore for temp files and archives

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

**Step 4: Verify commit**

```bash
git log -1 --oneline
# Expected: Commit hash with message
```

---

## Task 12: Push Changes

**Files:**
- Remote repository

**Step 1: Push to main**

```bash
git push origin main
```

**Step 2: Verify push succeeded**

```bash
git status
# Expected: "Your branch is up to date with 'origin/main'"
```

---

## Summary

After completing all tasks:
- ✅ 60+ temp files deleted
- ✅ PDB/JSON outputs archived
- ✅ Demo/test scripts in scripts/
- ✅ Outdated docs archived
- ✅ CLAUDE.md slimmed to ~100 lines
- ✅ Proper README.md created
- ✅ .gitignore updated
- ✅ Documentation references updated
- ✅ Changes committed and pushed
