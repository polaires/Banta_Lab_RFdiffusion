# RunPod Quick Setup Guide

## Overview
This guide enables fast GPU backend setup on RunPod for the Foundry Protein Design frontend.

## Prerequisites
- RunPod account with credits
- Network volume `foundry-checkpoints` (already created, persists `main.py`)

---

## Quick Start (< 2 minutes)

### Step 1: Deploy Pod
1. Go to [RunPod Console](https://console.runpod.io/deploy)
2. Select **A40** GPU ($0.40/hr) or any available GPU
3. Click **Edit** on the template section
4. Set **HTTP Ports**: `8888, 8000` (IMPORTANT: add port 8000!)
5. Ensure network volume `foundry-checkpoints` is attached
6. Click **Set Overrides** → **Deploy On-Demand**

### Step 2: Start FastAPI Server
1. Wait for pod to show "Ready" status
2. Click **Jupyter Lab** link
3. Create new Python notebook and run these cells:

**Cell 1 - Install dependencies:**
```python
!pip install fastapi uvicorn python-multipart -q
```

**Cell 2 - Download latest main.py from GitHub:**
```python
!curl -o main.py https://raw.githubusercontent.com/polaires/Banta_Lab_RFdiffusion/main/backend/main.py
```

**Cell 3 - Start server:**
```python
import subprocess
subprocess.Popen(["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"])
```

**Cell 4 - Verify health:**
```python
!curl -s http://localhost:8000/health
```

Expected output:
```json
{"status":"healthy","gpu_available":true,"models_loaded":["rfd3","rf3","proteinmpnn"]}
```

**Cell 5 - Test RFD3 endpoint:**
```python
!curl -s -X POST http://localhost:8000/api/rfd3/design -H "Content-Type: application/json" -d '{"contig": "100"}'
```

Expected output:
```json
{"job_id":"<uuid>","status":"pending"}
```

### Step 3: Connect Frontend
1. Copy your pod's API URL: `https://{POD_ID}-8000.proxy.runpod.net`
2. Update `frontend/.env.local`:
```
NEXT_PUBLIC_API_URL=https://{POD_ID}-8000.proxy.runpod.net
```
3. Restart frontend dev server or refresh browser

---

## Real Foundry Mode Setup (Advanced)

For real protein design (not mock data), you need to install `rc-foundry` and download model checkpoints.

### Step 1: Create Python 3.12 Virtual Environment
```python
# rc-foundry requires Python 3.12+
!apt update && apt install -y python3.12 python3.12-venv python3.12-dev -qq

# Create and activate venv
!python3.12 -m venv /workspace/foundry_env
!/workspace/foundry_env/bin/pip install --upgrade pip -q
```

### Step 2: Install rc-foundry
```python
!/workspace/foundry_env/bin/pip install "rc-foundry[all]" -q
```

### Step 3: Configure Checkpoint Directory
```python
import os
os.environ["FOUNDRY_CHECKPOINT_DIRS"] = "/workspace/checkpoints"
!mkdir -p /workspace/checkpoints
```

### Step 4: Start Backend with Foundry
```python
import subprocess
env = os.environ.copy()
env["FOUNDRY_CHECKPOINT_DIRS"] = "/workspace/checkpoints"
subprocess.Popen(
    ["/workspace/foundry_env/bin/python", "-m", "uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"],
    env=env,
    cwd="/workspace"
)
```

### Current Limitations
- **Checkpoint Download**: Model checkpoints (~4GB each) need to be downloaded from IPD servers
- The rc-foundry package should auto-download on first inference, but this may require specific setup
- For now, use **mock mode** for testing the frontend workflow

---

## Backend Features

The unified backend (`main.py`) includes:

### Dual Mode Operation
- **Auto Mode** (default): Detects if Foundry CLI is available
- **Mock Mode**: Returns synthetic data for testing
- **Real Mode**: Uses actual Foundry models

Set via environment variable:
```bash
export MOCK_MODE=auto   # Auto-detect (default)
export MOCK_MODE=true   # Force mock mode
export MOCK_MODE=false  # Force real mode
```

### Enhanced Health Check
```json
{
  "status": "healthy",
  "mode": "mock",
  "gpu_available": true,
  "gpu_name": "NVIDIA A40",
  "gpu_memory_gb": 48.0,
  "models": {
    "rfd3": {"available": false, "checkpoint_exists": false, "checkpoint_size_gb": 0},
    "rf3": {"available": false, "checkpoint_exists": false, "checkpoint_size_gb": 0},
    "proteinmpnn": {"available": false, "checkpoint_exists": false, "checkpoint_size_gb": 0}
  }
}
```

### API Compatibility
- Accepts both `contig` and `contigs` field names for RFD3 requests
- Backward compatible with existing frontend

---

## Troubleshooting

### Port 8000 not accessible externally
**Cause:** Port 8000 wasn't added to HTTP services when creating the pod.
**Fix:** Terminate pod and redeploy with `8888, 8000` in HTTP Ports field.

### "Connection refused" or CORS errors
**Cause:** Server not running or wrong URL.
**Fix:**
1. Check server is running in Jupyter: `!curl localhost:8000/health`
2. Verify URL matches pod ID exactly

### FastAPI/uvicorn not found
**Cause:** New pod doesn't have packages installed.
**Fix:** Run `!pip install fastapi uvicorn python-multipart -q` first.

### main.py not found
**Cause:** Network volume not attached or file deleted.
**Fix:** Create `main.py` using the code above, or use base64 method:
```python
import base64
code = """<paste main.py content here>"""
open("main.py", "w").write(code)
```

---

## Cost Management

| Action | Cost |
|--------|------|
| A40 GPU running | $0.40/hr |
| Network volume (50GB) | $0.005/hr |
| Pod stopped | $0 (only volume cost) |

**Tip:** Stop pod when not in use. Network volume persists all files.

---

## Pod Configuration Summary

| Setting | Value |
|---------|-------|
| GPU | A40 (48GB VRAM) recommended |
| Template | Runpod Pytorch 2.4.0 |
| HTTP Ports | `8888, 8000` |
| Network Volume | `foundry-checkpoints` |
| Mount Path | `/workspace` |

---

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/health` | GET | Health check with GPU, mode, and model status |
| `/api/rfd3/design` | POST | RFdiffusion3 structure design |
| `/api/rf3/predict` | POST | RosettaFold3 structure prediction |
| `/api/mpnn/design` | POST | ProteinMPNN sequence design |
| `/api/jobs/{job_id}` | GET | Get job status and results |
| `/api/jobs` | GET | List all jobs |
| `/api/jobs/{job_id}` | DELETE | Delete a job |

All endpoints work in both mock and real mode. Jobs are processed asynchronously with polling support.
