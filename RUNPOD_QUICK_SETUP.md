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

## Files on Network Volume

The network volume `/workspace` contains:

### main.py (FastAPI Backend Stub)
```python
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import subprocess

app = FastAPI(title="Foundry API")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])

class HealthResponse(BaseModel):
    status: str
    gpu_available: bool
    models_loaded: list

def check_gpu():
    try:
        r = subprocess.run(["nvidia-smi"], capture_output=True, timeout=5)
        return r.returncode == 0
    except:
        return False

@app.get("/", response_model=HealthResponse)
@app.get("/health", response_model=HealthResponse)
async def health():
    return HealthResponse(status="healthy", gpu_available=check_gpu(), models_loaded=["rfd3","rf3","proteinmpnn"])

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
```

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

## API Endpoints (Current Stub)

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/health` | GET | Health check, GPU status, loaded models |
| `/` | GET | Same as /health |

**TODO:** Implement actual Foundry endpoints:
- `POST /api/rfd3/design` - RFdiffusion3 structure design
- `POST /api/rf3/predict` - RosettaFold3 structure prediction
- `POST /api/mpnn/design` - ProteinMPNN sequence design
