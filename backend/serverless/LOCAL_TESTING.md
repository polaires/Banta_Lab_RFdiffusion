# Local Testing Guide for RFdiffusion Serverless

Test the exact same Docker container locally before deploying to RunPod.

## Prerequisites

### 1. Docker Desktop with WSL2 (Windows)

```powershell
# Install Docker Desktop from https://www.docker.com/products/docker-desktop/
# Enable WSL2 backend in Docker Desktop settings
```

### 2. NVIDIA Container Toolkit

For GPU support in Docker containers:

**Windows (WSL2):**
```bash
# In WSL2 terminal:
distribution=$(. /etc/os-release;echo $ID$VERSION_ID)
curl -s -L https://nvidia.github.io/nvidia-docker/gpgkey | sudo apt-key add -
curl -s -L https://nvidia.github.io/nvidia-docker/$distribution/nvidia-docker.list | sudo tee /etc/apt/sources.list.d/nvidia-docker.list

sudo apt-get update
sudo apt-get install -y nvidia-docker2
sudo systemctl restart docker
```

**Verify GPU access:**
```bash
docker run --rm --gpus all nvidia/cuda:12.4.1-base-ubuntu22.04 nvidia-smi
```

### 3. Model Checkpoints

The models are large (~15GB total). You need to download them once:

**Option A: Download via Foundry CLI (if installed locally)**
```bash
pip install rc-foundry[all]
export FOUNDRY_CHECKPOINT_DIRS=./checkpoints
foundry install base-models
```

**Option B: Copy from RunPod Network Volume**
If you already have checkpoints on RunPod, you can download them via the RunPod UI or SCP.

**Option C: Let container download on first run**
The container will use mock mode if checkpoints aren't available - useful for testing the API without GPU.

## Quick Start

### 1. Build and Run

```bash
cd backend/serverless

# Create checkpoints directory
mkdir -p checkpoints

# Build and start the container
docker-compose -f docker-compose.local.yml up --build
```

### 2. Test the API

In another terminal:

```bash
# Health check
python test_local.py health

# Test RFD3 design
python test_local.py rfd3

# Test RF3 prediction
python test_local.py rf3

# Test MPNN
python test_local.py mpnn

# Run all tests
python test_local.py all
```

Or use curl directly:

```bash
# Health check
curl http://localhost:8000/runsync -X POST \
  -H "Content-Type: application/json" \
  -d '{"input": {"task": "health"}}'

# RFD3 design
curl http://localhost:8000/runsync -X POST \
  -H "Content-Type: application/json" \
  -d '{"input": {"task": "rfd3", "contig": "100", "num_designs": 1}}'
```

## Development Workflow

### Live Code Editing

The docker-compose mounts `handler.py` and `inference_utils.py` from your local directory. Changes are reflected immediately - just restart the container:

```bash
# After editing handler.py or inference_utils.py:
docker-compose -f docker-compose.local.yml restart
```

### View Logs

```bash
# Follow logs
docker-compose -f docker-compose.local.yml logs -f

# View last 100 lines
docker-compose -f docker-compose.local.yml logs --tail=100
```

### Debug Inside Container

```bash
# Get a shell inside the running container
docker-compose -f docker-compose.local.yml exec rfdiffusion bash

# Check Python environment
python -c "from rfd3.engine import RFD3InferenceEngine; print('RFD3 OK')"

# Check checkpoints
ls -la /runpod-volume/checkpoints/
```

### Run Without GPU (Mock Mode)

If you don't have a GPU, the container will automatically run in mock mode, returning sample data. This is useful for testing the API structure:

```bash
# Run without GPU
docker-compose -f docker-compose.local.yml up --build
# The health check will show: "mode": "mock"
```

## Connecting Frontend to Local Backend

To test the frontend with your local container:

1. Edit `frontend/.env.local`:
```
RUNPOD_API_KEY=local-testing
RUNPOD_ENDPOINT_ID=local
LOCAL_BACKEND_URL=http://localhost:8000
```

2. Update the API proxy to use local backend (or create a local-only mode).

## Troubleshooting

### "CUDA out of memory"
- RFD3/RF3 require ~16GB VRAM. Reduce batch size or use a smaller GPU.

### "No GPU detected"
- Ensure NVIDIA Container Toolkit is installed
- Check `docker run --gpus all nvidia/cuda:12.4.1-base-ubuntu22.04 nvidia-smi`

### "Checkpoints not found"
- Place checkpoints in `./checkpoints/` directory
- Or run with `force=True` to download: `{"input": {"task": "download_checkpoints", "force": true}}`

### Container exits immediately
- Check logs: `docker-compose -f docker-compose.local.yml logs`
- Common issue: missing dependencies or Python version mismatch

### Slow first request
- Models are loaded on first request (~30-60s)
- Subsequent requests are fast
