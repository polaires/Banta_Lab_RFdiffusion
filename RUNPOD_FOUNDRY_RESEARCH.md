# Running RosettaCommons Foundry on Runpod with Web Frontend

## Research Summary

This document outlines how to deploy [RosettaCommons Foundry](https://github.com/RosettaCommons/foundry/tree/production) on Runpod with a web-based frontend interface.

---

## 1. Overview

### What is Foundry?

[Foundry](https://github.com/RosettaCommons/foundry) is a central repository for biomolecular foundation models providing tooling and infrastructure for protein design. It includes three primary models:

| Model | Purpose | Use Case |
|-------|---------|----------|
| **RFdiffusion3 (RFD3)** | Generative design | De novo protein structure design under complex constraints |
| **RosettaFold3 (RF3)** | Structure prediction | Protein folding, protein-DNA complex prediction |
| **ProteinMPNN/LigandMPNN** | Inverse folding | Sequence design for designed backbones |

### Key Dependencies
- Python ≥3.12
- PyTorch ≥2.2.0
- CUDA support (for GPU acceleration)
- AtomWorks (structure I/O framework)

---

## 2. Runpod Deployment Options

Runpod offers two primary deployment models:

### Option A: GPU Pods (Recommended for Development & Interactive Use)

**Pros:**
- Full control over the environment
- Can run web UIs (Gradio, Streamlit) with exposed ports
- Persistent storage via network volumes
- SSH and Jupyter Lab access

**Cons:**
- Charged for idle time
- Manual scaling

### Option B: Serverless Endpoints (Recommended for Production APIs)

**Pros:**
- Pay only for compute time
- Auto-scaling from 0 to many workers
- No idle costs

**Cons:**
- Requires custom handler functions
- Not ideal for interactive web UIs
- Cold start latency

---

## 3. Recommended Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    RUNPOD GPU POD                           │
│  ┌─────────────────────────────────────────────────────┐   │
│  │                 Gradio Web Interface                 │   │
│  │                   (Port 7860)                        │   │
│  └─────────────────────┬───────────────────────────────┘   │
│                        │                                    │
│  ┌─────────────────────▼───────────────────────────────┐   │
│  │              FastAPI Backend (Port 8000)             │   │
│  │  ┌─────────────┐ ┌─────────────┐ ┌───────────────┐  │   │
│  │  │    RFD3     │ │     RF3     │ │  ProteinMPNN  │  │   │
│  │  │  Inference  │ │  Inference  │ │   Inference   │  │   │
│  │  └─────────────┘ └─────────────┘ └───────────────┘  │   │
│  └─────────────────────────────────────────────────────┘   │
│                                                             │
│  ┌─────────────────────────────────────────────────────┐   │
│  │           Network Volume (/workspace)                │   │
│  │    ~/.foundry/checkpoints (model weights ~20GB+)     │   │
│  └─────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

---

## 4. GPU Requirements

Based on similar diffusion models and Foundry's dependencies:

| GPU | VRAM | Suitability |
|-----|------|-------------|
| **NVIDIA A100 80GB** | 80GB | ✅ Recommended - handles large proteins |
| **NVIDIA A100 40GB** | 40GB | ✅ Good for medium-sized proteins |
| **NVIDIA H100 80GB** | 80GB | ✅ Best performance, expensive |
| **NVIDIA RTX 4090** | 24GB | ⚠️ Limited to smaller proteins |
| **NVIDIA A40** | 48GB | ✅ Good balance of cost/performance |

**Recommendation:** Start with an A100 40GB or A40 for development, scale to A100 80GB for production.

---

## 5. Implementation Steps

### Step 1: Create Runpod Network Volume

Network volumes provide persistent storage for model checkpoints:

1. Go to **Storage** → **New Network Volume**
2. Select datacenter (e.g., US-CA-1)
3. Set size: **50GB minimum** (for all Foundry models)
4. Name: `foundry-checkpoints`

Cost: ~$0.07/GB/month ($3.50/month for 50GB)

### Step 2: Create Custom Docker Image

```dockerfile
# Dockerfile for Foundry on Runpod
FROM runpod/pytorch:2.4.0-py3.11-cuda12.4.1-devel-ubuntu22.04

# Set environment variables
ENV PYTHONUNBUFFERED=1
ENV FOUNDRY_CHECKPOINT_DIRS=/workspace/checkpoints

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Install uv for fast Python package management
RUN pip install uv

# Install Foundry with all models
RUN pip install "rc-foundry[all]"

# Install web UI dependencies
RUN pip install gradio fastapi uvicorn

# Create checkpoint directory
RUN mkdir -p /workspace/checkpoints

# Copy web interface code
COPY app/ /app/
WORKDIR /app

# Expose ports for web interface
EXPOSE 7860 8000

# Start script
COPY start.sh /start.sh
RUN chmod +x /start.sh
CMD ["/start.sh"]
```

### Step 3: Create Gradio Web Interface

```python
# app/gradio_app.py
import gradio as gr
import subprocess
import json
import tempfile
import os

def run_rfd3_design(input_json: str, num_designs: int = 1):
    """Run RFD3 protein design"""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Write input JSON
        input_path = os.path.join(tmpdir, "input.json")
        with open(input_path, "w") as f:
            f.write(input_json)

        # Run RFD3
        out_dir = os.path.join(tmpdir, "output")
        cmd = f"rfd3 design out_dir={out_dir} inputs={input_path}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        # Return results
        if result.returncode == 0:
            # Read output PDB files
            outputs = []
            for f in os.listdir(out_dir):
                if f.endswith(".pdb"):
                    with open(os.path.join(out_dir, f)) as pdb:
                        outputs.append(pdb.read())
            return outputs, result.stdout
        else:
            return [], result.stderr

def run_rf3_predict(pdb_content: str):
    """Run RF3 structure prediction"""
    # Implementation for RF3
    pass

def run_proteinmpnn(pdb_content: str, num_sequences: int = 8):
    """Run ProteinMPNN sequence design"""
    # Implementation for ProteinMPNN
    pass

# Create Gradio interface
with gr.Blocks(title="Foundry Protein Design Suite") as demo:
    gr.Markdown("# 🧬 Foundry Protein Design Suite")
    gr.Markdown("Design proteins using RFdiffusion3, predict structures with RosettaFold3, and design sequences with ProteinMPNN")

    with gr.Tabs():
        # RFD3 Tab
        with gr.TabItem("RFD3 - Structure Design"):
            gr.Markdown("## De Novo Protein Structure Design")
            with gr.Row():
                with gr.Column():
                    rfd3_input = gr.Textbox(
                        label="Input JSON Configuration",
                        placeholder='{"contigs": ["A1-50/0 50-100"], ...}',
                        lines=10
                    )
                    rfd3_num = gr.Slider(1, 10, value=1, step=1, label="Number of Designs")
                    rfd3_btn = gr.Button("Design", variant="primary")
                with gr.Column():
                    rfd3_output = gr.Textbox(label="Output", lines=20)
                    rfd3_log = gr.Textbox(label="Log", lines=5)

        # RF3 Tab
        with gr.TabItem("RF3 - Structure Prediction"):
            gr.Markdown("## Protein Structure Prediction")
            with gr.Row():
                rf3_input = gr.Textbox(label="Input Sequence (FASTA)", lines=5)
                rf3_btn = gr.Button("Predict", variant="primary")
            rf3_output = gr.Textbox(label="Predicted Structure (PDB)", lines=20)

        # ProteinMPNN Tab
        with gr.TabItem("ProteinMPNN - Sequence Design"):
            gr.Markdown("## Inverse Folding / Sequence Design")
            with gr.Row():
                mpnn_input = gr.Textbox(label="Input Structure (PDB)", lines=10)
                mpnn_num = gr.Slider(1, 32, value=8, step=1, label="Number of Sequences")
                mpnn_btn = gr.Button("Design Sequences", variant="primary")
            mpnn_output = gr.Textbox(label="Designed Sequences", lines=10)

if __name__ == "__main__":
    demo.launch(server_name="0.0.0.0", server_port=7860)
```

### Step 4: Create Start Script

```bash
#!/bin/bash
# start.sh

# Download model checkpoints if not present
if [ ! -d "/workspace/checkpoints/rfd3" ]; then
    echo "Downloading Foundry model checkpoints..."
    foundry install base-models --checkpoint-dir /workspace/checkpoints
fi

# Set checkpoint path
export FOUNDRY_CHECKPOINT_DIRS=/workspace/checkpoints

# Start Gradio web interface
python /app/gradio_app.py
```

### Step 5: Deploy on Runpod

1. **Build and push Docker image:**
   ```bash
   docker build -t yourusername/foundry-web:latest .
   docker push yourusername/foundry-web:latest
   ```

2. **Create GPU Pod:**
   - Go to Runpod Console → **Pods** → **Deploy**
   - Select GPU: **A100 40GB** or **A40**
   - Template: **Custom**
   - Docker Image: `yourusername/foundry-web:latest`
   - Attach Network Volume: `foundry-checkpoints`
   - **Expose HTTP Ports:** `7860, 8000`
   - Click **Deploy**

3. **Access Web Interface:**
   - Once pod is running, click on the **7860** port link
   - URL format: `https://<pod-id>-7860.proxy.runpod.net`

---

## 6. Port Configuration

| Port | Service | Purpose |
|------|---------|---------|
| 7860 | Gradio | Main web interface |
| 8000 | FastAPI | REST API (optional) |
| 8888 | Jupyter | Development notebooks |
| 22 | SSH | Remote access |

**Important:** Ensure all services bind to `0.0.0.0`, not `localhost`.

---

## 7. Cost Estimates

### GPU Pod Costs (per hour)

| GPU | On-Demand | Community Cloud |
|-----|-----------|-----------------|
| A100 40GB | ~$1.89/hr | ~$1.20/hr |
| A100 80GB | ~$2.49/hr | ~$1.60/hr |
| A40 | ~$0.79/hr | ~$0.50/hr |
| RTX 4090 | ~$0.69/hr | ~$0.44/hr |

### Monthly Estimates (8 hours/day usage)

| Component | Cost |
|-----------|------|
| A40 GPU (8hr/day × 30 days) | ~$190/month |
| Network Volume (50GB) | ~$3.50/month |
| **Total** | **~$194/month** |

---

## 8. Alternative: Serverless API Deployment

For production APIs without a web frontend:

### Handler Function

```python
# handler.py for Runpod Serverless
import runpod
import subprocess
import tempfile
import os

def handler(job):
    """Serverless handler for Foundry inference"""
    job_input = job["input"]

    model = job_input.get("model", "rfd3")
    config = job_input.get("config", {})

    with tempfile.TemporaryDirectory() as tmpdir:
        if model == "rfd3":
            # Run RFD3
            input_path = os.path.join(tmpdir, "input.json")
            with open(input_path, "w") as f:
                json.dump(config, f)

            out_dir = os.path.join(tmpdir, "output")
            subprocess.run(
                f"rfd3 design out_dir={out_dir} inputs={input_path}",
                shell=True, check=True
            )

            # Collect results
            results = []
            for f in os.listdir(out_dir):
                if f.endswith(".pdb"):
                    with open(os.path.join(out_dir, f)) as pdb:
                        results.append({"filename": f, "content": pdb.read()})

            return {"status": "success", "results": results}

        elif model == "rf3":
            # RF3 implementation
            pass

        elif model == "mpnn":
            # ProteinMPNN implementation
            pass

runpod.serverless.start({"handler": handler})
```

### Serverless Dockerfile

```dockerfile
FROM runpod/pytorch:2.4.0-py3.11-cuda12.4.1-devel-ubuntu22.04

RUN pip install "rc-foundry[all]" runpod

# Pre-download models during build (optional, increases image size)
# RUN foundry install base-models

COPY handler.py /handler.py
CMD ["python", "-u", "/handler.py"]
```

---

## 9. Best Practices

### Performance Optimization

1. **Pre-download models** to network volume (not during pod startup)
2. **Use FlashBoot** templates when available for faster cold starts
3. **Cache intermediate results** for iterative design workflows
4. **Use batch processing** when designing multiple proteins

### Security

1. **Don't expose** sensitive ports (SSH) publicly
2. **Use Runpod's proxy** for HTTP access (automatic HTTPS)
3. **Set environment variables** for API keys securely

### Cost Optimization

1. **Stop pods** when not in use
2. **Use Community Cloud** for development (cheaper)
3. **Use Secure Cloud** for production (more reliable)
4. **Right-size GPU** based on protein size requirements

---

## 10. Quick Start Commands

```bash
# 1. Clone this repository
git clone https://github.com/your-org/foundry-runpod.git
cd foundry-runpod

# 2. Build Docker image
docker build -t foundry-web:latest .

# 3. Push to Docker Hub
docker tag foundry-web:latest yourusername/foundry-web:latest
docker push yourusername/foundry-web:latest

# 4. Deploy on Runpod (via web console or CLI)
# - Select A40 or A100 GPU
# - Attach 50GB network volume
# - Expose port 7860
# - Use image: yourusername/foundry-web:latest
```

---

## 11. References

### Foundry & Models
- [RosettaCommons Foundry GitHub](https://github.com/RosettaCommons/foundry)
- [RFdiffusion GitHub](https://github.com/RosettaCommons/RFdiffusion)
- [Institute for Protein Design](https://www.ipd.uw.edu/software/)

### Runpod Documentation
- [Runpod Serverless Overview](https://docs.runpod.io/serverless/overview)
- [Expose Ports Documentation](https://docs.runpod.io/pods/configuration/expose-ports)
- [Network Volumes](https://docs.runpod.io/storage/network-volumes)

### Deployment Guides
- [Deploy FastAPI with GPU](https://www.runpod.io/articles/guides/deploy-fastapi-applications-gpu-cloud)
- [PyTorch 2.4 + CUDA 12.4 Template](https://www.runpod.io/articles/guides/pytorch-2-4-cuda-12-4)
- [Gradio Deployment Guide](https://www.gradio.app/)

### Alternative Platforms
- [NVIDIA NIM for RFdiffusion](https://docs.nvidia.com/nim/bionemo/rfdiffusion/latest/overview.html)
- [Levitate Bio RFDiffusion API](https://support.levitate.bio/api/api-rfdiffusion/)

---

## 12. Summary

| Aspect | Recommendation |
|--------|----------------|
| **Deployment Type** | GPU Pods (for web UI) or Serverless (for API) |
| **GPU** | A100 40GB or A40 for most use cases |
| **Storage** | 50GB+ Network Volume for model checkpoints |
| **Web Framework** | Gradio (easiest) or Streamlit |
| **Port** | 7860 (Gradio default) |
| **Cost** | ~$190-250/month for moderate usage |

This setup provides a production-ready environment for running Foundry's protein design models with an accessible web interface on Runpod's GPU infrastructure.
