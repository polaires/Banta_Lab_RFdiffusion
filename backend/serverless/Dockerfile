# RunPod Serverless Dockerfile for Foundry Protein Design
#
# Build:
#   docker build -t your-username/foundry-serverless:latest .
#   docker push your-username/foundry-serverless:latest
#
# Note: Model checkpoints should be stored on a Network Volume, not in the image

FROM runpod/pytorch:2.4.0-py3.11-cuda12.4.1-devel-ubuntu22.04

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    wget \
    curl \
    libboost-all-dev \
    libopenbabel-dev \
    && rm -rf /var/lib/apt/lists/*

# Install GNINA for protein-ligand docking and CNN scoring
# GNINA provides ~20-25% better scoring than AutoDock Vina
RUN wget -q https://github.com/gnina/gnina/releases/download/v1.1/gnina -O /usr/local/bin/gnina && \
    chmod +x /usr/local/bin/gnina

# Install xTB for high-quality conformer generation (semi-empirical quantum chemistry)
# GFN2-xTB provides excellent geometries for drug-like molecules and photoswitches
RUN mkdir -p /opt/xtb && \
    wget -q https://github.com/grimme-lab/xtb/releases/download/v6.7.1/xtb-6.7.1-linux-x86_64.tar.xz -O /tmp/xtb.tar.xz && \
    tar -xf /tmp/xtb.tar.xz -C /opt/xtb --strip-components=1 && \
    rm /tmp/xtb.tar.xz
ENV PATH="/opt/xtb/bin:$PATH"
ENV XTBHOME="/opt/xtb"

# Install CREST for conformer ensemble generation (uses xTB internally)
RUN mkdir -p /opt/crest && \
    wget -q https://github.com/crest-lab/crest/releases/download/v3.0.2/crest-gnu-12-ubuntu-latest.tar.xz -O /tmp/crest.tar.xz && \
    tar -xf /tmp/crest.tar.xz -C /opt/crest --strip-components=1 && \
    rm /tmp/crest.tar.xz && \
    chmod +x /opt/crest/crest
ENV PATH="/opt/crest:$PATH"

# Install uv (fast Python package manager that can install Python versions)
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.local/bin:$PATH"

# Install Python 3.12 using uv and create venv
RUN uv python install 3.12 && \
    uv venv /opt/venv --python 3.12
ENV PATH="/opt/venv/bin:$PATH"
ENV VIRTUAL_ENV="/opt/venv"

# Upgrade pip and install RunPod SDK
RUN uv pip install --upgrade pip runpod

# Install Foundry and dependencies (requires Python 3.12)
RUN uv pip install "rc-foundry[all]" biotite

# Install RDKit for SMILES-to-SDF conversion (needed for GNINA)
RUN uv pip install rdkit

# Install ESM3 for sequence scoring and generation
# Note: User must set HF_TOKEN env var and accept ESM3 license on HuggingFace
RUN uv pip install esm

# Copy handler code
COPY handler.py /app/
COPY inference_utils.py /app/
COPY binding_analysis.py /app/
COPY conformer_utils.py /app/

# Environment variables
# These can be overridden at runtime or in RunPod console
ENV CHECKPOINT_DIR=/runpod-volume/checkpoints
ENV FOUNDRY_CHECKPOINT_DIRS=/runpod-volume/checkpoints
ENV PYTHONUNBUFFERED=1

# HuggingFace cache for ESM3 model weights (stored on network volume for persistence)
ENV HF_HOME=/runpod-volume/hf_cache
ENV TRANSFORMERS_CACHE=/runpod-volume/hf_cache
ENV HF_HUB_CACHE=/runpod-volume/hf_cache

# Create workspace directory (for CLI tools that need it)
RUN mkdir -p /workspace && touch /workspace/.project-root

# RunPod handler entrypoint
CMD ["python", "-u", "handler.py"]
