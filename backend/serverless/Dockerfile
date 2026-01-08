# RunPod Serverless Dockerfile for Foundry Protein Design
#
# Build:
#   docker build -t your-username/foundry-serverless:latest .
#   docker push your-username/foundry-serverless:latest
#
# Note: Model checkpoints should be stored on a Network Volume, not in the image

FROM runpod/pytorch:2.4.0-py3.11-cuda12.4.1-devel-ubuntu22.04

WORKDIR /app

# Install system dependencies including Python 3.12 (required by rc-foundry)
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    wget \
    software-properties-common \
    && add-apt-repository ppa:deadsnakes/ppa -y \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
    python3.12 \
    python3.12-venv \
    python3.12-dev \
    && rm -rf /var/lib/apt/lists/*

# Create Python 3.12 virtual environment and use it
RUN python3.12 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Upgrade pip in venv
RUN pip install --upgrade pip

# Install RunPod SDK
RUN pip install --no-cache-dir runpod

# Install Foundry and dependencies (requires Python 3.12)
RUN pip install --no-cache-dir "rc-foundry[all]" && \
    pip install --no-cache-dir biotite

# Copy handler code
COPY handler.py /app/
COPY inference_utils.py /app/

# Environment variables
# These can be overridden at runtime or in RunPod console
ENV CHECKPOINT_DIR=/runpod-volume/checkpoints
ENV FOUNDRY_CHECKPOINT_DIRS=/runpod-volume/checkpoints
ENV PYTHONUNBUFFERED=1

# Create workspace directory (for CLI tools that need it)
RUN mkdir -p /workspace && touch /workspace/.project-root

# RunPod handler entrypoint
CMD ["python", "-u", "handler.py"]
