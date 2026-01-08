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
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip
RUN pip install --upgrade pip

# Install RunPod SDK
RUN pip install --no-cache-dir runpod

# Install Foundry and dependencies
# Note: This is a large installation (~5GB), consider using a pre-built base image
# Install rc-foundry separately to better handle potential dependency issues
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
