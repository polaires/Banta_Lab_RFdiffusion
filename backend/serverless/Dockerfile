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
    && rm -rf /var/lib/apt/lists/*

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
