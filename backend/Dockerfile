# Dockerfile for Foundry Backend on RunPod
FROM runpod/pytorch:2.4.0-py3.11-cuda12.4.1-devel-ubuntu22.04

# Set environment variables
ENV PYTHONUNBUFFERED=1
ENV FOUNDRY_CHECKPOINT_DIRS=/workspace/checkpoints

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

# Install Foundry with all models
RUN pip install "rc-foundry[all]"

# Create directories
RUN mkdir -p /workspace/checkpoints /app

# Copy application code
COPY app/ /app/

WORKDIR /app

# Expose API port
EXPOSE 8000

# Start FastAPI server
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
