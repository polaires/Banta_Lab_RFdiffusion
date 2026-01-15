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
# Note: python3-openbabel provides Python bindings for OpenBabel (required by PLIP)
# gfortran, BLAS/LAPACK needed to build xtb Python bindings (required by Architector)
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    wget \
    curl \
    libboost-all-dev \
    libopenbabel-dev \
    openbabel \
    python3-openbabel \
    swig \
    gfortran \
    pkg-config \
    libblas-dev \
    liblapack-dev \
    cmake \
    ninja-build \
    && rm -rf /var/lib/apt/lists/*

# Install GNINA for protein-ligand docking and CNN scoring
# GNINA provides ~20-25% better scoring than AutoDock Vina
RUN wget -q https://github.com/gnina/gnina/releases/download/v1.1/gnina -O /usr/local/bin/gnina && \
    chmod +x /usr/local/bin/gnina

# HBPLUS for hydrogen bond analysis (used by RFD3 for H-bond conditioning)
# Source code from EBI (requires signed confidentiality agreement)
# Note: Only compiling main hbplus binary; auxiliary Fortran programs (clean, accall)
# use ancient DEC format specifiers not supported by modern gfortran
COPY hbplus.zip /tmp/hbplus.zip
RUN apt-get update && apt-get install -y --no-install-recommends unzip make && \
    cd /tmp && unzip hbplus.zip && \
    cd /tmp/hbplus && \
    # Compile main HBPLUS program using Makefile
    make && \
    # Install binary and data files
    cp hbplus /usr/local/bin/ && \
    cp brkcln.par standard.data vdw.radii /usr/local/share/ && \
    chmod +x /usr/local/bin/hbplus && \
    # Cleanup - keep make installed for pynauty (Architector dependency)
    rm -rf /tmp/hbplus /tmp/hbplus.zip && \
    rm -rf /var/lib/apt/lists/*
ENV HBPLUS_PATH=/usr/local/bin/hbplus

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

# Install Miniforge3 for conda-based packages (xtb-python has no pip wheel)
# Using official miniforge installer - mambaforge is deprecated as of July 2024
ARG MINIFORGE_VERSION=24.7.1-2
ENV CONDA_DIR=/opt/conda
RUN wget --quiet https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/Miniforge3-${MINIFORGE_VERSION}-Linux-x86_64.sh -O /tmp/miniforge.sh && \
    /bin/bash /tmp/miniforge.sh -b -p ${CONDA_DIR} && \
    rm /tmp/miniforge.sh && \
    ${CONDA_DIR}/bin/conda clean --all --yes

# Install Architector for lanthanide coordination geometry via conda
# Architector provides validated 3D structure generation for f-block metal complexes
# with ~0.5 Å RMSD accuracy (Nandy et al., Nature Communications 2023)
# Installing in a separate environment to avoid conflicts with main Python
RUN ${CONDA_DIR}/bin/mamba create -y -n architector -c conda-forge python=3.11 architector && \
    ${CONDA_DIR}/bin/conda clean --all --yes && \
    echo "Architector installed via conda-forge"
ENV ARCHITECTOR_PYTHON="${CONDA_DIR}/envs/architector/bin/python"

# Install Foundry and dependencies (requires Python 3.12)
RUN uv pip install "rc-foundry[all]" biotite

# Install RDKit for SMILES-to-SDF conversion (needed for GNINA)
RUN uv pip install rdkit

# Install PyRosetta for structure refinement (FastRelax)
# Requires academic license - register at https://www.pyrosetta.org/home/licensing-pyrosetta
# Using pyrosetta-installer for automatic installation
RUN uv pip install pyrosetta-installer && \
    python -c "import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()"

# Install rdkit-to-params for generating PyRosetta ligand params from SMILES
# This allows FastRelax to properly handle ligands without pre-computed params files
RUN uv pip install rdkit-to-params

# Install ESM3 for sequence scoring and generation
# Note: User must set HF_TOKEN env var and accept ESM3 license on HuggingFace
RUN uv pip install esm

# Install PLIP for comprehensive protein-ligand interaction profiling
# Use openbabel-wheel (not openbabel) for Python 3.12 pre-built wheels
# --no-build-isolation allows plip to find installed openbabel during build
# --no-deps prevents pip from trying to install the unbuildable openbabel package
RUN pip install --only-binary=:all: "openbabel-wheel>=3.1.1.22" lxml && \
    pip install --no-build-isolation --no-deps plip

# Copy handler code
COPY handler.py /app/
COPY inference_utils.py /app/
COPY binding_analysis.py /app/
COPY conformer_utils.py /app/
COPY cleavage_utils.py /app/
COPY rosetta_utils.py /app/
COPY topology_validation.py /app/
COPY esm_utils.py /app/
COPY esmfold_utils.py /app/
COPY hotspot_detection.py /app/

# Copy lanthanide design modules (Caldwell et al. 2020 methodology)
COPY lanthanide_templates.py /app/
COPY metal_validation.py /app/
COPY tebl_analysis.py /app/
COPY hbplus_utils.py /app/
# Architector integration for physics-based coordination geometry
COPY architector_integration.py /app/
COPY geometry_validation.py /app/

# Copy utils package (SASA calculator, DSSP for secondary structure)
COPY utils/ /app/utils/

# Copy shared package (interaction analysis module)
COPY shared/ /app/shared/

# Copy test input for local testing mode
COPY test_input.json /app/test_input.json

# Environment variables
# These can be overridden at runtime or in RunPod console
ENV CHECKPOINT_DIR=/runpod-volume/checkpoints
ENV FOUNDRY_CHECKPOINT_DIRS=/runpod-volume/checkpoints
ENV PYTHONUNBUFFERED=1
ENV PYTHONPATH="/app"

# HuggingFace cache for ESM3 model weights (stored on network volume for persistence)
ENV HF_HOME=/runpod-volume/hf_cache
ENV TRANSFORMERS_CACHE=/runpod-volume/hf_cache
ENV HF_HUB_CACHE=/runpod-volume/hf_cache

# Create workspace directory (for CLI tools that need it)
RUN mkdir -p /workspace && touch /workspace/.project-root

# RunPod handler entrypoint
CMD ["python", "-u", "handler.py"]
