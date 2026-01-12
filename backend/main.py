"""
Foundry Protein Design API - Unified Backend (v2.0)
FastAPI backend for RFdiffusion3, RosettaFold3, and ProteinMPNN

Features:
- Python API integration for all models (matching IPD official pipeline)
- RF3 confidence metrics (pLDDT, PAE, pTM, ranking_score)
- MPNN model type selection (ligand_mpnn/protein_mpnn)
- RMSD validation endpoint for designability assessment
- Seed control for reproducible generation
- CIF export support
- Falls back to mock mode when models unavailable

Deploy to RunPod:
1. Copy this file to /workspace/main.py on your RunPod pod
2. Run: uvicorn main:app --host 0.0.0.0 --port 8000
"""

from fastapi import FastAPI, HTTPException, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, model_validator, Field
from typing import Optional, List, Dict, Any, Literal
import subprocess
import tempfile
import uuid
import os
import json
import asyncio
from datetime import datetime
from contextlib import asynccontextmanager

# ============== Configuration ==============

MOCK_MODE = os.environ.get("MOCK_MODE", "auto").lower()
CHECKPOINT_DIR = os.environ.get("FOUNDRY_CHECKPOINT_DIRS", "/workspace/checkpoints")

import sys
VENV_BIN_DIR = os.path.dirname(sys.executable)

KNOWN_VENV_PATHS = [
    VENV_BIN_DIR,
    "/workspace/foundry_env/bin",
    os.path.expanduser("~/foundry_env/bin"),
]

CLI_PREFIX = ""
for venv_path in KNOWN_VENV_PATHS:
    if os.path.exists(os.path.join(venv_path, "rfd3")):
        CLI_PREFIX = venv_path
        break

# Global state
FOUNDRY_AVAILABLE = False
GPU_INFO: Dict[str, Any] = {}

# In-memory job storage
jobs: Dict[str, Dict[str, Any]] = {}

# Store atom arrays for cross-panel data flow
stored_structures: Dict[str, Any] = {}


def get_cli_path(cmd: str) -> str:
    if CLI_PREFIX:
        return os.path.join(CLI_PREFIX, cmd)
    return cmd


def check_foundry_available() -> bool:
    """Check if Foundry Python modules are importable"""
    try:
        from rfd3.engine import RFD3InferenceEngine
        print("[STARTUP] Found rfd3.engine module")
        return True
    except ImportError:
        pass

    # Fallback: check CLI tools
    for cmd in ["rfd3", "rf3", "mpnn"]:
        try:
            cli_path = get_cli_path(cmd)
            result = subprocess.run([cli_path, "--help"], capture_output=True, timeout=10)
            if result.returncode == 0:
                print(f"[STARTUP] Found CLI tool: {cli_path}")
                return True
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
    return False


def get_gpu_info() -> Dict[str, Any]:
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=name,memory.total", "--format=csv,noheader,nounits"],
            capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0:
            parts = result.stdout.strip().split(", ")
            return {
                "available": True,
                "name": parts[0] if parts else "Unknown",
                "memory_gb": float(parts[1]) / 1024 if len(parts) > 1 else 0,
            }
    except Exception:
        pass
    return {"available": False, "name": None, "memory_gb": None}


def get_model_status() -> Dict[str, Dict[str, Any]]:
    status = {}
    for model, subdir in [("rfd3", "rfd3"), ("rf3", "rf3"), ("proteinmpnn", "mpnn")]:
        path = os.path.join(CHECKPOINT_DIR, subdir)
        exists = os.path.exists(path)
        size_gb = 0
        if exists:
            try:
                total = sum(
                    os.path.getsize(os.path.join(dp, f))
                    for dp, _, filenames in os.walk(path)
                    for f in filenames
                )
                size_gb = total / (1024**3)
            except Exception:
                pass
        status[model] = {
            "available": FOUNDRY_AVAILABLE,
            "checkpoint_exists": exists,
            "checkpoint_size_gb": round(size_gb, 2),
        }
    return status


@asynccontextmanager
async def lifespan(app: FastAPI):
    global FOUNDRY_AVAILABLE, GPU_INFO

    if MOCK_MODE == "true":
        FOUNDRY_AVAILABLE = False
        print("[STARTUP] Mock mode forced via MOCK_MODE=true")
    elif MOCK_MODE == "false":
        FOUNDRY_AVAILABLE = True
        print("[STARTUP] Real mode forced via MOCK_MODE=false")
    else:
        FOUNDRY_AVAILABLE = check_foundry_available()
        print(f"[STARTUP] Auto-detected Foundry available: {FOUNDRY_AVAILABLE}")

    GPU_INFO = get_gpu_info()
    print(f"[STARTUP] GPU: {GPU_INFO}")
    print(f"[STARTUP] Checkpoint dir: {CHECKPOINT_DIR}")
    print(f"[STARTUP] CLI prefix: {CLI_PREFIX or '(system PATH)'}")

    yield


app = FastAPI(
    title="Foundry Protein Design API",
    description="Backend API for RFdiffusion3, RosettaFold3, and ProteinMPNN (v2.0 - Python API)",
    version="2.0.0",
    lifespan=lifespan,
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# ============== Request/Response Models ==============

class HealthResponse(BaseModel):
    status: str
    mode: str
    gpu_available: bool
    gpu_name: Optional[str] = None
    gpu_memory_gb: Optional[float] = None
    models: Dict[str, Dict[str, Any]]


class SymmetryConfig(BaseModel):
    """Symmetry configuration for oligomeric design"""
    id: str  # C2, C3, C4, D2, D3, T, O, I
    is_unsym_motif: Optional[str] = None  # Chains to exclude from symmetry
    is_symmetric_motif: bool = True


class CFGConfig(BaseModel):
    """Classifier-Free Guidance configuration"""
    enabled: bool = False
    scale: float = Field(default=1.5, ge=0.5, le=3.0)
    t_max: float = Field(default=0.8, ge=0.0, le=1.0)
    features: Optional[List[str]] = None  # active_donor, active_acceptor, ref_atomwise_rasa


class RFD3Request(BaseModel):
    """RFdiffusion3 design request with advanced options"""
    contig: Optional[str] = None
    contigs: Optional[str] = None
    num_designs: int = Field(default=1, ge=1, le=10)
    pdb_content: Optional[str] = None
    seed: Optional[int] = None  # For reproducibility
    config: Optional[Dict[str, Any]] = None

    # Phase 1: Hotspot residues - dict mapping residue ranges to atom selection
    # e.g., {"A15-20": "ALL", "A25": "TIP", "A30-35": "BKBN"}
    select_hotspots: Optional[Dict[str, str]] = None

    # Phase 1: Partial diffusion - noise level in Angstroms (5-20 recommended)
    partial_t: Optional[float] = Field(default=None, ge=0.0, le=30.0)

    # Phase 1: Diffusion sampling parameters
    num_timesteps: Optional[int] = Field(default=None, ge=50, le=500)
    step_scale: Optional[float] = Field(default=None, ge=0.5, le=3.0)
    noise_scale: Optional[float] = Field(default=None, ge=0.9, le=1.1)
    gamma_0: Optional[float] = Field(default=None, ge=0.1, le=1.0)

    # Phase 2: Symmetry design - for oligomeric proteins
    symmetry: Optional[SymmetryConfig] = None

    # Phase 2: Ligand binding - chemical component ID (e.g., "ATP", "NAD", "ZN")
    ligand: Optional[str] = None

    # Phase 2: Classifier-Free Guidance
    cfg: Optional[CFGConfig] = None

    @model_validator(mode="after")
    def normalize_contig(self):
        if not self.contig and not self.contigs:
            raise ValueError("Either 'contig' or 'contigs' must be provided")
        if not self.contig and self.contigs:
            self.contig = self.contigs
        return self


class RF3Request(BaseModel):
    """RosettaFold3 prediction request"""
    sequence: str
    name: Optional[str] = "prediction"
    pdb_content: Optional[str] = None  # For structure-based input
    config: Optional[Dict[str, Any]] = None


class MPNNRequest(BaseModel):
    """ProteinMPNN sequence design request with model type selection"""
    pdb_content: str
    num_sequences: int = Field(default=8, ge=1, le=100)
    temperature: float = Field(default=0.1, ge=0.01, le=2.0)
    model_type: Literal["ligand_mpnn", "protein_mpnn"] = "ligand_mpnn"
    remove_waters: bool = True
    config: Optional[Dict[str, Any]] = None


class RMSDRequest(BaseModel):
    """RMSD validation request"""
    pdb_content_1: str  # Reference structure (e.g., RFD3 output)
    pdb_content_2: str  # Comparison structure (e.g., RF3 refolded)
    backbone_only: bool = True


class ExportRequest(BaseModel):
    """Structure export request"""
    pdb_content: str
    format: Literal["pdb", "cif"] = "cif"
    include_confidences: bool = False
    confidences: Optional[Dict[str, Any]] = None


class FetchPDBRequest(BaseModel):
    """Request to fetch PDB from RCSB"""
    pdb_id: str
    format: Literal["pdb", "cif"] = "pdb"
    use_cache: bool = True


class MetalBindingAnalysisRequest(BaseModel):
    """Request to analyze metal binding site"""
    pdb_content: Optional[str] = None
    pdb_id: Optional[str] = None  # Alternative: fetch from RCSB
    metal_chain: str = "A"
    metal_residue: str  # e.g., "FE", "ZN", "CA"
    metal_resnum: int
    distance_cutoff: float = Field(default=3.0, ge=2.0, le=5.0)


class AIRecommendRequest(BaseModel):
    """Request for AI-assisted parameter recommendation"""
    pdb_content: Optional[str] = None
    pdb_id: Optional[str] = None
    metal_chain: str = "A"
    metal_residue: str  # Current metal (e.g., "FE")
    metal_resnum: int
    target_metal: str  # Target metal (e.g., "TB")
    user_description: Optional[str] = None  # Natural language goal
    include_sasa: bool = True
    include_secondary_structure: bool = True


class JobResponse(BaseModel):
    job_id: str
    status: str
    message: str


class JobStatus(BaseModel):
    job_id: str
    status: str
    created_at: str
    completed_at: Optional[str] = None
    result: Optional[Dict[str, Any]] = None
    error: Optional[str] = None


# ============== Utility Functions ==============

def atom_array_to_pdb(atom_array) -> str:
    """Convert biotite AtomArray to PDB string"""
    from biotite.structure.io.pdb import PDBFile
    import io

    pdb_file = PDBFile()
    pdb_file.set_structure(atom_array)
    buf = io.StringIO()
    pdb_file.write(buf)
    return buf.getvalue()


def atom_array_to_cif(atom_array) -> str:
    """Convert biotite AtomArray to CIF string"""
    try:
        from atomworks.io.utils.io_utils import to_cif_file
        import tempfile
        import os

        with tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as f:
            temp_path = f.name

        to_cif_file(atom_array, temp_path)
        with open(temp_path) as f:
            content = f.read()
        os.unlink(temp_path)
        return content
    except ImportError:
        # Fallback to biotite CIF writer
        from biotite.structure.io.pdbx import PDBxFile, set_structure
        import io

        pdbx_file = PDBxFile()
        set_structure(pdbx_file, atom_array)
        buf = io.StringIO()
        pdbx_file.write(buf)
        return buf.getvalue()


def pdb_to_atom_array(content: str):
    """Convert PDB or CIF string to biotite AtomArray.

    Automatically detects format based on content:
    - CIF files start with 'data_'
    - PDB files start with HEADER, ATOM, etc.
    """
    import io

    # Detect format
    content_stripped = content.strip()
    if content_stripped.startswith('data_'):
        # CIF format (mmCIF)
        from biotite.structure.io.pdbx import CIFFile
        cif_file = CIFFile.read(io.StringIO(content))
        # Get the first block
        block = list(cif_file.values())[0] if cif_file else None
        if block is None:
            raise ValueError("Empty CIF file")
        from biotite.structure.io.pdbx import get_structure
        return get_structure(block, model=1)
    else:
        # PDB format
        from biotite.structure.io.pdb import PDBFile
        pdb_file = PDBFile.read(io.StringIO(content))
        return pdb_file.get_structure(model=1)


def extract_sequence_from_atom_array(atom_array) -> str:
    """Extract amino acid sequence from AtomArray"""
    from biotite.structure import get_residue_starts
    from biotite.sequence import ProteinSequence

    res_starts = get_residue_starts(atom_array)
    seq = ''.join(
        ProteinSequence.convert_letter_3to1(res_name)
        for res_name in atom_array.res_name[res_starts]
        if res_name in ProteinSequence._dict_3to1
    )
    return seq


# ============== Mock Implementations ==============

def generate_mock_pdb(length: int = 100) -> str:
    lines = ["HEADER    MOCK PROTEIN STRUCTURE"]
    for i in range(1, length + 1):
        x, y, z = i * 0.5, i * 0.3, i * 0.2
        lines.append(
            f"ATOM  {i:5d}  CA  ALA A{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
        )
    lines.append("END")
    return "\n".join(lines)


def generate_mock_fasta(num_seq: int = 8, length: int = 100) -> str:
    import random
    aa = "ACDEFGHIKLMNPQRSTVWY"
    lines = []
    for i in range(num_seq):
        seq = "".join(random.choices(aa, k=length))
        lines.append(f">design_{i+1}")
        lines.append(seq)
    return "\n".join(lines)


def generate_mock_confidences(length: int = 100) -> Dict[str, Any]:
    """Generate mock confidence metrics"""
    import random
    return {
        "summary_confidences": {
            "overall_plddt": round(random.uniform(0.7, 0.95), 3),
            "overall_pae": round(random.uniform(2.0, 8.0), 2),
            "overall_pde": round(random.uniform(0.1, 0.3), 3),
            "ptm": round(random.uniform(0.6, 0.9), 3),
            "iptm": None,
            "ranking_score": round(random.uniform(0.5, 0.85), 3),
            "has_clash": False,
        },
        "per_residue_plddt": [round(random.uniform(0.6, 1.0), 2) for _ in range(length)],
    }


# ============== RFD3 Job Processing ==============

async def process_rfd3_job(job_id: str, request: RFD3Request):
    jobs[job_id]["status"] = "running"

    if not FOUNDRY_AVAILABLE:
        await process_rfd3_mock(job_id, request)
    else:
        await process_rfd3_real(job_id, request)


async def process_rfd3_mock(job_id: str, request: RFD3Request):
    await asyncio.sleep(3)

    try:
        length = int(request.contig) if request.contig and request.contig.isdigit() else 100
    except Exception:
        length = 100

    pdb_content = generate_mock_pdb(length)

    jobs[job_id]["status"] = "completed"
    jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
    jobs[job_id]["result"] = {
        "designs": [{"filename": "design_1.pdb", "content": pdb_content}],
        "mode": "mock",
    }

    # Store for data flow
    stored_structures[job_id] = {"pdb": pdb_content, "type": "rfd3"}


async def process_rfd3_real(job_id: str, request: RFD3Request):
    try:
        from rfd3.engine import RFD3InferenceConfig, RFD3InferenceEngine

        # Set seed for reproducibility
        if request.seed is not None:
            from lightning.fabric import seed_everything
            seed_everything(request.seed)
            print(f"[RFD3] Set seed to {request.seed}")

        contig_str = request.contig or "100"
        is_simple_length = contig_str.replace("-", "").isdigit()

        # Build specification dict
        if is_simple_length:
            if "-" in contig_str:
                spec = {"length": contig_str}
            else:
                spec = {"length": int(contig_str)}
        else:
            spec = {"contig": contig_str}

        # Phase 1: Add hotspot residues to specification
        if request.select_hotspots:
            spec["select_hotspots"] = request.select_hotspots
            print(f"[RFD3] Added hotspots: {request.select_hotspots}")

        # Phase 1: Add partial diffusion to specification
        if request.partial_t is not None:
            spec["partial_t"] = request.partial_t
            print(f"[RFD3] Partial diffusion enabled: {request.partial_t}Å")

        # Phase 2: Add symmetry configuration
        if request.symmetry:
            symmetry_dict = {"id": request.symmetry.id}
            if request.symmetry.is_unsym_motif:
                symmetry_dict["is_unsym_motif"] = request.symmetry.is_unsym_motif
            symmetry_dict["is_symmetric_motif"] = request.symmetry.is_symmetric_motif
            spec["symmetry"] = symmetry_dict
            print(f"[RFD3] Symmetry enabled: {request.symmetry.id}")

        # Phase 2: Add ligand binding
        if request.ligand:
            spec["ligand"] = request.ligand
            print(f"[RFD3] Ligand binding enabled: {request.ligand}")

        # Build inference sampler config for diffusion parameters
        sampler_config = {}
        if request.num_timesteps is not None:
            sampler_config["num_timesteps"] = request.num_timesteps
        if request.step_scale is not None:
            sampler_config["step_scale"] = request.step_scale
        if request.noise_scale is not None:
            sampler_config["noise_scale"] = request.noise_scale
        if request.gamma_0 is not None:
            sampler_config["gamma_0"] = request.gamma_0

        # Phase 2: Add symmetry sampler kind
        if request.symmetry:
            sampler_config["kind"] = "symmetry"

        # Phase 2: Add Classifier-Free Guidance
        if request.cfg and request.cfg.enabled:
            sampler_config["use_classifier_free_guidance"] = True
            sampler_config["cfg_scale"] = request.cfg.scale
            sampler_config["cfg_t_max"] = request.cfg.t_max
            if request.cfg.features:
                sampler_config["cfg_features"] = request.cfg.features
            print(f"[RFD3] CFG enabled: scale={request.cfg.scale}, features={request.cfg.features}")

        if sampler_config:
            print(f"[RFD3] Sampler config: {sampler_config}")

        # Build the inference config
        # Symmetry requires batch_size=1 due to memory constraints
        batch_size = 1 if request.symmetry else request.num_designs

        config_kwargs = {
            "specification": spec,
            "diffusion_batch_size": batch_size,
        }

        # High-order symmetry requires low memory mode
        if request.symmetry and request.symmetry.id in ["T", "O", "I"]:
            config_kwargs["low_memory_mode"] = True
            print(f"[RFD3] Low memory mode enabled for {request.symmetry.id} symmetry")

        # Add sampler config if any parameters were set
        if sampler_config:
            config_kwargs["inference_sampler"] = sampler_config

        config = RFD3InferenceConfig(**config_kwargs)

        model = RFD3InferenceEngine(**config)
        outputs_dict = model.run(
            inputs=None,
            out_dir=None,
            n_batches=1,
        )

        outputs = []
        first_atom_array = None

        if outputs_dict:
            for key, result_list in outputs_dict.items():
                items = result_list if isinstance(result_list, list) else [result_list]

                for idx, item in enumerate(items):
                    pdb_content = None
                    cif_content = None
                    filename = f"{key}_{idx}.pdb" if len(items) > 1 else f"{key}.pdb"

                    if hasattr(item, 'atom_array'):
                        try:
                            atom_array = item.atom_array
                            if first_atom_array is None:
                                first_atom_array = atom_array

                            pdb_content = atom_array_to_pdb(atom_array)
                            try:
                                cif_content = atom_array_to_cif(atom_array)
                            except Exception:
                                pass
                        except Exception as e:
                            print(f"Error converting atom_array: {e}")

                    if pdb_content:
                        output_item = {"filename": filename, "content": pdb_content}
                        if cif_content:
                            output_item["cif_content"] = cif_content
                        outputs.append(output_item)

        jobs[job_id]["status"] = "completed"
        jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
        jobs[job_id]["result"] = {
            "designs": outputs,
            "mode": "real",
            "seed": request.seed,
        }

        # Store first structure for data flow
        if outputs:
            stored_structures[job_id] = {
                "pdb": outputs[0]["content"],
                "cif": outputs[0].get("cif_content"),
                "type": "rfd3",
            }

    except ImportError as e:
        jobs[job_id]["status"] = "failed"
        jobs[job_id]["error"] = f"Foundry not installed: {e}"
    except Exception as e:
        import traceback
        jobs[job_id]["status"] = "failed"
        jobs[job_id]["error"] = f"{str(e)}\n{traceback.format_exc()}"


# ============== RF3 Job Processing (Python API with Confidences) ==============

async def process_rf3_job(job_id: str, request: RF3Request):
    jobs[job_id]["status"] = "running"

    if not FOUNDRY_AVAILABLE:
        await process_rf3_mock(job_id, request)
    else:
        await process_rf3_real(job_id, request)


async def process_rf3_mock(job_id: str, request: RF3Request):
    await asyncio.sleep(3)

    length = len(request.sequence)
    pdb_content = generate_mock_pdb(length)
    confidences = generate_mock_confidences(length)

    jobs[job_id]["status"] = "completed"
    jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
    jobs[job_id]["result"] = {
        "predictions": [{
            "filename": f"{request.name}.pdb",
            "content": pdb_content,
        }],
        "confidences": confidences,
        "mode": "mock",
    }

    stored_structures[job_id] = {"pdb": pdb_content, "type": "rf3"}


async def process_rf3_real(job_id: str, request: RF3Request):
    """Real RF3 implementation using Python API with confidence metrics"""
    try:
        from rf3.inference_engines.rf3 import RF3InferenceEngine
        from rf3.utils.inference import InferenceInput

        # For sequence-only input, use CLI approach (Python API needs atom_array)
        if not request.pdb_content:
            print("[RF3] Sequence-only input, using CLI approach")
            await process_rf3_cli(job_id, request)
            return

        # Initialize RF3 engine for structure input
        inference_engine = RF3InferenceEngine(ckpt_path='rf3', verbose=False)

        # Structure-based input
        atom_array = pdb_to_atom_array(request.pdb_content)
        input_data = InferenceInput.from_atom_array(
            atom_array,
            example_id=request.name or "prediction"
        )

        # Run inference
        rf3_outputs = inference_engine.run(inputs=input_data)

        # Process outputs
        outputs = []
        confidences = None
        first_pdb = None

        example_id = request.name or "prediction"
        if example_id in rf3_outputs:
            results = rf3_outputs[example_id]

            for idx, rf3_output in enumerate(results):
                # Extract structure
                if hasattr(rf3_output, 'atom_array'):
                    pdb_content = atom_array_to_pdb(rf3_output.atom_array)
                    if first_pdb is None:
                        first_pdb = pdb_content

                    output_item = {
                        "filename": f"{example_id}_{idx}.pdb" if len(results) > 1 else f"{example_id}.pdb",
                        "content": pdb_content,
                    }

                    # Try to add CIF
                    try:
                        output_item["cif_content"] = atom_array_to_cif(rf3_output.atom_array)
                    except Exception:
                        pass

                    outputs.append(output_item)

                # Extract confidence metrics (from first output)
                if idx == 0 and confidences is None:
                    confidences = {}

                    # Summary confidences
                    if hasattr(rf3_output, 'summary_confidences') and rf3_output.summary_confidences:
                        sc = rf3_output.summary_confidences
                        confidences["summary_confidences"] = {
                            "overall_plddt": float(sc.get('overall_plddt', 0)),
                            "overall_pae": float(sc.get('overall_pae', 0)),
                            "overall_pde": float(sc.get('overall_pde', 0)) if 'overall_pde' in sc else None,
                            "ptm": float(sc.get('ptm', 0)),
                            "iptm": float(sc.get('iptm')) if sc.get('iptm') is not None else None,
                            "ranking_score": float(sc.get('ranking_score', 0)),
                            "has_clash": bool(sc.get('has_clash', False)),
                        }

                    # Per-residue confidences
                    if hasattr(rf3_output, 'confidences') and rf3_output.confidences:
                        conf = rf3_output.confidences
                        if 'atom_plddts' in conf:
                            # Convert to per-residue by averaging
                            atom_plddts = list(conf['atom_plddts'])
                            confidences["per_residue_plddt"] = [round(float(p), 3) for p in atom_plddts[:500]]  # Limit size

                        if 'pae' in conf:
                            # PAE matrix - limit size for JSON
                            pae = conf['pae']
                            pae_list = [[round(float(x), 2) for x in row[:100]] for row in pae[:100]]
                            confidences["pae_matrix"] = pae_list

        jobs[job_id]["status"] = "completed"
        jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
        jobs[job_id]["result"] = {
            "predictions": outputs,
            "confidences": confidences,
            "mode": "real",
        }

        if first_pdb:
            stored_structures[job_id] = {"pdb": first_pdb, "type": "rf3", "confidences": confidences}

    except ImportError as e:
        # Fallback to CLI if Python API not available
        print(f"[RF3] Python API not available ({e}), falling back to CLI")
        await process_rf3_cli(job_id, request)
    except Exception as e:
        import traceback
        jobs[job_id]["status"] = "failed"
        jobs[job_id]["error"] = f"{str(e)}\n{traceback.format_exc()}"


def ensure_rf3_checkpoint():
    """Ensure RF3 checkpoint is available, attempting download if needed"""
    rf3_ckpt_name = "rf3_foundry_01_24_latest.ckpt"
    rf3_ckpt_target = f"/workspace/{rf3_ckpt_name}"

    # Check if checkpoint already exists
    if os.path.exists(rf3_ckpt_target):
        print(f"[RF3] Checkpoint exists at {rf3_ckpt_target}")
        return True

    # Check common locations and create symlink
    # Note: The actual checkpoint might have "_remapped" suffix
    possible_locations = [
        f"/workspace/checkpoints/rf3_foundry_01_24_latest_remapped.ckpt",  # Most common
        f"/workspace/checkpoints/{rf3_ckpt_name}",
        f"/workspace/checkpoints/rf3/{rf3_ckpt_name}",
        f"/workspace/checkpoints/rf3/rf3_foundry_01_24_latest_remapped.ckpt",
        f"/workspace/foundry_checkpoints/rf3/{rf3_ckpt_name}",
        f"/root/.cache/foundry/rf3/{rf3_ckpt_name}",
        f"/root/.foundry/checkpoints/rf3/{rf3_ckpt_name}",
    ]

    for ckpt_path in possible_locations:
        if os.path.exists(ckpt_path):
            try:
                os.symlink(ckpt_path, rf3_ckpt_target)
                print(f"[RF3] Created symlink: {rf3_ckpt_target} -> {ckpt_path}")
                return True
            except FileExistsError:
                return True
            except Exception as e:
                print(f"[RF3] Could not create symlink: {e}")

    # Try to trigger checkpoint download via Python API initialization
    # This works because RFD3InferenceEngine auto-downloads, RF3 might too
    print("[RF3] Checkpoint not found, attempting to trigger download via Python API...")
    try:
        from rf3.inference_engines.rf3 import RF3InferenceEngine
        # Just initializing the engine should trigger checkpoint download
        # Use a try block since this might fail or take time
        engine = RF3InferenceEngine(ckpt_path='rf3', verbose=True)
        print("[RF3] Successfully initialized RF3 engine (checkpoint should be downloaded)")
        return True
    except Exception as e:
        print(f"[RF3] Python API initialization failed: {e}")

    # Try using Foundry CLI to install checkpoint
    print("[RF3] Attempting to install via Foundry CLI...")
    try:
        foundry_cli = get_cli_path("foundry")
        result = subprocess.run(
            f"{foundry_cli} install rf3",
            shell=True, capture_output=True, text=True, timeout=600,
            cwd="/workspace"
        )
        if result.returncode == 0:
            print(f"[RF3] Foundry install succeeded: {result.stdout}")
            return True
        else:
            print(f"[RF3] Foundry install failed: {result.stderr}")
    except Exception as e:
        print(f"[RF3] Foundry CLI install failed: {e}")

    # Debug: list available directories
    print("[RF3] Checkpoint not found. Debugging directory contents...")
    for check_dir in ["/workspace/checkpoints", "/workspace", "/root/.cache/foundry"]:
        if os.path.exists(check_dir):
            try:
                contents = os.listdir(check_dir)
                print(f"[RF3] {check_dir}: {contents[:20]}")  # Limit output
            except Exception:
                pass

    return False


async def process_rf3_cli(job_id: str, request: RF3Request):
    """CLI implementation for RF3 using JSON config for sequence input"""
    try:
        # Ensure checkpoint is available before running
        if not ensure_rf3_checkpoint():
            jobs[job_id]["status"] = "failed"
            jobs[job_id]["error"] = (
                "RF3 checkpoint not found and could not be downloaded. "
                "Please run 'foundry install rf3' or 'foundry install base-models' "
                "on your RunPod pod to install the required checkpoint."
            )
            return

        with tempfile.TemporaryDirectory() as tmpdir:
            # RF3 CLI uses JSON config with 'components' array
            # Format discovered from atomworks InferenceInput.from_json_dict()
            config_path = os.path.join(tmpdir, "input.json")
            example_id = request.name or "prediction"

            # RF3 JSON format for sequence prediction
            # Uses 'components' array with sequence and chain_id
            json_content = {
                "name": example_id,
                "components": [
                    {
                        "sequence": request.sequence,
                        "chain_id": "A"
                    }
                ]
            }
            with open(config_path, "w") as f:
                json.dump(json_content, f, indent=2)

            print(f"[RF3] Created JSON config: {json.dumps(json_content, indent=2)}")

            out_dir = os.path.join(tmpdir, "output")
            os.makedirs(out_dir, exist_ok=True)

            # Create .project-root file in /workspace for rootutils
            workspace_root = "/workspace/.project-root"
            if not os.path.exists(workspace_root):
                try:
                    with open(workspace_root, "w") as f:
                        f.write("")
                    print(f"[RF3] Created {workspace_root}")
                except Exception as e:
                    print(f"[RF3] Could not create {workspace_root}: {e}")

            # Set up environment with checkpoint directory
            env = os.environ.copy()
            env["FOUNDRY_CHECKPOINT_DIRS"] = "/workspace/checkpoints"

            rf3_cli = get_cli_path("rf3")
            cmd = f"{rf3_cli} predict out_dir={out_dir} inputs={config_path}"
            print(f"[RF3] Running: {cmd}")
            print(f"[RF3] FOUNDRY_CHECKPOINT_DIRS={env.get('FOUNDRY_CHECKPOINT_DIRS')}")
            result = subprocess.run(
                cmd, shell=True, capture_output=True, text=True, timeout=3600,
                cwd="/workspace", env=env
            )

            if result.returncode == 0:
                outputs = []
                # Check output directory for results
                for filename in os.listdir(out_dir):
                    if filename.endswith((".pdb", ".cif")):
                        with open(os.path.join(out_dir, filename)) as f:
                            outputs.append({"filename": filename, "content": f.read()})

                # Also check subdirectories (RF3 may create example_id subdirectory)
                for subdir in os.listdir(out_dir):
                    subpath = os.path.join(out_dir, subdir)
                    if os.path.isdir(subpath):
                        for filename in os.listdir(subpath):
                            if filename.endswith((".pdb", ".cif")):
                                with open(os.path.join(subpath, filename)) as f:
                                    outputs.append({"filename": filename, "content": f.read()})

                jobs[job_id]["status"] = "completed"
                jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
                jobs[job_id]["result"] = {
                    "predictions": outputs,
                    "confidences": None,  # Not available from CLI
                    "mode": "real (cli)",
                    "stdout": result.stdout[-2000:] if result.stdout else None,  # Last 2000 chars
                }
            else:
                jobs[job_id]["status"] = "failed"
                error_msg = result.stderr or result.stdout or "Unknown error"
                jobs[job_id]["error"] = error_msg

    except Exception as e:
        import traceback
        jobs[job_id]["status"] = "failed"
        jobs[job_id]["error"] = f"{str(e)}\n{traceback.format_exc()}"


# ============== MPNN Job Processing (Python API with model type) ==============

async def process_mpnn_job(job_id: str, request: MPNNRequest):
    jobs[job_id]["status"] = "running"

    if not FOUNDRY_AVAILABLE:
        await process_mpnn_mock(job_id, request)
    else:
        await process_mpnn_real(job_id, request)


async def process_mpnn_mock(job_id: str, request: MPNNRequest):
    await asyncio.sleep(2)

    # Parse PDB to get approximate length
    lines = request.pdb_content.split('\n')
    ca_count = sum(1 for l in lines if l.startswith('ATOM') and ' CA ' in l)
    length = ca_count if ca_count > 0 else 100

    jobs[job_id]["status"] = "completed"
    jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
    jobs[job_id]["result"] = {
        "sequences": [{
            "filename": "sequences.fasta",
            "content": generate_mock_fasta(request.num_sequences, length),
        }],
        "model_type": request.model_type,
        "mode": "mock",
    }


async def process_mpnn_real(job_id: str, request: MPNNRequest):
    """Real MPNN implementation using Python API"""
    try:
        from mpnn.inference_engines.mpnn import MPNNInferenceEngine
        from biotite.structure import get_residue_starts
        from biotite.sequence import ProteinSequence

        # Parse input PDB to atom array
        atom_array = pdb_to_atom_array(request.pdb_content)

        # Configure MPNN engine (per IPD official example)
        engine_config = {
            "model_type": request.model_type,  # "ligand_mpnn" or "protein_mpnn"
            "is_legacy_weights": True,  # Required for current models
            "out_directory": None,  # Return in memory
            "write_structures": False,
            "write_fasta": False,
        }

        # Per-input config
        input_configs = [{
            "batch_size": request.num_sequences,
            "remove_waters": request.remove_waters,
        }]

        # Run MPNN
        model = MPNNInferenceEngine(**engine_config)
        mpnn_outputs = model.run(input_dicts=input_configs, atom_arrays=[atom_array])

        # Extract sequences from outputs
        sequences = []
        for i, item in enumerate(mpnn_outputs):
            if hasattr(item, 'atom_array'):
                # Extract sequence from atom array
                res_starts = get_residue_starts(item.atom_array)
                seq = ''.join(
                    ProteinSequence.convert_letter_3to1(res_name)
                    for res_name in item.atom_array.res_name[res_starts]
                    if res_name in ProteinSequence._dict_3to1
                )
                sequences.append(f">design_{i+1}\n{seq}")

        fasta_content = "\n".join(sequences)

        jobs[job_id]["status"] = "completed"
        jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
        jobs[job_id]["result"] = {
            "sequences": [{
                "filename": "sequences.fasta",
                "content": fasta_content,
            }],
            "num_sequences": len(sequences),
            "model_type": request.model_type,
            "mode": "real",
        }

    except ImportError as e:
        print(f"[MPNN] Python API not available ({e}), falling back to CLI")
        await process_mpnn_cli(job_id, request)
    except Exception as e:
        # Try CLI fallback on any Python API error
        import traceback
        print(f"[MPNN] Python API failed ({e}), falling back to CLI\n{traceback.format_exc()}")
        try:
            await process_mpnn_cli(job_id, request)
        except Exception as cli_error:
            jobs[job_id]["status"] = "failed"
            jobs[job_id]["error"] = f"Python API: {str(e)}\nCLI fallback: {str(cli_error)}"


async def process_mpnn_cli(job_id: str, request: MPNNRequest):
    """Fallback CLI implementation for MPNN"""
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            pdb_path = os.path.join(tmpdir, "input.pdb")
            with open(pdb_path, "w") as f:
                f.write(request.pdb_content)

            out_dir = os.path.join(tmpdir, "output")
            os.makedirs(out_dir, exist_ok=True)

            mpnn_cli = get_cli_path("mpnn")

            # Build command with correct MPNN CLI arguments
            # CLI uses: --structure_path, --out_directory, --batch_size, --temperature
            # --model_type (ligand_mpnn or protein_mpnn), --is_legacy_weights True
            cmd_parts = [
                mpnn_cli,
                f"--structure_path {pdb_path}",
                f"--out_directory {out_dir}",
                "--name design",
                f"--model_type {request.model_type}",
                f"--batch_size {request.num_sequences}",
                f"--temperature {request.temperature}",
                "--is_legacy_weights True",
                "--write_fasta True",
                "--write_structures False",
            ]

            # Add remove_waters if specified (use True/False/None format)
            if request.remove_waters:
                cmd_parts.append("--remove_waters True")
            else:
                cmd_parts.append("--remove_waters None")

            cmd = " ".join(cmd_parts)
            print(f"[MPNN CLI] Running: {cmd}")

            result = subprocess.run(
                cmd, shell=True, capture_output=True, text=True, timeout=1800
            )

            if result.returncode == 0:
                sequences = []
                # MPNN outputs FASTA files
                for filename in sorted(os.listdir(out_dir)):
                    if filename.endswith((".fa", ".fasta")):
                        with open(os.path.join(out_dir, filename)) as f:
                            sequences.append({"filename": filename, "content": f.read()})

                # Also check for sequences in subdirectories
                if not sequences:
                    for root, dirs, files in os.walk(out_dir):
                        for filename in files:
                            if filename.endswith((".fa", ".fasta")):
                                with open(os.path.join(root, filename)) as f:
                                    sequences.append({"filename": filename, "content": f.read()})

                if sequences:
                    jobs[job_id]["status"] = "completed"
                    jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
                    jobs[job_id]["result"] = {
                        "sequences": sequences,
                        "model_type": request.model_type,
                        "mode": "real (cli)",
                    }
                else:
                    jobs[job_id]["status"] = "failed"
                    jobs[job_id]["error"] = f"No output sequences found. stdout: {result.stdout}"
            else:
                jobs[job_id]["status"] = "failed"
                jobs[job_id]["error"] = result.stderr or result.stdout or "Unknown error"

    except Exception as e:
        import traceback
        jobs[job_id]["status"] = "failed"
        jobs[job_id]["error"] = f"{str(e)}\n{traceback.format_exc()}"


# ============== RMSD Validation ==============

def calculate_rmsd(pdb1: str, pdb2: str, backbone_only: bool = True) -> Dict[str, Any]:
    """Calculate RMSD between two structures"""
    try:
        from biotite.structure import rmsd, superimpose
        import numpy as np

        # Backbone atom names
        BACKBONE_ATOMS = ['N', 'CA', 'C', 'O']

        # Parse structures
        aa1 = pdb_to_atom_array(pdb1)
        aa2 = pdb_to_atom_array(pdb2)

        if backbone_only:
            # Filter to backbone atoms
            mask1 = np.isin(aa1.atom_name, BACKBONE_ATOMS)
            mask2 = np.isin(aa2.atom_name, BACKBONE_ATOMS)
            aa1 = aa1[mask1]
            aa2 = aa2[mask2]

        # Check if structures have same number of atoms
        if len(aa1) != len(aa2):
            # Try to align by residue
            min_len = min(len(aa1), len(aa2))
            aa1 = aa1[:min_len]
            aa2 = aa2[:min_len]

        # Superimpose and calculate RMSD
        aa2_fitted, _ = superimpose(aa1, aa2)
        rmsd_value = float(rmsd(aa1, aa2_fitted))

        # Interpret RMSD
        if rmsd_value < 1.0:
            interpretation = "Excellent"
            description = "Very high designability - structure is highly likely to fold as designed"
        elif rmsd_value < 2.0:
            interpretation = "Good"
            description = "Good designability - structure will likely fold correctly"
        elif rmsd_value < 3.0:
            interpretation = "Moderate"
            description = "Moderate designability - some structural deviation expected"
        else:
            interpretation = "Poor"
            description = "Low designability - significant structural deviation"

        return {
            "rmsd": round(rmsd_value, 3),
            "interpretation": interpretation,
            "description": description,
            "backbone_only": backbone_only,
            "num_atoms_compared": len(aa1),
        }

    except Exception as e:
        return {
            "error": str(e),
            "rmsd": None,
        }


# ============== Endpoints ==============

@app.get("/", response_model=HealthResponse)
@app.get("/health", response_model=HealthResponse)
async def health():
    return HealthResponse(
        status="healthy",
        mode="real" if FOUNDRY_AVAILABLE else "mock",
        gpu_available=GPU_INFO.get("available", False),
        gpu_name=GPU_INFO.get("name"),
        gpu_memory_gb=GPU_INFO.get("memory_gb"),
        models=get_model_status(),
    )


@app.post("/api/rfd3/design", response_model=JobResponse)
async def submit_rfd3_design(request: RFD3Request, background_tasks: BackgroundTasks):
    job_id = str(uuid.uuid4())

    # Build request summary for logging
    request_summary = {
        "contig": request.contig,
        "num_designs": request.num_designs,
        "seed": request.seed,
    }
    # Add Phase 1 parameters if set
    if request.select_hotspots:
        request_summary["select_hotspots"] = request.select_hotspots
    if request.partial_t is not None:
        request_summary["partial_t"] = request.partial_t
    if request.num_timesteps is not None:
        request_summary["num_timesteps"] = request.num_timesteps
    if request.step_scale is not None:
        request_summary["step_scale"] = request.step_scale
    if request.gamma_0 is not None:
        request_summary["gamma_0"] = request.gamma_0
    # Add Phase 2 parameters if set
    if request.symmetry:
        request_summary["symmetry"] = request.symmetry.id
    if request.ligand:
        request_summary["ligand"] = request.ligand
    if request.cfg and request.cfg.enabled:
        request_summary["cfg"] = {"scale": request.cfg.scale, "features": request.cfg.features}

    jobs[job_id] = {
        "job_id": job_id,
        "status": "pending",
        "created_at": datetime.utcnow().isoformat(),
        "type": "rfd3",
        "request": request_summary,
    }

    background_tasks.add_task(process_rfd3_job, job_id, request)

    return JobResponse(
        job_id=job_id,
        status="pending",
        message="RFD3 design job submitted successfully",
    )


@app.post("/api/rf3/predict", response_model=JobResponse)
async def submit_rf3_prediction(request: RF3Request, background_tasks: BackgroundTasks):
    job_id = str(uuid.uuid4())

    jobs[job_id] = {
        "job_id": job_id,
        "status": "pending",
        "created_at": datetime.utcnow().isoformat(),
        "type": "rf3",
        "request": {"sequence_length": len(request.sequence)},
    }

    background_tasks.add_task(process_rf3_job, job_id, request)

    return JobResponse(
        job_id=job_id,
        status="pending",
        message="RF3 prediction job submitted successfully",
    )


@app.post("/api/mpnn/design", response_model=JobResponse)
async def submit_mpnn_design(request: MPNNRequest, background_tasks: BackgroundTasks):
    job_id = str(uuid.uuid4())

    jobs[job_id] = {
        "job_id": job_id,
        "status": "pending",
        "created_at": datetime.utcnow().isoformat(),
        "type": "mpnn",
        "request": {
            "num_sequences": request.num_sequences,
            "temperature": request.temperature,
            "model_type": request.model_type,
        },
    }

    background_tasks.add_task(process_mpnn_job, job_id, request)

    return JobResponse(
        job_id=job_id,
        status="pending",
        message="ProteinMPNN design job submitted successfully",
    )


@app.post("/api/validate/rmsd")
async def validate_rmsd(request: RMSDRequest):
    """Calculate RMSD between two structures for designability validation"""
    result = calculate_rmsd(
        request.pdb_content_1,
        request.pdb_content_2,
        request.backbone_only
    )

    if "error" in result:
        raise HTTPException(status_code=400, detail=result["error"])

    return result


@app.post("/api/export/structure")
async def export_structure(request: ExportRequest):
    """Export structure in different formats"""
    try:
        atom_array = pdb_to_atom_array(request.pdb_content)

        if request.format == "cif":
            content = atom_array_to_cif(atom_array)
            filename = "structure.cif"
        else:
            content = request.pdb_content
            filename = "structure.pdb"

        result = {
            "filename": filename,
            "content": content,
            "format": request.format,
        }

        if request.include_confidences and request.confidences:
            result["confidences_json"] = json.dumps(request.confidences, indent=2)

        return result

    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.get("/api/stored/{job_id}")
async def get_stored_structure(job_id: str):
    """Get stored structure from a completed job for cross-panel data flow"""
    if job_id not in stored_structures:
        raise HTTPException(status_code=404, detail="Structure not found")

    return stored_structures[job_id]


@app.get("/api/jobs/{job_id}", response_model=JobStatus)
async def get_job_status(job_id: str):
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")

    job = jobs[job_id]
    return JobStatus(
        job_id=job["job_id"],
        status=job["status"],
        created_at=job["created_at"],
        completed_at=job.get("completed_at"),
        result=job.get("result"),
        error=job.get("error"),
    )


@app.get("/api/jobs")
async def list_jobs():
    return {
        job_id: {
            "status": job["status"],
            "type": job.get("type"),
            "created_at": job["created_at"],
        }
        for job_id, job in jobs.items()
    }


@app.delete("/api/jobs/{job_id}")
async def delete_job(job_id: str):
    if job_id in jobs:
        del jobs[job_id]
    if job_id in stored_structures:
        del stored_structures[job_id]
    return {"status": "deleted"}


# ============== AI Protein Engineering Endpoints ==============

@app.post("/api/fetch/pdb")
async def fetch_pdb_from_rcsb(request: FetchPDBRequest):
    """
    Fetch PDB file from RCSB Protein Data Bank.

    Supports both PDB and mmCIF formats with local caching.
    """
    try:
        from utils.pdb_fetch import fetch_pdb, fetch_cif, parse_pdb_content

        if request.format == "cif":
            result = fetch_cif(request.pdb_id, use_cache=request.use_cache)
        else:
            result = fetch_pdb(request.pdb_id, use_cache=request.use_cache)

        if not result["success"]:
            raise HTTPException(status_code=404, detail=result["error"])

        # Parse basic info from content
        if request.format == "pdb":
            info = parse_pdb_content(result["content"])
        else:
            info = {"format": "cif"}

        return {
            "pdb_id": result["pdb_id"],
            "content": result["content"],
            "source": result["source"],
            "format": request.format,
            "info": info,
        }

    except ImportError:
        raise HTTPException(
            status_code=501,
            detail="PDB fetch module not available. Install required dependencies."
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/analyze/metal-binding")
async def analyze_metal_binding(request: MetalBindingAnalysisRequest):
    """
    Analyze metal binding site in a protein structure.

    Returns comprehensive analysis including:
    - Coordination number and geometry
    - Coordinating residues and atoms
    - Donor type classification
    - Bond distances and angles
    - Suggestions for optimization
    """
    try:
        from utils.coordination import analyze_coordination_geometry, suggest_lanthanide_conversion
        from utils.pdb_fetch import fetch_pdb

        # Get PDB content
        pdb_content = request.pdb_content
        if pdb_content is None and request.pdb_id:
            result = fetch_pdb(request.pdb_id)
            if not result["success"]:
                raise HTTPException(status_code=404, detail=result["error"])
            pdb_content = result["content"]

        if pdb_content is None:
            raise HTTPException(
                status_code=400,
                detail="Either pdb_content or pdb_id must be provided"
            )

        # Analyze coordination geometry
        analysis = analyze_coordination_geometry(
            pdb_content,
            request.metal_chain,
            request.metal_residue,
            request.metal_resnum,
            request.distance_cutoff
        )

        if not analysis.get("success"):
            raise HTTPException(status_code=400, detail=analysis.get("error"))

        return analysis

    except ImportError as e:
        raise HTTPException(
            status_code=501,
            detail=f"Analysis module not available: {str(e)}"
        )
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/ai/recommend-parameters")
async def recommend_rfd3_parameters(request: AIRecommendRequest):
    """
    AI-assisted RFD3 parameter recommendation.

    Uses enhanced hybrid approach:
    1. Chemistry rules from research/literature
    2. Optional LLM refinement for edge cases

    Returns:
    - Recommended RFD3 parameters
    - Design strategy and reasoning
    - Evaluation criteria for output assessment
    """
    try:
        from utils.coordination import analyze_coordination_geometry
        from utils.pdb_fetch import fetch_pdb
        from ai.recommender import generate_rfd3_parameters, AIRecommender

        # Get PDB content
        pdb_content = request.pdb_content
        if pdb_content is None and request.pdb_id:
            result = fetch_pdb(request.pdb_id)
            if not result["success"]:
                raise HTTPException(status_code=404, detail=result["error"])
            pdb_content = result["content"]

        if pdb_content is None:
            raise HTTPException(
                status_code=400,
                detail="Either pdb_content or pdb_id must be provided"
            )

        # Analyze current metal binding site
        analysis = analyze_coordination_geometry(
            pdb_content,
            request.metal_chain,
            request.metal_residue,
            request.metal_resnum,
        )

        if not analysis.get("success"):
            raise HTTPException(status_code=400, detail=analysis.get("error"))

        # Generate RFD3 parameters using AI recommender
        recommendation = generate_rfd3_parameters(
            current_analysis=analysis,
            target_metal=request.target_metal,
            user_description=request.user_description,
        )

        if not recommendation.get("success"):
            raise HTTPException(status_code=400, detail=recommendation.get("error"))

        # Add optional SASA analysis
        if request.include_sasa:
            try:
                from utils.sasa import calculate_sasa, get_binding_pocket_sasa

                coord_atoms = analysis["coordination"]["coordinating_atoms"]
                pocket_residues = [
                    f"{a['chain']}:{a['residue']}{a['residue_number']}"
                    for a in coord_atoms
                ]

                sasa_result = get_binding_pocket_sasa(pdb_content, pocket_residues)
                recommendation["sasa_analysis"] = sasa_result.get("pocket_analysis")
            except ImportError:
                recommendation["sasa_analysis"] = None
                recommendation["sasa_warning"] = "SASA module not available"

        # Add optional secondary structure context
        if request.include_secondary_structure:
            try:
                from utils.dssp import assign_secondary_structure, get_ss_context

                ss_result = assign_secondary_structure(pdb_content)
                coord_atoms = analysis["coordination"]["coordinating_atoms"]
                residue_ids = [
                    f"{a['chain']}:{a['residue']}{a['residue_number']}"
                    for a in coord_atoms
                ]

                ss_context = get_ss_context(ss_result, residue_ids)
                recommendation["secondary_structure"] = {
                    "content": ss_result.get("content"),
                    "binding_site_context": ss_context.get("context"),
                }
            except ImportError:
                recommendation["secondary_structure"] = None
                recommendation["ss_warning"] = "DSSP module not available"

        return {
            "success": True,
            "current_analysis": analysis,
            "recommendation": recommendation,
            "pdb_id": request.pdb_id,
        }

    except ImportError as e:
        raise HTTPException(
            status_code=501,
            detail=f"Required module not available: {str(e)}"
        )
    except HTTPException:
        raise
    except Exception as e:
        import traceback
        raise HTTPException(status_code=500, detail=f"{str(e)}\n{traceback.format_exc()}")


@app.post("/api/evaluate/design")
async def evaluate_design(
    pdb_content: str,
    target_metal: str,
    target_coordination: Optional[int] = None,
    metal_chain: str = "A",
    metal_resnum: int = 1,
):
    """
    Evaluate a designed structure against target specifications.

    Checks:
    - Coordination geometry
    - Donor types
    - Bond distances
    - Geometry quality
    """
    try:
        from utils.coordination import analyze_coordination_geometry
        from ai.recommender import MetalChemistryRules

        # Analyze the design
        analysis = analyze_coordination_geometry(
            pdb_content,
            metal_chain,
            target_metal,
            metal_resnum,
        )

        if not analysis.get("success"):
            return {
                "success": False,
                "error": "Could not find metal in structure",
                "suggestion": "The design may not have placed the metal. Try with higher partial_t or check output manually.",
            }

        # Get target properties
        target_props = MetalChemistryRules.get_metal_properties(target_metal)

        if target_props is None:
            return {
                "success": True,
                "analysis": analysis,
                "evaluation": "Unknown metal - cannot evaluate against target",
            }

        # Evaluate against targets
        evaluation = {
            "coordination": {
                "actual": analysis["coordination"]["number"],
                "target_range": target_props.typical_coordination,
                "status": "PASS" if analysis["coordination"]["number"] in target_props.typical_coordination else "FAIL",
            },
            "geometry": {
                "actual": analysis["coordination"]["geometry"],
                "target": target_props.preferred_geometries,
                "rmsd": analysis["coordination"]["geometry_rmsd"],
                "status": "PASS" if analysis["coordination"]["geometry"] in target_props.preferred_geometries else "WARN",
            },
            "bond_distances": {
                "average": analysis["bond_analysis"]["average_distance"],
                "target_range": target_props.bond_distance_range,
                "status": "PASS" if (
                    target_props.bond_distance_range[0] <= analysis["bond_analysis"]["average_distance"] <= target_props.bond_distance_range[1]
                ) else "WARN",
            },
            "donor_types": {
                "actual": analysis["donor_analysis"]["types"],
                "preferred": target_props.preferred_donors,
                "avoid": target_props.avoid_donors,
            },
        }

        # Calculate overall score
        pass_count = sum(1 for v in evaluation.values() if isinstance(v, dict) and v.get("status") == "PASS")
        total_checks = sum(1 for v in evaluation.values() if isinstance(v, dict) and "status" in v)
        score = (pass_count / total_checks * 100) if total_checks > 0 else 0

        # Generate recommendations
        recommendations = []
        if evaluation["coordination"]["status"] != "PASS":
            coord_diff = min(target_props.typical_coordination) - analysis["coordination"]["number"]
            if coord_diff > 0:
                recommendations.append(f"Increase coordination by {coord_diff}. Try higher partial_t or add more Asp/Glu.")
            else:
                recommendations.append(f"Decrease coordination. Current site may be too crowded.")

        if evaluation["bond_distances"]["status"] != "PASS":
            recommendations.append("Bond distances outside optimal range. Consider refining with lower partial_t.")

        return {
            "success": True,
            "analysis": analysis,
            "evaluation": evaluation,
            "score": round(score, 1),
            "recommendations": recommendations,
            "verdict": "ACCEPT" if score >= 70 else "NEEDS_IMPROVEMENT",
        }

    except ImportError as e:
        raise HTTPException(
            status_code=501,
            detail=f"Evaluation module not available: {str(e)}"
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# ============== AI-Driven Analysis Endpoint ==============

class AIAnalyzeRequest(BaseModel):
    """Request for AI-driven analysis of user intent."""
    user_input: str
    pdb_content: Optional[str] = None
    pdb_id: Optional[str] = None
    structure_info: Optional[Dict[str, Any]] = None
    conversation_history: Optional[List[Dict[str, str]]] = None


@app.post("/api/ai/analyze")
async def ai_analyze(request: AIAnalyzeRequest):
    """
    AI-driven analysis of user request for any RFdiffusion task.

    This endpoint uses a hybrid approach:
    1. Handbook knowledge base for parameter guidance
    2. Optional LLM (Claude) for natural language interpretation

    Supports ALL RFdiffusion tasks:
    - de_novo: Design new proteins from scratch
    - binder: Design protein binders
    - metal_redesign: Modify metal binding sites
    - enzyme: Design enzyme active sites
    - refinement: Improve existing structures
    - symmetric: Design symmetric oligomers
    - scaffold: Graft functional motifs

    Returns:
    - task_type: Detected design task
    - params: Suggested RFD3 parameters
    - reasoning: Plain English explanation
    - form_config: Dynamic form fields for UI
    - clarifying_questions: Questions if intent unclear
    """
    try:
        from ai.ai_engine import analyze_user_request
        import os

        # Get PDB content if ID provided
        pdb_content = request.pdb_content
        if pdb_content is None and request.pdb_id:
            try:
                from utils.pdb_fetch import fetch_pdb
                result = fetch_pdb(request.pdb_id)
                if result.get("success"):
                    pdb_content = result["content"]
            except ImportError:
                pass

        # Get API key from environment
        api_key = os.environ.get("ANTHROPIC_API_KEY")

        # Run AI analysis
        result = analyze_user_request(
            user_input=request.user_input,
            structure_analysis=request.structure_info,
            pdb_content=pdb_content,
            conversation_history=request.conversation_history,
            anthropic_api_key=api_key,
        )

        return result

    except ImportError as e:
        raise HTTPException(
            status_code=501,
            detail=f"AI engine module not available: {str(e)}"
        )
    except Exception as e:
        import traceback
        raise HTTPException(
            status_code=500,
            detail=f"AI analysis failed: {str(e)}\n{traceback.format_exc()}"
        )


@app.post("/api/ai/refine")
async def ai_refine(
    current_params: Dict[str, Any],
    user_feedback: str,
    task_type: str = "unknown",
):
    """
    Refine AI-suggested parameters based on user feedback.

    Args:
        current_params: Current parameter values
        user_feedback: User's adjustment request (e.g., "more aggressive")
        task_type: Current task type

    Returns:
        Updated parameters with reasoning
    """
    try:
        from ai.ai_engine import AIEngine, TaskType
        import os

        api_key = os.environ.get("ANTHROPIC_API_KEY")
        engine = AIEngine(api_key)

        # Parse task type
        try:
            task = TaskType(task_type)
        except ValueError:
            task = TaskType.UNKNOWN

        result = engine.refine_parameters(
            current_params=current_params,
            user_feedback=user_feedback,
            task_type=task,
        )

        return {
            "success": result.success,
            "params": result.params,
            "reasoning": result.reasoning,
            "confidence": result.confidence,
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# ============== AI Task Planning ==============

class TaskPlanRequest(BaseModel):
    """Request for AI task planning."""
    goal: str
    template: Optional[str] = None
    params: Optional[Dict[str, Any]] = None
    custom_steps: Optional[List[Dict[str, Any]]] = None


@app.post("/api/ai/plan")
async def ai_plan_task(request: TaskPlanRequest):
    """
    Plan a complex protein engineering task.

    Uses AI to break down natural language goals into
    executable step-by-step plans.

    Templates available:
    - metal_replacement: Convert metal binding sites
    - symmetric_binder: Design symmetric protein with ligand
    - de_novo: Design protein from scratch

    Returns:
        TaskPlan with steps, parameters, and AI context
    """
    try:
        from ai.task_planner import plan_task

        result = plan_task(
            goal=request.goal,
            template=request.template,
            params=request.params,
            custom_steps=request.custom_steps,
        )

        return result

    except Exception as e:
        import traceback
        raise HTTPException(
            status_code=500,
            detail=f"Task planning failed: {str(e)}\n{traceback.format_exc()}"
        )


@app.post("/api/ai/execute-step")
async def ai_execute_step(
    step_type: str,
    params: Dict[str, Any],
    artifacts: Optional[Dict[str, str]] = None,
):
    """
    Execute a single step from a task plan.

    This allows the AI to control step-by-step execution
    and review results between steps.

    Args:
        step_type: Type of step (rfd3_design, mpnn_design, etc.)
        params: Parameters for the step
        artifacts: Previous artifacts (pdb_content, sequences, etc.)

    Returns:
        Step result with metrics for AI review
    """
    try:
        from ai.task_planner import TaskExecutor, TaskStep, TaskStepType

        executor = TaskExecutor()

        # Restore artifacts from previous steps
        if artifacts:
            executor.artifacts = artifacts

        step = TaskStep(
            id="single-step",
            type=TaskStepType(step_type),
            name=f"Execute {step_type}",
            description=f"Single step execution: {step_type}",
            params=params,
        )

        # Create a minimal plan for single step execution
        from ai.task_planner import TaskPlan
        plan = TaskPlan(
            id="single",
            name="Single step",
            description="Single step execution",
            goal="Execute step",
            steps=[step],
        )

        result = await executor._execute_step(step, plan)

        return {
            "status": "completed",
            "result": result,
            "artifacts": executor.artifacts,
        }

    except Exception as e:
        import traceback
        return {
            "status": "failed",
            "error": str(e),
            "traceback": traceback.format_exc(),
        }


# ============== Main ==============

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
