"""
Inference Utilities for RunPod Serverless Handler

Shared functions for running RFD3, RF3, MPNN inference and utilities.
Extracted from main.py for use in serverless environment.
"""

import os
import io
import json
import random
import subprocess
import tempfile
from typing import Dict, Any, Optional, List


# ============== Utility Functions ==============

def get_gpu_info() -> Dict[str, Any]:
    """Get GPU information via nvidia-smi"""
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


def check_foundry_available() -> bool:
    """Check if Foundry Python modules are importable"""
    try:
        from rfd3.engine import RFD3InferenceEngine
        return True
    except ImportError:
        pass

    try:
        from rf3.inference_engines.rf3 import RF3InferenceEngine
        return True
    except ImportError:
        pass

    return False


def atom_array_to_pdb(atom_array) -> str:
    """Convert biotite AtomArray to PDB string"""
    from biotite.structure.io.pdb import PDBFile

    pdb_file = PDBFile()
    pdb_file.set_structure(atom_array)
    buf = io.StringIO()
    pdb_file.write(buf)
    return buf.getvalue()


def atom_array_to_cif(atom_array) -> str:
    """Convert biotite AtomArray to CIF string"""
    try:
        from atomworks.io.utils.io_utils import to_cif_file

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

        pdbx_file = PDBxFile()
        set_structure(pdbx_file, atom_array)
        buf = io.StringIO()
        pdbx_file.write(buf)
        return buf.getvalue()


def pdb_to_atom_array(content: str):
    """Convert PDB or CIF string to biotite AtomArray"""
    content_stripped = content.strip()

    if content_stripped.startswith('data_'):
        # CIF format
        from biotite.structure.io.pdbx import CIFFile, get_structure

        cif_file = CIFFile.read(io.StringIO(content))
        block = list(cif_file.values())[0] if cif_file else None
        if block is None:
            raise ValueError("Empty CIF file")
        return get_structure(block, model=1)
    else:
        # PDB format
        from biotite.structure.io.pdb import PDBFile

        pdb_file = PDBFile.read(io.StringIO(content))
        return pdb_file.get_structure(model=1)


# ============== Mock Implementations ==============

def generate_mock_pdb(length: int = 100) -> str:
    """Generate mock PDB content"""
    lines = ["HEADER    MOCK PROTEIN STRUCTURE"]
    for i in range(1, length + 1):
        x, y, z = i * 0.5, i * 0.3, i * 0.2
        lines.append(
            f"ATOM  {i:5d}  CA  ALA A{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
        )
    lines.append("END")
    return "\n".join(lines)


def generate_mock_fasta(num_seq: int = 8, length: int = 100) -> str:
    """Generate mock FASTA sequences"""
    aa = "ACDEFGHIKLMNPQRSTVWY"
    lines = []
    for i in range(num_seq):
        seq = "".join(random.choices(aa, k=length))
        lines.append(f">design_{i+1}")
        lines.append(seq)
    return "\n".join(lines)


def generate_mock_confidences(length: int = 100) -> Dict[str, Any]:
    """Generate mock confidence metrics"""
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


# ============== RFD3 Inference ==============

def run_rfd3_inference(
    contig: Optional[str] = None,
    length: Optional[str] = None,
    num_designs: int = 1,
    seed: Optional[int] = None,
    pdb_content: Optional[str] = None,
    # Quality settings
    num_timesteps: Optional[int] = None,
    step_scale: Optional[float] = None,
    gamma_0: Optional[float] = None,
    is_non_loopy: Optional[bool] = None,
    # Partial diffusion (refinement)
    partial_t: Optional[int] = None,
    # Symmetry
    symmetry: Optional[Dict[str, Any]] = None,
    # Small molecule / enzyme design
    ligand: Optional[str] = None,
    select_fixed_atoms: Optional[Dict[str, str]] = None,
    unindex: Optional[str] = None,
    # RASA conditioning (binding pocket design)
    select_buried: Optional[Dict[str, str]] = None,
    select_exposed: Optional[Dict[str, str]] = None,
    select_partially_buried: Optional[Dict[str, str]] = None,
    # Protein binder design
    hotspots: Optional[List[str]] = None,
    infer_ori_strategy: Optional[str] = None,
    # Nucleic acid binder design
    na_chains: Optional[str] = None,
    ori_token: Optional[List[float]] = None,
    select_hbond_donor: Optional[Dict[str, str]] = None,
    select_hbond_acceptor: Optional[Dict[str, str]] = None,
    # Mock mode
    use_mock: bool = False
) -> Dict[str, Any]:
    """
    Run RFD3 inference with full parameter support.

    Supports multiple design tasks:
    - De novo protein design (length only)
    - Protein binder design (contig + hotspots)
    - Small molecule binder design (ligand + RASA conditioning)
    - Nucleic acid binder design (na_chains + ori_token + H-bond conditioning)
    - Enzyme scaffold design (ligand + unindex + fixed atoms)
    - Symmetric oligomer design (symmetry config)
    - Structure refinement (partial_t)
    """

    if use_mock:
        mock_length = length or contig or "100"
        return run_rfd3_mock(mock_length, num_designs)

    try:
        from rfd3.engine import RFD3InferenceConfig, RFD3InferenceEngine

        # Set seed for reproducibility
        if seed is not None:
            try:
                from lightning.fabric import seed_everything
                seed_everything(seed)
                print(f"[RFD3] Set seed to {seed}")
            except ImportError:
                pass

        # Build specification dictionary
        spec: Dict[str, Any] = {}

        # Length or contig
        if contig:
            contig_str = contig
            is_simple_length = contig_str.replace("-", "").isdigit()
            if is_simple_length:
                if "-" in contig_str:
                    spec["length"] = contig_str
                else:
                    spec["length"] = int(contig_str)
            else:
                spec["contig"] = contig_str
        elif length:
            length_str = str(length)
            if "-" in length_str:
                spec["length"] = length_str
            elif length_str.isdigit():
                spec["length"] = int(length_str)
            else:
                spec["length"] = length_str

        # Ligand for small molecule / enzyme design
        if ligand:
            spec["ligand"] = ligand

        # Unindex for enzyme design (residues with inferred positions)
        if unindex:
            spec["unindex"] = unindex

        # Symmetry configuration (symmetry ID goes in spec)
        if symmetry:
            sym_id = symmetry.get("id") if isinstance(symmetry, dict) else symmetry
            if sym_id:
                spec["symmetry"] = sym_id
                # Symmetric motif flags
                if isinstance(symmetry, dict):
                    if symmetry.get("is_symmetric_motif"):
                        spec["is_symmetric_motif"] = True
                    if symmetry.get("is_unsym_motif"):
                        spec["is_unsym_motif"] = symmetry["is_unsym_motif"]

        # is_non_loopy goes in SPEC (not top-level)
        if is_non_loopy is not None:
            spec["is_non_loopy"] = is_non_loopy

        # partial_t goes in SPEC (not top-level)
        if partial_t is not None:
            spec["partial_t"] = partial_t

        # Hotspots -> select_hotspots in SPEC
        if hotspots:
            # Convert list to dict format: ["A15", "A20"] -> {"A15": "ALL", "A20": "ALL"}
            spec["select_hotspots"] = {h: "ALL" for h in hotspots}

        # Origin and strategy go in SPEC
        if ori_token:
            spec["ori_token"] = ori_token
        if infer_ori_strategy:
            spec["infer_ori_strategy"] = infer_ori_strategy

        # RASA conditioning goes in SPEC
        if select_buried:
            spec["select_buried"] = select_buried
        if select_exposed:
            spec["select_exposed"] = select_exposed
        if select_partially_buried:
            spec["select_partially_buried"] = select_partially_buried

        # Fixed atoms go in SPEC
        if select_fixed_atoms:
            spec["select_fixed_atoms"] = select_fixed_atoms

        # H-bond conditioning goes in SPEC
        if select_hbond_donor:
            spec["select_hbond_donor"] = select_hbond_donor
        if select_hbond_acceptor:
            spec["select_hbond_acceptor"] = select_hbond_acceptor

        # Build inference sampler config for DIFFUSION parameters ONLY
        sampler_config: Dict[str, Any] = {}
        if num_timesteps is not None:
            sampler_config["num_timesteps"] = num_timesteps
        if step_scale is not None:
            sampler_config["step_scale"] = step_scale
        if gamma_0 is not None:
            sampler_config["gamma_0"] = gamma_0

        # Symmetry requires special sampler kind
        if symmetry:
            sampler_config["kind"] = "symmetry"

        # Build config kwargs - ONLY these top-level params are allowed!
        config_kwargs: Dict[str, Any] = {
            "specification": spec,
            "diffusion_batch_size": num_designs,
        }

        # Add sampler config if any diffusion parameters were set
        if sampler_config:
            config_kwargs["inference_sampler"] = sampler_config

        # High-order symmetry requires low memory mode
        if symmetry:
            sym_id = symmetry.get("id") if isinstance(symmetry, dict) else symmetry
            if sym_id in ["T", "O", "I"]:
                config_kwargs["low_memory_mode"] = True

        # Handle PDB input - write to temp file and pass path via spec["input"]
        temp_pdb_path = None
        if pdb_content:
            # Write PDB to temp file - RFD3 expects file path, not AtomArray
            with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode='w') as f:
                f.write(pdb_content)
                temp_pdb_path = f.name
            spec["input"] = temp_pdb_path
            print(f"[RFD3] Wrote input structure to {temp_pdb_path}")

        print(f"[RFD3] Spec: {json.dumps({k: str(v)[:100] for k, v in spec.items()}, indent=2)}")
        print(f"[RFD3] Config: {json.dumps({k: str(v)[:100] for k, v in config_kwargs.items()}, indent=2)}")

        # Create config and run
        config = RFD3InferenceConfig(**config_kwargs)
        model = RFD3InferenceEngine(**config)

        # Run inference - inputs are passed via spec["input"], not model.run()
        outputs_dict = model.run(inputs=None, out_dir=None, n_batches=1)

        # Process outputs
        designs = []
        if outputs_dict:
            for key, result_list in outputs_dict.items():
                items = result_list if isinstance(result_list, list) else [result_list]

                for idx, item in enumerate(items):
                    if hasattr(item, 'atom_array'):
                        try:
                            output_pdb = atom_array_to_pdb(item.atom_array)
                            filename = f"{key}_{idx}.pdb" if len(items) > 1 else f"{key}.pdb"

                            design: Dict[str, Any] = {"filename": filename, "content": output_pdb}

                            try:
                                design["cif_content"] = atom_array_to_cif(item.atom_array)
                            except Exception:
                                pass

                            designs.append(design)
                        except Exception as e:
                            print(f"[RFD3] Error converting output: {e}")

        # Cleanup temp file
        if temp_pdb_path and os.path.exists(temp_pdb_path):
            os.unlink(temp_pdb_path)

        return {
            "status": "completed",
            "result": {
                "designs": designs,
                "mode": "real",
                "seed": seed,
            }
        }

    except Exception as e:
        import traceback
        # Cleanup temp file on error
        if 'temp_pdb_path' in locals() and temp_pdb_path and os.path.exists(temp_pdb_path):
            os.unlink(temp_pdb_path)
        return {
            "status": "failed",
            "error": str(e),
            "traceback": traceback.format_exc()
        }


def run_rfd3_mock(contig: str, num_designs: int) -> Dict[str, Any]:
    """Mock RFD3 inference"""
    try:
        length = int(contig) if contig.isdigit() else 100
    except Exception:
        length = 100

    designs = []
    for i in range(num_designs):
        designs.append({
            "filename": f"design_{i+1}.pdb",
            "content": generate_mock_pdb(length)
        })

    return {
        "status": "completed",
        "result": {
            "designs": designs,
            "mode": "mock",
        }
    }


# ============== RF3 Inference ==============

def run_rf3_inference(
    sequence: str,
    name: str = "prediction",
    pdb_content: Optional[str] = None,
    use_mock: bool = False
) -> Dict[str, Any]:
    """Run RF3 structure prediction"""

    if use_mock:
        return run_rf3_mock(sequence, name)

    try:
        # Try Python API first
        return run_rf3_python_api(sequence, name, pdb_content)
    except ImportError as e:
        print(f"[RF3] Python API not available: {e}")
        # Fallback to CLI
        return run_rf3_cli(sequence, name)
    except Exception as e:
        import traceback
        return {
            "status": "failed",
            "error": str(e),
            "traceback": traceback.format_exc()
        }


def run_rf3_python_api(sequence: str, name: str, pdb_content: Optional[str]) -> Dict[str, Any]:
    """RF3 via Python API"""
    from rf3.inference_engines.rf3 import RF3InferenceEngine
    from rf3.utils.inference import InferenceInput

    # For sequence-only, we need to use CLI approach
    if not pdb_content:
        return run_rf3_cli(sequence, name)

    # Initialize engine
    inference_engine = RF3InferenceEngine(ckpt_path='rf3', verbose=False)

    # Structure-based input
    atom_array = pdb_to_atom_array(pdb_content)
    input_data = InferenceInput.from_atom_array(atom_array, example_id=name)

    # Run inference
    rf3_outputs = inference_engine.run(inputs=input_data)

    # Process outputs
    predictions = []
    confidences = None

    if name in rf3_outputs:
        results = rf3_outputs[name]

        for idx, rf3_output in enumerate(results):
            if hasattr(rf3_output, 'atom_array'):
                pdb_out = atom_array_to_pdb(rf3_output.atom_array)

                pred = {
                    "filename": f"{name}_{idx}.pdb" if len(results) > 1 else f"{name}.pdb",
                    "content": pdb_out,
                }

                try:
                    pred["cif_content"] = atom_array_to_cif(rf3_output.atom_array)
                except Exception:
                    pass

                predictions.append(pred)

            # Extract confidences from first output
            if idx == 0 and confidences is None:
                confidences = extract_rf3_confidences(rf3_output)

    return {
        "status": "completed",
        "result": {
            "predictions": predictions,
            "confidences": confidences,
            "mode": "real",
        }
    }


def run_rf3_cli(sequence: str, name: str) -> Dict[str, Any]:
    """RF3 via CLI for sequence-only input"""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create JSON config
        config_path = os.path.join(tmpdir, "input.json")
        json_content = {
            "name": name,
            "components": [{"sequence": sequence, "chain_id": "A"}]
        }
        with open(config_path, "w") as f:
            json.dump(json_content, f)

        out_dir = os.path.join(tmpdir, "output")
        os.makedirs(out_dir, exist_ok=True)

        # Create .project-root for rootutils
        workspace_root = "/workspace/.project-root"
        if not os.path.exists(workspace_root):
            try:
                with open(workspace_root, "w") as f:
                    f.write("")
            except Exception:
                pass

        # Run CLI
        env = os.environ.copy()
        cmd = f"rf3 predict out_dir={out_dir} inputs={config_path}"

        result = subprocess.run(
            cmd, shell=True, capture_output=True, text=True,
            timeout=3600, cwd="/workspace", env=env
        )

        if result.returncode == 0:
            predictions = []
            # Find output files
            for root, _, files in os.walk(out_dir):
                for filename in files:
                    if filename.endswith((".pdb", ".cif")):
                        with open(os.path.join(root, filename)) as f:
                            predictions.append({"filename": filename, "content": f.read()})

            return {
                "status": "completed",
                "result": {
                    "predictions": predictions,
                    "confidences": None,
                    "mode": "real (cli)",
                }
            }
        else:
            return {
                "status": "failed",
                "error": result.stderr or result.stdout or "Unknown error"
            }


def run_rf3_mock(sequence: str, name: str) -> Dict[str, Any]:
    """Mock RF3 inference"""
    length = len(sequence)
    pdb_content = generate_mock_pdb(length)
    confidences = generate_mock_confidences(length)

    return {
        "status": "completed",
        "result": {
            "predictions": [{
                "filename": f"{name}.pdb",
                "content": pdb_content,
            }],
            "confidences": confidences,
            "mode": "mock",
        }
    }


def extract_rf3_confidences(rf3_output) -> Optional[Dict[str, Any]]:
    """Extract confidence metrics from RF3 output"""
    confidences = {}

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

    if hasattr(rf3_output, 'confidences') and rf3_output.confidences:
        conf = rf3_output.confidences
        if 'atom_plddts' in conf:
            atom_plddts = list(conf['atom_plddts'])
            confidences["per_residue_plddt"] = [round(float(p), 3) for p in atom_plddts[:500]]

        if 'pae' in conf:
            pae = conf['pae']
            pae_list = [[round(float(x), 2) for x in row[:100]] for row in pae[:100]]
            confidences["pae_matrix"] = pae_list

    return confidences if confidences else None


# ============== MPNN Inference ==============

def run_mpnn_inference(
    pdb_content: str,
    num_sequences: int = 8,
    temperature: float = 0.1,
    model_type: str = "ligand_mpnn",
    remove_waters: bool = True,
    use_mock: bool = False
) -> Dict[str, Any]:
    """Run MPNN sequence design"""

    if use_mock:
        return run_mpnn_mock(pdb_content, num_sequences)

    try:
        return run_mpnn_python_api(
            pdb_content, num_sequences, temperature, model_type, remove_waters
        )
    except ImportError as e:
        print(f"[MPNN] Python API not available: {e}")
        return run_mpnn_cli(pdb_content, num_sequences, temperature, model_type, remove_waters)
    except Exception as e:
        import traceback
        # Try CLI fallback
        try:
            return run_mpnn_cli(pdb_content, num_sequences, temperature, model_type, remove_waters)
        except Exception as cli_error:
            return {
                "status": "failed",
                "error": f"Python API: {str(e)}\nCLI fallback: {str(cli_error)}",
                "traceback": traceback.format_exc()
            }


def run_mpnn_python_api(
    pdb_content: str,
    num_sequences: int,
    temperature: float,
    model_type: str,
    remove_waters: bool
) -> Dict[str, Any]:
    """MPNN via Python API"""
    from mpnn.inference_engines.mpnn import MPNNInferenceEngine
    from biotite.structure import get_residue_starts
    from biotite.sequence import ProteinSequence

    # Parse input
    atom_array = pdb_to_atom_array(pdb_content)

    # Configure engine
    engine_config = {
        "model_type": model_type,
        "is_legacy_weights": True,
        "out_directory": None,
        "write_structures": False,
        "write_fasta": False,
    }

    input_configs = [{
        "batch_size": num_sequences,
        "remove_waters": remove_waters,
    }]

    # Run
    model = MPNNInferenceEngine(**engine_config)
    mpnn_outputs = model.run(input_dicts=input_configs, atom_arrays=[atom_array])

    # Extract sequences
    sequences = []
    for i, item in enumerate(mpnn_outputs):
        if hasattr(item, 'atom_array'):
            res_starts = get_residue_starts(item.atom_array)
            seq = ''.join(
                ProteinSequence.convert_letter_3to1(res_name)
                for res_name in item.atom_array.res_name[res_starts]
                if res_name in ProteinSequence._dict_3to1
            )
            sequences.append(f">design_{i+1}\n{seq}")

    fasta_content = "\n".join(sequences)

    return {
        "status": "completed",
        "result": {
            "sequences": [{"filename": "sequences.fasta", "content": fasta_content}],
            "num_sequences": len(sequences),
            "model_type": model_type,
            "mode": "real",
        }
    }


def run_mpnn_cli(
    pdb_content: str,
    num_sequences: int,
    temperature: float,
    model_type: str,
    remove_waters: bool
) -> Dict[str, Any]:
    """MPNN via CLI"""
    with tempfile.TemporaryDirectory() as tmpdir:
        pdb_path = os.path.join(tmpdir, "input.pdb")
        with open(pdb_path, "w") as f:
            f.write(pdb_content)

        out_dir = os.path.join(tmpdir, "output")
        os.makedirs(out_dir, exist_ok=True)

        cmd_parts = [
            "mpnn",
            f"--structure_path {pdb_path}",
            f"--out_directory {out_dir}",
            "--name design",
            f"--model_type {model_type}",
            f"--batch_size {num_sequences}",
            f"--temperature {temperature}",
            "--is_legacy_weights True",
            "--write_fasta True",
            "--write_structures False",
            f"--remove_waters {'True' if remove_waters else 'None'}",
        ]

        cmd = " ".join(cmd_parts)
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=1800)

        if result.returncode == 0:
            sequences = []
            for root, _, files in os.walk(out_dir):
                for filename in files:
                    if filename.endswith((".fa", ".fasta")):
                        with open(os.path.join(root, filename)) as f:
                            sequences.append({"filename": filename, "content": f.read()})

            if sequences:
                return {
                    "status": "completed",
                    "result": {
                        "sequences": sequences,
                        "model_type": model_type,
                        "mode": "real (cli)",
                    }
                }
            else:
                return {
                    "status": "failed",
                    "error": f"No output sequences found. stdout: {result.stdout}"
                }
        else:
            return {
                "status": "failed",
                "error": result.stderr or result.stdout or "Unknown error"
            }


def run_mpnn_mock(pdb_content: str, num_sequences: int) -> Dict[str, Any]:
    """Mock MPNN inference"""
    lines = pdb_content.split('\n')
    ca_count = sum(1 for l in lines if l.startswith('ATOM') and ' CA ' in l)
    length = ca_count if ca_count > 0 else 100

    return {
        "status": "completed",
        "result": {
            "sequences": [{
                "filename": "sequences.fasta",
                "content": generate_mock_fasta(num_sequences, length),
            }],
            "model_type": "mock",
            "mode": "mock",
        }
    }


# ============== RMSD Calculation ==============

def calculate_rmsd(pdb1: str, pdb2: str, backbone_only: bool = True) -> Dict[str, Any]:
    """Calculate RMSD between two structures"""
    try:
        from biotite.structure import rmsd, superimpose
        import numpy as np

        BACKBONE_ATOMS = ['N', 'CA', 'C', 'O']

        aa1 = pdb_to_atom_array(pdb1)
        aa2 = pdb_to_atom_array(pdb2)

        if backbone_only:
            mask1 = np.isin(aa1.atom_name, BACKBONE_ATOMS)
            mask2 = np.isin(aa2.atom_name, BACKBONE_ATOMS)
            aa1 = aa1[mask1]
            aa2 = aa2[mask2]

        if len(aa1) != len(aa2):
            min_len = min(len(aa1), len(aa2))
            aa1 = aa1[:min_len]
            aa2 = aa2[:min_len]

        aa2_fitted, _ = superimpose(aa1, aa2)
        rmsd_value = float(rmsd(aa1, aa2_fitted))

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
