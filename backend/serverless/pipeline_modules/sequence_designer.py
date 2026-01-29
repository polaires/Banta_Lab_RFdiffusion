"""
Sequence Designer Module (M5)

Designs sequences for backbones via LigandMPNN.
Wraps inference_utils.run_mpnn_inference() — does NOT rewrite it.

Reads: backbone_pdbs, mpnn_config (or params) from context
Writes: sequences to context
"""

import logging
import re
from typing import Any, Dict, List, Optional

from pipeline_types import StepContext, SequenceResult

logger = logging.getLogger(__name__)


class SequenceDesigner:
    """Design sequences for protein backbones using LigandMPNN."""

    name = "sequence_designer"
    description = "Design amino acid sequences for backbone structures using LigandMPNN"
    input_keys = ["backbone_pdbs"]
    output_keys = ["sequences"]
    optional_keys = ["mpnn_config"]

    def __init__(self, **kwargs):
        self.overrides = kwargs

    def run(self, context: StepContext) -> StepContext:
        """Execute MPNN sequence design for all backbones."""
        backend = context.backend
        if backend is None:
            raise RuntimeError("No inference backend configured on context")

        backbone_pdbs = context.backbone_pdbs
        if not backbone_pdbs:
            raise ValueError("No backbone PDBs in context — run backbone_generator first")

        num_sequences = self.overrides.get(
            "num_sequences",
            context.get_param("num_sequences", 8),
        )
        temperature = self.overrides.get(
            "temperature",
            context.get_param("temperature", 0.1),
        )

        logger.info(
            f"Designing {num_sequences} sequence(s) each for "
            f"{len(backbone_pdbs)} backbone(s)"
        )

        all_sequences: List[SequenceResult] = []

        for bb_idx, backbone_pdb in enumerate(backbone_pdbs):
            api_params = self._build_params(context, backbone_pdb, num_sequences, temperature)
            result = backend.run_mpnn(api_params)

            # Parse FASTA output from MPNN
            seqs = _parse_mpnn_result(result, bb_idx)
            all_sequences.extend(seqs)

        logger.info(f"Designed {len(all_sequences)} total sequence(s)")
        context.sequences = all_sequences
        context.metadata["sequence_designer"] = {
            "num_sequences_per_backbone": num_sequences,
            "total_sequences": len(all_sequences),
            "num_backbones": len(backbone_pdbs),
        }
        return context

    def _build_params(
        self,
        context: StepContext,
        backbone_pdb: str,
        num_sequences: int,
        temperature: float,
    ) -> Dict[str, Any]:
        """Build MPNN API params."""
        params: Dict[str, Any] = {
            "pdb_content": backbone_pdb,
            "num_sequences": num_sequences,
            "temperature": temperature,
            "model_type": "ligand_mpnn",
        }

        # Merge mpnn_config from context
        if context.mpnn_config:
            if isinstance(context.mpnn_config, dict):
                cfg = context.mpnn_config
            elif hasattr(context.mpnn_config, "to_api_params"):
                cfg = context.mpnn_config.to_api_params()
            else:
                cfg = {}
            # Config overrides defaults, but not pdb_content
            for k, v in cfg.items():
                if k not in ("pdb_content", "task") and v is not None:
                    params[k] = v

        # Context-level params
        for key in ("bias_AA", "omit_AA", "fixed_positions",
                     "pack_side_chains", "use_side_chain_context"):
            val = context.get_param(key)
            if val is not None:
                params[key] = val

        # Constructor overrides
        for k, v in self.overrides.items():
            if k not in ("pdb_content",) and v is not None:
                params[k] = v

        return params


def _parse_mpnn_result(
    result: Dict[str, Any],
    backbone_index: int,
) -> List[SequenceResult]:
    """Parse MPNN result into SequenceResult list.

    MPNN returns either:
    - {"sequences": [{"sequence": "...", "score": ...}, ...]}
    - {"fasta": ">seq_0001\\nMEVTI...\\n>seq_0002\\nMKILP..."}
    """
    sequences: List[SequenceResult] = []

    # Try structured sequences first
    seq_list = result.get("sequences", [])
    if seq_list:
        for i, entry in enumerate(seq_list):
            if isinstance(entry, dict):
                seq = entry.get("sequence", "")
                score = entry.get("score", entry.get("global_score", 0.0))
                pdb = entry.get("pdb_content")
            elif isinstance(entry, str):
                seq = entry
                score = 0.0
                pdb = None
            else:
                continue

            if seq:
                sequences.append(SequenceResult(
                    sequence=seq,
                    score=float(score),
                    backbone_index=backbone_index,
                    filename=f"seq_{backbone_index:03d}_{i:04d}",
                    pdb_content=pdb,
                ))
        return sequences

    # Fallback: parse FASTA string
    fasta = result.get("fasta", "")
    if fasta:
        return _parse_fasta(fasta, backbone_index)

    logger.warning(f"No sequences found in MPNN result for backbone {backbone_index}")
    return []


def _parse_fasta(fasta: str, backbone_index: int) -> List[SequenceResult]:
    """Parse a FASTA string into SequenceResult list."""
    sequences: List[SequenceResult] = []
    current_header = ""
    current_seq_lines: List[str] = []

    for line in fasta.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            # Save previous sequence
            if current_seq_lines:
                seq = "".join(current_seq_lines)
                score = _extract_score_from_header(current_header)
                sequences.append(SequenceResult(
                    sequence=seq,
                    score=score,
                    backbone_index=backbone_index,
                    filename=f"seq_{backbone_index:03d}_{len(sequences):04d}",
                ))
            current_header = line[1:]
            current_seq_lines = []
        else:
            current_seq_lines.append(line)

    # Don't forget last entry
    if current_seq_lines:
        seq = "".join(current_seq_lines)
        score = _extract_score_from_header(current_header)
        sequences.append(SequenceResult(
            sequence=seq,
            score=score,
            backbone_index=backbone_index,
            filename=f"seq_{backbone_index:03d}_{len(sequences):04d}",
        ))

    return sequences


def _extract_score_from_header(header: str) -> float:
    """Extract score from FASTA header (e.g., 'seq_0001, score=-1.234')."""
    match = re.search(r"score[=:]\s*([-\d.]+)", header, re.IGNORECASE)
    if match:
        try:
            return float(match.group(1))
        except ValueError:
            pass
    return 0.0
