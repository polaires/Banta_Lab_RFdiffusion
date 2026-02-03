"""
Structure Predictor Module (M6)

Predicts structures from designed sequences via RF3 or ESMFold.
Wraps inference_utils.run_rf3_inference() — does NOT rewrite it.

Reads: sequences (and optionally backbone_pdbs) from context
Writes: predictions to context
"""

import logging
from typing import Any, Dict, List, Optional

from pipeline_types import StepContext, PredictionResult

logger = logging.getLogger(__name__)


class StructurePredictor:
    """Predict structures for designed sequences using RF3 or ESMFold."""

    name = "structure_predictor"
    description = "Predict 3D structures from sequences using RoseTTAFold3 or ESMFold"
    input_keys = ["sequences"]
    output_keys = ["predictions"]
    optional_keys = ["backbone_pdbs"]

    def __init__(self, predictor: str = "rf3", **kwargs):
        """
        Args:
            predictor: "rf3" or "esmfold"
            **kwargs: Additional params (e.g., ligand_smiles for RF3)
        """
        self.predictor = predictor
        self.overrides = kwargs

    def run(self, context: StepContext) -> StepContext:
        """Run structure prediction for all sequences."""
        backend = context.backend
        if backend is None:
            raise RuntimeError("No inference backend configured on context")

        sequences = context.sequences
        if not sequences:
            raise ValueError("No sequences in context — run sequence_designer first")

        logger.info(f"Predicting structures for {len(sequences)} sequence(s) via {self.predictor}")

        predictions: List[PredictionResult] = []

        for seq_result in sequences:
            pred = self._predict_one(backend, context, seq_result)
            if pred:
                predictions.append(pred)

        # Sort by ptm descending for easy ranking
        predictions.sort(key=lambda p: p.ptm, reverse=True)

        logger.info(
            f"Completed {len(predictions)} prediction(s), "
            f"{sum(1 for p in predictions if p.passed)} passed"
        )
        context.predictions = predictions
        context.metadata["structure_predictor"] = {
            "predictor": self.predictor,
            "total_predictions": len(predictions),
            "passed": sum(1 for p in predictions if p.passed),
        }
        return context

    def _predict_one(
        self,
        backend: Any,
        context: StepContext,
        seq_result: Any,
    ) -> Optional[PredictionResult]:
        """Predict structure for a single sequence."""
        try:
            params = self._build_params(context, seq_result)
            result = backend.run_rf3(params)
            return _parse_rf3_result(result, seq_result, self.predictor)
        except Exception as e:
            logger.error(f"Prediction failed for seq {seq_result.filename}: {e}")
            return PredictionResult(
                sequence=seq_result.sequence,
                backbone_index=seq_result.backbone_index,
                predictor=self.predictor,
                metadata={"error": str(e)},
            )

    def _build_params(self, context: StepContext, seq_result: Any) -> Dict[str, Any]:
        """Build RF3 API params for a single sequence."""
        params: Dict[str, Any] = {
            "sequence": seq_result.sequence,
            "name": seq_result.filename or "prediction",
        }

        # Add ligand SMILES if available (for protein-ligand ipTM scoring)
        ligand_smiles = self.overrides.get(
            "ligand_smiles",
            context.get_param("ligand_smiles"),
        )
        if ligand_smiles:
            params["ligand_smiles"] = ligand_smiles

        # Add reference PDB for RMSD comparison if available
        if seq_result.pdb_content:
            params["pdb_content"] = seq_result.pdb_content

        return params


def _parse_rf3_result(
    result: Dict[str, Any],
    seq_result: Any,
    predictor: str,
) -> PredictionResult:
    """Parse RF3/ESMFold result into PredictionResult."""
    predictions = result.get("predictions", [{}])
    pred = predictions[0] if predictions else {}

    # RF3 stores confidences in result["confidences"], ESMFold may inline them
    conf = result.get("confidences") or {}
    summary = conf.get("summary_confidences", conf)

    plddt = float(
        summary.get("overall_plddt")
        or pred.get("mean_plddt")
        or result.get("mean_plddt")
        or 0.0
    )
    ptm = float(
        summary.get("ptm")
        or pred.get("ptm")
        or result.get("ptm")
        or 0.0
    )
    iptm = summary.get("iptm") or pred.get("iptm") or result.get("iptm")

    # RF3 returns "content" key; some paths use "pdb_content" — check both
    pdb_content = (
        pred.get("content")
        or pred.get("pdb_content")
        or result.get("pdb_content")
        or result.get("content")
    )

    # Determine pass/fail (basic thresholds, analyzer module does detailed filtering)
    passed = ptm >= 0.6 and plddt >= 0.65

    return PredictionResult(
        sequence=seq_result.sequence,
        predicted_pdb=pdb_content,
        plddt=plddt,
        ptm=ptm,
        iptm=float(iptm) if iptm is not None else None,
        passed=passed,
        predictor=predictor,
        backbone_index=seq_result.backbone_index,
        sequence_index=getattr(seq_result, "sequence_index", 0),
        metadata={
            "plddt_per_residue": pred.get("plddt_per_residue"),
            "chain_pair_pae": pred.get("chain_pair_pae"),
        },
    )
