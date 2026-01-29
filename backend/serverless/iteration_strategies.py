"""
Iteration Strategies — Scout and Sweep patterns for pipeline workflows.

These are composable strategies that the WorkflowRunner inserts between
pipeline steps to implement multi-pass design patterns.
"""

import logging
from typing import Any, Dict, List, Optional

from pipeline_types import StepContext, SequenceResult, PredictionResult

logger = logging.getLogger(__name__)


class ScoutStrategy:
    """
    Scout strategy: generate 1 sequence per backbone, validate,
    skip backbones that fail, then let downstream steps run on
    the filtered set.

    Inserted between backbone_generator and sequence_designer by
    WorkflowRunner when strategy="scout".
    """

    name = "scout_filter"
    description = "Scout: validate 1 sequence per backbone and filter bad backbones"
    input_keys = ["backbone_pdbs"]
    output_keys = ["backbone_pdbs"]
    optional_keys = []

    def __init__(self, ptm_threshold: float = 0.6, plddt_threshold: float = 0.65):
        self.ptm_threshold = ptm_threshold
        self.plddt_threshold = plddt_threshold

    def run(self, context: StepContext) -> StepContext:
        """Run scout round: 1 seq per backbone → validate → filter."""
        backend = context.backend
        if backend is None:
            raise RuntimeError("No inference backend configured")

        backbone_pdbs = context.backbone_pdbs
        if not backbone_pdbs:
            return context

        logger.info(f"Scout: evaluating {len(backbone_pdbs)} backbone(s) "
                    f"(ptm≥{self.ptm_threshold}, plddt≥{self.plddt_threshold})")

        passing_pdbs: List[str] = []
        scout_sequences: List[SequenceResult] = []
        scout_predictions: List[PredictionResult] = []

        for i, backbone_pdb in enumerate(backbone_pdbs):
            # Generate 1 scout sequence
            try:
                mpnn_params = self._build_scout_mpnn_params(context, backbone_pdb)
                mpnn_result = backend.run_mpnn(mpnn_params)
                seqs = mpnn_result.get("sequences", [])
                if not seqs:
                    fasta = mpnn_result.get("fasta", "")
                    if fasta:
                        # Quick FASTA parse for scout
                        seq_text = ""
                        for line in fasta.splitlines():
                            if not line.startswith(">") and line.strip():
                                seq_text += line.strip()
                                break
                        seqs = [{"sequence": seq_text}] if seq_text else []

                if not seqs:
                    logger.debug(f"Scout: backbone {i} — no sequence generated, skipping")
                    continue

                sequence = seqs[0].get("sequence", seqs[0]) if isinstance(seqs[0], dict) else seqs[0]
            except Exception as e:
                logger.warning(f"Scout: backbone {i} — MPNN failed: {e}")
                continue

            # Validate scout sequence
            try:
                rf3_params = {"sequence": sequence, "name": f"scout_{i:03d}"}
                ligand_smiles = context.get_param("ligand_smiles")
                if ligand_smiles:
                    rf3_params["ligand_smiles"] = ligand_smiles

                rf3_result = backend.run_rf3(rf3_params)
                preds = rf3_result.get("predictions", [{}])
                pred = preds[0] if preds else rf3_result
                ptm = float(pred.get("ptm", rf3_result.get("ptm", 0.0)))
                plddt = float(pred.get("mean_plddt", rf3_result.get("mean_plddt", 0.0)))
            except Exception as e:
                logger.warning(f"Scout: backbone {i} — RF3 failed: {e}")
                continue

            # Check thresholds
            if ptm >= self.ptm_threshold and plddt >= self.plddt_threshold:
                passing_pdbs.append(backbone_pdb)
                scout_sequences.append(SequenceResult(
                    sequence=sequence,
                    confidence=ptm,
                    backbone_index=i,
                ))
                scout_predictions.append(PredictionResult(
                    sequence=sequence,
                    ptm=ptm,
                    plddt=plddt,
                    passed=True,
                    predictor="rf3",
                    backbone_index=i,
                ))
                logger.debug(f"Scout: backbone {i} PASSED (ptm={ptm:.3f}, plddt={plddt:.3f})")
            else:
                logger.debug(f"Scout: backbone {i} FAILED (ptm={ptm:.3f}, plddt={plddt:.3f})")

        logger.info(
            f"Scout complete: {len(passing_pdbs)}/{len(backbone_pdbs)} backbones passed"
        )

        # Update context with filtered backbones
        context.backbone_pdbs = passing_pdbs
        context.metadata["scout_filter"] = {
            "original_backbones": len(backbone_pdbs),
            "passing_backbones": len(passing_pdbs),
            "ptm_threshold": self.ptm_threshold,
            "plddt_threshold": self.plddt_threshold,
        }

        return context

    def _build_scout_mpnn_params(self, context: StepContext, backbone_pdb: str) -> Dict[str, Any]:
        """Build MPNN params for a single scout sequence."""
        params: Dict[str, Any] = {
            "pdb_content": backbone_pdb,
            "num_sequences": 1,
            "temperature": 0.1,
            "model_type": "ligand_mpnn",
        }
        # Copy bias/omit from context mpnn_config
        if context.mpnn_config:
            cfg = context.mpnn_config
            if isinstance(cfg, dict):
                for k in ("bias_AA", "omit_AA", "fixed_positions"):
                    if k in cfg:
                        params[k] = cfg[k]
            elif hasattr(cfg, "to_api_params"):
                api = cfg.to_api_params()
                for k in ("bias_AA", "omit_AA", "fixed_positions"):
                    if k in api:
                        params[k] = api[k]
        return params


class SweepStrategy:
    """
    Sweep strategy: try multiple configurations, pick best, then
    optionally do a production run with the winner.

    Inserted before backbone_generator by WorkflowRunner when
    strategy="sweep".
    """

    name = "sweep_configurator"
    description = "Sweep: try multiple RFD3 configs and select the best"
    input_keys = ["rfd3_config"]
    output_keys = ["rfd3_config"]
    optional_keys = []

    def __init__(
        self,
        grid: Optional[Dict[str, List[Any]]] = None,
        designs_per_config: int = 10,
    ):
        """
        Args:
            grid: Dict of param_name → list of values to sweep.
                  e.g., {"cfg_scale": [1.5, 2.0, 2.5], "scaffold_size": ["120-150", "150-200"]}
            designs_per_config: Number of designs per sweep config.
        """
        self.grid = grid or {}
        self.designs_per_config = designs_per_config

    def run(self, context: StepContext) -> StepContext:
        """Generate sweep configs and store them for backbone_generator to iterate."""
        configs = self.generate_grid(context)

        if not configs:
            logger.warning("Sweep: no configs generated — using defaults")
            return context

        logger.info(f"Sweep: {len(configs)} configuration(s) to evaluate")

        # Store sweep configs in metadata for backbone_generator
        context.metadata["sweep_configs"] = configs
        context.metadata["sweep_designs_per_config"] = self.designs_per_config

        # Set first config as active (backbone_generator will iterate if sweep-aware)
        if context.rfd3_config and isinstance(context.rfd3_config, dict):
            context.rfd3_config.update(configs[0])
        elif hasattr(context.rfd3_config, "__dict__"):
            for k, v in configs[0].items():
                if hasattr(context.rfd3_config, k):
                    setattr(context.rfd3_config, k, v)

        context.metadata["sweep_configurator"] = {
            "num_configs": len(configs),
            "grid": self.grid,
            "designs_per_config": self.designs_per_config,
        }
        return context

    def generate_grid(self, context: StepContext = None) -> List[Dict[str, Any]]:
        """Generate all combinations from the grid."""
        if not self.grid:
            return []

        keys = list(self.grid.keys())
        values = list(self.grid.values())

        # Cartesian product
        configs: List[Dict[str, Any]] = [{}]
        for key, vals in zip(keys, values):
            new_configs = []
            for config in configs:
                for val in vals:
                    new_config = dict(config)
                    new_config[key] = val
                    new_configs.append(new_config)
            configs = new_configs

        return configs

    @staticmethod
    def generate_sweep_grid(
        sizes: List[str],
        cfg_scales: List[float],
        base_params: Optional[Dict[str, Any]] = None,
    ) -> List[Dict[str, Any]]:
        """
        Convenience: generate a standard scaffold size × cfg_scale grid.

        Args:
            sizes: List of scaffold size ranges (e.g., ["120-150", "150-200"])
            cfg_scales: List of CFG scale values (e.g., [1.5, 2.0, 2.5])
            base_params: Base params to include in every config

        Returns:
            List of config dicts
        """
        base = base_params or {}
        configs = []
        for size in sizes:
            for cfg in cfg_scales:
                config = dict(base)
                config["length"] = size
                config["cfg_scale"] = cfg
                configs.append(config)
        return configs
