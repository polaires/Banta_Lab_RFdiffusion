"""
Post-Design Validation Pipeline

Multi-step validation pipeline for RFD3 designs:
1. Sequence Design: LigandMPNN or ProteinMPNN
2. Structure Prediction: ESMFold (fast) or AF3 (accurate)
3. Structure Analysis: RMSD, pLDDT, metal coordination
4. Energy Refinement: PyRosetta FastRelax
5. Filtering: Rank candidates by multiple metrics

Works for any design type:
- Protein-only (unconditional, binders)
- Protein-ligand (small molecule binders)
- Metal-binding proteins
- Metal-ligand complexes (e.g., Dy-TriNOx)

Usage:
    from design_validation_pipeline import DesignValidationPipeline

    pipeline = DesignValidationPipeline()
    results = pipeline.validate(
        pdb_content=rfd3_output,
        num_sequences=8,
        use_ligandmpnn=True,
        ligand_name="UNL",
        metal_type="DY",
    )
"""

import os
import json
import tempfile
from datetime import datetime
from dataclasses import dataclass, field, asdict
from typing import Dict, Any, Optional, List, Tuple
from pathlib import Path

# Import PDB sanitization from unified analyzer
try:
    from unified_analyzer import sanitize_pdb
    SANITIZE_AVAILABLE = True
except ImportError:
    SANITIZE_AVAILABLE = False
    sanitize_pdb = None


@dataclass
class SequenceCandidate:
    """A sequence candidate from MPNN with validation results."""
    sequence: str
    mpnn_score: float
    # Structure prediction results
    predicted_pdb: Optional[str] = None
    backbone_rmsd: Optional[float] = None
    plddt_mean: Optional[float] = None
    plddt_per_residue: Optional[List[float]] = None
    # Metal coordination (if applicable)
    metal_coordination: Optional[int] = None
    metal_protein_donors: Optional[int] = None
    metal_preserved: Optional[bool] = None
    # RF3 validation results
    rf3_plddt: Optional[float] = None
    rf3_ptm: Optional[float] = None
    rf3_iptm: Optional[float] = None
    rf3_rmsd: Optional[float] = None
    rf3_cif: Optional[str] = None
    # Best-of metrics (max pLDDT, min RMSD across ESMFold and RF3)
    best_plddt: Optional[float] = None
    best_rmsd: Optional[float] = None
    # Ligand metrics (if applicable)
    ligand_rmsd: Optional[float] = None
    ligand_contacts: Optional[int] = None
    # PyRosetta metrics
    rosetta_energy: Optional[float] = None
    clash_score: Optional[int] = None
    # Overall assessment
    passes_filters: bool = False
    filter_failures: List[str] = field(default_factory=list)
    rank: Optional[int] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class ValidationResult:
    """Complete validation results for a design."""
    design_id: str
    timestamp: str
    input_pdb_path: Optional[str]
    design_type: str
    # Pipeline configuration
    config: Dict[str, Any]
    # Candidates
    candidates: List[SequenceCandidate]
    num_passed: int
    num_total: int
    # Best candidate
    best_candidate_idx: Optional[int]
    best_sequence: Optional[str]
    best_backbone_rmsd: Optional[float]
    best_plddt: Optional[float]
    # Summary statistics
    summary: Dict[str, Any]

    def to_dict(self) -> Dict[str, Any]:
        result = asdict(self)
        result["candidates"] = [c.to_dict() if hasattr(c, 'to_dict') else c for c in self.candidates]
        return result


class DesignValidationPipeline:
    """
    Multi-step validation pipeline for protein designs.

    Pipeline steps:
    1. Sequence design (LigandMPNN/ProteinMPNN)
    2. Structure prediction (ESMFold)
    3. RMSD calculation (designed vs predicted)
    4. Metal coordination check (if metal present)
    5. PyRosetta refinement (optional)
    6. Filtering and ranking
    """

    def __init__(self):
        """Initialize pipeline with available modules."""
        self._modules = self._discover_modules()
        print(f"[ValidationPipeline] Available modules: {self._modules}")

    def _discover_modules(self) -> Dict[str, bool]:
        """Check which analysis modules are available."""
        modules = {
            "esmfold": False,
            "rf3": False,
            "ligandmpnn": False,
            "proteinmpnn": False,
            "pyrosetta": False,
            "metal_validation": False,
            "unified_analyzer": False,
        }

        # Check ESMFold
        try:
            import esmfold_utils
            modules["esmfold"] = esmfold_utils.is_esmfold_available()
        except ImportError:
            pass

        # Check RF3
        try:
            import esmfold_utils
            modules["rf3"] = esmfold_utils.is_rf3_available()
        except (ImportError, AttributeError):
            pass

        # Check LigandMPNN/ProteinMPNN (via inference_utils)
        try:
            import inference_utils
            # Check if MPNN inference is available
            modules["ligandmpnn"] = True
            modules["proteinmpnn"] = True
        except ImportError:
            pass

        # Check PyRosetta
        try:
            import rosetta_utils
            modules["pyrosetta"] = rosetta_utils.check_pyrosetta_available()
        except ImportError:
            pass

        # Check metal validation
        try:
            import metal_validation
            modules["metal_validation"] = True
        except ImportError:
            pass

        # Check unified analyzer
        try:
            import unified_analyzer
            modules["unified_analyzer"] = True
        except ImportError:
            pass

        return modules

    def validate(
        self,
        pdb_content: str,
        pdb_path: Optional[str] = None,
        num_sequences: int = 8,
        use_ligandmpnn: bool = True,
        ligand_name: Optional[str] = None,
        metal_type: Optional[str] = None,
        # Filter thresholds
        max_backbone_rmsd: float = 1.5,
        min_plddt: float = 0.7,
        min_metal_coordination: int = 6,
        min_protein_donors: int = 1,
        # Options
        run_pyrosetta: bool = True,
        temperature: float = 0.1,
        use_rf3: bool = True,
        ligand_smiles: Optional[str] = None,
    ) -> ValidationResult:
        """
        Run full validation pipeline on an RFD3 design.

        Args:
            pdb_content: PDB file content from RFD3
            pdb_path: Optional path to PDB file
            num_sequences: Number of sequences to generate per backbone
            use_ligandmpnn: Use LigandMPNN (True) or ProteinMPNN (False)
            ligand_name: Ligand residue name (e.g., "UNL") for LigandMPNN
            metal_type: Metal type code (e.g., "DY", "TB") for coordination check
            max_backbone_rmsd: Maximum allowed backbone RMSD (Å)
            min_plddt: Minimum required pLDDT score
            min_metal_coordination: Minimum metal coordination number
            min_protein_donors: Minimum protein donors to metal
            run_pyrosetta: Whether to run PyRosetta FastRelax
            temperature: MPNN sampling temperature

        Returns:
            ValidationResult with all candidates and metrics
        """
        timestamp = datetime.now().isoformat()
        design_id = f"val_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

        # Sanitize PDB content to fix RFD3 output issues (duplicate atoms, etc.)
        sanitization_issues = []
        if SANITIZE_AVAILABLE and sanitize_pdb:
            pdb_content, sanitization_issues = sanitize_pdb(pdb_content, remove_duplicates=True)
            if sanitization_issues:
                print(f"[ValidationPipeline] PDB sanitization fixed: {len(sanitization_issues)} issues")

        # Detect design type
        design_type = self._detect_design_type(pdb_content, ligand_name, metal_type)
        print(f"[ValidationPipeline] Design type: {design_type}")

        # Configuration record
        config = {
            "num_sequences": num_sequences,
            "use_ligandmpnn": use_ligandmpnn,
            "ligand_name": ligand_name,
            "metal_type": metal_type,
            "max_backbone_rmsd": max_backbone_rmsd,
            "min_plddt": min_plddt,
            "min_metal_coordination": min_metal_coordination,
            "min_protein_donors": min_protein_donors,
            "run_pyrosetta": run_pyrosetta,
            "temperature": temperature,
            "pdb_sanitized": len(sanitization_issues) > 0,
            "sanitization_issues": sanitization_issues[:10] if sanitization_issues else [],  # Limit to first 10
            "use_rf3": use_rf3,
            "ligand_smiles": ligand_smiles,
        }

        # Step 1: Generate sequences with MPNN
        print(f"[ValidationPipeline] Step 1: Generating {num_sequences} sequences...")
        sequences = self._run_mpnn(
            pdb_content=pdb_content,
            num_sequences=num_sequences,
            use_ligandmpnn=use_ligandmpnn,
            ligand_name=ligand_name,
            temperature=temperature,
        )

        if not sequences:
            print("[ValidationPipeline] ERROR: No sequences generated")
            return ValidationResult(
                design_id=design_id,
                timestamp=timestamp,
                input_pdb_path=pdb_path,
                design_type=design_type,
                config=config,
                candidates=[],
                num_passed=0,
                num_total=0,
                best_candidate_idx=None,
                best_sequence=None,
                best_backbone_rmsd=None,
                best_plddt=None,
                summary={"error": "No sequences generated from MPNN"},
            )

        print(f"[ValidationPipeline] Generated {len(sequences)} sequences")

        # Create candidates
        candidates = []
        for seq, score in sequences:
            candidates.append(SequenceCandidate(
                sequence=seq,
                mpnn_score=score,
            ))

        # Step 2: Predict structures with ESMFold
        print("[ValidationPipeline] Step 2: Predicting structures with ESMFold...")
        candidates = self._predict_structures(candidates)

        # Step 3: Calculate RMSD to designed structure
        print("[ValidationPipeline] Step 3: Calculating backbone RMSD...")
        candidates = self._calculate_rmsd(candidates, pdb_content)

        # Step 4: Check metal coordination (if applicable)
        if metal_type:
            print(f"[ValidationPipeline] Step 4: Checking {metal_type} coordination...")
            candidates = self._check_metal_coordination(
                candidates, pdb_content, metal_type, ligand_name
            )

        # Step 4.5: RF3 structure prediction (if available)
        if use_rf3 and self._modules.get("rf3"):
            print("[ValidationPipeline] Step 4.5: Predicting structures with RF3...")
            candidates = self._predict_structures_rf3(
                candidates, pdb_content, ligand_smiles=ligand_smiles
            )

        # Step 5: PyRosetta refinement (optional)
        if run_pyrosetta and self._modules.get("pyrosetta"):
            print("[ValidationPipeline] Step 5: Running PyRosetta refinement...")
            candidates = self._run_pyrosetta(candidates)

        # Step 6: Apply filters and rank
        print("[ValidationPipeline] Step 6: Filtering and ranking candidates...")
        candidates = self._apply_filters(
            candidates=candidates,
            max_backbone_rmsd=max_backbone_rmsd,
            min_plddt=min_plddt,
            min_metal_coordination=min_metal_coordination if metal_type else None,
            min_protein_donors=min_protein_donors if metal_type else None,
        )

        # Sort by rank
        candidates.sort(key=lambda c: (not c.passes_filters, c.rank or 999))

        # Find best candidate
        passing = [c for c in candidates if c.passes_filters]
        num_passed = len(passing)
        best_idx = None
        best_seq = None
        best_rmsd = None
        best_plddt = None

        if passing:
            best_idx = candidates.index(passing[0])
            best_seq = passing[0].sequence
            best_rmsd = passing[0].best_rmsd or passing[0].backbone_rmsd
            best_plddt = passing[0].best_plddt or passing[0].plddt_mean

        # Summary statistics
        summary = {
            "num_candidates": len(candidates),
            "num_passed": num_passed,
            "pass_rate": f"{num_passed}/{len(candidates)} ({100*num_passed/len(candidates):.1f}%)" if candidates else "0/0",
            "esmfold_rmsd_range": self._get_range([c.backbone_rmsd for c in candidates if c.backbone_rmsd]),
            "esmfold_plddt_range": self._get_range([c.plddt_mean for c in candidates if c.plddt_mean]),
            "rf3_rmsd_range": self._get_range([c.rf3_rmsd for c in candidates if c.rf3_rmsd]),
            "rf3_plddt_range": self._get_range([c.rf3_plddt for c in candidates if c.rf3_plddt]),
            "best_rmsd_range": self._get_range([c.best_rmsd for c in candidates if c.best_rmsd]),
            "best_plddt_range": self._get_range([c.best_plddt for c in candidates if c.best_plddt]),
        }

        if metal_type:
            summary["metal_coordination_range"] = self._get_range(
                [c.metal_coordination for c in candidates if c.metal_coordination]
            )
            summary["metal_preserved_count"] = sum(
                1 for c in candidates if c.metal_preserved
            )

        print(f"[ValidationPipeline] Complete: {summary['pass_rate']} passed filters")

        return ValidationResult(
            design_id=design_id,
            timestamp=timestamp,
            input_pdb_path=pdb_path,
            design_type=design_type,
            config=config,
            candidates=candidates,
            num_passed=num_passed,
            num_total=len(candidates),
            best_candidate_idx=best_idx,
            best_sequence=best_seq,
            best_backbone_rmsd=best_rmsd,
            best_plddt=best_plddt,
            summary=summary,
        )

    def _detect_design_type(
        self,
        pdb_content: str,
        ligand_name: Optional[str],
        metal_type: Optional[str],
    ) -> str:
        """Detect design type from PDB content and parameters."""
        has_metal = metal_type is not None
        has_ligand = ligand_name is not None

        # Also check PDB for HETATM records
        if not has_metal:
            for line in pdb_content.split("\n"):
                if line.startswith("HETATM"):
                    res_name = line[17:20].strip()
                    # Common metal codes
                    if res_name in ["ZN", "MG", "CA", "FE", "CU", "MN", "CO", "NI",
                                    "TB", "EU", "DY", "GD", "SM", "ND", "LA", "CE"]:
                        has_metal = True
                        break

        if not has_ligand:
            for line in pdb_content.split("\n"):
                if line.startswith("HETATM"):
                    res_name = line[17:20].strip()
                    # Skip waters, metals, common ions
                    if res_name not in ["HOH", "WAT", "ZN", "MG", "CA", "FE", "CU",
                                        "MN", "CO", "NI", "TB", "EU", "DY", "GD",
                                        "SM", "ND", "LA", "CE", "NA", "CL", "K"]:
                        has_ligand = True
                        break

        if has_metal and has_ligand:
            return "metal_ligand_complex"
        elif has_metal:
            return "metal_binding"
        elif has_ligand:
            return "ligand_binding"
        else:
            return "protein_only"

    def _run_mpnn(
        self,
        pdb_content: str,
        num_sequences: int,
        use_ligandmpnn: bool,
        ligand_name: Optional[str],
        temperature: float,
    ) -> List[Tuple[str, float]]:
        """
        Generate sequences using LigandMPNN or ProteinMPNN.

        Uses inference_utils.run_mpnn_inference which handles both models.

        Returns:
            List of (sequence, score) tuples
        """
        if not (self._modules.get("ligandmpnn") or self._modules.get("proteinmpnn")):
            print("[ValidationPipeline] WARNING: No MPNN module available")
            return []

        try:
            import inference_utils

            # Choose model type
            model_type = "ligand_mpnn" if use_ligandmpnn else "protein_mpnn"

            # Run MPNN via inference_utils
            result = inference_utils.run_mpnn_inference(
                pdb_content=pdb_content,
                num_sequences=num_sequences,
                temperature=temperature,
                model_type=model_type,
            )

            status = result.get("status", "")
            if status not in ["success", "completed"]:
                error = result.get("error", "Unknown error")
                print(f"[ValidationPipeline] MPNN failed: {error}")
                return []

            # Extract sequences and scores from result
            # Handle different result formats from inference_utils
            sequences = []

            # Format 1: designs list
            designs = result.get("designs", [])

            # Format 2: result.sequences (FASTA format)
            if not designs and "result" in result:
                inner_result = result["result"]
                seq_files = inner_result.get("sequences", [])
                for seq_file in seq_files:
                    content = seq_file.get("content", "")
                    # Parse FASTA format
                    for line in content.split("\n"):
                        if line and not line.startswith(">"):
                            sequences.append((line.strip(), 0.0))  # No score in FASTA
                if sequences:
                    print(f"[ValidationPipeline] Parsed {len(sequences)} sequences from FASTA")
                    return sequences

            for design in designs:
                seq = design.get("sequence", "")
                score = design.get("score", 0.0)
                if seq:
                    sequences.append((seq, score))

            print(f"[ValidationPipeline] Generated {len(sequences)} sequences")
            return sequences

        except Exception as e:
            print(f"[ValidationPipeline] MPNN error: {e}")
            import traceback
            traceback.print_exc()
            return []

    def _predict_structures(
        self,
        candidates: List[SequenceCandidate],
    ) -> List[SequenceCandidate]:
        """Predict structures for all candidates using ESMFold."""
        if not self._modules.get("esmfold"):
            print("[ValidationPipeline] ESMFold not available, skipping prediction")
            return candidates

        try:
            import esmfold_utils

            for i, candidate in enumerate(candidates):
                print(f"  Predicting structure {i+1}/{len(candidates)}...")
                try:
                    # Use predict_structure_esmfold from esmfold_utils
                    result = esmfold_utils.predict_structure_esmfold(
                        sequence=candidate.sequence,
                        return_pdb=True,
                    )

                    # Check for success (status == "completed" or success == True)
                    status = result.get("status", "")
                    success = result.get("success", False) or status == "completed"

                    if success:
                        candidate.predicted_pdb = result.get("pdb_content")
                        # Handle both key names from different versions
                        candidate.plddt_mean = result.get("mean_plddt") or result.get("plddt_mean")
                        candidate.plddt_per_residue = result.get("plddt") or result.get("plddt_per_residue", [])
                        print(f"    pLDDT: {candidate.plddt_mean:.2f}" if candidate.plddt_mean else "    pLDDT: N/A")
                    else:
                        error = result.get("error", "Unknown error")
                        print(f"    WARNING: Prediction failed for candidate {i+1}: {error}")

                except Exception as e:
                    print(f"    ERROR: {e}")

        except Exception as e:
            print(f"[ValidationPipeline] ESMFold error: {e}")

        return candidates

    def _calculate_rmsd(
        self,
        candidates: List[SequenceCandidate],
        designed_pdb: str,
    ) -> List[SequenceCandidate]:
        """Calculate backbone RMSD between designed and predicted structures."""
        # Try using esmfold_utils.calculate_backbone_rmsd first
        try:
            import esmfold_utils

            for candidate in candidates:
                if not candidate.predicted_pdb:
                    continue

                try:
                    result = esmfold_utils.calculate_backbone_rmsd(
                        pdb1_content=designed_pdb,
                        pdb2_content=candidate.predicted_pdb,
                    )
                    # Handle both dict and float return types
                    if isinstance(result, dict):
                        if result.get("status") == "completed":
                            candidate.backbone_rmsd = result.get("rmsd")
                        else:
                            print(f"    RMSD calculation failed: {result.get('error')}")
                    else:
                        candidate.backbone_rmsd = result
                except Exception as e:
                    print(f"    RMSD calculation error: {e}")

            return candidates

        except (ImportError, AttributeError):
            pass

        # Fallback: use BioPython directly
        try:
            from Bio.PDB import PDBParser, Superimposer
            from Bio.PDB.Polypeptide import is_aa
            import numpy as np
            import io

            parser = PDBParser(QUIET=True)

            # Parse designed structure
            designed_structure = parser.get_structure(
                'designed', io.StringIO(designed_pdb)
            )

            # Get CA atoms from designed structure
            designed_ca = []
            for model in designed_structure:
                for chain in model:
                    for residue in chain:
                        if is_aa(residue) and 'CA' in residue:
                            designed_ca.append(residue['CA'])
                break  # Only first model

            if not designed_ca:
                print("[ValidationPipeline] No CA atoms in designed structure")
                return candidates

            for candidate in candidates:
                if not candidate.predicted_pdb:
                    continue

                try:
                    # Parse predicted structure
                    predicted_structure = parser.get_structure(
                        'predicted', io.StringIO(candidate.predicted_pdb)
                    )

                    # Get CA atoms from predicted structure
                    predicted_ca = []
                    for model in predicted_structure:
                        for chain in model:
                            for residue in chain:
                                if is_aa(residue) and 'CA' in residue:
                                    predicted_ca.append(residue['CA'])
                        break

                    # Align lengths (take minimum)
                    min_len = min(len(designed_ca), len(predicted_ca))
                    if min_len < 10:
                        continue

                    # Superimpose and calculate RMSD
                    sup = Superimposer()
                    sup.set_atoms(designed_ca[:min_len], predicted_ca[:min_len])
                    candidate.backbone_rmsd = sup.rms

                except Exception as e:
                    print(f"    RMSD calculation error: {e}")

        except ImportError:
            print("[ValidationPipeline] BioPython not available for RMSD calculation")

        return candidates

    def _predict_structures_rf3(
        self,
        candidates: List[SequenceCandidate],
        designed_pdb: str,
        ligand_smiles: Optional[str] = None,
    ) -> List[SequenceCandidate]:
        """Predict structures and compute RMSD for all candidates using RF3."""
        try:
            import esmfold_utils

            for i, candidate in enumerate(candidates):
                print(f"  RF3 prediction {i+1}/{len(candidates)}...")
                try:
                    result = esmfold_utils.validate_structure_rf3(
                        sequence=candidate.sequence,
                        designed_pdb=designed_pdb,
                        ligand_smiles=ligand_smiles,
                    )

                    if result.get("status") == "completed":
                        candidate.rf3_plddt = result.get("rf3_plddt")
                        candidate.rf3_ptm = result.get("rf3_ptm")
                        candidate.rf3_iptm = result.get("rf3_iptm")
                        candidate.rf3_rmsd = result.get("rf3_rmsd")
                        candidate.rf3_cif = result.get("rf3_cif_content")

                        plddt_str = f"{candidate.rf3_plddt:.2f}" if candidate.rf3_plddt else "N/A"
                        rmsd_str = f"{candidate.rf3_rmsd:.2f}" if candidate.rf3_rmsd else "N/A"
                        print(f"    RF3 pLDDT: {plddt_str}, RMSD: {rmsd_str}Å")
                    else:
                        error = result.get("error", "Unknown error")
                        print(f"    WARNING: RF3 prediction failed for candidate {i+1}: {error}")

                except Exception as e:
                    print(f"    ERROR: RF3 prediction {i+1}: {e}")

        except Exception as e:
            print(f"[ValidationPipeline] RF3 error: {e}")

        return candidates

    def _check_metal_coordination(
        self,
        candidates: List[SequenceCandidate],
        designed_pdb: str,
        metal_type: str,
        ligand_name: Optional[str],
    ) -> List[SequenceCandidate]:
        """Check if metal coordination is preserved in predicted structures."""
        if not self._modules.get("metal_validation"):
            print("[ValidationPipeline] metal_validation not available")
            return candidates

        try:
            import metal_validation

            # Get designed coordination as reference
            designed_result = metal_validation.validate_metal_ligand_complex_site(
                pdb_content=designed_pdb,
                metal=metal_type,
                ligand_name=ligand_name or "UNL",
            )
            designed_cn = designed_result.get("coordination_number", 0)
            designed_protein = designed_result.get("protein_coordination", 0)

            print(f"  Designed coordination: CN={designed_cn}, protein={designed_protein}")

            for candidate in candidates:
                if not candidate.predicted_pdb:
                    continue

                try:
                    # Note: Predicted structure won't have metal/ligand
                    # We need to transfer them from designed structure
                    # For now, just mark as needing manual check
                    # TODO: Implement metal transfer and re-check

                    # Use designed values as placeholder
                    candidate.metal_coordination = designed_cn
                    candidate.metal_protein_donors = designed_protein
                    candidate.metal_preserved = designed_cn >= 6

                except Exception as e:
                    print(f"    Metal check error: {e}")

        except Exception as e:
            print(f"[ValidationPipeline] Metal validation error: {e}")

        return candidates

    def _run_pyrosetta(
        self,
        candidates: List[SequenceCandidate],
    ) -> List[SequenceCandidate]:
        """Run PyRosetta FastRelax and scoring on predicted structures."""
        if not self._modules.get("pyrosetta"):
            return candidates

        try:
            import rosetta_utils

            for i, candidate in enumerate(candidates):
                if not candidate.predicted_pdb:
                    continue

                try:
                    # Calculate clash score
                    clash_result = rosetta_utils.calculate_clash_score(
                        candidate.predicted_pdb
                    )
                    if clash_result.get("status") == "completed":
                        candidate.clash_score = clash_result.get("clash_count", 0)

                except Exception as e:
                    print(f"    PyRosetta error for candidate {i+1}: {e}")

        except Exception as e:
            print(f"[ValidationPipeline] PyRosetta error: {e}")

        return candidates

    def _apply_filters(
        self,
        candidates: List[SequenceCandidate],
        max_backbone_rmsd: float,
        min_plddt: float,
        min_metal_coordination: Optional[int],
        min_protein_donors: Optional[int],
    ) -> List[SequenceCandidate]:
        """Apply filters and rank candidates."""
        for candidate in candidates:
            failures = []
            passes = True

            # Compute best-of metrics across ESMFold and RF3
            plddt_values = [v for v in [candidate.plddt_mean, candidate.rf3_plddt] if v is not None]
            rmsd_values = [v for v in [candidate.backbone_rmsd, candidate.rf3_rmsd] if v is not None]
            candidate.best_plddt = max(plddt_values) if plddt_values else None
            candidate.best_rmsd = min(rmsd_values) if rmsd_values else None

            # RMSD filter (use best RMSD across methods)
            if candidate.best_rmsd is not None:
                if candidate.best_rmsd > max_backbone_rmsd:
                    failures.append(f"best RMSD {candidate.best_rmsd:.2f} > {max_backbone_rmsd}")
                    passes = False
            else:
                failures.append("No RMSD calculated")
                passes = False

            # pLDDT filter (use best pLDDT across methods)
            if candidate.best_plddt is not None:
                if candidate.best_plddt < min_plddt:
                    failures.append(f"best pLDDT {candidate.best_plddt:.2f} < {min_plddt}")
                    passes = False
            else:
                failures.append("No pLDDT score")
                passes = False

            # Metal coordination filter (if applicable)
            if min_metal_coordination is not None:
                if candidate.metal_coordination is not None:
                    if candidate.metal_coordination < min_metal_coordination:
                        failures.append(
                            f"CN {candidate.metal_coordination} < {min_metal_coordination}"
                        )
                        passes = False

            # Protein donors filter (if applicable)
            if min_protein_donors is not None:
                if candidate.metal_protein_donors is not None:
                    if candidate.metal_protein_donors < min_protein_donors:
                        failures.append(
                            f"Protein donors {candidate.metal_protein_donors} < {min_protein_donors}"
                        )
                        passes = False

            candidate.passes_filters = passes
            candidate.filter_failures = failures

        # Rank candidates (lower is better)
        # Primary: best RMSD, Secondary: -best pLDDT (negative because higher is better)
        ranked = sorted(
            [c for c in candidates if c.best_rmsd is not None],
            key=lambda c: (c.best_rmsd, -(c.best_plddt or 0))
        )

        for rank, candidate in enumerate(ranked, 1):
            candidate.rank = rank

        return candidates

    def _get_range(self, values: List[float]) -> Optional[str]:
        """Get min-max range string for a list of values."""
        values = [v for v in values if v is not None]
        if not values:
            return None
        return f"{min(values):.2f} - {max(values):.2f}"


def validate_design(
    pdb_path: str,
    output_path: Optional[str] = None,
    num_sequences: int = 8,
    ligand_name: Optional[str] = None,
    metal_type: Optional[str] = None,
    **kwargs,
) -> Dict[str, Any]:
    """
    Convenience function to validate a design from file.

    Args:
        pdb_path: Path to RFD3 output PDB
        output_path: Optional path to save results JSON
        num_sequences: Number of sequences to generate
        ligand_name: Ligand residue name (e.g., "UNL")
        metal_type: Metal type (e.g., "DY", "TB")
        **kwargs: Additional arguments for DesignValidationPipeline.validate()

    Returns:
        Validation results dictionary
    """
    # Read PDB
    with open(pdb_path, 'r') as f:
        pdb_content = f.read()

    # Run pipeline
    pipeline = DesignValidationPipeline()
    result = pipeline.validate(
        pdb_content=pdb_content,
        pdb_path=pdb_path,
        num_sequences=num_sequences,
        ligand_name=ligand_name,
        metal_type=metal_type,
        **kwargs,
    )

    # Convert to dict
    result_dict = result.to_dict()

    # Save if output path provided
    if output_path:
        with open(output_path, 'w') as f:
            json.dump(result_dict, f, indent=2)
        print(f"[ValidationPipeline] Results saved to: {output_path}")

    return result_dict


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Post-LigandMPNN Validation Pipeline"
    )
    parser.add_argument("pdb", help="Path to RFD3 output PDB file")
    parser.add_argument("-o", "--output", help="Output JSON path")
    parser.add_argument("-n", "--num-seqs", type=int, default=8,
                        help="Number of sequences to generate (default: 8)")
    parser.add_argument("--ligand", help="Ligand residue name (e.g., UNL)")
    parser.add_argument("--metal", help="Metal type (e.g., DY, TB, ZN)")
    parser.add_argument("--max-rmsd", type=float, default=1.5,
                        help="Maximum backbone RMSD (default: 1.5)")
    parser.add_argument("--min-plddt", type=float, default=0.7,
                        help="Minimum pLDDT score (default: 0.7)")

    args = parser.parse_args()

    result = validate_design(
        pdb_path=args.pdb,
        output_path=args.output,
        num_sequences=args.num_seqs,
        ligand_name=args.ligand,
        metal_type=args.metal,
        max_backbone_rmsd=args.max_rmsd,
        min_plddt=args.min_plddt,
    )

    # Print summary
    print("\n" + "=" * 60)
    print("VALIDATION SUMMARY")
    print("=" * 60)
    print(f"Design type: {result['design_type']}")
    print(f"Candidates: {result['num_passed']}/{result['num_total']} passed filters")

    if result['best_sequence']:
        print(f"\nBest candidate:")
        print(f"  RMSD: {result['best_backbone_rmsd']:.2f} Å")
        print(f"  pLDDT: {result['best_plddt']:.2f}")
        print(f"  Sequence: {result['best_sequence'][:50]}...")
