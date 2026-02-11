/**
 * Shared step factory functions reused across pipeline definitions.
 * Each factory returns a PipelineStepDefinition with a real execute function
 * that calls the backend API.
 */

import { Dna, Shield, BarChart3, Search, Database, History, Lightbulb } from 'lucide-react';
import api from '@/lib/api';
import { extractPdbContent } from '@/lib/workflowHandlers';
import { getMetalMpnnConfig, getCofactorMpnnConfig } from './metal-mpnn-bias';
import type { MetalMpnnConfig } from './metal-mpnn-bias';
import type {
  PipelineStepDefinition,
  StepExecutionContext,
  StepResult,
  StepParameterField,
  PdbOutput,
  SequenceOutput,
} from '@/lib/pipeline-types';
import { ScaffoldSearchResultPreview } from '@/components/pipeline/ScaffoldSearchResultPreview';
import { ScoutResultPreview } from '@/components/pipeline/ScoutResultPreview';
import { LessonResultPreview } from '@/components/pipeline/LessonResultPreview';

// ---- Helpers ----

// FIX #12: Check abort signal and throw if cancelled
function checkAborted(signal: AbortSignal) {
  if (signal.aborted) {
    throw new DOMException('Step cancelled', 'AbortError');
  }
}

/** Submit an RFD3 job and wait for result, reporting progress. */
async function submitRfd3AndWait(
  ctx: StepExecutionContext,
  request: Record<string, unknown>,
): Promise<{ result: Record<string, unknown>; pdbContent?: string }> {
  ctx.onProgress(5, 'Submitting RFD3 design job...');
  checkAborted(ctx.abortSignal);

  const response = await api.submitRFD3Design(request as any);
  checkAborted(ctx.abortSignal);

  const jobId = response.job_id;
  ctx.onJobCreated?.(jobId, 'rfd3');

  ctx.onProgress(10, 'Job submitted, waiting for result...');

  let jobResult: { status: string; result?: Record<string, unknown>; error?: string };

  if (response.syncCompleted) {
    jobResult = {
      status: response.status,
      result: response.result as Record<string, unknown>,
      error: response.error,
    };
    ctx.onProgress(90, 'Design completed');
  } else {
    // FIX #12: Interpolate progress during polling
    let pollCount = 0;
    const polled = await api.waitForJob(jobId, (status) => {
      checkAborted(ctx.abortSignal);
      pollCount++;
      // Smooth progress: 10% → 85% over polling period
      const interpolated = Math.min(85, 10 + pollCount * 3);
      if (status.status === 'running') {
        ctx.onProgress(interpolated, 'Generating backbone structures...');
      }
    });
    jobResult = {
      status: polled.status,
      result: polled.result as Record<string, unknown>,
      error: polled.error,
    };
  }

  checkAborted(ctx.abortSignal);

  if (jobResult.status !== 'completed' || !jobResult.result) {
    throw new Error(jobResult.error || 'RFD3 design failed');
  }

  ctx.onProgress(95, 'Processing results...');

  const pdbContent = extractPdbContent(jobResult.result);
  return { result: jobResult.result, pdbContent };
}

/** Parse RFD3 designs array into PdbOutput items. */
function parseDesignsToPdbOutputs(result: Record<string, unknown>): PdbOutput[] {
  const designs = result.designs as Array<Record<string, unknown>> | undefined;
  if (!Array.isArray(designs)) {
    // Fallback: single PDB
    const pdb = extractPdbContent(result);
    if (pdb) {
      return [{ id: 'design-1', label: 'Design 1', pdbContent: pdb }];
    }
    return [];
  }

  return designs.map((d, idx) => {
    const pdb = (d.content as string) || (d.pdb_content as string) || '';
    const metrics: Record<string, number | string> = {};

    if (d.seed !== undefined) metrics.seed = d.seed as number;
    if (d.filename) metrics.file = d.filename as string;

    return {
      id: `design-${idx + 1}`,
      label: (d.filename as string)?.replace('.pdb', '') || `Design ${idx + 1}`,
      pdbContent: pdb,
      metrics,
    };
  });
}

// ---- Shared Step Factories ----

/**
 * Create an RFD3 backbone generation step.
 */
export function createRfd3Step(overrides?: Partial<PipelineStepDefinition>): PipelineStepDefinition {
  return {
    id: 'rfd3_backbone',
    name: 'RFD3 Backbone Design',
    description: 'Generate protein backbone structures using RFdiffusion3',
    icon: Shield,
    requiresReview: true,
    supportsSelection: true,
    optional: false,
    defaultParams: {
      num_designs: 4,
      num_timesteps: 200,
      step_scale: 1.5,
      gamma_0: 0.6,
    },
    parameterSchema: [
      { id: 'num_designs', label: 'Number of Designs', type: 'slider', required: false, defaultValue: 4, range: { min: 1, max: 10, step: 1 } },
      { id: 'num_timesteps', label: 'Timesteps', type: 'slider', required: false, defaultValue: 200, range: { min: 50, max: 500, step: 50 } },
      { id: 'step_scale', label: 'Step Scale', type: 'slider', required: false, defaultValue: 1.5, range: { min: 0.5, max: 3.0, step: 0.1 }, helpText: 'Higher = more designable, less diverse' },
      { id: 'gamma_0', label: 'Gamma', type: 'slider', required: false, defaultValue: 0.6, range: { min: 0.1, max: 1.0, step: 0.1 }, helpText: 'Lower = more designable' },
    ],

    async execute(ctx: StepExecutionContext): Promise<StepResult> {
      const { initialParams, params, previousResults } = ctx;

      // Extract rfd3_params from any analysis step's result
      let analysisRfd3Params: Record<string, unknown> = {};
      for (const result of Object.values(previousResults)) {
        const rfd3Data = (result.data as Record<string, unknown> | undefined)?.rfd3_params;
        if (rfd3Data && typeof rfd3Data === 'object') {
          analysisRfd3Params = { ...analysisRfd3Params, ...rfd3Data as Record<string, unknown> };
        }
      }

      const request: Record<string, unknown> = {
        ...initialParams,
        ...analysisRfd3Params,
        num_designs: params.num_designs ?? 4,
        num_timesteps: params.num_timesteps ?? analysisRfd3Params.num_timesteps ?? 200,
        step_scale: params.step_scale ?? 1.5,
        gamma_0: params.gamma_0 ?? 0.6,
      };

      // Ensure at least contig or length is present
      if (!request.contig && !request.length) {
        throw new Error('Missing contig or length parameter. The analysis step may not have produced structure information.');
      }

      const { result } = await submitRfd3AndWait(ctx, request);
      const pdbOutputs = parseDesignsToPdbOutputs(result);

      ctx.onProgress(100, 'Backbone generation complete');

      return {
        id: `rfd3-${Date.now()}`,
        summary: `Generated ${pdbOutputs.length} backbone design(s)`,
        pdbOutputs,
        data: { num_designs: pdbOutputs.length },
      };
    },

    ...overrides,
  };
}

/**
 * Create an MPNN sequence design step.
 */
export function createMpnnStep(overrides?: Partial<PipelineStepDefinition>): PipelineStepDefinition {
  return {
    id: 'mpnn_sequence',
    name: 'LigandMPNN Sequence Design',
    description: 'Design amino acid sequences for backbone structures',
    icon: Dna,
    requiresReview: true,
    supportsSelection: true,
    optional: false,
    defaultParams: {
      temperature: 0.1,
      temperature_high: 0.2,
      num_sequences: 8,
      model_type: 'ligand_mpnn',
    },
    parameterSchema: [
      { id: 'temperature', label: 'Low Temperature', type: 'slider', required: false, defaultValue: 0.1, range: { min: 0.01, max: 1.0, step: 0.01 }, helpText: 'Conservative sampling (half of sequences)' },
      { id: 'temperature_high', label: 'High Temperature', type: 'slider', required: false, defaultValue: 0.2, range: { min: 0.01, max: 1.0, step: 0.01 }, helpText: 'Exploratory sampling (half of sequences)' },
      { id: 'num_sequences', label: 'Sequences per Design (total)', type: 'slider', required: false, defaultValue: 8, range: { min: 2, max: 16, step: 2 }, helpText: 'Split evenly between low and high temperature' },
      {
        id: 'model_type', label: 'Model', type: 'select', required: false, defaultValue: 'ligand_mpnn',
        options: [
          { value: 'ligand_mpnn', label: 'LigandMPNN' },
          { value: 'protein_mpnn', label: 'ProteinMPNN' },
        ],
      },
    ],

    async execute(ctx: StepExecutionContext): Promise<StepResult> {
      const { previousResults, selectedItems, params, initialParams } = ctx;

      // Sweep pass-through: if a prior step ran metal binding sweep (session_id),
      // sequences are already designed — return them directly without re-running MPNN.
      const sweepResult = Object.values(previousResults).find(
        r => r.data?.session_id && r.sequences && r.sequences.length > 0
      );
      if (sweepResult) {
        ctx.onProgress(100, 'Sequences already designed by metal binding sweep');
        return {
          id: `mpnn-passthrough-${Date.now()}`,
          summary: `${sweepResult.sequences!.length} sequences from metal binding sweep (already designed)`,
          sequences: sweepResult.sequences,
          data: {
            total_sequences: sweepResult.sequences!.length,
            passthrough: true,
            source: 'metal_binding_sweep',
          },
        };
      }

      // Find PDB outputs from previous steps, excluding scaffold search results
      const pdbOutputs = findPdbOutputsExcluding(
        previousResults,
        selectedItems,
        (_stepId, result) => !!(result.data?.searched),
      );
      if (pdbOutputs.length === 0) {
        throw new Error('No backbone structures available for sequence design');
      }

      // Resolve metal-aware MPNN config.
      // Priority: previous step's mpnn_params > derive from target_metal > empty.
      let metalConfig: MetalMpnnConfig = {};
      for (const result of Object.values(previousResults)) {
        const mp = (result.data as Record<string, unknown> | undefined)?.mpnn_params;
        if (mp && typeof mp === 'object') {
          metalConfig = { ...metalConfig, ...mp as MetalMpnnConfig };
        }
      }
      if (!metalConfig.bias_AA && !metalConfig.omit_AA) {
        // Check initialParams first, then fall back to intent result from previous steps
        let targetMetal = (initialParams.target_metal as string) || (initialParams.targetMetal as string);
        let ligandName = '';
        if (!targetMetal) {
          const intentResult = Object.values(previousResults).find(r => r.data?.design_type);
          targetMetal = (intentResult?.data?.target_metal as string) || '';
          ligandName = (intentResult?.data?.ligand_name as string) || '';
        }
        if (!ligandName) {
          ligandName = (initialParams.ligand_name as string) || '';
        }
        if (targetMetal) {
          // Look for backend-provided cofactor bias from analyze_ligand_features
          let cofactorBias: Record<string, number> | undefined;
          for (const result of Object.values(previousResults)) {
            const biasData = (result.data as Record<string, unknown> | undefined)?.mpnn_bias_adjustments;
            if (biasData && typeof biasData === 'object') {
              cofactorBias = biasData as Record<string, number>;
              break;
            }
          }
          metalConfig = getCofactorMpnnConfig(targetMetal, ligandName.toUpperCase() || undefined, cofactorBias);
        }
      }

      const allSequences: SequenceOutput[] = [];
      const totalDesigns = pdbOutputs.length;

      for (let i = 0; i < totalDesigns; i++) {
        checkAborted(ctx.abortSignal);

        const pdb = pdbOutputs[i];
        const progressBase = (i / totalDesigns) * 90;
        ctx.onProgress(progressBase + 5, `Designing sequences for ${pdb.label} (${i + 1}/${totalDesigns})...`);

        // For metal-containing backbones, fix residues that coordinate the metal
        // so MPNN doesn't redesign the binding site.
        let perBackboneFixed: string[] | undefined;
        if (metalConfig.omit_AA || metalConfig.bias_AA) {
          const coordResidues = findMetalCoordinatingResidues(pdb.pdbContent);
          if (coordResidues) {
            // Format: ["A1", "A2", ...] — chain A + residue numbers
            perBackboneFixed = coordResidues.split(',').map(r => `A${r}`);
            console.log(`[MPNN] Design ${i + 1}: fixing ${perBackboneFixed.length} metal-coordinating residues: ${coordResidues}`);
          }
        }

        // Multi-temperature MPNN: split sequences between low (0.1) and high (0.2)
        // temperatures for better sequence diversity. RFD3 paper uses 8 seqs/backbone
        // for enzyme/small molecule designs (bioRxiv 2025.09.18.676967v2).
        const totalSeqsPerBackbone = (params.num_sequences as number) ?? 8;
        const tempLow = (params.temperature as number) ?? 0.1;
        const tempHigh = (params.temperature_high as number) ?? 0.2;
        const seqsLow = Math.ceil(totalSeqsPerBackbone / 2);
        const seqsHigh = totalSeqsPerBackbone - seqsLow;
        const temperatures = tempLow === tempHigh
          ? [{ temp: tempLow, count: totalSeqsPerBackbone }]
          : [{ temp: tempLow, count: seqsLow }, { temp: tempHigh, count: seqsHigh }];

        const baseRequest = {
          pdb_content: pdb.pdbContent,
          model_type: (params.model_type as 'ligand_mpnn' | 'protein_mpnn') ?? 'ligand_mpnn',
          ...metalConfig,
          ...(perBackboneFixed ? {
            fixed_positions: [
              ...(metalConfig.fixed_positions || []),
              ...perBackboneFixed,
            ],
          } : {}),
        };

        // Run MPNN at each temperature and merge sequences
        const tempResults: Array<{ status: string; result?: Record<string, unknown>; error?: string }> = [];
        for (const { temp, count } of temperatures) {
          const response = await api.submitMPNNDesign({
            ...baseRequest,
            temperature: temp,
            num_sequences: count,
          });

          checkAborted(ctx.abortSignal);
          ctx.onJobCreated?.(response.job_id, 'mpnn');

          let jobResult: { status: string; result?: Record<string, unknown>; error?: string };
          if (response.status === 'completed') {
            jobResult = { status: response.status, result: response.result as Record<string, unknown>, error: response.error };
          } else {
            const polled = await api.waitForJob(response.job_id, () => { checkAborted(ctx.abortSignal); });
            jobResult = { status: polled.status, result: polled.result as Record<string, unknown>, error: polled.error };
          }
          if (jobResult.status === 'completed' && jobResult.result) {
            tempResults.push(jobResult);
          }
        }

        // Merge results from multiple temperature runs
        let jobResult: { status: string; result?: Record<string, unknown>; error?: string };
        if (tempResults.length === 0) {
          continue; // All temperature runs failed
        } else if (tempResults.length === 1) {
          jobResult = tempResults[0];
        } else {
          // Merge sequences from both temperature runs
          const mergedResult = { ...tempResults[0].result! };
          const sequences1 = (mergedResult.sequences || []) as Array<Record<string, unknown>>;
          const sequences2 = ((tempResults[1].result?.sequences || []) as Array<Record<string, unknown>>);
          mergedResult.sequences = [...sequences1, ...sequences2];
          jobResult = { status: 'completed', result: mergedResult };
        }

        checkAborted(ctx.abortSignal);

        if (jobResult.status !== 'completed' || !jobResult.result) {
          throw new Error(jobResult.error || `MPNN design failed for ${pdb.label}`);
        }

        // Parse MPNN sequences from result
        const sequences = parseMpnnSequences(jobResult.result, pdb.id, i);
        allSequences.push(...sequences);

        ctx.onProgress(progressBase + (90 / totalDesigns), `${pdb.label}: ${sequences.length} sequences`);
      }

      ctx.onProgress(100, 'Sequence design complete');

      return {
        id: `mpnn-${Date.now()}`,
        summary: `Designed ${allSequences.length} sequences for ${totalDesigns} backbone(s)`,
        sequences: allSequences,
        data: { total_sequences: allSequences.length, backbones_processed: totalDesigns },
      };
    },

    ...overrides,
  };
}

/**
 * Create an RF3 structure validation step.
 */
export function createRf3Step(overrides?: Partial<PipelineStepDefinition>): PipelineStepDefinition {
  return {
    id: 'rf3_validation',
    name: 'RF3 Structure Validation',
    description: 'Predict 3D structures from sequences to validate designability',
    icon: Shield,
    requiresReview: true,
    supportsSelection: true,
    optional: false,
    defaultParams: {},
    parameterSchema: [],

    async execute(ctx: StepExecutionContext): Promise<StepResult> {
      const { previousResults, selectedItems } = ctx;

      // Sweep pass-through: if a prior step ran metal binding sweep (session_id),
      // structures are already validated — return them directly without re-running RF3.
      const sweepResult = Object.values(previousResults).find(
        r => r.data?.session_id && (r.pdbOutputs?.length ?? 0) > 0
      );
      if (sweepResult && sweepResult.pdbOutputs) {
        ctx.onProgress(100, 'Structures already validated by metal binding sweep');
        return {
          id: `rf3-passthrough-${Date.now()}`,
          summary: `${sweepResult.pdbOutputs.length} structures validated by sweep`,
          pdbOutputs: sweepResult.pdbOutputs,
          data: {
            total_validated: sweepResult.pdbOutputs.length,
            total_attempted: sweepResult.pdbOutputs.length,
            success_rate: '100.0%',
            passthrough: true,
            source: 'metal_binding_sweep',
          },
        };
      }
      // Also handle sweep with session_id but only sequences (no PDB content embedded)
      const sweepSeqOnly = Object.values(previousResults).find(
        r => r.data?.session_id && r.sequences && r.sequences.length > 0 && !(r.pdbOutputs?.length)
      );
      if (sweepSeqOnly) {
        ctx.onProgress(100, 'Sweep completed validation internally — skipping redundant RF3');
        return {
          id: `rf3-passthrough-${Date.now()}`,
          summary: `${sweepSeqOnly.sequences!.length} sequences validated by sweep (no PDB re-prediction needed)`,
          sequences: sweepSeqOnly.sequences,
          data: {
            total_validated: sweepSeqOnly.sequences!.length,
            total_attempted: sweepSeqOnly.sequences!.length,
            success_rate: '100.0%',
            passthrough: true,
            source: 'metal_binding_sweep',
          },
        };
      }

      // Find sequences from previous MPNN step
      const sequences = findSequenceOutputs(previousResults, selectedItems);
      if (sequences.length === 0) {
        throw new Error('No sequences available for validation');
      }

      // Extract ligand SMILES and metal type from previous results for RF3 prediction
      let ligandSmiles: string | undefined;
      let targetMetal: string | undefined;
      for (const result of Object.values(previousResults)) {
        if (!ligandSmiles) {
          const smiles = result.data?.ligand_smiles as string;
          if (smiles) ligandSmiles = smiles;
        }
        if (!targetMetal) {
          const metal = result.data?.target_metal as string;
          if (metal) targetMetal = metal;
        }
      }

      const pdbOutputs: PdbOutput[] = [];
      const totalSeqs = sequences.length;

      for (let i = 0; i < totalSeqs; i++) {
        checkAborted(ctx.abortSignal);

        const seq = sequences[i];
        const progressBase = (i / totalSeqs) * 90;
        ctx.onProgress(progressBase + 5, `Validating sequence ${i + 1}/${totalSeqs}...`);

        const response = await api.submitRF3Prediction({
          sequence: seq.sequence,
          name: seq.label ?? seq.id,
          ligand_smiles: ligandSmiles,
          metal: targetMetal,
        });

        checkAborted(ctx.abortSignal);
        ctx.onJobCreated?.(response.job_id, 'rf3');

        let jobResult: { status: string; result?: Record<string, unknown>; error?: string };

        if ((response as any).syncCompleted) {
          jobResult = { status: response.status, result: response.result as Record<string, unknown>, error: response.error };
        } else {
          const polled = await api.waitForJob(response.job_id, () => { checkAborted(ctx.abortSignal); });
          jobResult = { status: polled.status, result: polled.result as Record<string, unknown>, error: polled.error };
        }

        checkAborted(ctx.abortSignal);

        if (jobResult.status !== 'completed' || !jobResult.result) {
          // Don't fail entirely for individual failures
          continue;
        }

        const pdbContent = extractRf3PdbContent(jobResult.result);
        const confidences = extractConfidences(jobResult.result);

        if (pdbContent) {
          const metrics: Record<string, number | string> = {};
          if (confidences?.overall_plddt) metrics.pLDDT = confidences.overall_plddt;
          if (confidences?.ptm) metrics.pTM = confidences.ptm;
          if (confidences?.overall_pae) metrics.PAE = confidences.overall_pae;
          // AF3 paper ligand-interface metrics (bioRxiv 2025.09.18.676967v2)
          if (confidences?.iptm != null) metrics.iPTM = confidences.iptm;
          if (confidences?.min_chain_pair_pae != null) metrics.min_chain_pair_pae = confidences.min_chain_pair_pae;

          pdbOutputs.push({
            id: `rf3-${seq.id}`,
            label: seq.label ?? `Validation ${i + 1}`,
            pdbContent,
            metrics,
            sequence: seq.sequence,
          });
        }

        ctx.onProgress(progressBase + (90 / totalSeqs), `Validated ${i + 1}/${totalSeqs}`);
      }

      ctx.onProgress(100, 'Validation complete');

      // Compute average quality metrics for summary
      const ptmValues = pdbOutputs.map(p => p.metrics?.pTM as number).filter(v => typeof v === 'number');
      const plddtValues = pdbOutputs.map(p => p.metrics?.pLDDT as number).filter(v => typeof v === 'number');
      const avgPtm = ptmValues.length > 0 ? ptmValues.reduce((a, b) => a + b, 0) / ptmValues.length : undefined;
      const avgPlddt = plddtValues.length > 0 ? plddtValues.reduce((a, b) => a + b, 0) / plddtValues.length : undefined;

      const qualitySuffix = avgPtm !== undefined
        ? ` (avg pTM: ${avgPtm.toFixed(3)}, pLDDT: ${(avgPlddt ?? 0).toFixed(3)})`
        : '';

      return {
        id: `rf3-${Date.now()}`,
        summary: `Validated ${pdbOutputs.length}/${totalSeqs} sequences${qualitySuffix}`,
        pdbOutputs,
        data: {
          total_validated: pdbOutputs.length,
          total_attempted: totalSeqs,
          success_rate: totalSeqs > 0 ? ((pdbOutputs.length / totalSeqs) * 100).toFixed(1) + '%' : 'N/A',
          ...(avgPtm !== undefined ? { avg_pTM: avgPtm } : {}),
          ...(avgPlddt !== undefined ? { avg_pLDDT: avgPlddt } : {}),
        },
      };
    },

    ...overrides,
  };
}

/**
 * Create a design evaluation step (metal binding specific).
 */
export function createEvaluateStep(overrides?: Partial<PipelineStepDefinition>): PipelineStepDefinition {
  return {
    id: 'evaluate',
    name: 'Design Evaluation',
    description: 'Evaluate and score the final designs',
    icon: BarChart3,
    requiresReview: true,
    supportsSelection: false,
    optional: false,
    defaultParams: {},
    parameterSchema: [],

    async execute(ctx: StepExecutionContext): Promise<StepResult> {
      const { previousResults, selectedItems, initialParams } = ctx;

      // Find PDB structures from previous step (RF3 validation outputs)
      const pdbOutputs = findPdbOutputs(previousResults, selectedItems);
      if (pdbOutputs.length === 0) {
        throw new Error('No structures available for evaluation');
      }

      const evaluations: Array<{ pdb: PdbOutput; eval: Record<string, unknown> }> = [];

      for (let i = 0; i < pdbOutputs.length; i++) {
        checkAborted(ctx.abortSignal);

        const pdb = pdbOutputs[i];
        ctx.onProgress((i / pdbOutputs.length) * 90 + 5, `Evaluating ${pdb.label} (${i + 1}/${pdbOutputs.length})...`);

        try {
          const evaluation = await api.evaluateDesign({
            pdb_content: pdb.pdbContent,
            target_metal: (initialParams.target_metal as string) || (initialParams.targetMetal as string) || 'ZN',
            metal_chain: initialParams.metal_chain as string,
            metal_resnum: initialParams.metal_resnum as number,
          });

          evaluations.push({ pdb, eval: evaluation as unknown as Record<string, unknown> });
        } catch {
          // Continue evaluating other designs
        }
      }

      ctx.onProgress(100, 'Evaluation complete');

      // Build enriched PDB outputs with evaluation metrics
      const evaluatedPdbs: PdbOutput[] = evaluations.map(({ pdb, eval: ev }) => ({
        ...pdb,
        metrics: {
          ...pdb.metrics,
          coordination: ev.coordination_number as number,
          geometry: ev.geometry_type as string,
          bond_dist: ev.avg_bond_distance as number,
          pass: ev.overall_pass ? 'Yes' : 'No',
        },
      }));

      const passCount = evaluations.filter(e => e.eval.overall_pass).length;

      return {
        id: `eval-${Date.now()}`,
        summary: `${passCount}/${evaluations.length} designs passed evaluation criteria`,
        pdbOutputs: evaluatedPdbs,
        data: {
          total_evaluated: evaluations.length,
          passed: passCount,
          failed: evaluations.length - passCount,
        },
      };
    },

    ...overrides,
  };
}

/**
 * Create a structure analysis step.
 */
export function createAnalyzeStep(overrides?: Partial<PipelineStepDefinition>): PipelineStepDefinition {
  return {
    id: 'analyze',
    name: 'Structure Analysis',
    description: 'Analyze the input structure for binding sites and design opportunities',
    icon: Search,
    requiresReview: true,
    supportsSelection: false,
    optional: false,
    defaultParams: {},
    parameterSchema: [],

    async execute(ctx: StepExecutionContext): Promise<StepResult> {
      const { initialParams } = ctx;
      const pdbContent = initialParams.pdb_content as string;

      if (!pdbContent) {
        throw new Error('No PDB content provided for analysis');
      }

      ctx.onProgress(10, 'Analyzing structure...');

      const analysis = await api.analyzeStructure(pdbContent);

      ctx.onProgress(100, 'Analysis complete');

      return {
        id: `analysis-${Date.now()}`,
        summary: `Found ${analysis.binding_sites?.length ?? 0} binding site(s) and ${analysis.ligands?.length ?? 0} ligand(s)`,
        data: {
          num_residues: analysis.num_residues,
          num_ligands: analysis.ligands?.length ?? 0,
          num_binding_sites: analysis.binding_sites?.length ?? 0,
          num_suggestions: analysis.suggestions?.length ?? 0,
          ligands: analysis.ligands?.map(l => l.name).join(', ') || 'None',
        },
      };
    },

    ...overrides,
  };
}

/**
 * Create a scaffold search step that queries RCSB PDB for existing structures
 * containing a target metal-ligand complex. If a good scaffold is found,
 * its PDB content is passed forward for downstream steps to use.
 *
 * Reusable across any pipeline that designs around metal-ligand binding sites.
 */
export function createScaffoldSearchStep(overrides?: Partial<PipelineStepDefinition>): PipelineStepDefinition {
  return {
    id: 'scaffold_search',
    name: 'Scaffold Search',
    description: 'Search RCSB PDB for existing structures with the target metal-ligand complex',
    icon: Database,
    requiresReview: true,
    supportsSelection: true,
    optional: true,
    defaultParams: {
      resolution_max: 3.0,
      limit: 10,
    },
    parameterSchema: [
      { id: 'resolution_max', label: 'Max Resolution (\u00C5)', type: 'slider', required: false, defaultValue: 3.0, range: { min: 1.0, max: 5.0, step: 0.5 }, helpText: 'Maximum crystallographic resolution for PDB search' },
      { id: 'limit', label: 'Max Candidates', type: 'slider', required: false, defaultValue: 10, range: { min: 1, max: 25, step: 1 }, helpText: 'Maximum number of PDB hits to validate' },
    ],

    async execute(ctx: StepExecutionContext): Promise<StepResult> {
      const { previousResults, params, initialParams } = ctx;

      // Extract metal and ligand from previous parse step or initial params
      const intentResult = Object.values(previousResults).find(r => r.data?.design_type);
      const targetMetal = (intentResult?.data?.target_metal as string) || (initialParams.target_metal as string) || '';
      const ligandName = (intentResult?.data?.ligand_name as string) || (initialParams.ligand_name as string) || '';
      const ligandSmiles = (intentResult?.data?.ligand_smiles as string) || (initialParams.ligand_smiles as string) || '';

      // Guard: skip if no metal+ligand or already have a PDB
      if (!targetMetal || !ligandName) {
        return {
          id: `scaffold-skip-${Date.now()}`,
          summary: 'Skipped: metal and ligand required for scaffold search',
          data: { skipped: true, reason: 'No metal or ligand specified' },
        };
      }

      const structureResult = Object.values(previousResults).find(r => r.pdbOutputs && r.pdbOutputs.length > 0);
      if (structureResult) {
        return {
          id: `scaffold-skip-${Date.now()}`,
          summary: 'Skipped: input structure already provided',
          data: { skipped: true, reason: 'PDB structure already available from previous step' },
        };
      }

      ctx.onProgress(10, `Searching RCSB for ${targetMetal} + ${ligandName} structures...`);

      checkAborted(ctx.abortSignal);

      try {
        const searchResult = await api.searchScaffold({
          metal: targetMetal,
          ligand_name: ligandName,
          ligand_smiles: ligandSmiles || undefined,
          resolution_max: (params.resolution_max as number) ?? 3.0,
          limit: (params.limit as number) ?? 10,
          fetch_pdb: true,
        });

        checkAborted(ctx.abortSignal);

        ctx.onProgress(80, 'Processing search results...');

        const pdbOutputs: PdbOutput[] = [];
        const candidatePdbs = searchResult.candidate_pdbs || {};

        // Create PDB outputs for ALL candidates with fetched content
        if (searchResult.candidates) {
          for (const candidate of searchResult.candidates) {
            const pdbContent = candidatePdbs[candidate.pdb_id];
            if (pdbContent) {
              pdbOutputs.push({
                id: `scaffold-${candidate.pdb_id}`,
                label: `${candidate.pdb_id} (score: ${candidate.total_score})`,
                pdbContent,
                metrics: {
                  score: candidate.total_score,
                  coordination: candidate.coordination_number,
                  source_metal: candidate.source_metal,
                  needs_substitution: candidate.needs_substitution ? 'Yes' : 'No',
                },
              });
            }
          }
        }
        // Fallback: use best_pdb_content if candidate_pdbs wasn't populated
        if (pdbOutputs.length === 0 && searchResult.best_pdb_content && searchResult.best_candidate) {
          const best = searchResult.best_candidate;
          pdbOutputs.push({
            id: `scaffold-${best.pdb_id}`,
            label: `${best.pdb_id} (score: ${best.total_score})`,
            pdbContent: searchResult.best_pdb_content,
            metrics: {
              score: best.total_score,
              coordination: best.coordination_number,
              source_metal: best.source_metal,
              needs_substitution: best.needs_substitution ? 'Yes' : 'No',
            },
          });
        }

        ctx.onProgress(100, searchResult.recommended_action === 'scaffold'
          ? `Found scaffold: ${searchResult.best_candidate?.pdb_id}`
          : 'No suitable scaffold found — using de novo');

        return {
          id: `scaffold-${Date.now()}`,
          summary: searchResult.recommended_action === 'scaffold'
            ? `Found ${searchResult.best_candidate!.pdb_id} (score: ${searchResult.best_candidate!.total_score}, ${searchResult.num_pdb_hits} hits, ${searchResult.num_validated} validated)`
            : `No scaffold found (${searchResult.num_pdb_hits} hits, ${searchResult.num_validated} validated). ${searchResult.reason}`,
          pdbOutputs,
          data: {
            searched: true,
            recommended_action: searchResult.recommended_action,
            query_metal: searchResult.query_metal,
            query_ligand: searchResult.query_ligand,
            ligand_code: searchResult.ligand_code,
            num_pdb_hits: searchResult.num_pdb_hits,
            num_validated: searchResult.num_validated,
            candidates: searchResult.candidates,
            best_candidate: searchResult.best_candidate,
            reason: searchResult.reason,
            // Pass scaffold PDB ID forward for downstream steps
            scaffold_pdb_id: searchResult.best_candidate?.pdb_id,
          },
        };
      } catch (err) {
        // Non-fatal: scaffold search failure should not block the pipeline
        const message = err instanceof Error ? err.message : String(err);
        ctx.onProgress(100, `Scaffold search failed (non-fatal): ${message}`);

        return {
          id: `scaffold-fail-${Date.now()}`,
          summary: `Scaffold search failed: ${message}. Continuing with de novo.`,
          data: { searched: false, search_error: message, recommended_action: 'de_novo' },
        };
      }
    },

    ResultPreview: ScaffoldSearchResultPreview,

    ...overrides,
  };
}

/**
 * Create a scout filter step that validates 1 sequence per backbone
 * and removes backbones that fail pTM/pLDDT thresholds.
 *
 * Applies only to the general RFD3 path. Skips automatically when:
 * - Previous step was metal binding sweep (already validated)
 * - No pdbOutputs from previous step
 * - Only 1 backbone (not worth filtering)
 */
export function createScoutFilterStep(overrides?: Partial<PipelineStepDefinition>): PipelineStepDefinition {
  return {
    id: 'scout_filter',
    name: 'Scout Filter',
    description: 'Validate 1 sequence per backbone and filter out weak backbones',
    icon: Shield,
    requiresReview: true,
    supportsSelection: false,
    optional: true,
    defaultParams: {
      ptm_threshold: 0.6,
      plddt_threshold: 0.65,
    },
    parameterSchema: [
      { id: 'ptm_threshold', label: 'pTM Threshold', type: 'slider', required: false, defaultValue: 0.6, range: { min: 0.3, max: 0.9, step: 0.05 }, helpText: 'Minimum pTM to pass scout validation' },
      { id: 'plddt_threshold', label: 'pLDDT Threshold', type: 'slider', required: false, defaultValue: 0.65, range: { min: 0.3, max: 0.9, step: 0.05 }, helpText: 'Minimum pLDDT to pass scout validation' },
    ],

    ResultPreview: ScoutResultPreview,

    async execute(ctx: StepExecutionContext): Promise<StepResult> {
      const { previousResults, selectedItems, params, initialParams } = ctx;

      // Guard: skip if metal sweep (already validated internally)
      const rfd3Result = previousResults['rfd3_nl'];
      if (rfd3Result?.data?.session_id) {
        return {
          id: `scout-skip-${Date.now()}`,
          summary: 'Skipped: metal sweep already validated',
          data: { skipped: true, reason: 'Metal binding sweep performs internal validation' },
        };
      }

      // Metal single mode: backbones contain HETATM (metal/ligand) which scout's
      // simple RF3 call can't handle. But the CPU pre-filter CAN check these backbones
      // for broken chains, missing metal/ligand, etc. Run pre-filter only (no GPU).
      // NOTE: Don't re-emit pdbOutputs here — MPNN finds them from rfd3_nl directly.
      // Re-emitting would cause duplicates (MPNN collects from ALL previous steps).
      if (rfd3Result?.data?.metal_single_mode) {
        const metalPdbOutputs = findPdbOutputs(previousResults, selectedItems);
        if (metalPdbOutputs.length <= 1) {
          return {
            id: `scout-skip-${Date.now()}`,
            summary: 'Skipped: single metal backbone',
            data: { skipped: true, reason: 'Only 1 backbone — pre-filter not needed' },
          };
        }

        ctx.onProgress(5, `Pre-filtering ${metalPdbOutputs.length} metal backbones (CPU only)...`);
        checkAborted(ctx.abortSignal);

        const intentResult = Object.values(previousResults).find(r => r.data?.design_type);
        const scoutTargetMetal = (intentResult?.data?.target_metal as string)
          || (rfd3Result?.data?.target_metal as string)
          || (initialParams.target_metal as string)
          || '';
        const scoutLigandName = (intentResult?.data?.ligand_name as string)
          || (rfd3Result?.data?.ligand_name as string)
          || (initialParams.ligand_name as string)
          || '';

        // Call scout filter with pre-filter only flag — backend runs CPU checks, skips MPNN+RF3
        const pfResult = await api.scoutFilter({
          backbone_pdbs: metalPdbOutputs.map(p => p.pdbContent),
          ptm_threshold: (params.ptm_threshold as number) ?? 0.6,
          plddt_threshold: (params.plddt_threshold as number) ?? 0.65,
          target_metal: scoutTargetMetal,
          ligand_name: scoutLigandName || undefined,
          pre_filter_only: true,
        });

        checkAborted(ctx.abortSignal);

        const pfPassed = pfResult.pre_filter_passed ?? metalPdbOutputs.length;
        const pfFailed = pfResult.pre_filter_failed ?? 0;

        ctx.onProgress(100, pfFailed > 0
          ? `Pre-filter: ${pfFailed} backbone(s) eliminated`
          : 'Pre-filter: all backbones passed');

        return {
          id: `scout-prefilter-${Date.now()}`,
          summary: pfFailed > 0
            ? `Pre-filter: ${pfPassed}/${metalPdbOutputs.length} metal backbones passed (${pfFailed} eliminated)`
            : `Pre-filter: all ${metalPdbOutputs.length} metal backbones passed`,
          // Don't re-emit pdbOutputs — MPNN step finds them from rfd3_nl
          data: {
            metal_pre_filter_only: true,
            original_count: metalPdbOutputs.length,
            passing_count: pfPassed,
            scout_results: pfResult.scout_results,
            pre_filter_results: pfResult.pre_filter_results,
            pre_filter_passed: pfPassed,
            pre_filter_failed: pfFailed,
            ptm_threshold: pfResult.ptm_threshold,
            plddt_threshold: pfResult.plddt_threshold,
          },
        };
      }

      // Gather PDB outputs from previous steps
      const pdbOutputs = findPdbOutputs(previousResults, selectedItems);

      // Guard: skip if no backbones
      if (pdbOutputs.length === 0) {
        return {
          id: `scout-skip-${Date.now()}`,
          summary: 'Skipped: no backbone structures',
          data: { skipped: true, reason: 'No backbone PDB outputs from previous step' },
        };
      }

      // Guard: skip if only 1 backbone (not worth filtering)
      // NOTE: Don't re-emit pdbOutputs — MPNN finds them from rfd3_nl directly.
      if (pdbOutputs.length <= 1) {
        return {
          id: `scout-skip-${Date.now()}`,
          summary: 'Skipped: single backbone',
          data: { skipped: true, reason: 'Only 1 backbone — scout filtering not needed' },
        };
      }

      ctx.onProgress(5, `Scout filtering ${pdbOutputs.length} backbone(s): 1 seq each → MPNN → RF3...`);
      checkAborted(ctx.abortSignal);

      // Extract metal/ligand from intent result (NL pipeline) or initialParams (direct pipeline)
      const intentResult = Object.values(previousResults).find(r => r.data?.design_type);
      const scoutTargetMetal = (intentResult?.data?.target_metal as string)
        || (initialParams.target_metal as string)
        || '';
      const scoutLigandSmiles = (intentResult?.data?.ligand_smiles as string)
        || (initialParams.ligand_smiles as string)
        || '';
      const scoutLigandName = (intentResult?.data?.ligand_name as string)
        || (initialParams.ligand_name as string)
        || '';

      ctx.onProgress(10, `Validating ${pdbOutputs.length} backbone(s) (MPNN + RF3 per backbone)...`);

      const result = await api.scoutFilter({
        backbone_pdbs: pdbOutputs.map(p => p.pdbContent),
        ptm_threshold: (params.ptm_threshold as number) ?? 0.6,
        plddt_threshold: (params.plddt_threshold as number) ?? 0.65,
        target_metal: scoutTargetMetal,
        ligand_smiles: scoutLigandSmiles,
        ligand_name: scoutLigandName || undefined,
      });

      checkAborted(ctx.abortSignal);

      // Debug: log raw scout results to help diagnose display issues
      console.log('[Scout] Raw API result:', JSON.stringify({
        passing_count: result.passing_count,
        original_count: result.original_count,
        scout_results: result.scout_results?.map((sr: any) => ({
          idx: sr.backbone_index, passed: sr.passed, ptm: sr.ptm, plddt: sr.plddt,
          pre_filter_failed: sr.pre_filter_failed,
        })),
      }));

      ctx.onProgress(90, `Processing results: ${result.passing_count}/${result.original_count} passed...`);

      // Build filtered PDB outputs — only passing backbones
      const filteredOutputs: PdbOutput[] = [];
      for (const sr of result.scout_results) {
        if (sr.passed) {
          const originalPdb = pdbOutputs[sr.backbone_index];
          if (originalPdb) {
            filteredOutputs.push({
              id: `scout-pass-${sr.backbone_index}`,
              label: `${originalPdb.label} (scout: pTM=${sr.ptm.toFixed(2)}, pLDDT=${sr.plddt.toFixed(2)})`,
              pdbContent: originalPdb.pdbContent,
              metrics: {
                ...originalPdb.metrics,
                scout_ptm: sr.ptm,
                scout_plddt: sr.plddt,
              },
            });
          }
        }
      }

      ctx.onProgress(100, `Scout: ${result.passing_count}/${result.original_count} passed`);

      return {
        id: `scout-${Date.now()}`,
        summary: `Scout filter: ${result.passing_count}/${result.original_count} backbones passed`,
        pdbOutputs: filteredOutputs,
        data: {
          original_count: result.original_count,
          passing_count: result.passing_count,
          scout_results: result.scout_results,
          ptm_threshold: result.ptm_threshold,
          plddt_threshold: result.plddt_threshold,
          pre_filter_results: result.pre_filter_results,
          pre_filter_passed: result.pre_filter_passed,
          pre_filter_failed: result.pre_filter_failed,
        },
      };
    },

    ...overrides,
  };
}

/**
 * Create a save history step that persists pipeline results to the backend
 * design history store. Runs automatically, no user interaction needed.
 */
export function createSaveHistoryStep(overrides?: Partial<PipelineStepDefinition>): PipelineStepDefinition {
  return {
    id: 'save_history',
    name: 'Save History',
    description: 'Save design results to persistent history for pattern tracking',
    icon: History,
    requiresReview: false,
    supportsSelection: false,
    optional: true,
    defaultParams: {},
    parameterSchema: [],

    async execute(ctx: StepExecutionContext): Promise<StepResult> {
      const { previousResults, initialParams } = ctx;

      ctx.onProgress(10, 'Gathering pipeline results...');

      // Gather intent/config data from earlier steps
      const intentResult = Object.values(previousResults).find(r => r.data?.design_type);
      const designType = intentResult?.data?.design_type as string || 'general';
      const targetMetal = intentResult?.data?.target_metal as string || '';
      const ligandName = intentResult?.data?.ligand_name as string || '';

      // Gather analysis metrics from latest analysis step
      const analysisResult = previousResults['analysis'];
      const totalAnalyzed = analysisResult?.data?.total_analyzed as number || 0;
      const topScore = analysisResult?.data?.top_score as number || 0;

      // Count all PDBs and sequences across the pipeline
      let totalPdbs = 0;
      let totalSequences = 0;
      for (const r of Object.values(previousResults)) {
        if (r.pdbOutputs) totalPdbs += r.pdbOutputs.length;
        if (r.sequences) totalSequences += r.sequences.length;
      }

      // Prefer contig from configure step (actual RFD3 config), fall back to intent
      const configureResult = previousResults['configure'];
      const contig = (configureResult?.data?.contig as string)
        || (intentResult?.data?.contig as string)
        || '';

      const designParams: Record<string, unknown> = {
        design_type: designType,
        target_metal: targetMetal,
        ligand_name: ligandName,
        user_prompt: initialParams.user_prompt,
        contig,
      };

      const designOutputs: Record<string, unknown> = {
        backbone_count: totalPdbs,
        sequence_count: totalSequences,
      };

      const designMetrics: Record<string, unknown> = {
        design_type: designType,
        total_analyzed: totalAnalyzed,
        top_plddt: topScore,
      };

      ctx.onProgress(50, 'Saving to design history...');

      try {
        const saved = await api.saveDesignHistory({
          session_name: `nl_${designType}`,
          design_params: designParams,
          design_outputs: designOutputs,
          design_metrics: designMetrics,
        });

        ctx.onProgress(100, 'History saved');

        return {
          id: `history-${Date.now()}`,
          summary: `Saved to history (run: ${saved.run_id})`,
          data: {
            run_id: saved.run_id,
            session_id: saved.session_id,
          },
        };
      } catch (err) {
        // Non-fatal — don't block the pipeline
        const message = err instanceof Error ? err.message : String(err);
        ctx.onProgress(100, `History save failed (non-fatal): ${message}`);

        return {
          id: `history-fail-${Date.now()}`,
          summary: `History save failed: ${message}`,
          data: { save_error: message },
        };
      }
    },

    ...overrides,
  };
}

/**
 * Create a check lessons step that detects failure patterns,
 * breakthroughs, and improvements in the design history.
 */
export function createCheckLessonsStep(overrides?: Partial<PipelineStepDefinition>): PipelineStepDefinition {
  return {
    id: 'check_lessons',
    name: 'Check Lessons',
    description: 'Detect failure patterns, breakthroughs, or improvements in design history',
    icon: Lightbulb,
    requiresReview: true,
    supportsSelection: false,
    optional: true,
    defaultParams: {},
    parameterSchema: [],

    ResultPreview: LessonResultPreview,

    async execute(ctx: StepExecutionContext): Promise<StepResult> {
      const { previousResults } = ctx;

      ctx.onProgress(10, 'Gathering metrics for lesson detection...');

      // Build a result summary from analysis step
      const analysisResult = previousResults['analysis'];
      const intentResult = Object.values(previousResults).find(r => r.data?.design_type);

      const metricsResult: Record<string, unknown> = {
        design_type: intentResult?.data?.design_type || 'general',
        total_analyzed: analysisResult?.data?.total_analyzed || 0,
        top_plddt: analysisResult?.data?.top_score || 0,
      };

      // Determine outcome (pass/failure) from analysis
      const topScore = analysisResult?.data?.top_score as number || 0;
      if (topScore < 0.5) {
        metricsResult.outcome = 'failure';
      }

      // Include metrics from analysis PDB outputs
      if (analysisResult?.pdbOutputs?.length) {
        const best = analysisResult.pdbOutputs[0];
        if (best.metrics) {
          metricsResult.metrics = best.metrics;
        }
      }

      ctx.onProgress(50, 'Checking for lesson triggers...');

      try {
        const lessons = await api.checkLessons({ result: metricsResult });

        ctx.onProgress(100, lessons.trigger_detected
          ? `Lesson detected: ${lessons.trigger?.type}`
          : 'No patterns detected');

        return {
          id: `lessons-${Date.now()}`,
          summary: lessons.trigger_detected
            ? `Lesson: ${lessons.trigger?.type} — ${lessons.trigger?.description}`
            : `No significant patterns (${lessons.history_count} designs in history)`,
          data: {
            trigger_detected: lessons.trigger_detected,
            trigger: lessons.trigger,
            history_count: lessons.history_count,
          },
        };
      } catch (err) {
        // Non-fatal
        const message = err instanceof Error ? err.message : String(err);
        ctx.onProgress(100, `Lesson check failed (non-fatal): ${message}`);

        return {
          id: `lessons-fail-${Date.now()}`,
          summary: `Lesson check failed: ${message}`,
          data: { check_error: message },
        };
      }
    },

    ...overrides,
  };
}

// ---- Internal Helpers ----

/** Find PDB outputs from previous results, filtered by selected IDs. */
function findPdbOutputs(
  previousResults: Record<string, StepResult>,
  selectedItems: string[],
): PdbOutput[] {
  const allPdbs: PdbOutput[] = [];

  // Iterate in insertion order (which matches step order)
  for (const result of Object.values(previousResults)) {
    if (result.pdbOutputs) {
      allPdbs.push(...result.pdbOutputs);
    }
  }

  if (selectedItems.length > 0) {
    return allPdbs.filter(p => selectedItems.includes(p.id));
  }

  return allPdbs;
}

/**
 * Find PDB outputs from previous results, excluding steps whose data
 * matches a predicate (e.g. scaffold search results that aren't designed backbones).
 */
function findPdbOutputsExcluding(
  previousResults: Record<string, StepResult>,
  selectedItems: string[],
  excludeStepPredicate: (stepId: string, result: StepResult) => boolean,
): PdbOutput[] {
  const allPdbs: PdbOutput[] = [];

  for (const [stepId, result] of Object.entries(previousResults)) {
    if (excludeStepPredicate(stepId, result)) continue;
    if (result.pdbOutputs) {
      allPdbs.push(...result.pdbOutputs);
    }
  }

  if (selectedItems.length > 0) {
    return allPdbs.filter(p => selectedItems.includes(p.id));
  }

  return allPdbs;
}

/** Find sequence outputs from previous results, filtered by selected IDs. */
function findSequenceOutputs(
  previousResults: Record<string, StepResult>,
  selectedItems: string[],
): SequenceOutput[] {
  const allSeqs: SequenceOutput[] = [];

  for (const result of Object.values(previousResults)) {
    if (result.sequences) {
      allSeqs.push(...result.sequences);
    }
  }

  if (selectedItems.length > 0) {
    return allSeqs.filter(s => selectedItems.includes(s.id));
  }

  return allSeqs;
}

/** Parse a FASTA string into individual sequence entries. */
function parseFastaContent(fastaContent: string): Array<{ header: string; sequence: string }> {
  const entries: Array<{ header: string; sequence: string }> = [];
  const lines = fastaContent.trim().split('\n');
  let currentHeader = '';
  let currentSeq = '';

  for (const line of lines) {
    const trimmed = line.trim();
    if (trimmed.startsWith('>')) {
      if (currentHeader && currentSeq) {
        entries.push({ header: currentHeader, sequence: currentSeq });
      }
      currentHeader = trimmed.slice(1).trim();
      currentSeq = '';
    } else if (trimmed) {
      currentSeq += trimmed;
    }
  }
  if (currentHeader && currentSeq) {
    entries.push({ header: currentHeader, sequence: currentSeq });
  }

  return entries;
}

/** Parse MPNN result into SequenceOutput items. */
function parseMpnnSequences(
  result: Record<string, unknown>,
  parentDesignId: string,
  designIndex: number,
): SequenceOutput[] {
  const sequences = result.sequences as Array<Record<string, unknown>> | undefined;
  if (!Array.isArray(sequences)) {
    // Fallback: check for single sequence in result
    if (typeof result.sequence === 'string') {
      return [{
        id: `mpnn-${parentDesignId}-1`,
        sequence: result.sequence as string,
        score: result.score as number | undefined,
        label: `Design ${designIndex + 1} - Seq 1`,
      }];
    }
    return [];
  }

  // Each entry in sequences may contain multi-sequence FASTA content.
  // Parse each entry and split into individual sequences.
  const allParsed: SequenceOutput[] = [];
  for (const s of sequences) {
    const content = (s.content as string) || (s.sequence as string) || '';
    const globalScore = s.global_score as number | undefined;

    // Check if content is multi-sequence FASTA (contains '>' headers)
    if (content.includes('>')) {
      const fastaEntries = parseFastaContent(content);
      for (let fi = 0; fi < fastaEntries.length; fi++) {
        const entry = fastaEntries[fi];
        // Extract score from FASTA header if present (e.g., "score=0.1234" or "global_score=0.1234")
        const scoreMatch = entry.header.match(/(?:global_)?score[=:]\s*([\d.]+)/i);
        const entryScore = scoreMatch ? parseFloat(scoreMatch[1]) : globalScore;

        allParsed.push({
          id: `mpnn-${parentDesignId}-${allParsed.length + 1}`,
          sequence: entry.sequence,
          score: entryScore,
          label: `Design ${designIndex + 1} - Seq ${allParsed.length + 1}`,
          metrics: entryScore !== undefined ? { global_score: entryScore } : undefined,
        });
      }
    } else {
      // Plain sequence string (no FASTA headers)
      allParsed.push({
        id: `mpnn-${parentDesignId}-${allParsed.length + 1}`,
        sequence: content,
        score: globalScore,
        label: `Design ${designIndex + 1} - Seq ${allParsed.length + 1}`,
        metrics: globalScore !== undefined ? { global_score: globalScore } : undefined,
      });
    }
  }

  return allParsed;
}

/** Extract PDB content from RF3 prediction result. */
function extractRf3PdbContent(result: Record<string, unknown>): string | undefined {
  // RF3 returns predictions array
  const predictions = result.predictions as Array<Record<string, unknown>> | undefined;
  if (Array.isArray(predictions) && predictions.length > 0) {
    return (predictions[0].content as string) || (predictions[0].pdb_content as string);
  }
  // Fallback
  return extractPdbContent(result);
}

/** Extract confidence metrics from RF3 result.
 *
 * RF3 returns confidences in varying formats:
 * - Direct: { confidences: { overall_plddt, chain_ptm, overall_pae } }
 * - Nested: { confidences: { summary_confidences: { ptm, overall_plddt } } }
 * - Per-prediction: { predictions: [{ confidences: { ... } }] }
 */
function extractConfidences(result: Record<string, unknown>): {
  overall_plddt?: number;
  ptm?: number;
  overall_pae?: number;
  // AF3 paper ligand-interface metrics (bioRxiv 2025.09.18.676967v2):
  // iPTM > 0.8 and min chain-pair PAE < 1.5 for ligand binding success.
  iptm?: number;
  min_chain_pair_pae?: number;
} | undefined {
  // Try per-prediction confidences first
  const predictions = result.predictions as Array<Record<string, unknown>> | undefined;
  const predConf = predictions?.[0]?.confidences as Record<string, unknown> | undefined;

  // Then top-level confidences
  const topConf = result.confidences as Record<string, unknown> | undefined;
  const confidences = predConf || topConf;
  if (!confidences) return undefined;

  // Try summary_confidences (nested format), then fall back to direct fields
  const summary = confidences.summary_confidences as Record<string, unknown> | undefined;
  const src = summary || confidences;

  // pTM: try ptm, then chain_ptm[0]
  let ptm = src.ptm as number | undefined;
  if (ptm === undefined) {
    const chainPtm = (confidences.chain_ptm || src.chain_ptm) as number[] | undefined;
    if (Array.isArray(chainPtm) && chainPtm.length > 0) {
      ptm = chainPtm[0];
    }
  }

  // iPTM: interface pTM for protein-ligand/multi-chain complexes
  const iptm = (src.iptm as number | undefined) ?? (confidences.iptm as number | undefined);

  // min_chain_pair_pae: pre-computed in backend extract_rf3_confidences from full PAE matrix
  const min_chain_pair_pae = (confidences.min_chain_pair_pae as number | undefined)
    ?? (src.min_chain_pair_pae as number | undefined);

  return {
    overall_plddt: (src.overall_plddt as number | undefined) ?? (src.mean_plddt as number | undefined),
    ptm,
    overall_pae: (src.overall_pae as number | undefined),
    iptm,
    min_chain_pair_pae,
  };
}

/**
 * Find protein residues that coordinate a metal ion in a PDB.
 * Returns a comma-separated string of residue numbers (chain A)
 * that have any heavy atom within `cutoff` Å of a HETATM metal atom.
 *
 * These residues should be fixed during MPNN to preserve the binding site.
 */
function findMetalCoordinatingResidues(pdbContent: string, cutoff = 3.5): string {
  // Parse HETATM positions (metals are single-atom HETATM with element like TB, MG, etc.)
  const metalAtoms: Array<[number, number, number]> = [];
  const lines = pdbContent.split('\n');

  for (const line of lines) {
    if (!line.startsWith('HETATM')) continue;
    const resName = line.substring(17, 20).trim();
    // Metal residues are typically 1-3 letter elements (TB, MG, CA, ZN, FE, etc.)
    // Skip known ligands (CIT, ATP, NAD, HEM, PQQ, etc.)
    if (resName.length > 3 || resName.length === 0) continue;
    // If residue name is all uppercase letters and ≤2 chars, likely a metal
    if (!/^[A-Z]{1,2}$/.test(resName)) continue;

    const x = parseFloat(line.substring(30, 38));
    const y = parseFloat(line.substring(38, 46));
    const z = parseFloat(line.substring(46, 54));
    if (!isNaN(x) && !isNaN(y) && !isNaN(z)) {
      metalAtoms.push([x, y, z]);
    }
  }

  if (metalAtoms.length === 0) return '';

  // Find ATOM residues near metals
  const coordResidues = new Set<number>();
  const cutoffSq = cutoff * cutoff;

  for (const line of lines) {
    if (!line.startsWith('ATOM')) continue;
    const x = parseFloat(line.substring(30, 38));
    const y = parseFloat(line.substring(38, 46));
    const z = parseFloat(line.substring(46, 54));
    if (isNaN(x) || isNaN(y) || isNaN(z)) continue;

    for (const [mx, my, mz] of metalAtoms) {
      const dx = x - mx, dy = y - my, dz = z - mz;
      if (dx * dx + dy * dy + dz * dz <= cutoffSq) {
        const resNum = parseInt(line.substring(22, 26).trim(), 10);
        if (!isNaN(resNum)) coordResidues.add(resNum);
        break;
      }
    }
  }

  return [...coordResidues].sort((a, b) => a - b).join(',');
}
