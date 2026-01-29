/**
 * Shared step factory functions reused across pipeline definitions.
 * Each factory returns a PipelineStepDefinition with a real execute function
 * that calls the backend API.
 */

import { Dna, Shield, BarChart3, Search } from 'lucide-react';
import api from '@/lib/api';
import { extractPdbContent } from '@/lib/workflowHandlers';
import { getMetalMpnnConfig } from './metal-mpnn-bias';
import type { MetalMpnnConfig } from './metal-mpnn-bias';
import type {
  PipelineStepDefinition,
  StepExecutionContext,
  StepResult,
  StepParameterField,
  PdbOutput,
  SequenceOutput,
} from '@/lib/pipeline-types';

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
      num_sequences: 4,
      model_type: 'ligand_mpnn',
    },
    parameterSchema: [
      { id: 'temperature', label: 'Temperature', type: 'slider', required: false, defaultValue: 0.1, range: { min: 0.01, max: 1.0, step: 0.01 }, helpText: 'Lower = more confident sequences' },
      { id: 'num_sequences', label: 'Sequences per Design', type: 'slider', required: false, defaultValue: 4, range: { min: 1, max: 16, step: 1 } },
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

      // Find PDB outputs from previous steps
      const pdbOutputs = findPdbOutputs(previousResults, selectedItems);
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
        const targetMetal = (initialParams.target_metal as string) || (initialParams.targetMetal as string);
        if (targetMetal) {
          metalConfig = getMetalMpnnConfig(targetMetal);
        }
      }

      const allSequences: SequenceOutput[] = [];
      const totalDesigns = pdbOutputs.length;

      for (let i = 0; i < totalDesigns; i++) {
        checkAborted(ctx.abortSignal);

        const pdb = pdbOutputs[i];
        const progressBase = (i / totalDesigns) * 90;
        ctx.onProgress(progressBase + 5, `Designing sequences for ${pdb.label} (${i + 1}/${totalDesigns})...`);

        const response = await api.submitMPNNDesign({
          pdb_content: pdb.pdbContent,
          temperature: (params.temperature as number) ?? 0.1,
          num_sequences: (params.num_sequences as number) ?? 4,
          model_type: (params.model_type as 'ligand_mpnn' | 'protein_mpnn') ?? 'ligand_mpnn',
          ...metalConfig,
        });

        checkAborted(ctx.abortSignal);
        ctx.onJobCreated?.(response.job_id, 'mpnn');

        let jobResult: { status: string; result?: Record<string, unknown>; error?: string };

        if (response.syncCompleted) {
          jobResult = { status: response.status, result: response.result as Record<string, unknown>, error: response.error };
        } else {
          const polled = await api.waitForJob(response.job_id, () => { checkAborted(ctx.abortSignal); });
          jobResult = { status: polled.status, result: polled.result as Record<string, unknown>, error: polled.error };
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

      // Find sequences from previous MPNN step
      const sequences = findSequenceOutputs(previousResults, selectedItems);
      if (sequences.length === 0) {
        throw new Error('No sequences available for validation');
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

      return {
        id: `rf3-${Date.now()}`,
        summary: `Validated ${pdbOutputs.length}/${totalSeqs} sequences successfully`,
        pdbOutputs,
        data: {
          total_validated: pdbOutputs.length,
          total_attempted: totalSeqs,
          success_rate: totalSeqs > 0 ? ((pdbOutputs.length / totalSeqs) * 100).toFixed(1) + '%' : 'N/A',
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

/** Extract confidence metrics from RF3 result. */
function extractConfidences(result: Record<string, unknown>): {
  overall_plddt?: number;
  ptm?: number;
  overall_pae?: number;
} | undefined {
  const confidences = result.confidences as Record<string, unknown> | undefined;
  if (!confidences) return undefined;

  const summary = confidences.summary_confidences as Record<string, unknown> | undefined;
  if (!summary) return undefined;

  return {
    overall_plddt: summary.overall_plddt as number | undefined,
    ptm: summary.ptm as number | undefined,
    overall_pae: summary.overall_pae as number | undefined,
  };
}
