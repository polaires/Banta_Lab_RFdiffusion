/**
 * Workflow handlers for the AI Design Assistant
 * Consolidates common patterns across Metal, Ligand, and Binder workflows
 */

import api from './api';
import {
  METAL_CASE_STUDY,
  AZOBENZENE_CASE_STUDY,
  BINDER_CASE_STUDY,
  type DesignResult,
  type BinderDesign,
  type BinderStatistics,
  type BinderEvaluation,
  type LigandEvaluation,
} from './demoData';
import type { DesignEvaluation, RFD3Request } from './api';

// Common delay utility
export function delay(ms: number): Promise<void> {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

// Types for workflow results
export interface WorkflowJobResult<T> {
  success: boolean;
  result?: T;
  error?: string;
  jobId?: string;
  pdbContent?: string;
}

export interface MetalWorkflowResult {
  evaluation: DesignEvaluation;
  pdbContent?: string;
}

export interface LigandWorkflowResult {
  evaluation: LigandEvaluation;
  designs: DesignResult[];
  pdbContent?: string;
}

export interface BinderWorkflowResult {
  evaluation: BinderEvaluation;
  designs: BinderDesign[];
  statistics: BinderStatistics;
  pdbContent?: string;
}

/**
 * Execute a design workflow with automatic demo mode fallback
 */
export async function executeWorkflow<T>(
  isBackendConnected: boolean,
  executeBackend: () => Promise<WorkflowJobResult<T>>,
  executeDemoMode: () => Promise<T>,
  options: {
    onError?: (error: string) => void;
    fallbackOnError?: boolean;
  } = {}
): Promise<WorkflowJobResult<T>> {
  const { onError, fallbackOnError = true } = options;

  if (!isBackendConnected) {
    const result = await executeDemoMode();
    return { success: true, result };
  }

  try {
    const result = await executeBackend();
    if (result.success) {
      return result;
    }

    // Backend returned error - fall back to demo mode if enabled
    if (fallbackOnError) {
      onError?.(result.error || 'Backend error, using demo mode');
      const demoResult = await executeDemoMode();
      return { success: true, result: demoResult };
    }

    return result;
  } catch (err) {
    const errorMsg = err instanceof Error ? err.message : 'Unknown error';
    onError?.(errorMsg);

    if (fallbackOnError) {
      const demoResult = await executeDemoMode();
      return { success: true, result: demoResult };
    }

    return { success: false, error: errorMsg };
  }
}

/**
 * Submit a design job and wait for result
 * Handles both sync and async response modes
 */
export async function submitAndWaitForJob(
  request: RFD3Request,
  callbacks: {
    onJobCreated?: (jobId: string) => void;
    onStatusUpdate?: (status: string) => void;
  } = {}
): Promise<WorkflowJobResult<Record<string, unknown>>> {
  const { onJobCreated, onStatusUpdate } = callbacks;

  try {
    const response = await api.submitRFD3Design(request);
    const jobId = response.job_id;

    onJobCreated?.(jobId);

    let jobResult: { status: string; result?: Record<string, unknown>; error?: string };

    if (response.syncCompleted) {
      jobResult = {
        status: response.status,
        result: response.result as Record<string, unknown>,
        error: response.error,
      };
    } else {
      const polledResult = await api.waitForJob(jobId, (status) => {
        onStatusUpdate?.(status.status);
      });
      jobResult = {
        status: polledResult.status,
        result: polledResult.result as Record<string, unknown>,
        error: polledResult.error,
      };
    }

    if (jobResult.status === 'completed' && jobResult.result) {
      return {
        success: true,
        result: jobResult.result,
        jobId,
        pdbContent: extractPdbContent(jobResult.result),
      };
    }

    return {
      success: false,
      error: jobResult.error || 'Job failed',
      jobId,
    };
  } catch (err) {
    const errorMsg = err instanceof Error ? err.message : 'Unknown error';
    return { success: false, error: errorMsg };
  }
}

/**
 * Extract PDB content from various result formats
 */
export function extractPdbContent(result: unknown): string | undefined {
  if (!result || typeof result !== 'object') return undefined;

  const r = result as Record<string, unknown>;

  // Check designs array format
  if (Array.isArray(r.designs) && r.designs.length > 0) {
    const firstDesign = r.designs[0] as Record<string, unknown>;
    if (typeof firstDesign.content === 'string') return firstDesign.content;
    if (typeof firstDesign.pdb_content === 'string') return firstDesign.pdb_content;
  }

  // Check direct pdb field
  if (typeof r.pdb === 'string') return r.pdb;

  // Check dimer format (ligand design)
  if (r.dimer && typeof r.dimer === 'object') {
    const dimer = r.dimer as Record<string, unknown>;
    if (typeof dimer.pdb_content === 'string') return dimer.pdb_content;
  }

  // Check if result itself is a string
  if (typeof result === 'string') return result;

  return undefined;
}

/**
 * Parse binder designs from backend response
 */
export function parseBinderDesigns(result: Record<string, unknown>): BinderDesign[] {
  const rawDesigns = result.designs as Array<Record<string, unknown>> | undefined;
  if (!Array.isArray(rawDesigns)) return [];

  return rawDesigns.map((d, idx) => ({
    rank: (d.rank as number) || idx + 1,
    binder_sequence: (d.binder_sequence as string) || (d.sequence as string) || '',
    mpnn_score: d.mpnn_score as number,
    esm_perplexity: d.esm_perplexity as number,
    esm_confidence: d.esm_confidence as number,
    interface_contacts: d.interface_contacts as number,
    interface_hbonds: d.interface_hbonds as number,
    buried_sasa: d.buried_sasa as number,
    packstat: d.packstat as number,
    pdb_content: d.pdb_content as string | undefined,
    shape_complementarity: d.shape_complementarity as number | undefined,
    surface_hydrophobicity: d.surface_hydrophobicity as number | undefined,
    unsaturated_hbonds: d.unsaturated_hbonds as number | undefined,
    interface_residue_count: d.interface_residue_count as number | undefined,
    esmfold_plddt: d.esmfold_plddt as number | undefined,
    esmfold_rmsd: d.esmfold_rmsd as number | undefined,
    esmfold_validation_passed: d.esmfold_validation_passed as boolean | undefined,
  }));
}

/**
 * Parse binder statistics from backend response
 */
export function parseBinderStatistics(
  result: Record<string, unknown>,
  designCount: number
): BinderStatistics {
  const stats = result.statistics as Record<string, number> | undefined;

  if (stats) {
    return {
      generated: stats.generated || designCount,
      mpnn_designed: stats.mpnn_designed || designCount,
      esm_passed: stats.esm_passed || designCount,
      relaxed: stats.relaxed || designCount,
      interface_analyzed: stats.interface_analyzed || designCount,
      passed_filters: stats.passed_filters || designCount,
      returned: stats.returned || designCount,
    };
  }

  return {
    generated: designCount,
    mpnn_designed: designCount,
    esm_passed: designCount,
    relaxed: designCount,
    interface_analyzed: designCount,
    passed_filters: designCount,
    returned: designCount,
  };
}

/**
 * Create binder evaluation from designs and statistics
 */
export function createBinderEvaluation(
  designs: BinderDesign[],
  statistics: BinderStatistics
): BinderEvaluation {
  const bestDesign = designs[0];

  return {
    success: true,
    statistics,
    best_design: bestDesign
      ? {
          rank: bestDesign.rank,
          binder_sequence: bestDesign.binder_sequence,
          mpnn_score: bestDesign.mpnn_score,
          esm_perplexity: bestDesign.esm_perplexity,
          esm_confidence: bestDesign.esm_confidence,
          interface_contacts: bestDesign.interface_contacts,
          interface_hbonds: bestDesign.interface_hbonds,
          buried_sasa: bestDesign.buried_sasa,
          packstat: bestDesign.packstat,
        }
      : undefined,
    overall_pass: designs.length > 0,
  };
}

/**
 * Parse ligand designs from backend response
 */
export function parseLigandDesigns(result: Record<string, unknown>): DesignResult[] {
  const rawDesigns = result.designs as Array<Record<string, unknown>> | undefined;
  if (!Array.isArray(rawDesigns)) return [];

  return rawDesigns.map((d, idx) => ({
    id: `design-${idx + 1}`,
    rank: idx + 1,
    pdbContent: (d.pdb_content as string) || (d.content as string) || '',
    metrics: {
      affinity: d.metrics
        ? (d.metrics as Record<string, unknown>).affinity as number | undefined
        : undefined,
      contacts_a: d.metrics
        ? (d.metrics as Record<string, unknown>).contacts_a as number | undefined
        : undefined,
      contacts_b: d.metrics
        ? (d.metrics as Record<string, unknown>).contacts_b as number | undefined
        : undefined,
      has_clashes: d.metrics
        ? ((d.metrics as Record<string, unknown>).has_clashes as boolean) || false
        : false,
      separable: d.metrics
        ? ((d.metrics as Record<string, unknown>).separable as boolean) ?? true
        : true,
      interface_area: d.metrics
        ? (d.metrics as Record<string, unknown>).interface_area as number | undefined
        : undefined,
    },
  }));
}

/**
 * Create ligand evaluation from dimer metrics
 */
export function createLigandEvaluation(
  dimerMetrics: Record<string, unknown>
): LigandEvaluation {
  return {
    success: true,
    approach: 'full',
    dimer: {
      affinity: (dimerMetrics.affinity as number) ?? -3.5,
      contacts_a: (dimerMetrics.contacts_a as number) ?? 6,
      contacts_b: (dimerMetrics.contacts_b as number) ?? 5,
      has_clashes: (dimerMetrics.has_clashes as boolean) || false,
      separable: (dimerMetrics.separable as boolean) ?? true,
      interface_area: dimerMetrics.interface_area as number | undefined,
    },
    overall_pass:
      !(dimerMetrics.has_clashes as boolean) &&
      ((dimerMetrics.separable as boolean) ?? true),
  };
}

/**
 * Format evaluation message for metal binding design
 */
export function formatMetalEvaluationMessage(
  evaluation: DesignEvaluation,
  targetMetalLabel: string
): string {
  const passCount = evaluation.criteria_passed;
  if (evaluation.overall_pass) {
    return `Excellent results! The design meets ${passCount}/${evaluation.criteria_total} target criteria for ${targetMetalLabel} binding. The coordination geometry looks ideal.`;
  }
  return `The design meets ${passCount}/${evaluation.criteria_total} criteria for ${targetMetalLabel} binding. Consider adjusting settings.`;
}

/**
 * Format evaluation message for ligand design
 */
export function formatLigandEvaluationMessage(
  evaluation: LigandEvaluation,
  designCount: number,
  passingCount: number
): string {
  const dimer = evaluation.dimer;
  if (!dimer) {
    return `Design complete! Generated ${designCount} designs (${passingCount} pass criteria).`;
  }

  return `Design complete! Generated ${designCount} design${designCount > 1 ? 's' : ''} (${passingCount} pass criteria).

**Best Design:**
- Affinity: ${dimer.affinity.toFixed(2)} kcal/mol
- Separable: ${dimer.separable ? 'Yes' : 'No'}
- Clashes: ${dimer.has_clashes ? 'Detected' : 'None'}

The structure is now visible in the viewer below!`;
}

/**
 * Format evaluation message for binder design
 */
export function formatBinderEvaluationMessage(
  statistics: BinderStatistics,
  bestDesign: BinderDesign | undefined
): string {
  if (!bestDesign) {
    return `Design complete! Generated ${statistics.generated} designs, but none passed all filters.`;
  }

  return `Design complete! Generated ${statistics.generated} designs, ${statistics.returned} passed filters.

**Best Design (#${bestDesign.rank}):**
- ESM Confidence: ${bestDesign.esm_confidence ? `${(bestDesign.esm_confidence * 100).toFixed(0)}%` : 'N/A'}
- Interface Contacts: ${bestDesign.interface_contacts || 'N/A'}
- H-Bonds: ${bestDesign.interface_hbonds || 'N/A'}

Results are displayed below.`;
}

// Demo mode execution functions

/**
 * Execute metal design in demo mode
 */
export async function executeDemoMetalDesign(
  onProgress: (progress: number) => void
): Promise<MetalWorkflowResult> {
  for (let i = 0; i <= 100; i += 20) {
    await delay(1000);
    onProgress(i);
  }

  return {
    evaluation: METAL_CASE_STUDY.simulatedEvaluation,
  };
}

/**
 * Execute ligand design in demo mode
 */
export async function executeDemoLigandDesign(
  onProgress: (progress: number) => void,
  onStage?: (stage: string) => void
): Promise<LigandWorkflowResult> {
  onStage?.('Step 1: Designing Chain A...');

  for (let i = 0; i <= 50; i += 10) {
    await delay(800);
    onProgress(i);
  }

  onStage?.('Step 2: Designing Chain B...');

  for (let i = 50; i <= 100; i += 10) {
    await delay(800);
    onProgress(i);
  }

  return {
    evaluation: AZOBENZENE_CASE_STUDY.simulatedEvaluation,
    designs: AZOBENZENE_CASE_STUDY.simulatedDesigns,
  };
}

/**
 * Execute binder design in demo mode
 */
export async function executeDemoBinderDesign(
  onProgress: (progress: number) => void,
  onStage?: (stage: string) => void
): Promise<BinderWorkflowResult> {
  onStage?.('Stage 1: Generating backbones with RFdiffusion...');
  for (let i = 0; i <= 25; i += 5) {
    await delay(400);
    onProgress(i);
  }

  onStage?.('Stage 2: Designing sequences with ProteinMPNN...');
  for (let i = 25; i <= 50; i += 5) {
    await delay(400);
    onProgress(i);
  }

  onStage?.('Stage 3: Scoring with ESM-3...');
  for (let i = 50; i <= 75; i += 5) {
    await delay(400);
    onProgress(i);
  }

  onStage?.('Stage 4: Analyzing interfaces...');
  for (let i = 75; i <= 100; i += 5) {
    await delay(400);
    onProgress(i);
  }

  return {
    evaluation: BINDER_CASE_STUDY.simulatedEvaluation,
    designs: BINDER_CASE_STUDY.simulatedDesigns,
    statistics: BINDER_CASE_STUDY.simulatedStatistics,
    pdbContent: BINDER_CASE_STUDY.simulatedDesigns[0]?.pdb_content,
  };
}
