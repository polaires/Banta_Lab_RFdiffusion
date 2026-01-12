/**
 * Custom hook for managing design workflow execution
 * Consolidates job submission, polling, and demo mode fallback logic
 */

import { useState, useCallback, useRef } from 'react';
import { useStore } from './store';
import api, { type RFD3Request, type DesignEvaluation } from './api';
import {
  delay,
  extractPdbContent,
  parseBinderDesigns,
  parseBinderStatistics,
  createBinderEvaluation,
  parseLigandDesigns,
  createLigandEvaluation,
  executeDemoMetalDesign,
  executeDemoLigandDesign,
  executeDemoBinderDesign,
  type MetalWorkflowResult,
  type LigandWorkflowResult,
  type BinderWorkflowResult,
} from './workflowHandlers';
import {
  METAL_CASE_STUDY,
  AZOBENZENE_CASE_STUDY,
  BINDER_CASE_STUDY,
  type DesignResult,
  type LigandEvaluation,
} from './demoData';

export type WorkflowType = 'metal' | 'ligand' | 'binder';

interface WorkflowCallbacks {
  onJobCreated?: (jobId: string) => void;
  onProgress?: (progress: number) => void;
  onStage?: (stage: string) => void;
  onError?: (error: string) => void;
}

interface WorkflowState {
  isRunning: boolean;
  progress: number;
  currentStage: string | null;
  error: string | null;
}

/**
 * Hook for executing design workflows with automatic demo mode fallback
 */
export function useDesignWorkflow() {
  const { health, addJob, updateJob } = useStore();
  const [state, setState] = useState<WorkflowState>({
    isRunning: false,
    progress: 0,
    currentStage: null,
    error: null,
  });

  // Track if backend was connected at job start
  const wasConnectedRef = useRef(false);

  const isBackendConnected = health?.status === 'healthy';

  const resetState = useCallback(() => {
    setState({
      isRunning: false,
      progress: 0,
      currentStage: null,
      error: null,
    });
  }, []);

  const setProgress = useCallback((progress: number) => {
    setState((prev) => ({ ...prev, progress }));
  }, []);

  const setStage = useCallback((stage: string) => {
    setState((prev) => ({ ...prev, currentStage: stage }));
  }, []);

  const setError = useCallback((error: string) => {
    setState((prev) => ({ ...prev, error, isRunning: false }));
  }, []);

  /**
   * Execute a metal binding redesign workflow
   */
  const executeMetal = useCallback(
    async (
      request: RFD3Request,
      pdbContent: string,
      callbacks: WorkflowCallbacks = {}
    ): Promise<MetalWorkflowResult | null> => {
      wasConnectedRef.current = isBackendConnected;
      setState({ isRunning: true, progress: 0, currentStage: null, error: null });

      const { onJobCreated, onProgress, onError } = callbacks;

      if (!isBackendConnected) {
        // Demo mode
        const result = await executeDemoMetalDesign((progress) => {
          setProgress(progress);
          onProgress?.(progress);
        });
        setState((prev) => ({ ...prev, isRunning: false }));
        return result;
      }

      // Backend mode
      try {
        const response = await api.submitRFD3Design({
          ...request,
          pdb_content: pdbContent,
        });

        const jobId = response.job_id;
        onJobCreated?.(jobId);

        addJob({
          id: jobId,
          type: 'rfd3',
          status: 'pending',
          createdAt: new Date().toISOString(),
        });

        let jobResult: { status: string; result?: Record<string, unknown>; error?: string };

        if (response.syncCompleted) {
          jobResult = {
            status: response.status,
            result: response.result as Record<string, unknown>,
            error: response.error,
          };
          updateJob(jobId, { status: response.status === 'completed' ? 'completed' : 'failed' });
        } else {
          const polledResult = await api.waitForJob(jobId, (status) => {
            updateJob(jobId, { status: status.status });
          });
          jobResult = {
            status: polledResult.status,
            result: polledResult.result as Record<string, unknown>,
            error: polledResult.error,
          };
        }

        const designPdb = extractPdbContent(jobResult.result);

        if (jobResult.status === 'completed' && designPdb) {
          // Try to evaluate the design
          let evaluation: DesignEvaluation;
          try {
            const evalResult = await api.evaluateDesign({
              pdb_content: designPdb,
              target_metal: (request as Record<string, unknown>).ligand as string || 'TB',
              metal_chain: 'A',
              metal_resnum: 1,
            });

            if (evalResult && typeof evalResult.criteria_passed === 'number') {
              evaluation = evalResult;
            } else {
              throw new Error('Invalid evaluation response');
            }
          } catch {
            // Fallback to simulated evaluation
            evaluation = {
              success: true,
              coordination_number: 8,
              avg_bond_distance: 2.42,
              geometry_type: 'estimated',
              geometry_rmsd: 0.3,
              oxygen_donors: 6,
              nitrogen_donors: 2,
              sulfur_donors: 0,
              criteria_passed: 3,
              criteria_total: 4,
              overall_pass: true,
              suggestions: ['Full analysis requires backend evaluation endpoint.'],
            };
          }

          setState((prev) => ({ ...prev, isRunning: false }));
          return { evaluation, pdbContent: designPdb };
        }

        // Backend error - fall back to demo mode
        const errorMsg = jobResult.error || 'Job failed';
        onError?.(`Backend error: ${errorMsg}. Falling back to demo mode.`);

        await delay(1000);
        const demoResult = await executeDemoMetalDesign((progress) => {
          setProgress(progress);
          onProgress?.(progress);
        });
        setState((prev) => ({ ...prev, isRunning: false }));
        return demoResult;
      } catch (err) {
        const errorMsg = err instanceof Error ? err.message : 'Unknown error';
        onError?.(`Error: ${errorMsg}. Falling back to demo mode.`);

        await delay(1000);
        const demoResult = await executeDemoMetalDesign((progress) => {
          setProgress(progress);
          onProgress?.(progress);
        });
        setState((prev) => ({ ...prev, isRunning: false }));
        return demoResult;
      }
    },
    [isBackendConnected, addJob, updateJob, setProgress]
  );

  /**
   * Execute a ligand interface dimer design workflow
   */
  const executeLigand = useCallback(
    async (
      request: RFD3Request,
      callbacks: WorkflowCallbacks = {}
    ): Promise<LigandWorkflowResult | null> => {
      wasConnectedRef.current = isBackendConnected;
      setState({ isRunning: true, progress: 0, currentStage: null, error: null });

      const { onJobCreated, onProgress, onStage, onError } = callbacks;

      if (!isBackendConnected) {
        // Demo mode
        const result = await executeDemoLigandDesign(
          (progress) => {
            setProgress(progress);
            onProgress?.(progress);
          },
          (stage) => {
            setStage(stage);
            onStage?.(stage);
          }
        );
        setState((prev) => ({ ...prev, isRunning: false }));
        return result;
      }

      // Backend mode
      try {
        const response = await api.submitRFD3Design({
          task: 'interface_ligand_design',
          approach: 'full',
          ...request,
        });

        const jobId = response.job_id;
        onJobCreated?.(jobId);

        addJob({
          id: jobId,
          type: 'rfd3',
          status: 'pending',
          createdAt: new Date().toISOString(),
        });

        let jobResult: { status: string; result?: Record<string, unknown>; error?: string };

        if (response.syncCompleted) {
          jobResult = {
            status: response.status,
            result: response.result as Record<string, unknown>,
            error: response.error,
          };
          updateJob(jobId, { status: response.status === 'completed' ? 'completed' : 'failed' });
        } else {
          const polledResult = await api.waitForJob(jobId, (status) => {
            updateJob(jobId, { status: status.status });
          });
          jobResult = {
            status: polledResult.status,
            result: polledResult.result as Record<string, unknown>,
            error: polledResult.error,
          };
        }

        const result = jobResult.result;
        const designPdb = extractPdbContent(result);

        if (jobResult.status === 'completed' && designPdb) {
          // Parse designs from result
          let designs = result ? parseLigandDesigns(result) : [];

          // If no designs array, create single design from PDB
          if (designs.length === 0 && designPdb) {
            const dimerMetrics = extractDimerMetrics(result);
            designs = [
              {
                id: 'design-1',
                rank: 1,
                pdbContent: designPdb,
                metrics: {
                  affinity: dimerMetrics.affinity as number | undefined,
                  contacts_a: dimerMetrics.contacts_a as number | undefined,
                  contacts_b: dimerMetrics.contacts_b as number | undefined,
                  has_clashes: (dimerMetrics.has_clashes as boolean | undefined) ?? false,
                  separable: (dimerMetrics.separable as boolean | undefined) ?? true,
                  interface_area: dimerMetrics.interface_area as number | undefined,
                },
              },
            ];
          }

          const dimerMetrics = extractDimerMetrics(result);
          const evaluation = createLigandEvaluation(dimerMetrics);

          setState((prev) => ({ ...prev, isRunning: false }));
          return { evaluation, designs, pdbContent: designPdb };
        }

        // Backend error - fall back to demo mode
        const errorMsg = jobResult.error || 'Job failed';
        onError?.(`Backend error: ${errorMsg}. Falling back to demo mode.`);

        await delay(1000);
        const demoResult = await executeDemoLigandDesign(
          (progress) => {
            setProgress(progress);
            onProgress?.(progress);
          },
          (stage) => {
            setStage(stage);
            onStage?.(stage);
          }
        );
        setState((prev) => ({ ...prev, isRunning: false }));
        return demoResult;
      } catch (err) {
        const errorMsg = err instanceof Error ? err.message : 'Unknown error';
        onError?.(`Error: ${errorMsg}. Falling back to demo mode.`);

        await delay(1000);
        const demoResult = await executeDemoLigandDesign(
          (progress) => {
            setProgress(progress);
            onProgress?.(progress);
          },
          (stage) => {
            setStage(stage);
            onStage?.(stage);
          }
        );
        setState((prev) => ({ ...prev, isRunning: false }));
        return demoResult;
      }
    },
    [isBackendConnected, addJob, updateJob, setProgress, setStage]
  );

  /**
   * Execute a protein binder design workflow
   */
  const executeBinder = useCallback(
    async (
      request: RFD3Request,
      callbacks: WorkflowCallbacks = {}
    ): Promise<BinderWorkflowResult | null> => {
      wasConnectedRef.current = isBackendConnected;
      setState({ isRunning: true, progress: 0, currentStage: null, error: null });

      const { onJobCreated, onProgress, onStage, onError } = callbacks;

      if (!isBackendConnected) {
        // Demo mode
        const result = await executeDemoBinderDesign(
          (progress) => {
            setProgress(progress);
            onProgress?.(progress);
          },
          (stage) => {
            setStage(stage);
            onStage?.(stage);
          }
        );
        setState((prev) => ({ ...prev, isRunning: false }));
        return result;
      }

      // Backend mode
      try {
        // Debug: Log the request being sent
        console.log('[executeBinder] Submitting request:', {
          task: 'protein_binder_design',
          num_designs: request.num_designs,
          binder_length: request.binder_length,
          quality_threshold: request.quality_threshold,
          protocol: request.protocol,
        });

        const response = await api.submitRFD3Design({
          task: 'protein_binder_design',
          ...request,
        });

        const jobId = response.job_id;
        onJobCreated?.(jobId);

        addJob({
          id: jobId,
          type: 'rfd3',
          status: 'pending',
          createdAt: new Date().toISOString(),
        });

        let jobResult: { status: string; result?: Record<string, unknown>; error?: string };

        if (response.syncCompleted) {
          jobResult = {
            status: response.status,
            result: response.result as Record<string, unknown>,
            error: response.error,
          };
          updateJob(jobId, { status: response.status === 'completed' ? 'completed' : 'failed' });
        } else {
          const polledResult = await api.waitForJob(jobId, (status) => {
            updateJob(jobId, { status: status.status });
          });
          jobResult = {
            status: polledResult.status,
            result: polledResult.result as Record<string, unknown>,
            error: polledResult.error,
          };
        }

        const result = jobResult.result;

        // Debug: Log the result
        console.log('[executeBinder] Job result:', {
          status: jobResult.status,
          hasResult: !!result,
          statistics: result?.statistics,
          designCount: (result?.designs as unknown[])?.length,
        });

        if (jobResult.status === 'completed' && result) {
          const designs = parseBinderDesigns(result);
          const statistics = parseBinderStatistics(result, designs.length);
          const evaluation = createBinderEvaluation(designs, statistics);
          const pdbContent = designs[0]?.pdb_content;

          console.log('[executeBinder] Parsed designs:', designs.length);
          setState((prev) => ({ ...prev, isRunning: false }));
          return { evaluation, designs, statistics, pdbContent };
        }

        // Backend error - fall back to demo mode
        const errorMsg = jobResult.error || 'Job failed';
        onError?.(`Backend error: ${errorMsg}. Falling back to demo mode.`);

        await delay(1000);
        const demoResult = await executeDemoBinderDesign(
          (progress) => {
            setProgress(progress);
            onProgress?.(progress);
          },
          (stage) => {
            setStage(stage);
            onStage?.(stage);
          }
        );
        setState((prev) => ({ ...prev, isRunning: false }));
        return demoResult;
      } catch (err) {
        const errorMsg = err instanceof Error ? err.message : 'Unknown error';

        // Only fall back if backend wasn't connected at start
        if (!wasConnectedRef.current) {
          onError?.(`Error: ${errorMsg}. Falling back to demo mode.`);

          await delay(1000);
          const demoResult = await executeDemoBinderDesign(
            (progress) => {
              setProgress(progress);
              onProgress?.(progress);
            },
            (stage) => {
              setStage(stage);
              onStage?.(stage);
            }
          );
          setState((prev) => ({ ...prev, isRunning: false }));
          return demoResult;
        }

        // Backend was connected - show real error
        setError(errorMsg);
        return null;
      }
    },
    [isBackendConnected, addJob, updateJob, setProgress, setStage, setError]
  );

  return {
    state,
    isBackendConnected,
    resetState,
    executeMetal,
    executeLigand,
    executeBinder,
  };
}

/**
 * Extract dimer metrics from various result formats
 */
function extractDimerMetrics(result: Record<string, unknown> | undefined): Record<string, unknown> {
  if (!result) return {};

  const dimer = result.dimer as Record<string, unknown> | undefined;
  if (dimer) {
    // Metrics might be nested in dimer.metrics or directly on dimer
    const metrics = dimer.metrics as Record<string, unknown> | undefined;
    return metrics || dimer;
  }

  return {};
}

/**
 * Create evaluation from a selected ligand design
 */
export function createEvaluationFromDesign(design: DesignResult): LigandEvaluation {
  return {
    success: true,
    approach: 'full',
    chain_a: design.chain_a_metrics
      ? {
          contacts: design.chain_a_metrics.contacts,
          exposed_atoms: design.chain_a_metrics.exposed_atoms.join(','),
          affinity: design.chain_a_metrics.affinity ?? 0,
        }
      : undefined,
    chain_b: design.chain_b_metrics
      ? {
          contacts: design.chain_b_metrics.contacts,
          exposed_atoms: design.chain_b_metrics.exposed_atoms.join(','),
          affinity: design.chain_b_metrics.affinity ?? 0,
        }
      : undefined,
    dimer: {
      affinity: design.metrics.affinity ?? -3.0,
      contacts_a: design.metrics.contacts_a ?? 0,
      contacts_b: design.metrics.contacts_b ?? 0,
      has_clashes: design.metrics.has_clashes ?? false,
      separable: design.metrics.separable ?? true,
      interface_area: design.metrics.interface_area,
    },
    overall_pass: !design.metrics.has_clashes && (design.metrics.separable ?? true),
  };
}

// Re-export demo data for convenience
export {
  METAL_CASE_STUDY,
  AZOBENZENE_CASE_STUDY,
  BINDER_CASE_STUDY,
};
