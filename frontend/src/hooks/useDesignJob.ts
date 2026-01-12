/**
 * Shared hook for design job submission and polling
 * Consolidates common patterns across Metal, Ligand, and Binder workflows
 */

import { useState, useCallback, useRef } from 'react';
import api, { type JobStatus, type RFD3Request } from '@/lib/api';
import { useStore } from '@/lib/store';

export type DesignJobStatus = 'idle' | 'submitting' | 'polling' | 'completed' | 'failed';

export interface DesignJobResult<T = unknown> {
  status: 'completed' | 'failed';
  result?: T;
  error?: string;
  jobId: string;
  syncCompleted?: boolean;
}

export interface UseDesignJobOptions {
  onProgress?: (progress: number) => void;
  onStatusChange?: (status: JobStatus) => void;
  pollInterval?: number;
}

export interface UseDesignJobReturn<T = unknown> {
  status: DesignJobStatus;
  isRunning: boolean;
  error: string | null;
  result: T | null;
  jobId: string | null;
  submitJob: (request: RFD3Request) => Promise<DesignJobResult<T>>;
  reset: () => void;
}

/**
 * Hook for managing design job submission and polling
 * Supports both synchronous (traditional) and asynchronous (serverless) modes
 */
export function useDesignJob<T = unknown>(
  options: UseDesignJobOptions = {}
): UseDesignJobReturn<T> {
  const { onProgress, onStatusChange, pollInterval = 2000 } = options;

  const [status, setStatus] = useState<DesignJobStatus>('idle');
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<T | null>(null);
  const [jobId, setJobId] = useState<string | null>(null);

  const { addJob, updateJob } = useStore();

  const abortControllerRef = useRef<AbortController | null>(null);

  const isRunning = status === 'submitting' || status === 'polling';

  const reset = useCallback(() => {
    if (abortControllerRef.current) {
      abortControllerRef.current.abort();
      abortControllerRef.current = null;
    }
    setStatus('idle');
    setError(null);
    setResult(null);
    setJobId(null);
  }, []);

  const submitJob = useCallback(
    async (request: RFD3Request): Promise<DesignJobResult<T>> => {
      setStatus('submitting');
      setError(null);
      setResult(null);

      abortControllerRef.current = new AbortController();

      try {
        const response = await api.submitRFD3Design(request);
        const newJobId = response.job_id;
        setJobId(newJobId);

        addJob({
          id: newJobId,
          type: 'rfd3',
          status: 'pending',
          createdAt: new Date().toISOString(),
        });

        let jobResult: { status: string; result?: T; error?: string };

        if (response.syncCompleted) {
          // Synchronous mode - result already available
          jobResult = {
            status: response.status,
            result: response.result as T,
            error: response.error,
          };
          updateJob(newJobId, {
            status: response.status === 'completed' ? 'completed' : 'failed',
          });
        } else {
          // Asynchronous mode - poll for completion
          setStatus('polling');
          const polledResult = await api.waitForJob(
            newJobId,
            (pollStatus) => {
              updateJob(newJobId, { status: pollStatus.status });
              onStatusChange?.(pollStatus);
            },
            pollInterval
          );
          jobResult = {
            status: polledResult.status,
            result: polledResult.result as T,
            error: polledResult.error,
          };
        }

        if (jobResult.status === 'completed' && jobResult.result) {
          setResult(jobResult.result);
          setStatus('completed');
          return {
            status: 'completed',
            result: jobResult.result,
            jobId: newJobId,
            syncCompleted: response.syncCompleted,
          };
        } else {
          const errorMsg = jobResult.error || 'Job failed with unknown error';
          setError(errorMsg);
          setStatus('failed');
          return {
            status: 'failed',
            error: errorMsg,
            jobId: newJobId,
          };
        }
      } catch (err) {
        const errorMsg = err instanceof Error ? err.message : 'Unknown error';
        setError(errorMsg);
        setStatus('failed');
        return {
          status: 'failed',
          error: errorMsg,
          jobId: jobId || 'unknown',
        };
      }
    },
    [addJob, updateJob, onStatusChange, pollInterval, jobId]
  );

  return {
    status,
    isRunning,
    error,
    result,
    jobId,
    submitJob,
    reset,
  };
}

/**
 * Simulates job progress for demo mode
 * Returns a cleanup function to cancel the simulation
 */
export function simulateProgress(
  onProgress: (progress: number) => void,
  stages: { name: string; duration: number; startProgress: number; endProgress: number }[],
  onStageChange?: (stageName: string) => void
): () => void {
  let cancelled = false;
  let currentTimeout: ReturnType<typeof setTimeout> | null = null;

  async function runSimulation(): Promise<void> {
    for (const stage of stages) {
      if (cancelled) break;

      onStageChange?.(stage.name);

      const steps = Math.ceil((stage.endProgress - stage.startProgress) / 5);
      const stepDuration = stage.duration / steps;

      for (let i = 0; i <= steps; i++) {
        if (cancelled) break;

        const progress = stage.startProgress + (i / steps) * (stage.endProgress - stage.startProgress);
        onProgress(Math.round(progress));

        if (i < steps) {
          await new Promise<void>((resolve) => {
            currentTimeout = setTimeout(resolve, stepDuration);
          });
        }
      }
    }
  }

  runSimulation();

  return () => {
    cancelled = true;
    if (currentTimeout) {
      clearTimeout(currentTimeout);
    }
  };
}

/**
 * Helper to extract PDB content from various response formats
 */
export function extractPdbFromResult(result: unknown): string | undefined {
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
