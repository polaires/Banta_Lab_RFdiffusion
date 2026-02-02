'use client';

import { useReducer, useCallback, useRef, useEffect } from 'react';
import type {
  PipelineDefinition,
  PipelineStepDefinition,
  PipelineRuntimeState,
  StepRuntimeState,
  StepStatus,
  StepResult,
  StepExecutionContext,
  SequenceOutput,
} from '@/lib/pipeline-types';
import { computeProtParam } from '@/lib/protparam';

// ---- FIX #14: Pipeline State Persistence ----

const PIPELINE_STORAGE_KEY = 'pipeline-runtime-state';

/** Serializable subset of runtime state (excludes non-serializable fields). */
function savePipelineState(state: PipelineRuntimeState): void {
  try {
    // Only persist at stable points (paused, completed, failed)
    if (state.status === 'running' || state.status === 'idle') return;
    // PDB content can be large — truncate if needed
    const serializable = {
      ...state,
      steps: state.steps.map(s => ({
        ...s,
        // Cap PDB content size to prevent storage quota issues
        result: s.result ? {
          ...s.result,
          pdbOutputs: s.result.pdbOutputs?.map(p => ({
            ...p,
            pdbContent: p.pdbContent.length > 100_000 ? '' : p.pdbContent,
          })),
        } : undefined,
      })),
    };
    sessionStorage.setItem(PIPELINE_STORAGE_KEY, JSON.stringify(serializable));
  } catch {
    // Storage full or unavailable — silently ignore
  }
}

function loadPipelineState(): PipelineRuntimeState | null {
  try {
    const stored = sessionStorage.getItem(PIPELINE_STORAGE_KEY);
    if (!stored) return null;
    return JSON.parse(stored) as PipelineRuntimeState;
  } catch {
    return null;
  }
}

function clearPipelineState(): void {
  try {
    sessionStorage.removeItem(PIPELINE_STORAGE_KEY);
  } catch {
    // Ignore
  }
}

// ---- Reducer Actions ----

type PipelineAction =
  | { type: 'INIT'; pipelineId: string; sessionId: string; steps: StepRuntimeState[]; initialParams: Record<string, unknown> }
  | { type: 'STEP_START'; stepIndex: number }
  | { type: 'STEP_PROGRESS'; stepIndex: number; progress: number; message?: string }
  | { type: 'STEP_COMPLETE'; stepIndex: number; result: StepResult }
  | { type: 'STEP_PAUSE'; stepIndex: number }
  | { type: 'STEP_FAIL'; stepIndex: number; error: string }
  | { type: 'SET_PARAMS'; stepId: string; params: Record<string, unknown> }
  | { type: 'SET_SELECTION'; stepId: string; outputIds: string[] }
  | { type: 'CONFIRM'; stepIndex: number }
  | { type: 'SKIP'; stepIndex: number }
  | { type: 'CANCEL' }
  | { type: 'RESET' };

const INITIAL_STATE: PipelineRuntimeState = {
  pipelineId: '',
  sessionId: '',
  status: 'idle',
  activeStepIndex: 0,
  steps: [],
  initialParams: {},
};

function pipelineReducer(state: PipelineRuntimeState, action: PipelineAction): PipelineRuntimeState {
  switch (action.type) {
    case 'INIT':
      return {
        pipelineId: action.pipelineId,
        sessionId: action.sessionId,
        status: 'idle',
        activeStepIndex: 0,
        steps: action.steps,
        initialParams: action.initialParams,
        startedAt: Date.now(),
      };

    case 'STEP_START': {
      const steps = state.steps.map((s, i) =>
        i === action.stepIndex
          ? { ...s, status: 'running' as StepStatus, progress: 0, progressMessage: undefined, error: undefined, startedAt: Date.now() }
          : s
      );
      return { ...state, status: 'running', activeStepIndex: action.stepIndex, steps };
    }

    case 'STEP_PROGRESS': {
      const steps = state.steps.map((s, i) =>
        i === action.stepIndex
          ? { ...s, progress: action.progress, progressMessage: action.message }
          : s
      );
      return { ...state, steps };
    }

    case 'STEP_COMPLETE': {
      const steps = state.steps.map((s, i) =>
        i === action.stepIndex
          ? { ...s, status: 'completed' as StepStatus, progress: 100, result: action.result, completedAt: Date.now() }
          : s
      );
      const allDone = steps.every(s => s.status === 'completed' || s.status === 'skipped');
      return {
        ...state,
        status: allDone ? 'completed' : 'running',
        steps,
        completedAt: allDone ? Date.now() : undefined,
      };
    }

    case 'STEP_PAUSE': {
      const steps = state.steps.map((s, i) =>
        i === action.stepIndex
          ? { ...s, status: 'paused' as StepStatus, progress: 100 }
          : s
      );
      return { ...state, status: 'paused', activeStepIndex: action.stepIndex, steps };
    }

    case 'STEP_FAIL': {
      const steps = state.steps.map((s, i) =>
        i === action.stepIndex
          ? { ...s, status: 'failed' as StepStatus, error: action.error, completedAt: Date.now() }
          : s
      );
      return { ...state, status: 'failed', steps };
    }

    case 'SET_PARAMS': {
      const steps = state.steps.map(s =>
        s.stepId === action.stepId ? { ...s, params: { ...s.params, ...action.params } } : s
      );
      return { ...state, steps };
    }

    case 'SET_SELECTION': {
      const steps = state.steps.map(s =>
        s.stepId === action.stepId ? { ...s, selectedOutputIds: action.outputIds } : s
      );
      return { ...state, steps };
    }

    case 'CONFIRM': {
      const steps = state.steps.map((s, i) =>
        i === action.stepIndex ? { ...s, status: 'completed' as StepStatus } : s
      );
      const nextIndex = action.stepIndex + 1;
      const allDone = nextIndex >= steps.length;
      return {
        ...state,
        status: allDone ? 'completed' : 'running',
        activeStepIndex: allDone ? action.stepIndex : nextIndex,
        steps,
        completedAt: allDone ? Date.now() : undefined,
      };
    }

    case 'SKIP': {
      const steps = state.steps.map((s, i) =>
        i === action.stepIndex ? { ...s, status: 'skipped' as StepStatus, completedAt: Date.now() } : s
      );
      const nextIndex = action.stepIndex + 1;
      const allDone = nextIndex >= steps.length;
      return {
        ...state,
        status: allDone ? 'completed' : 'running',
        activeStepIndex: allDone ? action.stepIndex : nextIndex,
        steps,
        completedAt: allDone ? Date.now() : undefined,
      };
    }

    case 'CANCEL':
      return { ...state, status: 'cancelled', completedAt: Date.now() };

    case 'RESET':
      return INITIAL_STATE;

    default:
      return state;
  }
}

// ---- Callbacks interface (decoupled from store) ----

export interface PipelineCallbacks {
  onPdbUpdate?: (pdbContent: string) => void;
  onJobCreated?: (job: { id: string; type: string; status: string; createdAt: string }) => void;
}

// ---- Hook ----

export interface UsePipelineReturn {
  definition: PipelineDefinition | null;
  runtime: PipelineRuntimeState;
  activeStep: PipelineStepDefinition | null;
  activeStepState: StepRuntimeState | null;
  isPaused: boolean;
  isExecuting: boolean;

  initialize: (definition: PipelineDefinition, initialParams: Record<string, unknown>) => void;
  /** FIX #14: Restore pipeline state from session storage. Returns true if restored. */
  restore: (definition: PipelineDefinition) => boolean;
  runNextStep: () => Promise<void>;
  retryStep: (newParams?: Record<string, unknown>) => Promise<void>;
  skipStep: () => void;
  updateStepParams: (stepId: string, params: Record<string, unknown>) => void;
  setSelectedOutputs: (stepId: string, outputIds: string[]) => void;
  confirmAndContinue: () => void;
  cancel: () => void;
  reset: () => void;
}

export function usePipeline(callbacks?: PipelineCallbacks): UsePipelineReturn {
  const [runtime, dispatch] = useReducer(pipelineReducer, INITIAL_STATE);
  const definitionRef = useRef<PipelineDefinition | null>(null);
  const abortRef = useRef<AbortController | null>(null);
  // FIX #2: Use ref to always read fresh state (avoids stale closures)
  const runtimeRef = useRef(runtime);
  runtimeRef.current = runtime;
  // FIX #17: Guard against duplicate step execution
  const executingRef = useRef(false);
  // Store callbacks in ref to avoid re-creating executeStep
  const callbacksRef = useRef(callbacks);
  callbacksRef.current = callbacks;

  const definition = definitionRef.current;
  const activeStep = definition?.steps[runtime.activeStepIndex] ?? null;
  const activeStepState = runtime.steps[runtime.activeStepIndex] ?? null;
  const isPaused = runtime.status === 'paused';
  const isExecuting = runtime.status === 'running' && activeStepState?.status === 'running';

  const initialize = useCallback((def: PipelineDefinition, initialParams: Record<string, unknown>) => {
    definitionRef.current = def;
    abortRef.current = new AbortController();
    executingRef.current = false;

    const steps: StepRuntimeState[] = def.steps.map(step => ({
      stepId: step.id,
      status: 'pending',
      progress: 0,
      selectedOutputIds: [],
      params: { ...step.defaultParams },
    }));

    dispatch({
      type: 'INIT',
      pipelineId: def.id,
      sessionId: `${def.id}-${Date.now()}`,
      steps,
      initialParams,
    });
  }, []);

  const executeStep = useCallback(async (stepIndex: number) => {
    const def = definitionRef.current;
    if (!def) return;
    // FIX #17: Prevent concurrent execution
    if (executingRef.current) return;
    executingRef.current = true;

    const stepDef = def.steps[stepIndex];
    // FIX #2: Read fresh state from ref instead of closure
    const currentRuntime = runtimeRef.current;
    const stepState = currentRuntime.steps[stepIndex];
    if (!stepDef || !stepState) {
      executingRef.current = false;
      return;
    }

    dispatch({ type: 'STEP_START', stepIndex });

    // Build previous results map from ALL prior steps
    const previousResults: Record<string, StepResult> = {};
    for (let i = 0; i < stepIndex; i++) {
      const prev = currentRuntime.steps[i];
      if (prev.result) {
        previousResults[prev.stepId] = prev.result;
      }
    }

    // FIX #5: Aggregate selected items from ALL previous steps, not just immediate predecessor
    const selectedItems: string[] = [];
    for (let i = 0; i < stepIndex; i++) {
      const prev = currentRuntime.steps[i];
      if (prev.selectedOutputIds.length > 0) {
        selectedItems.push(...prev.selectedOutputIds);
      }
    }

    const ctx: StepExecutionContext = {
      previousResults,
      selectedItems,
      // FIX #3: Read fresh params from ref (supports retry with new params)
      params: runtimeRef.current.steps[stepIndex]?.params ?? stepState.params,
      initialParams: currentRuntime.initialParams,
      abortSignal: abortRef.current?.signal ?? new AbortController().signal,
      onProgress: (percent, message) => {
        dispatch({ type: 'STEP_PROGRESS', stepIndex, progress: percent, message });
      },
      onJobCreated: (jobId, type) => {
        callbacksRef.current?.onJobCreated?.({
          id: jobId,
          type,
          status: 'running',
          createdAt: new Date().toISOString(),
        });
      },
    };

    try {
      const result = await stepDef.execute(ctx);

      // Check if cancelled during execution
      if (abortRef.current?.signal.aborted) {
        executingRef.current = false;
        return;
      }

      // FIX #12: Notify viewer via callback instead of direct store access
      if (result.pdbOutputs && result.pdbOutputs.length > 0) {
        callbacksRef.current?.onPdbUpdate?.(result.pdbOutputs[0].pdbContent);
      }

      dispatch({ type: 'STEP_COMPLETE', stepIndex, result });

      if (stepDef.requiresReview) {
        dispatch({ type: 'STEP_PAUSE', stepIndex });

        // Auto-select outputs when supportsSelection is true
        if (stepDef.supportsSelection) {
          if (result.sequences && result.sequences.length > 0) {
            // For sequence steps: auto-select only stable sequences (II < 40)
            const stableIds = result.sequences
              .filter((seq: SequenceOutput) => {
                const pp = computeProtParam(seq.sequence);
                return pp.isStable;
              })
              .map((seq: SequenceOutput) => seq.id);
            // Fall back to all if none are stable
            const autoIds = stableIds.length > 0 ? stableIds : result.sequences.map((s: SequenceOutput) => s.id);
            dispatch({ type: 'SET_SELECTION', stepId: stepDef.id, outputIds: autoIds });
          } else if (result.pdbOutputs && result.pdbOutputs.length > 0) {
            // For backbone generation (rfd3_nl): auto-select ALL designs since there's no pre-filtering
            // For other steps (scaffold search): auto-select only the first (best) candidate
            const isBackboneStep = stepDef.id === 'rfd3_nl';
            const autoIds = isBackboneStep
              ? result.pdbOutputs.map(p => p.id)
              : [result.pdbOutputs[0].id];
            dispatch({ type: 'SET_SELECTION', stepId: stepDef.id, outputIds: autoIds });
          }
        }
      }
    } catch (err) {
      // FIX #12: Handle both abort signal and AbortError from shared-steps
      if (abortRef.current?.signal.aborted || (err instanceof DOMException && err.name === 'AbortError')) {
        dispatch({ type: 'CANCEL' });
        executingRef.current = false;
        return;
      }
      const errorMsg = err instanceof Error ? err.message : 'Unknown error';
      dispatch({ type: 'STEP_FAIL', stepIndex, error: errorMsg });
    } finally {
      executingRef.current = false;
    }
  }, []); // No deps needed — reads from refs

  const runNextStep = useCallback(async () => {
    const def = definitionRef.current;
    if (!def || executingRef.current) return;

    // FIX #2: Read from ref
    const current = runtimeRef.current;
    let nextIndex = current.activeStepIndex;
    while (nextIndex < current.steps.length) {
      const step = current.steps[nextIndex];
      if (step.status === 'pending' || step.status === 'ready') {
        await executeStep(nextIndex);
        return;
      }
      nextIndex++;
    }
  }, [executeStep]);

  const retryStep = useCallback(async (newParams?: Record<string, unknown>) => {
    // FIX #3: Apply params synchronously via ref read, then execute
    const current = runtimeRef.current;
    const stepIndex = current.activeStepIndex;
    if (newParams) {
      const stepId = current.steps[stepIndex]?.stepId;
      if (stepId) {
        dispatch({ type: 'SET_PARAMS', stepId, params: newParams });
      }
    }
    // Use requestAnimationFrame to ensure params are applied before execution
    await new Promise(resolve => requestAnimationFrame(resolve));
    await executeStep(stepIndex);
  }, [executeStep]);

  const skipStep = useCallback(() => {
    const def = definitionRef.current;
    const stepIndex = runtimeRef.current.activeStepIndex;
    if (def?.steps[stepIndex]?.optional) {
      dispatch({ type: 'SKIP', stepIndex });
    }
  }, []);

  const updateStepParams = useCallback((stepId: string, params: Record<string, unknown>) => {
    dispatch({ type: 'SET_PARAMS', stepId, params });
  }, []);

  const setSelectedOutputs = useCallback((stepId: string, outputIds: string[]) => {
    dispatch({ type: 'SET_SELECTION', stepId, outputIds });
  }, []);

  const confirmAndContinue = useCallback(() => {
    const stepIndex = runtimeRef.current.activeStepIndex;
    dispatch({ type: 'CONFIRM', stepIndex });
  }, []);

  const cancel = useCallback(() => {
    abortRef.current?.abort();
    executingRef.current = false;
    dispatch({ type: 'CANCEL' });
  }, []);

  const resetPipeline = useCallback(() => {
    abortRef.current?.abort();
    executingRef.current = false;
    definitionRef.current = null;
    dispatch({ type: 'RESET' });
    clearPipelineState();
  }, []);

  // FIX #14: Persist state at stable points (paused, completed, failed)
  useEffect(() => {
    savePipelineState(runtime);
  }, [runtime.status, runtime.activeStepIndex]);

  // FIX #14: Restore from session storage — returns saved state if available
  const restore = useCallback((def: PipelineDefinition): boolean => {
    const saved = loadPipelineState();
    if (!saved || saved.pipelineId !== def.id) return false;

    // Only restore paused or completed states
    if (saved.status !== 'paused' && saved.status !== 'completed' && saved.status !== 'failed') {
      clearPipelineState();
      return false;
    }

    definitionRef.current = def;
    abortRef.current = new AbortController();
    executingRef.current = false;

    // Replay the saved state via INIT + individual step updates
    dispatch({
      type: 'INIT',
      pipelineId: saved.pipelineId,
      sessionId: saved.sessionId,
      steps: saved.steps,
      initialParams: saved.initialParams,
    });

    // The INIT sets status to 'idle' — we need to restore the actual status.
    // Dispatch step-level actions to reconstruct the state.
    for (let i = 0; i < saved.steps.length; i++) {
      const step = saved.steps[i];
      if (step.status === 'completed' && step.result) {
        dispatch({ type: 'STEP_COMPLETE', stepIndex: i, result: step.result });
      } else if (step.status === 'paused' && step.result) {
        dispatch({ type: 'STEP_COMPLETE', stepIndex: i, result: step.result });
        dispatch({ type: 'STEP_PAUSE', stepIndex: i });
      } else if (step.status === 'failed') {
        dispatch({ type: 'STEP_FAIL', stepIndex: i, error: step.error ?? 'Unknown error' });
      } else if (step.status === 'skipped') {
        dispatch({ type: 'SKIP', stepIndex: i });
      }
      // Restore selections
      if (step.selectedOutputIds.length > 0) {
        dispatch({ type: 'SET_SELECTION', stepId: step.stepId, outputIds: step.selectedOutputIds });
      }
      // Restore params
      if (Object.keys(step.params).length > 0) {
        dispatch({ type: 'SET_PARAMS', stepId: step.stepId, params: step.params });
      }
    }

    return true;
  }, []);

  return {
    definition,
    runtime,
    activeStep,
    activeStepState,
    isPaused,
    isExecuting,
    initialize,
    restore,
    runNextStep,
    retryStep,
    skipStep,
    updateStepParams,
    setSelectedOutputs,
    confirmAndContinue,
    cancel,
    reset: resetPipeline,
  };
}
