'use client';

import { useEffect, useCallback, useRef, useState, Component, type ReactNode } from 'react';
import { Card, CardContent, CardHeader } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { Separator } from '@/components/ui/separator';
// ScrollArea removed — parent panel handles scrolling
import { XCircle, AlertTriangle, Zap, PartyPopper, Trophy, TrendingUp } from 'lucide-react';
import { cn } from '@/lib/utils';
import { usePipeline, type PipelineCallbacks } from '@/hooks/usePipeline';
import { useStore } from '@/lib/store';
import type { PipelineDefinition, StepResult } from '@/lib/pipeline-types';
import { PipelineStepper } from './PipelineStepper';
import { StepCard } from './StepCard';
import { pipelineCompleteMessages, getRandomMessage, trackPipelineRun } from '@/lib/step-personality';
import { Confetti } from './Confetti';

// ---- Error Boundary (FIX #8) ----

interface ErrorBoundaryProps {
  children: ReactNode;
  fallback: (error: Error, reset: () => void) => ReactNode;
}

interface ErrorBoundaryState {
  error: Error | null;
}

class PipelineErrorBoundary extends Component<ErrorBoundaryProps, ErrorBoundaryState> {
  state: ErrorBoundaryState = { error: null };

  static getDerivedStateFromError(error: Error): ErrorBoundaryState {
    return { error };
  }

  reset = () => this.setState({ error: null });

  render() {
    if (this.state.error) {
      return this.props.fallback(this.state.error, this.reset);
    }
    return this.props.children;
  }
}

// ---- Cancel Confirmation Dialog (FIX #16) ----

function CancelConfirmation({ onConfirm, onDismiss }: { onConfirm: () => void; onDismiss: () => void }) {
  return (
    <div className="flex items-center gap-2 p-2 bg-destructive/5 border border-destructive/20 rounded-md">
      <AlertTriangle className="h-3.5 w-3.5 text-destructive shrink-0" />
      <span className="text-xs text-destructive">Cancel pipeline? Running jobs will continue on the server but results will be lost.</span>
      <div className="flex gap-1 shrink-0 ml-auto">
        <Button variant="destructive" size="sm" className="h-6 px-2 text-xs" onClick={onConfirm}>
          Cancel Pipeline
        </Button>
        <Button variant="ghost" size="sm" className="h-6 px-2 text-xs" onClick={onDismiss}>
          Keep Running
        </Button>
      </div>
    </div>
  );
}

// ---- Main Component ----

interface PipelineRunnerProps {
  definition: PipelineDefinition;
  initialParams: Record<string, unknown>;
  onStepComplete?: (stepId: string, result: StepResult) => void;
  onPipelineComplete?: (results: Record<string, StepResult>) => void;
  onPipelineFailed?: (stepId: string, error: string) => void;
  onDesignSelected?: (pdbContent: string) => void;
  onCancel?: () => void;
}

function PipelineRunnerInner({
  definition,
  initialParams,
  onStepComplete,
  onPipelineComplete,
  onPipelineFailed,
  onDesignSelected,
  onCancel,
}: PipelineRunnerProps) {
  const { setSelectedPdb, addJob, updateJob } = useStore();
  const [showCancelConfirm, setShowCancelConfirm] = useState(false);
  const [selectedStepIndex, setSelectedStepIndex] = useState<number | null>(null);
  const [autoRun, setAutoRun] = useState(false);

  // FIX #12: Pass store callbacks via the decoupled interface
  const pipelineCallbacks: PipelineCallbacks = {
    onPdbUpdate: (pdb) => setSelectedPdb(pdb),
    onJobCreated: (job) => addJob(job as any),
    onJobCompleted: (jobId, status) => {
      updateJob(jobId, {
        status,
        completedAt: new Date().toISOString(),
      });
    },
  };

  const pipeline = usePipeline(pipelineCallbacks);

  // FIX #6: Track which step completions have already been notified
  const notifiedStepsRef = useRef<Set<string>>(new Set());
  const notifiedCompleteRef = useRef(false);
  const notifiedFailedRef = useRef(false);
  const restoringRef = useRef(false);

  // Try to restore from session storage first, otherwise initialize fresh
  // restore() checks user_prompt to ensure we don't restore old state for new designs
  useEffect(() => {
    restoringRef.current = true;
    notifiedStepsRef.current = new Set();
    notifiedCompleteRef.current = false;
    notifiedFailedRef.current = false;

    // Try to restore existing pipeline state (same user_prompt = same design session)
    // restore() returns false if: no saved state, different pipeline, or different user_prompt
    const restored = pipeline.restore(definition, initialParams);

    if (!restored) {
      // No valid restore - initialize fresh
      pipeline.initialize(definition, initialParams);
    }

    // Allow one render cycle for restore dispatches to settle before enabling notifications
    requestAnimationFrame(() => {
      // After restore settles, mark all currently completed/paused steps as notified
      for (const s of pipeline.runtime.steps) {
        if (s.result && (s.status === 'completed' || s.status === 'paused' || s.status === 'skipped')) {
          notifiedStepsRef.current.add(s.stepId);
        }
      }
      restoringRef.current = false;
    });
    // Re-run when definition changes OR when user_prompt changes (new design request)
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [definition.id, initialParams.user_prompt]);

  // FIX #7: Cancel on unmount to prevent memory leaks
  useEffect(() => {
    return () => {
      pipeline.cancel();
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  // Auto-run the first step after initialization
  useEffect(() => {
    if (pipeline.runtime.status === 'idle' && pipeline.runtime.steps.length > 0) {
      pipeline.runNextStep();
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [pipeline.runtime.status, pipeline.runtime.steps.length]);

  // Auto-run next step after confirm (when status transitions from paused to running)
  useEffect(() => {
    if (pipeline.runtime.status === 'running') {
      const activeState = pipeline.runtime.steps[pipeline.runtime.activeStepIndex];
      if (activeState && activeState.status === 'pending') {
        pipeline.runNextStep();
      }
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [pipeline.runtime.status, pipeline.runtime.activeStepIndex]);

  // Auto-confirm paused steps when Quick Run is active
  // Exception: don't auto-confirm when all sweep designs failed — user must review
  useEffect(() => {
    if (autoRun && pipeline.runtime.status === 'paused') {
      const activeState = pipeline.runtime.steps[pipeline.runtime.activeStepIndex];
      if (activeState?.result?.data?.sweep_all_failed) {
        return; // Stop auto-run — user must review failure and decide to retry or continue
      }
      pipeline.confirmAndContinue();
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [autoRun, pipeline.runtime.status]);

  // FIX #6: Notify on step completion — only once per step, skip during restore
  useEffect(() => {
    if (restoringRef.current) return;
    for (const stepState of pipeline.runtime.steps) {
      if (stepState.status === 'paused' && stepState.result && !notifiedStepsRef.current.has(stepState.stepId)) {
        notifiedStepsRef.current.add(stepState.stepId);
        onStepComplete?.(stepState.stepId, stepState.result);
      }
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [pipeline.runtime.steps]);

  // Notify on pipeline completion — only once
  useEffect(() => {
    if (pipeline.runtime.status === 'completed' && !notifiedCompleteRef.current) {
      notifiedCompleteRef.current = true;
      const results: Record<string, StepResult> = {};
      for (const stepState of pipeline.runtime.steps) {
        if (stepState.result) {
          results[stepState.stepId] = stepState.result;
        }
      }
      onPipelineComplete?.(results);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [pipeline.runtime.status]);

  // Notify on failure — only once
  useEffect(() => {
    if (pipeline.runtime.status === 'failed' && !notifiedFailedRef.current) {
      notifiedFailedRef.current = true;
      const failedStep = pipeline.runtime.steps.find(s => s.status === 'failed');
      if (failedStep) {
        onPipelineFailed?.(failedStep.stepId, failedStep.error ?? 'Unknown error');
      }
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [pipeline.runtime.status]);

  // Auto-select the active step when it changes (running/paused/failed)
  useEffect(() => {
    const activeState = pipeline.runtime.steps[pipeline.runtime.activeStepIndex];
    if (activeState && (activeState.status === 'running' || activeState.status === 'paused' || activeState.status === 'failed')) {
      setSelectedStepIndex(pipeline.runtime.activeStepIndex);
    }
  }, [pipeline.runtime.activeStepIndex, pipeline.runtime.steps]);

  const handleStepClick = useCallback((index: number) => {
    setSelectedStepIndex(prev => prev === index ? null : index);
  }, []);

  // FIX #16: Cancel with confirmation
  const handleCancelClick = useCallback(() => {
    setShowCancelConfirm(true);
  }, []);

  const handleCancelConfirm = useCallback(() => {
    setShowCancelConfirm(false);
    pipeline.cancel();
    onCancel?.();
  }, [pipeline, onCancel]);

  const handleCancelDismiss = useCallback(() => {
    setShowCancelConfirm(false);
  }, []);

  const isTerminal = pipeline.runtime.status === 'completed' || pipeline.runtime.status === 'cancelled' || pipeline.runtime.status === 'failed';

  // Completion stats
  const completionMessage = useRef<string | null>(null);
  const milestoneMessage = useRef<string | null>(null);
  const [showConfetti, setShowConfetti] = useState(false);

  useEffect(() => {
    if (pipeline.runtime.status === 'completed' && !completionMessage.current) {
      completionMessage.current = getRandomMessage(pipelineCompleteMessages);
      milestoneMessage.current = trackPipelineRun();
      setShowConfetti(true);
      const timer = setTimeout(() => setShowConfetti(false), 2500);
      return () => clearTimeout(timer);
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [pipeline.runtime.status]);

  // Compute stats for the completion card
  const completionStats = pipeline.runtime.status === 'completed' ? (() => {
    const steps = pipeline.runtime.steps;
    const completedCount = steps.filter(s => s.status === 'completed').length;
    const skippedCount = steps.filter(s => s.status === 'skipped').length;
    return { completedCount, skippedCount, totalCount: steps.length };
  })() : null;

  return (
    <Card className="border-border relative overflow-hidden">
      {showConfetti && <Confetti />}
      <CardHeader className="px-4 py-3 space-y-3">
        {/* Controls row — no duplicate name/description */}
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-2">
            {pipeline.runtime.status === 'completed' && (
              <Badge variant="secondary" className="text-xs gap-1">
                <PartyPopper className="h-3 w-3" />
                Complete
              </Badge>
            )}
            {pipeline.runtime.status === 'failed' && (
              <Badge variant="destructive" className="text-xs">Failed</Badge>
            )}
            {pipeline.runtime.status === 'cancelled' && (
              <Badge variant="secondary" className="text-xs">Cancelled</Badge>
            )}
          </div>
          <div className="flex items-center gap-2">
            {!isTerminal && (
              <div className="flex items-center gap-1">
                <Button
                  variant={autoRun ? 'default' : 'ghost'}
                  size="sm"
                  onClick={() => setAutoRun(prev => !prev)}
                  className={cn(
                    'h-7 px-2 text-xs',
                    autoRun
                      ? 'bg-amber-500 hover:bg-amber-600 text-white'
                      : 'text-muted-foreground hover:text-amber-600',
                  )}
                  title={autoRun ? 'Auto-proceeding through all steps' : 'Run all steps without pausing for review'}
                >
                  <Zap className="h-3.5 w-3.5 mr-1" />
                  {autoRun ? 'Auto-Proceed On' : 'Auto-Proceed'}
                </Button>
                <Button
                  variant="ghost"
                  size="sm"
                  onClick={handleCancelClick}
                  className="h-7 px-2 text-xs text-muted-foreground hover:text-destructive"
                >
                  <XCircle className="h-3.5 w-3.5 mr-1" />
                  Cancel
                </Button>
              </div>
            )}
          </div>
        </div>

        {showCancelConfirm && (
          <CancelConfirmation onConfirm={handleCancelConfirm} onDismiss={handleCancelDismiss} />
        )}

        <PipelineStepper
          steps={definition.steps}
          stepStates={pipeline.runtime.steps}
          activeStepIndex={pipeline.runtime.activeStepIndex}
          selectedStepIndex={selectedStepIndex}
          onStepClick={handleStepClick}
        />
      </CardHeader>

      <Separator />

      <CardContent className="p-3">
          {/* Completion celebration card */}
          {pipeline.runtime.status === 'completed' && completionStats && (
            <div className="mb-3 rounded-lg bg-success/5 border border-success/20 p-4 animate-step-enter">
              <div className="flex items-center gap-2 mb-2">
                <Trophy className="h-4 w-4 text-amber-500" />
                <span className="text-sm font-medium text-foreground">Design Complete</span>
              </div>
              <div className="flex items-center gap-4 text-xs text-muted-foreground mb-2">
                <span className="flex items-center gap-1">
                  <TrendingUp className="h-3 w-3 text-success" />
                  {completionStats.completedCount} steps completed
                </span>
                {completionStats.skippedCount > 0 && (
                  <span>{completionStats.skippedCount} skipped</span>
                )}
              </div>
              {completionMessage.current && (
                <p className="text-xs italic text-muted-foreground/80">
                  {completionMessage.current}
                </p>
              )}
              {milestoneMessage.current && (
                <p className="text-xs font-medium text-amber-600 mt-1.5">
                  {milestoneMessage.current}
                </p>
              )}
            </div>
          )}

          {selectedStepIndex !== null && (() => {
            const step = definition.steps[selectedStepIndex];
            const state = pipeline.runtime.steps[selectedStepIndex];
            if (!step || !state) return null;

            const isActive = selectedStepIndex === pipeline.runtime.activeStepIndex;
            const nextStep = definition.steps[selectedStepIndex + 1];
            const nextState = pipeline.runtime.steps[selectedStepIndex + 1];

            // Allow going back to completed/skipped steps that aren't the current active step
            // Also allow when pipeline is paused (user is reviewing) or failed
            const canGoBackToStep = (state.status === 'completed' || state.status === 'skipped')
              && !isActive
              && !pipeline.isExecuting;

            return (
              <StepCard
                key={step.id}
                step={step}
                state={state}
                isActive={isActive}
                allowExpand
                canGoBack={canGoBackToStep}
                nextStepSchema={
                  // Show next step's params for review — but skip optional steps
                  // (e.g. scout filter thresholds shouldn't appear in backbone gen
                  // since scout is auto-skipped for metal designs).
                  isActive && state.status === 'paused' && nextStep && !nextStep.optional
                    ? nextStep.parameterSchema
                    : undefined
                }
                nextStepParams={nextState?.params}
                onNextStepParamsChange={
                  nextStep
                    ? (params) => pipeline.updateStepParams(nextStep.id, params)
                    : undefined
                }
                onConfirm={() => pipeline.confirmAndContinue()}
                onRetry={() => pipeline.retryStep()}
                onSkip={() => pipeline.skipStep()}
                onGoBack={() => pipeline.goBack(selectedStepIndex)}
                onSelectionChange={(ids) => pipeline.setSelectedOutputs(step.id, ids)}
                onViewDesign={onDesignSelected}
              />
            );
          })()}
          {selectedStepIndex === null && (
            <p className="text-sm text-muted-foreground text-center py-4">
              Click a step icon above to view its details
            </p>
          )}
      </CardContent>
    </Card>
  );
}

// Wrap in error boundary
export function PipelineRunner(props: PipelineRunnerProps) {
  return (
    <PipelineErrorBoundary
      fallback={(error, reset) => (
        <Card className="border-destructive/30">
          <CardContent className="p-4">
            <div className="flex items-start gap-3">
              <AlertTriangle className="h-5 w-5 text-destructive shrink-0 mt-0.5" />
              <div className="space-y-2">
                <p className="text-sm font-medium text-destructive">Pipeline Error</p>
                <p className="text-xs text-muted-foreground">{error.message}</p>
                <div className="flex gap-2">
                  <Button variant="outline" size="sm" className="h-7 text-xs" onClick={reset}>
                    Try Again
                  </Button>
                  {props.onCancel && (
                    <Button variant="ghost" size="sm" className="h-7 text-xs" onClick={props.onCancel}>
                      Cancel
                    </Button>
                  )}
                </div>
              </div>
            </div>
          </CardContent>
        </Card>
      )}
    >
      <PipelineRunnerInner {...props} />
    </PipelineErrorBoundary>
  );
}
