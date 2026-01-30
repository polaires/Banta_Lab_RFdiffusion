'use client';

import { useState, useCallback, useRef } from 'react';
import { Check, X, Loader2, Pause, ChevronDown, ChevronRight, RotateCcw, SkipForward, ArrowRight, Info } from 'lucide-react';
import { Button } from '@/components/ui/button';
import { Card, CardContent } from '@/components/ui/card';
import { Progress } from '@/components/ui/progress';
import { Badge } from '@/components/ui/badge';
import { Separator } from '@/components/ui/separator';
import {
  Collapsible,
  CollapsibleContent,
  CollapsibleTrigger,
} from '@/components/ui/collapsible';
import { cn } from '@/lib/utils';
import type { PipelineStepDefinition, StepRuntimeState } from '@/lib/pipeline-types';
import { StepResultPreview } from './StepResultPreview';
import { StepParameterEditor } from './StepParameterEditor';
import { DesignSelector } from './DesignSelector';

interface StepCardProps {
  step: PipelineStepDefinition;
  state: StepRuntimeState;
  isActive: boolean;
  /** Allow expanding completed steps to review results (FIX #18) */
  allowExpand?: boolean;
  /** Next step's parameter schema (to show editor for upcoming step) */
  nextStepSchema?: PipelineStepDefinition['parameterSchema'];
  nextStepParams?: Record<string, unknown>;
  onNextStepParamsChange?: (params: Record<string, unknown>) => void;
  onConfirm: () => void;
  onRetry: () => void;
  onSkip: () => void;
  onSelectionChange: (ids: string[]) => void;
  onViewDesign?: (pdbContent: string) => void;
}

function getStatusBadge(status: StepRuntimeState['status']) {
  switch (status) {
    case 'completed':
      return <Badge variant="secondary" className="text-[10px] h-5"><Check className="h-3 w-3 mr-1" />Done</Badge>;
    case 'running':
      return <Badge variant="default" className="text-[10px] h-5"><Loader2 className="h-3 w-3 mr-1 animate-spin" />Running</Badge>;
    case 'paused':
      return <Badge variant="outline" className="text-[10px] h-5 border-primary text-primary"><Pause className="h-3 w-3 mr-1" />Review</Badge>;
    case 'failed':
      return <Badge variant="destructive" className="text-[10px] h-5"><X className="h-3 w-3 mr-1" />Failed</Badge>;
    case 'skipped':
      return <Badge variant="secondary" className="text-[10px] h-5 line-through">Skipped</Badge>;
    default:
      return <Badge variant="secondary" className="text-[10px] h-5 text-muted-foreground">Pending</Badge>;
  }
}

// FIX #30: Map common errors to user-friendly messages
function getActionableError(error: string): { message: string; suggestion?: string } {
  const lower = error.toLowerCase();

  if (lower.includes('no backbone') || lower.includes('no structures')) {
    return {
      message: error,
      suggestion: 'Try adjusting the design parameters (more designs, different timesteps) and retry.',
    };
  }
  if (lower.includes('no sequences')) {
    return {
      message: error,
      suggestion: 'Ensure the previous step produced structures, then retry.',
    };
  }
  if (lower.includes('backend') || lower.includes('network') || lower.includes('fetch')) {
    return {
      message: error,
      suggestion: 'Check your backend connection and try again.',
    };
  }
  if (lower.includes('timeout') || lower.includes('timed out')) {
    return {
      message: error,
      suggestion: 'The job may still be running on the server. Try retrying or check the Jobs tab.',
    };
  }
  if (lower.includes('pdb') || lower.includes('structure')) {
    return {
      message: error,
      suggestion: 'The input structure may be invalid. Try a different PDB or re-upload.',
    };
  }

  return { message: error };
}

export function StepCard({
  step,
  state,
  isActive,
  allowExpand,
  nextStepSchema,
  nextStepParams,
  onNextStepParamsChange,
  onConfirm,
  onRetry,
  onSkip,
  onSelectionChange,
  onViewDesign,
}: StepCardProps) {
  const Icon = step.icon;
  const isPaused = state.status === 'paused';
  const isFailed = state.status === 'failed';
  const isRunning = state.status === 'running';
  const isDone = state.status === 'completed' || state.status === 'skipped';

  // FIX #18: Allow toggling completed steps to review results
  const [manualOpen, setManualOpen] = useState<boolean | null>(null);
  const isExpandedByState = isActive && (isRunning || isPaused || isFailed);
  // When allowExpand is set, default to expanded (icon-toggled view)
  const isExpanded = manualOpen !== null ? manualOpen : (allowExpand || isExpandedByState);

  // FIX #17: Debounce action buttons to prevent rapid clicking
  const actionPendingRef = useRef(false);
  const debounceAction = useCallback((action: () => void) => {
    if (actionPendingRef.current) return;
    actionPendingRef.current = true;
    action();
    setTimeout(() => { actionPendingRef.current = false; }, 500);
  }, []);

  // FIX #4: Determine if continue should be disabled
  // Only require selection if step supports it AND has selectable outputs
  const hasSelectableOutputs = step.supportsSelection && (
    (state.result?.pdbOutputs && state.result.pdbOutputs.length > 0) ||
    (state.result?.sequences && state.result.sequences.length > 0)
  );
  const continueDisabled = hasSelectableOutputs && state.selectedOutputIds.length === 0;

  const handleToggle = () => {
    if (isDone && allowExpand) {
      setManualOpen(prev => prev === null ? true : !prev);
    }
  };

  return (
    <Collapsible open={isExpanded}>
      <Card className={cn(
        'transition-all',
        isActive && !isDone && 'border-primary/30',
        isFailed && 'border-destructive/30',
      )}>
        <CollapsibleTrigger asChild>
          <div
            className={cn(
              'flex items-center gap-3 px-3 py-2.5 cursor-pointer hover:bg-accent/50 transition-colors rounded-t-lg',
              isDone && !isExpanded && 'rounded-b-lg',
            )}
            onClick={handleToggle}
          >
            <div className={cn(
              'flex items-center justify-center w-7 h-7 rounded-full shrink-0',
              isRunning ? 'bg-primary/10 text-primary' :
              isPaused ? 'bg-primary/10 text-primary' :
              isFailed ? 'bg-destructive/10 text-destructive' :
              isDone ? 'bg-primary/10 text-primary' :
              'bg-muted text-muted-foreground'
            )}>
              {isRunning ? <Loader2 className="h-3.5 w-3.5 animate-spin" /> :
               isDone ? <Check className="h-3.5 w-3.5" /> :
               isFailed ? <X className="h-3.5 w-3.5" /> :
               isPaused ? <Pause className="h-3.5 w-3.5" /> :
               <Icon className="h-3.5 w-3.5" />}
            </div>

            <div className="flex-1 min-w-0">
              <div className="flex items-center gap-2">
                <span className={cn(
                  'text-sm font-medium',
                  isDone || isRunning || isPaused ? 'text-foreground' : 'text-muted-foreground'
                )}>
                  {step.name}
                </span>
                {getStatusBadge(state.status)}
              </div>
              {!isExpanded && state.result?.summary && (
                <p className="text-xs text-muted-foreground truncate mt-0.5">
                  {state.result.summary}
                </p>
              )}
            </div>

            {isExpanded ? (
              <ChevronDown className="h-4 w-4 text-muted-foreground shrink-0" />
            ) : (
              <ChevronRight className="h-4 w-4 text-muted-foreground shrink-0" />
            )}
          </div>
        </CollapsibleTrigger>

        <CollapsibleContent>
          <Separator />
          <CardContent className="p-3 space-y-3">
            {/* Running: progress bar */}
            {isRunning && (
              <div className="space-y-1.5">
                <div className="flex items-center justify-between">
                  <span className="text-xs text-muted-foreground">
                    {state.progressMessage ?? step.description}
                  </span>
                  <span className="text-xs font-mono text-muted-foreground">
                    {state.progress}%
                  </span>
                </div>
                <Progress value={state.progress} className="h-1.5" />
              </div>
            )}

            {/* FIX #30: Failed: actionable error message */}
            {isFailed && state.error && (() => {
              const { message, suggestion } = getActionableError(state.error);
              return (
                <div className="text-xs bg-destructive/5 rounded-md p-2.5 border border-destructive/20 space-y-1.5">
                  <p className="text-destructive">{message}</p>
                  {suggestion && (
                    <p className="text-muted-foreground flex items-start gap-1">
                      <Info className="h-3 w-3 shrink-0 mt-0.5" />
                      {suggestion}
                    </p>
                  )}
                </div>
              );
            })()}

            {/* Paused or completed with results: result preview */}
            {(isPaused || (isDone && isExpanded)) && state.result && (
              <>
                {step.ResultPreview ? (
                  <step.ResultPreview result={state.result} onSelectDesign={onViewDesign} />
                ) : (
                  <StepResultPreview result={state.result} onSelectDesign={onViewDesign} />
                )}

                {/* Selection grid — only on paused (active review) */}
                {isPaused && step.supportsSelection && (state.result.pdbOutputs || state.result.sequences) && (
                  <>
                    <Separator />
                    <div>
                      {/* FIX #20: Explain what selection does */}
                      <p className="text-xs font-medium text-foreground mb-1">
                        Select outputs to carry forward
                      </p>
                      <p className="text-[10px] text-muted-foreground mb-2">
                        Only selected items will be used by the next step. Unselected items are discarded.
                      </p>
                      <DesignSelector
                        pdbOutputs={state.result.pdbOutputs}
                        sequences={state.result.sequences}
                        selectedIds={state.selectedOutputIds}
                        onSelectionChange={onSelectionChange}
                        onViewDesign={onViewDesign}
                      />
                    </div>
                  </>
                )}

                {/* Next step parameter editor — only on paused */}
                {isPaused && nextStepSchema && nextStepSchema.length > 0 && nextStepParams && onNextStepParamsChange && (
                  <>
                    <Separator />
                    <div>
                      <p className="text-xs font-medium text-foreground mb-2">
                        Parameters for next step
                      </p>
                      <StepParameterEditor
                        schema={nextStepSchema}
                        params={nextStepParams}
                        onChange={onNextStepParamsChange}
                      />
                    </div>
                  </>
                )}
              </>
            )}

            {/* FIX #19: Empty results message */}
            {isPaused && state.result && !state.result.pdbOutputs?.length && !state.result.sequences?.length && !state.result.data && (
              <div className="text-xs text-muted-foreground bg-muted/50 rounded-md p-2.5 border border-border">
                This step completed but produced no structures or sequences. Check the summary above for details.
              </div>
            )}

            {/* Action buttons */}
            {(isPaused || isFailed) && (
              <div className="flex items-center gap-2 pt-1">
                {isPaused && (
                  <Button
                    size="sm"
                    onClick={() => debounceAction(onConfirm)}
                    className="h-8 text-xs"
                    disabled={continueDisabled}
                    title={continueDisabled ? 'Select at least one output to continue' : undefined}
                  >
                    <ArrowRight className="h-3.5 w-3.5 mr-1" />
                    Continue
                  </Button>
                )}
                {isFailed && (
                  <Button
                    variant="outline"
                    size="sm"
                    onClick={() => debounceAction(onRetry)}
                    className="h-8 text-xs"
                  >
                    <RotateCcw className="h-3.5 w-3.5 mr-1" />
                    Retry
                  </Button>
                )}
                {step.optional && (isPaused || isFailed) && (
                  <Button
                    variant="ghost"
                    size="sm"
                    onClick={() => debounceAction(onSkip)}
                    className="h-8 text-xs"
                  >
                    <SkipForward className="h-3.5 w-3.5 mr-1" />
                    Skip
                  </Button>
                )}
                {/* FIX #4: Show hint when selection is required */}
                {continueDisabled && (
                  <span className="text-[10px] text-muted-foreground">
                    Select at least one output above
                  </span>
                )}
              </div>
            )}
          </CardContent>
        </CollapsibleContent>
      </Card>
    </Collapsible>
  );
}
