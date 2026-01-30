'use client';

import { Check, X, Loader2, Pause } from 'lucide-react';
import { cn } from '@/lib/utils';
import type { PipelineStepDefinition, StepRuntimeState } from '@/lib/pipeline-types';

interface PipelineStepperProps {
  steps: PipelineStepDefinition[];
  stepStates: StepRuntimeState[];
  activeStepIndex: number;
  /** Index of the step whose detail card is currently shown */
  selectedStepIndex?: number | null;
  /** Called when a step icon is clicked */
  onStepClick?: (index: number) => void;
}

function getStatusIcon(status: StepRuntimeState['status']) {
  switch (status) {
    case 'completed':
      return <Check className="h-5 w-5" />;
    case 'running':
      return <Loader2 className="h-5 w-5 animate-spin" />;
    case 'paused':
      return <Pause className="h-5 w-5" />;
    case 'failed':
      return <X className="h-5 w-5" />;
    default:
      return null;
  }
}

function getStatusStyles(status: StepRuntimeState['status']) {
  if (status === 'completed') {
    return { iconBg: 'bg-primary/10 text-primary', line: 'bg-primary' };
  }
  if (status === 'running' || status === 'paused') {
    return { iconBg: 'bg-primary/10 text-primary ring-2 ring-primary/30', line: 'bg-border' };
  }
  if (status === 'failed') {
    return { iconBg: 'bg-destructive/10 text-destructive', line: 'bg-destructive/30' };
  }
  if (status === 'skipped') {
    return { iconBg: 'bg-muted text-muted-foreground', line: 'bg-border' };
  }
  return { iconBg: 'bg-muted text-muted-foreground', line: 'bg-border' };
}

export function PipelineStepper({
  steps,
  stepStates,
  activeStepIndex,
  selectedStepIndex,
  onStepClick,
}: PipelineStepperProps) {
  return (
    <div className="flex items-center w-full">
      {steps.map((step, index) => {
        const state = stepStates[index];
        const isActive = index === activeStepIndex;
        const isSelected = selectedStepIndex === index;
        const styles = getStatusStyles(state?.status ?? 'pending');
        const Icon = step.icon;
        const statusIcon = getStatusIcon(state?.status ?? 'pending');

        return (
          <div key={step.id} className="flex items-center flex-1 last:flex-none">
            {/* Step circle + label */}
            <div className="flex flex-col items-center gap-1 group relative">
              <button
                type="button"
                onClick={() => onStepClick?.(index)}
                title={step.name}
                className={cn(
                  'flex items-center justify-center w-11 h-11 rounded-full transition-all shrink-0',
                  styles.iconBg,
                  'cursor-pointer hover:scale-110 hover:shadow-md',
                  isSelected && 'ring-2 ring-offset-2 ring-primary',
                )}
              >
                {statusIcon ?? <Icon className="h-5 w-5" />}
              </button>
              <span
                className={cn(
                  'text-[10px] leading-tight text-center whitespace-nowrap',
                  isActive ? 'font-semibold text-foreground' : 'text-muted-foreground',
                  state?.status === 'failed' && 'text-destructive',
                  state?.status === 'completed' && 'text-primary',
                  isSelected && 'font-semibold',
                )}
              >
                {step.name.length > 12 ? step.name.slice(0, 11) + '…' : step.name}
              </span>
              {/* Progress bar under active running step */}
              {state?.status === 'running' && state.progress > 0 && (
                <div className="w-11 h-1 bg-muted rounded-full">
                  <div
                    className="h-full bg-primary rounded-full transition-all duration-300"
                    style={{ width: `${state.progress}%` }}
                  />
                </div>
              )}
            </div>

            {/* Connector line */}
            {index < steps.length - 1 && (
              <div className={cn('flex-1 h-px mx-1', styles.line)} />
            )}
          </div>
        );
      })}
    </div>
  );
}
