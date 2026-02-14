'use client';

import { useEffect, useRef } from 'react';
import { Check, X, Pause } from 'lucide-react';
import { cn } from '@/lib/utils';
import type { PipelineStepDefinition, StepRuntimeState } from '@/lib/pipeline-types';
import { getStepPersonality } from '@/lib/step-personality';

interface PipelineStepperProps {
  steps: PipelineStepDefinition[];
  stepStates: StepRuntimeState[];
  activeStepIndex: number;
  /** Index of the step whose detail card is currently shown */
  selectedStepIndex?: number | null;
  /** Called when a step icon is clicked */
  onStepClick?: (index: number) => void;
}

/** Auto-step IDs that get smaller dots in the minimap */
const AUTO_STEP_IDS = new Set(['save_history_nl', 'check_lessons_nl']);

function getStatusStyles(status: StepRuntimeState['status']) {
  if (status === 'completed') {
    return { iconBg: 'bg-primary/10 text-primary', line: 'bg-primary' };
  }
  if (status === 'running' || status === 'paused') {
    return { iconBg: 'bg-primary/10 text-primary', line: 'bg-border' };
  }
  if (status === 'failed') {
    return { iconBg: 'bg-destructive/10 text-destructive animate-gentle-shake', line: 'bg-destructive/30' };
  }
  if (status === 'skipped') {
    return { iconBg: 'bg-muted text-muted-foreground', line: 'bg-border' };
  }
  return { iconBg: 'bg-muted/60 text-muted-foreground/50', line: 'bg-border' };
}

export function PipelineStepper({
  steps,
  stepStates,
  activeStepIndex,
  selectedStepIndex,
  onStepClick,
}: PipelineStepperProps) {
  const activeRef = useRef<HTMLButtonElement>(null);

  // Auto-scroll to keep active step visible
  useEffect(() => {
    if (activeRef.current) {
      activeRef.current.scrollIntoView({
        behavior: 'smooth',
        block: 'nearest',
        inline: 'center',
      });
    }
  }, [activeStepIndex]);

  return (
    <div className="flex items-center w-full overflow-x-auto overflow-y-hidden scrollbar-hide px-1">
      {steps.map((step, index) => {
        const state = stepStates[index];
        const status = state?.status ?? 'pending';
        const isActive = index === activeStepIndex;
        const isSelected = selectedStepIndex === index;
        const isRunning = status === 'running';
        const isCompleted = status === 'completed';
        const isFailed = status === 'failed';
        const isPaused = status === 'paused';
        const isAutoStep = AUTO_STEP_IDS.has(step.id);
        const styles = getStatusStyles(status);
        const personality = getStepPersonality(step.id);

        // Choose which icon to show
        const Icon = step.icon;
        const RunningIcon = personality?.runningIcon ?? Icon;
        const animClass = personality?.animationClass ?? 'animate-spin';

        // Node sizing
        const nodeSize = isAutoStep ? 'w-8 h-8' : 'w-11 h-11';
        const iconSize = isAutoStep ? 'h-3.5 w-3.5' : 'h-5 w-5';

        return (
          <div key={step.id} className="flex items-center flex-1 last:flex-none">
            {/* Step circle + label */}
            <div className="flex flex-col items-center gap-1 group relative">
              <button
                ref={isActive ? activeRef : undefined}
                type="button"
                onClick={() => onStepClick?.(index)}
                title={isRunning && personality ? personality.verbPhrase : step.name}
                className={cn(
                  'flex items-center justify-center rounded-full transition-all shrink-0 relative',
                  nodeSize,
                  styles.iconBg,
                  'cursor-pointer hover:scale-110 hover:shadow-md',
                  isSelected && 'ring-2 ring-offset-2 ring-primary',
                  isRunning && 'animate-active-glow',
                )}
              >
                {/* Running: personality icon + animation */}
                {isRunning && (
                  <div className={animClass}>
                    <RunningIcon className={iconSize} />
                  </div>
                )}
                {/* Completed: checkmark with bloom */}
                {isCompleted && (
                  <div className="animate-check-bloom">
                    <Check className={iconSize} />
                  </div>
                )}
                {/* Failed */}
                {isFailed && <X className={iconSize} />}
                {/* Paused */}
                {isPaused && <Pause className={iconSize} />}
                {/* Pending: default step icon */}
                {status === 'pending' && <Icon className={iconSize} />}
                {/* Skipped */}
                {status === 'skipped' && <Icon className={cn(iconSize, 'opacity-50')} />}

                {/* Progress ring (SVG) for running steps */}
                {isRunning && state.progress > 0 && !isAutoStep && (
                  <svg className="absolute inset-0" width="100%" height="100%" viewBox="0 0 44 44">
                    <circle
                      cx="22" cy="22" r="20"
                      fill="none"
                      stroke="currentColor"
                      strokeWidth="2"
                      strokeDasharray={2 * Math.PI * 20}
                      strokeDashoffset={2 * Math.PI * 20 * (1 - state.progress / 100)}
                      strokeLinecap="round"
                      className="text-primary/40 -rotate-90 origin-center transition-all duration-500"
                    />
                  </svg>
                )}
              </button>

              {/* Label */}
              {!isAutoStep && (
                <span
                  className={cn(
                    'text-[10px] leading-tight text-center max-w-[72px] truncate transition-colors',
                    isActive ? 'font-semibold text-foreground' : 'text-muted-foreground',
                    isFailed && 'text-destructive',
                    isCompleted && 'text-primary',
                    isSelected && 'font-semibold',
                  )}
                  title={step.name}
                >
                  {step.name}
                </span>
              )}

              {/* Caret indicator for selected step */}
              {isSelected && (
                <div className="absolute -bottom-2.5 left-1/2 -translate-x-1/2 w-0 h-0 border-l-[5px] border-l-transparent border-r-[5px] border-r-transparent border-t-[5px] border-t-primary/40" />
              )}
            </div>

            {/* Connector line */}
            {index < steps.length - 1 && (
              <div
                className={cn(
                  'flex-1 h-0.5 mx-1 rounded-full transition-all duration-500',
                  isCompleted ? styles.line : 'bg-border',
                  status === 'pending' && 'bg-border/50',
                )}
              />
            )}
          </div>
        );
      })}
    </div>
  );
}
