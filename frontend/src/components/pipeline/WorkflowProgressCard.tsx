'use client';

/**
 * WorkflowProgressCard — Displays step-by-step progress for workflow_run jobs.
 *
 * Shows each workflow step with status icon, name, summary, and timing.
 * Reuses existing Card styling from the design system.
 */

import React from 'react';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { CheckCircle2, XCircle, Loader2, Circle } from 'lucide-react';
import type { WorkflowProgress, WorkflowStepResult } from '@/lib/workflow-types';
import { MODULE_DISPLAY_NAMES } from '@/lib/workflow-types';

interface WorkflowProgressCardProps {
  progress: WorkflowProgress;
  /** Module names in order from the workflow spec */
  stepNames: string[];
}

export function WorkflowProgressCard({ progress, stepNames }: WorkflowProgressCardProps) {
  const { current_step, step_results, status } = progress;

  return (
    <Card className="w-full">
      <CardHeader className="pb-3">
        <CardTitle className="text-sm font-medium flex items-center justify-between">
          <span>{progress.workflow_name || 'Workflow'}</span>
          <StatusBadge status={status} />
        </CardTitle>
      </CardHeader>
      <CardContent className="space-y-2">
        {stepNames.map((stepName, idx) => {
          const result = step_results[stepName];
          const isCurrent = stepName === current_step && status === 'running';
          const isCompleted = result?.status === 'completed';
          const isFailed = result?.status === 'failed';
          const isPending = !result && !isCurrent;

          return (
            <div
              key={stepName}
              className={`flex items-center gap-3 rounded-md px-3 py-2 text-sm transition-colors ${
                isCurrent
                  ? 'bg-blue-50 dark:bg-blue-950 border border-blue-200 dark:border-blue-800'
                  : isCompleted
                  ? 'bg-green-50/50 dark:bg-green-950/30'
                  : isFailed
                  ? 'bg-red-50/50 dark:bg-red-950/30'
                  : 'opacity-50'
              }`}
            >
              {/* Step icon */}
              <StepIcon
                isCurrent={isCurrent}
                isCompleted={isCompleted}
                isFailed={isFailed}
              />

              {/* Step info */}
              <div className="flex-1 min-w-0">
                <div className="font-medium truncate">
                  {MODULE_DISPLAY_NAMES[stepName] || stepName}
                </div>
                {result?.summary && (
                  <div className="text-xs text-muted-foreground truncate">
                    {result.summary}
                  </div>
                )}
              </div>

              {/* Timing */}
              {result?.timing != null && (
                <span className="text-xs text-muted-foreground whitespace-nowrap">
                  {result.timing < 1
                    ? `${(result.timing * 1000).toFixed(0)}ms`
                    : `${result.timing.toFixed(1)}s`}
                </span>
              )}

              {/* Step number */}
              <span className="text-xs text-muted-foreground">
                {idx + 1}/{stepNames.length}
              </span>
            </div>
          );
        })}
      </CardContent>
    </Card>
  );
}

function StepIcon({
  isCurrent,
  isCompleted,
  isFailed,
}: {
  isCurrent: boolean;
  isCompleted: boolean;
  isFailed: boolean;
}) {
  if (isCurrent) {
    return <Loader2 className="h-4 w-4 text-blue-500 animate-spin flex-shrink-0" />;
  }
  if (isCompleted) {
    return <CheckCircle2 className="h-4 w-4 text-green-500 flex-shrink-0" />;
  }
  if (isFailed) {
    return <XCircle className="h-4 w-4 text-red-500 flex-shrink-0" />;
  }
  return <Circle className="h-4 w-4 text-muted-foreground/40 flex-shrink-0" />;
}

function StatusBadge({ status }: { status: string }) {
  const colors: Record<string, string> = {
    running: 'bg-blue-100 text-blue-700 dark:bg-blue-900 dark:text-blue-300',
    completed: 'bg-green-100 text-green-700 dark:bg-green-900 dark:text-green-300',
    failed: 'bg-red-100 text-red-700 dark:bg-red-900 dark:text-red-300',
  };

  return (
    <span
      className={`px-2 py-0.5 rounded-full text-xs font-medium ${
        colors[status] || 'bg-gray-100 text-gray-700'
      }`}
    >
      {status}
    </span>
  );
}
