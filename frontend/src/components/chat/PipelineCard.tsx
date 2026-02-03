'use client';

import { useState, useMemo } from 'react';
import { ChevronDown, ChevronUp } from 'lucide-react';
import { Badge } from '@/components/ui/badge';
import { Progress } from '@/components/ui/progress';
import { useStore } from '@/lib/store';
import { PipelineRunner } from '@/components/pipeline/PipelineRunner';
import { getPipeline } from '@/lib/pipelines';
import type { StepResult } from '@/lib/pipeline-types';
import { cn } from '@/lib/utils';

interface PipelineCardProps {
  pipelineId: string;
  onSelectDesign?: (pdbContent: string) => void;
  onStepComplete?: (stepId: string, result: StepResult) => void;
  onPipelineComplete?: (results: Record<string, StepResult>) => void;
  onPipelineFailed?: (stepId: string, error: string) => void;
  onCancel?: () => void;
  initialParams?: Record<string, unknown>;
}

export function PipelineCard({
  pipelineId,
  onSelectDesign,
  onStepComplete,
  onPipelineComplete,
  onPipelineFailed,
  onCancel,
  initialParams = {},
}: PipelineCardProps) {
  const [expanded, setExpanded] = useState(true);
  const definition = useMemo(() => getPipeline(pipelineId), [pipelineId]);

  if (!definition) {
    return (
      <div className="bg-destructive/5 border border-destructive/20 rounded-xl p-3 text-sm text-destructive">
        Unknown pipeline: {pipelineId}
      </div>
    );
  }

  const Icon = definition.icon;

  return (
    <div className="border border-border rounded-xl overflow-hidden bg-card">
      {/* Collapsible header */}
      <button
        onClick={() => setExpanded((prev) => !prev)}
        className="w-full flex items-center gap-3 px-4 py-3 hover:bg-muted/50 transition-colors"
      >
        <div className="flex items-center justify-center w-7 h-7 rounded-full bg-primary/10">
          <Icon className="h-4 w-4 text-primary" />
        </div>
        <div className="flex-1 text-left">
          <div className="text-sm font-medium text-foreground">{definition.name}</div>
          <div className="text-[10px] text-muted-foreground">{definition.description}</div>
        </div>
        {expanded ? (
          <ChevronUp className="h-4 w-4 text-muted-foreground" />
        ) : (
          <ChevronDown className="h-4 w-4 text-muted-foreground" />
        )}
      </button>

      {/* Expanded pipeline runner */}
      {expanded && (
        <div className="border-t border-border">
          <PipelineRunner
            definition={definition}
            initialParams={initialParams}
            onStepComplete={onStepComplete}
            onPipelineComplete={onPipelineComplete}
            onPipelineFailed={onPipelineFailed}
            onDesignSelected={onSelectDesign}
            onCancel={onCancel}
          />
        </div>
      )}
    </div>
  );
}
