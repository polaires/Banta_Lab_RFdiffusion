'use client';

import { Badge } from '@/components/ui/badge';
import {
  AlertTriangle,
  CheckCircle,
  TrendingUp,
  Lightbulb,
  History,
} from 'lucide-react';
import { cn } from '@/lib/utils';
import type { StepResult } from '@/lib/pipeline-types';

interface LessonResultPreviewProps {
  result: StepResult;
  onSelectDesign?: (pdbContent: string) => void;
}

export function LessonResultPreview({
  result,
}: LessonResultPreviewProps) {
  const data = result.data || {};

  const triggerDetected = data.trigger_detected as boolean;
  const trigger = data.trigger as {
    type: 'failure_pattern' | 'breakthrough' | 'improvement';
    description: string;
    relevant_designs: string[];
    metrics_involved: string[];
  } | undefined;
  const historyCount = data.history_count as number || 0;

  // No trigger detected
  if (!triggerDetected || !trigger) {
    return (
      <div className="flex items-center gap-2 text-xs text-muted-foreground">
        <CheckCircle className="h-3.5 w-3.5 text-muted-foreground" />
        <span>No significant patterns detected</span>
        <Badge variant="outline" className="text-[10px] h-5 gap-1">
          <History className="h-3 w-3" />
          {historyCount} designs in history
        </Badge>
      </div>
    );
  }

  // Trigger detected — style by type
  const typeConfig = {
    failure_pattern: {
      icon: AlertTriangle,
      bgClass: 'bg-red-500/10 border-red-500/30',
      textClass: 'text-red-700 dark:text-red-400',
      badgeVariant: 'destructive' as const,
      label: 'Failure Pattern',
    },
    breakthrough: {
      icon: Lightbulb,
      bgClass: 'bg-emerald-500/10 border-emerald-500/30',
      textClass: 'text-emerald-700 dark:text-emerald-400',
      badgeVariant: 'default' as const,
      label: 'Breakthrough',
    },
    improvement: {
      icon: TrendingUp,
      bgClass: 'bg-blue-500/10 border-blue-500/30',
      textClass: 'text-blue-700 dark:text-blue-400',
      badgeVariant: 'secondary' as const,
      label: 'Improvement',
    },
  };

  const config = typeConfig[trigger.type] || typeConfig.improvement;
  const Icon = config.icon;

  return (
    <div className="space-y-2">
      {/* Trigger banner */}
      <div className={cn(
        'flex items-start gap-2 rounded-md border px-3 py-2',
        config.bgClass,
      )}>
        <Icon className={cn('h-4 w-4 mt-0.5 shrink-0', config.textClass)} />
        <div className="space-y-1 flex-1">
          <div className="flex items-center gap-2">
            <Badge variant={config.badgeVariant} className="text-[10px] h-5">
              {config.label}
            </Badge>
            <Badge variant="outline" className="text-[10px] h-5 gap-1">
              <History className="h-3 w-3" />
              {historyCount} in history
            </Badge>
          </div>
          <p className={cn('text-xs', config.textClass)}>
            {trigger.description}
          </p>
        </div>
      </div>

      {/* Metrics involved */}
      {trigger.metrics_involved.length > 0 && (
        <div className="flex flex-wrap gap-1.5">
          <span className="text-[10px] text-muted-foreground">Metrics:</span>
          {trigger.metrics_involved.map((metric) => (
            <Badge key={metric} variant="outline" className="text-[10px] h-5 font-mono">
              {metric}
            </Badge>
          ))}
        </div>
      )}

      {/* Related designs */}
      {trigger.relevant_designs.length > 0 && (
        <div className="flex flex-wrap gap-1.5">
          <span className="text-[10px] text-muted-foreground">Related:</span>
          {trigger.relevant_designs.slice(0, 5).map((design) => (
            <Badge key={design} variant="outline" className="text-[10px] h-5 font-mono">
              {design}
            </Badge>
          ))}
          {trigger.relevant_designs.length > 5 && (
            <span className="text-[10px] text-muted-foreground">
              +{trigger.relevant_designs.length - 5} more
            </span>
          )}
        </div>
      )}
    </div>
  );
}
