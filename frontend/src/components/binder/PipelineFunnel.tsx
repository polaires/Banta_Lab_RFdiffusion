'use client';

import {
  Filter,
  Blocks,
  FileEdit,
  Brain,
  Network,
  FilterX,
  CheckCircle2,
  type LucideIcon,
} from 'lucide-react';

// Map stage names to Lucide icons
const STAGE_ICONS: Record<string, LucideIcon> = {
  'Generated': Blocks,
  'MPNN Designed': FileEdit,
  'ESM Passed': Brain,
  'Interface Analyzed': Network,
  'Passed Filters': FilterX,
  'Returned': CheckCircle2,
};

interface PipelineStage {
  name: string;
  count: number;
  color: string;
  icon: string;
}

interface PipelineFunnelProps {
  statistics: {
    generated: number;
    mpnn_designed: number;
    esm_passed: number;
    relaxed?: number;
    interface_analyzed: number;
    passed_filters: number;
    returned: number;
  };
}

export function PipelineFunnel({ statistics }: PipelineFunnelProps) {
  const maxCount = statistics.generated || 1;

  const stages: PipelineStage[] = [
    { name: 'Generated', count: statistics.generated, color: 'bg-muted/500', icon: 'architecture' },
    { name: 'MPNN Designed', count: statistics.mpnn_designed, color: 'bg-blue-500', icon: 'edit_note' },
    { name: 'ESM Passed', count: statistics.esm_passed, color: 'bg-violet-500', icon: 'psychology' },
    { name: 'Interface Analyzed', count: statistics.interface_analyzed, color: 'bg-cyan-500', icon: 'hub' },
    { name: 'Passed Filters', count: statistics.passed_filters, color: 'bg-emerald-500', icon: 'filter_alt' },
    { name: 'Returned', count: statistics.returned, color: 'bg-green-600', icon: 'check_circle' },
  ];

  // Calculate pass rate
  const passRate = statistics.generated > 0
    ? ((statistics.returned / statistics.generated) * 100).toFixed(0)
    : '0';

  return (
    <div className="bg-card rounded-xl border border-border overflow-hidden">
      {/* Header */}
      <div className="bg-muted px-4 py-3 border-b border-border">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-2">
            <Filter className="h-5 w-5 text-muted-foreground" />
            <h4 className="font-semibold text-foreground text-sm">Pipeline Progress</h4>
          </div>
          <div className="flex items-center gap-2">
            <span className="text-xs text-muted-foreground">Pass Rate:</span>
            <span className={`text-sm font-bold ${
              parseInt(passRate) >= 50 ? 'text-green-600' :
              parseInt(passRate) >= 25 ? 'text-muted-foreground' : 'text-red-600'
            }`}>
              {passRate}%
            </span>
          </div>
        </div>
      </div>

      {/* Funnel Visualization */}
      <div className="p-4 space-y-3">
        {stages.map((stage, index) => {
          const widthPercent = (stage.count / maxCount) * 100;
          const dropoff = index > 0 ? stages[index - 1].count - stage.count : 0;

          return (
            <div key={stage.name} className="relative">
              {/* Stage Row */}
              <div className="flex items-center gap-3">
                {/* Icon */}
                <div className={`w-8 h-8 rounded-lg ${stage.color} bg-opacity-10 flex items-center justify-center`}>
                  {(() => {
                    const IconComponent = STAGE_ICONS[stage.name];
                    return IconComponent ? (
                      <IconComponent className={`h-4 w-4 ${stage.color.replace('bg-', 'text-')}`} />
                    ) : null;
                  })()}
                </div>

                {/* Bar Container */}
                <div className="flex-1">
                  <div className="flex items-center justify-between mb-1">
                    <span className="text-xs text-muted-foreground">{stage.name}</span>
                    <div className="flex items-center gap-2">
                      <span className="text-sm font-semibold text-foreground">{stage.count}</span>
                      {dropoff > 0 && (
                        <span className="text-xs text-red-500">-{dropoff}</span>
                      )}
                    </div>
                  </div>

                  {/* Progress Bar */}
                  <div className="h-2 bg-muted rounded-full overflow-hidden">
                    <div
                      className={`h-full ${stage.color} transition-all duration-500 ease-out rounded-full`}
                      style={{ width: `${widthPercent}%` }}
                    />
                  </div>
                </div>
              </div>

              {/* Connector Arrow (except last) */}
              {index < stages.length - 1 && (
                <div className="ml-4 h-2 border-l-2 border-dashed border-border" />
              )}
            </div>
          );
        })}
      </div>

      {/* Summary Footer */}
      <div className="bg-muted/50 px-4 py-2 border-t border-border">
        <div className="flex items-center justify-between text-xs">
          <span className="text-muted-foreground">
            {statistics.generated - statistics.returned} designs filtered out
          </span>
          <span className="text-muted-foreground">
            {statistics.returned} high-quality designs available
          </span>
        </div>
      </div>
    </div>
  );
}
