'use client';

import { Badge } from '@/components/ui/badge';
import { Progress } from '@/components/ui/progress';
import {
  CheckCircle,
  XCircle,
  AlertTriangle,
  Shield,
} from 'lucide-react';
import { cn } from '@/lib/utils';
import type { StepResult } from '@/lib/pipeline-types';

interface ScoutResultPreviewProps {
  result: StepResult;
  onSelectDesign?: (pdbContent: string) => void;
}

export function ScoutResultPreview({
  result,
}: ScoutResultPreviewProps) {
  const data = result.data || {};

  // Handle skipped state
  if (data.skipped) {
    return (
      <div className="flex items-center gap-2 text-xs text-muted-foreground">
        <AlertTriangle className="h-3.5 w-3.5" />
        <span>{data.reason as string || data.summary as string || 'Scout filter skipped'}</span>
      </div>
    );
  }

  const originalCount = data.original_count as number || 0;
  const passingCount = data.passing_count as number || 0;
  const scoutResults = (data.scout_results as Array<{
    backbone_index: number;
    ptm: number;
    plddt: number;
    passed: boolean;
    sequence: string;
  }>) || [];
  const ptmThreshold = data.ptm_threshold as number || 0.6;
  const plddtThreshold = data.plddt_threshold as number || 0.65;

  const passRate = originalCount > 0 ? (passingCount / originalCount) * 100 : 0;

  return (
    <div className="space-y-3">
      {/* Pass/fail count bar */}
      <div className="space-y-1.5">
        <div className="flex items-center justify-between text-xs">
          <div className="flex items-center gap-1.5">
            <Shield className="h-3.5 w-3.5 text-primary" />
            <span className="font-medium">
              {passingCount}/{originalCount} backbones passed scout validation
            </span>
          </div>
          <Badge variant={passingCount > 0 ? 'default' : 'destructive'} className="text-[10px] h-5">
            {passRate.toFixed(0)}% pass rate
          </Badge>
        </div>
        <Progress value={passRate} className="h-2" />
      </div>

      {/* Threshold indicators */}
      <div className="flex gap-3 text-[10px] text-muted-foreground">
        <span>pTM threshold: <span className="font-mono font-medium text-foreground">{ptmThreshold}</span></span>
        <span>pLDDT threshold: <span className="font-mono font-medium text-foreground">{plddtThreshold}</span></span>
      </div>

      {/* Per-backbone results table */}
      {scoutResults.length > 0 && (
        <div className="border rounded-md overflow-hidden">
          <table className="w-full text-[10px]">
            <thead>
              <tr className="bg-muted/50 border-b">
                <th className="px-2 py-1.5 text-left font-medium text-muted-foreground">Backbone</th>
                <th className="px-2 py-1.5 text-right font-medium text-muted-foreground">pTM</th>
                <th className="px-2 py-1.5 text-right font-medium text-muted-foreground">pLDDT</th>
                <th className="px-2 py-1.5 text-center font-medium text-muted-foreground">Status</th>
              </tr>
            </thead>
            <tbody>
              {scoutResults.map((sr) => (
                <tr
                  key={sr.backbone_index}
                  className={cn(
                    'border-b last:border-b-0',
                    sr.passed ? 'bg-emerald-500/5' : 'bg-red-500/5',
                  )}
                >
                  <td className="px-2 py-1.5 font-mono">#{sr.backbone_index + 1}</td>
                  <td className={cn(
                    'px-2 py-1.5 text-right font-mono',
                    sr.ptm >= ptmThreshold ? 'text-emerald-600 dark:text-emerald-400' : 'text-red-500 dark:text-red-400',
                  )}>
                    {sr.ptm.toFixed(3)}
                  </td>
                  <td className={cn(
                    'px-2 py-1.5 text-right font-mono',
                    sr.plddt >= plddtThreshold ? 'text-emerald-600 dark:text-emerald-400' : 'text-red-500 dark:text-red-400',
                  )}>
                    {sr.plddt.toFixed(3)}
                  </td>
                  <td className="px-2 py-1.5 text-center">
                    {sr.passed ? (
                      <CheckCircle className="h-3.5 w-3.5 text-emerald-500 inline" />
                    ) : (
                      <XCircle className="h-3.5 w-3.5 text-red-500 inline" />
                    )}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}
    </div>
  );
}
