'use client';

import { Badge } from '@/components/ui/badge';
import { Progress } from '@/components/ui/progress';
import {
  CheckCircle,
  XCircle,
  AlertTriangle,
  Shield,
  Filter,
} from 'lucide-react';
import { cn } from '@/lib/utils';
import type { StepResult } from '@/lib/pipeline-types';

interface ScoutResultPreviewProps {
  result: StepResult;
  onSelectDesign?: (pdbContent: string) => void;
}

interface PreFilterCheck {
  passed: boolean;
  reasons?: string[];
  chain_breaks?: Array<{ chain: string; res_i: number; res_j: number; distance: number }>;
  coordination_number?: number;
  geometry_rmsd?: number;
  ligand_name?: string;
  atom_count?: number;
  error?: string;
}

interface PreFilterResult {
  passed: boolean;
  checks: Record<string, PreFilterCheck>;
  failed_checks: string[];
  skipped_checks: string[];
}

/** Format pre-filter check name for display. */
function formatCheckName(check: string): string {
  return check.replace(/_/g, ' ').replace(/\b\w/g, c => c.toUpperCase());
}

/** Build tooltip text summarizing why pre-filter failed. */
function preFilterTooltip(pfResult: PreFilterResult): string {
  if (pfResult.passed) return 'Pre-filter: passed all checks';

  const lines: string[] = ['Pre-filter failed:'];
  for (const checkName of pfResult.failed_checks) {
    const check = pfResult.checks[checkName];
    if (!check) continue;

    const name = formatCheckName(checkName);
    if (check.reasons && check.reasons.length > 0) {
      lines.push(`• ${name}: ${check.reasons.join('; ')}`);
    } else if (check.chain_breaks && check.chain_breaks.length > 0) {
      const breaks = check.chain_breaks
        .map(b => `${b.chain}:${b.res_i}-${b.res_j} (${b.distance}Å)`)
        .join(', ');
      lines.push(`• ${name}: breaks at ${breaks}`);
    } else {
      lines.push(`• ${name}`);
    }
  }
  return lines.join('\n');
}

/** Extract a short metric summary from pre-filter check data. */
function getCheckMetric(check: PreFilterCheck, checkName: string): string | null {
  if (checkName === 'metal_coordination' && check.coordination_number !== undefined) {
    const parts = [`CN: ${check.coordination_number}`];
    if (check.geometry_rmsd !== undefined) parts.push(`geom: ${check.geometry_rmsd.toFixed(2)}Å`);
    return parts.join(', ');
  }
  if (checkName === 'chain_continuity' && check.chain_breaks) {
    return check.chain_breaks.length === 0 ? 'no breaks' : `${check.chain_breaks.length} break(s)`;
  }
  if (checkName === 'ligand_presence' && check.atom_count !== undefined) {
    return `${check.atom_count} atoms`;
  }
  return null;
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
  const rawScoutResults = (data.scout_results as Array<Record<string, unknown>>) || [];

  // Normalize scout results — coerce types from JSON
  const scoutResults = rawScoutResults.map(sr => ({
    backbone_index: Number(sr.backbone_index ?? 0),
    ptm: Number(sr.ptm ?? 0),
    plddt: Number(sr.plddt ?? 0),
    passed: Boolean(sr.passed),
    sequence: String(sr.sequence ?? ''),
    pre_filter_failed: Boolean(sr.pre_filter_failed),
  }));

  // Fix data inconsistency: if passing_count disagrees with individual entries,
  // trust passing_count and derive per-entry status from it
  const individualPassCount = scoutResults.filter(sr => sr.passed).length;
  if (passingCount > 0 && individualPassCount === 0 && scoutResults.length > 0) {
    console.warn('[ScoutPreview] Data inconsistency: passing_count=', passingCount,
      'but individual passed count=', individualPassCount, '— deriving from counts');
    // All passed if passing_count === original_count, otherwise can't determine which
    if (passingCount === originalCount) {
      for (const sr of scoutResults) {
        if (!sr.pre_filter_failed) {
          sr.passed = true;
        }
      }
    }
  }

  const ptmThreshold = data.ptm_threshold as number || 0.6;
  const plddtThreshold = data.plddt_threshold as number || 0.65;

  const preFilterResults = (data.pre_filter_results as PreFilterResult[] | undefined) || [];
  const preFilterFailed = data.pre_filter_failed as number || 0;
  const hasPreFilter = preFilterResults.length > 0;
  const isPreFilterOnly = !!data.metal_pre_filter_only;

  // Detect whether any scout result has actual RF3 metrics (non-zero ptm/plddt)
  const hasRf3Metrics = !isPreFilterOnly && scoutResults.some(
    sr => !sr.pre_filter_failed && (sr.ptm > 0 || sr.plddt > 0)
  );

  // Collect unique check names from pre-filter results for column headers
  const checkNames = hasPreFilter
    ? [...new Set(preFilterResults.flatMap(pf => pf.checks ? Object.keys(pf.checks) : []))]
    : [];

  const passRate = originalCount > 0 ? (passingCount / originalCount) * 100 : 0;

  return (
    <div className="space-y-3">
      {/* Pre-filter summary banner (amber) — only when some failed */}
      {hasPreFilter && preFilterFailed > 0 && (
        <div className="flex items-center gap-2 rounded-md bg-amber-500/10 border border-amber-500/20 px-3 py-2 text-xs">
          <Filter className="h-3.5 w-3.5 text-amber-500 shrink-0" />
          <span className="text-amber-700 dark:text-amber-400">
            Pre-filter eliminated <span className="font-semibold">{preFilterFailed}</span> backbone{preFilterFailed !== 1 ? 's' : ''} (broken chain, missing metal/ligand) before GPU validation
          </span>
        </div>
      )}

      {/* Pass/fail count bar */}
      <div className="space-y-1.5">
        <div className="flex items-center justify-between text-xs">
          <div className="flex items-center gap-1.5">
            <Shield className="h-3.5 w-3.5 text-primary" />
            <span className="font-medium">
              {passingCount}/{originalCount} backbones passed {isPreFilterOnly ? 'pre-filter' : 'scout validation'}
            </span>
          </div>
          <Badge variant={passingCount > 0 ? 'default' : 'destructive'} className="text-[10px] h-5">
            {passRate.toFixed(0)}% pass rate
          </Badge>
        </div>
        <Progress value={passRate} className="h-2" />
      </div>

      {/* Threshold indicators (shown when RF3 metrics exist) */}
      {hasRf3Metrics && (
        <div className="flex gap-3 text-[10px] text-muted-foreground">
          <span>pTM threshold: <span className="font-mono font-medium text-foreground">{ptmThreshold}</span></span>
          <span>pLDDT threshold: <span className="font-mono font-medium text-foreground">{plddtThreshold}</span></span>
        </div>
      )}

      {/* Per-backbone results table */}
      {scoutResults.length > 0 && (
        <div className="border rounded-md overflow-hidden">
          <table className="w-full text-[10px]">
            <thead>
              <tr className="bg-muted/50 border-b">
                <th className="px-2 py-1.5 text-left font-medium text-muted-foreground">Backbone</th>
                {/* Pre-filter check columns — show actual metrics */}
                {isPreFilterOnly && checkNames.map(name => (
                  <th key={name} className="px-2 py-1.5 text-center font-medium text-muted-foreground">
                    {formatCheckName(name)}
                  </th>
                ))}
                {/* Non-pre-filter: show pre-filter column if available */}
                {!isPreFilterOnly && hasPreFilter && (
                  <th className="px-2 py-1.5 text-center font-medium text-muted-foreground">Pre-filter</th>
                )}
                {hasRf3Metrics && (
                  <>
                    <th className="px-2 py-1.5 text-right font-medium text-muted-foreground">pTM</th>
                    <th className="px-2 py-1.5 text-right font-medium text-muted-foreground">pLDDT</th>
                  </>
                )}
                <th className="px-2 py-1.5 text-center font-medium text-muted-foreground">Status</th>
              </tr>
            </thead>
            <tbody>
              {scoutResults.map((sr) => {
                const pfResult = preFilterResults[sr.backbone_index] as PreFilterResult | undefined;
                const isPreFilterFail = sr.pre_filter_failed;
                const isPassed = sr.passed;

                return (
                  <tr
                    key={sr.backbone_index}
                    className={cn(
                      'border-b last:border-b-0',
                      isPreFilterFail
                        ? 'bg-amber-500/5'
                        : isPassed
                          ? 'bg-emerald-500/5'
                          : 'bg-red-500/5',
                    )}
                  >
                    <td className="px-2 py-1.5 font-mono">#{sr.backbone_index + 1}</td>

                    {/* Pre-filter-only: show per-check metrics */}
                    {isPreFilterOnly && checkNames.map(checkName => {
                      const check = pfResult?.checks?.[checkName];
                      if (!check) {
                        return <td key={checkName} className="px-2 py-1.5 text-center text-muted-foreground">—</td>;
                      }
                      const metric = getCheckMetric(check, checkName);
                      return (
                        <td
                          key={checkName}
                          className={cn(
                            'px-2 py-1.5 text-center font-mono',
                            check.passed
                              ? 'text-emerald-600 dark:text-emerald-400'
                              : 'text-red-500 dark:text-red-400',
                          )}
                          title={check.reasons?.join('; ') || (check.passed ? 'Passed' : 'Failed')}
                        >
                          {metric || (check.passed ? '✓' : '✗')}
                        </td>
                      );
                    })}

                    {/* Non-pre-filter: single pre-filter status column */}
                    {!isPreFilterOnly && hasPreFilter && (
                      <td className="px-2 py-1.5 text-center">
                        {pfResult ? (
                          <span title={preFilterTooltip(pfResult)}>
                            {pfResult.passed ? (
                              <CheckCircle className="h-3 w-3 text-emerald-500 inline" />
                            ) : (
                              <AlertTriangle className="h-3 w-3 text-amber-500 inline" />
                            )}
                          </span>
                        ) : (
                          <span className="text-muted-foreground">—</span>
                        )}
                      </td>
                    )}

                    {hasRf3Metrics && (
                      <>
                        <td className={cn(
                          'px-2 py-1.5 text-right font-mono',
                          isPreFilterFail
                            ? 'text-muted-foreground'
                            : sr.ptm >= ptmThreshold
                              ? 'text-emerald-600 dark:text-emerald-400'
                              : 'text-red-500 dark:text-red-400',
                        )}>
                          {isPreFilterFail ? '—' : sr.ptm.toFixed(3)}
                        </td>
                        <td className={cn(
                          'px-2 py-1.5 text-right font-mono',
                          isPreFilterFail
                            ? 'text-muted-foreground'
                            : sr.plddt >= plddtThreshold
                              ? 'text-emerald-600 dark:text-emerald-400'
                              : 'text-red-500 dark:text-red-400',
                        )}>
                          {isPreFilterFail ? '—' : sr.plddt.toFixed(3)}
                        </td>
                      </>
                    )}

                    <td className="px-2 py-1.5 text-center">
                      {isPreFilterFail ? (
                        <span title={pfResult ? preFilterTooltip(pfResult) : 'Pre-filter failed'}>
                          <AlertTriangle className="h-3.5 w-3.5 text-amber-500 inline" />
                        </span>
                      ) : isPassed ? (
                        <CheckCircle className="h-3.5 w-3.5 text-emerald-500 inline" />
                      ) : (
                        <XCircle className="h-3.5 w-3.5 text-red-500 inline" />
                      )}
                    </td>
                  </tr>
                );
              })}
            </tbody>
          </table>
        </div>
      )}

      {/* Fallback: no per-backbone data but we have summary counts */}
      {scoutResults.length === 0 && originalCount > 0 && (
        <div className="text-xs text-muted-foreground">
          {passingCount === originalCount
            ? 'All backbones passed scout validation.'
            : `${originalCount - passingCount} backbone(s) were filtered out.`}
        </div>
      )}
    </div>
  );
}
