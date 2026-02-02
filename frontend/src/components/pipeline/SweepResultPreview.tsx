'use client';

import { Card } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { Progress } from '@/components/ui/progress';
import { Button } from '@/components/ui/button';
import {
  CheckCircle,
  Eye,
  XCircle,
  TrendingUp,
  Check,
} from 'lucide-react';
import { cn } from '@/lib/utils';
import type { StepResult } from '@/lib/pipeline-types';

interface SweepResultPreviewProps {
  result: StepResult;
  onSelectDesign?: (pdbContent: string) => void;
  // Selection props (passed by StepCard when supportsSelection=true)
  selectedIds?: string[];
  onSelectionChange?: (ids: string[]) => void;
  isPaused?: boolean;
}

interface FailedDesign {
  name: string;
  config_name: string;
  plddt: number;
  ptm: number;
  pae: number;
  tier: string;
  status: string;
  sequence?: string;
  pdb_content?: string;
}

function getTierBadgeVariant(tier: string): 'default' | 'secondary' | 'destructive' {
  if (tier === 'S' || tier === 'A') return 'default';
  if (tier === 'B' || tier === 'C') return 'secondary';
  return 'destructive';
}

function MetricBar({ label, value, threshold, max }: {
  label: string;
  value: number;
  threshold: number;
  max: number;
}) {
  const pct = Math.min(100, (value / max) * 100);
  const thresholdPct = (threshold / max) * 100;
  const passing = value >= threshold;

  return (
    <div className="space-y-0.5">
      <div className="flex items-center justify-between text-[10px]">
        <span className="text-muted-foreground">{label}</span>
        <span className={cn('font-mono font-medium', passing ? 'text-emerald-600 dark:text-emerald-400' : 'text-red-500 dark:text-red-400')}>
          {value.toFixed(3)}
        </span>
      </div>
      <div className="relative">
        <Progress value={pct} className="h-1.5" />
        {/* Threshold marker */}
        <div
          className="absolute top-0 h-1.5 w-px bg-foreground/50"
          style={{ left: `${thresholdPct}%` }}
        />
      </div>
      <div className="flex justify-between text-[9px] text-muted-foreground">
        <span>0</span>
        <span>threshold: {threshold}</span>
      </div>
    </div>
  );
}

/** Count residues in a PDB string (unique residue numbers in ATOM lines). */
function countResidues(pdbContent: string): number {
  const residues = new Set<string>();
  for (const line of pdbContent.split('\n')) {
    if (line.startsWith('ATOM')) {
      residues.add(line.substring(22, 26).trim());
    }
  }
  return residues.size;
}

/** Count HETATM lines (metal/ligand atoms) in a PDB string. */
function countHetatm(pdbContent: string): number {
  let count = 0;
  for (const line of pdbContent.split('\n')) {
    if (line.startsWith('HETATM')) count++;
  }
  return count;
}

export function SweepResultPreview({
  result,
  onSelectDesign,
  selectedIds,
  onSelectionChange,
  isPaused,
}: SweepResultPreviewProps) {
  const data = result.data || {};
  const allFailed = data.sweep_all_failed as boolean;
  const failedDesigns = (data.failed_designs as FailedDesign[]) || [];
  const totalGenerated = data.total_generated as number || 0;
  const totalPassing = data.total_passing as number || 0;
  const totalReview = data.total_review as number || 0;
  const passRate = data.pass_rate as number || 0;
  const bestDesign = data.best_design as FailedDesign | null;
  const configRankings = data.config_rankings as Array<{ config: string; avg_ptm: number; avg_plddt: number; count: number }> | null;
  const configErrors = data.config_errors as Record<string, string[]> | null;

  // Non-sweep result (general RFD3 or metal single) — show backbone cards with selection
  if (!data.session_id) {
    const pdbs = result.pdbOutputs || [];
    const isMetal = !!data.metal_single_mode;
    const selectedSet = new Set(selectedIds || []);

    const toggleSelection = (id: string) => {
      if (!onSelectionChange) return;
      const newIds = selectedSet.has(id)
        ? [...selectedSet].filter(x => x !== id)
        : [...selectedSet, id];
      onSelectionChange(newIds);
    };

    const selectAll = () => {
      if (!onSelectionChange) return;
      onSelectionChange(pdbs.map(p => p.id));
    };

    return (
      <div className="space-y-2">
        <div className="flex items-center justify-between">
          <p className="text-xs text-muted-foreground">{result.summary}</p>
          {isPaused && pdbs.length > 1 && onSelectionChange && (
            <Button
              variant="ghost"
              size="sm"
              className="h-5 text-[10px] px-1.5"
              onClick={selectAll}
            >
              Select all
            </Button>
          )}
        </div>

        {pdbs.length > 0 && (
          <div className="grid gap-1.5">
            {pdbs.map((pdb) => {
              const residueCount = countResidues(pdb.pdbContent);
              const hetatmCount = countHetatm(pdb.pdbContent);
              const isSelected = selectedSet.has(pdb.id);

              return (
                <Card
                  key={pdb.id}
                  className={cn(
                    'p-2 flex items-center gap-2 cursor-pointer transition-colors',
                    isSelected ? 'border-primary/50 bg-primary/5' : 'hover:bg-accent/50',
                  )}
                  onClick={() => isPaused && toggleSelection(pdb.id)}
                >
                  {/* Selection checkbox (only when paused for review) */}
                  {isPaused && (
                    <div className={cn(
                      'w-4 h-4 rounded border flex items-center justify-center shrink-0',
                      isSelected ? 'bg-primary border-primary' : 'border-muted-foreground/40',
                    )}>
                      {isSelected && <Check className="h-3 w-3 text-primary-foreground" />}
                    </div>
                  )}

                  <div className="flex-1 min-w-0">
                    <div className="flex items-center gap-1.5">
                      <span className="text-xs font-medium">{pdb.label}</span>
                      <Badge variant="secondary" className="text-[9px] h-4">
                        {residueCount} res
                      </Badge>
                      {isMetal && hetatmCount > 0 && (
                        <Badge variant="outline" className="text-[9px] h-4">
                          {hetatmCount} HETATM
                        </Badge>
                      )}
                    </div>
                  </div>

                  {onSelectDesign && (
                    <Button
                      variant="ghost"
                      size="sm"
                      className="h-6 w-6 p-0 shrink-0"
                      onClick={(e) => {
                        e.stopPropagation();
                        onSelectDesign(pdb.pdbContent);
                      }}
                    >
                      <Eye className="h-3 w-3" />
                    </Button>
                  )}
                </Card>
              );
            })}
          </div>
        )}

        {isPaused && pdbs.length > 0 && (
          <p className="text-[10px] text-muted-foreground">
            Select backbone(s) to carry forward to sequence design.
          </p>
        )}
      </div>
    );
  }

  // Successful sweep — show passing designs summary
  if (!allFailed && totalPassing > 0) {
    return (
      <div className="space-y-2">
        <div className="flex items-center gap-2 text-xs">
          <CheckCircle className="h-3.5 w-3.5 text-emerald-500" />
          <span className="font-medium">{result.summary}</span>
        </div>
        {result.pdbOutputs && result.pdbOutputs.length > 0 && (
          <div className="text-[10px] text-muted-foreground">
            {result.pdbOutputs.length} structure(s) ready for sequence design
          </div>
        )}
      </div>
    );
  }

  // All designs failed — show detailed failure analysis
  return (
    <div className="space-y-3">
      {/* Failure banner */}
      <div className="flex items-start gap-2 p-2.5 rounded-md bg-destructive/5 border border-destructive/20">
        <XCircle className="h-4 w-4 text-destructive shrink-0 mt-0.5" />
        <div className="space-y-1">
          <p className="text-xs font-medium text-destructive">
            {totalGenerated === 0
              ? 'No designs were generated — all attempts failed'
              : `All ${totalGenerated} design(s) failed validation`}
          </p>
          {totalGenerated > 0 && (
            <p className="text-[10px] text-muted-foreground">
              No designs met the minimum quality thresholds (pTM &ge; 0.5 for B-tier).
            </p>
          )}
        </div>
      </div>

      {/* Per-config errors (when generation itself failed) */}
      {configErrors && Object.keys(configErrors).length > 0 && (
        <div className="space-y-1.5">
          <p className="text-[10px] font-medium text-muted-foreground">Error Details</p>
          <div className="space-y-1">
            {Object.entries(configErrors).map(([cfgName, errors]) => (
              <div key={cfgName} className="p-2 rounded bg-muted/30 border border-border">
                <span className="text-[10px] font-mono font-medium">{cfgName}</span>
                {errors.map((err, i) => (
                  <p key={i} className="text-[10px] text-destructive/80 mt-0.5 break-all">{err}</p>
                ))}
              </div>
            ))}
          </div>
        </div>
      )}

      {/* Per-design breakdown */}
      {failedDesigns.length > 0 && (
        <div className="space-y-1.5">
          <p className="text-[10px] font-medium text-muted-foreground">Design Results</p>
          <div className="space-y-2">
            {failedDesigns.map((d) => (
              <Card
                key={d.name}
                className={cn(
                  'p-2.5 border-destructive/15',
                  d.pdb_content && onSelectDesign && 'cursor-pointer hover:bg-muted/50 transition-colors',
                )}
                onClick={d.pdb_content && onSelectDesign ? () => onSelectDesign(d.pdb_content!) : undefined}
              >
                <div className="flex items-center justify-between mb-1.5">
                  <div className="flex items-center gap-2">
                    <span className="text-xs font-mono font-medium">{d.name}</span>
                    <Badge variant={getTierBadgeVariant(d.tier)} className="text-[10px] h-4 px-1.5">
                      {d.tier}-tier
                    </Badge>
                    <span className="text-[10px] text-muted-foreground">{d.config_name}</span>
                  </div>
                  {d.pdb_content && onSelectDesign && (
                    <Button
                      variant="ghost"
                      size="sm"
                      className="h-5 px-1.5 text-[10px]"
                      onClick={(e) => { e.stopPropagation(); onSelectDesign(d.pdb_content!); }}
                    >
                      <Eye className="h-3 w-3 mr-1" />
                      View
                    </Button>
                  )}
                </div>
                <div className="grid grid-cols-2 gap-x-4 gap-y-1">
                  <MetricBar label="pTM" value={d.ptm} threshold={0.5} max={1.0} />
                  <MetricBar label="pLDDT" value={d.plddt} threshold={0.7} max={1.0} />
                </div>
                {d.pae > 0 && (
                  <div className="mt-1 text-[10px] text-muted-foreground">
                    pAE: <span className="font-mono">{d.pae.toFixed(1)}</span>
                  </div>
                )}
              </Card>
            ))}
          </div>
        </div>
      )}

      {/* Config rankings */}
      {configRankings && configRankings.length > 0 && (
        <div className="space-y-1">
          <p className="text-[10px] font-medium text-muted-foreground flex items-center gap-1">
            <TrendingUp className="h-3 w-3" />
            Config Performance
          </p>
          <div className="grid gap-1">
            {configRankings.map((cr) => (
              <div key={cr.config} className="flex items-center justify-between text-[10px] px-2 py-1 rounded bg-muted/30">
                <span className="font-mono">{cr.config}</span>
                <div className="flex gap-3 text-muted-foreground">
                  <span>pTM: <span className="font-mono font-medium text-foreground">{cr.avg_ptm.toFixed(3)}</span></span>
                  <span>pLDDT: <span className="font-mono font-medium text-foreground">{cr.avg_plddt.toFixed(3)}</span></span>
                </div>
              </div>
            ))}
          </div>
        </div>
      )}
    </div>
  );
}
