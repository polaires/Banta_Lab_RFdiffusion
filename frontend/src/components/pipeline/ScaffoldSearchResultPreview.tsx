'use client';

import { Card } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { Button } from '@/components/ui/button';
import { Progress } from '@/components/ui/progress';
import {
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger,
} from '@/components/ui/tooltip';
import { ScrollArea } from '@/components/ui/scroll-area';
import {
  Eye,
  Database,
  ArrowRight,
  CheckCircle,
  XCircle,
  AlertTriangle,
  ExternalLink,
  Atom,
  FlaskConical,
} from 'lucide-react';
import { cn } from '@/lib/utils';
import type { StepResult } from '@/lib/pipeline-types';

interface ScaffoldSearchResultPreviewProps {
  result: StepResult;
  onSelectDesign?: (pdbContent: string) => void;
}

interface CandidateData {
  pdb_id: string;
  source_metal: string;
  target_metal: string;
  ligand_code: string;
  needs_substitution: boolean;
  coordination_number: number;
  protein_donors?: string[];
  ligand_donors?: string[];
  ligand_distance?: number;
  resolution?: number;
  score_cn?: number;
  score_hsab?: number;
  score_lig_donors?: number;
  score_prot_donors?: number;
  score_distance?: number;
  score_resolution?: number;
  total_score: number;
}

// Score thresholds for color coding
function getScoreColor(score: number): string {
  if (score >= 70) return 'text-emerald-600 dark:text-emerald-400';
  if (score >= 50) return 'text-amber-600 dark:text-amber-400';
  return 'text-red-500 dark:text-red-400';
}

function getScoreBadgeVariant(score: number): 'default' | 'secondary' | 'destructive' {
  if (score >= 70) return 'default';
  if (score >= 50) return 'secondary';
  return 'destructive';
}

function getScoreLabel(score: number): string {
  if (score >= 80) return 'Excellent';
  if (score >= 70) return 'Good';
  if (score >= 50) return 'Acceptable';
  if (score >= 30) return 'Marginal';
  return 'Poor';
}

// Score breakdown bar component
function ScoreBar({ label, value, max, icon }: {
  label: string;
  value: number;
  max: number;
  icon?: string;
}) {
  const pct = Math.min(100, (value / max) * 100);
  return (
    <div className="space-y-0.5">
      <div className="flex items-center justify-between text-[10px]">
        <span className="text-muted-foreground">{icon ? `${icon} ` : ''}{label}</span>
        <span className="font-mono font-medium">{value.toFixed(1)}/{max}</span>
      </div>
      <Progress value={pct} className="h-1.5" />
    </div>
  );
}

export function ScaffoldSearchResultPreview({
  result,
  onSelectDesign,
}: ScaffoldSearchResultPreviewProps) {
  const data = result.data || {};

  // Handle skipped or error states
  if (data.skipped) {
    return (
      <div className="flex items-center gap-2 text-xs text-muted-foreground">
        <AlertTriangle className="h-3.5 w-3.5" />
        <span>{data.reason as string || 'Scaffold search skipped'}</span>
      </div>
    );
  }

  if (data.search_error) {
    return (
      <div className="flex items-center gap-2 text-xs text-muted-foreground">
        <XCircle className="h-3.5 w-3.5 text-destructive" />
        <span>Search failed: {data.search_error as string}. Continuing with de novo design.</span>
      </div>
    );
  }

  const queryMetal = data.query_metal as string || '?';
  const queryLigand = data.query_ligand as string || '?';
  const ligandCode = data.ligand_code as string || '';
  const numHits = data.num_pdb_hits as number || 0;
  const numValidated = data.num_validated as number || 0;
  const action = data.recommended_action as string || 'de_novo';
  const reason = data.reason as string || '';
  const bestCandidate = data.best_candidate as CandidateData | null;
  const candidates = (data.candidates as CandidateData[]) || [];

  // Get best PDB content from pdbOutputs
  const bestPdbOutput = result.pdbOutputs?.[0];

  const isScaffold = action === 'scaffold';

  return (
    <div className="space-y-3">
      {/* Query summary */}
      <div className="flex items-center gap-2 flex-wrap">
        <Badge variant="outline" className="text-[10px] h-5 gap-1">
          <Atom className="h-3 w-3" />
          {queryMetal}
        </Badge>
        <span className="text-[10px] text-muted-foreground">+</span>
        <Badge variant="outline" className="text-[10px] h-5 gap-1">
          <FlaskConical className="h-3 w-3" />
          {queryLigand}{ligandCode ? ` (${ligandCode})` : ''}
        </Badge>
        <ArrowRight className="h-3 w-3 text-muted-foreground" />
        <Badge
          variant={isScaffold ? 'default' : 'secondary'}
          className="text-[10px] h-5 gap-1"
        >
          {isScaffold ? (
            <><CheckCircle className="h-3 w-3" /> Scaffold found</>
          ) : (
            <><Database className="h-3 w-3" /> De novo</>
          )}
        </Badge>
      </div>

      {/* Search stats */}
      <div className="flex items-center gap-4 text-[10px] text-muted-foreground">
        <span>{numHits} PDB hit{numHits !== 1 ? 's' : ''}</span>
        <span>{numValidated} validated</span>
        {candidates.length > 0 && <span>{candidates.length} scored</span>}
      </div>

      {/* Best candidate card */}
      {bestCandidate && (
        <Card className={cn(
          'p-3 border-2',
          isScaffold ? 'border-primary/30 bg-primary/5' : 'border-border',
        )}>
          <div className="flex items-start justify-between gap-2">
            <div className="space-y-1.5 flex-1">
              {/* PDB ID and score */}
              <div className="flex items-center gap-2">
                <TooltipProvider>
                  <Tooltip>
                    <TooltipTrigger asChild>
                      <a
                        href={`https://www.rcsb.org/structure/${bestCandidate.pdb_id}`}
                        target="_blank"
                        rel="noopener noreferrer"
                        className="text-sm font-bold font-mono text-primary hover:underline flex items-center gap-1"
                      >
                        {bestCandidate.pdb_id}
                        <ExternalLink className="h-3 w-3" />
                      </a>
                    </TooltipTrigger>
                    <TooltipContent>
                      <p>View on RCSB PDB</p>
                    </TooltipContent>
                  </Tooltip>
                </TooltipProvider>

                <Badge variant={getScoreBadgeVariant(bestCandidate.total_score)} className="text-[10px] h-5">
                  Score: {bestCandidate.total_score.toFixed(1)}/100 — {getScoreLabel(bestCandidate.total_score)}
                </Badge>
              </div>

              {/* Metadata row */}
              <div className="flex flex-wrap gap-x-3 gap-y-1 text-[10px] text-muted-foreground">
                <span>
                  Metal: <span className="font-mono font-medium text-foreground">{bestCandidate.source_metal}</span>
                  {bestCandidate.needs_substitution && (
                    <span className="text-amber-500"> → {bestCandidate.target_metal}</span>
                  )}
                </span>
                <span>CN: <span className="font-mono font-medium text-foreground">{bestCandidate.coordination_number}</span></span>
                {bestCandidate.resolution != null && bestCandidate.resolution > 0 && (
                  <span>Res: <span className="font-mono font-medium text-foreground">{bestCandidate.resolution.toFixed(1)}\u00C5</span></span>
                )}
                {bestCandidate.ligand_distance != null && bestCandidate.ligand_distance > 0 && (
                  <span>M-L dist: <span className="font-mono font-medium text-foreground">{bestCandidate.ligand_distance.toFixed(2)}\u00C5</span></span>
                )}
                {bestCandidate.protein_donors && bestCandidate.protein_donors.length > 0 && (
                  <span>Protein donors: <span className="font-mono font-medium text-foreground">{bestCandidate.protein_donors.join(', ')}</span></span>
                )}
                {bestCandidate.ligand_donors && bestCandidate.ligand_donors.length > 0 && (
                  <span>Ligand donors: <span className="font-mono font-medium text-foreground">{bestCandidate.ligand_donors.join(', ')}</span></span>
                )}
              </div>

              {/* Score breakdown */}
              {(bestCandidate.score_cn != null) && (
                <div className="grid grid-cols-2 gap-x-4 gap-y-1 mt-2">
                  <ScoreBar label="CN match" value={bestCandidate.score_cn ?? 0} max={25} />
                  <ScoreBar label="HSAB donors" value={bestCandidate.score_hsab ?? 0} max={20} />
                  <ScoreBar label="Ligand donors" value={bestCandidate.score_lig_donors ?? 0} max={15} />
                  <ScoreBar label="Protein donors" value={bestCandidate.score_prot_donors ?? 0} max={15} />
                  <ScoreBar label="Bond distance" value={bestCandidate.score_distance ?? 0} max={15} />
                  <ScoreBar label="Resolution" value={bestCandidate.score_resolution ?? 0} max={10} />
                </div>
              )}
            </div>

            {/* View button */}
            {bestPdbOutput && onSelectDesign && (
              <Button
                variant="outline"
                size="sm"
                className="shrink-0 gap-1.5"
                onClick={() => onSelectDesign(bestPdbOutput.pdbContent)}
              >
                <Eye className="h-3.5 w-3.5" />
                View 3D
              </Button>
            )}
          </div>

          {bestCandidate.needs_substitution && (
            <div className="mt-2 flex items-center gap-1.5 text-[10px] text-amber-600 dark:text-amber-400 bg-amber-500/10 rounded px-2 py-1">
              <AlertTriangle className="h-3 w-3 shrink-0" />
              Metal substitution needed: {bestCandidate.source_metal} → {bestCandidate.target_metal} (HSAB-compatible)
            </div>
          )}
        </Card>
      )}

      {/* Other candidates */}
      {candidates.length > 1 && (
        <div className="space-y-1.5">
          <p className="text-[10px] font-medium text-muted-foreground">Other candidates</p>
          <ScrollArea className="max-h-32">
            <div className="space-y-1">
              {candidates.slice(1).map((c, i) => {
                // Find matching pdbOutput for this candidate
                const pdbOut = result.pdbOutputs?.find(p => p.id === `scaffold-${c.pdb_id}`);
                return (
                  <div
                    key={c.pdb_id}
                    className={cn(
                      'flex items-center justify-between px-2.5 py-1.5 rounded-md text-[10px] border border-border',
                      pdbOut && onSelectDesign ? 'cursor-pointer hover:bg-muted/50 transition-colors' : '',
                    )}
                    onClick={pdbOut && onSelectDesign ? () => onSelectDesign(pdbOut.pdbContent) : undefined}
                  >
                    <div className="flex items-center gap-3">
                      <span className="text-muted-foreground w-4">{i + 2}.</span>
                      <a
                        href={`https://www.rcsb.org/structure/${c.pdb_id}`}
                        target="_blank"
                        rel="noopener noreferrer"
                        className="font-mono font-medium text-primary hover:underline"
                        onClick={e => e.stopPropagation()}
                      >
                        {c.pdb_id}
                      </a>
                      <span className="text-muted-foreground">
                        {c.source_metal}{c.needs_substitution ? ` → ${c.target_metal}` : ''}
                      </span>
                      <span className="text-muted-foreground">CN {c.coordination_number}</span>
                    </div>
                    <div className="flex items-center gap-2">
                      <span className={cn('font-mono font-medium', getScoreColor(c.total_score))}>
                        {c.total_score.toFixed(1)}
                      </span>
                      {pdbOut && onSelectDesign && (
                        <Eye className="h-3 w-3 text-muted-foreground" />
                      )}
                    </div>
                  </div>
                );
              })}
            </div>
          </ScrollArea>
        </div>
      )}

      {/* Decision reason */}
      {reason && (
        <p className="text-[10px] text-muted-foreground italic">{reason}</p>
      )}
    </div>
  );
}
