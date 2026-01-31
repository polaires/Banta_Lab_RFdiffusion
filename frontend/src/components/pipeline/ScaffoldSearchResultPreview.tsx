'use client';

import { useMemo } from 'react';
import { Card } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { Button } from '@/components/ui/button';
import { Checkbox } from '@/components/ui/checkbox';
import { Progress } from '@/components/ui/progress';
import { Separator } from '@/components/ui/separator';
import {
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger,
} from '@/components/ui/tooltip';
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
  // Selection props (passed from StepCard when supportsSelection + isPaused)
  selectedIds?: string[];
  onSelectionChange?: (ids: string[]) => void;
  isPaused?: boolean;
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

// --- Helpers ---

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

function ScoreBar({ label, value, max }: { label: string; value: number; max: number }) {
  const pct = Math.min(100, (value / max) * 100);
  return (
    <div className="space-y-0.5">
      <div className="flex items-center justify-between text-[10px]">
        <span className="text-muted-foreground">{label}</span>
        <span className="font-mono font-medium">{value.toFixed(1)}/{max}</span>
      </div>
      <Progress value={pct} className="h-1.5" />
    </div>
  );
}

/** Parse a protein donor string like "C26:SER O@2.29\u00C5" into parts */
function parseProteinDonor(s: string): { residue: string; atom: string; distance: string } | null {
  // Format: "C26:SER O@2.29\u00C5" or "C28:ASP OD1@2.28\u00C5"
  const match = s.match(/^([A-Z]?\d+:\w+)\s+(\w+)@([\d.]+)/);
  if (match) return { residue: match[1], atom: match[2], distance: match[3] };
  return null;
}

/** Parse a ligand donor string like "O1@2.26\u00C5" into parts */
function parseLigandDonor(s: string): { atom: string; distance: string } | null {
  const match = s.match(/^(\w+)@([\d.]+)/);
  if (match) return { atom: match[1], distance: match[2] };
  return null;
}

/** Classify donor distance for color */
function getDonorDistColor(dist: number): string {
  if (dist <= 2.5) return 'bg-emerald-100 text-emerald-800 dark:bg-emerald-900/40 dark:text-emerald-300 border-emerald-200 dark:border-emerald-800';
  if (dist <= 3.0) return 'bg-amber-100 text-amber-800 dark:bg-amber-900/40 dark:text-amber-300 border-amber-200 dark:border-amber-800';
  return 'bg-muted text-muted-foreground border-border';
}

// --- Detail card for selected candidate ---

function CandidateDetailCard({ candidate, pdbContent, onSelectDesign }: {
  candidate: CandidateData;
  pdbContent?: string;
  onSelectDesign?: (pdb: string) => void;
}) {
  const proteinDonors = useMemo(() =>
    (candidate.protein_donors || []).map(parseProteinDonor).filter(Boolean),
    [candidate.protein_donors]
  );
  const ligandDonors = useMemo(() =>
    (candidate.ligand_donors || []).map(parseLigandDonor).filter(Boolean),
    [candidate.ligand_donors]
  );

  return (
    <div className="space-y-3">
      {/* Header: PDB ID + Score + View 3D */}
      <div className="flex items-start justify-between gap-2">
        <div className="flex items-center gap-2 flex-wrap">
          <TooltipProvider>
            <Tooltip>
              <TooltipTrigger asChild>
                <a
                  href={`https://www.rcsb.org/structure/${candidate.pdb_id}`}
                  target="_blank"
                  rel="noopener noreferrer"
                  className="text-sm font-bold font-mono text-primary hover:underline flex items-center gap-1"
                >
                  {candidate.pdb_id}
                  <ExternalLink className="h-3 w-3" />
                </a>
              </TooltipTrigger>
              <TooltipContent><p>View on RCSB PDB</p></TooltipContent>
            </Tooltip>
          </TooltipProvider>
          <Badge variant={getScoreBadgeVariant(candidate.total_score)} className="text-[10px] h-5">
            {candidate.total_score.toFixed(1)}/100 — {getScoreLabel(candidate.total_score)}
          </Badge>
        </div>
        {pdbContent && onSelectDesign && (
          <Button variant="outline" size="sm" className="shrink-0 gap-1.5" onClick={() => onSelectDesign(pdbContent)}>
            <Eye className="h-3.5 w-3.5" />
            View 3D
          </Button>
        )}
      </div>

      {/* Key metrics grid */}
      <div className="grid grid-cols-2 sm:grid-cols-4 gap-2">
        <div className="rounded-md border border-border px-2.5 py-1.5">
          <div className="text-[9px] uppercase tracking-wider text-muted-foreground">Metal</div>
          <div className="text-xs font-mono font-semibold mt-0.5">
            {candidate.source_metal}
            {candidate.needs_substitution && (
              <span className="text-amber-500"> {'\u2192'} {candidate.target_metal}</span>
            )}
          </div>
        </div>
        <div className="rounded-md border border-border px-2.5 py-1.5">
          <div className="text-[9px] uppercase tracking-wider text-muted-foreground">Coordination</div>
          <div className="text-xs font-mono font-semibold mt-0.5">CN {candidate.coordination_number}</div>
        </div>
        {candidate.resolution != null && candidate.resolution > 0 && (
          <div className="rounded-md border border-border px-2.5 py-1.5">
            <div className="text-[9px] uppercase tracking-wider text-muted-foreground">Resolution</div>
            <div className="text-xs font-mono font-semibold mt-0.5">{candidate.resolution.toFixed(1)} {'\u00C5'}</div>
          </div>
        )}
        {candidate.ligand_distance != null && candidate.ligand_distance > 0 && (
          <div className="rounded-md border border-border px-2.5 py-1.5">
            <div className="text-[9px] uppercase tracking-wider text-muted-foreground">M-L Distance</div>
            <div className="text-xs font-mono font-semibold mt-0.5">{candidate.ligand_distance.toFixed(2)} {'\u00C5'}</div>
          </div>
        )}
      </div>

      {/* Protein donors */}
      {proteinDonors.length > 0 && (
        <div>
          <div className="text-[10px] font-medium text-muted-foreground mb-1">Protein Donors</div>
          <div className="flex flex-wrap gap-1">
            {proteinDonors.map((d, i) => {
              const dist = parseFloat(d!.distance);
              return (
                <span
                  key={i}
                  className={cn('inline-flex items-center gap-1 rounded-md border px-1.5 py-0.5 text-[10px] font-mono', getDonorDistColor(dist))}
                >
                  <span className="font-semibold">{d!.residue}</span>
                  <span className="opacity-70">{d!.atom}</span>
                  <span className="opacity-60">{d!.distance}{'\u00C5'}</span>
                </span>
              );
            })}
          </div>
        </div>
      )}

      {/* Ligand donors */}
      {ligandDonors.length > 0 && (
        <div>
          <div className="text-[10px] font-medium text-muted-foreground mb-1">Ligand Donors</div>
          <div className="flex flex-wrap gap-1">
            {ligandDonors.map((d, i) => {
              const dist = parseFloat(d!.distance);
              return (
                <span
                  key={i}
                  className={cn('inline-flex items-center gap-1 rounded-md border px-1.5 py-0.5 text-[10px] font-mono', getDonorDistColor(dist))}
                >
                  <span className="font-semibold">{d!.atom}</span>
                  <span className="opacity-60">{d!.distance}{'\u00C5'}</span>
                </span>
              );
            })}
          </div>
        </div>
      )}

      {/* Metal substitution warning */}
      {candidate.needs_substitution && (
        <div className="flex items-center gap-1.5 text-[10px] text-amber-600 dark:text-amber-400 bg-amber-500/10 rounded px-2 py-1">
          <AlertTriangle className="h-3 w-3 shrink-0" />
          Metal substitution: {candidate.source_metal} {'\u2192'} {candidate.target_metal} (HSAB-compatible)
        </div>
      )}

      {/* Score breakdown */}
      {candidate.score_cn != null && (
        <div className="grid grid-cols-2 gap-x-4 gap-y-1">
          <ScoreBar label="CN match" value={candidate.score_cn ?? 0} max={25} />
          <ScoreBar label="HSAB donors" value={candidate.score_hsab ?? 0} max={20} />
          <ScoreBar label="Ligand donors" value={candidate.score_lig_donors ?? 0} max={20} />
          <ScoreBar label="Protein donors" value={candidate.score_prot_donors ?? 0} max={10} />
          <ScoreBar label="Bond distance" value={candidate.score_distance ?? 0} max={15} />
          <ScoreBar label="Resolution" value={candidate.score_resolution ?? 0} max={10} />
        </div>
      )}
    </div>
  );
}

// --- Main component ---

export function ScaffoldSearchResultPreview({
  result,
  onSelectDesign,
  selectedIds,
  onSelectionChange,
  isPaused,
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
  const candidates = (data.candidates as CandidateData[]) || [];
  const isScaffold = action === 'scaffold';

  // Determine which candidate is selected (for detail display)
  const selectionActive = isPaused && selectedIds && onSelectionChange;
  const selectedPdbId = selectedIds && selectedIds.length > 0
    ? selectedIds[0].replace('scaffold-', '')
    : null;
  const displayCandidate = candidates.find(c => c.pdb_id === selectedPdbId) || candidates[0] || null;
  const displayPdbOutput = displayCandidate
    ? result.pdbOutputs?.find(p => p.id === `scaffold-${displayCandidate.pdb_id}`)
    : result.pdbOutputs?.[0];

  const handleToggle = (pdbId: string) => {
    if (!onSelectionChange) return;
    const outputId = `scaffold-${pdbId}`;
    // Radio-style: select only this one
    onSelectionChange([outputId]);
  };

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

      {/* Selected candidate detail card */}
      {displayCandidate && (
        <Card className={cn(
          'p-3 border-2',
          isScaffold ? 'border-primary/30 bg-primary/5' : 'border-border',
        )}>
          <CandidateDetailCard
            candidate={displayCandidate}
            pdbContent={displayPdbOutput?.pdbContent}
            onSelectDesign={onSelectDesign}
          />
        </Card>
      )}

      {/* Candidate selection list (consolidates "Other candidates" + "Select outputs") */}
      {candidates.length > 0 && (
        <>
          <Separator />
          <div className="space-y-1.5">
            <div className="flex items-center justify-between">
              <p className="text-xs font-medium text-foreground">
                {selectionActive ? 'Select scaffold to use' : `Candidates (${candidates.length})`}
              </p>
              {selectionActive && (
                <p className="text-[10px] text-muted-foreground">
                  Click to select a different scaffold
                </p>
              )}
            </div>
            <div className="max-h-52 overflow-y-auto rounded-md space-y-1">
              {candidates.map((c, i) => {
                const pdbOut = result.pdbOutputs?.find(p => p.id === `scaffold-${c.pdb_id}`);
                const isSelected = selectedPdbId === c.pdb_id;

                return (
                  <div
                    key={c.pdb_id}
                    className={cn(
                      'flex items-center gap-2 px-2.5 py-2 rounded-md text-[11px] border transition-all',
                      isSelected
                        ? 'border-primary bg-primary/5'
                        : 'border-border hover:border-primary/50 hover:bg-muted/50',
                      selectionActive ? 'cursor-pointer' : (pdbOut && onSelectDesign ? 'cursor-pointer' : ''),
                    )}
                    onClick={() => {
                      if (selectionActive) {
                        handleToggle(c.pdb_id);
                      } else if (pdbOut && onSelectDesign) {
                        onSelectDesign(pdbOut.pdbContent);
                      }
                    }}
                  >
                    {/* Selection indicator */}
                    {selectionActive ? (
                      <Checkbox
                        checked={isSelected}
                        onCheckedChange={() => handleToggle(c.pdb_id)}
                        className="shrink-0"
                      />
                    ) : (
                      <span className="text-muted-foreground w-4 text-center shrink-0 text-[10px]">{i + 1}.</span>
                    )}

                    {/* PDB ID */}
                    <a
                      href={`https://www.rcsb.org/structure/${c.pdb_id}`}
                      target="_blank"
                      rel="noopener noreferrer"
                      className="font-mono font-semibold text-primary hover:underline shrink-0"
                      onClick={e => e.stopPropagation()}
                    >
                      {c.pdb_id}
                    </a>

                    {/* Metal info */}
                    <span className="text-muted-foreground font-mono shrink-0">
                      {c.source_metal}{c.needs_substitution ? ` \u2192 ${c.target_metal}` : ''}
                    </span>

                    {/* CN */}
                    <Badge variant="outline" className="text-[9px] h-4 px-1.5 shrink-0">
                      CN {c.coordination_number}
                    </Badge>

                    {/* Resolution */}
                    {c.resolution != null && c.resolution > 0 && (
                      <span className="text-[10px] text-muted-foreground shrink-0">
                        {c.resolution.toFixed(1)}{'\u00C5'}
                      </span>
                    )}

                    {/* Spacer */}
                    <span className="flex-1" />

                    {/* Score */}
                    <span className={cn('font-mono font-semibold shrink-0', getScoreColor(c.total_score))}>
                      {c.total_score.toFixed(1)}
                    </span>

                    {/* View icon */}
                    {pdbOut && onSelectDesign && !selectionActive && (
                      <Eye className="h-3 w-3 text-muted-foreground shrink-0" />
                    )}
                  </div>
                );
              })}
            </div>
          </div>
        </>
      )}

      {/* Decision reason */}
      {reason && (
        <p className="text-[10px] text-muted-foreground italic">{reason}</p>
      )}
    </div>
  );
}
