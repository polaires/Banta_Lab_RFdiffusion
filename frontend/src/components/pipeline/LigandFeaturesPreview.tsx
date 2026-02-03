'use client';

import { useState, useCallback, useMemo } from 'react';
import { Card } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { Button } from '@/components/ui/button';
import { Checkbox } from '@/components/ui/checkbox';
import { Progress } from '@/components/ui/progress';
import { Separator } from '@/components/ui/separator';
import { ScrollArea } from '@/components/ui/scroll-area';
import {
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger,
} from '@/components/ui/tooltip';
import {
  Eye,
  Database,
  FlaskConical,
  Atom,
  AlertTriangle,
  Beaker,
  Layers,
} from 'lucide-react';
import { cn } from '@/lib/utils';
import type { StepResult } from '@/lib/pipeline-types';

interface LigandFeaturesPreviewProps {
  result: StepResult;
  onSelectDesign?: (pdbContent: string) => void;
  selectedIds?: string[];
  onSelectionChange?: (ids: string[]) => void;
  isPaused?: boolean;
}

interface FeatureData {
  atom_idx: number;
  atom_name: string;
  element: string;
  type: string;
  is_coordination_donor: boolean;
  coords: [number, number, number] | null;
  hsab: string | null;
  enabled: boolean;
}

// --- Helpers ---

function getSourceBadge(source: string) {
  switch (source) {
    case 'knowledge_base':
      return { label: 'Knowledge Base', variant: 'default' as const, icon: Database };
    case 'chemicalfeatures':
      return { label: 'Predicted', variant: 'secondary' as const, icon: Beaker };
    case 'geometry_filter':
      return { label: 'Geometry', variant: 'outline' as const, icon: Layers };
    default:
      return { label: 'Unknown', variant: 'outline' as const, icon: AlertTriangle };
  }
}

function getTypeBadgeClass(type: string): string {
  switch (type) {
    case 'donor':
      return 'bg-blue-100 text-blue-800 dark:bg-blue-900/40 dark:text-blue-300 border-blue-200 dark:border-blue-800';
    case 'acceptor':
      return 'bg-red-100 text-red-800 dark:bg-red-900/40 dark:text-red-300 border-red-200 dark:border-red-800';
    case 'aromatic':
      return 'bg-purple-100 text-purple-800 dark:bg-purple-900/40 dark:text-purple-300 border-purple-200 dark:border-purple-800';
    case 'hydrophobic':
      return 'bg-green-100 text-green-800 dark:bg-green-900/40 dark:text-green-300 border-green-200 dark:border-green-800';
    default:
      return 'bg-muted text-muted-foreground border-border';
  }
}

function getHsabLabel(hsab: string | null): string {
  if (!hsab || hsab === 'unknown') return '-';
  return hsab.charAt(0).toUpperCase() + hsab.slice(1);
}

// --- Main Component ---

export function LigandFeaturesPreview({
  result,
  onSelectDesign,
}: LigandFeaturesPreviewProps) {
  const data = result.data || {};

  // Handle skipped states
  if (data.skipped) {
    return (
      <div className="flex items-center gap-2 text-xs text-muted-foreground">
        <AlertTriangle className="h-3.5 w-3.5" />
        <span>{data.reason as string || 'Ligand feature analysis skipped'}</span>
      </div>
    );
  }

  const ligandName = data.ligand_name as string || '?';
  const smiles = data.smiles as string || '';
  const source = data.source as string || 'unknown';
  const metal = data.metal as string || null;
  const features = (data.features as FeatureData[]) || [];
  const coordinationDonors = (data.coordination_donors as string[]) || [];
  const maxDenticity = data.max_denticity as number || 0;
  const evidenceCount = data.evidence_count as number || 0;
  const compatibilityScore = data.compatibility_score as number || 0;
  const coordinationMode = data.coordination_mode as string || 'unknown';
  const notes = data.notes as string || '';
  const pdbEvidence = (data.pdb_evidence as string[]) || [];
  const ligandPdbContent = data.ligand_pdb_content as string | null;

  const sourceBadge = getSourceBadge(source);
  const SourceIcon = sourceBadge.icon;

  // Track user overrides for coordination donors
  const [overrides, setOverrides] = useState<Record<string, boolean>>(() => {
    const initial: Record<string, boolean> = {};
    features.forEach(f => {
      if (f.is_coordination_donor) {
        initial[f.atom_name] = f.enabled;
      }
    });
    return initial;
  });

  const handleToggleDonor = useCallback((atomName: string) => {
    setOverrides(prev => {
      const next = { ...prev, [atomName]: !prev[atomName] };
      // Store overrides in result data for downstream steps
      const activeDonors = Object.entries(next)
        .filter(([, enabled]) => enabled)
        .map(([name]) => name);
      if (result.data) {
        (result.data as Record<string, unknown>).user_overrides = {
          active_donors: activeDonors,
          modified: true,
        };
      }
      return next;
    });
  }, [result]);

  const activeDonors = useMemo(() =>
    Object.entries(overrides)
      .filter(([, enabled]) => enabled)
      .map(([name]) => name),
    [overrides],
  );

  const hasOverride = useMemo(() => {
    return features.some(f =>
      f.is_coordination_donor && overrides[f.atom_name] !== f.enabled
    );
  }, [features, overrides]);

  // Separate coordination donors and other features
  const coordFeatures = features.filter(f => f.is_coordination_donor);
  const otherFeatures = features.filter(f => !f.is_coordination_donor);

  return (
    <div className="space-y-3">
      {/* Header row */}
      <div className="flex items-center gap-2 flex-wrap">
        <Badge variant={sourceBadge.variant} className="text-[10px] h-5 gap-1">
          <SourceIcon className="h-3 w-3" />
          {sourceBadge.label}
        </Badge>
        {source === 'knowledge_base' && evidenceCount > 0 && (
          <span className="text-[10px] text-muted-foreground">
            {evidenceCount} PDB structure{evidenceCount !== 1 ? 's' : ''}
          </span>
        )}
        {metal && (
          <Badge variant="outline" className="text-[10px] h-5 gap-1">
            <Atom className="h-3 w-3" />
            {metal}
          </Badge>
        )}
        {compatibilityScore > 0 && (
          <div className="flex items-center gap-1.5 ml-auto">
            <span className="text-[10px] text-muted-foreground">Compatibility</span>
            <Progress value={compatibilityScore * 100} className="h-1.5 w-16" />
            <span className="text-[10px] font-mono font-medium">
              {(compatibilityScore * 100).toFixed(0)}%
            </span>
          </div>
        )}
      </div>

      {/* Feature table */}
      {features.length > 0 && (
        <Card className="p-0 overflow-hidden">
          <ScrollArea className="max-h-64">
            <table className="w-full text-[11px]">
              <thead>
                <tr className="border-b bg-muted/50">
                  <th className="px-2 py-1.5 text-left font-medium w-8"></th>
                  <th className="px-2 py-1.5 text-left font-medium">Atom</th>
                  <th className="px-2 py-1.5 text-left font-medium">Element</th>
                  <th className="px-2 py-1.5 text-left font-medium">Type</th>
                  <th className="px-2 py-1.5 text-left font-medium">HSAB</th>
                  <th className="px-2 py-1.5 text-center font-medium">Coord?</th>
                </tr>
              </thead>
              <tbody>
                {/* Coordination donors first */}
                {coordFeatures.map((f) => (
                  <tr
                    key={`${f.atom_idx}-${f.type}`}
                    className={cn(
                      'border-b transition-colors',
                      overrides[f.atom_name] !== false
                        ? 'bg-primary/5'
                        : 'opacity-50',
                    )}
                  >
                    <td className="px-2 py-1">
                      <Checkbox
                        checked={overrides[f.atom_name] !== false}
                        onCheckedChange={() => handleToggleDonor(f.atom_name)}
                        className="h-3.5 w-3.5"
                      />
                    </td>
                    <td className="px-2 py-1 font-mono font-semibold">{f.atom_name}</td>
                    <td className="px-2 py-1 font-mono">{f.element}</td>
                    <td className="px-2 py-1">
                      <span className={cn(
                        'inline-flex items-center rounded border px-1.5 py-0.5 text-[9px] font-medium',
                        getTypeBadgeClass(f.type),
                      )}>
                        {f.type}
                      </span>
                    </td>
                    <td className="px-2 py-1 text-muted-foreground">{getHsabLabel(f.hsab)}</td>
                    <td className="px-2 py-1 text-center">
                      <span className="text-emerald-600 dark:text-emerald-400 font-medium">Yes</span>
                    </td>
                  </tr>
                ))}
                {/* Non-coordination features */}
                {otherFeatures.map((f) => (
                  <tr
                    key={`${f.atom_idx}-${f.type}`}
                    className="border-b opacity-60"
                  >
                    <td className="px-2 py-1"></td>
                    <td className="px-2 py-1 font-mono">{f.atom_name}</td>
                    <td className="px-2 py-1 font-mono">{f.element}</td>
                    <td className="px-2 py-1">
                      <span className={cn(
                        'inline-flex items-center rounded border px-1.5 py-0.5 text-[9px] font-medium',
                        getTypeBadgeClass(f.type),
                      )}>
                        {f.type}
                      </span>
                    </td>
                    <td className="px-2 py-1 text-muted-foreground">{getHsabLabel(f.hsab)}</td>
                    <td className="px-2 py-1 text-center text-muted-foreground">-</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </ScrollArea>
        </Card>
      )}

      {/* Metadata row */}
      <div className="flex items-center gap-3 flex-wrap text-[10px]">
        {coordinationMode !== 'unknown' && (
          <Badge variant="outline" className="text-[10px] h-5">
            {coordinationMode}
          </Badge>
        )}
        {maxDenticity > 0 && (
          <span className="text-muted-foreground">
            Max denticity: <span className="font-mono font-medium">{maxDenticity}</span>
          </span>
        )}
        {activeDonors.length > 0 && (
          <span className="text-muted-foreground">
            Active donors: <span className="font-mono font-medium">{activeDonors.length}</span>
          </span>
        )}
      </div>

      {/* PDB evidence list */}
      {pdbEvidence.length > 0 && (
        <div className="text-[10px] text-muted-foreground">
          <span className="font-medium">PDB evidence: </span>
          {pdbEvidence.map((ev, i) => (
            <span key={ev}>
              {i > 0 && ', '}
              <span className="font-mono">{ev}</span>
            </span>
          ))}
        </div>
      )}

      {/* Notes */}
      {notes && (
        <p className="text-[10px] text-muted-foreground italic">{notes}</p>
      )}

      {/* View 3D button */}
      {ligandPdbContent && onSelectDesign && (
        <>
          <Separator />
          <Button
            variant="outline"
            size="sm"
            className="gap-1.5"
            onClick={() => onSelectDesign(ligandPdbContent)}
          >
            <Eye className="h-3.5 w-3.5" />
            View 3D
          </Button>
        </>
      )}

      {/* Override indicator */}
      {hasOverride && (
        <div className="flex items-center gap-1.5 text-[10px] text-amber-600 dark:text-amber-400 bg-amber-500/10 rounded px-2 py-1">
          <AlertTriangle className="h-3 w-3 shrink-0" />
          Override applied: {activeDonors.length} active donor{activeDonors.length !== 1 ? 's' : ''}
        </div>
      )}
    </div>
  );
}
