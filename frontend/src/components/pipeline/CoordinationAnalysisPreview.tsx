'use client';

import { useState, useMemo } from 'react';
import { Card } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { Button } from '@/components/ui/button';
import { Checkbox } from '@/components/ui/checkbox';
import { Separator } from '@/components/ui/separator';
import { ScrollArea } from '@/components/ui/scroll-area';
import { Atom, Eye, AlertTriangle, Database, Beaker } from 'lucide-react';
import { cn } from '@/lib/utils';
import { useStore } from '@/lib/store';
import type { StepResult } from '@/lib/pipeline-types';
import type { LigandCoordinationFeature, ProteinDonor } from '@/lib/store';

interface CoordinationAnalysisPreviewProps {
  result: StepResult;
  onSelectDesign?: (pdbContent: string) => void;
}

function getHsabLabel(hsab: string | null | undefined): string {
  if (!hsab || hsab === 'unknown') return '-';
  return hsab.charAt(0).toUpperCase() + hsab.slice(1);
}

export function CoordinationAnalysisPreview({ result, onSelectDesign }: CoordinationAnalysisPreviewProps) {
  const setCoordinationSphereData = useStore(s => s.setCoordinationSphereData);
  const setLigandCoordinationFeatures = useStore(s => s.setLigandCoordinationFeatures);
  const setShouldFocusCoordinationSphere = useStore(s => s.setShouldFocusCoordinationSphere);

  const data = result.data || {};

  console.log('[CoordinationAnalysisPreview] result.data:', data);
  console.log('[CoordinationAnalysisPreview] protein_donors:', data.protein_donors);
  console.log('[CoordinationAnalysisPreview] ligand_donors:', data.ligand_donors);

  // Handle skipped state
  if (data.skipped) {
    return (
      <div className="flex items-center gap-2 text-xs text-muted-foreground">
        <AlertTriangle className="h-3.5 w-3.5" />
        <span>{(data.reason as string) || 'Coordination analysis skipped'}</span>
      </div>
    );
  }

  const proteinDonors = (data.protein_donors as ProteinDonor[]) || [];
  const ligandDonors = (data.ligand_donors as LigandCoordinationFeature[]) || [];
  const metal = (data.metal as string) || '';
  const ligandName = (data.ligand_name as string) || '';
  const scaffoldPdbId = data.scaffold_pdb_id as string | undefined;
  const coordinationNumber = (data.coordination_number as number) || 0;
  const hasScaffoldData = data.has_scaffold_data as boolean;
  const hasLigandKbData = data.has_ligand_kb_data as boolean;

  // Track enabled state for protein donors (key includes atomName for uniqueness)
  const [enabledProtein, setEnabledProtein] = useState<Record<string, boolean>>(() => {
    const init: Record<string, boolean> = {};
    proteinDonors.forEach(d => {
      init[`${d.chain}${d.residue}_${d.atomName}`] = d.enabled !== false;
    });
    return init;
  });

  // Track enabled state for ligand donors
  const [enabledLigand, setEnabledLigand] = useState<Record<string, boolean>>(() => {
    const init: Record<string, boolean> = {};
    ligandDonors.forEach(d => {
      init[d.atom_name] = d.enabled !== false;
    });
    return init;
  });

  const activeProteinCount = useMemo(() =>
    Object.values(enabledProtein).filter(Boolean).length,
    [enabledProtein]);

  const activeLigandCount = useMemo(() =>
    Object.values(enabledLigand).filter(Boolean).length,
    [enabledLigand]);

  const handleView3D = () => {
    console.log('[CoordinationAnalysisPreview] handleView3D clicked');

    // Filter enabled donors
    const activeLigandDonors = ligandDonors
      .filter(d => enabledLigand[d.atom_name])
      .map(d => ({ ...d, enabled: true }));

    const activeProteinDonors = proteinDonors
      .filter(d => enabledProtein[`${d.chain}${d.residue}_${d.atomName}`]);

    // Set ligand features in store (used by focusOnCoordinationSphere)
    setLigandCoordinationFeatures(activeLigandDonors);

    // Set full coordination data in store (used by focusOnCoordinationSphere)
    setCoordinationSphereData({
      metal,
      ligandName,
      coordinationNumber,
      scaffoldPdbId,
      ligandDonors: activeLigandDonors,
      proteinDonors: activeProteinDonors,
    });

    // Store user overrides in result data for configure step
    if (result.data) {
      (result.data as Record<string, unknown>).user_overrides = {
        active_protein_donors: activeProteinDonors.map(d => `${d.chain}${d.residue}`),
        active_ligand_donors: activeLigandDonors.map(d => d.atom_name),
        modified: true,
      };
    }

    // Load scaffold PDB and trigger focus view
    const scaffoldPdb = result.pdbOutputs?.find(p => p.id === 'coordination-scaffold')?.pdbContent;
    if (scaffoldPdb && onSelectDesign) {
      console.log('[CoordinationAnalysisPreview] Loading scaffold PDB into viewer');
      onSelectDesign(scaffoldPdb);

      // Trigger proper focus view after structure loads
      setTimeout(() => {
        setShouldFocusCoordinationSphere(true);
      }, 300);
    } else {
      console.warn('[CoordinationAnalysisPreview] Cannot load scaffold - no PDB content or onSelectDesign');
    }
  };

  return (
    <div className="space-y-4">
      {/* Header */}
      <div className="flex items-center gap-2 flex-wrap">
        {metal && (
          <Badge variant="default" className="gap-1">
            <Atom className="h-3 w-3" />
            {metal}
          </Badge>
        )}
        {scaffoldPdbId && (
          <Badge variant="outline" className="gap-1">
            <Database className="h-3 w-3" />
            {scaffoldPdbId}
          </Badge>
        )}
        {hasLigandKbData && (
          <Badge variant="secondary" className="gap-1 text-[10px]">
            <Beaker className="h-3 w-3" />
            KB
          </Badge>
        )}
        {coordinationNumber > 0 && (
          <span className="text-xs text-muted-foreground">
            CN: <span className="font-mono font-medium">{coordinationNumber}</span>
          </span>
        )}
      </div>

      {/* Two-panel layout */}
      <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
        {/* Protein Donors - residues coordinating to metal */}
        <Card className="p-0 overflow-hidden">
          <div className="px-3 py-2 bg-muted/50 border-b">
            <h4 className="text-xs font-medium">Metal Coordination</h4>
            <p className="text-[10px] text-muted-foreground">
              {hasScaffoldData ? 'Protein residues → metal' : 'Coordinating residues'}
            </p>
          </div>
          <ScrollArea className="max-h-48">
            {proteinDonors.length > 0 ? (
              <table className="w-full text-[11px]">
                <thead>
                  <tr className="border-b">
                    <th className="px-2 py-1 w-6"></th>
                    <th className="px-2 py-1 text-left">Residue</th>
                    <th className="px-2 py-1 text-left">Atom</th>
                    <th className="px-2 py-1 text-right">Dist</th>
                  </tr>
                </thead>
                <tbody>
                  {proteinDonors.map(d => {
                    const key = `${d.chain}${d.residue}_${d.atomName}`;
                    return (
                      <tr key={key} className={cn('border-b', !enabledProtein[key] && 'opacity-50')}>
                        <td className="px-2 py-1">
                          <Checkbox
                            checked={enabledProtein[key]}
                            onCheckedChange={() => setEnabledProtein(prev => ({ ...prev, [key]: !prev[key] }))}
                            className="h-3 w-3"
                          />
                        </td>
                        <td className="px-2 py-1 font-mono">{d.name} {d.chain}{d.residue}</td>
                        <td className="px-2 py-1 font-mono">{d.atomName}</td>
                        <td className="px-2 py-1 text-right font-mono">
                          {d.distance?.toFixed(2) || '-'}
                        </td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            ) : (
              <p className="px-3 py-4 text-xs text-muted-foreground text-center">
                No scaffold data
              </p>
            )}
          </ScrollArea>
        </Card>

        {/* Ligand Donors - ligand atoms coordinating to metal */}
        <Card className="p-0 overflow-hidden">
          <div className="px-3 py-2 bg-muted/50 border-b">
            <h4 className="text-xs font-medium">Ligand → Metal</h4>
            <p className="text-[10px] text-muted-foreground">{ligandName || 'Coordinating atoms'}</p>
          </div>
          <ScrollArea className="max-h-48">
            {ligandDonors.length > 0 ? (
              <table className="w-full text-[11px]">
                <thead>
                  <tr className="border-b">
                    <th className="px-2 py-1 w-6"></th>
                    <th className="px-2 py-1 text-left">Atom</th>
                    <th className="px-2 py-1 text-left">Type</th>
                    <th className="px-2 py-1 text-left">HSAB</th>
                  </tr>
                </thead>
                <tbody>
                  {ligandDonors.map(d => (
                    <tr key={d.atom_name} className={cn('border-b', !enabledLigand[d.atom_name] && 'opacity-50')}>
                      <td className="px-2 py-1">
                        <Checkbox
                          checked={enabledLigand[d.atom_name]}
                          onCheckedChange={() => setEnabledLigand(prev => ({ ...prev, [d.atom_name]: !prev[d.atom_name] }))}
                          className="h-3 w-3"
                        />
                      </td>
                      <td className="px-2 py-1 font-mono font-semibold">{d.atom_name}</td>
                      <td className="px-2 py-1">{d.type}</td>
                      <td className="px-2 py-1 text-muted-foreground">{getHsabLabel(d.hsab)}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            ) : (
              <p className="px-3 py-4 text-xs text-muted-foreground text-center">
                No ligand features
              </p>
            )}
          </ScrollArea>
        </Card>
      </div>

      {/* Selection summary */}
      <div className="flex items-center gap-3 text-[10px] text-muted-foreground">
        <span>
          Protein: <span className="font-mono font-medium">{activeProteinCount}</span> selected
        </span>
        <span>
          Ligand: <span className="font-mono font-medium">{activeLigandCount}</span> selected
        </span>
      </div>

      {/* View 3D button */}
      {(proteinDonors.length > 0 || ligandDonors.length > 0) && onSelectDesign && (
        <>
          <Separator />
          <Button variant="outline" size="sm" className="gap-1.5" onClick={handleView3D}>
            <Eye className="h-3.5 w-3.5" />
            View Coordination Sphere
          </Button>
        </>
      )}
    </div>
  );
}
