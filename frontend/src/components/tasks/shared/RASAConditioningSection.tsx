'use client';

import { useMemo, useCallback } from 'react';
import { AlertCircle } from 'lucide-react';
import { FormSection } from './FormSection';
import { AtomCheckboxGrid } from './AtomCheckboxGrid';
import type { EnzymeAnalysisResult } from '../../../lib/store';

interface RASAConditioningSectionProps {
  /** Whether RASA conditioning is enabled */
  enabled: boolean;
  /** Toggle RASA conditioning on/off */
  onEnabledChange: (enabled: boolean) => void;
  /** Selected buried atoms: { ligandName: "O1,O2,O3" } */
  buriedAtoms: Record<string, string>;
  /** Callback when buried atoms change */
  onBuriedAtomsChange: (key: string, atoms: string) => void;
  /** Selected exposed atoms: { ligandName: "O4,O5,O6" } */
  exposedAtoms: Record<string, string>;
  /** Callback when exposed atoms change */
  onExposedAtomsChange: (key: string, atoms: string) => void;
  /** Whether buried selection has been manually overridden */
  buriedOverridden: boolean;
  /** Whether exposed selection has been manually overridden */
  exposedOverridden: boolean;
  /** Enzyme analysis result for ligand info */
  enzymeAnalysis: EnzymeAnalysisResult | null;
  /** Target metal (if any) for burial suggestion */
  targetMetal?: string | null;
  /** Apply auto suggestions */
  onApplySuggestions: () => void;
  /** Additional CSS classes */
  className?: string;
}

export function RASAConditioningSection({
  enabled,
  onEnabledChange,
  buriedAtoms,
  onBuriedAtomsChange,
  exposedAtoms,
  onExposedAtomsChange,
  buriedOverridden,
  exposedOverridden,
  enzymeAnalysis,
  targetMetal,
  onApplySuggestions,
  className = '',
}: RASAConditioningSectionProps) {
  // Get ligand info
  const ligandInfo = useMemo(() => {
    if (!enzymeAnalysis || enzymeAnalysis.ligands.length === 0) return null;
    return enzymeAnalysis.ligands[0]; // Use first detected ligand
  }, [enzymeAnalysis]);

  // Get metal info
  const metalInfo = useMemo(() => {
    if (!enzymeAnalysis || enzymeAnalysis.metals.length === 0) return null;
    return enzymeAnalysis.metals[0];
  }, [enzymeAnalysis]);

  // Convert comma-separated string to array
  const parseAtomString = (str: string): string[] => {
    return str ? str.split(',').map(s => s.trim()).filter(Boolean) : [];
  };

  // Convert array to comma-separated string
  const toAtomString = (atoms: string[]): string => {
    return atoms.join(',');
  };

  // Get buried atoms for ligand as array
  const ligandBuriedAtoms = useMemo(() => {
    if (!ligandInfo) return [];
    return parseAtomString(buriedAtoms[ligandInfo.name] || '');
  }, [ligandInfo, buriedAtoms]);

  // Get exposed atoms for ligand as array
  const ligandExposedAtoms = useMemo(() => {
    if (!ligandInfo) return [];
    return parseAtomString(exposedAtoms[ligandInfo.name] || '');
  }, [ligandInfo, exposedAtoms]);

  // Handle buried atom selection change
  const handleBuriedChange = useCallback((atoms: string[]) => {
    if (!ligandInfo) return;
    onBuriedAtomsChange(ligandInfo.name, toAtomString(atoms));
  }, [ligandInfo, onBuriedAtomsChange]);

  // Handle exposed atom selection change
  const handleExposedChange = useCallback((atoms: string[]) => {
    if (!ligandInfo) return;
    onExposedAtomsChange(ligandInfo.name, toAtomString(atoms));
  }, [ligandInfo, onExposedAtomsChange]);

  // Apply buried suggestions
  const applyBuriedSuggestions = useCallback(() => {
    if (!ligandInfo) return;
    onBuriedAtomsChange(ligandInfo.name, toAtomString(ligandInfo.suggestedBuried));
  }, [ligandInfo, onBuriedAtomsChange]);

  // Apply exposed suggestions
  const applyExposedSuggestions = useCallback(() => {
    if (!ligandInfo) return;
    onExposedAtomsChange(ligandInfo.name, toAtomString(ligandInfo.suggestedExposed));
  }, [ligandInfo, onExposedAtomsChange]);

  // Check for overlap between buried and exposed
  const overlapError = useMemo(() => {
    const buriedSet = new Set(ligandBuriedAtoms);
    const overlapping = ligandExposedAtoms.filter(a => buriedSet.has(a));
    return overlapping.length > 0 ? overlapping : null;
  }, [ligandBuriedAtoms, ligandExposedAtoms]);

  // Metal display name
  const metalDisplay = targetMetal || metalInfo?.element || null;

  if (!ligandInfo && !metalInfo) {
    return null; // Don't render if no ligand or metal detected
  }

  return (
    <FormSection
      title="Active Site Accessibility (RASA)"
      description="Control burial and exposure of active site atoms"
      className={className}
    >
      {/* Enable toggle */}
      <label className="flex items-center gap-2 cursor-pointer">
        <input
          type="checkbox"
          checked={enabled}
          onChange={(e) => onEnabledChange(e.target.checked)}
          className="w-4 h-4 rounded border-border text-primary focus:ring-primary"
        />
        <span className="text-sm text-foreground">Enable RASA conditioning</span>
      </label>

      {enabled && (
        <div className="space-y-4 mt-3">
          {/* Ligand header */}
          {ligandInfo && (
            <div className="flex items-center gap-2 text-sm">
              <span className="font-medium text-foreground">Ligand:</span>
              <span className="font-mono bg-muted px-2 py-0.5 rounded text-foreground">
                {ligandInfo.name}
              </span>
              <span className="text-muted-foreground">
                ({ligandInfo.atoms.length} atoms detected)
              </span>
            </div>
          )}

          {/* Overlap warning */}
          {overlapError && (
            <div className="flex items-start gap-2 p-3 bg-red-50 text-red-800 dark:bg-red-950/30 dark:text-red-300 rounded-lg">
              <AlertCircle className="h-4 w-4 mt-0.5 flex-shrink-0" />
              <div className="text-xs">
                <p className="font-medium">Conflicting selections</p>
                <p>
                  {overlapError.join(', ')} cannot be both buried and exposed.
                  Please remove from one category.
                </p>
              </div>
            </div>
          )}

          {/* Buried atoms section */}
          <div className="p-3 bg-muted/30 rounded-lg border border-border space-y-3">
            <div className="flex items-center justify-between">
              <h4 className="text-sm font-medium text-foreground">
                Buried (enclosed pocket)
              </h4>
            </div>

            {/* Metal burial */}
            {metalDisplay && (
              <label className="flex items-center gap-2 cursor-pointer">
                <input
                  type="checkbox"
                  checked={buriedAtoms[metalDisplay] === 'ALL'}
                  onChange={(e) => {
                    if (e.target.checked) {
                      onBuriedAtomsChange(metalDisplay, 'ALL');
                    } else {
                      onBuriedAtomsChange(metalDisplay, '');
                    }
                  }}
                  className="w-4 h-4 rounded border-border text-primary focus:ring-primary"
                />
                <span className="text-sm text-foreground">
                  Metal: <span className="font-mono font-semibold">{metalDisplay}</span>
                </span>
                <span className="text-xs text-muted-foreground">(always bury the metal ion)</span>
              </label>
            )}

            {/* Ligand buried atoms */}
            {ligandInfo && (
              <AtomCheckboxGrid
                label="Ligand atoms:"
                atoms={ligandInfo.atoms}
                selectedAtoms={ligandBuriedAtoms}
                onChange={handleBuriedChange}
                suggestedAtoms={ligandInfo.suggestedBuried}
                suggestionButtonLabel="Apply suggested burial"
                onApplySuggestions={applyBuriedSuggestions}
                isOverridden={buriedOverridden}
                // hint="Coordination atoms (near metal) should be buried" - commented out for cleaner UI
              />
            )}
          </div>

          {/* Exposed atoms section */}
          <div className="p-3 bg-muted/30 rounded-lg border border-border space-y-3">
            <div className="flex items-center justify-between">
              <h4 className="text-sm font-medium text-foreground">
                Exposed (substrate entry/exit)
              </h4>
            </div>

            {/* Ligand exposed atoms */}
            {ligandInfo && (
              <AtomCheckboxGrid
                label="Ligand atoms:"
                atoms={ligandInfo.atoms}
                selectedAtoms={ligandExposedAtoms}
                onChange={handleExposedChange}
                suggestedAtoms={ligandInfo.suggestedExposed}
                suggestionButtonLabel="Apply suggested exposure"
                onApplySuggestions={applyExposedSuggestions}
                isOverridden={exposedOverridden}
                // hint="Terminal atoms allow substrate access" - commented out for cleaner UI
              />
            )}
          </div>

          {/* Apply all suggestions button - removed, using inline Auto buttons instead */}

          {/* Info box - commented out for cleaner UI
          <div className="p-3 bg-blue-50 dark:bg-blue-950/30 rounded-lg border border-blue-200 dark:border-blue-800">
            <p className="text-xs text-blue-800 dark:text-blue-300">
              <strong>How RASA works:</strong> Buried atoms will be enclosed by protein (low solvent accessibility),
              while exposed atoms will remain accessible to solvent. This shapes the active site pocket geometry.
            </p>
          </div>
          */}
        </div>
      )}
    </FormSection>
  );
}

export default RASAConditioningSection;
