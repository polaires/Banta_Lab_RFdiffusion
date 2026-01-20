'use client';

import { useMemo, useEffect } from 'react';
import { ArrowRight, AlertTriangle } from 'lucide-react';
// Unused icons after UI simplification: Info, Sparkles, Zap, Settings2
import { FormSection } from './FormSection';
import {
  METAL_COORDINATION_INFO,
  getCoordinationWarning,
  isLanthanide,
  getRecommendedDonors,
  recommendPreset,
  getPresetConfig,
  METAL_REPLACEMENT_PRESETS,
  TARGET_METAL_OPTIONS,
} from '../../../lib/enzymeAnalysis';
import type { MetalReplacementPreset } from '../../../lib/enzymeAnalysis';
import type { CoordinationMode, EnzymeAnalysisResult } from '../../../lib/store';

interface MetalReplacementSectionProps {
  /** Whether metal replacement is enabled */
  enabled: boolean;
  /** Toggle metal replacement on/off */
  onEnabledChange: (enabled: boolean) => void;
  /** Source metal from PDB analysis */
  sourceMetal: string | null;
  /** Target metal for replacement */
  targetMetal: string | null;
  /** Callback when target metal changes */
  onTargetMetalChange: (metal: string | null) => void;
  /** Current coordination mode */
  coordinationMode: CoordinationMode;
  /** Callback when coordination mode changes */
  onCoordinationModeChange: (mode: CoordinationMode) => void;
  /** Whether to fix the ligand position */
  fixLigandPosition: boolean;
  /** Callback when fix ligand changes */
  onFixLigandPositionChange: (fix: boolean) => void;
  /** Current preset */
  preset: MetalReplacementPreset;
  /** Callback when preset changes (applies all settings) */
  onPresetChange: (preset: MetalReplacementPreset) => void;
  /** Excluded catalytic residues (metal coordinators in expand/lanthanide modes) */
  excludedResidues: Set<string>;
  /** Enzyme analysis result for metal info */
  enzymeAnalysis: EnzymeAnalysisResult | null;
  /** Detected ligand name (if any) */
  ligandName?: string;
  /** Additional CSS classes */
  className?: string;
}

// Preset display info
// Preset display icons - commented out after UI simplification
// const PRESET_DISPLAY: Record<MetalReplacementPreset, { icon: React.ReactNode; color: string }> = {
//   keep_coordination: { icon: <Settings2 className="h-4 w-4" />, color: 'bg-blue-500' },
//   expand_coordination: { icon: <Zap className="h-4 w-4" />, color: 'bg-orange-500' },
//   reduce_coordination: { icon: <Settings2 className="h-4 w-4" />, color: 'bg-gray-500' },
//   lanthanide_design: { icon: <Sparkles className="h-4 w-4" />, color: 'bg-purple-500' },
//   custom: { icon: <Settings2 className="h-4 w-4" />, color: 'bg-gray-400' },
// };

export function MetalReplacementSection({
  enabled,
  onEnabledChange,
  sourceMetal,
  targetMetal,
  onTargetMetalChange,
  coordinationMode,
  onCoordinationModeChange,
  fixLigandPosition,
  onFixLigandPositionChange,
  preset,
  onPresetChange,
  excludedResidues,
  enzymeAnalysis,
  ligandName,
  className = '',
}: MetalReplacementSectionProps) {
  // Get source metal info
  const sourceMetalInfo = useMemo(() => {
    if (!enzymeAnalysis || enzymeAnalysis.metals.length === 0) return null;
    const metal = enzymeAnalysis.metals[0]; // Use first detected metal
    return {
      element: metal.element,
      cn: metal.coordinationNumber,
      chain: metal.chain,
      residue: metal.residueNumber,
      coordinators: metal.coordinatingAtoms,
    };
  }, [enzymeAnalysis]);

  // Get coordination warning
  const warning = useMemo(() => {
    if (!sourceMetal || !targetMetal || !sourceMetalInfo) return null;
    return getCoordinationWarning(sourceMetal, targetMetal, sourceMetalInfo.cn);
  }, [sourceMetal, targetMetal, sourceMetalInfo]);

  // Get recommended preset when target changes
  const recommendedPreset = useMemo(() => {
    if (!sourceMetal || !targetMetal) return null;
    return recommendPreset(sourceMetal, targetMetal);
  }, [sourceMetal, targetMetal]);

  // Get current preset config
  const presetConfig = useMemo(() => {
    return getPresetConfig(preset, sourceMetal || undefined, targetMetal || undefined);
  }, [preset, sourceMetal, targetMetal]);

  // Get recommended donors for target
  const recommendedDonors = useMemo(() => {
    if (!targetMetal) return [];
    return getRecommendedDonors(targetMetal);
  }, [targetMetal]);

  // Check if target is a lanthanide
  const targetIsLanthanide = targetMetal ? isLanthanide(targetMetal) : false;

  // Calculate CN difference
  const cnDifference = useMemo(() => {
    if (!sourceMetalInfo || !targetMetal) return 0;
    const targetInfo = METAL_COORDINATION_INFO[targetMetal] || METAL_COORDINATION_INFO['DEFAULT'];
    return targetInfo.typicalCN - sourceMetalInfo.cn;
  }, [sourceMetalInfo, targetMetal]);

  // Auto-apply recommended preset when target metal changes
  useEffect(() => {
    if (recommendedPreset && targetMetal && preset === 'custom') {
      // Only auto-apply if currently on custom preset
      onPresetChange(recommendedPreset.preset);
    }
  }, [targetMetal]); // Intentionally not including recommendedPreset/preset/onPresetChange to avoid loops

  if (!sourceMetalInfo) {
    return null; // Don't render if no metal detected
  }

  return (
    <FormSection
      title="Metal Replacement"
      description="Replace the detected metal with a different element"
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
        <span className="text-sm text-foreground">Enable metal replacement</span>
      </label>

      {enabled && (
        <div className="space-y-4 mt-3 pl-6 border-l-2 border-primary/20">
          {/* Source → Target display */}
          <div className="flex items-center gap-4">
            {/* Source metal (read-only) */}
            <div className="flex-1">
              <label className="block text-xs font-medium text-muted-foreground mb-1">
                Source (detected)
              </label>
              <div className="px-3 py-2 bg-muted rounded-lg border border-border">
                <span className="font-mono text-sm font-semibold text-foreground">
                  {sourceMetalInfo.element}
                </span>
                <span className="text-xs text-muted-foreground ml-2">
                  CN {sourceMetalInfo.cn}
                </span>
                <span className="text-xs text-muted-foreground ml-2">
                  ({sourceMetalInfo.chain}{sourceMetalInfo.residue})
                </span>
              </div>
            </div>

            {/* Arrow */}
            <ArrowRight className="h-5 w-5 text-muted-foreground mt-5" />

            {/* Target metal (dropdown) */}
            <div className="flex-1">
              <label className="block text-xs font-medium text-muted-foreground mb-1">
                Target
              </label>
              <select
                value={targetMetal || ''}
                onChange={(e) => onTargetMetalChange(e.target.value || null)}
                className="w-full px-3 py-2 border border-border rounded-lg bg-card text-sm focus:ring-2 focus:ring-ring focus:border-transparent"
              >
                <option value="">Select target metal...</option>
                {Object.entries(TARGET_METAL_OPTIONS).map(([category, metals]) => (
                  <optgroup key={category} label={category}>
                    {metals.map((metal) => (
                      <option key={metal.symbol} value={metal.symbol}>
                        {metal.symbol} - {metal.name} (CN {metal.cn})
                      </option>
                    ))}
                  </optgroup>
                ))}
              </select>
            </div>
          </div>

          {/* Preset Selector - Show when target is selected */}
          {targetMetal && (
            <div className="space-y-2">
              <div className="flex items-center justify-between">
                <label className="block text-xs font-medium text-muted-foreground">
                  Design Preset
                </label>
                {recommendedPreset && preset !== recommendedPreset.preset && (
                  <button
                    type="button"
                    onClick={() => onPresetChange(recommendedPreset.preset)}
                    className="text-xs text-primary hover:underline"
                  >
                    Use recommended
                  </button>
                )}
              </div>

              {/* Preset Cards - shadcn style */}
              <div className="grid grid-cols-2 gap-2">
                {(['keep_coordination', 'expand_coordination', 'lanthanide_design', 'custom'] as MetalReplacementPreset[]).map((p) => {
                  const config = METAL_REPLACEMENT_PRESETS[p];
                  const isRecommended = recommendedPreset?.preset === p;
                  const isSelected = preset === p;

                  return (
                    <button
                      key={p}
                      type="button"
                      onClick={() => onPresetChange(p)}
                      className={`
                        relative px-3 py-2 rounded-xl border-2 text-left transition-all
                        ${isSelected
                          ? 'border-primary bg-primary/5'
                          : 'border-border hover:border-primary/50 bg-card'
                        }
                      `}
                    >
                      {isRecommended && (
                        <span className="absolute -top-2 -right-2 px-1.5 py-0.5 bg-primary text-primary-foreground text-[10px] font-medium rounded-full">
                          ✓
                        </span>
                      )}
                      <span className="text-sm font-medium text-foreground">{config.name}</span>
                    </button>
                  );
                })}
              </div>

              {/* Preset Explanation - commented out for cleaner UI
              {preset !== 'custom' && (
                <div className={`p-3 rounded-lg ${
                  preset === 'lanthanide_design' || preset === 'expand_coordination'
                    ? 'bg-orange-50 dark:bg-orange-950/30 border border-orange-200 dark:border-orange-800'
                    : 'bg-blue-50 dark:bg-blue-950/30 border border-blue-200 dark:border-blue-800'
                }`}>
                  <div className="text-xs space-y-2">
                    {presetConfig.warning && (
                      <div className="flex items-start gap-2 text-orange-800 dark:text-orange-300">
                        <AlertTriangle className="h-3.5 w-3.5 mt-0.5 flex-shrink-0" />
                        <span>{presetConfig.warning}</span>
                      </div>
                    )}
                    <div className="space-y-1 text-muted-foreground">
                      <p><strong>This preset will:</strong></p>
                      <ul className="list-disc list-inside space-y-0.5 ml-2">
                        <li>Coordination mode: <span className="font-mono">{presetConfig.coordinationMode}</span></li>
                        <li>{presetConfig.fixLigandPosition ? 'Fix' : 'Allow movement of'} ligand position</li>
                        <li>{presetConfig.fixMetalCoordinators ? 'Keep' : 'Redesign'} original metal coordinators</li>
                        {presetConfig.buryMetal && <li>Bury {targetMetal} in the binding pocket</li>}
                        {presetConfig.enableHBonds && <li>Enable H-bond conditioning on ligand</li>}
                        {presetConfig.extraCoordinatorsNeeded > 0 && (
                          <li className="text-primary font-medium">
                            RFD3 will design ~{presetConfig.extraCoordinatorsNeeded} additional Asp/Glu coordinators
                          </li>
                        )}
                      </ul>
                    </div>
                  </div>
                </div>
              )}
              */}

              {/* Excluded Coordinators Warning - commented out for cleaner UI
              {excludedResidues.size > 0 && (
                <div className="flex items-start gap-2 p-2 bg-yellow-50 dark:bg-yellow-950/30 rounded-lg border border-yellow-200 dark:border-yellow-800">
                  <Info className="h-4 w-4 text-yellow-600 dark:text-yellow-400 mt-0.5 flex-shrink-0" />
                  <div className="text-xs text-yellow-800 dark:text-yellow-300">
                    <p className="font-medium">Metal coordinators will NOT be fixed:</p>
                    <p className="font-mono mt-1">
                      {Array.from(excludedResidues).slice(0, 6).join(', ')}
                      {excludedResidues.size > 6 && ` +${excludedResidues.size - 6} more`}
                    </p>
                    <p className="mt-1 text-yellow-700 dark:text-yellow-400">
                      These residues will be redesigned for {targetMetal}'s larger coordination sphere.
                    </p>
                  </div>
                </div>
              )}
              */}
            </div>
          )}

          {/* CN Warning (only for custom mode) */}
          {warning && preset === 'custom' && (
            <div className={`flex items-start gap-2 p-3 rounded-lg ${
              warning.severity === 'error'
                ? 'bg-red-50 text-red-800 dark:bg-red-950/30 dark:text-red-300'
                : warning.severity === 'warning'
                ? 'bg-yellow-50 text-yellow-800 dark:bg-yellow-950/30 dark:text-yellow-300'
                : 'bg-blue-50 text-blue-800 dark:bg-blue-950/30 dark:text-blue-300'
            }`}>
              <AlertTriangle className="h-4 w-4 mt-0.5 flex-shrink-0" />
              <p className="text-xs">{warning.warning}</p>
            </div>
          )}

          {/* Manual Coordination Mode (only for custom preset) */}
          {targetMetal && preset === 'custom' && (
            <div className="space-y-2">
              <label className="block text-xs font-medium text-muted-foreground">
                Coordination Mode (Custom)
              </label>
              <div className="flex gap-2">
                {[
                  { value: 'keep', label: 'Keep existing', hint: 'Use current coordinators only' },
                  { value: 'explore', label: 'Explore new', hint: cnDifference > 0 ? `Find ${cnDifference} more donors` : 'Let RFD3 optimize' },
                  { value: 'hybrid', label: 'Hybrid', hint: 'Keep some, explore some' },
                ].map((option) => (
                  <label
                    key={option.value}
                    className={`
                      flex-1 px-3 py-2 rounded-lg border cursor-pointer transition-all
                      ${coordinationMode === option.value
                        ? 'bg-primary/10 border-primary text-primary'
                        : 'bg-card border-border text-foreground hover:border-primary/50'
                      }
                    `}
                  >
                    <input
                      type="radio"
                      name="coordinationMode"
                      value={option.value}
                      checked={coordinationMode === option.value}
                      onChange={(e) => onCoordinationModeChange(e.target.value as CoordinationMode)}
                      className="sr-only"
                    />
                    <div className="text-center">
                      <span className="block text-sm font-medium">{option.label}</span>
                      <span className="block text-[10px] text-muted-foreground mt-0.5">{option.hint}</span>
                    </div>
                  </label>
                ))}
              </div>
            </div>
          )}

          {/* Fix Ligand Position */}
          {ligandName && targetMetal && (
            <label className="flex items-center gap-2 cursor-pointer">
              <input
                type="checkbox"
                checked={fixLigandPosition}
                onChange={(e) => onFixLigandPositionChange(e.target.checked)}
                className="w-4 h-4 rounded border-border text-primary focus:ring-primary"
              />
              <span className="text-sm text-foreground">
                Fix {ligandName} position
              </span>
              <span className="text-xs text-muted-foreground">
                (ligand atoms will not move during diffusion)
              </span>
            </label>
          )}

          {/* Preset Tips - commented out for cleaner UI
          {preset !== 'custom' && presetConfig.tips.length > 0 && (
            <div className="flex items-start gap-2 p-3 bg-primary/5 rounded-lg border border-primary/20">
              <Sparkles className="h-4 w-4 text-primary mt-0.5 flex-shrink-0" />
              <div className="text-xs">
                <p className="font-medium text-foreground">{presetConfig.name} Tips</p>
                <ul className="mt-1 text-muted-foreground space-y-0.5">
                  {presetConfig.tips.map((tip, i) => (
                    <li key={i}>• {tip}</li>
                  ))}
                </ul>
              </div>
            </div>
          )}
          */}
        </div>
      )}
    </FormSection>
  );
}

export default MetalReplacementSection;
