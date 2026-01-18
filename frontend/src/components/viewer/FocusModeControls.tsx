'use client';

import { useState } from 'react';
import {
  Droplet,
  Link2,
  CircleDot,
  SlidersHorizontal,
  X,
  ChevronDown,
} from 'lucide-react';
import { useStore } from '@/lib/store';
import { Switch } from '@/components/ui/switch';
import { Slider } from '@/components/ui/slider';

interface FocusModeControlsProps {
  focusType: 'metal' | 'ligand';
  hasWaters: boolean;
  hasInteractions: boolean;
  hasPharmacophores: boolean;
  onRadiusChange: (radius: number) => void;
  onReset: () => void;
}

export function FocusModeControls({
  focusType,
  hasWaters,
  hasInteractions,
  hasPharmacophores,
  onRadiusChange,
  onReset,
}: FocusModeControlsProps) {
  const { focusSettings, setFocusSettings } = useStore();
  const [expanded, setExpanded] = useState(true);

  const currentRadius = focusType === 'metal'
    ? focusSettings.coordinationRadius
    : focusSettings.bindingPocketRadius;

  const radiusLabel = focusType === 'metal'
    ? 'Coordination Radius'
    : 'Pocket Radius';

  // Check if any controls should be shown
  const hasAnyControls = hasWaters || hasInteractions || hasPharmacophores;

  const handleRadiusChange = (values: number[]) => {
    const value = values[0];
    if (focusType === 'metal') {
      setFocusSettings({ coordinationRadius: value });
    } else {
      setFocusSettings({ bindingPocketRadius: value });
    }
    onRadiusChange(value);
  };

  return (
    <div className="absolute top-2 left-2 z-10">
      {/* Focus Mode Badge */}
      <div className="flex items-center gap-2 mb-2">
        <div className={`px-3 py-1.5 rounded-lg shadow-lg text-xs font-semibold text-white ${
          focusType === 'metal' ? 'bg-purple-600' : 'bg-emerald-600'
        }`}>
          {focusType === 'metal' ? 'Metal Focus' : 'Ligand Focus'}
        </div>
        <button
          onClick={onReset}
          className="p-1.5 rounded-lg bg-card/90 shadow border border-border hover:bg-muted transition-colors"
          title="Reset View"
        >
          <X className="w-4 h-4 text-muted-foreground" />
        </button>
      </div>

      {/* Controls Panel - only if there are controls to show */}
      {hasAnyControls && (
        <div className="bg-card/95 backdrop-blur-sm rounded-lg shadow-lg border border-border overflow-hidden min-w-[200px]">
          {/* Header */}
          <button
            onClick={() => setExpanded(!expanded)}
            className="w-full px-3 py-2 flex items-center justify-between text-xs font-medium text-foreground hover:bg-muted"
          >
            <span className="flex items-center gap-1.5">
              <SlidersHorizontal className="w-3.5 h-3.5" />
              View Options
            </span>
            <ChevronDown className={`w-4 h-4 transition-transform ${expanded ? '' : '-rotate-90'}`} />
          </button>

          {expanded && (
            <div className="px-3 pb-3 space-y-3 border-t border-border">
              {/* Radius Slider */}
              <div className="pt-3">
                <div className="flex items-center justify-between mb-2">
                  <label className="text-xs text-muted-foreground">{radiusLabel}</label>
                  <span className="text-xs font-mono text-muted-foreground">{currentRadius.toFixed(1)}Å</span>
                </div>
                <Slider
                  min={2.0}
                  max={8.0}
                  step={0.5}
                  value={[currentRadius]}
                  onValueChange={handleRadiusChange}
                />
              </div>

              {/* Water Toggle - only if waters exist */}
              {hasWaters && (
                <div className="flex items-center justify-between">
                  <span className="flex items-center gap-1.5 text-xs text-muted-foreground">
                    <Droplet className="w-3.5 h-3.5 text-cyan-500" />
                    Show Waters
                  </span>
                  <Switch
                    checked={focusSettings.showWaters}
                    onCheckedChange={(checked) => setFocusSettings({ showWaters: checked })}
                  />
                </div>
              )}

              {/* Interaction Lines Toggle - only if interactions exist */}
              {hasInteractions && (
                <div className="flex items-center justify-between">
                  <span className="flex items-center gap-1.5 text-xs text-muted-foreground">
                    <Link2 className="w-3.5 h-3.5 text-blue-500" />
                    Interaction Lines
                  </span>
                  <Switch
                    checked={focusSettings.showInteractionLines}
                    onCheckedChange={(checked) => setFocusSettings({ showInteractionLines: checked })}
                  />
                </div>
              )}

              {/* Pharmacophores Toggle - only if pharmacophores exist */}
              {hasPharmacophores && (
                <div className="flex items-center justify-between">
                  <span className="flex items-center gap-1.5 text-xs text-muted-foreground">
                    <CircleDot className="w-3.5 h-3.5 text-purple-500" />
                    Pharmacophores
                  </span>
                  <Switch
                    checked={focusSettings.showPharmacophores}
                    onCheckedChange={(checked) => setFocusSettings({ showPharmacophores: checked })}
                  />
                </div>
              )}
            </div>
          )}
        </div>
      )}
    </div>
  );
}

export default FocusModeControls;
