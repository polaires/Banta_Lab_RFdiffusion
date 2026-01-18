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
          className="p-1.5 rounded-lg bg-white/90 shadow hover:bg-white transition-colors"
          title="Reset View"
        >
          <X className="w-4 h-4 text-gray-600" />
        </button>
      </div>

      {/* Controls Panel - only if there are controls to show */}
      {hasAnyControls && (
        <div className="bg-white/95 rounded-lg shadow-lg border border-gray-200 overflow-hidden min-w-[200px]">
          {/* Header */}
          <button
            onClick={() => setExpanded(!expanded)}
            className="w-full px-3 py-2 flex items-center justify-between text-xs font-medium text-gray-700 hover:bg-gray-50"
          >
            <span className="flex items-center gap-1.5">
              <SlidersHorizontal className="w-3.5 h-3.5" />
              View Options
            </span>
            <ChevronDown className={`w-4 h-4 transition-transform ${expanded ? '' : '-rotate-90'}`} />
          </button>

          {expanded && (
            <div className="px-3 pb-3 space-y-3 border-t border-gray-100">
              {/* Radius Slider */}
              <div className="pt-3">
                <div className="flex items-center justify-between mb-1.5">
                  <label className="text-xs text-gray-600">{radiusLabel}</label>
                  <span className="text-xs font-mono text-gray-500">{currentRadius.toFixed(1)}Å</span>
                </div>
                <input
                  type="range"
                  min={2.0}
                  max={8.0}
                  step={0.5}
                  value={currentRadius}
                  onChange={(e) => {
                    const value = parseFloat(e.target.value);
                    if (focusType === 'metal') {
                      setFocusSettings({ coordinationRadius: value });
                    } else {
                      setFocusSettings({ bindingPocketRadius: value });
                    }
                    onRadiusChange(value);
                  }}
                  className="w-full h-1.5 bg-gray-200 rounded-lg appearance-none cursor-pointer accent-blue-600"
                />
              </div>

              {/* Water Toggle - only if waters exist */}
              {hasWaters && (
                <label className="flex items-center justify-between cursor-pointer group">
                  <span className="flex items-center gap-1.5 text-xs text-gray-600">
                    <Droplet className="w-3.5 h-3.5 text-cyan-500" />
                    Show Waters
                  </span>
                  <div className="relative">
                    <input
                      type="checkbox"
                      checked={focusSettings.showWaters}
                      onChange={(e) => setFocusSettings({ showWaters: e.target.checked })}
                      className="sr-only peer"
                    />
                    <div className="w-8 h-4 bg-gray-200 rounded-full peer peer-checked:bg-blue-600 transition-colors" />
                    <div className="absolute top-0.5 left-0.5 w-3 h-3 bg-white rounded-full shadow peer-checked:translate-x-4 transition-transform" />
                  </div>
                </label>
              )}

              {/* Interaction Lines Toggle - only if interactions exist */}
              {hasInteractions && (
                <label className="flex items-center justify-between cursor-pointer group">
                  <span className="flex items-center gap-1.5 text-xs text-gray-600">
                    <Link2 className="w-3.5 h-3.5 text-blue-500" />
                    Interaction Lines
                  </span>
                  <div className="relative">
                    <input
                      type="checkbox"
                      checked={focusSettings.showInteractionLines}
                      onChange={(e) => setFocusSettings({ showInteractionLines: e.target.checked })}
                      className="sr-only peer"
                    />
                    <div className="w-8 h-4 bg-gray-200 rounded-full peer peer-checked:bg-blue-600 transition-colors" />
                    <div className="absolute top-0.5 left-0.5 w-3 h-3 bg-white rounded-full shadow peer-checked:translate-x-4 transition-transform" />
                  </div>
                </label>
              )}

              {/* Pharmacophores Toggle - only if pharmacophores exist */}
              {hasPharmacophores && (
                <label className="flex items-center justify-between cursor-pointer group">
                  <span className="flex items-center gap-1.5 text-xs text-gray-600">
                    <CircleDot className="w-3.5 h-3.5 text-purple-500" />
                    Pharmacophores
                  </span>
                  <div className="relative">
                    <input
                      type="checkbox"
                      checked={focusSettings.showPharmacophores}
                      onChange={(e) => setFocusSettings({ showPharmacophores: e.target.checked })}
                      className="sr-only peer"
                    />
                    <div className="w-8 h-4 bg-gray-200 rounded-full peer peer-checked:bg-blue-600 transition-colors" />
                    <div className="absolute top-0.5 left-0.5 w-3 h-3 bg-white rounded-full shadow peer-checked:translate-x-4 transition-transform" />
                  </div>
                </label>
              )}
            </div>
          )}
        </div>
      )}
    </div>
  );
}

export default FocusModeControls;
