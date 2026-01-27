'use client';

import { useState } from 'react';
import {
  Atom,
  Pill,
  Palette,
  RotateCcw,
  Layers,
  ChevronDown,
  Activity,
  GitCompare,
  Search,
  Flame,
  Dna
} from 'lucide-react';
import type { ViewerMode, RepresentationStyle, ColorScheme } from '@/lib/store';

interface ViewerControlsProps {
  viewerMode: ViewerMode;
  representationStyle: RepresentationStyle;
  colorScheme: ColorScheme;
  onViewerModeChange: (mode: ViewerMode) => void;
  onRepresentationChange: (style: RepresentationStyle) => void;
  onColorSchemeChange: (scheme: ColorScheme) => void;
  onResetView: () => void;
  onRunAnalysis: () => void;
  onFocusFirstMetal: () => void;
  onFocusFirstLigand: () => void;
  hasStructure: boolean;
  hasConfidences: boolean;
  hasComparison: boolean;
  hasMetals: boolean;
  hasLigands: boolean;
  hasHotspots: boolean;
  showHotspots3D: boolean;
  onToggleHotspots: () => void;
  hasConservation: boolean;
  conservationLoading: boolean;
  showConservation3D: boolean;
  onToggleConservation: () => void;
  analysisLoading?: boolean;
}

// Representation options
const REPRESENTATIONS: { value: RepresentationStyle; label: string }[] = [
  { value: 'cartoon', label: 'Cartoon' },
  { value: 'ball-and-stick', label: 'Ball & Stick' },
  { value: 'spacefill', label: 'Spacefill' },
  { value: 'surface', label: 'Surface' }
];

// Color scheme options
const COLOR_SCHEMES: { value: ColorScheme; label: string; requiresConfidence?: boolean }[] = [
  { value: 'default', label: 'Default' },
  { value: 'chain', label: 'By Chain' },
  { value: 'residue-type', label: 'Residue Type' },
  { value: 'secondary-structure', label: 'Secondary Structure' },
  { value: 'confidence', label: 'Confidence (pLDDT)', requiresConfidence: true },
  { value: 'hydrophobicity', label: 'Hydrophobicity' }
];

export function ViewerControls({
  viewerMode,
  representationStyle,
  colorScheme,
  onViewerModeChange,
  onRepresentationChange,
  onColorSchemeChange,
  onResetView,
  onRunAnalysis,
  onFocusFirstMetal,
  onFocusFirstLigand,
  hasStructure,
  hasConfidences,
  hasComparison,
  hasMetals,
  hasLigands,
  hasHotspots,
  showHotspots3D,
  onToggleHotspots,
  hasConservation,
  conservationLoading,
  showConservation3D,
  onToggleConservation,
  analysisLoading = false
}: ViewerControlsProps) {
  const [showRepDropdown, setShowRepDropdown] = useState(false);
  const [showColorDropdown, setShowColorDropdown] = useState(false);

  return (
    <div className="flex items-center gap-2 px-4 py-2 bg-muted/80 border-b border-border">
      {/* Analysis mode buttons */}
      <div className="flex items-center gap-1 border-r border-border pr-3 mr-1">
        <button
          onClick={() => {
            if (viewerMode === 'metal') {
              // Already in metal mode, reset view
              onResetView();
              onViewerModeChange('default');
            } else if (hasMetals) {
              // Focus on first metal
              onFocusFirstMetal();
              onViewerModeChange('metal');
            } else {
              // No metals, just toggle mode (will run analysis if needed)
              onViewerModeChange('metal');
            }
          }}
          className={`px-2.5 py-1.5 rounded-lg text-xs font-medium flex items-center gap-1.5 transition-colors ${
            viewerMode === 'metal'
              ? 'bg-purple-100 text-purple-700 border border-purple-300'
              : 'hover:bg-muted text-muted-foreground'
          }`}
          disabled={!hasStructure}
          title={hasMetals ? "Focus on metal binding site" : "Show metal coordination"}
        >
          <Atom className="w-4 h-4" />
          <span>Metals</span>
          {hasMetals && viewerMode !== 'metal' && (
            <span className="w-1.5 h-1.5 rounded-full bg-purple-500" />
          )}
        </button>

        <button
          onClick={() => {
            if (viewerMode === 'ligand') {
              // Already in ligand mode, reset view
              onResetView();
              onViewerModeChange('default');
            } else if (hasLigands) {
              // Focus on first ligand
              onFocusFirstLigand();
              onViewerModeChange('ligand');
            } else {
              // No ligands, just toggle mode
              onViewerModeChange('ligand');
            }
          }}
          className={`px-2.5 py-1.5 rounded-lg text-xs font-medium flex items-center gap-1.5 transition-colors ${
            viewerMode === 'ligand'
              ? 'bg-emerald-100 text-emerald-700 border border-emerald-300'
              : 'hover:bg-muted text-muted-foreground'
          }`}
          disabled={!hasStructure}
          title={hasLigands ? "Focus on ligand binding site" : "Show ligand contacts"}
        >
          <Pill className="w-4 h-4" />
          <span>Ligands</span>
          {hasLigands && viewerMode !== 'ligand' && (
            <span className="w-1.5 h-1.5 rounded-full bg-emerald-500" />
          )}
        </button>

        {hasHotspots && (
          <button
            onClick={onToggleHotspots}
            className={`px-2.5 py-1.5 rounded-lg text-xs font-medium flex items-center gap-1.5 transition-colors ${
              showHotspots3D
                ? 'bg-rose-100 text-rose-700 border border-rose-300'
                : 'hover:bg-muted text-muted-foreground'
            }`}
            disabled={!hasStructure}
            title="Show detected binding hotspots"
          >
            <Flame className="w-4 h-4" />
            <span>Hotspots</span>
            {!showHotspots3D && (
              <span className="w-1.5 h-1.5 rounded-full bg-rose-500" />
            )}
          </button>
        )}

        {(hasConservation || conservationLoading) && (
          <button
            onClick={onToggleConservation}
            className={`px-2.5 py-1.5 rounded-lg text-xs font-medium flex items-center gap-1.5 transition-colors ${
              showConservation3D
                ? 'bg-violet-100 text-violet-700 border border-violet-300'
                : 'hover:bg-muted text-muted-foreground'
            }`}
            disabled={!hasStructure || conservationLoading}
            title="Color by evolutionary conservation (ConSurf)"
          >
            {conservationLoading ? (
              <div className="w-4 h-4 border-2 border-violet-300 border-t-transparent rounded-full animate-spin" />
            ) : (
              <Dna className="w-4 h-4" />
            )}
            <span>ConSurf</span>
            {hasConservation && !showConservation3D && !conservationLoading && (
              <span className="w-1.5 h-1.5 rounded-full bg-violet-500" />
            )}
          </button>
        )}

        {hasConfidences && (
          <button
            onClick={() => onViewerModeChange(viewerMode === 'confidence' ? 'default' : 'confidence')}
            className={`px-2.5 py-1.5 rounded-lg text-xs font-medium flex items-center gap-1.5 transition-colors ${
              viewerMode === 'confidence'
                ? 'bg-blue-100 text-blue-700 border border-blue-300'
                : 'hover:bg-muted text-muted-foreground'
            }`}
            title="Color by pLDDT confidence"
          >
            <Activity className="w-4 h-4" />
            <span>Confidence</span>
          </button>
        )}

        {hasComparison && (
          <button
            onClick={() => onViewerModeChange(viewerMode === 'comparison' ? 'default' : 'comparison')}
            className={`px-2.5 py-1.5 rounded-lg text-xs font-medium flex items-center gap-1.5 transition-colors ${
              viewerMode === 'comparison'
                ? 'bg-amber-100 text-amber-700 border border-amber-300'
                : 'hover:bg-muted text-muted-foreground'
            }`}
            title="Compare structures"
          >
            <GitCompare className="w-4 h-4" />
            <span>Compare</span>
          </button>
        )}
      </div>

      {/* Run analysis button */}
      <button
        onClick={onRunAnalysis}
        disabled={!hasStructure || analysisLoading}
        className={`px-2.5 py-1.5 rounded-lg text-xs font-medium flex items-center gap-1.5 transition-colors ${
          analysisLoading
            ? 'bg-muted text-muted-foreground cursor-wait'
            : 'bg-blue-50 hover:bg-blue-100 text-blue-700 border border-blue-200'
        }`}
        title="Run binding site analysis"
      >
        {analysisLoading ? (
          <div className="w-4 h-4 border-2 border-blue-300 border-t-transparent rounded-full animate-spin" />
        ) : (
          <Search className="w-4 h-4" />
        )}
        <span>Analyze</span>
      </button>

      {/* Spacer */}
      <div className="flex-1" />

      {/* Representation dropdown */}
      <div className="relative">
        <button
          onClick={() => {
            setShowRepDropdown(!showRepDropdown);
            setShowColorDropdown(false);
          }}
          className="px-2.5 py-1.5 rounded-lg text-xs font-medium flex items-center gap-1.5 hover:bg-muted text-muted-foreground transition-colors"
          disabled={!hasStructure}
        >
          <Layers className="w-4 h-4" />
          <span className="hidden md:inline">{REPRESENTATIONS.find(r => r.value === representationStyle)?.label}</span>
          <ChevronDown className="w-3 h-3" />
        </button>

        {showRepDropdown && (
          <div className="absolute right-0 mt-1 w-40 bg-card border border-border rounded-lg shadow-lg z-10">
            {REPRESENTATIONS.map((rep) => (
              <button
                key={rep.value}
                onClick={() => {
                  onRepresentationChange(rep.value);
                  setShowRepDropdown(false);
                }}
                className={`w-full px-3 py-2 text-left text-sm hover:bg-muted/50 first:rounded-t-lg last:rounded-b-lg ${
                  representationStyle === rep.value ? 'text-blue-600 bg-blue-50 font-medium' : 'text-foreground'
                }`}
              >
                {rep.label}
              </button>
            ))}
          </div>
        )}
      </div>

      {/* Color scheme dropdown */}
      <div className="relative">
        <button
          onClick={() => {
            setShowColorDropdown(!showColorDropdown);
            setShowRepDropdown(false);
          }}
          className="px-2.5 py-1.5 rounded-lg text-xs font-medium flex items-center gap-1.5 hover:bg-muted text-muted-foreground transition-colors"
          disabled={!hasStructure}
        >
          <Palette className="w-4 h-4" />
          <span className="hidden md:inline">{COLOR_SCHEMES.find(c => c.value === colorScheme)?.label}</span>
          <ChevronDown className="w-3 h-3" />
        </button>

        {showColorDropdown && (
          <div className="absolute right-0 mt-1 w-48 bg-card border border-border rounded-lg shadow-lg z-10">
            {COLOR_SCHEMES.map((scheme) => {
              const disabled = scheme.requiresConfidence && !hasConfidences;
              return (
                <button
                  key={scheme.value}
                  onClick={() => {
                    if (!disabled) {
                      onColorSchemeChange(scheme.value);
                      setShowColorDropdown(false);
                    }
                  }}
                  disabled={disabled}
                  className={`w-full px-3 py-2 text-left text-sm first:rounded-t-lg last:rounded-b-lg ${
                    disabled
                      ? 'text-muted-foreground/50 cursor-not-allowed'
                      : colorScheme === scheme.value
                        ? 'text-blue-600 bg-blue-50 font-medium'
                        : 'text-foreground hover:bg-muted/50'
                  }`}
                >
                  {scheme.label}
                  {disabled && <span className="text-xs text-muted-foreground ml-1">(no data)</span>}
                </button>
              );
            })}
          </div>
        )}
      </div>

      {/* Reset view button */}
      <button
        onClick={onResetView}
        className="p-1.5 rounded-lg hover:bg-muted text-muted-foreground transition-colors"
        disabled={!hasStructure}
        title="Reset view"
      >
        <RotateCcw className="w-4 h-4" />
      </button>

      {/* Click outside handler */}
      {(showRepDropdown || showColorDropdown) && (
        <div
          className="fixed inset-0 z-0"
          onClick={() => {
            setShowRepDropdown(false);
            setShowColorDropdown(false);
          }}
        />
      )}
    </div>
  );
}

export default ViewerControls;
