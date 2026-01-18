'use client';

import { useState, useCallback, useEffect } from 'react';
import { useStore } from '@/lib/store';
import { GitCompare, Layers, SplitSquareHorizontal, Eye, EyeOff, RefreshCw, Info } from 'lucide-react';

interface ComparisonViewProps {
  onLoadReference: (pdb: string, label: string, source: 'rfd3' | 'rf3') => void;
  onClearReference: () => void;
  onToggleOverlay: (visible: boolean) => void;
}

export function ComparisonView({
  onLoadReference,
  onClearReference,
  onToggleOverlay,
}: ComparisonViewProps) {
  const {
    latestRfd3Design,
    latestRf3Prediction,
    latestRmsdResult,
    referenceStructure,
    comparisonEnabled,
    setComparisonEnabled,
    comparisonMode,
    setComparisonMode,
    setReferenceStructure,
  } = useStore();

  const [overlayVisible, setOverlayVisible] = useState(true);
  const [selectedReference, setSelectedReference] = useState<'rfd3' | 'rf3' | null>(null);

  // Auto-select RFD3 as reference when both structures are available
  useEffect(() => {
    if (latestRfd3Design && latestRf3Prediction && !referenceStructure) {
      // Default: RFD3 design is the reference (what we designed)
      // RF3 prediction is shown as the main structure (what it folds to)
      setSelectedReference('rfd3');
    }
  }, [latestRfd3Design, latestRf3Prediction, referenceStructure]);

  const handleSetReference = useCallback((source: 'rfd3' | 'rf3') => {
    setSelectedReference(source);

    if (source === 'rfd3' && latestRfd3Design) {
      setReferenceStructure({
        pdb: latestRfd3Design.pdbContent,
        label: 'RFD3 Design',
        source: 'rfd3',
      });
      onLoadReference(latestRfd3Design.pdbContent, 'RFD3 Design', 'rfd3');
    } else if (source === 'rf3' && latestRf3Prediction) {
      setReferenceStructure({
        pdb: latestRf3Prediction.pdbContent,
        label: 'RF3 Prediction',
        source: 'rf3',
      });
      onLoadReference(latestRf3Prediction.pdbContent, 'RF3 Prediction', 'rf3');
    }

    setComparisonEnabled(true);
  }, [latestRfd3Design, latestRf3Prediction, setReferenceStructure, setComparisonEnabled, onLoadReference]);

  const handleClearComparison = useCallback(() => {
    setSelectedReference(null);
    setReferenceStructure(null);
    setComparisonEnabled(false);
    onClearReference();
  }, [setReferenceStructure, setComparisonEnabled, onClearReference]);

  const handleToggleOverlay = useCallback(() => {
    const newVisible = !overlayVisible;
    setOverlayVisible(newVisible);
    onToggleOverlay(newVisible);
  }, [overlayVisible, onToggleOverlay]);

  // Get RMSD color based on value
  const getRmsdColor = (rmsd: number | null | undefined): string => {
    if (rmsd === null || rmsd === undefined) return 'text-muted-foreground';
    if (rmsd < 1.0) return 'text-green-600';
    if (rmsd < 2.0) return 'text-blue-600';
    if (rmsd < 3.0) return 'text-amber-600';
    return 'text-red-600';
  };

  const getRmsdBgColor = (rmsd: number | null | undefined): string => {
    if (rmsd === null || rmsd === undefined) return 'bg-muted';
    if (rmsd < 1.0) return 'bg-green-50 border-green-200';
    if (rmsd < 2.0) return 'bg-blue-50 border-blue-200';
    if (rmsd < 3.0) return 'bg-amber-50 border-amber-200';
    return 'bg-red-50 border-red-200';
  };

  // Check if we have structures to compare
  const hasRfd3 = !!latestRfd3Design?.pdbContent;
  const hasRf3 = !!latestRf3Prediction?.pdbContent;
  const canCompare = hasRfd3 && hasRf3;

  if (!canCompare) {
    return (
      <div className="p-4 bg-muted/50 rounded-xl border border-border">
        <div className="flex items-start gap-3">
          <div className="p-2 rounded-lg bg-muted">
            <GitCompare className="w-5 h-5 text-muted-foreground" />
          </div>
          <div>
            <h3 className="font-medium text-foreground text-sm">Structure Comparison</h3>
            <p className="text-xs text-muted-foreground mt-1">
              {!hasRfd3 && !hasRf3
                ? 'Run RFD3 design and RF3 prediction to enable comparison'
                : !hasRfd3
                  ? 'Run RFD3 design to compare with RF3 prediction'
                  : 'Run RF3 prediction to compare with RFD3 design'
              }
            </p>
          </div>
        </div>
      </div>
    );
  }

  return (
    <div className="space-y-4">
      {/* Header */}
      <div className="flex items-center justify-between">
        <div className="flex items-center gap-2">
          <GitCompare className="w-5 h-5 text-violet-600" />
          <h3 className="font-semibold text-foreground text-sm">Structure Comparison</h3>
        </div>

        {comparisonEnabled && (
          <button
            onClick={handleClearComparison}
            className="text-xs text-muted-foreground hover:text-foreground flex items-center gap-1"
          >
            <RefreshCw className="w-3 h-3" />
            Reset
          </button>
        )}
      </div>

      {/* RMSD Result Banner */}
      {latestRmsdResult && (
        <div className={`p-3 rounded-xl border ${getRmsdBgColor(latestRmsdResult.rmsd)}`}>
          <div className="flex items-center justify-between">
            <div className="flex items-center gap-2">
              <span className={`text-lg font-bold ${getRmsdColor(latestRmsdResult.rmsd)}`}>
                {latestRmsdResult.rmsd?.toFixed(2)} Å
              </span>
              <span className="text-xs font-medium text-muted-foreground">
                {latestRmsdResult.interpretation}
              </span>
            </div>
            <div className="group relative">
              <Info className="w-4 h-4 text-muted-foreground cursor-help" />
              <div className="absolute right-0 bottom-6 w-64 p-3 bg-card rounded-lg shadow-lg border border-border hidden group-hover:block z-10">
                <p className="text-xs text-muted-foreground">{latestRmsdResult.description}</p>
                <div className="mt-2 space-y-1 text-xs">
                  <div className="flex items-center gap-2">
                    <span className="w-2 h-2 rounded-full bg-green-500" />
                    <span className="text-muted-foreground">&lt; 1.0 Å: Excellent</span>
                  </div>
                  <div className="flex items-center gap-2">
                    <span className="w-2 h-2 rounded-full bg-blue-500" />
                    <span className="text-muted-foreground">1.0-2.0 Å: Good</span>
                  </div>
                  <div className="flex items-center gap-2">
                    <span className="w-2 h-2 rounded-full bg-amber-500" />
                    <span className="text-muted-foreground">2.0-3.0 Å: Moderate</span>
                  </div>
                  <div className="flex items-center gap-2">
                    <span className="w-2 h-2 rounded-full bg-red-500" />
                    <span className="text-muted-foreground">&gt; 3.0 Å: Poor</span>
                  </div>
                </div>
              </div>
            </div>
          </div>
          <p className="text-xs text-muted-foreground mt-1">
            Backbone RMSD between RFD3 design and RF3 prediction
          </p>
        </div>
      )}

      {/* Mode selector */}
      <div className="flex gap-2 p-1 bg-muted rounded-lg">
        <button
          onClick={() => setComparisonMode('overlay')}
          className={`flex-1 flex items-center justify-center gap-1.5 py-2 rounded-md text-xs font-medium transition-colors ${
            comparisonMode === 'overlay'
              ? 'bg-card text-violet-700 shadow-sm'
              : 'text-muted-foreground hover:text-foreground'
          }`}
        >
          <Layers className="w-3.5 h-3.5" />
          Overlay
        </button>
        <button
          onClick={() => setComparisonMode('side-by-side')}
          className={`flex-1 flex items-center justify-center gap-1.5 py-2 rounded-md text-xs font-medium transition-colors ${
            comparisonMode === 'side-by-side'
              ? 'bg-card text-violet-700 shadow-sm'
              : 'text-muted-foreground hover:text-foreground'
          }`}
        >
          <SplitSquareHorizontal className="w-3.5 h-3.5" />
          Side by Side
        </button>
      </div>

      {/* Structure selection */}
      <div className="space-y-2">
        <p className="text-xs font-medium text-muted-foreground">Select reference structure:</p>

        <div className="grid grid-cols-2 gap-2">
          {/* RFD3 Design */}
          <button
            onClick={() => handleSetReference('rfd3')}
            disabled={!hasRfd3}
            className={`p-3 rounded-xl border-2 text-left transition-all ${
              selectedReference === 'rfd3'
                ? 'border-violet-400 bg-violet-50'
                : 'border-border hover:border-border'
            } ${!hasRfd3 ? 'opacity-50 cursor-not-allowed' : ''}`}
          >
            <div className="flex items-center gap-2">
              <div className={`w-3 h-3 rounded-full ${selectedReference === 'rfd3' ? 'bg-violet-500' : 'bg-muted'}`} />
              <span className="text-sm font-medium text-foreground">RFD3 Design</span>
            </div>
            <p className="text-xs text-muted-foreground mt-1">Original backbone design</p>
          </button>

          {/* RF3 Prediction */}
          <button
            onClick={() => handleSetReference('rf3')}
            disabled={!hasRf3}
            className={`p-3 rounded-xl border-2 text-left transition-all ${
              selectedReference === 'rf3'
                ? 'border-emerald-400 bg-emerald-50'
                : 'border-border hover:border-border'
            } ${!hasRf3 ? 'opacity-50 cursor-not-allowed' : ''}`}
          >
            <div className="flex items-center gap-2">
              <div className={`w-3 h-3 rounded-full ${selectedReference === 'rf3' ? 'bg-emerald-500' : 'bg-muted'}`} />
              <span className="text-sm font-medium text-foreground">RF3 Prediction</span>
            </div>
            <p className="text-xs text-muted-foreground mt-1">Structure prediction</p>
          </button>
        </div>
      </div>

      {/* Overlay controls (when in overlay mode and comparison is active) */}
      {comparisonMode === 'overlay' && comparisonEnabled && (
        <div className="flex items-center justify-between p-3 bg-muted/50 rounded-xl">
          <div className="flex items-center gap-2">
            <span className="text-xs font-medium text-muted-foreground">Reference overlay:</span>
            <span className="text-xs text-muted-foreground">
              {referenceStructure?.label || 'None'}
            </span>
          </div>
          <button
            onClick={handleToggleOverlay}
            className={`p-2 rounded-lg transition-colors ${
              overlayVisible
                ? 'bg-violet-100 text-violet-700'
                : 'bg-muted text-muted-foreground'
            }`}
            title={overlayVisible ? 'Hide reference' : 'Show reference'}
          >
            {overlayVisible ? <Eye className="w-4 h-4" /> : <EyeOff className="w-4 h-4" />}
          </button>
        </div>
      )}

      {/* Color legend */}
      {comparisonEnabled && (
        <div className="p-3 bg-muted/50 rounded-xl">
          <p className="text-xs font-medium text-muted-foreground mb-2">Structure colors:</p>
          <div className="flex flex-wrap gap-3">
            <div className="flex items-center gap-1.5">
              <div className="w-3 h-3 rounded-full bg-violet-500" />
              <span className="text-xs text-muted-foreground">RFD3 Design</span>
            </div>
            <div className="flex items-center gap-1.5">
              <div className="w-3 h-3 rounded-full bg-emerald-500" />
              <span className="text-xs text-muted-foreground">RF3 Prediction</span>
            </div>
          </div>
        </div>
      )}

      {/* Info about comparison */}
      {!comparisonEnabled && (
        <div className="p-3 bg-blue-50 rounded-xl border border-blue-200">
          <div className="flex items-start gap-2">
            <Info className="w-4 h-4 text-blue-500 mt-0.5 flex-shrink-0" />
            <div className="text-xs text-blue-700">
              <p className="font-medium">How comparison works:</p>
              <ul className="mt-1 space-y-0.5 text-blue-600">
                <li>• Select a reference structure above</li>
                <li>• The reference will be shown as a translucent overlay</li>
                <li>• Compare backbone alignment to assess designability</li>
                <li>• Lower RMSD = better fold prediction</li>
              </ul>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}

export default ComparisonView;
