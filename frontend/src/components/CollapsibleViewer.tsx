'use client';

import { useCallback, useEffect, useState } from 'react';
import dynamic from 'next/dynamic';
import { useStore } from '@/lib/store';
import { findMetalCoordinationFromPDB } from '@/lib/metalAnalysis';
import { findLigandContactsFromPDB, extractPharmacophoreFeatures } from '@/lib/ligandAnalysis';
import api from '@/lib/api';
import { Eye, ChevronDown, ChevronUp, Box } from 'lucide-react';

// Dynamic imports for components that use Molstar (SSR incompatible)
const ComparisonView = dynamic(
  () => import('@/components/viewer/ComparisonView').then((mod) => mod.ComparisonView),
  { ssr: false }
);

const ProteinViewer = dynamic(
  () => import('@/components/ProteinViewer').then((mod) => mod.ProteinViewer),
  {
    ssr: false,
    loading: () => (
      <div className="h-[520px] flex items-center justify-center bg-muted rounded-xl">
        <div className="flex items-center gap-3 text-muted-foreground">
          <div className="w-5 h-5 border-2 border-border border-t-transparent rounded-full animate-spin" />
          <span className="text-sm font-medium">Loading 3D viewer...</span>
        </div>
      </div>
    ),
  }
);

const ViewerControls = dynamic(
  () => import('@/components/viewer/ViewerControls').then((mod) => mod.ViewerControls),
  { ssr: false }
);

const BindingSitePanel = dynamic(
  () => import('@/components/viewer/BindingSitePanel').then((mod) => mod.BindingSitePanel),
  { ssr: false }
);

const LigandBindingPanel = dynamic(
  () => import('@/components/LigandBindingPanel').then((mod) => mod.LigandBindingPanel),
  { ssr: false }
);

export function CollapsibleViewer() {
  const {
    selectedPdb,
    viewerCollapsed,
    setViewerCollapsed,
    viewerMode,
    setViewerMode,
    representationStyle,
    setRepresentationStyle,
    colorScheme,
    setColorScheme,
    metalCoordination,
    setMetalCoordination,
    ligandData,
    setLigandData,
    focusedMetalIndex,
    setFocusedMetalIndex,
    focusedLigandIndex,
    setFocusedLigandIndex,
    analysisLoading,
    setAnalysisLoading,
    latestConfidences,
    latestRfd3Design,
    latestRf3Prediction,
    comparisonEnabled,
    comparisonMode,
    // Pharmacophore visualization state
    showPharmacophores3D,
    setShowPharmacophores3D,
    pharmacophoreFeatures,
    setPharmacophoreFeatures,
    // Hotspot visualization state
    hotspotsData,
    setHotspotsData,
    showHotspots3D,
    setShowHotspots3D,
    // Conservation visualization state
    conservationData,
    setConservationData,
    conservationLoading,
    setConservationLoading,
    showConservation3D,
    setShowConservation3D,
    // Catalytic suggestions (for conservation merge)
    catalyticSuggestions,
    setCatalyticSuggestions,
    suggestionsSource,
  } = useStore();

  // Track if analysis has been run for current structure
  const [analyzedPdb, setAnalyzedPdb] = useState<string | null>(null);

  // Show binding site panel when there's analysis data
  const showBindingSitePanel = (metalCoordination && metalCoordination.length > 0) ||
    (ligandData && ligandData.ligandCount > 0);

  // Run structure analysis
  const runAnalysis = useCallback(async () => {
    if (!selectedPdb || analysisLoading) return;

    setAnalysisLoading(true);

    try {
      // Run metal coordination analysis
      const metals = findMetalCoordinationFromPDB(selectedPdb, 3.0);
      setMetalCoordination(metals);

      // Run ligand contact analysis
      const ligands = findLigandContactsFromPDB(selectedPdb, 4.0);
      setLigandData(ligands);

      setAnalyzedPdb(selectedPdb);

      // Auto-switch to metal mode if metals found
      if (metals.length > 0) {
        setViewerMode('metal');
      } else if (ligands.ligandCount > 0) {
        setViewerMode('ligand');
      }

      // Run hotspot detection in background (non-blocking)
      api.detectHotspots({
        target_pdb: selectedPdb,
        target_chain: 'A',
        method: 'exposed_clustered',
        max_hotspots: 5,
        prefer_hydrophobic: true,
      }).then((result) => {
        if (result.status === 'completed' && result.residue_details?.length > 0) {
          setHotspotsData(result);
          console.log(`[Analysis] Detected ${result.hotspots.length} hotspots`);
        }
      }).catch((err) => {
        console.warn('[Analysis] Hotspot detection failed (non-critical):', err);
      });

      // Run conservation analysis in background (non-blocking)
      setConservationLoading(true);
      api.analyzeConservation({ pdb_content: selectedPdb, chain: 'A' })
        .then((result) => {
          if (result.status === 'success' && result.grades?.length > 0) {
            setConservationData(result);
            console.log(`[Analysis] Conservation: ${result.grades.length} residues, MSA depth=${result.msa_depth}`);

            // Merge highly conserved residues into catalytic suggestions
            const conservedSuggestions = result.grades
              .filter((g) => g.grade <= 3)
              .map((g) => ({
                chain: 'A',
                residue: g.position,
                name: g.residue,
                role: `Conserved (grade ${g.grade})`,
                confidence: (4 - g.grade) / 3,
                source: 'conservation' as const,
              }));

            const existing = useStore.getState().catalyticSuggestions;
            const existingKeys = new Set(existing.map((s) => `${s.chain}${s.residue}`));
            const newSuggestions = conservedSuggestions.filter(
              (s) => !existingKeys.has(`${s.chain}${s.residue}`)
            );
            if (newSuggestions.length > 0) {
              setCatalyticSuggestions(
                [...existing, ...newSuggestions],
                existing.length > 0 ? suggestionsSource : 'conservation'
              );
              console.log(`[Analysis] Added ${newSuggestions.length} conserved residue suggestions`);
            }
          }
        })
        .catch((err) => {
          console.warn('[Analysis] Conservation analysis failed (non-critical):', err);
        })
        .finally(() => {
          setConservationLoading(false);
        });
    } catch (error) {
      console.error('Analysis failed:', error);
    } finally {
      setAnalysisLoading(false);
    }
  }, [selectedPdb, analysisLoading, setAnalysisLoading, setMetalCoordination, setLigandData, setViewerMode, setHotspotsData, setConservationData, setConservationLoading, setCatalyticSuggestions, suggestionsSource]);

  // Clear analysis when structure changes
  useEffect(() => {
    if (selectedPdb !== analyzedPdb) {
      // Structure changed, clear old analysis
      setMetalCoordination(null);
      setLigandData(null);
      setFocusedMetalIndex(null);
      setFocusedLigandIndex(null);
      setViewerMode('default');
      setPharmacophoreFeatures(null);
      setShowPharmacophores3D(false);
      setHotspotsData(null);
      setShowHotspots3D(false);
      setConservationData(null);
      setShowConservation3D(false);
    }
  }, [selectedPdb, analyzedPdb, setMetalCoordination, setLigandData, setFocusedMetalIndex, setFocusedLigandIndex, setViewerMode, setPharmacophoreFeatures, setShowPharmacophores3D, setHotspotsData, setShowHotspots3D, setConservationData, setShowConservation3D]);

  // Compute pharmacophore features when a ligand is focused
  useEffect(() => {
    if (ligandData && focusedLigandIndex !== null) {
      const ligand = ligandData.ligandDetails[focusedLigandIndex];
      if (ligand) {
        const features = extractPharmacophoreFeatures(ligand.contacts);
        setPharmacophoreFeatures(features);
      }
    } else {
      setPharmacophoreFeatures(null);
    }
  }, [ligandData, focusedLigandIndex, setPharmacophoreFeatures]);

  // Reset view handler
  const handleResetView = useCallback(() => {
    setFocusedMetalIndex(null);
    setFocusedLigandIndex(null);
    // The viewer will reset camera on focus change
  }, [setFocusedMetalIndex, setFocusedLigandIndex]);

  // Focus on first metal handler
  const handleFocusFirstMetal = useCallback(() => {
    if (metalCoordination && metalCoordination.length > 0) {
      setFocusedLigandIndex(null); // Clear ligand focus
      setFocusedMetalIndex(0);
    }
  }, [metalCoordination, setFocusedMetalIndex, setFocusedLigandIndex]);

  // Focus on first ligand handler
  const handleFocusFirstLigand = useCallback(() => {
    if (ligandData && ligandData.ligandCount > 0) {
      setFocusedMetalIndex(null); // Clear metal focus
      setFocusedLigandIndex(0);
    }
  }, [ligandData, setFocusedMetalIndex, setFocusedLigandIndex]);

  // Check if we have comparison structures
  const hasComparison = !!latestRfd3Design && !!latestRf3Prediction;

  // Check if we have metals/ligands
  const hasMetals = !!(metalCoordination && metalCoordination.length > 0);
  const hasLigands = !!(ligandData && ligandData.ligandCount > 0);

  // State for reference structure overlay visibility
  const [overlayVisible, setOverlayVisible] = useState(true);

  // Handler for loading reference structure into viewer
  const handleLoadReference = useCallback((pdb: string, label: string, source: 'rfd3' | 'rf3') => {
    // The ProteinViewer will handle the reference structure via the store
    console.log(`[ComparisonView] Loading reference: ${label} (${source})`);
    // Reference is already set in store by ComparisonView
    setOverlayVisible(true);
  }, []);

  // Handler for clearing reference structure
  const handleClearReference = useCallback(() => {
    console.log('[ComparisonView] Clearing reference');
    setOverlayVisible(false);
  }, []);

  // Handler for toggling overlay visibility
  const handleToggleOverlay = useCallback((visible: boolean) => {
    setOverlayVisible(visible);
  }, []);

  return (
    <section className="bg-card rounded-2xl shadow-card overflow-hidden mt-8">
      {/* Header */}
      <button
        onClick={() => setViewerCollapsed(!viewerCollapsed)}
        className="w-full px-6 py-4 border-b border-border flex justify-between items-center bg-muted/50 hover:bg-muted/50 transition-colors"
      >
        <div className="flex items-center gap-3">
          <Eye className="w-5 h-5 text-blue-600" />
          <h3 className="text-xs font-bold text-foreground uppercase tracking-widest">Structure Viewer</h3>
          {selectedPdb ? (
            <span className="text-[10px] font-semibold text-emerald-700 bg-emerald-50 px-2.5 py-1 rounded-full border border-emerald-200">
              Structure loaded
            </span>
          ) : (
            <span className="text-[10px] font-medium text-muted-foreground bg-card px-2 py-0.5 rounded-full border border-border shadow-sm">
              No structure
            </span>
          )}
          {metalCoordination && metalCoordination.length > 0 && (
            <span className="text-[10px] font-semibold text-purple-700 bg-purple-50 px-2 py-0.5 rounded-full border border-purple-200">
              {metalCoordination.length} metal{metalCoordination.length !== 1 ? 's' : ''}
            </span>
          )}
          {ligandData && ligandData.ligandCount > 0 && (
            <span className="text-[10px] font-semibold text-green-700 bg-green-50 px-2 py-0.5 rounded-full border border-green-200">
              {ligandData.ligandCount} ligand{ligandData.ligandCount !== 1 ? 's' : ''}
            </span>
          )}
          {comparisonEnabled && (
            <span className="text-[10px] font-semibold text-violet-700 bg-violet-50 px-2 py-0.5 rounded-full border border-violet-200">
              Comparing
            </span>
          )}
        </div>
        {viewerCollapsed ? (
          <ChevronDown className="w-5 h-5 text-muted-foreground hover:text-muted-foreground transition-colors" />
        ) : (
          <ChevronUp className="w-5 h-5 text-muted-foreground hover:text-muted-foreground transition-colors" />
        )}
      </button>

      {/* Viewer content */}
      {!viewerCollapsed && (
        <div className="flex">
          {/* Main viewer area */}
          <div className="flex-1">
            {/* Controls toolbar */}
            {selectedPdb && (
              <ViewerControls
                viewerMode={viewerMode}
                representationStyle={representationStyle}
                colorScheme={colorScheme}
                onViewerModeChange={setViewerMode}
                onRepresentationChange={setRepresentationStyle}
                onColorSchemeChange={setColorScheme}
                onResetView={handleResetView}
                onRunAnalysis={runAnalysis}
                onFocusFirstMetal={handleFocusFirstMetal}
                onFocusFirstLigand={handleFocusFirstLigand}
                hasStructure={!!selectedPdb}
                hasConfidences={!!latestConfidences?.per_residue_plddt}
                hasComparison={hasComparison}
                hasMetals={hasMetals}
                hasLigands={hasLigands}
                hasHotspots={!!(hotspotsData?.residue_details && hotspotsData.residue_details.length > 0)}
                showHotspots3D={showHotspots3D}
                onToggleHotspots={() => setShowHotspots3D(!showHotspots3D)}
                hasConservation={!!(conservationData?.grades && conservationData.grades.length > 0)}
                conservationLoading={conservationLoading}
                showConservation3D={showConservation3D}
                onToggleConservation={() => setShowConservation3D(!showConservation3D)}
                analysisLoading={analysisLoading}
              />
            )}

            {/* Viewer */}
            <div className="p-4">
              {selectedPdb ? (
                comparisonMode === 'side-by-side' && comparisonEnabled && latestRfd3Design && latestRf3Prediction ? (
                  // Side-by-side view: Two viewers
                  <div className="grid grid-cols-2 gap-4">
                    <div className="space-y-2">
                      <div className="flex items-center gap-2">
                        <div className="w-3 h-3 rounded-full bg-violet-500" />
                        <span className="text-sm font-medium text-foreground">RFD3 Design</span>
                      </div>
                      <ProteinViewer
                        pdbContent={latestRfd3Design.pdbContent}
                        className="h-[480px] rounded-xl overflow-hidden"
                        focusedMetalIndex={focusedMetalIndex}
                        focusedLigandIndex={focusedLigandIndex}
                        metalCoordination={metalCoordination}
                        ligandData={ligandData}
                        pharmacophoreFeatures={pharmacophoreFeatures ?? undefined}
                        showPharmacophores={showPharmacophores3D}
                      />
                    </div>
                    <div className="space-y-2">
                      <div className="flex items-center gap-2">
                        <div className="w-3 h-3 rounded-full bg-emerald-500" />
                        <span className="text-sm font-medium text-foreground">RF3 Prediction</span>
                      </div>
                      <ProteinViewer
                        pdbContent={latestRf3Prediction.pdbContent}
                        className="h-[480px] rounded-xl overflow-hidden"
                      />
                    </div>
                  </div>
                ) : (
                  // Single viewer (overlay mode or no comparison)
                  <ProteinViewer
                    pdbContent={selectedPdb}
                    className="h-[520px] rounded-xl overflow-hidden"
                    focusedMetalIndex={focusedMetalIndex}
                    focusedLigandIndex={focusedLigandIndex}
                    metalCoordination={metalCoordination}
                    ligandData={ligandData}
                    pharmacophoreFeatures={pharmacophoreFeatures ?? undefined}
                    showPharmacophores={showPharmacophores3D}
                  />
                )
              ) : (
                <div className="h-[520px] flex flex-col items-center justify-center bg-muted rounded-xl border-2 border-dashed border-border">
                  <Box className="w-10 h-10 text-muted-foreground mb-3" />
                  <p className="text-sm font-medium text-muted-foreground">No structure to display</p>
                  <p className="text-xs text-muted-foreground mt-1">Run a design job to visualize structures</p>
                </div>
              )}
            </div>
          </div>

          {/* Side panels container */}
          {(showBindingSitePanel || hasComparison) && (
            <div className="w-80 border-l border-border bg-muted/50 overflow-y-auto max-h-[640px]">
              {/* Comparison view panel */}
              {hasComparison && (
                <div className="p-4 border-b border-border">
                  <ComparisonView
                    onLoadReference={handleLoadReference}
                    onClearReference={handleClearReference}
                    onToggleOverlay={handleToggleOverlay}
                  />
                </div>
              )}

              {/* Binding site analysis panel */}
              {showBindingSitePanel && (
                <BindingSitePanel
                  metalCoordination={metalCoordination}
                  ligandData={ligandData}
                  focusedMetalIndex={focusedMetalIndex}
                  focusedLigandIndex={focusedLigandIndex}
                  onFocusMetal={setFocusedMetalIndex}
                  onFocusLigand={setFocusedLigandIndex}
                  loading={analysisLoading}
                />
              )}

              {/* Ligand binding analysis panel with pharmacophore controls */}
              {hasLigands && (
                <div className="p-4 border-t border-border">
                  <LigandBindingPanel
                    show3D={showPharmacophores3D}
                    onToggle3D={setShowPharmacophores3D}
                  />
                </div>
              )}
            </div>
          )}
        </div>
      )}
    </section>
  );
}
