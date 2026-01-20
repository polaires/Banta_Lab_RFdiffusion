'use client';

import { useState, useCallback, useRef, useEffect, useMemo } from 'react';
import {
  RotateCcw,
  Maximize2,
  Play,
  Pause,
  Microscope,
  GripHorizontal,
  Atom,
  Pill,
  ChevronDown,
  ChevronRight,
  Focus,
  Loader2,
} from 'lucide-react';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { Separator } from '@/components/ui/separator';
import {
  Collapsible,
  CollapsibleContent,
  CollapsibleTrigger,
} from '@/components/ui/collapsible';
import { cn } from '@/lib/utils';
import { useStore } from '@/lib/store';
import { FocusModeControls } from '@/components/viewer/FocusModeControls';
import { CatalyticSuggestionsPanel } from '@/components/CatalyticSuggestionsPanel';
import { findMetalCoordinationFromPDB, type MetalCoordination } from '@/lib/metalAnalysis';
import { findLigandContactsFromPDB, getInteractionLabel, type LigandData } from '@/lib/ligandAnalysis';

interface ViewerPanelProps {
  children: React.ReactNode;
  isSpinning?: boolean;
  onToggleSpin?: () => void;
  onReset?: () => void;
  onExpand?: () => void;
  showControls?: boolean;
}

// Metal card component - minimal theme
function MetalCard({
  metal,
  index,
  isFocused,
  onFocus,
}: {
  metal: MetalCoordination;
  index: number;
  isFocused: boolean;
  onFocus: (index: number | null) => void;
}) {
  const [expanded, setExpanded] = useState(false);

  return (
    <div
      className={cn(
        'rounded-md border transition-all',
        isFocused
          ? 'border-foreground/30 bg-accent'
          : 'border-border bg-card hover:bg-accent/50'
      )}
    >
      <div className="p-2">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-1.5 min-w-0">
            <button
              onClick={() => setExpanded(!expanded)}
              className="p-0.5 hover:bg-muted rounded text-muted-foreground shrink-0"
            >
              {expanded ? <ChevronDown className="w-3 h-3" /> : <ChevronRight className="w-3 h-3" />}
            </button>
            <Atom className="w-3.5 h-3.5 text-muted-foreground shrink-0" />
            <span className="font-medium text-xs text-foreground">{metal.element}</span>
            <span className="text-[10px] text-muted-foreground truncate">
              {metal.chainId}:{metal.resSeq}
            </span>
          </div>
          <button
            onClick={() => onFocus(isFocused ? null : index)}
            className={cn(
              'p-1 rounded transition-colors shrink-0',
              isFocused ? 'bg-primary text-primary-foreground' : 'hover:bg-muted text-muted-foreground'
            )}
            title={isFocused ? 'Clear focus' : 'Focus on metal'}
          >
            <Focus className="w-3 h-3" />
          </button>
        </div>

        <div className="flex flex-wrap gap-1 mt-1.5">
          <Badge variant="outline" className="text-[10px] h-5 px-1.5">
            {metal.bindingSiteType === 'functional' ? 'Functional' :
             metal.bindingSiteType === 'crystal_artifact' ? 'Artifact' : 'Uncertain'}
          </Badge>
          {metal.geometry && (
            <Badge variant="secondary" className="text-[10px] h-5 px-1.5">
              {metal.geometry.geometryType}
            </Badge>
          )}
          <Badge variant="secondary" className="text-[10px] h-5 px-1.5">
            CN: {metal.coordinating.length}
          </Badge>
        </div>
      </div>

      {expanded && (
        <div className="border-t border-border p-2 bg-muted/50">
          <div className="text-[10px] font-medium text-muted-foreground uppercase tracking-wide mb-1.5">
            Coordinating Atoms
          </div>
          <div className="space-y-0.5 max-h-20 overflow-y-auto">
            {metal.coordinating.slice(0, 6).map((coord, i) => (
              <div key={i} className="flex items-center justify-between text-[10px] px-1 py-0.5 rounded hover:bg-muted">
                <span className="text-foreground">
                  {coord.residue}
                  {coord.isWater && <span className="text-muted-foreground ml-1">(water)</span>}
                </span>
                <span className="font-mono text-muted-foreground">{coord.distance.toFixed(2)}Å</span>
              </div>
            ))}
            {metal.coordinating.length > 6 && (
              <div className="text-[10px] text-muted-foreground text-center py-0.5">
                +{metal.coordinating.length - 6} more
              </div>
            )}
          </div>
        </div>
      )}
    </div>
  );
}

// Ligand card component - minimal theme
function LigandCard({
  ligand,
  index,
  isFocused,
  onFocus,
}: {
  ligand: LigandData;
  index: number;
  isFocused: boolean;
  onFocus: (index: number | null) => void;
}) {
  const [expanded, setExpanded] = useState(false);

  return (
    <div
      className={cn(
        'rounded-md border transition-all',
        isFocused
          ? 'border-foreground/30 bg-accent'
          : 'border-border bg-card hover:bg-accent/50'
      )}
    >
      <div className="p-2">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-1.5 min-w-0">
            <button
              onClick={() => setExpanded(!expanded)}
              className="p-0.5 hover:bg-muted rounded text-muted-foreground shrink-0"
            >
              {expanded ? <ChevronDown className="w-3 h-3" /> : <ChevronRight className="w-3 h-3" />}
            </button>
            <Pill className="w-3.5 h-3.5 text-muted-foreground shrink-0" />
            <span className="font-medium text-xs text-foreground">{ligand.name}</span>
            <span className="text-[10px] text-muted-foreground truncate">
              {ligand.chainId}:{ligand.resSeq}
            </span>
          </div>
          <button
            onClick={() => onFocus(isFocused ? null : index)}
            className={cn(
              'p-1 rounded transition-colors shrink-0',
              isFocused ? 'bg-primary text-primary-foreground' : 'hover:bg-muted text-muted-foreground'
            )}
            title={isFocused ? 'Clear focus' : 'Focus on ligand'}
          >
            <Focus className="w-3 h-3" />
          </button>
        </div>

        <div className="flex flex-wrap gap-1 mt-1.5">
          <Badge variant="outline" className="text-[10px] h-5 px-1.5">
            {ligand.bindingSiteType === 'functional' ? 'Functional' :
             ligand.bindingSiteType === 'crystal_artifact' ? 'Artifact' : 'Uncertain'}
          </Badge>
          <Badge variant="secondary" className="text-[10px] h-5 px-1.5">
            {ligand.contacts.length} contacts
          </Badge>
        </div>
      </div>

      {expanded && (
        <div className="border-t border-border p-2 bg-muted/50">
          <div className="text-[10px] font-medium text-muted-foreground uppercase tracking-wide mb-1.5">
            Protein Contacts
          </div>
          <div className="space-y-0.5 max-h-20 overflow-y-auto">
            {ligand.contacts.slice(0, 6).map((contact, i) => (
              <div key={i} className="flex items-center justify-between text-[10px] px-1 py-0.5 rounded hover:bg-muted">
                <div className="flex items-center gap-1.5 min-w-0">
                  <span className="text-foreground">{contact.residue}</span>
                  <span className="text-muted-foreground truncate">{getInteractionLabel(contact.interactionType)}</span>
                </div>
                <span className="font-mono text-muted-foreground shrink-0">{contact.distance.toFixed(2)}Å</span>
              </div>
            ))}
            {ligand.contacts.length > 6 && (
              <div className="text-[10px] text-muted-foreground text-center py-0.5">
                +{ligand.contacts.length - 6} more
              </div>
            )}
          </div>
        </div>
      )}
    </div>
  );
}

export function ViewerPanel({
  children,
  isSpinning = false,
  onToggleSpin,
  onReset,
  onExpand,
  showControls = true,
}: ViewerPanelProps) {
  // Get state from store
  const {
    selectedPdb,
    focusedMetalIndex,
    setFocusedMetalIndex,
    focusedLigandIndex,
    setFocusedLigandIndex,
    metalCoordination,
    setMetalCoordination,
    ligandData,
    setLigandData,
    pharmacophoreFeatures,
    analysisLoading,
    setAnalysisLoading,
    catalyticSuggestions,
    suggestionsSource,
    suggestionsLoading,
    suggestionsError,
    bottomPanelMode,
    setBottomPanelMode,
    enzymeCatalyticResidues,
    enzymeLigandCodes,
    addEnzymeCatalyticResidue,
    removeEnzymeCatalyticResidue,
    setHoveredCatalyticSuggestion,
  } = useStore();

  // Filter catalytic suggestions based on current ligand codes
  const filteredCatalyticSuggestions = useMemo(() => {
    if (catalyticSuggestions.length === 0) return [];

    // Parse ligand codes into array (uppercase, trimmed)
    const codes = enzymeLigandCodes
      .split(',')
      .map(c => c.trim().toUpperCase())
      .filter(c => c.length > 0);

    // If no ligand codes specified, show all suggestions
    if (codes.length === 0) return catalyticSuggestions;

    // Filter suggestions to only those matching current ligand codes
    return catalyticSuggestions.filter(s => {
      // If suggestion has no ligandCode, include it (legacy/M-CSA data)
      if (!s.ligandCode) return true;
      // Include if suggestion's ligandCode matches any of the user's codes
      return codes.includes(s.ligandCode.toUpperCase());
    });
  }, [catalyticSuggestions, enzymeLigandCodes]);

  // Resizing state
  const [viewerHeight, setViewerHeight] = useState(65);
  const [isResizing, setIsResizing] = useState(false);
  const containerRef = useRef<HTMLDivElement>(null);

  // Section expanded state
  const [metalsOpen, setMetalsOpen] = useState(true);
  const [ligandsOpen, setLigandsOpen] = useState(true);

  // Focus mode
  const isInFocusMode = focusedMetalIndex !== null || focusedLigandIndex !== null;
  const focusType = focusedMetalIndex !== null ? 'metal' : 'ligand';

  const hasWaters = focusedMetalIndex !== null
    ? (metalCoordination?.[focusedMetalIndex]?.hydrationAnalysis?.waterCount ?? 0) > 0
    : (ligandData?.ligandDetails[focusedLigandIndex ?? 0]?.waterContactCount ?? 0) > 0;

  const hasInteractions = focusedLigandIndex !== null
    ? (ligandData?.ligandDetails[focusedLigandIndex]?.contacts?.length ?? 0) > 0
    : (metalCoordination?.[focusedMetalIndex ?? 0]?.coordinating?.length ?? 0) > 0;

  const hasPharmacophores = focusedLigandIndex !== null && (pharmacophoreFeatures?.length ?? 0) > 0;

  // Analysis data checks
  const hasMetals = metalCoordination && metalCoordination.length > 0;
  const hasLigands = ligandData && ligandData.ligandCount > 0;
  const hasAnalysis = hasMetals || hasLigands;

  // Panel priority: suggestions > metals > ligands > none
  const showSuggestions = bottomPanelMode === 'suggestions' || suggestionsLoading;
  const showAnalysis = !showSuggestions;

  // Run analysis
  const runAnalysis = useCallback(async () => {
    if (!selectedPdb || analysisLoading) return;

    setAnalysisLoading(true);

    try {
      const metals = findMetalCoordinationFromPDB(selectedPdb, 3.0);
      setMetalCoordination(metals);

      const ligands = findLigandContactsFromPDB(selectedPdb, 4.0);
      setLigandData(ligands);
    } catch (error) {
      console.error('Analysis failed:', error);
    } finally {
      setAnalysisLoading(false);
    }
  }, [selectedPdb, analysisLoading, setAnalysisLoading, setMetalCoordination, setLigandData]);

  // Resize handling
  const handleMouseDown = useCallback((e: React.MouseEvent) => {
    e.preventDefault();
    setIsResizing(true);
  }, []);

  useEffect(() => {
    const handleMouseMove = (e: MouseEvent) => {
      if (!isResizing || !containerRef.current) return;

      const containerRect = containerRef.current.getBoundingClientRect();
      const relativeY = e.clientY - containerRect.top;
      const percentage = (relativeY / containerRect.height) * 100;

      setViewerHeight(Math.min(85, Math.max(30, percentage)));
    };

    const handleMouseUp = () => {
      setIsResizing(false);
    };

    if (isResizing) {
      document.addEventListener('mousemove', handleMouseMove);
      document.addEventListener('mouseup', handleMouseUp);
      document.body.style.cursor = 'row-resize';
      document.body.style.userSelect = 'none';
    }

    return () => {
      document.removeEventListener('mousemove', handleMouseMove);
      document.removeEventListener('mouseup', handleMouseUp);
      document.body.style.cursor = '';
      document.body.style.userSelect = '';
    };
  }, [isResizing]);

  return (
    <div ref={containerRef} className="h-full flex flex-col bg-card">
      {/* Viewer Area */}
      <div
        className="relative bg-muted/30 overflow-hidden [&>*]:absolute [&>*]:inset-0"
        style={{ height: `${viewerHeight}%` }}
      >
        {children}

        {isInFocusMode && (
          <FocusModeControls
            focusType={focusType}
            hasWaters={hasWaters}
            hasInteractions={hasInteractions}
            hasPharmacophores={hasPharmacophores}
            onRadiusChange={() => {
              if (focusedMetalIndex !== null) {
                setFocusedMetalIndex(focusedMetalIndex);
              } else if (focusedLigandIndex !== null) {
                setFocusedLigandIndex(focusedLigandIndex);
              }
            }}
            onReset={() => {
              setFocusedMetalIndex(null);
              setFocusedLigandIndex(null);
            }}
          />
        )}
      </div>

      {/* Resize Handle */}
      <div
        onMouseDown={handleMouseDown}
        className={cn(
          'h-1.5 flex items-center justify-center cursor-row-resize group transition-colors border-y border-border',
          isResizing ? 'bg-accent' : 'hover:bg-accent/50'
        )}
      >
        <GripHorizontal className="h-3 w-3 text-muted-foreground/40 group-hover:text-muted-foreground" />
      </div>

      {/* Catalytic Suggestions Panel */}
      {showSuggestions && (
        <CatalyticSuggestionsPanel
          suggestions={filteredCatalyticSuggestions}
          source={suggestionsSource}
          loading={suggestionsLoading}
          error={suggestionsError}
          existingResidues={enzymeCatalyticResidues}
          onAdd={(s, atomType) => addEnzymeCatalyticResidue(s.chain, s.residue, s.name, atomType)}
          onRemove={(s) => removeEnzymeCatalyticResidue(s.chain, s.residue)}
          onAddAll={(atomType) => {
            filteredCatalyticSuggestions
              .filter((s) => !enzymeCatalyticResidues.some((r) => r.chain === s.chain && r.residue === s.residue))
              .forEach((s) => addEnzymeCatalyticResidue(s.chain, s.residue, s.name, atomType));
          }}
          onClose={() => setBottomPanelMode('none')}
          onHoverResidue={setHoveredCatalyticSuggestion}
        />
      )}

      {/* Analysis Panel */}
      {showAnalysis && (
      <div className="flex-1 flex flex-col overflow-hidden">
        {/* Header */}
        <div className="px-3 py-2 border-b border-border flex items-center justify-between shrink-0 bg-muted/30">
          <div className="flex items-center gap-2">
            <Microscope className="w-3.5 h-3.5 text-muted-foreground" />
            <span className="text-xs font-medium text-foreground">Analysis</span>
            {hasMetals && (
              <Badge variant="secondary" className="text-[10px] h-5 px-1.5">
                {metalCoordination!.length} metal{metalCoordination!.length !== 1 ? 's' : ''}
              </Badge>
            )}
            {hasLigands && (
              <Badge variant="secondary" className="text-[10px] h-5 px-1.5">
                {ligandData!.ligandCount} ligand{ligandData!.ligandCount !== 1 ? 's' : ''}
              </Badge>
            )}
          </div>
          <Button
            variant="ghost"
            size="sm"
            onClick={runAnalysis}
            disabled={!selectedPdb || analysisLoading}
            className="h-6 px-2 text-xs"
          >
            {analysisLoading ? (
              <Loader2 className="w-3 h-3 animate-spin" />
            ) : (
              'Analyze'
            )}
          </Button>
        </div>

        {/* Content */}
        <div className="flex-1 overflow-y-auto p-2 space-y-2">
          {!selectedPdb ? (
            <div className="flex flex-col items-center justify-center h-full text-muted-foreground">
              <Microscope className="w-6 h-6 mb-2 opacity-30" />
              <p className="text-xs">Load a structure to analyze</p>
            </div>
          ) : !hasAnalysis && !analysisLoading ? (
            <div className="flex flex-col items-center justify-center h-full text-muted-foreground">
              <p className="text-xs">Click Analyze to detect metals & ligands</p>
            </div>
          ) : analysisLoading ? (
            <div className="flex flex-col items-center justify-center h-full text-muted-foreground">
              <Loader2 className="w-5 h-5 mb-2 animate-spin opacity-50" />
              <p className="text-xs">Analyzing...</p>
            </div>
          ) : (
            <>
              {/* Metals */}
              {hasMetals && (
                <Collapsible open={metalsOpen} onOpenChange={setMetalsOpen}>
                  <CollapsibleTrigger className="flex items-center gap-1.5 w-full text-left py-1 group">
                    <ChevronRight className={cn(
                      "w-3 h-3 text-muted-foreground transition-transform",
                      metalsOpen && "rotate-90"
                    )} />
                    <Atom className="w-3.5 h-3.5 text-muted-foreground" />
                    <span className="text-xs font-medium text-foreground">Metals</span>
                    <span className="text-[10px] text-muted-foreground">({metalCoordination!.length})</span>
                  </CollapsibleTrigger>
                  <CollapsibleContent className="space-y-1.5 pt-1.5">
                    {metalCoordination!.map((metal, i) => (
                      <MetalCard
                        key={`${metal.chainId}-${metal.resSeq}-${i}`}
                        metal={metal}
                        index={i}
                        isFocused={focusedMetalIndex === i}
                        onFocus={(idx) => {
                          setFocusedLigandIndex(null);
                          setFocusedMetalIndex(idx);
                        }}
                      />
                    ))}
                  </CollapsibleContent>
                </Collapsible>
              )}

              {/* Ligands */}
              {hasLigands && (
                <Collapsible open={ligandsOpen} onOpenChange={setLigandsOpen}>
                  <CollapsibleTrigger className="flex items-center gap-1.5 w-full text-left py-1 group">
                    <ChevronRight className={cn(
                      "w-3 h-3 text-muted-foreground transition-transform",
                      ligandsOpen && "rotate-90"
                    )} />
                    <Pill className="w-3.5 h-3.5 text-muted-foreground" />
                    <span className="text-xs font-medium text-foreground">Ligands</span>
                    <span className="text-[10px] text-muted-foreground">({ligandData!.ligandCount})</span>
                  </CollapsibleTrigger>
                  <CollapsibleContent className="space-y-1.5 pt-1.5">
                    {ligandData!.ligandDetails.map((ligand, i) => (
                      <LigandCard
                        key={`${ligand.chainId}-${ligand.resSeq}-${i}`}
                        ligand={ligand}
                        index={i}
                        isFocused={focusedLigandIndex === i}
                        onFocus={(idx) => {
                          setFocusedMetalIndex(null);
                          setFocusedLigandIndex(idx);
                        }}
                      />
                    ))}
                  </CollapsibleContent>
                </Collapsible>
              )}

              {!hasMetals && !hasLigands && (
                <div className="flex items-center justify-center h-full text-muted-foreground">
                  <p className="text-xs">No metals or ligands detected</p>
                </div>
              )}
            </>
          )}
        </div>
      </div>
      )}

      {/* Bottom Controls */}
      {showControls && (
        <>
          <Separator />
          <div className="p-1.5 flex items-center gap-1 shrink-0 bg-muted/30">
            <Button variant="ghost" size="sm" onClick={onToggleSpin} className="h-7 px-2">
              {isSpinning ? <Pause className="h-3.5 w-3.5" /> : <Play className="h-3.5 w-3.5" />}
              <span className="ml-1 text-xs">Spin</span>
            </Button>
            <Button variant="ghost" size="sm" onClick={onReset} className="h-7 px-2">
              <RotateCcw className="h-3.5 w-3.5" />
              <span className="ml-1 text-xs">Reset</span>
            </Button>
            <div className="flex-1" />
            <Button variant="ghost" size="sm" onClick={onExpand} className="h-7 px-2">
              <Maximize2 className="h-3.5 w-3.5" />
            </Button>
          </div>
        </>
      )}
    </div>
  );
}

export default ViewerPanel;
