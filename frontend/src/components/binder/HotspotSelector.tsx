'use client';

import { useState, useRef, useEffect, useCallback } from 'react';
import { Crosshair, Sparkles, Check, X } from 'lucide-react';
import api, { HotspotResidue } from '@/lib/api';

// Import Molstar CSS
import 'molstar/lib/mol-plugin-ui/skin/light.scss';

// Force-include Molstar modules to prevent tree-shaking in production builds
import 'molstar/lib/mol-model-formats/structure/pdb';
import 'molstar/lib/mol-model-formats/structure/mmcif';

// Types
type PluginUIContext = any;

interface HotspotSelectorProps {
  pdbContent: string;
  targetChain?: string;
  onConfirm: (hotspots: string[]) => void;
  onCancel: () => void;
  initialHotspots?: string[];
  maxHotspots?: number;
}

// Amino acid properties for sequence display
const AA_PROPERTIES: Record<string, { color: string; type: string }> = {
  // Hydrophobic
  A: { color: 'bg-amber-200', type: 'hydrophobic' },
  V: { color: 'bg-amber-200', type: 'hydrophobic' },
  L: { color: 'bg-amber-200', type: 'hydrophobic' },
  I: { color: 'bg-amber-200', type: 'hydrophobic' },
  M: { color: 'bg-amber-200', type: 'hydrophobic' },
  F: { color: 'bg-amber-300', type: 'hydrophobic' },
  Y: { color: 'bg-amber-300', type: 'hydrophobic' },
  W: { color: 'bg-amber-300', type: 'hydrophobic' },
  P: { color: 'bg-amber-200', type: 'hydrophobic' },
  // Polar
  S: { color: 'bg-blue-200', type: 'polar' },
  T: { color: 'bg-blue-200', type: 'polar' },
  N: { color: 'bg-blue-200', type: 'polar' },
  Q: { color: 'bg-blue-200', type: 'polar' },
  C: { color: 'bg-blue-200', type: 'polar' },
  G: { color: 'bg-muted', type: 'special' },
  // Positive
  K: { color: 'bg-red-200', type: 'positive' },
  R: { color: 'bg-red-200', type: 'positive' },
  H: { color: 'bg-red-200', type: 'positive' },
  // Negative
  D: { color: 'bg-purple-200', type: 'negative' },
  E: { color: 'bg-purple-200', type: 'negative' },
};

// Parse sequence from PDB
function parseSequenceFromPDB(pdbContent: string, chain: string): { seq: string; residues: { num: number; aa: string }[] } {
  const residues: { num: number; aa: string }[] = [];
  const seen = new Set<number>();

  const AA_MAP: Record<string, string> = {
    ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C',
    GLN: 'Q', GLU: 'E', GLY: 'G', HIS: 'H', ILE: 'I',
    LEU: 'L', LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P',
    SER: 'S', THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V',
  };

  for (const line of pdbContent.split('\n')) {
    if (!line.startsWith('ATOM')) continue;
    const lineChain = line.substring(21, 22).trim();
    if (lineChain !== chain) continue;

    const resNum = parseInt(line.substring(22, 26).trim());
    if (seen.has(resNum)) continue;
    seen.add(resNum);

    const resName = line.substring(17, 20).trim();
    const aa = AA_MAP[resName] || 'X';
    residues.push({ num: resNum, aa });
  }

  residues.sort((a, b) => a.num - b.num);
  return { seq: residues.map(r => r.aa).join(''), residues };
}

// Colors
const COLORS = {
  target: 0x94A3B8,      // slate-400 (gray)
  autoDetected: 0x14B8A6, // teal-500
  selected: 0x8B5CF6,    // violet-500
  hovered: 0xF59E0B,     // amber-500
};

// Mol* imports - loaded dynamically
let MS: any = null;
let Color: any = null;

async function loadMolstarModules() {
  if (MS) return;

  const [msModule, colorModule] = await Promise.all([
    import('molstar/lib/mol-script/language/builder'),
    import('molstar/lib/mol-util/color'),
  ]);

  MS = msModule.MolScriptBuilder;
  Color = colorModule.Color;
}

export function HotspotSelector({
  pdbContent,
  targetChain = 'A',
  onConfirm,
  onCancel,
  initialHotspots = [],
  maxHotspots = 3,  // BindCraft recommends 3-5 hotspots for focused binding
}: HotspotSelectorProps) {
  // State
  const [selectedHotspots, setSelectedHotspots] = useState<string[]>(initialHotspots);
  const [autoDetectedHotspots, setAutoDetectedHotspots] = useState<HotspotResidue[]>([]);
  const [hoveredResidue, setHoveredResidue] = useState<string | null>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [isDetecting, setIsDetecting] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [detectionMethod, setDetectionMethod] = useState<'exposed_clustered' | 'exposed' | 'interface_like'>('exposed_clustered');
  const [viewMode, setViewMode] = useState<'cartoon' | 'surface'>('surface');
  const [parsedSequence, setParsedSequence] = useState<{ seq: string; residues: { num: number; aa: string }[] }>({ seq: '', residues: [] });

  // Parse sequence on mount
  useEffect(() => {
    if (pdbContent) {
      const parsed = parseSequenceFromPDB(pdbContent, targetChain);
      setParsedSequence(parsed);
    }
  }, [pdbContent, targetChain]);

  // Molstar refs
  const containerRef = useRef<HTMLDivElement>(null);
  const pluginRef = useRef<PluginUIContext | null>(null);
  const structureRef = useRef<any>(null);
  const initializingRef = useRef(false);  // Track initialization state
  const [isReady, setIsReady] = useState(false);

  // Initialize Molstar
  useEffect(() => {
    const container = containerRef.current;
    // Prevent double initialization (React StrictMode)
    if (!container || pluginRef.current || initializingRef.current) return;

    initializingRef.current = true;

    (async () => {
      try {
        const molstarUI = await import('molstar/lib/mol-plugin-ui');
        const molstarReact = await import('molstar/lib/mol-plugin-ui/react18');
        const molstarSpec = await import('molstar/lib/mol-plugin-ui/spec');

        // Clear any existing content
        while (container.firstChild) {
          container.removeChild(container.firstChild);
        }

        const plugin = await molstarUI.createPluginUI({
          target: container,
          render: molstarReact.renderReact18,
          spec: {
            ...molstarSpec.DefaultPluginUISpec(),
            layout: {
              initial: {
                isExpanded: false,
                showControls: false,
                controlsDisplay: 'reactive' as const,
              },
            },
          },
        });

        pluginRef.current = plugin;
        await loadMolstarModules();
        setIsReady(true);
      } catch (err) {
        console.error('[HotspotSelector] Failed to initialize:', err);
        setError('Failed to initialize 3D viewer');
        initializingRef.current = false;
      }
    })();

    return () => {
      // Cleanup on unmount
      if (pluginRef.current) {
        try {
          pluginRef.current.dispose();
        } catch (e) {
          // Ignore disposal errors
        }
        pluginRef.current = null;
      }
      initializingRef.current = false;
      setIsReady(false);
    };
  }, []);

  // Load structure
  const loadStructure = useCallback(async () => {
    const plugin = pluginRef.current;
    if (!plugin || !pdbContent || !isReady) return;

    setIsLoading(true);
    setError(null);

    try {
      // Ensure plugin canvas is fully initialized (critical for serverless/production)
      if (!plugin.canvas3d) {
        console.warn('[HotspotSelector] Waiting for canvas3d initialization...');
        await new Promise(resolve => setTimeout(resolve, 100));
        if (!plugin.canvas3d) {
          throw new Error('Molstar canvas3d not initialized - WebGL may not be available');
        }
      }

      await loadMolstarModules();

      const isCif = pdbContent.trimStart().startsWith('data_');
      const format = isCif ? 'mmcif' : 'pdb';

      // Wrap in dataTransaction to prevent Molstar hierarchy builder from crashing
      // on partially-formed state cells in production mode
      await plugin.dataTransaction(async () => {
        await plugin.clear();
        const data = await plugin.builders.data.rawData(
          { data: pdbContent, label: 'target.pdb' },
          { state: { isGhost: true } }
        );
        const trajectory = await plugin.builders.structure.parseTrajectory(data, format);
        if (!trajectory) {
          throw new Error('Failed to parse trajectory - PDB data may be malformed');
        }
        await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', {
          representationPreset: 'empty',
        });
      });

      // Get the structure reference from the hierarchy
      const structures = plugin.managers.structure.hierarchy.current.structures;
      if (!structures || structures.length === 0) {
        throw new Error('No structure found after preset application');
      }
      const structureCell = structures[0]?.cell;
      if (!structureCell?.transform?.ref) {
        throw new Error('Structure reference not available');
      }
      const structureRefValue = structureCell.transform.ref;
      structureRef.current = { ref: structureRefValue };

      // All residues as gray base
      const allExpr = MS.struct.generator.atomGroups({
        'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_asym_id(), targetChain])
      });

      const allComp = await plugin.builders.structure.tryCreateComponentFromExpression(
        structureRefValue,
        allExpr,
        'all-residues'
      );

      if (allComp) {
        await plugin.builders.structure.representation.addRepresentation(allComp, {
          type: viewMode === 'surface' ? 'molecular-surface' : 'cartoon',
          color: 'uniform',
          colorParams: { value: Color(COLORS.target) },
          typeParams: viewMode === 'surface' ? { alpha: 0.7 } : undefined,
        });
      }

      plugin.canvas3d?.requestCameraReset();
      console.log('[HotspotSelector] Structure loaded');
    } catch (err) {
      console.error('[HotspotSelector] Failed to load structure:', err);
      setError(`Failed to load structure: ${err instanceof Error ? err.message : 'Unknown error'}`);
    } finally {
      setIsLoading(false);
    }
  }, [pdbContent, isReady, targetChain, viewMode]);

  // Update highlights when selection changes
  const updateHighlights = useCallback(async () => {
    const plugin = pluginRef.current;
    const structure = structureRef.current;

    if (!plugin || !structure || !isReady) {
      console.log('[HotspotSelector] Skipping highlights - not ready:', { plugin: !!plugin, structure: !!structure, isReady });
      return;
    }

    // Ensure Molstar modules are loaded
    if (!MS || !Color) {
      await loadMolstarModules();
      if (!MS || !Color) {
        console.warn('[HotspotSelector] Molstar modules not loaded yet');
        return;
      }
    }

    console.log('[HotspotSelector] Updating highlights:', {
      autoDetected: autoDetectedHotspots.length,
      selected: selectedHotspots.length,
      hovered: hoveredResidue
    });

    try {
      // Remove old highlights by iterating over cells
      const state = plugin.state.data;
      const toRemove: string[] = [];
      state.cells.forEach((cell: any, ref: string) => {
        const label = cell.obj?.label || '';
        if (label.startsWith('hotspot-') || label.startsWith('hovered-')) {
          toRemove.push(ref);
        }
      });
      for (const ref of toRemove) {
        try {
          await plugin.build().delete(ref).commit();
        } catch (e) {
          // Ignore deletion errors
        }
      }

      // Add auto-detected hotspot highlights (teal) - use spacefill for visibility
      const autoDetectedIds = autoDetectedHotspots.map(h => h.residue);
      const unselectedAuto = autoDetectedIds.filter(id => !selectedHotspots.includes(id));

      if (unselectedAuto.length > 0) {
        const resNums = unselectedAuto.map(id => parseInt(id.slice(1)));
        const autoExpr = MS.struct.generator.atomGroups({
          'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_asym_id(), targetChain]),
          'residue-test': MS.core.set.has([
            MS.set(...resNums),
            MS.struct.atomProperty.macromolecular.auth_seq_id()
          ])
        });

        const autoComp = await plugin.builders.structure.tryCreateComponentFromExpression(
          structure.ref,
          autoExpr,
          'hotspot-auto'
        );

        if (autoComp) {
          // Use spacefill with slight transparency for auto-detected
          await plugin.builders.structure.representation.addRepresentation(autoComp, {
            type: 'spacefill',
            color: 'uniform',
            colorParams: { value: Color(COLORS.autoDetected) },
            typeParams: { sizeFactor: 1.0 },
          });
        }
      }

      // Add selected hotspot highlights (violet) - larger and more visible
      if (selectedHotspots.length > 0) {
        const resNums = selectedHotspots.map(id => parseInt(id.slice(1)));
        const selectedExpr = MS.struct.generator.atomGroups({
          'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_asym_id(), targetChain]),
          'residue-test': MS.core.set.has([
            MS.set(...resNums),
            MS.struct.atomProperty.macromolecular.auth_seq_id()
          ])
        });

        const selectedComp = await plugin.builders.structure.tryCreateComponentFromExpression(
          structure.ref,
          selectedExpr,
          'hotspot-selected'
        );

        if (selectedComp) {
          // Use spacefill for selected residues - most prominent
          await plugin.builders.structure.representation.addRepresentation(selectedComp, {
            type: 'spacefill',
            color: 'uniform',
            colorParams: { value: Color(COLORS.selected) },
            typeParams: { sizeFactor: 1.1 },  // Slightly larger
          });
        }
      }

      // Add hovered residue highlight (amber)
      if (hoveredResidue) {
        const resNum = parseInt(hoveredResidue.slice(1));
        const hoverExpr = MS.struct.generator.atomGroups({
          'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_asym_id(), targetChain]),
          'residue-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_seq_id(), resNum])
        });

        const hoverComp = await plugin.builders.structure.tryCreateComponentFromExpression(
          structure.ref,
          hoverExpr,
          'hovered-residue'
        );

        if (hoverComp) {
          await plugin.builders.structure.representation.addRepresentation(hoverComp, {
            type: 'spacefill',
            color: 'uniform',
            colorParams: { value: Color(COLORS.hovered) },
            typeParams: { sizeFactor: 0.4 },
          });
        }
      }
    } catch (err) {
      const errMsg = err instanceof Error ? err.message : String(err);
      console.error('[HotspotSelector] Failed to update highlights:', errMsg, err);
    }
  }, [autoDetectedHotspots, selectedHotspots, hoveredResidue, targetChain, isReady]);

  // Load structure on mount
  useEffect(() => {
    loadStructure();
  }, [loadStructure]);

  // Update highlights when selection changes
  useEffect(() => {
    updateHighlights();
  }, [updateHighlights]);

  // Auto-detect hotspots
  const detectHotspots = async () => {
    setIsDetecting(true);
    setError(null);

    try {
      const result = await api.detectHotspots({
        target_pdb: pdbContent,
        target_chain: targetChain,
        method: detectionMethod,
        max_hotspots: maxHotspots,
        prefer_hydrophobic: true,
      });

      if (result.status === 'completed' && result.residue_details) {
        // Store all detected for display, but only select up to maxHotspots
        setAutoDetectedHotspots(result.residue_details);
        // Auto-select the detected hotspots (limited by maxHotspots from backend)
        const hotspotsToSelect = result.hotspots.slice(0, maxHotspots);
        setSelectedHotspots(hotspotsToSelect);
        console.log('[HotspotSelector] Detected hotspots:', hotspotsToSelect);

        // Force highlight update after short delay to ensure state is updated
        setTimeout(() => {
          updateHighlights();
        }, 100);
      } else {
        throw new Error(result.error || 'Detection failed');
      }
    } catch (err) {
      console.error('[HotspotSelector] Detection failed:', err);
      setError(`Detection failed: ${err instanceof Error ? err.message : 'Unknown error'}`);
    } finally {
      setIsDetecting(false);
    }
  };

  // Toggle hotspot selection
  const toggleHotspot = (residueId: string) => {
    setSelectedHotspots(prev => {
      if (prev.includes(residueId)) {
        return prev.filter(id => id !== residueId);
      } else if (prev.length < maxHotspots) {
        return [...prev, residueId];
      }
      return prev;
    });
  };

  // Get property color
  const getPropertyColor = (property: string) => {
    switch (property) {
      case 'hydrophobic': return 'bg-muted text-foreground';
      case 'polar': return 'bg-blue-100 text-blue-700';
      case 'charged': return 'bg-red-100 text-red-700';
      default: return 'bg-muted text-foreground';
    }
  };

  return (
    <div className="bg-card rounded-xl border border-border overflow-hidden shadow-lg">
      {/* Header */}
      <div className="bg-muted px-4 py-3 border-b border-border">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-2">
            <Crosshair className="h-5 w-5 text-primary" />
            <h3 className="font-semibold text-foreground">Select Binding Hotspots</h3>
          </div>
          <span className="text-xs text-muted-foreground">
            {selectedHotspots.length}/{maxHotspots} selected
          </span>
        </div>
        <p className="text-sm text-muted-foreground mt-1">
          Choose residues on the target where the binder should attach
        </p>
      </div>

      <div className="grid grid-cols-2 gap-0">
        {/* Left: 3D Viewer */}
        <div className="border-r border-border">
          {/* View controls */}
          <div className="px-3 py-2 border-b border-border flex items-center justify-between bg-muted/50">
            <div className="flex items-center gap-2">
              <span className="text-xs text-muted-foreground">View:</span>
              <div className="flex rounded-md border border-border overflow-hidden">
                <button
                  onClick={() => { setViewMode('cartoon'); loadStructure(); }}
                  className={`px-2 py-1 text-xs ${viewMode === 'cartoon' ? 'bg-teal-500 text-white' : 'bg-card text-muted-foreground'}`}
                >
                  Cartoon
                </button>
                <button
                  onClick={() => { setViewMode('surface'); loadStructure(); }}
                  className={`px-2 py-1 text-xs border-l border-border ${viewMode === 'surface' ? 'bg-teal-500 text-white' : 'bg-card text-muted-foreground'}`}
                >
                  Surface
                </button>
              </div>
            </div>
          </div>

          {/* Molstar container */}
          <div className="relative" style={{ height: '350px' }}>
            <div ref={containerRef} className="w-full h-full" />

            {isLoading && (
              <div className="absolute inset-0 flex items-center justify-center bg-card/80">
                <div className="flex items-center gap-2 text-muted-foreground">
                  <div className="w-5 h-5 border-2 border-border border-t-transparent rounded-full animate-spin" />
                  Loading...
                </div>
              </div>
            )}
          </div>

          {/* Legend */}
          <div className="px-3 py-2 border-t border-border bg-muted/50">
            <div className="flex items-center gap-4 text-xs">
              <div className="flex items-center gap-1">
                <div className="w-3 h-3 rounded" style={{ backgroundColor: '#14B8A6' }} />
                <span className="text-muted-foreground">Auto-detected</span>
              </div>
              <div className="flex items-center gap-1">
                <div className="w-3 h-3 rounded" style={{ backgroundColor: '#8B5CF6' }} />
                <span className="text-muted-foreground">Selected</span>
              </div>
            </div>
          </div>
        </div>

        {/* Right: Selection panel */}
        <div className="flex flex-col">
          {/* Auto-detect controls */}
          <div className="px-4 py-3 border-b border-border">
            <div className="flex items-center gap-2 mb-2">
              <select
                value={detectionMethod}
                onChange={(e) => setDetectionMethod(e.target.value as any)}
                className="text-sm border border-border rounded-md px-2 py-1.5 flex-1"
              >
                <option value="exposed_clustered">Clustered Surface (Recommended)</option>
                <option value="exposed">All Exposed Residues</option>
                <option value="interface_like">Interface-like Residues</option>
              </select>
              <button
                onClick={detectHotspots}
                disabled={isDetecting}
                className="px-3 py-1.5 bg-teal-500 text-white text-sm font-medium rounded-md hover:bg-teal-600 disabled:opacity-50 flex items-center gap-1"
              >
                {isDetecting ? (
                  <>
                    <div className="w-4 h-4 border-2 border-white border-t-transparent rounded-full animate-spin" />
                    Detecting...
                  </>
                ) : (
                  <>
                    <Sparkles className="h-4 w-4" />
                    Detect
                  </>
                )}
              </button>
            </div>

            {error && (
              <div className="text-xs text-red-600 bg-red-50 px-2 py-1 rounded">
                {error}
              </div>
            )}
          </div>

          {/* Sequence Minimap for Manual Selection */}
          <div className="px-3 py-2 border-b border-border">
            <div className="flex items-center justify-between mb-2">
              <span className="text-xs font-medium text-muted-foreground">Sequence (click to select)</span>
              <span className="text-xs text-muted-foreground">{parsedSequence.residues.length} residues</span>
            </div>
            <div className="flex flex-wrap gap-0.5 max-h-24 overflow-y-auto p-1 bg-muted/50 rounded">
              {parsedSequence.residues.map((res) => {
                const residueId = `${targetChain}${res.num}`;
                const isSelected = selectedHotspots.includes(residueId);
                const isAutoDetected = autoDetectedHotspots.some(h => h.residue === residueId);
                const isHovered = hoveredResidue === residueId;
                const aaProps = AA_PROPERTIES[res.aa] || { color: 'bg-muted', type: 'unknown' };

                return (
                  <button
                    key={residueId}
                    onClick={() => toggleHotspot(residueId)}
                    onMouseEnter={() => setHoveredResidue(residueId)}
                    onMouseLeave={() => setHoveredResidue(null)}
                    title={`${residueId} (${res.aa}) - ${aaProps.type}`}
                    className={`w-4 h-4 text-[8px] font-mono font-bold rounded-sm transition-all ${
                      isSelected
                        ? 'bg-violet-500 text-white ring-2 ring-violet-300'
                        : isAutoDetected
                        ? 'bg-teal-400 text-white ring-1 ring-teal-300'
                        : isHovered
                        ? 'bg-amber-400 text-white'
                        : aaProps.color + ' text-muted-foreground hover:ring-1 hover:ring-border'
                    }`}
                  >
                    {res.aa}
                  </button>
                );
              })}
            </div>
            {/* Sequence legend */}
            <div className="flex items-center gap-2 mt-2 text-[10px] text-muted-foreground">
              <span className="flex items-center gap-0.5"><span className="w-2 h-2 bg-amber-200 rounded-sm" /> Hydrophobic</span>
              <span className="flex items-center gap-0.5"><span className="w-2 h-2 bg-blue-200 rounded-sm" /> Polar</span>
              <span className="flex items-center gap-0.5"><span className="w-2 h-2 bg-red-200 rounded-sm" /> +</span>
              <span className="flex items-center gap-0.5"><span className="w-2 h-2 bg-purple-200 rounded-sm" /> -</span>
            </div>
          </div>

          {/* Detected Hotspots list */}
          <div className="flex-1 overflow-y-auto" style={{ maxHeight: '180px' }}>
            {autoDetectedHotspots.length === 0 ? (
              <div className="flex flex-col items-center justify-center h-full text-muted-foreground p-4">
                <Sparkles className="h-8 w-8 mb-2" />
                <p className="text-xs text-center">
                  Click "Detect" for AI suggestions<br />or select from sequence above
                </p>
              </div>
            ) : (
              <div className="divide-y divide-border">
                {autoDetectedHotspots.map((hotspot) => {
                  const isSelected = selectedHotspots.includes(hotspot.residue);
                  const isHovered = hoveredResidue === hotspot.residue;
                  // SASA is relative (0=buried, 1=avg exposed, >1=very exposed)
                  // Always multiply by 100 for percentage display
                  const sasaPercent = hotspot.relative_sasa * 100;

                  return (
                    <button
                      key={hotspot.residue}
                      onClick={() => toggleHotspot(hotspot.residue)}
                      onMouseEnter={() => setHoveredResidue(hotspot.residue)}
                      onMouseLeave={() => setHoveredResidue(null)}
                      className={`w-full px-3 py-1.5 text-left flex items-center justify-between transition-colors ${
                        isSelected
                          ? 'bg-violet-50'
                          : isHovered
                          ? 'bg-amber-50'
                          : 'hover:bg-muted/50'
                      }`}
                    >
                      <div className="flex items-center gap-2">
                        <div className={`w-4 h-4 rounded border-2 flex items-center justify-center ${
                          isSelected ? 'border-primary bg-primary' : 'border-border'
                        }`}>
                          {isSelected && (
                            <Check className="h-2.5 w-2.5 text-white" />
                          )}
                        </div>
                        <div>
                          <div className="flex items-center gap-1.5">
                            <span className="font-mono font-medium text-sm text-foreground">{hotspot.residue}</span>
                            <span className="text-xs text-muted-foreground">{hotspot.restype}</span>
                          </div>
                          <div className="text-[10px] text-muted-foreground">
                            SASA: {sasaPercent.toFixed(0)}%
                          </div>
                        </div>
                      </div>
                      <span className={`text-[10px] px-1.5 py-0.5 rounded-full ${getPropertyColor(hotspot.property)}`}>
                        {hotspot.property}
                      </span>
                    </button>
                  );
                })}
              </div>
            )}
          </div>

          {/* Selected summary */}
          <div className="px-3 py-2 border-t border-border bg-muted/50">
            <div className="flex items-center justify-between mb-1">
              <span className="text-xs font-medium text-muted-foreground">Selected ({selectedHotspots.length}/{maxHotspots})</span>
              {selectedHotspots.length > 0 && (
                <button
                  onClick={() => setSelectedHotspots([])}
                  className="text-[10px] text-muted-foreground hover:text-muted-foreground"
                >
                  Clear all
                </button>
              )}
            </div>
            <div className="flex flex-wrap gap-1 min-h-[24px]">
              {selectedHotspots.length === 0 ? (
                <span className="text-xs text-muted-foreground italic">No hotspots selected</span>
              ) : (
                selectedHotspots.map(id => (
                  <span
                    key={id}
                    className="text-[10px] px-1.5 py-0.5 bg-primary/10 text-primary rounded flex items-center gap-0.5"
                  >
                    {id}
                    <button
                      onClick={() => toggleHotspot(id)}
                      className="hover:text-primary"
                    >
                      <X className="h-2.5 w-2.5" />
                    </button>
                  </span>
                ))
              )}
            </div>
          </div>
        </div>
      </div>

      {/* Footer */}
      <div className="px-4 py-3 border-t border-border bg-muted/50 flex items-center justify-between">
        <button
          onClick={onCancel}
          className="px-4 py-2 text-muted-foreground text-sm font-medium hover:text-foreground"
        >
          Cancel
        </button>
        <button
          onClick={() => onConfirm(selectedHotspots)}
          disabled={selectedHotspots.length === 0}
          className="px-5 py-2 bg-teal-600 text-white text-sm font-medium rounded-lg hover:bg-teal-700 disabled:opacity-50 flex items-center gap-2"
        >
          Confirm {selectedHotspots.length} Hotspot{selectedHotspots.length !== 1 ? 's' : ''}
          <Check className="h-4 w-4" />
        </button>
      </div>
    </div>
  );
}
