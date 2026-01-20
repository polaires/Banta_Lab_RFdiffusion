'use client';

import { useRef, useEffect, useState, useCallback } from 'react';
import { Box, AlertCircle } from 'lucide-react';

// Import Molstar CSS
import 'molstar/lib/mol-plugin-ui/skin/light.scss';

// Force-include Molstar modules to prevent tree-shaking in production builds
import 'molstar/lib/mol-model-formats/structure/pdb';
import 'molstar/lib/mol-model-formats/structure/mmcif';
import 'molstar/lib/mol-io/reader/cif';
import 'molstar/lib/mol-io/reader/pdb';

// Types
type PluginUIContext = any;

interface BinderViewerInnerProps {
  pdbContent: string | null;
  targetChain?: string;
  binderChain?: string;
  interfaceResidues?: {
    target: number[];
    binder: number[];
  };
  highlightedResidues?: { chain: string; residue: number }[];
  className?: string;
  colorMode?: 'chain' | 'interface' | 'confidence';
  onResidueClick?: (chain: string, residue: number) => void;
}

// Colors
const COLORS = {
  target: 0x94A3B8,      // slate-400 (gray)
  binder: 0x8B5CF6,      // violet-500
  interface: 0xF59E0B,   // amber-500
  highlight: 0xEF4444,   // red-500
};

// Global state for Molstar plugin
let globalPlugin: PluginUIContext | null = null;
let globalInitPromise: Promise<any> | null = null;
let globalContainer: HTMLDivElement | null = null;

// Mol* imports
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

export function BinderViewerInner({
  pdbContent,
  targetChain = 'A',
  binderChain = 'B',
  interfaceResidues,
  highlightedResidues,
  className = '',
  colorMode = 'chain',
  onResidueClick,
}: BinderViewerInnerProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [isReady, setIsReady] = useState(false);
  const [viewMode, setViewMode] = useState<'cartoon' | 'surface' | 'both'>('cartoon');
  const lastPdbContentRef = useRef<string | null>(null);

  // Initialize Molstar
  useEffect(() => {
    const container = containerRef.current;
    if (!container) return;

    if (globalPlugin) {
      if (globalContainer !== container) {
        container.innerHTML = '';
        if (globalContainer?.firstChild) {
          while (globalContainer.firstChild) {
            container.appendChild(globalContainer.firstChild);
          }
        }
        globalContainer = container;
      }
      setIsReady(true);
      return;
    }

    if (globalInitPromise) {
      globalInitPromise.then(() => {
        if (globalPlugin) {
          globalContainer = container;
          setIsReady(true);
        }
      }).catch(() => setError('Failed to initialize viewer'));
      return;
    }

    globalContainer = container;

    globalInitPromise = (async () => {
      try {
        const molstarUI = await import('molstar/lib/mol-plugin-ui');
        const molstarReact = await import('molstar/lib/mol-plugin-ui/react18');
        const molstarSpec = await import('molstar/lib/mol-plugin-ui/spec');

        if (!globalContainer) return null;

        globalContainer.innerHTML = '';

        const plugin = await molstarUI.createPluginUI({
          target: globalContainer,
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

        globalPlugin = plugin;
        await loadMolstarModules();
        return plugin;
      } catch (err) {
        console.error('[BinderViewer] Failed to initialize:', err);
        globalInitPromise = null;
        throw err;
      }
    })();

    globalInitPromise.then((plugin) => {
      if (plugin) setIsReady(true);
    }).catch(() => setError('Failed to initialize viewer'));
  }, []);

  // Load structure with binder-specific coloring
  const loadBinderStructure = useCallback(async () => {
    if (!pdbContent || !globalPlugin || !isReady) return;

    setLoading(true);
    setError(null);

    try {
      // Ensure plugin canvas is fully initialized (critical for serverless/production)
      if (!globalPlugin.canvas3d) {
        console.warn('[BinderViewer] Waiting for canvas3d initialization...');
        await new Promise(resolve => setTimeout(resolve, 100));
        if (!globalPlugin.canvas3d) {
          throw new Error('Molstar canvas3d not initialized - WebGL may not be available');
        }
      }

      await globalPlugin.clear();
      await loadMolstarModules();

      const isCif = pdbContent.trimStart().startsWith('data_');
      const format = isCif ? 'mmcif' : 'pdb';

      const data = await globalPlugin.builders.data.rawData(
        { data: pdbContent, label: 'binder-complex.pdb' },
        { state: { isGhost: true } }
      );

      const trajectory = await globalPlugin.builders.structure.parseTrajectory(data, format);
      if (!trajectory) {
        throw new Error('Failed to parse trajectory - PDB data may be malformed');
      }

      // Use applyPreset with empty representation - this is bundler-safe
      // Then we'll add our custom representations
      await globalPlugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', {
        representationPreset: 'empty',
      });

      // Get the structure reference from the hierarchy
      const structures = globalPlugin.managers.structure.hierarchy.current.structures;
      if (!structures || structures.length === 0) {
        throw new Error('No structure found after preset application');
      }
      const structureCell = structures[0]?.cell;
      if (!structureCell?.transform?.ref) {
        throw new Error('Structure reference not available');
      }
      const structureRef = structureCell.transform.ref;

      // Target chain expression
      const targetExpression = MS.struct.generator.atomGroups({
        'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_asym_id(), targetChain])
      });

      // Binder chain expression
      const binderExpression = MS.struct.generator.atomGroups({
        'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_asym_id(), binderChain])
      });

      // 1. Target chain - gray cartoon
      const targetComp = await globalPlugin.builders.structure.tryCreateComponentFromExpression(
        structureRef,
        targetExpression,
        'target-chain'
      );

      if (targetComp) {
        if (viewMode === 'cartoon' || viewMode === 'both') {
          await globalPlugin.builders.structure.representation.addRepresentation(targetComp, {
            type: 'cartoon',
            color: 'uniform',
            colorParams: { value: Color(COLORS.target) },
          });
        }
        if (viewMode === 'surface' || viewMode === 'both') {
          await globalPlugin.builders.structure.representation.addRepresentation(targetComp, {
            type: 'molecular-surface',
            color: 'uniform',
            colorParams: { value: Color(COLORS.target) },
            typeParams: { alpha: viewMode === 'both' ? 0.3 : 0.7 },
          });
        }
      }

      // 2. Binder chain - violet with spectrum coloring
      const binderComp = await globalPlugin.builders.structure.tryCreateComponentFromExpression(
        structureRef,
        binderExpression,
        'binder-chain'
      );

      if (binderComp) {
        if (viewMode === 'cartoon' || viewMode === 'both') {
          await globalPlugin.builders.structure.representation.addRepresentation(binderComp, {
            type: 'cartoon',
            color: colorMode === 'chain' ? 'uniform' : 'residue-index',
            colorParams: colorMode === 'chain' ? { value: Color(COLORS.binder) } : undefined,
          });
        }
        if (viewMode === 'surface' || viewMode === 'both') {
          await globalPlugin.builders.structure.representation.addRepresentation(binderComp, {
            type: 'molecular-surface',
            color: colorMode === 'chain' ? 'uniform' : 'residue-index',
            colorParams: colorMode === 'chain' ? { value: Color(COLORS.binder) } : undefined,
            typeParams: { alpha: viewMode === 'both' ? 0.3 : 0.7 },
          });
        }
      }

      // 3. Interface residues highlighting (if provided)
      if (colorMode === 'interface' && interfaceResidues) {
        // Highlight target interface
        if (interfaceResidues.target.length > 0) {
          const targetInterfaceExpr = MS.struct.generator.atomGroups({
            'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_asym_id(), targetChain]),
            'residue-test': MS.core.set.has([
              MS.set(...interfaceResidues.target),
              MS.struct.atomProperty.macromolecular.auth_seq_id()
            ])
          });

          const targetIntComp = await globalPlugin.builders.structure.tryCreateComponentFromExpression(
            structureRef,
            targetInterfaceExpr,
            'target-interface'
          );

          if (targetIntComp) {
            await globalPlugin.builders.structure.representation.addRepresentation(targetIntComp, {
              type: 'ball-and-stick',
              color: 'uniform',
              colorParams: { value: Color(COLORS.interface) },
              typeParams: { sizeFactor: 0.3 },
            });
          }
        }

        // Highlight binder interface
        if (interfaceResidues.binder.length > 0) {
          const binderInterfaceExpr = MS.struct.generator.atomGroups({
            'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_asym_id(), binderChain]),
            'residue-test': MS.core.set.has([
              MS.set(...interfaceResidues.binder),
              MS.struct.atomProperty.macromolecular.auth_seq_id()
            ])
          });

          const binderIntComp = await globalPlugin.builders.structure.tryCreateComponentFromExpression(
            structureRef,
            binderInterfaceExpr,
            'binder-interface'
          );

          if (binderIntComp) {
            await globalPlugin.builders.structure.representation.addRepresentation(binderIntComp, {
              type: 'ball-and-stick',
              color: 'uniform',
              colorParams: { value: Color(COLORS.interface) },
              typeParams: { sizeFactor: 0.3 },
            });
          }
        }
      }

      // 4. Highlighted residues (from sequence viewer click, etc.)
      if (highlightedResidues && highlightedResidues.length > 0) {
        for (const hr of highlightedResidues) {
          const highlightExpr = MS.struct.generator.atomGroups({
            'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_asym_id(), hr.chain]),
            'residue-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_seq_id(), hr.residue])
          });

          const highlightComp = await globalPlugin.builders.structure.tryCreateComponentFromExpression(
            structureRef,
            highlightExpr,
            `highlight-${hr.chain}-${hr.residue}`
          );

          if (highlightComp) {
            await globalPlugin.builders.structure.representation.addRepresentation(highlightComp, {
              type: 'ball-and-stick',
              color: 'uniform',
              colorParams: { value: Color(COLORS.highlight) },
              typeParams: { sizeFactor: 0.5 },
            });
          }
        }
      }

      globalPlugin.canvas3d?.requestCameraReset();
      console.log('[BinderViewer] Binder structure loaded with chain coloring');
    } catch (err) {
      console.error('[BinderViewer] Failed to load structure:', err);
      setError(`Failed to load structure: ${err instanceof Error ? err.message : 'Unknown error'}`);
    } finally {
      setLoading(false);
    }
  }, [pdbContent, isReady, targetChain, binderChain, colorMode, viewMode, interfaceResidues, highlightedResidues]);

  // Reload when content or settings change
  useEffect(() => {
    if (pdbContent !== lastPdbContentRef.current || isReady) {
      lastPdbContentRef.current = pdbContent;
      loadBinderStructure();
    }
  }, [pdbContent, loadBinderStructure, isReady, colorMode, viewMode, highlightedResidues]);

  return (
    <div className={`bg-card rounded-xl border border-border overflow-hidden ${className}`}>
      {/* Header with controls */}
      <div className="bg-muted px-4 py-3 border-b border-border">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-2">
            <Box className="h-5 w-5 text-primary" />
            <h4 className="font-semibold text-foreground text-sm">3D Structure</h4>
          </div>

          {/* View mode controls */}
          <div className="flex items-center gap-2">
            <span className="text-xs text-muted-foreground">View:</span>
            <div className="flex rounded-lg border border-border overflow-hidden">
              <button
                onClick={() => setViewMode('cartoon')}
                className={`px-2 py-1 text-xs ${viewMode === 'cartoon' ? 'bg-violet-500 text-white' : 'bg-card text-muted-foreground hover:bg-muted/50'}`}
              >
                Cartoon
              </button>
              <button
                onClick={() => setViewMode('surface')}
                className={`px-2 py-1 text-xs border-l border-border ${viewMode === 'surface' ? 'bg-violet-500 text-white' : 'bg-card text-muted-foreground hover:bg-muted/50'}`}
              >
                Surface
              </button>
              <button
                onClick={() => setViewMode('both')}
                className={`px-2 py-1 text-xs border-l border-border ${viewMode === 'both' ? 'bg-violet-500 text-white' : 'bg-card text-muted-foreground hover:bg-muted/50'}`}
              >
                Both
              </button>
            </div>
          </div>
        </div>
      </div>

      {/* Viewer container */}
      <div className="relative" style={{ minHeight: '350px' }}>
        <div ref={containerRef} className="w-full h-full" style={{ minHeight: '350px' }} />

        {loading && (
          <div className="absolute inset-0 flex items-center justify-center bg-card/80">
            <div className="flex items-center gap-2 text-muted-foreground">
              <div className="w-5 h-5 border-2 border-border border-t-transparent rounded-full animate-spin" />
              Loading structure...
            </div>
          </div>
        )}

        {error && (
          <div className="absolute inset-0 flex items-center justify-center bg-red-50">
            <div className="text-center p-4">
              <AlertCircle className="h-8 w-8 text-red-500 mx-auto mb-2" />
              <p className="text-red-600 text-sm">{error}</p>
            </div>
          </div>
        )}

        {!pdbContent && !error && isReady && (
          <div className="absolute inset-0 flex items-center justify-center bg-muted/50">
            <p className="text-muted-foreground">No structure to display</p>
          </div>
        )}
      </div>

      {/* Legend */}
      <div className="px-4 py-2 border-t border-border bg-muted/50">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-4 text-xs">
            <div className="flex items-center gap-1">
              <div className="w-3 h-3 rounded" style={{ backgroundColor: `#${COLORS.target.toString(16).padStart(6, '0')}` }} />
              <span className="text-muted-foreground">Target ({targetChain})</span>
            </div>
            <div className="flex items-center gap-1">
              <div className="w-3 h-3 rounded" style={{ backgroundColor: `#${COLORS.binder.toString(16).padStart(6, '0')}` }} />
              <span className="text-muted-foreground">Binder ({binderChain})</span>
            </div>
            {colorMode === 'interface' && (
              <div className="flex items-center gap-1">
                <div className="w-3 h-3 rounded" style={{ backgroundColor: `#${COLORS.interface.toString(16).padStart(6, '0')}` }} />
                <span className="text-muted-foreground">Interface</span>
              </div>
            )}
          </div>
          <span className="text-xs text-muted-foreground">Drag to rotate • Scroll to zoom</span>
        </div>
      </div>
    </div>
  );
}
