'use client';

import { useEffect, useRef, useState, useCallback, useImperativeHandle, forwardRef } from 'react';
import type { MetalCoordination } from '@/lib/metalAnalysis';
import type { LigandData, PharmacophoreFeature } from '@/lib/ligandAnalysis';
import { computeInteractionLines, type InteractionLine } from '@/lib/interactionGeometry';
import { useStore } from '@/lib/store';

// Import expression helpers from new modules
import {
  WATER_EXPRESSION,
  POLYMER_EXPRESSION,
  NON_POLYMER_EXPRESSION,
  residueExpression,
  surroundingsExpression,
  exceptExpression,
  intersectExpression,
} from '@/lib/molstar-expressions';
import { renderInteractionLinesBatched } from '@/lib/molstar-shapes';

// Mol* imports - static imports for core utilities
import type { PluginUIContext } from 'molstar/lib/mol-plugin-ui/context';
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder';
import { Script } from 'molstar/lib/mol-script/script';
import { StructureSelection } from 'molstar/lib/mol-model/structure';
import { Color } from 'molstar/lib/mol-util/color';
type StateObjectRef = any;

interface ProteinViewerClientProps {
  pdbContent: string | null;
  className?: string;
  // Focus props
  focusedMetalIndex?: number | null;
  focusedLigandIndex?: number | null;
  metalCoordination?: MetalCoordination[] | null;
  ligandData?: { ligandDetails: LigandData[] } | null;
  // Pharmacophore visualization
  pharmacophoreFeatures?: PharmacophoreFeature[];
  showPharmacophores?: boolean;
  // Callbacks
  onReady?: () => void;
}

// Pharmacophore feature colors
const PHARMACOPHORE_COLORS: Record<string, number> = {
  donor: 0x3B82F6,      // Blue
  acceptor: 0xEF4444,   // Red
  aromatic: 0xA855F7,   // Purple
  hydrophobic: 0x22C55E, // Green
  positive: 0x0EA5E9,   // Cyan
  negative: 0xF97316,   // Orange
};

export interface ProteinViewerHandle {
  resetView: () => void;
  focusOnMetal: (index: number) => Promise<void>;
  focusOnLigand: (index: number) => Promise<void>;
}

// Global state for Molstar plugin - persists across component remounts
// These are kept as globals to support the singleton plugin pattern
let globalPlugin: PluginUIContext | null = null;
let globalInitPromise: Promise<any> | null = null;
let globalContainer: HTMLDivElement | null = null;
let globalStructureRef: StateObjectRef | null = null;
let globalErrorHandler: ((event: ErrorEvent) => void) | null = null;
let globalUnhandledRejectionHandler: ((event: PromiseRejectionEvent) => void) | null = null;

export const ProteinViewerClient = forwardRef<ProteinViewerHandle, ProteinViewerClientProps>(
  function ProteinViewerClient(
    {
      pdbContent,
      className = '',
      focusedMetalIndex,
      focusedLigandIndex,
      metalCoordination,
      ligandData,
      pharmacophoreFeatures,
      showPharmacophores,
      onReady,
    },
    ref
  ) {
    const containerRef = useRef<HTMLDivElement>(null);
    const [loading, setLoading] = useState(false);
    const [error, setError] = useState<string | null>(null);
    const [isReady, setIsReady] = useState(false);
    const [focusMode, setFocusMode] = useState<'none' | 'metal' | 'ligand'>('none');
    const lastPdbContentRef = useRef<string | null>(null);
    const prevFocusRef = useRef<{ metalIndex: number | null | undefined; ligandIndex: number | null | undefined }>({
      metalIndex: undefined,
      ligandIndex: undefined,
    });

    // Get focus settings from store
    const { focusSettings } = useStore();

    // Focus on a metal site
    const focusOnMetal = useCallback(async (index: number) => {
      if (!globalPlugin || !metalCoordination || !metalCoordination[index] || !pdbContent) return;

      const plugin = globalPlugin;
      const metal = metalCoordination[index];

      try {
        console.log(`[ProteinViewer] Focusing on metal: ${metal.element} at ${metal.chainId}:${metal.resSeq}`);

        // Clear current view
        await plugin.clear();

        // Reload structure
        const isCif = pdbContent.trimStart().startsWith('data_');
        const format = isCif ? 'mmcif' : 'pdb';

        const data = await plugin.builders.data.rawData(
          { data: pdbContent, label: 'structure.pdb' },
          { state: { isGhost: true } }
        );

        const trajectory = await plugin.builders.structure.parseTrajectory(data, format);
        const model = await plugin.builders.structure.createModel(trajectory);
        const structure = await plugin.builders.structure.createStructure(model);
        globalStructureRef = structure.ref;

        // Build metal selection expression using helper
        const metalExpression = residueExpression(metal.resName, metal.chainId, metal.resSeq);

        // Coordination sphere (residues within radius + 2A margin)
        const coordSphereExpression = surroundingsExpression(
          metalExpression,
          focusSettings.coordinationRadius + 2.0,
          true
        );

        // Rest of protein (excluding coordination sphere)
        const proteinExpression = exceptExpression(POLYMER_EXPRESSION, coordSphereExpression);

        // 1. Metal ion - large purple spacefill
        const metalComp = await plugin.builders.structure.tryCreateComponentFromExpression(
          structure.ref,
          metalExpression,
          `metal-${metal.element}`
        );
        if (metalComp) {
          await plugin.builders.structure.representation.addRepresentation(metalComp, {
            type: 'spacefill',
            color: 'uniform',
            colorParams: { value: Color(0x7C3AED) }, // Purple
            typeParams: { sizeFactor: 1.5 },
          });
        }

        // 2. Coordination sphere - ball-and-stick with element colors
        const coordWithoutMetal = exceptExpression(coordSphereExpression, metalExpression);

        // Exclude waters for cleaner view
        const coordProteinExpression = exceptExpression(coordWithoutMetal, WATER_EXPRESSION);

        const coordComp = await plugin.builders.structure.tryCreateComponentFromExpression(
          structure.ref,
          coordProteinExpression,
          'coordination-sphere'
        );
        if (coordComp) {
          // Ball and stick
          await plugin.builders.structure.representation.addRepresentation(coordComp, {
            type: 'ball-and-stick',
            color: 'element-symbol',
            typeParams: { sizeFactor: 0.3 },
          });
          // Labels
          await plugin.builders.structure.representation.addRepresentation(coordComp, {
            type: 'label',
            typeParams: { level: 'residue' },
            color: 'uniform',
            colorParams: { value: Color(0x333333) },
          });
        }

        // 3. Coordinating waters - conditionally rendered
        if (focusSettings.showWaters) {
          const waterExpression = intersectExpression(coordWithoutMetal, WATER_EXPRESSION);
          const waterComp = await plugin.builders.structure.tryCreateComponentFromExpression(
            structure.ref,
            waterExpression,
            'coord-water'
          );
          if (waterComp) {
            await plugin.builders.structure.representation.addRepresentation(waterComp, {
              type: 'ball-and-stick',
              color: 'element-symbol',
              typeParams: { sizeFactor: 0.2, alpha: 0.7 },
            });
          }
        }

        // 4. Rest of protein - semi-transparent cartoon
        const proteinComp = await plugin.builders.structure.tryCreateComponentFromExpression(
          structure.ref,
          proteinExpression,
          'protein-background'
        );
        if (proteinComp) {
          await plugin.builders.structure.representation.addRepresentation(proteinComp, {
            type: 'cartoon',
            color: 'chain-id',
            typeParams: { alpha: 0.25 },
          });
        }

        // 5. Focus camera on coordination sphere
        const state = plugin.state.data;
        const cell = state.cells.get(structure.ref);
        if (cell?.obj?.data) {
          const focusStructure = cell.obj.data;
          const focusSelection = Script.getStructureSelection(coordSphereExpression, focusStructure);
          if (!StructureSelection.isEmpty(focusSelection)) {
            const loci = StructureSelection.toLociWithSourceUnits(focusSelection);
            plugin.managers.camera.focusLoci(loci, { durationMs: 400 });
          }
        }

        setFocusMode('metal');
        console.log(`[ProteinViewer] Focused on ${metal.element} coordination sphere`);

      } catch (err) {
        console.error('[ProteinViewer] Failed to focus on metal:', err);
      }
    }, [pdbContent, metalCoordination, focusSettings]);

    // Render interaction lines as 3D cylinders using batched renderer
    const renderInteractionLines = useCallback(async (lines: InteractionLine[]) => {
      if (!globalPlugin || lines.length === 0) return;
      await renderInteractionLinesBatched(globalPlugin, lines, { radius: 0.08 });
    }, []);

    // Focus on a ligand site
    const focusOnLigand = useCallback(async (index: number) => {
      if (!globalPlugin || !ligandData || !ligandData.ligandDetails[index] || !pdbContent) return;

      const plugin = globalPlugin;
      const ligand = ligandData.ligandDetails[index];

      try {
        console.log(`[ProteinViewer] Focusing on ligand: ${ligand.name} at ${ligand.chainId}:${ligand.resSeq}`);

        // Clear current view
        await plugin.clear();

        // Reload structure
        const isCif = pdbContent.trimStart().startsWith('data_');
        const format = isCif ? 'mmcif' : 'pdb';

        const data = await plugin.builders.data.rawData(
          { data: pdbContent, label: 'structure.pdb' },
          { state: { isGhost: true } }
        );

        const trajectory = await plugin.builders.structure.parseTrajectory(data, format);
        const model = await plugin.builders.structure.createModel(trajectory);
        const structure = await plugin.builders.structure.createStructure(model);
        globalStructureRef = structure.ref;

        // Build ligand selection expression using helper
        const ligandExpression = residueExpression(ligand.name, ligand.chainId, ligand.resSeq);

        // Binding site (residues within configurable radius)
        const bindingSiteExpression = surroundingsExpression(
          ligandExpression,
          focusSettings.bindingPocketRadius,
          true
        );

        // Rest of protein
        const proteinExpression = exceptExpression(POLYMER_EXPRESSION, bindingSiteExpression);

        // 1. Ligand - ball-and-stick with green carbons for visibility
        const ligandComp = await plugin.builders.structure.tryCreateComponentFromExpression(
          structure.ref,
          ligandExpression,
          `ligand-${ligand.name}`
        );
        if (ligandComp) {
          await plugin.builders.structure.representation.addRepresentation(ligandComp, {
            type: 'ball-and-stick',
            color: 'uniform',
            colorParams: { value: Color(focusSettings.ligandCarbonColor) },
            typeParams: { sizeFactor: 0.4 },
          });
        }

        // 2. Binding site residues
        const siteWithoutLigand = exceptExpression(bindingSiteExpression, ligandExpression);

        // Exclude waters
        const siteProteinExpression = exceptExpression(siteWithoutLigand, WATER_EXPRESSION);

        const siteComp = await plugin.builders.structure.tryCreateComponentFromExpression(
          structure.ref,
          siteProteinExpression,
          'binding-site'
        );
        if (siteComp) {
          await plugin.builders.structure.representation.addRepresentation(siteComp, {
            type: 'ball-and-stick',
            color: 'element-symbol',
            typeParams: { sizeFactor: 0.25 },
          });
          await plugin.builders.structure.representation.addRepresentation(siteComp, {
            type: 'label',
            typeParams: { level: 'residue' },
            color: 'uniform',
            colorParams: { value: Color(0x333333) },
          });
        }

        // 3. Binding site waters - conditionally rendered
        if (focusSettings.showWaters) {
          const waterExpression = intersectExpression(siteWithoutLigand, WATER_EXPRESSION);
          const waterComp = await plugin.builders.structure.tryCreateComponentFromExpression(
            structure.ref,
            waterExpression,
            'binding-water'
          );
          if (waterComp) {
            await plugin.builders.structure.representation.addRepresentation(waterComp, {
              type: 'ball-and-stick',
              color: 'element-symbol',
              typeParams: { sizeFactor: 0.15, alpha: 0.7 },
            });
          }
        }

        // 4. Rest of protein - semi-transparent cartoon
        const proteinComp = await plugin.builders.structure.tryCreateComponentFromExpression(
          structure.ref,
          proteinExpression,
          'protein-background'
        );
        if (proteinComp) {
          await plugin.builders.structure.representation.addRepresentation(proteinComp, {
            type: 'cartoon',
            color: 'chain-id',
            typeParams: { alpha: 0.25 },
          });
        }

        // 5. Render interaction lines if enabled
        if (focusSettings.showInteractionLines && ligand.contacts.length > 0) {
          const lineContacts = ligand.contacts.map(c => ({
            ligandChain: ligand.chainId,
            ligandResSeq: ligand.resSeq,
            ligandAtom: c.atom,
            proteinResidue: c.residue,
            proteinChain: c.chain,
            proteinAtom: c.atom,
            interactionType: c.interactionType,
            distance: c.distance,
          }));

          const lines = computeInteractionLines(pdbContent || '', lineContacts);
          await renderInteractionLines(lines);
        }

        // 6. Focus camera
        const state = plugin.state.data;
        const cell = state.cells.get(structure.ref);
        if (cell?.obj?.data) {
          const focusStructure = cell.obj.data;
          const focusSelection = Script.getStructureSelection(bindingSiteExpression, focusStructure);
          if (!StructureSelection.isEmpty(focusSelection)) {
            const loci = StructureSelection.toLociWithSourceUnits(focusSelection);
            plugin.managers.camera.focusLoci(loci, { durationMs: 400 });
          }
        }

        setFocusMode('ligand');
        console.log(`[ProteinViewer] Focused on ${ligand.name} binding site`);

      } catch (err) {
        console.error('[ProteinViewer] Failed to focus on ligand:', err);
      }
    }, [pdbContent, ligandData, focusSettings, renderInteractionLines]);

    // Clear all pharmacophore representations
    const clearPharmacophores = useCallback(async () => {
      if (!globalPlugin || !globalStructureRef) return;

      try {
        const state = globalPlugin.state.data;
        // Find and remove all pharmacophore components
        const toRemove: string[] = [];
        state.cells.forEach((cell: any, ref: string) => {
          if (cell.obj?.label?.startsWith('pharmacophore-')) {
            toRemove.push(ref);
          }
        });

        for (const ref of toRemove) {
          const cell = state.cells.get(ref);
          if (cell) {
            await globalPlugin.build().delete(ref).commit();
          }
        }

        if (toRemove.length > 0) {
          console.log(`[ProteinViewer] Cleared ${toRemove.length} pharmacophore components`);
        }
      } catch (err) {
        console.error('[ProteinViewer] Failed to clear pharmacophores:', err);
      }
    }, []);

    // Render pharmacophore features as colored spheres
    const renderPharmacophores = useCallback(async (features: PharmacophoreFeature[]) => {
      if (!globalPlugin || !globalStructureRef) {
        console.log('[ProteinViewer] Cannot render pharmacophores: plugin or structure not ready');
        return;
      }

      const plugin = globalPlugin;

      try {
        // Clear existing pharmacophores first
        await clearPharmacophores();

        console.log(`[ProteinViewer] Rendering ${features.length} pharmacophore features`);

        for (const feature of features) {
          // Extract residue number from feature.residue (e.g., "ASP45" -> 45)
          const resNumMatch = feature.residue.match(/\d+/);
          if (!resNumMatch) {
            console.warn(`[ProteinViewer] Could not extract residue number from ${feature.residue}`);
            continue;
          }
          const resNum = parseInt(resNumMatch[0], 10);

          // Extract residue name (e.g., "ASP45" -> "ASP")
          const resName = feature.residue.replace(/\d+/g, '');

          // Create Mol* selection expression for this residue/atom
          const pharmacophoreExpression = MS.struct.generator.atomGroups({
            'residue-test': MS.core.logic.and([
              MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_comp_id(), resName]),
              MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_asym_id(), feature.chain]),
              MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_seq_id(), resNum]),
            ]),
            'atom-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_atom_id(), feature.atom]),
          });

          // Get color for this pharmacophore type
          const color = PHARMACOPHORE_COLORS[feature.type] || 0x808080;

          // Create component and add spacefill representation
          const pharmacophoreComp = await plugin.builders.structure.tryCreateComponentFromExpression(
            globalStructureRef,
            pharmacophoreExpression,
            `pharmacophore-${feature.type}-${feature.residue}-${feature.atom}`
          );

          if (pharmacophoreComp) {
            await plugin.builders.structure.representation.addRepresentation(pharmacophoreComp, {
              type: 'spacefill',
              color: 'uniform',
              colorParams: { value: Color(color) },
              typeParams: { sizeFactor: 0.8, alpha: 0.7 },
            });
          }
        }

        console.log('[ProteinViewer] Pharmacophore features rendered');
      } catch (err) {
        console.error('[ProteinViewer] Failed to render pharmacophores:', err);
      }
    }, [clearPharmacophores]);

    // Reset view to default
    const resetView = useCallback(async () => {
      if (!globalPlugin || !pdbContent) return;

      try {
        await globalPlugin.clear();

        const isCif = pdbContent.trimStart().startsWith('data_');
        const format = isCif ? 'mmcif' : 'pdb';

        const data = await globalPlugin.builders.data.rawData(
          { data: pdbContent, label: 'structure.pdb' },
          { state: { isGhost: true } }
        );

        const trajectory = await globalPlugin.builders.structure.parseTrajectory(data, format);
        await globalPlugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');

        globalPlugin.canvas3d?.requestCameraReset();
        setFocusMode('none');
        console.log('[ProteinViewer] View reset to default');
      } catch (err) {
        console.error('[ProteinViewer] Failed to reset view:', err);
      }
    }, [pdbContent]);

    // Expose methods via ref
    useImperativeHandle(ref, () => ({
      resetView,
      focusOnMetal,
      focusOnLigand,
    }), [resetView, focusOnMetal, focusOnLigand]);

    // Handle container resize - notify Molstar to update canvas
    useEffect(() => {
      const container = containerRef.current;
      if (!container || !globalPlugin) return;

      const resizeObserver = new ResizeObserver(() => {
        // Debounce resize handling
        requestAnimationFrame(() => {
          if (globalPlugin?.canvas3d) {
            globalPlugin.canvas3d.handleResize();
          }
        });
      });

      resizeObserver.observe(container);

      return () => {
        resizeObserver.disconnect();
      };
    }, [isReady]);

    // Initialize Molstar on mount
    useEffect(() => {
      const container = containerRef.current;
      if (!container) return;

      // If plugin already exists and initialized, just mark ready
      if (globalPlugin) {
        console.log('[ProteinViewer] Reusing existing plugin');
        if (globalContainer !== container) {
          console.log('[ProteinViewer] Re-attaching to new container');
          if (globalContainer && globalContainer.firstChild) {
            container.innerHTML = '';
            while (globalContainer.firstChild) {
              container.appendChild(globalContainer.firstChild);
            }
          }
          globalContainer = container;
          // Trigger resize after re-attach
          requestAnimationFrame(() => {
            globalPlugin?.canvas3d?.handleResize();
          });
        }
        setIsReady(true);
        onReady?.();
        return;
      }

      // If already initializing, wait for completion
      if (globalInitPromise) {
        console.log('[ProteinViewer] Waiting for existing initialization...');
        globalInitPromise.then(() => {
          if (globalPlugin) {
            globalContainer = container;
            setIsReady(true);
            onReady?.();
          }
        }).catch(() => {
          setError('Failed to initialize 3D viewer');
        });
        return;
      }

      console.log('[ProteinViewer] Starting Molstar initialization...');
      globalContainer = container;

      globalInitPromise = (async () => {
        try {
          const molstarUI = await import('molstar/lib/mol-plugin-ui');
          const molstarReact = await import('molstar/lib/mol-plugin-ui/react18');
          const molstarSpec = await import('molstar/lib/mol-plugin-ui/spec');

          if (!globalContainer) {
            console.log('[ProteinViewer] No container available');
            return null;
          }

          globalContainer.innerHTML = '';

          // Set up global error handler to catch Molstar internal errors (e.g., LociHighlightManager)
          if (globalErrorHandler) {
            window.removeEventListener('error', globalErrorHandler);
          }
          globalErrorHandler = (event: ErrorEvent) => {
            // Suppress Molstar highlight manager errors that occur during hover
            const isMolstarError = event.filename?.includes('molstar') ||
                                   event.error?.stack?.includes('molstar') ||
                                   event.error?.stack?.includes('LociHighlightManager') ||
                                   event.error?.stack?.includes('normalizedLoci');
            if (event.message?.includes('Cannot read properties of undefined') && isMolstarError) {
              console.warn('[ProteinViewer] Suppressed Molstar internal error:', event.message);
              event.preventDefault();
              return true;
            }
          };
          window.addEventListener('error', globalErrorHandler);

          // Also handle unhandled promise rejections from Molstar
          if (globalUnhandledRejectionHandler) {
            window.removeEventListener('unhandledrejection', globalUnhandledRejectionHandler);
          }
          globalUnhandledRejectionHandler = (event: PromiseRejectionEvent) => {
            const errorString = String(event.reason?.stack || event.reason?.message || event.reason);
            if (errorString.includes('molstar') ||
                errorString.includes('LociHighlightManager') ||
                errorString.includes('normalizedLoci')) {
              console.warn('[ProteinViewer] Suppressed Molstar promise rejection:', event.reason);
              event.preventDefault();
            }
          };
          window.addEventListener('unhandledrejection', globalUnhandledRejectionHandler);

          console.log('[ProteinViewer] Creating plugin UI...');
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
              behaviors: [
                // Keep default behaviors but we'll wrap them for error handling
                ...(molstarSpec.DefaultPluginUISpec().behaviors || []),
              ],
            },
          });

          console.log('[ProteinViewer] Plugin created successfully');
          globalPlugin = plugin;

          // Wrap the highlight and interactivity managers to prevent crashes on hover
          const wrapInteractivityMethod = (obj: any, methodName: string, label: string) => {
            if (obj && typeof obj[methodName] === 'function') {
              const original = obj[methodName].bind(obj);
              obj[methodName] = (...args: any[]) => {
                try {
                  return original(...args);
                } catch (err) {
                  console.warn(`[ProteinViewer] ${label} error suppressed:`, err);
                }
              };
            }
          };

          // Wrap highlight manager methods
          if (plugin.managers?.interactivity?.lociHighlights) {
            const highlights = plugin.managers.interactivity.lociHighlights;
            wrapInteractivityMethod(highlights, 'highlightOnly', 'Highlight');
            wrapInteractivityMethod(highlights, 'highlight', 'Highlight');
            wrapInteractivityMethod(highlights, 'highlightOnlyExtend', 'HighlightExtend');
          }

          // Wrap selection manager methods
          if (plugin.managers?.interactivity?.lociSelects) {
            const selects = plugin.managers.interactivity.lociSelects;
            wrapInteractivityMethod(selects, 'select', 'Select');
            wrapInteractivityMethod(selects, 'selectOnly', 'SelectOnly');
            wrapInteractivityMethod(selects, 'selectToggle', 'SelectToggle');
          }

          // Wrap label manager methods (if available - may not exist in all Molstar versions)
          const interactivity = plugin.managers?.interactivity as any;
          if (interactivity?.lociLabels) {
            const labels = interactivity.lociLabels;
            wrapInteractivityMethod(labels, 'mark', 'Label');
            wrapInteractivityMethod(labels, 'markOnlyExtend', 'LabelExtend');
          }

          return plugin;
        } catch (err) {
          console.error('[ProteinViewer] Failed to initialize Molstar:', err);
          globalInitPromise = null;
          throw err;
        }
      })();

      globalInitPromise.then((plugin) => {
        if (plugin) {
          setIsReady(true);
          onReady?.();
        }
      }).catch(() => {
        setError('Failed to initialize 3D viewer');
      });
    }, [onReady]);

    // Load structure when pdbContent changes
    useEffect(() => {
      if (!pdbContent || !globalPlugin || !isReady) return;

      // Skip if same content
      if (lastPdbContentRef.current === pdbContent) return;
      lastPdbContentRef.current = pdbContent;

      // Capture plugin reference for async function (TypeScript narrowing)
      const plugin = globalPlugin;

      const loadStructure = async () => {
        setLoading(true);
        setError(null);
        setFocusMode('none');

        try {
          console.log('[ProteinViewer] Clearing existing structures...');
          await plugin.clear();

          const isCif = pdbContent.trimStart().startsWith('data_');
          const format = isCif ? 'mmcif' : 'pdb';
          const extension = isCif ? 'cif' : 'pdb';

          console.log(`[ProteinViewer] Detected format: ${format}`);

          const data = await plugin.builders.data.rawData(
            { data: pdbContent, label: `structure.${extension}` },
            { state: { isGhost: true } }
          );

          const trajectory = await plugin.builders.structure.parseTrajectory(data, format);
          const model = await plugin.builders.structure.createModel(trajectory);
          const structure = await plugin.builders.structure.createStructure(model);

          // Store structure reference
          globalStructureRef = structure.ref;

          console.log('[ProteinViewer] Structure created, adding cartoon representation...');

          // Create polymer component and add cartoon representation
          try {
            // Get polymer component using expression constant
            const polymerComp = await plugin.builders.structure.tryCreateComponentFromExpression(
              structure.ref,
              POLYMER_EXPRESSION,
              'polymer'
            );

            if (polymerComp) {
              await plugin.builders.structure.representation.addRepresentation(polymerComp, {
                type: 'cartoon',
                color: 'chain-id',
              });
              console.log('[ProteinViewer] Added cartoon representation with chain-id coloring');
            }

            // Add ball-and-stick for ligands/heteroatoms using expression constant
            const hetComp = await plugin.builders.structure.tryCreateComponentFromExpression(
              structure.ref,
              NON_POLYMER_EXPRESSION,
              'ligand'
            );

            if (hetComp) {
              await plugin.builders.structure.representation.addRepresentation(hetComp, {
                type: 'ball-and-stick',
                color: 'element-symbol',
              });
              console.log('[ProteinViewer] Added ball-and-stick for ligands');
            }
          } catch (reprErr) {
            console.log('[ProteinViewer] Custom representation failed, trying preset fallback:', reprErr);
            // Fallback to preset
            await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', {
              representationPreset: 'polymer-cartoon',
            });
          }

          console.log('[ProteinViewer] Representation applied');

          plugin.canvas3d?.requestCameraReset();
          console.log('[ProteinViewer] Structure loaded successfully');
        } catch (err) {
          console.error('[ProteinViewer] Failed to load structure:', err);
          setError(`Failed to load structure: ${err instanceof Error ? err.message : 'Unknown error'}`);
        } finally {
          setLoading(false);
        }
      };

      loadStructure();
    }, [pdbContent, isReady]);

    // Handle focus changes from parent - auto-focus when index changes
    useEffect(() => {
      if (!isReady || !globalPlugin || !pdbContent) return;

      const prev = prevFocusRef.current;
      const metalChanged = focusedMetalIndex !== prev.metalIndex;
      const ligandChanged = focusedLigandIndex !== prev.ligandIndex;

      // Update ref
      prevFocusRef.current = { metalIndex: focusedMetalIndex, ligandIndex: focusedLigandIndex };

      // Only act if something changed
      if (!metalChanged && !ligandChanged) return;

      const handleFocusChange = async () => {
        // Metal focus requested
        if (focusedMetalIndex !== null && focusedMetalIndex !== undefined && metalCoordination) {
          console.log(`[ProteinViewer] Auto-focusing on metal index ${focusedMetalIndex}`);
          await focusOnMetal(focusedMetalIndex);
        }
        // Ligand focus requested
        else if (focusedLigandIndex !== null && focusedLigandIndex !== undefined && ligandData) {
          console.log(`[ProteinViewer] Auto-focusing on ligand index ${focusedLigandIndex}`);
          await focusOnLigand(focusedLigandIndex);
        }
        // Focus cleared - reset to default view
        else if (focusedMetalIndex === null || focusedLigandIndex === null) {
          console.log('[ProteinViewer] Focus cleared, resetting view');
          await resetView();
        }
      };

      handleFocusChange();
    }, [focusedMetalIndex, focusedLigandIndex, isReady, pdbContent, metalCoordination, ligandData, focusOnMetal, focusOnLigand, resetView]);

    // Render pharmacophore features when enabled, clear when disabled
    useEffect(() => {
      if (showPharmacophores && pharmacophoreFeatures && pharmacophoreFeatures.length > 0) {
        renderPharmacophores(pharmacophoreFeatures);
      } else if (!showPharmacophores) {
        clearPharmacophores();
      }
    }, [showPharmacophores, pharmacophoreFeatures, renderPharmacophores, clearPharmacophores]);

    return (
      <div className={`relative bg-gray-100 rounded-lg overflow-hidden ${className}`}>
        <div ref={containerRef} className="w-full h-full" style={{ minHeight: '300px' }} />

        {!pdbContent && !error && !isReady && (
          <div className="absolute inset-0 flex items-center justify-center bg-gray-100 pointer-events-none">
            <div className="flex items-center gap-2 text-gray-600">
              <div className="w-5 h-5 border-2 border-gray-400 border-t-transparent rounded-full animate-spin" />
              Initializing viewer...
            </div>
          </div>
        )}

        {!pdbContent && !error && isReady && (
          <div className="absolute inset-0 flex items-center justify-center bg-gray-100 pointer-events-none">
            <p className="text-gray-500">No structure to display</p>
          </div>
        )}

        {loading && (
          <div className="absolute inset-0 flex items-center justify-center bg-white/80 pointer-events-none">
            <div className="flex items-center gap-2 text-gray-600">
              <div className="w-5 h-5 border-2 border-gray-400 border-t-transparent rounded-full animate-spin" />
              Loading structure...
            </div>
          </div>
        )}

        {error && (
          <div className="absolute inset-0 flex items-center justify-center bg-white/95">
            <div className="text-center p-4">
              <p className="text-red-600 mb-2">{error}</p>
              <p className="text-xs text-gray-500">Fallback: Showing raw PDB</p>
              <pre className="mt-2 text-xs text-gray-700 bg-gray-100 p-2 rounded max-h-48 overflow-auto text-left">
                {pdbContent?.slice(0, 2000)}...
              </pre>
            </div>
          </div>
        )}

        {pdbContent && !error && isReady && (
          <div className="absolute bottom-2 left-2 text-xs text-gray-600 bg-white/70 px-2 py-1 rounded pointer-events-none">
            Drag to rotate • Scroll to zoom • Shift+drag to pan
          </div>
        )}

        {/* Focus mode indicator */}
        {focusMode !== 'none' && (
          <div className="absolute top-2 left-2 text-xs font-medium text-white bg-purple-600 px-2 py-1 rounded shadow">
            {focusMode === 'metal' ? 'Metal Focus' : 'Ligand Focus'}
          </div>
        )}
      </div>
    );
  }
);

// Default export for dynamic import compatibility
export default ProteinViewerClient;
