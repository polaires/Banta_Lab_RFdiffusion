'use client';

import { useEffect, useRef, useCallback } from 'react';
import { useStore, type HotspotResidue } from '@/lib/store';
import { Color } from 'molstar/lib/mol-util/color';
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder';
import { Script } from 'molstar/lib/mol-script/script';
import { StructureSelection } from 'molstar/lib/mol-model/structure';
import type { PluginUIContext } from 'molstar/lib/mol-plugin-ui/context';
import {
  POLYMER_EXPRESSION,
  surroundingsExpression,
  exceptExpression,
} from '@/lib/molstar-expressions';

// Colors for hotspot visualization (warm colors to indicate "hot" spots)
const HOTSPOT_COLOR = Color(0xe11d48);       // Rose-600 - primary hotspot color

// Component labels for cleanup
const HOTSPOT_COMPONENT_LABEL = 'hotspot-residues';
const HOTSPOT_BACKGROUND_LABEL = 'hotspot-background-cartoon';

/**
 * Parse hotspot residue ID (e.g., "A25") into chain and number
 */
function parseResidueId(residueId: string): { chain: string; residue: number } | null {
  // Format: "A25" - single letter chain followed by residue number
  if (!residueId || residueId.length < 2) return null;

  const chain = residueId[0];
  const residueStr = residueId.slice(1);
  const residue = parseInt(residueStr, 10);

  if (isNaN(residue)) return null;

  return { chain, residue };
}

/**
 * Hook to visualize detected hotspots in Molstar viewer.
 * - Shows hotspot residues as rose-colored ball-and-stick when showHotspots3D is true
 * - Handles hover highlighting when user hovers over a hotspot in a panel
 * - Adds transparent background cartoon for context
 * - Retries if structure is not yet loaded (up to 5 attempts)
 */
export function useHotspotVisualization(plugin: PluginUIContext | null) {
  const hotspotsData = useStore((s) => s.hotspotsData);
  const showHotspots3D = useStore((s) => s.showHotspots3D);
  const hoveredHotspot = useStore((s) => s.hoveredHotspot);
  const previousShowRef = useRef<boolean>(false);
  const previousDataRef = useRef<string>('');
  const retryTimeoutRef = useRef<NodeJS.Timeout | null>(null);

  // Clear existing hotspot components
  const clearHotspotComponents = useCallback(async () => {
    if (!plugin) return;

    try {
      const state = plugin.state.data;
      const toRemove: string[] = [];

      state.cells.forEach((cell: any, ref: string) => {
        if (cell.obj?.label === HOTSPOT_COMPONENT_LABEL ||
            cell.obj?.label === HOTSPOT_BACKGROUND_LABEL) {
          toRemove.push(ref);
        }
      });

      for (const ref of toRemove) {
        await plugin.build().delete(ref).commit();
      }
    } catch (err) {
      // Ignore cleanup errors
    }
  }, [plugin]);

  // Create visualization for hotspot residues
  const highlightHotspots = useCallback(async (retryCount = 0) => {
    if (!plugin || !showHotspots3D) return;

    const structure = plugin.managers.structure.hierarchy.current.structures[0];
    if (!structure?.cell?.obj?.data) {
      // Structure not ready, retry after delay (max 5 retries)
      if (retryCount < 5 && hotspotsData?.residue_details && hotspotsData.residue_details.length > 0) {
        console.log(`[HotspotVisualization] Structure not ready, retrying in 500ms (attempt ${retryCount + 1})`);
        retryTimeoutRef.current = setTimeout(() => highlightHotspots(retryCount + 1), 500);
      }
      return;
    }

    // Clear existing highlights first
    await clearHotspotComponents();

    // Check if we have valid data
    if (!hotspotsData?.residue_details || hotspotsData.residue_details.length === 0) {
      console.log('[HotspotVisualization] No hotspot data to visualize');
      return;
    }

    try {
      // Get the state reference from the structure
      const structureRef = (structure as any).cell?.transform?.ref ?? (structure as any).ref;
      if (!structureRef) return;

      console.log(`[HotspotVisualization] Highlighting ${hotspotsData.residue_details.length} hotspot residues`);

      // Build combined expression for all hotspot residues
      const groups = hotspotsData.residue_details
        .map((r) => {
          const parsed = parseResidueId(r.residue);
          if (!parsed) return null;

          return MS.struct.generator.atomGroups({
            'chain-test': MS.core.rel.eq([
              MS.struct.atomProperty.macromolecular.auth_asym_id(),
              parsed.chain,
            ]),
            'residue-test': MS.core.rel.eq([
              MS.struct.atomProperty.macromolecular.auth_seq_id(),
              parsed.residue,
            ]),
          });
        })
        .filter((g): g is NonNullable<typeof g> => g !== null);

      if (groups.length === 0) {
        console.log('[HotspotVisualization] No valid residue IDs to visualize');
        return;
      }

      const expression = MS.struct.combinator.merge(groups);

      // 1. Create ball-and-stick for hotspot residues
      const component = await plugin.builders.structure.tryCreateComponentFromExpression(
        structureRef,
        expression,
        HOTSPOT_COMPONENT_LABEL
      );

      if (component) {
        // Use uniform color for now - could enhance to color by property
        await plugin.builders.structure.representation.addRepresentation(component, {
          type: 'ball-and-stick',
          color: 'uniform',
          colorParams: { value: HOTSPOT_COLOR },
          typeParams: { sizeFactor: 0.35 },
        });
      }

      // 2. Create transparent background cartoon for context
      const surroundingsExpr = surroundingsExpression(expression, 8.0);
      const backgroundExpr = exceptExpression(POLYMER_EXPRESSION, surroundingsExpr);

      const backgroundComponent = await plugin.builders.structure.tryCreateComponentFromExpression(
        structureRef,
        backgroundExpr,
        HOTSPOT_BACKGROUND_LABEL
      );

      if (backgroundComponent) {
        await plugin.builders.structure.representation.addRepresentation(backgroundComponent, {
          type: 'cartoon',
          color: 'chain-id',
          typeParams: { alpha: 0.3 },
        });
      }

      // 3. Focus camera on the hotspot cluster center if available
      if (hotspotsData.cluster_center) {
        // Auto-zoom to show all hotspots
        plugin.managers.camera.reset();
      }
    } catch (error) {
      console.error('[HotspotVisualization] Failed to highlight hotspots:', error);
    }
  }, [plugin, showHotspots3D, hotspotsData, clearHotspotComponents]);

  // Update highlights when showHotspots3D or data changes
  useEffect(() => {
    // Clean up any pending retry
    if (retryTimeoutRef.current) {
      clearTimeout(retryTimeoutRef.current);
      retryTimeoutRef.current = null;
    }

    if (!plugin) return;

    // Create a stable string representation to check for changes
    const currentData = JSON.stringify(hotspotsData?.hotspots?.sort() ?? []);
    const dataChanged = currentData !== previousDataRef.current;
    const showChanged = showHotspots3D !== previousShowRef.current;

    // Skip if nothing changed
    if (!dataChanged && !showChanged) return;

    previousDataRef.current = currentData;
    previousShowRef.current = showHotspots3D;

    if (showHotspots3D) {
      highlightHotspots();
    } else {
      // Clear visualization when toggled off
      clearHotspotComponents();
    }

    return () => {
      if (retryTimeoutRef.current) {
        clearTimeout(retryTimeoutRef.current);
      }
    };
  }, [plugin, showHotspots3D, hotspotsData, highlightHotspots, clearHotspotComponents]);

  // Handle hover highlighting - when user hovers over a hotspot in a panel
  useEffect(() => {
    if (!plugin) return;

    const structure = plugin.managers.structure.hierarchy.current.structures[0];
    if (!structure?.cell?.obj?.data) return;

    if (!hoveredHotspot) {
      // Clear highlight when nothing is hovered
      try {
        plugin.managers.interactivity.lociHighlights.clearHighlights();
      } catch (err) {
        // Ignore errors when clearing highlights
      }
      return;
    }

    try {
      // Parse the residue ID
      const parsed = parseResidueId(hoveredHotspot.residue);
      if (!parsed) return;

      // Build expression for the hovered residue
      const expression = MS.struct.generator.atomGroups({
        'chain-test': MS.core.rel.eq([
          MS.struct.atomProperty.macromolecular.auth_asym_id(),
          parsed.chain,
        ]),
        'residue-test': MS.core.rel.eq([
          MS.struct.atomProperty.macromolecular.auth_seq_id(),
          parsed.residue,
        ]),
      });

      // Get selection from structure
      const structureData = structure.cell.obj.data;
      const selection = Script.getStructureSelection(expression, structureData);

      if (!StructureSelection.isEmpty(selection)) {
        // Convert to loci and highlight
        const loci = StructureSelection.toLociWithSourceUnits(selection);
        plugin.managers.interactivity.lociHighlights.highlightOnly({ loci });

        // Focus camera on the residue for better visibility
        plugin.managers.camera.focusLoci(loci, { durationMs: 200 });
      }
    } catch (err) {
      console.warn('[HotspotVisualization] Hover highlight error:', err);
    }
  }, [plugin, hoveredHotspot]);

  return { hotspotsData, showHotspots3D, hoveredHotspot };
}
