'use client';

import { useEffect, useRef, useCallback } from 'react';
import { useStore, type CatalyticSuggestion } from '@/lib/store';
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

// Colors for highlighting
const SELECTED_COLOR = Color(0xf97316);  // Orange - for selected catalytic residues
const HOVER_COLOR = Color(0xfbbf24);     // Yellow/amber - for hover highlight

// Component labels for cleanup
const SELECTED_COMPONENT_LABEL = 'selected-catalytic-residues';
const BACKGROUND_COMPONENT_LABEL = 'catalytic-background-cartoon';

/**
 * Hook to highlight catalytic residues in Molstar viewer.
 * - Highlights SELECTED residues (enzymeCatalyticResidues) with orange ball-and-stick
 * - Adds transparent background cartoon for context (like focus mode)
 * - Handles hover highlighting when user hovers over a suggestion in the panel
 * Call this from ProteinViewerClient or a parent component.
 */
export function useCatalyticVisualization(plugin: PluginUIContext | null) {
  const catalyticSuggestions = useStore((s) => s.catalyticSuggestions);
  const enzymeCatalyticResidues = useStore((s) => s.enzymeCatalyticResidues);
  const hoveredCatalyticSuggestion = useStore((s) => s.hoveredCatalyticSuggestion);
  const previousSelectedRef = useRef<string>('');
  const retryTimeoutRef = useRef<NodeJS.Timeout | null>(null);

  // Clear existing highlight components
  const clearHighlightComponents = useCallback(async () => {
    if (!plugin) return;

    try {
      const state = plugin.state.data;
      const toRemove: string[] = [];

      state.cells.forEach((cell: any, ref: string) => {
        if (cell.obj?.label === SELECTED_COMPONENT_LABEL ||
            cell.obj?.label === BACKGROUND_COMPONENT_LABEL) {
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

  // Create highlight for selected residues
  const highlightSelectedResidues = useCallback(async (retryCount = 0) => {
    if (!plugin) return;

    const structure = plugin.managers.structure.hierarchy.current.structures[0];
    if (!structure?.cell?.obj?.data) {
      // Structure not ready, retry after delay (max 5 retries)
      if (retryCount < 5 && enzymeCatalyticResidues.length > 0) {
        console.log(`[CatalyticVisualization] Structure not ready, retrying in 500ms (attempt ${retryCount + 1})`);
        retryTimeoutRef.current = setTimeout(() => highlightSelectedResidues(retryCount + 1), 500);
      }
      return;
    }

    // Clear existing highlights first
    await clearHighlightComponents();

    // If no selected residues, we're done
    if (enzymeCatalyticResidues.length === 0) {
      console.log('[CatalyticVisualization] No selected residues to highlight');
      return;
    }

    try {
      // Get the state reference from the structure
      const structureRef = (structure as any).cell?.transform?.ref ?? (structure as any).ref;
      if (!structureRef) return;

      console.log(`[CatalyticVisualization] Highlighting ${enzymeCatalyticResidues.length} selected residues`);

      // Build combined expression for all selected residues
      const groups = enzymeCatalyticResidues.map((r) =>
        MS.struct.generator.atomGroups({
          'chain-test': MS.core.rel.eq([
            MS.struct.atomProperty.macromolecular.auth_asym_id(),
            r.chain,
          ]),
          'residue-test': MS.core.rel.eq([
            MS.struct.atomProperty.macromolecular.auth_seq_id(),
            r.residue,
          ]),
        })
      );

      const expression = MS.struct.combinator.merge(groups);

      // 1. Create ball-and-stick for catalytic residues (orange)
      const component = await plugin.builders.structure.tryCreateComponentFromExpression(
        structureRef,
        expression,
        SELECTED_COMPONENT_LABEL
      );

      if (component) {
        await plugin.builders.structure.representation.addRepresentation(component, {
          type: 'ball-and-stick',
          color: 'uniform',
          colorParams: { value: SELECTED_COLOR },
          typeParams: { sizeFactor: 0.35 },
        });
      }

      // 2. Create transparent background cartoon for context (like focus mode)
      // Get surroundings of the selected residues and exclude them to show background
      const surroundingsExpr = surroundingsExpression(expression, 8.0);
      const backgroundExpr = exceptExpression(POLYMER_EXPRESSION, surroundingsExpr);

      const backgroundComponent = await plugin.builders.structure.tryCreateComponentFromExpression(
        structureRef,
        backgroundExpr,
        BACKGROUND_COMPONENT_LABEL
      );

      if (backgroundComponent) {
        await plugin.builders.structure.representation.addRepresentation(backgroundComponent, {
          type: 'cartoon',
          color: 'chain-id',
          typeParams: { alpha: 0.3 },
        });
      }
    } catch (error) {
      console.error('[CatalyticVisualization] Failed to highlight selected residues:', error);
    }
  }, [plugin, enzymeCatalyticResidues, clearHighlightComponents]);

  // Update highlights when selected residues change
  useEffect(() => {
    // Clean up any pending retry
    if (retryTimeoutRef.current) {
      clearTimeout(retryTimeoutRef.current);
      retryTimeoutRef.current = null;
    }

    if (!plugin) return;

    // Create a stable string representation to check for changes
    const currentSelected = JSON.stringify(
      enzymeCatalyticResidues.map((r) => `${r.chain}${r.residue}`).sort()
    );

    // Skip if selection hasn't changed
    if (currentSelected === previousSelectedRef.current) return;
    previousSelectedRef.current = currentSelected;

    highlightSelectedResidues();

    return () => {
      if (retryTimeoutRef.current) {
        clearTimeout(retryTimeoutRef.current);
      }
    };
  }, [plugin, enzymeCatalyticResidues, highlightSelectedResidues]);

  // Handle hover highlighting - when user hovers over a suggestion in the panel
  useEffect(() => {
    if (!plugin) return;

    const structure = plugin.managers.structure.hierarchy.current.structures[0];
    if (!structure?.cell?.obj?.data) return;

    if (!hoveredCatalyticSuggestion) {
      // Clear highlight when nothing is hovered
      try {
        plugin.managers.interactivity.lociHighlights.clearHighlights();
      } catch (err) {
        // Ignore errors when clearing highlights
      }
      return;
    }

    try {
      // Build expression for the hovered residue
      const expression = MS.struct.generator.atomGroups({
        'chain-test': MS.core.rel.eq([
          MS.struct.atomProperty.macromolecular.auth_asym_id(),
          hoveredCatalyticSuggestion.chain,
        ]),
        'residue-test': MS.core.rel.eq([
          MS.struct.atomProperty.macromolecular.auth_seq_id(),
          hoveredCatalyticSuggestion.residue,
        ]),
      });

      // Get selection from structure
      const structureData = structure.cell.obj.data;
      const selection = Script.getStructureSelection(expression, structureData);

      if (!StructureSelection.isEmpty(selection)) {
        // Convert to loci and highlight
        const loci = StructureSelection.toLociWithSourceUnits(selection);
        plugin.managers.interactivity.lociHighlights.highlightOnly({ loci });

        // Also focus camera on the residue for better visibility
        plugin.managers.camera.focusLoci(loci, { durationMs: 200 });
      }
    } catch (err) {
      console.warn('[CatalyticVisualization] Hover highlight error:', err);
    }
  }, [plugin, hoveredCatalyticSuggestion]);

  return { catalyticSuggestions, enzymeCatalyticResidues, hoveredCatalyticSuggestion };
}
