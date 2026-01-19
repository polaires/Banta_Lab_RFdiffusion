'use client';

import { useEffect, useRef, useCallback } from 'react';
import { useStore, type CatalyticSuggestion } from '@/lib/store';
import { Color } from 'molstar/lib/mol-util/color';
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder';
import type { PluginUIContext } from 'molstar/lib/mol-plugin-ui/context';

// Colors for suggestion highlighting
const MCSA_COLOR = Color(0x2563eb);  // Blue - for curated M-CSA catalytic residues
const LOCAL_COLOR = Color(0xf97316); // Orange - for local binding pocket detection

/**
 * Hook to highlight catalytic suggestions in Molstar viewer.
 * Call this from ProteinViewerClient or a parent component.
 */
export function useCatalyticVisualization(plugin: PluginUIContext | null) {
  const catalyticSuggestions = useStore((s) => s.catalyticSuggestions);
  const previousSuggestionsRef = useRef<CatalyticSuggestion[]>([]);
  const retryTimeoutRef = useRef<NodeJS.Timeout | null>(null);

  const attemptHighlight = useCallback(async (retryCount = 0) => {
    if (!plugin || catalyticSuggestions.length === 0) return;

    const structure = plugin.managers.structure.hierarchy.current.structures[0];
    if (!structure?.cell?.obj?.data) {
      // Structure not ready, retry after delay (max 5 retries)
      if (retryCount < 5) {
        console.log(`[CatalyticVisualization] Structure not ready, retrying in 500ms (attempt ${retryCount + 1})`);
        retryTimeoutRef.current = setTimeout(() => attemptHighlight(retryCount + 1), 500);
      }
      return;
    }

    await highlightSuggestions(plugin, catalyticSuggestions);
  }, [plugin, catalyticSuggestions]);

  useEffect(() => {
    // Clean up any pending retry
    if (retryTimeoutRef.current) {
      clearTimeout(retryTimeoutRef.current);
      retryTimeoutRef.current = null;
    }

    if (!plugin) return;

    // Skip if suggestions haven't changed
    const current = JSON.stringify(catalyticSuggestions);
    const previous = JSON.stringify(previousSuggestionsRef.current);
    if (current === previous) return;

    previousSuggestionsRef.current = catalyticSuggestions;

    attemptHighlight();

    return () => {
      if (retryTimeoutRef.current) {
        clearTimeout(retryTimeoutRef.current);
      }
    };
  }, [plugin, catalyticSuggestions, attemptHighlight]);

  return { catalyticSuggestions };
}

async function highlightSuggestions(
  plugin: PluginUIContext,
  suggestions: CatalyticSuggestion[]
) {
  if (suggestions.length === 0) {
    return;
  }

  try {
    const structure = plugin.managers.structure.hierarchy.current.structures[0];
    if (!structure?.cell?.obj?.data) {
      return;
    }

    // Group suggestions by source for batch processing
    const mcsaResidues = suggestions.filter((s) => s.source === 'mcsa');
    const localResidues = suggestions.filter((s) => s.source === 'local');

    // Get the state reference from the structure
    const structureRef = (structure as any).cell?.transform?.ref ?? (structure as any).ref;
    if (!structureRef) {
      return;
    }

    console.log(`[CatalyticVisualization] Creating highlights: ${mcsaResidues.length} mcsa, ${localResidues.length} local`);

    // Create highlight components for each group
    if (mcsaResidues.length > 0) {
      await createHighlightGroup(plugin, structureRef, mcsaResidues, MCSA_COLOR, 'mcsa-suggestions');
    }

    if (localResidues.length > 0) {
      await createHighlightGroup(plugin, structureRef, localResidues, LOCAL_COLOR, 'local-suggestions');
    }
  } catch (error) {
    console.error('[CatalyticVisualization] Failed to highlight:', error);
  }
}

async function createHighlightGroup(
  plugin: PluginUIContext,
  structureRef: any,
  residues: CatalyticSuggestion[],
  color: Color,
  label: string
) {
  // Build combined expression for all residues in this group
  const groups = residues.map((r) =>
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

  try {
    const component = await plugin.builders.structure.tryCreateComponentFromExpression(
      structureRef,
      expression,
      label
    );

    if (component) {
      await plugin.builders.structure.representation.addRepresentation(component, {
        type: 'ball-and-stick',
        color: 'uniform',
        colorParams: { value: color },
        typeParams: { sizeFactor: 0.3 },
      });
    }
  } catch (error) {
    console.error(`[CatalyticVisualization] Failed to create ${label}:`, error);
  }
}
