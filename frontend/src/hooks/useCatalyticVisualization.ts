'use client';

import { useEffect, useRef } from 'react';
import { useStore, type CatalyticSuggestion } from '@/lib/store';
import { Color } from 'molstar/lib/mol-util/color';
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder';
import type { PluginUIContext } from 'molstar/lib/mol-plugin-ui/context';

// Colors for suggestion highlighting
const MCSA_COLOR = Color(0x2563eb);  // Blue
const P2RANK_COLOR = Color(0xf97316); // Orange

/**
 * Hook to highlight catalytic suggestions in Molstar viewer.
 * Call this from ProteinViewerClient or a parent component.
 */
export function useCatalyticVisualization(plugin: PluginUIContext | null) {
  const catalyticSuggestions = useStore((s) => s.catalyticSuggestions);
  const previousSuggestionsRef = useRef<CatalyticSuggestion[]>([]);

  useEffect(() => {
    if (!plugin) return;

    // Skip if suggestions haven't changed
    const current = JSON.stringify(catalyticSuggestions);
    const previous = JSON.stringify(previousSuggestionsRef.current);
    if (current === previous) return;

    previousSuggestionsRef.current = catalyticSuggestions;

    highlightSuggestions(plugin, catalyticSuggestions);
  }, [plugin, catalyticSuggestions]);

  return { catalyticSuggestions };
}

async function highlightSuggestions(
  plugin: PluginUIContext,
  suggestions: CatalyticSuggestion[]
) {
  if (suggestions.length === 0) {
    // Clear any existing suggestion highlights
    // This will be handled by clearing custom components
    return;
  }

  try {
    const structure = plugin.managers.structure.hierarchy.current.structures[0];
    if (!structure?.cell?.obj?.data) return;

    // Group suggestions by source for batch processing
    const mcsaResidues = suggestions.filter((s) => s.source === 'mcsa');
    const p2rankResidues = suggestions.filter((s) => s.source === 'p2rank');

    // Get the state reference from the structure
    // StructureRef from hierarchy uses cell.transform.ref, not .ref directly
    const structureRef = (structure as any).cell?.transform?.ref ?? (structure as any).ref;
    if (!structureRef) return;

    // Create highlight components for each group
    if (mcsaResidues.length > 0) {
      await createHighlightGroup(plugin, structureRef, mcsaResidues, MCSA_COLOR, 'mcsa-suggestions');
    }

    if (p2rankResidues.length > 0) {
      await createHighlightGroup(plugin, structureRef, p2rankResidues, P2RANK_COLOR, 'p2rank-suggestions');
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
