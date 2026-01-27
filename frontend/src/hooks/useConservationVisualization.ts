'use client';

import { useEffect, useRef, useCallback } from 'react';
import { useStore, type ConservationGrade } from '@/lib/store';
import { Color } from 'molstar/lib/mol-util/color';
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder';
import { Script } from 'molstar/lib/mol-script/script';
import { StructureSelection } from 'molstar/lib/mol-model/structure';
import type { PluginUIContext } from 'molstar/lib/mol-plugin-ui/context';
import { POLYMER_EXPRESSION } from '@/lib/molstar-expressions';

// ConSurf color scheme: 9-grade scale from conserved (burgundy) to variable (cyan)
const CONSURF_COLORS: Record<number, number> = {
  1: 0x10004D,  // Most conserved - dark burgundy
  2: 0x4C006B,
  3: 0x8B0057,
  4: 0xD2004B,
  5: 0xFFFFFF,  // Average - white
  6: 0x98D6EA,
  7: 0x61D2DA,
  8: 0x35CBCB,
  9: 0x00FFFF,  // Most variable - cyan
};

// Component labels for cleanup
const GRADE_LABEL_PREFIX = 'conservation-grade-';
const CONSERVATION_BG_LABEL = 'conservation-background';

function gradeLabel(grade: number): string {
  return `${GRADE_LABEL_PREFIX}${grade}`;
}

/**
 * Hook to visualize conservation data using ConSurf-style coloring in Molstar.
 * Groups residues by grade (1-9), creates a cartoon component per grade
 * with ConSurf colors (burgundy=conserved -> white -> cyan=variable).
 *
 * Strategy: Create a transparent background cartoon for context, then
 * overlay opaque grade-colored cartoon segments on top.
 */
export function useConservationVisualization(plugin: PluginUIContext | null) {
  const conservationData = useStore((s) => s.conservationData);
  const showConservation3D = useStore((s) => s.showConservation3D);
  const hoveredResidue = useStore((s) => s.hoveredConservationResidue);
  const previousShowRef = useRef<boolean>(false);
  const previousDataRef = useRef<string>('');
  const retryTimeoutRef = useRef<NodeJS.Timeout | null>(null);

  // Clear existing conservation components
  const clearConservationComponents = useCallback(async () => {
    if (!plugin) return;

    try {
      const state = plugin.state.data;
      const toRemove: string[] = [];

      state.cells.forEach((cell: any, ref: string) => {
        if (cell.obj?.label?.startsWith(GRADE_LABEL_PREFIX) ||
            cell.obj?.label === CONSERVATION_BG_LABEL) {
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

  // Create ConSurf-style visualization
  const visualizeConservation = useCallback(async (retryCount = 0) => {
    if (!plugin || !showConservation3D) return;

    const structure = plugin.managers.structure.hierarchy.current.structures[0];
    if (!structure?.cell?.obj?.data) {
      if (retryCount < 5 && conservationData?.grades && conservationData.grades.length > 0) {
        console.log(`[ConservationVis] Structure not ready, retrying (attempt ${retryCount + 1})`);
        retryTimeoutRef.current = setTimeout(() => visualizeConservation(retryCount + 1), 500);
      }
      return;
    }

    // Clear existing conservation components
    await clearConservationComponents();

    if (!conservationData?.grades || conservationData.grades.length === 0) {
      console.log('[ConservationVis] No conservation data to visualize');
      return;
    }

    try {
      const structureRef = (structure as any).cell?.transform?.ref ?? (structure as any).ref;
      if (!structureRef) return;

      console.log(`[ConservationVis] Coloring ${conservationData.grades.length} residues by conservation grade`);

      // Group residues by grade
      const byGrade = new Map<number, ConservationGrade[]>();
      for (const g of conservationData.grades) {
        if (!byGrade.has(g.grade)) byGrade.set(g.grade, []);
        byGrade.get(g.grade)!.push(g);
      }

      // 1. Create a semi-transparent background cartoon for structural context
      const bgComponent = await plugin.builders.structure.tryCreateComponentFromExpression(
        structureRef,
        POLYMER_EXPRESSION,
        CONSERVATION_BG_LABEL
      );

      if (bgComponent) {
        await plugin.builders.structure.representation.addRepresentation(bgComponent, {
          type: 'cartoon',
          color: 'uniform',
          colorParams: { value: Color(0xDDDDDD) },
          typeParams: { alpha: 0.15 },
        });
      }

      // 2. Create opaque cartoon component per grade with ConSurf color
      for (const [grade, residues] of byGrade) {
        const colorHex = CONSURF_COLORS[grade] ?? 0xCCCCCC;

        const groups = residues.map((r) =>
          MS.struct.generator.atomGroups({
            'chain-test': MS.core.rel.eq([
              MS.struct.atomProperty.macromolecular.auth_asym_id(),
              'A',
            ]),
            'residue-test': MS.core.rel.eq([
              MS.struct.atomProperty.macromolecular.auth_seq_id(),
              r.position,
            ]),
          })
        );

        if (groups.length === 0) continue;

        const expression = groups.length === 1 ? groups[0] : MS.struct.combinator.merge(groups);

        const component = await plugin.builders.structure.tryCreateComponentFromExpression(
          structureRef,
          expression,
          gradeLabel(grade)
        );

        if (component) {
          await plugin.builders.structure.representation.addRepresentation(component, {
            type: 'cartoon',
            color: 'uniform',
            colorParams: { value: Color(colorHex) },
          });
        }
      }

      console.log(`[ConservationVis] Created ${byGrade.size} grade groups`);
    } catch (error) {
      console.error('[ConservationVis] Failed to visualize conservation:', error);
    }
  }, [plugin, showConservation3D, conservationData, clearConservationComponents]);

  // Update visualization when toggle or data changes
  useEffect(() => {
    if (retryTimeoutRef.current) {
      clearTimeout(retryTimeoutRef.current);
      retryTimeoutRef.current = null;
    }

    if (!plugin) return;

    const currentData = JSON.stringify(conservationData?.grades?.map(g => g.position).sort() ?? []);
    const dataChanged = currentData !== previousDataRef.current;
    const showChanged = showConservation3D !== previousShowRef.current;

    if (!dataChanged && !showChanged) return;

    previousDataRef.current = currentData;
    previousShowRef.current = showConservation3D;

    if (showConservation3D) {
      visualizeConservation();
    } else {
      clearConservationComponents();
    }

    return () => {
      if (retryTimeoutRef.current) {
        clearTimeout(retryTimeoutRef.current);
      }
    };
  }, [plugin, showConservation3D, conservationData, visualizeConservation, clearConservationComponents]);

  // Handle hover highlighting
  useEffect(() => {
    if (!plugin) return;

    const structure = plugin.managers.structure.hierarchy.current.structures[0];
    if (!structure?.cell?.obj?.data) return;

    if (!hoveredResidue) {
      try {
        plugin.managers.interactivity.lociHighlights.clearHighlights();
      } catch {
        // Ignore
      }
      return;
    }

    try {
      const expression = MS.struct.generator.atomGroups({
        'chain-test': MS.core.rel.eq([
          MS.struct.atomProperty.macromolecular.auth_asym_id(),
          'A',
        ]),
        'residue-test': MS.core.rel.eq([
          MS.struct.atomProperty.macromolecular.auth_seq_id(),
          hoveredResidue.position,
        ]),
      });

      const structureData = structure.cell.obj.data;
      const selection = Script.getStructureSelection(expression, structureData);

      if (!StructureSelection.isEmpty(selection)) {
        const loci = StructureSelection.toLociWithSourceUnits(selection);
        plugin.managers.interactivity.lociHighlights.highlightOnly({ loci });
        plugin.managers.camera.focusLoci(loci, { durationMs: 200 });
      }
    } catch (err) {
      console.warn('[ConservationVis] Hover highlight error:', err);
    }
  }, [plugin, hoveredResidue]);

  return { conservationData, showConservation3D };
}
