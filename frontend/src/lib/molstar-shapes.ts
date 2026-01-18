/**
 * Molstar Shape Utilities
 * Efficient batched rendering for custom 3D shapes
 */

import type { PluginUIContext } from 'molstar/lib/mol-plugin-ui/context';
import type { InteractionLine } from './interactionGeometry';
import { Color } from 'molstar/lib/mol-util/color';

/**
 * Render interaction lines as batched cylinders
 * Creates a single mesh for all lines instead of individual shapes
 */
export async function renderInteractionLinesBatched(
  plugin: PluginUIContext,
  lines: InteractionLine[],
  options: {
    radius?: number;
    alpha?: number;
    label?: string;
  } = {}
): Promise<string | null> {
  if (!plugin || lines.length === 0) return null;

  const { radius = 0.08, alpha = 1.0, label = 'interaction-lines' } = options;

  try {
    // Clear existing interaction lines
    const state = plugin.state.data;
    const toRemove: string[] = [];
    state.cells.forEach((cell: any, ref: string) => {
      if (cell.obj?.label === label || cell.obj?.label?.startsWith('interaction-line')) {
        toRemove.push(ref);
      }
    });
    for (const ref of toRemove) {
      await plugin.build().delete(ref).commit();
    }

    // Import Molstar geometry modules
    const [
      { Shape },
      { MeshBuilder },
      { addCylinder },
      { Vec3 },
      { StateTransforms },
    ] = await Promise.all([
      import('molstar/lib/mol-model/shape'),
      import('molstar/lib/mol-geo/geometry/mesh/mesh-builder'),
      import('molstar/lib/mol-geo/geometry/mesh/builder/cylinder'),
      import('molstar/lib/mol-math/linear-algebra'),
      import('molstar/lib/mol-plugin-state/transforms'),
    ]);

    // Group lines by color for efficient rendering
    const linesByColor = new Map<number, InteractionLine[]>();
    for (const line of lines) {
      const existing = linesByColor.get(line.color) || [];
      existing.push(line);
      linesByColor.set(line.color, existing);
    }

    let totalRef: string | null = null;

    // Create one mesh per color group
    for (const [color, colorLines] of linesByColor) {
      // Estimate buffer sizes
      const vertexCount = colorLines.length * 64;  // ~64 vertices per cylinder
      const indexCount = colorLines.length * 128;  // ~128 indices per cylinder

      const builderState = MeshBuilder.createState(vertexCount, indexCount);

      for (const line of colorLines) {
        const startVec = Vec3.create(line.start[0], line.start[1], line.start[2]);
        const endVec = Vec3.create(line.end[0], line.end[1], line.end[2]);

        addCylinder(builderState, startVec, endVec, radius, {
          radiusTop: radius,
          radiusBottom: radius,
        });
      }

      const mesh = MeshBuilder.getMesh(builderState);

      const shape = Shape.create(
        `${label}-${color.toString(16)}`,
        {},
        mesh,
        () => Color(color),
        () => 1,
        () => `Interactions (${colorLines.length})`
      );

      // Note: Molstar's rawData accepts shapes at runtime despite TypeScript types
      // This pattern is used throughout the codebase (see ProteinViewerClient.tsx)
      const shapeData = await plugin.builders.data.rawData(
        { data: shape as any, label: `${label}-${color.toString(16)}` },
        { state: { isGhost: true } }
      );

      const result = await plugin.build()
        .to(shapeData)
        .apply(StateTransforms.Representation.ShapeRepresentation3D as any, { alpha })
        .commit();

      if (!totalRef) totalRef = result.ref;
    }

    console.log(`[MolstarShapes] Rendered ${lines.length} interaction lines in ${linesByColor.size} batches`);
    return totalRef;
  } catch (err) {
    console.error('[MolstarShapes] Failed to render interaction lines:', err);
    return null;
  }
}

/**
 * Clear all interaction line shapes
 */
export async function clearInteractionLines(
  plugin: PluginUIContext,
  label = 'interaction-lines'
): Promise<void> {
  if (!plugin) return;

  const state = plugin.state.data;
  const toRemove: string[] = [];

  state.cells.forEach((cell: any, ref: string) => {
    if (cell.obj?.label?.includes(label) || cell.obj?.label?.startsWith('interaction-line')) {
      toRemove.push(ref);
    }
  });

  for (const ref of toRemove) {
    await plugin.build().delete(ref).commit();
  }
}
