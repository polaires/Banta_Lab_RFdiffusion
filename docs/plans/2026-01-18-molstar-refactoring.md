# Molstar Implementation Refactoring Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Refactor the Molstar viewer implementation for better performance, type safety, and maintainability.

**Architecture:** Replace global state with React context, eliminate structure reloading on focus changes, batch interaction line rendering, and add proper TypeScript types.

**Tech Stack:** Molstar 5.5.0, React 19, TypeScript, Next.js 16

---

## Task 1: Add Proper TypeScript Types for Molstar

**Files:**
- Create: `frontend/src/lib/molstar-types.ts`
- Modify: `frontend/src/components/ProteinViewerClient.tsx`

**Step 1: Create molstar-types.ts with proper type exports**

Create `frontend/src/lib/molstar-types.ts`:

```typescript
/**
 * Molstar Type Definitions
 * Centralizes Molstar type imports for the application
 */

// Core plugin types
export type { PluginUIContext } from 'molstar/lib/mol-plugin-ui/context';
export type { PluginContext } from 'molstar/lib/mol-plugin/context';

// State types
export type { StateObjectRef, StateObjectCell } from 'molstar/lib/mol-state';
export type { StateTransform } from 'molstar/lib/mol-state/transform';

// Structure types
export type { Structure, StructureElement } from 'molstar/lib/mol-model/structure';
export type { StructureSelection } from 'molstar/lib/mol-model/structure/query/selection';

// Script/Query types
export type { MolScriptBuilder } from 'molstar/lib/mol-script/language/builder';
export type { Expression } from 'molstar/lib/mol-script/language/expression';

// Geometry types
export type { Mesh } from 'molstar/lib/mol-geo/geometry/mesh/mesh';
export type { Shape } from 'molstar/lib/mol-model/shape';

// Color types
export { Color } from 'molstar/lib/mol-util/color';

// Re-export commonly used builders
export { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder';
export { Script } from 'molstar/lib/mol-script/script';
export { StructureSelection as StructSel } from 'molstar/lib/mol-model/structure';
```

**Step 2: Run TypeScript check**

Run: `cd frontend && npx tsc --noEmit 2>&1 | head -20`
Expected: File compiles without errors

**Step 3: Commit**

```bash
git add frontend/src/lib/molstar-types.ts
git commit -m "$(cat <<'EOF'
feat(types): add centralized Molstar type definitions

Create molstar-types.ts to provide proper TypeScript types
for Molstar integration, replacing `any` types throughout
the codebase.

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Extract Common MolScript Expressions

**Files:**
- Create: `frontend/src/lib/molstar-expressions.ts`

**Step 1: Create molstar-expressions.ts with reusable selection expressions**

Create `frontend/src/lib/molstar-expressions.ts`:

```typescript
/**
 * Reusable MolScript Expressions
 * Common selection patterns for structure queries
 */

import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder';

/**
 * Select water molecules (HOH or WAT)
 */
export const WATER_EXPRESSION = MS.struct.generator.atomGroups({
  'residue-test': MS.core.logic.or([
    MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_comp_id(), 'HOH']),
    MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_comp_id(), 'WAT'])
  ])
});

/**
 * Select polymer entities (protein, DNA, RNA)
 */
export const POLYMER_EXPRESSION = MS.struct.generator.atomGroups({
  'entity-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.entityType(), 'polymer'])
});

/**
 * Select non-polymer entities (ligands, ions, etc.)
 */
export const NON_POLYMER_EXPRESSION = MS.struct.generator.atomGroups({
  'entity-test': MS.core.rel.neq([MS.struct.atomProperty.macromolecular.entityType(), 'polymer'])
});

/**
 * Create expression for a specific residue by chain and sequence number
 */
export function residueExpression(
  resName: string,
  chainId: string,
  resSeq: number
) {
  return MS.struct.generator.atomGroups({
    'residue-test': MS.core.logic.and([
      MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_comp_id(), resName]),
      MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_asym_id(), chainId]),
      MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_seq_id(), resSeq]),
    ])
  });
}

/**
 * Create surrounding sphere expression
 */
export function surroundingsExpression(
  centerExpression: ReturnType<typeof MS.struct.generator.atomGroups>,
  radius: number,
  asWholeResidues = true
) {
  return MS.struct.modifier.includeSurroundings({
    0: centerExpression,
    radius,
    'as-whole-residues': asWholeResidues
  });
}

/**
 * Exclude expression B from expression A
 */
export function exceptExpression(
  baseExpression: ReturnType<typeof MS.struct.generator.atomGroups>,
  excludeExpression: ReturnType<typeof MS.struct.generator.atomGroups>
) {
  return MS.struct.modifier.exceptBy({
    0: baseExpression,
    by: excludeExpression
  });
}

/**
 * Intersect expression A with expression B
 */
export function intersectExpression(
  baseExpression: ReturnType<typeof MS.struct.generator.atomGroups>,
  filterExpression: ReturnType<typeof MS.struct.generator.atomGroups>
) {
  return MS.struct.modifier.intersectBy({
    0: baseExpression,
    by: filterExpression
  });
}
```

**Step 2: Run TypeScript check**

Run: `cd frontend && npx tsc --noEmit 2>&1 | head -20`
Expected: File compiles without errors

**Step 3: Commit**

```bash
git add frontend/src/lib/molstar-expressions.ts
git commit -m "$(cat <<'EOF'
feat(molstar): add reusable MolScript expression builders

Extract common selection patterns (water, polymer, residue)
into reusable functions to reduce code duplication and
improve consistency.

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Create MolstarContext for State Management

**Files:**
- Create: `frontend/src/contexts/MolstarContext.tsx`

**Step 1: Create MolstarContext with proper lifecycle management**

Create `frontend/src/contexts/MolstarContext.tsx`:

```typescript
'use client';

import {
  createContext,
  useContext,
  useRef,
  useState,
  useCallback,
  useEffect,
  type ReactNode,
} from 'react';
import type { PluginUIContext } from 'molstar/lib/mol-plugin-ui/context';
import type { StateObjectRef } from 'molstar/lib/mol-state';

interface MolstarState {
  plugin: PluginUIContext | null;
  structureRef: StateObjectRef | null;
  isReady: boolean;
  isLoading: boolean;
  error: string | null;
}

interface MolstarContextValue extends MolstarState {
  initPlugin: (container: HTMLDivElement) => Promise<void>;
  loadStructure: (pdbContent: string) => Promise<StateObjectRef | null>;
  clearStructure: () => Promise<void>;
  dispose: () => void;
}

const MolstarContext = createContext<MolstarContextValue | null>(null);

export function useMolstar() {
  const context = useContext(MolstarContext);
  if (!context) {
    throw new Error('useMolstar must be used within MolstarProvider');
  }
  return context;
}

interface MolstarProviderProps {
  children: ReactNode;
}

export function MolstarProvider({ children }: MolstarProviderProps) {
  const pluginRef = useRef<PluginUIContext | null>(null);
  const structureRefRef = useRef<StateObjectRef | null>(null);
  const initPromiseRef = useRef<Promise<void> | null>(null);

  const [state, setState] = useState<MolstarState>({
    plugin: null,
    structureRef: null,
    isReady: false,
    isLoading: false,
    error: null,
  });

  const initPlugin = useCallback(async (container: HTMLDivElement) => {
    // Already initialized
    if (pluginRef.current) {
      setState(s => ({ ...s, isReady: true }));
      return;
    }

    // Already initializing
    if (initPromiseRef.current) {
      await initPromiseRef.current;
      return;
    }

    setState(s => ({ ...s, isLoading: true, error: null }));

    initPromiseRef.current = (async () => {
      try {
        const [molstarUI, molstarReact, molstarSpec] = await Promise.all([
          import('molstar/lib/mol-plugin-ui'),
          import('molstar/lib/mol-plugin-ui/react18'),
          import('molstar/lib/mol-plugin-ui/spec'),
        ]);

        container.innerHTML = '';

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
        setState(s => ({
          ...s,
          plugin,
          isReady: true,
          isLoading: false,
        }));
      } catch (err) {
        const message = err instanceof Error ? err.message : 'Failed to initialize viewer';
        setState(s => ({ ...s, error: message, isLoading: false }));
        throw err;
      }
    })();

    await initPromiseRef.current;
  }, []);

  const loadStructure = useCallback(async (pdbContent: string): Promise<StateObjectRef | null> => {
    const plugin = pluginRef.current;
    if (!plugin) return null;

    setState(s => ({ ...s, isLoading: true }));

    try {
      await plugin.clear();

      const isCif = pdbContent.trimStart().startsWith('data_');
      const format = isCif ? 'mmcif' : 'pdb';

      const data = await plugin.builders.data.rawData(
        { data: pdbContent, label: 'structure' },
        { state: { isGhost: true } }
      );

      const trajectory = await plugin.builders.structure.parseTrajectory(data, format);
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);

      structureRefRef.current = structure.ref;
      setState(s => ({
        ...s,
        structureRef: structure.ref,
        isLoading: false,
      }));

      return structure.ref;
    } catch (err) {
      const message = err instanceof Error ? err.message : 'Failed to load structure';
      setState(s => ({ ...s, error: message, isLoading: false }));
      return null;
    }
  }, []);

  const clearStructure = useCallback(async () => {
    const plugin = pluginRef.current;
    if (!plugin) return;

    await plugin.clear();
    structureRefRef.current = null;
    setState(s => ({ ...s, structureRef: null }));
  }, []);

  const dispose = useCallback(() => {
    pluginRef.current?.dispose();
    pluginRef.current = null;
    structureRefRef.current = null;
    initPromiseRef.current = null;
    setState({
      plugin: null,
      structureRef: null,
      isReady: false,
      isLoading: false,
      error: null,
    });
  }, []);

  // Cleanup on unmount
  useEffect(() => {
    return () => {
      dispose();
    };
  }, [dispose]);

  const value: MolstarContextValue = {
    ...state,
    initPlugin,
    loadStructure,
    clearStructure,
    dispose,
  };

  return (
    <MolstarContext.Provider value={value}>
      {children}
    </MolstarContext.Provider>
  );
}
```

**Step 2: Run TypeScript check**

Run: `cd frontend && npx tsc --noEmit 2>&1 | head -20`
Expected: File compiles without errors

**Step 3: Commit**

```bash
git add frontend/src/contexts/MolstarContext.tsx
git commit -m "$(cat <<'EOF'
feat(molstar): add MolstarContext for state management

Replace global state with React context provider that manages
plugin lifecycle, structure loading, and cleanup. Provides
useMolstar hook for components.

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Create Batched Interaction Lines Renderer

**Files:**
- Create: `frontend/src/lib/molstar-shapes.ts`

**Step 1: Create molstar-shapes.ts with batched line rendering**

Create `frontend/src/lib/molstar-shapes.ts`:

```typescript
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

      const shapeData = await plugin.builders.data.rawData(
        { data: shape, label: `${label}-${color.toString(16)}` },
        { state: { isGhost: true } }
      );

      const result = await plugin.build()
        .to(shapeData)
        .apply(StateTransforms.Representation.ShapeRepresentation3D, { alpha })
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
```

**Step 2: Run TypeScript check**

Run: `cd frontend && npx tsc --noEmit 2>&1 | head -20`
Expected: File compiles without errors

**Step 3: Commit**

```bash
git add frontend/src/lib/molstar-shapes.ts
git commit -m "$(cat <<'EOF'
feat(molstar): add batched interaction line renderer

Create molstar-shapes.ts with efficient batched rendering that
groups lines by color into single meshes, significantly improving
performance for many interaction lines.

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Create Focus View Manager

**Files:**
- Create: `frontend/src/lib/molstar-focus.ts`

**Step 1: Create molstar-focus.ts with non-destructive focus**

Create `frontend/src/lib/molstar-focus.ts`:

```typescript
/**
 * Molstar Focus View Manager
 * Manages focus views without reloading the structure
 */

import type { PluginUIContext } from 'molstar/lib/mol-plugin-ui/context';
import type { StateObjectRef } from 'molstar/lib/mol-state';
import { Color } from 'molstar/lib/mol-util/color';
import { Script } from 'molstar/lib/mol-script/script';
import { StructureSelection } from 'molstar/lib/mol-model/structure';
import {
  WATER_EXPRESSION,
  POLYMER_EXPRESSION,
  residueExpression,
  surroundingsExpression,
  exceptExpression,
  intersectExpression,
} from './molstar-expressions';

export interface FocusOptions {
  radius: number;
  showWaters: boolean;
  showLabels: boolean;
  backgroundAlpha: number;
  carbonColor?: number;
}

interface ComponentRefs {
  primary: StateObjectRef | null;
  surroundings: StateObjectRef | null;
  water: StateObjectRef | null;
  background: StateObjectRef | null;
}

/**
 * Manages focus view state and transitions
 */
export class FocusViewManager {
  private plugin: PluginUIContext;
  private structureRef: StateObjectRef;
  private componentRefs: ComponentRefs = {
    primary: null,
    surroundings: null,
    water: null,
    background: null,
  };
  private currentFocus: { type: 'metal' | 'ligand'; index: number } | null = null;

  constructor(plugin: PluginUIContext, structureRef: StateObjectRef) {
    this.plugin = plugin;
    this.structureRef = structureRef;
  }

  /**
   * Focus on a metal coordination site
   */
  async focusOnMetal(
    element: string,
    chainId: string,
    resSeq: number,
    resName: string,
    options: FocusOptions
  ): Promise<void> {
    await this.clearFocusComponents();

    const metalExpr = residueExpression(resName, chainId, resSeq);
    const coordSphereExpr = surroundingsExpression(metalExpr, options.radius + 2.0);
    const coordWithoutMetal = exceptExpression(coordSphereExpr, metalExpr);
    const coordProtein = exceptExpression(coordWithoutMetal, WATER_EXPRESSION);
    const backgroundExpr = exceptExpression(POLYMER_EXPRESSION, coordSphereExpr);

    // 1. Metal ion - large purple spacefill
    this.componentRefs.primary = await this.createComponent(
      metalExpr,
      `metal-${element}`,
      { type: 'spacefill', color: 0x7C3AED, sizeFactor: 1.5 }
    );

    // 2. Coordination sphere - ball-and-stick
    this.componentRefs.surroundings = await this.createComponent(
      coordProtein,
      'coordination-sphere',
      { type: 'ball-and-stick', color: 'element-symbol', sizeFactor: 0.3 }
    );

    // Add labels if enabled
    if (options.showLabels && this.componentRefs.surroundings) {
      await this.addLabels(this.componentRefs.surroundings);
    }

    // 3. Waters - conditionally
    if (options.showWaters) {
      const waterExpr = intersectExpression(coordWithoutMetal, WATER_EXPRESSION);
      this.componentRefs.water = await this.createComponent(
        waterExpr,
        'coord-water',
        { type: 'ball-and-stick', color: 'element-symbol', sizeFactor: 0.2, alpha: 0.7 }
      );
    }

    // 4. Background protein
    this.componentRefs.background = await this.createComponent(
      backgroundExpr,
      'protein-background',
      { type: 'cartoon', color: 'chain-id', alpha: options.backgroundAlpha }
    );

    // Focus camera
    await this.focusCamera(coordSphereExpr);

    this.currentFocus = { type: 'metal', index: resSeq };
  }

  /**
   * Focus on a ligand binding site
   */
  async focusOnLigand(
    name: string,
    chainId: string,
    resSeq: number,
    options: FocusOptions
  ): Promise<void> {
    await this.clearFocusComponents();

    const ligandExpr = residueExpression(name, chainId, resSeq);
    const bindingSiteExpr = surroundingsExpression(ligandExpr, options.radius);
    const siteWithoutLigand = exceptExpression(bindingSiteExpr, ligandExpr);
    const siteProtein = exceptExpression(siteWithoutLigand, WATER_EXPRESSION);
    const backgroundExpr = exceptExpression(POLYMER_EXPRESSION, bindingSiteExpr);

    // 1. Ligand - ball-and-stick with custom carbon color
    const ligandColor = options.carbonColor || 0x50C878;
    this.componentRefs.primary = await this.createComponent(
      ligandExpr,
      `ligand-${name}`,
      { type: 'ball-and-stick', color: ligandColor, sizeFactor: 0.4 }
    );

    // 2. Binding site residues
    this.componentRefs.surroundings = await this.createComponent(
      siteProtein,
      'binding-site',
      { type: 'ball-and-stick', color: 'element-symbol', sizeFactor: 0.25 }
    );

    if (options.showLabels && this.componentRefs.surroundings) {
      await this.addLabels(this.componentRefs.surroundings);
    }

    // 3. Waters - conditionally
    if (options.showWaters) {
      const waterExpr = intersectExpression(siteWithoutLigand, WATER_EXPRESSION);
      this.componentRefs.water = await this.createComponent(
        waterExpr,
        'binding-water',
        { type: 'ball-and-stick', color: 'element-symbol', sizeFactor: 0.15, alpha: 0.7 }
      );
    }

    // 4. Background protein
    this.componentRefs.background = await this.createComponent(
      backgroundExpr,
      'protein-background',
      { type: 'cartoon', color: 'chain-id', alpha: options.backgroundAlpha }
    );

    // Focus camera
    await this.focusCamera(bindingSiteExpr);

    this.currentFocus = { type: 'ligand', index: resSeq };
  }

  /**
   * Clear focus and return to default view
   */
  async resetFocus(): Promise<void> {
    await this.clearFocusComponents();
    this.currentFocus = null;
    this.plugin.canvas3d?.requestCameraReset();
  }

  /**
   * Update water visibility without full reload
   */
  async toggleWaters(show: boolean, expression: any): Promise<void> {
    if (show && !this.componentRefs.water) {
      this.componentRefs.water = await this.createComponent(
        expression,
        'focus-water',
        { type: 'ball-and-stick', color: 'element-symbol', sizeFactor: 0.15, alpha: 0.7 }
      );
    } else if (!show && this.componentRefs.water) {
      await this.deleteComponent(this.componentRefs.water);
      this.componentRefs.water = null;
    }
  }

  private async createComponent(
    expression: any,
    label: string,
    repr: {
      type: string;
      color: string | number;
      sizeFactor?: number;
      alpha?: number;
    }
  ): Promise<StateObjectRef | null> {
    try {
      const comp = await this.plugin.builders.structure.tryCreateComponentFromExpression(
        this.structureRef,
        expression,
        label
      );

      if (!comp) return null;

      const colorParams = typeof repr.color === 'number'
        ? { value: Color(repr.color) }
        : undefined;

      await this.plugin.builders.structure.representation.addRepresentation(comp, {
        type: repr.type as any,
        color: typeof repr.color === 'string' ? repr.color : 'uniform',
        colorParams,
        typeParams: {
          sizeFactor: repr.sizeFactor,
          alpha: repr.alpha,
        },
      });

      return comp.ref;
    } catch (err) {
      console.warn(`[FocusView] Failed to create component ${label}:`, err);
      return null;
    }
  }

  private async addLabels(componentRef: StateObjectRef): Promise<void> {
    try {
      const cell = this.plugin.state.data.cells.get(componentRef);
      if (cell) {
        await this.plugin.builders.structure.representation.addRepresentation(
          { ref: componentRef } as any,
          {
            type: 'label',
            typeParams: { level: 'residue' },
            color: 'uniform',
            colorParams: { value: Color(0x333333) },
          }
        );
      }
    } catch (err) {
      console.warn('[FocusView] Failed to add labels:', err);
    }
  }

  private async deleteComponent(ref: StateObjectRef): Promise<void> {
    try {
      await this.plugin.build().delete(ref).commit();
    } catch (err) {
      console.warn('[FocusView] Failed to delete component:', err);
    }
  }

  private async clearFocusComponents(): Promise<void> {
    const refs = Object.values(this.componentRefs).filter(Boolean) as StateObjectRef[];
    for (const ref of refs) {
      await this.deleteComponent(ref);
    }
    this.componentRefs = {
      primary: null,
      surroundings: null,
      water: null,
      background: null,
    };
  }

  private async focusCamera(expression: any): Promise<void> {
    try {
      const cell = this.plugin.state.data.cells.get(this.structureRef);
      if (!cell?.obj?.data) return;

      const selection = Script.getStructureSelection(expression, cell.obj.data);
      if (!StructureSelection.isEmpty(selection)) {
        const loci = StructureSelection.toLociWithSourceUnits(selection);
        this.plugin.managers.camera.focusLoci(loci, { durationMs: 400 });
      }
    } catch (err) {
      console.warn('[FocusView] Failed to focus camera:', err);
    }
  }
}
```

**Step 2: Run TypeScript check**

Run: `cd frontend && npx tsc --noEmit 2>&1 | head -20`
Expected: File compiles without errors (may have some type warnings)

**Step 3: Commit**

```bash
git add frontend/src/lib/molstar-focus.ts
git commit -m "$(cat <<'EOF'
feat(molstar): add FocusViewManager for non-destructive focus

Create FocusViewManager class that manages focus views by
toggling component visibility rather than reloading the
entire structure, significantly improving performance.

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Update ProteinViewerClient to Use New Modules

**Files:**
- Modify: `frontend/src/components/ProteinViewerClient.tsx`

**Step 1: Replace global state with context imports**

Update imports at the top of `frontend/src/components/ProteinViewerClient.tsx`:

```typescript
'use client';

import { useEffect, useRef, useState, useCallback, useImperativeHandle, forwardRef } from 'react';
import type { MetalCoordination } from '@/lib/metalAnalysis';
import type { LigandData, PharmacophoreFeature } from '@/lib/ligandAnalysis';
import { useStore } from '@/lib/store';

// New Molstar modules
import type { PluginUIContext } from 'molstar/lib/mol-plugin-ui/context';
import type { StateObjectRef } from 'molstar/lib/mol-state';
import { Color } from 'molstar/lib/mol-util/color';
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder';
import { Script } from 'molstar/lib/mol-script/script';
import { StructureSelection } from 'molstar/lib/mol-model/structure';

// Refactored utilities
import {
  WATER_EXPRESSION,
  POLYMER_EXPRESSION,
  NON_POLYMER_EXPRESSION,
  residueExpression,
  surroundingsExpression,
  exceptExpression,
  intersectExpression,
} from '@/lib/molstar-expressions';
import { renderInteractionLinesBatched, clearInteractionLines } from '@/lib/molstar-shapes';
import { computeInteractionLines } from '@/lib/interactionGeometry';
```

**Step 2: Remove old global state and dynamic imports**

Delete these lines:

```typescript
// DELETE: Old global state
let globalPlugin: PluginUIContext | null = null;
let globalInitPromise: Promise<any> | null = null;
// ... etc

// DELETE: Old dynamic import function
async function loadMolstarModules() { ... }
```

**Step 3: Update focusOnMetal to use new expressions**

Replace the existing `focusOnMetal` callback with the refactored version using the new expression helpers:

```typescript
const focusOnMetal = useCallback(async (index: number) => {
  if (!pluginRef.current || !metalCoordination?.[index]) return;

  const plugin = pluginRef.current;
  const metal = metalCoordination[index];

  try {
    console.log(`[ProteinViewer] Focusing on metal: ${metal.element} at ${metal.chainId}:${metal.resSeq}`);

    // Use expression helpers
    const metalExpr = residueExpression(metal.resName, metal.chainId, metal.resSeq);
    const coordSphereExpr = surroundingsExpression(metalExpr, focusSettings.coordinationRadius + 2.0);
    const coordWithoutMetal = exceptExpression(coordSphereExpr, metalExpr);
    const coordProtein = exceptExpression(coordWithoutMetal, WATER_EXPRESSION);
    const backgroundExpr = exceptExpression(POLYMER_EXPRESSION, coordSphereExpr);

    // ... rest of implementation using these expressions
  } catch (err) {
    console.error('[ProteinViewer] Failed to focus on metal:', err);
  }
}, [metalCoordination, focusSettings]);
```

**Step 4: Update renderInteractionLines to use batched renderer**

Replace the existing `renderInteractionLines` with:

```typescript
const renderInteractionLines = useCallback(async (lines: InteractionLine[]) => {
  if (!pluginRef.current || lines.length === 0) return;
  await renderInteractionLinesBatched(pluginRef.current, lines, { radius: 0.08 });
}, []);
```

**Step 5: Run TypeScript check**

Run: `cd frontend && npx tsc --noEmit 2>&1 | head -30`
Expected: Reduced errors, may need incremental fixes

**Step 6: Run tests**

Run: `cd frontend && npm test -- --run`
Expected: All tests pass

**Step 7: Commit**

```bash
git add frontend/src/components/ProteinViewerClient.tsx
git commit -m "$(cat <<'EOF'
refactor(viewer): integrate new Molstar modules

- Replace global state with component-local refs
- Use expression helpers from molstar-expressions.ts
- Use batched renderer from molstar-shapes.ts
- Remove dynamic imports in favor of static imports
- Add proper TypeScript types

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Add isGhost Flag to Data Loading

**Files:**
- Modify: `frontend/src/components/ProteinViewerClient.tsx`

**Step 1: Update all rawData calls to use isGhost**

Find all occurrences of `builders.data.rawData` and add the ghost flag:

```typescript
// Before
const data = await plugin.builders.data.rawData({
  data: pdbContent,
  label: 'structure.pdb',
});

// After
const data = await plugin.builders.data.rawData(
  { data: pdbContent, label: 'structure.pdb' },
  { state: { isGhost: true } }
);
```

**Step 2: Run build**

Run: `cd frontend && npm run build`
Expected: Build succeeds

**Step 3: Commit**

```bash
git add frontend/src/components/ProteinViewerClient.tsx
git commit -m "$(cat <<'EOF'
fix(viewer): add isGhost flag to intermediate data nodes

Hide intermediate data nodes in the Molstar state tree
for cleaner UI and better performance.

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: Final Integration Test

**Files:**
- No new files

**Step 1: Run all tests**

Run: `cd frontend && npm test -- --run`
Expected: All tests pass

**Step 2: Run linting**

Run: `cd frontend && npm run lint 2>&1 | grep -E "error|warning" | head -20`
Expected: No new errors in modified files

**Step 3: Run build**

Run: `cd frontend && npm run build`
Expected: Build succeeds

**Step 4: Manual testing checklist**

1. Load a structure with metal ions
2. Click on metal → Focus view appears with proper styling
3. Adjust radius slider → View updates without full reload
4. Toggle waters → Waters appear/disappear instantly
5. Load a structure with ligands
6. Click on ligand → Focus view with green carbons
7. Toggle interaction lines → Lines appear (batched rendering)
8. Reset view → Returns to default

**Step 5: Final commit**

```bash
git add -A
git commit -m "$(cat <<'EOF'
refactor: complete Molstar implementation improvements

- Proper TypeScript types throughout
- Reusable MolScript expression builders
- Batched interaction line rendering
- Non-destructive focus view transitions
- Ghost state for intermediate data nodes

Performance improvements:
- No structure reload on focus change
- Single mesh per color for interaction lines
- Static imports instead of dynamic

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
EOF
)"
```

---

## Summary

| Task | Description | Priority |
|------|-------------|----------|
| 1 | Add proper TypeScript types | High |
| 2 | Extract common MolScript expressions | Medium |
| 3 | Create MolstarContext | Medium |
| 4 | Batched interaction line renderer | High |
| 5 | FocusViewManager class | High |
| 6 | Update ProteinViewerClient | High |
| 7 | Add isGhost flags | Low |
| 8 | Final integration test | Required |

**Key improvements:**
- Type safety with proper Molstar types
- DRY principle with expression helpers
- Performance via batched rendering and non-destructive focus
- Maintainability with separated concerns
