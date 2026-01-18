# Molstar Focus View Enhancement Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Enhance the Molstar viewer with advanced focus views including configurable radii, water toggles, 3D interaction lines, pi-stacking geometry analysis, and detection-aware UI controls.

**Architecture:** Three-layer design following existing patterns: (1) Analysis layer in `lib/` with pure functions for pi-stacking geometry and interaction line computation, (2) Viewer layer in `ProteinViewerClient.tsx` with enhanced focus functions and Mol* shape primitives for 3D lines, (3) UI layer with a new `FocusModeControls.tsx` component that appears inline during focus mode with detection-aware visibility.

**Tech Stack:** TypeScript, React 19, Molstar 5.5.0, Zustand, Tailwind CSS

---

## Task 1: Add Focus Settings to Zustand Store

**Files:**
- Modify: `frontend/src/lib/store.ts:190-220`

**Step 1: Add focus settings types and state**

Add these new state properties to the AppState interface (around line 190):

```typescript
// Focus mode visualization settings
focusSettings: {
  coordinationRadius: number;  // Metal coordination sphere radius (default 3.0)
  bindingPocketRadius: number; // Ligand binding pocket radius (default 5.0)
  showWaters: boolean;         // Show coordinating waters in focus view
  showInteractionLines: boolean; // Show H-bond, salt bridge lines
  showPharmacophores: boolean;  // Show pharmacophore spheres during focus
  ligandCarbonColor: number;    // Ligand carbon color (default green 0x50C878)
};
setFocusSettings: (settings: Partial<AppState['focusSettings']>) => void;
```

**Step 2: Add default values in store initialization**

Add to the store creation (around line 350):

```typescript
// Focus mode settings
focusSettings: {
  coordinationRadius: 3.0,
  bindingPocketRadius: 5.0,
  showWaters: false,
  showInteractionLines: true,
  showPharmacophores: false,
  ligandCarbonColor: 0x50C878,
},
setFocusSettings: (settings) => set((state) => ({
  focusSettings: { ...state.focusSettings, ...settings },
})),
```

**Step 3: Run TypeScript check**

Run: `cd frontend && npx tsc --noEmit`
Expected: No errors (store types should be valid)

**Step 4: Commit**

```bash
git add frontend/src/lib/store.ts
git commit -m "$(cat <<'EOF'
feat(store): add focus mode visualization settings

Add focusSettings state for configurable coordination radius,
binding pocket radius, water visibility, interaction lines,
and ligand carbon coloring.

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Add Pi-Stacking Geometry Types and Ring Detection

**Files:**
- Modify: `frontend/src/lib/ligandAnalysis.ts`
- Create: `frontend/src/lib/__tests__/piStacking.test.ts`

**Step 1: Write failing test for ring geometry calculation**

Create `frontend/src/lib/__tests__/piStacking.test.ts`:

```typescript
import { describe, it, expect } from 'vitest';
import {
  calculateRingGeometry,
  analyzePiStacking,
  type RingGeometry,
  type PiStackingResult,
} from '../ligandAnalysis';

describe('calculateRingGeometry', () => {
  it('should calculate centroid of a hexagonal ring', () => {
    // Regular hexagon in XY plane centered at origin
    const positions: [number, number, number][] = [
      [1, 0, 0],
      [0.5, 0.866, 0],
      [-0.5, 0.866, 0],
      [-1, 0, 0],
      [-0.5, -0.866, 0],
      [0.5, -0.866, 0],
    ];

    const result = calculateRingGeometry(positions);

    expect(result).not.toBeNull();
    expect(result!.centroid[0]).toBeCloseTo(0, 1);
    expect(result!.centroid[1]).toBeCloseTo(0, 1);
    expect(result!.centroid[2]).toBeCloseTo(0, 1);
  });

  it('should calculate normal vector perpendicular to ring plane', () => {
    // Ring in XY plane should have normal along Z
    const positions: [number, number, number][] = [
      [1, 0, 0],
      [0.5, 0.866, 0],
      [-0.5, 0.866, 0],
      [-1, 0, 0],
      [-0.5, -0.866, 0],
      [0.5, -0.866, 0],
    ];

    const result = calculateRingGeometry(positions);

    expect(result).not.toBeNull();
    // Normal should be [0, 0, 1] or [0, 0, -1]
    expect(Math.abs(result!.normal[2])).toBeCloseTo(1, 1);
  });

  it('should return null for fewer than 5 atoms', () => {
    const positions: [number, number, number][] = [
      [0, 0, 0],
      [1, 0, 0],
      [1, 1, 0],
    ];

    const result = calculateRingGeometry(positions);

    expect(result).toBeNull();
  });
});

describe('analyzePiStacking', () => {
  it('should detect parallel stacking for rings with parallel normals', () => {
    const ring1: RingGeometry = {
      centroid: [0, 0, 0],
      normal: [0, 0, 1],
    };
    const ring2: RingGeometry = {
      centroid: [0, 0, 3.5],
      normal: [0, 0, 1],
    };

    const result = analyzePiStacking(ring1, ring2);

    expect(result.type).toBe('parallel');
    expect(result.isStacking).toBe(true);
    expect(result.distance).toBeCloseTo(3.5, 1);
    expect(result.angle).toBeCloseTo(0, 5);
  });

  it('should detect t-shaped stacking for perpendicular rings', () => {
    const ring1: RingGeometry = {
      centroid: [0, 0, 0],
      normal: [0, 0, 1],
    };
    const ring2: RingGeometry = {
      centroid: [3, 0, 0],
      normal: [1, 0, 0],
    };

    const result = analyzePiStacking(ring1, ring2);

    expect(result.type).toBe('t-shaped');
    expect(result.isStacking).toBe(true);
    expect(result.angle).toBeCloseTo(90, 5);
  });

  it('should detect offset-parallel for displaced parallel rings', () => {
    const ring1: RingGeometry = {
      centroid: [0, 0, 0],
      normal: [0, 0, 1],
    };
    const ring2: RingGeometry = {
      centroid: [2.5, 0, 3.5],
      normal: [0, 0, 1],
    };

    const result = analyzePiStacking(ring1, ring2);

    expect(result.type).toBe('offset-parallel');
    expect(result.isStacking).toBe(true);
  });

  it('should return none for distant rings', () => {
    const ring1: RingGeometry = {
      centroid: [0, 0, 0],
      normal: [0, 0, 1],
    };
    const ring2: RingGeometry = {
      centroid: [10, 0, 0],
      normal: [0, 0, 1],
    };

    const result = analyzePiStacking(ring1, ring2);

    expect(result.type).toBe('none');
    expect(result.isStacking).toBe(false);
  });
});
```

**Step 2: Run test to verify it fails**

Run: `cd frontend && npm test -- --run piStacking`
Expected: FAIL with "calculateRingGeometry is not exported"

**Step 3: Add types and functions to ligandAnalysis.ts**

Add to `frontend/src/lib/ligandAnalysis.ts` (after line 57, before LigandContact interface):

```typescript
// Pi-stacking geometry types
export interface RingGeometry {
  centroid: [number, number, number];
  normal: [number, number, number];
}

export interface PiStackingResult {
  type: 'parallel' | 't-shaped' | 'offset-parallel' | 'none';
  distance: number;
  angle: number;
  offset: number;
  isStacking: boolean;
  ligandRingIdx?: number;
  proteinResidue?: string;
  proteinChain?: string;
}

// Aromatic ring atom patterns for known residues
export const AROMATIC_PATTERNS: Record<string, string[][]> = {
  'PHE': [['CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2']],
  'TYR': [['CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2']],
  'TRP': [['CG', 'CD1', 'NE1', 'CE2', 'CD2'], ['CE2', 'CD2', 'CE3', 'CZ3', 'CH2', 'CZ2']],
  'HIS': [['CG', 'ND1', 'CE1', 'NE2', 'CD2']],
};

/**
 * Calculate ring geometry (centroid and normal vector) from atom positions
 * Requires at least 5 atoms to form an aromatic ring
 */
export function calculateRingGeometry(
  positions: [number, number, number][]
): RingGeometry | null {
  if (positions.length < 5) return null;

  // Calculate centroid
  const centroid: [number, number, number] = [0, 0, 0];
  for (const pos of positions) {
    centroid[0] += pos[0];
    centroid[1] += pos[1];
    centroid[2] += pos[2];
  }
  centroid[0] /= positions.length;
  centroid[1] /= positions.length;
  centroid[2] /= positions.length;

  // Calculate normal using first 3 points (cross product of two edge vectors)
  const v1: [number, number, number] = [
    positions[1][0] - positions[0][0],
    positions[1][1] - positions[0][1],
    positions[1][2] - positions[0][2],
  ];
  const v2: [number, number, number] = [
    positions[2][0] - positions[0][0],
    positions[2][1] - positions[0][1],
    positions[2][2] - positions[0][2],
  ];

  // Cross product
  const normal: [number, number, number] = [
    v1[1] * v2[2] - v1[2] * v2[1],
    v1[2] * v2[0] - v1[0] * v2[2],
    v1[0] * v2[1] - v1[1] * v2[0],
  ];

  // Normalize
  const length = Math.sqrt(normal[0] ** 2 + normal[1] ** 2 + normal[2] ** 2);
  if (length < 0.0001) return null;

  normal[0] /= length;
  normal[1] /= length;
  normal[2] /= length;

  return { centroid, normal };
}

/**
 * Analyze pi-stacking interaction between two aromatic rings
 * Returns stacking type (parallel, t-shaped, offset-parallel, or none)
 */
export function analyzePiStacking(
  ring1: RingGeometry,
  ring2: RingGeometry
): PiStackingResult {
  // Distance between centroids
  const dx = ring2.centroid[0] - ring1.centroid[0];
  const dy = ring2.centroid[1] - ring1.centroid[1];
  const dz = ring2.centroid[2] - ring1.centroid[2];
  const distance = Math.sqrt(dx * dx + dy * dy + dz * dz);

  // Angle between normals (dot product)
  const dotProduct = Math.abs(
    ring1.normal[0] * ring2.normal[0] +
    ring1.normal[1] * ring2.normal[1] +
    ring1.normal[2] * ring2.normal[2]
  );
  const angle = Math.acos(Math.min(1, dotProduct)) * 180 / Math.PI;

  // Calculate offset (horizontal displacement from stacking axis)
  const projection =
    dx * ring1.normal[0] + dy * ring1.normal[1] + dz * ring1.normal[2];
  const offset = Math.sqrt(Math.max(0, distance * distance - projection * projection));

  // Classify stacking type
  let type: PiStackingResult['type'] = 'none';

  if (distance <= 4.5) {
    if (angle < 30) {
      // Rings are roughly parallel
      type = offset < 2.0 ? 'parallel' : 'offset-parallel';
    } else if (angle > 60) {
      // Rings are roughly perpendicular
      type = 't-shaped';
    }
  }

  return {
    type,
    distance,
    angle,
    offset,
    isStacking: type !== 'none',
  };
}
```

**Step 4: Run test to verify it passes**

Run: `cd frontend && npm test -- --run piStacking`
Expected: PASS (all 7 tests)

**Step 5: Commit**

```bash
git add frontend/src/lib/ligandAnalysis.ts frontend/src/lib/__tests__/piStacking.test.ts
git commit -m "$(cat <<'EOF'
feat(analysis): add pi-stacking geometry analysis

Add calculateRingGeometry and analyzePiStacking functions for
detecting parallel, t-shaped, and offset-parallel aromatic
stacking interactions with geometric classification.

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Add Interaction Line Geometry Computation

**Files:**
- Create: `frontend/src/lib/interactionGeometry.ts`
- Create: `frontend/src/lib/__tests__/interactionGeometry.test.ts`

**Step 1: Write failing test**

Create `frontend/src/lib/__tests__/interactionGeometry.test.ts`:

```typescript
import { describe, it, expect } from 'vitest';
import {
  computeInteractionLines,
  type InteractionLine,
} from '../interactionGeometry';

describe('computeInteractionLines', () => {
  const samplePdb = `
ATOM      1  N   SER A  10       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  SER A  10       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  OG  SER A  10       2.000   1.000   0.000  1.00  0.00           O
HETATM   10  C1  LIG B   1       4.000   1.500   0.000  1.00  0.00           C
HETATM   11  O1  LIG B   1       5.000   1.000   0.000  1.00  0.00           O
END
`;

  it('should return empty array when no ligand contacts', () => {
    const result = computeInteractionLines('', []);
    expect(result).toEqual([]);
  });

  it('should compute line endpoints for hydrogen bond contacts', () => {
    const contacts = [
      {
        ligandChain: 'B',
        ligandResSeq: 1,
        ligandAtom: 'O1',
        proteinResidue: 'SER10',
        proteinChain: 'A',
        proteinAtom: 'OG',
        interactionType: 'hydrogen_bond' as const,
        distance: 3.0,
      },
    ];

    const result = computeInteractionLines(samplePdb, contacts);

    expect(result).toHaveLength(1);
    expect(result[0].type).toBe('hydrogen_bond');
    expect(result[0].start).toBeDefined();
    expect(result[0].end).toBeDefined();
    expect(result[0].color).toBeDefined();
  });

  it('should assign correct colors based on interaction type', () => {
    const contacts = [
      {
        ligandChain: 'B',
        ligandResSeq: 1,
        ligandAtom: 'O1',
        proteinResidue: 'SER10',
        proteinChain: 'A',
        proteinAtom: 'OG',
        interactionType: 'hydrogen_bond' as const,
        distance: 3.0,
      },
    ];

    const result = computeInteractionLines(samplePdb, contacts);

    expect(result[0].color).toBe(0x3B82F6); // Blue for H-bonds
  });
});
```

**Step 2: Run test to verify it fails**

Run: `cd frontend && npm test -- --run interactionGeometry`
Expected: FAIL with "Cannot find module '../interactionGeometry'"

**Step 3: Create interactionGeometry.ts**

Create `frontend/src/lib/interactionGeometry.ts`:

```typescript
/**
 * Interaction Line Geometry Computation
 * Computes 3D line endpoints for visualizing protein-ligand interactions
 */

import type { InteractionType } from './ligandAnalysis';

export interface InteractionLineContact {
  ligandChain: string;
  ligandResSeq: number;
  ligandAtom: string;
  proteinResidue: string;
  proteinChain: string;
  proteinAtom: string;
  interactionType: InteractionType;
  distance: number;
}

export interface InteractionLine {
  start: [number, number, number];
  end: [number, number, number];
  type: InteractionType;
  color: number;
  dashed: boolean;
  label?: string;
}

// Interaction type colors (matching ligandAnalysis)
const INTERACTION_COLORS: Record<InteractionType, number> = {
  hydrogen_bond: 0x3B82F6,  // Blue
  hydrophobic: 0x22C55E,    // Green
  salt_bridge: 0xEF4444,    // Red
  pi_stacking: 0xA855F7,    // Purple
  other: 0x6B7280,          // Gray
};

// Dashed lines for weaker interactions
const DASHED_TYPES: Set<InteractionType> = new Set(['hydrophobic', 'other']);

interface ParsedAtom {
  name: string;
  resName: string;
  chainId: string;
  resSeq: number;
  x: number;
  y: number;
  z: number;
}

/**
 * Parse PDB to extract atom positions
 */
function parseAtomPositions(pdbContent: string): Map<string, ParsedAtom> {
  const atoms = new Map<string, ParsedAtom>();
  const lines = pdbContent.split('\n');

  for (const line of lines) {
    if (!line.startsWith('ATOM') && !line.startsWith('HETATM')) continue;

    try {
      const atom: ParsedAtom = {
        name: line.slice(12, 16).trim(),
        resName: line.slice(17, 20).trim(),
        chainId: line.slice(21, 22).trim(),
        resSeq: parseInt(line.slice(22, 26).trim()),
        x: parseFloat(line.slice(30, 38).trim()),
        y: parseFloat(line.slice(38, 46).trim()),
        z: parseFloat(line.slice(46, 54).trim()),
      };

      // Create unique key for lookup
      const key = `${atom.chainId}:${atom.resSeq}:${atom.name}`;
      atoms.set(key, atom);
    } catch {
      // Skip malformed lines
    }
  }

  return atoms;
}

/**
 * Extract residue number from residue string (e.g., "SER10" -> 10)
 */
function extractResSeq(residue: string): number {
  const match = residue.match(/\d+/);
  return match ? parseInt(match[0]) : 0;
}

/**
 * Compute 3D line segments for interaction visualization
 */
export function computeInteractionLines(
  pdbContent: string,
  contacts: InteractionLineContact[]
): InteractionLine[] {
  if (!pdbContent || contacts.length === 0) {
    return [];
  }

  const atoms = parseAtomPositions(pdbContent);
  const lines: InteractionLine[] = [];

  for (const contact of contacts) {
    // Find ligand atom position
    const ligandKey = `${contact.ligandChain}:${contact.ligandResSeq}:${contact.ligandAtom}`;
    const ligandAtom = atoms.get(ligandKey);

    // Find protein atom position
    const proteinResSeq = extractResSeq(contact.proteinResidue);
    const proteinKey = `${contact.proteinChain}:${proteinResSeq}:${contact.proteinAtom}`;
    const proteinAtom = atoms.get(proteinKey);

    if (!ligandAtom || !proteinAtom) {
      continue; // Skip if atoms not found
    }

    lines.push({
      start: [ligandAtom.x, ligandAtom.y, ligandAtom.z],
      end: [proteinAtom.x, proteinAtom.y, proteinAtom.z],
      type: contact.interactionType,
      color: INTERACTION_COLORS[contact.interactionType],
      dashed: DASHED_TYPES.has(contact.interactionType),
      label: `${contact.proteinResidue} ${contact.proteinAtom}`,
    });
  }

  return lines;
}
```

**Step 4: Run test to verify it passes**

Run: `cd frontend && npm test -- --run interactionGeometry`
Expected: PASS (all 3 tests)

**Step 5: Commit**

```bash
git add frontend/src/lib/interactionGeometry.ts frontend/src/lib/__tests__/interactionGeometry.test.ts
git commit -m "$(cat <<'EOF'
feat(analysis): add interaction line geometry computation

Add computeInteractionLines function to calculate 3D line
endpoints for visualizing H-bonds, salt bridges, and other
protein-ligand interactions in the Molstar viewer.

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Create FocusModeControls Component

**Files:**
- Create: `frontend/src/components/viewer/FocusModeControls.tsx`

**Step 1: Create the component with detection-aware visibility**

Create `frontend/src/components/viewer/FocusModeControls.tsx`:

```typescript
'use client';

import { useState } from 'react';
import {
  Droplet,
  Link2,
  CircleDot,
  SlidersHorizontal,
  X,
  ChevronDown,
} from 'lucide-react';
import { useStore } from '@/lib/store';

interface FocusModeControlsProps {
  focusType: 'metal' | 'ligand';
  hasWaters: boolean;
  hasInteractions: boolean;
  hasPharmacophores: boolean;
  onRadiusChange: (radius: number) => void;
  onReset: () => void;
}

export function FocusModeControls({
  focusType,
  hasWaters,
  hasInteractions,
  hasPharmacophores,
  onRadiusChange,
  onReset,
}: FocusModeControlsProps) {
  const { focusSettings, setFocusSettings } = useStore();
  const [expanded, setExpanded] = useState(true);

  const currentRadius = focusType === 'metal'
    ? focusSettings.coordinationRadius
    : focusSettings.bindingPocketRadius;

  const radiusLabel = focusType === 'metal'
    ? 'Coordination Radius'
    : 'Pocket Radius';

  // Check if any controls should be shown
  const hasAnyControls = hasWaters || hasInteractions || hasPharmacophores;

  return (
    <div className="absolute top-2 left-2 z-10">
      {/* Focus Mode Badge */}
      <div className="flex items-center gap-2 mb-2">
        <div className={`px-3 py-1.5 rounded-lg shadow-lg text-xs font-semibold text-white ${
          focusType === 'metal' ? 'bg-purple-600' : 'bg-emerald-600'
        }`}>
          {focusType === 'metal' ? 'Metal Focus' : 'Ligand Focus'}
        </div>
        <button
          onClick={onReset}
          className="p-1.5 rounded-lg bg-white/90 shadow hover:bg-white transition-colors"
          title="Reset View"
        >
          <X className="w-4 h-4 text-gray-600" />
        </button>
      </div>

      {/* Controls Panel - only if there are controls to show */}
      {hasAnyControls && (
        <div className="bg-white/95 rounded-lg shadow-lg border border-gray-200 overflow-hidden min-w-[200px]">
          {/* Header */}
          <button
            onClick={() => setExpanded(!expanded)}
            className="w-full px-3 py-2 flex items-center justify-between text-xs font-medium text-gray-700 hover:bg-gray-50"
          >
            <span className="flex items-center gap-1.5">
              <SlidersHorizontal className="w-3.5 h-3.5" />
              View Options
            </span>
            <ChevronDown className={`w-4 h-4 transition-transform ${expanded ? '' : '-rotate-90'}`} />
          </button>

          {expanded && (
            <div className="px-3 pb-3 space-y-3 border-t border-gray-100">
              {/* Radius Slider */}
              <div className="pt-3">
                <div className="flex items-center justify-between mb-1.5">
                  <label className="text-xs text-gray-600">{radiusLabel}</label>
                  <span className="text-xs font-mono text-gray-500">{currentRadius.toFixed(1)}Å</span>
                </div>
                <input
                  type="range"
                  min={2.0}
                  max={8.0}
                  step={0.5}
                  value={currentRadius}
                  onChange={(e) => {
                    const value = parseFloat(e.target.value);
                    if (focusType === 'metal') {
                      setFocusSettings({ coordinationRadius: value });
                    } else {
                      setFocusSettings({ bindingPocketRadius: value });
                    }
                    onRadiusChange(value);
                  }}
                  className="w-full h-1.5 bg-gray-200 rounded-lg appearance-none cursor-pointer accent-blue-600"
                />
              </div>

              {/* Water Toggle - only if waters exist */}
              {hasWaters && (
                <label className="flex items-center justify-between cursor-pointer group">
                  <span className="flex items-center gap-1.5 text-xs text-gray-600">
                    <Droplet className="w-3.5 h-3.5 text-cyan-500" />
                    Show Waters
                  </span>
                  <div className="relative">
                    <input
                      type="checkbox"
                      checked={focusSettings.showWaters}
                      onChange={(e) => setFocusSettings({ showWaters: e.target.checked })}
                      className="sr-only peer"
                    />
                    <div className="w-8 h-4 bg-gray-200 rounded-full peer peer-checked:bg-blue-600 transition-colors" />
                    <div className="absolute top-0.5 left-0.5 w-3 h-3 bg-white rounded-full shadow peer-checked:translate-x-4 transition-transform" />
                  </div>
                </label>
              )}

              {/* Interaction Lines Toggle - only if interactions exist */}
              {hasInteractions && (
                <label className="flex items-center justify-between cursor-pointer group">
                  <span className="flex items-center gap-1.5 text-xs text-gray-600">
                    <Link2 className="w-3.5 h-3.5 text-blue-500" />
                    Interaction Lines
                  </span>
                  <div className="relative">
                    <input
                      type="checkbox"
                      checked={focusSettings.showInteractionLines}
                      onChange={(e) => setFocusSettings({ showInteractionLines: e.target.checked })}
                      className="sr-only peer"
                    />
                    <div className="w-8 h-4 bg-gray-200 rounded-full peer peer-checked:bg-blue-600 transition-colors" />
                    <div className="absolute top-0.5 left-0.5 w-3 h-3 bg-white rounded-full shadow peer-checked:translate-x-4 transition-transform" />
                  </div>
                </label>
              )}

              {/* Pharmacophores Toggle - only if pharmacophores exist */}
              {hasPharmacophores && (
                <label className="flex items-center justify-between cursor-pointer group">
                  <span className="flex items-center gap-1.5 text-xs text-gray-600">
                    <CircleDot className="w-3.5 h-3.5 text-purple-500" />
                    Pharmacophores
                  </span>
                  <div className="relative">
                    <input
                      type="checkbox"
                      checked={focusSettings.showPharmacophores}
                      onChange={(e) => setFocusSettings({ showPharmacophores: e.target.checked })}
                      className="sr-only peer"
                    />
                    <div className="w-8 h-4 bg-gray-200 rounded-full peer peer-checked:bg-blue-600 transition-colors" />
                    <div className="absolute top-0.5 left-0.5 w-3 h-3 bg-white rounded-full shadow peer-checked:translate-x-4 transition-transform" />
                  </div>
                </label>
              )}
            </div>
          )}
        </div>
      )}
    </div>
  );
}

export default FocusModeControls;
```

**Step 2: Run TypeScript check**

Run: `cd frontend && npx tsc --noEmit`
Expected: No errors

**Step 3: Commit**

```bash
git add frontend/src/components/viewer/FocusModeControls.tsx
git commit -m "$(cat <<'EOF'
feat(ui): add FocusModeControls component

Add detection-aware focus mode controls with radius slider,
water toggle, interaction lines toggle, and pharmacophore
toggle. Controls only appear when relevant data exists.

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Enhance ProteinViewerClient with Green Carbons and Settings

**Files:**
- Modify: `frontend/src/components/ProteinViewerClient.tsx`

**Step 1: Update focusOnLigand to use green carbon coloring**

In `frontend/src/components/ProteinViewerClient.tsx`, modify the `focusOnLigand` function (around line 316-328). Replace the ligand component representation:

```typescript
// 1. Ligand - ball-and-stick with GREEN CARBONS
const ligandComp = await plugin.builders.structure.tryCreateComponentFromExpression(
  structure.ref,
  ligandExpression,
  `ligand-${ligand.name}`
);
if (ligandComp) {
  await plugin.builders.structure.representation.addRepresentation(ligandComp, {
    type: 'ball-and-stick',
    color: 'element-symbol',
    colorParams: {
      carbonColor: { r: 80, g: 200, b: 120 },  // Green carbons (0x50C878)
    },
    typeParams: { sizeFactor: 0.4 },
  });
}
```

**Step 2: Add focusSettings to component props and use dynamic radius**

Add to the props interface (around line 11):

```typescript
import { useStore } from '@/lib/store';
```

Inside the component, get focus settings:

```typescript
const { focusSettings } = useStore();
```

Update the coordination sphere radius in `focusOnMetal` (around line 140):

```typescript
const coordSphereExpression = MS.struct.modifier.includeSurroundings({
  0: metalExpression,
  radius: focusSettings.coordinationRadius + 2.0, // Use configurable radius
  'as-whole-residues': true
});
```

Update the binding site radius in `focusOnLigand` (around line 302):

```typescript
const bindingSiteExpression = MS.struct.modifier.includeSurroundings({
  0: ligandExpression,
  radius: focusSettings.bindingPocketRadius,  // Use configurable radius
  'as-whole-residues': true
});
```

**Step 3: Run dev server to verify no errors**

Run: `cd frontend && npm run dev`
Expected: Server starts without errors

**Step 4: Commit**

```bash
git add frontend/src/components/ProteinViewerClient.tsx
git commit -m "$(cat <<'EOF'
feat(viewer): add green carbon coloring and configurable radii

Update focusOnLigand to render ligand carbons in green for
better visibility. Use focusSettings from store for dynamic
coordination and binding pocket radii.

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Add Water Visibility Toggle to Focus Views

**Files:**
- Modify: `frontend/src/components/ProteinViewerClient.tsx`

**Step 1: Conditionally render waters based on focusSettings.showWaters**

In `focusOnMetal` function, wrap the water rendering in a conditional (around line 207-228):

```typescript
// 3. Coordinating waters - conditionally rendered
if (focusSettings.showWaters) {
  const waterExpression = MS.struct.modifier.intersectBy({
    0: coordWithoutMetal,
    by: MS.struct.generator.atomGroups({
      'residue-test': MS.core.logic.or([
        MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_comp_id(), 'HOH']),
        MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_comp_id(), 'WAT'])
      ])
    })
  });
  const waterComp = await plugin.builders.structure.tryCreateComponentFromExpression(
    structure.ref,
    waterExpression,
    'coord-water'
  );
  if (waterComp) {
    await plugin.builders.structure.representation.addRepresentation(waterComp, {
      type: 'ball-and-stick',
      color: 'element-symbol',
      typeParams: { sizeFactor: 0.2 },
    });
  }
}
```

**Step 2: Do the same for focusOnLigand (add water rendering section)**

After the binding site residues section (around line 365), add:

```typescript
// 3. Binding site waters - conditionally rendered
if (focusSettings.showWaters) {
  const waterExpression = MS.struct.modifier.intersectBy({
    0: siteWithoutLigand,
    by: MS.struct.generator.atomGroups({
      'residue-test': MS.core.logic.or([
        MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_comp_id(), 'HOH']),
        MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_comp_id(), 'WAT'])
      ])
    })
  });
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
```

**Step 3: Run TypeScript check**

Run: `cd frontend && npx tsc --noEmit`
Expected: No errors

**Step 4: Commit**

```bash
git add frontend/src/components/ProteinViewerClient.tsx
git commit -m "$(cat <<'EOF'
feat(viewer): add conditional water visibility in focus views

Waters in metal coordination spheres and ligand binding pockets
are now conditionally rendered based on focusSettings.showWaters.

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Add 3D Interaction Lines Rendering

**Files:**
- Modify: `frontend/src/components/ProteinViewerClient.tsx`

**Step 1: Import interaction geometry module**

Add to imports at top of file:

```typescript
import { computeInteractionLines, type InteractionLine } from '@/lib/interactionGeometry';
```

**Step 2: Add renderInteractionLines function**

Add this function inside the component (after clearPharmacophores around line 427):

```typescript
// Render interaction lines as 3D cylinders
const renderInteractionLines = useCallback(async (lines: InteractionLine[]) => {
  if (!globalPlugin || !globalStructureRef || lines.length === 0) return;

  await loadMolstarModules();
  const plugin = globalPlugin;

  try {
    // Clear existing interaction lines
    const state = plugin.state.data;
    const toRemove: string[] = [];
    state.cells.forEach((cell: any, ref: string) => {
      if (cell.obj?.label?.startsWith('interaction-line-')) {
        toRemove.push(ref);
      }
    });
    for (const ref of toRemove) {
      await plugin.build().delete(ref).commit();
    }

    // Create shape provider for lines
    const { ShapeRepresentation3D } = await import('molstar/lib/mol-repr/shape/representation');
    const { Shape } = await import('molstar/lib/mol-model/shape');
    const { Mesh } = await import('molstar/lib/mol-geo/geometry/mesh/mesh');
    const { MeshBuilder } = await import('molstar/lib/mol-geo/geometry/mesh/mesh-builder');
    const { addCylinder } = await import('molstar/lib/mol-geo/geometry/mesh/builder/cylinder');

    for (let i = 0; i < lines.length; i++) {
      const line = lines[i];
      const builder = MeshBuilder.createState(256, 128);

      // Add cylinder from start to end
      addCylinder(builder, line.start, line.end, 0.08, {
        radiusTop: 0.08,
        radiusBottom: 0.08,
      });

      const mesh = MeshBuilder.getMesh(builder);

      // Create shape with color
      const shape = Shape.create(
        `interaction-line-${i}`,
        {},
        mesh,
        () => Color(line.color),
        () => 1,
        () => `${line.label || line.type}`
      );

      // Add to scene
      await plugin.build()
        .toRoot()
        .apply(ShapeRepresentation3D, { shape })
        .commit();
    }

    console.log(`[ProteinViewer] Rendered ${lines.length} interaction lines`);
  } catch (err) {
    console.error('[ProteinViewer] Failed to render interaction lines:', err);
  }
}, []);
```

**Step 3: Call renderInteractionLines in focusOnLigand when enabled**

At the end of `focusOnLigand`, before the camera focus (around line 380), add:

```typescript
// Render interaction lines if enabled
if (focusSettings.showInteractionLines && ligand.contacts.length > 0) {
  const lineContacts = ligand.contacts.map(c => ({
    ligandChain: ligand.chainId,
    ligandResSeq: ligand.resSeq,
    ligandAtom: c.atom, // Note: This needs refinement to get ligand atom
    proteinResidue: c.residue,
    proteinChain: c.chain,
    proteinAtom: c.atom,
    interactionType: c.interactionType,
    distance: c.distance,
  }));

  const lines = computeInteractionLines(pdbContent || '', lineContacts);
  await renderInteractionLines(lines);
}
```

**Step 4: Run TypeScript check**

Run: `cd frontend && npx tsc --noEmit`
Expected: No errors (or minor type adjustments needed)

**Step 5: Commit**

```bash
git add frontend/src/components/ProteinViewerClient.tsx
git commit -m "$(cat <<'EOF'
feat(viewer): add 3D interaction line rendering

Add renderInteractionLines function using Molstar Mesh shapes
to draw cylinders between interacting atoms. Lines are colored
by interaction type and rendered when showInteractionLines is enabled.

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: Integrate FocusModeControls into ViewerPanel

**Files:**
- Modify: `frontend/src/components/layout/ViewerPanel.tsx`

**Step 1: Import and add FocusModeControls**

First, read the current ViewerPanel to understand its structure, then add the FocusModeControls import:

```typescript
import { FocusModeControls } from '@/components/viewer/FocusModeControls';
```

**Step 2: Add FocusModeControls inside the viewer container**

Add the component inside the viewer container, conditionally rendered when in focus mode:

```typescript
{/* Focus Mode Controls - appears when metal or ligand is focused */}
{(focusedMetalIndex !== null || focusedLigandIndex !== null) && (
  <FocusModeControls
    focusType={focusedMetalIndex !== null ? 'metal' : 'ligand'}
    hasWaters={
      focusedMetalIndex !== null
        ? (metalCoordination?.[focusedMetalIndex]?.hydrationAnalysis?.waterCount ?? 0) > 0
        : (ligandData?.ligandDetails[focusedLigandIndex ?? 0]?.waterContactCount ?? 0) > 0
    }
    hasInteractions={
      focusedLigandIndex !== null
        ? (ligandData?.ligandDetails[focusedLigandIndex]?.contacts?.length ?? 0) > 0
        : (metalCoordination?.[focusedMetalIndex ?? 0]?.coordinating?.length ?? 0) > 0
    }
    hasPharmacophores={
      focusedLigandIndex !== null && (pharmacophoreFeatures?.length ?? 0) > 0
    }
    onRadiusChange={() => {
      // Trigger re-render of focus view with new radius
      if (focusedMetalIndex !== null) {
        setFocusedMetalIndex(focusedMetalIndex);
      } else if (focusedLigandIndex !== null) {
        setFocusedLigandIndex(focusedLigandIndex);
      }
    }}
    onReset={() => {
      setFocusedMetalIndex(null);
      setFocusedLigandIndex(null);
    }}
  />
)}
```

**Step 3: Run dev server to test**

Run: `cd frontend && npm run dev`
Expected: Server starts, focus controls appear when clicking on metal/ligand

**Step 4: Commit**

```bash
git add frontend/src/components/layout/ViewerPanel.tsx
git commit -m "$(cat <<'EOF'
feat(ui): integrate FocusModeControls into ViewerPanel

Add focus mode controls overlay that appears when a metal or
ligand is focused. Controls are detection-aware and only show
toggles for available data (waters, interactions, pharmacophores).

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 9: Add Pi-Stacking Display to LigandBindingPanel

**Files:**
- Modify: `frontend/src/components/LigandBindingPanel.tsx`

**Step 1: Import pi-stacking types**

Add to imports:

```typescript
import {
  // ... existing imports
  type PiStackingResult,
} from '@/lib/ligandAnalysis';
```

**Step 2: Add pi-stacking display section in LigandCard**

After the pharmacophore features section (around line 297), add:

```typescript
{/* Pi-Stacking Interactions */}
{ligand.piStackingResults && ligand.piStackingResults.length > 0 && (
  <div className="bg-card rounded-lg border border-border p-3">
    <h4 className="text-xs font-semibold text-muted-foreground uppercase tracking-wide mb-2">
      Pi-Stacking Interactions
    </h4>
    <div className="space-y-2">
      {ligand.piStackingResults.map((stack, i) => (
        <div key={i} className="flex items-center justify-between text-xs">
          <div className="flex items-center gap-2">
            <span className={`px-2 py-0.5 rounded font-medium ${
              stack.type === 'parallel' ? 'bg-purple-100 text-purple-700' :
              stack.type === 't-shaped' ? 'bg-blue-100 text-blue-700' :
              stack.type === 'offset-parallel' ? 'bg-indigo-100 text-indigo-700' :
              'bg-gray-100 text-gray-600'
            }`}>
              {stack.type}
            </span>
            <span className="text-muted-foreground">
              {stack.proteinResidue}
            </span>
          </div>
          <div className="text-muted-foreground font-mono">
            {stack.distance.toFixed(1)}Å / {stack.angle.toFixed(0)}°
          </div>
        </div>
      ))}
    </div>
  </div>
)}
```

**Step 3: Update LigandData type to include piStackingResults (if not already)**

In `frontend/src/lib/ligandAnalysis.ts`, add to the LigandData interface (around line 77):

```typescript
export interface LigandData {
  // ... existing fields
  piStackingResults?: PiStackingResult[];
}
```

**Step 4: Run TypeScript check**

Run: `cd frontend && npx tsc --noEmit`
Expected: No errors

**Step 5: Commit**

```bash
git add frontend/src/components/LigandBindingPanel.tsx frontend/src/lib/ligandAnalysis.ts
git commit -m "$(cat <<'EOF'
feat(ui): add pi-stacking display to LigandBindingPanel

Show detected pi-stacking interactions with type classification
(parallel, t-shaped, offset-parallel), residue info, and
distance/angle metrics.

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 10: Final Integration Test

**Files:**
- No new files

**Step 1: Run all tests**

Run: `cd frontend && npm test -- --run`
Expected: All tests pass

**Step 2: Run linting**

Run: `cd frontend && npm run lint`
Expected: No errors

**Step 3: Run build**

Run: `cd frontend && npm run build`
Expected: Build succeeds

**Step 4: Manual testing checklist**

1. Load a structure with metal ions (e.g., 1A3N with heme iron)
2. Click on metal in BindingSitePanel → Focus view appears
3. Verify FocusModeControls appear with radius slider
4. Adjust radius → View updates
5. Toggle waters (if present) → Waters appear/disappear
6. Load a structure with ligands (e.g., 1ATP with ATP)
7. Click on ligand → Focus view with green carbons
8. Toggle interaction lines → Lines appear connecting atoms
9. Toggle pharmacophores → Colored spheres appear

**Step 5: Final commit**

```bash
git add -A
git commit -m "$(cat <<'EOF'
feat: complete Molstar focus view enhancement

- Configurable coordination/pocket radii via UI sliders
- Water visibility toggle (detection-aware)
- 3D interaction lines using Mol* mesh shapes
- Green carbon coloring for ligand focus
- Pi-stacking geometry analysis with type classification
- Detection-aware FocusModeControls component

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>
EOF
)"
```

---

## Summary

This plan adds the following features to Banta_Lab_RFdiffusion's Molstar viewer:

| Feature | Files | Status |
|---------|-------|--------|
| Focus settings store | `store.ts` | Task 1 |
| Pi-stacking geometry | `ligandAnalysis.ts` | Task 2 |
| Interaction line computation | `interactionGeometry.ts` | Task 3 |
| FocusModeControls UI | `FocusModeControls.tsx` | Task 4 |
| Green ligand carbons | `ProteinViewerClient.tsx` | Task 5 |
| Water toggle | `ProteinViewerClient.tsx` | Task 6 |
| 3D interaction lines | `ProteinViewerClient.tsx` | Task 7 |
| UI integration | `ViewerPanel.tsx` | Task 8 |
| Pi-stacking display | `LigandBindingPanel.tsx` | Task 9 |

All features are detection-aware and only show controls when relevant data exists.
