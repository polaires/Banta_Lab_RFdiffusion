/**
 * Metal Coordination Analysis for MolStar
 * Extracted and adapted from Protein_engineering_tools/src/components/ProteinViewer.tsx
 *
 * Provides:
 * - Metal ion detection (lanthanides, transition metals, etc.)
 * - Coordination geometry analysis using Continuous Shape Measures (CShM)
 * - Binding site classification (functional vs crystal artifact)
 * - Hydration analysis
 */

import type { Vec3 as Vec3Type } from 'molstar/lib/mol-math/linear-algebra';

// Types
export interface CoordinatingAtom {
  atom: string;
  residue: string;
  chain: string;
  distance: number;
  isWater: boolean;
  position: [number, number, number];
}

export interface GeometryAnalysis {
  coordinationNumber: number;
  geometryType: string;
  idealGeometry: string;
  angles: { atom1: string; atom2: string; angle: number }[];
  avgAngle: number;
  rmsd: number;
  distortion: 'ideal' | 'low' | 'moderate' | 'high' | 'severe';
}

export interface HydrationAnalysis {
  waterCount: number;
  proteinLigandCount: number;
  hydrationState: 'fully_hydrated' | 'partially_hydrated' | 'dehydrated';
  expectedHydration: number;
  waterDisplacement: number;
  hydrationNote: string;
}

export interface MetalCoordination {
  element: string;
  info: string;
  chainId: string;
  resSeq: number;
  resName: string;
  pos: [number, number, number];
  coordinating: CoordinatingAtom[];
  geometry: GeometryAnalysis | null;
  bindingSiteType: 'functional' | 'crystal_artifact' | 'uncertain';
  bindingSiteReason: string;
  hydrationAnalysis: HydrationAnalysis;
}

// Comprehensive list of ALL metal elements
export const METAL_ELEMENTS = new Set([
  // Alkali metals
  'LI', 'NA', 'K', 'RB', 'CS', 'FR',
  // Alkaline earth metals
  'BE', 'MG', 'CA', 'SR', 'BA', 'RA',
  // Transition metals (3d, 4d, 5d)
  'SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN',
  'Y', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD',
  'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG',
  // Post-transition metals
  'AL', 'GA', 'IN', 'SN', 'TL', 'PB', 'BI', 'PO',
  // Lanthanides (rare earths)
  'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU',
  // Actinides
  'AC', 'TH', 'PA', 'U', 'NP', 'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM', 'MD', 'NO', 'LR'
]);

// Expected hydration numbers based on metal type (from literature)
const EXPECTED_HYDRATION_BY_METAL: Record<string, number> = {
  // Alkali metals - high hydration
  'NA': 6, 'K': 6, 'LI': 4, 'RB': 8, 'CS': 8,
  // Alkaline earth - variable
  'MG': 6, 'CA': 7, 'SR': 8, 'BA': 8,
  // Common transition metals
  'ZN': 4, 'CU': 4, 'FE': 6, 'MN': 6, 'CO': 6, 'NI': 6,
  // Lanthanides - high coordination
  'LA': 9, 'CE': 9, 'PR': 9, 'ND': 9, 'EU': 9, 'GD': 9, 'TB': 9, 'DY': 9,
  // Default
  'DEFAULT': 6
};

// Reference polyhedra coordinates (normalized, centered at origin)
// These are ideal vertex positions for each geometry type
const REFERENCE_POLYHEDRA: {
  [cn: number]: { name: string; vertices: number[][] }[];
} = {
  2: [
    { name: 'Linear', vertices: [[0, 0, 1], [0, 0, -1]] }
  ],
  3: [
    { name: 'Trigonal Planar', vertices: [
      [1, 0, 0], [-0.5, 0.866, 0], [-0.5, -0.866, 0]
    ]},
    { name: 'T-shaped', vertices: [
      [1, 0, 0], [-1, 0, 0], [0, 1, 0]
    ]},
    { name: 'Trigonal Pyramidal', vertices: [
      [0.943, 0, -0.333], [-0.471, 0.816, -0.333], [-0.471, -0.816, -0.333]
    ]}
  ],
  4: [
    { name: 'Tetrahedral', vertices: [
      [1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]
    ].map(v => v.map(c => c / Math.sqrt(3)))},
    { name: 'Square Planar', vertices: [
      [1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]
    ]},
    { name: 'See-saw', vertices: [
      [0, 0, 1], [0, 0, -1], [1, 0, 0], [-0.5, 0.866, 0]
    ]}
  ],
  5: [
    { name: 'Trigonal Bipyramidal', vertices: [
      [0, 0, 1], [0, 0, -1], // axial
      [1, 0, 0], [-0.5, 0.866, 0], [-0.5, -0.866, 0] // equatorial
    ]},
    { name: 'Square Pyramidal', vertices: [
      [0, 0, 1], // apex
      [1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0] // base
    ]},
    { name: 'Pentagonal Planar', vertices: [
      [1, 0, 0], [0.309, 0.951, 0], [-0.809, 0.588, 0], [-0.809, -0.588, 0], [0.309, -0.951, 0]
    ]}
  ],
  6: [
    { name: 'Octahedral', vertices: [
      [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]
    ]},
    { name: 'Trigonal Prismatic', vertices: [
      [1, 0, 0.5], [-0.5, 0.866, 0.5], [-0.5, -0.866, 0.5],
      [1, 0, -0.5], [-0.5, 0.866, -0.5], [-0.5, -0.866, -0.5]
    ]},
    { name: 'Pentagonal Pyramidal', vertices: [
      [0, 0, 1], // apex
      [1, 0, 0], [0.309, 0.951, 0], [-0.809, 0.588, 0], [-0.809, -0.588, 0], [0.309, -0.951, 0]
    ]}
  ],
  7: [
    { name: 'Pentagonal Bipyramidal', vertices: [
      [0, 0, 1], [0, 0, -1], // axial
      [1, 0, 0], [0.309, 0.951, 0], [-0.809, 0.588, 0], [-0.809, -0.588, 0], [0.309, -0.951, 0]
    ]},
    { name: 'Capped Octahedral', vertices: [
      [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1],
      [0.577, 0.577, 0.577] // cap
    ]},
    { name: 'Capped Trigonal Prismatic', vertices: [
      [1, 0, 0.5], [-0.5, 0.866, 0.5], [-0.5, -0.866, 0.5],
      [1, 0, -0.5], [-0.5, 0.866, -0.5], [-0.5, -0.866, -0.5],
      [0, 0, 0] // face cap
    ]}
  ],
  8: [
    { name: 'Square Antiprismatic', vertices: (() => {
      const h = 0.5;
      const r = 1;
      const top = [[r, 0, h], [0, r, h], [-r, 0, h], [0, -r, h]];
      const a = Math.PI / 4;
      const bottom = [
        [r * Math.cos(a), r * Math.sin(a), -h],
        [-r * Math.sin(a), r * Math.cos(a), -h],
        [-r * Math.cos(a), -r * Math.sin(a), -h],
        [r * Math.sin(a), -r * Math.cos(a), -h]
      ];
      return [...top, ...bottom];
    })()},
    { name: 'Dodecahedral', vertices: [
      [1, 0, 0.5], [-1, 0, 0.5], [0, 1, -0.5], [0, -1, -0.5],
      [0.707, 0.707, 0], [-0.707, 0.707, 0], [-0.707, -0.707, 0], [0.707, -0.707, 0]
    ]},
    { name: 'Bicapped Trigonal Prismatic', vertices: [
      [1, 0, 0.577], [-0.5, 0.866, 0.577], [-0.5, -0.866, 0.577],
      [1, 0, -0.577], [-0.5, 0.866, -0.577], [-0.5, -0.866, -0.577],
      [0, 0, 1], [0, 0, -1]
    ]},
    { name: 'Cubic', vertices: [
      [1, 1, 1], [1, 1, -1], [1, -1, 1], [1, -1, -1],
      [-1, 1, 1], [-1, 1, -1], [-1, -1, 1], [-1, -1, -1]
    ].map(v => v.map(c => c / Math.sqrt(3)))}
  ],
  9: [
    { name: 'Tricapped Trigonal Prismatic', vertices: [
      [1, 0, 0.816], [-0.5, 0.866, 0.816], [-0.5, -0.866, 0.816],
      [1, 0, -0.816], [-0.5, 0.866, -0.816], [-0.5, -0.866, -0.816],
      [0.25, 0.433, 0], [-0.5, 0, 0], [0.25, -0.433, 0]
    ]},
    { name: 'Capped Square Antiprismatic', vertices: (() => {
      const h = 0.5;
      const r = 1;
      const top = [[r, 0, h], [0, r, h], [-r, 0, h], [0, -r, h]];
      const a = Math.PI / 4;
      const bottom = [
        [r * Math.cos(a), r * Math.sin(a), -h],
        [-r * Math.sin(a), r * Math.cos(a), -h],
        [-r * Math.cos(a), -r * Math.sin(a), -h],
        [r * Math.sin(a), -r * Math.cos(a), -h]
      ];
      return [...top, ...bottom, [0, 0, 1]];
    })()},
    { name: 'Muffin', vertices: [
      [1, 0, 0], [0.5, 0.866, 0], [-0.5, 0.866, 0], [-1, 0, 0], [-0.5, -0.866, 0], [0.5, -0.866, 0],
      [0.5, 0.289, 0.816], [-0.5, 0.289, 0.816], [0, -0.577, 0.816]
    ]}
  ]
};

// Vector math helpers (for use without full Vec3 dependency)
function vec3Distance(a: [number, number, number], b: [number, number, number]): number {
  const dx = a[0] - b[0];
  const dy = a[1] - b[1];
  const dz = a[2] - b[2];
  return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

function vec3Sub(a: [number, number, number], b: [number, number, number]): [number, number, number] {
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
}

function vec3Dot(a: [number, number, number], b: [number, number, number]): number {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

function vec3Magnitude(v: [number, number, number]): number {
  return Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

// Calculate angle between three points (in degrees) - angle at point B
function calculateAngle(
  a: [number, number, number],
  b: [number, number, number],
  c: [number, number, number]
): number {
  const ba = vec3Sub(a, b);
  const bc = vec3Sub(c, b);

  const dot = vec3Dot(ba, bc);
  const magBA = vec3Magnitude(ba);
  const magBC = vec3Magnitude(bc);

  if (magBA === 0 || magBC === 0) return 0;

  const cosAngle = Math.max(-1, Math.min(1, dot / (magBA * magBC)));
  return Math.acos(cosAngle) * (180 / Math.PI);
}

// Calculate centroid of a set of points
function calculateCentroid(points: number[][]): number[] {
  const n = points.length;
  const centroid = [0, 0, 0];
  for (const p of points) {
    centroid[0] += p[0];
    centroid[1] += p[1];
    centroid[2] += p[2];
  }
  return centroid.map(c => c / n);
}

// Center points at origin
function centerPoints(points: number[][]): number[][] {
  const centroid = calculateCentroid(points);
  return points.map(p => [p[0] - centroid[0], p[1] - centroid[1], p[2] - centroid[2]]);
}

// Scale points to unit RMS distance from origin
function normalizeScale(points: number[][]): number[][] {
  let sumSq = 0;
  for (const p of points) {
    sumSq += p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
  }
  const rms = Math.sqrt(sumSq / points.length);
  if (rms === 0) return points;
  return points.map(p => [p[0] / rms, p[1] / rms, p[2] / rms]);
}

// Calculate CShM for a specific permutation
function calculateCShMForPermutation(
  observed: number[][],
  reference: number[][],
  permutation: number[]
): number {
  let denominator = 0;
  for (const p of observed) {
    denominator += p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
  }
  if (denominator === 0) return 100;

  let numerator = 0;
  for (let i = 0; i < observed.length; i++) {
    const obs = observed[i];
    const ref = reference[permutation[i]];
    const dx = obs[0] - ref[0];
    const dy = obs[1] - ref[1];
    const dz = obs[2] - ref[2];
    numerator += dx * dx + dy * dy + dz * dz;
  }

  return 100 * numerator / denominator;
}

// Simple rotation matrix application
function rotatePoint(p: number[], rotMatrix: number[][]): number[] {
  return [
    rotMatrix[0][0] * p[0] + rotMatrix[0][1] * p[1] + rotMatrix[0][2] * p[2],
    rotMatrix[1][0] * p[0] + rotMatrix[1][1] * p[1] + rotMatrix[1][2] * p[2],
    rotMatrix[2][0] * p[0] + rotMatrix[2][1] * p[1] + rotMatrix[2][2] * p[2]
  ];
}

// Generate rotation matrix from Euler angles
function eulerRotationMatrix(alpha: number, beta: number, gamma: number): number[][] {
  const ca = Math.cos(alpha), sa = Math.sin(alpha);
  const cb = Math.cos(beta), sb = Math.sin(beta);
  const cg = Math.cos(gamma), sg = Math.sin(gamma);
  return [
    [ca * cb * cg - sa * sg, -ca * cb * sg - sa * cg, ca * sb],
    [sa * cb * cg + ca * sg, -sa * cb * sg + ca * cg, sa * sb],
    [-sb * cg, sb * sg, cb]
  ];
}

// Calculate CShM with optimal rotation
function calculateCShMWithRotation(
  observed: number[][],
  reference: number[][],
  permutation: number[]
): number {
  let minCShM = Infinity;
  const nTrials = 50;

  for (let t = 0; t < nTrials; t++) {
    const alpha = (t === 0) ? 0 : Math.random() * 2 * Math.PI;
    const beta = (t === 0) ? 0 : Math.random() * Math.PI;
    const gamma = (t === 0) ? 0 : Math.random() * 2 * Math.PI;

    const rotMatrix = eulerRotationMatrix(alpha, beta, gamma);
    const rotatedRef = reference.map(p => rotatePoint(p, rotMatrix));

    const cshm = calculateCShMForPermutation(observed, rotatedRef, permutation);
    if (cshm < minCShM) {
      minCShM = cshm;
    }
  }

  return minCShM;
}

// Generate permutations (for small n) or sample permutations (for large n)
function generatePermutations(n: number, maxPerm: number = 5000): number[][] {
  if (n <= 6) {
    const result: number[][] = [];
    const permute = (arr: number[], start: number) => {
      if (start === arr.length) {
        result.push([...arr]);
        return;
      }
      for (let i = start; i < arr.length; i++) {
        [arr[start], arr[i]] = [arr[i], arr[start]];
        permute(arr, start + 1);
        [arr[start], arr[i]] = [arr[i], arr[start]];
      }
    };
    permute(Array.from({ length: n }, (_, i) => i), 0);
    return result;
  } else {
    const result: number[][] = [];
    const identity = Array.from({ length: n }, (_, i) => i);
    result.push([...identity]);

    for (let i = 0; i < maxPerm - 1; i++) {
      const perm = [...identity];
      for (let j = n - 1; j > 0; j--) {
        const k = Math.floor(Math.random() * (j + 1));
        [perm[j], perm[k]] = [perm[k], perm[j]];
      }
      result.push(perm);
    }
    return result;
  }
}

// Main CShM calculation function
function calculateCShM(
  observedPositions: [number, number, number][],
  referenceGeometry: { name: string; vertices: number[][] }
): { cshm: number; geometryName: string } {
  const n = observedPositions.length;
  if (n !== referenceGeometry.vertices.length) {
    return { cshm: 100, geometryName: referenceGeometry.name };
  }

  let observed = observedPositions.map(v => [v[0], v[1], v[2]]);
  let reference = referenceGeometry.vertices.map(v => [...v]);

  observed = normalizeScale(centerPoints(observed));
  reference = normalizeScale(centerPoints(reference));

  const permutations = generatePermutations(n, n <= 8 ? 5000 : 1000);
  let minCShM = Infinity;

  for (const perm of permutations) {
    const cshm = calculateCShMWithRotation(observed, reference, perm);
    if (cshm < minCShM) {
      minCShM = cshm;
    }
    if (minCShM < 0.1) break;
  }

  return { cshm: minCShM, geometryName: referenceGeometry.name };
}

// Analyze coordination geometry for a metal center using CShM
export function analyzeCoordinationGeometry(
  metalPos: [number, number, number],
  ligandPositions: { atom: string; pos: [number, number, number] }[]
): GeometryAnalysis | null {
  const cn = ligandPositions.length;

  if (cn < 2) {
    return {
      coordinationNumber: cn,
      geometryType: cn === 1 ? 'Monocoordinate' : 'None',
      idealGeometry: cn === 1 ? 'Terminal' : 'N/A',
      angles: [],
      avgAngle: 0,
      rmsd: 0,
      distortion: 'ideal'
    };
  }

  // Calculate all L-M-L angles
  const angles: { atom1: string; atom2: string; angle: number }[] = [];
  for (let i = 0; i < cn; i++) {
    for (let j = i + 1; j < cn; j++) {
      const angle = calculateAngle(ligandPositions[i].pos, metalPos, ligandPositions[j].pos);
      angles.push({
        atom1: ligandPositions[i].atom,
        atom2: ligandPositions[j].atom,
        angle: angle
      });
    }
  }

  const avgAngle = angles.reduce((sum, a) => sum + a.angle, 0) / angles.length;

  // Get reference geometries for this CN
  const referenceGeometries = REFERENCE_POLYHEDRA[cn];

  if (!referenceGeometries || referenceGeometries.length === 0) {
    return {
      coordinationNumber: cn,
      geometryType: `CN-${cn} (no reference)`,
      idealGeometry: `CN-${cn}`,
      angles: angles.sort((a, b) => a.angle - b.angle),
      avgAngle,
      rmsd: 0,
      distortion: 'moderate'
    };
  }

  // Calculate CShM for each reference geometry
  const ligandVectors = ligandPositions.map(lp => lp.pos);
  let bestGeometry = '';
  let bestCShM = Infinity;

  for (const refGeom of referenceGeometries) {
    const result = calculateCShM(ligandVectors, refGeom);
    if (result.cshm < bestCShM) {
      bestCShM = result.cshm;
      bestGeometry = result.geometryName;
    }
  }

  // Determine distortion level based on CShM value
  let distortion: 'ideal' | 'low' | 'moderate' | 'high' | 'severe';
  if (bestCShM < 0.5) distortion = 'ideal';
  else if (bestCShM < 1.0) distortion = 'low';
  else if (bestCShM < 3.0) distortion = 'moderate';
  else if (bestCShM < 5.0) distortion = 'high';
  else distortion = 'severe';

  return {
    coordinationNumber: cn,
    geometryType: bestGeometry,
    idealGeometry: bestGeometry,
    angles: angles.sort((a, b) => a.angle - b.angle),
    avgAngle,
    rmsd: bestCShM,
    distortion
  };
}

// Classify binding site as functional vs crystal artifact
export function classifyBindingSite(
  proteinContacts: number,
  waterContacts: number,
  geometry: GeometryAnalysis | null
): { type: 'functional' | 'crystal_artifact' | 'uncertain'; reason: string } {
  const totalContacts = proteinContacts + waterContacts;

  if (proteinContacts >= 3) {
    return {
      type: 'functional',
      reason: `${proteinContacts} protein contacts - likely functional binding site`
    };
  } else if (proteinContacts === 0 && waterContacts > 0) {
    return {
      type: 'crystal_artifact',
      reason: `Only water coordination (${waterContacts}) - likely crystal additive`
    };
  } else if (proteinContacts === 1 && waterContacts >= 3) {
    return {
      type: 'crystal_artifact',
      reason: `Minimal protein contact (1), mostly water (${waterContacts}) - likely adventitious`
    };
  } else if (proteinContacts >= 2 && geometry && geometry.distortion !== 'severe') {
    return {
      type: 'functional',
      reason: `${proteinContacts} protein contacts with ${geometry.geometryType} geometry`
    };
  } else if (proteinContacts === 2 && waterContacts >= 2) {
    return {
      type: 'uncertain',
      reason: `Mixed coordination (${proteinContacts} protein, ${waterContacts} water) - needs manual review`
    };
  } else if (totalContacts <= 2) {
    return {
      type: 'uncertain',
      reason: `Low coordination number (${totalContacts}) - incomplete or artifact`
    };
  } else {
    return {
      type: 'uncertain',
      reason: `${proteinContacts} protein, ${waterContacts} water contacts`
    };
  }
}

// Analyze hydration state
export function analyzeHydration(
  element: string,
  proteinContacts: number,
  waterContacts: number
): HydrationAnalysis {
  const expectedHydration = EXPECTED_HYDRATION_BY_METAL[element] || EXPECTED_HYDRATION_BY_METAL['DEFAULT'];
  const waterDisplacement = expectedHydration > 0
    ? Math.round(((expectedHydration - waterContacts) / expectedHydration) * 100)
    : 0;

  let hydrationState: 'fully_hydrated' | 'partially_hydrated' | 'dehydrated';
  let hydrationNote: string;

  if (proteinContacts === 0 && waterContacts > 0) {
    hydrationState = 'fully_hydrated';
    hydrationNote = `All ${waterContacts} coordination sites occupied by water - metal is solvated, not protein-bound`;
  } else if (waterContacts === 0 && proteinContacts > 0) {
    hydrationState = 'dehydrated';
    hydrationNote = `Complete water displacement by ${proteinContacts} protein ligands - mature binding site`;
  } else if (proteinContacts > 0 && waterContacts > 0) {
    hydrationState = 'partially_hydrated';
    const ratio = proteinContacts / (proteinContacts + waterContacts);
    if (ratio >= 0.7) {
      hydrationNote = `High protein occupancy (${Math.round(ratio * 100)}%) - ${Math.max(0, waterDisplacement)}% water displaced`;
    } else if (ratio >= 0.4) {
      hydrationNote = `Mixed coordination (${proteinContacts} protein, ${waterContacts} water) - binding site may be incomplete`;
    } else {
      hydrationNote = `Low protein occupancy (${Math.round(ratio * 100)}%) - mostly hydrated, weak binding`;
    }
  } else {
    hydrationState = 'dehydrated';
    hydrationNote = 'No coordination detected';
  }

  return {
    waterCount: waterContacts,
    proteinLigandCount: proteinContacts,
    hydrationState,
    expectedHydration,
    waterDisplacement: Math.max(0, waterDisplacement),
    hydrationNote
  };
}

/**
 * Parse PDB content and find metal ions with their coordination
 * This is a standalone function that works without Mol* structure
 */
export function findMetalCoordinationFromPDB(
  pdbContent: string,
  coordinationRadius: number = 3.0
): MetalCoordination[] {
  const lines = pdbContent.split('\n');

  // Parse all atoms
  const atoms: {
    serial: number;
    name: string;
    resName: string;
    chainId: string;
    resSeq: number;
    x: number;
    y: number;
    z: number;
    element: string;
    isHETATM: boolean;
  }[] = [];

  for (const line of lines) {
    if (line.startsWith('ATOM') || line.startsWith('HETATM')) {
      const isHETATM = line.startsWith('HETATM');
      try {
        atoms.push({
          serial: parseInt(line.slice(6, 11).trim()),
          name: line.slice(12, 16).trim(),
          resName: line.slice(17, 20).trim(),
          chainId: line.slice(21, 22).trim(),
          resSeq: parseInt(line.slice(22, 26).trim()),
          x: parseFloat(line.slice(30, 38).trim()),
          y: parseFloat(line.slice(38, 46).trim()),
          z: parseFloat(line.slice(46, 54).trim()),
          element: line.slice(76, 78).trim().toUpperCase() || line.slice(12, 14).trim().toUpperCase(),
          isHETATM
        });
      } catch {
        // Skip malformed lines
      }
    }
  }

  // Find metal ions
  const metals = atoms.filter(a => METAL_ELEMENTS.has(a.element));
  const results: MetalCoordination[] = [];

  for (const metal of metals) {
    const metalPos: [number, number, number] = [metal.x, metal.y, metal.z];
    const coordinating: CoordinatingAtom[] = [];

    // Find coordinating atoms
    for (const atom of atoms) {
      if (atom.serial === metal.serial) continue;

      const atomPos: [number, number, number] = [atom.x, atom.y, atom.z];
      const dist = vec3Distance(metalPos, atomPos);

      if (dist <= coordinationRadius) {
        const isWater = atom.resName === 'HOH' || atom.resName === 'WAT';
        coordinating.push({
          atom: atom.name,
          residue: `${atom.resName}${atom.resSeq}`,
          chain: atom.chainId,
          distance: dist,
          isWater,
          position: atomPos
        });
      }
    }

    // Sort by distance
    coordinating.sort((a, b) => a.distance - b.distance);

    // Analyze geometry
    const ligandPositions = coordinating.map(c => ({
      atom: `${c.atom}(${c.residue})`,
      pos: c.position
    }));
    const geometry = analyzeCoordinationGeometry(metalPos, ligandPositions);

    // Classify binding site
    const proteinContacts = coordinating.filter(c => !c.isWater).length;
    const waterContacts = coordinating.filter(c => c.isWater).length;
    const classification = classifyBindingSite(proteinContacts, waterContacts, geometry);

    // Hydration analysis
    const hydrationAnalysis = analyzeHydration(metal.element, proteinContacts, waterContacts);

    results.push({
      element: metal.element,
      info: `${metal.element} (${metal.resName}${metal.resSeq}, Chain ${metal.chainId})`,
      chainId: metal.chainId,
      resSeq: metal.resSeq,
      resName: metal.resName,
      pos: metalPos,
      coordinating,
      geometry,
      bindingSiteType: classification.type,
      bindingSiteReason: classification.reason,
      hydrationAnalysis
    });
  }

  return results;
}
