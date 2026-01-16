# Ligand-Dimer Enhancements V2 Port to Banta_Lab_RFdiffusion

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Port the v2 ligand analysis enhancements from Protein_engineering_tools to Banta_Lab_RFdiffusion's Next.js architecture.

**Architecture:** Banta_Lab_RFdiffusion uses Next.js 16 with:
- Server/Client component separation (ProteinViewer.tsx wrapper, ProteinViewerClient.tsx implementation)
- Zustand store for state management (store.ts)
- Tailwind CSS for styling
- Modular lib files (ligandAnalysis.ts, metalAnalysis.ts)

**Tech Stack:** Next.js 16, React 19, Zustand 5, Mol*, Tailwind CSS 4

---

## Context: What Exists

### In Banta_Lab_RFdiffusion (v1 - already present):
- `frontend/src/lib/ligandAnalysis.ts` (359 lines) - Basic contact analysis, interaction types, binding site classification
- `frontend/src/components/ProteinViewerClient.tsx` (701 lines) - Mol* viewer with focusOnLigand capability
- `frontend/src/lib/store.ts` (436 lines) - Zustand store with ligandData state

### What v2 Adds (from Protein_engineering_tools):
1. Enhanced pharmacophore extraction (donor/acceptor/aromatic/hydrophobic classification)
2. Shell analysis (primary/secondary shell residue classification)
3. Binder Quality Score (PLIP-inspired scoring algorithm)
4. Export functionality (JSON, CSV, RFdiffusion hotspot format)
5. 3D pharmacophore visualization (spheres in Mol* viewer)
6. Pi-stacking detection with proper ring center calculations

---

## Task 1: Enhance Pharmacophore Feature Extraction

**Files:**
- Modify: `frontend/src/lib/ligandAnalysis.ts`

**Step 1: Write the failing test**

Create test file first:

```typescript
// frontend/src/lib/__tests__/ligandAnalysis.test.ts
import { describe, it, expect } from 'vitest';
import { extractPharmacophoreFeatures, PharmacophoreFeature } from '../ligandAnalysis';

describe('extractPharmacophoreFeatures', () => {
  it('should identify H-bond donors from NH and OH groups', () => {
    const contacts = [
      { residue: 'SER45', chain: 'A', atom: 'OG', distance: 2.8, interactionType: 'hydrogen_bond' as const },
      { residue: 'LYS72', chain: 'A', atom: 'NZ', distance: 3.1, interactionType: 'hydrogen_bond' as const },
    ];
    const features = extractPharmacophoreFeatures(contacts);

    const donors = features.filter(f => f.type === 'donor');
    expect(donors.length).toBeGreaterThan(0);
  });

  it('should identify aromatic features from PHE, TYR, TRP', () => {
    const contacts = [
      { residue: 'PHE123', chain: 'A', atom: 'CZ', distance: 3.8, interactionType: 'pi_stacking' as const },
    ];
    const features = extractPharmacophoreFeatures(contacts);

    const aromatics = features.filter(f => f.type === 'aromatic');
    expect(aromatics.length).toBe(1);
    expect(aromatics[0].residue).toBe('PHE123');
  });

  it('should identify hydrophobic features from carbon contacts', () => {
    const contacts = [
      { residue: 'LEU45', chain: 'A', atom: 'CD1', distance: 3.5, interactionType: 'hydrophobic' as const },
    ];
    const features = extractPharmacophoreFeatures(contacts);

    const hydrophobic = features.filter(f => f.type === 'hydrophobic');
    expect(hydrophobic.length).toBe(1);
  });

  it('should identify charged features from salt bridges', () => {
    const contacts = [
      { residue: 'ARG89', chain: 'A', atom: 'NH1', distance: 3.2, interactionType: 'salt_bridge' as const },
    ];
    const features = extractPharmacophoreFeatures(contacts);

    const positive = features.filter(f => f.type === 'positive');
    expect(positive.length).toBe(1);
  });
});
```

**Step 2: Run test to verify it fails**

Run: `cd frontend && npm test -- --run src/lib/__tests__/ligandAnalysis.test.ts`
Expected: FAIL with "extractPharmacophoreFeatures is not exported"

**Step 3: Add PharmacophoreFeature type and extraction function**

Add to `frontend/src/lib/ligandAnalysis.ts`:

```typescript
// Pharmacophore feature types
export type PharmacophoreType = 'donor' | 'acceptor' | 'aromatic' | 'hydrophobic' | 'positive' | 'negative';

export interface PharmacophoreFeature {
  type: PharmacophoreType;
  residue: string;
  chain: string;
  atom: string;
  distance: number;
  // Position for 3D rendering (optional, populated during visualization)
  position?: [number, number, number];
}

// Atom classification for pharmacophore features
const DONOR_ATOMS = new Set(['N', 'O']); // When bonded to H
const ACCEPTOR_ATOMS = new Set(['O', 'N', 'S', 'F']);
const AROMATIC_RESIDUES = new Set(['PHE', 'TYR', 'TRP', 'HIS']);

/**
 * Extract pharmacophore features from ligand contacts
 * Classifies each contact into donor, acceptor, aromatic, hydrophobic, or charged
 */
export function extractPharmacophoreFeatures(contacts: LigandContact[]): PharmacophoreFeature[] {
  const features: PharmacophoreFeature[] = [];

  for (const contact of contacts) {
    const resName = contact.residue.replace(/[0-9]/g, ''); // Extract residue name
    const atomElement = contact.atom.charAt(0); // First char is usually element

    let type: PharmacophoreType;

    // Determine feature type based on interaction and atom
    if (contact.interactionType === 'pi_stacking' || AROMATIC_RESIDUES.has(resName)) {
      type = 'aromatic';
    } else if (contact.interactionType === 'salt_bridge') {
      type = POSITIVE_RESIDUES.has(resName) ? 'positive' : 'negative';
    } else if (contact.interactionType === 'hydrogen_bond') {
      // Donors have N or O with attached H; acceptors have lone pairs
      // Simplification: N atoms are donors, O atoms can be either
      type = atomElement === 'N' ? 'donor' : 'acceptor';
    } else if (contact.interactionType === 'hydrophobic') {
      type = 'hydrophobic';
    } else {
      // Default: classify by atom type
      if (DONOR_ATOMS.has(atomElement)) {
        type = 'acceptor'; // Conservative: lone pair acceptors
      } else {
        type = 'hydrophobic';
      }
    }

    features.push({
      type,
      residue: contact.residue,
      chain: contact.chain,
      atom: contact.atom,
      distance: contact.distance,
    });
  }

  return features;
}
```

**Step 4: Run test to verify it passes**

Run: `cd frontend && npm test -- --run src/lib/__tests__/ligandAnalysis.test.ts`
Expected: PASS (4 tests)

**Step 5: Commit**

```bash
git add frontend/src/lib/ligandAnalysis.ts frontend/src/lib/__tests__/ligandAnalysis.test.ts
git commit -m "feat(ligand): add pharmacophore feature extraction

- Add PharmacophoreFeature type and extractPharmacophoreFeatures function
- Classify contacts as donor/acceptor/aromatic/hydrophobic/charged
- Add comprehensive test suite"
```

---

## Task 2: Add Shell Analysis (Primary/Secondary Classification)

**Files:**
- Modify: `frontend/src/lib/ligandAnalysis.ts`

**Step 1: Write the failing test**

Add to test file:

```typescript
describe('classifyShellResidues', () => {
  it('should classify residues within 4A as primary shell', () => {
    const contacts = [
      { residue: 'ASP45', chain: 'A', atom: 'OD1', distance: 2.8, interactionType: 'hydrogen_bond' as const },
      { residue: 'LEU72', chain: 'A', atom: 'CD1', distance: 3.5, interactionType: 'hydrophobic' as const },
    ];
    const shells = classifyShellResidues(contacts);

    expect(shells.primary.length).toBe(2);
    expect(shells.secondary.length).toBe(0);
  });

  it('should classify residues 4-6A as secondary shell', () => {
    const contacts = [
      { residue: 'PHE89', chain: 'A', atom: 'CZ', distance: 5.2, interactionType: 'hydrophobic' as const },
    ];
    const shells = classifyShellResidues(contacts);

    expect(shells.primary.length).toBe(0);
    expect(shells.secondary.length).toBe(1);
  });

  it('should provide shell statistics', () => {
    const contacts = [
      { residue: 'ASP45', chain: 'A', atom: 'OD1', distance: 2.8, interactionType: 'hydrogen_bond' as const },
      { residue: 'LEU72', chain: 'A', atom: 'CD1', distance: 3.5, interactionType: 'hydrophobic' as const },
      { residue: 'PHE89', chain: 'A', atom: 'CZ', distance: 5.2, interactionType: 'hydrophobic' as const },
    ];
    const shells = classifyShellResidues(contacts);

    expect(shells.stats.primaryCount).toBe(2);
    expect(shells.stats.secondaryCount).toBe(1);
    expect(shells.stats.totalCount).toBe(3);
  });
});
```

**Step 2: Run test to verify it fails**

Run: `cd frontend && npm test -- --run src/lib/__tests__/ligandAnalysis.test.ts`
Expected: FAIL with "classifyShellResidues is not exported"

**Step 3: Implement shell classification**

Add to `frontend/src/lib/ligandAnalysis.ts`:

```typescript
// Shell classification thresholds
const PRIMARY_SHELL_CUTOFF = 4.0;  // Direct contacts
const SECONDARY_SHELL_CUTOFF = 6.0; // Extended environment

export interface ShellAnalysis {
  primary: LigandContact[];   // < 4A: direct binding
  secondary: LigandContact[]; // 4-6A: extended shell
  stats: {
    primaryCount: number;
    secondaryCount: number;
    totalCount: number;
    avgPrimaryDistance: number;
    avgSecondaryDistance: number;
  };
}

/**
 * Classify contacts into primary (direct) and secondary (extended) shells
 * Primary shell: < 4A - direct binding interactions
 * Secondary shell: 4-6A - supporting/stabilizing residues
 */
export function classifyShellResidues(contacts: LigandContact[]): ShellAnalysis {
  const primary: LigandContact[] = [];
  const secondary: LigandContact[] = [];

  for (const contact of contacts) {
    if (contact.distance < PRIMARY_SHELL_CUTOFF) {
      primary.push(contact);
    } else if (contact.distance < SECONDARY_SHELL_CUTOFF) {
      secondary.push(contact);
    }
    // Beyond 6A: not included in shell analysis
  }

  // Calculate averages
  const avgPrimary = primary.length > 0
    ? primary.reduce((sum, c) => sum + c.distance, 0) / primary.length
    : 0;
  const avgSecondary = secondary.length > 0
    ? secondary.reduce((sum, c) => sum + c.distance, 0) / secondary.length
    : 0;

  return {
    primary,
    secondary,
    stats: {
      primaryCount: primary.length,
      secondaryCount: secondary.length,
      totalCount: primary.length + secondary.length,
      avgPrimaryDistance: avgPrimary,
      avgSecondaryDistance: avgSecondary,
    },
  };
}
```

**Step 4: Run test to verify it passes**

Run: `cd frontend && npm test -- --run src/lib/__tests__/ligandAnalysis.test.ts`
Expected: PASS

**Step 5: Commit**

```bash
git add frontend/src/lib/ligandAnalysis.ts frontend/src/lib/__tests__/ligandAnalysis.test.ts
git commit -m "feat(ligand): add shell analysis for binding site characterization

- Add ShellAnalysis interface and classifyShellResidues function
- Primary shell: direct contacts within 4A
- Secondary shell: supporting residues 4-6A
- Include statistics (counts, average distances)"
```

---

## Task 3: Add Binder Quality Score

**Files:**
- Modify: `frontend/src/lib/ligandAnalysis.ts`

**Step 1: Write the failing test**

Add to test file:

```typescript
describe('calculateBinderQualityScore', () => {
  it('should return high score for well-formed binding site', () => {
    const contacts = [
      { residue: 'ASP45', chain: 'A', atom: 'OD1', distance: 2.8, interactionType: 'hydrogen_bond' as const },
      { residue: 'ASN48', chain: 'A', atom: 'ND2', distance: 3.0, interactionType: 'hydrogen_bond' as const },
      { residue: 'LEU72', chain: 'A', atom: 'CD1', distance: 3.5, interactionType: 'hydrophobic' as const },
      { residue: 'PHE89', chain: 'A', atom: 'CZ', distance: 3.8, interactionType: 'pi_stacking' as const },
      { residue: 'ARG92', chain: 'A', atom: 'NH1', distance: 3.2, interactionType: 'salt_bridge' as const },
    ];
    const result = calculateBinderQualityScore(contacts);

    expect(result.score).toBeGreaterThan(70);
    expect(result.rating).toBe('EXCELLENT');
  });

  it('should return low score for sparse contacts', () => {
    const contacts = [
      { residue: 'ALA45', chain: 'A', atom: 'CB', distance: 3.9, interactionType: 'hydrophobic' as const },
    ];
    const result = calculateBinderQualityScore(contacts);

    expect(result.score).toBeLessThan(40);
    expect(result.rating).toBe('POOR');
  });

  it('should provide detailed breakdown', () => {
    const contacts = [
      { residue: 'ASP45', chain: 'A', atom: 'OD1', distance: 2.8, interactionType: 'hydrogen_bond' as const },
    ];
    const result = calculateBinderQualityScore(contacts);

    expect(result.breakdown).toBeDefined();
    expect(result.breakdown.hbondScore).toBeDefined();
    expect(result.breakdown.hydrophobicScore).toBeDefined();
  });
});
```

**Step 2: Run test to verify it fails**

Run: `cd frontend && npm test -- --run src/lib/__tests__/ligandAnalysis.test.ts`
Expected: FAIL

**Step 3: Implement Binder Quality Score**

Add to `frontend/src/lib/ligandAnalysis.ts`:

```typescript
export type QualityRating = 'EXCELLENT' | 'GOOD' | 'MODERATE' | 'POOR';

export interface BinderQualityResult {
  score: number;      // 0-100
  rating: QualityRating;
  breakdown: {
    hbondScore: number;       // H-bond contribution (0-30)
    hydrophobicScore: number; // Hydrophobic contribution (0-25)
    piStackScore: number;     // Pi-stacking contribution (0-15)
    saltBridgeScore: number;  // Salt bridge contribution (0-15)
    burialnessScore: number;  // Contact density (0-15)
  };
  suggestions: string[];
}

/**
 * Calculate a PLIP-inspired binder quality score
 * Evaluates the binding site based on interaction diversity and strength
 */
export function calculateBinderQualityScore(contacts: LigandContact[]): BinderQualityResult {
  // Count interactions by type
  const hbonds = contacts.filter(c => c.interactionType === 'hydrogen_bond').length;
  const hydrophobic = contacts.filter(c => c.interactionType === 'hydrophobic').length;
  const piStack = contacts.filter(c => c.interactionType === 'pi_stacking').length;
  const saltBridge = contacts.filter(c => c.interactionType === 'salt_bridge').length;

  // Score each category (with diminishing returns)
  const hbondScore = Math.min(30, hbonds * 10);          // Max 30 pts, 3+ H-bonds
  const hydrophobicScore = Math.min(25, hydrophobic * 5); // Max 25 pts, 5+ contacts
  const piStackScore = Math.min(15, piStack * 7.5);       // Max 15 pts, 2+ stacks
  const saltBridgeScore = Math.min(15, saltBridge * 7.5); // Max 15 pts, 2+ bridges

  // Burialness: contact density bonus
  const totalContacts = contacts.length;
  const burialnessScore = Math.min(15, totalContacts * 1.5); // Max 15 pts, 10+ contacts

  // Total score
  const score = Math.round(
    hbondScore + hydrophobicScore + piStackScore + saltBridgeScore + burialnessScore
  );

  // Rating
  let rating: QualityRating;
  if (score >= 70) rating = 'EXCELLENT';
  else if (score >= 50) rating = 'GOOD';
  else if (score >= 30) rating = 'MODERATE';
  else rating = 'POOR';

  // Suggestions for improvement
  const suggestions: string[] = [];
  if (hbonds < 2) suggestions.push('Add H-bond capable residues (Ser, Thr, Asn, Gln)');
  if (hydrophobic < 3) suggestions.push('Increase hydrophobic contacts (Leu, Ile, Val, Phe)');
  if (piStack === 0) suggestions.push('Consider aromatic residues for pi-stacking (Phe, Tyr, Trp)');
  if (saltBridge === 0 && totalContacts > 3) suggestions.push('Salt bridges could enhance specificity (Arg, Lys, Asp, Glu)');

  return {
    score,
    rating,
    breakdown: {
      hbondScore,
      hydrophobicScore,
      piStackScore,
      saltBridgeScore,
      burialnessScore,
    },
    suggestions,
  };
}
```

**Step 4: Run test to verify it passes**

Run: `cd frontend && npm test -- --run src/lib/__tests__/ligandAnalysis.test.ts`
Expected: PASS

**Step 5: Commit**

```bash
git add frontend/src/lib/ligandAnalysis.ts frontend/src/lib/__tests__/ligandAnalysis.test.ts
git commit -m "feat(ligand): add PLIP-inspired binder quality score

- Score 0-100 based on interaction diversity
- Categories: H-bonds (30), hydrophobic (25), pi-stacking (15), salt bridges (15), burial (15)
- Quality ratings: EXCELLENT/GOOD/MODERATE/POOR
- Actionable suggestions for improvement"
```

---

## Task 4: Add Export Functionality

**Files:**
- Create: `frontend/src/lib/ligandExport.ts`

**Step 1: Write the failing test**

```typescript
// frontend/src/lib/__tests__/ligandExport.test.ts
import { describe, it, expect } from 'vitest';
import { exportToJSON, exportToCSV, exportToRFdiffusionHotspots } from '../ligandExport';

const mockLigandData = {
  ligandCount: 1,
  ligandDetails: [{
    name: 'ATP',
    info: 'ATP (Chain A, 501)',
    atoms: 31,
    chainId: 'A',
    resSeq: 501,
    contacts: [
      { residue: 'ASP45', chain: 'A', atom: 'OD1', distance: 2.8, interactionType: 'hydrogen_bond' as const },
      { residue: 'LEU72', chain: 'A', atom: 'CD1', distance: 3.5, interactionType: 'hydrophobic' as const },
    ],
    proteinContactCount: 5,
    waterContactCount: 2,
    bindingSiteType: 'functional' as const,
    bindingSiteReason: '5 protein residue contacts',
  }],
};

describe('exportToJSON', () => {
  it('should export complete analysis as JSON', () => {
    const json = exportToJSON(mockLigandData);
    const parsed = JSON.parse(json);

    expect(parsed.ligandCount).toBe(1);
    expect(parsed.ligandDetails[0].name).toBe('ATP');
  });
});

describe('exportToCSV', () => {
  it('should export contacts as CSV', () => {
    const csv = exportToCSV(mockLigandData);

    expect(csv).toContain('Ligand,Chain,ResSeq,ContactResidue');
    expect(csv).toContain('ATP,A,501,ASP45');
  });
});

describe('exportToRFdiffusionHotspots', () => {
  it('should generate RFdiffusion hotspot config', () => {
    const config = exportToRFdiffusionHotspots(mockLigandData);

    expect(config).toContain('hotspot_residues');
    expect(config).toContain('A45'); // ASP45 on chain A
  });
});
```

**Step 2: Run test to verify it fails**

Run: `cd frontend && npm test -- --run src/lib/__tests__/ligandExport.test.ts`
Expected: FAIL - module not found

**Step 3: Implement export functions**

Create `frontend/src/lib/ligandExport.ts`:

```typescript
/**
 * Export utilities for ligand binding analysis
 * Supports JSON, CSV, and RFdiffusion hotspot format
 */

import type { LigandAnalysisResult, LigandData, LigandContact } from './ligandAnalysis';

/**
 * Export complete analysis as JSON
 */
export function exportToJSON(data: LigandAnalysisResult): string {
  return JSON.stringify(data, null, 2);
}

/**
 * Export contacts as CSV for spreadsheet analysis
 */
export function exportToCSV(data: LigandAnalysisResult): string {
  const headers = [
    'Ligand', 'Chain', 'ResSeq', 'ContactResidue', 'ContactChain',
    'ContactAtom', 'Distance', 'InteractionType', 'BindingSiteType'
  ];

  const rows: string[][] = [headers];

  for (const ligand of data.ligandDetails) {
    for (const contact of ligand.contacts) {
      rows.push([
        ligand.name,
        ligand.chainId,
        String(ligand.resSeq),
        contact.residue,
        contact.chain,
        contact.atom,
        contact.distance.toFixed(2),
        contact.interactionType,
        ligand.bindingSiteType,
      ]);
    }
  }

  return rows.map(row => row.join(',')).join('\n');
}

/**
 * Export as RFdiffusion hotspot configuration
 * Generates hotspot_residues format for binder design
 */
export function exportToRFdiffusionHotspots(data: LigandAnalysisResult): string {
  const hotspots: string[] = [];

  for (const ligand of data.ligandDetails) {
    // Only include functional binding sites
    if (ligand.bindingSiteType !== 'functional') continue;

    // Extract unique residue positions
    const seenResidues = new Set<string>();
    for (const contact of ligand.contacts) {
      // Extract residue number from "ASP45" format
      const resNum = contact.residue.replace(/[A-Z]/g, '');
      const hotspot = `${contact.chain}${resNum}`;

      if (!seenResidues.has(hotspot)) {
        seenResidues.add(hotspot);
        hotspots.push(hotspot);
      }
    }
  }

  // Format as RFdiffusion config
  return `# RFdiffusion Hotspot Configuration
# Generated from ligand binding analysis
# Use with: --hotspot_residues

hotspot_residues: [${hotspots.join(', ')}]

# Example usage:
# python run_inference.py \\
#   --config-name base \\
#   contigmap.contigs="[A1-100/0 B1-50]" \\
#   ppi.hotspot_residues=[${hotspots.join(',')}]
`;
}

/**
 * Trigger browser download of exported data
 */
export function downloadFile(content: string, filename: string, mimeType: string): void {
  const blob = new Blob([content], { type: mimeType });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
}
```

**Step 4: Run test to verify it passes**

Run: `cd frontend && npm test -- --run src/lib/__tests__/ligandExport.test.ts`
Expected: PASS

**Step 5: Commit**

```bash
git add frontend/src/lib/ligandExport.ts frontend/src/lib/__tests__/ligandExport.test.ts
git commit -m "feat(ligand): add export functionality for binding analysis

- JSON export for complete analysis data
- CSV export for spreadsheet analysis
- RFdiffusion hotspot config generation
- Browser download utility"
```

---

## Task 5: Add LigandBindingPanel UI Component

**Files:**
- Create: `frontend/src/components/LigandBindingPanel.tsx`

**Step 1: Create the UI component**

This is a UI task - no failing test needed. Create `frontend/src/components/LigandBindingPanel.tsx`:

```typescript
'use client';

import { useState, useMemo } from 'react';
import { useStore } from '@/lib/store';
import {
  extractPharmacophoreFeatures,
  classifyShellResidues,
  calculateBinderQualityScore,
  getInteractionColor,
  getInteractionLabel,
  type PharmacophoreFeature,
  type ShellAnalysis,
  type BinderQualityResult,
} from '@/lib/ligandAnalysis';
import { exportToJSON, exportToCSV, exportToRFdiffusionHotspots, downloadFile } from '@/lib/ligandExport';

interface LigandBindingPanelProps {
  onToggle3D?: (show: boolean) => void;
  show3D?: boolean;
}

export function LigandBindingPanel({ onToggle3D, show3D = false }: LigandBindingPanelProps) {
  const { ligandData, focusedLigandIndex, setFocusedLigandIndex } = useStore();
  const [expandedLigand, setExpandedLigand] = useState<number | null>(null);

  // Compute analysis for focused ligand
  const analysis = useMemo(() => {
    if (!ligandData || focusedLigandIndex === null) return null;

    const ligand = ligandData.ligandDetails[focusedLigandIndex];
    if (!ligand) return null;

    const pharmacophores = extractPharmacophoreFeatures(ligand.contacts);
    const shells = classifyShellResidues(ligand.contacts);
    const quality = calculateBinderQualityScore(ligand.contacts);

    return { ligand, pharmacophores, shells, quality };
  }, [ligandData, focusedLigandIndex]);

  if (!ligandData || ligandData.ligandCount === 0) {
    return (
      <div className="p-4 text-center text-gray-500">
        <p>No ligands detected in structure</p>
      </div>
    );
  }

  const handleExport = (format: 'json' | 'csv' | 'rfdiffusion') => {
    if (!ligandData) return;

    switch (format) {
      case 'json':
        downloadFile(exportToJSON(ligandData), 'ligand-analysis.json', 'application/json');
        break;
      case 'csv':
        downloadFile(exportToCSV(ligandData), 'ligand-contacts.csv', 'text/csv');
        break;
      case 'rfdiffusion':
        downloadFile(exportToRFdiffusionHotspots(ligandData), 'hotspots.txt', 'text/plain');
        break;
    }
  };

  const getRatingColor = (rating: string) => {
    switch (rating) {
      case 'EXCELLENT': return 'bg-green-100 text-green-800 border-green-200';
      case 'GOOD': return 'bg-blue-100 text-blue-800 border-blue-200';
      case 'MODERATE': return 'bg-yellow-100 text-yellow-800 border-yellow-200';
      case 'POOR': return 'bg-red-100 text-red-800 border-red-200';
      default: return 'bg-gray-100 text-gray-800 border-gray-200';
    }
  };

  return (
    <div className="space-y-4">
      {/* Header with Export Buttons */}
      <div className="flex items-center justify-between">
        <h3 className="font-semibold text-gray-900">
          Ligand Binding Analysis
          <span className="ml-2 text-sm font-normal text-gray-500">
            ({ligandData.ligandCount} ligand{ligandData.ligandCount > 1 ? 's' : ''})
          </span>
        </h3>
        <div className="flex gap-2">
          <button
            onClick={() => handleExport('json')}
            className="px-2 py-1 text-xs bg-gray-100 hover:bg-gray-200 rounded"
            title="Export as JSON"
          >
            JSON
          </button>
          <button
            onClick={() => handleExport('csv')}
            className="px-2 py-1 text-xs bg-gray-100 hover:bg-gray-200 rounded"
            title="Export as CSV"
          >
            CSV
          </button>
          <button
            onClick={() => handleExport('rfdiffusion')}
            className="px-2 py-1 text-xs bg-purple-100 hover:bg-purple-200 text-purple-700 rounded"
            title="Export RFdiffusion hotspots"
          >
            RFdiffusion
          </button>
        </div>
      </div>

      {/* Ligand List */}
      <div className="space-y-2">
        {ligandData.ligandDetails.map((ligand, index) => (
          <div
            key={`${ligand.chainId}-${ligand.resSeq}`}
            className={`border rounded-lg overflow-hidden ${
              focusedLigandIndex === index ? 'border-purple-500 ring-2 ring-purple-200' : 'border-gray-200'
            }`}
          >
            {/* Ligand Header */}
            <button
              onClick={() => {
                setFocusedLigandIndex(focusedLigandIndex === index ? null : index);
                setExpandedLigand(expandedLigand === index ? null : index);
              }}
              className="w-full px-4 py-3 flex items-center justify-between bg-white hover:bg-gray-50"
            >
              <div className="flex items-center gap-3">
                <span className={`px-2 py-0.5 text-xs rounded ${
                  ligand.bindingSiteType === 'functional'
                    ? 'bg-green-100 text-green-700'
                    : ligand.bindingSiteType === 'crystal_artifact'
                    ? 'bg-gray-100 text-gray-600'
                    : 'bg-yellow-100 text-yellow-700'
                }`}>
                  {ligand.bindingSiteType}
                </span>
                <span className="font-medium">{ligand.name}</span>
                <span className="text-sm text-gray-500">
                  Chain {ligand.chainId}:{ligand.resSeq}
                </span>
              </div>
              <span className="text-sm text-gray-500">
                {ligand.contacts.length} contacts
              </span>
            </button>

            {/* Expanded Analysis */}
            {expandedLigand === index && analysis && focusedLigandIndex === index && (
              <div className="px-4 pb-4 space-y-4 bg-gray-50">
                {/* Binder Quality Score */}
                <div className="p-3 bg-white rounded-lg border">
                  <div className="flex items-center justify-between mb-2">
                    <span className="text-sm font-medium">Binder Quality Score</span>
                    <span className={`px-2 py-0.5 text-xs font-bold rounded border ${getRatingColor(analysis.quality.rating)}`}>
                      {analysis.quality.score}/100 - {analysis.quality.rating}
                    </span>
                  </div>

                  {/* Score Breakdown */}
                  <div className="grid grid-cols-5 gap-2 text-xs">
                    <div className="text-center p-1 bg-blue-50 rounded">
                      <div className="font-bold text-blue-700">{analysis.quality.breakdown.hbondScore}</div>
                      <div className="text-blue-500">H-bonds</div>
                    </div>
                    <div className="text-center p-1 bg-green-50 rounded">
                      <div className="font-bold text-green-700">{analysis.quality.breakdown.hydrophobicScore}</div>
                      <div className="text-green-500">Hydrophobic</div>
                    </div>
                    <div className="text-center p-1 bg-purple-50 rounded">
                      <div className="font-bold text-purple-700">{analysis.quality.breakdown.piStackScore}</div>
                      <div className="text-purple-500">Pi-stack</div>
                    </div>
                    <div className="text-center p-1 bg-red-50 rounded">
                      <div className="font-bold text-red-700">{analysis.quality.breakdown.saltBridgeScore}</div>
                      <div className="text-red-500">Salt bridge</div>
                    </div>
                    <div className="text-center p-1 bg-gray-100 rounded">
                      <div className="font-bold text-gray-700">{analysis.quality.breakdown.burialnessScore}</div>
                      <div className="text-gray-500">Burial</div>
                    </div>
                  </div>
                </div>

                {/* Shell Analysis */}
                <div className="p-3 bg-white rounded-lg border">
                  <span className="text-sm font-medium">Shell Analysis</span>
                  <div className="mt-2 grid grid-cols-2 gap-2 text-xs">
                    <div className="p-2 bg-purple-50 rounded">
                      <div className="font-bold text-purple-700">{analysis.shells.stats.primaryCount}</div>
                      <div className="text-purple-500">Primary Shell (&lt;4Å)</div>
                    </div>
                    <div className="p-2 bg-blue-50 rounded">
                      <div className="font-bold text-blue-700">{analysis.shells.stats.secondaryCount}</div>
                      <div className="text-blue-500">Secondary Shell (4-6Å)</div>
                    </div>
                  </div>
                </div>

                {/* Pharmacophore Features */}
                <div className="p-3 bg-white rounded-lg border">
                  <div className="flex items-center justify-between mb-2">
                    <span className="text-sm font-medium">
                      Pharmacophore Features ({analysis.pharmacophores.length})
                    </span>
                    {onToggle3D && (
                      <button
                        onClick={() => onToggle3D(!show3D)}
                        className={`px-2 py-1 text-xs rounded ${
                          show3D
                            ? 'bg-purple-600 text-white'
                            : 'bg-purple-100 text-purple-700 hover:bg-purple-200'
                        }`}
                      >
                        {show3D ? 'Hide 3D' : 'Show 3D'}
                      </button>
                    )}
                  </div>

                  {/* Feature counts by type */}
                  <div className="flex flex-wrap gap-1">
                    {['donor', 'acceptor', 'aromatic', 'hydrophobic', 'positive', 'negative'].map(type => {
                      const count = analysis.pharmacophores.filter(f => f.type === type).length;
                      if (count === 0) return null;
                      return (
                        <span key={type} className="px-2 py-0.5 text-xs bg-gray-100 rounded">
                          {type}: {count}
                        </span>
                      );
                    })}
                  </div>
                </div>

                {/* Suggestions */}
                {analysis.quality.suggestions.length > 0 && (
                  <div className="p-3 bg-amber-50 rounded-lg border border-amber-200">
                    <span className="text-sm font-medium text-amber-800">Suggestions</span>
                    <ul className="mt-1 text-xs text-amber-700 space-y-1">
                      {analysis.quality.suggestions.map((s, i) => (
                        <li key={i}>• {s}</li>
                      ))}
                    </ul>
                  </div>
                )}
              </div>
            )}
          </div>
        ))}
      </div>
    </div>
  );
}

export default LigandBindingPanel;
```

**Step 2: Commit**

```bash
git add frontend/src/components/LigandBindingPanel.tsx
git commit -m "feat(ui): add LigandBindingPanel component

- Display ligand list with binding site classification
- Show binder quality score with breakdown
- Shell analysis visualization (primary/secondary)
- Pharmacophore feature summary
- Export buttons (JSON, CSV, RFdiffusion)
- 3D toggle for pharmacophore visualization"
```

---

## Task 6: Add 3D Pharmacophore Visualization to ProteinViewerClient

**Files:**
- Modify: `frontend/src/components/ProteinViewerClient.tsx`

**Step 1: Add pharmacophore rendering function**

This integrates with Mol* to render colored spheres for pharmacophore features. Add to `ProteinViewerClient.tsx`:

```typescript
// Add to imports at top
import type { PharmacophoreFeature } from '@/lib/ligandAnalysis';

// Add new prop to interface
interface ProteinViewerClientProps {
  // ... existing props
  pharmacophoreFeatures?: PharmacophoreFeature[];
  showPharmacophores?: boolean;
}

// Add pharmacophore colors
const PHARMACOPHORE_COLORS: Record<string, number> = {
  donor: 0x3B82F6,      // Blue
  acceptor: 0xEF4444,   // Red
  aromatic: 0xA855F7,   // Purple
  hydrophobic: 0x22C55E, // Green
  positive: 0x0EA5E9,   // Cyan
  negative: 0xF97316,   // Orange
};

// Add renderPharmacophores function inside component
const renderPharmacophores = useCallback(async (features: PharmacophoreFeature[]) => {
  if (!globalPlugin || !globalStructureRef) return;

  await loadMolstarModules();

  try {
    // Get current structure data
    const state = globalPlugin.state.data;
    const structureCell = state.cells.get(globalStructureRef);
    if (!structureCell?.obj?.data) return;

    const structure = structureCell.obj.data;

    // For each pharmacophore feature, create a sphere representation
    for (let i = 0; i < features.length; i++) {
      const feature = features[i];
      const resNum = parseInt(feature.residue.replace(/[A-Z]/g, ''));

      // Create selection for this residue's atom
      const expression = MS.struct.generator.atomGroups({
        'residue-test': MS.core.logic.and([
          MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_seq_id(), resNum]),
          MS.core.rel.eq([MS.struct.atomProperty.macromolecular.auth_asym_id(), feature.chain]),
        ]),
        'atom-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_atom_id(), feature.atom]),
      });

      const comp = await globalPlugin.builders.structure.tryCreateComponentFromExpression(
        globalStructureRef,
        expression,
        `pharmacophore-${feature.type}-${i}`
      );

      if (comp) {
        const color = PHARMACOPHORE_COLORS[feature.type] || 0x888888;
        await globalPlugin.builders.structure.representation.addRepresentation(comp, {
          type: 'spacefill',
          color: 'uniform',
          colorParams: { value: Color(color) },
          typeParams: { sizeFactor: 0.8, alpha: 0.7 },
        });
      }
    }

    console.log(`[ProteinViewer] Rendered ${features.length} pharmacophore features`);
  } catch (err) {
    console.error('[ProteinViewer] Failed to render pharmacophores:', err);
  }
}, []);

// Add effect to handle pharmacophore prop changes
useEffect(() => {
  if (showPharmacophores && pharmacophoreFeatures && pharmacophoreFeatures.length > 0) {
    renderPharmacophores(pharmacophoreFeatures);
  }
}, [showPharmacophores, pharmacophoreFeatures, renderPharmacophores]);
```

**Step 2: Commit**

```bash
git add frontend/src/components/ProteinViewerClient.tsx
git commit -m "feat(viewer): add 3D pharmacophore visualization

- Render pharmacophore features as colored spheres
- Color coding: donor (blue), acceptor (red), aromatic (purple), etc.
- Toggle via showPharmacophores prop
- Semi-transparent spheres (alpha 0.7) for clarity"
```

---

## Task 7: Update Zustand Store for Pharmacophore State

**Files:**
- Modify: `frontend/src/lib/store.ts`

**Step 1: Add pharmacophore state to store**

Add to store interface and implementation:

```typescript
// Add to imports
import type { PharmacophoreFeature } from './ligandAnalysis';

// Add to AppState interface
pharmacophoreFeatures: PharmacophoreFeature[] | null;
setPharmacophoreFeatures: (features: PharmacophoreFeature[] | null) => void;
showPharmacophores3D: boolean;
setShowPharmacophores3D: (show: boolean) => void;

// Add to store implementation
pharmacophoreFeatures: null,
setPharmacophoreFeatures: (features) => set({ pharmacophoreFeatures: features }),
showPharmacophores3D: false,
setShowPharmacophores3D: (show) => set({ showPharmacophores3D: show }),
```

**Step 2: Commit**

```bash
git add frontend/src/lib/store.ts
git commit -m "feat(store): add pharmacophore state management

- pharmacophoreFeatures: array of extracted features
- showPharmacophores3D: toggle for 3D visualization"
```

---

## Task 8: Integration - Wire Components Together

**Files:**
- Modify: `frontend/src/app/page.tsx` or relevant layout component

**Step 1: Find and update the main page component**

Locate where ProteinViewerClient and analysis panels are rendered. Add LigandBindingPanel integration:

```typescript
import { LigandBindingPanel } from '@/components/LigandBindingPanel';
import { extractPharmacophoreFeatures } from '@/lib/ligandAnalysis';

// In the component:
const {
  ligandData,
  focusedLigandIndex,
  showPharmacophores3D,
  setShowPharmacophores3D,
  pharmacophoreFeatures,
  setPharmacophoreFeatures,
} = useStore();

// Compute pharmacophores when ligand is focused
useEffect(() => {
  if (ligandData && focusedLigandIndex !== null) {
    const ligand = ligandData.ligandDetails[focusedLigandIndex];
    if (ligand) {
      const features = extractPharmacophoreFeatures(ligand.contacts);
      setPharmacophoreFeatures(features);
    }
  } else {
    setPharmacophoreFeatures(null);
  }
}, [ligandData, focusedLigandIndex, setPharmacophoreFeatures]);

// In JSX:
<LigandBindingPanel
  show3D={showPharmacophores3D}
  onToggle3D={setShowPharmacophores3D}
/>

<ProteinViewerClient
  // ... existing props
  pharmacophoreFeatures={pharmacophoreFeatures}
  showPharmacophores={showPharmacophores3D}
/>
```

**Step 2: Commit**

```bash
git add frontend/src/app/page.tsx
git commit -m "feat: integrate LigandBindingPanel with ProteinViewer

- Wire pharmacophore state between panel and viewer
- Auto-compute features when ligand is focused
- Pass 3D toggle state to viewer"
```

---

## Task 9: Add Vitest Configuration and Run All Tests

**Files:**
- Create: `frontend/vitest.config.ts` (if not exists)

**Step 1: Ensure vitest is configured**

Check if vitest.config.ts exists, create if needed:

```typescript
// frontend/vitest.config.ts
import { defineConfig } from 'vitest/config';
import react from '@vitejs/plugin-react';
import path from 'path';

export default defineConfig({
  plugins: [react()],
  test: {
    environment: 'jsdom',
    globals: true,
  },
  resolve: {
    alias: {
      '@': path.resolve(__dirname, './src'),
    },
  },
});
```

**Step 2: Run all tests**

```bash
cd frontend && npm test -- --run
```

Expected: All tests pass

**Step 3: Final commit**

```bash
git add -A
git commit -m "chore: finalize ligand v2 enhancements port

- All tests passing
- TypeScript builds without errors
- ESLint clean"
```

---

## Summary

**Total Tasks:** 9

**New Files Created:**
- `frontend/src/lib/__tests__/ligandAnalysis.test.ts`
- `frontend/src/lib/__tests__/ligandExport.test.ts`
- `frontend/src/lib/ligandExport.ts`
- `frontend/src/components/LigandBindingPanel.tsx`
- `frontend/vitest.config.ts` (if needed)

**Files Modified:**
- `frontend/src/lib/ligandAnalysis.ts` - Added pharmacophore extraction, shell analysis, binder quality score
- `frontend/src/lib/store.ts` - Added pharmacophore state
- `frontend/src/components/ProteinViewerClient.tsx` - Added 3D pharmacophore rendering
- `frontend/src/app/page.tsx` - Integration wiring

**Key Differences from Protein_engineering_tools Implementation:**
1. Modular file structure (separate ligandExport.ts, LigandBindingPanel.tsx)
2. Next.js client/server component patterns
3. Zustand instead of React state
4. Tailwind CSS instead of inline styles
5. TypeScript strict mode
