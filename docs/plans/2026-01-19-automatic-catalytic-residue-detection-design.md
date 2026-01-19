# Automatic Catalytic Residue Detection & Visual Selection

**Date:** 2026-01-19
**Status:** Approved
**Author:** Claude + User

## Overview

Enhance the Enzyme Scaffold design workflow with automatic catalytic residue detection from curated databases (M-CSA) and ML-based prediction (P2Rank), combined with visual selection via right-click context menu in Molstar.

## Goals

1. **Auto-detect on load**: Query M-CSA/P2Rank immediately when a PDB is loaded or uploaded
2. **Visual feedback**: Highlight suggested residues in Molstar (blue = curated, orange = predicted)
3. **Easy addition**: Right-click context menu with atom type selection (Backbone/All/Tip/Flexible)
4. **User control**: Show suggestions as recommendations, not auto-added to the list

## Architecture

### Data Flow

```
PDB Upload/Load
    ↓
┌─────────────────────────────────────────┐
│  Backend: Catalytic Residue Detection   │
├─────────────────────────────────────────┤
│  1. Extract PDB ID from structure       │
│  2. Query M-CSA API by PDB ID           │
│  3. If no results → Query PrankWeb API  │
│  4. Return unified suggestion list      │
└─────────────────────────────────────────┘
    ↓
┌─────────────────────────────────────────┐
│  Frontend: Suggestion Display           │
├─────────────────────────────────────────┤
│  1. Highlight residues in Molstar       │
│     - Blue = M-CSA curated              │
│     - Orange = P2Rank predicted         │
│  2. Show suggestion panel below viewer  │
│  3. Enable right-click context menu     │
└─────────────────────────────────────────┘
    ↓
User interaction (click suggestions or right-click)
    ↓
Catalytic residues list populated
```

### Detection Triggers

| Trigger | Action |
|---------|--------|
| **Load by PDB ID** | Query M-CSA directly with PDB ID (fastest path) |
| **Upload PDB file** | Parse HEADER for PDB ID → M-CSA, fallback to P2Rank |
| **Upload without ID** | Skip M-CSA, go directly to P2Rank |

## External APIs

### M-CSA (Mechanism and Catalytic Site Atlas)

- **Endpoint**: `https://www.ebi.ac.uk/thornton-srv/m-csa/api/residues/?format=json`
- **Coverage**: 1003 curated entries, ~15,500 PDB structures
- **Data**: Chain, residue number, name, catalytic role, mechanism
- **Auth**: None required

### PrankWeb / P2Rank

- **Endpoint**: `https://prankweb.cz/api/v2/prediction`
- **Method**: Submit structure, poll for results (<5 seconds typical)
- **Data**: Binding pocket residues with confidence scores
- **Auth**: None required

## API Design

### New Endpoint

```
POST /api/catalytic-suggestions
Content-Type: application/json

Request:
{
  "pdb_id": "1TRZ",           // Optional, if known
  "pdb_content": "ATOM..."    // Required, full PDB content
}

Response:
{
  "source": "mcsa" | "p2rank" | "none",
  "residues": [
    {
      "chain": "A",
      "residue": 57,
      "name": "HIS",
      "role": "Proton shuttle",    // M-CSA only
      "confidence": 1.0,
      "source": "mcsa"
    }
  ]
}
```

## UI Design

### Layout

```
┌────────────────────────────────────────────────────────────────────┐
│                        HEADER BAR                                  │
├──────────────────────┬─────────────────────────────────────────────┤
│                      │                                             │
│   TASK PANEL         │         MOLSTAR VIEWER                      │
│                      │                                             │
│  ┌────────────────┐  │  ┌─────────────────────────────────────┐   │
│  │ Enzyme Form    │  │  │                                     │   │
│  │                │  │  │     3D Structure                    │   │
│  │ Structure...   │  │  │     (with blue/orange highlights)   │   │
│  │                │  │  │                                     │   │
│  │ Catalytic      │  │  └─────────────────────────────────────┘   │
│  │ Residues:      │  │                                             │
│  │ [HIS A57 ×]    │  │  ┌─────────────────────────────────────┐   │
│  │ [Show Panel]   │  │  │ 🔬 CATALYTIC SUGGESTIONS            │   │
│  │                │  │  │                                     │   │
│  │ + Add manually │  │  │ 🔵 HIS 57 (A) - Proton shuttle [+] │   │
│  │                │  │  │ 🔵 ASP 102 (A) - Electrophile  [+] │   │
│  └────────────────┘  │  │ 🔵 SER 195 (A) - Nucleophile   [+] │   │
│                      │  │                                     │   │
│                      │  │ Source: M-CSA       [Add All ▼]     │   │
│                      │  └─────────────────────────────────────┘   │
└──────────────────────┴─────────────────────────────────────────────┘
```

### Panel Priority (Below Viewer)

The bottom panel shows in order of priority:
1. **Suggestions** - if available and no analysis results yet
2. **Metal Analysis** - if metal coordination data exists
3. **Ligand Analysis** - if ligand binding data exists
4. **None** - collapsed

### Right-Click Context Menu

```
┌─────────────────────────────┐
│ HIS A57                     │
├─────────────────────────────┤
│ Add as catalytic residue  → │──┐
├─────────────────────────────┤  │  ┌──────────────────┐
│ Remove from catalytic       │  └──│ Backbone only    │
└─────────────────────────────┘     │ All atoms        │
                                    │ Functional tip   │
                                    │ Flexible         │
                                    └──────────────────┘
```

### Color Scheme

| Source | Color | Hex |
|--------|-------|-----|
| M-CSA (curated) | Blue | `#2563eb` |
| P2Rank (predicted) | Orange | `#f97316` |

## State Management

### Zustand Store Additions

```typescript
interface CatalyticSuggestion {
  chain: string;
  residue: number;
  name: string;
  role?: string;
  confidence: number;
  source: 'mcsa' | 'p2rank';
}

interface CatalyticSuggestionsSlice {
  catalyticSuggestions: CatalyticSuggestion[];
  suggestionsSource: 'mcsa' | 'p2rank' | 'none';
  suggestionsLoading: boolean;
  suggestionsError: string | null;
  bottomPanelMode: 'suggestions' | 'metal-analysis' | 'ligand-analysis' | 'none';

  fetchCatalyticSuggestions: (pdbId: string | null, pdbContent: string) => Promise<void>;
  clearSuggestions: () => void;
  setBottomPanelMode: (mode: BottomPanelMode) => void;
}
```

## File Changes

### New Files

| File | Purpose |
|------|---------|
| `frontend/src/components/CatalyticSuggestionsPanel.tsx` | Panel below viewer showing suggestions |
| `frontend/src/components/MolstarContextMenu.tsx` | Right-click menu with atom selection |
| `frontend/src/components/AtomTypeDropdown.tsx` | Reusable dropdown for BKBN/ALL/TIP/Flexible |
| `frontend/src/hooks/useCatalyticVisualization.ts` | Syncs suggestions with Molstar highlights |
| `frontend/src/lib/catalyticDetection.ts` | API client for suggestion fetching |
| `backend/app/routers/catalytic_suggestions.py` | New API endpoint |
| `backend/app/services/mcsa_client.py` | M-CSA API integration |
| `backend/app/services/prankweb_client.py` | PrankWeb/P2Rank API integration |

### Modified Files

| File | Changes |
|------|---------|
| `frontend/src/lib/store.ts` | Add suggestions slice, bottomPanelMode |
| `frontend/src/components/tasks/EnzymeForm.tsx` | Add "Show Suggestions" button, collapsible manual entry |
| `frontend/src/components/tasks/shared/PdbUploader.tsx` | Trigger suggestion fetch on load |
| `frontend/src/components/ProteinViewerClient.tsx` | Add right-click handler, suggestion highlighting |
| `frontend/src/components/layout/ViewerPanel.tsx` | Add bottom panel switcher for suggestions |
| `backend/app/main.py` | Register new router |

## Implementation Phases

### Phase 1: Backend API
1. Create M-CSA client with PDB ID lookup
2. Create PrankWeb client with structure submission
3. Create unified `/api/catalytic-suggestions` endpoint
4. Test with sample PDB IDs

### Phase 2: Frontend State & Data Flow
1. Add suggestions slice to Zustand store
2. Create `catalyticDetection.ts` API client
3. Wire up auto-fetch on structure load
4. Test data flow end-to-end

### Phase 3: Suggestion Panel UI
1. Create `CatalyticSuggestionsPanel` component
2. Create `AtomTypeDropdown` component
3. Add bottom panel switching to `ViewerPanel`
4. Add "Show Suggestions" button to `EnzymeForm`

### Phase 4: Molstar Integration
1. Add suggestion highlighting (blue/orange)
2. Implement right-click event subscription
3. Create `MolstarContextMenu` component
4. Wire up "Add as catalytic" → form state

### Phase 5: Polish
1. Loading states and error handling
2. Hover-to-highlight in panel ↔ viewer
3. Clear suggestions on structure change
4. Test full workflow

## Dependencies

No new package installations required:
- M-CSA API: Public REST, no auth
- PrankWeb API: Public REST, no auth
- Molstar: Already installed (v5.5.0)

## References

- [M-CSA Database](https://www.ebi.ac.uk/thornton-srv/m-csa/)
- [M-CSA API Documentation](https://www.ebi.ac.uk/thornton-srv/m-csa/api/residues/?format=json)
- [P2Rank GitHub](https://github.com/rdk/p2rank)
- [PrankWeb](https://prankweb.cz/)
- [Molstar Selection Documentation](https://molstar.org/docs/plugin/selections/)
- [Molstar Click Events](https://github.com/molstar/molstar/discussions/998)
