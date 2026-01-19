# Automatic Catalytic Residue Detection - Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add automatic catalytic residue detection from M-CSA/P2Rank APIs with visual selection via Molstar right-click menu.

**Architecture:** Backend queries M-CSA by PDB ID (fallback to P2Rank for unknown structures). Frontend displays suggestions below viewer with blue (curated) / orange (predicted) highlighting. Right-click context menu enables adding residues with atom type selection.

**Tech Stack:** FastAPI (backend), React/Zustand (frontend), Molstar (3D viewer), M-CSA API, PrankWeb API

**Design Doc:** `docs/plans/2026-01-19-automatic-catalytic-residue-detection-design.md`

---

## Phase 1: Backend API

### Task 1.1: Create M-CSA Client

**Files:**
- Create: `backend/app/services/mcsa_client.py`

**Step 1: Write the M-CSA client module**

```python
"""
M-CSA (Mechanism and Catalytic Site Atlas) API client.
Queries curated catalytic residue data by PDB ID.
"""

import httpx
from typing import List, Optional
from pydantic import BaseModel


class CatalyticResidue(BaseModel):
    """A catalytic residue from M-CSA."""
    chain: str
    residue: int
    name: str
    role: Optional[str] = None
    confidence: float = 1.0
    source: str = "mcsa"


MCSA_API_BASE = "https://www.ebi.ac.uk/thornton-srv/m-csa/api"


async def query_mcsa_by_pdb(pdb_id: str) -> List[CatalyticResidue]:
    """
    Query M-CSA for catalytic residues by PDB ID.

    Args:
        pdb_id: 4-character PDB ID (e.g., "1TRZ")

    Returns:
        List of catalytic residues with chain, position, name, and role
    """
    if not pdb_id or len(pdb_id) != 4:
        return []

    pdb_id = pdb_id.upper()
    url = f"{MCSA_API_BASE}/residues/?pdb_id={pdb_id}&format=json"

    try:
        async with httpx.AsyncClient(timeout=10.0) as client:
            response = await client.get(url)

            if response.status_code != 200:
                return []

            data = response.json()

            if not data:
                return []

            residues = []
            seen = set()  # Deduplicate by chain+residue

            for entry in data:
                # Extract residue info from M-CSA response
                chains = entry.get("residue_chains", {})
                chain = chains.get("chain_name", "A")
                resid = chains.get("resid")
                code = chains.get("code", "UNK")
                role = entry.get("roles_summary", "")

                if resid is None:
                    continue

                key = f"{chain}{resid}"
                if key in seen:
                    continue
                seen.add(key)

                residues.append(CatalyticResidue(
                    chain=chain,
                    residue=int(resid),
                    name=code,
                    role=role if role else None,
                    confidence=1.0,
                    source="mcsa"
                ))

            return residues

    except Exception as e:
        print(f"[M-CSA] Query failed: {e}")
        return []
```

**Step 2: Verify the module is syntactically correct**

Run: `cd "G:\Github_local_repo\Banta_Lab_RFdiffusion\backend" && python -c "from app.services.mcsa_client import query_mcsa_by_pdb; print('OK')"`
Expected: `OK`

**Step 3: Commit**

```bash
git add backend/app/services/mcsa_client.py
git commit -m "feat(backend): add M-CSA API client for catalytic residue lookup"
```

---

### Task 1.2: Create PrankWeb Client

**Files:**
- Create: `backend/app/services/prankweb_client.py`

**Step 1: Write the PrankWeb client module**

```python
"""
PrankWeb/P2Rank API client for binding site prediction.
Used as fallback when M-CSA has no data for a structure.
"""

import httpx
import asyncio
from typing import List, Optional
from pydantic import BaseModel


class PredictedResidue(BaseModel):
    """A predicted binding site residue from P2Rank."""
    chain: str
    residue: int
    name: str
    confidence: float
    source: str = "p2rank"


PRANKWEB_API = "https://prankweb.cz/api/v2"


async def query_prankweb(pdb_content: str, max_wait: int = 30) -> List[PredictedResidue]:
    """
    Submit structure to PrankWeb for binding site prediction.

    Args:
        pdb_content: PDB file content as string
        max_wait: Maximum seconds to wait for results

    Returns:
        List of predicted binding site residues from top pocket
    """
    if not pdb_content or len(pdb_content) < 100:
        return []

    try:
        async with httpx.AsyncClient(timeout=60.0) as client:
            # Submit prediction job
            submit_response = await client.post(
                f"{PRANKWEB_API}/predict/structure",
                json={"structure": pdb_content},
                headers={"Content-Type": "application/json"}
            )

            if submit_response.status_code not in (200, 201, 202):
                print(f"[PrankWeb] Submit failed: {submit_response.status_code}")
                return []

            job_data = submit_response.json()
            job_id = job_data.get("id") or job_data.get("taskId")

            if not job_id:
                # Some versions return results directly
                return _parse_prankweb_results(job_data)

            # Poll for results
            for _ in range(max_wait):
                await asyncio.sleep(1)

                status_response = await client.get(f"{PRANKWEB_API}/predict/{job_id}")

                if status_response.status_code != 200:
                    continue

                status_data = status_response.json()
                status = status_data.get("status", "").lower()

                if status in ("completed", "finished", "done"):
                    return _parse_prankweb_results(status_data)
                elif status in ("failed", "error"):
                    print(f"[PrankWeb] Job failed: {status_data.get('error')}")
                    return []

            print("[PrankWeb] Timeout waiting for results")
            return []

    except Exception as e:
        print(f"[PrankWeb] Query failed: {e}")
        return []


def _parse_prankweb_results(data: dict) -> List[PredictedResidue]:
    """Parse PrankWeb response into residue list."""
    residues = []
    seen = set()

    # Get pockets from response
    pockets = data.get("pockets", [])
    if not pockets:
        pockets = data.get("predictions", {}).get("pockets", [])

    if not pockets:
        return []

    # Use top pocket (highest score)
    top_pocket = pockets[0]
    pocket_residues = top_pocket.get("residues", [])
    pocket_score = top_pocket.get("score", 0.5)

    for res in pocket_residues:
        chain = res.get("chain", "A")
        resnum = res.get("residue_number") or res.get("resNum")
        resname = res.get("residue_name") or res.get("resName", "UNK")
        score = res.get("score", pocket_score)

        if resnum is None:
            continue

        key = f"{chain}{resnum}"
        if key in seen:
            continue
        seen.add(key)

        residues.append(PredictedResidue(
            chain=chain,
            residue=int(resnum),
            name=resname,
            confidence=float(score) if score else 0.5,
            source="p2rank"
        ))

    # Sort by confidence, take top 10
    residues.sort(key=lambda r: r.confidence, reverse=True)
    return residues[:10]
```

**Step 2: Verify the module is syntactically correct**

Run: `cd "G:\Github_local_repo\Banta_Lab_RFdiffusion\backend" && python -c "from app.services.prankweb_client import query_prankweb; print('OK')"`
Expected: `OK`

**Step 3: Commit**

```bash
git add backend/app/services/prankweb_client.py
git commit -m "feat(backend): add PrankWeb/P2Rank client for binding site prediction"
```

---

### Task 1.3: Create Catalytic Suggestions Endpoint

**Files:**
- Create: `backend/app/routers/catalytic_suggestions.py`
- Modify: `backend/app/main.py`

**Step 1: Write the router module**

```python
"""
Catalytic residue suggestions API endpoint.
Queries M-CSA first, falls back to P2Rank for unknown structures.
"""

import re
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Optional, Literal

from app.services.mcsa_client import query_mcsa_by_pdb, CatalyticResidue
from app.services.prankweb_client import query_prankweb, PredictedResidue


router = APIRouter(prefix="/api", tags=["catalytic"])


class SuggestionRequest(BaseModel):
    """Request for catalytic residue suggestions."""
    pdb_id: Optional[str] = None
    pdb_content: str


class SuggestionResidue(BaseModel):
    """A suggested catalytic residue."""
    chain: str
    residue: int
    name: str
    role: Optional[str] = None
    confidence: float
    source: Literal["mcsa", "p2rank"]


class SuggestionResponse(BaseModel):
    """Response containing catalytic residue suggestions."""
    source: Literal["mcsa", "p2rank", "none"]
    residues: List[SuggestionResidue]


def extract_pdb_id_from_content(pdb_content: str) -> Optional[str]:
    """Extract PDB ID from HEADER line of PDB file."""
    if not pdb_content:
        return None

    # Look for HEADER line with PDB ID
    # Format: HEADER    HYDROLASE                               01-JAN-00   1ABC
    for line in pdb_content.split('\n')[:20]:
        if line.startswith('HEADER'):
            # PDB ID is typically at positions 62-66
            if len(line) >= 66:
                pdb_id = line[62:66].strip()
                if len(pdb_id) == 4 and pdb_id.isalnum():
                    return pdb_id.upper()
            # Also try regex for flexibility
            match = re.search(r'\b([0-9][A-Za-z0-9]{3})\s*$', line)
            if match:
                return match.group(1).upper()

    return None


@router.post("/catalytic-suggestions", response_model=SuggestionResponse)
async def get_catalytic_suggestions(request: SuggestionRequest) -> SuggestionResponse:
    """
    Get catalytic residue suggestions for a protein structure.

    First queries M-CSA for curated catalytic residues.
    If no results, falls back to P2Rank binding site prediction.
    """
    # Determine PDB ID
    pdb_id = request.pdb_id
    if not pdb_id:
        pdb_id = extract_pdb_id_from_content(request.pdb_content)

    # Try M-CSA first if we have a PDB ID
    if pdb_id:
        mcsa_residues = await query_mcsa_by_pdb(pdb_id)

        if mcsa_residues:
            return SuggestionResponse(
                source="mcsa",
                residues=[
                    SuggestionResidue(
                        chain=r.chain,
                        residue=r.residue,
                        name=r.name,
                        role=r.role,
                        confidence=r.confidence,
                        source="mcsa"
                    )
                    for r in mcsa_residues
                ]
            )

    # Fallback to P2Rank
    p2rank_residues = await query_prankweb(request.pdb_content)

    if p2rank_residues:
        return SuggestionResponse(
            source="p2rank",
            residues=[
                SuggestionResidue(
                    chain=r.chain,
                    residue=r.residue,
                    name=r.name,
                    role=None,
                    confidence=r.confidence,
                    source="p2rank"
                )
                for r in p2rank_residues
            ]
        )

    # No suggestions found
    return SuggestionResponse(source="none", residues=[])
```

**Step 2: Register router in main.py**

Add to `backend/app/main.py` after line 54 (after CORS middleware):

```python
# Import and register catalytic suggestions router
from app.routers.catalytic_suggestions import router as catalytic_router
app.include_router(catalytic_router)
```

**Step 3: Verify the endpoint loads**

Run: `cd "G:\Github_local_repo\Banta_Lab_RFdiffusion\backend" && python -c "from app.main import app; print('Routes:', [r.path for r in app.routes if hasattr(r, 'path')])"`
Expected: Should include `/api/catalytic-suggestions`

**Step 4: Commit**

```bash
git add backend/app/routers/catalytic_suggestions.py backend/app/main.py
git commit -m "feat(api): add /api/catalytic-suggestions endpoint with M-CSA + P2Rank"
```

---

## Phase 2: Frontend State & Data Flow

### Task 2.1: Add Suggestions State to Zustand Store

**Files:**
- Modify: `frontend/src/lib/store.ts`

**Step 1: Add type definitions after line 17**

```typescript
// Catalytic residue suggestion from API
export interface CatalyticSuggestion {
  chain: string;
  residue: number;
  name: string;
  role?: string;
  confidence: number;
  source: 'mcsa' | 'p2rank';
}

// Bottom panel mode
export type BottomPanelMode = 'suggestions' | 'metal-analysis' | 'ligand-analysis' | 'none';
```

**Step 2: Add state fields to AppState interface (after line 241)**

```typescript
  // Catalytic residue suggestions
  catalyticSuggestions: CatalyticSuggestion[];
  suggestionsSource: 'mcsa' | 'p2rank' | 'none';
  suggestionsLoading: boolean;
  suggestionsError: string | null;
  bottomPanelMode: BottomPanelMode;

  // Actions
  setCatalyticSuggestions: (suggestions: CatalyticSuggestion[], source: 'mcsa' | 'p2rank' | 'none') => void;
  setSuggestionsLoading: (loading: boolean) => void;
  setSuggestionsError: (error: string | null) => void;
  setBottomPanelMode: (mode: BottomPanelMode) => void;
  clearSuggestions: () => void;
```

**Step 3: Add initial state and actions in create() (after line 446)**

```typescript
  // Catalytic residue suggestions
  catalyticSuggestions: [],
  suggestionsSource: 'none',
  suggestionsLoading: false,
  suggestionsError: null,
  bottomPanelMode: 'none',

  setCatalyticSuggestions: (suggestions, source) => set({
    catalyticSuggestions: suggestions,
    suggestionsSource: source,
    suggestionsError: null,
    // Auto-show panel if suggestions found
    bottomPanelMode: suggestions.length > 0 ? 'suggestions' : 'none',
  }),
  setSuggestionsLoading: (loading) => set({ suggestionsLoading: loading }),
  setSuggestionsError: (error) => set({ suggestionsError: error, suggestionsLoading: false }),
  setBottomPanelMode: (mode) => set({ bottomPanelMode: mode }),
  clearSuggestions: () => set({
    catalyticSuggestions: [],
    suggestionsSource: 'none',
    suggestionsError: null,
    bottomPanelMode: 'none',
  }),
```

**Step 4: Verify TypeScript compiles**

Run: `cd "G:\Github_local_repo\Banta_Lab_RFdiffusion\frontend" && npx tsc --noEmit`
Expected: No errors related to store.ts

**Step 5: Commit**

```bash
git add frontend/src/lib/store.ts
git commit -m "feat(store): add catalytic suggestions state slice"
```

---

### Task 2.2: Create API Client for Suggestions

**Files:**
- Create: `frontend/src/lib/catalyticDetection.ts`

**Step 1: Write the API client**

```typescript
/**
 * API client for catalytic residue detection.
 * Fetches suggestions from backend which queries M-CSA and P2Rank.
 */

import type { CatalyticSuggestion } from './store';

interface SuggestionResponse {
  source: 'mcsa' | 'p2rank' | 'none';
  residues: CatalyticSuggestion[];
}

/**
 * Fetch catalytic residue suggestions for a PDB structure.
 *
 * @param pdbContent - Full PDB file content
 * @param pdbId - Optional PDB ID if known (speeds up M-CSA lookup)
 * @param backendUrl - Backend API URL
 */
export async function fetchCatalyticSuggestions(
  pdbContent: string,
  pdbId: string | null,
  backendUrl: string
): Promise<SuggestionResponse> {
  const response = await fetch(`${backendUrl}/api/catalytic-suggestions`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      pdb_id: pdbId,
      pdb_content: pdbContent,
    }),
  });

  if (!response.ok) {
    throw new Error(`Failed to fetch suggestions: ${response.statusText}`);
  }

  return response.json();
}

/**
 * Extract PDB ID from structure content.
 * Checks HEADER line for 4-character PDB code.
 */
export function extractPdbId(pdbContent: string): string | null {
  if (!pdbContent) return null;

  const lines = pdbContent.split('\n').slice(0, 20);
  for (const line of lines) {
    if (line.startsWith('HEADER')) {
      // PDB ID at positions 62-66
      const pdbId = line.slice(62, 66).trim();
      if (pdbId.length === 4 && /^[0-9A-Z]{4}$/i.test(pdbId)) {
        return pdbId.toUpperCase();
      }
      // Try regex fallback
      const match = line.match(/\b([0-9][A-Za-z0-9]{3})\s*$/);
      if (match) {
        return match[1].toUpperCase();
      }
    }
  }

  return null;
}
```

**Step 2: Verify TypeScript compiles**

Run: `cd "G:\Github_local_repo\Banta_Lab_RFdiffusion\frontend" && npx tsc --noEmit`
Expected: No errors

**Step 3: Commit**

```bash
git add frontend/src/lib/catalyticDetection.ts
git commit -m "feat(lib): add catalytic detection API client"
```

---

### Task 2.3: Wire Up Auto-Fetch on Structure Load

**Files:**
- Modify: `frontend/src/components/tasks/shared/PdbUploader.tsx`

**Step 1: Add imports at top of file**

```typescript
import { useStore } from '@/lib/store';
import { fetchCatalyticSuggestions, extractPdbId } from '@/lib/catalyticDetection';
```

**Step 2: Add store hooks inside PdbUploader function (after line 27)**

```typescript
  const {
    backendUrl,
    setCatalyticSuggestions,
    setSuggestionsLoading,
    setSuggestionsError,
    clearSuggestions,
    setSelectedPdb,
  } = useStore();
```

**Step 3: Create fetchSuggestions helper (after store hooks)**

```typescript
  // Fetch catalytic suggestions when structure is loaded
  const fetchSuggestions = async (content: string) => {
    setSuggestionsLoading(true);
    try {
      const pdbId = extractPdbId(content);
      const result = await fetchCatalyticSuggestions(content, pdbId, backendUrl);
      setCatalyticSuggestions(result.residues, result.source);
    } catch (error) {
      console.error('Failed to fetch catalytic suggestions:', error);
      setSuggestionsError(error instanceof Error ? error.message : 'Unknown error');
    }
  };
```

**Step 4: Modify handleFileUpload to trigger fetch (replace lines 30-38)**

```typescript
  const handleFileUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (file) {
      const reader = new FileReader();
      reader.onload = (e) => {
        const content = e.target?.result as string;
        onChange(content, file.name);
        setSelectedPdb(content);
        // Auto-fetch catalytic suggestions
        fetchSuggestions(content);
      };
      reader.readAsText(file);
    }
  };
```

**Step 5: Modify handleDrop similarly (replace lines 41-51)**

```typescript
  const handleDrop = (e: React.DragEvent) => {
    e.preventDefault();
    setIsDragOver(false);
    const file = e.dataTransfer.files[0];
    if (file) {
      const reader = new FileReader();
      reader.onload = (ev) => {
        const content = ev.target?.result as string;
        onChange(content, file.name);
        setSelectedPdb(content);
        // Auto-fetch catalytic suggestions
        fetchSuggestions(content);
      };
      reader.readAsText(file);
    }
  };
```

**Step 6: Modify handleClear to clear suggestions (replace lines 54-59)**

```typescript
  const handleClear = (e: React.MouseEvent) => {
    e.stopPropagation();
    onChange(null, null);
    clearSuggestions();
    if (fileInputRef.current) {
      fileInputRef.current.value = '';
    }
  };
```

**Step 7: Verify TypeScript compiles**

Run: `cd "G:\Github_local_repo\Banta_Lab_RFdiffusion\frontend" && npx tsc --noEmit`
Expected: No errors

**Step 8: Commit**

```bash
git add frontend/src/components/tasks/shared/PdbUploader.tsx
git commit -m "feat(uploader): auto-fetch catalytic suggestions on structure load"
```

---

## Phase 3: Suggestion Panel UI

### Task 3.1: Create AtomTypeDropdown Component

**Files:**
- Create: `frontend/src/components/AtomTypeDropdown.tsx`

**Step 1: Write the component**

```typescript
'use client';

import { useState } from 'react';
import { ChevronDown } from 'lucide-react';
import { Button } from '@/components/ui/button';
import {
  Popover,
  PopoverContent,
  PopoverTrigger,
} from '@/components/ui/popover';

interface AtomTypeDropdownProps {
  onSelect: (atomType: string) => void;
  children: React.ReactNode;
  align?: 'start' | 'center' | 'end';
}

const ATOM_TYPES = [
  { value: 'BKBN', label: 'Backbone only', desc: 'N, CA, C, O' },
  { value: 'ALL', label: 'All atoms', desc: 'Entire residue' },
  { value: 'TIP', label: 'Functional tip', desc: 'Side chain tip' },
  { value: '', label: 'Flexible', desc: 'No constraints' },
];

export function AtomTypeDropdown({ onSelect, children, align = 'end' }: AtomTypeDropdownProps) {
  const [open, setOpen] = useState(false);

  const handleSelect = (atomType: string) => {
    onSelect(atomType);
    setOpen(false);
  };

  return (
    <Popover open={open} onOpenChange={setOpen}>
      <PopoverTrigger asChild>
        {children}
      </PopoverTrigger>
      <PopoverContent align={align} className="w-48 p-1">
        <div className="space-y-0.5">
          {ATOM_TYPES.map((type) => (
            <button
              key={type.value}
              onClick={() => handleSelect(type.value)}
              className="w-full text-left px-3 py-2 rounded-md hover:bg-accent transition-colors"
            >
              <div className="text-sm font-medium">{type.label}</div>
              <div className="text-xs text-muted-foreground">{type.desc}</div>
            </button>
          ))}
        </div>
      </PopoverContent>
    </Popover>
  );
}
```

**Step 2: Verify TypeScript compiles**

Run: `cd "G:\Github_local_repo\Banta_Lab_RFdiffusion\frontend" && npx tsc --noEmit`
Expected: No errors

**Step 3: Commit**

```bash
git add frontend/src/components/AtomTypeDropdown.tsx
git commit -m "feat(ui): add AtomTypeDropdown component for fixed atom selection"
```

---

### Task 3.2: Create CatalyticSuggestionsPanel Component

**Files:**
- Create: `frontend/src/components/CatalyticSuggestionsPanel.tsx`

**Step 1: Write the component**

```typescript
'use client';

import { FlaskConical, Plus, Check, X, Loader2, Sparkles } from 'lucide-react';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { AtomTypeDropdown } from '@/components/AtomTypeDropdown';
import type { CatalyticSuggestion } from '@/lib/store';

interface CatalyticSuggestionsPanelProps {
  suggestions: CatalyticSuggestion[];
  source: 'mcsa' | 'p2rank' | 'none';
  loading: boolean;
  error: string | null;
  existingResidues: Array<{ chain: string; residue: number }>;
  onAdd: (suggestion: CatalyticSuggestion, atomType: string) => void;
  onAddAll: (atomType: string) => void;
  onClose: () => void;
  onHoverResidue?: (suggestion: CatalyticSuggestion | null) => void;
}

export function CatalyticSuggestionsPanel({
  suggestions,
  source,
  loading,
  error,
  existingResidues,
  onAdd,
  onAddAll,
  onClose,
  onHoverResidue,
}: CatalyticSuggestionsPanelProps) {
  const isAdded = (s: CatalyticSuggestion) =>
    existingResidues.some((r) => r.chain === s.chain && r.residue === s.residue);

  const unadded = suggestions.filter((s) => !isAdded(s));

  if (loading) {
    return (
      <div className="border-t bg-card p-4">
        <div className="flex items-center justify-center gap-2 text-muted-foreground">
          <Loader2 className="h-4 w-4 animate-spin" />
          <span className="text-sm">Detecting catalytic residues...</span>
        </div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="border-t bg-card p-4">
        <div className="text-sm text-destructive text-center">{error}</div>
      </div>
    );
  }

  if (suggestions.length === 0) {
    return (
      <div className="border-t bg-card p-4">
        <div className="text-sm text-muted-foreground text-center">
          No catalytic residues detected. Add manually or right-click in viewer.
        </div>
      </div>
    );
  }

  return (
    <div className="border-t bg-card">
      {/* Header */}
      <div className="flex items-center justify-between px-4 py-2 border-b">
        <h3 className="font-medium text-sm flex items-center gap-2">
          <FlaskConical className="h-4 w-4" />
          Catalytic Residue Suggestions
          <Badge variant="secondary" className="text-xs">
            {suggestions.length}
          </Badge>
        </h3>
        <div className="flex items-center gap-2">
          {unadded.length > 0 && (
            <AtomTypeDropdown onSelect={onAddAll}>
              <Button variant="outline" size="sm" className="h-7 text-xs">
                <Plus className="h-3 w-3 mr-1" />
                Add All ({unadded.length})
              </Button>
            </AtomTypeDropdown>
          )}
          <Button variant="ghost" size="icon" className="h-7 w-7" onClick={onClose}>
            <X className="h-4 w-4" />
          </Button>
        </div>
      </div>

      {/* Content */}
      <div className="p-3 max-h-40 overflow-y-auto">
        <div className="grid grid-cols-2 md:grid-cols-3 gap-2">
          {suggestions.map((s) => {
            const added = isAdded(s);
            return (
              <div
                key={`${s.chain}${s.residue}`}
                className={`flex items-center justify-between p-2 rounded-md border transition-colors ${
                  added ? 'bg-muted/50 border-border' : 'bg-card border-border hover:bg-accent/50'
                }`}
                onMouseEnter={() => onHoverResidue?.(s)}
                onMouseLeave={() => onHoverResidue?.(null)}
              >
                <div className="flex items-center gap-2 min-w-0">
                  <div
                    className={`w-2 h-2 rounded-full shrink-0 ${
                      s.source === 'mcsa' ? 'bg-blue-500' : 'bg-orange-500'
                    }`}
                  />
                  <div className="min-w-0">
                    <div className="text-xs font-medium truncate">
                      {s.name} {s.chain}{s.residue}
                    </div>
                    {s.role && (
                      <div className="text-[10px] text-muted-foreground truncate">
                        {s.role}
                      </div>
                    )}
                  </div>
                </div>
                {added ? (
                  <Check className="h-4 w-4 text-green-500 shrink-0" />
                ) : (
                  <AtomTypeDropdown onSelect={(type) => onAdd(s, type)} align="end">
                    <Button variant="ghost" size="icon" className="h-6 w-6 shrink-0">
                      <Plus className="h-3 w-3" />
                    </Button>
                  </AtomTypeDropdown>
                )}
              </div>
            );
          })}
        </div>
      </div>

      {/* Footer */}
      <div className="px-4 py-2 border-t text-xs text-muted-foreground flex justify-between">
        <span>
          Source: {source === 'mcsa' ? 'M-CSA (curated)' : source === 'p2rank' ? 'P2Rank (predicted)' : 'None'}
        </span>
        <span className="flex items-center gap-3">
          <span className="flex items-center gap-1">
            <div className="w-2 h-2 rounded-full bg-blue-500" /> Curated
          </span>
          <span className="flex items-center gap-1">
            <div className="w-2 h-2 rounded-full bg-orange-500" /> Predicted
          </span>
        </span>
      </div>
    </div>
  );
}
```

**Step 2: Verify TypeScript compiles**

Run: `cd "G:\Github_local_repo\Banta_Lab_RFdiffusion\frontend" && npx tsc --noEmit`
Expected: No errors

**Step 3: Commit**

```bash
git add frontend/src/components/CatalyticSuggestionsPanel.tsx
git commit -m "feat(ui): add CatalyticSuggestionsPanel component"
```

---

### Task 3.3: Integrate Panel into ViewerPanel

**Files:**
- Modify: `frontend/src/components/layout/ViewerPanel.tsx`

**Step 1: Add imports at top**

```typescript
import { CatalyticSuggestionsPanel } from '@/components/CatalyticSuggestionsPanel';
```

**Step 2: Add store selectors in component (after line 244)**

```typescript
  const {
    catalyticSuggestions,
    suggestionsSource,
    suggestionsLoading,
    suggestionsError,
    bottomPanelMode,
    setBottomPanelMode,
  } = useStore();
```

**Step 3: Add panel mode logic (after line 273)**

```typescript
  // Panel priority: suggestions > metals > ligands > none
  const showSuggestions = bottomPanelMode === 'suggestions' && catalyticSuggestions.length > 0;
  const showAnalysis = bottomPanelMode !== 'suggestions' || catalyticSuggestions.length === 0;
```

**Step 4: Add suggestions panel rendering before analysis panel (after line 370, before "Analysis Panel" comment)**

```typescript
      {/* Catalytic Suggestions Panel */}
      {showSuggestions && (
        <CatalyticSuggestionsPanel
          suggestions={catalyticSuggestions}
          source={suggestionsSource}
          loading={suggestionsLoading}
          error={suggestionsError}
          existingResidues={[]} // Will be connected in Task 3.4
          onAdd={(s, atomType) => {
            // Will be connected in Task 3.4
            console.log('Add residue:', s, atomType);
          }}
          onAddAll={(atomType) => {
            // Will be connected in Task 3.4
            console.log('Add all with:', atomType);
          }}
          onClose={() => setBottomPanelMode('none')}
        />
      )}
```

**Step 5: Wrap analysis panel with condition**

Change the Analysis Panel div (around line 371) to:

```typescript
      {/* Analysis Panel */}
      {showAnalysis && (
        <div className="flex-1 flex flex-col overflow-hidden">
          {/* ... existing analysis panel content ... */}
        </div>
      )}
```

**Step 6: Verify TypeScript compiles**

Run: `cd "G:\Github_local_repo\Banta_Lab_RFdiffusion\frontend" && npx tsc --noEmit`
Expected: No errors

**Step 7: Commit**

```bash
git add frontend/src/components/layout/ViewerPanel.tsx
git commit -m "feat(viewer): integrate CatalyticSuggestionsPanel below viewer"
```

---

### Task 3.4: Add "Show Suggestions" Button to EnzymeForm

**Files:**
- Modify: `frontend/src/components/tasks/EnzymeForm.tsx`

**Step 1: Add imports at top**

```typescript
import { Sparkles } from 'lucide-react';
import { useStore } from '@/lib/store';
```

**Step 2: Add store hooks (after line 67)**

```typescript
  const {
    catalyticSuggestions,
    suggestionsSource,
    bottomPanelMode,
    setBottomPanelMode,
  } = useStore();
```

**Step 3: Add callback to add residue from suggestion (after addCatalyticResidue function)**

```typescript
  const addFromSuggestion = (chain: string, residue: number, name: string, atomType: string) => {
    const exists = catalyticResidues.some(
      (r) => r.chain === chain && r.residue === residue
    );
    if (!exists) {
      setCatalyticResidues([
        ...catalyticResidues,
        { chain, residue, name },
      ]);
      // Update fixed atoms for this residue
      if (atomType) {
        // Store will be updated through fixedAtomType state
      }
    }
  };
```

**Step 4: Add "Show Suggestions" button after badges (after line 280)**

```typescript
          {/* Show Suggestions button */}
          {catalyticSuggestions.length > 0 && (
            <button
              onClick={() => setBottomPanelMode('suggestions')}
              className="flex items-center gap-1.5 px-3 py-1.5 rounded-lg bg-primary/10 text-primary text-sm font-medium hover:bg-primary/20 transition-colors"
            >
              <Sparkles className="w-4 h-4" />
              Show {catalyticSuggestions.length} Suggestions
              <span className="flex gap-0.5 ml-1">
                {catalyticSuggestions.some(s => s.source === 'mcsa') && (
                  <span className="w-2 h-2 rounded-full bg-blue-500" />
                )}
                {catalyticSuggestions.some(s => s.source === 'p2rank') && (
                  <span className="w-2 h-2 rounded-full bg-orange-500" />
                )}
              </span>
            </button>
          )}
```

**Step 5: Verify TypeScript compiles**

Run: `cd "G:\Github_local_repo\Banta_Lab_RFdiffusion\frontend" && npx tsc --noEmit`
Expected: No errors

**Step 6: Commit**

```bash
git add frontend/src/components/tasks/EnzymeForm.tsx
git commit -m "feat(enzyme-form): add Show Suggestions button for catalytic residues"
```

---

## Phase 4: Molstar Integration

### Task 4.1: Add Suggestion Highlighting

**Files:**
- Create: `frontend/src/hooks/useCatalyticVisualization.ts`

**Step 1: Write the hook**

```typescript
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

    // Create highlight components for each group
    if (mcsaResidues.length > 0) {
      await createHighlightGroup(plugin, structure.ref, mcsaResidues, MCSA_COLOR, 'mcsa-suggestions');
    }

    if (p2rankResidues.length > 0) {
      await createHighlightGroup(plugin, structure.ref, p2rankResidues, P2RANK_COLOR, 'p2rank-suggestions');
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
```

**Step 2: Verify TypeScript compiles**

Run: `cd "G:\Github_local_repo\Banta_Lab_RFdiffusion\frontend" && npx tsc --noEmit`
Expected: No errors

**Step 3: Commit**

```bash
git add frontend/src/hooks/useCatalyticVisualization.ts
git commit -m "feat(hooks): add useCatalyticVisualization for Molstar highlighting"
```

---

### Task 4.2: Create MolstarContextMenu Component

**Files:**
- Create: `frontend/src/components/MolstarContextMenu.tsx`

**Step 1: Write the component**

```typescript
'use client';

import { useState, useEffect, useCallback } from 'react';
import { AtomTypeDropdown } from '@/components/AtomTypeDropdown';

interface ResidueInfo {
  chain: string;
  residue: number;
  name: string;
}

interface ContextMenuState {
  visible: boolean;
  x: number;
  y: number;
  residue: ResidueInfo | null;
}

interface MolstarContextMenuProps {
  onAddResidue: (residue: ResidueInfo, atomType: string) => void;
  onRemoveResidue?: (residue: ResidueInfo) => void;
  isResidueSelected?: (residue: ResidueInfo) => boolean;
}

export function useMolstarContextMenu({
  onAddResidue,
  onRemoveResidue,
  isResidueSelected,
}: MolstarContextMenuProps) {
  const [menuState, setMenuState] = useState<ContextMenuState>({
    visible: false,
    x: 0,
    y: 0,
    residue: null,
  });

  const showMenu = useCallback((x: number, y: number, residue: ResidueInfo) => {
    setMenuState({ visible: true, x, y, residue });
  }, []);

  const hideMenu = useCallback(() => {
    setMenuState((prev) => ({ ...prev, visible: false }));
  }, []);

  // Close menu on click outside
  useEffect(() => {
    if (!menuState.visible) return;

    const handleClick = () => hideMenu();
    const handleEscape = (e: KeyboardEvent) => {
      if (e.key === 'Escape') hideMenu();
    };

    document.addEventListener('click', handleClick);
    document.addEventListener('keydown', handleEscape);

    return () => {
      document.removeEventListener('click', handleClick);
      document.removeEventListener('keydown', handleEscape);
    };
  }, [menuState.visible, hideMenu]);

  const MenuComponent = menuState.visible && menuState.residue && (
    <div
      className="fixed z-50 bg-popover border rounded-md shadow-lg py-1 min-w-[180px]"
      style={{ left: menuState.x, top: menuState.y }}
      onClick={(e) => e.stopPropagation()}
    >
      {/* Header */}
      <div className="px-3 py-1.5 text-sm font-medium border-b">
        {menuState.residue.name} {menuState.residue.chain}{menuState.residue.residue}
      </div>

      {/* Add as catalytic */}
      {!isResidueSelected?.(menuState.residue) && (
        <div className="relative group">
          <button className="w-full text-left px-3 py-2 text-sm hover:bg-accent flex items-center justify-between">
            Add as catalytic residue
            <span className="text-muted-foreground">→</span>
          </button>
          <div className="absolute left-full top-0 ml-1 hidden group-hover:block">
            <div className="bg-popover border rounded-md shadow-lg py-1 min-w-[140px]">
              {[
                { value: 'BKBN', label: 'Backbone only' },
                { value: 'ALL', label: 'All atoms' },
                { value: 'TIP', label: 'Functional tip' },
                { value: '', label: 'Flexible' },
              ].map((opt) => (
                <button
                  key={opt.value}
                  className="w-full text-left px-3 py-1.5 text-sm hover:bg-accent"
                  onClick={() => {
                    if (menuState.residue) {
                      onAddResidue(menuState.residue, opt.value);
                    }
                    hideMenu();
                  }}
                >
                  {opt.label}
                </button>
              ))}
            </div>
          </div>
        </div>
      )}

      {/* Remove from catalytic */}
      {isResidueSelected?.(menuState.residue) && onRemoveResidue && (
        <button
          className="w-full text-left px-3 py-2 text-sm hover:bg-accent text-destructive"
          onClick={() => {
            if (menuState.residue) {
              onRemoveResidue(menuState.residue);
            }
            hideMenu();
          }}
        >
          Remove from catalytic residues
        </button>
      )}
    </div>
  );

  return { showMenu, hideMenu, MenuComponent };
}
```

**Step 2: Verify TypeScript compiles**

Run: `cd "G:\Github_local_repo\Banta_Lab_RFdiffusion\frontend" && npx tsc --noEmit`
Expected: No errors

**Step 3: Commit**

```bash
git add frontend/src/components/MolstarContextMenu.tsx
git commit -m "feat(ui): add MolstarContextMenu for right-click residue selection"
```

---

### Task 4.3: Integrate Context Menu with ProteinViewerClient

**Files:**
- Modify: `frontend/src/components/ProteinViewerClient.tsx`

**Step 1: Add imports**

```typescript
import { useMolstarContextMenu } from '@/components/MolstarContextMenu';
```

**Step 2: Add context menu props to interface (after line 42)**

```typescript
  // Context menu callbacks
  onAddCatalyticResidue?: (chain: string, residue: number, name: string, atomType: string) => void;
  onRemoveCatalyticResidue?: (chain: string, residue: number) => void;
  isCatalyticResidue?: (chain: string, residue: number) => boolean;
```

**Step 3: Set up context menu hook in component (after line 96)**

```typescript
    // Context menu for right-click residue selection
    const { showMenu, hideMenu, MenuComponent } = useMolstarContextMenu({
      onAddResidue: (residue, atomType) => {
        onAddCatalyticResidue?.(residue.chain, residue.residue, residue.name, atomType);
      },
      onRemoveResidue: (residue) => {
        onRemoveCatalyticResidue?.(residue.chain, residue.residue);
      },
      isResidueSelected: (residue) =>
        isCatalyticResidue?.(residue.chain, residue.residue) ?? false,
    });
```

**Step 4: Add right-click handler after plugin initialization (in the useEffect that initializes Molstar)**

After the plugin is ready (around line 200), add:

```typescript
        // Subscribe to right-click events for context menu
        globalPlugin.behaviors.interaction.click.subscribe((event) => {
          // Check for right-click (button 2)
          if (event.button !== 2) return;

          const loci = event.current.loci;
          if (loci.kind !== 'element-loci') return;

          try {
            const { StructureElement, StructureProperties } = require('molstar/lib/mol-model/structure');
            const loc = StructureElement.Location.create(loci.structure);
            StructureElement.Location.set(loc, loci.structure, loci.elements[0].unit, loci.elements[0].indices[0]);

            const chain = StructureProperties.chain.auth_asym_id(loc);
            const residue = StructureProperties.residue.auth_seq_id(loc);
            const name = StructureProperties.atom.label_comp_id(loc);

            showMenu(event.pageX ?? event.x, event.pageY ?? event.y, { chain, residue, name });
          } catch (error) {
            console.error('[ProteinViewer] Failed to get residue info:', error);
          }
        });
```

**Step 5: Render context menu in component return (before closing div)**

```typescript
      {MenuComponent}
```

**Step 6: Prevent default context menu on viewer container**

Add to the container div:

```typescript
      onContextMenu={(e) => e.preventDefault()}
```

**Step 7: Verify TypeScript compiles**

Run: `cd "G:\Github_local_repo\Banta_Lab_RFdiffusion\frontend" && npx tsc --noEmit`
Expected: No errors

**Step 8: Commit**

```bash
git add frontend/src/components/ProteinViewerClient.tsx
git commit -m "feat(viewer): integrate right-click context menu for residue selection"
```

---

## Phase 5: Wire Up Complete Data Flow

### Task 5.1: Connect EnzymeForm to ViewerPanel Actions

**Files:**
- Modify: `frontend/src/components/layout/ViewerPanel.tsx`
- Modify: `frontend/src/components/tasks/EnzymeForm.tsx`

This task connects all the pieces:
1. ViewerPanel passes add/remove callbacks to CatalyticSuggestionsPanel
2. EnzymeForm exposes its catalytic residue state and setters
3. Create a shared context or use store for cross-component communication

**Step 1: Create a shared store slice for enzyme form state**

Add to `frontend/src/lib/store.ts` (in AppState interface):

```typescript
  // Enzyme form shared state (for cross-component communication)
  enzymeCatalyticResidues: Array<{ chain: string; residue: number; name: string }>;
  enzymeFixedAtomTypes: Record<string, string>;
  addEnzymeCatalyticResidue: (chain: string, residue: number, name: string, atomType: string) => void;
  removeEnzymeCatalyticResidue: (chain: string, residue: number) => void;
  clearEnzymeCatalyticResidues: () => void;
```

**Step 2: Add implementations in store create()**

```typescript
  enzymeCatalyticResidues: [],
  enzymeFixedAtomTypes: {},
  addEnzymeCatalyticResidue: (chain, residue, name, atomType) => set((state) => {
    const key = `${chain}${residue}`;
    if (state.enzymeCatalyticResidues.some((r) => r.chain === chain && r.residue === residue)) {
      return state; // Already exists
    }
    return {
      enzymeCatalyticResidues: [...state.enzymeCatalyticResidues, { chain, residue, name }],
      enzymeFixedAtomTypes: { ...state.enzymeFixedAtomTypes, [key]: atomType },
    };
  }),
  removeEnzymeCatalyticResidue: (chain, residue) => set((state) => {
    const key = `${chain}${residue}`;
    const { [key]: _, ...restAtomTypes } = state.enzymeFixedAtomTypes;
    return {
      enzymeCatalyticResidues: state.enzymeCatalyticResidues.filter(
        (r) => !(r.chain === chain && r.residue === residue)
      ),
      enzymeFixedAtomTypes: restAtomTypes,
    };
  }),
  clearEnzymeCatalyticResidues: () => set({
    enzymeCatalyticResidues: [],
    enzymeFixedAtomTypes: {},
  }),
```

**Step 3: Update ViewerPanel to use store actions**

In ViewerPanel.tsx, update the CatalyticSuggestionsPanel usage:

```typescript
  const {
    // ... existing
    enzymeCatalyticResidues,
    addEnzymeCatalyticResidue,
  } = useStore();

  // In the JSX:
  <CatalyticSuggestionsPanel
    suggestions={catalyticSuggestions}
    source={suggestionsSource}
    loading={suggestionsLoading}
    error={suggestionsError}
    existingResidues={enzymeCatalyticResidues}
    onAdd={(s, atomType) => addEnzymeCatalyticResidue(s.chain, s.residue, s.name, atomType)}
    onAddAll={(atomType) => {
      catalyticSuggestions
        .filter((s) => !enzymeCatalyticResidues.some((r) => r.chain === s.chain && r.residue === s.residue))
        .forEach((s) => addEnzymeCatalyticResidue(s.chain, s.residue, s.name, atomType));
    }}
    onClose={() => setBottomPanelMode('none')}
  />
```

**Step 4: Update EnzymeForm to sync with store**

Replace the local `catalyticResidues` state with store state:

```typescript
  const {
    enzymeCatalyticResidues,
    enzymeFixedAtomTypes,
    addEnzymeCatalyticResidue,
    removeEnzymeCatalyticResidue,
    clearEnzymeCatalyticResidues,
  } = useStore();

  // Use enzymeCatalyticResidues instead of local catalyticResidues state
  // Remove the local useState for catalyticResidues
```

**Step 5: Verify TypeScript compiles**

Run: `cd "G:\Github_local_repo\Banta_Lab_RFdiffusion\frontend" && npx tsc --noEmit`
Expected: No errors

**Step 6: Commit**

```bash
git add frontend/src/lib/store.ts frontend/src/components/layout/ViewerPanel.tsx frontend/src/components/tasks/EnzymeForm.tsx
git commit -m "feat: wire up complete data flow between enzyme form and viewer panel"
```

---

## Phase 6: Testing & Polish

### Task 6.1: Manual Integration Test

**Steps:**
1. Start backend: `cd backend && uvicorn app.main:app --reload`
2. Start frontend: `cd frontend && npm run dev`
3. Navigate to Enzyme Scaffold design
4. Upload a PDB file (e.g., 1TRZ - trypsin)
5. Verify:
   - Loading spinner shows during detection
   - Suggestions appear in panel below viewer
   - Blue/orange colors match source (M-CSA vs P2Rank)
   - "Show Suggestions" button appears in form
   - Right-click on residue shows context menu
   - Adding residues updates both form and viewer
   - "Add All" works correctly

### Task 6.2: Error Handling Polish

**Files:**
- Modify: `frontend/src/components/CatalyticSuggestionsPanel.tsx`

Add retry button for errors:

```typescript
  if (error) {
    return (
      <div className="border-t bg-card p-4">
        <div className="flex items-center justify-center gap-3">
          <span className="text-sm text-destructive">{error}</span>
          <Button variant="outline" size="sm" onClick={onRetry}>
            Retry
          </Button>
        </div>
      </div>
    );
  }
```

---

## Summary

**New Files Created:**
- `backend/app/services/mcsa_client.py` - M-CSA API client
- `backend/app/services/prankweb_client.py` - PrankWeb/P2Rank client
- `backend/app/routers/catalytic_suggestions.py` - API endpoint
- `frontend/src/lib/catalyticDetection.ts` - Frontend API client
- `frontend/src/components/AtomTypeDropdown.tsx` - Atom type selector
- `frontend/src/components/CatalyticSuggestionsPanel.tsx` - Suggestions UI
- `frontend/src/components/MolstarContextMenu.tsx` - Right-click menu
- `frontend/src/hooks/useCatalyticVisualization.ts` - Molstar highlighting

**Files Modified:**
- `backend/app/main.py` - Register new router
- `frontend/src/lib/store.ts` - Add suggestions + enzyme state
- `frontend/src/components/tasks/shared/PdbUploader.tsx` - Auto-fetch on load
- `frontend/src/components/layout/ViewerPanel.tsx` - Show suggestions panel
- `frontend/src/components/tasks/EnzymeForm.tsx` - Show Suggestions button
- `frontend/src/components/ProteinViewerClient.tsx` - Right-click handling

**Total: 8 new files, 6 modified files**

---

Plan complete and saved to `docs/plans/2026-01-19-catalytic-residue-detection-implementation.md`.

**Two execution options:**

**1. Subagent-Driven (this session)** - I dispatch fresh subagent per task, review between tasks, fast iteration

**2. Parallel Session (separate)** - Open new session with executing-plans, batch execution with checkpoints

**Which approach?**
