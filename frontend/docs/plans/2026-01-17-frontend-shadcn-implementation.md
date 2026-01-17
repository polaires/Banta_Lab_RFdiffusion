# Frontend shadcn/ui Redesign Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Migrate the entire frontend from custom components to shadcn/ui with warm professional aesthetics, AI-first interface, and sidebar navigation.

**Architecture:** Three-column layout (Sidebar 240px | Main Content flex-1 | Protein Viewer 400px collapsible). AI chat as primary interface with manual mode toggle. React Hook Form + Zod for form validation.

**Tech Stack:** Next.js 16, React 19, shadcn/ui, Tailwind CSS v4, React Hook Form, Zod, Lucide React, Zustand v5

**Design Document:** `docs/plans/2026-01-17-frontend-redesign-design.md`

**Working Directory:** `G:\Github_local_repo\Banta_Lab_RFdiffusion\.worktrees\frontend-shadcn-redesign\frontend`

---

## Phase 1: Foundation

### Task 1.1: Initialize shadcn/ui

**Files:**
- Create: `components.json`
- Create: `src/lib/utils.ts`
- Modify: `tailwind.config.ts` (if needed)

**Step 1: Install shadcn/ui CLI dependencies**

Run:
```bash
cd frontend && npm install clsx tailwind-merge class-variance-authority
```

Expected: Dependencies added to package.json

**Step 2: Create lib/utils.ts with cn() helper**

Create file `src/lib/utils.ts`:
```typescript
import { type ClassValue, clsx } from 'clsx';
import { twMerge } from 'tailwind-merge';

export function cn(...inputs: ClassValue[]) {
  return twMerge(clsx(inputs));
}
```

**Step 3: Create components.json for shadcn/ui**

Create file `components.json`:
```json
{
  "$schema": "https://ui.shadcn.com/schema.json",
  "style": "new-york",
  "rsc": true,
  "tsx": true,
  "tailwind": {
    "config": "",
    "css": "src/app/globals.css",
    "baseColor": "stone",
    "cssVariables": true,
    "prefix": ""
  },
  "aliases": {
    "components": "@/components",
    "utils": "@/lib/utils",
    "ui": "@/components/ui",
    "lib": "@/lib",
    "hooks": "@/hooks"
  },
  "iconLibrary": "lucide"
}
```

**Step 4: Verify setup**

Run:
```bash
npm run build
```

Expected: Build succeeds without errors

**Step 5: Commit**

```bash
git add src/lib/utils.ts components.json package.json package-lock.json
git commit -m "feat: initialize shadcn/ui foundation with cn() utility"
```

---

### Task 1.2: Configure Warm Professional Theme

**Files:**
- Modify: `src/app/globals.css`

**Step 1: Backup existing globals.css**

Run:
```bash
cp src/app/globals.css src/app/globals.css.backup
```

**Step 2: Add shadcn/ui CSS variables with warm stone palette**

Add to `src/app/globals.css` after the @import statements:

```css
@layer base {
  :root {
    /* Backgrounds */
    --background: 30 6% 98%;
    --foreground: 24 10% 10%;
    --card: 0 0% 100%;
    --card-foreground: 24 10% 10%;
    --popover: 0 0% 100%;
    --popover-foreground: 24 10% 10%;

    /* Primary - stone-800 */
    --primary: 24 10% 15%;
    --primary-foreground: 30 6% 98%;

    /* Secondary - stone-100 */
    --secondary: 60 5% 96%;
    --secondary-foreground: 24 10% 10%;

    /* Muted - stone-100/500 */
    --muted: 60 5% 96%;
    --muted-foreground: 25 6% 45%;

    /* Accent - stone-100 */
    --accent: 60 5% 96%;
    --accent-foreground: 24 10% 10%;

    /* Semantic */
    --destructive: 0 84% 60%;
    --destructive-foreground: 0 0% 100%;

    /* Borders & Rings */
    --border: 24 6% 90%;
    --input: 24 6% 90%;
    --ring: 24 6% 65%;

    /* Chart colors */
    --chart-1: 24 10% 15%;
    --chart-2: 25 6% 45%;
    --chart-3: 24 6% 65%;
    --chart-4: 24 6% 80%;
    --chart-5: 60 5% 96%;

    /* Radius */
    --radius: 0.5rem;

    /* Sidebar specific */
    --sidebar-background: 0 0% 100%;
    --sidebar-foreground: 24 10% 10%;
    --sidebar-primary: 24 10% 15%;
    --sidebar-primary-foreground: 30 6% 98%;
    --sidebar-accent: 60 5% 96%;
    --sidebar-accent-foreground: 24 10% 10%;
    --sidebar-border: 24 6% 90%;
    --sidebar-ring: 24 6% 65%;
  }
}

@layer base {
  * {
    @apply border-border;
  }
  body {
    @apply bg-background text-foreground;
  }
}
```

**Step 3: Verify theme applies**

Run:
```bash
npm run dev
```

Expected: App loads with warm stone background (#fafaf9)

**Step 4: Commit**

```bash
git add src/app/globals.css
git commit -m "feat: configure warm professional theme with stone palette"
```

---

### Task 1.3: Install React Hook Form and Zod

**Files:**
- Modify: `package.json`

**Step 1: Install dependencies**

Run:
```bash
npm install react-hook-form @hookform/resolvers zod
```

Expected: Dependencies added

**Step 2: Verify installation**

Run:
```bash
npm ls react-hook-form zod
```

Expected: Shows installed versions

**Step 3: Commit**

```bash
git add package.json package-lock.json
git commit -m "feat: add React Hook Form and Zod for form validation"
```

---

### Task 1.4: Install Core shadcn/ui Components

**Files:**
- Create: `src/components/ui/button.tsx`
- Create: `src/components/ui/card.tsx`
- Create: `src/components/ui/separator.tsx`
- Create: `src/components/ui/scroll-area.tsx`
- Create: `src/components/ui/tooltip.tsx`

**Step 1: Install button component**

Run:
```bash
npx shadcn@latest add button -y
```

**Step 2: Install card component**

Run:
```bash
npx shadcn@latest add card -y
```

**Step 3: Install separator component**

Run:
```bash
npx shadcn@latest add separator -y
```

**Step 4: Install scroll-area component**

Run:
```bash
npx shadcn@latest add scroll-area -y
```

**Step 5: Install tooltip component**

Run:
```bash
npx shadcn@latest add tooltip -y
```

**Step 6: Verify components created**

Run:
```bash
ls src/components/ui/
```

Expected: button.tsx, card.tsx, separator.tsx, scroll-area.tsx, tooltip.tsx

**Step 7: Commit**

```bash
git add src/components/ui/
git commit -m "feat: add core shadcn/ui components (button, card, separator, scroll-area, tooltip)"
```

---

### Task 1.5: Install Form Components

**Files:**
- Create: `src/components/ui/form.tsx`
- Create: `src/components/ui/input.tsx`
- Create: `src/components/ui/label.tsx`
- Create: `src/components/ui/select.tsx`
- Create: `src/components/ui/slider.tsx`
- Create: `src/components/ui/switch.tsx`
- Create: `src/components/ui/checkbox.tsx`
- Create: `src/components/ui/radio-group.tsx`

**Step 1: Install form components**

Run:
```bash
npx shadcn@latest add form input label select slider switch checkbox radio-group -y
```

**Step 2: Verify installation**

Run:
```bash
ls src/components/ui/
```

Expected: All form components present

**Step 3: Run build to verify no errors**

Run:
```bash
npm run build
```

Expected: Build succeeds

**Step 4: Commit**

```bash
git add src/components/ui/
git commit -m "feat: add shadcn/ui form components (form, input, select, slider, switch, checkbox, radio-group)"
```

---

### Task 1.6: Install Overlay Components

**Files:**
- Create: `src/components/ui/dialog.tsx`
- Create: `src/components/ui/sheet.tsx`
- Create: `src/components/ui/alert-dialog.tsx`
- Create: `src/components/ui/popover.tsx`

**Step 1: Install overlay components**

Run:
```bash
npx shadcn@latest add dialog sheet alert-dialog popover -y
```

**Step 2: Verify installation**

Run:
```bash
ls src/components/ui/
```

Expected: dialog.tsx, sheet.tsx, alert-dialog.tsx, popover.tsx present

**Step 3: Commit**

```bash
git add src/components/ui/
git commit -m "feat: add shadcn/ui overlay components (dialog, sheet, alert-dialog, popover)"
```

---

### Task 1.7: Install Feedback Components

**Files:**
- Create: `src/components/ui/toast.tsx`
- Create: `src/components/ui/toaster.tsx` (or sonner)
- Create: `src/components/ui/alert.tsx`
- Create: `src/components/ui/progress.tsx`
- Create: `src/components/ui/badge.tsx`

**Step 1: Install toast (using sonner)**

Run:
```bash
npx shadcn@latest add sonner -y
```

**Step 2: Install alert, progress, badge**

Run:
```bash
npx shadcn@latest add alert progress badge -y
```

**Step 3: Verify installation**

Run:
```bash
ls src/components/ui/
```

Expected: sonner.tsx (or toast files), alert.tsx, progress.tsx, badge.tsx present

**Step 4: Commit**

```bash
git add src/components/ui/
git commit -m "feat: add shadcn/ui feedback components (sonner toast, alert, progress, badge)"
```

---

### Task 1.8: Install Navigation Components

**Files:**
- Create: `src/components/ui/collapsible.tsx`
- Create: `src/components/ui/toggle.tsx`
- Create: `src/components/ui/tabs.tsx`

**Step 1: Install navigation components**

Run:
```bash
npx shadcn@latest add collapsible toggle tabs -y
```

**Step 2: Verify installation**

Run:
```bash
ls src/components/ui/
```

Expected: collapsible.tsx, toggle.tsx, tabs.tsx present

**Step 3: Run full build**

Run:
```bash
npm run build
```

Expected: Build succeeds with all shadcn/ui components

**Step 4: Commit**

```bash
git add src/components/ui/
git commit -m "feat: add shadcn/ui navigation components (collapsible, toggle, tabs)"
```

---

### Task 1.9: Install Data Display Components

**Files:**
- Create: `src/components/ui/table.tsx`
- Create: `src/components/ui/avatar.tsx`

**Step 1: Install data display components**

Run:
```bash
npx shadcn@latest add table avatar -y
```

**Step 2: Verify all components**

Run:
```bash
ls -la src/components/ui/ | wc -l
```

Expected: ~20+ files (all shadcn components)

**Step 3: Run tests**

Run:
```bash
npm test
```

Expected: All tests pass

**Step 4: Commit**

```bash
git add src/components/ui/
git commit -m "feat: add shadcn/ui data display components (table, avatar)"
```

---

## Phase 2: Layout Shell

### Task 2.1: Create MainLayout Component

**Files:**
- Create: `src/components/layout/MainLayout.tsx`

**Step 1: Create layout directory**

Run:
```bash
mkdir -p src/components/layout
```

**Step 2: Create MainLayout.tsx**

Create file `src/components/layout/MainLayout.tsx`:
```tsx
'use client';

import { useState } from 'react';
import { cn } from '@/lib/utils';

interface MainLayoutProps {
  sidebar: React.ReactNode;
  main: React.ReactNode;
  viewer?: React.ReactNode;
  header?: React.ReactNode;
}

export function MainLayout({ sidebar, main, viewer, header }: MainLayoutProps) {
  const [viewerCollapsed, setViewerCollapsed] = useState(false);

  return (
    <div className="h-screen flex flex-col bg-background">
      {/* Header Bar */}
      {header && (
        <header className="h-14 border-b border-border flex items-center px-4 shrink-0">
          {header}
        </header>
      )}

      {/* Main Content Area */}
      <div className="flex-1 flex overflow-hidden">
        {/* Sidebar */}
        <aside className="w-60 border-r border-border flex flex-col shrink-0 bg-sidebar-background">
          {sidebar}
        </aside>

        {/* Main Content */}
        <main className="flex-1 overflow-auto">
          {main}
        </main>

        {/* Protein Viewer Panel */}
        {viewer && (
          <aside
            className={cn(
              'border-l border-border shrink-0 transition-all duration-200 bg-card',
              viewerCollapsed ? 'w-12' : 'w-[400px]'
            )}
          >
            <div className="h-full flex flex-col">
              <div className="h-10 border-b border-border flex items-center justify-between px-3">
                <span className={cn('text-sm font-medium', viewerCollapsed && 'hidden')}>
                  Structure Viewer
                </span>
                <button
                  onClick={() => setViewerCollapsed(!viewerCollapsed)}
                  className="p-1 hover:bg-accent rounded"
                >
                  {viewerCollapsed ? '◀' : '▶'}
                </button>
              </div>
              <div className={cn('flex-1 overflow-hidden', viewerCollapsed && 'hidden')}>
                {viewer}
              </div>
            </div>
          </aside>
        )}
      </div>
    </div>
  );
}

export default MainLayout;
```

**Step 3: Verify file created**

Run:
```bash
cat src/components/layout/MainLayout.tsx | head -20
```

Expected: Shows component code

**Step 4: Commit**

```bash
git add src/components/layout/MainLayout.tsx
git commit -m "feat: create MainLayout three-column component"
```

---

### Task 2.2: Create Sidebar Component

**Files:**
- Create: `src/components/layout/Sidebar.tsx`

**Step 1: Create Sidebar.tsx**

Create file `src/components/layout/Sidebar.tsx`:
```tsx
'use client';

import { useState } from 'react';
import {
  Plus,
  Settings,
  Command,
  Check,
  Circle,
  RefreshCw,
  FlaskConical
} from 'lucide-react';
import { cn } from '@/lib/utils';
import { Button } from '@/components/ui/button';
import { Switch } from '@/components/ui/switch';
import { ScrollArea } from '@/components/ui/scroll-area';
import { Separator } from '@/components/ui/separator';
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from '@/components/ui/tooltip';

type WorkflowStage = 'task' | 'rfd3' | 'mpnn' | 'rf3';

interface WorkflowStep {
  id: WorkflowStage;
  label: string;
  status: 'pending' | 'active' | 'completed';
  iterations?: number;
}

interface DesignHistoryItem {
  id: string;
  name: string;
  timestamp: string;
}

interface SidebarProps {
  currentStage: WorkflowStage;
  workflowSteps: WorkflowStep[];
  history: DesignHistoryItem[];
  manualMode: boolean;
  onManualModeChange: (enabled: boolean) => void;
  onStageClick: (stage: WorkflowStage) => void;
  onNewDesign: () => void;
  onHistoryClick: (id: string) => void;
  onSettingsClick: () => void;
}

export function Sidebar({
  currentStage,
  workflowSteps,
  history,
  manualMode,
  onManualModeChange,
  onStageClick,
  onNewDesign,
  onHistoryClick,
  onSettingsClick,
}: SidebarProps) {
  return (
    <TooltipProvider>
      <div className="h-full flex flex-col">
        {/* Logo */}
        <div className="p-4 border-b border-sidebar-border">
          <div className="flex items-center gap-2">
            <FlaskConical className="h-6 w-6 text-primary" />
            <div>
              <div className="font-semibold text-sm">Banta Lab</div>
              <div className="text-xs text-muted-foreground">RFdiffusion</div>
            </div>
          </div>
        </div>

        {/* New Design Button */}
        <div className="p-3">
          <Button onClick={onNewDesign} className="w-full" size="sm">
            <Plus className="h-4 w-4 mr-2" />
            New Design
          </Button>
        </div>

        <Separator />

        {/* Workflow Steps */}
        <div className="p-3">
          <div className="text-xs font-medium text-muted-foreground uppercase tracking-wider mb-2">
            Current Workflow
          </div>
          <div className="space-y-1">
            {workflowSteps.map((step) => (
              <WorkflowStepItem
                key={step.id}
                step={step}
                isActive={currentStage === step.id}
                onClick={() => onStageClick(step.id)}
              />
            ))}
          </div>
        </div>

        <Separator />

        {/* History */}
        <div className="flex-1 p-3 overflow-hidden">
          <div className="text-xs font-medium text-muted-foreground uppercase tracking-wider mb-2">
            History
          </div>
          <ScrollArea className="h-[calc(100%-24px)]">
            <div className="space-y-1">
              {history.map((item) => (
                <button
                  key={item.id}
                  onClick={() => onHistoryClick(item.id)}
                  className="w-full flex items-center gap-2 px-2 py-1.5 rounded text-sm hover:bg-sidebar-accent text-left"
                >
                  <Check className="h-3.5 w-3.5 text-muted-foreground" />
                  <span className="flex-1 truncate">{item.name}</span>
                  <span className="text-xs text-muted-foreground">{item.timestamp}</span>
                </button>
              ))}
            </div>
          </ScrollArea>
        </div>

        <Separator />

        {/* Bottom Actions */}
        <div className="p-3 space-y-2">
          <button
            onClick={onSettingsClick}
            className="w-full flex items-center gap-2 px-2 py-1.5 rounded text-sm hover:bg-sidebar-accent"
          >
            <Settings className="h-4 w-4" />
            Settings
          </button>

          <div className="flex items-center gap-2 px-2 py-1.5 text-sm text-muted-foreground">
            <Command className="h-3.5 w-3.5" />
            <span className="text-xs">K</span>
            <span className="flex-1">AI Assistant</span>
          </div>

          <div className="flex items-center justify-between px-2 py-1.5">
            <span className="text-sm">Manual Mode</span>
            <Switch
              checked={manualMode}
              onCheckedChange={onManualModeChange}
            />
          </div>
        </div>
      </div>
    </TooltipProvider>
  );
}

function WorkflowStepItem({
  step,
  isActive,
  onClick,
}: {
  step: WorkflowStep;
  isActive: boolean;
  onClick: () => void;
}) {
  const isClickable = step.status !== 'pending';

  const content = (
    <button
      onClick={isClickable ? onClick : undefined}
      disabled={!isClickable}
      className={cn(
        'w-full flex items-center gap-2 px-2 py-1.5 rounded text-sm transition-colors',
        isActive && 'bg-sidebar-accent',
        isClickable ? 'hover:bg-sidebar-accent cursor-pointer' : 'cursor-not-allowed opacity-60'
      )}
    >
      {step.status === 'completed' ? (
        <Check className="h-3.5 w-3.5 text-primary" />
      ) : step.status === 'active' ? (
        <Circle className="h-3.5 w-3.5 fill-primary text-primary" />
      ) : (
        <Circle className="h-3.5 w-3.5 text-muted-foreground" />
      )}
      <span className="flex-1 text-left">{step.label}</span>
      {step.iterations && step.iterations > 1 && (
        <span className="flex items-center gap-0.5 text-xs text-muted-foreground">
          <RefreshCw className="h-3 w-3" />
          {step.iterations}
        </span>
      )}
    </button>
  );

  if (!isClickable) {
    return (
      <Tooltip>
        <TooltipTrigger asChild>{content}</TooltipTrigger>
        <TooltipContent side="right">
          <p>Complete previous stage first</p>
        </TooltipContent>
      </Tooltip>
    );
  }

  return content;
}

export default Sidebar;
```

**Step 2: Verify file created**

Run:
```bash
wc -l src/components/layout/Sidebar.tsx
```

Expected: ~180 lines

**Step 3: Commit**

```bash
git add src/components/layout/Sidebar.tsx
git commit -m "feat: create Sidebar navigation component with workflow tracking"
```

---

### Task 2.3: Create HeaderBar Component

**Files:**
- Create: `src/components/layout/HeaderBar.tsx`

**Step 1: Create HeaderBar.tsx**

Create file `src/components/layout/HeaderBar.tsx`:
```tsx
'use client';

import { Wifi, WifiOff, User, ChevronDown } from 'lucide-react';
import { cn } from '@/lib/utils';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';

interface HeaderBarProps {
  connectionStatus: 'connected' | 'disconnected' | 'connecting';
  userName?: string;
  onConnectionClick: () => void;
  onUserClick: () => void;
}

export function HeaderBar({
  connectionStatus,
  userName,
  onConnectionClick,
  onUserClick,
}: HeaderBarProps) {
  return (
    <div className="flex-1 flex items-center justify-between">
      {/* Left side - can be empty or have breadcrumbs */}
      <div />

      {/* Right side - connection status and user */}
      <div className="flex items-center gap-3">
        {/* Connection Status */}
        <button
          onClick={onConnectionClick}
          className={cn(
            'flex items-center gap-2 px-3 py-1.5 rounded-md text-sm transition-colors',
            'hover:bg-accent'
          )}
        >
          {connectionStatus === 'connected' ? (
            <>
              <Wifi className="h-4 w-4 text-green-600" />
              <span className="text-muted-foreground">Connected</span>
            </>
          ) : connectionStatus === 'connecting' ? (
            <>
              <Wifi className="h-4 w-4 text-amber-500 animate-pulse" />
              <span className="text-muted-foreground">Connecting...</span>
            </>
          ) : (
            <>
              <WifiOff className="h-4 w-4 text-destructive" />
              <span className="text-muted-foreground">Disconnected</span>
            </>
          )}
        </button>

        {/* User Menu */}
        <button
          onClick={onUserClick}
          className="flex items-center gap-2 px-3 py-1.5 rounded-md text-sm hover:bg-accent transition-colors"
        >
          <div className="h-6 w-6 rounded-full bg-primary flex items-center justify-center">
            <User className="h-3.5 w-3.5 text-primary-foreground" />
          </div>
          {userName && <span>{userName}</span>}
          <ChevronDown className="h-4 w-4 text-muted-foreground" />
        </button>
      </div>
    </div>
  );
}

export default HeaderBar;
```

**Step 2: Commit**

```bash
git add src/components/layout/HeaderBar.tsx
git commit -m "feat: create HeaderBar component with connection status"
```

---

### Task 2.4: Create ViewerPanel Component

**Files:**
- Create: `src/components/layout/ViewerPanel.tsx`

**Step 1: Create ViewerPanel.tsx**

Create file `src/components/layout/ViewerPanel.tsx`:
```tsx
'use client';

import { RotateCcw, Maximize2, Play, Pause } from 'lucide-react';
import { Button } from '@/components/ui/button';
import { Separator } from '@/components/ui/separator';
import { cn } from '@/lib/utils';

interface StructureInfo {
  residues?: number;
  plddt?: number;
  rmsd?: number;
}

interface ViewerPanelProps {
  children: React.ReactNode;
  structureInfo?: StructureInfo;
  isSpinning?: boolean;
  onToggleSpin?: () => void;
  onReset?: () => void;
  onExpand?: () => void;
  showControls?: boolean;
}

export function ViewerPanel({
  children,
  structureInfo,
  isSpinning = false,
  onToggleSpin,
  onReset,
  onExpand,
  showControls = true,
}: ViewerPanelProps) {
  return (
    <div className="h-full flex flex-col">
      {/* Viewer Area */}
      <div className="flex-1 relative bg-muted/30">
        {children}
      </div>

      {/* Controls */}
      {showControls && (
        <>
          <Separator />
          <div className="p-2 flex items-center gap-1">
            <Button
              variant="ghost"
              size="sm"
              onClick={onToggleSpin}
              className="h-8 px-2"
            >
              {isSpinning ? (
                <Pause className="h-4 w-4" />
              ) : (
                <Play className="h-4 w-4" />
              )}
              <span className="ml-1 text-xs">Spin</span>
            </Button>
            <Button
              variant="ghost"
              size="sm"
              onClick={onReset}
              className="h-8 px-2"
            >
              <RotateCcw className="h-4 w-4" />
              <span className="ml-1 text-xs">Reset</span>
            </Button>
            <div className="flex-1" />
            <Button
              variant="ghost"
              size="sm"
              onClick={onExpand}
              className="h-8 px-2"
            >
              <Maximize2 className="h-4 w-4" />
            </Button>
          </div>
        </>
      )}

      {/* Structure Info */}
      {structureInfo && (
        <>
          <Separator />
          <div className="p-3 space-y-2">
            <div className="text-xs font-medium text-muted-foreground uppercase tracking-wider">
              Structure Info
            </div>
            <div className="grid grid-cols-3 gap-2 text-sm">
              {structureInfo.residues !== undefined && (
                <div>
                  <div className="text-muted-foreground text-xs">Residues</div>
                  <div className="font-mono">{structureInfo.residues}</div>
                </div>
              )}
              {structureInfo.plddt !== undefined && (
                <div>
                  <div className="text-muted-foreground text-xs">pLDDT</div>
                  <div className={cn(
                    'font-mono',
                    structureInfo.plddt >= 90 ? 'text-green-600' :
                    structureInfo.plddt >= 70 ? 'text-amber-600' :
                    'text-red-600'
                  )}>
                    {structureInfo.plddt.toFixed(1)}
                  </div>
                </div>
              )}
              {structureInfo.rmsd !== undefined && (
                <div>
                  <div className="text-muted-foreground text-xs">RMSD</div>
                  <div className="font-mono">{structureInfo.rmsd.toFixed(2)} Å</div>
                </div>
              )}
            </div>
          </div>
        </>
      )}

      {/* Confidence Legend */}
      <Separator />
      <div className="p-3">
        <div className="text-xs font-medium text-muted-foreground uppercase tracking-wider mb-2">
          Confidence Coloring
        </div>
        <div className="grid grid-cols-2 gap-1 text-xs">
          <div className="flex items-center gap-1.5">
            <div className="w-3 h-3 rounded bg-blue-600" />
            <span>Very High (&gt;90)</span>
          </div>
          <div className="flex items-center gap-1.5">
            <div className="w-3 h-3 rounded bg-cyan-500" />
            <span>High (70-90)</span>
          </div>
          <div className="flex items-center gap-1.5">
            <div className="w-3 h-3 rounded bg-yellow-500" />
            <span>Medium (50-70)</span>
          </div>
          <div className="flex items-center gap-1.5">
            <div className="w-3 h-3 rounded bg-orange-500" />
            <span>Low (&lt;50)</span>
          </div>
        </div>
      </div>
    </div>
  );
}

export default ViewerPanel;
```

**Step 2: Commit**

```bash
git add src/components/layout/ViewerPanel.tsx
git commit -m "feat: create ViewerPanel component with controls and structure info"
```

---

### Task 2.5: Create Layout Index Export

**Files:**
- Create: `src/components/layout/index.ts`

**Step 1: Create index.ts**

Create file `src/components/layout/index.ts`:
```typescript
export { MainLayout } from './MainLayout';
export { Sidebar } from './Sidebar';
export { HeaderBar } from './HeaderBar';
export { ViewerPanel } from './ViewerPanel';
```

**Step 2: Verify build**

Run:
```bash
npm run build
```

Expected: Build succeeds

**Step 3: Commit**

```bash
git add src/components/layout/index.ts
git commit -m "feat: add layout components index export"
```

---

## Phase 3: Connection & Notifications

### Task 3.1: Create ConnectionSheet Component

**Files:**
- Create: `src/components/connection/ConnectionSheet.tsx`

**Step 1: Create connection directory**

Run:
```bash
mkdir -p src/components/connection
```

**Step 2: Create ConnectionSheet.tsx**

Create file `src/components/connection/ConnectionSheet.tsx`:
```tsx
'use client';

import { useState } from 'react';
import { Cloud, Server, Wrench, Loader2, CheckCircle, XCircle } from 'lucide-react';
import {
  Sheet,
  SheetContent,
  SheetDescription,
  SheetHeader,
  SheetTitle,
} from '@/components/ui/sheet';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Alert, AlertDescription } from '@/components/ui/alert';
import { cn } from '@/lib/utils';

type ConnectionMode = 'runpod' | 'traditional' | 'local';
type ConnectionStatus = 'idle' | 'connecting' | 'connected' | 'failed';

interface ConnectionSheetProps {
  open: boolean;
  onOpenChange: (open: boolean) => void;
  onConnect: (mode: ConnectionMode, url: string) => Promise<boolean>;
  defaultUrl?: string;
  defaultMode?: ConnectionMode;
}

const modes = [
  {
    id: 'runpod' as const,
    label: 'RunPod Serverless',
    description: 'Recommended for most users',
    icon: Cloud,
  },
  {
    id: 'traditional' as const,
    label: 'Traditional Server',
    description: 'Self-hosted backend',
    icon: Server,
  },
  {
    id: 'local' as const,
    label: 'Local Development',
    description: 'For contributors',
    icon: Wrench,
  },
];

export function ConnectionSheet({
  open,
  onOpenChange,
  onConnect,
  defaultUrl = '',
  defaultMode = 'runpod',
}: ConnectionSheetProps) {
  const [mode, setMode] = useState<ConnectionMode>(defaultMode);
  const [url, setUrl] = useState(defaultUrl);
  const [status, setStatus] = useState<ConnectionStatus>('idle');
  const [error, setError] = useState<string | null>(null);

  const handleConnect = async () => {
    setStatus('connecting');
    setError(null);

    try {
      const success = await onConnect(mode, url);
      if (success) {
        setStatus('connected');
        setTimeout(() => onOpenChange(false), 1000);
      } else {
        setStatus('failed');
        setError('Connection failed. Please check the URL and try again.');
      }
    } catch (err) {
      setStatus('failed');
      setError(err instanceof Error ? err.message : 'Connection failed');
    }
  };

  return (
    <Sheet open={open} onOpenChange={onOpenChange}>
      <SheetContent className="w-[400px] sm:w-[540px]">
        <SheetHeader>
          <SheetTitle>Backend Connection</SheetTitle>
          <SheetDescription>
            Configure how to connect to the RFdiffusion backend
          </SheetDescription>
        </SheetHeader>

        <div className="mt-6 space-y-6">
          {/* Mode Selection */}
          <div className="space-y-3">
            <Label>Choose your setup</Label>
            <div className="space-y-2">
              {modes.map((m) => (
                <button
                  key={m.id}
                  onClick={() => setMode(m.id)}
                  className={cn(
                    'w-full flex items-center gap-3 p-3 rounded-lg border transition-colors text-left',
                    mode === m.id
                      ? 'border-primary bg-accent'
                      : 'border-border hover:bg-accent/50'
                  )}
                >
                  <m.icon className="h-5 w-5 text-muted-foreground" />
                  <div>
                    <div className="font-medium text-sm">{m.label}</div>
                    <div className="text-xs text-muted-foreground">{m.description}</div>
                  </div>
                </button>
              ))}
            </div>
          </div>

          {/* URL Input */}
          <div className="space-y-2">
            <Label htmlFor="server-url">Server URL</Label>
            <Input
              id="server-url"
              placeholder={
                mode === 'runpod'
                  ? 'https://api.runpod.ai/v2/...'
                  : mode === 'local'
                  ? 'http://localhost:8000'
                  : 'https://your-server.com/api'
              }
              value={url}
              onChange={(e) => setUrl(e.target.value)}
            />
          </div>

          {/* Error Message */}
          {error && (
            <Alert variant="destructive">
              <XCircle className="h-4 w-4" />
              <AlertDescription>{error}</AlertDescription>
            </Alert>
          )}

          {/* Success Message */}
          {status === 'connected' && (
            <Alert>
              <CheckCircle className="h-4 w-4 text-green-600" />
              <AlertDescription>Connected successfully!</AlertDescription>
            </Alert>
          )}

          {/* Connect Button */}
          <Button
            onClick={handleConnect}
            disabled={!url || status === 'connecting'}
            className="w-full"
          >
            {status === 'connecting' ? (
              <>
                <Loader2 className="h-4 w-4 mr-2 animate-spin" />
                Connecting...
              </>
            ) : (
              'Test Connection'
            )}
          </Button>

          {/* Status */}
          <div className="flex items-center gap-2 text-sm">
            <span className="text-muted-foreground">Status:</span>
            {status === 'connected' ? (
              <span className="flex items-center gap-1 text-green-600">
                <span className="w-2 h-2 rounded-full bg-green-600" />
                Connected
              </span>
            ) : status === 'connecting' ? (
              <span className="flex items-center gap-1 text-amber-600">
                <span className="w-2 h-2 rounded-full bg-amber-600 animate-pulse" />
                Connecting
              </span>
            ) : status === 'failed' ? (
              <span className="flex items-center gap-1 text-destructive">
                <span className="w-2 h-2 rounded-full bg-destructive" />
                Failed
              </span>
            ) : (
              <span className="flex items-center gap-1 text-muted-foreground">
                <span className="w-2 h-2 rounded-full bg-muted-foreground" />
                Not connected
              </span>
            )}
          </div>
        </div>
      </SheetContent>
    </Sheet>
  );
}

export default ConnectionSheet;
```

**Step 3: Create index export**

Create file `src/components/connection/index.ts`:
```typescript
export { ConnectionSheet } from './ConnectionSheet';
```

**Step 4: Commit**

```bash
git add src/components/connection/
git commit -m "feat: create ConnectionSheet component replacing ConnectionModal"
```

---

### Task 3.2: Add Toaster to Layout

**Files:**
- Modify: `src/app/layout.tsx`

**Step 1: Read current layout.tsx**

Run:
```bash
cat src/app/layout.tsx
```

**Step 2: Add Toaster import and component**

Add to imports:
```typescript
import { Toaster } from '@/components/ui/sonner';
```

Add inside the body, after children:
```tsx
<Toaster position="bottom-right" />
```

**Step 3: Remove Material Symbols font link if present**

Remove any `<link>` tag referencing Material Symbols from the layout.

**Step 4: Verify build**

Run:
```bash
npm run build
```

Expected: Build succeeds

**Step 5: Commit**

```bash
git add src/app/layout.tsx
git commit -m "feat: add Toaster provider to layout, remove Material Symbols"
```

---

## Phase 4: Form Validation Schemas

### Task 4.1: Create Metal Form Schema

**Files:**
- Create: `src/lib/validations/metal-form.ts`

**Step 1: Create validations directory**

Run:
```bash
mkdir -p src/lib/validations
```

**Step 2: Create metal-form.ts**

Create file `src/lib/validations/metal-form.ts`:
```typescript
import { z } from 'zod';

export const metalFormSchema = z.object({
  // Metal Selection
  metal: z.enum([
    'Cu', 'Zn', 'Fe', 'Mn', 'Co', 'Ni', 'Ca', 'Mg',
    'Na', 'K', 'Cd', 'Hg', 'Pb', 'Mo', 'W'
  ]),
  coordination: z.number().min(2).max(8),
  geometry: z.enum([
    'tetrahedral',
    'square_planar',
    'octahedral',
    'trigonal_bipyramidal',
    'square_pyramidal',
    'linear',
    'trigonal_planar'
  ]).optional(),

  // Binding Site
  bindingSiteResidues: z.string().optional(),
  distanceConstraint: z.number().min(1.5).max(4.0).default(2.3),

  // Quality Settings
  qualityPreset: z.enum(['Quick', 'Balanced', 'High', 'Custom']).default('Balanced'),
  numTimesteps: z.number().min(10).max(500).default(50),
  stepScale: z.number().min(0.1).max(2.0).default(1.0),
  gamma0: z.number().min(0).max(1).default(0.0),

  // Output
  numDesigns: z.number().min(1).max(10).default(4),

  // Advanced
  usePotentials: z.boolean().default(true),
  potentialScale: z.number().min(0.1).max(10).default(1.0),
});

export type MetalFormValues = z.infer<typeof metalFormSchema>;

export const metalFormDefaults: MetalFormValues = {
  metal: 'Cu',
  coordination: 4,
  geometry: 'tetrahedral',
  bindingSiteResidues: '',
  distanceConstraint: 2.3,
  qualityPreset: 'Balanced',
  numTimesteps: 50,
  stepScale: 1.0,
  gamma0: 0.0,
  numDesigns: 4,
  usePotentials: true,
  potentialScale: 1.0,
};
```

**Step 3: Commit**

```bash
git add src/lib/validations/metal-form.ts
git commit -m "feat: create Zod schema for metal form validation"
```

---

### Task 4.2: Create Ligand Form Schema

**Files:**
- Create: `src/lib/validations/ligand-form.ts`

**Step 1: Create ligand-form.ts**

Create file `src/lib/validations/ligand-form.ts`:
```typescript
import { z } from 'zod';

export const ligandFormSchema = z.object({
  // Ligand Selection
  ligand: z.string().min(1, 'Ligand is required'),
  customSmiles: z.string().optional(),

  // Binding Site
  bindingPocket: z.enum(['auto', 'manual']).default('auto'),
  pocketResidues: z.string().optional(),
  pocketRadius: z.number().min(3).max(15).default(6),

  // Design Parameters
  scaffoldLength: z.object({
    min: z.number().min(30).max(500),
    max: z.number().min(30).max(500),
  }).refine(data => data.max >= data.min, {
    message: 'Max must be greater than or equal to min',
  }),

  // Quality Settings
  qualityPreset: z.enum(['Quick', 'Balanced', 'High', 'Custom']).default('Balanced'),
  numTimesteps: z.number().min(10).max(500).default(50),
  stepScale: z.number().min(0.1).max(2.0).default(1.0),
  gamma0: z.number().min(0).max(1).default(0.0),

  // Output
  numDesigns: z.number().min(1).max(10).default(4),

  // Advanced
  useConstraints: z.boolean().default(true),
  constraintWeight: z.number().min(0.1).max(10).default(1.0),
});

export type LigandFormValues = z.infer<typeof ligandFormSchema>;

export const ligandFormDefaults: LigandFormValues = {
  ligand: 'ATP',
  customSmiles: '',
  bindingPocket: 'auto',
  pocketResidues: '',
  pocketRadius: 6,
  scaffoldLength: { min: 80, max: 150 },
  qualityPreset: 'Balanced',
  numTimesteps: 50,
  stepScale: 1.0,
  gamma0: 0.0,
  numDesigns: 4,
  useConstraints: true,
  constraintWeight: 1.0,
};
```

**Step 2: Create validations index**

Create file `src/lib/validations/index.ts`:
```typescript
export * from './metal-form';
export * from './ligand-form';
```

**Step 3: Commit**

```bash
git add src/lib/validations/
git commit -m "feat: create Zod schema for ligand form validation"
```

---

## Phase 5: Integration (Summary Tasks)

The remaining phases involve integrating the new components with existing code. Each task follows the same pattern:

### Task 5.1: Update page.tsx to Use MainLayout
### Task 5.2: Connect Sidebar to Zustand Store
### Task 5.3: Migrate InterfaceMetalForm to React Hook Form
### Task 5.4: Migrate InterfaceLigandForm to React Hook Form
### Task 5.5: Migrate remaining task forms
### Task 5.6: Update AI components styling
### Task 5.7: Update binder components styling
### Task 5.8: Remove deprecated components
### Task 5.9: Final cleanup and testing

---

## Verification Checklist

After completing all phases, verify:

- [ ] `npm run build` succeeds
- [ ] `npm test` passes (70+ tests)
- [ ] `npm run dev` shows warm stone theme
- [ ] Sidebar navigation works
- [ ] Manual mode toggle switches between AI and forms
- [ ] Connection sheet opens and closes
- [ ] Toast notifications appear
- [ ] All forms validate correctly
- [ ] No Material Symbols icons remain (all Lucide)
- [ ] No console errors

---

## Rollback Plan

If issues arise:

1. Keep backup: `src/app/globals.css.backup`
2. Git history preserves all changes
3. Can revert individual commits
4. Original components not deleted until Phase 7

---

**End of Implementation Plan**
