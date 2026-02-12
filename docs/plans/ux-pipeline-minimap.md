# Pipeline Minimap UX Design

## Status: Proposal (Feb 2026)

## Context

The Natural Language design pipeline has 11 steps displayed via `PipelineStepper` (a horizontal row of 44px icon circles connected by lines) inside `PipelineRunner`. When a user clicks a step icon, a single `StepCard` renders below the stepper with full details (progress, results, actions).

**Current pain points:**

1. The stepper is purely horizontal with 11 circles -- it overflows or compresses heavily on narrow panels (the AI assistant panel is ~380px wide).
2. Step labels are truncated to 11 characters and rendered at 10px -- hard to read.
3. Only one StepCard is visible at a time. Users lose context of what happened in earlier steps.
4. No animation or visual feedback beyond a spinning Loader2 icon.
5. Completed steps collapse to a single line -- reviewing past results requires clicking back through each step.
6. The separator between stepper and detail card is a hard visual break -- no sense of the minimap "zooming into" a step.

## 1. Pipeline Minimap Overview

### Option A: Horizontal Track (Subway Map) -- RECOMMENDED

A compact horizontal rail where each step is a small node (icon + short label) connected by a colored progress line. This is the closest evolution of the existing `PipelineStepper` but redesigned for density and polish.

```
 +---------------------------------------------------------------------------+
 |  [Parse]---[Resolve]---[Scaffold]---[Coord]---[Config]---[RFD3]---[...]  |
 |    (1)        (2)         (3)        (4)       (5)        (6)             |
 +---------------------------------------------------------------------------+
```

**Layout:** Flex row with `overflow-x-auto` and `scrollbar-hide`. Nodes are 32px circles (down from 44px) with 6px connecting lines. Labels sit below each node at 9px font. The active step auto-scrolls into the center of the viewport.

**Progress visualization:**
- Completed nodes: filled primary color, white check icon
- Active node: pulsing ring animation (`animate-pulse ring-2 ring-primary/40`), step icon visible
- Pending nodes: muted background, dimmed icon
- Connector lines: solid primary up to the active step, dashed muted after

**Handling 11 steps:**
- Nodes shrink to 28px on containers under 400px wide
- Horizontal scroll with momentum scrolling (`scroll-smooth`)
- A subtle gradient fade on left/right edges indicates scrollable content
- The last 2 steps (save_history, check_lessons) are automatic/non-interactive -- render them as smaller 20px dots without labels to reduce visual weight

**Tailwind sketch:**
```
// Container
<div className="flex items-center gap-0 overflow-x-auto scrollbar-hide px-2 py-1.5">
  {steps.map((step, i) => (
    <Fragment key={step.id}>
      <MinimapNode step={step} state={states[i]} isActive={i === activeIdx} />
      {i < steps.length - 1 && <MinimapConnector completed={states[i].status === 'completed'} />}
    </Fragment>
  ))}
</div>
```

### Option B: Vertical Sidebar

A compact vertical list that stays visible alongside the step detail. Better for the AI panel layout where vertical space is more abundant than horizontal.

```
  +--[ Minimap ]--+--[ Step Detail ]-----------+
  | (1) Parse     |                             |
  | (2) Resolve   |  [Backbone Generation]      |
  | (3) Scaffold  |  Status: Running (47%)      |
  |>(4) Config    |  ========================   |
  | (5) RFD3      |  Generating 4 designs...    |
  | (6) Scout     |                             |
  | (7) MPNN      |                             |
  | (8) RF3       |                             |
  | (9) Analysis  |                             |
  | (.) Save      |                             |
  | (.) Lessons   |                             |
  +---------------+-----------------------------+
```

**Layout:** 140px fixed-width sidebar with vertically stacked rows. Each row: 24px icon + truncated label. Connector is a 2px vertical line in the left margin. Step detail panel fills remaining width.

**Progress:** Vertical line fills downward with primary color as steps complete. Active row has a left border accent (`border-l-2 border-primary`).

**Handling 11 steps:** Vertical layout handles any count naturally. At 24px per row + 4px gap, 11 steps = 308px -- fits in most viewport heights. Automatic steps (save/lessons) get a smaller font and reduced opacity.

**Trade-off:** Requires a layout refactor from stacked (stepper above, detail below) to side-by-side. This changes the PipelineRunner component structure significantly.

### Option C: Compact Arc (Not Recommended)

Steps arranged in a semicircle above the detail area. Visually interesting but wastes space in a narrow panel and makes labels awkward (rotated text or radial labels).

**Verdict:** Rejected. The 380px panel width cannot accommodate a readable arc with 11 nodes. This option works better for 5-7 steps on wider layouts.

### Recommendation

**Option A (Horizontal Track)** as the primary layout, with a future path to Option B for users who pin the pipeline panel wider. Option A requires minimal structural changes (evolve `PipelineStepper` in place), while Option B requires a layout rewrite.

---

## 2. Step Icons

Each pipeline step gets a unique Lucide icon that communicates its purpose at a glance. The icons below are chosen from the existing Lucide set already imported in the project.

| # | Step ID | Name | Icon | Lucide Name | Rationale |
|---|---------|------|------|-------------|-----------|
| 1 | `parse_intent` | Parse Intent | Brain sparkle | `MessageSquare` (current) | NL parsing -- speech bubble fits |
| 2 | `resolve_structure` | Resolve Structure | Database box | `Search` (current) | Fetching/resolving a resource |
| 3 | `scaffold_search_nl` | Scaffold Search | Library search | `Database` (current) | Querying PDB database |
| 4 | `coordination_analysis_nl` | Coordination Analysis | Atom orbital | `Atom` (current) | Chemistry/coordination |
| 5 | `configure` | Configure | Gear sliders | `Settings` (current) | Parameter tuning |
| 6 | `rfd3_nl` | Backbone Generation | Building blocks | `Building2` (current) | Constructing structure |
| 7 | `scout_filter_nl` | Scout Filter | Shield with check | `ShieldCheck` | Filtering/validation gate |
| 8 | `mpnn_nl` | MPNN Sequence | DNA strand | `Dna` (current) | Sequence design |
| 9 | `rf3_nl` | RF3 Validation | Microscope | `Microscope` | Structure prediction/validation |
| 10 | `analysis` | Final Analysis | Bar chart | `BarChart3` (current) | Scoring and metrics |
| 11 | `save_history_nl` | Save History | Clock/history | `History` (current) | Persistence |
| 12 | `check_lessons_nl` | Check Lessons | Light bulb | `Lightbulb` (current) | Learning/insights |

**Changes from current:** Only two icons change:
- Scout Filter: `Shield` -> `ShieldCheck` (adds the checkmark to distinguish from RF3 validation, which also used `Shield`)
- RF3 Validation: `Shield` -> `Microscope` (validation is closer to "examining under a microscope")

**Icon rendering in minimap node:**
```tsx
// MinimapNode icon sizing
const iconSize = isAutoStep ? 'h-3 w-3' : 'h-4 w-4';
// When completed, replace step icon with Check
{state.status === 'completed' ? <Check className={iconSize} /> : <Icon className={iconSize} />}
```

---

## 3. Zoom Interaction

The core interaction model: clicking a minimap node "zooms into" that step, showing its full detail card. The minimap stays visible above as a navigation rail.

### 3.1 Minimap Behavior on Zoom

When a step is selected (clicked or auto-focused):

1. **Minimap stays at full width** -- it does not shrink or move. It remains the top navigation bar.
2. **Selected node gets a ring indicator** -- `ring-2 ring-offset-1 ring-primary` on the clicked node, like the current `isSelected` style.
3. **Detail card appears below** with a slide-down animation (`animate-in slide-in-from-top-2 duration-200`).
4. **Connector visual** -- a small downward-pointing triangle/caret appears below the selected minimap node, pointing at the detail card. This creates a visual "zoom path" from minimap to detail.

```
  [1]---[2]---[3]---[>4<]---[5]---[6]---[7]---[8]---[9]
                      |
                      V  (caret indicator)
  +--------------------------------------------------+
  |  Configure Parameters                             |
  |  Status: Review                                   |
  |  ...                                              |
  +--------------------------------------------------+
```

### 3.2 Expanded Step View

The detail card (currently `StepCard`) gets a refined layout:

```
  +--------------------------------------------------+
  |  [Icon]  Step Name              [Status Badge]    |
  |          Step description text                    |
  +--------------------------------------------------+
  |                                                    |
  |  [ Progress bar or Result preview ]                |
  |                                                    |
  |  [ Parameter editor (if paused + next step) ]      |
  |                                                    |
  +--------------------------------------------------+
  |  [<- Back]   [Continue ->]   [Skip]   [Retry]     |
  +--------------------------------------------------+
```

This is essentially the current `StepCard` layout, which works well. The key change is how it connects to the minimap above.

### 3.3 Navigation

- **Click any minimap node** to jump to that step's detail (already implemented in current code).
- **Keyboard:** Left/Right arrow keys move between steps when minimap is focused. Enter opens the detail card. Escape returns focus to the minimap.
- **Back to overview:** Click the same minimap node again to deselect (toggle behavior, already implemented).

### 3.4 Animation Transitions

When switching between steps:

1. Current detail card fades out (`animate-out fade-out duration-150`)
2. Minimap selection ring moves to new node (CSS transition on ring position, `transition-all duration-200`)
3. New detail card slides in (`animate-in fade-in slide-in-from-top-2 duration-200`)

For step completion:
1. Active node's pulsing ring stops
2. Checkmark icon replaces step icon with a scale-in animation (`animate-in zoom-in duration-300`)
3. Connector line to the next step fills with primary color (CSS transition `transition-colors duration-500`)
4. Next node's ring begins pulsing

These use Tailwind's built-in animation utilities from `tailwindcss-animate` (already in the project via shadcn/ui).

---

## 4. Active Step Focus

### 4.1 Auto-Scroll to Active Step

When the pipeline advances to a new step, the minimap auto-scrolls to center the active node:

```tsx
// Inside MinimapNode, when this node becomes active:
useEffect(() => {
  if (isActive && nodeRef.current) {
    nodeRef.current.scrollIntoView({
      behavior: 'smooth',
      block: 'nearest',
      inline: 'center',
    });
  }
}, [isActive]);
```

The detail card also auto-expands for the new active step (current behavior via `setSelectedStepIndex` in PipelineRunner).

### 4.2 Active Step Pulse

The running step gets a subtle pulse animation on its ring:

```tsx
// MinimapNode for running status
<div className={cn(
  'rounded-full flex items-center justify-center',
  isRunning && 'animate-pulse ring-2 ring-primary/30',
)}>
```

Additionally, a small progress arc can be rendered around the node circle using an SVG overlay:

```
  Before (0%):     During (60%):     Complete (100%):
    .---.             .---.              .---.
   |     |           |#####|            |#####|
   |  >  |           |## > |            |  v  |
   |     |           |     |            |#####|
    '---'             '---'              '---'
```

Implementation: An SVG `<circle>` with `stroke-dasharray` and `stroke-dashoffset` driven by `state.progress`. This provides a progress ring similar to a loading spinner but showing actual completion percentage.

```tsx
// ProgressRing component (SVG overlay on minimap node)
function ProgressRing({ progress, size = 32 }: { progress: number; size?: number }) {
  const r = (size - 4) / 2;
  const circumference = 2 * Math.PI * r;
  const offset = circumference - (progress / 100) * circumference;

  return (
    <svg className="absolute inset-0" width={size} height={size}>
      <circle
        cx={size / 2} cy={size / 2} r={r}
        fill="none"
        stroke="currentColor"
        strokeWidth={2}
        strokeDasharray={circumference}
        strokeDashoffset={offset}
        className="text-primary transition-all duration-300 -rotate-90 origin-center"
      />
    </svg>
  );
}
```

### 4.3 Pending Step Dimming

Steps that haven't started yet are visually suppressed:

```tsx
// MinimapNode for pending status
<div className={cn(
  'rounded-full flex items-center justify-center',
  'bg-muted/50 text-muted-foreground/50',
)}>
```

The connector lines after the active step use dashed styling:
```tsx
// MinimapConnector
<div className={cn(
  'flex-1 h-px mx-1',
  completed ? 'bg-primary' : 'bg-border border-dashed',
)} />
```

### 4.4 Completion Animation

When a step completes:

1. The node icon crossfades from the step icon to a Check icon
2. The node background transitions from `bg-primary/10` (active) to `bg-primary/15` (completed)
3. A brief scale bounce: `scale-110` for 200ms then back to `scale-100`
4. The connector line to the next step transitions from dashed to solid primary

For failed steps, the node turns red (`bg-destructive/10 text-destructive`) with a brief shake animation:
```css
@keyframes shake {
  0%, 100% { transform: translateX(0); }
  25% { transform: translateX(-2px); }
  75% { transform: translateX(2px); }
}
```

---

## 5. Component Architecture

### 5.1 Component Tree

```
PipelineRunner (orchestration, state management)
  |
  +-- PipelineMinimap (overview rail)
  |     |
  |     +-- MinimapNode (per step: icon, status, progress ring)
  |     +-- MinimapConnector (line between nodes)
  |     +-- MinimapCaret (indicator pointing to selected step's detail)
  |
  +-- PipelineStepDetail (zoomed-in view of selected step)
        |
        +-- StepHeader (icon, name, status badge, timing)
        +-- StepContent (context-dependent body)
        |     +-- StepProgress (progress bar + message, when running)
        |     +-- StepResultPreview / custom ResultPreview (when paused/completed)
        |     +-- DesignSelector (when supportsSelection + paused)
        |     +-- StepParameterEditor (next step params, when paused)
        |
        +-- StepActions (action buttons: continue, retry, skip, edit)
```

### 5.2 File Layout

```
frontend/src/components/pipeline/
  PipelineRunner.tsx          # Existing -- update to use new components
  PipelineMinimap.tsx         # NEW -- replaces PipelineStepper
  MinimapNode.tsx             # NEW -- single step node with progress ring
  MinimapConnector.tsx        # NEW -- connector line segment
  PipelineStepDetail.tsx      # NEW -- refactored from StepCard
  StepCard.tsx                # DEPRECATED -- replaced by PipelineStepDetail

  # Existing, unchanged:
  StepResultPreview.tsx
  StepParameterEditor.tsx
  DesignSelector.tsx
  IntentResultPreview.tsx
  ScaffoldSearchResultPreview.tsx
  CoordinationAnalysisPreview.tsx
  ScoutResultPreview.tsx
  SweepResultPreview.tsx
  LessonResultPreview.tsx
  WorkflowProgressCard.tsx
```

### 5.3 Key Interfaces

```tsx
// MinimapNode props
interface MinimapNodeProps {
  step: PipelineStepDefinition;
  state: StepRuntimeState;
  index: number;
  isActive: boolean;
  isSelected: boolean;
  isAutoStep: boolean;  // save_history, check_lessons
  onClick: () => void;
}

// MinimapConnector props
interface MinimapConnectorProps {
  status: 'completed' | 'active' | 'pending';
}

// PipelineMinimap props (replaces PipelineStepperProps)
interface PipelineMinimapProps {
  steps: PipelineStepDefinition[];
  stepStates: StepRuntimeState[];
  activeStepIndex: number;
  selectedStepIndex: number | null;
  onStepClick: (index: number) => void;
  /** Step IDs that are automatic (no user interaction) */
  autoStepIds?: Set<string>;
}

// PipelineStepDetail props (replaces StepCardProps)
interface PipelineStepDetailProps {
  step: PipelineStepDefinition;
  state: StepRuntimeState;
  isActive: boolean;
  canGoBack: boolean;
  nextStepSchema?: StepParameterField[];
  nextStepParams?: Record<string, unknown>;
  onNextStepParamsChange?: (params: Record<string, unknown>) => void;
  onConfirm: () => void;
  onRetry: () => void;
  onSkip: () => void;
  onGoBack?: () => void;
  onSelectionChange: (ids: string[]) => void;
  onViewDesign?: (pdbContent: string) => void;
}
```

### 5.4 State Management

No new state is required. The existing `usePipeline` hook provides all runtime state (`PipelineRuntimeState`). The `selectedStepIndex` state in `PipelineRunner` controls which step detail is shown.

The only new local state is inside `MinimapNode`:
- `nodeRef` for scroll-into-view
- No additional state -- all appearance is derived from `StepRuntimeState.status` and `StepRuntimeState.progress`

---

## 6. Responsive Behavior

### 6.1 Narrow Screens (< 400px panel width)

The AI assistant panel can be as narrow as 360px. At this width:

- **Minimap nodes shrink** from 32px to 24px circles
- **Labels are hidden** entirely (only icons visible) -- tooltip on hover shows the step name
- **Auto-step dots** (save/lessons) shrink to 16px
- **Horizontal scroll** with momentum, gradient fade hints on edges
- **Detail card** uses full width below the minimap

```tsx
// Responsive node sizing
const nodeSize = containerWidth < 400 ? 'w-6 h-6' : 'w-8 h-8';
const showLabels = containerWidth >= 400;
```

### 6.2 Medium Screens (400-600px)

- Full minimap with 32px nodes and labels
- Labels truncated to 8 characters
- Horizontal scroll for overflow

### 6.3 Wide Screens (> 600px)

- Full minimap with 36px nodes and full labels
- No horizontal scroll needed (11 nodes at 36px + gaps fit in ~600px)
- Consider Option B (vertical sidebar) as an alternative layout toggle

### 6.4 Mobile (< 360px)

- Minimap collapses to a single line showing `Step 5/11: Configure` with left/right navigation arrows
- Tapping the collapsed bar expands the full minimap as an overlay
- Detail card fills the full viewport

### 6.5 Keyboard Navigation

| Key | Action |
|-----|--------|
| `Left Arrow` | Select previous step (when minimap focused) |
| `Right Arrow` | Select next step |
| `Enter` / `Space` | Open selected step's detail card |
| `Escape` | Close detail card, return focus to minimap |
| `Tab` | Move focus between minimap and detail card actions |
| `C` | Continue (when detail card is focused and step is paused) |
| `S` | Skip (when step is optional and paused) |
| `R` | Retry (when step has failed) |

Implementation: A `useEffect` keydown listener on the minimap container with `tabIndex={0}`:

```tsx
function PipelineMinimap({ ... }: PipelineMinimapProps) {
  const containerRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    const el = containerRef.current;
    if (!el) return;

    const handleKeyDown = (e: KeyboardEvent) => {
      if (e.key === 'ArrowLeft' && selectedStepIndex !== null && selectedStepIndex > 0) {
        e.preventDefault();
        onStepClick(selectedStepIndex - 1);
      } else if (e.key === 'ArrowRight' && selectedStepIndex !== null && selectedStepIndex < steps.length - 1) {
        e.preventDefault();
        onStepClick(selectedStepIndex + 1);
      } else if (e.key === 'Escape') {
        onStepClick(-1); // deselect
      }
    };

    el.addEventListener('keydown', handleKeyDown);
    return () => el.removeEventListener('keydown', handleKeyDown);
  }, [selectedStepIndex, steps.length, onStepClick]);

  return (
    <div ref={containerRef} tabIndex={0} role="tablist" aria-label="Pipeline steps" ...>
      ...
    </div>
  );
}
```

Each `MinimapNode` renders as `role="tab"` with `aria-selected`, `aria-label`, and `aria-controls` pointing to the detail panel.

---

## 7. Migration Path

### Phase 1: MinimapNode + Progress Ring (replaces PipelineStepper)
1. Create `PipelineMinimap.tsx` with `MinimapNode` and `MinimapConnector` inline
2. Add SVG progress ring to running nodes
3. Replace `PipelineStepper` import in `PipelineRunner.tsx`
4. Add responsive sizing based on container width (ResizeObserver)

### Phase 2: StepCard Refinement (rename to PipelineStepDetail)
1. Rename `StepCard` to `PipelineStepDetail`
2. Add entrance/exit animations using `tailwindcss-animate`
3. Add the caret indicator connecting minimap to detail

### Phase 3: Polish
1. Completion animation (scale bounce + checkmark crossfade)
2. Failure shake animation
3. Keyboard navigation
4. ARIA attributes for accessibility
5. Mobile collapsed mode

---

## 8. Rejected Alternatives

| Alternative | Reason for Rejection |
|-------------|---------------------|
| Circular/arc layout | Too wide for 380px panel; labels unreadable |
| Tabbed interface | Loses at-a-glance pipeline overview |
| Accordion (all steps stacked) | Too much vertical space; loses minimap benefit |
| Fixed sidebar (always visible) | Requires layout restructure; panel too narrow for side-by-side |
| Breadcrumb-style | Only shows current position, not full pipeline |
| Timeline with timestamps | Adds clutter; timestamps not useful during pipeline execution |
