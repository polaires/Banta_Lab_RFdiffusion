# Frontend UI Redesign with shadcn/ui

**Date:** 2026-01-17
**Status:** Approved
**Author:** AI Design Assistant + User collaboration

---

## Executive Summary

Complete frontend redesign migrating from custom components to shadcn/ui, implementing an AI-first interface with warm professional aesthetics. The redesign addresses consistency, visual design, and accessibility while supporting both AI-assisted and manual workflows.

---

## Goals

1. **Consistency & Maintainability** - Standardize all 66 components on shadcn/ui
2. **Visual Refresh** - Warm professional aesthetic (stone/zinc palette)
3. **Accessibility & UX** - Improved keyboard navigation, ARIA attributes, form validation
4. **AI-First Interface** - Chat-centric design with AI assistant as primary workflow

---

## Design Decisions

| Decision | Choice |
|----------|--------|
| Goal | Complete overhaul (consistency + visual + accessibility) |
| Aesthetic | Warm Professional (stone/zinc, Notion/Cal.com-like) |
| Icons | Lucide React (full migration from Material Symbols) |
| Forms | shadcn/ui Form + React Hook Form + Zod validation |
| Navigation | Sidebar with workflow stages and iteration tracking |
| AI Assistant | Main content panel (AI-first, chat-centric) |
| Manual Mode | Toggle in sidebar for power users |
| Protein Viewer | Collapsible right panel (400px default) |
| Connection Setup | Sheet component (slide-in from right) |
| Responsive | Desktop only (1200px minimum) |

---

## Layout Architecture

### Three-Column Layout (1200px+)

```
┌─────────────────────────────────────────────────────────────────────────┐
│  Logo                                        Connection Status  User ▾  │
├────────────┬────────────────────────────────────────┬───────────────────┤
│            │                                        │                   │
│  SIDEBAR   │           MAIN CONTENT                 │  PROTEIN VIEWER   │
│   240px    │           flex-1                       │   400px           │
│   fixed    │                                        │   collapsible     │
│            │                                        │                   │
└────────────┴────────────────────────────────────────┴───────────────────┘
```

### Sidebar Navigation

```
┌────────────────────────┐
│  Banta Lab             │  ← Logo + brand
│  RFdiffusion           │
├────────────────────────┤
│  + New Design          │  ← Primary action
├────────────────────────┤
│  CURRENT WORKFLOW      │
│  ○ Task Selection      │  ← Pending (stone-400)
│  ● RFdiffusion    ⟳2   │  ← Active + iteration badge
│  ○ ProteinMPNN         │
│  ○ Validation          │
├────────────────────────┤
│  HISTORY               │
│  ✓ Design #47   2h ago │
│  ✓ Design #46   1d ago │
├────────────────────────┤
│  ⚙ Settings            │
│  ⌘K AI Assistant       │
│  Manual Mode ○         │  ← Toggle switch
└────────────────────────┘
```

**Behaviors:**
- Completed/active stages are clickable
- Pending stages show tooltip "Complete previous stage first"
- Iteration badge (⟳2) expands to show run history
- Manual Mode toggle switches main content between AI chat and forms
- `⌘K` keyboard shortcut opens AI assistant

---

## Theme Configuration

### Warm Professional Color Palette

```css
:root {
  /* Backgrounds */
  --background: 28 25% 98%;           /* stone-50 #fafaf9 */
  --card: 0 0% 100%;                  /* white */
  --muted: 60 5% 96%;                 /* stone-100 */

  /* Foregrounds */
  --foreground: 24 10% 10%;           /* stone-900 #1c1917 */
  --muted-foreground: 25 6% 45%;      /* stone-500 #78716c */
  --card-foreground: 24 10% 10%;      /* stone-900 */

  /* Primary (buttons, active states) */
  --primary: 24 10% 15%;              /* stone-800 #292524 */
  --primary-foreground: 28 25% 98%;   /* stone-50 */

  /* Borders & Rings */
  --border: 24 6% 90%;                /* stone-200 #e7e5e4 */
  --ring: 24 6% 65%;                  /* stone-400 #a8a29e */

  /* Accent (hover states) */
  --accent: 60 5% 96%;                /* stone-100 #f5f5f4 */
  --accent-foreground: 24 10% 10%;    /* stone-900 */

  /* Semantic colors */
  --destructive: 0 84% 60%;           /* red-500 */
  --warning: 38 92% 50%;              /* amber-500 */
  --success: 142 76% 36%;             /* green-600 */

  /* Radius */
  --radius: 0.5rem;
}
```

### Typography

- **UI Font:** Inter (current, keep)
- **Monospace:** JetBrains Mono (for scientific values, parameters)
- **Scale:** 12px (xs), 14px (sm/base), 16px (lg), 18px (xl), 24px (2xl)

---

## Component Specifications

### AI Assistant Panel (Main Content - AI Mode)

**Structure:**
- Header with title and "Clear Chat" action
- Scrollable chat history
- Quick-start cards for common design tasks
- Multiple-choice option buttons for AI questions
- Evaluation cards for generated designs
- "Review Parameters" expandable section before submission
- Input field with send button

**Quick-Start Cards:**
- Protein Binder
- Enzyme Scaffold
- Binder Interface
- De Novo Design
- (More based on task types)

### Manual Mode Forms

**Form Architecture:**
- React Hook Form for state management
- Zod schemas for validation
- shadcn/ui Form components
- Grouped sections with clear hierarchy
- Advanced Options in collapsible sections

**Validation Approach:**
```typescript
const metalFormSchema = z.object({
  metal: z.enum(['Cu', 'Zn', 'Fe', 'Mn', ...]),
  coordination: z.number().min(2).max(6),
  geometry: z.enum(['tetrahedral', 'octahedral', ...]),
  num_designs: z.number().min(1).max(10),
  // ... other fields
});
```

### Protein Viewer Panel

**Specifications:**
- Default width: 400px
- Resizable via drag handle
- Collapsed width: 48px strip
- Auto-expand on new structure generation
- Keyboard shortcut: `V` to toggle

**Sections:**
1. Mol* 3D viewer (main area)
2. View controls (Spin, Reset, pLDDT coloring, etc.)
3. Structure info (residues, pLDDT avg, RMSD)
4. Confidence coloring legend

### Connection Sheet

**Trigger:** First launch, connection failure, or manual open

**Flow:**
1. Mode selection (RunPod Serverless / Traditional / Local Dev)
2. Server URL input
3. Test Connection button
4. Status display (Connected/Failed/Connecting)

### Toast Notifications

**Variants:**

| Type | Icon | Usage |
|------|------|-------|
| Success | `CheckCircle` | Job completed, connection established |
| Warning | `AlertTriangle` | Queue busy, slow response |
| Error | `XCircle` | Connection failed, job error |
| Info | `Info` | General notifications |

**Position:** Bottom-right, stacked

---

## shadcn/ui Components Required

```bash
# Core layout
npx shadcn@latest add card button separator scroll-area resizable

# Forms
npx shadcn@latest add form input label select slider switch checkbox radio-group

# Overlays
npx shadcn@latest add dialog sheet alert-dialog tooltip popover

# Feedback
npx shadcn@latest add toast alert progress badge

# Navigation
npx shadcn@latest add collapsible toggle tabs

# Data display
npx shadcn@latest add table avatar
```

---

## File Structure Changes

### New Files

```
src/
├── components/
│   ├── ui/                    # shadcn/ui components (auto-generated)
│   │   ├── button.tsx
│   │   ├── card.tsx
│   │   ├── form.tsx
│   │   ├── input.tsx
│   │   ├── select.tsx
│   │   ├── sheet.tsx
│   │   ├── slider.tsx
│   │   ├── toast.tsx
│   │   └── ...
│   ├── layout/
│   │   ├── MainLayout.tsx     # Three-column layout wrapper
│   │   ├── Sidebar.tsx        # Navigation sidebar
│   │   ├── ViewerPanel.tsx    # Collapsible protein viewer
│   │   └── HeaderBar.tsx      # Top bar (connection, user)
│   └── connection/
│       └── ConnectionSheet.tsx # Replaces ConnectionModal
├── lib/
│   ├── utils.ts               # cn() helper for class merging
│   └── validations/
│       ├── metal-form.ts      # Zod schema for metal form
│       ├── ligand-form.ts     # Zod schema for ligand form
│       └── ...                # Other form schemas
```

### Modified Files

| File | Changes |
|------|---------|
| `app/globals.css` | shadcn/ui CSS variables, warm theme |
| `app/layout.tsx` | Remove Material Symbols, add Toaster |
| `app/page.tsx` | Use MainLayout, remove old Header |
| `components/ai/AIChatPanel.tsx` | Redesign with shadcn/ui |
| `components/ai/InterviewMode.tsx` | Update styling |
| `components/tasks/*.tsx` | Migrate to React Hook Form |
| `components/tasks/shared/*.tsx` | Update shared components |

### Deprecated Files

| File | Replacement |
|------|-------------|
| `Header.tsx` | `layout/HeaderBar.tsx` + `layout/Sidebar.tsx` |
| `TaskPanel.tsx` | AI quick-start cards |
| `TaskSelector.tsx` | AI quick-start cards |
| `ConnectionModal.tsx` | `ConnectionSheet.tsx` |
| `NotificationToast.tsx` | shadcn/ui Toast |

---

## Implementation Phases

### Phase 1: Foundation
- [ ] Install shadcn/ui CLI and initialize
- [ ] Configure warm professional theme in globals.css
- [ ] Create `lib/utils.ts` with `cn()` helper
- [ ] Install React Hook Form + Zod
- [ ] Replace Material Symbols with Lucide React icons

### Phase 2: Layout Shell
- [ ] Create `MainLayout.tsx` (three-column structure)
- [ ] Build `Sidebar.tsx` with navigation
- [ ] Build `ViewerPanel.tsx` (collapsible)
- [ ] Create `HeaderBar.tsx`
- [ ] Integrate with existing Zustand store

### Phase 3: Connection & Notifications
- [ ] Create `ConnectionSheet.tsx`
- [ ] Migrate to shadcn/ui Toast
- [ ] Add Toaster provider to layout
- [ ] Update connection status display

### Phase 4: AI Assistant Panel
- [ ] Redesign `AIChatPanel.tsx` with shadcn/ui
- [ ] Create quick-start cards component
- [ ] Update interview mode components
- [ ] Add "Review Parameters" expandable section
- [ ] Style evaluation cards

### Phase 5: Manual Mode Forms
- [ ] Create Zod schemas for each form type
- [ ] Migrate `InterfaceMetalForm.tsx`
- [ ] Migrate `InterfaceLigandForm.tsx`
- [ ] Migrate `DeNovoForm.tsx`
- [ ] Migrate `ProteinBinderForm.tsx`
- [ ] Migrate `EnzymeForm.tsx`
- [ ] Migrate `SymmetricForm.tsx`
- [ ] Migrate `SmallMoleculeForm.tsx`
- [ ] Migrate `NucleicAcidForm.tsx`
- [ ] Migrate `RefinementForm.tsx`
- [ ] Update shared form components

### Phase 6: Data Visualization & Binder
- [ ] Update binder analysis components styling
- [ ] Migrate metrics displays
- [ ] Update heatmaps and charts
- [ ] Style results panels

### Phase 7: Cleanup & Polish
- [ ] Remove deprecated components
- [ ] Remove unused CSS classes from globals.css
- [ ] Final visual QA pass
- [ ] Accessibility audit (keyboard nav, ARIA)
- [ ] Performance check

---

## Success Criteria

1. **Visual Consistency:** All components use shadcn/ui with warm theme
2. **No Material Symbols:** All icons are Lucide React
3. **Form Validation:** All forms use React Hook Form + Zod
4. **Accessibility:** Full keyboard navigation, proper ARIA labels
5. **AI-First:** Default view is AI chat, manual mode is optional
6. **Performance:** No regression in load time or interaction speed

---

## Technical Notes

### React Hook Form Integration

```typescript
// Example form setup
import { useForm } from 'react-hook-form';
import { zodResolver } from '@hookform/resolvers/zod';
import { metalFormSchema, type MetalFormValues } from '@/lib/validations/metal-form';

export function InterfaceMetalForm() {
  const form = useForm<MetalFormValues>({
    resolver: zodResolver(metalFormSchema),
    defaultValues: {
      metal: 'Cu',
      coordination: 4,
      geometry: 'tetrahedral',
      num_designs: 4,
    },
  });

  const onSubmit = (data: MetalFormValues) => {
    // Submit to API
  };

  return (
    <Form {...form}>
      <form onSubmit={form.handleSubmit(onSubmit)}>
        {/* Form fields */}
      </form>
    </Form>
  );
}
```

### Class Merging Utility

```typescript
// lib/utils.ts
import { type ClassValue, clsx } from 'clsx';
import { twMerge } from 'tailwind-merge';

export function cn(...inputs: ClassValue[]) {
  return twMerge(clsx(inputs));
}
```

### Icon Migration Map

| Material Symbol | Lucide Replacement |
|-----------------|-------------------|
| `tune` | `SlidersHorizontal` |
| `expand_more` | `ChevronDown` |
| `info` | `Info` |
| `check_circle` | `CheckCircle` |
| `error` | `XCircle` |
| `warning` | `AlertTriangle` |
| `settings` | `Settings` |
| `science` | `FlaskConical` |
| `upload_file` | `Upload` |
| `download` | `Download` |

---

## References

- [shadcn/ui Documentation](https://ui.shadcn.com/)
- [React Hook Form](https://react-hook-form.com/)
- [Zod Validation](https://zod.dev/)
- [Lucide Icons](https://lucide.dev/)
- [Cradle.bio](https://www.cradle.bio/) - Workflow inspiration
- [Notion](https://notion.so/) - Sidebar navigation reference
- [Cal.com](https://cal.com/) - Warm professional aesthetic reference
