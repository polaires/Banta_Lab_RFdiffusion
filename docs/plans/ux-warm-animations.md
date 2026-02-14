# UX Warm Animations & Personality Design

**Date:** 2026-02-11
**Status:** Brainstorm / proposal
**Scope:** Frontend pipeline UI (`PipelineRunner`, `StepCard`, `PipelineStepper`, step previews)

---

## Problem Statement

The current pipeline UI is functional but feels clinical. Steps show generic "Running" badges with a spinning `Loader2` icon. Progress messages are technical (`Generating 4 backbone(s) via RFdiffusion3...`). The color palette is neutral stone/gray. There is no celebration on success, no personality between steps, and error states feel alarming rather than helpful.

The tool builds *proteins* — molecular machines of extraordinary beauty. The UI should reflect that wonder.

---

## 1. Step-Specific "Thinking" Animations

Each pipeline step gets a unique **verb phrase**, **icon**, and **CSS animation** that plays while the step is running. These replace the current generic `Loader2 animate-spin` treatment.

### Current state (all steps)

```tsx
// PipelineStepper.tsx line 22
<Loader2 className="h-5 w-5 animate-spin" />

// StepCard.tsx line 170
<Loader2 className="h-3.5 w-3.5 animate-spin" />
```

### Proposed step personalities

| Step | Verb Phrase | Lucide Icon | Animation | Duration | Notes |
|------|------------|-------------|-----------|----------|-------|
| **Parse Intent** | "Deciphering your vision..." | `Brain` | Gentle pulse + sparkle particles | 1.5s loop | Brain icon scales 1.0-1.08, tiny dots orbit around it |
| **Resolve Structure** | "Summoning molecular blueprints..." | `Telescope` | Slow zoom-in + subtle glow | 2s loop | Icon drifts inward, outer ring glows amber |
| **Scaffold Search** | "Scouting nature's library..." | `BookOpen` | Page-flip rotation on Y axis | 1.2s loop | Book flips open/closed, subtle paper-shuffle feel |
| **Coordination Analysis** | "Mapping the binding pocket..." | `Atom` | Orbital spin (electrons circling nucleus) | 2.5s loop | Central atom stays, 2-3 small dots orbit on ellipses |
| **Configure** | "Brewing the perfect recipe..." | `FlaskConical` | Bubbles rising from flask | 1.8s loop | 3-4 small circles float upward, fade at top |
| **Backbone Generation** | "Sculpting molecular architecture..." | `Hammer` | Building-block stack (translate Y) | 2s loop | Small blocks slide in from sides, stack upward |
| **Scout Filter** | "Separating wheat from chaff..." | `Filter` | Particles falling through funnel | 1.5s loop | Dots enter top, some pass through, some deflect |
| **MPNN Sequence** | "Spelling out the genetic code..." | `Dna` | Helix rotation + typewriter letters | 2s loop | Double helix rotates, A/T/C/G letters fade in below |
| **RF3 Validation** | "Stress-testing your creation..." | `Microscope` | Slow scan line + pulse | 2.2s loop | Horizontal line sweeps down icon, green pulse on pass |
| **Analysis** | "Crunching the numbers..." | `BarChart3` | Bars growing upward staggered | 1.5s loop | 3-4 bars animate height 0-100% with 100ms stagger |
| **Save History** | "Filing away wisdom..." | `BookMarked` | Book slides into shelf (translate X) | 1s once | Subtle rightward slide + opacity fade-in |
| **Check Lessons** | "Searching for patterns..." | `Lightbulb` | Slow brightening glow | 1.8s loop | Radial gradient opacity oscillates, warm amber color |

### CSS Keyframes (examples)

```css
/* Gentle pulse for Parse Intent (Brain) */
@keyframes thinking-pulse {
  0%, 100% { transform: scale(1); opacity: 0.9; }
  50% { transform: scale(1.08); opacity: 1; }
}

/* Bubble animation for Configure (Flask) */
@keyframes bubble-rise {
  0% { transform: translateY(0); opacity: 0.7; }
  50% { opacity: 1; }
  100% { transform: translateY(-12px); opacity: 0; }
}

/* Bar chart grow for Analysis */
@keyframes bar-grow {
  0% { transform: scaleY(0.2); }
  60% { transform: scaleY(1.05); }
  100% { transform: scaleY(1); }
}

/* Orbital spin for Coordination Analysis (Atom) */
@keyframes orbital-spin {
  0% { transform: rotate(0deg) translateX(8px) rotate(0deg); }
  100% { transform: rotate(360deg) translateX(8px) rotate(-360deg); }
}

/* Helix rotation for MPNN (DNA) */
@keyframes helix-rotate {
  0% { transform: rotateY(0deg); }
  100% { transform: rotateY(360deg); }
}

/* Scan line for RF3 Validation (Microscope) */
@keyframes scan-line {
  0% { transform: translateY(-100%); opacity: 0.3; }
  50% { opacity: 0.8; }
  100% { transform: translateY(100%); opacity: 0.3; }
}

/* Glow for Check Lessons (Lightbulb) */
@keyframes glow-breathe {
  0%, 100% { filter: drop-shadow(0 0 2px rgba(245, 158, 11, 0.3)); }
  50% { filter: drop-shadow(0 0 8px rgba(245, 158, 11, 0.6)); }
}
```

### Implementation approach

**Option A: Pure CSS (recommended for v1)**
Add animation keyframes to `globals.css` and apply them via Tailwind arbitrary classes or utility classes. Each step icon wrapper gets a conditional class when `status === 'running'`.

```tsx
// In PipelineStepper.tsx
const runningAnimation: Record<string, string> = {
  parse_intent: 'animate-thinking-pulse',
  resolve_structure: 'animate-zoom-glow',
  scaffold_search_nl: 'animate-page-flip',
  coordination_analysis_nl: 'animate-orbital',
  configure: 'animate-bubble',
  rfd3_nl: 'animate-build-stack',
  scout_filter_nl: 'animate-filter-fall',
  mpnn_nl: 'animate-helix-spin',
  rf3_nl: 'animate-scan-pulse',
  analysis: 'animate-bar-grow',
  save_history_nl: 'animate-file-slide',
  check_lessons_nl: 'animate-glow-breathe',
};
```

**Option B: Framer Motion (richer, but adds dependency)**
Each animation is a `<motion.div>` with step-specific `animate` props. Better for complex multi-element animations (orbital dots, bubbles, bars) but adds ~40KB to bundle.

**Option C: Lottie micro-animations (highest quality)**
Pre-made JSON animations from LottieFiles. Highest visual quality, zero CSS complexity, but requires `lottie-react` dependency and curating/commissioning 12 animation files.

**Recommendation:** Start with Option A for most steps. Use Option B (Framer Motion) selectively for the 3-4 most complex animations (orbital, bubbles, bars, filter particles) if the team already plans to use Framer Motion elsewhere.

---

## 2. Warm Color Palette

### Current palette analysis

The existing palette is stone/warm-gray (`hsl(24, 10%, 15%)` primary, `hsl(30, 6%, 98%)` background). This is already slightly warm — the hue is in the 24-30 range (orange-ish neutral). But the saturation is very low (6-10%), so it reads as gray.

### Proposed: "Warm Lab" palette

The idea: a university lab at golden hour. Amber light streaming through flasks. Sage green chalkboards. Warm wood surfaces. Rose chalk dust.

| Token | Current | Proposed | Hex | Usage |
|-------|---------|----------|-----|-------|
| `--background` | `30 6% 98%` | `36 30% 98%` | `#FDFBF7` | Page bg — very subtle cream |
| `--foreground` | `24 10% 10%` | `24 15% 12%` | `#1F1A16` | Text — warm near-black |
| `--card` | `0 0% 100%` | `40 30% 99%` | `#FEFCF9` | Card bg — warm white |
| `--primary` | `24 10% 15%` | `28 60% 30%` | `#7A4A1F` | Warm umber — primary actions |
| `--primary/10` (running) | stone-800/10 | `35 80% 55%` | `#D4922A` | Amber glow for active states |
| `--success` | `142 76% 36%` | `152 50% 40%` | `#339966` | Sage green — softer than pure green |
| `--destructive` | `0 84% 60%` | `4 60% 55%` | `#C45544` | Warm rose — less alarming than pure red |
| `--warning` | `38 92% 50%` | `38 85% 52%` | `#E09520` | Rich amber (nearly unchanged) |
| `--info` | `217 91% 60%` | `210 50% 55%` | `#5588AA` | Dusty blue — calmer |
| `--border` | `24 6% 90%` | `30 15% 88%` | `#E5DDD4` | Warm border — hint of parchment |
| `--muted` | `60 5% 96%` | `36 20% 95%` | `#F5F0E8` | Warm muted bg |
| `--muted-foreground` | `25 6% 45%` | `25 12% 45%` | `#7A6E63` | Warm gray text |
| `--ring` | `24 6% 65%` | `30 30% 55%` | `#A88A66` | Focus ring — warm gold |

### Semantic step colors (new tokens)

For the pipeline stepper, each step gets a subtle color accent in its "running" state:

| Step Category | Color Token | Hex | Used For |
|--------------|-------------|-----|----------|
| AI/Parse | `--step-ai` | `#8B6FC0` | Soft purple — intelligence |
| Structure | `--step-structure` | `#5B8A72` | Sage green — molecular |
| Configure | `--step-config` | `#C4883A` | Amber — brewing/tuning |
| Generate | `--step-generate` | `#D4922A` | Rich gold — creation |
| Validate | `--step-validate` | `#5588AA` | Dusty blue — analysis |
| Learn | `--step-learn` | `#A68B5B` | Warm bronze — wisdom |

### Alternative palette: "Deep Ocean"

If warm amber feels too earthy, an alternative is a deep ocean palette — indigo/teal primary, soft coral accents, pearl white backgrounds. This would feel more "high-tech" while still being warm.

| Token | Value | Hex |
|-------|-------|-----|
| `--background` | `200 20% 98%` | `#F7FAFB` |
| `--primary` | `200 60% 35%` | `#234E6C` |
| `--success` | `160 50% 40%` | `#338866` |
| `--accent` | `15 80% 60%` | `#E87850` |

### Alternative palette: "Greenhouse"

Botanical garden vibes. Deep forest green primary, warm terracotta accents, cream paper texture.

| Token | Value | Hex |
|-------|-------|-----|
| `--background` | `42 30% 97%` | `#FBF8F2` |
| `--primary` | `150 40% 25%` | `#264D3B` |
| `--success` | `95 50% 45%` | `#6B9944` |
| `--accent` | `18 70% 55%` | `#CC7744` |

**Recommendation:** "Warm Lab" palette. It builds on the existing stone hue (24-30) so the change is evolutionary, not revolutionary. The saturation bumps (6% to 15-30%) add warmth without becoming garish.

---

## 3. Micro-Interactions

### 3.1 Step Transitions

**Current:** Steps transition instantly. One step's card disappears, the next appears.

**Proposed:** Smooth handoff animation.

**Option A: Slide + Fade**
```css
@keyframes step-enter {
  from { opacity: 0; transform: translateY(8px); }
  to { opacity: 1; transform: translateY(0); }
}

@keyframes step-exit {
  from { opacity: 1; transform: translateY(0); }
  to { opacity: 0; transform: translateY(-8px); }
}

.step-entering { animation: step-enter 0.3s ease-out; }
.step-exiting { animation: step-exit 0.2s ease-in; }
```

**Option B: Crossfade**
Both cards visible briefly, outgoing fades to 0 while incoming fades from 0. Smoother but requires absolute positioning.

**Option C: Accordion (recommended)**
The StepCard already uses `<Collapsible>`. Enhance the collapse/expand with:
- Outgoing step collapses smoothly (300ms ease-out)
- Connector line between stepper dots animates from gray to primary color (200ms)
- Incoming step expands with a subtle spring (transform: scaleY, 400ms cubic-bezier)

### 3.2 Success Celebrations

When a step completes successfully:

**Option A: Checkmark bloom (recommended)**
The step's circle in the PipelineStepper transitions from running state to a brief "bloom":
```css
@keyframes check-bloom {
  0% { transform: scale(1); background: hsl(var(--primary)/0.1); }
  30% { transform: scale(1.3); background: hsl(var(--success)/0.2); }
  100% { transform: scale(1); background: hsl(var(--primary)/0.1); }
}
```
Duration: 500ms. The check icon fades in during the bloom.

**Option B: Confetti burst (reserved for pipeline completion)**
When the entire pipeline finishes, a subtle confetti burst:
- 20-30 small colored dots (palette colors) burst from the center
- Fall with gravity + slight horizontal drift
- Fade out after 1.5s
- Canvas-based or CSS pseudo-elements
- Only trigger on `pipeline.runtime.status === 'completed'`

**Option C: Ripple ring**
A single expanding ring (like a water ripple) emanates from the completed step's icon and fades:
```css
@keyframes ripple-success {
  0% { transform: scale(1); opacity: 0.4; }
  100% { transform: scale(2.5); opacity: 0; }
}
```

**Recommendation:** Checkmark bloom for individual steps (subtle, fast). Confetti burst for full pipeline completion only (reserved celebration).

### 3.3 Error States

**Current:** Red background, X icon, red text. Feels alarming.

**Proposed:** Warm, helpful error treatment.

- Replace pure red (`destructive`) with warm rose (`4 60% 55%`)
- Add a "What went wrong" collapsible with the technical error
- Show the suggestion text (already exists in `getActionableError()`) more prominently
- Icon: replace `X` with `AlertCircle` (less hostile)
- Add a subtle "don't worry" preamble: "This step hit a snag. Here's what happened..."
- Warm rose gradient background instead of flat red:

```css
.error-card {
  background: linear-gradient(135deg, hsl(4 60% 98%), hsl(4 60% 95%));
  border-color: hsl(4 60% 80%);
}
```

### 3.4 Progress Indicators

**Current:** A flat `<Progress>` bar (`h-1.5`) with instant updates.

**Proposed: "Breathing" progress**

**Option A: Pulsing glow**
The progress bar's filled portion has a subtle pulsing glow:
```css
.progress-alive {
  box-shadow: 0 0 8px 0 hsl(var(--primary)/0.3);
  animation: progress-breathe 2s ease-in-out infinite;
}

@keyframes progress-breathe {
  0%, 100% { box-shadow: 0 0 4px 0 hsl(var(--primary)/0.2); }
  50% { box-shadow: 0 0 12px 0 hsl(var(--primary)/0.4); }
}
```

**Option B: Shimmer sweep**
A highlight sweeps across the filled portion of the bar:
```css
.progress-shimmer::after {
  content: '';
  position: absolute;
  top: 0; left: -50%; width: 50%; height: 100%;
  background: linear-gradient(90deg, transparent, rgba(255,255,255,0.3), transparent);
  animation: shimmer 2s infinite;
}

@keyframes shimmer {
  0% { left: -50%; }
  100% { left: 150%; }
}
```

**Option C: Gradient bar**
The progress bar color transitions from amber (early) through gold (middle) to sage green (near completion):
```css
.progress-gradient {
  background: linear-gradient(90deg,
    hsl(38 85% 52%) 0%,
    hsl(45 80% 50%) 50%,
    hsl(152 50% 40%) 100%
  );
  background-size: 200% 100%;
  background-position: calc(100% - var(--progress) * 100%) 0;
}
```

**Recommendation:** Option B (shimmer) for the main progress bar, Option C (gradient) for the stepper's per-step mini progress bar. Both are achievable in pure CSS.

---

## 4. Personality Touches

### 4.1 Randomized Encouraging Messages

Between steps (during the "paused / review" state), show a brief encouraging message. These rotate randomly.

**Message pools by context:**

**After Parse Intent:**
- "Great choice. Let's bring this protein to life."
- "Understood. The molecular blueprint is forming."
- "An interesting challenge. Here's what I found."

**After Backbone Generation (success):**
- "Beautiful backbones. Nature would be proud."
- "The scaffolds are ready. Time to add the amino acids."
- "Solid architecture. Let's dress these in sequences."

**After Backbone Generation (mixed results):**
- "Some winners, some learning opportunities. That's science."
- "Not every backbone folds well — that's why we validate."

**After MPNN Sequence:**
- "The genetic code is written. Now for the final test."
- "Sequences designed. Let's see how they fold."

**After Analysis (high pass rate):**
- "Excellent results. You've got strong candidates here."
- "Multiple designs passed. Time to pick your favorites."

**After Analysis (low pass rate):**
- "Tough filtering, but quality over quantity."
- "A few strong candidates emerged. Sometimes one is all you need."

**Pipeline complete:**
- "Design complete. From idea to protein in minutes."
- "Your protein designs are ready for the wet lab."
- "Mission accomplished. Go make some proteins."

**Implementation:** A simple function that picks randomly from a pool:

```tsx
const encouragingMessages: Record<string, string[]> = {
  parse_intent: [
    'Great choice. Let\'s bring this protein to life.',
    'Understood. The molecular blueprint is forming.',
  ],
  rfd3_nl: [
    'Beautiful backbones. Nature would be proud.',
    'The scaffolds are ready. Time to add the amino acids.',
  ],
  // ...
};

function getEncouragingMessage(stepId: string): string {
  const pool = encouragingMessages[stepId];
  if (!pool || pool.length === 0) return '';
  return pool[Math.floor(Math.random() * pool.length)];
}
```

Display as a subtle italic line below the step summary in the paused state:

```tsx
{isPaused && state.result && (
  <p className="text-xs italic text-muted-foreground/70 mt-1">
    {getEncouragingMessage(step.id)}
  </p>
)}
```

### 4.2 Completion Celebration

When the pipeline finishes (`status === 'completed'`), show a **stats summary card** with personality:

```
+------------------------------------------+
|  Design Complete                         |
|                                          |
|  12 backbones generated                  |
|  8 passed scout filter (67%)             |
|  32 sequences designed                   |
|  5 designs passed analysis               |
|                                          |
|  Best: Design #3 — pLDDT 0.89, pTM 0.91 |
|                                          |
|  "Your protein designs are ready for     |
|   the wet lab."                          |
+------------------------------------------+
```

The card appears with the confetti animation described above.

### 4.3 Easter Eggs for Repeat Users

Track usage in `localStorage`:

- **5th pipeline run:** "You're getting the hang of this! 5 designs and counting."
- **10th run:** A subtle golden border on the pipeline card.
- **25th run:** "Power user detected. You've designed 25 proteins!"
- **50th run:** A special lab coat icon appears next to the pipeline name.
- **100th run:** "You've designed 100 proteins. Professor Banta would be proud."
- **Friday runs:** "Friday protein design session? Excellent choice."

**Implementation:** Increment a counter in localStorage on pipeline completion. Display Easter egg messages inline, never as popups.

---

## 5. Typography and Spacing

### Current

- Font: `Inter` with `cv11` and `ss01` features
- Very tight spacing: `py-2.5`, `gap-2`, `text-xs` / `text-[10px]`
- Cards are dense with minimal padding

### Proposed Changes

**Font:** Keep Inter — it is warm for a sans-serif (the `cv11` alternate makes the 'a' single-story, which is friendlier). Consider enabling `cv08` (old-style numerals) for metric displays to add character.

**Spacing — breathe more:**

| Element | Current | Proposed | Rationale |
|---------|---------|----------|-----------|
| StepCard padding | `px-3 py-2.5` | `px-4 py-3` | More breathing room |
| Step header gap | `gap-3` | `gap-3.5` | Slightly roomier |
| Progress bar height | `h-1.5` | `h-2` | More visible, less fragile |
| Badge text | `text-[10px]` | `text-[11px]` | Slightly more readable |
| Stepper icon size | `w-11 h-11` | `w-12 h-12` | Bigger tap target, more room for animation |
| Inter-step line | `h-px` (1px) | `h-0.5` (2px) | More visible connection |
| Card border radius | `0.5rem` | `0.625rem` (10px) | Slightly softer |

**Visual hierarchy improvements:**

1. **Step name in stepper:** Increase from `text-[10px]` to `text-[11px]`, add `tracking-wide` for the active step.

2. **Running step:** Add a subtle warm background to the entire StepCard when running:
```tsx
isRunning && 'bg-amber-50/50 dark:bg-amber-950/10'
```

3. **Completed step summary:** Currently truncated single line. Consider showing 2 lines max with a slightly bolder weight for key metrics.

4. **Section dividers:** Replace flat `<Separator />` with a gentler gradient line:
```css
.separator-warm {
  background: linear-gradient(90deg, transparent, hsl(var(--border)), transparent);
  height: 1px;
}
```

---

## 6. Pipeline Stepper Redesign

The stepper (horizontal icon row) is the pipeline's signature visual. Currently it is functional but flat.

### Proposed enhancements

**Animated connector lines:**
When a step completes, the connector line to the next step fills from left to right (like liquid flowing through a tube):
```css
@keyframes line-fill {
  from { background-size: 0% 100%; }
  to { background-size: 100% 100%; }
}

.connector-filling {
  background: linear-gradient(90deg, hsl(var(--primary)) 0%, hsl(var(--primary)) 100%);
  background-repeat: no-repeat;
  animation: line-fill 0.5s ease-out forwards;
}
```

**Step number badges:**
Small step number (1-11) shown as a tiny badge on each icon. This helps orientation in a long pipeline.

**Active step glow:**
The running step's circle gets a subtle amber glow ring that pulses:
```css
@keyframes active-glow {
  0%, 100% { box-shadow: 0 0 0 3px hsl(35 80% 55% / 0.15); }
  50% { box-shadow: 0 0 0 6px hsl(35 80% 55% / 0.08); }
}
```

**Tooltip on hover:**
Currently the step name is truncated to 12 chars. On hover, show a tooltip with the full name + status + elapsed time.

---

## 7. Implementation Plan

### Phase 1: Quick wins (1-2 days)
- Add warm color palette CSS variables to `globals.css`
- Add `thinking-pulse` and `glow-breathe` animations to `globals.css`
- Update `PipelineStepper` to use step-specific icons during running state
- Add encouraging messages to paused state in `StepCard`
- Bump spacing values (padding, icon sizes)

### Phase 2: Animations (2-3 days)
- Implement all 12 step-specific running animations
- Add checkmark bloom animation for step completion
- Add shimmer/gradient to progress bars
- Animate connector lines in stepper

### Phase 3: Celebrations (1 day)
- Add confetti animation for pipeline completion
- Add stats summary card on completion
- Add Easter egg counter in localStorage

### Phase 4: Polish (1 day)
- Error state redesign (warm rose, helpful messaging)
- Step transition animations (accordion spring)
- Tooltip hover on stepper icons
- Dark mode adjustments for all new colors

---

## 8. Technical Considerations

### Bundle size
- Pure CSS animations: 0 KB added to JS bundle
- Framer Motion: ~40KB gzipped (only if not already used)
- Lottie: ~15KB for player + per-animation JSON files
- Recommendation: pure CSS for v1

### Performance
- All animations use `transform` and `opacity` only (GPU-composited, no layout thrash)
- `will-change: transform` on animated elements
- Animations only render when steps are `running` (not on idle/completed steps)
- No `setInterval` — all CSS keyframe based

### Accessibility
- All animations respect `prefers-reduced-motion`:
```css
@media (prefers-reduced-motion: reduce) {
  .animate-thinking-pulse,
  .animate-bubble,
  .animate-orbital,
  /* ... all custom animations ... */ {
    animation: none !important;
  }
}
```
- Encouraging messages are decorative — screen readers can skip via `aria-hidden`
- Color changes maintain WCAG AA contrast ratios (all text > 4.5:1)
- Confetti uses `aria-hidden` and is non-interactive

### Dark mode
- All proposed colors include dark mode variants
- Warm cream backgrounds become warm dark browns (`hsl(30 15% 8%)`)
- Amber accents remain amber (high saturation reads well on dark)
- Rose error states lighten slightly for dark backgrounds

---

## 9. Files to Modify

| File | Changes |
|------|---------|
| `frontend/src/app/globals.css` | New CSS variables, animation keyframes, prefers-reduced-motion |
| `frontend/src/components/pipeline/PipelineStepper.tsx` | Step-specific icons when running, connector animation, glow ring |
| `frontend/src/components/pipeline/StepCard.tsx` | Encouraging messages, warm error styling, spacing bumps |
| `frontend/src/components/pipeline/PipelineRunner.tsx` | Completion celebration card, confetti, Easter eggs |
| `frontend/src/lib/pipelines/natural-language.ts` | (No changes — step IDs already stable) |
| `frontend/src/lib/step-personality.ts` | **NEW** — Map of step IDs to verb phrases, icons, animation classes, message pools |
