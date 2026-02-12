/**
 * Step personality — unique verb phrases, icons, animations, and encouraging messages
 * for each pipeline step. Replaces the generic Loader2 spinner with character.
 */

import {
  Brain,
  Telescope,
  BookOpen,
  Atom,
  FlaskConical,
  Hammer,
  Filter,
  Dna,
  Microscope,
  BarChart3,
  BookMarked,
  Lightbulb,
  type LucideIcon,
} from 'lucide-react';

export interface StepPersonality {
  /** Playful verb phrase shown while running */
  verbPhrase: string;
  /** Icon to show while running (overrides step.icon) */
  runningIcon: LucideIcon;
  /** CSS animation class applied to the icon wrapper when running */
  animationClass: string;
  /** Encouraging messages shown after step completes (randomly selected) */
  encouragingMessages: string[];
}

const personalities: Record<string, StepPersonality> = {
  parse_intent: {
    verbPhrase: 'Deciphering your vision...',
    runningIcon: Brain,
    animationClass: 'animate-thinking-pulse',
    encouragingMessages: [
      'Great choice. Let\'s bring this protein to life.',
      'Understood. The molecular blueprint is forming.',
      'An interesting challenge. Here\'s what I found.',
    ],
  },
  resolve_structure: {
    verbPhrase: 'Summoning molecular blueprints...',
    runningIcon: Telescope,
    animationClass: 'animate-zoom-glow',
    encouragingMessages: [
      'Structure locked in. Ready to build on it.',
      'Got the blueprint. Time to get creative.',
    ],
  },
  scaffold_search_nl: {
    verbPhrase: 'Scouting nature\'s library...',
    runningIcon: BookOpen,
    animationClass: 'animate-page-flip',
    encouragingMessages: [
      'Found some promising templates from nature.',
      'Nature\'s designs are a great starting point.',
    ],
  },
  coordination_analysis_nl: {
    verbPhrase: 'Mapping the binding pocket...',
    runningIcon: Atom,
    animationClass: 'animate-orbital',
    encouragingMessages: [
      'The coordination sphere is taking shape.',
      'Metal binding site mapped out nicely.',
    ],
  },
  configure: {
    verbPhrase: 'Brewing the perfect recipe...',
    runningIcon: FlaskConical,
    animationClass: 'animate-bubble',
    encouragingMessages: [
      'Parameters dialed in. Let\'s cook.',
      'Recipe locked and loaded.',
    ],
  },
  rfd3_nl: {
    verbPhrase: 'Sculpting molecular architecture...',
    runningIcon: Hammer,
    animationClass: 'animate-build',
    encouragingMessages: [
      'Beautiful backbones. Nature would be proud.',
      'The scaffolds are ready. Time to add amino acids.',
      'Solid architecture. Let\'s dress these in sequences.',
    ],
  },
  scout_filter_nl: {
    verbPhrase: 'Separating wheat from chaff...',
    runningIcon: Filter,
    animationClass: 'animate-filter-fall',
    encouragingMessages: [
      'The best candidates made it through.',
      'Quality over quantity. Smart filtering.',
    ],
  },
  mpnn_nl: {
    verbPhrase: 'Spelling out the genetic code...',
    runningIcon: Dna,
    animationClass: 'animate-helix-spin',
    encouragingMessages: [
      'The genetic code is written. Now for the final test.',
      'Sequences designed. Let\'s see how they fold.',
    ],
  },
  rf3_nl: {
    verbPhrase: 'Stress-testing your creation...',
    runningIcon: Microscope,
    animationClass: 'animate-scan-pulse',
    encouragingMessages: [
      'Validation complete. Let\'s see what made it.',
      'Structures predicted. Time to analyze.',
    ],
  },
  analysis: {
    verbPhrase: 'Crunching the numbers...',
    runningIcon: BarChart3,
    animationClass: 'animate-bar-grow',
    encouragingMessages: [
      'Excellent results. You\'ve got strong candidates here.',
      'Multiple designs passed. Time to pick your favorites.',
      'The numbers are in. Looking good.',
    ],
  },
  save_history_nl: {
    verbPhrase: 'Filing away wisdom...',
    runningIcon: BookMarked,
    animationClass: 'animate-file-slide',
    encouragingMessages: [
      'Results saved for future reference.',
    ],
  },
  check_lessons_nl: {
    verbPhrase: 'Searching for patterns...',
    runningIcon: Lightbulb,
    animationClass: 'animate-glow-breathe',
    encouragingMessages: [
      'Always learning, always improving.',
    ],
  },
};

/**
 * Pipeline completion messages — randomly selected when the full pipeline finishes.
 */
export const pipelineCompleteMessages = [
  'Design complete. From idea to protein in minutes.',
  'Your protein designs are ready for the wet lab.',
  'Mission accomplished. Go make some proteins.',
  'Another day, another protein designed.',
];

/**
 * Messages for mixed results (some passed, some didn't)
 */
export const mixedResultMessages = [
  'Some winners, some learning opportunities. That\'s science.',
  'Not every backbone folds well \u2014 that\u2019s why we validate.',
  'Tough filtering, but quality over quantity.',
  'A few strong candidates emerged. Sometimes one is all you need.',
];

export function getStepPersonality(stepId: string): StepPersonality | undefined {
  return personalities[stepId];
}

export function getEncouragingMessage(stepId: string): string {
  const personality = personalities[stepId];
  if (!personality?.encouragingMessages.length) return '';
  return personality.encouragingMessages[Math.floor(Math.random() * personality.encouragingMessages.length)];
}

export function getRandomMessage(pool: string[]): string {
  return pool[Math.floor(Math.random() * pool.length)];
}

/** Friendly step names for chat messages (avoids raw IDs like "parse_intent") */
const friendlyStepNames: Record<string, string> = {
  parse_intent: 'understanding your design',
  resolve_structure: 'finding your starting point',
  scaffold_search_nl: 'searching nature\'s library',
  coordination_analysis_nl: 'analyzing binding sites',
  configure: 'tuning the recipe',
  rfd3_nl: 'building protein structures',
  scout_filter_nl: 'finding the strongest',
  mpnn_nl: 'writing the genetic code',
  rf3_nl: 'testing your designs',
  analysis: 'evaluating results',
  save_history_nl: 'saving your work',
  check_lessons_nl: 'learning what works',
};

/**
 * Build a warm chat message for step completion.
 * Replaces cold "Step "parse_intent" completed: ..." with personality.
 */
export function getStepCompletionMessage(stepId: string, summary: string): string {
  const personality = personalities[stepId];
  const friendlyName = friendlyStepNames[stepId] ?? stepId.replace(/_/g, ' ');

  // Use an encouraging message if available, otherwise a generic warm phrase
  if (personality?.encouragingMessages.length) {
    const msg = personality.encouragingMessages[Math.floor(Math.random() * personality.encouragingMessages.length)];
    return `${msg} ${summary}`;
  }

  return `Done ${friendlyName} \u2014 ${summary}`;
}

/** Warm kickoff messages — randomly selected when a pipeline starts. */
export const pipelineStartMessages = [
  'On it. Let me work through this step by step.',
  'Great prompt. Firing up the design engine.',
  'Challenge accepted. Here we go.',
  'Interesting request \u2014 let me put the pieces together.',
];

/** Warm failure messages — randomly selected when a step fails. */
export function getStepFailureMessage(stepId: string, error: string): string {
  const friendlyName = friendlyStepNames[stepId] ?? stepId.replace(/_/g, ' ');
  const prefixes = [
    `Hit a snag while ${friendlyName}.`,
    `Something went wrong during ${friendlyName}.`,
    `Ran into trouble while ${friendlyName}.`,
  ];
  const prefix = prefixes[Math.floor(Math.random() * prefixes.length)];
  return `${prefix} ${error}`;
}

/**
 * Track pipeline usage for Easter eggs. Returns milestone message if applicable.
 */
export function trackPipelineRun(): string | null {
  if (typeof window === 'undefined') return null;
  const key = 'pipeline_run_count';
  const count = parseInt(localStorage.getItem(key) ?? '0', 10) + 1;
  localStorage.setItem(key, String(count));

  const day = new Date().getDay();
  if (day === 5 && count > 1) return 'Friday protein design session? Excellent choice.';

  switch (count) {
    case 5: return 'You\'re getting the hang of this! 5 designs and counting.';
    case 10: return '10 designs in the books. You\'re on a roll.';
    case 25: return 'Power user detected. 25 proteins designed!';
    case 50: return '50 proteins. You might need a bigger lab.';
    case 100: return '100 proteins designed. Professor Banta would be proud.';
    default: return null;
  }
}
