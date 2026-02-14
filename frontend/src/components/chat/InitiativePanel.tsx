'use client';

import { useState, useEffect, useCallback, useRef } from 'react';
import { Sparkles, Dna, ArrowRight, Send } from 'lucide-react';
import { Button } from '@/components/ui/button';

// ---- Greetings (rotate on each mount) ----

const GREETINGS = [
  'What protein shall we create today?',
  'Ready to engineer something new?',
  'What binding challenge can I help with?',
  "Let's design a protein together.",
  'What would you like to build?',
];

// ---- Suggestion pills ----

const SUGGESTIONS = [
  { text: 'Redesign rubredoxin for terbium', prompt: 'Redesign 1BRF to bind terbium instead of iron' },
  { text: 'Terbium-citrate metalloenzyme', prompt: 'Design a scaffold for terbium with citrate' },
  { text: 'Light-switchable azobenzene dimer', prompt: 'Design an azobenzene interface dimer' },
  { text: 'Calcium sensor protein', prompt: 'Design a calcium sensing protein with exposed binding site' },
  { text: 'PQQ-binding catalytic site', prompt: 'Design a protein with PQQ cofactor for catalysis' },
  { text: 'Europium fluorescence probe', prompt: 'Design a europium binding protein for fluorescence sensing' },
  { text: 'Zinc-finger redesign', prompt: 'Redesign a zinc finger domain for improved binding' },
  { text: 'Gadolinium MRI contrast agent', prompt: 'Design a gadolinium binding protein for MRI contrast' },
  { text: 'De novo iron-sulfur cluster', prompt: 'Design a de novo protein with an iron-sulfur cluster' },
  { text: 'Lanthanide luminescence tag', prompt: 'Design a lanthanide binding tag for luminescence applications' },
];

function shuffle<T>(arr: T[]): T[] {
  const a = [...arr];
  for (let i = a.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [a[i], a[j]] = [a[j], a[i]];
  }
  return a;
}

interface InitiativePanelProps {
  onSuggestionClick: (prompt: string) => void;
  onGuidedDesign: () => void;
  onSend: (text: string) => void;
}

export function InitiativePanel({ onSuggestionClick, onGuidedDesign, onSend }: InitiativePanelProps) {
  // Defer randomization to client-side only to avoid SSR hydration mismatch
  const [greeting, setGreeting] = useState(GREETINGS[0]);
  const [pills, setPills] = useState(SUGGESTIONS.slice(0, 4));
  const [inputValue, setInputValue] = useState('');
  const textareaRef = useRef<HTMLTextAreaElement>(null);

  useEffect(() => {
    setGreeting(GREETINGS[Math.floor(Math.random() * GREETINGS.length)]);
    setPills(shuffle(SUGGESTIONS).slice(0, 4));
  }, []);

  // Auto-resize textarea
  const adjustHeight = useCallback(() => {
    const el = textareaRef.current;
    if (!el) return;
    el.style.height = 'auto';
    el.style.height = Math.min(el.scrollHeight, 150) + 'px';
  }, []);

  useEffect(() => {
    adjustHeight();
  }, [inputValue, adjustHeight]);

  const handleSubmit = useCallback(() => {
    const trimmed = inputValue.trim();
    if (!trimmed) return;
    setInputValue('');
    onSend(trimmed);
  }, [inputValue, onSend]);

  const handleKeyDown = useCallback(
    (e: React.KeyboardEvent) => {
      if (e.key === 'Enter' && !e.shiftKey) {
        e.preventDefault();
        handleSubmit();
      }
    },
    [handleSubmit]
  );

  return (
    <div className="flex-1 flex flex-col items-center justify-center px-6 py-12">
      <div className="w-full max-w-2xl space-y-8">
        {/* Greeting — large centered text like Claude */}
        <div className="text-center space-y-2">
          <div className="flex items-center justify-center gap-3 mb-4">
            <div className="p-2.5 bg-gradient-to-br from-primary to-primary/80 rounded-xl">
              <Sparkles className="h-6 w-6 text-primary-foreground" />
            </div>
          </div>
          <h1 className="text-3xl font-semibold text-foreground tracking-tight">
            {greeting}
          </h1>
        </div>

        {/* Input box — centered, prominent */}
        <div className="bg-card border border-border rounded-2xl shadow-sm overflow-hidden">
          <div className="px-4 pt-4 pb-2">
            <textarea
              ref={textareaRef}
              value={inputValue}
              onChange={(e) => setInputValue(e.target.value)}
              onKeyDown={handleKeyDown}
              placeholder="What would you like to design?"
              className="w-full bg-transparent text-foreground text-sm placeholder:text-muted-foreground focus:outline-none resize-none min-h-[44px] max-h-[150px]"
              rows={1}
            />
          </div>
          <div className="flex items-center justify-end px-4 pb-3">
            <Button
              onClick={handleSubmit}
              disabled={!inputValue.trim()}
              size="icon"
              className="h-9 w-9 rounded-xl"
            >
              <Send className="h-4 w-4" />
            </Button>
          </div>
        </div>

        {/* Suggestion pills — centered row */}
        <div className="flex flex-wrap justify-center gap-2">
          {pills.map((s) => (
            <button
              key={s.text}
              onClick={() => onSuggestionClick(s.prompt)}
              className="inline-flex items-center gap-1.5 px-3.5 py-2 bg-card hover:bg-muted border border-border hover:border-primary/30 rounded-full text-xs font-medium text-foreground transition-all"
            >
              {s.text}
            </button>
          ))}
        </div>

        {/* Guided Design button — card style like Claude's project cards */}
        <div className="flex justify-center">
          <button
            onClick={onGuidedDesign}
            className="flex items-center gap-4 px-6 py-4 bg-card hover:bg-muted border border-border hover:border-primary/30 rounded-2xl transition-all group max-w-sm w-full"
          >
            <div className="p-2.5 bg-primary/10 rounded-xl group-hover:bg-primary/15 transition-colors">
              <Dna className="h-5 w-5 text-primary" />
            </div>
            <div className="flex-1 text-left">
              <div className="font-semibold text-foreground text-sm">Guided Design</div>
              <div className="text-xs text-muted-foreground">Walk through design options step by step</div>
            </div>
            <ArrowRight className="h-4 w-4 text-muted-foreground group-hover:text-primary transition-colors" />
          </button>
        </div>
      </div>
    </div>
  );
}
