'use client';

import { FlaskConical, Plus, Check, X, Loader2 } from 'lucide-react';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { AtomTypeDropdown } from '@/components/AtomTypeDropdown';
import { cn } from '@/lib/utils';
import type { CatalyticSuggestion } from '@/lib/store';

// Color mapping for ligand/metal codes
const LIGAND_COLORS: Record<string, { bg: string; text: string }> = {
  // Metals - warm colors
  MG: { bg: 'bg-amber-100 dark:bg-amber-900/30', text: 'text-amber-700 dark:text-amber-400' },
  ZN: { bg: 'bg-slate-100 dark:bg-slate-800/50', text: 'text-slate-700 dark:text-slate-300' },
  CA: { bg: 'bg-stone-100 dark:bg-stone-800/50', text: 'text-stone-700 dark:text-stone-300' },
  FE: { bg: 'bg-red-100 dark:bg-red-900/30', text: 'text-red-700 dark:text-red-400' },
  MN: { bg: 'bg-purple-100 dark:bg-purple-900/30', text: 'text-purple-700 dark:text-purple-400' },
  CU: { bg: 'bg-orange-100 dark:bg-orange-900/30', text: 'text-orange-700 dark:text-orange-400' },
  TB: { bg: 'bg-indigo-100 dark:bg-indigo-900/30', text: 'text-indigo-700 dark:text-indigo-400' },
};

const DEFAULT_LIGAND_COLOR = { bg: 'bg-emerald-100 dark:bg-emerald-900/30', text: 'text-emerald-700 dark:text-emerald-400' };

function getLigandColor(code: string | undefined) {
  if (!code) return { bg: 'bg-muted', text: 'text-muted-foreground' };
  return LIGAND_COLORS[code.toUpperCase()] || DEFAULT_LIGAND_COLOR;
}

interface CatalyticSuggestionsPanelProps {
  suggestions: CatalyticSuggestion[];
  source: 'mcsa' | 'local' | 'none';
  loading: boolean;
  error: string | null;
  existingResidues: Array<{ chain: string; residue: number }>;
  onAdd: (suggestion: CatalyticSuggestion, atomType: string) => void;
  onRemove: (suggestion: CatalyticSuggestion) => void;
  onAddAll: (atomType: string) => void;
  onClose: () => void;
  onHoverResidue?: (suggestion: CatalyticSuggestion | null) => void;
  onRetry?: () => void;
  /** Default atom type for click-to-add (default: 'ALL') */
  defaultAtomType?: string;
}

export function CatalyticSuggestionsPanel({
  suggestions,
  source,
  loading,
  error,
  existingResidues,
  onAdd,
  onRemove,
  onAddAll,
  onClose,
  onHoverResidue,
  onRetry,
  defaultAtomType = 'ALL',
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
        <div className="flex items-center justify-center gap-3">
          <span className="text-sm text-destructive">{error}</span>
          {onRetry && (
            <Button variant="outline" size="sm" onClick={onRetry}>
              Retry
            </Button>
          )}
        </div>
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
              <button
                type="button"
                key={`${s.chain}${s.residue}`}
                className={`flex items-center justify-between p-2 rounded-md border transition-colors cursor-pointer text-left ${
                  added
                    ? 'bg-primary/10 border-primary/30 hover:bg-primary/20'
                    : 'bg-card border-border hover:bg-accent/50'
                }`}
                onClick={() => {
                  if (added) {
                    onRemove(s);
                  } else {
                    onAdd(s, defaultAtomType);
                  }
                }}
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
                    <div className="flex items-center gap-1.5">
                      <span className="text-xs font-medium truncate">
                        {s.name} {s.chain}{s.residue}
                      </span>
                      {s.ligandCode && (
                        <span className={cn(
                          "text-[9px] px-1 py-0.5 rounded font-medium",
                          getLigandColor(s.ligandCode).bg,
                          getLigandColor(s.ligandCode).text
                        )}>
                          {s.ligandCode}
                        </span>
                      )}
                    </div>
                    {s.role && (
                      <div className="text-[10px] text-muted-foreground truncate">
                        {s.role}
                      </div>
                    )}
                  </div>
                </div>
                {added ? (
                  <Check className="h-4 w-4 text-primary shrink-0" />
                ) : (
                  <Plus className="h-4 w-4 text-muted-foreground shrink-0" />
                )}
              </button>
            );
          })}
        </div>
      </div>

      {/* Footer */}
      <div className="px-4 py-2 border-t text-xs text-muted-foreground flex justify-between">
        <span>
          Source: {source === 'mcsa' ? 'M-CSA (curated)' : source === 'local' ? 'Binding Pocket (local)' : 'None'}
        </span>
        <span className="flex items-center gap-3">
          <span className="flex items-center gap-1">
            <div className="w-2 h-2 rounded-full bg-blue-500" /> M-CSA
          </span>
          <span className="flex items-center gap-1">
            <div className="w-2 h-2 rounded-full bg-orange-500" /> Local
          </span>
        </span>
      </div>
    </div>
  );
}
