'use client';

import { useCallback, useMemo } from 'react';
import { Sparkles } from 'lucide-react';

interface AtomInfo {
  name: string;
  element: string;
}

interface AtomCheckboxGridProps {
  /** Label for the grid section */
  label: string;
  /** Optional description text */
  description?: string;
  /** Available atoms to display as checkboxes */
  atoms: AtomInfo[];
  /** Currently selected atom names (comma-separated in state, converted here) */
  selectedAtoms: string[];
  /** Callback when selection changes */
  onChange: (atoms: string[]) => void;
  /** Suggested atoms to highlight */
  suggestedAtoms?: string[];
  /** Button label for applying suggestions */
  suggestionButtonLabel?: string;
  /** Callback when apply suggestions is clicked */
  onApplySuggestions?: () => void;
  /** Whether suggestions have been overridden by user */
  isOverridden?: boolean;
  /** Hint text shown below the grid */
  hint?: string;
  /** Additional CSS classes */
  className?: string;
  /** Whether the component is disabled */
  disabled?: boolean;
}

// Shadcn-style colors - neutral with primary accent for selected
const UNCHECKED_STYLE = 'text-foreground border-border bg-card hover:bg-muted/50';
const CHECKED_STYLE = 'bg-primary border-primary text-primary-foreground';

export function AtomCheckboxGrid({
  label,
  description,
  atoms,
  selectedAtoms,
  onChange,
  suggestedAtoms = [],
  suggestionButtonLabel = 'Apply suggestions',
  onApplySuggestions,
  isOverridden = false,
  hint,
  className = '',
  disabled = false,
}: AtomCheckboxGridProps) {
  // Convert selected atoms to a Set for O(1) lookup
  const selectedSet = useMemo(() => new Set(selectedAtoms), [selectedAtoms]);
  const suggestedSet = useMemo(() => new Set(suggestedAtoms), [suggestedAtoms]);

  // Toggle a single atom
  const toggleAtom = useCallback((atomName: string) => {
    if (disabled) return;

    if (selectedSet.has(atomName)) {
      onChange(selectedAtoms.filter(a => a !== atomName));
    } else {
      onChange([...selectedAtoms, atomName]);
    }
  }, [selectedAtoms, selectedSet, onChange, disabled]);

  // Select all atoms
  const selectAll = useCallback(() => {
    if (disabled) return;
    onChange(atoms.map(a => a.name));
  }, [atoms, onChange, disabled]);

  // Clear all selections
  const clearAll = useCallback(() => {
    if (disabled) return;
    onChange([]);
  }, [onChange, disabled]);

  // Get color classes for an atom - shadcn style
  const getColorClasses = (isChecked: boolean) => {
    return isChecked ? CHECKED_STYLE : UNCHECKED_STYLE;
  };

  // Check if atom is suggested but not selected
  const isSuggested = (atomName: string) => {
    return suggestedSet.has(atomName) && !selectedSet.has(atomName);
  };

  return (
    <div className={`space-y-2 ${className}`}>
      {/* Header with inline actions */}
      <div className="flex items-center justify-between gap-2">
        <div className="flex items-center gap-2 min-w-0">
          <span className="text-sm font-medium text-foreground whitespace-nowrap">{label}</span>
          {/* Apply suggestions button - inline */}
          {onApplySuggestions && suggestedAtoms.length > 0 && (
            <button
              type="button"
              onClick={onApplySuggestions}
              disabled={disabled}
              className={`
                inline-flex items-center gap-1 px-2 py-0.5 text-[10px] font-medium rounded-full
                transition-colors whitespace-nowrap
                ${isOverridden
                  ? 'bg-muted text-muted-foreground hover:bg-muted/80'
                  : 'bg-primary/10 text-primary hover:bg-primary/20'
                }
                disabled:opacity-50 disabled:cursor-not-allowed
              `}
            >
              <Sparkles className="h-3 w-3" />
              {isOverridden ? 'Reset' : 'Auto'}
            </button>
          )}
        </div>
        <div className="flex items-center gap-2 flex-shrink-0">
          <button
            type="button"
            onClick={selectAll}
            disabled={disabled || atoms.length === 0}
            className="text-xs text-primary hover:text-primary/80 disabled:text-muted-foreground disabled:cursor-not-allowed"
          >
            All
          </button>
          <span className="text-muted-foreground">|</span>
          <button
            type="button"
            onClick={clearAll}
            disabled={disabled || selectedAtoms.length === 0}
            className="text-xs text-primary hover:text-primary/80 disabled:text-muted-foreground disabled:cursor-not-allowed"
          >
            None
          </button>
        </div>
      </div>
      {description && (
        <p className="text-xs text-muted-foreground">{description}</p>
      )}

      {/* Atom grid */}
      {atoms.length > 0 ? (
        <div className="flex flex-wrap gap-1.5">
          {atoms.map((atom) => {
            const isChecked = selectedSet.has(atom.name);
            const suggested = isSuggested(atom.name);

            return (
              <button
                key={atom.name}
                type="button"
                onClick={() => toggleAtom(atom.name)}
                disabled={disabled}
                className={`
                  relative px-2 py-1 text-xs font-medium rounded-md border transition-all
                  ${getColorClasses(isChecked)}
                  ${disabled ? 'opacity-50 cursor-not-allowed' : 'cursor-pointer'}
                  ${suggested ? 'ring-2 ring-primary/50 ring-offset-1' : ''}
                `}
                title={`${atom.name} (${atom.element})${suggested ? ' - Suggested' : ''}`}
              >
                {isChecked && (
                  <span className="mr-1">✓</span>
                )}
                {atom.name}
                {suggested && !isChecked && (
                  <span className="absolute -top-1 -right-1 w-2 h-2 bg-primary rounded-full" />
                )}
              </button>
            );
          })}
        </div>
      ) : (
        <p className="text-xs text-muted-foreground italic">No atoms available</p>
      )}

      {/* Suggestion button moved to header for compactness */}

      {/* Hint text */}
      {hint && (
        <p className="text-xs text-muted-foreground flex items-center gap-1">
          <span className="text-primary">💡</span>
          {hint}
        </p>
      )}

      {/* Selected count */}
      {selectedAtoms.length > 0 && (
        <p className="text-xs text-muted-foreground">
          {selectedAtoms.length} atom{selectedAtoms.length !== 1 ? 's' : ''} selected: {selectedAtoms.join(', ')}
        </p>
      )}
    </div>
  );
}

export default AtomCheckboxGrid;
