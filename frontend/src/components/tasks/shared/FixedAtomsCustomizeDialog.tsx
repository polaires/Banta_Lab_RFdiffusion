'use client';

import { useMemo, useCallback } from 'react';
import { Sparkles } from 'lucide-react';
import {
  Dialog,
  DialogContent,
  DialogHeader,
  DialogTitle,
  DialogDescription,
  DialogFooter,
} from '@/components/ui/dialog';
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from '@/components/ui/select';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { ScrollArea } from '@/components/ui/scroll-area';
import { cn } from '@/lib/utils';

interface FixedAtomsCustomizeDialogProps {
  open: boolean;
  onOpenChange: (open: boolean) => void;
  catalyticResidues: Array<{ chain: string; residue: number; name: string }>;
  fixedAtomTypes: Record<string, string>;
  onAtomTypeChange: (key: string, atomType: string) => void;
  onBatchAtomTypeChange?: (types: Record<string, string>) => void;
  excludedResidues: Set<string>;
  metalCoordinators: Array<{ chain: string; residue: number }>;
}

const ATOM_TYPE_OPTIONS = [
  { value: 'ALL', label: 'All', description: 'Fix all atoms (full fixation)' },
  { value: 'BKBN', label: 'Backbone', description: 'Fix backbone only (N, CA, C, O)' },
  { value: 'TIP', label: 'Tips', description: 'Fix sidechain tips only' },
  { value: '', label: 'None', description: 'No atoms fixed' },
];

export function FixedAtomsCustomizeDialog({
  open,
  onOpenChange,
  catalyticResidues,
  fixedAtomTypes,
  onAtomTypeChange,
  onBatchAtomTypeChange,
  excludedResidues,
  metalCoordinators,
}: FixedAtomsCustomizeDialogProps) {
  // Build a set of metal coordinator keys for quick lookup
  const metalCoordinatorKeys = useMemo(() => {
    return new Set(metalCoordinators.map((c) => `${c.chain}${c.residue}`));
  }, [metalCoordinators]);

  // Get suggested atom type for a residue
  const getSuggestedAtomType = useCallback(
    (key: string): string => {
      const isCoordinator = metalCoordinatorKeys.has(key);
      const isExcluded = excludedResidues.has(key);

      if (isCoordinator) {
        // Metal coordinators need full fixation for coordination geometry
        // Unless they're excluded (redesign mode), then backbone only
        return isExcluded ? 'BKBN' : 'ALL';
      }

      // Regular catalytic residues - backbone is standard for enzymes
      return 'BKBN';
    },
    [metalCoordinatorKeys, excludedResidues]
  );

  // Get role for a residue
  const getResidueRole = useCallback(
    (key: string): 'coordinator' | 'catalytic' => {
      return metalCoordinatorKeys.has(key) ? 'coordinator' : 'catalytic';
    },
    [metalCoordinatorKeys]
  );

  // Apply all suggestions
  const handleApplySuggestions = useCallback(() => {
    // Build all changes first to avoid state batching issues
    const newTypes: Record<string, string> = { ...fixedAtomTypes };
    for (const res of catalyticResidues) {
      const key = `${res.chain}${res.residue}`;
      const suggested = getSuggestedAtomType(key);
      newTypes[key] = suggested;
    }

    // Use batch update if available, otherwise fall back to individual updates
    if (onBatchAtomTypeChange) {
      onBatchAtomTypeChange(newTypes);
    } else {
      // Fallback: set all at once by calling with the merged object
      for (const res of catalyticResidues) {
        const key = `${res.chain}${res.residue}`;
        onAtomTypeChange(key, newTypes[key]);
      }
    }
  }, [catalyticResidues, fixedAtomTypes, getSuggestedAtomType, onAtomTypeChange, onBatchAtomTypeChange]);

  // Check if current values match suggestions
  const allMatchSuggestions = useMemo(() => {
    return catalyticResidues.every((res) => {
      const key = `${res.chain}${res.residue}`;
      const current = fixedAtomTypes[key] || 'BKBN';
      const suggested = getSuggestedAtomType(key);
      return current === suggested;
    });
  }, [catalyticResidues, fixedAtomTypes, getSuggestedAtomType]);

  if (catalyticResidues.length === 0) {
    return null;
  }

  return (
    <Dialog open={open} onOpenChange={onOpenChange}>
      <DialogContent className="sm:max-w-lg">
        <DialogHeader>
          <DialogTitle>Customize Fixed Atoms</DialogTitle>
          <DialogDescription>
            Set atom fixation mode for each catalytic residue
          </DialogDescription>
        </DialogHeader>

        <ScrollArea className="max-h-[400px] pr-4">
          <div className="space-y-2">
            {catalyticResidues.map((res) => {
              const key = `${res.chain}${res.residue}`;
              const isExcluded = excludedResidues.has(key);
              const role = getResidueRole(key);
              const suggested = getSuggestedAtomType(key);
              const currentValue = fixedAtomTypes[key] || 'BKBN';

              return (
                <div
                  key={key}
                  className={cn(
                    'flex items-center gap-3 p-3 rounded-lg border',
                    isExcluded
                      ? 'bg-orange-50/50 dark:bg-orange-950/10 border-orange-200/50 dark:border-orange-800/50 opacity-60'
                      : 'bg-card border-border'
                  )}
                >
                  {/* Residue label */}
                  <div className="flex-1 min-w-0">
                    <div className="flex items-center gap-2">
                      <span className="font-mono text-sm font-medium">
                        {res.chain}
                        {res.residue}
                        {res.name ? `:${res.name}` : ''}
                      </span>
                      <Badge
                        variant="default"
                        className={cn(
                          'text-[10px] px-1.5 py-0',
                          role === 'coordinator'
                            ? 'bg-orange-500/90 hover:bg-orange-500/90 text-white'
                            : 'bg-sky-500/90 hover:bg-sky-500/90 text-white'
                        )}
                      >
                        {role}
                      </Badge>
                      {isExcluded && (
                        <Badge
                          variant="outline"
                          className="text-[10px] px-1.5 py-0 border-yellow-400 text-yellow-600 dark:text-yellow-400"
                        >
                          redesign
                        </Badge>
                      )}
                    </div>
                  </div>

                  {/* Atom type selector */}
                  <Select
                    value={currentValue || 'NONE'}
                    onValueChange={(value) => onAtomTypeChange(key, value === 'NONE' ? '' : value)}
                    disabled={isExcluded}
                  >
                    <SelectTrigger className="w-28 h-8 text-xs">
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      {ATOM_TYPE_OPTIONS.map((opt) => {
                        // GLY doesn't have sidechain atoms, so 'TIP' is invalid
                        const isGly = res.name === 'GLY';
                        const isTipDisabled = isGly && opt.value === 'TIP';
                        return (
                          <SelectItem
                            key={opt.value || 'NONE'}
                            value={opt.value || 'NONE'}
                            className="text-xs"
                            disabled={isTipDisabled}
                          >
                            {opt.label}{isTipDisabled ? ' (N/A)' : ''}
                          </SelectItem>
                        );
                      })}
                    </SelectContent>
                  </Select>

                  {/* Suggested value hint */}
                  <div className="text-[10px] text-muted-foreground w-24 text-right whitespace-nowrap shrink-0">
                    Suggested: {suggested || 'None'}
                  </div>
                </div>
              );
            })}
          </div>
        </ScrollArea>

        <DialogFooter className="flex-col sm:flex-row gap-2 sm:gap-0">
          <Button
            type="button"
            variant="outline"
            size="sm"
            onClick={handleApplySuggestions}
            disabled={allMatchSuggestions}
            className="gap-1.5"
          >
            <Sparkles className="w-3.5 h-3.5" />
            Apply Suggestions
          </Button>

          <div className="flex-1" />

          <div className="flex gap-2">
            <Button
              type="button"
              variant="ghost"
              size="sm"
              onClick={() => onOpenChange(false)}
            >
              Cancel
            </Button>

            <Button
              type="button"
              size="sm"
              onClick={() => onOpenChange(false)}
            >
              Save
            </Button>
          </div>
        </DialogFooter>
      </DialogContent>
    </Dialog>
  );
}

export default FixedAtomsCustomizeDialog;
