'use client';

import { useMemo } from 'react';
import { Checkbox } from '@/components/ui/checkbox';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { Card } from '@/components/ui/card';
import { ScrollArea } from '@/components/ui/scroll-area';
import { Eye, Copy } from 'lucide-react';
import type { PdbOutput, SequenceOutput } from '@/lib/pipeline-types';
import { computeProtParam, type ProtParamResult } from '@/lib/protparam';

interface DesignSelectorProps {
  pdbOutputs?: PdbOutput[];
  sequences?: SequenceOutput[];
  selectedIds: string[];
  onSelectionChange: (ids: string[]) => void;
  onViewDesign?: (pdbContent: string) => void;
}

export function DesignSelector({
  pdbOutputs,
  sequences,
  selectedIds,
  onSelectionChange,
  onViewDesign,
}: DesignSelectorProps) {
  const items = pdbOutputs ?? sequences ?? [];
  if (items.length === 0) return null;

  // Compute ProtParam metrics for sequence items (memoized)
  const protparamCache = useMemo(() => {
    const cache = new Map<string, ProtParamResult>();
    for (const item of items) {
      if ('sequence' in item && (item as SequenceOutput).sequence) {
        cache.set(item.id, computeProtParam((item as SequenceOutput).sequence));
      }
    }
    return cache;
  }, [items]);

  const allIds = items.map(item => item.id);
  const allSelected = allIds.every(id => selectedIds.includes(id));
  const noneSelected = selectedIds.length === 0;

  // Count stable sequences for the badge
  const stableCount = Array.from(protparamCache.values()).filter(p => p.isStable).length;

  const toggleItem = (id: string) => {
    if (selectedIds.includes(id)) {
      onSelectionChange(selectedIds.filter(s => s !== id));
    } else {
      onSelectionChange([...selectedIds, id]);
    }
  };

  const selectAll = () => onSelectionChange(allIds);
  const selectNone = () => onSelectionChange([]);
  const selectStable = () => {
    const stableIds = items
      .filter(item => {
        const pp = protparamCache.get(item.id);
        return !pp || pp.isStable; // keep PDB items + stable sequences
      })
      .map(item => item.id);
    onSelectionChange(stableIds);
  };

  return (
    <div className="space-y-2">
      <div className="flex items-center justify-between">
        <div className="flex items-center gap-2">
          <Badge variant="secondary" className="text-xs">
            {selectedIds.length} of {items.length} selected
          </Badge>
          {protparamCache.size > 0 && (
            <span className="text-[10px] text-muted-foreground">
              {stableCount} stable (II&lt;40)
            </span>
          )}
        </div>
        <div className="flex items-center gap-1">
          {protparamCache.size > 0 && stableCount < items.length && (
            <Button
              variant="ghost"
              size="sm"
              onClick={selectStable}
              className="h-6 px-2 text-xs"
            >
              Stable Only
            </Button>
          )}
          <Button
            variant="ghost"
            size="sm"
            onClick={selectAll}
            disabled={allSelected}
            className="h-6 px-2 text-xs"
          >
            Select All
          </Button>
          <Button
            variant="ghost"
            size="sm"
            onClick={selectNone}
            disabled={noneSelected}
            className="h-6 px-2 text-xs"
          >
            Clear
          </Button>
        </div>
      </div>

      <ScrollArea className="max-h-[50vh]">
        <div className="grid grid-cols-2 gap-2">
          {items.map(item => {
            const isSelected = selectedIds.includes(item.id);
            const isPdb = 'pdbContent' in item;
            const metrics = (item as PdbOutput).metrics;
            const isSeq = !isPdb && 'sequence' in item;
            // PDB items from RF3 validation also carry a sequence
            const itemSequence = (item as any).sequence as string | undefined;
            const pp = protparamCache.get(item.id);

            return (
              <Card
                key={item.id}
                className={`p-2.5 cursor-pointer transition-all overflow-hidden ${
                  isSelected
                    ? 'border-primary bg-primary/5'
                    : pp && !pp.isStable
                      ? 'border-destructive/30 hover:border-destructive/50'
                      : 'hover:border-primary/50'
                }`}
                onClick={() => toggleItem(item.id)}
              >
                <div className="flex items-start gap-2 min-w-0">
                  <Checkbox
                    checked={isSelected}
                    onCheckedChange={() => toggleItem(item.id)}
                    className="mt-0.5 shrink-0"
                  />
                  <div className="flex-1 min-w-0">
                    <div className="flex items-center justify-between gap-1 min-w-0">
                      <span className="text-xs font-medium text-foreground truncate min-w-0">
                        {item.label ?? item.id}
                      </span>
                      <div className="flex items-center gap-0.5 shrink-0">
                        {pp && !pp.isStable && (
                          <Badge variant="destructive" className="text-[9px] h-4 px-1">
                            unstable
                          </Badge>
                        )}
                        {itemSequence && (
                          <Button
                            variant="ghost"
                            size="sm"
                            className="h-5 w-5 p-0"
                            title="Copy sequence"
                            onClick={e => {
                              e.stopPropagation();
                              navigator.clipboard.writeText(itemSequence);
                            }}
                          >
                            <Copy className="h-3 w-3" />
                          </Button>
                        )}
                        {isPdb && onViewDesign && (
                          <Button
                            variant="ghost"
                            size="sm"
                            className="h-5 w-5 p-0"
                            onClick={e => {
                              e.stopPropagation();
                              onViewDesign((item as PdbOutput).pdbContent);
                            }}
                          >
                            <Eye className="h-3 w-3" />
                          </Button>
                        )}
                      </div>
                    </div>
                    {/* PDB metrics */}
                    {metrics && Object.keys(metrics).length > 0 && (
                      <div className="flex flex-wrap gap-x-2 gap-y-0.5 mt-1">
                        {Object.entries(metrics).slice(0, 3).map(([key, val]) => (
                          <span key={key} className="text-[10px] text-muted-foreground">
                            {key}: <span className="font-mono">{typeof val === 'number' ? val.toFixed(2) : val}</span>
                          </span>
                        ))}
                      </div>
                    )}
                    {/* ProtParam metrics for sequences */}
                    {pp && (
                      <div className="flex flex-wrap gap-x-2 gap-y-0.5 mt-1">
                        <span className={`text-[10px] ${pp.isStable ? 'text-muted-foreground' : 'text-destructive'}`}>
                          II: <span className="font-mono">{pp.instabilityIndex.toFixed(1)}</span>
                        </span>
                        <span className="text-[10px] text-muted-foreground">
                          pI: <span className="font-mono">{pp.theoreticalPI.toFixed(1)}</span>
                        </span>
                        <span className="text-[10px] text-muted-foreground">
                          AI: <span className="font-mono">{pp.aliphaticIndex.toFixed(0)}</span>
                        </span>
                        <span className="text-[10px] text-muted-foreground">
                          <span className="font-mono">{pp.length}</span> aa
                        </span>
                      </div>
                    )}
                    {/* Sequence preview */}
                    {itemSequence && (
                      <p className="text-[10px] text-muted-foreground font-mono truncate mt-0.5 overflow-hidden">
                        {itemSequence.slice(0, 30)}...
                      </p>
                    )}
                  </div>
                </div>
              </Card>
            );
          })}
        </div>
      </ScrollArea>
    </div>
  );
}
