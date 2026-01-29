'use client';

import { Checkbox } from '@/components/ui/checkbox';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { Card } from '@/components/ui/card';
import { Eye } from 'lucide-react';
import type { PdbOutput, SequenceOutput } from '@/lib/pipeline-types';

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

  const allIds = items.map(item => item.id);
  const allSelected = allIds.every(id => selectedIds.includes(id));
  const noneSelected = selectedIds.length === 0;

  const toggleItem = (id: string) => {
    if (selectedIds.includes(id)) {
      onSelectionChange(selectedIds.filter(s => s !== id));
    } else {
      onSelectionChange([...selectedIds, id]);
    }
  };

  const selectAll = () => onSelectionChange(allIds);
  const selectNone = () => onSelectionChange([]);

  return (
    <div className="space-y-2">
      <div className="flex items-center justify-between">
        <Badge variant="secondary" className="text-xs">
          {selectedIds.length} of {items.length} selected
        </Badge>
        <div className="flex items-center gap-1">
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

      <div className="grid grid-cols-2 gap-2">
        {items.map(item => {
          const isSelected = selectedIds.includes(item.id);
          const isPdb = 'pdbContent' in item;
          const metrics = (item as PdbOutput).metrics;

          return (
            <Card
              key={item.id}
              className={`p-2.5 cursor-pointer transition-all ${
                isSelected
                  ? 'border-primary bg-primary/5'
                  : 'hover:border-primary/50'
              }`}
              onClick={() => toggleItem(item.id)}
            >
              <div className="flex items-start gap-2">
                <Checkbox
                  checked={isSelected}
                  onCheckedChange={() => toggleItem(item.id)}
                  className="mt-0.5"
                />
                <div className="flex-1 min-w-0">
                  <div className="flex items-center justify-between">
                    <span className="text-xs font-medium text-foreground truncate">
                      {item.label ?? item.id}
                    </span>
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
                  {metrics && Object.keys(metrics).length > 0 && (
                    <div className="flex flex-wrap gap-x-2 gap-y-0.5 mt-1">
                      {Object.entries(metrics).slice(0, 3).map(([key, val]) => (
                        <span key={key} className="text-[10px] text-muted-foreground">
                          {key}: <span className="font-mono">{typeof val === 'number' ? val.toFixed(2) : val}</span>
                        </span>
                      ))}
                    </div>
                  )}
                  {!isPdb && 'sequence' in item && (
                    <p className="text-[10px] text-muted-foreground font-mono truncate mt-0.5">
                      {(item as SequenceOutput).sequence.slice(0, 30)}...
                    </p>
                  )}
                </div>
              </div>
            </Card>
          );
        })}
      </div>
    </div>
  );
}
