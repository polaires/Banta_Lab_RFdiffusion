'use client';

import { Card } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { Button } from '@/components/ui/button';
import { ScrollArea } from '@/components/ui/scroll-area';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Eye, Copy } from 'lucide-react';
import type { StepResult } from '@/lib/pipeline-types';

interface StepResultPreviewProps {
  result: StepResult;
  onSelectDesign?: (pdbContent: string) => void;
}

export function StepResultPreview({ result, onSelectDesign }: StepResultPreviewProps) {
  const hasStructures = result.pdbOutputs && result.pdbOutputs.length > 0;
  const hasSequences = result.sequences && result.sequences.length > 0;

  // Only show Details tab when data has meaningful content the user cares about.
  // Hide it when it's just a few scalar metadata fields duplicating the summary.
  const dataEntries = result.data ? Object.entries(result.data).filter(
    ([, v]) => typeof v !== 'object' || v === null
  ) : [];
  const hasData = dataEntries.length > 0 && (
    // Always show if it's the only content
    (!hasStructures && !hasSequences) ||
    // Show if there are more than 3 meaningful entries (worth a tab)
    dataEntries.length > 3
  );

  const tabCount = [hasStructures, hasSequences, hasData].filter(Boolean).length;

  // If only a summary, show it directly
  if (tabCount === 0) {
    return (
      <div className="text-sm text-foreground">{result.summary}</div>
    );
  }

  // Determine default tab
  const defaultTab = hasStructures ? 'structures' : hasSequences ? 'sequences' : 'data';

  return (
    <div className="space-y-2 min-w-0 overflow-hidden">
      <p className="text-xs text-muted-foreground">{result.summary}</p>

      {tabCount === 1 && hasStructures && (
        <StructuresGrid
          outputs={result.pdbOutputs!}
          onSelectDesign={onSelectDesign}
        />
      )}

      {tabCount === 1 && hasSequences && (
        <SequencesList sequences={result.sequences!} />
      )}

      {tabCount === 1 && hasData && (
        <DataTable data={result.data!} />
      )}

      {tabCount > 1 && (
        <Tabs defaultValue={defaultTab}>
          <TabsList className="h-8">
            {hasStructures && (
              <TabsTrigger value="structures" className="text-xs h-7">
                Structures ({result.pdbOutputs!.length})
              </TabsTrigger>
            )}
            {hasSequences && (
              <TabsTrigger value="sequences" className="text-xs h-7">
                Sequences ({result.sequences!.length})
              </TabsTrigger>
            )}
            {hasData && (
              <TabsTrigger value="data" className="text-xs h-7">
                Details
              </TabsTrigger>
            )}
          </TabsList>
          {hasStructures && (
            <TabsContent value="structures" className="mt-2">
              <StructuresGrid
                outputs={result.pdbOutputs!}
                onSelectDesign={onSelectDesign}
              />
            </TabsContent>
          )}
          {hasSequences && (
            <TabsContent value="sequences" className="mt-2">
              <SequencesList sequences={result.sequences!} />
            </TabsContent>
          )}
          {hasData && (
            <TabsContent value="data" className="mt-2">
              <DataTable data={result.data!} />
            </TabsContent>
          )}
        </Tabs>
      )}
    </div>
  );
}

function StructuresGrid({
  outputs,
  onSelectDesign,
}: {
  outputs: StepResult['pdbOutputs'] & object;
  onSelectDesign?: (pdbContent: string) => void;
}) {
  return (
    <div className="grid grid-cols-2 gap-2">
      {outputs.map(out => (
        <Card key={out.id} className="p-2.5">
          <div className="flex items-center justify-between mb-1">
            <span className="text-xs font-medium text-foreground truncate">
              {out.label}
            </span>
            {onSelectDesign && (
              <Button
                variant="ghost"
                size="sm"
                className="h-6 px-1.5"
                onClick={() => onSelectDesign(out.pdbContent)}
              >
                <Eye className="h-3 w-3 mr-1" />
                <span className="text-[10px]">View</span>
              </Button>
            )}
          </div>
          {out.metrics && Object.keys(out.metrics).length > 0 && (
            <div className="flex flex-wrap gap-x-2 gap-y-0.5">
              {Object.entries(out.metrics).slice(0, 4).map(([key, val]) => (
                <span key={key} className="text-[10px] text-muted-foreground">
                  {key}: <span className="font-mono">{typeof val === 'number' ? val.toFixed(2) : val}</span>
                </span>
              ))}
            </div>
          )}
        </Card>
      ))}
    </div>
  );
}

function SequencesList({ sequences }: { sequences: StepResult['sequences'] & object }) {
  return (
    <ScrollArea className="max-h-[40vh]">
      <div className="space-y-1.5 min-w-0">
        {sequences.map(seq => (
          <Card key={seq.id} className="p-2.5 overflow-hidden">
            <div className="flex items-center justify-between mb-1 min-w-0">
              <span className="text-xs font-medium text-foreground truncate min-w-0">
                {seq.label ?? seq.id}
              </span>
              <div className="flex items-center gap-1 shrink-0">
                {seq.score !== undefined && (
                  <Badge variant="secondary" className="text-[10px] h-5 px-1.5">
                    {seq.score.toFixed(2)}
                  </Badge>
                )}
                <Button
                  variant="ghost"
                  size="sm"
                  className="h-5 w-5 p-0"
                  onClick={() => navigator.clipboard.writeText(seq.sequence)}
                >
                  <Copy className="h-3 w-3" />
                </Button>
              </div>
            </div>
            <p className="text-[10px] text-muted-foreground font-mono break-all leading-relaxed overflow-hidden">
              {seq.sequence.length > 120
                ? `${seq.sequence.slice(0, 120)}...`
                : seq.sequence
              }
            </p>
          </Card>
        ))}
      </div>
    </ScrollArea>
  );
}

function DataTable({ data }: { data: Record<string, unknown> }) {
  const entries = Object.entries(data).filter(
    ([, v]) => typeof v !== 'object' || v === null
  );

  if (entries.length === 0) return null;

  return (
    <div className="space-y-1">
      {entries.map(([key, val]) => (
        <div key={key} className="flex items-center justify-between text-xs py-1 border-b border-border last:border-0">
          <span className="text-muted-foreground">{key}</span>
          <span className="font-medium text-foreground font-mono">
            {typeof val === 'number' ? val.toFixed(3) : String(val ?? 'N/A')}
          </span>
        </div>
      ))}
    </div>
  );
}
