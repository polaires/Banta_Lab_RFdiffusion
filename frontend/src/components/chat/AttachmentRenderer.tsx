'use client';

import type { ChatAttachment } from '@/lib/chat-types';
import { PipelineCard } from './PipelineCard';

interface AttachmentRendererProps {
  attachment: ChatAttachment;
  onSelectDesign?: (pdbContent: string) => void;
}

/**
 * Dispatcher that renders inline components based on attachment type.
 * Currently supports pipeline_progress and quick_start.
 * Other types render a placeholder until the legacy components are adapted.
 */
export function AttachmentRenderer({ attachment, onSelectDesign }: AttachmentRendererProps) {
  switch (attachment.type) {
    case 'pipeline_progress':
      return <PipelineCard pipelineId={attachment.pipelineId} onSelectDesign={onSelectDesign} />;

    case 'quick_start':
      // Quick start cards are rendered by ChatPanel directly, not here
      return null;

    case 'step_result':
      return (
        <div className="bg-muted/50 rounded-xl p-3 border border-border text-xs text-muted-foreground">
          Step <span className="font-medium text-foreground">{attachment.stepId}</span> completed: {attachment.result.summary}
        </div>
      );

    case 'interview':
      // Interviews are rendered inline by ChatPanel which manages their state
      return null;

    case 'preference_summary':
      return (
        <div className="bg-muted rounded-xl p-4 border border-border space-y-2">
          <div className="text-xs font-medium text-muted-foreground uppercase tracking-wider">Design Configuration</div>
          {Object.entries(attachment.preferences).map(([key, val]) => (
            <div key={key} className="flex justify-between text-sm">
              <span className="text-muted-foreground capitalize">{key.replace(/_/g, ' ')}</span>
              <span className="font-medium text-foreground">{String(val)}</span>
            </div>
          ))}
        </div>
      );

    case 'design_gallery':
      return (
        <div className="bg-muted rounded-xl p-4 border border-border">
          <div className="text-xs font-medium text-muted-foreground uppercase tracking-wider mb-2">
            {attachment.designs.length} Design{attachment.designs.length !== 1 ? 's' : ''} Generated
          </div>
          <div className="space-y-1">
            {attachment.designs.slice(0, 5).map((d) => (
              <button
                key={d.id}
                onClick={() => d.pdbContent && onSelectDesign?.(d.pdbContent)}
                className="w-full text-left px-3 py-2 rounded-lg hover:bg-card text-sm transition-colors border border-transparent hover:border-border"
              >
                <span className="font-medium text-foreground">{d.name}</span>
                {d.metrics && (
                  <span className="text-xs text-muted-foreground ml-2">
                    {Object.entries(d.metrics).slice(0, 3).map(([k, v]) => `${k}: ${v}`).join(' | ')}
                  </span>
                )}
              </button>
            ))}
            {attachment.designs.length > 5 && (
              <div className="text-xs text-muted-foreground text-center py-1">
                +{attachment.designs.length - 5} more
              </div>
            )}
          </div>
        </div>
      );

    case 'evaluation_card':
      return (
        <div className="bg-muted rounded-xl p-4 border border-border text-sm text-muted-foreground">
          Evaluation complete
        </div>
      );

    default:
      return null;
  }
}
