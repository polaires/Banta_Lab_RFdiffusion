'use client';

import { Check, Circle, XCircle, Archive, Trash2 } from 'lucide-react';
import { useStore } from '@/lib/store';
import type { ProjectStatus } from '@/lib/chat-types';
import { cn } from '@/lib/utils';
import { ScrollArea } from '@/components/ui/scroll-area';

function formatRelativeTime(ts: number): string {
  const diffMs = Date.now() - ts;
  const diffMins = Math.floor(diffMs / 60000);
  const diffHours = Math.floor(diffMins / 60);
  const diffDays = Math.floor(diffHours / 24);

  if (diffMins < 1) return 'now';
  if (diffMins < 60) return `${diffMins}m`;
  if (diffHours < 24) return `${diffHours}h`;
  return `${diffDays}d`;
}

function StatusIcon({ status }: { status: ProjectStatus }) {
  switch (status) {
    case 'completed':
      return <Check className="h-3 w-3 text-green-600" />;
    case 'failed':
      return <XCircle className="h-3 w-3 text-destructive" />;
    case 'archived':
      return <Archive className="h-3 w-3 text-muted-foreground" />;
    case 'active':
    default:
      return <Circle className="h-3 w-3 text-primary fill-primary/20" />;
  }
}

export function ProjectSidebar() {
  const { projects, activeProjectId, setActiveProject, deleteProject } = useStore();

  // Sort by updatedAt descending
  const sorted = [...projects].sort((a, b) => b.updatedAt - a.updatedAt);

  if (sorted.length === 0) {
    return (
      <div className="text-xs text-muted-foreground px-2 py-4 text-center">
        No projects yet. Click &quot;New Design&quot; to start.
      </div>
    );
  }

  return (
    <ScrollArea className="h-full">
      <div className="space-y-0.5">
        {sorted.map((project) => {
          const isActive = project.id === activeProjectId;
          return (
            <button
              key={project.id}
              onClick={() => setActiveProject(project.id)}
              className={cn(
                'w-full flex items-center gap-2 px-2 py-2 rounded text-sm hover:bg-sidebar-accent text-left group transition-colors',
                isActive && 'bg-sidebar-accent'
              )}
            >
              <StatusIcon status={project.status} />
              <span className="flex-1 truncate text-foreground">{project.name}</span>
              <span className="text-[10px] text-muted-foreground shrink-0">
                {formatRelativeTime(project.updatedAt)}
              </span>
              <button
                onClick={(e) => {
                  e.stopPropagation();
                  deleteProject(project.id);
                }}
                className="opacity-0 group-hover:opacity-100 transition-opacity p-0.5 hover:text-destructive"
                title="Delete project"
              >
                <Trash2 className="h-3 w-3" />
              </button>
            </button>
          );
        })}
      </div>
    </ScrollArea>
  );
}
