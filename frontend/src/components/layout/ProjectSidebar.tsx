'use client';

import { useState } from 'react';
import { Check, Circle, XCircle, Archive, Trash2, MoreHorizontal } from 'lucide-react';
import { useStore } from '@/lib/store';
import type { ProjectStatus } from '@/lib/chat-types';
import { cn } from '@/lib/utils';
import { ScrollArea } from '@/components/ui/scroll-area';
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuItem,
  DropdownMenuTrigger,
} from '@/components/ui/dropdown-menu';

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
            <div
              key={project.id}
              role="button"
              tabIndex={0}
              onClick={() => setActiveProject(project.id)}
              onKeyDown={(e) => { if (e.key === 'Enter' || e.key === ' ') setActiveProject(project.id); }}
              onContextMenu={(e) => {
                e.preventDefault();
                if (window.confirm(`Delete project "${project.name}"?`)) {
                  deleteProject(project.id);
                }
              }}
              className={cn(
                'w-full flex items-center gap-2 px-2 py-2 rounded text-sm hover:bg-sidebar-accent text-left group transition-colors cursor-pointer',
                isActive && 'bg-sidebar-accent'
              )}
            >
              <StatusIcon status={project.status} />
              <span className="flex-1 truncate text-foreground">{project.name}</span>
              <span className="text-[10px] text-muted-foreground shrink-0 group-hover:hidden">
                {formatRelativeTime(project.updatedAt)}
              </span>
              <DropdownMenu>
                <DropdownMenuTrigger asChild>
                  <button
                    onClick={(e) => e.stopPropagation()}
                    className="hidden group-hover:flex shrink-0 p-0.5 rounded hover:bg-sidebar-accent/80 text-muted-foreground hover:text-foreground"
                  >
                    <MoreHorizontal className="h-3.5 w-3.5" />
                  </button>
                </DropdownMenuTrigger>
                <DropdownMenuContent side="right" align="start" className="w-36">
                  <DropdownMenuItem
                    onClick={(e) => {
                      e.stopPropagation();
                      deleteProject(project.id);
                    }}
                    className="text-destructive focus:text-destructive"
                  >
                    <Trash2 className="h-3.5 w-3.5 mr-2" />
                    Delete
                  </DropdownMenuItem>
                </DropdownMenuContent>
              </DropdownMenu>
            </div>
          );
        })}
      </div>
    </ScrollArea>
  );
}
