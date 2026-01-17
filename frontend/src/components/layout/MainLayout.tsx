'use client';

import { useState } from 'react';
import { cn } from '@/lib/utils';

interface MainLayoutProps {
  sidebar: React.ReactNode;
  main: React.ReactNode;
  viewer?: React.ReactNode;
  header?: React.ReactNode;
}

export function MainLayout({ sidebar, main, viewer, header }: MainLayoutProps) {
  const [viewerCollapsed, setViewerCollapsed] = useState(false);

  return (
    <div className="h-screen flex flex-col bg-background">
      {/* Header Bar */}
      {header && (
        <header className="h-14 border-b border-border flex items-center px-4 shrink-0">
          {header}
        </header>
      )}

      {/* Main Content Area */}
      <div className="flex-1 flex overflow-hidden">
        {/* Sidebar */}
        <aside className="w-60 border-r border-border flex flex-col shrink-0 bg-sidebar-background">
          {sidebar}
        </aside>

        {/* Main Content */}
        <main className="flex-1 overflow-auto">
          {main}
        </main>

        {/* Protein Viewer Panel */}
        {viewer && (
          <aside
            className={cn(
              'border-l border-border shrink-0 transition-all duration-200 bg-card',
              viewerCollapsed ? 'w-12' : 'w-[400px]'
            )}
          >
            <div className="h-full flex flex-col">
              <div className="h-10 border-b border-border flex items-center justify-between px-3">
                <span className={cn('text-sm font-medium', viewerCollapsed && 'hidden')}>
                  Structure Viewer
                </span>
                <button
                  onClick={() => setViewerCollapsed(!viewerCollapsed)}
                  className="p-1 hover:bg-accent rounded"
                >
                  {viewerCollapsed ? '◀' : '▶'}
                </button>
              </div>
              <div className={cn('flex-1 overflow-hidden', viewerCollapsed && 'hidden')}>
                {viewer}
              </div>
            </div>
          </aside>
        )}
      </div>
    </div>
  );
}

export default MainLayout;
