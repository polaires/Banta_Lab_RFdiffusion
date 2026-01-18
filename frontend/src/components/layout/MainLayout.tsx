'use client';

import { useState, useCallback, useRef, useEffect, createContext, useContext } from 'react';
import { cn } from '@/lib/utils';
import { GripVertical, PanelRightClose, PanelRight } from 'lucide-react';

// Context to notify children of viewer resize
export const ViewerResizeContext = createContext<{ width: number; collapsed: boolean }>({ width: 400, collapsed: false });
export const useViewerResize = () => useContext(ViewerResizeContext);

// Layout constants
const MIN_VIEWER_WIDTH = 300;

interface MainLayoutProps {
  sidebar: React.ReactNode;
  main: React.ReactNode;
  viewer?: React.ReactNode;
  header?: React.ReactNode;
}

export function MainLayout({ sidebar, main, viewer, header }: MainLayoutProps) {
  const [viewerCollapsed, setViewerCollapsed] = useState(false);
  const [viewerWidth, setViewerWidth] = useState(400);
  const [isResizing, setIsResizing] = useState(false);
  const resizeRef = useRef<HTMLDivElement>(null);

  const handleMouseDown = useCallback((e: React.MouseEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setIsResizing(true);
  }, []);

  // Handle resize drag - no max limit
  useEffect(() => {
    const handleMouseMove = (e: MouseEvent) => {
      if (!isResizing) return;

      // Calculate new width based on mouse position from right edge
      const newWidth = window.innerWidth - e.clientX;
      // Only enforce minimum width, no maximum
      setViewerWidth(Math.max(MIN_VIEWER_WIDTH, newWidth));
    };

    const handleMouseUp = () => {
      setIsResizing(false);
    };

    if (isResizing) {
      document.addEventListener('mousemove', handleMouseMove);
      document.addEventListener('mouseup', handleMouseUp);
      document.body.style.cursor = 'col-resize';
      document.body.style.userSelect = 'none';
    }

    return () => {
      document.removeEventListener('mousemove', handleMouseMove);
      document.removeEventListener('mouseup', handleMouseUp);
      document.body.style.cursor = '';
      document.body.style.userSelect = '';
    };
  }, [isResizing]);

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
        {/* Sidebar - handles overflow gracefully */}
        <aside className="w-60 border-r border-border flex flex-col shrink-0 overflow-hidden bg-sidebar-background">
          <div className="flex-1 overflow-y-auto overflow-x-hidden">
            {sidebar}
          </div>
        </aside>

        {/* Main Content - can shrink but handles overflow */}
        <main className="flex-1 min-w-0 overflow-hidden">
          {main}
        </main>

        {/* Protein Viewer Panel */}
        {viewer && (
          <aside
            className={cn(
              'border-l border-border shrink-0 bg-card relative',
              viewerCollapsed ? 'w-12' : '',
              isResizing ? 'transition-none' : 'transition-[width] duration-150'
            )}
            style={{ width: viewerCollapsed ? 48 : viewerWidth }}
          >
            {/* Resize Handle - Claude artifact style */}
            {!viewerCollapsed && (
              <div
                ref={resizeRef}
                onMouseDown={handleMouseDown}
                className={cn(
                  'absolute left-0 top-0 bottom-0 w-2 cursor-col-resize z-20 group',
                  'hover:bg-primary/10',
                  isResizing && 'bg-primary/20'
                )}
              >
                {/* Visible drag indicator line */}
                <div className={cn(
                  'absolute left-0 top-0 bottom-0 w-[3px]',
                  'bg-transparent group-hover:bg-primary/40 transition-colors',
                  isResizing && 'bg-primary/60'
                )} />
                {/* Grip icon on hover */}
                <div className={cn(
                  'absolute left-1/2 top-1/2 -translate-x-1/2 -translate-y-1/2',
                  'opacity-0 group-hover:opacity-100 transition-opacity pointer-events-none',
                  isResizing && 'opacity-100'
                )}>
                  <GripVertical className="h-8 w-8 text-muted-foreground drop-shadow-sm" />
                </div>
              </div>
            )}

            <div className="h-full flex flex-col">
              <div className="h-10 border-b border-border flex items-center justify-between px-3">
                <span className={cn('text-sm font-medium', viewerCollapsed && 'hidden')}>
                  Structure Viewer
                </span>
                <button
                  onClick={() => setViewerCollapsed(!viewerCollapsed)}
                  className="p-1.5 hover:bg-accent rounded-md transition-colors"
                  title={viewerCollapsed ? 'Expand panel' : 'Collapse panel'}
                >
                  {viewerCollapsed ? (
                    <PanelRight className="h-4 w-4 text-muted-foreground" />
                  ) : (
                    <PanelRightClose className="h-4 w-4 text-muted-foreground" />
                  )}
                </button>
              </div>
              <div className={cn('flex-1 overflow-hidden', viewerCollapsed && 'hidden')}>
                <ViewerResizeContext.Provider value={{ width: viewerWidth, collapsed: viewerCollapsed }}>
                  {viewer}
                </ViewerResizeContext.Provider>
              </div>
            </div>
          </aside>
        )}
      </div>
    </div>
  );
}

export default MainLayout;
