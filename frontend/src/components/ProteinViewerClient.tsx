'use client';

import { useEffect, useRef, useState } from 'react';

interface ProteinViewerClientProps {
  pdbContent: string | null;
  className?: string;
}

// Global state for Molstar plugin - persists across component remounts
let globalPlugin: any = null;
let globalInitPromise: Promise<any> | null = null;
let globalContainer: HTMLDivElement | null = null;

export function ProteinViewerClient({ pdbContent, className = '' }: ProteinViewerClientProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [isReady, setIsReady] = useState(false);

  // Initialize Molstar on mount
  useEffect(() => {
    const container = containerRef.current;
    if (!container) return;

    // If plugin already exists and initialized, just mark ready
    if (globalPlugin) {
      console.log('[ProteinViewer] Reusing existing plugin');
      // Re-attach plugin to new container if needed
      if (globalContainer !== container) {
        console.log('[ProteinViewer] Re-attaching to new container');
        // Move the canvas to the new container
        if (globalContainer && globalContainer.firstChild) {
          container.innerHTML = '';
          while (globalContainer.firstChild) {
            container.appendChild(globalContainer.firstChild);
          }
        }
        globalContainer = container;
      }
      setIsReady(true);
      return;
    }

    // If already initializing, wait for completion
    if (globalInitPromise) {
      console.log('[ProteinViewer] Waiting for existing initialization...');
      globalInitPromise.then(() => {
        if (globalPlugin) {
          globalContainer = container;
          setIsReady(true);
        }
      }).catch(() => {
        setError('Failed to initialize 3D viewer');
      });
      return;
    }

    console.log('[ProteinViewer] Starting Molstar initialization...');
    globalContainer = container;

    globalInitPromise = (async () => {
      try {
        // Dynamic imports
        const molstarUI = await import('molstar/lib/mol-plugin-ui');
        const molstarReact = await import('molstar/lib/mol-plugin-ui/react18');
        const molstarSpec = await import('molstar/lib/mol-plugin-ui/spec');

        // Use the stored global container
        if (!globalContainer) {
          console.log('[ProteinViewer] No container available');
          return null;
        }

        // Clear container
        globalContainer.innerHTML = '';

        console.log('[ProteinViewer] Creating plugin UI...');
        const plugin = await molstarUI.createPluginUI({
          target: globalContainer,
          render: molstarReact.renderReact18,
          spec: {
            ...molstarSpec.DefaultPluginUISpec(),
            layout: {
              initial: {
                isExpanded: false,
                showControls: false,
                controlsDisplay: 'reactive' as const,
              },
            },
          },
        });

        console.log('[ProteinViewer] Plugin created successfully');
        globalPlugin = plugin;
        return plugin;
      } catch (err) {
        console.error('[ProteinViewer] Failed to initialize Molstar:', err);
        globalInitPromise = null; // Allow retry
        throw err;
      }
    })();

    globalInitPromise.then((plugin) => {
      if (plugin) {
        setIsReady(true);
      }
    }).catch(() => {
      setError('Failed to initialize 3D viewer');
    });
  }, []);

  // Load structure when pdbContent changes
  useEffect(() => {
    if (!pdbContent || !globalPlugin || !isReady) return;

    const loadStructure = async () => {
      setLoading(true);
      setError(null);

      try {
        console.log('[ProteinViewer] Clearing existing structures...');
        await globalPlugin.clear();

        const isCif = pdbContent.trimStart().startsWith('data_');
        const format = isCif ? 'mmcif' : 'pdb';
        const extension = isCif ? 'cif' : 'pdb';

        console.log(`[ProteinViewer] Detected format: ${format}`);

        const data = await globalPlugin.builders.data.rawData({
          data: pdbContent,
          label: `structure.${extension}`,
        });

        const trajectory = await globalPlugin.builders.structure.parseTrajectory(data, format);
        await globalPlugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');

        globalPlugin.canvas3d?.requestCameraReset();
        console.log('[ProteinViewer] Structure loaded successfully');
      } catch (err) {
        console.error('[ProteinViewer] Failed to load structure:', err);
        setError(`Failed to load structure: ${err instanceof Error ? err.message : 'Unknown error'}`);
      } finally {
        setLoading(false);
      }
    };

    loadStructure();
  }, [pdbContent, isReady]);

  return (
    <div className={`relative bg-gray-100 rounded-lg overflow-hidden ${className}`}>
      <div ref={containerRef} className="w-full h-full" style={{ minHeight: '300px' }} />

      {!pdbContent && !error && !isReady && (
        <div className="absolute inset-0 flex items-center justify-center bg-gray-100 pointer-events-none">
          <div className="flex items-center gap-2 text-gray-600">
            <div className="w-5 h-5 border-2 border-gray-400 border-t-transparent rounded-full animate-spin" />
            Initializing viewer...
          </div>
        </div>
      )}

      {!pdbContent && !error && isReady && (
        <div className="absolute inset-0 flex items-center justify-center bg-gray-100 pointer-events-none">
          <p className="text-gray-500">No structure to display</p>
        </div>
      )}

      {loading && (
        <div className="absolute inset-0 flex items-center justify-center bg-white/80 pointer-events-none">
          <div className="flex items-center gap-2 text-gray-600">
            <div className="w-5 h-5 border-2 border-gray-400 border-t-transparent rounded-full animate-spin" />
            Loading structure...
          </div>
        </div>
      )}

      {error && (
        <div className="absolute inset-0 flex items-center justify-center bg-white/95">
          <div className="text-center p-4">
            <p className="text-red-600 mb-2">{error}</p>
            <p className="text-xs text-gray-500">Fallback: Showing raw PDB</p>
            <pre className="mt-2 text-xs text-gray-700 bg-gray-100 p-2 rounded max-h-48 overflow-auto text-left">
              {pdbContent?.slice(0, 2000)}...
            </pre>
          </div>
        </div>
      )}

      {pdbContent && !error && isReady && (
        <div className="absolute bottom-2 left-2 text-xs text-gray-600 bg-white/70 px-2 py-1 rounded pointer-events-none">
          Drag to rotate • Scroll to zoom • Shift+drag to pan
        </div>
      )}
    </div>
  );
}
