'use client';

import { useEffect, useRef, useState } from 'react';

interface ProteinViewerClientProps {
  pdbContent: string | null;
  className?: string;
}

// Track initialization globally to handle StrictMode
let globalPlugin: any = null;
let globalInitPromise: Promise<any> | null = null;

export function ProteinViewerClient({ pdbContent, className = '' }: ProteinViewerClientProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [isReady, setIsReady] = useState(false);

  // Initialize on mount
  useEffect(() => {
    let mounted = true;

    const initPlugin = async () => {
      if (!containerRef.current) return;

      // If already initialized, reuse the plugin
      if (globalPlugin) {
        console.log('[ProteinViewer] Reusing existing plugin');
        setIsReady(true);
        return;
      }

      // If already initializing, wait for that to complete
      if (globalInitPromise) {
        console.log('[ProteinViewer] Waiting for existing initialization...');
        try {
          await globalInitPromise;
          if (mounted && globalPlugin) {
            setIsReady(true);
          }
        } catch (err) {
          if (mounted) {
            setError('Failed to initialize 3D viewer');
          }
        }
        return;
      }

      console.log('[ProteinViewer] Starting Molstar initialization...');

      globalInitPromise = (async () => {
        try {
          // Dynamic imports to avoid bundler issues
          const molstarUI = await import('molstar/lib/mol-plugin-ui');
          const molstarReact = await import('molstar/lib/mol-plugin-ui/react18');
          const molstarSpec = await import('molstar/lib/mol-plugin-ui/spec');

          if (!containerRef.current || !mounted) {
            console.log('[ProteinViewer] Aborted - unmounted');
            return null;
          }

          // Clear any existing content
          containerRef.current.innerHTML = '';

          console.log('[ProteinViewer] Creating plugin UI...');
          const plugin = await molstarUI.createPluginUI({
            target: containerRef.current,
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
          throw err;
        }
      })();

      try {
        const plugin = await globalInitPromise;
        if (plugin && mounted) {
          setIsReady(true);
        }
      } catch (err) {
        if (mounted) {
          setError('Failed to initialize 3D viewer');
        }
      }
    };

    initPlugin();

    return () => {
      mounted = false;
    };
  }, []);

  // Load structure content when it changes
  useEffect(() => {
    if (!pdbContent || !globalPlugin || !isReady) return;

    const loadStructure = async () => {
      setLoading(true);
      setError(null);

      try {
        console.log('[ProteinViewer] Clearing existing structures...');
        await globalPlugin.clear();

        // Detect format: CIF files start with "data_", PDB files start with HEADER/ATOM/etc
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

      {/* No structure placeholder */}
      {!pdbContent && !error && (
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

      {/* Viewer controls hint */}
      {pdbContent && !error && (
        <div className="absolute bottom-2 left-2 text-xs text-gray-600 bg-white/70 px-2 py-1 rounded pointer-events-none">
          Drag to rotate • Scroll to zoom • Shift+drag to pan
        </div>
      )}
    </div>
  );
}
