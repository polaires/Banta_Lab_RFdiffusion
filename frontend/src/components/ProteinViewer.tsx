'use client';

import { useEffect, useRef, useState, useCallback } from 'react';
// Import Molstar CSS
import 'molstar/lib/mol-plugin-ui/skin/dark.scss';

interface ProteinViewerProps {
  pdbContent: string | null;
  className?: string;
}

// Track initialization globally to handle StrictMode
let globalInitPromise: Promise<any> | null = null;

export function ProteinViewer({ pdbContent, className = '' }: ProteinViewerProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const pluginRef = useRef<any>(null);
  const mountedRef = useRef(true);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [isReady, setIsReady] = useState(false);

  // Initialize on mount
  useEffect(() => {
    mountedRef.current = true;

    const initPlugin = async () => {
      if (!containerRef.current) return;

      // If already initializing, wait for that to complete
      if (globalInitPromise) {
        console.log('[ProteinViewer] Waiting for existing initialization...');
        await globalInitPromise;
        if (pluginRef.current && mountedRef.current) {
          setIsReady(true);
        }
        return;
      }

      console.log('[ProteinViewer] Starting Molstar initialization...');

      globalInitPromise = (async () => {
        try {
          const { createPluginUI } = await import('molstar/lib/mol-plugin-ui');
          const { renderReact18 } = await import('molstar/lib/mol-plugin-ui/react18');
          const { DefaultPluginUISpec } = await import('molstar/lib/mol-plugin-ui/spec');

          if (!containerRef.current || !mountedRef.current) {
            console.log('[ProteinViewer] Aborted - unmounted');
            return null;
          }

          // Clear any existing content
          containerRef.current.innerHTML = '';

          console.log('[ProteinViewer] Creating plugin UI...');
          const plugin = await createPluginUI({
            target: containerRef.current,
            render: renderReact18,
            spec: {
              ...DefaultPluginUISpec(),
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
          return plugin;
        } catch (err) {
          console.error('[ProteinViewer] Failed to initialize Molstar:', err);
          throw err;
        }
      })();

      try {
        const plugin = await globalInitPromise;
        if (plugin && mountedRef.current) {
          pluginRef.current = plugin;
          (window as any).__molstarPlugin = plugin;
          setIsReady(true);
        }
      } catch (err) {
        if (mountedRef.current) {
          setError('Failed to initialize 3D viewer');
        }
      }
    };

    initPlugin();

    return () => {
      console.log('[ProteinViewer] Cleanup');
      mountedRef.current = false;
      // Don't dispose the plugin on StrictMode cleanup - let it persist
    };
  }, []);

  // Load PDB content when it changes
  useEffect(() => {
    const plugin = pluginRef.current;
    console.log('[ProteinViewer] Load PDB effect - pdbContent:', pdbContent ? `${pdbContent.length} chars` : null, 'isReady:', isReady);

    if (!pdbContent || !plugin || !isReady) return;

    // Check if WebGL context is lost
    if (plugin.canvas3d?.webgl?.isContextLost) {
      console.log('[ProteinViewer] WebGL context lost, cannot load PDB');
      setError('3D viewer context lost. Please refresh the page.');
      return;
    }

    const loadPdb = async () => {
      setLoading(true);
      setError(null);

      try {
        console.log('[ProteinViewer] Clearing existing structures...');
        await plugin.clear();

        console.log('[ProteinViewer] Loading PDB data...');
        const data = await plugin.builders.data.rawData({
          data: pdbContent,
          label: 'structure.pdb',
        });

        console.log('[ProteinViewer] Parsing trajectory...');
        const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');

        console.log('[ProteinViewer] Applying preset...');
        await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');

        console.log('[ProteinViewer] Resetting camera...');
        plugin.canvas3d?.requestCameraReset();
        console.log('[ProteinViewer] PDB loaded successfully');
      } catch (err) {
        console.error('[ProteinViewer] Failed to load PDB:', err);
        setError('Failed to load structure');
      } finally {
        setLoading(false);
      }
    };

    loadPdb();
  }, [pdbContent, isReady]);

  return (
    <div className={`relative bg-gray-900 rounded-lg overflow-hidden ${className}`}>
      <div ref={containerRef} className="w-full h-full" style={{ minHeight: '300px' }} />

      {/* No structure placeholder */}
      {!pdbContent && !error && (
        <div className="absolute inset-0 flex items-center justify-center bg-gray-800/90 pointer-events-none">
          <p className="text-gray-400">No structure to display</p>
        </div>
      )}

      {loading && (
        <div className="absolute inset-0 flex items-center justify-center bg-gray-900/80 pointer-events-none">
          <div className="flex items-center gap-2 text-gray-300">
            <div className="w-5 h-5 border-2 border-gray-300 border-t-transparent rounded-full animate-spin" />
            Loading structure...
          </div>
        </div>
      )}

      {error && (
        <div className="absolute inset-0 flex items-center justify-center bg-gray-900/90">
          <div className="text-center p-4">
            <p className="text-red-400 mb-2">{error}</p>
            <p className="text-xs text-gray-500">Fallback: Showing raw PDB</p>
            <pre className="mt-2 text-xs text-green-400 max-h-48 overflow-auto text-left">
              {pdbContent?.slice(0, 2000)}...
            </pre>
          </div>
        </div>
      )}

      {/* Viewer controls hint */}
      {pdbContent && !error && (
        <div className="absolute bottom-2 left-2 text-xs text-gray-500 pointer-events-none">
          Drag to rotate • Scroll to zoom • Shift+drag to pan
        </div>
      )}
    </div>
  );
}
