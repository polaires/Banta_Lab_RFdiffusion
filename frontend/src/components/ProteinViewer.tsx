'use client';

import { useEffect, useRef, useState } from 'react';

interface ProteinViewerProps {
  pdbContent: string | null;
  className?: string;
}

export function ProteinViewer({ pdbContent, className = '' }: ProteinViewerProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<{ plugin: unknown; dispose: () => void } | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [initialized, setInitialized] = useState(false);

  // Initialize Molstar viewer
  useEffect(() => {
    if (!containerRef.current || initialized) return;

    const initViewer = async () => {
      try {
        // Dynamically import Molstar to avoid SSR issues
        const { createPluginUI } = await import('molstar/lib/mol-plugin-ui');
        const { renderReact18 } = await import('molstar/lib/mol-plugin-ui/react18');
        const { DefaultPluginUISpec } = await import('molstar/lib/mol-plugin-ui/spec');

        const spec = DefaultPluginUISpec();
        spec.layout = {
          initial: {
            isExpanded: false,
            showControls: false,
            controlsDisplay: 'reactive',
          },
        };

        const plugin = await createPluginUI({
          target: containerRef.current!,
          spec,
          render: renderReact18,
        });

        viewerRef.current = {
          plugin,
          dispose: () => plugin.dispose(),
        };
        setInitialized(true);
      } catch (err) {
        console.error('Failed to initialize Molstar:', err);
        setError('Failed to initialize 3D viewer');
      }
    };

    initViewer();

    return () => {
      viewerRef.current?.dispose();
      viewerRef.current = null;
    };
  }, [initialized]);

  // Load PDB content when it changes
  useEffect(() => {
    if (!pdbContent || !viewerRef.current || !initialized) return;

    const loadPdb = async () => {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      const plugin = viewerRef.current?.plugin as any;
      if (!plugin) return;

      setLoading(true);
      setError(null);

      try {
        // Clear existing structures
        await plugin.clear();

        // Load PDB from string data
        const data = await plugin.builders.data.rawData({
          data: pdbContent,
          label: 'structure.pdb',
        });

        const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
        await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');

        // Reset camera to fit structure
        plugin.canvas3d?.requestCameraReset();
      } catch (err) {
        console.error('Failed to load PDB:', err);
        setError('Failed to load structure');
      } finally {
        setLoading(false);
      }
    };

    loadPdb();
  }, [pdbContent, initialized]);

  if (!pdbContent) {
    return (
      <div className={`flex items-center justify-center bg-gray-800 rounded-lg ${className}`}>
        <p className="text-gray-400">No structure to display</p>
      </div>
    );
  }

  return (
    <div className={`relative bg-gray-900 rounded-lg overflow-hidden ${className}`}>
      <div ref={containerRef} className="w-full h-full" style={{ minHeight: '300px' }} />

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
              {pdbContent.slice(0, 2000)}...
            </pre>
          </div>
        </div>
      )}

      {/* Viewer controls hint */}
      <div className="absolute bottom-2 left-2 text-xs text-gray-500 pointer-events-none">
        Drag to rotate • Scroll to zoom • Shift+drag to pan
      </div>
    </div>
  );
}
