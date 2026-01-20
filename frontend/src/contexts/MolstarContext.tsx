'use client';

import {
  createContext,
  useContext,
  useRef,
  useState,
  useCallback,
  useEffect,
  type ReactNode,
} from 'react';
import type { PluginUIContext } from 'molstar/lib/mol-plugin-ui/context';
import type { StateObjectRef } from 'molstar/lib/mol-state';

interface MolstarState {
  plugin: PluginUIContext | null;
  structureRef: StateObjectRef | null;
  isReady: boolean;
  isLoading: boolean;
  error: string | null;
}

interface MolstarContextValue extends MolstarState {
  initPlugin: (container: HTMLDivElement) => Promise<void>;
  loadStructure: (pdbContent: string) => Promise<StateObjectRef | null>;
  clearStructure: () => Promise<void>;
  dispose: () => void;
}

const MolstarContext = createContext<MolstarContextValue | null>(null);

export function useMolstar() {
  const context = useContext(MolstarContext);
  if (!context) {
    throw new Error('useMolstar must be used within MolstarProvider');
  }
  return context;
}

interface MolstarProviderProps {
  children: ReactNode;
}

export function MolstarProvider({ children }: MolstarProviderProps) {
  const pluginRef = useRef<PluginUIContext | null>(null);
  const structureRefRef = useRef<StateObjectRef | null>(null);
  const initPromiseRef = useRef<Promise<void> | null>(null);

  const [state, setState] = useState<MolstarState>({
    plugin: null,
    structureRef: null,
    isReady: false,
    isLoading: false,
    error: null,
  });

  const initPlugin = useCallback(async (container: HTMLDivElement) => {
    // Already initialized
    if (pluginRef.current) {
      setState(s => ({ ...s, isReady: true }));
      return;
    }

    // Already initializing
    if (initPromiseRef.current) {
      await initPromiseRef.current;
      return;
    }

    setState(s => ({ ...s, isLoading: true, error: null }));

    initPromiseRef.current = (async () => {
      try {
        const [molstarUI, molstarReact, molstarSpec] = await Promise.all([
          import('molstar/lib/mol-plugin-ui'),
          import('molstar/lib/mol-plugin-ui/react18'),
          import('molstar/lib/mol-plugin-ui/spec'),
        ]);

        container.innerHTML = '';

        const plugin = await molstarUI.createPluginUI({
          target: container,
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

        pluginRef.current = plugin;
        setState(s => ({
          ...s,
          plugin,
          isReady: true,
          isLoading: false,
        }));
      } catch (err) {
        const message = err instanceof Error ? err.message : 'Failed to initialize viewer';
        setState(s => ({ ...s, error: message, isLoading: false }));
        throw err;
      }
    })();

    await initPromiseRef.current;
  }, []);

  const loadStructure = useCallback(async (pdbContent: string): Promise<StateObjectRef | null> => {
    const plugin = pluginRef.current;
    if (!plugin) return null;

    setState(s => ({ ...s, isLoading: true, error: null }));

    try {
      // Ensure plugin canvas is fully initialized (critical for serverless/production)
      if (!plugin.canvas3d) {
        console.warn('[MolstarContext] Waiting for canvas3d initialization...');
        await new Promise(resolve => setTimeout(resolve, 100));
        if (!plugin.canvas3d) {
          throw new Error('Molstar canvas3d not initialized - WebGL may not be available');
        }
      }

      await plugin.clear();

      const isCif = pdbContent.trimStart().startsWith('data_');
      const format = isCif ? 'mmcif' : 'pdb';

      const data = await plugin.builders.data.rawData(
        { data: pdbContent, label: 'structure' },
        { state: { isGhost: true } }
      );

      const trajectory = await plugin.builders.structure.parseTrajectory(data, format);
      if (!trajectory) {
        throw new Error('Failed to parse trajectory - PDB data may be malformed');
      }
      const model = await plugin.builders.structure.createModel(trajectory);
      if (!model) {
        throw new Error('Failed to create model from trajectory');
      }
      const structure = await plugin.builders.structure.createStructure(model);
      if (!structure || !structure.ref) {
        throw new Error('Failed to create structure from model - this may be a WebGL or Molstar initialization issue');
      }

      structureRefRef.current = structure.ref;
      setState(s => ({
        ...s,
        structureRef: structure.ref,
        isLoading: false,
      }));

      return structure.ref;
    } catch (err) {
      const message = err instanceof Error ? err.message : 'Failed to load structure';
      setState(s => ({ ...s, error: message, isLoading: false }));
      return null;
    }
  }, []);

  const clearStructure = useCallback(async () => {
    const plugin = pluginRef.current;
    if (!plugin) return;

    await plugin.clear();
    structureRefRef.current = null;
    setState(s => ({ ...s, structureRef: null, error: null }));
  }, []);

  const dispose = useCallback(() => {
    pluginRef.current?.dispose();
    pluginRef.current = null;
    structureRefRef.current = null;
    initPromiseRef.current = null;
    setState({
      plugin: null,
      structureRef: null,
      isReady: false,
      isLoading: false,
      error: null,
    });
  }, []);

  // Cleanup on unmount
  useEffect(() => {
    return () => {
      dispose();
    };
  }, [dispose]);

  const value: MolstarContextValue = {
    ...state,
    initPlugin,
    loadStructure,
    clearStructure,
    dispose,
  };

  return (
    <MolstarContext.Provider value={value}>
      {children}
    </MolstarContext.Provider>
  );
}
