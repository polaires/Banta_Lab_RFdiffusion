'use client';

import { useEffect, useState } from 'react';
import { useStore } from '@/lib/store';
import api, { ApiMode } from '@/lib/api';

export function ConnectionModal() {
  const {
    backendUrl,
    setBackendUrl,
    health,
    setHealth,
    connectionModalOpen,
    setConnectionModalOpen,
  } = useStore();

  const [inputUrl, setInputUrl] = useState(backendUrl);
  const [checking, setChecking] = useState(false);
  const [showHelp, setShowHelp] = useState(false);
  const [connectionMode, setConnectionMode] = useState<ApiMode>(api.getMode());

  const checkConnection = async () => {
    setChecking(true);
    try {
      api.setMode(connectionMode);
      if (connectionMode === 'traditional') {
        api.setBaseUrl(inputUrl);
        setBackendUrl(inputUrl);
      }
      const healthResponse = await api.checkHealth();
      setHealth(healthResponse);
    } catch {
      setHealth(null);
    } finally {
      setChecking(false);
    }
  };

  const handleModeChange = (mode: ApiMode) => {
    setConnectionMode(mode);
    api.setMode(mode);
    setHealth(null); // Reset health when mode changes
  };

  // Initial connection check when modal opens
  useEffect(() => {
    if (connectionModalOpen && !health) {
      checkConnection();
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [connectionModalOpen]);

  // Sync inputUrl when backendUrl changes
  useEffect(() => {
    setInputUrl(backendUrl);
  }, [backendUrl]);

  // Close on Escape key
  useEffect(() => {
    const handleEscape = (e: KeyboardEvent) => {
      if (e.key === 'Escape') setConnectionModalOpen(false);
    };
    if (connectionModalOpen) {
      document.addEventListener('keydown', handleEscape);
      return () => document.removeEventListener('keydown', handleEscape);
    }
  }, [connectionModalOpen, setConnectionModalOpen]);

  if (!connectionModalOpen) return null;

  const isConnected = health?.status === 'healthy';
  const isMockMode = health?.mode === 'mock';

  // Handle both old format (models_loaded array) and new format (models object)
  const getModelsDisplay = () => {
    if (!health) return 'None';

    if (health.models && typeof health.models === 'object') {
      const availableModels = Object.entries(health.models)
        .filter(([, info]) => (info as { available?: boolean })?.available)
        .map(([name]) => name.toUpperCase());
      return availableModels.length > 0 ? availableModels.join(', ') : 'None (mock mode)';
    }

    if (health.models_loaded && Array.isArray(health.models_loaded)) {
      return health.models_loaded.join(', ') || 'None';
    }

    return 'None';
  };

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center p-4">
      {/* Backdrop */}
      <div
        className="absolute inset-0 bg-black/40 backdrop-blur-sm"
        onClick={() => setConnectionModalOpen(false)}
      />

      {/* Modal */}
      <div className="relative bg-white rounded-2xl shadow-2xl w-full max-w-md border border-slate-200">
        {/* Header */}
        <div className="flex items-center justify-between px-6 py-4 border-b border-slate-200">
          <h2 className="text-base font-bold text-slate-900 flex items-center gap-2">
            <span className={`material-symbols-outlined ${isConnected ? 'text-emerald-600' : 'text-red-500'}`}>
              {isConnected ? 'wifi' : 'wifi_off'}
            </span>
            Backend Connection
          </h2>
          <button
            onClick={() => setConnectionModalOpen(false)}
            className="p-1.5 rounded-lg text-slate-400 hover:text-slate-600 hover:bg-slate-100 transition-colors"
          >
            <span className="material-symbols-outlined text-xl">close</span>
          </button>
        </div>

        {/* Content */}
        <div className="px-6 py-5 space-y-5">
          {/* Mode Selection */}
          <div className="space-y-3">
            <label className="text-xs font-bold text-slate-700 uppercase tracking-wider">Connection Mode</label>
            <div className="grid grid-cols-3 gap-2">
              <button
                onClick={() => {
                  handleModeChange('traditional');
                  setInputUrl('http://localhost:8000');
                }}
                className={`p-3 rounded-xl border-2 transition-all text-left ${
                  connectionMode === 'traditional' && inputUrl === 'http://localhost:8000'
                    ? 'border-amber-500 bg-amber-50 shadow-sm'
                    : 'border-slate-200 hover:border-slate-300 bg-white'
                }`}
              >
                <div className="flex items-center gap-1.5 mb-1">
                  <span className="material-symbols-outlined text-amber-600 text-lg">developer_mode</span>
                  <span className="font-semibold text-slate-900 text-xs">Local Dev</span>
                </div>
                <p className="text-[10px] text-slate-500">Docker on localhost</p>
              </button>
              <button
                onClick={() => handleModeChange('traditional')}
                className={`p-3 rounded-xl border-2 transition-all text-left ${
                  connectionMode === 'traditional' && inputUrl !== 'http://localhost:8000'
                    ? 'border-blue-500 bg-blue-50 shadow-sm'
                    : 'border-slate-200 hover:border-slate-300 bg-white'
                }`}
              >
                <div className="flex items-center gap-1.5 mb-1">
                  <span className="material-symbols-outlined text-slate-600 text-lg">dns</span>
                  <span className="font-semibold text-slate-900 text-xs">Traditional</span>
                </div>
                <p className="text-[10px] text-slate-500">Always-on GPU Pod</p>
              </button>
              <button
                onClick={() => handleModeChange('serverless')}
                className={`p-3 rounded-xl border-2 transition-all text-left ${
                  connectionMode === 'serverless'
                    ? 'border-blue-500 bg-blue-50 shadow-sm'
                    : 'border-slate-200 hover:border-slate-300 bg-white'
                }`}
              >
                <div className="flex items-center gap-1.5 mb-1">
                  <span className="material-symbols-outlined text-slate-600 text-lg">cloud</span>
                  <span className="font-semibold text-slate-900 text-xs">Serverless</span>
                </div>
                <p className="text-[10px] text-slate-500">Pay-per-use</p>
              </button>
            </div>
          </div>

          {/* URL Input (Traditional Mode Only) */}
          {connectionMode === 'traditional' && (
            <div className="space-y-3">
              <label className="text-xs font-bold text-slate-700 uppercase tracking-wider">RunPod API URL</label>
              <div className="flex gap-2">
                <input
                  type="text"
                  value={inputUrl}
                  onChange={(e) => setInputUrl(e.target.value)}
                  placeholder="https://your-runpod-endpoint:8000"
                  className="flex-1 px-4 py-3 bg-slate-50 rounded-xl border border-slate-200 focus:border-blue-500 focus:ring-2 focus:ring-blue-500/20 focus:outline-none text-sm text-slate-900 transition-all"
                  onKeyDown={(e) => e.key === 'Enter' && checkConnection()}
                />
                <button
                  onClick={checkConnection}
                  disabled={checking}
                  className="px-5 py-3 bg-blue-600 hover:bg-blue-700 disabled:bg-slate-300 text-white rounded-xl text-sm font-semibold transition-all flex items-center gap-2 shadow-sm hover:shadow-md"
                >
                  {checking ? (
                    <>
                      <span className="material-symbols-outlined animate-spin text-lg">progress_activity</span>
                      <span className="hidden sm:inline">Checking</span>
                    </>
                  ) : (
                    'Connect'
                  )}
                </button>
              </div>
            </div>
          )}

          {/* Serverless Info */}
          {connectionMode === 'serverless' && (
            <div className="space-y-4">
              <div className="p-4 bg-emerald-50 border border-emerald-200 rounded-xl">
                <div className="flex items-center gap-2 mb-1">
                  <span className="material-symbols-outlined text-emerald-600">bolt</span>
                  <p className="text-sm font-semibold text-emerald-700">Serverless Mode</p>
                </div>
                <p className="text-xs text-emerald-600">
                  Using RunPod Serverless via Vercel Edge Function. Cold start: ~60-90 seconds.
                </p>
              </div>
              <button
                onClick={checkConnection}
                disabled={checking}
                className="w-full py-3 bg-blue-600 hover:bg-blue-700 disabled:bg-slate-300 text-white rounded-xl text-sm font-semibold transition-all flex items-center justify-center gap-2 shadow-sm hover:shadow-md"
              >
                {checking ? (
                  <>
                    <span className="material-symbols-outlined animate-spin text-lg">progress_activity</span>
                    Checking Endpoint...
                  </>
                ) : (
                  'Test Connection'
                )}
              </button>
            </div>
          )}

          {/* Status Display */}
          {health && (
            <div className="bg-slate-50 rounded-xl p-4 space-y-3 border border-slate-200">
              {/* Connection status */}
              <div className="flex items-center justify-between">
                <span className="text-sm text-slate-600">Status</span>
                <span className={`text-sm font-semibold flex items-center gap-1.5 ${isConnected ? 'text-emerald-600' : 'text-red-600'}`}>
                  <span className={`w-2 h-2 rounded-full ${isConnected ? 'bg-emerald-500' : 'bg-red-500'}`} />
                  {isConnected ? 'Connected' : 'Disconnected'}
                </span>
              </div>

              {/* Mode */}
              {health.mode && (
                <div className="flex items-center justify-between">
                  <span className="text-sm text-slate-600">Mode</span>
                  <span className={`text-sm font-semibold flex items-center gap-1.5 ${isMockMode ? 'text-amber-600' : 'text-emerald-600'}`}>
                    <span className="material-symbols-outlined text-base">
                      {isMockMode ? 'warning' : 'auto_awesome'}
                    </span>
                    {isMockMode ? 'Mock (Demo)' : 'Real (Foundry)'}
                  </span>
                </div>
              )}

              {/* GPU */}
              <div className="flex items-center justify-between">
                <span className="text-sm text-slate-600">GPU</span>
                {health.gpu_available ? (
                  <span className="text-sm font-semibold text-emerald-600 flex items-center gap-1.5">
                    <span className="material-symbols-outlined text-base">memory</span>
                    {health.gpu_name || 'Available'}
                    {health.gpu_memory_gb && ` (${Math.round(health.gpu_memory_gb)}GB)`}
                  </span>
                ) : (
                  <span className="text-sm text-red-600 font-medium">Not Available</span>
                )}
              </div>

              {/* Models */}
              <div className="flex items-center justify-between">
                <span className="text-sm text-slate-600">Models</span>
                <span className="text-sm font-semibold text-slate-900">{getModelsDisplay()}</span>
              </div>
            </div>
          )}

          {/* Not connected warning */}
          {!health && !checking && (
            <div className="bg-red-50 border border-red-200 rounded-xl p-4 flex items-start gap-3">
              <span className="material-symbols-outlined text-red-500">error</span>
              <p className="text-sm text-red-700">
                Unable to connect to backend. Make sure your RunPod pod is running.
              </p>
            </div>
          )}

          {/* Setup Help Accordion */}
          <div className="border-t border-slate-200 pt-4">
            <button
              onClick={() => setShowHelp(!showHelp)}
              className="flex items-center gap-2 text-sm text-slate-600 hover:text-slate-800 w-full transition-colors"
            >
              <span className="material-symbols-outlined text-lg">help</span>
              <span className="font-medium">Setup Instructions</span>
              <span className="material-symbols-outlined text-lg ml-auto">
                {showHelp ? 'expand_less' : 'expand_more'}
              </span>
            </button>

            {showHelp && (
              <div className="mt-4 space-y-3 text-sm">
                <div className="bg-slate-50 rounded-xl p-4 space-y-2 border border-slate-200">
                  <h3 className="font-semibold text-slate-900">Quick Setup (RunPod)</h3>
                  <ol className="list-decimal list-inside space-y-1.5 text-slate-600 text-xs">
                    <li>
                      Deploy a GPU pod (A40 recommended) on{' '}
                      <a
                        href="https://console.runpod.io/deploy"
                        target="_blank"
                        rel="noopener noreferrer"
                        className="text-blue-600 hover:underline font-medium"
                      >
                        RunPod
                        <span className="material-symbols-outlined text-xs align-middle ml-0.5">open_in_new</span>
                      </a>
                    </li>
                    <li>
                      Add <code className="bg-slate-200 px-1.5 py-0.5 rounded text-slate-800">8000</code> to HTTP Ports
                    </li>
                    <li>
                      Open Jupyter Lab and run the setup cell from{' '}
                      <a
                        href="https://github.com/polaires/Banta_Lab_RFdiffusion/blob/main/RUNPOD_QUICK_SETUP.md"
                        target="_blank"
                        rel="noopener noreferrer"
                        className="text-blue-600 hover:underline font-medium"
                      >
                        RUNPOD_QUICK_SETUP.md
                        <span className="material-symbols-outlined text-xs align-middle ml-0.5">open_in_new</span>
                      </a>
                    </li>
                    <li>Copy your API URL and paste above</li>
                  </ol>
                </div>

                <div className="bg-slate-50 rounded-xl p-4 space-y-2 border border-slate-200">
                  <h3 className="font-semibold text-slate-900">URL Format</h3>
                  <code className="block bg-slate-200 px-3 py-2 rounded-lg text-xs text-slate-800 font-mono">
                    https://abc123xyz-8000.proxy.runpod.net
                  </code>
                </div>

                <div className="bg-slate-50 rounded-xl p-4 space-y-2 border border-slate-200">
                  <h3 className="font-semibold text-slate-900">Troubleshooting</h3>
                  <ul className="space-y-1.5 text-slate-600 text-xs">
                    <li>
                      <span className="text-amber-600 font-medium">Connection refused:</span> Re-run setup cell
                    </li>
                    <li>
                      <span className="text-amber-600 font-medium">CORS error:</span> Check port 8000 is exposed
                    </li>
                    <li>
                      <span className="text-amber-600 font-medium">Mock mode:</span> Checkpoints still downloading
                    </li>
                  </ul>
                </div>

                <p className="text-xs text-slate-500">
                  Your URL is automatically saved and persists across sessions.
                </p>
              </div>
            )}
          </div>
        </div>

        {/* Footer */}
        <div className="px-6 py-3 bg-slate-50 border-t border-slate-200 rounded-b-2xl">
          <p className="text-xs text-slate-500 text-center">
            Connection status auto-refreshes every 30 seconds
          </p>
        </div>
      </div>
    </div>
  );
}
