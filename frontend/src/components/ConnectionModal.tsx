'use client';

import { useEffect, useState } from 'react';
import { useStore } from '@/lib/store';
import api from '@/lib/api';
import {
  X,
  Wifi,
  WifiOff,
  Cpu,
  AlertCircle,
  Sparkles,
  HelpCircle,
  ChevronDown,
  ChevronUp,
  ExternalLink,
  Loader2,
} from 'lucide-react';

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

  const checkConnection = async () => {
    setChecking(true);
    try {
      api.setBaseUrl(inputUrl);
      const healthResponse = await api.checkHealth();
      setHealth(healthResponse);
      setBackendUrl(inputUrl);
    } catch {
      setHealth(null);
    } finally {
      setChecking(false);
    }
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
      <div className="relative bg-white rounded-xl shadow-2xl w-full max-w-md border border-gray-200">
        {/* Header */}
        <div className="flex items-center justify-between px-5 py-4 border-b border-gray-200">
          <h2 className="text-lg font-semibold text-gray-900 flex items-center gap-2">
            {isConnected ? (
              <Wifi className="w-5 h-5 text-green-600" />
            ) : (
              <WifiOff className="w-5 h-5 text-red-500" />
            )}
            Backend Connection
          </h2>
          <button
            onClick={() => setConnectionModalOpen(false)}
            className="p-1 rounded-lg text-gray-400 hover:text-gray-600 hover:bg-gray-100 transition"
          >
            <X className="w-5 h-5" />
          </button>
        </div>

        {/* Content */}
        <div className="px-5 py-4 space-y-4">
          {/* URL Input */}
          <div className="space-y-2">
            <label className="text-sm text-gray-600">RunPod API URL</label>
            <div className="flex gap-2">
              <input
                type="text"
                value={inputUrl}
                onChange={(e) => setInputUrl(e.target.value)}
                placeholder="https://your-runpod-endpoint:8000"
                className="flex-1 px-3 py-2.5 bg-gray-50 rounded-lg border border-gray-300 focus:border-blue-500 focus:outline-none text-sm text-gray-900"
                onKeyDown={(e) => e.key === 'Enter' && checkConnection()}
              />
              <button
                onClick={checkConnection}
                disabled={checking}
                className="px-4 py-2.5 bg-blue-600 hover:bg-blue-700 disabled:bg-gray-400 text-white rounded-lg text-sm font-medium transition flex items-center gap-2"
              >
                {checking ? (
                  <>
                    <Loader2 className="w-4 h-4 animate-spin" />
                    <span className="hidden sm:inline">Checking</span>
                  </>
                ) : (
                  'Connect'
                )}
              </button>
            </div>
          </div>

          {/* Status Display */}
          {health && (
            <div className="bg-gray-50 rounded-lg p-4 space-y-3 border border-gray-200">
              {/* Connection status */}
              <div className="flex items-center justify-between">
                <span className="text-sm text-gray-600">Status</span>
                <span className={`text-sm font-medium flex items-center gap-1.5 ${isConnected ? 'text-green-600' : 'text-red-600'}`}>
                  <span className={`w-2 h-2 rounded-full ${isConnected ? 'bg-green-500' : 'bg-red-500'}`} />
                  {isConnected ? 'Connected' : 'Disconnected'}
                </span>
              </div>

              {/* Mode */}
              {health.mode && (
                <div className="flex items-center justify-between">
                  <span className="text-sm text-gray-600">Mode</span>
                  <span className={`text-sm font-medium flex items-center gap-1.5 ${isMockMode ? 'text-amber-600' : 'text-green-600'}`}>
                    {isMockMode ? (
                      <AlertCircle className="w-4 h-4" />
                    ) : (
                      <Sparkles className="w-4 h-4" />
                    )}
                    {isMockMode ? 'Mock (Demo)' : 'Real (Foundry)'}
                  </span>
                </div>
              )}

              {/* GPU */}
              <div className="flex items-center justify-between">
                <span className="text-sm text-gray-600">GPU</span>
                {health.gpu_available ? (
                  <span className="text-sm font-medium text-green-600 flex items-center gap-1.5">
                    <Cpu className="w-4 h-4" />
                    {health.gpu_name || 'Available'}
                    {health.gpu_memory_gb && ` (${Math.round(health.gpu_memory_gb)}GB)`}
                  </span>
                ) : (
                  <span className="text-sm text-red-600">Not Available</span>
                )}
              </div>

              {/* Models */}
              <div className="flex items-center justify-between">
                <span className="text-sm text-gray-600">Models</span>
                <span className="text-sm font-medium text-gray-900">{getModelsDisplay()}</span>
              </div>
            </div>
          )}

          {/* Not connected warning */}
          {!health && !checking && (
            <div className="bg-red-50 border border-red-200 rounded-lg p-4">
              <p className="text-sm text-red-700">
                Unable to connect to backend. Make sure your RunPod pod is running.
              </p>
            </div>
          )}

          {/* Setup Help Accordion */}
          <div className="border-t border-gray-200 pt-4">
            <button
              onClick={() => setShowHelp(!showHelp)}
              className="flex items-center gap-2 text-sm text-gray-600 hover:text-gray-800 w-full"
            >
              <HelpCircle className="w-4 h-4" />
              <span>Setup Instructions</span>
              {showHelp ? (
                <ChevronUp className="w-4 h-4 ml-auto" />
              ) : (
                <ChevronDown className="w-4 h-4 ml-auto" />
              )}
            </button>

            {showHelp && (
              <div className="mt-4 space-y-3 text-sm">
                <div className="bg-gray-50 rounded-lg p-3 space-y-2 border border-gray-200">
                  <h3 className="font-medium text-gray-900">Quick Setup (RunPod)</h3>
                  <ol className="list-decimal list-inside space-y-1 text-gray-600 text-xs">
                    <li>
                      Deploy a GPU pod (A40 recommended) on{' '}
                      <a
                        href="https://console.runpod.io/deploy"
                        target="_blank"
                        rel="noopener noreferrer"
                        className="text-blue-600 hover:underline"
                      >
                        RunPod <ExternalLink className="w-3 h-3 inline" />
                      </a>
                    </li>
                    <li>
                      Add <code className="bg-gray-200 px-1 rounded">8000</code> to HTTP Ports
                    </li>
                    <li>
                      Open Jupyter Lab and run the setup cell from{' '}
                      <a
                        href="https://github.com/polaires/Banta_Lab_RFdiffusion/blob/main/RUNPOD_QUICK_SETUP.md"
                        target="_blank"
                        rel="noopener noreferrer"
                        className="text-blue-600 hover:underline"
                      >
                        RUNPOD_QUICK_SETUP.md <ExternalLink className="w-3 h-3 inline" />
                      </a>
                    </li>
                    <li>Copy your API URL and paste above</li>
                  </ol>
                </div>

                <div className="bg-gray-50 rounded-lg p-3 space-y-2 border border-gray-200">
                  <h3 className="font-medium text-gray-900">URL Format</h3>
                  <code className="block bg-gray-200 px-2 py-1 rounded text-xs text-gray-800">
                    https://abc123xyz-8000.proxy.runpod.net
                  </code>
                </div>

                <div className="bg-gray-50 rounded-lg p-3 space-y-2 border border-gray-200">
                  <h3 className="font-medium text-gray-900">Troubleshooting</h3>
                  <ul className="space-y-1 text-gray-600 text-xs">
                    <li>
                      <span className="text-amber-600">Connection refused:</span> Re-run setup cell
                    </li>
                    <li>
                      <span className="text-amber-600">CORS error:</span> Check port 8000 is exposed
                    </li>
                    <li>
                      <span className="text-amber-600">Mock mode:</span> Checkpoints still downloading
                    </li>
                  </ul>
                </div>

                <p className="text-xs text-gray-500">
                  Your URL is automatically saved and persists across sessions.
                </p>
              </div>
            )}
          </div>
        </div>

        {/* Footer */}
        <div className="px-5 py-3 bg-gray-50 border-t border-gray-200 rounded-b-xl">
          <p className="text-xs text-gray-500 text-center">
            Connection status auto-refreshes every 30 seconds
          </p>
        </div>
      </div>
    </div>
  );
}
