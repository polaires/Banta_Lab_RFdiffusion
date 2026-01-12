'use client';

import { useEffect, useState } from 'react';
import { useStore } from '@/lib/store';
import api from '@/lib/api';
import { Wifi, WifiOff, Cpu, Check, X, AlertCircle, Sparkles, HelpCircle, ChevronDown, ChevronUp, ExternalLink } from 'lucide-react';

export function ConnectionStatus() {
  const { backendUrl, setBackendUrl, health, setHealth, aiCaseStudy } = useStore();
  const [inputUrl, setInputUrl] = useState(backendUrl);
  const [checking, setChecking] = useState(false);
  const [showHelp, setShowHelp] = useState(false);

  // Check if a design job is currently running (either by job ID or workflow phase)
  const isJobRunning =
    (aiCaseStudy?.pendingJobId !== null && aiCaseStudy?.pendingJobId !== undefined) ||
    aiCaseStudy?.workflowPhase === 'running';

  const checkConnection = async () => {
    setChecking(true);
    try {
      api.setBaseUrl(inputUrl);
      const healthResponse = await api.checkHealth();
      setHealth(healthResponse);
      setBackendUrl(inputUrl);
    } catch {
      // Don't mark as disconnected if a job is running - backend may just be busy
      // Only set health to null if no job is running
      if (!isJobRunning) {
        setHealth(null);
      }
      // If job is running, keep previous health status (assume still connected)
    } finally {
      setChecking(false);
    }
  };

  useEffect(() => {
    checkConnection();
    const interval = setInterval(checkConnection, 30000); // Check every 30s
    return () => clearInterval(interval);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  // Sync inputUrl when backendUrl changes (e.g., from localStorage)
  useEffect(() => {
    setInputUrl(backendUrl);
  }, [backendUrl]);

  const isConnected = health?.status === 'healthy';

  // Handle both old format (models_loaded array) and new format (models object)
  const getModelsDisplay = () => {
    if (!health) return 'None';

    // New format: models object with availability info
    if (health.models && typeof health.models === 'object') {
      const availableModels = Object.entries(health.models)
        .filter(([, info]) => (info as { available?: boolean })?.available)
        .map(([name]) => name.toUpperCase());
      return availableModels.length > 0 ? availableModels.join(', ') : 'None (mock mode)';
    }

    // Old format: models_loaded array
    if (health.models_loaded && Array.isArray(health.models_loaded)) {
      return health.models_loaded.join(', ') || 'None';
    }

    return 'None';
  };

  const isMockMode = health?.mode === 'mock';

  return (
    <div className="bg-gray-800 rounded-lg p-4 space-y-4">
      <h2 className="text-lg font-semibold flex items-center gap-2">
        {isConnected ? (
          <Wifi className="w-5 h-5 text-green-400" />
        ) : (
          <WifiOff className="w-5 h-5 text-red-400" />
        )}
        Backend Connection
      </h2>

      <div className="flex gap-2">
        <input
          type="text"
          value={inputUrl}
          onChange={(e) => setInputUrl(e.target.value)}
          placeholder="https://your-runpod-endpoint:8000"
          className="flex-1 px-3 py-2 bg-gray-700 rounded border border-gray-600 focus:border-blue-500 focus:outline-none text-sm"
        />
        <button
          onClick={checkConnection}
          disabled={checking}
          className="px-4 py-2 bg-blue-600 hover:bg-blue-700 disabled:bg-gray-600 rounded text-sm font-medium transition"
        >
          {checking ? 'Checking...' : 'Connect'}
        </button>
      </div>

      {health && (
        <div className="space-y-2 text-sm">
          {/* Mode indicator */}
          {health.mode && (
            <div className="flex items-center gap-2">
              {isMockMode ? (
                <AlertCircle className="w-4 h-4 text-yellow-400" />
              ) : (
                <Sparkles className="w-4 h-4 text-green-400" />
              )}
              <span>Mode:</span>
              <span className={isMockMode ? 'text-yellow-400' : 'text-green-400'}>
                {isMockMode ? 'Mock (Demo)' : 'Real (Foundry)'}
              </span>
            </div>
          )}

          {/* GPU info */}
          <div className="flex items-center gap-2">
            <Cpu className="w-4 h-4" />
            <span>GPU:</span>
            {health.gpu_available ? (
              <span className="text-green-400 flex items-center gap-1">
                <Check className="w-4 h-4" />
                {health.gpu_name || 'Available'}
                {health.gpu_memory_gb && ` (${Math.round(health.gpu_memory_gb)}GB)`}
              </span>
            ) : (
              <span className="text-red-400 flex items-center gap-1">
                <X className="w-4 h-4" /> Not Available
              </span>
            )}
          </div>

          {/* Models */}
          <div>
            <span className="text-gray-400">Models: </span>
            <span>{getModelsDisplay()}</span>
          </div>
        </div>
      )}

      {!health && !checking && (
        <div className="space-y-2">
          <p className="text-sm text-red-400">
            Unable to connect to backend. Make sure your RunPod pod is running.
          </p>
          <button
            onClick={() => setShowHelp(true)}
            className="text-sm text-blue-400 hover:text-blue-300 flex items-center gap-1"
          >
            <HelpCircle className="w-4 h-4" />
            Need help setting up?
          </button>
        </div>
      )}

      {/* Setup Help Section */}
      <div className="border-t border-gray-700 pt-3 mt-3">
        <button
          onClick={() => setShowHelp(!showHelp)}
          className="flex items-center gap-2 text-sm text-gray-400 hover:text-gray-300 w-full"
        >
          <HelpCircle className="w-4 h-4" />
          <span>Setup Instructions</span>
          {showHelp ? <ChevronUp className="w-4 h-4 ml-auto" /> : <ChevronDown className="w-4 h-4 ml-auto" />}
        </button>

        {showHelp && (
          <div className="mt-3 space-y-3 text-sm">
            <div className="bg-gray-700/50 rounded p-3 space-y-2">
              <h3 className="font-medium text-gray-200">Quick Setup (RunPod)</h3>
              <ol className="list-decimal list-inside space-y-1 text-gray-400">
                <li>Deploy a GPU pod (A40 recommended) on <a href="https://console.runpod.io/deploy" target="_blank" rel="noopener noreferrer" className="text-blue-400 hover:underline">RunPod <ExternalLink className="w-3 h-3 inline" /></a></li>
                <li>Add <code className="bg-gray-800 px-1 rounded">8000</code> to HTTP Ports</li>
                <li>Open Jupyter Lab and run the setup cell from <a href="https://github.com/polaires/Banta_Lab_RFdiffusion/blob/main/RUNPOD_QUICK_SETUP.md" target="_blank" rel="noopener noreferrer" className="text-blue-400 hover:underline">RUNPOD_QUICK_SETUP.md <ExternalLink className="w-3 h-3 inline" /></a></li>
                <li>Copy your API URL: <code className="bg-gray-800 px-1 rounded text-xs">https://POD_ID-8000.proxy.runpod.net</code></li>
                <li>Paste above and click Connect</li>
              </ol>
            </div>

            <div className="bg-gray-700/50 rounded p-3 space-y-2">
              <h3 className="font-medium text-gray-200">URL Format</h3>
              <p className="text-gray-400">Your backend URL should look like:</p>
              <code className="block bg-gray-800 px-2 py-1 rounded text-xs text-green-400">
                https://abc123xyz-8000.proxy.runpod.net
              </code>
              <p className="text-xs text-gray-500">Find this in your RunPod pod&apos;s Connect menu</p>
            </div>

            <div className="bg-gray-700/50 rounded p-3 space-y-2">
              <h3 className="font-medium text-gray-200">Troubleshooting</h3>
              <ul className="space-y-1 text-gray-400">
                <li><span className="text-yellow-400">Connection refused:</span> Re-run the setup cell in Jupyter</li>
                <li><span className="text-yellow-400">CORS error:</span> Check that port 8000 is in HTTP Ports</li>
                <li><span className="text-yellow-400">Mock mode:</span> Checkpoints may still be downloading</li>
              </ul>
            </div>

            <p className="text-xs text-gray-500">
              Your URL is automatically saved and will persist across browser sessions.
            </p>
          </div>
        )}
      </div>
    </div>
  );
}
