'use client';

import { useEffect, useState } from 'react';
import { useStore } from '@/lib/store';
import api from '@/lib/api';
import { Wifi, WifiOff, Cpu, Check, X, AlertCircle, Sparkles } from 'lucide-react';

export function ConnectionStatus() {
  const { backendUrl, setBackendUrl, health, setHealth } = useStore();
  const [inputUrl, setInputUrl] = useState(backendUrl);
  const [checking, setChecking] = useState(false);

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

  useEffect(() => {
    checkConnection();
    const interval = setInterval(checkConnection, 30000); // Check every 30s
    return () => clearInterval(interval);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

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
        <p className="text-sm text-red-400">
          Unable to connect to backend. Make sure your RunPod pod is running.
        </p>
      )}
    </div>
  );
}
