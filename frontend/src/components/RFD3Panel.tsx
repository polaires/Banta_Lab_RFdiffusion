'use client';

import { useState } from 'react';
import { useStore } from '@/lib/store';
import api from '@/lib/api';
import { Play, Loader2, Info } from 'lucide-react';

const EXAMPLE_CONFIGS = {
  'De novo helix bundle': {
    contig: '100',
    description: 'Design a 100-residue de novo protein',
  },
  'Binder design': {
    contig: 'A1-100/0 100-150',
    description: 'Design a binder for chain A residues 1-100',
  },
  'Scaffold design': {
    contig: 'A10-20/0 30/0 A40-50',
    description: 'Design scaffold connecting fixed regions',
  },
};

export function RFD3Panel() {
  const { health, addJob, updateJob, setSelectedPdb } = useStore();
  const [contig, setContig] = useState('100');
  const [numDesigns, setNumDesigns] = useState(1);
  const [submitting, setSubmitting] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const handleSubmit = async () => {
    if (!health) {
      setError('Backend not connected');
      return;
    }

    setError(null);
    setSubmitting(true);

    try {
      const response = await api.submitRFD3Design({
        contig,
        num_designs: numDesigns,
      });

      addJob({
        id: response.job_id,
        type: 'rfd3',
        status: 'pending',
        createdAt: new Date().toISOString(),
      });

      // Poll for completion
      const result = await api.waitForJob(response.job_id, (status) => {
        updateJob(response.job_id, {
          status: status.status,
          completedAt: status.completed_at,
          result: status.result,
          error: status.error,
        });
      });

      if (result.status === 'completed' && result.result?.designs?.[0]) {
        setSelectedPdb(result.result.designs[0].content);
      }
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to submit job');
    } finally {
      setSubmitting(false);
    }
  };

  return (
    <div className="space-y-6">
      <div>
        <h2 className="text-xl font-bold mb-2">RFdiffusion3 - Structure Design</h2>
        <p className="text-gray-400 text-sm">
          Generate de novo protein structures using diffusion-based generative AI
        </p>
      </div>

      {/* Quick Examples */}
      <div className="space-y-2">
        <label className="text-sm font-medium text-gray-300">Quick Examples</label>
        <div className="flex flex-wrap gap-2">
          {Object.entries(EXAMPLE_CONFIGS).map(([name, config]) => (
            <button
              key={name}
              onClick={() => setContig(config.contig)}
              className="px-3 py-1 text-xs bg-gray-700 hover:bg-gray-600 rounded transition"
              title={config.description}
            >
              {name}
            </button>
          ))}
        </div>
      </div>

      {/* Contig Specification */}
      <div className="space-y-2">
        <label className="text-sm font-medium text-gray-300 flex items-center gap-2">
          Contig Specification
          <Info className="w-4 h-4 text-gray-500" title="Define regions to design" />
        </label>
        <input
          type="text"
          value={contig}
          onChange={(e) => setContig(e.target.value)}
          placeholder="e.g., A1-50/0 50-100"
          className="w-full px-4 py-2 bg-gray-700 rounded border border-gray-600 focus:border-blue-500 focus:outline-none"
        />
        <p className="text-xs text-gray-500">
          Examples: &quot;100&quot; for 100 residues, &quot;A1-50/0 50-100&quot; for binder with 50-100 new residues
        </p>
      </div>

      {/* Number of Designs */}
      <div className="space-y-2">
        <label className="text-sm font-medium text-gray-300">Number of Designs</label>
        <input
          type="number"
          value={numDesigns}
          onChange={(e) => setNumDesigns(Math.max(1, Math.min(10, parseInt(e.target.value) || 1)))}
          min={1}
          max={10}
          className="w-24 px-4 py-2 bg-gray-700 rounded border border-gray-600 focus:border-blue-500 focus:outline-none"
        />
      </div>

      {/* Error Display */}
      {error && (
        <div className="p-3 bg-red-900/50 border border-red-700 rounded text-sm text-red-200">
          {error}
        </div>
      )}

      {/* Submit Button */}
      <button
        onClick={handleSubmit}
        disabled={!health || submitting || !contig}
        className="w-full py-3 bg-blue-600 hover:bg-blue-700 disabled:bg-gray-600 disabled:cursor-not-allowed rounded font-medium flex items-center justify-center gap-2 transition"
      >
        {submitting ? (
          <>
            <Loader2 className="w-5 h-5 animate-spin" />
            Running RFD3...
          </>
        ) : (
          <>
            <Play className="w-5 h-5" />
            Design Structure
          </>
        )}
      </button>

      {!health && (
        <p className="text-sm text-yellow-400 text-center">
          Connect to backend to enable design
        </p>
      )}
    </div>
  );
}
