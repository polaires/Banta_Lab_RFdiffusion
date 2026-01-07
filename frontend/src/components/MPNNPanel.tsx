'use client';

import { useState } from 'react';
import { useStore } from '@/lib/store';
import api from '@/lib/api';
import { Play, Loader2, Upload } from 'lucide-react';

export function MPNNPanel() {
  const { health, addJob, updateJob } = useStore();
  const [pdbContent, setPdbContent] = useState('');
  const [numSequences, setNumSequences] = useState(8);
  const [temperature, setTemperature] = useState(0.1);
  const [submitting, setSubmitting] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<string | null>(null);

  const handleFileUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (file) {
      const reader = new FileReader();
      reader.onload = (e) => {
        setPdbContent(e.target?.result as string);
      };
      reader.readAsText(file);
    }
  };

  const handleSubmit = async () => {
    if (!health) {
      setError('Backend not connected');
      return;
    }

    if (!pdbContent || pdbContent.length < 100) {
      setError('Please upload or paste a valid PDB file');
      return;
    }

    setError(null);
    setResult(null);
    setSubmitting(true);

    try {
      const response = await api.submitMPNNDesign({
        pdb_content: pdbContent,
        num_sequences: numSequences,
        temperature,
      });

      addJob({
        id: response.job_id,
        type: 'mpnn',
        status: 'pending',
        createdAt: new Date().toISOString(),
      });

      const jobResult = await api.waitForJob(response.job_id, (status) => {
        updateJob(response.job_id, {
          status: status.status,
          completedAt: status.completed_at,
          result: status.result,
          error: status.error,
        });
      });

      if (jobResult.status === 'completed' && jobResult.result?.sequences?.[0]) {
        setResult(jobResult.result.sequences[0].content);
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
        <h2 className="text-xl font-bold mb-2">ProteinMPNN - Sequence Design</h2>
        <p className="text-gray-400 text-sm">
          Design amino acid sequences for a given protein backbone structure
        </p>
      </div>

      {/* PDB Upload */}
      <div className="space-y-2">
        <label className="text-sm font-medium text-gray-300">Input Structure (PDB)</label>
        <div className="flex gap-2">
          <label className="flex-1 py-8 border-2 border-dashed border-gray-600 rounded-lg hover:border-gray-500 cursor-pointer transition flex flex-col items-center justify-center gap-2">
            <Upload className="w-8 h-8 text-gray-500" />
            <span className="text-sm text-gray-400">
              {pdbContent ? 'File loaded' : 'Click to upload PDB'}
            </span>
            <input
              type="file"
              accept=".pdb"
              onChange={handleFileUpload}
              className="hidden"
            />
          </label>
        </div>
      </div>

      {/* Or paste PDB */}
      <div className="space-y-2">
        <label className="text-sm font-medium text-gray-300">Or Paste PDB Content</label>
        <textarea
          value={pdbContent}
          onChange={(e) => setPdbContent(e.target.value)}
          placeholder="ATOM      1  N   ALA A   1      ..."
          rows={6}
          className="w-full px-4 py-2 bg-gray-700 rounded border border-gray-600 focus:border-blue-500 focus:outline-none font-mono text-xs"
        />
      </div>

      {/* Parameters */}
      <div className="grid grid-cols-2 gap-4">
        <div className="space-y-2">
          <label className="text-sm font-medium text-gray-300">Number of Sequences</label>
          <input
            type="number"
            value={numSequences}
            onChange={(e) => setNumSequences(Math.max(1, Math.min(32, parseInt(e.target.value) || 8)))}
            min={1}
            max={32}
            className="w-full px-4 py-2 bg-gray-700 rounded border border-gray-600 focus:border-blue-500 focus:outline-none"
          />
        </div>
        <div className="space-y-2">
          <label className="text-sm font-medium text-gray-300">Sampling Temperature</label>
          <input
            type="number"
            value={temperature}
            onChange={(e) => setTemperature(Math.max(0.01, Math.min(2, parseFloat(e.target.value) || 0.1)))}
            min={0.01}
            max={2}
            step={0.01}
            className="w-full px-4 py-2 bg-gray-700 rounded border border-gray-600 focus:border-blue-500 focus:outline-none"
          />
          <p className="text-xs text-gray-500">Lower = more conservative, higher = more diverse</p>
        </div>
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
        disabled={!health || submitting || !pdbContent}
        className="w-full py-3 bg-purple-600 hover:bg-purple-700 disabled:bg-gray-600 disabled:cursor-not-allowed rounded font-medium flex items-center justify-center gap-2 transition"
      >
        {submitting ? (
          <>
            <Loader2 className="w-5 h-5 animate-spin" />
            Designing Sequences...
          </>
        ) : (
          <>
            <Play className="w-5 h-5" />
            Design Sequences
          </>
        )}
      </button>

      {/* Results */}
      {result && (
        <div className="space-y-2">
          <label className="text-sm font-medium text-gray-300">Designed Sequences (FASTA)</label>
          <pre className="p-4 bg-gray-900 rounded text-xs font-mono overflow-auto max-h-64 text-green-400">
            {result}
          </pre>
        </div>
      )}
    </div>
  );
}
