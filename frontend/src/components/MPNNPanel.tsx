'use client';

import { useState } from 'react';
import { useStore } from '@/lib/store';
import api from '@/lib/api';
import { Play, Loader2, Upload, Sparkles, Info, Download } from 'lucide-react';

type ModelType = 'ligand_mpnn' | 'protein_mpnn';

export function MPNNPanel() {
  const {
    health,
    addJob,
    updateJob,
    addNotification,
    latestDesignPdb,
    latestRfd3Design,
    setLastCompletedJobType
  } = useStore();
  const [pdbContent, setPdbContent] = useState('');
  const [numSequences, setNumSequences] = useState(8);
  const [temperature, setTemperature] = useState(0.1);
  const [modelType, setModelType] = useState<ModelType>('ligand_mpnn');
  const [removeWaters, setRemoveWaters] = useState(true);
  const [submitting, setSubmitting] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<string | null>(null);
  const [resultModelType, setResultModelType] = useState<string | null>(null);

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
        model_type: modelType,
        remove_waters: removeWaters,
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
        setResultModelType(jobResult.result.model_type || modelType);
        setLastCompletedJobType('mpnn');

        // Notify user with next step suggestion
        addNotification({
          type: 'success',
          title: 'Sequences designed!',
          message: `Generated ${numSequences} sequences using ${modelType}. Predict structure with RF3 to validate.`,
          action: {
            label: 'Predict Structure',
            tab: 'rf3',
          },
        });
      } else if (jobResult.status === 'failed') {
        addNotification({
          type: 'error',
          title: 'Sequence design failed',
          message: jobResult.error || 'Unknown error occurred',
        });
      }
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to submit job');
      addNotification({
        type: 'error',
        title: 'Submission failed',
        message: err instanceof Error ? err.message : 'Failed to submit job',
      });
    } finally {
      setSubmitting(false);
    }
  };

  const handleUseLatestDesign = () => {
    if (latestRfd3Design?.pdbContent) {
      setPdbContent(latestRfd3Design.pdbContent);
      addNotification({
        type: 'info',
        title: 'RFD3 design loaded',
        message: 'Latest RFD3 backbone has been loaded. Adjust parameters and submit.',
      });
    } else if (latestDesignPdb) {
      setPdbContent(latestDesignPdb);
      addNotification({
        type: 'info',
        title: 'Structure loaded',
        message: 'Latest design has been loaded. Adjust parameters and submit.',
      });
    }
  };

  const downloadFasta = () => {
    if (!result) return;
    const blob = new Blob([result], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'designed_sequences.fasta';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  const hasDesignAvailable = latestRfd3Design?.pdbContent || latestDesignPdb;

  return (
    <div className="space-y-6">
      <div>
        <div className="flex items-center gap-2 mb-2">
          <span className="text-xs font-medium px-2 py-0.5 bg-purple-100 text-purple-700 rounded-full">Step 2</span>
          <h2 className="text-xl font-bold">ProteinMPNN - Sequence Design</h2>
        </div>
        <p className="text-gray-600 text-sm">
          Design amino acid sequences for your RFD3 backbone. After design, validate with RF3 structure prediction.
        </p>
      </div>

      {/* Use Latest Design button - shown when a design is available */}
      {hasDesignAvailable && (
        <button
          onClick={handleUseLatestDesign}
          className="w-full py-3 bg-gradient-to-r from-blue-600 to-purple-600 hover:from-blue-700 hover:to-purple-700 rounded-lg font-medium flex items-center justify-center gap-2 transition"
        >
          <Sparkles className="w-5 h-5" />
          Use Latest Design from RFD3/RF3
        </button>
      )}

      {/* PDB Upload */}
      <div className="space-y-2">
        <label className="text-sm font-medium text-gray-700">Input Structure (PDB)</label>

        <div className="flex gap-2">
          <label className="flex-1 py-8 border-2 border-dashed border-gray-300 rounded-lg hover:border-gray-400 cursor-pointer transition flex flex-col items-center justify-center gap-2 bg-gray-50">
            <Upload className="w-8 h-8 text-gray-400" />
            <span className="text-sm text-gray-500">
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
        <label className="text-sm font-medium text-gray-700">Or Paste PDB Content</label>
        <textarea
          value={pdbContent}
          onChange={(e) => setPdbContent(e.target.value)}
          placeholder="ATOM      1  N   ALA A   1      ..."
          rows={6}
          className="w-full px-4 py-2 bg-gray-50 rounded border border-gray-300 focus:border-blue-500 focus:outline-none font-mono text-xs text-gray-900"
        />
      </div>

      {/* Model Type Selection */}
      <div className="space-y-2">
        <label className="text-sm font-medium text-gray-700">Model Type</label>
        <div className="grid grid-cols-2 gap-2">
          <button
            onClick={() => setModelType('ligand_mpnn')}
            className={`p-3 rounded-lg border-2 transition text-left ${
              modelType === 'ligand_mpnn'
                ? 'border-purple-500 bg-purple-50'
                : 'border-gray-200 hover:border-gray-300 bg-white'
            }`}
          >
            <div className="font-medium text-gray-900">LigandMPNN</div>
            <div className="text-xs text-gray-500">Ligand-aware design (recommended)</div>
          </button>
          <button
            onClick={() => setModelType('protein_mpnn')}
            className={`p-3 rounded-lg border-2 transition text-left ${
              modelType === 'protein_mpnn'
                ? 'border-purple-500 bg-purple-50'
                : 'border-gray-200 hover:border-gray-300 bg-white'
            }`}
          >
            <div className="font-medium text-gray-900">ProteinMPNN</div>
            <div className="text-xs text-gray-500">Original protein-only model</div>
          </button>
        </div>
        <div className="flex items-start gap-2 text-xs text-gray-500">
          <Info className="w-3 h-3 mt-0.5 flex-shrink-0" />
          <span>LigandMPNN is recommended for most use cases. Use ProteinMPNN for pure protein design without ligands.</span>
        </div>
      </div>

      {/* Parameters */}
      <div className="grid grid-cols-2 gap-4">
        <div className="space-y-2">
          <label className="text-sm font-medium text-gray-700">Number of Sequences</label>
          <input
            type="number"
            value={numSequences}
            onChange={(e) => setNumSequences(Math.max(1, Math.min(100, parseInt(e.target.value) || 8)))}
            min={1}
            max={100}
            className="w-full px-4 py-2 bg-gray-50 rounded border border-gray-300 focus:border-blue-500 focus:outline-none text-gray-900"
          />
        </div>
        <div className="space-y-2">
          <label className="text-sm font-medium text-gray-700">Sampling Temperature</label>
          <input
            type="number"
            value={temperature}
            onChange={(e) => setTemperature(Math.max(0.01, Math.min(2, parseFloat(e.target.value) || 0.1)))}
            min={0.01}
            max={2}
            step={0.01}
            className="w-full px-4 py-2 bg-gray-50 rounded border border-gray-300 focus:border-blue-500 focus:outline-none text-gray-900"
          />
          <p className="text-xs text-gray-500">Lower = more conservative, higher = more diverse</p>
        </div>
      </div>

      {/* Remove Waters Toggle */}
      <div className="flex items-center gap-2">
        <input
          type="checkbox"
          id="removeWaters"
          checked={removeWaters}
          onChange={(e) => setRemoveWaters(e.target.checked)}
          className="w-4 h-4 rounded border-gray-300 bg-white text-purple-600 focus:ring-purple-500"
        />
        <label htmlFor="removeWaters" className="text-sm text-gray-700">
          Remove water molecules from context
        </label>
      </div>

      {/* Error Display */}
      {error && (
        <div className="p-3 bg-red-50 border border-red-200 rounded text-sm text-red-700">
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
          <div className="flex items-center justify-between">
            <label className="text-sm font-medium text-gray-700">
              Designed Sequences (FASTA)
              {resultModelType && (
                <span className="ml-2 text-xs text-gray-500">via {resultModelType}</span>
              )}
            </label>
            <button
              onClick={downloadFasta}
              className="px-2 py-1 text-xs bg-gray-100 hover:bg-gray-200 text-gray-700 rounded flex items-center gap-1 transition"
            >
              <Download className="w-3 h-3" />
              Download
            </button>
          </div>
          <pre className="p-4 bg-gray-900 rounded text-xs font-mono overflow-auto max-h-64 text-green-400">
            {result}
          </pre>
        </div>
      )}
    </div>
  );
}
