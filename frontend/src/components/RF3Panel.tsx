'use client';

import { useState } from 'react';
import { useStore } from '@/lib/store';
import api from '@/lib/api';
import { Play, Loader2 } from 'lucide-react';

const EXAMPLE_SEQUENCES = {
  'GFP (partial)': 'MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK',
  'Ubiquitin': 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG',
  'Insulin B chain': 'FVNQHLCGSHLVEALYLVCGERGFFYTPKT',
};

export function RF3Panel() {
  const {
    health,
    addJob,
    updateJob,
    setSelectedPdb,
    addNotification,
    setLatestDesignPdb,
    setLastCompletedJobType
  } = useStore();
  const [sequence, setSequence] = useState('');
  const [submitting, setSubmitting] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const handleSubmit = async () => {
    if (!health) {
      setError('Backend not connected');
      return;
    }

    // Clean sequence
    const cleanSeq = sequence.replace(/[^A-Za-z]/g, '').toUpperCase();
    if (cleanSeq.length < 10) {
      setError('Sequence must be at least 10 residues');
      return;
    }

    setError(null);
    setSubmitting(true);

    try {
      const response = await api.submitRF3Prediction({ sequence: cleanSeq });

      addJob({
        id: response.job_id,
        type: 'rf3',
        status: 'pending',
        createdAt: new Date().toISOString(),
      });

      const result = await api.waitForJob(response.job_id, (status) => {
        updateJob(response.job_id, {
          status: status.status,
          completedAt: status.completed_at,
          result: status.result,
          error: status.error,
        });
      });

      if (result.status === 'completed' && result.result?.predictions?.[0]) {
        const pdbContent = result.result.predictions[0].content;
        setSelectedPdb(pdbContent);
        setLatestDesignPdb(pdbContent);
        setLastCompletedJobType('rf3');

        // Notify user with next step suggestion
        addNotification({
          type: 'success',
          title: 'Structure predicted!',
          message: 'Fold validated. Design sequences with MPNN or redesign with RFD3.',
          action: {
            label: 'Design Sequences',
            tab: 'mpnn',
          },
        });
      } else if (result.status === 'failed') {
        addNotification({
          type: 'error',
          title: 'Prediction failed',
          message: result.error || 'Unknown error occurred',
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

  return (
    <div className="space-y-6">
      <div>
        <h2 className="text-xl font-bold mb-2">RosettaFold3 - Structure Prediction</h2>
        <p className="text-gray-400 text-sm">
          Predict 3D protein structure from amino acid sequence
        </p>
      </div>

      {/* Quick Examples */}
      <div className="space-y-2">
        <label className="text-sm font-medium text-gray-300">Example Sequences</label>
        <div className="flex flex-wrap gap-2">
          {Object.entries(EXAMPLE_SEQUENCES).map(([name, seq]) => (
            <button
              key={name}
              onClick={() => setSequence(seq)}
              className="px-3 py-1 text-xs bg-gray-700 hover:bg-gray-600 rounded transition"
            >
              {name}
            </button>
          ))}
        </div>
      </div>

      {/* Sequence Input */}
      <div className="space-y-2">
        <label className="text-sm font-medium text-gray-300">
          Protein Sequence (FASTA or plain)
        </label>
        <textarea
          value={sequence}
          onChange={(e) => setSequence(e.target.value)}
          placeholder="MSKGEELFTGVVPILVELDGDVNGHKFSVSG..."
          rows={6}
          className="w-full px-4 py-2 bg-gray-700 rounded border border-gray-600 focus:border-blue-500 focus:outline-none font-mono text-sm"
        />
        <p className="text-xs text-gray-500">
          {sequence.replace(/[^A-Za-z]/g, '').length} residues
        </p>
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
        disabled={!health || submitting || !sequence}
        className="w-full py-3 bg-green-600 hover:bg-green-700 disabled:bg-gray-600 disabled:cursor-not-allowed rounded font-medium flex items-center justify-center gap-2 transition"
      >
        {submitting ? (
          <>
            <Loader2 className="w-5 h-5 animate-spin" />
            Predicting Structure...
          </>
        ) : (
          <>
            <Play className="w-5 h-5" />
            Predict Structure
          </>
        )}
      </button>
    </div>
  );
}
