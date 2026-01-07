'use client';

import { useState } from 'react';
import { useStore } from '@/lib/store';
import api from '@/lib/api';
import { Play, Loader2, ArrowRight, Calculator } from 'lucide-react';
import { ConfidenceMetricsDisplay } from './ConfidenceMetrics';
import { ExportPanel } from './ExportPanel';
import type { ConfidenceMetrics, RMSDResult } from '@/lib/api';

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
    setLastCompletedJobType,
    setLatestRf3Prediction,
    setLatestConfidences,
    latestRfd3Design,
    setLatestRmsdResult,
  } = useStore();
  const [sequence, setSequence] = useState('');
  const [submitting, setSubmitting] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [confidences, setConfidences] = useState<ConfidenceMetrics | null>(null);
  const [rmsdResult, setRmsdResult] = useState<RMSDResult | null>(null);
  const [calculatingRmsd, setCalculatingRmsd] = useState(false);
  const [lastPrediction, setLastPrediction] = useState<string | null>(null);

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
    setConfidences(null);
    setRmsdResult(null);

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
        setLastPrediction(pdbContent);

        // Store confidences if available
        if (result.result.confidences) {
          setConfidences(result.result.confidences);
          setLatestConfidences(result.result.confidences);
        }

        // Store for cross-panel data flow
        setLatestRf3Prediction({
          jobId: response.job_id,
          pdbContent,
          cifContent: result.result.predictions[0].cif_content,
          source: 'rf3',
          confidences: result.result.confidences,
          timestamp: Date.now(),
        });

        // Notify user with next step suggestion
        addNotification({
          type: 'success',
          title: 'Structure predicted!',
          message: result.result.confidences
            ? `pLDDT: ${(result.result.confidences.summary_confidences?.overall_plddt || 0 * 100).toFixed(1)}. Design sequences with MPNN.`
            : 'Fold validated. Design sequences with MPNN or redesign with RFD3.',
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

  const handleCalculateRMSD = async () => {
    if (!latestRfd3Design?.pdbContent || !lastPrediction) {
      setError('Need both RFD3 design and RF3 prediction to calculate RMSD');
      return;
    }

    setCalculatingRmsd(true);
    setError(null);

    try {
      const result = await api.calculateRMSD({
        pdb_content_1: latestRfd3Design.pdbContent,
        pdb_content_2: lastPrediction,
        backbone_only: true,
      });

      setRmsdResult(result);
      setLatestRmsdResult(result);

      addNotification({
        type: result.interpretation === 'Excellent' || result.interpretation === 'Good' ? 'success' : 'info',
        title: `RMSD: ${result.rmsd.toFixed(2)}A - ${result.interpretation}`,
        message: result.description,
      });
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to calculate RMSD');
    } finally {
      setCalculatingRmsd(false);
    }
  };

  // Extract sequence from latest RFD3 design
  const useRfd3Sequence = () => {
    if (latestRfd3Design?.sequence) {
      setSequence(latestRfd3Design.sequence);
    } else if (latestRfd3Design?.pdbContent) {
      // Extract sequence from PDB
      const seqFromPdb = extractSequenceFromPdb(latestRfd3Design.pdbContent);
      setSequence(seqFromPdb);
    }
  };

  return (
    <div className="space-y-6">
      <div>
        <h2 className="text-xl font-bold mb-2">RosettaFold3 - Structure Prediction</h2>
        <p className="text-gray-400 text-sm">
          Predict 3D protein structure from amino acid sequence. Returns confidence metrics (pLDDT, PAE, pTM).
        </p>
      </div>

      {/* Use from RFD3 button */}
      {latestRfd3Design && (
        <div className="p-3 bg-blue-900/30 border border-blue-700 rounded">
          <div className="flex items-center justify-between">
            <div>
              <p className="text-sm font-medium text-blue-300">RFD3 Design Available</p>
              <p className="text-xs text-gray-400">
                Use sequence from your latest RFD3 design for validation
              </p>
            </div>
            <button
              onClick={useRfd3Sequence}
              className="px-3 py-1.5 bg-blue-600 hover:bg-blue-700 rounded text-sm flex items-center gap-1 transition"
            >
              <ArrowRight className="w-4 h-4" />
              Use Sequence
            </button>
          </div>
        </div>
      )}

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

      {/* RMSD Validation Button */}
      {lastPrediction && latestRfd3Design && (
        <button
          onClick={handleCalculateRMSD}
          disabled={calculatingRmsd}
          className="w-full py-2 bg-purple-600 hover:bg-purple-700 disabled:bg-gray-600 disabled:cursor-not-allowed rounded font-medium flex items-center justify-center gap-2 transition text-sm"
        >
          {calculatingRmsd ? (
            <>
              <Loader2 className="w-4 h-4 animate-spin" />
              Calculating RMSD...
            </>
          ) : (
            <>
              <Calculator className="w-4 h-4" />
              Validate Design (Calculate RMSD vs RFD3)
            </>
          )}
        </button>
      )}

      {/* Confidence Metrics Display */}
      <ConfidenceMetricsDisplay confidences={confidences} rmsdResult={rmsdResult} />

      {/* Export Panel */}
      {lastPrediction && (
        <ExportPanel
          pdbContent={lastPrediction}
          source="rf3"
          filename="rf3_prediction"
        />
      )}
    </div>
  );
}

// Helper to extract sequence from PDB content
function extractSequenceFromPdb(pdbContent: string): string {
  const threeToOne: Record<string, string> = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
  };

  const residues: { resSeq: number; resName: string }[] = [];
  const seenResidues = new Set<string>();

  for (const line of pdbContent.split('\n')) {
    if (line.startsWith('ATOM') && line.includes(' CA ')) {
      const resName = line.substring(17, 20).trim();
      const resSeq = parseInt(line.substring(22, 26).trim());
      const key = `${resSeq}`;

      if (!seenResidues.has(key)) {
        seenResidues.add(key);
        residues.push({ resSeq, resName });
      }
    }
  }

  residues.sort((a, b) => a.resSeq - b.resSeq);

  return residues
    .map(r => threeToOne[r.resName] || 'X')
    .join('');
}
