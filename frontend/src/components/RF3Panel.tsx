'use client';

import { useState, useRef } from 'react';
import { useStore } from '@/lib/store';
import api from '@/lib/api';
import { ConfidenceMetricsDisplay } from './ConfidenceMetrics';
import { ExportPanel } from './ExportPanel';
import type { ConfidenceMetrics, RMSDResult } from '@/lib/api';
import { saveJob as saveJobToSupabase, updateJob as updateJobInSupabase } from '@/lib/supabase';
import { FlaskConical, ArrowRight, ChevronUp, ChevronDown, Library, CheckCircle, X, Upload, Info, AlertCircle, Loader2, Play, Calculator } from 'lucide-react';
import { useAuth } from '@/contexts/AuthContext';

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
  const { user, isConfigured: authConfigured, signInWithGoogle } = useAuth();
  const [sequence, setSequence] = useState('');
  const [submitting, setSubmitting] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [confidences, setConfidences] = useState<ConfidenceMetrics | null>(null);
  const [rmsdResult, setRmsdResult] = useState<RMSDResult | null>(null);
  const [calculatingRmsd, setCalculatingRmsd] = useState(false);
  const [lastPrediction, setLastPrediction] = useState<string | null>(null);

  // MSA support
  const [msaContent, setMsaContent] = useState<string | null>(null);
  const [msaFileName, setMsaFileName] = useState<string | null>(null);
  const [showAdvanced, setShowAdvanced] = useState(false);
  const msaInputRef = useRef<HTMLInputElement>(null);

  // Handle MSA file upload
  const handleMsaUpload = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;

    // Validate file extension
    const ext = file.name.toLowerCase().split('.').pop();
    if (!['a3m', 'fasta', 'fa'].includes(ext || '')) {
      setError('MSA file must be .a3m, .fasta, or .fa format');
      return;
    }

    const reader = new FileReader();
    reader.onload = (event) => {
      const content = event.target?.result as string;
      setMsaContent(content);
      setMsaFileName(file.name);
      setError(null);
    };
    reader.onerror = () => {
      setError('Failed to read MSA file');
    };
    reader.readAsText(file);
  };

  const clearMsa = () => {
    setMsaContent(null);
    setMsaFileName(null);
    if (msaInputRef.current) {
      msaInputRef.current.value = '';
    }
  };

  const handleSubmit = async () => {
    if (authConfigured && !user && process.env.NODE_ENV !== 'development') {
      addNotification({ type: 'error', title: 'Sign in required', message: 'Please sign in to submit designs' });
      signInWithGoogle();
      return;
    }

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
      const response = await api.submitRF3Prediction({
        sequence: cleanSeq,
        msa_content: msaContent || undefined,
      });

      // Add to local store (UI state)
      addJob({
        id: response.job_id,
        type: 'rf3',
        status: 'pending',
        createdAt: new Date().toISOString(),
      });

      // Persist to Supabase (async, non-blocking)
      saveJobToSupabase({
        runpod_id: response.job_id,
        type: 'rf3',
        request: { sequence: cleanSeq },
        user_id: user?.id ?? null,
      });

      const result = await api.waitForJob(response.job_id, (status) => {
        // Update local store
        updateJob(response.job_id, {
          status: status.status,
          completedAt: status.completed_at,
          result: status.result,
          error: status.error,
        });
        // Update Supabase (async, non-blocking)
        updateJobInSupabase(response.job_id, {
          status: status.status,
          completed_at: status.completed_at || null,
          result: status.result || null,
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
        const hasRfd3ForComparison = !!latestRfd3Design?.pdbContent;
        addNotification({
          type: 'success',
          title: 'Structure predicted!',
          message: result.result.confidences
            ? `pLDDT: ${(result.result.confidences.summary_confidences?.overall_plddt || 0 * 100).toFixed(1)}${hasRfd3ForComparison ? '. Compare structures in viewer below.' : '. Design sequences with MPNN.'}`
            : `Fold validated.${hasRfd3ForComparison ? ' Compare structures in viewer below.' : ' Design sequences with MPNN or redesign with RFD3.'}`,
          action: hasRfd3ForComparison ? undefined : {
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
    <div className="space-y-8">
      {/* Header */}
      <div>
        <div className="flex items-center gap-3 mb-3">
          <span className="inline-flex items-center justify-center w-7 h-7 rounded-full bg-emerald-100 text-emerald-700 text-xs font-bold">3</span>
          <h2 className="text-lg font-bold text-foreground">RosettaFold3 — Structure Validation</h2>
        </div>
        <p className="text-muted-foreground text-sm leading-relaxed">
          Validate designability by predicting the structure from MPNN sequences. Compare to RFD3 design via RMSD.
        </p>
      </div>

      {/* Use from RFD3 button */}
      {latestRfd3Design && (
        <div className="p-4 bg-primary/10 border border-primary/20 rounded-xl">
          <div className="flex items-center justify-between">
            <div className="flex items-center gap-3">
              <div className="w-10 h-10 rounded-xl bg-primary/20 flex items-center justify-center">
                <FlaskConical className="w-5 h-5 text-primary" />
              </div>
              <div>
                <p className="text-sm font-semibold text-primary">RFD3 Design Available</p>
                <p className="text-xs text-muted-foreground">
                  Use sequence from your latest RFD3 design for validation
                </p>
              </div>
            </div>
            <button
              onClick={useRfd3Sequence}
              className="px-4 py-2 bg-primary hover:bg-primary/90 text-primary-foreground rounded-lg text-sm font-medium flex items-center gap-2 transition-all shadow-sm hover:shadow-md"
            >
              <ArrowRight className="w-4 h-4" />
              Use Sequence
            </button>
          </div>
        </div>
      )}

      {/* Quick Examples */}
      <div className="space-y-3">
        <label className="text-xs font-bold text-muted-foreground uppercase tracking-wider">Example Sequences</label>
        <div className="flex flex-wrap gap-2">
          {Object.entries(EXAMPLE_SEQUENCES).map(([name, seq]) => (
            <button
              key={name}
              onClick={() => setSequence(seq)}
              className="px-4 py-2 text-xs font-medium bg-muted hover:bg-muted/80 text-foreground rounded-lg transition-colors"
            >
              {name}
            </button>
          ))}
        </div>
      </div>

      {/* Sequence Input */}
      <div className="space-y-3">
        <label className="text-xs font-bold text-muted-foreground uppercase tracking-wider">
          Protein Sequence
        </label>
        <textarea
          value={sequence}
          onChange={(e) => setSequence(e.target.value)}
          placeholder="MSKGEELFTGVVPILVELDGDVNGHKFSVSG..."
          rows={6}
          className="w-full px-4 py-3 bg-muted rounded-xl border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 focus:outline-none font-mono text-sm text-foreground transition-all"
        />
        <div className="flex items-center justify-between">
          <p className="text-xs text-muted-foreground">
            Accepts FASTA or plain sequence
          </p>
          <span className="text-xs font-medium text-foreground bg-muted px-2 py-1 rounded-lg">
            {sequence.replace(/[^A-Za-z]/g, '').length} residues
          </span>
        </div>
      </div>

      {/* Advanced Options Toggle */}
      <button
        onClick={() => setShowAdvanced(!showAdvanced)}
        className="flex items-center gap-2 text-sm text-muted-foreground hover:text-foreground transition-colors"
      >
        {showAdvanced ? <ChevronUp className="w-4 h-4" /> : <ChevronDown className="w-4 h-4" />}
        Advanced Options
      </button>

      {/* Advanced Options (MSA Upload) */}
      {showAdvanced && (
        <div className="space-y-4 p-4 bg-muted rounded-xl border border-border">
          {/* MSA Upload */}
          <div className="space-y-2">
            <label className="text-xs font-bold text-muted-foreground uppercase tracking-wider flex items-center gap-2">
              <Library className="w-4 h-4" />
              Multiple Sequence Alignment (Optional)
            </label>
            <p className="text-xs text-muted-foreground">
              Upload an .a3m or .fasta MSA file to improve prediction quality for natural proteins.
            </p>

            {msaFileName ? (
              <div className="flex items-center justify-between p-3 bg-card rounded-lg border border-emerald-200">
                <div className="flex items-center gap-2">
                  <CheckCircle className="w-5 h-5 text-emerald-600" />
                  <div>
                    <p className="text-sm font-medium text-foreground">{msaFileName}</p>
                    <p className="text-xs text-muted-foreground">
                      {msaContent ? `${msaContent.split('\n').filter(l => l.startsWith('>')).length} sequences` : ''}
                    </p>
                  </div>
                </div>
                <button
                  onClick={clearMsa}
                  className="p-1.5 rounded-lg hover:bg-muted text-muted-foreground hover:text-foreground transition-colors"
                >
                  <X className="w-4 h-4" />
                </button>
              </div>
            ) : (
              <div className="relative">
                <input
                  ref={msaInputRef}
                  type="file"
                  accept=".a3m,.fasta,.fa"
                  onChange={handleMsaUpload}
                  className="absolute inset-0 opacity-0 cursor-pointer"
                />
                <div className="flex items-center justify-center gap-2 p-4 border-2 border-dashed border-border rounded-xl hover:border-primary hover:bg-primary/5 transition-colors cursor-pointer">
                  <Upload className="w-5 h-5 text-muted-foreground" />
                  <span className="text-sm text-muted-foreground">Upload MSA file (.a3m, .fasta)</span>
                </div>
              </div>
            )}

            <div className="p-3 bg-primary/10 rounded-lg border border-primary/20">
              <div className="flex items-start gap-2">
                <Info className="w-4 h-4 text-primary mt-0.5 flex-shrink-0" />
                <div className="text-xs text-primary">
                  <p className="font-medium">When to use MSA:</p>
                  <ul className="mt-1 space-y-0.5 text-primary/80">
                    <li>• Predicting natural protein structures</li>
                    <li>• Improving pLDDT confidence scores</li>
                    <li>• Getting more accurate PAE matrices</li>
                  </ul>
                  <p className="mt-2 text-primary/80">
                    For de novo designed proteins, MSA is typically not needed.
                  </p>
                </div>
              </div>
            </div>
          </div>
        </div>
      )}

      {/* Error Display */}
      {error && (
        <div className="p-4 bg-red-50 border border-red-200 rounded-xl text-sm text-red-700 flex items-start gap-3">
          <AlertCircle className="w-5 h-5 text-red-500 flex-shrink-0" />
          {error}
        </div>
      )}

      {/* Submit Button */}
      <button
        onClick={handleSubmit}
        disabled={!health || submitting || !sequence}
        className="w-full py-4 bg-emerald-600 hover:bg-emerald-700 disabled:bg-muted disabled:text-muted-foreground disabled:cursor-not-allowed text-white rounded-xl font-semibold flex items-center justify-center gap-2 transition-all shadow-sm hover:shadow-md"
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
          className="w-full py-3 bg-primary hover:bg-primary/90 disabled:bg-muted disabled:text-muted-foreground disabled:cursor-not-allowed text-primary-foreground rounded-xl font-medium flex items-center justify-center gap-2 transition-all text-sm"
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
