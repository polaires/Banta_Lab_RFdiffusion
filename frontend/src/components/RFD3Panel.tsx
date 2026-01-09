'use client';

import { useState } from 'react';
import { useStore } from '@/lib/store';
import api from '@/lib/api';
import { DesignTask } from './TaskSelector';
import { ErrorDetails } from './ErrorDetails';
import { saveJob as saveJobToSupabase, updateJob as updateJobInSupabase } from '@/lib/supabase';

// Import task forms
import {
  DeNovoForm,
  ProteinBinderForm,
  SmallMoleculeForm,
  NucleicAcidForm,
  EnzymeForm,
  SymmetricForm,
  RefinementForm,
  RFD3Request,
} from './tasks';

// Task display names
const TASK_NAMES: Record<DesignTask, string> = {
  denovo: 'De Novo Protein',
  protein_binder: 'Protein Binder',
  small_molecule: 'Small Molecule Binder',
  nucleic_acid: 'Nucleic Acid Binder',
  enzyme: 'Enzyme Scaffold',
  symmetric: 'Symmetric Oligomer',
  refinement: 'Structure Refinement',
};

export function RFD3Panel() {
  const {
    health,
    addJob,
    updateJob,
    setSelectedPdb,
    addNotification,
    setLatestDesignPdb,
    setLastCompletedJobType,
    setLatestRfd3Design,
    selectedDesignTask,
    setSelectedDesignTask,
    setActiveTab,
  } = useStore();

  const [submitting, setSubmitting] = useState(false);
  const [error, setError] = useState<{
    message: string;
    errorType?: string;
    traceback?: string;
    context?: Record<string, any>;
  } | null>(null);
  const [lastDesignPdb, setLastDesignPdb] = useState<string | null>(null);

  const handleSubmit = async (request: RFD3Request) => {
    if (!health) {
      setError({ message: 'Backend not connected' });
      return;
    }

    setError(null);
    setSubmitting(true);

    try {
      const response = await api.submitRFD3Design(request as Parameters<typeof api.submitRFD3Design>[0]);

      // Add to local store (UI state)
      addJob({
        id: response.job_id,
        type: 'rfd3',
        status: 'pending',
        createdAt: new Date().toISOString(),
      });

      // Persist to Supabase (async, non-blocking)
      saveJobToSupabase({
        runpod_id: response.job_id,
        type: 'rfd3',
        request: request as Record<string, any>,
      });

      const result = await api.waitForJob(response.job_id, (status) => {
        // Update local store with full error details
        updateJob(response.job_id, {
          status: status.status,
          completedAt: status.completed_at,
          result: status.result,
          error: status.error,
          errorType: status.error_type,
          traceback: status.traceback,
          errorContext: status.context,
        });
        // Update Supabase (async, non-blocking)
        updateJobInSupabase(response.job_id, {
          status: status.status,
          completed_at: status.completed_at || null,
          result: status.result || null,
        });
      });

      if (result.status === 'completed' && result.result?.designs?.[0]) {
        const pdbContent = result.result.designs[0].content;
        const cifContent = result.result.designs[0].cif_content;
        setSelectedPdb(pdbContent);
        setLatestDesignPdb(pdbContent);
        setLastCompletedJobType('rfd3');
        setLastDesignPdb(pdbContent);

        const sequence = extractSequenceFromPdb(pdbContent);

        setLatestRfd3Design({
          jobId: response.job_id,
          pdbContent,
          cifContent,
          source: 'rfd3',
          sequence,
          timestamp: Date.now(),
        });

        addNotification({
          type: 'success',
          title: 'Structure designed!',
          message: result.result.seed
            ? `Backbone ready (seed: ${result.result.seed}). Design sequences with MPNN next.`
            : 'Your protein backbone is ready. Design sequences with MPNN next.',
          action: {
            label: 'Design Sequences',
            tab: 'mpnn',
          },
        });
      } else if (result.status === 'failed') {
        // Set detailed error for display
        setError({
          message: result.error || 'Unknown error occurred',
          errorType: result.error_type,
          traceback: result.traceback,
          context: result.context,
        });
        addNotification({
          type: 'error',
          title: 'Design failed',
          message: result.error || 'Unknown error occurred',
        });
      }
    } catch (err) {
      setError({ message: err instanceof Error ? err.message : 'Failed to submit job' });
      addNotification({
        type: 'error',
        title: 'Submission failed',
        message: err instanceof Error ? err.message : 'Failed to submit job',
      });
    } finally {
      setSubmitting(false);
    }
  };

  const downloadPdb = () => {
    if (!lastDesignPdb) return;
    const blob = new Blob([lastDesignPdb], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'rfd3_design.pdb';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  const formProps = {
    onSubmit: handleSubmit,
    isSubmitting: submitting,
    health,
  };

  // Task form components
  const TaskForms: Record<DesignTask, React.ComponentType<typeof formProps>> = {
    denovo: DeNovoForm,
    protein_binder: ProteinBinderForm,
    small_molecule: SmallMoleculeForm,
    nucleic_acid: NucleicAcidForm,
    enzyme: EnzymeForm,
    symmetric: SymmetricForm,
    refinement: RefinementForm,
  };

  const CurrentForm = selectedDesignTask ? TaskForms[selectedDesignTask as DesignTask] : null;
  const taskName = selectedDesignTask ? TASK_NAMES[selectedDesignTask as DesignTask] : '';

  // If no task selected, redirect to task selection
  if (!selectedDesignTask) {
    return (
      <div className="p-8 space-y-8">
        <div className="border-b border-slate-100 pb-6">
          <h2 className="text-xl font-bold text-slate-900">No Task Selected</h2>
          <p className="text-slate-500 text-sm mt-2">
            Please select a design task first.
          </p>
        </div>
        <button
          onClick={() => setActiveTab('task')}
          className="flex items-center gap-2 px-4 py-2 bg-blue-600 text-white rounded-xl hover:bg-blue-700 transition-colors"
        >
          <span className="material-symbols-outlined">arrow_back</span>
          Choose Design Task
        </button>
      </div>
    );
  }

  return (
    <div className="p-8 space-y-8">
      {/* Header */}
      <div className="border-b border-slate-100 pb-6">
        <div className="flex items-center gap-3 mb-2">
          <span className="bg-blue-50 text-blue-700 border border-blue-100 text-[11px] font-bold px-2.5 py-0.5 rounded-full uppercase tracking-wide">
            Step 1
          </span>
          <h2 className="text-xl font-bold text-slate-900">RFdiffusion3 - {taskName}</h2>
        </div>
        <p className="text-slate-500 text-sm leading-relaxed max-w-3xl pl-1">
          Configure your design parameters. After design, proceed to MPNN for sequence design.
        </p>
      </div>

      <div className="space-y-6">
        {/* Back Button */}
        <button
          onClick={() => {
            setSelectedDesignTask(null);
            setActiveTab('task');
            setError(null);
          }}
          className="flex items-center gap-2 text-slate-600 hover:text-blue-600 transition-colors group"
        >
          <span className="material-symbols-outlined text-lg group-hover:-translate-x-1 transition-transform">
            arrow_back
          </span>
          <span className="text-sm font-medium">Change Design Task</span>
        </button>

        {/* Current Task Form */}
        {CurrentForm && <CurrentForm {...formProps} />}

        {/* Error display */}
        {error && (
          <ErrorDetails
            error={error.message}
            errorType={error.errorType}
            traceback={error.traceback}
            context={error.context}
          />
        )}

        {/* Download button */}
        {lastDesignPdb && (
          <button
            onClick={downloadPdb}
            className="w-full flex justify-center items-center gap-2 py-2.5 px-4 rounded-lg text-sm font-medium text-slate-600 bg-slate-100 hover:bg-slate-200 transition-all"
          >
            <span className="material-symbols-outlined text-lg">download</span>
            Download Design (PDB)
          </button>
        )}
      </div>
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
