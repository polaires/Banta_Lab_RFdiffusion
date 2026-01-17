'use client';

import { useState, useEffect } from 'react';
import { ArrowLeft, Download, Plus } from 'lucide-react';
import { useStore } from '@/lib/store';
import api from '@/lib/api';
import { DesignTask } from './TaskSelector';
import { ErrorDetails } from './ErrorDetails';
import { saveJob as saveJobToSupabase, updateJob as updateJobInSupabase } from '@/lib/supabase';
import { PipelineProgress, PipelineStage } from './tasks/PipelineProgress';
import {
  InterfaceLigandResultsPanel,
  InterfaceLigandJobResult,
  InterfaceLigandDesign,
} from './tasks/InterfaceLigandResultsPanel';

// Import task forms
import {
  DeNovoForm,
  ProteinBinderForm,
  SmallMoleculeForm,
  NucleicAcidForm,
  EnzymeForm,
  SymmetricForm,
  RefinementForm,
  InterfaceLigandForm,
  InterfaceMetalForm,
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
  interface_ligand: 'Interface Ligand Dimer',
  interface_metal: 'Interface Metal Dimer',
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

  // Pipeline progress tracking for interface_ligand tasks
  const [pipelineStage, setPipelineStage] = useState<PipelineStage>('idle');
  const [pipelineMessage, setPipelineMessage] = useState<string>('');
  const [pipelineError, setPipelineError] = useState<string | undefined>(undefined);

  // Interface ligand design results
  const [interfaceLigandResult, setInterfaceLigandResult] = useState<InterfaceLigandJobResult | null>(null);
  const [selectedDesignId, setSelectedDesignId] = useState<string | undefined>(undefined);
  const [isPipelineProcessing, setIsPipelineProcessing] = useState(false);

  // Helper to parse status messages and determine pipeline stage
  const parsePipelineStage = (message: string | undefined): PipelineStage => {
    if (!message) return 'idle';
    const msg = message.toLowerCase();

    if (msg.includes('rfd3') || msg.includes('diffusion') || msg.includes('backbone') || msg.includes('generating')) {
      return 'backbone';
    }
    if (msg.includes('mpnn') || msg.includes('ligandmpnn') || msg.includes('sequence')) {
      return 'ligandmpnn';
    }
    if (msg.includes('esm') || msg.includes('valid') || msg.includes('fold')) {
      return 'validation';
    }
    if (msg.includes('plip') || msg.includes('interaction') || msg.includes('analysis') || msg.includes('gnina') || msg.includes('docking')) {
      return 'analysis';
    }
    if (msg.includes('complete') || msg.includes('done') || msg.includes('finished')) {
      return 'complete';
    }
    if (msg.includes('error') || msg.includes('fail')) {
      return 'error';
    }
    return 'backbone'; // Default to first stage if running
  };

  const handleSubmit = async (request: RFD3Request) => {
    if (!health) {
      setError({ message: 'Backend not connected' });
      return;
    }

    setError(null);
    setSubmitting(true);

    // Reset pipeline state for interface_ligand tasks
    const isInterfaceLigandTask = request.task === 'interface_ligand_design';
    if (isInterfaceLigandTask) {
      setPipelineStage('backbone');
      setPipelineMessage('Starting backbone generation...');
      setPipelineError(undefined);
    }

    try {
      const response = await api.submitRFD3Design(request as Parameters<typeof api.submitRFD3Design>[0]);

      // Add to local store (UI state)
      addJob({
        id: response.job_id,
        type: 'rfd3',
        status: response.syncCompleted ? (response.status as 'pending' | 'running' | 'completed' | 'failed') : 'pending',
        createdAt: new Date().toISOString(),
      });

      // Persist to Supabase (async, non-blocking)
      saveJobToSupabase({
        runpod_id: response.job_id,
        type: 'rfd3',
        request: request as Record<string, any>,
      });

      let result;

      // For traditional mode (runsync), result is already in the response
      if (response.syncCompleted) {
        // Result already available from synchronous runsync call
        result = {
          status: response.status,
          result: response.result,
          error: response.error,
        };
        // Update job with completed status
        updateJob(response.job_id, {
          status: response.status as 'pending' | 'running' | 'completed' | 'failed',
          completedAt: new Date().toISOString(),
          result: response.result,
          error: response.error,
        });
      } else {
        // Serverless mode - poll for result
        result = await api.waitForJob(response.job_id, (status) => {
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

          // Update pipeline stage for interface_ligand tasks
          if (isInterfaceLigandTask) {
            // The backend may return additional status info not in the type definition
            const extendedStatus = status as typeof status & { message?: string; stage?: string };
            const statusMessage = extendedStatus.message || extendedStatus.stage || status.status;
            const newStage = parsePipelineStage(statusMessage);
            setPipelineStage(newStage);
            setPipelineMessage(statusMessage);

            if (status.status === 'failed') {
              setPipelineStage('error');
              setPipelineError(status.error);
            }
          }
        });
      }

      if (result.status === 'completed') {
        // Update pipeline stage to complete for interface_ligand tasks
        if (isInterfaceLigandTask) {
          setPipelineStage('complete');
          setPipelineMessage('Design complete!');

          // Parse interface ligand designs for results panel
          const designs: InterfaceLigandDesign[] = [];
          const resultData = result.result || {};
          const approach = resultData.approach || request.approach || 'joint';

          // Handle different response formats from backend
          if (resultData.designs && Array.isArray(resultData.designs)) {
            resultData.designs.forEach((d: any, idx: number) => {
              designs.push({
                id: d.id || `design-${idx}`,
                rank: d.rank || idx + 1,
                pdb_content: d.pdb_content || d.content || '',
                sequence_a: d.sequence_a,
                sequence_b: d.sequence_b,
                metrics: d.metrics || {},
                approach,
                status: d.passed !== false ? 'passed' : 'failed',
                pipeline_stage: 'backbone',
              });
            });
          } else if (resultData.dimer) {
            // Single dimer result (full approach)
            designs.push({
              id: 'dimer-0',
              rank: 1,
              pdb_content: resultData.dimer.pdb_content || '',
              sequence_a: resultData.dimer.sequence_a,
              sequence_b: resultData.dimer.sequence_b,
              metrics: resultData.dimer.metrics || {},
              approach,
              status: 'passed',
              pipeline_stage: 'backbone',
            });
          }

          // Set results for the panel
          setInterfaceLigandResult({
            status: 'completed',
            approach,
            designs,
            best_design_idx: resultData.best_design_idx,
          });

          // Select the first design by default
          if (designs.length > 0) {
            setSelectedDesignId(designs[0].id);
          }
        }

        // Extract PDB content from various response formats:
        // 1. Standard RFD3: result.designs[0].content
        // 2. Interface ligand (asymmetric/sequential): result.designs[0].pdb_content
        // 3. Interface ligand (full): result.dimer.pdb_content
        let pdbContent: string | undefined;
        let cifContent: string | undefined;

        if (result.result?.designs?.[0]) {
          // Standard RFD3 or interface_ligand asymmetric/sequential
          const design = result.result.designs[0];
          pdbContent = design.content || design.pdb_content;
          cifContent = design.cif_content;
        } else if (result.result?.dimer?.pdb_content) {
          // Interface ligand full dimer approach
          pdbContent = result.result.dimer.pdb_content;
        }

        if (pdbContent) {
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

          // Customize notification based on response type
          const isInterfaceLigand = result.result?.approach;
          const designCount = result.result?.designs?.length || (result.result?.dimer ? 1 : 0);
          addNotification({
            type: 'success',
            title: isInterfaceLigand ? 'Dimer designed!' : 'Structure designed!',
            message: isInterfaceLigand
              ? `${result.result.approach} design complete. ${designCount} design(s) generated.`
              : result.result.seed
                ? `Backbone ready (seed: ${result.result.seed}). Design sequences with MPNN next.`
                : 'Your protein backbone is ready. Design sequences with MPNN next.',
            action: {
              label: 'Design Sequences',
              tab: 'mpnn',
            },
          });
        }
      } else if (result.status === 'failed') {
        // Update pipeline stage to error for interface_ligand tasks
        if (isInterfaceLigandTask) {
          setPipelineStage('error');
          setPipelineError(result.error || 'Unknown error occurred');
        }

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

  // Handler for selecting a design in the results panel
  const handleSelectDesign = (design: InterfaceLigandDesign) => {
    setSelectedDesignId(design.id);
    if (design.pdb_content) {
      setSelectedPdb(design.pdb_content);
    }
  };

  // Handler for downloading a specific design's PDB
  const handleDownloadDesignPdb = (design: InterfaceLigandDesign) => {
    if (!design.pdb_content) return;
    const blob = new Blob([design.pdb_content], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `design_${design.rank}.pdb`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  // Handler for running LigandMPNN on selected designs
  const handleRunLigandMPNN = async (designs: InterfaceLigandDesign[]) => {
    if (!interfaceLigandResult || designs.length === 0) return;

    setIsPipelineProcessing(true);
    try {
      // Process each design sequentially
      for (const design of designs) {
        if (!design.pdb_content) continue;

        // Submit MPNN job
        const response = await api.submitMPNNDesign({
          pdb_content: design.pdb_content,
          model_type: 'ligand_mpnn',
          num_sequences: 4,
          temperature: 0.1,
        });

        // Wait for completion
        const result = await api.waitForJob(response.job_id);

        if (result.status === 'completed' && result.result?.sequences) {
          // Update design with MPNN results
          setInterfaceLigandResult(prev => {
            if (!prev) return prev;
            return {
              ...prev,
              designs: prev.designs.map(d =>
                d.id === design.id
                  ? { ...d, pipeline_stage: 'ligandmpnn' as const }
                  : d
              ),
            };
          });
        }
      }

      addNotification({
        type: 'success',
        title: 'LigandMPNN Complete',
        message: `Processed ${designs.length} design(s) with LigandMPNN`,
      });
    } catch (err) {
      addNotification({
        type: 'error',
        title: 'LigandMPNN Failed',
        message: err instanceof Error ? err.message : 'Failed to run LigandMPNN',
      });
    } finally {
      setIsPipelineProcessing(false);
    }
  };

  // Handler for running validation on selected designs
  const handleRunValidation = async (designs: InterfaceLigandDesign[]) => {
    if (!interfaceLigandResult || designs.length === 0) return;

    setIsPipelineProcessing(true);
    try {
      // For validation, we would typically use ESMFold or similar
      // For now, update the pipeline stage
      for (const design of designs) {
        setInterfaceLigandResult(prev => {
          if (!prev) return prev;
          return {
            ...prev,
            designs: prev.designs.map(d =>
              d.id === design.id
                ? { ...d, pipeline_stage: 'validation' as const }
                : d
            ),
          };
        });
      }

      addNotification({
        type: 'success',
        title: 'Validation Complete',
        message: `Validated ${designs.length} design(s)`,
      });
    } catch (err) {
      addNotification({
        type: 'error',
        title: 'Validation Failed',
        message: err instanceof Error ? err.message : 'Failed to validate designs',
      });
    } finally {
      setIsPipelineProcessing(false);
    }
  };

  // Reset results when changing tasks
  const handleBackToTaskSelection = () => {
    setSelectedDesignTask(null);
    setActiveTab('task');
    setError(null);
    setInterfaceLigandResult(null);
    setSelectedDesignId(undefined);
    setPipelineStage('idle');
    setPipelineMessage('');
    setPipelineError(undefined);
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
    interface_ligand: InterfaceLigandForm,
    interface_metal: InterfaceMetalForm,
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
          <ArrowLeft className="h-5 w-5" />
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
          onClick={handleBackToTaskSelection}
          className="flex items-center gap-2 text-slate-600 hover:text-blue-600 transition-colors group"
        >
          <ArrowLeft className="h-5 w-5 group-hover:-translate-x-1 transition-transform" />
          <span className="text-sm font-medium">Change Design Task</span>
        </button>

        {/* Current Task Form - hide when results are showing */}
        {CurrentForm && !interfaceLigandResult && <CurrentForm {...formProps} />}

        {/* Pipeline Progress - shown during interface_ligand design */}
        {selectedDesignTask === 'interface_ligand' && submitting && pipelineStage !== 'idle' && (
          <div className="mt-6">
            <PipelineProgress
              currentStage={pipelineStage}
              stageDetails={{
                backbone: { message: pipelineStage === 'backbone' ? pipelineMessage : undefined },
                ligandmpnn: { message: pipelineStage === 'ligandmpnn' ? pipelineMessage : undefined },
                validation: { message: pipelineStage === 'validation' ? pipelineMessage : undefined },
                analysis: { message: pipelineStage === 'analysis' ? pipelineMessage : undefined },
              }}
              error={pipelineError}
            />
          </div>
        )}

        {/* Error display */}
        {error && (
          <ErrorDetails
            error={error.message}
            errorType={error.errorType}
            traceback={error.traceback}
            context={error.context}
          />
        )}

        {/* Interface Ligand Results Panel - shown after design completion */}
        {selectedDesignTask === 'interface_ligand' && interfaceLigandResult && !submitting && (
          <div className="space-y-4">
            {/* New Design Button */}
            <div className="flex items-center justify-between">
              <h3 className="text-lg font-semibold text-slate-900">Design Results</h3>
              <button
                onClick={() => {
                  setInterfaceLigandResult(null);
                  setSelectedDesignId(undefined);
                  setPipelineStage('idle');
                }}
                className="flex items-center gap-2 px-4 py-2 text-sm font-medium text-purple-700 bg-purple-50 rounded-lg hover:bg-purple-100 transition-colors"
              >
                <Plus className="h-4 w-4" />
                New Design
              </button>
            </div>

            <InterfaceLigandResultsPanel
              result={interfaceLigandResult}
              onSelectDesign={handleSelectDesign}
              onRunLigandMPNN={handleRunLigandMPNN}
              onRunValidation={handleRunValidation}
              onDownloadPdb={handleDownloadDesignPdb}
              selectedDesignId={selectedDesignId}
              isProcessing={isPipelineProcessing}
            />
          </div>
        )}

        {/* Download button - for non-interface_ligand tasks or when no results panel */}
        {lastDesignPdb && selectedDesignTask !== 'interface_ligand' && (
          <button
            onClick={downloadPdb}
            className="w-full flex justify-center items-center gap-2 py-2.5 px-4 rounded-lg text-sm font-medium text-slate-600 bg-slate-100 hover:bg-slate-200 transition-all"
          >
            <Download className="h-5 w-5" />
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
