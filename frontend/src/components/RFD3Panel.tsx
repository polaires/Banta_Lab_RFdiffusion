'use client';

import { useState, useRef } from 'react';
import { useStore } from '@/lib/store';
import api from '@/lib/api';
import { ContigBuilder } from './ContigBuilder';

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
  const {
    health,
    addJob,
    updateJob,
    setSelectedPdb,
    addNotification,
    setLatestDesignPdb,
    setLastCompletedJobType,
    setLatestRfd3Design
  } = useStore();
  const [contig, setContig] = useState('100');
  const [numDesigns, setNumDesigns] = useState(1);
  const [submitting, setSubmitting] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [showBuilder, setShowBuilder] = useState(false);
  const [seed, setSeed] = useState<number | null>(null);
  const [useSeed, setUseSeed] = useState(false);
  const [inputPdb, setInputPdb] = useState<string | null>(null);
  const [inputPdbName, setInputPdbName] = useState<string | null>(null);
  const [lastDesignPdb, setLastDesignPdb] = useState<string | null>(null);
  const [isDragOver, setIsDragOver] = useState(false);
  const fileInputRef = useRef<HTMLInputElement>(null);

  const handleFileUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (file) {
      const reader = new FileReader();
      reader.onload = (e) => {
        setInputPdb(e.target?.result as string);
        setInputPdbName(file.name);
      };
      reader.readAsText(file);
    }
  };

  const handleDrop = (e: React.DragEvent) => {
    e.preventDefault();
    setIsDragOver(false);
    const file = e.dataTransfer.files[0];
    if (file && file.name.endsWith('.pdb')) {
      const reader = new FileReader();
      reader.onload = (ev) => {
        setInputPdb(ev.target?.result as string);
        setInputPdbName(file.name);
      };
      reader.readAsText(file);
    }
  };

  const generateRandomSeed = () => {
    setSeed(Math.floor(Math.random() * 1000000));
    setUseSeed(true);
  };

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
        ...(useSeed && seed !== null && { seed }),
        ...(inputPdb && { pdb_content: inputPdb }),
      });

      addJob({
        id: response.job_id,
        type: 'rfd3',
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
        addNotification({
          type: 'error',
          title: 'Design failed',
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

  return (
    <div className="p-8 space-y-8">
      {/* Header */}
      <div className="border-b border-slate-100 pb-6">
        <div className="flex items-center gap-3 mb-2">
          <span className="bg-blue-50 text-blue-700 border border-blue-100 text-[11px] font-bold px-2.5 py-0.5 rounded-full uppercase tracking-wide">
            Step 1
          </span>
          <h2 className="text-xl font-bold text-slate-900">RFdiffusion3 - Structure Design</h2>
        </div>
        <p className="text-slate-500 text-sm leading-relaxed max-w-3xl pl-1">
          Generate de novo protein backbone structures. After design, proceed to MPNN for sequence design.
        </p>
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-12 gap-10">
        {/* Left column - Main controls */}
        <div className="lg:col-span-7 space-y-8">
          {/* Quick Examples */}
          <div className="space-y-3">
            <label className="text-xs font-semibold text-slate-400 uppercase tracking-wider pl-1">
              Quick Examples
            </label>
            <div className="flex flex-wrap gap-2">
              {Object.entries(EXAMPLE_CONFIGS).map(([name, config]) => (
                <button
                  key={name}
                  onClick={() => setContig(config.contig)}
                  className="group bg-white border border-slate-200 hover:border-blue-300 hover:bg-blue-50/50 text-slate-600 hover:text-blue-700 px-4 py-2 rounded-lg text-xs font-medium transition-all duration-200 shadow-sm"
                  title={config.description}
                >
                  {name}
                </button>
              ))}
            </div>
          </div>

          {/* Contig Specification */}
          <div className="space-y-2">
            <div className="flex justify-between items-end px-1">
              <label className="text-sm font-medium text-slate-700 flex items-center gap-1.5">
                Contig Specification
                <span
                  className="material-symbols-outlined text-slate-300 text-sm cursor-help hover:text-slate-500 transition-colors"
                  title="Define the length or structure of regions to design"
                >
                  info
                </span>
              </label>
              <button
                onClick={() => setShowBuilder(!showBuilder)}
                className="text-blue-600 hover:text-blue-800 text-xs font-semibold flex items-center gap-1 group transition-colors"
              >
                <span className="material-symbols-outlined text-base group-hover:rotate-12 transition-transform">
                  design_services
                </span>
                {showBuilder ? 'Hide Builder' : 'Visual Builder'}
              </button>
            </div>
            <div>
              <input
                type="text"
                value={contig}
                onChange={(e) => setContig(e.target.value)}
                className="block w-full rounded-lg border-slate-200 bg-slate-50 focus:bg-white text-slate-900 shadow-sm focus:border-blue-500 focus:ring-blue-500 sm:text-sm py-2.5 px-3 transition-colors placeholder:text-slate-400"
                placeholder="e.g., 100 or A1-50/0 50-100"
              />
              <p className="text-xs text-slate-400 mt-2 pl-1">
                Format: &quot;100&quot; for 100 residues, &quot;A1-50/0 50-100&quot; for binders.
              </p>
            </div>
          </div>

          {/* Visual Builder */}
          {showBuilder && (
            <ContigBuilder
              onContigChange={(newContig) => setContig(newContig)}
              initialContig={contig}
            />
          )}

          <div className="border-t border-slate-100 my-4" />

          {/* Number of Designs & Seed */}
          <div className="grid grid-cols-1 sm:grid-cols-2 gap-6">
            <div className="space-y-2">
              <label className="block text-sm font-medium text-slate-700 pl-1">
                Number of Designs
              </label>
              <input
                type="number"
                value={numDesigns}
                onChange={(e) => setNumDesigns(Math.max(1, Math.min(10, parseInt(e.target.value) || 1)))}
                min={1}
                max={10}
                className="block w-full rounded-lg border-slate-200 bg-slate-50 focus:bg-white text-slate-900 shadow-sm focus:border-blue-500 focus:ring-blue-500 sm:text-sm py-2.5 px-3"
              />
            </div>

            <div className="space-y-2">
              <div className="flex items-center justify-between pl-1">
                <label className="block text-sm font-medium text-slate-700">Seed</label>
                <div className="flex items-center gap-2">
                  <input
                    type="checkbox"
                    id="use-seed"
                    checked={useSeed}
                    onChange={(e) => setUseSeed(e.target.checked)}
                    className="h-3.5 w-3.5 text-blue-600 focus:ring-blue-500 border-slate-300 rounded cursor-pointer"
                  />
                  <label htmlFor="use-seed" className="text-xs text-slate-500 cursor-pointer select-none">
                    Fix seed
                  </label>
                </div>
              </div>
              <div className="flex rounded-lg shadow-sm">
                <input
                  type="text"
                  value={useSeed && seed !== null ? seed : ''}
                  onChange={(e) => setSeed(parseInt(e.target.value) || null)}
                  disabled={!useSeed}
                  placeholder="Random"
                  className="flex-1 min-w-0 block w-full px-3 py-2.5 rounded-l-lg border border-slate-200 bg-slate-100 text-slate-500 sm:text-sm focus:ring-blue-500 focus:border-blue-500 disabled:opacity-70 disabled:cursor-not-allowed"
                />
                <button
                  type="button"
                  onClick={generateRandomSeed}
                  className="inline-flex items-center px-3 py-2 border border-l-0 border-slate-200 rounded-r-lg bg-white text-slate-500 hover:text-blue-600 hover:bg-slate-50 focus:outline-none focus:ring-1 focus:ring-blue-500 transition-colors"
                  title="Generate random seed"
                >
                  <span className="material-symbols-outlined text-lg">shuffle</span>
                </button>
              </div>
            </div>
          </div>
        </div>

        {/* Right column - File upload */}
        <div className="lg:col-span-5 flex flex-col">
          <div className="flex-grow space-y-2 flex flex-col">
            <label className="block text-sm font-medium text-slate-700 pl-1">
              Input Structure <span className="text-slate-400 font-normal">(Optional)</span>
            </label>
            <div
              onDragOver={(e) => { e.preventDefault(); setIsDragOver(true); }}
              onDragLeave={() => setIsDragOver(false)}
              onDrop={handleDrop}
              onClick={() => fileInputRef.current?.click()}
              className={`flex-grow min-h-[200px] border-2 border-dashed rounded-xl transition-all cursor-pointer group flex flex-col items-center justify-center p-6 text-center relative ${
                isDragOver
                  ? 'border-blue-400 bg-blue-50/50'
                  : inputPdb
                  ? 'border-emerald-300 bg-emerald-50/30'
                  : 'border-slate-200 hover:border-blue-400 hover:bg-blue-50/30 bg-slate-50/50'
              }`}
            >
              <input
                ref={fileInputRef}
                type="file"
                accept=".pdb"
                onChange={handleFileUpload}
                className="absolute inset-0 w-full h-full opacity-0 cursor-pointer z-10"
              />

              {inputPdb ? (
                <>
                  <div className="w-14 h-14 bg-emerald-100 shadow-sm rounded-full flex items-center justify-center mb-4">
                    <span className="material-symbols-outlined text-emerald-600 text-2xl">check_circle</span>
                  </div>
                  <div className="text-sm text-emerald-700 font-medium">{inputPdbName}</div>
                  <button
                    onClick={(e) => {
                      e.stopPropagation();
                      setInputPdb(null);
                      setInputPdbName(null);
                    }}
                    className="mt-2 text-xs text-slate-500 hover:text-red-600 transition-colors"
                  >
                    Remove file
                  </button>
                </>
              ) : (
                <>
                  <div className="w-14 h-14 bg-white shadow-sm rounded-full flex items-center justify-center mb-4 group-hover:scale-110 transition-transform duration-300">
                    <span className="material-symbols-outlined text-slate-400 group-hover:text-blue-600 text-2xl">
                      upload_file
                    </span>
                  </div>
                  <div className="text-sm text-slate-700 font-medium">
                    Drop PDB file here or <span className="text-blue-600 hover:underline">browse</span>
                  </div>
                  <p className="text-xs text-slate-400 mt-2 max-w-[200px]">
                    Required for binder or scaffold design workflows.
                  </p>
                </>
              )}
            </div>
          </div>
        </div>
      </div>

      {/* Error display */}
      {error && (
        <div className="p-4 bg-red-50 border border-red-200 rounded-lg text-sm text-red-700 flex items-center gap-2">
          <span className="material-symbols-outlined text-red-500">error</span>
          {error}
        </div>
      )}

      {/* Submit button */}
      <div className="pt-6 mt-2 border-t border-slate-100">
        <button
          onClick={handleSubmit}
          disabled={!health || submitting || !contig}
          className="w-full flex justify-center items-center gap-2 py-3.5 px-4 rounded-xl shadow-lg shadow-blue-500/20 text-sm font-semibold text-white bg-blue-600 hover:bg-blue-700 hover:shadow-blue-600/30 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-blue-500 transition-all transform active:scale-[0.99] disabled:opacity-50 disabled:cursor-not-allowed disabled:shadow-none"
        >
          {submitting ? (
            <>
              <span className="material-symbols-outlined text-xl animate-spin">progress_activity</span>
              Running RFdiffusion3...
            </>
          ) : (
            <>
              <span className="material-symbols-outlined text-xl">play_circle</span>
              Design Structure
            </>
          )}
        </button>

        {!health && (
          <p className="text-sm text-amber-600 text-center mt-3">
            Connect to backend to enable design
          </p>
        )}
      </div>

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
