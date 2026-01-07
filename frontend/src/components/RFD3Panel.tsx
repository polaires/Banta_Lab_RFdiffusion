'use client';

import { useState } from 'react';
import { useStore } from '@/lib/store';
import api from '@/lib/api';
import { Play, Loader2, Info, Wand2, Upload, Shuffle, Download } from 'lucide-react';
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
        const pdbContent = result.result.designs[0].content;
        const cifContent = result.result.designs[0].cif_content;
        setSelectedPdb(pdbContent);
        setLatestDesignPdb(pdbContent);
        setLastCompletedJobType('rfd3');
        setLastDesignPdb(pdbContent);

        // Extract sequence from PDB for cross-panel use
        const sequence = extractSequenceFromPdb(pdbContent);

        // Store for cross-panel data flow
        setLatestRfd3Design({
          jobId: response.job_id,
          pdbContent,
          cifContent,
          source: 'rfd3',
          sequence,
          timestamp: Date.now(),
        });

        // Notify user with next step suggestion
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
        <div className="flex items-center justify-between">
          <label className="text-sm font-medium text-gray-300 flex items-center gap-2">
            Contig Specification
            <span title="Define regions to design">
              <Info className="w-4 h-4 text-gray-500" />
            </span>
          </label>
          <button
            onClick={() => setShowBuilder(!showBuilder)}
            className={`text-xs px-2 py-1 rounded flex items-center gap-1 transition ${
              showBuilder
                ? 'bg-blue-600 text-white'
                : 'bg-gray-700 hover:bg-gray-600 text-gray-300'
            }`}
          >
            <Wand2 className="w-3 h-3" />
            {showBuilder ? 'Hide Builder' : 'Visual Builder'}
          </button>
        </div>
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

      {/* Visual Contig Builder */}
      {showBuilder && (
        <ContigBuilder
          onContigChange={(newContig) => setContig(newContig)}
          initialContig={contig}
        />
      )}

      {/* Input PDB for Conditional Design */}
      <div className="space-y-2">
        <label className="text-sm font-medium text-gray-300">
          Input Structure (Optional - for binder/scaffold design)
        </label>
        <label className="flex py-4 border-2 border-dashed border-gray-600 rounded-lg hover:border-gray-500 cursor-pointer transition items-center justify-center gap-2">
          <Upload className="w-5 h-5 text-gray-500" />
          <span className="text-sm text-gray-400">
            {inputPdbName || 'Upload target PDB'}
          </span>
          <input
            type="file"
            accept=".pdb"
            onChange={handleFileUpload}
            className="hidden"
          />
        </label>
        {inputPdb && (
          <div className="flex items-center justify-between text-xs text-gray-500">
            <span>✓ {inputPdbName} loaded</span>
            <button
              onClick={() => {
                setInputPdb(null);
                setInputPdbName(null);
              }}
              className="text-red-400 hover:text-red-300"
            >
              Remove
            </button>
          </div>
        )}
      </div>

      {/* Number of Designs & Seed */}
      <div className="grid grid-cols-2 gap-4">
        <div className="space-y-2">
          <label className="text-sm font-medium text-gray-300">Number of Designs</label>
          <input
            type="number"
            value={numDesigns}
            onChange={(e) => setNumDesigns(Math.max(1, Math.min(10, parseInt(e.target.value) || 1)))}
            min={1}
            max={10}
            className="w-full px-4 py-2 bg-gray-700 rounded border border-gray-600 focus:border-blue-500 focus:outline-none"
          />
        </div>

        <div className="space-y-2">
          <div className="flex items-center gap-2">
            <input
              type="checkbox"
              id="useSeed"
              checked={useSeed}
              onChange={(e) => setUseSeed(e.target.checked)}
              className="w-4 h-4 rounded border-gray-600 bg-gray-700 text-blue-600 focus:ring-blue-500"
            />
            <label htmlFor="useSeed" className="text-sm font-medium text-gray-300">
              Use Seed (reproducibility)
            </label>
          </div>
          <div className="flex gap-2">
            <input
              type="number"
              value={seed ?? ''}
              onChange={(e) => setSeed(parseInt(e.target.value) || null)}
              placeholder="Random seed"
              disabled={!useSeed}
              className="flex-1 px-4 py-2 bg-gray-700 rounded border border-gray-600 focus:border-blue-500 focus:outline-none disabled:opacity-50"
            />
            <button
              onClick={generateRandomSeed}
              disabled={!useSeed}
              className="p-2 bg-gray-600 hover:bg-gray-500 rounded disabled:opacity-50"
              title="Generate random seed"
            >
              <Shuffle className="w-4 h-4" />
            </button>
          </div>
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

      {/* Download Design Button */}
      {lastDesignPdb && (
        <button
          onClick={downloadPdb}
          className="w-full py-2 bg-gray-700 hover:bg-gray-600 rounded font-medium flex items-center justify-center gap-2 transition text-sm"
        >
          <Download className="w-4 h-4" />
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
