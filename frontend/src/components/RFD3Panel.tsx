'use client';

import { useState, useRef } from 'react';
import { useStore } from '@/lib/store';
import api from '@/lib/api';
import { ContigBuilder } from './ContigBuilder';
import { ErrorDetails } from './ErrorDetails';
import { saveJob as saveJobToSupabase, updateJob as updateJobInSupabase } from '@/lib/supabase';

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

// Diffusion quality presets
const DIFFUSION_PRESETS = {
  'Quick': {
    num_timesteps: 50,
    step_scale: 1.8,
    gamma_0: 0.5,
    description: 'Fast generation, good for exploration',
  },
  'Balanced': {
    num_timesteps: 200,
    step_scale: 1.5,
    gamma_0: 0.6,
    description: 'Default settings, balanced quality',
  },
  'High Quality': {
    num_timesteps: 500,
    step_scale: 1.2,
    gamma_0: 0.7,
    description: 'Slower but higher quality designs',
  },
  'Diverse': {
    num_timesteps: 200,
    step_scale: 1.0,
    gamma_0: 0.8,
    description: 'More diverse, creative designs',
  },
};

// Hotspot atom selection options
const HOTSPOT_ATOM_OPTIONS = [
  { value: 'ALL', label: 'All atoms', description: 'Contact any atom in residue' },
  { value: 'BKBN', label: 'Backbone', description: 'Contact backbone (N, CA, C, O)' },
  { value: 'TIP', label: 'Tip atom', description: 'Contact functional tip only' },
];

// Phase 2: Symmetry options
const SYMMETRY_OPTIONS = [
  { id: 'C2', label: 'C2 (Dimer)', subunits: 2, description: '2-fold rotational symmetry' },
  { id: 'C3', label: 'C3 (Trimer)', subunits: 3, description: '3-fold rotational symmetry' },
  { id: 'C4', label: 'C4 (Tetramer)', subunits: 4, description: '4-fold rotational symmetry' },
  { id: 'C5', label: 'C5 (Pentamer)', subunits: 5, description: '5-fold rotational symmetry' },
  { id: 'C6', label: 'C6 (Hexamer)', subunits: 6, description: '6-fold rotational symmetry' },
  { id: 'D2', label: 'D2 (4 chains)', subunits: 4, description: 'Dihedral with 2-fold axes' },
  { id: 'D3', label: 'D3 (6 chains)', subunits: 6, description: 'Dihedral with 3-fold axis' },
  { id: 'D4', label: 'D4 (8 chains)', subunits: 8, description: 'Dihedral with 4-fold axis' },
  { id: 'T', label: 'Tetrahedral', subunits: 12, description: '12 subunits, cage-like' },
  { id: 'O', label: 'Octahedral', subunits: 24, description: '24 subunits, cubic' },
  { id: 'I', label: 'Icosahedral', subunits: 60, description: '60 subunits, viral capsid' },
];

// Phase 2: Common ligands
const COMMON_LIGANDS = [
  { id: 'ATP', label: 'ATP', category: 'Nucleotide', description: 'Adenosine triphosphate' },
  { id: 'ADP', label: 'ADP', category: 'Nucleotide', description: 'Adenosine diphosphate' },
  { id: 'NAD', label: 'NAD+', category: 'Cofactor', description: 'Nicotinamide adenine dinucleotide' },
  { id: 'FAD', label: 'FAD', category: 'Cofactor', description: 'Flavin adenine dinucleotide' },
  { id: 'HEM', label: 'Heme', category: 'Cofactor', description: 'Iron porphyrin' },
  { id: 'ZN', label: 'Zinc', category: 'Metal', description: 'Zinc ion' },
  { id: 'MG', label: 'Magnesium', category: 'Metal', description: 'Magnesium ion' },
  { id: 'CA', label: 'Calcium', category: 'Metal', description: 'Calcium ion' },
  { id: 'FE', label: 'Iron', category: 'Metal', description: 'Iron ion' },
  { id: 'SAM', label: 'SAM', category: 'Cofactor', description: 'S-adenosyl methionine' },
];

// Phase 2: CFG feature options
const CFG_FEATURES = [
  { id: 'active_donor', label: 'H-bond Donors', description: 'Favor hydrogen bond donor surfaces' },
  { id: 'active_acceptor', label: 'H-bond Acceptors', description: 'Favor hydrogen bond acceptor surfaces' },
  { id: 'ref_atomwise_rasa', label: 'Surface Accessibility', description: 'Guide based on RASA patterns' },
];

interface Hotspot {
  id: string;
  chain: string;
  residues: string;
  atoms: string;
}

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
  const [error, setError] = useState<{
    message: string;
    errorType?: string;
    traceback?: string;
    context?: Record<string, any>;
  } | null>(null);
  const [showBuilder, setShowBuilder] = useState(false);
  const [seed, setSeed] = useState<number | null>(null);
  const [useSeed, setUseSeed] = useState(false);
  const [inputPdb, setInputPdb] = useState<string | null>(null);
  const [inputPdbName, setInputPdbName] = useState<string | null>(null);
  const [lastDesignPdb, setLastDesignPdb] = useState<string | null>(null);
  const [isDragOver, setIsDragOver] = useState(false);
  const fileInputRef = useRef<HTMLInputElement>(null);

  // Phase 1: Advanced options state
  const [showAdvanced, setShowAdvanced] = useState(false);
  const [selectedPreset, setSelectedPreset] = useState<string>('Balanced');
  const [numTimesteps, setNumTimesteps] = useState<number | null>(null);
  const [stepScale, setStepScale] = useState<number | null>(null);
  const [gamma0, setGamma0] = useState<number | null>(null);

  // Phase 1: Partial diffusion (only when PDB uploaded)
  const [usePartialDiffusion, setUsePartialDiffusion] = useState(false);
  const [partialT, setPartialT] = useState(10);

  // Phase 1: Hotspot residues
  const [hotspots, setHotspots] = useState<Hotspot[]>([]);
  const [newHotspotChain, setNewHotspotChain] = useState('A');
  const [newHotspotResidues, setNewHotspotResidues] = useState('');
  const [newHotspotAtoms, setNewHotspotAtoms] = useState('ALL');

  // Phase 2: Symmetry design
  const [useSymmetry, setUseSymmetry] = useState(false);
  const [symmetryType, setSymmetryType] = useState('C2');
  const [asymmetricChains, setAsymmetricChains] = useState('');

  // Phase 2: Ligand binding
  const [useLigand, setUseLigand] = useState(false);
  const [selectedLigand, setSelectedLigand] = useState('');
  const [customLigand, setCustomLigand] = useState('');

  // Phase 2: Classifier-Free Guidance
  const [useCFG, setUseCFG] = useState(false);
  const [cfgScale, setCfgScale] = useState(1.5);
  const [cfgFeatures, setCfgFeatures] = useState<string[]>([]);

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

  // Apply diffusion preset
  const applyPreset = (presetName: string) => {
    const preset = DIFFUSION_PRESETS[presetName as keyof typeof DIFFUSION_PRESETS];
    if (preset) {
      setSelectedPreset(presetName);
      setNumTimesteps(preset.num_timesteps);
      setStepScale(preset.step_scale);
      setGamma0(preset.gamma_0);
    }
  };

  // Hotspot management
  const addHotspot = () => {
    if (!newHotspotResidues.trim()) return;
    const newHotspot: Hotspot = {
      id: `${Date.now()}`,
      chain: newHotspotChain,
      residues: newHotspotResidues.trim(),
      atoms: newHotspotAtoms,
    };
    setHotspots([...hotspots, newHotspot]);
    setNewHotspotResidues('');
  };

  const removeHotspot = (id: string) => {
    setHotspots(hotspots.filter((h: Hotspot) => h.id !== id));
  };

  // Convert hotspots to API format
  const buildHotspotsDict = (): Record<string, string> | undefined => {
    if (hotspots.length === 0) return undefined;
    const dict: Record<string, string> = {};
    for (const h of hotspots) {
      const key = `${h.chain}${h.residues}`;
      dict[key] = h.atoms;
    }
    return dict;
  };

  // Phase 2: Toggle CFG feature
  const toggleCfgFeature = (featureId: string) => {
    setCfgFeatures(prev =>
      prev.includes(featureId)
        ? prev.filter(f => f !== featureId)
        : [...prev, featureId]
    );
  };

  // Phase 2: Get effective ligand (custom or selected)
  const getEffectiveLigand = (): string | undefined => {
    if (!useLigand) return undefined;
    if (customLigand.trim()) return customLigand.trim().toUpperCase();
    if (selectedLigand) return selectedLigand;
    return undefined;
  };

  // Count active advanced features for badge
  const countActiveFeatures = (): number => {
    let count = hotspots.length;
    if (usePartialDiffusion) count++;
    if (useSymmetry) count++;
    if (useLigand && getEffectiveLigand()) count++;
    if (useCFG) count++;
    return count;
  };

  const handleSubmit = async () => {
    if (!health) {
      setError({ message: 'Backend not connected' });
      return;
    }

    setError(null);
    setSubmitting(true);

    try {
      // Build request with Phase 1 advanced options
      const request: Parameters<typeof api.submitRFD3Design>[0] = {
        contig,
        num_designs: numDesigns,
        ...(useSeed && seed !== null && { seed }),
        ...(inputPdb && { pdb_content: inputPdb }),
      };

      // Add hotspots if defined
      const hotspotsDict = buildHotspotsDict();
      if (hotspotsDict) {
        request.select_hotspots = hotspotsDict;
      }

      // Add partial diffusion if enabled and PDB is provided
      if (usePartialDiffusion && inputPdb && partialT > 0) {
        request.partial_t = partialT;
      }

      // Add diffusion parameters if not using defaults
      if (numTimesteps !== null && numTimesteps !== 200) {
        request.num_timesteps = numTimesteps;
      }
      if (stepScale !== null && stepScale !== 1.5) {
        request.step_scale = stepScale;
      }
      if (gamma0 !== null && gamma0 !== 0.6) {
        request.gamma_0 = gamma0;
      }

      // Phase 2: Add symmetry configuration
      if (useSymmetry) {
        request.symmetry = {
          id: symmetryType,
          is_symmetric_motif: true,
        };
        if (asymmetricChains.trim()) {
          request.symmetry.is_unsym_motif = asymmetricChains.trim();
        }
      }

      // Phase 2: Add ligand binding
      const effectiveLigand = getEffectiveLigand();
      if (effectiveLigand) {
        request.ligand = effectiveLigand;
      }

      // Phase 2: Add Classifier-Free Guidance
      if (useCFG) {
        request.cfg = {
          enabled: true,
          scale: cfgScale,
          features: cfgFeatures.length > 0 ? cfgFeatures : undefined,
        };
      }

      const response = await api.submitRFD3Design(request);

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

          {/* Diffusion Quality Presets */}
          <div className="space-y-3">
            <label className="text-xs font-semibold text-slate-400 uppercase tracking-wider pl-1 flex items-center gap-1.5">
              Quality Preset
              <span
                className="material-symbols-outlined text-slate-300 text-sm cursor-help hover:text-slate-500 transition-colors"
                title="Control speed vs quality tradeoff"
              >
                info
              </span>
            </label>
            <div className="flex flex-wrap gap-2">
              {Object.entries(DIFFUSION_PRESETS).map(([name, preset]) => (
                <button
                  key={name}
                  onClick={() => applyPreset(name)}
                  className={`group px-4 py-2 rounded-lg text-xs font-medium transition-all duration-200 shadow-sm ${
                    selectedPreset === name
                      ? 'bg-blue-600 text-white border border-blue-600 shadow-blue-200'
                      : 'bg-white border border-slate-200 hover:border-blue-300 hover:bg-blue-50/50 text-slate-600 hover:text-blue-700'
                  }`}
                  title={preset.description}
                >
                  {name}
                </button>
              ))}
            </div>
            <p className="text-xs text-slate-400 pl-1">
              {DIFFUSION_PRESETS[selectedPreset as keyof typeof DIFFUSION_PRESETS]?.description || 'Select a preset'}
            </p>
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

          {/* Advanced Options Toggle */}
          <div className="border-t border-slate-100 pt-4">
            <button
              onClick={() => setShowAdvanced(!showAdvanced)}
              className="flex items-center gap-2 text-sm font-medium text-slate-600 hover:text-blue-600 transition-colors"
            >
              <span className={`material-symbols-outlined text-lg transition-transform ${showAdvanced ? 'rotate-90' : ''}`}>
                chevron_right
              </span>
              Advanced Options
              {countActiveFeatures() > 0 && (
                <span className="bg-blue-100 text-blue-700 text-[10px] px-1.5 py-0.5 rounded-full font-semibold">
                  {countActiveFeatures()}
                </span>
              )}
            </button>
          </div>

          {/* Advanced Options Panel */}
          {showAdvanced && (
            <div className="space-y-6 pl-2 border-l-2 border-slate-100">
              {/* Hotspot Residues */}
              <div className="space-y-3">
                <label className="text-xs font-semibold text-slate-400 uppercase tracking-wider pl-1 flex items-center gap-1.5">
                  Hotspot Residues
                  <span
                    className="material-symbols-outlined text-slate-300 text-sm cursor-help hover:text-slate-500 transition-colors"
                    title="Specify residues the designed protein must contact (~4.5Å)"
                  >
                    info
                  </span>
                </label>

                {/* Hotspot input row */}
                <div className="flex gap-2 items-end">
                  <div className="w-16">
                    <label className="block text-xs text-slate-500 mb-1">Chain</label>
                    <input
                      type="text"
                      value={newHotspotChain}
                      onChange={(e) => setNewHotspotChain(e.target.value.toUpperCase().slice(0, 1))}
                      className="block w-full rounded-lg border-slate-200 bg-slate-50 focus:bg-white text-slate-900 shadow-sm focus:border-blue-500 focus:ring-blue-500 text-sm py-2 px-2.5 text-center"
                      placeholder="A"
                    />
                  </div>
                  <div className="flex-1">
                    <label className="block text-xs text-slate-500 mb-1">Residues</label>
                    <input
                      type="text"
                      value={newHotspotResidues}
                      onChange={(e) => setNewHotspotResidues(e.target.value)}
                      className="block w-full rounded-lg border-slate-200 bg-slate-50 focus:bg-white text-slate-900 shadow-sm focus:border-blue-500 focus:ring-blue-500 text-sm py-2 px-2.5"
                      placeholder="15-20 or 15,16,17"
                    />
                  </div>
                  <div className="w-28">
                    <label className="block text-xs text-slate-500 mb-1">Atoms</label>
                    <select
                      value={newHotspotAtoms}
                      onChange={(e) => setNewHotspotAtoms(e.target.value)}
                      className="block w-full rounded-lg border-slate-200 bg-slate-50 focus:bg-white text-slate-900 shadow-sm focus:border-blue-500 focus:ring-blue-500 text-sm py-2 px-2"
                    >
                      {HOTSPOT_ATOM_OPTIONS.map(opt => (
                        <option key={opt.value} value={opt.value}>{opt.label}</option>
                      ))}
                    </select>
                  </div>
                  <button
                    onClick={addHotspot}
                    disabled={!newHotspotResidues.trim()}
                    className="px-3 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700 disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
                  >
                    <span className="material-symbols-outlined text-lg">add</span>
                  </button>
                </div>

                {/* Hotspot list */}
                {hotspots.length > 0 && (
                  <div className="space-y-2">
                    {hotspots.map(h => (
                      <div key={h.id} className="flex items-center justify-between bg-slate-50 rounded-lg px-3 py-2 text-sm">
                        <div className="flex items-center gap-2">
                          <span className="font-mono text-blue-700">{h.chain}{h.residues}</span>
                          <span className="text-slate-400">·</span>
                          <span className="text-slate-600">{HOTSPOT_ATOM_OPTIONS.find(o => o.value === h.atoms)?.label}</span>
                        </div>
                        <button
                          onClick={() => removeHotspot(h.id)}
                          className="text-slate-400 hover:text-red-500 transition-colors"
                        >
                          <span className="material-symbols-outlined text-lg">close</span>
                        </button>
                      </div>
                    ))}
                  </div>
                )}

                <p className="text-xs text-slate-400 pl-1">
                  Selected residues will be within ~4.5Å of the designed protein for functional binding.
                </p>
              </div>

              {/* Partial Diffusion - only show when PDB is uploaded */}
              {inputPdb && (
                <div className="space-y-3">
                  <div className="flex items-center justify-between">
                    <label className="text-xs font-semibold text-slate-400 uppercase tracking-wider pl-1 flex items-center gap-1.5">
                      Partial Diffusion
                      <span
                        className="material-symbols-outlined text-slate-300 text-sm cursor-help hover:text-slate-500 transition-colors"
                        title="Refine existing structure instead of generating from scratch"
                      >
                        info
                      </span>
                    </label>
                    <div className="flex items-center gap-2">
                      <input
                        type="checkbox"
                        id="use-partial"
                        checked={usePartialDiffusion}
                        onChange={(e) => setUsePartialDiffusion(e.target.checked)}
                        className="h-3.5 w-3.5 text-blue-600 focus:ring-blue-500 border-slate-300 rounded cursor-pointer"
                      />
                      <label htmlFor="use-partial" className="text-xs text-slate-500 cursor-pointer select-none">
                        Enable
                      </label>
                    </div>
                  </div>

                  {usePartialDiffusion && (
                    <div className="space-y-2">
                      <div className="flex items-center gap-4">
                        <input
                          type="range"
                          min={5}
                          max={20}
                          step={1}
                          value={partialT}
                          onChange={(e) => setPartialT(parseInt(e.target.value))}
                          className="flex-1 h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-blue-600"
                        />
                        <span className="text-sm font-mono text-slate-700 w-12 text-right">{partialT}Å</span>
                      </div>
                      <div className="flex justify-between text-xs text-slate-400 px-1">
                        <span>Minor refinement</span>
                        <span>Major redesign</span>
                      </div>
                      <p className="text-xs text-slate-400 pl-1">
                        Lower values preserve more of the input structure. Higher values allow more changes.
                      </p>
                    </div>
                  )}
                </div>
              )}

              {/* Detailed Parameters */}
              <div className="space-y-3">
                <label className="text-xs font-semibold text-slate-400 uppercase tracking-wider pl-1">
                  Detailed Parameters
                </label>
                <div className="grid grid-cols-3 gap-3">
                  <div>
                    <label className="block text-xs text-slate-500 mb-1">Timesteps</label>
                    <input
                      type="number"
                      value={numTimesteps ?? 200}
                      onChange={(e) => setNumTimesteps(parseInt(e.target.value) || 200)}
                      min={50}
                      max={500}
                      className="block w-full rounded-lg border-slate-200 bg-slate-50 focus:bg-white text-slate-900 shadow-sm focus:border-blue-500 focus:ring-blue-500 text-sm py-2 px-2.5"
                    />
                  </div>
                  <div>
                    <label className="block text-xs text-slate-500 mb-1">Step Scale</label>
                    <input
                      type="number"
                      value={stepScale ?? 1.5}
                      onChange={(e) => setStepScale(parseFloat(e.target.value) || 1.5)}
                      min={0.5}
                      max={3.0}
                      step={0.1}
                      className="block w-full rounded-lg border-slate-200 bg-slate-50 focus:bg-white text-slate-900 shadow-sm focus:border-blue-500 focus:ring-blue-500 text-sm py-2 px-2.5"
                    />
                  </div>
                  <div>
                    <label className="block text-xs text-slate-500 mb-1">Gamma</label>
                    <input
                      type="number"
                      value={gamma0 ?? 0.6}
                      onChange={(e) => setGamma0(parseFloat(e.target.value) || 0.6)}
                      min={0.1}
                      max={1.0}
                      step={0.1}
                      className="block w-full rounded-lg border-slate-200 bg-slate-50 focus:bg-white text-slate-900 shadow-sm focus:border-blue-500 focus:ring-blue-500 text-sm py-2 px-2.5"
                    />
                  </div>
                </div>
                <p className="text-xs text-slate-400 pl-1">
                  Timesteps: more = higher quality. Step Scale: higher = more designable. Gamma: lower = more designable.
                </p>
              </div>

              {/* Phase 2: Symmetry Design */}
              <div className="space-y-3 border-t border-slate-100 pt-4">
                <div className="flex items-center justify-between">
                  <label className="text-xs font-semibold text-slate-400 uppercase tracking-wider pl-1 flex items-center gap-1.5">
                    <span className="material-symbols-outlined text-amber-500 text-sm">hub</span>
                    Symmetry Design
                    <span
                      className="material-symbols-outlined text-slate-300 text-sm cursor-help hover:text-slate-500 transition-colors"
                      title="Design symmetric oligomeric proteins (dimers, trimers, etc.)"
                    >
                      info
                    </span>
                  </label>
                  <div className="flex items-center gap-2">
                    <input
                      type="checkbox"
                      id="use-symmetry"
                      checked={useSymmetry}
                      onChange={(e) => setUseSymmetry(e.target.checked)}
                      className="h-3.5 w-3.5 text-blue-600 focus:ring-blue-500 border-slate-300 rounded cursor-pointer"
                    />
                    <label htmlFor="use-symmetry" className="text-xs text-slate-500 cursor-pointer select-none">
                      Enable
                    </label>
                  </div>
                </div>

                {useSymmetry && (
                  <div className="space-y-3">
                    <div>
                      <label className="block text-xs text-slate-500 mb-2">Symmetry Type</label>
                      <div className="grid grid-cols-4 gap-2">
                        {SYMMETRY_OPTIONS.slice(0, 8).map(sym => (
                          <button
                            key={sym.id}
                            onClick={() => setSymmetryType(sym.id)}
                            className={`px-2 py-1.5 rounded-lg text-xs font-medium transition-all ${
                              symmetryType === sym.id
                                ? 'bg-amber-500 text-white border border-amber-500'
                                : 'bg-white border border-slate-200 text-slate-600 hover:border-amber-300 hover:bg-amber-50'
                            }`}
                            title={sym.description}
                          >
                            {sym.label}
                          </button>
                        ))}
                      </div>
                      <div className="grid grid-cols-3 gap-2 mt-2">
                        {SYMMETRY_OPTIONS.slice(8).map(sym => (
                          <button
                            key={sym.id}
                            onClick={() => setSymmetryType(sym.id)}
                            className={`px-2 py-1.5 rounded-lg text-xs font-medium transition-all ${
                              symmetryType === sym.id
                                ? 'bg-amber-500 text-white border border-amber-500'
                                : 'bg-white border border-slate-200 text-slate-600 hover:border-amber-300 hover:bg-amber-50'
                            }`}
                            title={sym.description}
                          >
                            {sym.label}
                          </button>
                        ))}
                      </div>
                    </div>

                    <div className="bg-amber-50 rounded-lg p-3 text-xs">
                      <div className="flex items-center gap-2 text-amber-700 font-medium mb-1">
                        <span className="material-symbols-outlined text-sm">info</span>
                        {SYMMETRY_OPTIONS.find(s => s.id === symmetryType)?.label}
                      </div>
                      <p className="text-amber-600">
                        {SYMMETRY_OPTIONS.find(s => s.id === symmetryType)?.subunits} identical subunits.
                        {['T', 'O', 'I'].includes(symmetryType) && ' High memory usage - uses low memory mode.'}
                      </p>
                    </div>

                    <div>
                      <label className="block text-xs text-slate-500 mb-1">Asymmetric Chains (optional)</label>
                      <input
                        type="text"
                        value={asymmetricChains}
                        onChange={(e) => setAsymmetricChains(e.target.value)}
                        className="block w-full rounded-lg border-slate-200 bg-slate-50 focus:bg-white text-slate-900 shadow-sm focus:border-blue-500 focus:ring-blue-500 text-sm py-2 px-2.5"
                        placeholder="e.g., Y1-11,Z16-25"
                      />
                      <p className="text-xs text-slate-400 mt-1">
                        Chains/residues that should NOT be symmetrized (e.g., DNA, ligand).
                      </p>
                    </div>
                  </div>
                )}
              </div>

              {/* Phase 2: Ligand Binding */}
              <div className="space-y-3 border-t border-slate-100 pt-4">
                <div className="flex items-center justify-between">
                  <label className="text-xs font-semibold text-slate-400 uppercase tracking-wider pl-1 flex items-center gap-1.5">
                    <span className="material-symbols-outlined text-emerald-500 text-sm">science</span>
                    Ligand Binding
                    <span
                      className="material-symbols-outlined text-slate-300 text-sm cursor-help hover:text-slate-500 transition-colors"
                      title="Design protein with binding site for specific ligand"
                    >
                      info
                    </span>
                  </label>
                  <div className="flex items-center gap-2">
                    <input
                      type="checkbox"
                      id="use-ligand"
                      checked={useLigand}
                      onChange={(e) => setUseLigand(e.target.checked)}
                      className="h-3.5 w-3.5 text-blue-600 focus:ring-blue-500 border-slate-300 rounded cursor-pointer"
                    />
                    <label htmlFor="use-ligand" className="text-xs text-slate-500 cursor-pointer select-none">
                      Enable
                    </label>
                  </div>
                </div>

                {useLigand && (
                  <div className="space-y-3">
                    <div>
                      <label className="block text-xs text-slate-500 mb-2">Common Ligands</label>
                      <div className="flex flex-wrap gap-2">
                        {COMMON_LIGANDS.map(lig => (
                          <button
                            key={lig.id}
                            onClick={() => { setSelectedLigand(lig.id); setCustomLigand(''); }}
                            className={`px-3 py-1.5 rounded-lg text-xs font-medium transition-all ${
                              selectedLigand === lig.id && !customLigand
                                ? 'bg-emerald-500 text-white border border-emerald-500'
                                : 'bg-white border border-slate-200 text-slate-600 hover:border-emerald-300 hover:bg-emerald-50'
                            }`}
                            title={`${lig.description} (${lig.category})`}
                          >
                            {lig.label}
                          </button>
                        ))}
                      </div>
                    </div>

                    <div>
                      <label className="block text-xs text-slate-500 mb-1">Or Custom Ligand ID</label>
                      <input
                        type="text"
                        value={customLigand}
                        onChange={(e) => { setCustomLigand(e.target.value); if (e.target.value) setSelectedLigand(''); }}
                        className="block w-full rounded-lg border-slate-200 bg-slate-50 focus:bg-white text-slate-900 shadow-sm focus:border-blue-500 focus:ring-blue-500 text-sm py-2 px-2.5"
                        placeholder="e.g., GTP, FMN, or custom 3-letter code"
                      />
                    </div>

                    {(selectedLigand || customLigand) && (
                      <div className="bg-emerald-50 rounded-lg p-3 text-xs">
                        <div className="flex items-center gap-2 text-emerald-700 font-medium">
                          <span className="material-symbols-outlined text-sm">check_circle</span>
                          Designing binding site for: <span className="font-mono">{getEffectiveLigand()}</span>
                        </div>
                      </div>
                    )}
                  </div>
                )}
              </div>

              {/* Phase 2: Classifier-Free Guidance */}
              <div className="space-y-3 border-t border-slate-100 pt-4">
                <div className="flex items-center justify-between">
                  <label className="text-xs font-semibold text-slate-400 uppercase tracking-wider pl-1 flex items-center gap-1.5">
                    <span className="material-symbols-outlined text-purple-500 text-sm">tune</span>
                    Classifier-Free Guidance
                    <span
                      className="material-symbols-outlined text-slate-300 text-sm cursor-help hover:text-slate-500 transition-colors"
                      title="Guide designs toward specific surface properties"
                    >
                      info
                    </span>
                  </label>
                  <div className="flex items-center gap-2">
                    <input
                      type="checkbox"
                      id="use-cfg"
                      checked={useCFG}
                      onChange={(e) => setUseCFG(e.target.checked)}
                      className="h-3.5 w-3.5 text-blue-600 focus:ring-blue-500 border-slate-300 rounded cursor-pointer"
                    />
                    <label htmlFor="use-cfg" className="text-xs text-slate-500 cursor-pointer select-none">
                      Enable
                    </label>
                  </div>
                </div>

                {useCFG && (
                  <div className="space-y-3">
                    <div>
                      <label className="block text-xs text-slate-500 mb-1">Guidance Strength</label>
                      <div className="flex items-center gap-4">
                        <input
                          type="range"
                          min={0.5}
                          max={3.0}
                          step={0.1}
                          value={cfgScale}
                          onChange={(e) => setCfgScale(parseFloat(e.target.value))}
                          className="flex-1 h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-purple-600"
                        />
                        <span className="text-sm font-mono text-slate-700 w-10 text-right">{cfgScale.toFixed(1)}</span>
                      </div>
                      <div className="flex justify-between text-xs text-slate-400 px-1 mt-1">
                        <span>Subtle</span>
                        <span>Strong</span>
                      </div>
                    </div>

                    <div>
                      <label className="block text-xs text-slate-500 mb-2">Guide Toward Features</label>
                      <div className="space-y-2">
                        {CFG_FEATURES.map(feature => (
                          <label
                            key={feature.id}
                            className={`flex items-center gap-3 p-2 rounded-lg cursor-pointer transition-colors ${
                              cfgFeatures.includes(feature.id)
                                ? 'bg-purple-50 border border-purple-200'
                                : 'bg-slate-50 border border-transparent hover:bg-slate-100'
                            }`}
                          >
                            <input
                              type="checkbox"
                              checked={cfgFeatures.includes(feature.id)}
                              onChange={() => toggleCfgFeature(feature.id)}
                              className="h-3.5 w-3.5 text-purple-600 focus:ring-purple-500 border-slate-300 rounded"
                            />
                            <div>
                              <span className="text-sm font-medium text-slate-700">{feature.label}</span>
                              <p className="text-xs text-slate-500">{feature.description}</p>
                            </div>
                          </label>
                        ))}
                      </div>
                    </div>

                    <p className="text-xs text-slate-400 pl-1">
                      CFG biases the diffusion process to favor designs with selected surface properties.
                    </p>
                  </div>
                )}
              </div>
            </div>
          )}
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
        <ErrorDetails
          error={error.message}
          errorType={error.errorType}
          traceback={error.traceback}
          context={error.context}
        />
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
