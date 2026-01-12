'use client';

import { useState } from 'react';
import { FormSection, FormField } from './shared/FormSection';
import { QUALITY_PRESETS, RFD3Request, TaskFormProps } from './shared/types';

type QualityPreset = keyof typeof QUALITY_PRESETS;

// Common ligands with SMILES
const COMMON_LIGANDS_SMILES = [
  { id: 'azobenzene', name: 'Azobenzene', smiles: 'c1ccc(cc1)N=Nc2ccccc2', description: 'Light-switchable dye' },
  { id: 'caffeine', name: 'Caffeine', smiles: 'Cn1cnc2c1c(=O)n(c(=O)n2C)C', description: 'Stimulant alkaloid' },
  { id: 'benzene', name: 'Benzene', smiles: 'c1ccccc1', description: 'Simple aromatic ring' },
  { id: 'naphthalene', name: 'Naphthalene', smiles: 'c1ccc2ccccc2c1', description: 'Fused bicyclic aromatic' },
  { id: 'biphenyl', name: 'Biphenyl', smiles: 'c1ccc(cc1)c2ccccc2', description: 'Two connected phenyl rings' },
];

// Design approaches
const APPROACHES = [
  {
    id: 'asymmetric',
    name: 'Asymmetric Binder',
    description: 'Design one-sided binder (Chain A only)',
    details: 'Creates a protein that binds to one side of the ligand, leaving the other side exposed for a second chain.',
  },
  {
    id: 'full',
    name: 'Full Dimer',
    description: 'Design complete separable dimer (A + B)',
    details: 'First designs Chain A, then designs Chain B using A as context. Creates two chains that can physically separate.',
  },
];

export function InterfaceLigandForm({ onSubmit, isSubmitting, health }: TaskFormProps) {
  // Ligand selection
  const [ligandSmiles, setLigandSmiles] = useState('c1ccc(cc1)N=Nc2ccccc2'); // Azobenzene default
  const [selectedPreset, setSelectedPreset] = useState('azobenzene');
  const [customSmiles, setCustomSmiles] = useState(false);

  // Approach
  const [approach, setApproach] = useState<'asymmetric' | 'full'>('full');
  const [side, setSide] = useState<'left' | 'right'>('left');

  // Design parameters
  const [chainLength, setChainLength] = useState('60-80');
  const [numDesigns, setNumDesigns] = useState(3);
  const [seed, setSeed] = useState<string>('');
  const [qualityPreset, setQualityPreset] = useState<QualityPreset>('Balanced');

  // Advanced options
  const [showAdvanced, setShowAdvanced] = useState(false);
  const [useOriToken, setUseOriToken] = useState(false);
  const [oriOffset, setOriOffset] = useState('12.0, 0.0, 0.0');

  const handlePresetSelect = (presetId: string) => {
    setSelectedPreset(presetId);
    const preset = COMMON_LIGANDS_SMILES.find(l => l.id === presetId);
    if (preset) {
      setLigandSmiles(preset.smiles);
      setCustomSmiles(false);
    }
  };

  const handleSubmit = async () => {
    const preset = QUALITY_PRESETS[qualityPreset];

    // Parse ori_offset if provided
    let parsedOriOffset: number[] | undefined;
    if (useOriToken && oriOffset) {
      parsedOriOffset = oriOffset.split(',').map(v => parseFloat(v.trim()));
      if (parsedOriOffset.length !== 3 || parsedOriOffset.some(isNaN)) {
        parsedOriOffset = undefined;
      }
    }

    // Build the request - using task field for interface_ligand_design
    const request: RFD3Request = {
      task: 'interface_ligand_design',
      approach,
      ligand_smiles: ligandSmiles,
      chain_length: chainLength,
      num_designs: numDesigns,
      side,
      num_timesteps: preset.num_timesteps,
      step_scale: preset.step_scale,
      gamma_0: preset.gamma_0,
    };

    if (seed) {
      request.seed = parseInt(seed, 10);
    }

    if (parsedOriOffset) {
      request.ori_offset = parsedOriOffset;
      request.use_ori_token = true;
    }

    await onSubmit(request);
  };

  const isValid = ligandSmiles.trim() !== '' && chainLength.trim() !== '';
  const selectedLigandName = COMMON_LIGANDS_SMILES.find(l => l.id === selectedPreset)?.name || 'Custom';

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex items-center gap-3 pb-4 border-b border-slate-200">
        <div className="w-10 h-10 rounded-lg bg-gradient-to-br from-purple-500 to-pink-500 flex items-center justify-center">
          <span className="material-symbols-outlined text-white">link</span>
        </div>
        <div>
          <h2 className="font-semibold text-slate-900">Interface Ligand Dimer Design</h2>
          <p className="text-sm text-slate-500">Design protein dimers with ligand at the interface (separable topology)</p>
        </div>
      </div>

      {/* Info Banner */}
      <div className="p-4 rounded-xl bg-gradient-to-r from-purple-50 to-pink-50 border border-purple-100">
        <div className="flex items-start gap-3">
          <span className="material-symbols-outlined text-purple-600 mt-0.5">info</span>
          <div className="text-sm text-purple-900">
            <p className="font-medium mb-1">Separable Dimer Design</p>
            <p className="text-purple-700">
              Unlike buried ligand designs, this creates dimers where chains A and B can physically separate.
              The ligand sits at the interface between the two chains.
            </p>
          </div>
        </div>
      </div>

      {/* Ligand Selection */}
      <FormSection
        title="Ligand"
        description="Select a ligand or enter custom SMILES"
        required
      >
        <div className="space-y-4">
          {/* Preset buttons */}
          <div className="grid grid-cols-2 sm:grid-cols-3 gap-2">
            {COMMON_LIGANDS_SMILES.map((lig) => (
              <button
                key={lig.id}
                onClick={() => handlePresetSelect(lig.id)}
                className={`p-3 rounded-lg border text-left transition-all ${
                  selectedPreset === lig.id && !customSmiles
                    ? 'border-purple-400 bg-purple-50 ring-1 ring-purple-200'
                    : 'border-slate-200 hover:border-slate-300 hover:bg-slate-50'
                }`}
              >
                <div className="font-medium text-sm text-slate-800">{lig.name}</div>
                <div className="text-xs text-slate-500 mt-0.5">{lig.description}</div>
              </button>
            ))}
            <button
              onClick={() => setCustomSmiles(true)}
              className={`p-3 rounded-lg border text-left transition-all ${
                customSmiles
                  ? 'border-purple-400 bg-purple-50 ring-1 ring-purple-200'
                  : 'border-slate-200 hover:border-slate-300 hover:bg-slate-50'
              }`}
            >
              <div className="font-medium text-sm text-slate-800">Custom SMILES</div>
              <div className="text-xs text-slate-500 mt-0.5">Enter your own molecule</div>
            </button>
          </div>

          {/* SMILES input */}
          <div>
            <label className="text-xs text-slate-500 font-medium">SMILES String</label>
            <input
              type="text"
              value={ligandSmiles}
              onChange={(e) => {
                setLigandSmiles(e.target.value);
                setCustomSmiles(true);
              }}
              placeholder="e.g., c1ccc(cc1)N=Nc2ccccc2"
              className="mt-1 w-full px-3 py-2 rounded-lg border border-slate-200 focus:border-purple-400 focus:ring-1 focus:ring-purple-200 outline-none font-mono text-sm"
            />
          </div>
        </div>
      </FormSection>

      {/* Design Approach */}
      <FormSection
        title="Design Approach"
        description="Choose how to design the dimer"
        required
      >
        <div className="grid grid-cols-1 sm:grid-cols-2 gap-3">
          {APPROACHES.map((app) => (
            <button
              key={app.id}
              onClick={() => setApproach(app.id as 'asymmetric' | 'full')}
              className={`p-4 rounded-xl border text-left transition-all ${
                approach === app.id
                  ? 'border-purple-400 bg-purple-50 ring-1 ring-purple-200'
                  : 'border-slate-200 hover:border-slate-300 hover:bg-slate-50'
              }`}
            >
              <div className="flex items-center gap-2 mb-1">
                <span className={`w-4 h-4 rounded-full border-2 flex items-center justify-center ${
                  approach === app.id ? 'border-purple-600' : 'border-slate-300'
                }`}>
                  {approach === app.id && <span className="w-2 h-2 rounded-full bg-purple-600" />}
                </span>
                <span className="font-medium text-slate-900">{app.name}</span>
              </div>
              <p className="text-sm text-slate-600 ml-6">{app.description}</p>
              <p className="text-xs text-slate-400 ml-6 mt-1">{app.details}</p>
            </button>
          ))}
        </div>

        {/* Side selection for asymmetric */}
        {approach === 'asymmetric' && (
          <div className="mt-4 p-3 rounded-lg bg-slate-50 border border-slate-200">
            <label className="text-xs text-slate-500 font-medium">Binding Side</label>
            <div className="mt-2 flex gap-3">
              <button
                onClick={() => setSide('left')}
                className={`flex-1 py-2 px-4 rounded-lg border text-sm font-medium transition-all ${
                  side === 'left'
                    ? 'border-purple-400 bg-purple-100 text-purple-700'
                    : 'border-slate-200 text-slate-600 hover:bg-slate-100'
                }`}
              >
                Left side (expose right)
              </button>
              <button
                onClick={() => setSide('right')}
                className={`flex-1 py-2 px-4 rounded-lg border text-sm font-medium transition-all ${
                  side === 'right'
                    ? 'border-purple-400 bg-purple-100 text-purple-700'
                    : 'border-slate-200 text-slate-600 hover:bg-slate-100'
                }`}
              >
                Right side (expose left)
              </button>
            </div>
            <p className="text-xs text-slate-500 mt-2">
              For azobenzene: left = first phenyl ring, right = second phenyl ring
            </p>
          </div>
        )}
      </FormSection>

      {/* Design Parameters */}
      <FormSection title="Design Parameters" required>
        <div className="space-y-4">
          <div className="grid grid-cols-3 gap-4">
            <FormField label="Chain Length" hint="e.g., 60-80">
              <input
                type="text"
                value={chainLength}
                onChange={(e) => setChainLength(e.target.value)}
                placeholder="60-80"
                className="w-full px-3 py-2 rounded-lg border border-slate-200 focus:border-purple-400 focus:ring-1 focus:ring-purple-200 outline-none"
              />
            </FormField>
            <FormField label="# Designs">
              <input
                type="number"
                value={numDesigns}
                onChange={(e) => setNumDesigns(Math.max(1, parseInt(e.target.value) || 1))}
                min={1}
                max={10}
                className="w-full px-3 py-2 rounded-lg border border-slate-200 focus:border-purple-400 focus:ring-1 focus:ring-purple-200 outline-none"
              />
            </FormField>
            <FormField label="Seed" hint="Optional">
              <input
                type="number"
                value={seed}
                onChange={(e) => setSeed(e.target.value)}
                placeholder="Random"
                className="w-full px-3 py-2 rounded-lg border border-slate-200 focus:border-purple-400 focus:ring-1 focus:ring-purple-200 outline-none"
              />
            </FormField>
          </div>

          {/* Quality preset */}
          <div>
            <label className="text-xs text-slate-500 font-medium">Quality</label>
            <div className="mt-1 grid grid-cols-2 sm:grid-cols-4 gap-2">
              {(Object.keys(QUALITY_PRESETS) as QualityPreset[]).map((preset) => (
                <button
                  key={preset}
                  onClick={() => setQualityPreset(preset)}
                  className={`p-2.5 rounded-lg border text-left transition-all ${
                    qualityPreset === preset
                      ? 'border-purple-400 bg-purple-50'
                      : 'border-slate-200 hover:border-slate-300 hover:bg-slate-50'
                  }`}
                >
                  <div className="font-medium text-sm text-slate-800">{preset}</div>
                  <div className="text-xs text-slate-500 mt-0.5 line-clamp-1">
                    {QUALITY_PRESETS[preset].description}
                  </div>
                </button>
              ))}
            </div>
          </div>
        </div>
      </FormSection>

      {/* Advanced Options */}
      <div className="border border-slate-200 rounded-xl overflow-hidden">
        <button
          onClick={() => setShowAdvanced(!showAdvanced)}
          className="w-full px-4 py-3 flex items-center justify-between bg-slate-50 hover:bg-slate-100 transition-colors"
        >
          <span className="text-sm font-medium text-slate-700">Advanced Options</span>
          <span className={`material-symbols-outlined text-slate-400 transition-transform ${showAdvanced ? 'rotate-180' : ''}`}>
            expand_more
          </span>
        </button>

        {showAdvanced && (
          <div className="p-4 space-y-4 border-t border-slate-200">
            <label className="flex items-center gap-3 cursor-pointer">
              <input
                type="checkbox"
                checked={useOriToken}
                onChange={(e) => setUseOriToken(e.target.checked)}
                className="w-4 h-4 rounded border-slate-300 text-purple-600"
              />
              <span className="text-sm text-slate-600">Use ori_token offset (experimental)</span>
            </label>

            {useOriToken && (
              <div>
                <label className="text-xs text-slate-500 font-medium">Offset from ligand (x, y, z)</label>
                <input
                  type="text"
                  value={oriOffset}
                  onChange={(e) => setOriOffset(e.target.value)}
                  placeholder="12.0, 0.0, 0.0"
                  className="mt-1 w-full px-3 py-2 rounded-lg border border-slate-200 text-sm font-mono"
                />
                <p className="text-xs text-slate-500 mt-1">
                  Positions protein origin away from ligand center. Use cautiously.
                </p>
              </div>
            )}
          </div>
        )}
      </div>

      {/* Submit Button */}
      <div className="pt-4 border-t border-slate-200">
        <button
          onClick={handleSubmit}
          disabled={!isValid || isSubmitting || !health}
          className={`w-full py-3.5 px-6 rounded-xl font-semibold text-white transition-all flex items-center justify-center gap-2 ${
            isValid && !isSubmitting && !!health
              ? 'bg-gradient-to-r from-purple-600 to-pink-600 hover:from-purple-700 hover:to-pink-700 shadow-lg shadow-purple-600/20'
              : 'bg-slate-300 cursor-not-allowed'
          }`}
        >
          {isSubmitting ? (
            <>
              <span className="material-symbols-outlined animate-spin">progress_activity</span>
              Designing {approach === 'full' ? 'Dimer' : 'Binder'}...
            </>
          ) : (
            <>
              <span className="material-symbols-outlined">link</span>
              Design {selectedLigandName} {approach === 'full' ? 'Dimer' : 'Binder'}
            </>
          )}
        </button>
        {!health && (
          <p className="text-center text-sm text-slate-500 mt-2">
            Backend service unavailable
          </p>
        )}

        {/* Expected Results */}
        <div className="mt-4 p-3 rounded-lg bg-slate-50 border border-slate-200">
          <p className="text-xs text-slate-500 font-medium mb-2">Expected Results</p>
          <div className="grid grid-cols-2 gap-2 text-xs text-slate-600">
            <div>Affinity: -3 to -5 kcal/mol</div>
            <div>Topology: Separable</div>
            <div>Chain contacts: ≥3 each</div>
            <div>Clashes: None</div>
          </div>
        </div>
      </div>
    </div>
  );
}
