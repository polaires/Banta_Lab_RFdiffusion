'use client';

import { useState } from 'react';
import { FormSection, FormField } from './shared/FormSection';
import { LengthRangeInput } from './shared/LengthRangeInput';
import { QualityPresetSelector, QualityPreset, QualityParams } from './shared/QualityPresetSelector';
import { AdvancedOptionsWrapper } from './shared/AdvancedOptionsWrapper';
import { QUALITY_PRESETS, RFD3Request, TaskFormProps } from './shared/types';

// Common ligands with SMILES
const COMMON_LIGANDS_SMILES = [
  { id: 'azobenzene', name: 'Azobenzene', smiles: 'c1ccc(cc1)N=Nc2ccccc2', description: 'Light-switchable dye' },
  { id: 'caffeine', name: 'Caffeine', smiles: 'Cn1cnc2c1c(=O)n(c(=O)n2C)C', description: 'Stimulant alkaloid' },
  { id: 'benzene', name: 'Benzene', smiles: 'c1ccccc1', description: 'Simple aromatic ring' },
  { id: 'naphthalene', name: 'Naphthalene', smiles: 'c1ccc2ccccc2c1', description: 'Fused bicyclic aromatic' },
  { id: 'biphenyl', name: 'Biphenyl', smiles: 'c1ccc(cc1)c2ccccc2', description: 'Two connected phenyl rings' },
];

// Design approaches - based on tested workflows in docs/
// See: joint_heterodimer_workflow.md and heterodimer_context_protein_workflow.md
const APPROACHES = [
  {
    id: 'joint',
    name: 'Joint Heterodimer',
    description: 'Both chains simultaneously',
    details: 'Multi-chain RFD3 diffusion. Both chains co-evolve around the ligand in same coordinate frame. Best for true heterodimers.',
    icon: 'hub',
    recommended: true,
  },
  {
    id: 'asymmetric_rasa',
    name: 'Asymmetric RASA',
    description: 'RASA-conditioned dimer',
    details: 'Complementary burial/exposure for each chain. Chain A buries one side, Chain B buries the other.',
    icon: 'tune',
    recommended: false,
  },
  {
    id: 'induced',
    name: 'Induced Dimerization',
    description: 'Context protein workflow',
    details: 'Design Chain A first, then Chain B using A + ligand as fixed context. Good for induced dimerization designs.',
    icon: 'account_tree',
    recommended: false,
  },
  {
    id: 'asymmetric',
    name: 'Single Chain',
    description: 'One-sided binder only',
    details: 'Design a protein binding one side of ligand. Select which side to bind (left/right).',
    icon: 'link',
    recommended: false,
  },
];

// Expected Results thresholds matching backend handler.py
const EXPECTED_RESULTS = {
  excellent: {
    affinity: '< -7 kcal/mol',
    contacts: '≥ 10 per chain',
    identity: '< 40%',
    antiHomo: '> 80/100',
    hbonds: 'N5 ≥ 2, N6 ≥ 2',
  },
  good: {
    affinity: '-5 to -7 kcal/mol',
    contacts: '5-10 per chain',
    identity: '< 70%',
    antiHomo: '> 60/100',
    hbonds: 'N5 ≥ 1, N6 ≥ 1',
  },
  minimum: {
    affinity: '< -2 kcal/mol',
    contacts: '≥ 2 per chain',
    identity: '< 90%',
    antiHomo: '> 40/100',
    hbonds: 'At least 1 total',
  },
};

export function InterfaceLigandForm({ onSubmit, isSubmitting, health }: TaskFormProps) {
  // Ligand selection
  const [ligandSmiles, setLigandSmiles] = useState('c1ccc(cc1)N=Nc2ccccc2'); // Azobenzene default
  const [selectedPreset, setSelectedPreset] = useState('azobenzene');
  const [customSmiles, setCustomSmiles] = useState(false);

  // Approach - default to joint (best results per docs)
  const [approach, setApproach] = useState<'joint' | 'asymmetric_rasa' | 'induced' | 'asymmetric'>('joint');
  const [side, setSide] = useState<'left' | 'right'>('left');

  // Design parameters
  const [chainLength, setChainLength] = useState('60-80');
  const [numDesigns, setNumDesigns] = useState(3);
  const [seed, setSeed] = useState<string>('');
  const [qualityPreset, setQualityPreset] = useState<QualityPreset>('Balanced');
  const [qualityParams, setQualityParams] = useState<QualityParams>(QUALITY_PRESETS.Balanced);

  // Advanced options
  const [useOriToken, setUseOriToken] = useState(false);
  const [oriOffset, setOriOffset] = useState('12.0, 0.0, 0.0');

  const handleQualityChange = (preset: QualityPreset, params: QualityParams) => {
    setQualityPreset(preset);
    setQualityParams(params);
  };

  const handlePresetSelect = (presetId: string) => {
    setSelectedPreset(presetId);
    const preset = COMMON_LIGANDS_SMILES.find(l => l.id === presetId);
    if (preset) {
      setLigandSmiles(preset.smiles);
      setCustomSmiles(false);
    }
  };

  const handleSubmit = async () => {
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
      num_timesteps: qualityParams.num_timesteps,
      step_scale: qualityParams.step_scale,
      gamma_0: qualityParams.gamma_0,
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
              className="mt-1 w-full px-3 py-2 rounded-lg border border-slate-200 focus:border-purple-400 focus:ring-2 focus:ring-purple-100 outline-none font-mono text-sm"
            />
          </div>
        </div>
      </FormSection>

      {/* Design Approach */}
      <FormSection
        title="Design Approach"
        description="Choose how to design the heterodimer"
        required
      >
        <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-3">
          {APPROACHES.map((app) => (
            <button
              key={app.id}
              onClick={() => setApproach(app.id as 'joint' | 'asymmetric_rasa' | 'induced' | 'asymmetric')}
              className={`p-4 rounded-xl border text-left transition-all relative ${
                approach === app.id
                  ? 'border-purple-400 bg-purple-50 ring-1 ring-purple-200'
                  : 'border-slate-200 hover:border-slate-300 hover:bg-slate-50'
              }`}
            >
              {app.recommended && (
                <span className="absolute -top-2 -right-2 px-2 py-0.5 bg-emerald-500 text-white text-xs font-medium rounded-full">
                  Recommended
                </span>
              )}
              <div className="flex items-center gap-2 mb-2">
                <span className={`material-symbols-outlined text-lg ${
                  approach === app.id ? 'text-purple-600' : 'text-slate-400'
                }`}>
                  {app.icon}
                </span>
                <span className="font-medium text-slate-900">{app.name}</span>
              </div>
              <p className="text-sm text-slate-600">{app.description}</p>
              <p className="text-xs text-slate-400 mt-1">{app.details}</p>
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
          <div className="flex gap-4 items-start">
            <div className="flex-1">
              <LengthRangeInput
                value={chainLength}
                onChange={setChainLength}
                label="Chain Length"
                placeholder="60-80"
                hint="Per chain length range"
              />
            </div>
            <FormField label="# Designs" className="w-24">
              <input
                type="number"
                value={numDesigns}
                onChange={(e) => setNumDesigns(Math.max(1, parseInt(e.target.value) || 1))}
                min={1}
                max={10}
                className="w-full px-3 py-2 rounded-lg border border-slate-200 focus:border-purple-400 focus:ring-2 focus:ring-purple-100 outline-none text-sm"
              />
            </FormField>
            <FormField label="Seed" hint="Optional" className="w-24">
              <input
                type="number"
                value={seed}
                onChange={(e) => setSeed(e.target.value)}
                placeholder="Random"
                className="w-full px-3 py-2 rounded-lg border border-slate-200 focus:border-purple-400 focus:ring-2 focus:ring-purple-100 outline-none text-sm"
              />
            </FormField>
          </div>

          {/* Quality preset */}
          <QualityPresetSelector
            value={qualityPreset}
            onChange={handleQualityChange}
            showDescription
          />
        </div>
      </FormSection>

      {/* Advanced Options */}
      <AdvancedOptionsWrapper title="Advanced Options">
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
              className="mt-1 w-full px-3 py-2 rounded-lg border border-slate-200 focus:border-purple-400 focus:ring-2 focus:ring-purple-100 outline-none text-sm font-mono"
            />
            <p className="text-xs text-slate-500 mt-1">
              Positions protein origin away from ligand center. Use cautiously.
            </p>
          </div>
        )}
      </AdvancedOptionsWrapper>

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
              Designing {approach === 'asymmetric' ? 'Binder' : 'Dimer'}...
            </>
          ) : (
            <>
              <span className="material-symbols-outlined">link</span>
              Design {selectedLigandName} {approach === 'asymmetric' ? 'Binder' : 'Dimer'}
            </>
          )}
        </button>
        {!health && (
          <p className="text-center text-sm text-slate-500 mt-2">
            Backend service unavailable
          </p>
        )}

        {/* Expected Results - Quality Tiers */}
        <div className="mt-4 p-4 rounded-xl bg-gradient-to-br from-slate-50 to-slate-100 border border-slate-200">
          <div className="flex items-center gap-2 mb-3">
            <span className="material-symbols-outlined text-slate-500 text-sm">checklist</span>
            <p className="text-xs text-slate-600 font-semibold uppercase tracking-wide">Success Criteria</p>
          </div>

          {/* Quality tiers */}
          <div className="space-y-3">
            {/* Excellent */}
            <div className="p-2 rounded-lg bg-emerald-50 border border-emerald-100">
              <div className="flex items-center gap-2 mb-1">
                <span className="w-2 h-2 rounded-full bg-emerald-500" />
                <span className="text-xs font-medium text-emerald-800">Excellent</span>
              </div>
              <div className="grid grid-cols-2 gap-x-4 gap-y-0.5 text-xs text-emerald-700 ml-4">
                <div>GNINA: {EXPECTED_RESULTS.excellent.affinity}</div>
                <div>Contacts: {EXPECTED_RESULTS.excellent.contacts}</div>
                <div>Identity: {EXPECTED_RESULTS.excellent.identity}</div>
                <div>Anti-homo: {EXPECTED_RESULTS.excellent.antiHomo}</div>
              </div>
            </div>

            {/* Good */}
            <div className="p-2 rounded-lg bg-blue-50 border border-blue-100">
              <div className="flex items-center gap-2 mb-1">
                <span className="w-2 h-2 rounded-full bg-blue-500" />
                <span className="text-xs font-medium text-blue-800">Good</span>
              </div>
              <div className="grid grid-cols-2 gap-x-4 gap-y-0.5 text-xs text-blue-700 ml-4">
                <div>GNINA: {EXPECTED_RESULTS.good.affinity}</div>
                <div>Contacts: {EXPECTED_RESULTS.good.contacts}</div>
                <div>Identity: {EXPECTED_RESULTS.good.identity}</div>
                <div>Anti-homo: {EXPECTED_RESULTS.good.antiHomo}</div>
              </div>
            </div>

            {/* Minimum */}
            <div className="p-2 rounded-lg bg-amber-50 border border-amber-100">
              <div className="flex items-center gap-2 mb-1">
                <span className="w-2 h-2 rounded-full bg-amber-500" />
                <span className="text-xs font-medium text-amber-800">Minimum Viable</span>
              </div>
              <div className="grid grid-cols-2 gap-x-4 gap-y-0.5 text-xs text-amber-700 ml-4">
                <div>GNINA: {EXPECTED_RESULTS.minimum.affinity}</div>
                <div>Contacts: {EXPECTED_RESULTS.minimum.contacts}</div>
                <div>Identity: {EXPECTED_RESULTS.minimum.identity}</div>
                <div>H-bonds: {EXPECTED_RESULTS.minimum.hbonds}</div>
              </div>
            </div>
          </div>

          <p className="text-xs text-slate-500 mt-3 text-center">
            Pipeline: Backbone → LigandMPNN → Validation → PLIP Analysis
          </p>
        </div>
      </div>
    </div>
  );
}
