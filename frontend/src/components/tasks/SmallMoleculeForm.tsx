'use client';

import { useState } from 'react';
import { FormSection, FormField, FormRow } from './shared/FormSection';
import { PdbUploader } from './shared/PdbUploader';
import {
  QUALITY_PRESETS,
  COMMON_LIGANDS,
  ATOM_SELECTION_OPTIONS,
  RFD3Request,
  TaskFormProps,
} from './shared/types';

type QualityPreset = keyof typeof QUALITY_PRESETS;

export function SmallMoleculeForm({ onSubmit, isSubmitting, health }: TaskFormProps) {
  // Required
  const [pdbContent, setPdbContent] = useState<string | null>(null);
  const [pdbFileName, setPdbFileName] = useState<string | null>(null);
  const [ligandCode, setLigandCode] = useState('');
  const [proteinLength, setProteinLength] = useState('100');

  // RASA conditioning (optional)
  const [useBuried, setUseBuried] = useState(false);
  const [buriedAtoms, setBuriedAtoms] = useState('');
  const [useExposed, setUseExposed] = useState(false);
  const [exposedAtoms, setExposedAtoms] = useState('');
  const [usePartiallyBuried, setUsePartiallyBuried] = useState(false);
  const [partiallyBuriedAtoms, setPartiallyBuriedAtoms] = useState('');

  // Fixed atoms (optional)
  const [fixedAtomSelection, setFixedAtomSelection] = useState('');

  // Options
  const [qualityPreset, setQualityPreset] = useState<QualityPreset>('Balanced');
  const [isNonLoopy, setIsNonLoopy] = useState(true);
  const [numDesigns, setNumDesigns] = useState(1);
  const [seed, setSeed] = useState<string>('');

  const handleSubmit = async () => {
    const preset = QUALITY_PRESETS[qualityPreset];

    const request: RFD3Request = {
      pdb_content: pdbContent || undefined,
      ligand: ligandCode.toUpperCase(),
      length: proteinLength,
      num_designs: numDesigns,
      is_non_loopy: isNonLoopy,
      num_timesteps: preset.num_timesteps,
      step_scale: preset.step_scale,
      gamma_0: preset.gamma_0,
    };

    // RASA conditioning
    if (useBuried && buriedAtoms) {
      request.select_buried = { [ligandCode.toUpperCase()]: buriedAtoms };
    }
    if (useExposed && exposedAtoms) {
      request.select_exposed = { [ligandCode.toUpperCase()]: exposedAtoms };
    }
    if (usePartiallyBuried && partiallyBuriedAtoms) {
      request.select_partially_buried = { [ligandCode.toUpperCase()]: partiallyBuriedAtoms };
    }

    // Fixed atoms
    if (fixedAtomSelection) {
      request.select_fixed_atoms = { [ligandCode.toUpperCase()]: fixedAtomSelection };
    }

    if (seed) {
      request.seed = parseInt(seed, 10);
    }

    await onSubmit(request);
  };

  const isValid =
    pdbContent !== null &&
    ligandCode.trim() !== '' &&
    proteinLength.trim() !== '';

  // Group ligands by category
  const ligandsByCategory = COMMON_LIGANDS.reduce((acc, lig) => {
    if (!acc[lig.category]) acc[lig.category] = [];
    acc[lig.category].push(lig);
    return acc;
  }, {} as Record<string, typeof COMMON_LIGANDS[number][]>);

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex items-center gap-3 pb-4 border-b border-slate-200">
        <div className="w-10 h-10 rounded-lg bg-slate-100 flex items-center justify-center">
          <span className="material-symbols-outlined text-slate-600">science</span>
        </div>
        <div>
          <h2 className="font-semibold text-slate-900">Small Molecule Binder Design</h2>
          <p className="text-sm text-slate-500">Design a protein with a binding pocket for a ligand</p>
        </div>
      </div>

      {/* Info Banner */}
      <div className="p-4 rounded-xl bg-slate-50 border border-slate-200">
        <div className="flex items-start gap-3">
          <span className="material-symbols-outlined text-slate-600 text-xl">info</span>
          <div className="text-sm text-slate-700">
            <strong>Requirements:</strong> Your PDB must contain the ligand you want to bind.
            The model will design a protein scaffold around the ligand, creating a binding pocket.
          </div>
        </div>
      </div>

      {/* Input PDB - Required */}
      <FormSection
        title="Structure with Ligand"
        description="Upload PDB containing the ligand you want to design a binder for"
        required
      >
        <PdbUploader
          label="PDB with Ligand"
          description="Must contain the target ligand molecule"
          required
          value={pdbContent}
          fileName={pdbFileName}
          onChange={(content, name) => {
            setPdbContent(content);
            setPdbFileName(name);
          }}
        />
      </FormSection>

      {/* Ligand Selection - Required */}
      <FormSection
        title="Ligand Code"
        description="3-letter code identifying the ligand in the PDB file"
        required
      >
        <div className="space-y-4">
          <FormField label="Ligand Code" required>
            <input
              type="text"
              value={ligandCode}
              onChange={(e) => setLigandCode(e.target.value.toUpperCase())}
              placeholder="e.g., ATP, NAD, HEM"
              maxLength={4}
              className="w-full px-4 py-2.5 rounded-xl border border-slate-200 focus:border-blue-400 focus:ring-2 focus:ring-blue-100 outline-none transition-all font-mono uppercase"
            />
          </FormField>

          {/* Common ligand shortcuts */}
          <div className="space-y-3">
            <p className="text-xs text-slate-500 font-medium">Common Ligands:</p>
            {Object.entries(ligandsByCategory).map(([category, ligands]) => (
              <div key={category}>
                <p className="text-xs text-slate-400 uppercase tracking-wider mb-1">{category}</p>
                <div className="flex flex-wrap gap-1.5">
                  {ligands.map((lig) => (
                    <button
                      key={lig.id}
                      onClick={() => setLigandCode(lig.id)}
                      className={`px-2.5 py-1 text-xs rounded-lg font-medium transition-all ${
                        ligandCode === lig.id
                          ? 'bg-blue-600 text-white'
                          : 'bg-slate-100 text-slate-700 hover:bg-slate-200'
                      }`}
                    >
                      {lig.label}
                    </button>
                  ))}
                </div>
              </div>
            ))}
          </div>
        </div>
      </FormSection>

      {/* Protein Length - Required */}
      <FormSection
        title="Protein Length"
        description="Length of the binding protein to design"
        required
      >
        <FormRow>
          <FormField label="Length" required hint="e.g., 100 or 80-120">
            <input
              type="text"
              value={proteinLength}
              onChange={(e) => setProteinLength(e.target.value)}
              placeholder="100"
              className="w-full px-4 py-2.5 rounded-xl border border-slate-200 focus:border-blue-400 focus:ring-2 focus:ring-blue-100 outline-none transition-all"
            />
          </FormField>
          <FormField label="# Designs">
            <input
              type="number"
              value={numDesigns}
              onChange={(e) => setNumDesigns(Math.max(1, parseInt(e.target.value) || 1))}
              min={1}
              max={10}
              className="w-full px-4 py-2.5 rounded-xl border border-slate-200 focus:border-blue-400 focus:ring-2 focus:ring-blue-100 outline-none transition-all"
            />
          </FormField>
        </FormRow>
      </FormSection>

      {/* RASA Conditioning - Optional but powerful */}
      <FormSection
        title="Binding Pocket Design (RASA)"
        description="Control which ligand atoms should be buried vs exposed in the binding pocket"
      >
        <div className="space-y-4">
          {/* Buried atoms */}
          <div className={`p-4 rounded-xl border transition-all ${
            useBuried ? 'bg-slate-50 border-slate-200' : 'bg-slate-50 border-slate-200'
          }`}>
            <label className="flex items-center gap-3 cursor-pointer mb-3">
              <input
                type="checkbox"
                checked={useBuried}
                onChange={(e) => setUseBuried(e.target.checked)}
                className="w-5 h-5 rounded border-slate-300 text-blue-600 focus:ring-blue-500"
              />
              <div>
                <div className="font-medium text-sm text-slate-900">Buried Atoms</div>
                <div className="text-xs text-slate-500">
                  Atoms that should be fully enclosed in the protein
                </div>
              </div>
            </label>
            {useBuried && (
              <input
                type="text"
                value={buriedAtoms}
                onChange={(e) => setBuriedAtoms(e.target.value)}
                placeholder="Atom names (e.g., C1,C2,N1) or ALL"
                className="w-full px-4 py-2.5 rounded-xl border border-slate-200 focus:border-blue-400 focus:ring-2 focus:ring-blue-100 outline-none transition-all font-mono text-sm"
              />
            )}
          </div>

          {/* Exposed atoms */}
          <div className={`p-4 rounded-xl border transition-all ${
            useExposed ? 'bg-slate-50 border-slate-200' : 'bg-slate-50 border-slate-200'
          }`}>
            <label className="flex items-center gap-3 cursor-pointer mb-3">
              <input
                type="checkbox"
                checked={useExposed}
                onChange={(e) => setUseExposed(e.target.checked)}
                className="w-5 h-5 rounded border-slate-300 text-blue-600 focus:ring-blue-500"
              />
              <div>
                <div className="font-medium text-sm text-slate-900">Exposed Atoms</div>
                <div className="text-xs text-slate-500">
                  Atoms that should remain solvent accessible
                </div>
              </div>
            </label>
            {useExposed && (
              <input
                type="text"
                value={exposedAtoms}
                onChange={(e) => setExposedAtoms(e.target.value)}
                placeholder="Atom names (e.g., O1,O2) or ALL"
                className="w-full px-4 py-2.5 rounded-xl border border-slate-200 focus:border-blue-400 focus:ring-2 focus:ring-blue-100 outline-none transition-all font-mono text-sm"
              />
            )}
          </div>

          {/* Partially buried */}
          <div className={`p-4 rounded-xl border transition-all ${
            usePartiallyBuried ? 'bg-slate-50 border-slate-200' : 'bg-slate-50 border-slate-200'
          }`}>
            <label className="flex items-center gap-3 cursor-pointer mb-3">
              <input
                type="checkbox"
                checked={usePartiallyBuried}
                onChange={(e) => setUsePartiallyBuried(e.target.checked)}
                className="w-5 h-5 rounded border-slate-300 text-blue-600 focus:ring-blue-500"
              />
              <div>
                <div className="font-medium text-sm text-slate-900">Partially Buried Atoms</div>
                <div className="text-xs text-slate-500">
                  Atoms at the protein surface edge
                </div>
              </div>
            </label>
            {usePartiallyBuried && (
              <input
                type="text"
                value={partiallyBuriedAtoms}
                onChange={(e) => setPartiallyBuriedAtoms(e.target.value)}
                placeholder="Atom names"
                className="w-full px-4 py-2.5 rounded-xl border border-slate-200 focus:border-blue-400 focus:ring-2 focus:ring-blue-100 outline-none transition-all font-mono text-sm"
              />
            )}
          </div>
        </div>
      </FormSection>

      {/* Fixed Atoms - Optional */}
      <FormSection
        title="Fixed Atoms (Optional)"
        description="Which ligand atoms should stay fixed in space during design"
      >
        <FormField label="Fixed Atoms" hint="Leave empty for default (backbone atoms)">
          <select
            value={fixedAtomSelection}
            onChange={(e) => setFixedAtomSelection(e.target.value)}
            className="w-full px-4 py-2.5 rounded-xl border border-slate-200 focus:border-blue-400 focus:ring-2 focus:ring-blue-100 outline-none transition-all bg-white"
          >
            {ATOM_SELECTION_OPTIONS.map((opt) => (
              <option key={opt.value} value={opt.value}>
                {opt.label} - {opt.description}
              </option>
            ))}
          </select>
        </FormField>
      </FormSection>

      {/* Quality Settings */}
      <FormSection
        title="Quality Settings"
        description="Higher quality takes longer but produces better designs"
      >
        <div className="grid grid-cols-2 sm:grid-cols-4 gap-2">
          {(Object.keys(QUALITY_PRESETS) as QualityPreset[]).map((preset) => (
            <button
              key={preset}
              onClick={() => setQualityPreset(preset)}
              className={`p-3 rounded-xl border-2 text-left transition-all ${
                qualityPreset === preset
                  ? 'border-blue-400 bg-blue-50'
                  : 'border-slate-200 hover:border-slate-300'
              }`}
            >
              <div className="font-medium text-sm text-slate-900">{preset}</div>
              <div className="text-xs text-slate-500 mt-0.5">
                {QUALITY_PRESETS[preset].description}
              </div>
            </button>
          ))}
        </div>
      </FormSection>

      {/* Structure Options */}
      <FormSection title="Structure Options">
        <label className="flex items-center gap-3 p-3 rounded-xl bg-slate-50 hover:bg-slate-100 cursor-pointer transition-colors">
          <input
            type="checkbox"
            checked={isNonLoopy}
            onChange={(e) => setIsNonLoopy(e.target.checked)}
            className="w-5 h-5 rounded border-slate-300 text-blue-600 focus:ring-blue-500"
          />
          <div>
            <div className="font-medium text-sm text-slate-900">Non-loopy Mode</div>
            <div className="text-xs text-slate-500">
              Produces cleaner secondary structures (recommended)
            </div>
          </div>
        </label>
      </FormSection>

      {/* Advanced Options */}
      <FormSection title="Advanced" description="Additional options for fine-tuning">
        <FormRow>
          <FormField label="Random Seed" hint="For reproducible results">
            <input
              type="number"
              value={seed}
              onChange={(e) => setSeed(e.target.value)}
              placeholder="Optional"
              className="w-full px-4 py-2.5 rounded-xl border border-slate-200 focus:border-blue-400 focus:ring-2 focus:ring-blue-100 outline-none transition-all"
            />
          </FormField>
        </FormRow>
      </FormSection>

      {/* Submit Button */}
      <div className="pt-4 border-t border-slate-200">
        <button
          onClick={handleSubmit}
          disabled={!isValid || isSubmitting || !health}
          className={`w-full py-3 px-6 rounded-xl font-semibold text-white transition-all flex items-center justify-center gap-2 ${
            isValid && !isSubmitting && !!health
              ? 'bg-blue-600 hover:bg-blue-700 shadow-lg shadow-blue-600/20'
              : 'bg-slate-300 cursor-not-allowed'
          }`}
        >
          {isSubmitting ? (
            <>
              <span className="material-symbols-outlined animate-spin">progress_activity</span>
              Submitting...
            </>
          ) : (
            <>
              <span className="material-symbols-outlined">rocket_launch</span>
              Design {numDesigns} {ligandCode || 'Ligand'} Binder{numDesigns > 1 ? 's' : ''}
            </>
          )}
        </button>
        {!health && (
          <p className="text-center text-sm text-amber-600 mt-2">
            Backend service unavailable. Please check connection.
          </p>
        )}
      </div>
    </div>
  );
}
