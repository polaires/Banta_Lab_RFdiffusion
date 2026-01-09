'use client';

import { useState } from 'react';
import { FormSection, FormField, FormRow } from './shared/FormSection';
import { PdbUploader } from './shared/PdbUploader';
import { QUALITY_PRESETS, RFD3Request, TaskFormProps } from './shared/types';

type QualityPreset = keyof typeof QUALITY_PRESETS;

interface CatalyticResidue {
  chain: string;
  residue: number;
  name: string;
}

export function EnzymeForm({ onSubmit, isSubmitting, health }: TaskFormProps) {
  // Required
  const [pdbContent, setPdbContent] = useState<string | null>(null);
  const [pdbFileName, setPdbFileName] = useState<string | null>(null);
  const [ligandCodes, setLigandCodes] = useState('');
  const [proteinLength, setProteinLength] = useState('150');

  // Catalytic residues (unindex)
  const [catalyticResidues, setCatalyticResidues] = useState<CatalyticResidue[]>([]);
  const [newCatChain, setNewCatChain] = useState('A');
  const [newCatResidue, setNewCatResidue] = useState('');
  const [newCatName, setNewCatName] = useState('');

  // Fixed atoms on catalytic residues
  const [fixedAtomType, setFixedAtomType] = useState('BKBN');

  // Options
  const [qualityPreset, setQualityPreset] = useState<QualityPreset>('Balanced');
  const [isNonLoopy, setIsNonLoopy] = useState(true);
  const [numDesigns, setNumDesigns] = useState(1);
  const [seed, setSeed] = useState<string>('');

  const addCatalyticResidue = () => {
    if (newCatResidue) {
      const residueNum = parseInt(newCatResidue, 10);
      if (!isNaN(residueNum)) {
        const exists = catalyticResidues.some(
          (r) => r.chain === newCatChain && r.residue === residueNum
        );
        if (!exists) {
          setCatalyticResidues([
            ...catalyticResidues,
            { chain: newCatChain, residue: residueNum, name: newCatName },
          ]);
          setNewCatResidue('');
          setNewCatName('');
        }
      }
    }
  };

  const removeCatalyticResidue = (index: number) => {
    setCatalyticResidues(catalyticResidues.filter((_, i) => i !== index));
  };

  const handleSubmit = async () => {
    const preset = QUALITY_PRESETS[qualityPreset];

    // Build unindex string: comma-separated residue positions
    // Format: Chain:Residue[ResName] -> simplified to just "chain residue"
    const unindexParts = catalyticResidues.map((r) => `${r.chain}${r.residue}`);
    const unindex = unindexParts.join(',');

    // Build fixed atoms selection for catalytic residues
    const selectFixedAtoms: Record<string, string> = {};
    for (const res of catalyticResidues) {
      const key = `${res.chain}${res.residue}`;
      selectFixedAtoms[key] = fixedAtomType;
    }

    const request: RFD3Request = {
      pdb_content: pdbContent || undefined,
      ligand: ligandCodes.toUpperCase(),
      length: proteinLength,
      unindex,
      num_designs: numDesigns,
      is_non_loopy: isNonLoopy,
      num_timesteps: preset.num_timesteps,
      step_scale: preset.step_scale,
      gamma_0: preset.gamma_0,
    };

    if (Object.keys(selectFixedAtoms).length > 0 && fixedAtomType) {
      request.select_fixed_atoms = selectFixedAtoms;
    }

    if (seed) {
      request.seed = parseInt(seed, 10);
    }

    await onSubmit(request);
  };

  const isValid =
    pdbContent !== null &&
    ligandCodes.trim() !== '' &&
    catalyticResidues.length > 0 &&
    proteinLength.trim() !== '';

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex items-center gap-3 pb-4 border-b border-slate-200">
        <div className="w-10 h-10 rounded-lg bg-slate-100 flex items-center justify-center">
          <span className="material-symbols-outlined text-slate-600">biotech</span>
        </div>
        <div>
          <h2 className="font-semibold text-slate-900">Enzyme Scaffold Design</h2>
          <p className="text-sm text-slate-500">Build a protein scaffold around an active site (theozyme)</p>
        </div>
      </div>

      {/* Info Banner */}
      <div className="p-4 rounded-xl bg-slate-50 border border-slate-200">
        <div className="flex items-start gap-3">
          <span className="material-symbols-outlined text-slate-600 text-xl">info</span>
          <div className="text-sm text-slate-700">
            <strong>Theozyme scaffolding:</strong> Your PDB should contain the active site model
            with catalytic residues and substrate/ligand. The model will design a protein scaffold
            that positions these elements correctly for catalysis.
          </div>
        </div>
      </div>

      {/* Theozyme PDB - Required */}
      <FormSection
        title="Theozyme Structure"
        description="Upload the PDB with catalytic residues and ligand(s)"
        required
      >
        <PdbUploader
          label="Theozyme PDB"
          description="Active site model with catalytic residues and substrate"
          required
          value={pdbContent}
          fileName={pdbFileName}
          onChange={(content, name) => {
            setPdbContent(content);
            setPdbFileName(name);
          }}
        />
      </FormSection>

      {/* Ligand Codes - Required */}
      <FormSection
        title="Ligand/Substrate Codes"
        description="3-letter codes for substrates and cofactors in the theozyme"
        required
      >
        <FormField label="Ligand Codes" required hint="Comma-separated (e.g., SUB,NAD)">
          <input
            type="text"
            value={ligandCodes}
            onChange={(e) => setLigandCodes(e.target.value.toUpperCase())}
            placeholder="SUB,NAD,MG"
            className="w-full px-4 py-2.5 rounded-xl border border-slate-200 focus:border-blue-400 focus:ring-2 focus:ring-blue-100 outline-none transition-all font-mono uppercase"
          />
        </FormField>
      </FormSection>

      {/* Catalytic Residues - Required */}
      <FormSection
        title="Catalytic Residues"
        description="Residues that form the active site. Their positions will be inferred by the model (unindex)."
        required
      >
        <div className="space-y-4">
          {/* Current residues */}
          {catalyticResidues.length > 0 && (
            <div className="flex flex-wrap gap-2 mb-3">
              {catalyticResidues.map((res, i) => (
                <div
                  key={`${res.chain}-${res.residue}`}
                  className="flex items-center gap-1 px-3 py-1.5 rounded-lg bg-slate-100 text-slate-700 text-sm font-medium"
                >
                  <span>
                    {res.chain}:{res.residue}
                    {res.name && <span className="text-slate-600 ml-1">({res.name})</span>}
                  </span>
                  <button
                    onClick={() => removeCatalyticResidue(i)}
                    className="ml-1 hover:text-red-600 transition-colors"
                  >
                    <span className="material-symbols-outlined text-sm">close</span>
                  </button>
                </div>
              ))}
            </div>
          )}

          {/* Add new residue */}
          <div className="flex gap-2">
            <input
              type="text"
              value={newCatChain}
              onChange={(e) => setNewCatChain(e.target.value.toUpperCase())}
              placeholder="Chain"
              maxLength={1}
              className="w-20 px-3 py-2 rounded-xl border border-slate-200 focus:border-blue-400 focus:ring-2 focus:ring-blue-100 outline-none transition-all text-center"
            />
            <input
              type="number"
              value={newCatResidue}
              onChange={(e) => setNewCatResidue(e.target.value)}
              placeholder="Residue #"
              onKeyDown={(e) => e.key === 'Enter' && addCatalyticResidue()}
              className="w-28 px-3 py-2 rounded-xl border border-slate-200 focus:border-blue-400 focus:ring-2 focus:ring-blue-100 outline-none transition-all"
            />
            <input
              type="text"
              value={newCatName}
              onChange={(e) => setNewCatName(e.target.value.toUpperCase())}
              placeholder="Name (opt)"
              maxLength={3}
              className="w-24 px-3 py-2 rounded-xl border border-slate-200 focus:border-blue-400 focus:ring-2 focus:ring-blue-100 outline-none transition-all font-mono"
            />
            <button
              onClick={addCatalyticResidue}
              disabled={!newCatResidue}
              className={`px-4 py-2 rounded-xl font-medium transition-all ${
                newCatResidue
                  ? 'bg-blue-600 text-white hover:bg-blue-700'
                  : 'bg-slate-100 text-slate-400 cursor-not-allowed'
              }`}
            >
              Add
            </button>
          </div>

          {catalyticResidues.length === 0 && (
            <p className="text-xs text-amber-600 flex items-center gap-1">
              <span className="material-symbols-outlined text-sm">warning</span>
              Add at least one catalytic residue for enzyme scaffolding
            </p>
          )}

          {/* Info about unindex */}
          <div className="p-3 rounded-xl bg-slate-50 border border-slate-200 text-sm text-slate-600">
            <strong>Note:</strong> "Unindex" means the exact positions of these residues are not
            fixed — the model will determine their optimal placement within the scaffold.
          </div>
        </div>
      </FormSection>

      {/* Fixed Atom Selection */}
      <FormSection
        title="Fixed Atoms on Catalytic Residues"
        description="Which atoms of the catalytic residues should stay fixed during design"
      >
        <div className="grid grid-cols-2 sm:grid-cols-4 gap-2">
          {[
            { value: 'BKBN', label: 'Backbone', desc: 'N, CA, C, O' },
            { value: 'ALL', label: 'All atoms', desc: 'Entire residue' },
            { value: 'TIP', label: 'Tip atoms', desc: 'Functional groups' },
            { value: '', label: 'None', desc: 'Fully flexible' },
          ].map((opt) => (
            <button
              key={opt.value}
              onClick={() => setFixedAtomType(opt.value)}
              className={`p-3 rounded-xl border-2 text-left transition-all ${
                fixedAtomType === opt.value
                  ? 'border-blue-400 bg-blue-50'
                  : 'border-slate-200 hover:border-slate-300'
              }`}
            >
              <div className="font-medium text-sm text-slate-900">{opt.label}</div>
              <div className="text-xs text-slate-500 mt-0.5">{opt.desc}</div>
            </button>
          ))}
        </div>
      </FormSection>

      {/* Scaffold Length - Required */}
      <FormSection
        title="Scaffold Length"
        description="Total length of the enzyme scaffold to design"
        required
      >
        <FormRow>
          <FormField label="Length" required hint="e.g., 150 or 120-180">
            <input
              type="text"
              value={proteinLength}
              onChange={(e) => setProteinLength(e.target.value)}
              placeholder="150"
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
              Design {numDesigns} Enzyme Scaffold{numDesigns > 1 ? 's' : ''}
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
