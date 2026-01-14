'use client';

import { useState } from 'react';
import { FormSection, FormField, FormRow } from './shared/FormSection';
import { PdbUploader } from './shared/PdbUploader';
import { LengthRangeInput } from './shared/LengthRangeInput';
import { QualityPresetSelector, QualityPreset, QualityParams } from './shared/QualityPresetSelector';
import { AdvancedOptionsWrapper } from './shared/AdvancedOptionsWrapper';
import { QUALITY_PRESETS, SYMMETRY_OPTIONS, RFD3Request, TaskFormProps } from './shared/types';

export function SymmetricForm({ onSubmit, isSubmitting, health }: TaskFormProps) {
  // Required
  const [symmetry, setSymmetry] = useState('C2');
  const [lengthPerSubunit, setLengthPerSubunit] = useState('80');

  // Optional motif scaffolding
  const [useMotif, setUseMotif] = useState(false);
  const [pdbContent, setPdbContent] = useState<string | null>(null);
  const [pdbFileName, setPdbFileName] = useState<string | null>(null);
  const [motifChain, setMotifChain] = useState('A');
  const [motifResidues, setMotifResidues] = useState('');
  const [isSymmetricMotif, setIsSymmetricMotif] = useState(false);
  const [isUnsymMotif, setIsUnsymMotif] = useState('');

  // Options
  const [qualityPreset, setQualityPreset] = useState<QualityPreset>('Balanced');
  const [qualityParams, setQualityParams] = useState<QualityParams>(QUALITY_PRESETS.Balanced);
  const [isNonLoopy, setIsNonLoopy] = useState(true);
  const [numDesigns, setNumDesigns] = useState(1);
  const [seed, setSeed] = useState<string>('');

  const handleQualityChange = (preset: QualityPreset, params: QualityParams) => {
    setQualityPreset(preset);
    setQualityParams(params);
  };

  const handleSubmit = async () => {
    const selectedSymmetry = SYMMETRY_OPTIONS.find((s) => s.id === symmetry);

    let contig: string | undefined;

    if (useMotif && pdbContent && motifResidues) {
      // Motif scaffolding mode
      // Format: scaffolding around a motif with symmetry
      // e.g., {length_before},{motif_chain}{residues},{length_after}
      contig = `${lengthPerSubunit}`;
    } else {
      // Simple symmetric design - just length
      contig = undefined; // Use length parameter instead
    }

    const request: RFD3Request = {
      length: useMotif ? undefined : lengthPerSubunit,
      contig: useMotif ? contig : undefined,
      pdb_content: useMotif ? pdbContent || undefined : undefined,
      num_designs: numDesigns,
      is_non_loopy: isNonLoopy,
      num_timesteps: qualityParams.num_timesteps,
      step_scale: qualityParams.step_scale,
      gamma_0: qualityParams.gamma_0,
      symmetry: {
        id: symmetry,
        is_symmetric_motif: useMotif ? isSymmetricMotif : undefined,
        is_unsym_motif: useMotif && isUnsymMotif ? isUnsymMotif : undefined,
      },
    };

    if (seed) {
      request.seed = parseInt(seed, 10);
    }

    await onSubmit(request);
  };

  const isValid = symmetry && lengthPerSubunit.trim() !== '' && (!useMotif || pdbContent !== null);
  const selectedSymmetry = SYMMETRY_OPTIONS.find((s) => s.id === symmetry);
  const totalSubunits = selectedSymmetry?.subunits || 1;

  // Group symmetry options by category
  const symmetryByCategory = SYMMETRY_OPTIONS.reduce((acc, opt) => {
    if (!acc[opt.category]) acc[opt.category] = [];
    acc[opt.category].push(opt);
    return acc;
  }, {} as Record<string, typeof SYMMETRY_OPTIONS[number][]>);

  // Memory warnings for high symmetries
  const isHighSymmetry = totalSubunits >= 12;
  const estimatedLength = parseInt(lengthPerSubunit, 10) || 0;
  const totalLength = estimatedLength * totalSubunits;

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex items-center gap-3 pb-4 border-b border-slate-200">
        <div className="w-10 h-10 rounded-lg bg-slate-100 flex items-center justify-center">
          <span className="material-symbols-outlined text-slate-600">hexagon</span>
        </div>
        <div>
          <h2 className="font-semibold text-slate-900">Symmetric Oligomer Design</h2>
          <p className="text-sm text-slate-500">Design homo-oligomeric proteins (dimers, trimers, etc.)</p>
        </div>
      </div>

      {/* Symmetry Selection - Required */}
      <FormSection
        title="Symmetry Type"
        description="Select the symmetry of your oligomer"
        required
      >
        <div className="space-y-4">
          {Object.entries(symmetryByCategory).map(([category, options]) => (
            <div key={category}>
              <h4 className="text-xs font-medium text-slate-500 uppercase tracking-wider mb-2">
                {category}
              </h4>
              <div className="grid grid-cols-3 sm:grid-cols-4 md:grid-cols-6 gap-2">
                {options.map((opt) => (
                  <button
                    key={opt.id}
                    onClick={() => setSymmetry(opt.id)}
                    className={`p-3 rounded-xl border-2 text-center transition-all ${
                      symmetry === opt.id
                        ? 'border-blue-400 bg-blue-50'
                        : 'border-slate-200 hover:border-slate-300'
                    }`}
                  >
                    <div className="font-bold text-lg text-slate-900">{opt.id}</div>
                    <div className="text-xs text-slate-500">{opt.subunits} units</div>
                  </button>
                ))}
              </div>
            </div>
          ))}
        </div>
      </FormSection>

      {/* Size Settings - Required */}
      <FormSection
        title="Subunit Size"
        description="Specify the length per subunit"
        required
      >
        <div className="flex gap-4 items-start">
          <div className="flex-1">
            <LengthRangeInput
              value={lengthPerSubunit}
              onChange={setLengthPerSubunit}
              label="Length per Subunit"
              placeholder="80 or 60-100"
              hint="Number of residues in each subunit"
            />
          </div>
          <FormField label="# Designs" className="w-28">
            <input
              type="number"
              value={numDesigns}
              onChange={(e) => setNumDesigns(Math.max(1, parseInt(e.target.value) || 1))}
              min={1}
              max={10}
              className="w-full px-3 py-2 rounded-lg border border-slate-200 focus:border-blue-400 focus:ring-2 focus:ring-blue-100 outline-none transition-all text-sm"
            />
          </FormField>
        </div>

        {/* Total size info */}
        <div className={`p-3 rounded-xl border mt-3 ${
          isHighSymmetry ? 'bg-amber-50 border-amber-200' : 'bg-slate-50 border-slate-200'
        }`}>
          <div className="flex items-start gap-2">
            <span className={`material-symbols-outlined text-lg ${
              isHighSymmetry ? 'text-amber-600' : 'text-slate-600'
            }`}>
              {isHighSymmetry ? 'warning' : 'info'}
            </span>
            <div className={`text-sm ${isHighSymmetry ? 'text-amber-800' : 'text-slate-700'}`}>
              <strong>{symmetry}</strong> with {estimatedLength} residues/subunit
              = <strong>{totalLength} total residues</strong> ({totalSubunits} subunits)
              {isHighSymmetry && (
                <div className="mt-1 text-amber-700">
                  High symmetry orders require significant GPU memory. Consider reducing subunit length if you encounter memory issues.
                </div>
              )}
            </div>
          </div>
        </div>
      </FormSection>

      {/* Motif Scaffolding - Optional */}
      <FormSection
        title="Motif Scaffolding (Optional)"
        description="Design symmetric oligomers around an existing motif"
      >
        <label className="flex items-center gap-3 p-3 rounded-xl bg-slate-50 hover:bg-slate-100 cursor-pointer transition-colors mb-4">
          <input
            type="checkbox"
            checked={useMotif}
            onChange={(e) => setUseMotif(e.target.checked)}
            className="w-5 h-5 rounded border-slate-300 text-blue-600 focus:ring-blue-500"
          />
          <div>
            <div className="font-medium text-sm text-slate-900">Use Motif Scaffolding</div>
            <div className="text-xs text-slate-500">
              Build symmetric structure around an existing motif
            </div>
          </div>
        </label>

        {useMotif && (
          <div className="space-y-4 p-4 rounded-xl bg-slate-50 border border-slate-200">
            <PdbUploader
              label="Motif PDB"
              description="PDB containing the motif to scaffold"
              required
              value={pdbContent}
              fileName={pdbFileName}
              onChange={(content, name) => {
                setPdbContent(content);
                setPdbFileName(name);
              }}
            />

            <FormRow>
              <FormField label="Motif Chain">
                <input
                  type="text"
                  value={motifChain}
                  onChange={(e) => setMotifChain(e.target.value.toUpperCase())}
                  placeholder="A"
                  maxLength={1}
                  className="w-full px-4 py-2.5 rounded-xl border border-slate-200 focus:border-blue-400 focus:ring-2 focus:ring-blue-100 outline-none transition-all"
                />
              </FormField>
              <FormField label="Motif Residues" hint="e.g., 1-20 or 5,10,15">
                <input
                  type="text"
                  value={motifResidues}
                  onChange={(e) => setMotifResidues(e.target.value)}
                  placeholder="1-20"
                  className="w-full px-4 py-2.5 rounded-xl border border-slate-200 focus:border-blue-400 focus:ring-2 focus:ring-blue-100 outline-none transition-all"
                />
              </FormField>
            </FormRow>

            <div className="space-y-2">
              <label className="flex items-center gap-3 p-3 rounded-xl bg-white hover:bg-slate-100 cursor-pointer transition-colors">
                <input
                  type="checkbox"
                  checked={isSymmetricMotif}
                  onChange={(e) => setIsSymmetricMotif(e.target.checked)}
                  className="w-5 h-5 rounded border-slate-300 text-blue-600 focus:ring-blue-500"
                />
                <div>
                  <div className="font-medium text-sm text-slate-900">Symmetric Motif</div>
                  <div className="text-xs text-slate-500">
                    Motif already has the target symmetry
                  </div>
                </div>
              </label>

              <FormField label="Asymmetric Chains" hint="Chains not to symmetrize (comma-separated)">
                <input
                  type="text"
                  value={isUnsymMotif}
                  onChange={(e) => setIsUnsymMotif(e.target.value.toUpperCase())}
                  placeholder="B,C"
                  className="w-full px-4 py-2.5 rounded-xl border border-slate-200 focus:border-blue-400 focus:ring-2 focus:ring-blue-100 outline-none transition-all"
                />
              </FormField>
            </div>
          </div>
        )}
      </FormSection>

      {/* Quality Settings */}
      <FormSection
        title="Quality Settings"
        description="Higher quality takes longer but produces better designs"
      >
        <QualityPresetSelector
          value={qualityPreset}
          onChange={handleQualityChange}
          showDescription
        />
      </FormSection>

      {/* Structure Options */}
      <FormSection title="Structure Options">
        <label className="flex items-center gap-3 p-3 rounded-lg bg-slate-50 hover:bg-slate-100 cursor-pointer transition-colors">
          <input
            type="checkbox"
            checked={isNonLoopy}
            onChange={(e) => setIsNonLoopy(e.target.checked)}
            className="w-4 h-4 rounded border-slate-300 text-blue-600 focus:ring-blue-500"
          />
          <div>
            <div className="font-medium text-sm text-slate-800">Non-loopy Mode</div>
            <div className="text-xs text-slate-500">
              Produces cleaner secondary structures (recommended)
            </div>
          </div>
        </label>
      </FormSection>

      {/* Advanced Options */}
      <AdvancedOptionsWrapper title="Advanced Options">
        <FormRow>
          <FormField label="Random Seed" hint="For reproducible results">
            <input
              type="number"
              value={seed}
              onChange={(e) => setSeed(e.target.value)}
              placeholder="Optional"
              className="w-full px-3 py-2 rounded-lg border border-slate-200 focus:border-blue-400 focus:ring-2 focus:ring-blue-100 outline-none transition-all text-sm"
            />
          </FormField>
        </FormRow>
      </AdvancedOptionsWrapper>

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
              Generate {numDesigns} {symmetry} Oligomer{numDesigns > 1 ? 's' : ''}
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
