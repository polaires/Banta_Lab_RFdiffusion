'use client';

import { useState, useMemo } from 'react';
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

// Extract HETATM ligand codes from PDB content
function extractLigandsFromPdb(pdbContent: string): string[] {
  const ligands = new Set<string>();
  const lines = pdbContent.split('\n');
  for (const line of lines) {
    if (line.startsWith('HETATM')) {
      const resName = line.substring(17, 20).trim();
      if (resName && !['HOH', 'WAT', 'DOD', 'H2O'].includes(resName)) {
        ligands.add(resName);
      }
    }
  }
  return Array.from(ligands).sort();
}

// Metal/ion codes that are single atoms
const METAL_CODES = new Set([
  'ZN', 'MG', 'CA', 'FE', 'MN', 'CO', 'CU', 'NI',
  'LA', 'GD', 'EU', 'TB', 'YB', 'LU',
  'K', 'NA', 'CL',
]);

// Replace ligand code in PDB content
function replaceLigandInPdb(pdbContent: string, sourceCode: string, targetCode: string): string {
  if (!sourceCode || !targetCode || sourceCode === targetCode) return pdbContent;

  const isMetal = METAL_CODES.has(targetCode.toUpperCase());
  const lines = pdbContent.split('\n');

  return lines.map(line => {
    if (line.startsWith('HETATM') || line.startsWith('ATOM')) {
      const resName = line.substring(17, 20).trim();
      if (resName === sourceCode) {
        let paddedLine = line.padEnd(80);
        const paddedResName = targetCode.padStart(3);
        paddedLine = paddedLine.substring(0, 17) + paddedResName + paddedLine.substring(20);

        if (isMetal) {
          const atomName = targetCode.length === 1
            ? ' ' + targetCode + '  '
            : targetCode.padStart(2) + '  ';
          paddedLine = paddedLine.substring(0, 12) + atomName + paddedLine.substring(16);
          const element = targetCode.length === 1 ? ' ' + targetCode : targetCode.substring(0, 2);
          paddedLine = paddedLine.substring(0, 76) + element + paddedLine.substring(78);
        }
        return paddedLine.trimEnd();
      }
    }
    return line;
  }).join('\n');
}

export function SmallMoleculeForm({ onSubmit, isSubmitting, health }: TaskFormProps) {
  // Required
  const [pdbContent, setPdbContent] = useState<string | null>(null);
  const [pdbFileName, setPdbFileName] = useState<string | null>(null);
  const [ligandCode, setLigandCode] = useState('');
  const [proteinLength, setProteinLength] = useState('100');

  // Ligand replacement
  const [replaceLigand, setReplaceLigand] = useState(false);
  const [sourceLigand, setSourceLigand] = useState('');
  const [targetLigand, setTargetLigand] = useState('');

  // Template refinement (partial diffusion)
  const [useTemplateRefinement, setUseTemplateRefinement] = useState(false);
  const [partialT, setPartialT] = useState('10');

  // Advanced options - collapsed by default
  const [showAdvanced, setShowAdvanced] = useState(false);

  // RASA conditioning
  const [useBuried, setUseBuried] = useState(false);
  const [buriedAtoms, setBuriedAtoms] = useState('');
  const [useExposed, setUseExposed] = useState(false);
  const [exposedAtoms, setExposedAtoms] = useState('');
  const [usePartiallyBuried, setUsePartiallyBuried] = useState(false);
  const [partiallyBuriedAtoms, setPartiallyBuriedAtoms] = useState('');

  // Fixed atoms
  const [fixedAtomSelection, setFixedAtomSelection] = useState('');

  // Unindex (coordinating residues)
  const [useUnindex, setUseUnindex] = useState(false);
  const [unindexResidues, setUnindexResidues] = useState('');

  // Options
  const [qualityPreset, setQualityPreset] = useState<QualityPreset>('Balanced');
  const [isNonLoopy, setIsNonLoopy] = useState(true);
  const [numDesigns, setNumDesigns] = useState(1);
  const [seed, setSeed] = useState<string>('');

  // Detected ligands
  const detectedLigands = useMemo(() => {
    if (!pdbContent) return [];
    return extractLigandsFromPdb(pdbContent);
  }, [pdbContent]);

  // Current ligand code (for display and validation)
  const currentLigand = replaceLigand ? targetLigand : ligandCode;
  const isMetal = METAL_CODES.has(currentLigand.toUpperCase());

  const handleSubmit = async () => {
    const preset = QUALITY_PRESETS[qualityPreset];

    let finalPdbContent = pdbContent;
    let finalLigandCode = ligandCode.toUpperCase();

    if (replaceLigand && sourceLigand && targetLigand && pdbContent) {
      finalPdbContent = replaceLigandInPdb(pdbContent, sourceLigand, targetLigand);
      finalLigandCode = targetLigand.toUpperCase();
    }

    const request: RFD3Request = {
      pdb_content: finalPdbContent || undefined,
      ligand: finalLigandCode,
      length: proteinLength,
      num_designs: numDesigns,
      is_non_loopy: isNonLoopy,
      num_timesteps: preset.num_timesteps,
      step_scale: preset.step_scale,
      gamma_0: preset.gamma_0,
    };

    if (useBuried && buriedAtoms) {
      request.select_buried = { [finalLigandCode]: buriedAtoms };
    }
    if (useExposed && exposedAtoms) {
      request.select_exposed = { [finalLigandCode]: exposedAtoms };
    }
    if (usePartiallyBuried && partiallyBuriedAtoms) {
      request.select_partially_buried = { [finalLigandCode]: partiallyBuriedAtoms };
    }
    if (fixedAtomSelection) {
      request.select_fixed_atoms = { [finalLigandCode]: fixedAtomSelection };
    }
    if (useTemplateRefinement && partialT) {
      request.partial_t = parseFloat(partialT);
    }
    if (useUnindex && unindexResidues) {
      request.unindex = unindexResidues;
    }
    if (seed) {
      request.seed = parseInt(seed, 10);
    }

    await onSubmit(request);
  };

  const isValid =
    pdbContent !== null &&
    proteinLength.trim() !== '' &&
    (replaceLigand
      ? sourceLigand.trim() !== '' && targetLigand.trim() !== ''
      : ligandCode.trim() !== '');

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

      {/* Structure Input Section - Consolidated */}
      <FormSection
        title="Input Structure"
        description="Upload PDB containing the ligand"
        required
      >
        <div className="space-y-4">
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

          {/* Detected ligands - inline */}
          {detectedLigands.length > 0 && (
            <div className="flex items-center gap-3 p-3 rounded-lg bg-slate-50 border border-slate-200">
              <span className="text-xs text-slate-500 font-medium whitespace-nowrap">Detected:</span>
              <div className="flex flex-wrap gap-1.5">
                {detectedLigands.map((lig) => (
                  <button
                    key={lig}
                    onClick={() => {
                      if (!replaceLigand) setLigandCode(lig);
                      else setSourceLigand(lig);
                    }}
                    className={`px-2 py-0.5 text-xs rounded font-mono font-medium transition-all ${
                      (replaceLigand ? sourceLigand : ligandCode) === lig
                        ? 'bg-slate-700 text-white'
                        : 'bg-white text-slate-600 hover:bg-slate-100 border border-slate-300'
                    }`}
                  >
                    {lig}
                  </button>
                ))}
              </div>
            </div>
          )}

          {/* Replace ligand toggle - inline */}
          {detectedLigands.length > 0 && (
            <label className="flex items-center gap-3 cursor-pointer group">
              <input
                type="checkbox"
                checked={replaceLigand}
                onChange={(e) => {
                  setReplaceLigand(e.target.checked);
                  if (e.target.checked && detectedLigands.length > 0) {
                    setSourceLigand(detectedLigands[0]);
                  }
                }}
                className="w-4 h-4 rounded border-slate-300 text-slate-700 focus:ring-slate-500"
              />
              <span className="text-sm text-slate-600 group-hover:text-slate-900">
                Replace ligand in PDB (e.g., Ca → Gd)
              </span>
            </label>
          )}

          {/* Replacement controls - show when enabled */}
          {replaceLigand && (
            <div className="grid grid-cols-2 gap-3 p-3 rounded-lg bg-slate-50 border border-slate-200">
              <div>
                <label className="text-xs text-slate-500 font-medium">Source</label>
                <select
                  value={sourceLigand}
                  onChange={(e) => setSourceLigand(e.target.value)}
                  className="mt-1 w-full px-3 py-2 rounded-lg border border-slate-200 focus:border-slate-400 focus:ring-1 focus:ring-slate-200 outline-none bg-white font-mono text-sm"
                >
                  <option value="">Select...</option>
                  {detectedLigands.map((lig) => (
                    <option key={lig} value={lig}>{lig}</option>
                  ))}
                </select>
              </div>
              <div>
                <label className="text-xs text-slate-500 font-medium">Target</label>
                <input
                  type="text"
                  value={targetLigand}
                  onChange={(e) => setTargetLigand(e.target.value.toUpperCase())}
                  placeholder="GD, LA, EU..."
                  maxLength={4}
                  className="mt-1 w-full px-3 py-2 rounded-lg border border-slate-200 focus:border-slate-400 focus:ring-1 focus:ring-slate-200 outline-none font-mono text-sm uppercase"
                />
              </div>
            </div>
          )}

          {/* Template refinement toggle - inline */}
          <label className="flex items-center gap-3 cursor-pointer group">
            <input
              type="checkbox"
              checked={useTemplateRefinement}
              onChange={(e) => setUseTemplateRefinement(e.target.checked)}
              className="w-4 h-4 rounded border-slate-300 text-slate-700 focus:ring-slate-500"
            />
            <span className="text-sm text-slate-600 group-hover:text-slate-900">
              Refine from template (partial diffusion)
            </span>
            {useTemplateRefinement && (
              <input
                type="number"
                value={partialT}
                onChange={(e) => setPartialT(e.target.value)}
                min={1}
                max={50}
                className="w-16 px-2 py-1 rounded border border-slate-200 text-sm text-center"
                title="Noise level (lower = more similar to input)"
              />
            )}
          </label>
        </div>
      </FormSection>

      {/* Ligand Selection */}
      <FormSection
        title="Ligand"
        description={isMetal ? "Metal binding requires template with coordinating residues" : "Select the ligand to design around"}
        required
      >
        <div className="space-y-3">
          {!replaceLigand && (
            <div className="grid grid-cols-2 gap-3">
              {/* Dropdown with grouped options */}
              <div>
                <label className="text-xs text-slate-500 font-medium">Common Ligands</label>
                <select
                  value={ligandCode}
                  onChange={(e) => setLigandCode(e.target.value)}
                  className="mt-1 w-full px-3 py-2 rounded-lg border border-slate-200 focus:border-slate-400 focus:ring-1 focus:ring-slate-200 outline-none bg-white text-sm"
                >
                  <option value="">Select ligand...</option>
                  {Object.entries(ligandsByCategory).map(([category, ligands]) => (
                    <optgroup key={category} label={category}>
                      {ligands.map((lig) => (
                        <option key={lig.id} value={lig.id}>
                          {lig.id} - {lig.label}
                        </option>
                      ))}
                    </optgroup>
                  ))}
                </select>
              </div>

              {/* Custom code input */}
              <div>
                <label className="text-xs text-slate-500 font-medium">Or enter code</label>
                <input
                  type="text"
                  value={ligandCode}
                  onChange={(e) => setLigandCode(e.target.value.toUpperCase())}
                  placeholder="e.g., ATP"
                  maxLength={4}
                  className="mt-1 w-full px-3 py-2 rounded-lg border border-slate-200 focus:border-slate-400 focus:ring-1 focus:ring-slate-200 outline-none text-sm font-mono uppercase"
                />
              </div>
            </div>
          )}

          {replaceLigand && sourceLigand && targetLigand && (
            <div className="text-sm text-slate-600">
              Using <code className="font-mono bg-slate-100 px-1.5 py-0.5 rounded">{targetLigand}</code> (replacing {sourceLigand})
            </div>
          )}

          {/* Metal hint - inline */}
          {isMetal && (
            <p className="text-xs text-slate-500">
              For proper coordination, use a template PDB with Asp/Glu residues and enable refinement mode above.
            </p>
          )}
        </div>
      </FormSection>

      {/* Design Parameters */}
      <FormSection title="Design Parameters" required>
        <div className="space-y-4">
          <div className="grid grid-cols-2 gap-4">
            <FormField label="Length" hint="e.g., 100 or 80-120">
              <input
                type="text"
                value={proteinLength}
                onChange={(e) => setProteinLength(e.target.value)}
                placeholder="100"
                className="w-full px-3 py-2 rounded-lg border border-slate-200 focus:border-slate-400 focus:ring-1 focus:ring-slate-200 outline-none transition-all"
              />
            </FormField>
            <FormField label="# Designs">
              <input
                type="number"
                value={numDesigns}
                onChange={(e) => setNumDesigns(Math.max(1, parseInt(e.target.value) || 1))}
                min={1}
                max={10}
                className="w-full px-3 py-2 rounded-lg border border-slate-200 focus:border-slate-400 focus:ring-1 focus:ring-slate-200 outline-none transition-all"
              />
            </FormField>
          </div>

          {/* Quality preset buttons */}
          <div>
            <label className="text-xs text-slate-500 font-medium">Quality</label>
            <div className="mt-1 grid grid-cols-2 sm:grid-cols-4 gap-2">
              {(Object.keys(QUALITY_PRESETS) as QualityPreset[]).map((preset) => (
                <button
                  key={preset}
                  onClick={() => setQualityPreset(preset)}
                  className={`p-2.5 rounded-lg border text-left transition-all ${
                    qualityPreset === preset
                      ? 'border-slate-400 bg-slate-100'
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

      {/* Advanced Options - Collapsible */}
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
            {/* RASA Conditioning */}
            <div className="space-y-3">
              <p className="text-xs text-slate-500 font-medium uppercase tracking-wider">Binding Pocket (RASA)</p>

              <label className="flex items-center gap-3 cursor-pointer">
                <input
                  type="checkbox"
                  checked={useBuried}
                  onChange={(e) => setUseBuried(e.target.checked)}
                  className="w-4 h-4 rounded border-slate-300 text-slate-700"
                />
                <span className="text-sm text-slate-600">Buried atoms</span>
                {useBuried && (
                  <input
                    type="text"
                    value={buriedAtoms}
                    onChange={(e) => setBuriedAtoms(e.target.value)}
                    placeholder="ALL or atom names"
                    className="flex-1 px-3 py-1.5 rounded border border-slate-200 text-sm font-mono"
                  />
                )}
              </label>

              <label className="flex items-center gap-3 cursor-pointer">
                <input
                  type="checkbox"
                  checked={useExposed}
                  onChange={(e) => setUseExposed(e.target.checked)}
                  className="w-4 h-4 rounded border-slate-300 text-slate-700"
                />
                <span className="text-sm text-slate-600">Exposed atoms</span>
                {useExposed && (
                  <input
                    type="text"
                    value={exposedAtoms}
                    onChange={(e) => setExposedAtoms(e.target.value)}
                    placeholder="ALL or atom names"
                    className="flex-1 px-3 py-1.5 rounded border border-slate-200 text-sm font-mono"
                  />
                )}
              </label>

              <label className="flex items-center gap-3 cursor-pointer">
                <input
                  type="checkbox"
                  checked={usePartiallyBuried}
                  onChange={(e) => setUsePartiallyBuried(e.target.checked)}
                  className="w-4 h-4 rounded border-slate-300 text-slate-700"
                />
                <span className="text-sm text-slate-600">Partially buried</span>
                {usePartiallyBuried && (
                  <input
                    type="text"
                    value={partiallyBuriedAtoms}
                    onChange={(e) => setPartiallyBuriedAtoms(e.target.value)}
                    placeholder="Atom names"
                    className="flex-1 px-3 py-1.5 rounded border border-slate-200 text-sm font-mono"
                  />
                )}
              </label>
            </div>

            {/* Fixed Atoms */}
            <div className="space-y-2">
              <p className="text-xs text-slate-500 font-medium uppercase tracking-wider">Fixed Atoms</p>
              <select
                value={fixedAtomSelection}
                onChange={(e) => setFixedAtomSelection(e.target.value)}
                className="w-full px-3 py-2 rounded-lg border border-slate-200 bg-white text-sm"
              >
                {ATOM_SELECTION_OPTIONS.map((opt) => (
                  <option key={opt.value} value={opt.value}>
                    {opt.label} - {opt.description}
                  </option>
                ))}
              </select>
            </div>

            {/* Coordinating Residues */}
            <div className="space-y-2">
              <label className="flex items-center gap-3 cursor-pointer">
                <input
                  type="checkbox"
                  checked={useUnindex}
                  onChange={(e) => setUseUnindex(e.target.checked)}
                  className="w-4 h-4 rounded border-slate-300 text-slate-700"
                />
                <span className="text-sm text-slate-600">Coordinating residues (unindex)</span>
              </label>
              {useUnindex && (
                <input
                  type="text"
                  value={unindexResidues}
                  onChange={(e) => setUnindexResidues(e.target.value)}
                  placeholder="e.g., A10,A15,A20"
                  className="w-full px-3 py-2 rounded-lg border border-slate-200 text-sm font-mono"
                />
              )}
            </div>

            {/* Other Options */}
            <div className="flex items-center gap-6">
              <label className="flex items-center gap-2 cursor-pointer">
                <input
                  type="checkbox"
                  checked={isNonLoopy}
                  onChange={(e) => setIsNonLoopy(e.target.checked)}
                  className="w-4 h-4 rounded border-slate-300 text-slate-700"
                />
                <span className="text-sm text-slate-600">Non-loopy mode</span>
              </label>

              <div className="flex items-center gap-2">
                <span className="text-sm text-slate-600">Seed:</span>
                <input
                  type="number"
                  value={seed}
                  onChange={(e) => setSeed(e.target.value)}
                  placeholder="Random"
                  className="w-24 px-2 py-1 rounded border border-slate-200 text-sm"
                />
              </div>
            </div>
          </div>
        )}
      </div>

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
              Design {numDesigns} {currentLigand || 'Ligand'} Binder{numDesigns > 1 ? 's' : ''}
            </>
          )}
        </button>
        {!health && (
          <p className="text-center text-sm text-slate-500 mt-2">
            Backend service unavailable
          </p>
        )}
      </div>
    </div>
  );
}
