'use client';

import { useState } from 'react';
import { Dna, Info, Loader2, Rocket } from 'lucide-react';
import { FormSection, FormField, FormRow } from './shared/FormSection';
import { PdbUploader } from './shared/PdbUploader';
import { LengthRangeInput } from './shared/LengthRangeInput';
import { QualityPresetSelector, QualityPreset, QualityParams } from './shared/QualityPresetSelector';
import { AdvancedOptionsWrapper } from './shared/AdvancedOptionsWrapper';
import { QUALITY_PRESETS, RFD3Request, TaskFormProps } from './shared/types';

const NA_PRESETS = [
  { id: 'dsDNA', label: 'Double-stranded DNA', chains: 'A,B', description: 'Two complementary DNA chains' },
  { id: 'ssDNA', label: 'Single-stranded DNA', chains: 'A', description: 'One DNA chain' },
  { id: 'dsRNA', label: 'Double-stranded RNA', chains: 'A,B', description: 'Two complementary RNA chains' },
  { id: 'ssRNA', label: 'Single-stranded RNA', chains: 'A', description: 'One RNA chain' },
];

export function NucleicAcidForm({ onSubmit, isSubmitting, health }: TaskFormProps) {
  // Required
  const [pdbContent, setPdbContent] = useState<string | null>(null);
  const [pdbFileName, setPdbFileName] = useState<string | null>(null);
  const [naChains, setNaChains] = useState('A,B');
  const [binderLength, setBinderLength] = useState('80-120');

  // Origin token (positioning relative to NA)
  const [useOriToken, setUseOriToken] = useState(false);
  const [oriTokenX, setOriTokenX] = useState('0');
  const [oriTokenY, setOriTokenY] = useState('0');
  const [oriTokenZ, setOriTokenZ] = useState('20');

  // H-bond conditioning
  const [useHbonds, setUseHbonds] = useState(false);
  const [hbondDonors, setHbondDonors] = useState('');
  const [hbondAcceptors, setHbondAcceptors] = useState('');

  // Fixed atoms
  const [useFixedAtoms, setUseFixedAtoms] = useState(false);
  const [fixedAtoms, setFixedAtoms] = useState('');

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
    // Build contig: {binder_length},/0,{chain1}{residues},/0,{chain2}{residues}...
    // For now, simplify to just chains
    const chainList = naChains.split(',').map((c) => c.trim());
    const naContig = chainList.map((c) => `/0,${c}`).join('');
    const contig = `${binderLength}${naContig}`;

    const request: RFD3Request = {
      contig,
      pdb_content: pdbContent || undefined,
      num_designs: numDesigns,
      is_non_loopy: isNonLoopy,
      num_timesteps: qualityParams.num_timesteps,
      step_scale: qualityParams.step_scale,
      gamma_0: qualityParams.gamma_0,
    };

    // Origin token for positioning
    if (useOriToken) {
      request.ori_token = [
        parseFloat(oriTokenX) || 0,
        parseFloat(oriTokenY) || 0,
        parseFloat(oriTokenZ) || 0,
      ];
    }

    // H-bond conditioning
    if (useHbonds) {
      if (hbondDonors) {
        const donors: Record<string, string> = {};
        for (const chain of chainList) {
          donors[chain] = hbondDonors;
        }
        request.select_hbond_donor = donors;
      }
      if (hbondAcceptors) {
        const acceptors: Record<string, string> = {};
        for (const chain of chainList) {
          acceptors[chain] = hbondAcceptors;
        }
        request.select_hbond_acceptor = acceptors;
      }
    }

    // Fixed atoms
    if (useFixedAtoms && fixedAtoms) {
      const fixed: Record<string, string> = {};
      for (const chain of chainList) {
        fixed[chain] = fixedAtoms;
      }
      request.select_fixed_atoms = fixed;
    }

    if (seed) {
      request.seed = parseInt(seed, 10);
    }

    await onSubmit(request);
  };

  const isValid =
    pdbContent !== null &&
    naChains.trim() !== '' &&
    binderLength.trim() !== '';

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex items-center gap-3 pb-4 border-b border-border">
        <div className="w-10 h-10 rounded-lg bg-muted flex items-center justify-center">
          <Dna className="w-5 h-5 text-muted-foreground" />
        </div>
        <div>
          <h2 className="font-semibold text-foreground">Nucleic Acid Binder Design</h2>
          <p className="text-sm text-muted-foreground">Design a protein that binds DNA or RNA</p>
        </div>
      </div>

      {/* Info Banner */}
      <div className="p-4 rounded-xl bg-muted/50 border border-border">
        <div className="flex items-start gap-3">
          <Info className="w-5 h-5 text-muted-foreground flex-shrink-0 mt-0.5" />
          <div className="text-sm text-foreground">
            <strong>Requirements:</strong> Your PDB must contain the DNA or RNA structure you want
            to target. The model will design a protein that interacts with the nucleic acid.
          </div>
        </div>
      </div>

      {/* Input PDB - Required */}
      <FormSection
        title="Structure with Nucleic Acid"
        description="Upload PDB containing the DNA or RNA you want to target"
        required
      >
        <PdbUploader
          label="PDB with DNA/RNA"
          description="Must contain the nucleic acid chains"
          required
          value={pdbContent}
          fileName={pdbFileName}
          onChange={(content, name) => {
            setPdbContent(content);
            setPdbFileName(name);
          }}
        />
      </FormSection>

      {/* NA Chain Selection - Required */}
      <FormSection
        title="Nucleic Acid Chains"
        description="Specify which chains in the PDB are the target nucleic acid"
        required
      >
        <div className="space-y-4">
          {/* Presets */}
          <div className="grid grid-cols-2 gap-2">
            {NA_PRESETS.map((preset) => (
              <button
                key={preset.id}
                onClick={() => setNaChains(preset.chains)}
                className={`p-3 rounded-xl border-2 text-left transition-all ${
                  naChains === preset.chains
                    ? 'border-primary bg-primary/10'
                    : 'border-border hover:border-border'
                }`}
              >
                <div className="font-medium text-sm text-foreground">{preset.label}</div>
                <div className="text-xs text-muted-foreground mt-0.5">{preset.description}</div>
              </button>
            ))}
          </div>

          <FormField label="Chain IDs" required hint="Comma-separated chain letters">
            <input
              type="text"
              value={naChains}
              onChange={(e) => setNaChains(e.target.value.toUpperCase())}
              placeholder="A,B"
              className="w-full px-4 py-2.5 rounded-xl border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none transition-all font-mono"
            />
          </FormField>
        </div>
      </FormSection>

      {/* Binder Length - Required */}
      <FormSection
        title="Binder Length"
        description="Length of the DNA/RNA-binding protein to design"
        required
      >
        <div className="flex gap-4 items-start">
          <div className="flex-1">
            <LengthRangeInput
              value={binderLength}
              onChange={setBinderLength}
              label="Length Range"
              placeholder="80-120"
              hint="e.g., 80-120 for flexible length"
            />
          </div>
          <FormField label="# Designs" className="w-28">
            <input
              type="number"
              value={numDesigns}
              onChange={(e) => setNumDesigns(Math.max(1, parseInt(e.target.value) || 1))}
              min={1}
              max={10}
              className="w-full px-3 py-2 rounded-lg border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none transition-all text-sm"
            />
          </FormField>
        </div>
      </FormSection>

      {/* Origin Token - Optional */}
      <FormSection
        title="Origin Positioning (Optional)"
        description="Position the protein binder relative to the nucleic acid using an origin token"
      >
        <div className={`p-4 rounded-xl border transition-all ${
          useOriToken ? 'bg-muted/50 border-border' : 'bg-muted/50 border-border'
        }`}>
          <label className="flex items-center gap-3 cursor-pointer mb-4">
            <input
              type="checkbox"
              checked={useOriToken}
              onChange={(e) => setUseOriToken(e.target.checked)}
              className="w-5 h-5 rounded border-border text-primary focus:ring-primary"
            />
            <div>
              <div className="font-medium text-sm text-foreground">Use Origin Token</div>
              <div className="text-xs text-muted-foreground">
                Define a 3D position for the binder relative to NA center of mass
              </div>
            </div>
          </label>

          {useOriToken && (
            <div className="grid grid-cols-3 gap-3">
              <FormField label="X (Å)">
                <input
                  type="number"
                  value={oriTokenX}
                  onChange={(e) => setOriTokenX(e.target.value)}
                  className="w-full px-4 py-2.5 rounded-xl border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none transition-all"
                />
              </FormField>
              <FormField label="Y (Å)">
                <input
                  type="number"
                  value={oriTokenY}
                  onChange={(e) => setOriTokenY(e.target.value)}
                  className="w-full px-4 py-2.5 rounded-xl border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none transition-all"
                />
              </FormField>
              <FormField label="Z (Å)">
                <input
                  type="number"
                  value={oriTokenZ}
                  onChange={(e) => setOriTokenZ(e.target.value)}
                  className="w-full px-4 py-2.5 rounded-xl border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none transition-all"
                />
              </FormField>
            </div>
          )}
        </div>
      </FormSection>

      {/* H-bond Conditioning - Optional */}
      <FormSection
        title="Hydrogen Bond Conditioning (Optional)"
        description="Guide which NA atoms should form H-bonds with the designed protein"
      >
        <div className={`p-4 rounded-xl border transition-all ${
          useHbonds ? 'bg-muted/50 border-border' : 'bg-muted/50 border-border'
        }`}>
          <label className="flex items-center gap-3 cursor-pointer mb-4">
            <input
              type="checkbox"
              checked={useHbonds}
              onChange={(e) => setUseHbonds(e.target.checked)}
              className="w-5 h-5 rounded border-border text-primary focus:ring-primary"
            />
            <div>
              <div className="font-medium text-sm text-foreground">Enable H-bond Conditioning</div>
              <div className="text-xs text-muted-foreground">
                Specify which NA atoms should donate/accept H-bonds
              </div>
            </div>
          </label>

          {useHbonds && (
            <div className="space-y-3">
              <FormField label="H-bond Donor Atoms" hint="Atoms that should donate H-bonds (e.g., N1,N3)">
                <input
                  type="text"
                  value={hbondDonors}
                  onChange={(e) => setHbondDonors(e.target.value)}
                  placeholder="N1,N3,N6"
                  className="w-full px-4 py-2.5 rounded-xl border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none transition-all font-mono text-sm"
                />
              </FormField>
              <FormField label="H-bond Acceptor Atoms" hint="Atoms that should accept H-bonds (e.g., O2,O4)">
                <input
                  type="text"
                  value={hbondAcceptors}
                  onChange={(e) => setHbondAcceptors(e.target.value)}
                  placeholder="O2,O4,O6"
                  className="w-full px-4 py-2.5 rounded-xl border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none transition-all font-mono text-sm"
                />
              </FormField>
            </div>
          )}
        </div>
      </FormSection>

      {/* Fixed Atoms - Optional */}
      <FormSection
        title="Fixed Atoms (Optional)"
        description="Which nucleic acid atoms should stay fixed in space"
      >
        <div className={`p-4 rounded-xl border transition-all ${
          useFixedAtoms ? 'bg-muted/50 border-border' : 'bg-muted/50 border-border'
        }`}>
          <label className="flex items-center gap-3 cursor-pointer mb-3">
            <input
              type="checkbox"
              checked={useFixedAtoms}
              onChange={(e) => setUseFixedAtoms(e.target.checked)}
              className="w-5 h-5 rounded border-border text-primary focus:ring-primary"
            />
            <div>
              <div className="font-medium text-sm text-foreground">Specify Fixed Atoms</div>
              <div className="text-xs text-muted-foreground">
                By default, phosphate backbone is fixed
              </div>
            </div>
          </label>

          {useFixedAtoms && (
            <input
              type="text"
              value={fixedAtoms}
              onChange={(e) => setFixedAtoms(e.target.value)}
              placeholder="P,O5',O3' or ALL"
              className="w-full px-4 py-2.5 rounded-xl border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none transition-all font-mono text-sm"
            />
          )}
        </div>
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
        <label className="flex items-center gap-3 p-3 rounded-lg bg-muted/50 hover:bg-muted cursor-pointer transition-colors">
          <input
            type="checkbox"
            checked={isNonLoopy}
            onChange={(e) => setIsNonLoopy(e.target.checked)}
            className="w-4 h-4 rounded border-border text-primary focus:ring-primary"
          />
          <div>
            <div className="font-medium text-sm text-foreground">Non-loopy Mode</div>
            <div className="text-xs text-muted-foreground">
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
              className="w-full px-3 py-2 rounded-lg border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none transition-all text-sm"
            />
          </FormField>
        </FormRow>
      </AdvancedOptionsWrapper>

      {/* Submit Button */}
      <div className="pt-4 border-t border-border">
        <button
          onClick={handleSubmit}
          disabled={!isValid || isSubmitting || !health}
          className={`w-full py-3 px-6 rounded-xl font-semibold text-white transition-all flex items-center justify-center gap-2 ${
            isValid && !isSubmitting && !!health
              ? 'bg-primary hover:bg-primary/90 shadow-lg shadow-primary/20'
              : 'bg-muted cursor-not-allowed'
          }`}
        >
          {isSubmitting ? (
            <>
              <Loader2 className="w-5 h-5 animate-spin" />
              Submitting...
            </>
          ) : (
            <>
              <Rocket className="w-5 h-5" />
              Design {numDesigns} NA Binder{numDesigns > 1 ? 's' : ''}
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
