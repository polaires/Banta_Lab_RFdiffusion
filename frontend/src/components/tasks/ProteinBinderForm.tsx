'use client';

import { useState } from 'react';
import { Network, X, Lightbulb, Loader2, PlayCircle } from 'lucide-react';
import { FormSection, FormField, FormRow } from './shared/FormSection';
import { PdbUploader } from './shared/PdbUploader';
import { LengthRangeInput } from './shared/LengthRangeInput';
import { QualityPresetSelector, QualityPreset, QualityParams } from './shared/QualityPresetSelector';
import { AdvancedOptionsWrapper } from './shared/AdvancedOptionsWrapper';
import { QUALITY_PRESETS, RFD3Request, TaskFormProps } from './shared/types';

interface HotspotResidue {
  chain: string;
  residue: number;
}

// Extended props for binder form with results callback
interface ProteinBinderFormProps extends TaskFormProps {
  onBinderResult?: (result: BinderDesignResult) => void;
}

// Response from protein_binder_design API
export interface BinderDesignResult {
  status: 'completed' | 'error';
  statistics: {
    generated: number;
    mpnn_designed: number;
    esm_passed: number;
    relaxed: number;
    interface_analyzed: number;
    passed_filters: number;
    returned: number;
  };
  designs: Array<{
    binder_sequence: string;
    mpnn_score: number;
    esm_perplexity: number;
    esm_confidence: number;
    interface_contacts: number;
    interface_hbonds: number;
    buried_sasa: number;
    packstat: number;
    rank: number;
    pdb_content?: string;
  }>;
  error?: string;
}

export function ProteinBinderForm({ onSubmit, isSubmitting, health, onBinderResult }: ProteinBinderFormProps) {
  // PDB file
  const [pdbContent, setPdbContent] = useState<string | null>(null);
  const [pdbFileName, setPdbFileName] = useState<string | null>(null);

  // Target definition
  const [targetChain, setTargetChain] = useState('A');
  const [targetStart, setTargetStart] = useState('');
  const [targetEnd, setTargetEnd] = useState('');

  // Binder settings
  const [binderLength, setBinderLength] = useState('60-80');

  // Hotspots
  const [hotspots, setHotspots] = useState<HotspotResidue[]>([]);
  const [newHotspotChain, setNewHotspotChain] = useState('A');
  const [newHotspotResidue, setNewHotspotResidue] = useState('');

  // Quality & options
  const [qualityPreset, setQualityPreset] = useState<QualityPreset>('Binder Optimized');
  const [qualityParams, setQualityParams] = useState<QualityParams>(QUALITY_PRESETS['Binder Optimized']);
  const [isNonLoopy, setIsNonLoopy] = useState(true);
  const [numDesigns, setNumDesigns] = useState(4);  // Default to 4 for pipeline
  const [seed, setSeed] = useState<string>('');

  const handleQualityChange = (preset: QualityPreset, params: QualityParams) => {
    setQualityPreset(preset);
    setQualityParams(params);
  };

  // Pipeline mode toggle
  const [usePipeline, setUsePipeline] = useState(true);

  const addHotspot = () => {
    if (newHotspotResidue) {
      const residueNum = parseInt(newHotspotResidue, 10);
      if (!isNaN(residueNum)) {
        const exists = hotspots.some(
          (h) => h.chain === newHotspotChain && h.residue === residueNum
        );
        if (!exists) {
          setHotspots([...hotspots, { chain: newHotspotChain, residue: residueNum }]);
          setNewHotspotResidue('');
        }
      }
    }
  };

  const removeHotspot = (index: number) => {
    setHotspots(hotspots.filter((_, i) => i !== index));
  };

  const handleSubmit = async () => {
    // Build hotspot selections
    const selectHotspots: Record<string, string> = {};
    if (hotspots.length > 0) {
      const hotspotsByChain: Record<string, number[]> = {};
      for (const h of hotspots) {
        if (!hotspotsByChain[h.chain]) hotspotsByChain[h.chain] = [];
        hotspotsByChain[h.chain].push(h.residue);
      }
      for (const [chain, residues] of Object.entries(hotspotsByChain)) {
        selectHotspots[chain] = residues.join(',');
      }
    }

    if (usePipeline) {
      // Use the new protein_binder_design pipeline
      const request: RFD3Request = {
        task: 'protein_binder_design',
        pdb_content: pdbContent || undefined,
        // Target region
        contig: targetStart && targetEnd
          ? `${targetChain}${targetStart}-${targetEnd}`
          : targetChain,
        // Binder length
        chain_length: binderLength,
        num_designs: numDesigns,
        // Optional hotspots
        ...(Object.keys(selectHotspots).length > 0 && { select_hotspots: selectHotspots }),
        // Seed
        ...(seed && { seed: parseInt(seed, 10) }),
      };

      await onSubmit(request);
    } else {
      // Legacy RFD3-only mode (for testing/comparison)
      const targetRange = targetStart && targetEnd ? `${targetStart}-${targetEnd}` : '';
      const contig = `${binderLength},/0,${targetChain}${targetRange}`;

      const request: RFD3Request = {
        contig,
        pdb_content: pdbContent || undefined,
        num_designs: numDesigns,
        is_non_loopy: isNonLoopy,
        num_timesteps: qualityParams.num_timesteps,
        step_scale: qualityParams.step_scale,
        gamma_0: qualityParams.gamma_0,
        infer_ori_strategy: 'hotspots',
        ...(Object.keys(selectHotspots).length > 0 && { select_hotspots: selectHotspots }),
        ...(seed && { seed: parseInt(seed, 10) }),
      };

      await onSubmit(request);
    }
  };

  const isValid =
    pdbContent !== null &&
    targetChain.trim() !== '' &&
    binderLength.trim() !== '';

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex items-center gap-3 pb-4 border-b border-border">
        <div className="w-9 h-9 rounded-lg bg-muted flex items-center justify-center">
          <Network className="w-5 h-5 text-muted-foreground" />
        </div>
        <div>
          <h2 className="font-semibold text-foreground">Protein Binder Design</h2>
          <p className="text-sm text-muted-foreground">Design a protein that binds to a target protein</p>
        </div>
      </div>

      {/* Target PDB Upload - Required */}
      <FormSection
        title="Target Protein"
        description="Upload the PDB file containing your target protein"
        required
      >
        <PdbUploader
          label="Target PDB File"
          description="The protein you want to design a binder for"
          required
          value={pdbContent}
          fileName={pdbFileName}
          onChange={(content, name) => {
            setPdbContent(content);
            setPdbFileName(name);
          }}
        />
      </FormSection>

      {/* Target Definition - Required */}
      <FormSection
        title="Target Region"
        description="Specify which chain and residues to target for binding"
        required
      >
        <FormRow>
          <FormField label="Target Chain" required>
            <input
              type="text"
              value={targetChain}
              onChange={(e) => setTargetChain(e.target.value.toUpperCase())}
              placeholder="A"
              maxLength={1}
              className="w-full px-3 py-2 rounded-lg border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none transition-all text-sm"
            />
          </FormField>
          <FormField label="Residue Range" hint="Optional: limit to specific residues">
            <div className="flex items-center gap-2">
              <input
                type="number"
                value={targetStart}
                onChange={(e) => setTargetStart(e.target.value)}
                placeholder="Start"
                className="flex-1 px-3 py-2 rounded-lg border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none transition-all text-sm"
              />
              <span className="text-muted-foreground">-</span>
              <input
                type="number"
                value={targetEnd}
                onChange={(e) => setTargetEnd(e.target.value)}
                placeholder="End"
                className="flex-1 px-3 py-2 rounded-lg border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none transition-all text-sm"
              />
            </div>
          </FormField>
        </FormRow>
      </FormSection>

      {/* Binder Length - Required */}
      <FormSection
        title="Binder Length"
        description="Specify the length of the designed binder protein"
        required
      >
        <div className="flex gap-4 items-start">
          <div className="flex-1">
            <LengthRangeInput
              value={binderLength}
              onChange={setBinderLength}
              label="Length Range"
              placeholder="60-80"
              hint="e.g., 60-80 for flexible length"
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

      {/* Hotspots - Recommended */}
      <FormSection
        title="Hotspot Residues"
        description="Residues on the target that the binder should contact. Highly recommended for better designs."
      >
        {/* Hotspot list */}
        {hotspots.length > 0 && (
          <div className="flex flex-wrap gap-2 mb-3">
            {hotspots.map((h, i) => (
              <div
                key={`${h.chain}-${h.residue}`}
                className="flex items-center gap-1 px-2.5 py-1 rounded-lg bg-muted text-foreground text-sm font-medium"
              >
                <span>{h.chain}:{h.residue}</span>
                <button
                  onClick={() => removeHotspot(i)}
                  className="ml-1 hover:text-red-600 transition-colors"
                >
                  <X className="w-4 h-4" />
                </button>
              </div>
            ))}
          </div>
        )}

        {/* Add hotspot */}
        <div className="flex gap-2">
          <input
            type="text"
            value={newHotspotChain}
            onChange={(e) => setNewHotspotChain(e.target.value.toUpperCase())}
            placeholder="Chain"
            maxLength={1}
            className="w-16 px-3 py-2 rounded-lg border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none transition-all text-center text-sm"
          />
          <input
            type="number"
            value={newHotspotResidue}
            onChange={(e) => setNewHotspotResidue(e.target.value)}
            placeholder="Residue #"
            onKeyDown={(e) => e.key === 'Enter' && addHotspot()}
            className="flex-1 px-3 py-2 rounded-lg border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none transition-all text-sm"
          />
          <button
            onClick={addHotspot}
            disabled={!newHotspotResidue}
            className={`px-4 py-2 rounded-lg font-medium transition-all text-sm ${
              newHotspotResidue
                ? 'bg-primary text-white hover:bg-primary/90'
                : 'bg-muted text-muted-foreground cursor-not-allowed'
            }`}
          >
            Add
          </button>
        </div>

        {hotspots.length === 0 && (
          <p className="text-xs text-muted-foreground mt-2 flex items-center gap-1">
            <Lightbulb className="w-4 h-4 text-muted-foreground" />
            Adding 2-5 hotspot residues is strongly recommended for successful binder design
          </p>
        )}
      </FormSection>

      {/* Quality Settings (only for legacy mode) */}
      {!usePipeline && (
        <FormSection
          title="Quality Settings"
          description="Binder Optimized preset is recommended for most use cases"
        >
          <QualityPresetSelector
            value={qualityPreset}
            onChange={handleQualityChange}
            showDescription
          />
        </FormSection>
      )}

      {/* Pipeline Mode */}
      <FormSection
        title="Design Pipeline"
        description="Choose between full pipeline with validation or quick RFD3-only mode"
      >
        <div className="space-y-2">
          <label className={`flex items-center gap-3 p-3 rounded-lg cursor-pointer transition-all ${
            usePipeline ? 'bg-primary/10 border-2 border-primary' : 'bg-muted/50 border-2 border-transparent hover:bg-muted'
          }`}>
            <input
              type="radio"
              checked={usePipeline}
              onChange={() => setUsePipeline(true)}
              className="w-4 h-4 text-primary focus:ring-primary"
            />
            <div className="flex-1">
              <div className="font-medium text-sm text-foreground flex items-center gap-2">
                Full Pipeline
                <span className="text-xs px-2 py-0.5 rounded-full bg-primary/20 text-primary">Recommended</span>
              </div>
              <div className="text-xs text-muted-foreground">
                RFD3 + MPNN sequence design + ESM-3 validation + interface analysis
              </div>
            </div>
          </label>

          <label className={`flex items-center gap-3 p-3 rounded-lg cursor-pointer transition-all ${
            !usePipeline ? 'bg-secondary border-2 border-secondary-foreground/30' : 'bg-muted/50 border-2 border-transparent hover:bg-muted'
          }`}>
            <input
              type="radio"
              checked={!usePipeline}
              onChange={() => setUsePipeline(false)}
              className="w-4 h-4 text-primary focus:ring-primary"
            />
            <div className="flex-1">
              <div className="font-medium text-sm text-foreground">RFD3 Only (Quick)</div>
              <div className="text-xs text-muted-foreground">
                Fast backbone generation without sequence/validation (for testing)
              </div>
            </div>
          </label>
        </div>
      </FormSection>

      {/* Structure Options (only for legacy mode) */}
      {!usePipeline && (
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
      )}

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
              {usePipeline ? 'Running Pipeline...' : 'Submitting...'}
            </>
          ) : (
            <>
              {usePipeline ? <Network className="w-5 h-5" /> : <PlayCircle className="w-5 h-5" />}
              {usePipeline
                ? `Design & Validate ${numDesigns} Binder${numDesigns > 1 ? 's' : ''}`
                : `Design ${numDesigns} Binder${numDesigns > 1 ? 's' : ''}`
              }
            </>
          )}
        </button>
        {!health && (
          <p className="text-center text-sm text-muted-foreground mt-2">
            Connect to backend to enable design
          </p>
        )}
      </div>
    </div>
  );
}
