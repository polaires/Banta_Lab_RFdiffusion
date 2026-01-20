'use client';

import { useState, useEffect, useCallback } from 'react';
import { Beaker, Info, X, AlertTriangle, Plus, Loader2, Rocket, Sparkles } from 'lucide-react';
import { useStore } from '@/lib/store';
import { FormSection, FormField, FormRow } from './shared/FormSection';
import { PdbUploader } from './shared/PdbUploader';
import { LengthRangeInput } from './shared/LengthRangeInput';
import { QualityPresetSelector, QualityPreset, QualityParams } from './shared/QualityPresetSelector';
import { AdvancedOptionsWrapper } from './shared/AdvancedOptionsWrapper';
import { MetalReplacementSection } from './shared/MetalReplacementSection';
import { RASAConditioningSection } from './shared/RASAConditioningSection';
import { HBondConditioningSection } from './shared/HBondConditioningSection';
import { QUALITY_PRESETS, RFD3Request, TaskFormProps } from './shared/types';
import {
  analyzeEnzymeStructure,
  replaceMetalInPdb,
  buildBuriedSelection,
  buildExposedSelection,
  buildHBondAcceptorSelection,
  buildHBondDonorSelection,
} from '@/lib/enzymeAnalysis';

interface CovalentBond {
  // Protein side: chain/resName/resId/atomName
  proteinChain: string;
  proteinRes: number;
  proteinResName: string;
  proteinAtom: string;
  // Ligand side: chain/resName/resId/atomName
  ligandChain: string;
  ligandRes: number;
  ligandResName: string;
  ligandAtom: string;
}

export function EnzymeForm({ onSubmit, isSubmitting, health }: TaskFormProps) {
  // Store hooks for catalytic suggestions and shared residue state
  const {
    catalyticSuggestions,
    setBottomPanelMode,
    enzymeCatalyticResidues,
    enzymeFixedAtomTypes,
    addEnzymeCatalyticResidue,
    removeEnzymeCatalyticResidue,
    setEnzymeCatalyticResidues,
    // Enzyme analysis state
    enzymeAnalysis,
    enzymeAnalysisLoading,
    setEnzymeAnalysis,
    setEnzymeAnalysisLoading,
    // Metal replacement
    metalReplacementEnabled,
    sourceMetal,
    targetMetal,
    coordinationMode,
    fixLigandPosition,
    metalReplacementPreset,
    excludedCatalyticResidues,
    setMetalReplacementEnabled,
    setSourceMetal,
    setTargetMetal,
    setCoordinationMode,
    setFixLigandPosition,
    setMetalReplacementPreset,
    applyMetalReplacementPreset,
    // RASA conditioning
    selectedBuriedAtoms,
    selectedExposedAtoms,
    buriedOverridden,
    exposedOverridden,
    setBuriedAtoms,
    setExposedAtoms,
    clearRASAConditioning,
    // H-Bond conditioning
    selectedHBondAcceptors,
    selectedHBondDonors,
    hbondOverridden,
    setHBondAcceptors,
    setHBondDonors,
    removeHBondDonors,
    clearHBondConditioning,
    // Apply suggestions
    applyEnzymeSuggestions,
  } = useStore();

  // Required
  const [pdbContent, setPdbContent] = useState<string | null>(null);
  const [pdbFileName, setPdbFileName] = useState<string | null>(null);
  const [ligandCodes, setLigandCodes] = useState('');
  const [proteinLength, setProteinLength] = useState('150');

  // Catalytic residues input state (actual residues are in store: enzymeCatalyticResidues)
  const [newCatChain, setNewCatChain] = useState('A');
  const [newCatResidue, setNewCatResidue] = useState('');
  const [newCatName, setNewCatName] = useState('');

  // Fixed atoms on catalytic residues
  const [fixedAtomType, setFixedAtomType] = useState('BKBN');

  // Covalent bonds between protein and ligand
  const [covalentBonds, setCovalentBonds] = useState<CovalentBond[]>([]);
  const [showCovalentSection, setShowCovalentSection] = useState(false);
  // New covalent bond input state
  const [newBond, setNewBond] = useState<CovalentBond>({
    proteinChain: 'A',
    proteinRes: 0,
    proteinResName: '',
    proteinAtom: '',
    ligandChain: '',
    ligandRes: 1,
    ligandResName: '',
    ligandAtom: '',
  });

  // Options
  const [qualityPreset, setQualityPreset] = useState<QualityPreset>('Balanced');
  const [qualityParams, setQualityParams] = useState<QualityParams>(QUALITY_PRESETS.Balanced);
  const [isNonLoopy, setIsNonLoopy] = useState(true);
  const [numDesigns, setNumDesigns] = useState(1);
  const [seed, setSeed] = useState<string>('');

  // Conditioning section enables
  const [rasaEnabled, setRasaEnabled] = useState(false);
  const [hbondEnabled, setHbondEnabled] = useState(false);

  // Run enzyme analysis when PDB content changes
  useEffect(() => {
    if (pdbContent) {
      setEnzymeAnalysisLoading(true);
      try {
        const analysis = analyzeEnzymeStructure(pdbContent);
        setEnzymeAnalysis(analysis);

        // Auto-set source metal if detected
        if (analysis.metals.length > 0) {
          setSourceMetal(analysis.metals[0].element);
        }

        // Auto-fill ligand codes from detected ligands and metals
        const detectedCodes: string[] = [];
        for (const ligand of analysis.ligands) {
          if (!detectedCodes.includes(ligand.name)) {
            detectedCodes.push(ligand.name);
          }
        }
        for (const metal of analysis.metals) {
          if (!detectedCodes.includes(metal.element)) {
            detectedCodes.push(metal.element);
          }
        }
        if (detectedCodes.length > 0) {
          setLigandCodes(detectedCodes.join(','));
        }
      } catch (error) {
        console.error('Enzyme analysis failed:', error);
        setEnzymeAnalysis(null);
      } finally {
        setEnzymeAnalysisLoading(false);
      }
    } else {
      setEnzymeAnalysis(null);
      setSourceMetal(null);
    }
  }, [pdbContent, setEnzymeAnalysis, setEnzymeAnalysisLoading, setSourceMetal]);

  // Get first detected ligand name for display
  const detectedLigandName = enzymeAnalysis?.ligands[0]?.name || null;

  const handleQualityChange = (preset: QualityPreset, params: QualityParams) => {
    setQualityPreset(preset);
    setQualityParams(params);
  };

  const addCatalyticResidue = () => {
    if (newCatResidue) {
      const residueNum = parseInt(newCatResidue, 10);
      if (!isNaN(residueNum)) {
        addEnzymeCatalyticResidue(newCatChain, residueNum, newCatName, fixedAtomType);
        setNewCatResidue('');
        setNewCatName('');
      }
    }
  };

  const removeCatalyticResidue = (chain: string, residue: number) => {
    removeEnzymeCatalyticResidue(chain, residue);
  };

  // Apply all catalytic suggestions
  const applyCatalyticSuggestions = useCallback(() => {
    if (catalyticSuggestions.length === 0) return;

    const newResidues = catalyticSuggestions.map((s) => ({
      chain: s.chain,
      residue: s.residue,
      name: s.name || '',
      atomType: fixedAtomType,
    }));
    setEnzymeCatalyticResidues(newResidues);
  }, [catalyticSuggestions, fixedAtomType, setEnzymeCatalyticResidues]);

  // Covalent bond management
  const addCovalentBond = () => {
    if (
      newBond.proteinRes > 0 &&
      newBond.proteinResName &&
      newBond.proteinAtom &&
      newBond.ligandResName &&
      newBond.ligandAtom
    ) {
      setCovalentBonds([...covalentBonds, { ...newBond }]);
      // Reset form
      setNewBond({
        proteinChain: 'A',
        proteinRes: 0,
        proteinResName: '',
        proteinAtom: '',
        ligandChain: '',
        ligandRes: 1,
        ligandResName: '',
        ligandAtom: '',
      });
    }
  };

  const removeCovalentBond = (index: number) => {
    setCovalentBonds(covalentBonds.filter((_, i) => i !== index));
  };

  // Format covalent bond for display
  const formatBondString = (bond: CovalentBond): string => {
    const proteinStr = `${bond.proteinChain}/${bond.proteinResName}/${bond.proteinRes}/${bond.proteinAtom}`;
    const ligandStr = `${bond.ligandChain || 'L'}/${bond.ligandResName}/${bond.ligandRes}/${bond.ligandAtom}`;
    return `${proteinStr} ↔ ${ligandStr}`;
  };

  const handleSubmit = async () => {
    // Build unindex string: comma-separated residue positions
    // Format: Chain:Residue[ResName] -> simplified to just "chain residue"
    // IMPORTANT: Filter out excluded residues (metal coordinators in expand/lanthanide presets)
    const filteredResidues = enzymeCatalyticResidues.filter((r) => {
      const key = `${r.chain}${r.residue}`;
      return !excludedCatalyticResidues.has(key);
    });
    const unindexParts = filteredResidues.map((r) => `${r.chain}${r.residue}`);
    const unindex = unindexParts.join(',');

    // Build fixed atoms selection for catalytic residues (use store's atom types)
    // Also filter out excluded residues from fixed atoms
    const selectFixedAtoms: Record<string, string> = {};
    for (const [key, value] of Object.entries(enzymeFixedAtomTypes)) {
      // Check if this residue is excluded (format: "A123")
      if (!excludedCatalyticResidues.has(key)) {
        selectFixedAtoms[key] = value;
      }
    }

    // Handle metal replacement - modify PDB content if enabled
    let finalPdbContent = pdbContent;
    let finalLigandCodes = ligandCodes.toUpperCase();

    // DEBUG: Log metal replacement state
    console.log('[EnzymeForm] Metal replacement state:', {
      metalReplacementEnabled,
      sourceMetal,
      targetMetal,
      hasPdbContent: !!pdbContent,
      originalLigandCodes: ligandCodes,
    });

    if (metalReplacementEnabled && sourceMetal && targetMetal && pdbContent) {
      // Replace metal in PDB
      finalPdbContent = replaceMetalInPdb(pdbContent, sourceMetal, targetMetal);

      // Update ligand codes to use target metal instead of source
      const codes = finalLigandCodes.split(',').map(c => c.trim());
      const updatedCodes = codes.map(code =>
        code === sourceMetal.toUpperCase() ? targetMetal.toUpperCase() : code
      );
      // Ensure target metal is in the list
      if (!updatedCodes.includes(targetMetal.toUpperCase())) {
        updatedCodes.unshift(targetMetal.toUpperCase());
      }
      finalLigandCodes = updatedCodes.join(',');

      // Fix ligand position if requested - add ligand to select_fixed_atoms
      // NOTE: Use ligand NAME (e.g., "CIT") as key, not chain+residue,
      // because RFD3 parser separates polymer/non-polymer and renames chains
      if (fixLigandPosition && enzymeAnalysis && enzymeAnalysis.ligands.length > 0) {
        const ligand = enzymeAnalysis.ligands[0];
        selectFixedAtoms[ligand.name] = 'ALL';
      }
    }

    const request: RFD3Request = {
      pdb_content: finalPdbContent || undefined,
      ligand: finalLigandCodes,
      length: proteinLength,
      unindex,
      num_designs: numDesigns,
      is_non_loopy: isNonLoopy,
      num_timesteps: qualityParams.num_timesteps,
      step_scale: qualityParams.step_scale,
      gamma_0: qualityParams.gamma_0,
    };

    // Add fixed atoms if any are specified (catalytic residues or fixed ligand)
    if (Object.keys(selectFixedAtoms).length > 0) {
      request.select_fixed_atoms = selectFixedAtoms;
    }

    if (seed) {
      request.seed = parseInt(seed, 10);
    }

    // Add covalent bonds if specified
    if (covalentBonds.length > 0) {
      request.covalent_bonds = covalentBonds.map((bond) => ({
        protein: {
          chain: bond.proteinChain,
          res_name: bond.proteinResName,
          res_num: bond.proteinRes,
          atom_name: bond.proteinAtom,
        },
        ligand: {
          chain: bond.ligandChain || 'L',
          res_name: bond.ligandResName,
          res_num: bond.ligandRes,
          atom_name: bond.ligandAtom,
        },
      }));
    }

    // Add RASA conditioning if enabled
    // IMPORTANT: Always bury the target metal when metal replacement is enabled
    const finalBuriedAtoms = { ...selectedBuriedAtoms };
    if (metalReplacementEnabled && targetMetal) {
      // Ensure target metal is buried
      finalBuriedAtoms[targetMetal] = 'ALL';
    }

    if (rasaEnabled || (metalReplacementEnabled && targetMetal)) {
      const buried = buildBuriedSelection(finalBuriedAtoms);
      const exposed = buildExposedSelection(selectedExposedAtoms);

      if (buried) {
        request.select_buried = buried;
      }
      if (exposed) {
        request.select_exposed = exposed;
      }
    }

    // Add H-bond conditioning if enabled
    if (hbondEnabled) {
      const acceptors = buildHBondAcceptorSelection(selectedHBondAcceptors);
      const donors = buildHBondDonorSelection(selectedHBondDonors);

      if (acceptors) {
        request.select_hbond_acceptor = acceptors;
      }
      if (donors) {
        request.select_hbond_donor = donors;
      }
    }

    // DEBUG: Log final request
    console.log('[EnzymeForm] Final request:', JSON.stringify(request, null, 2));

    await onSubmit(request);
  };

  const isValid =
    pdbContent !== null &&
    ligandCodes.trim() !== '' &&
    enzymeCatalyticResidues.length > 0 &&
    proteinLength.trim() !== '';

  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="flex items-center gap-3 pb-4 border-b border-border">
        <div className="w-10 h-10 rounded-lg bg-muted flex items-center justify-center">
          <Beaker className="w-5 h-5 text-muted-foreground" />
        </div>
        <div>
          <h2 className="font-semibold text-foreground">Enzyme Scaffold Design</h2>
          <p className="text-sm text-muted-foreground">Build a protein scaffold around an active site (theozyme)</p>
        </div>
      </div>

      {/* Info Banner */}
      <div className="p-4 rounded-xl bg-muted/50 border border-border">
        <div className="flex items-start gap-3">
          <Info className="w-5 h-5 text-muted-foreground flex-shrink-0 mt-0.5" />
          <div className="text-sm text-foreground">
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
            className="w-full px-4 py-2.5 rounded-xl border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none transition-all font-mono uppercase"
          />
        </FormField>

        {/* Detected ligands display - commented out for cleaner UI
        {enzymeAnalysis && enzymeAnalysis.ligands.length > 0 && (
          <div className="mt-2 flex flex-wrap gap-2">
            <span className="text-xs text-muted-foreground">Detected:</span>
            {enzymeAnalysis.ligands.map((lig) => (
              <span
                key={`${lig.chain}_${lig.name}_${lig.residueNumber}`}
                className="px-2 py-0.5 text-xs font-mono bg-muted rounded"
              >
                {lig.name}
              </span>
            ))}
            {enzymeAnalysis.metals.length > 0 && enzymeAnalysis.metals.map((metal) => (
              <span
                key={`${metal.chain}_${metal.element}_${metal.residueNumber}`}
                className="px-2 py-0.5 text-xs font-mono bg-primary/10 text-primary rounded"
              >
                {metal.element} (CN{metal.coordinationNumber})
              </span>
            ))}
          </div>
        )}
        */}
      </FormSection>

      {/* Metal Replacement Section - shown when metal is detected */}
      {enzymeAnalysis && enzymeAnalysis.metals.length > 0 && (
        <MetalReplacementSection
          enabled={metalReplacementEnabled}
          onEnabledChange={setMetalReplacementEnabled}
          sourceMetal={sourceMetal}
          targetMetal={targetMetal}
          onTargetMetalChange={setTargetMetal}
          coordinationMode={coordinationMode}
          onCoordinationModeChange={setCoordinationMode}
          fixLigandPosition={fixLigandPosition}
          onFixLigandPositionChange={setFixLigandPosition}
          preset={metalReplacementPreset}
          onPresetChange={applyMetalReplacementPreset}
          excludedResidues={excludedCatalyticResidues}
          enzymeAnalysis={enzymeAnalysis}
          ligandName={detectedLigandName || undefined}
        />
      )}

      {/* Catalytic Residues - Required */}
      <FormSection
        title="Catalytic Residues"
        description="Residues that form the active site"
        required
      >
        <div className="space-y-3">
          {/* Current residues as inline chips */}
          <div className="flex flex-wrap gap-1.5">
            {enzymeCatalyticResidues.map((res) => (
              <span
                key={`${res.chain}-${res.residue}`}
                className="inline-flex items-center gap-1 px-2 py-1 rounded-full bg-primary/10 text-primary text-xs font-medium"
              >
                {res.chain}{res.residue}{res.name ? `:${res.name}` : ''}
                <button
                  onClick={() => removeCatalyticResidue(res.chain, res.residue)}
                  className="hover:text-red-600 transition-colors"
                >
                  <X className="w-3 h-3" />
                </button>
              </span>
            ))}
            {enzymeCatalyticResidues.length === 0 && (
              <span className="text-xs text-muted-foreground italic">No residues added</span>
            )}
          </div>

          {/* Smart input + suggestions */}
          <div className="flex gap-2 items-center">
            <div className="flex-1 relative">
              <input
                type="text"
                value={newCatResidue ? `${newCatChain}${newCatResidue}${newCatName ? ':' + newCatName : ''}` : ''}
                onChange={(e) => {
                  const val = e.target.value.toUpperCase();
                  // Parse format: A145 or A145:CYS
                  const match = val.match(/^([A-Z]?)(\d*)(?::([A-Z]{0,3}))?$/);
                  if (match) {
                    setNewCatChain(match[1] || 'A');
                    setNewCatResidue(match[2] || '');
                    setNewCatName(match[3] || '');
                  }
                }}
                onKeyDown={(e) => e.key === 'Enter' && newCatResidue && addCatalyticResidue()}
                placeholder="A145:CYS"
                className="w-full px-3 py-2 rounded-lg border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none transition-all text-sm font-mono"
              />
            </div>
            <button
              onClick={addCatalyticResidue}
              disabled={!newCatResidue}
              className={`px-3 py-2 rounded-lg font-medium text-sm transition-all ${
                newCatResidue
                  ? 'bg-primary text-primary-foreground hover:bg-primary/90'
                  : 'bg-muted text-muted-foreground cursor-not-allowed'
              }`}
            >
              <Plus className="w-4 h-4" />
            </button>
            {catalyticSuggestions.length > 0 && (
              <>
                <button
                  onClick={applyCatalyticSuggestions}
                  className="inline-flex items-center gap-1 px-2 py-1.5 rounded-lg bg-primary/10 text-primary text-[10px] font-medium hover:bg-primary/20 transition-colors"
                  title="Apply all suggestions"
                >
                  <Sparkles className="w-3 h-3" />
                  Auto
                </button>
                <button
                  onClick={() => setBottomPanelMode('suggestions')}
                  className="inline-flex items-center justify-center w-7 h-7 rounded-lg bg-muted text-muted-foreground text-xs font-medium hover:bg-muted/80 transition-colors"
                  title="View suggestions"
                >
                  {catalyticSuggestions.length}
                </button>
              </>
            )}
          </div>

          {enzymeCatalyticResidues.length === 0 && (
            <p className="text-xs text-amber-600 flex items-center gap-1">
              <AlertTriangle className="w-3.5 h-3.5" />
              Required for enzyme scaffolding
            </p>
          )}
        </div>
      </FormSection>

      {/* Fixed Atom Selection - Horizontal segmented control */}
      <FormSection
        title="Fixed Atoms"
        description="Which atoms stay fixed during design"
      >
        <div className="inline-flex rounded-lg border border-border p-1 bg-muted/30">
          {[
            { value: 'BKBN', label: 'Backbone' },
            { value: 'ALL', label: 'All' },
            { value: 'TIP', label: 'Tips' },
            { value: '', label: 'None' },
          ].map((opt, i, arr) => (
            <button
              key={opt.value}
              onClick={() => setFixedAtomType(opt.value)}
              className={`px-4 py-1.5 text-sm font-medium transition-all ${
                fixedAtomType === opt.value
                  ? 'bg-primary text-primary-foreground rounded-md shadow-sm'
                  : 'text-muted-foreground hover:text-foreground'
              } ${i === 0 ? 'rounded-l-md' : ''} ${i === arr.length - 1 ? 'rounded-r-md' : ''}`}
            >
              {opt.label}
            </button>
          ))}
        </div>
      </FormSection>

      {/* Scaffold Length - Required */}
      <FormSection
        title="Scaffold Length"
        description="Total length of the enzyme scaffold"
        required
      >
        <div className="space-y-3">
          {/* Quick presets */}
          <div className="flex flex-wrap gap-1.5">
            {[
              { label: 'Small', value: '100-120' },
              { label: 'Medium', value: '150-180' },
              { label: 'Large', value: '200-250' },
              { label: 'TIM barrel', value: '220-260' },
            ].map((preset) => (
              <button
                key={preset.label}
                type="button"
                onClick={() => setProteinLength(preset.value)}
                className={`px-2.5 py-1 text-xs rounded-full transition-colors ${
                  proteinLength === preset.value
                    ? 'bg-primary text-primary-foreground'
                    : 'bg-muted text-muted-foreground hover:bg-muted/80'
                }`}
              >
                {preset.label}
              </button>
            ))}
          </div>
          <div className="flex gap-4 items-start">
            <div className="flex-1">
              <LengthRangeInput
                value={proteinLength}
                onChange={setProteinLength}
                label="Length"
                placeholder="150 or 120-180"
                hint="Single number or range"
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

      {/* Advanced Options */}
      <AdvancedOptionsWrapper title="Advanced Options">
        {/* Active Site Conditioning Group */}
        {enzymeAnalysis && (enzymeAnalysis.ligands.length > 0 || enzymeAnalysis.metals.length > 0) && (
          <div className="space-y-4">
            <div className="text-xs font-medium text-muted-foreground uppercase tracking-wide">Active Site Conditioning</div>

            <RASAConditioningSection
              enabled={rasaEnabled}
              onEnabledChange={setRasaEnabled}
              buriedAtoms={selectedBuriedAtoms}
              onBuriedAtomsChange={setBuriedAtoms}
              exposedAtoms={selectedExposedAtoms}
              onExposedAtomsChange={setExposedAtoms}
              buriedOverridden={buriedOverridden}
              exposedOverridden={exposedOverridden}
              enzymeAnalysis={enzymeAnalysis}
              targetMetal={metalReplacementEnabled ? targetMetal : sourceMetal}
              onApplySuggestions={applyEnzymeSuggestions}
            />

            {enzymeAnalysis.ligands.length > 0 && (
              <HBondConditioningSection
                enabled={hbondEnabled}
                onEnabledChange={setHbondEnabled}
                acceptorAtoms={selectedHBondAcceptors}
                onAcceptorAtomsChange={setHBondAcceptors}
                donorAtoms={selectedHBondDonors}
                onDonorAtomsChange={setHBondDonors}
                onRemoveDonor={removeHBondDonors}
                acceptorOverridden={hbondOverridden}
                enzymeAnalysis={enzymeAnalysis}
                onApplySuggestions={applyEnzymeSuggestions}
              />
            )}
          </div>
        )}

        {/* Structure Options Group */}
        <div className="space-y-3">
          <div className="text-xs font-medium text-muted-foreground uppercase tracking-wide">Structure Options</div>

          <label className="flex items-center gap-3 cursor-pointer">
            <input
              type="checkbox"
              checked={isNonLoopy}
              onChange={(e) => setIsNonLoopy(e.target.checked)}
              className="w-4 h-4 rounded border-border text-primary focus:ring-primary"
            />
            <div>
              <span className="text-sm text-foreground">Non-loopy mode</span>
              <span className="text-xs text-muted-foreground ml-2">(cleaner secondary structures)</span>
            </div>
          </label>

          <label className="flex items-center gap-3 cursor-pointer">
            <input
              type="checkbox"
              checked={showCovalentSection}
              onChange={(e) => setShowCovalentSection(e.target.checked)}
              className="w-4 h-4 rounded border-border text-primary focus:ring-primary"
            />
            <div>
              <span className="text-sm text-foreground">Covalent modifications</span>
              <span className="text-xs text-muted-foreground ml-2">(for covalent inhibitors)</span>
            </div>
          </label>
        </div>

        {/* Covalent Bonds Form - shown when enabled */}
        {showCovalentSection && (
          <div className="p-4 rounded-xl bg-muted/30 border border-border space-y-4">
            <div className="text-sm font-medium text-foreground">Covalent Bonds</div>

            {/* Existing bonds */}
            {covalentBonds.length > 0 && (
              <div className="space-y-1.5">
                {covalentBonds.map((bond, i) => (
                  <div
                    key={i}
                    className="flex items-center justify-between px-2 py-1.5 rounded-lg bg-card border border-border"
                  >
                    <code className="text-xs font-mono text-foreground">
                      {formatBondString(bond)}
                    </code>
                    <button
                      onClick={() => removeCovalentBond(i)}
                      className="text-muted-foreground hover:text-red-600 transition-colors"
                    >
                      <X className="w-3.5 h-3.5" />
                    </button>
                  </div>
                ))}
              </div>
            )}

            {/* Add new bond - compact form */}
            <div className="grid grid-cols-2 gap-3">
              <div className="space-y-1.5">
                <div className="text-xs text-muted-foreground">Protein (chain/res/name/atom)</div>
                <div className="flex gap-1">
                  <input
                    type="text"
                    value={newBond.proteinChain}
                    onChange={(e) => setNewBond({ ...newBond, proteinChain: e.target.value.toUpperCase() })}
                    placeholder="A"
                    maxLength={1}
                    className="w-10 px-2 py-1.5 rounded border border-border text-xs text-center"
                  />
                  <input
                    type="number"
                    value={newBond.proteinRes || ''}
                    onChange={(e) => setNewBond({ ...newBond, proteinRes: parseInt(e.target.value) || 0 })}
                    placeholder="#"
                    className="w-14 px-2 py-1.5 rounded border border-border text-xs"
                  />
                  <input
                    type="text"
                    value={newBond.proteinResName}
                    onChange={(e) => setNewBond({ ...newBond, proteinResName: e.target.value.toUpperCase() })}
                    placeholder="CYS"
                    maxLength={3}
                    className="w-12 px-2 py-1.5 rounded border border-border text-xs font-mono"
                  />
                  <input
                    type="text"
                    value={newBond.proteinAtom}
                    onChange={(e) => setNewBond({ ...newBond, proteinAtom: e.target.value.toUpperCase() })}
                    placeholder="SG"
                    maxLength={4}
                    className="w-12 px-2 py-1.5 rounded border border-border text-xs font-mono"
                  />
                </div>
              </div>
              <div className="space-y-1.5">
                <div className="text-xs text-muted-foreground">Ligand (name/res/atom)</div>
                <div className="flex gap-1">
                  <input
                    type="text"
                    value={newBond.ligandResName}
                    onChange={(e) => setNewBond({ ...newBond, ligandResName: e.target.value.toUpperCase() })}
                    placeholder="LIG"
                    className="w-14 px-2 py-1.5 rounded border border-border text-xs font-mono"
                  />
                  <input
                    type="number"
                    value={newBond.ligandRes || ''}
                    onChange={(e) => setNewBond({ ...newBond, ligandRes: parseInt(e.target.value) || 1 })}
                    placeholder="#"
                    className="w-12 px-2 py-1.5 rounded border border-border text-xs"
                  />
                  <input
                    type="text"
                    value={newBond.ligandAtom}
                    onChange={(e) => setNewBond({ ...newBond, ligandAtom: e.target.value.toUpperCase() })}
                    placeholder="C1"
                    maxLength={4}
                    className="w-12 px-2 py-1.5 rounded border border-border text-xs font-mono"
                  />
                  <button
                    onClick={addCovalentBond}
                    disabled={!newBond.proteinRes || !newBond.proteinResName || !newBond.proteinAtom || !newBond.ligandResName || !newBond.ligandAtom}
                    className="px-2 py-1.5 rounded bg-primary text-primary-foreground text-xs disabled:bg-muted disabled:text-muted-foreground"
                  >
                    <Plus className="w-3.5 h-3.5" />
                  </button>
                </div>
              </div>
            </div>
          </div>
        )}

        {/* Other Options */}
        <div className="space-y-3">
          <div className="text-xs font-medium text-muted-foreground uppercase tracking-wide">Other</div>
          <FormField label="Random Seed" hint="For reproducible results">
            <input
              type="number"
              value={seed}
              onChange={(e) => setSeed(e.target.value)}
              placeholder="Optional"
              className="w-full px-3 py-2 rounded-lg border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none transition-all text-sm"
            />
          </FormField>
        </div>
      </AdvancedOptionsWrapper>

      {/* Submit Button */}
      <div className="pt-4 border-t border-border space-y-2">
        <button
          onClick={handleSubmit}
          disabled={!isValid || isSubmitting || !health}
          className={`w-full py-3 px-6 rounded-xl font-semibold transition-all flex items-center justify-center gap-2 ${
            isValid && !isSubmitting && !!health
              ? 'bg-primary text-primary-foreground hover:bg-primary/90 shadow-lg shadow-primary/20'
              : 'bg-muted text-muted-foreground cursor-not-allowed'
          }`}
        >
          {isSubmitting ? (
            <>
              <Loader2 className="w-5 h-5 animate-spin" />
              Designing...
            </>
          ) : (
            <>
              <Rocket className="w-5 h-5" />
              Design {numDesigns} Scaffold{numDesigns > 1 ? 's' : ''}
            </>
          )}
        </button>
        {/* Missing requirements hint */}
        {!isValid && !isSubmitting && (
          <p className="text-center text-xs text-muted-foreground">
            {!pdbContent && 'Upload PDB'}
            {pdbContent && !ligandCodes.trim() && 'Add ligand codes'}
            {pdbContent && ligandCodes.trim() && enzymeCatalyticResidues.length === 0 && 'Add catalytic residues'}
            {pdbContent && ligandCodes.trim() && enzymeCatalyticResidues.length > 0 && !proteinLength.trim() && 'Set scaffold length'}
          </p>
        )}
        {!health && (
          <p className="text-center text-xs text-amber-600">
            Backend unavailable
          </p>
        )}
      </div>
    </div>
  );
}
