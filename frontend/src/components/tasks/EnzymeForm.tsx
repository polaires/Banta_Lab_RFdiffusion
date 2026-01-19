'use client';

import { useState } from 'react';
import { Beaker, Info, X, AlertTriangle, Plus, Loader2, Rocket, Sparkles } from 'lucide-react';
import { useStore } from '@/lib/store';
import { FormSection, FormField, FormRow } from './shared/FormSection';
import { PdbUploader } from './shared/PdbUploader';
import { LengthRangeInput } from './shared/LengthRangeInput';
import { QualityPresetSelector, QualityPreset, QualityParams } from './shared/QualityPresetSelector';
import { AdvancedOptionsWrapper } from './shared/AdvancedOptionsWrapper';
import { QUALITY_PRESETS, RFD3Request, TaskFormProps } from './shared/types';

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
    const unindexParts = enzymeCatalyticResidues.map((r) => `${r.chain}${r.residue}`);
    const unindex = unindexParts.join(',');

    // Build fixed atoms selection for catalytic residues (use store's atom types)
    const selectFixedAtoms: Record<string, string> = { ...enzymeFixedAtomTypes };

    const request: RFD3Request = {
      pdb_content: pdbContent || undefined,
      ligand: ligandCodes.toUpperCase(),
      length: proteinLength,
      unindex,
      num_designs: numDesigns,
      is_non_loopy: isNonLoopy,
      num_timesteps: qualityParams.num_timesteps,
      step_scale: qualityParams.step_scale,
      gamma_0: qualityParams.gamma_0,
    };

    if (Object.keys(selectFixedAtoms).length > 0 && fixedAtomType) {
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
      </FormSection>

      {/* Catalytic Residues - Required */}
      <FormSection
        title="Catalytic Residues"
        description="Residues that form the active site. Their positions will be inferred by the model (unindex)."
        required
      >
        <div className="space-y-4">
          {/* Current residues */}
          {enzymeCatalyticResidues.length > 0 && (
            <div className="flex flex-wrap gap-2 mb-3">
              {enzymeCatalyticResidues.map((res) => (
                <div
                  key={`${res.chain}-${res.residue}`}
                  className="flex items-center gap-1 px-3 py-1.5 rounded-lg bg-muted text-foreground text-sm font-medium"
                >
                  <span>
                    {res.chain}:{res.residue}
                    {res.name && <span className="text-muted-foreground ml-1">({res.name})</span>}
                  </span>
                  <button
                    onClick={() => removeCatalyticResidue(res.chain, res.residue)}
                    className="ml-1 hover:text-red-600 transition-colors"
                  >
                    <X className="w-4 h-4" />
                  </button>
                </div>
              ))}
            </div>
          )}

          {/* Show Suggestions button */}
          {catalyticSuggestions.length > 0 && (
            <button
              onClick={() => setBottomPanelMode('suggestions')}
              className="flex items-center gap-1.5 px-3 py-1.5 rounded-lg bg-primary/10 text-primary text-sm font-medium hover:bg-primary/20 transition-colors"
            >
              <Sparkles className="w-4 h-4" />
              Show {catalyticSuggestions.length} Suggestions
              <span className="flex gap-0.5 ml-1">
                {catalyticSuggestions.some(s => s.source === 'mcsa') && (
                  <span className="w-2 h-2 rounded-full bg-blue-500" />
                )}
                {catalyticSuggestions.some(s => s.source === 'local') && (
                  <span className="w-2 h-2 rounded-full bg-orange-500" />
                )}
              </span>
            </button>
          )}

          {/* Add new residue */}
          <div className="flex gap-2">
            <input
              type="text"
              value={newCatChain}
              onChange={(e) => setNewCatChain(e.target.value.toUpperCase())}
              placeholder="Chain"
              maxLength={1}
              className="w-20 px-3 py-2 rounded-xl border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none transition-all text-center"
            />
            <input
              type="number"
              value={newCatResidue}
              onChange={(e) => setNewCatResidue(e.target.value)}
              placeholder="Residue #"
              onKeyDown={(e) => e.key === 'Enter' && addCatalyticResidue()}
              className="w-28 px-3 py-2 rounded-xl border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none transition-all"
            />
            <input
              type="text"
              value={newCatName}
              onChange={(e) => setNewCatName(e.target.value.toUpperCase())}
              placeholder="Name (opt)"
              maxLength={3}
              className="w-24 px-3 py-2 rounded-xl border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none transition-all font-mono"
            />
            <button
              onClick={addCatalyticResidue}
              disabled={!newCatResidue}
              className={`px-4 py-2 rounded-xl font-medium transition-all ${
                newCatResidue
                  ? 'bg-primary text-white hover:bg-primary/90'
                  : 'bg-muted text-muted-foreground cursor-not-allowed'
              }`}
            >
              Add
            </button>
          </div>

          {enzymeCatalyticResidues.length === 0 && (
            <p className="text-xs text-amber-600 flex items-center gap-1">
              <AlertTriangle className="w-4 h-4" />
              Add at least one catalytic residue for enzyme scaffolding
            </p>
          )}

          {/* Info about unindex */}
          <div className="p-3 rounded-xl bg-muted/50 border border-border text-sm text-muted-foreground">
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
                  ? 'border-primary bg-primary/10'
                  : 'border-border hover:border-border'
              }`}
            >
              <div className="font-medium text-sm text-foreground">{opt.label}</div>
              <div className="text-xs text-muted-foreground mt-0.5">{opt.desc}</div>
            </button>
          ))}
        </div>
      </FormSection>

      {/* Covalent Bonds Section */}
      <FormSection
        title="Covalent Modifications"
        description="Optional: Specify covalent bonds between protein and ligand (e.g., covalent inhibitors)"
      >
        <div className="space-y-4">
          {/* Toggle to show/hide covalent section */}
          <label className="flex items-center gap-3 p-3 rounded-xl bg-muted/50 hover:bg-muted cursor-pointer transition-colors">
            <input
              type="checkbox"
              checked={showCovalentSection}
              onChange={(e) => setShowCovalentSection(e.target.checked)}
              className="w-5 h-5 rounded border-border text-primary focus:ring-primary"
            />
            <div>
              <div className="font-medium text-sm text-foreground">Add Covalent Bonds</div>
              <div className="text-xs text-muted-foreground">
                For covalent inhibitors, suicide substrates, or covalently-linked cofactors
              </div>
            </div>
          </label>

          {showCovalentSection && (
            <div className="space-y-4 p-4 rounded-xl bg-muted/50 border border-border">
              {/* Existing covalent bonds */}
              {covalentBonds.length > 0 && (
                <div className="space-y-2">
                  <div className="text-sm font-medium text-foreground">Defined Bonds:</div>
                  {covalentBonds.map((bond, i) => (
                    <div
                      key={i}
                      className="flex items-center justify-between px-3 py-2 rounded-lg bg-card border border-border"
                    >
                      <code className="text-sm font-mono text-foreground">
                        {formatBondString(bond)}
                      </code>
                      <button
                        onClick={() => removeCovalentBond(i)}
                        className="text-muted-foreground hover:text-red-600 transition-colors"
                      >
                        <X className="w-4 h-4" />
                      </button>
                    </div>
                  ))}
                </div>
              )}

              {/* Add new covalent bond form */}
              <div className="space-y-3">
                <div className="text-sm font-medium text-foreground">Add New Bond:</div>

                {/* Protein side */}
                <div className="p-3 rounded-lg bg-card border border-border">
                  <div className="text-xs font-medium text-muted-foreground mb-2">Protein Residue</div>
                  <div className="flex gap-2 flex-wrap">
                    <input
                      type="text"
                      value={newBond.proteinChain}
                      onChange={(e) => setNewBond({ ...newBond, proteinChain: e.target.value.toUpperCase() })}
                      placeholder="Chain"
                      maxLength={1}
                      className="w-16 px-3 py-2 rounded-lg border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none text-sm text-center"
                    />
                    <input
                      type="text"
                      value={newBond.proteinResName}
                      onChange={(e) => setNewBond({ ...newBond, proteinResName: e.target.value.toUpperCase() })}
                      placeholder="ResName (e.g., CYS)"
                      maxLength={3}
                      className="w-28 px-3 py-2 rounded-lg border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none text-sm font-mono"
                    />
                    <input
                      type="number"
                      value={newBond.proteinRes || ''}
                      onChange={(e) => setNewBond({ ...newBond, proteinRes: parseInt(e.target.value) || 0 })}
                      placeholder="ResNum"
                      className="w-20 px-3 py-2 rounded-lg border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none text-sm"
                    />
                    <input
                      type="text"
                      value={newBond.proteinAtom}
                      onChange={(e) => setNewBond({ ...newBond, proteinAtom: e.target.value.toUpperCase() })}
                      placeholder="Atom (e.g., SG)"
                      maxLength={4}
                      className="w-24 px-3 py-2 rounded-lg border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none text-sm font-mono"
                    />
                  </div>
                </div>

                {/* Ligand side */}
                <div className="p-3 rounded-lg bg-card border border-border">
                  <div className="text-xs font-medium text-muted-foreground mb-2">Ligand Atom</div>
                  <div className="flex gap-2 flex-wrap">
                    <input
                      type="text"
                      value={newBond.ligandChain}
                      onChange={(e) => setNewBond({ ...newBond, ligandChain: e.target.value.toUpperCase() })}
                      placeholder="Chain (opt)"
                      maxLength={1}
                      className="w-20 px-2 py-1.5 rounded-lg border border-border focus:border-primary focus:ring-1 focus:ring-ring/20 outline-none text-sm text-center"
                    />
                    <input
                      type="text"
                      value={newBond.ligandResName}
                      onChange={(e) => setNewBond({ ...newBond, ligandResName: e.target.value.toUpperCase() })}
                      placeholder="LigandName (e.g., LIG)"
                      className="w-32 px-3 py-2 rounded-lg border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none text-sm font-mono"
                    />
                    <input
                      type="number"
                      value={newBond.ligandRes || ''}
                      onChange={(e) => setNewBond({ ...newBond, ligandRes: parseInt(e.target.value) || 1 })}
                      placeholder="ResNum"
                      className="w-20 px-3 py-2 rounded-lg border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none text-sm"
                    />
                    <input
                      type="text"
                      value={newBond.ligandAtom}
                      onChange={(e) => setNewBond({ ...newBond, ligandAtom: e.target.value.toUpperCase() })}
                      placeholder="Atom (e.g., C1)"
                      maxLength={4}
                      className="w-24 px-3 py-2 rounded-lg border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 outline-none text-sm font-mono"
                    />
                  </div>
                </div>

                {/* Add bond button */}
                <button
                  onClick={addCovalentBond}
                  disabled={
                    !newBond.proteinRes ||
                    !newBond.proteinResName ||
                    !newBond.proteinAtom ||
                    !newBond.ligandResName ||
                    !newBond.ligandAtom
                  }
                  className={`w-full py-2 px-4 rounded-lg font-medium text-sm transition-all ${
                    newBond.proteinRes &&
                    newBond.proteinResName &&
                    newBond.proteinAtom &&
                    newBond.ligandResName &&
                    newBond.ligandAtom
                      ? 'bg-primary text-white hover:bg-primary/90'
                      : 'bg-muted text-muted-foreground cursor-not-allowed'
                  }`}
                >
                  <span className="flex items-center justify-center gap-2">
                    <Plus className="w-4 h-4" />
                    Add Covalent Bond
                  </span>
                </button>
              </div>

              {/* Help text */}
              <div className="p-3 rounded-lg bg-primary/10 border border-primary/30 text-sm text-primary">
                <strong>Format:</strong> Chain/ResName/ResNum/AtomName (e.g., A/CYS/145/SG for cysteine sulfur)
              </div>
            </div>
          )}
        </div>
      </FormSection>

      {/* Scaffold Length - Required */}
      <FormSection
        title="Scaffold Length"
        description="Total length of the enzyme scaffold to design"
        required
      >
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
