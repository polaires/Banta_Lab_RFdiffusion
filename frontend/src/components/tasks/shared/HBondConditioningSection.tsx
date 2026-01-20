'use client';

import { useState, useMemo, useCallback } from 'react';
import { AlertTriangle, Plus, X } from 'lucide-react';
import { FormSection } from './FormSection';
import { AtomCheckboxGrid } from './AtomCheckboxGrid';
import { suggestProteinDonors } from '../../../lib/enzymeAnalysis';
import type { EnzymeAnalysisResult } from '../../../lib/store';

interface HBondConditioningSectionProps {
  /** Whether H-bond conditioning is enabled */
  enabled: boolean;
  /** Toggle H-bond conditioning on/off */
  onEnabledChange: (enabled: boolean) => void;
  /** Selected H-bond acceptor atoms: { ligandName: "O1,O2,O3" } */
  acceptorAtoms: Record<string, string>;
  /** Callback when acceptor atoms change */
  onAcceptorAtomsChange: (key: string, atoms: string) => void;
  /** Selected H-bond donor atoms: { "A193": "N", "A195": "N" } */
  donorAtoms: Record<string, string>;
  /** Callback when donor atoms change */
  onDonorAtomsChange: (key: string, atoms: string) => void;
  /** Remove a donor entry */
  onRemoveDonor: (key: string) => void;
  /** Whether acceptor selection has been manually overridden */
  acceptorOverridden: boolean;
  /** Enzyme analysis result for ligand info */
  enzymeAnalysis: EnzymeAnalysisResult | null;
  /** Apply auto suggestions */
  onApplySuggestions: () => void;
  /** Additional CSS classes */
  className?: string;
}

// Quick-add presets for common H-bond donors
const QUICK_ADD_PRESETS = [
  { label: 'Backbone NH', atom: 'N', residues: ['ALA', 'GLY', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP'] },
  { label: 'Ser OG', atom: 'OG', residues: ['SER'] },
  { label: 'Thr OG1', atom: 'OG1', residues: ['THR'] },
  { label: 'Tyr OH', atom: 'OH', residues: ['TYR'] },
  { label: 'His NE2', atom: 'NE2', residues: ['HIS'] },
  { label: 'Asn ND2', atom: 'ND2', residues: ['ASN'] },
  { label: 'Gln NE2', atom: 'NE2', residues: ['GLN'] },
  { label: 'Arg NH1', atom: 'NH1', residues: ['ARG'] },
  { label: 'Trp NE1', atom: 'NE1', residues: ['TRP'] },
];

export function HBondConditioningSection({
  enabled,
  onEnabledChange,
  acceptorAtoms,
  onAcceptorAtomsChange,
  donorAtoms,
  onDonorAtomsChange,
  onRemoveDonor,
  acceptorOverridden,
  enzymeAnalysis,
  onApplySuggestions,
  className = '',
}: HBondConditioningSectionProps) {
  // Local state for adding new donors
  const [newDonorChain, setNewDonorChain] = useState('A');
  const [newDonorResidue, setNewDonorResidue] = useState('');
  const [newDonorAtom, setNewDonorAtom] = useState('N');
  const [addError, setAddError] = useState<string | null>(null);

  // Get ligand info
  const ligandInfo = useMemo(() => {
    if (!enzymeAnalysis || enzymeAnalysis.ligands.length === 0) return null;
    return enzymeAnalysis.ligands[0];
  }, [enzymeAnalysis]);

  // Convert comma-separated string to array
  const parseAtomString = (str: string): string[] => {
    return str ? str.split(',').map(s => s.trim()).filter(Boolean) : [];
  };

  // Convert array to comma-separated string
  const toAtomString = (atoms: string[]): string => {
    return atoms.join(',');
  };

  // Get acceptor atoms for ligand as array
  const ligandAcceptorAtoms = useMemo(() => {
    if (!ligandInfo) return [];
    return parseAtomString(acceptorAtoms[ligandInfo.name] || '');
  }, [ligandInfo, acceptorAtoms]);

  // Handle acceptor atom selection change
  const handleAcceptorChange = useCallback((atoms: string[]) => {
    if (!ligandInfo) return;
    onAcceptorAtomsChange(ligandInfo.name, toAtomString(atoms));
  }, [ligandInfo, onAcceptorAtomsChange]);

  // Apply acceptor suggestions
  const applyAcceptorSuggestions = useCallback(() => {
    if (!ligandInfo) return;
    onAcceptorAtomsChange(ligandInfo.name, toAtomString(ligandInfo.suggestedHBondAcceptors));
  }, [ligandInfo, onAcceptorAtomsChange]);

  // Add a new protein donor
  const handleAddDonor = useCallback(() => {
    const residueNum = parseInt(newDonorResidue);
    if (isNaN(residueNum) || residueNum < 1) {
      setAddError('Invalid residue number');
      return;
    }

    const key = `${newDonorChain}${residueNum}`;

    // Check if already exists
    if (donorAtoms[key]) {
      // Append atom if different
      const existingAtoms = parseAtomString(donorAtoms[key]);
      if (!existingAtoms.includes(newDonorAtom)) {
        onDonorAtomsChange(key, toAtomString([...existingAtoms, newDonorAtom]));
      }
    } else {
      onDonorAtomsChange(key, newDonorAtom);
    }

    setNewDonorResidue('');
    setAddError(null);
  }, [newDonorChain, newDonorResidue, newDonorAtom, donorAtoms, onDonorAtomsChange]);

  // Get donor entries as array
  const donorEntries = useMemo(() => {
    return Object.entries(donorAtoms)
      .filter(([, atoms]) => atoms && atoms.trim())
      .map(([key, atoms]) => ({ key, atoms }));
  }, [donorAtoms]);

  // Count total acceptor selections
  const totalAcceptors = useMemo(() => {
    return Object.values(acceptorAtoms).reduce((sum, atoms) => {
      return sum + parseAtomString(atoms).length;
    }, 0);
  }, [acceptorAtoms]);

  if (!ligandInfo) {
    return null; // Don't render if no ligand detected
  }

  return (
    <FormSection
      title="H-Bond Conditioning"
      description="Specify hydrogen bond donors and acceptors for transition state stabilization"
      className={className}
    >
      {/* Enable toggle */}
      <label className="flex items-center gap-2 cursor-pointer">
        <input
          type="checkbox"
          checked={enabled}
          onChange={(e) => onEnabledChange(e.target.checked)}
          className="w-4 h-4 rounded border-border text-primary focus:ring-primary"
        />
        <span className="text-sm text-foreground">Enable H-bond conditioning</span>
      </label>

      {/* Success rate warning - commented out for cleaner UI
      {enabled && (
        <div className="flex items-start gap-2 p-3 bg-yellow-50 text-yellow-800 dark:bg-yellow-950/30 dark:text-yellow-300 rounded-lg mt-3">
          <AlertTriangle className="h-4 w-4 mt-0.5 flex-shrink-0" />
          <p className="text-xs">
            <strong>Expected success rate: ~37% per H-bond.</strong> Not all requested H-bonds
            will form. Generate multiple designs and filter for best results.
          </p>
        </div>
      )}
      */}

      {enabled && (
        <div className="space-y-4 mt-3">
          {/* Ligand Acceptors Section */}
          <div className="p-3 bg-muted/30 rounded-lg border border-border space-y-3">
            <h4 className="text-sm font-medium text-foreground">
              Ligand Acceptors
              <span className="text-xs text-muted-foreground ml-2">
                (receive H-bonds from protein)
              </span>
            </h4>

            <AtomCheckboxGrid
              label={`${ligandInfo.name} atoms:`}
              atoms={ligandInfo.atoms.filter(a => ['O', 'N', 'S'].includes(a.element))}
              selectedAtoms={ligandAcceptorAtoms}
              onChange={handleAcceptorChange}
              suggestedAtoms={ligandInfo.suggestedHBondAcceptors}
              suggestionButtonLabel="Select all O/N atoms"
              onApplySuggestions={applyAcceptorSuggestions}
              isOverridden={acceptorOverridden}
              // hint="O and N atoms can accept H-bonds from protein donors" - commented out for cleaner UI
            />
          </div>

          {/* Protein Donors Section */}
          <div className="p-3 bg-muted/30 rounded-lg border border-border space-y-3">
            <h4 className="text-sm font-medium text-foreground">
              Protein Donors
              <span className="text-xs text-muted-foreground ml-2">
                (optional: specify specific residues as donors)
              </span>
            </h4>

            {/* Current donors */}
            {donorEntries.length > 0 && (
              <div className="flex flex-wrap gap-1.5">
                {donorEntries.map(({ key, atoms }) => (
                  <span
                    key={key}
                    className="inline-flex items-center gap-1 px-2 py-1 bg-blue-50 text-blue-700 dark:bg-blue-950/30 dark:text-blue-300 text-xs font-medium rounded-full border border-blue-200 dark:border-blue-800"
                  >
                    {key}:{atoms}
                    <button
                      type="button"
                      onClick={() => onRemoveDonor(key)}
                      className="hover:text-blue-500 transition-colors"
                    >
                      <X className="h-3.5 w-3.5" />
                    </button>
                  </span>
                ))}
              </div>
            )}

            {donorEntries.length === 0 && (
              <p className="text-xs text-muted-foreground italic">
                No protein donors specified. RFD3 will auto-design appropriate donors.
              </p>
            )}

            {/* Add new donor */}
            <div className="flex gap-2 items-end">
              <div className="w-16">
                <label className="block text-xs text-muted-foreground mb-1">Chain</label>
                <input
                  type="text"
                  value={newDonorChain}
                  onChange={(e) => setNewDonorChain(e.target.value.toUpperCase())}
                  maxLength={1}
                  className="w-full px-2 py-1.5 border border-border rounded text-sm bg-card focus:ring-2 focus:ring-ring"
                />
              </div>
              <div className="w-20">
                <label className="block text-xs text-muted-foreground mb-1">Residue</label>
                <input
                  type="text"
                  value={newDonorResidue}
                  onChange={(e) => {
                    setNewDonorResidue(e.target.value);
                    setAddError(null);
                  }}
                  placeholder="e.g., 193"
                  className="w-full px-2 py-1.5 border border-border rounded text-sm bg-card focus:ring-2 focus:ring-ring"
                />
              </div>
              <div className="w-24">
                <label className="block text-xs text-muted-foreground mb-1">Atom</label>
                <select
                  value={newDonorAtom}
                  onChange={(e) => setNewDonorAtom(e.target.value)}
                  className="w-full px-2 py-1.5 border border-border rounded text-sm bg-card focus:ring-2 focus:ring-ring"
                >
                  <option value="N">N (backbone)</option>
                  <option value="OG">OG (Ser)</option>
                  <option value="OG1">OG1 (Thr)</option>
                  <option value="OH">OH (Tyr)</option>
                  <option value="ND1">ND1 (His)</option>
                  <option value="NE2">NE2 (His/Gln)</option>
                  <option value="ND2">ND2 (Asn)</option>
                  <option value="NH1">NH1 (Arg)</option>
                  <option value="NH2">NH2 (Arg)</option>
                  <option value="NE1">NE1 (Trp)</option>
                  <option value="NZ">NZ (Lys)</option>
                </select>
              </div>
              <button
                type="button"
                onClick={handleAddDonor}
                disabled={!newDonorResidue}
                className="px-3 py-1.5 bg-primary text-primary-foreground text-sm font-medium rounded hover:bg-primary/90 disabled:bg-muted disabled:text-muted-foreground disabled:cursor-not-allowed transition-colors"
              >
                <Plus className="h-4 w-4" />
              </button>
            </div>

            {addError && (
              <p className="text-xs text-red-500">{addError}</p>
            )}

            {/* Quick-add presets */}
            <div className="pt-2 border-t border-border">
              <p className="text-xs text-muted-foreground mb-2">Quick add:</p>
              <div className="flex flex-wrap gap-1">
                {QUICK_ADD_PRESETS.slice(0, 6).map((preset) => (
                  <button
                    key={preset.label}
                    type="button"
                    onClick={() => {
                      setNewDonorAtom(preset.atom);
                    }}
                    className="px-2 py-1 text-[10px] bg-muted text-muted-foreground rounded hover:bg-muted/80 transition-colors"
                  >
                    {preset.label}
                  </button>
                ))}
              </div>
            </div>
          </div>

          {/* Apply all suggestions button - removed, using inline Auto buttons instead */}

          {/* Summary - commented out for cleaner UI
          {totalAcceptors > 0 && (
            <div className="p-3 bg-blue-50 dark:bg-blue-950/30 rounded-lg border border-blue-200 dark:border-blue-800">
              <p className="text-xs text-blue-800 dark:text-blue-300">
                <strong>Summary:</strong> {totalAcceptors} acceptor atom{totalAcceptors !== 1 ? 's' : ''} selected
                {donorEntries.length > 0 && `, ${donorEntries.length} protein donor${donorEntries.length !== 1 ? 's' : ''} specified`}.
                Expected H-bonds formed: ~{Math.round(totalAcceptors * 0.37)} (at 37% success rate).
              </p>
            </div>
          )}
          */}
        </div>
      )}
    </FormSection>
  );
}

export default HBondConditioningSection;
