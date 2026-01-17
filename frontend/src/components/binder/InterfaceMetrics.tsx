'use client';

import {
  Network,
  MousePointerClick,
  Link2,
  Layers,
  Box,
  FlaskConical,
  BadgeCheck,
  AlertTriangle,
  Brain,
  FileEdit,
} from 'lucide-react';

// PLIP interaction profile from ligand analysis
export interface InteractionProfile {
  hydrogen_bonds: number;
  hydrophobic_contacts: number;
  pi_stacking: number;
  salt_bridges: number;
  halogen_bonds?: number;
  total: number;
}

interface InterfaceMetricsProps {
  metrics: {
    interface_contacts: number;
    interface_hbonds: number;
    buried_sasa: number;
    packstat: number;
    // New BindCraft-style metrics
    shape_complementarity?: number;
    surface_hydrophobicity?: number;
    unsaturated_hbonds?: number;
    interface_residue_count?: number;
  };
  esm?: {
    perplexity: number;
    confidence: number;
  };
  mpnn_score?: number;
  // ESMFold validation metrics
  esmfold?: {
    plddt: number;
    rmsd?: number;
    validation_passed: boolean;
  };
  // PLIP interaction profile
  interactions?: InteractionProfile;
  key_binding_residues?: string[];
  analysis_method?: 'plip' | 'distance_based';
}

function getContactsQuality(contacts: number): { label: string; color: string } {
  if (contacts >= 5000) return { label: 'Excellent', color: 'text-green-600' };
  if (contacts >= 3000) return { label: 'Good', color: 'text-blue-600' };
  if (contacts >= 1000) return { label: 'Moderate', color: 'text-amber-600' };
  return { label: 'Low', color: 'text-red-600' };
}

function getHbondsQuality(hbonds: number): { label: string; color: string } {
  if (hbonds >= 15) return { label: 'Excellent', color: 'text-green-600' };
  if (hbonds >= 10) return { label: 'Good', color: 'text-blue-600' };
  if (hbonds >= 5) return { label: 'Moderate', color: 'text-amber-600' };
  return { label: 'Low', color: 'text-red-600' };
}

function getPerplexityQuality(perplexity: number): { label: string; color: string } {
  if (perplexity < 5.0) return { label: 'Excellent', color: 'text-green-600' };
  if (perplexity < 8.0) return { label: 'Good', color: 'text-blue-600' };
  if (perplexity < 12.0) return { label: 'Moderate', color: 'text-amber-600' };
  return { label: 'Poor', color: 'text-red-600' };
}

function getPackstatQuality(packstat: number): { label: string; color: string } {
  if (packstat >= 0.8) return { label: 'Excellent', color: 'text-green-600' };
  if (packstat >= 0.6) return { label: 'Good', color: 'text-blue-600' };
  if (packstat >= 0.4) return { label: 'Moderate', color: 'text-amber-600' };
  return { label: 'Poor', color: 'text-red-600' };
}

// New BindCraft-style quality functions
function getSCQuality(sc: number): { label: string; color: string } {
  if (sc >= 0.65) return { label: 'Excellent', color: 'text-green-600' };
  if (sc >= 0.5) return { label: 'Good', color: 'text-blue-600' };
  if (sc >= 0.4) return { label: 'Moderate', color: 'text-amber-600' };
  return { label: 'Poor', color: 'text-red-600' };
}

function getHydrophobicityQuality(hydrophobicity: number): { label: string; color: string } {
  if (hydrophobicity < 0.25) return { label: 'Excellent', color: 'text-green-600' };
  if (hydrophobicity < 0.37) return { label: 'Good', color: 'text-blue-600' };
  if (hydrophobicity < 0.50) return { label: 'Moderate', color: 'text-amber-600' };
  return { label: 'High', color: 'text-red-600' };
}

function getPlddtQuality(plddt: number): { label: string; color: string } {
  if (plddt >= 0.8) return { label: 'High Confidence', color: 'text-green-600' };
  if (plddt >= 0.7) return { label: 'Good', color: 'text-blue-600' };
  if (plddt >= 0.5) return { label: 'Low', color: 'text-amber-600' };
  return { label: 'Very Low', color: 'text-red-600' };
}

export function InterfaceMetrics({
  metrics,
  esm,
  mpnn_score,
  esmfold,
  interactions,
  key_binding_residues,
  analysis_method
}: InterfaceMetricsProps) {
  const contactsQ = getContactsQuality(metrics.interface_contacts);
  const hbondsQ = getHbondsQuality(metrics.interface_hbonds);
  const scQ = metrics.shape_complementarity != null ? getSCQuality(metrics.shape_complementarity) : null;
  const hydrophobicityQ = metrics.surface_hydrophobicity != null ? getHydrophobicityQuality(metrics.surface_hydrophobicity) : null;
  const plddtQ = esmfold ? getPlddtQuality(esmfold.plddt) : null;
  const packstatQ = getPackstatQuality(metrics.packstat);
  const perplexityQ = esm ? getPerplexityQuality(esm.perplexity) : null;

  return (
    <div className="bg-white rounded-xl border border-slate-200 overflow-hidden">
      {/* Header */}
      <div className="bg-gradient-to-r from-cyan-50 to-blue-50 px-4 py-3 border-b border-cyan-100">
        <div className="flex items-center gap-2">
          <Network className="h-5 w-5 text-cyan-600" />
          <h4 className="font-semibold text-slate-900 text-sm">Interface Quality</h4>
        </div>
        <p className="text-xs text-slate-500 mt-1">
          Protein-protein binding interface metrics
        </p>
      </div>

      {/* Metrics Grid */}
      <div className="p-4 grid grid-cols-2 gap-4">
        {/* Interface Contacts */}
        <div className="bg-slate-50 rounded-lg p-3">
          <div className="flex items-center gap-2 mb-2">
            <MousePointerClick className="h-4 w-4 text-slate-400" />
            <span className="text-xs text-slate-500 uppercase tracking-wide">Contacts</span>
          </div>
          <div className="text-2xl font-bold text-slate-900">
            {metrics.interface_contacts.toLocaleString()}
          </div>
          <div className={`text-xs font-medium ${contactsQ.color}`}>
            {contactsQ.label}
          </div>
        </div>

        {/* Hydrogen Bonds */}
        <div className="bg-slate-50 rounded-lg p-3">
          <div className="flex items-center gap-2 mb-2">
            <Link2 className="h-4 w-4 text-slate-400" />
            <span className="text-xs text-slate-500 uppercase tracking-wide">H-Bonds</span>
          </div>
          <div className="text-2xl font-bold text-slate-900">
            {metrics.interface_hbonds}
          </div>
          <div className={`text-xs font-medium ${hbondsQ.color}`}>
            {hbondsQ.label}
          </div>
        </div>

        {/* Buried SASA */}
        <div className="bg-slate-50 rounded-lg p-3">
          <div className="flex items-center gap-2 mb-2">
            <Layers className="h-4 w-4 text-slate-400" />
            <span className="text-xs text-slate-500 uppercase tracking-wide">Buried SASA</span>
          </div>
          <div className="text-2xl font-bold text-slate-900">
            {metrics.buried_sasa.toFixed(1)}
            <span className="text-sm font-normal text-slate-500 ml-1">Å²</span>
          </div>
          <div className="text-xs text-slate-500">
            Interface area
          </div>
        </div>

        {/* Packstat */}
        <div className="bg-slate-50 rounded-lg p-3">
          <div className="flex items-center gap-2 mb-2">
            <Box className="h-4 w-4 text-slate-400" />
            <span className="text-xs text-slate-500 uppercase tracking-wide">Packstat</span>
          </div>
          <div className="text-2xl font-bold text-slate-900">
            {metrics.packstat.toFixed(2)}
          </div>
          <div className={`text-xs font-medium ${packstatQ.color}`}>
            {packstatQ.label}
          </div>
        </div>
      </div>

      {/* BindCraft-style metrics (if available) */}
      {(metrics.shape_complementarity != null || metrics.surface_hydrophobicity != null) && (
        <div className="px-4 pb-4">
          <div className="bg-emerald-50 rounded-lg p-3 border border-emerald-100">
            <div className="flex items-center gap-2 mb-3">
              <FlaskConical className="h-4 w-4 text-emerald-600" />
              <span className="text-xs text-emerald-700 font-medium uppercase tracking-wide">
                BindCraft Metrics
              </span>
            </div>
            <div className="grid grid-cols-2 gap-4">
              {/* Shape Complementarity */}
              {metrics.shape_complementarity != null && (
                <div>
                  <div className="text-xs text-slate-500 mb-1">Shape Complementarity</div>
                  <div className="text-xl font-bold text-slate-900">
                    {metrics.shape_complementarity.toFixed(2)}
                  </div>
                  <div className={`text-xs font-medium ${scQ?.color}`}>
                    {scQ?.label}
                  </div>
                </div>
              )}
              {/* Surface Hydrophobicity */}
              {metrics.surface_hydrophobicity != null && (
                <div>
                  <div className="text-xs text-slate-500 mb-1">Surface Hydrophobicity</div>
                  <div className="text-xl font-bold text-slate-900">
                    {(metrics.surface_hydrophobicity * 100).toFixed(0)}%
                  </div>
                  <div className={`text-xs font-medium ${hydrophobicityQ?.color}`}>
                    {hydrophobicityQ?.label}
                  </div>
                </div>
              )}
              {/* Unsaturated H-bonds */}
              {metrics.unsaturated_hbonds != null && (
                <div>
                  <div className="text-xs text-slate-500 mb-1">Unsat. H-bonds</div>
                  <div className="text-xl font-bold text-slate-900">
                    {metrics.unsaturated_hbonds}
                  </div>
                  <div className={`text-xs font-medium ${metrics.unsaturated_hbonds <= 6 ? 'text-green-600' : 'text-amber-600'}`}>
                    {metrics.unsaturated_hbonds <= 6 ? 'Good' : 'High'}
                  </div>
                </div>
              )}
              {/* Interface Residue Count */}
              {metrics.interface_residue_count != null && (
                <div>
                  <div className="text-xs text-slate-500 mb-1">Interface Residues</div>
                  <div className="text-xl font-bold text-slate-900">
                    {metrics.interface_residue_count}
                  </div>
                  <div className={`text-xs font-medium ${metrics.interface_residue_count >= 6 ? 'text-green-600' : 'text-amber-600'}`}>
                    {metrics.interface_residue_count >= 6 ? 'Good' : 'Small'}
                  </div>
                </div>
              )}
            </div>
          </div>
        </div>
      )}

      {/* PLIP Interaction Profile (if available) */}
      {interactions && (
        <div className="px-4 pb-4">
          <div className="bg-gradient-to-br from-blue-50 to-purple-50 rounded-lg p-4 border border-blue-100">
            <div className="flex items-center justify-between mb-3">
              <div className="flex items-center gap-2">
                <FlaskConical className="h-4 w-4 text-blue-600" />
                <span className="text-xs text-blue-700 font-medium uppercase tracking-wide">
                  Interaction Profile
                </span>
              </div>
              {analysis_method && (
                <span className={`text-xs px-2 py-0.5 rounded-full ${
                  analysis_method === 'plip'
                    ? 'bg-purple-100 text-purple-700'
                    : 'bg-slate-100 text-slate-600'
                }`}>
                  {analysis_method === 'plip' ? 'PLIP Analysis' : 'Distance-Based'}
                </span>
              )}
            </div>

            {/* Interaction counts grid */}
            <div className="grid grid-cols-5 gap-2 mb-4">
              {/* H-bonds */}
              <div className="text-center p-2 bg-white/60 rounded-lg">
                <div className="text-xl font-bold text-blue-700">
                  {interactions.hydrogen_bonds}
                </div>
                <div className="text-xs text-slate-500">H-bonds</div>
              </div>

              {/* Hydrophobic */}
              <div className="text-center p-2 bg-white/60 rounded-lg">
                <div className="text-xl font-bold text-amber-700">
                  {interactions.hydrophobic_contacts}
                </div>
                <div className="text-xs text-slate-500">Hydrophobic</div>
              </div>

              {/* Pi-stacking */}
              <div className="text-center p-2 bg-white/60 rounded-lg">
                <div className="text-xl font-bold text-purple-700">
                  {interactions.pi_stacking}
                </div>
                <div className="text-xs text-slate-500">Pi-Stack</div>
              </div>

              {/* Salt bridges */}
              <div className="text-center p-2 bg-white/60 rounded-lg">
                <div className="text-xl font-bold text-red-700">
                  {interactions.salt_bridges}
                </div>
                <div className="text-xs text-slate-500">Salt Bridges</div>
              </div>

              {/* Total */}
              <div className="text-center p-2 bg-white/60 rounded-lg border-2 border-blue-200">
                <div className="text-xl font-bold text-slate-900">
                  {interactions.total}
                </div>
                <div className="text-xs text-slate-500 font-medium">Total</div>
              </div>
            </div>

            {/* Halogen bonds (if present) */}
            {interactions.halogen_bonds !== undefined && interactions.halogen_bonds > 0 && (
              <div className="mb-3 text-center p-2 bg-white/60 rounded-lg">
                <span className="text-sm text-teal-700">
                  <span className="font-bold">{interactions.halogen_bonds}</span> halogen bonds
                </span>
              </div>
            )}

            {/* Key binding residues */}
            {key_binding_residues && key_binding_residues.length > 0 && (
              <div>
                <div className="text-xs text-slate-500 mb-2">Key Binding Residues:</div>
                <div className="flex flex-wrap gap-1">
                  {key_binding_residues.slice(0, 10).map((residue, idx) => (
                    <span
                      key={idx}
                      className="px-2 py-0.5 bg-white/80 text-xs font-mono text-slate-700 rounded border border-slate-200"
                    >
                      {residue}
                    </span>
                  ))}
                  {key_binding_residues.length > 10 && (
                    <span className="px-2 py-0.5 text-xs text-slate-500">
                      +{key_binding_residues.length - 10} more
                    </span>
                  )}
                </div>
              </div>
            )}
          </div>
        </div>
      )}

      {/* ESMFold Validation (if available) */}
      {esmfold && (
        <div className="px-4 pb-4">
          <div className={`rounded-lg p-3 border ${
            esmfold.validation_passed
              ? 'bg-green-50 border-green-100'
              : 'bg-amber-50 border-amber-100'
          }`}>
            <div className="flex items-center gap-2 mb-3">
              {esmfold.validation_passed ? (
                <BadgeCheck className="h-4 w-4 text-green-600" />
              ) : (
                <AlertTriangle className="h-4 w-4 text-amber-600" />
              )}
              <span className={`text-xs font-medium uppercase tracking-wide ${
                esmfold.validation_passed ? 'text-green-700' : 'text-amber-700'
              }`}>
                ESMFold Structure Validation
              </span>
            </div>
            <div className="grid grid-cols-2 gap-4">
              {/* pLDDT */}
              <div>
                <div className="text-xs text-slate-500 mb-1">pLDDT Score</div>
                <div className="text-xl font-bold text-slate-900">
                  {(esmfold.plddt * 100).toFixed(0)}%
                </div>
                <div className={`text-xs font-medium ${plddtQ?.color}`}>
                  {plddtQ?.label}
                </div>
              </div>
              {/* RMSD */}
              {esmfold.rmsd != null && (
                <div>
                  <div className="text-xs text-slate-500 mb-1">Backbone RMSD</div>
                  <div className="text-xl font-bold text-slate-900">
                    {esmfold.rmsd.toFixed(2)}
                    <span className="text-sm font-normal text-slate-500 ml-1">Å</span>
                  </div>
                  <div className={`text-xs font-medium ${
                    esmfold.rmsd <= 1.5 ? 'text-green-600' :
                    esmfold.rmsd <= 2.0 ? 'text-blue-600' : 'text-amber-600'
                  }`}>
                    {esmfold.rmsd <= 1.5 ? 'Excellent' : esmfold.rmsd <= 2.0 ? 'Good' : 'Moderate'}
                  </div>
                </div>
              )}
            </div>
          </div>
        </div>
      )}

      {/* ESM Metrics (if available) */}
      {esm && (
        <div className="px-4 pb-4">
          <div className="bg-violet-50 rounded-lg p-3 border border-violet-100">
            <div className="flex items-center gap-2 mb-3">
              <Brain className="h-4 w-4 text-violet-600" />
              <span className="text-xs text-violet-700 font-medium uppercase tracking-wide">
                ESM-3 Sequence Quality
              </span>
            </div>
            <div className="grid grid-cols-2 gap-4">
              {/* Perplexity */}
              <div>
                <div className="text-xs text-slate-500 mb-1">Perplexity</div>
                <div className="text-xl font-bold text-slate-900">
                  {esm.perplexity.toFixed(2)}
                </div>
                <div className={`text-xs font-medium ${perplexityQ?.color}`}>
                  {perplexityQ?.label}
                </div>
              </div>
              {/* Confidence */}
              <div>
                <div className="text-xs text-slate-500 mb-1">Confidence</div>
                <div className="flex items-center gap-2">
                  <div className="text-xl font-bold text-slate-900">
                    {(esm.confidence * 100).toFixed(0)}%
                  </div>
                  <div className="flex-1 h-2 bg-slate-200 rounded-full overflow-hidden">
                    <div
                      className={`h-full rounded-full ${
                        esm.confidence >= 0.7 ? 'bg-green-500' :
                        esm.confidence >= 0.5 ? 'bg-amber-500' : 'bg-red-500'
                      }`}
                      style={{ width: `${esm.confidence * 100}%` }}
                    />
                  </div>
                </div>
              </div>
            </div>
          </div>
        </div>
      )}

      {/* MPNN Score (if available) */}
      {mpnn_score !== undefined && (
        <div className="px-4 pb-4">
          <div className="bg-blue-50 rounded-lg p-3 border border-blue-100">
            <div className="flex items-center justify-between">
              <div className="flex items-center gap-2">
                <FileEdit className="h-4 w-4 text-blue-600" />
                <span className="text-xs text-blue-700 font-medium">MPNN Score</span>
              </div>
              <span className="text-lg font-bold text-slate-900">
                {mpnn_score.toFixed(2)}
              </span>
            </div>
          </div>
        </div>
      )}

      {/* Quality Indicator Footer */}
      <div className="bg-slate-50 px-4 py-2 border-t border-slate-200">
        <div className="flex items-center gap-4 text-xs">
          <div className="flex items-center gap-1">
            <div className="w-2 h-2 rounded-full bg-green-500" />
            <span className="text-slate-600">Excellent</span>
          </div>
          <div className="flex items-center gap-1">
            <div className="w-2 h-2 rounded-full bg-blue-500" />
            <span className="text-slate-600">Good</span>
          </div>
          <div className="flex items-center gap-1">
            <div className="w-2 h-2 rounded-full bg-amber-500" />
            <span className="text-slate-600">Moderate</span>
          </div>
          <div className="flex items-center gap-1">
            <div className="w-2 h-2 rounded-full bg-red-500" />
            <span className="text-slate-600">Low</span>
          </div>
        </div>
      </div>
    </div>
  );
}
