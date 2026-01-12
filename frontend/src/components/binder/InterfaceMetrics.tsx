'use client';

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

export function InterfaceMetrics({ metrics, esm, mpnn_score, esmfold }: InterfaceMetricsProps) {
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
          <span className="material-symbols-outlined text-cyan-600">hub</span>
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
            <span className="material-symbols-outlined text-sm text-slate-400">touch_app</span>
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
            <span className="material-symbols-outlined text-sm text-slate-400">link</span>
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
            <span className="material-symbols-outlined text-sm text-slate-400">layers</span>
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
            <span className="material-symbols-outlined text-sm text-slate-400">deployed_code</span>
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
              <span className="material-symbols-outlined text-sm text-emerald-600">biotech</span>
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

      {/* ESMFold Validation (if available) */}
      {esmfold && (
        <div className="px-4 pb-4">
          <div className={`rounded-lg p-3 border ${
            esmfold.validation_passed
              ? 'bg-green-50 border-green-100'
              : 'bg-amber-50 border-amber-100'
          }`}>
            <div className="flex items-center gap-2 mb-3">
              <span className={`material-symbols-outlined text-sm ${
                esmfold.validation_passed ? 'text-green-600' : 'text-amber-600'
              }`}>
                {esmfold.validation_passed ? 'verified' : 'warning'}
              </span>
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
              <span className="material-symbols-outlined text-sm text-violet-600">psychology</span>
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
                <span className="material-symbols-outlined text-sm text-blue-600">edit_note</span>
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
