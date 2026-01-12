'use client';

import { useState } from 'react';

// Ligand analysis data structure
export interface LigandAnalysis {
  success: boolean;
  ligand: {
    name: string;
    smiles: string;
    num_atoms?: number;
    num_heavy_atoms?: number;
    molecular_weight?: number;
    // Enhanced properties
    num_rings?: number;
    num_rotatable_bonds?: number;
    num_hbd?: number;  // H-bond donors
    num_hba?: number;  // H-bond acceptors
    logp?: number;
    tpsa?: number;  // Topological polar surface area
    symmetry?: 'symmetric' | 'asymmetric' | 'pseudo-symmetric';
    functional_groups?: string[];
  };
  binding_site: {
    type: 'interface' | 'buried' | 'surface';
    approach: string;
    description: string;
    recommended_chain_length?: string;
  };
  topology: {
    separable: boolean;
    description: string;
    ring_threading_risk?: 'low' | 'medium' | 'high';
  };
  // Atom exposure analysis
  atom_analysis?: {
    face_a_atoms: string[];  // Atoms exposed to Chain A
    face_b_atoms: string[];  // Atoms exposed to Chain B
    core_atoms: string[];    // Buried atoms
    contact_potential: number;  // Estimated contact score
  };
  // Design recommendations
  recommendations?: string[];
}

interface LigandAnalysisCardProps {
  analysis: LigandAnalysis;
  expanded?: boolean;
}

export function LigandAnalysisCard({ analysis, expanded = false }: LigandAnalysisCardProps) {
  const [isExpanded, setIsExpanded] = useState(expanded);

  const getRiskColor = (risk?: 'low' | 'medium' | 'high') => {
    if (risk === 'low') return 'bg-green-100 text-green-700';
    if (risk === 'medium') return 'bg-amber-100 text-amber-700';
    if (risk === 'high') return 'bg-red-100 text-red-700';
    return 'bg-slate-100 text-slate-600';
  };

  return (
    <div className="bg-gradient-to-br from-purple-50 to-pink-50 rounded-xl border border-purple-200 overflow-hidden">
      {/* Header */}
      <div className="px-4 py-3 flex items-center justify-between">
        <div className="flex items-center gap-2">
          <div className="w-8 h-8 rounded-lg bg-gradient-to-br from-purple-500 to-pink-500 flex items-center justify-center">
            <span className="material-symbols-outlined text-white text-sm">science</span>
          </div>
          <h4 className="font-semibold text-slate-900">Ligand Analysis</h4>
          {analysis.success && (
            <span className="text-xs bg-green-100 text-green-700 px-2 py-0.5 rounded-full">Complete</span>
          )}
        </div>
        <button
          onClick={() => setIsExpanded(!isExpanded)}
          className="text-sm text-purple-600 hover:text-purple-800 flex items-center gap-1"
        >
          <span className="material-symbols-outlined text-sm">
            {isExpanded ? 'expand_less' : 'expand_more'}
          </span>
          {isExpanded ? 'Less' : 'More'}
        </button>
      </div>

      <div className="px-4 pb-4 space-y-3">
        {/* Basic Ligand Info */}
        <div className="bg-white/60 rounded-lg p-3">
          <div className="flex items-start justify-between mb-2">
            <div>
              <div className="font-medium text-slate-900">{analysis.ligand.name}</div>
              <code className="text-xs bg-purple-100 text-purple-800 px-2 py-0.5 rounded inline-block mt-1 max-w-[280px] truncate">
                {analysis.ligand.smiles}
              </code>
            </div>
            {analysis.ligand.symmetry && (
              <span className={`text-xs px-2 py-0.5 rounded ${
                analysis.ligand.symmetry === 'symmetric'
                  ? 'bg-blue-100 text-blue-700'
                  : 'bg-purple-100 text-purple-700'
              }`}>
                {analysis.ligand.symmetry}
              </span>
            )}
          </div>

          {/* Key Properties Grid */}
          <div className="grid grid-cols-4 gap-2 mt-3">
            {analysis.ligand.molecular_weight && (
              <div className="text-center p-2 bg-slate-50 rounded">
                <div className="text-lg font-bold text-slate-900">{analysis.ligand.molecular_weight.toFixed(0)}</div>
                <div className="text-[10px] text-slate-500 uppercase">MW (Da)</div>
              </div>
            )}
            {analysis.ligand.num_heavy_atoms && (
              <div className="text-center p-2 bg-slate-50 rounded">
                <div className="text-lg font-bold text-slate-900">{analysis.ligand.num_heavy_atoms}</div>
                <div className="text-[10px] text-slate-500 uppercase">Heavy Atoms</div>
              </div>
            )}
            {analysis.ligand.num_rings !== undefined && (
              <div className="text-center p-2 bg-slate-50 rounded">
                <div className="text-lg font-bold text-slate-900">{analysis.ligand.num_rings}</div>
                <div className="text-[10px] text-slate-500 uppercase">Rings</div>
              </div>
            )}
            {analysis.ligand.num_rotatable_bonds !== undefined && (
              <div className="text-center p-2 bg-slate-50 rounded">
                <div className="text-lg font-bold text-slate-900">{analysis.ligand.num_rotatable_bonds}</div>
                <div className="text-[10px] text-slate-500 uppercase">Rotatable</div>
              </div>
            )}
          </div>
        </div>

        {/* Binding Site & Topology */}
        <div className="grid grid-cols-2 gap-3">
          {/* Binding Type */}
          <div className="bg-white/60 rounded-lg p-3">
            <div className="text-xs text-slate-500 uppercase mb-1">Binding Type</div>
            <div className={`inline-flex items-center gap-1 font-medium px-2 py-1 rounded text-sm ${
              analysis.binding_site.type === 'interface'
                ? 'bg-purple-100 text-purple-700'
                : analysis.binding_site.type === 'buried'
                  ? 'bg-amber-100 text-amber-700'
                  : 'bg-blue-100 text-blue-700'
            }`}>
              <span className="material-symbols-outlined text-sm">
                {analysis.binding_site.type === 'interface' ? 'call_split' :
                 analysis.binding_site.type === 'buried' ? 'deployed_code' : 'layers'}
              </span>
              {analysis.binding_site.type.charAt(0).toUpperCase() + analysis.binding_site.type.slice(1)}
            </div>
          </div>

          {/* Topology */}
          <div className="bg-white/60 rounded-lg p-3">
            <div className="text-xs text-slate-500 uppercase mb-1">Topology</div>
            <div className={`inline-flex items-center gap-1 font-medium px-2 py-1 rounded text-sm ${
              analysis.topology.separable
                ? 'bg-green-100 text-green-700'
                : 'bg-red-100 text-red-700'
            }`}>
              <span className="material-symbols-outlined text-sm">
                {analysis.topology.separable ? 'check_circle' : 'link'}
              </span>
              {analysis.topology.separable ? 'Separable' : 'Entangled'}
            </div>
          </div>
        </div>

        {/* Expanded Details */}
        {isExpanded && (
          <>
            {/* Additional Properties */}
            {(analysis.ligand.num_hbd !== undefined || analysis.ligand.logp !== undefined) && (
              <div className="bg-white/60 rounded-lg p-3">
                <div className="text-xs text-slate-500 uppercase mb-2">Drug-like Properties</div>
                <div className="grid grid-cols-4 gap-2">
                  {analysis.ligand.num_hbd !== undefined && (
                    <div className="text-center">
                      <div className="font-bold text-slate-900">{analysis.ligand.num_hbd}</div>
                      <div className="text-[10px] text-slate-500">HB Donors</div>
                    </div>
                  )}
                  {analysis.ligand.num_hba !== undefined && (
                    <div className="text-center">
                      <div className="font-bold text-slate-900">{analysis.ligand.num_hba}</div>
                      <div className="text-[10px] text-slate-500">HB Acceptors</div>
                    </div>
                  )}
                  {analysis.ligand.logp !== undefined && (
                    <div className="text-center">
                      <div className="font-bold text-slate-900">{analysis.ligand.logp.toFixed(1)}</div>
                      <div className="text-[10px] text-slate-500">LogP</div>
                    </div>
                  )}
                  {analysis.ligand.tpsa !== undefined && (
                    <div className="text-center">
                      <div className="font-bold text-slate-900">{analysis.ligand.tpsa.toFixed(0)}</div>
                      <div className="text-[10px] text-slate-500">TPSA (Å²)</div>
                    </div>
                  )}
                </div>
              </div>
            )}

            {/* Atom Exposure Analysis */}
            {analysis.atom_analysis && (
              <div className="bg-white/60 rounded-lg p-3">
                <div className="text-xs text-slate-500 uppercase mb-2">Atom Exposure Analysis</div>
                <div className="grid grid-cols-3 gap-2">
                  <div className="bg-purple-50 rounded p-2">
                    <div className="text-xs text-purple-600 font-medium mb-1">Face A</div>
                    <div className="text-sm font-bold text-purple-700">{analysis.atom_analysis.face_a_atoms.length} atoms</div>
                    <div className="text-[10px] text-purple-500 truncate">
                      {analysis.atom_analysis.face_a_atoms.slice(0, 4).join(', ')}
                      {analysis.atom_analysis.face_a_atoms.length > 4 && '...'}
                    </div>
                  </div>
                  <div className="bg-slate-100 rounded p-2">
                    <div className="text-xs text-slate-600 font-medium mb-1">Core</div>
                    <div className="text-sm font-bold text-slate-700">{analysis.atom_analysis.core_atoms.length} atoms</div>
                    <div className="text-[10px] text-slate-500">Buried/Central</div>
                  </div>
                  <div className="bg-cyan-50 rounded p-2">
                    <div className="text-xs text-cyan-600 font-medium mb-1">Face B</div>
                    <div className="text-sm font-bold text-cyan-700">{analysis.atom_analysis.face_b_atoms.length} atoms</div>
                    <div className="text-[10px] text-cyan-500 truncate">
                      {analysis.atom_analysis.face_b_atoms.slice(0, 4).join(', ')}
                      {analysis.atom_analysis.face_b_atoms.length > 4 && '...'}
                    </div>
                  </div>
                </div>
                {analysis.atom_analysis.contact_potential > 0 && (
                  <div className="mt-2 flex items-center justify-between text-xs">
                    <span className="text-slate-500">Contact Potential</span>
                    <span className="font-medium text-purple-600">{analysis.atom_analysis.contact_potential.toFixed(1)}</span>
                  </div>
                )}
              </div>
            )}

            {/* Risk Assessment */}
            {analysis.topology.ring_threading_risk && (
              <div className="bg-white/60 rounded-lg p-3">
                <div className="text-xs text-slate-500 uppercase mb-2">Risk Assessment</div>
                <div className="flex items-center justify-between">
                  <span className="text-sm text-slate-600">Ring Threading Risk</span>
                  <span className={`text-xs px-2 py-0.5 rounded font-medium ${getRiskColor(analysis.topology.ring_threading_risk)}`}>
                    {analysis.topology.ring_threading_risk.toUpperCase()}
                  </span>
                </div>
                {analysis.binding_site.recommended_chain_length && (
                  <div className="flex items-center justify-between mt-2">
                    <span className="text-sm text-slate-600">Recommended Chain Length</span>
                    <span className="text-sm font-medium text-slate-700">{analysis.binding_site.recommended_chain_length}</span>
                  </div>
                )}
              </div>
            )}

            {/* Functional Groups */}
            {analysis.ligand.functional_groups && analysis.ligand.functional_groups.length > 0 && (
              <div className="bg-white/60 rounded-lg p-3">
                <div className="text-xs text-slate-500 uppercase mb-2">Functional Groups</div>
                <div className="flex flex-wrap gap-1">
                  {analysis.ligand.functional_groups.map((group, idx) => (
                    <span key={idx} className="text-xs bg-slate-200 text-slate-700 px-2 py-0.5 rounded">
                      {group}
                    </span>
                  ))}
                </div>
              </div>
            )}

            {/* Recommendations */}
            {analysis.recommendations && analysis.recommendations.length > 0 && (
              <div className="bg-white/60 rounded-lg p-3">
                <div className="text-xs text-slate-500 uppercase mb-2">Recommendations</div>
                <ul className="space-y-1">
                  {analysis.recommendations.map((rec, idx) => (
                    <li key={idx} className="text-xs text-slate-600 flex items-start gap-1.5">
                      <span className="material-symbols-outlined text-purple-500 text-sm mt-0.5">arrow_right</span>
                      {rec}
                    </li>
                  ))}
                </ul>
              </div>
            )}
          </>
        )}

        {/* Design Approach Diagram */}
        <div className="bg-white/60 rounded-lg p-3">
          <div className="text-xs text-slate-500 uppercase mb-2">Design Workflow</div>
          <div className="flex items-center justify-center gap-1">
            {/* Step 1: Chain A */}
            <div className="flex flex-col items-center">
              <div className="w-10 h-10 rounded bg-purple-200 flex items-center justify-center text-purple-700 font-bold text-sm border-2 border-purple-300">
                A
              </div>
              <span className="text-[10px] text-slate-500 mt-1">Step 1</span>
            </div>

            {/* Arrow */}
            <span className="material-symbols-outlined text-slate-400 text-sm">arrow_forward</span>

            {/* Ligand */}
            <div className="flex flex-col items-center">
              <div className="w-10 h-10 rounded bg-gradient-to-br from-pink-200 to-purple-200 flex items-center justify-center border-2 border-pink-300">
                <span className="material-symbols-outlined text-pink-700 text-base">hexagon</span>
              </div>
              <span className="text-[10px] text-slate-500 mt-1">Interface</span>
            </div>

            {/* Arrow */}
            <span className="material-symbols-outlined text-slate-400 text-sm">arrow_forward</span>

            {/* Step 2: Chain B */}
            <div className="flex flex-col items-center">
              <div className="w-10 h-10 rounded bg-cyan-200 flex items-center justify-center text-cyan-700 font-bold text-sm border-2 border-cyan-300">
                B
              </div>
              <span className="text-[10px] text-slate-500 mt-1">Step 2</span>
            </div>

            {/* Arrow */}
            <span className="material-symbols-outlined text-slate-400 text-sm">arrow_forward</span>

            {/* Result: Dimer */}
            <div className="flex flex-col items-center">
              <div className="w-12 h-10 rounded bg-green-100 flex items-center justify-center gap-0.5 border-2 border-green-300">
                <span className="text-purple-600 font-bold text-xs">A</span>
                <span className="text-pink-500 text-xs">+</span>
                <span className="text-cyan-600 font-bold text-xs">B</span>
              </div>
              <span className="text-[10px] text-green-600 font-medium mt-1">Dimer</span>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
