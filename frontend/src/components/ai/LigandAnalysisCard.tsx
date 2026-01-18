'use client';

import { useState } from 'react';
import { FlaskConical, ChevronUp, ChevronDown, Unlink, CheckCircle, Link, Box, Layers, ArrowRight, Hexagon } from 'lucide-react';

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
  // Interaction profile (from PLIP analysis)
  interactions?: {
    hbonds: number;
    hydrophobic: number;
    pi_stacking: number;
    salt_bridges: number;
    total: number;
  };
  key_binding_residues?: string[];
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
    return 'bg-muted text-muted-foreground';
  };

  return (
    <div className="bg-muted rounded-xl border border-border overflow-hidden">
      {/* Header */}
      <div className="px-4 py-3 flex items-center justify-between">
        <div className="flex items-center gap-2">
          <div className="w-8 h-8 rounded-lg bg-primary flex items-center justify-center">
            <FlaskConical className="h-4 w-4 text-primary-foreground" />
          </div>
          <h4 className="font-semibold text-foreground">Ligand Analysis</h4>
          {analysis.success && (
            <span className="text-xs bg-green-100 text-green-700 px-2 py-0.5 rounded-full">Complete</span>
          )}
        </div>
        <button
          onClick={() => setIsExpanded(!isExpanded)}
          className="text-sm text-purple-600 hover:text-purple-800 flex items-center gap-1"
        >
          {isExpanded ? <ChevronUp className="h-4 w-4" /> : <ChevronDown className="h-4 w-4" />}
          {isExpanded ? 'Less' : 'More'}
        </button>
      </div>

      <div className="px-4 pb-4 space-y-3">
        {/* Basic Ligand Info */}
        <div className="bg-card/60 rounded-lg p-3">
          <div className="flex items-start justify-between mb-2">
            <div>
              <div className="font-medium text-foreground">{analysis.ligand.name}</div>
              <code className="text-xs bg-primary/10 text-primary px-2 py-0.5 rounded inline-block mt-1 max-w-[280px] truncate">
                {analysis.ligand.smiles}
              </code>
            </div>
            {analysis.ligand.symmetry && (
              <span className={`text-xs px-2 py-0.5 rounded ${
                analysis.ligand.symmetry === 'symmetric'
                  ? 'bg-muted text-muted-foreground'
                  : 'bg-primary/10 text-primary'
              }`}>
                {analysis.ligand.symmetry}
              </span>
            )}
          </div>

          {/* Key Properties Grid */}
          <div className="grid grid-cols-4 gap-2 mt-3">
            {analysis.ligand.molecular_weight && (
              <div className="text-center p-2 bg-muted/50 rounded">
                <div className="text-lg font-bold text-foreground">{analysis.ligand.molecular_weight.toFixed(0)}</div>
                <div className="text-[10px] text-muted-foreground uppercase">MW (Da)</div>
              </div>
            )}
            {analysis.ligand.num_heavy_atoms && (
              <div className="text-center p-2 bg-muted/50 rounded">
                <div className="text-lg font-bold text-foreground">{analysis.ligand.num_heavy_atoms}</div>
                <div className="text-[10px] text-muted-foreground uppercase">Heavy Atoms</div>
              </div>
            )}
            {analysis.ligand.num_rings !== undefined && (
              <div className="text-center p-2 bg-muted/50 rounded">
                <div className="text-lg font-bold text-foreground">{analysis.ligand.num_rings}</div>
                <div className="text-[10px] text-muted-foreground uppercase">Rings</div>
              </div>
            )}
            {analysis.ligand.num_rotatable_bonds !== undefined && (
              <div className="text-center p-2 bg-muted/50 rounded">
                <div className="text-lg font-bold text-foreground">{analysis.ligand.num_rotatable_bonds}</div>
                <div className="text-[10px] text-muted-foreground uppercase">Rotatable</div>
              </div>
            )}
          </div>
        </div>

        {/* Binding Site & Topology */}
        <div className="grid grid-cols-2 gap-3">
          {/* Binding Type */}
          <div className="bg-card/60 rounded-lg p-3">
            <div className="text-xs text-muted-foreground uppercase mb-1">Binding Type</div>
            <div className={`inline-flex items-center gap-1 font-medium px-2 py-1 rounded text-sm ${
              analysis.binding_site.type === 'interface'
                ? 'bg-primary/10 text-primary'
                : analysis.binding_site.type === 'buried'
                  ? 'bg-muted text-muted-foreground'
                  : 'bg-muted text-muted-foreground'
            }`}>
              {analysis.binding_site.type === 'interface' ? <Unlink className="h-4 w-4" /> :
               analysis.binding_site.type === 'buried' ? <Box className="h-4 w-4" /> : <Layers className="h-4 w-4" />}
              {analysis.binding_site.type.charAt(0).toUpperCase() + analysis.binding_site.type.slice(1)}
            </div>
          </div>

          {/* Topology */}
          <div className="bg-card/60 rounded-lg p-3">
            <div className="text-xs text-muted-foreground uppercase mb-1">Topology</div>
            <div className={`inline-flex items-center gap-1 font-medium px-2 py-1 rounded text-sm ${
              analysis.topology.separable
                ? 'bg-green-100 text-green-700'
                : 'bg-red-100 text-red-700'
            }`}>
              {analysis.topology.separable ? <CheckCircle className="h-4 w-4" /> : <Link className="h-4 w-4" />}
              {analysis.topology.separable ? 'Separable' : 'Entangled'}
            </div>
          </div>
        </div>

        {/* Expanded Details */}
        {isExpanded && (
          <>
            {/* Additional Properties */}
            {(analysis.ligand.num_hbd !== undefined || analysis.ligand.logp !== undefined) && (
              <div className="bg-card/60 rounded-lg p-3">
                <div className="text-xs text-muted-foreground uppercase mb-2">Drug-like Properties</div>
                <div className="grid grid-cols-4 gap-2">
                  {analysis.ligand.num_hbd !== undefined && (
                    <div className="text-center">
                      <div className="font-bold text-foreground">{analysis.ligand.num_hbd}</div>
                      <div className="text-[10px] text-muted-foreground">HB Donors</div>
                    </div>
                  )}
                  {analysis.ligand.num_hba !== undefined && (
                    <div className="text-center">
                      <div className="font-bold text-foreground">{analysis.ligand.num_hba}</div>
                      <div className="text-[10px] text-muted-foreground">HB Acceptors</div>
                    </div>
                  )}
                  {analysis.ligand.logp !== undefined && (
                    <div className="text-center">
                      <div className="font-bold text-foreground">{analysis.ligand.logp.toFixed(1)}</div>
                      <div className="text-[10px] text-muted-foreground">LogP</div>
                    </div>
                  )}
                  {analysis.ligand.tpsa !== undefined && (
                    <div className="text-center">
                      <div className="font-bold text-foreground">{analysis.ligand.tpsa.toFixed(0)}</div>
                      <div className="text-[10px] text-muted-foreground">TPSA (Å²)</div>
                    </div>
                  )}
                </div>
              </div>
            )}

            {/* Atom Exposure Analysis */}
            {analysis.atom_analysis && (
              <div className="bg-card/60 rounded-lg p-3">
                <div className="text-xs text-muted-foreground uppercase mb-2">Atom Exposure Analysis</div>
                <div className="grid grid-cols-3 gap-2">
                  <div className="bg-primary/5 rounded p-2">
                    <div className="text-xs text-primary font-medium mb-1">Face A</div>
                    <div className="text-sm font-bold text-primary">{analysis.atom_analysis.face_a_atoms.length} atoms</div>
                    <div className="text-[10px] text-primary/70 truncate">
                      {analysis.atom_analysis.face_a_atoms.slice(0, 4).join(', ')}
                      {analysis.atom_analysis.face_a_atoms.length > 4 && '...'}
                    </div>
                  </div>
                  <div className="bg-muted rounded p-2">
                    <div className="text-xs text-muted-foreground font-medium mb-1">Core</div>
                    <div className="text-sm font-bold text-foreground">{analysis.atom_analysis.core_atoms.length} atoms</div>
                    <div className="text-[10px] text-muted-foreground">Buried/Central</div>
                  </div>
                  <div className="bg-muted/50 rounded p-2">
                    <div className="text-xs text-muted-foreground font-medium mb-1">Face B</div>
                    <div className="text-sm font-bold text-foreground">{analysis.atom_analysis.face_b_atoms.length} atoms</div>
                    <div className="text-[10px] text-muted-foreground truncate">
                      {analysis.atom_analysis.face_b_atoms.slice(0, 4).join(', ')}
                      {analysis.atom_analysis.face_b_atoms.length > 4 && '...'}
                    </div>
                  </div>
                </div>
                {analysis.atom_analysis.contact_potential > 0 && (
                  <div className="mt-2 flex items-center justify-between text-xs">
                    <span className="text-muted-foreground">Contact Potential</span>
                    <span className="font-medium text-primary">{analysis.atom_analysis.contact_potential.toFixed(1)}</span>
                  </div>
                )}
              </div>
            )}

            {/* Interaction Profile (from PLIP analysis) */}
            {analysis.interactions && (
              <div className="bg-card/60 rounded-lg p-3">
                <div className="text-xs text-muted-foreground uppercase mb-2">Interaction Profile</div>
                <div className="grid grid-cols-5 gap-2">
                  <div className="text-center p-2 bg-muted rounded">
                    <div className="text-lg font-bold text-primary">{analysis.interactions.hbonds}</div>
                    <div className="text-[10px] text-muted-foreground">H-bonds</div>
                  </div>
                  <div className="text-center p-2 bg-muted rounded">
                    <div className="text-lg font-bold text-foreground">{analysis.interactions.hydrophobic}</div>
                    <div className="text-[10px] text-muted-foreground">Hydrophobic</div>
                  </div>
                  <div className="text-center p-2 bg-muted rounded">
                    <div className="text-lg font-bold text-foreground">{analysis.interactions.pi_stacking}</div>
                    <div className="text-[10px] text-muted-foreground">Pi-stack</div>
                  </div>
                  <div className="text-center p-2 bg-muted rounded">
                    <div className="text-lg font-bold text-foreground">{analysis.interactions.salt_bridges}</div>
                    <div className="text-[10px] text-muted-foreground">Salt Bridge</div>
                  </div>
                  <div className="text-center p-2 bg-muted rounded">
                    <div className="text-lg font-bold text-foreground">{analysis.interactions.total}</div>
                    <div className="text-[10px] text-muted-foreground">Total</div>
                  </div>
                </div>
                {analysis.key_binding_residues && analysis.key_binding_residues.length > 0 && (
                  <div className="mt-2">
                    <div className="text-xs text-muted-foreground mb-1">Key Binding Residues</div>
                    <div className="flex flex-wrap gap-1">
                      {analysis.key_binding_residues.slice(0, 6).map((residue, idx) => (
                        <span key={idx} className="text-xs bg-primary/10 text-primary px-1.5 py-0.5 rounded">
                          {residue}
                        </span>
                      ))}
                      {analysis.key_binding_residues.length > 6 && (
                        <span className="text-xs text-muted-foreground">+{analysis.key_binding_residues.length - 6} more</span>
                      )}
                    </div>
                  </div>
                )}
              </div>
            )}

            {/* Risk Assessment */}
            {analysis.topology.ring_threading_risk && (
              <div className="bg-card/60 rounded-lg p-3">
                <div className="text-xs text-muted-foreground uppercase mb-2">Risk Assessment</div>
                <div className="flex items-center justify-between">
                  <span className="text-sm text-muted-foreground">Ring Threading Risk</span>
                  <span className={`text-xs px-2 py-0.5 rounded font-medium ${getRiskColor(analysis.topology.ring_threading_risk)}`}>
                    {analysis.topology.ring_threading_risk.toUpperCase()}
                  </span>
                </div>
                {analysis.binding_site.recommended_chain_length && (
                  <div className="flex items-center justify-between mt-2">
                    <span className="text-sm text-muted-foreground">Recommended Chain Length</span>
                    <span className="text-sm font-medium text-foreground">{analysis.binding_site.recommended_chain_length}</span>
                  </div>
                )}
              </div>
            )}

            {/* Functional Groups */}
            {analysis.ligand.functional_groups && analysis.ligand.functional_groups.length > 0 && (
              <div className="bg-card/60 rounded-lg p-3">
                <div className="text-xs text-muted-foreground uppercase mb-2">Functional Groups</div>
                <div className="flex flex-wrap gap-1">
                  {analysis.ligand.functional_groups.map((group, idx) => (
                    <span key={idx} className="text-xs bg-muted text-foreground px-2 py-0.5 rounded">
                      {group}
                    </span>
                  ))}
                </div>
              </div>
            )}

            {/* Recommendations */}
            {analysis.recommendations && analysis.recommendations.length > 0 && (
              <div className="bg-card/60 rounded-lg p-3">
                <div className="text-xs text-muted-foreground uppercase mb-2">Recommendations</div>
                <ul className="space-y-1">
                  {analysis.recommendations.map((rec, idx) => (
                    <li key={idx} className="text-xs text-muted-foreground flex items-start gap-1.5">
                      <ArrowRight className="h-4 w-4 text-primary mt-0.5 flex-shrink-0" />
                      {rec}
                    </li>
                  ))}
                </ul>
              </div>
            )}
          </>
        )}

        {/* Design Approach Diagram */}
        <div className="bg-card/60 rounded-lg p-3">
          <div className="text-xs text-muted-foreground uppercase mb-2">Design Workflow</div>
          <div className="flex items-center justify-center gap-1">
            {/* Step 1: Chain A */}
            <div className="flex flex-col items-center">
              <div className="w-10 h-10 rounded bg-primary/20 flex items-center justify-center text-primary font-bold text-sm border-2 border-primary/30">
                A
              </div>
              <span className="text-[10px] text-muted-foreground mt-1">Step 1</span>
            </div>

            {/* Arrow */}
            <ArrowRight className="h-4 w-4 text-muted-foreground" />

            {/* Ligand */}
            <div className="flex flex-col items-center">
              <div className="w-10 h-10 rounded bg-muted flex items-center justify-center border-2 border-border">
                <Hexagon className="h-5 w-5 text-primary" />
              </div>
              <span className="text-[10px] text-muted-foreground mt-1">Interface</span>
            </div>

            {/* Arrow */}
            <ArrowRight className="h-4 w-4 text-muted-foreground" />

            {/* Step 2: Chain B */}
            <div className="flex flex-col items-center">
              <div className="w-10 h-10 rounded bg-muted flex items-center justify-center text-muted-foreground font-bold text-sm border-2 border-border">
                B
              </div>
              <span className="text-[10px] text-muted-foreground mt-1">Step 2</span>
            </div>

            {/* Arrow */}
            <ArrowRight className="h-4 w-4 text-muted-foreground" />

            {/* Result: Dimer */}
            <div className="flex flex-col items-center">
              <div className="w-12 h-10 rounded bg-green-100 flex items-center justify-center gap-0.5 border-2 border-green-300">
                <span className="text-primary font-bold text-xs">A</span>
                <span className="text-muted-foreground text-xs">+</span>
                <span className="text-foreground font-bold text-xs">B</span>
              </div>
              <span className="text-[10px] text-green-600 font-medium mt-1">Dimer</span>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
