'use client';

import { useState } from 'react';
import { BarChart3, ChevronUp, ChevronDown, XCircle, CheckCircle, Unlink, Link, Info, AlertTriangle, Hexagon } from 'lucide-react';

// Ligand evaluation data structure
export interface LigandEvaluation {
  success: boolean;
  approach: 'asymmetric' | 'sequential' | 'full';
  chain_a?: {
    contacts: number;
    exposed_atoms: string;
    affinity?: number;
    // Enhanced metrics
    residue_contacts?: { resnum: number; resname: string; distance: number }[];
    buried_sasa?: number;
    hbonds?: number;
    hydrophobic_contacts?: number;
  };
  chain_b?: {
    contacts: number;
    exposed_atoms: string;
    affinity?: number;
    // Enhanced metrics
    residue_contacts?: { resnum: number; resname: string; distance: number }[];
    buried_sasa?: number;
    hbonds?: number;
    hydrophobic_contacts?: number;
  };
  dimer?: {
    affinity: number;
    contacts_a: number;
    contacts_b: number;
    has_clashes: boolean;
    separable: boolean;
    interface_area?: number;
    // Enhanced metrics
    energy_breakdown?: {
      vdw: number;
      electrostatic: number;
      hbond: number;
      desolvation: number;
    };
    shape_complementarity?: number;
    buried_sasa?: number;
    crossing_count?: number;  // For topology check
  };
  // Additional scoring
  gnina_score?: {
    cnn_score: number;
    cnn_affinity: number;
    pose_rmsd?: number;
  };
  esm_score?: number;
  overall_pass: boolean;
  warnings?: string[];
}

interface LigandEvaluationCardProps {
  evaluation: LigandEvaluation;
  expanded?: boolean;
}

export function LigandEvaluationCard({ evaluation, expanded = false }: LigandEvaluationCardProps) {
  const [isExpanded, setIsExpanded] = useState(expanded);
  const [activeTab, setActiveTab] = useState<'overview' | 'chain_a' | 'chain_b' | 'energy'>('overview');

  const passColor = evaluation.overall_pass ? 'text-green-600' : 'text-red-600';
  const passBg = evaluation.overall_pass ? 'bg-green-50 border-green-200' : 'bg-red-50 border-red-200';

  const getAffinityQuality = (affinity: number) => {
    if (affinity < -5) return { label: 'Excellent', color: 'text-green-600 bg-green-100' };
    if (affinity < -3) return { label: 'Good', color: 'text-emerald-600 bg-emerald-100' };
    if (affinity < -1.5) return { label: 'Moderate', color: 'text-amber-600 bg-amber-100' };
    return { label: 'Weak', color: 'text-red-600 bg-red-100' };
  };

  return (
    <div className="bg-muted rounded-xl border border-border overflow-hidden">
      {/* Header */}
      <div className="px-4 py-3 flex items-center justify-between border-b border-border">
        <div className="flex items-center gap-2">
          <div className="w-8 h-8 rounded-lg bg-primary flex items-center justify-center">
            <BarChart3 className="h-4 w-4 text-primary-foreground" />
          </div>
          <h4 className="font-semibold text-foreground">Dimer Evaluation</h4>
        </div>
        <div className="flex items-center gap-2">
          <div className={`px-3 py-1 rounded-full text-sm font-medium border ${passBg} ${passColor}`}>
            {evaluation.overall_pass ? 'PASS' : 'FAIL'}
          </div>
          <button
            onClick={() => setIsExpanded(!isExpanded)}
            className="text-sm text-primary hover:text-primary/80 flex items-center gap-1"
          >
            {isExpanded ? <ChevronUp className="h-4 w-4" /> : <ChevronDown className="h-4 w-4" />}
          </button>
        </div>
      </div>

      <div className="p-4 space-y-4">
        {/* Approach Badge */}
        <div className="flex items-center gap-3">
          <span className="text-sm text-muted-foreground">Design Approach:</span>
          <span className="px-2 py-0.5 bg-primary/10 text-primary rounded text-sm font-medium capitalize">
            {evaluation.approach}
          </span>
          {evaluation.gnina_score && (
            <span className="px-2 py-0.5 bg-blue-100 text-blue-700 rounded text-xs font-medium">
              GNINA Scored
            </span>
          )}
        </div>

        {/* Quick Stats Bar */}
        {evaluation.dimer && (
          <div className="grid grid-cols-4 gap-2">
            {/* Affinity */}
            <div className="bg-white/70 rounded-lg p-2 text-center">
              <div className={`text-lg font-bold ${getAffinityQuality(evaluation.dimer.affinity).color.split(' ')[0]}`}>
                {evaluation.dimer.affinity.toFixed(1)}
              </div>
              <div className="text-[10px] text-muted-foreground uppercase">kcal/mol</div>
            </div>
            {/* Contacts */}
            <div className="bg-white/70 rounded-lg p-2 text-center">
              <div className="text-lg font-bold text-foreground">
                {evaluation.dimer.contacts_a + evaluation.dimer.contacts_b}
              </div>
              <div className="text-[10px] text-muted-foreground uppercase">Contacts</div>
            </div>
            {/* Interface Area */}
            {evaluation.dimer.interface_area && (
              <div className="bg-white/70 rounded-lg p-2 text-center">
                <div className="text-lg font-bold text-foreground">
                  {evaluation.dimer.interface_area.toFixed(0)}
                </div>
                <div className="text-[10px] text-muted-foreground uppercase">Å² Area</div>
              </div>
            )}
            {/* Status */}
            <div className="bg-white/70 rounded-lg p-2 text-center">
              <div className={`text-lg font-bold ${evaluation.dimer.separable && !evaluation.dimer.has_clashes ? 'text-green-600' : 'text-red-600'}`}>
                {evaluation.dimer.separable && !evaluation.dimer.has_clashes ? '✓' : '✗'}
              </div>
              <div className="text-[10px] text-muted-foreground uppercase">Valid</div>
            </div>
          </div>
        )}

        {/* Tab Navigation */}
        {isExpanded && (
          <div className="flex border-b border-border">
            {['overview', 'chain_a', 'chain_b', 'energy'].map((tab) => (
              <button
                key={tab}
                onClick={() => setActiveTab(tab as typeof activeTab)}
                className={`px-3 py-2 text-xs font-medium transition-colors ${
                  activeTab === tab
                    ? 'text-primary border-b-2 border-primary'
                    : 'text-muted-foreground hover:text-foreground'
                }`}
              >
                {tab === 'chain_a' ? 'Chain A' : tab === 'chain_b' ? 'Chain B' : tab.charAt(0).toUpperCase() + tab.slice(1)}
              </button>
            ))}
          </div>
        )}

        {/* Tab Content */}
        {isExpanded && (
          <div className="min-h-[200px]">
            {/* Overview Tab */}
            {activeTab === 'overview' && evaluation.dimer && (
              <div className="space-y-3">
                {/* Affinity Analysis */}
                <div className="bg-white/60 rounded-lg p-3">
                  <div className="flex items-center justify-between mb-2">
                    <span className="text-sm font-medium text-foreground">Binding Affinity</span>
                    <span className={`px-2 py-0.5 rounded text-xs font-medium ${getAffinityQuality(evaluation.dimer.affinity).color}`}>
                      {getAffinityQuality(evaluation.dimer.affinity).label}
                    </span>
                  </div>
                  <div className="flex items-end gap-2">
                    <span className={`text-2xl font-bold ${getAffinityQuality(evaluation.dimer.affinity).color.split(' ')[0]}`}>
                      {evaluation.dimer.affinity.toFixed(2)}
                    </span>
                    <span className="text-sm text-muted-foreground mb-1">kcal/mol</span>
                  </div>
                  {/* Affinity bar visualization */}
                  <div className="mt-2 h-2 bg-muted rounded-full overflow-hidden">
                    <div
                      className={`h-full rounded-full transition-all ${
                        evaluation.dimer.affinity < -5 ? 'bg-green-500' :
                        evaluation.dimer.affinity < -3 ? 'bg-emerald-500' :
                        evaluation.dimer.affinity < -1.5 ? 'bg-amber-500' : 'bg-red-500'
                      }`}
                      style={{ width: `${Math.min(100, Math.max(0, (Math.abs(evaluation.dimer.affinity) / 8) * 100))}%` }}
                    />
                  </div>
                  <div className="flex justify-between text-[10px] text-muted-foreground mt-1">
                    <span>Weak (0)</span>
                    <span>Strong (-8)</span>
                  </div>
                </div>

                {/* Contact Distribution */}
                <div className="bg-card/60 rounded-lg p-3">
                  <div className="text-sm font-medium text-foreground mb-2">Contact Distribution</div>
                  <div className="flex items-center gap-2">
                    <div className="flex-1">
                      <div className="flex items-center justify-between text-xs mb-1">
                        <span className="text-primary font-medium">Chain A</span>
                        <span>{evaluation.dimer.contacts_a}</span>
                      </div>
                      <div className="h-3 bg-primary/20 rounded-full overflow-hidden">
                        <div
                          className="h-full bg-primary rounded-full"
                          style={{ width: `${(evaluation.dimer.contacts_a / (evaluation.dimer.contacts_a + evaluation.dimer.contacts_b)) * 100}%` }}
                        />
                      </div>
                    </div>
                    <div className="flex-1">
                      <div className="flex items-center justify-between text-xs mb-1">
                        <span className="text-muted-foreground font-medium">Chain B</span>
                        <span>{evaluation.dimer.contacts_b}</span>
                      </div>
                      <div className="h-3 bg-muted rounded-full overflow-hidden">
                        <div
                          className="h-full bg-muted-foreground rounded-full"
                          style={{ width: `${(evaluation.dimer.contacts_b / (evaluation.dimer.contacts_a + evaluation.dimer.contacts_b)) * 100}%` }}
                        />
                      </div>
                    </div>
                  </div>
                </div>

                {/* Status Checks */}
                <div className="bg-card/60 rounded-lg p-3">
                  <div className="text-sm font-medium text-foreground mb-2">Quality Checks</div>
                  <div className="space-y-2">
                    <div className="flex items-center justify-between">
                      <span className="text-sm text-muted-foreground">Steric Clashes</span>
                      <span className={`flex items-center gap-1 text-sm font-medium ${
                        evaluation.dimer.has_clashes ? 'text-red-600' : 'text-green-600'
                      }`}>
                        {evaluation.dimer.has_clashes ? <XCircle className="h-4 w-4" /> : <CheckCircle className="h-4 w-4" />}
                        {evaluation.dimer.has_clashes ? 'Detected' : 'None'}
                      </span>
                    </div>
                    <div className="flex items-center justify-between">
                      <span className="text-sm text-muted-foreground">Chain Separability</span>
                      <span className={`flex items-center gap-1 text-sm font-medium ${
                        evaluation.dimer.separable ? 'text-green-600' : 'text-red-600'
                      }`}>
                        {evaluation.dimer.separable ? <Unlink className="h-4 w-4" /> : <Link className="h-4 w-4" />}
                        {evaluation.dimer.separable ? 'Separable' : 'Entangled'}
                      </span>
                    </div>
                    {evaluation.dimer.shape_complementarity !== undefined && (
                      <div className="flex items-center justify-between">
                        <span className="text-sm text-muted-foreground">Shape Complementarity</span>
                        <span className="text-sm font-medium text-foreground">
                          {(evaluation.dimer.shape_complementarity * 100).toFixed(0)}%
                        </span>
                      </div>
                    )}
                  </div>
                </div>

                {/* GNINA Score (if available) */}
                {evaluation.gnina_score && (
                  <div className="bg-card/60 rounded-lg p-3">
                    <div className="text-sm font-medium text-foreground mb-2">GNINA CNN Scoring</div>
                    <div className="grid grid-cols-2 gap-3">
                      <div>
                        <div className="text-xs text-muted-foreground">CNN Score</div>
                        <div className="text-lg font-bold text-primary">{evaluation.gnina_score.cnn_score.toFixed(2)}</div>
                      </div>
                      <div>
                        <div className="text-xs text-muted-foreground">CNN Affinity</div>
                        <div className="text-lg font-bold text-primary">{evaluation.gnina_score.cnn_affinity.toFixed(2)}</div>
                      </div>
                    </div>
                  </div>
                )}
              </div>
            )}

            {/* Chain A Tab */}
            {activeTab === 'chain_a' && evaluation.chain_a && (
              <div className="space-y-3">
                <div className="bg-primary/5 rounded-lg p-3">
                  <div className="flex items-center gap-2 mb-3">
                    <div className="w-8 h-8 rounded bg-primary text-primary-foreground flex items-center justify-center font-bold">A</div>
                    <div>
                      <div className="font-medium text-foreground">Chain A Metrics</div>
                      <div className="text-xs text-muted-foreground">Asymmetric binder (one-sided)</div>
                    </div>
                  </div>

                  <div className="grid grid-cols-2 gap-3">
                    <div className="bg-card/70 rounded p-2">
                      <div className="text-xs text-muted-foreground">Ligand Contacts</div>
                      <div className="text-xl font-bold text-primary">{evaluation.chain_a.contacts}</div>
                    </div>
                    {evaluation.chain_a.affinity !== undefined && (
                      <div className="bg-card/70 rounded p-2">
                        <div className="text-xs text-muted-foreground">Affinity</div>
                        <div className="text-xl font-bold text-primary">{evaluation.chain_a.affinity.toFixed(2)}</div>
                      </div>
                    )}
                    {evaluation.chain_a.hbonds !== undefined && (
                      <div className="bg-card/70 rounded p-2">
                        <div className="text-xs text-muted-foreground">H-Bonds</div>
                        <div className="text-xl font-bold text-primary">{evaluation.chain_a.hbonds}</div>
                      </div>
                    )}
                    {evaluation.chain_a.hydrophobic_contacts !== undefined && (
                      <div className="bg-card/70 rounded p-2">
                        <div className="text-xs text-muted-foreground">Hydrophobic</div>
                        <div className="text-xl font-bold text-primary">{evaluation.chain_a.hydrophobic_contacts}</div>
                      </div>
                    )}
                  </div>

                  {evaluation.chain_a.exposed_atoms && (
                    <div className="mt-3">
                      <div className="text-xs text-muted-foreground mb-1">Exposed Ligand Atoms</div>
                      <div className="flex flex-wrap gap-1">
                        {evaluation.chain_a.exposed_atoms.split(',').map((atom, idx) => (
                          <span key={idx} className="text-xs bg-primary/20 text-primary px-2 py-0.5 rounded">
                            {atom.trim()}
                          </span>
                        ))}
                      </div>
                    </div>
                  )}

                  {evaluation.chain_a.residue_contacts && evaluation.chain_a.residue_contacts.length > 0 && (
                    <div className="mt-3">
                      <div className="text-xs text-muted-foreground mb-1">Contacting Residues</div>
                      <div className="max-h-24 overflow-y-auto">
                        <table className="w-full text-xs">
                          <thead className="text-muted-foreground">
                            <tr>
                              <th className="text-left py-1">Residue</th>
                              <th className="text-right py-1">Distance (Å)</th>
                            </tr>
                          </thead>
                          <tbody>
                            {evaluation.chain_a.residue_contacts.map((contact, idx) => (
                              <tr key={idx} className="border-t border-border">
                                <td className="py-1">{contact.resname}{contact.resnum}</td>
                                <td className="text-right py-1">{contact.distance.toFixed(2)}</td>
                              </tr>
                            ))}
                          </tbody>
                        </table>
                      </div>
                    </div>
                  )}
                </div>
              </div>
            )}

            {/* Chain B Tab */}
            {activeTab === 'chain_b' && evaluation.chain_b && (
              <div className="space-y-3">
                <div className="bg-muted/50 rounded-lg p-3">
                  <div className="flex items-center gap-2 mb-3">
                    <div className="w-8 h-8 rounded bg-muted-foreground text-card flex items-center justify-center font-bold">B</div>
                    <div>
                      <div className="font-medium text-foreground">Chain B Metrics</div>
                      <div className="text-xs text-muted-foreground">Sequential binder (complementary)</div>
                    </div>
                  </div>

                  <div className="grid grid-cols-2 gap-3">
                    <div className="bg-card/70 rounded p-2">
                      <div className="text-xs text-muted-foreground">Ligand Contacts</div>
                      <div className="text-xl font-bold text-foreground">{evaluation.chain_b.contacts}</div>
                    </div>
                    {evaluation.chain_b.affinity !== undefined && (
                      <div className="bg-card/70 rounded p-2">
                        <div className="text-xs text-muted-foreground">Affinity</div>
                        <div className="text-xl font-bold text-foreground">{evaluation.chain_b.affinity.toFixed(2)}</div>
                      </div>
                    )}
                    {evaluation.chain_b.hbonds !== undefined && (
                      <div className="bg-card/70 rounded p-2">
                        <div className="text-xs text-muted-foreground">H-Bonds</div>
                        <div className="text-xl font-bold text-foreground">{evaluation.chain_b.hbonds}</div>
                      </div>
                    )}
                    {evaluation.chain_b.hydrophobic_contacts !== undefined && (
                      <div className="bg-card/70 rounded p-2">
                        <div className="text-xs text-muted-foreground">Hydrophobic</div>
                        <div className="text-xl font-bold text-foreground">{evaluation.chain_b.hydrophobic_contacts}</div>
                      </div>
                    )}
                  </div>

                  {evaluation.chain_b.exposed_atoms && (
                    <div className="mt-3">
                      <div className="text-xs text-muted-foreground mb-1">Exposed Ligand Atoms</div>
                      <div className="flex flex-wrap gap-1">
                        {evaluation.chain_b.exposed_atoms.split(',').map((atom, idx) => (
                          <span key={idx} className="text-xs bg-muted text-muted-foreground px-2 py-0.5 rounded">
                            {atom.trim()}
                          </span>
                        ))}
                      </div>
                    </div>
                  )}

                  {evaluation.chain_b.residue_contacts && evaluation.chain_b.residue_contacts.length > 0 && (
                    <div className="mt-3">
                      <div className="text-xs text-muted-foreground mb-1">Contacting Residues</div>
                      <div className="max-h-24 overflow-y-auto">
                        <table className="w-full text-xs">
                          <thead className="text-muted-foreground">
                            <tr>
                              <th className="text-left py-1">Residue</th>
                              <th className="text-right py-1">Distance (Å)</th>
                            </tr>
                          </thead>
                          <tbody>
                            {evaluation.chain_b.residue_contacts.map((contact, idx) => (
                              <tr key={idx} className="border-t border-border">
                                <td className="py-1">{contact.resname}{contact.resnum}</td>
                                <td className="text-right py-1">{contact.distance.toFixed(2)}</td>
                              </tr>
                            ))}
                          </tbody>
                        </table>
                      </div>
                    </div>
                  )}
                </div>
              </div>
            )}

            {/* Energy Tab */}
            {activeTab === 'energy' && evaluation.dimer?.energy_breakdown && (
              <div className="space-y-3">
                <div className="bg-card/60 rounded-lg p-3">
                  <div className="text-sm font-medium text-foreground mb-3">Energy Breakdown</div>

                  {/* Energy components */}
                  <div className="space-y-2">
                    {Object.entries(evaluation.dimer.energy_breakdown).map(([key, value]) => {
                      const labels: Record<string, string> = {
                        vdw: 'Van der Waals',
                        electrostatic: 'Electrostatic',
                        hbond: 'Hydrogen Bonds',
                        desolvation: 'Desolvation',
                      };
                      const isNegative = value < 0;
                      return (
                        <div key={key} className="flex items-center justify-between">
                          <span className="text-sm text-muted-foreground">{labels[key] || key}</span>
                          <div className="flex items-center gap-2">
                            <div className="w-24 h-2 bg-muted rounded-full overflow-hidden">
                              <div
                                className={`h-full rounded-full ${isNegative ? 'bg-green-500' : 'bg-red-500'}`}
                                style={{ width: `${Math.min(100, Math.abs(value) / 5 * 100)}%` }}
                              />
                            </div>
                            <span className={`text-sm font-medium w-16 text-right ${isNegative ? 'text-green-600' : 'text-red-600'}`}>
                              {value.toFixed(2)}
                            </span>
                          </div>
                        </div>
                      );
                    })}
                  </div>

                  {/* Total */}
                  <div className="mt-3 pt-3 border-t border-border flex items-center justify-between">
                    <span className="text-sm font-medium text-foreground">Total Binding Energy</span>
                    <span className={`text-lg font-bold ${evaluation.dimer.affinity < 0 ? 'text-green-600' : 'text-red-600'}`}>
                      {evaluation.dimer.affinity.toFixed(2)} kcal/mol
                    </span>
                  </div>
                </div>
              </div>
            )}

            {/* No energy data message */}
            {activeTab === 'energy' && !evaluation.dimer?.energy_breakdown && (
              <div className="flex items-center justify-center h-32 text-muted-foreground text-sm">
                <Info className="h-4 w-4 mr-2" />
                Energy breakdown not available for this design
              </div>
            )}
          </div>
        )}

        {/* Warnings */}
        {evaluation.warnings && evaluation.warnings.length > 0 && (
          <div className="bg-amber-50 border border-amber-200 rounded-lg p-3">
            <div className="flex items-center gap-2 text-amber-700 text-sm font-medium mb-1">
              <AlertTriangle className="h-4 w-4" />
              Warnings
            </div>
            <ul className="space-y-1">
              {evaluation.warnings.map((warning, idx) => (
                <li key={idx} className="text-xs text-amber-600 flex items-start gap-1">
                  <span className="text-amber-400">•</span>
                  {warning}
                </li>
              ))}
            </ul>
          </div>
        )}

        {/* Visual Summary (compact version for non-expanded) */}
        {!isExpanded && (
          <div className="bg-card/60 rounded-lg p-3">
            <div className="flex items-center justify-center gap-2">
              {/* Chain A */}
              <div className="flex flex-col items-center">
                <div className={`w-10 h-10 rounded flex items-center justify-center text-sm font-medium ${
                  evaluation.chain_a ? 'bg-primary/20 text-primary' : 'bg-muted text-muted-foreground'
                }`}>
                  A
                </div>
                <span className="text-[10px] text-muted-foreground mt-1">
                  {evaluation.dimer?.contacts_a || evaluation.chain_a?.contacts || 0}
                </span>
              </div>

              {/* Connection */}
              <div className="w-4 h-0.5 bg-primary/30"></div>

              {/* Ligand */}
              <div className="flex flex-col items-center">
                <div className="w-10 h-10 rounded bg-muted flex items-center justify-center">
                  <Hexagon className="h-5 w-5 text-primary" />
                </div>
                <span className="text-[10px] text-muted-foreground mt-1">Ligand</span>
              </div>

              {/* Connection */}
              <div className="w-4 h-0.5 bg-muted-foreground/30"></div>

              {/* Chain B */}
              <div className="flex flex-col items-center">
                <div className={`w-10 h-10 rounded flex items-center justify-center text-sm font-medium ${
                  evaluation.chain_b || evaluation.dimer ? 'bg-muted text-muted-foreground' : 'bg-muted text-muted-foreground'
                }`}>
                  B
                </div>
                <span className="text-[10px] text-muted-foreground mt-1">
                  {evaluation.dimer?.contacts_b || evaluation.chain_b?.contacts || 0}
                </span>
              </div>
            </div>

            {/* Separability indicator */}
            {evaluation.dimer && (
              <div className="mt-2 flex items-center justify-center gap-1 text-xs">
                {evaluation.dimer.separable
                  ? <CheckCircle className="h-4 w-4 text-green-500" />
                  : <XCircle className="h-4 w-4 text-red-500" />
                }
                <span className={evaluation.dimer.separable ? 'text-green-600' : 'text-red-600'}>
                  {evaluation.dimer.separable ? 'Chains can separate' : 'Chains entangled'}
                </span>
              </div>
            )}
          </div>
        )}
      </div>
    </div>
  );
}
