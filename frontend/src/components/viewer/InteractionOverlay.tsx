'use client';

import { useState, useMemo } from 'react';

/**
 * Types for protein-ligand interaction data from the backend.
 */
export interface HydrogenBond {
  ligand_atom: string;
  protein_residue: string;
  protein_chain: string;
  protein_atom: string;
  distance: number;
  angle?: number | null;
  donor_coords?: [number, number, number] | null;
  acceptor_coords?: [number, number, number] | null;
  type: 'ligand_donor' | 'ligand_acceptor' | string;
}

export interface HydrophobicContact {
  ligand_atom: string;
  protein_residue: string;
  protein_chain: string;
  protein_atom: string;
  distance: number;
  coords?: [number, number, number] | null;
}

export interface PiStacking {
  ligand_ring: string;
  protein_residue: string;
  protein_chain: string;
  type: 'face-to-face' | 'edge-to-face';
  distance: number;
  angle?: number | null;
}

export interface SaltBridge {
  ligand_group: string;
  protein_residue: string;
  protein_chain: string;
  distance: number;
  ligand_charge: 'positive' | 'negative';
  protein_charge: 'positive' | 'negative';
}

export interface InteractionData {
  interactions: {
    hbonds: number;
    hydrophobic: number;
    pi_stacking: number;
    salt_bridges: number;
    halogen_bonds: number;
    total: number;
  };
  key_residues: string[];
  details: {
    hydrogen_bonds: HydrogenBond[];
    hydrophobic_contacts: HydrophobicContact[];
    pi_stacking?: PiStacking[];
    salt_bridges?: SaltBridge[];
  };
  visualization?: {
    lines: Array<{
      start: [number, number, number];
      end: [number, number, number];
      color: string;
      dashed: boolean;
      label?: string;
    }>;
    spheres: Array<{
      center: [number, number, number];
      radius: number;
      color: string;
      opacity: number;
    }>;
    highlights: string[];
  };
  analysis_method: 'plip' | 'distance_based';
  status: string;
}

interface InteractionOverlayProps {
  interactions: InteractionData | null;
  className?: string;
  showDetails?: boolean;
  onResidueClick?: (residue: string) => void;
}

/**
 * Badge component for interaction counts.
 */
function InteractionBadge({
  count,
  label,
  color
}: {
  count: number;
  label: string;
  color: string;
}) {
  return (
    <div className="flex flex-col items-center p-2 rounded-lg bg-gray-50">
      <span className={`text-lg font-bold ${color}`}>{count}</span>
      <span className="text-xs text-gray-500">{label}</span>
    </div>
  );
}

/**
 * InteractionOverlay displays protein-ligand interaction analysis results.
 *
 * Shows:
 * - Interaction summary (counts by type)
 * - Key binding residues
 * - Detailed interaction list (expandable)
 */
export function InteractionOverlay({
  interactions,
  className = '',
  showDetails = false,
  onResidueClick,
}: InteractionOverlayProps) {
  const [expandedSection, setExpandedSection] = useState<string | null>(
    showDetails ? 'hbonds' : null
  );

  // Memoize sorted H-bonds by distance
  const sortedHBonds = useMemo(() => {
    if (!interactions?.details?.hydrogen_bonds) return [];
    return [...interactions.details.hydrogen_bonds].sort((a, b) => a.distance - b.distance);
  }, [interactions?.details?.hydrogen_bonds]);

  // Memoize sorted hydrophobic contacts
  const sortedHydrophobic = useMemo(() => {
    if (!interactions?.details?.hydrophobic_contacts) return [];
    return [...interactions.details.hydrophobic_contacts].sort((a, b) => a.distance - b.distance);
  }, [interactions?.details?.hydrophobic_contacts]);

  if (!interactions) {
    return (
      <div className={`p-4 bg-gray-50 rounded-lg ${className}`}>
        <p className="text-gray-500 text-sm">No interaction data available</p>
      </div>
    );
  }

  const { interactions: counts, key_residues, analysis_method } = interactions;

  return (
    <div className={`bg-white rounded-lg shadow-sm border border-gray-200 ${className}`}>
      {/* Header */}
      <div className="p-3 border-b border-gray-100">
        <div className="flex items-center justify-between">
          <h3 className="font-medium text-gray-800">Interaction Profile</h3>
          <span className="text-xs px-2 py-0.5 bg-gray-100 text-gray-600 rounded">
            {analysis_method === 'plip' ? 'PLIP' : 'Distance-based'}
          </span>
        </div>
      </div>

      {/* Summary counts */}
      <div className="p-3 border-b border-gray-100">
        <div className="grid grid-cols-5 gap-2">
          <InteractionBadge count={counts.hbonds} label="H-bonds" color="text-blue-600" />
          <InteractionBadge count={counts.hydrophobic} label="Hydrophobic" color="text-yellow-600" />
          <InteractionBadge count={counts.pi_stacking} label="Pi-stack" color="text-purple-600" />
          <InteractionBadge count={counts.salt_bridges} label="Salt bridge" color="text-red-600" />
          <InteractionBadge count={counts.total} label="Total" color="text-gray-800" />
        </div>
      </div>

      {/* Key residues */}
      {key_residues.length > 0 && (
        <div className="p-3 border-b border-gray-100">
          <h4 className="text-xs font-medium text-gray-500 mb-2">KEY BINDING RESIDUES</h4>
          <div className="flex flex-wrap gap-1">
            {key_residues.slice(0, 8).map((residue) => (
              <button
                key={residue}
                onClick={() => onResidueClick?.(residue)}
                className="px-2 py-0.5 text-xs bg-emerald-50 text-emerald-700 rounded hover:bg-emerald-100 transition-colors"
              >
                {residue}
              </button>
            ))}
            {key_residues.length > 8 && (
              <span className="px-2 py-0.5 text-xs text-gray-400">
                +{key_residues.length - 8} more
              </span>
            )}
          </div>
        </div>
      )}

      {/* Expandable details */}
      <div className="divide-y divide-gray-100">
        {/* Hydrogen Bonds */}
        {counts.hbonds > 0 && (
          <div>
            <button
              onClick={() => setExpandedSection(expandedSection === 'hbonds' ? null : 'hbonds')}
              className="w-full p-3 flex items-center justify-between hover:bg-gray-50 transition-colors"
            >
              <span className="text-sm font-medium text-gray-700">
                Hydrogen Bonds ({counts.hbonds})
              </span>
              <svg
                className={`w-4 h-4 text-gray-400 transition-transform ${
                  expandedSection === 'hbonds' ? 'rotate-180' : ''
                }`}
                fill="none"
                viewBox="0 0 24 24"
                stroke="currentColor"
              >
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
              </svg>
            </button>
            {expandedSection === 'hbonds' && (
              <div className="px-3 pb-3">
                <div className="bg-gray-50 rounded-lg overflow-hidden">
                  <table className="w-full text-xs">
                    <thead className="bg-gray-100">
                      <tr>
                        <th className="px-2 py-1 text-left text-gray-600">Protein</th>
                        <th className="px-2 py-1 text-left text-gray-600">Ligand</th>
                        <th className="px-2 py-1 text-right text-gray-600">Distance</th>
                        <th className="px-2 py-1 text-right text-gray-600">Angle</th>
                      </tr>
                    </thead>
                    <tbody className="divide-y divide-gray-100">
                      {sortedHBonds.slice(0, 10).map((hb, i) => (
                        <tr key={i} className="hover:bg-blue-50">
                          <td className="px-2 py-1 text-gray-700">
                            {hb.protein_chain}:{hb.protein_residue}.{hb.protein_atom}
                          </td>
                          <td className="px-2 py-1 text-gray-700">{hb.ligand_atom}</td>
                          <td className="px-2 py-1 text-right text-blue-600">
                            {hb.distance.toFixed(2)} A
                          </td>
                          <td className="px-2 py-1 text-right text-gray-500">
                            {hb.angle ? `${hb.angle.toFixed(0)}deg` : '-'}
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                  {sortedHBonds.length > 10 && (
                    <div className="px-2 py-1 text-xs text-gray-400 text-center bg-gray-100">
                      +{sortedHBonds.length - 10} more
                    </div>
                  )}
                </div>
              </div>
            )}
          </div>
        )}

        {/* Hydrophobic Contacts */}
        {counts.hydrophobic > 0 && (
          <div>
            <button
              onClick={() => setExpandedSection(expandedSection === 'hydrophobic' ? null : 'hydrophobic')}
              className="w-full p-3 flex items-center justify-between hover:bg-gray-50 transition-colors"
            >
              <span className="text-sm font-medium text-gray-700">
                Hydrophobic Contacts ({counts.hydrophobic})
              </span>
              <svg
                className={`w-4 h-4 text-gray-400 transition-transform ${
                  expandedSection === 'hydrophobic' ? 'rotate-180' : ''
                }`}
                fill="none"
                viewBox="0 0 24 24"
                stroke="currentColor"
              >
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
              </svg>
            </button>
            {expandedSection === 'hydrophobic' && (
              <div className="px-3 pb-3">
                <div className="bg-gray-50 rounded-lg overflow-hidden">
                  <table className="w-full text-xs">
                    <thead className="bg-gray-100">
                      <tr>
                        <th className="px-2 py-1 text-left text-gray-600">Protein</th>
                        <th className="px-2 py-1 text-left text-gray-600">Ligand</th>
                        <th className="px-2 py-1 text-right text-gray-600">Distance</th>
                      </tr>
                    </thead>
                    <tbody className="divide-y divide-gray-100">
                      {sortedHydrophobic.slice(0, 10).map((hc, i) => (
                        <tr key={i} className="hover:bg-yellow-50">
                          <td className="px-2 py-1 text-gray-700">
                            {hc.protein_chain}:{hc.protein_residue}.{hc.protein_atom}
                          </td>
                          <td className="px-2 py-1 text-gray-700">{hc.ligand_atom}</td>
                          <td className="px-2 py-1 text-right text-yellow-600">
                            {hc.distance.toFixed(2)} A
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                  {sortedHydrophobic.length > 10 && (
                    <div className="px-2 py-1 text-xs text-gray-400 text-center bg-gray-100">
                      +{sortedHydrophobic.length - 10} more
                    </div>
                  )}
                </div>
              </div>
            )}
          </div>
        )}

        {/* Pi-Stacking */}
        {counts.pi_stacking > 0 && interactions.details.pi_stacking && (
          <div>
            <button
              onClick={() => setExpandedSection(expandedSection === 'pi_stacking' ? null : 'pi_stacking')}
              className="w-full p-3 flex items-center justify-between hover:bg-gray-50 transition-colors"
            >
              <span className="text-sm font-medium text-gray-700">
                Pi-Stacking ({counts.pi_stacking})
              </span>
              <svg
                className={`w-4 h-4 text-gray-400 transition-transform ${
                  expandedSection === 'pi_stacking' ? 'rotate-180' : ''
                }`}
                fill="none"
                viewBox="0 0 24 24"
                stroke="currentColor"
              >
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
              </svg>
            </button>
            {expandedSection === 'pi_stacking' && (
              <div className="px-3 pb-3">
                <div className="bg-gray-50 rounded-lg p-2">
                  {interactions.details.pi_stacking.map((ps, i) => (
                    <div
                      key={i}
                      className="flex items-center justify-between py-1 text-xs"
                    >
                      <span className="text-gray-700">
                        {ps.protein_chain}:{ps.protein_residue}
                      </span>
                      <span className="text-purple-600">
                        {ps.type} ({ps.distance.toFixed(2)} A)
                      </span>
                    </div>
                  ))}
                </div>
              </div>
            )}
          </div>
        )}
      </div>
    </div>
  );
}

export default InteractionOverlay;
