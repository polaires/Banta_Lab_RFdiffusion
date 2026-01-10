'use client';

import { useState } from 'react';
import {
  ChevronDown,
  ChevronRight,
  Atom,
  Pill,
  Focus,
  Info,
  CheckCircle,
  AlertCircle,
  HelpCircle
} from 'lucide-react';
import type { MetalCoordination } from '@/lib/metalAnalysis';
import type { LigandAnalysisResult, LigandData } from '@/lib/ligandAnalysis';
import { getInteractionColor, getInteractionLabel } from '@/lib/ligandAnalysis';

interface BindingSitePanelProps {
  metalCoordination: MetalCoordination[] | null;
  ligandData: LigandAnalysisResult | null;
  focusedMetalIndex: number | null;
  focusedLigandIndex: number | null;
  onFocusMetal: (index: number | null) => void;
  onFocusLigand: (index: number | null) => void;
  loading?: boolean;
}

// Binding site type badge - light theme
function BindingSiteTypeBadge({ type }: { type: 'functional' | 'crystal_artifact' | 'uncertain' }) {
  switch (type) {
    case 'functional':
      return (
        <span className="inline-flex items-center gap-1 px-2 py-0.5 rounded text-xs bg-emerald-50 text-emerald-700 border border-emerald-200">
          <CheckCircle className="w-3 h-3" />
          Functional
        </span>
      );
    case 'crystal_artifact':
      return (
        <span className="inline-flex items-center gap-1 px-2 py-0.5 rounded text-xs bg-red-50 text-red-700 border border-red-200">
          <AlertCircle className="w-3 h-3" />
          Crystal Artifact
        </span>
      );
    default:
      return (
        <span className="inline-flex items-center gap-1 px-2 py-0.5 rounded text-xs bg-amber-50 text-amber-700 border border-amber-200">
          <HelpCircle className="w-3 h-3" />
          Uncertain
        </span>
      );
  }
}

// Distortion level badge for geometry - light theme
function DistortionBadge({ level }: { level: 'ideal' | 'low' | 'moderate' | 'high' | 'severe' }) {
  const colors = {
    ideal: 'bg-emerald-50 text-emerald-700 border-emerald-200',
    low: 'bg-blue-50 text-blue-700 border-blue-200',
    moderate: 'bg-amber-50 text-amber-700 border-amber-200',
    high: 'bg-orange-50 text-orange-700 border-orange-200',
    severe: 'bg-red-50 text-red-700 border-red-200'
  };

  return (
    <span className={`px-1.5 py-0.5 rounded text-xs border ${colors[level]}`}>
      {level}
    </span>
  );
}

// Metal card component - light theme with expanded details by default
function MetalCard({
  metal,
  index,
  isFocused,
}: {
  metal: MetalCoordination;
  index: number;
  isFocused: boolean;
  onFocus: (index: number | null) => void;
}) {
  const [expanded, setExpanded] = useState(true); // Expanded by default

  return (
    <div
      className={`rounded-lg border shadow-sm transition-all ${
        isFocused
          ? 'border-blue-400 bg-blue-50 ring-2 ring-blue-200'
          : 'border-slate-200 bg-white'
      }`}
    >
      {/* Header */}
      <div className="p-3">
        <div className="flex items-center justify-between mb-2">
          <div className="flex items-center gap-2">
            <button
              onClick={() => setExpanded(!expanded)}
              className="p-0.5 hover:bg-slate-200 rounded text-slate-500"
            >
              {expanded ? <ChevronDown className="w-4 h-4" /> : <ChevronRight className="w-4 h-4" />}
            </button>
            <Atom className="w-4 h-4 text-purple-600" />
            <span className="font-semibold text-slate-800">{metal.element}</span>
            <span className="text-xs text-slate-500">
              Chain {metal.chainId}, Res {metal.resSeq}
            </span>
          </div>
          {isFocused && (
            <div className="p-1.5 rounded bg-blue-600 text-white">
              <Focus className="w-4 h-4" />
            </div>
          )}
        </div>

        {/* Quick info */}
        <div className="flex flex-wrap gap-1.5 mb-2">
          <BindingSiteTypeBadge type={metal.bindingSiteType} />
          {metal.geometry && (
            <>
              <span className="px-2 py-0.5 rounded text-xs bg-slate-100 text-slate-600 border border-slate-200">
                {metal.geometry.geometryType}
              </span>
              <span className="px-2 py-0.5 rounded text-xs bg-purple-50 text-purple-700 border border-purple-200">
                CN: {metal.geometry.coordinationNumber}
              </span>
              <DistortionBadge level={metal.geometry.distortion} />
            </>
          )}
        </div>

        {/* Reason */}
        <p className="text-xs text-slate-500">{metal.bindingSiteReason}</p>
      </div>

      {/* Expanded content - coordinating atoms table */}
      {expanded && (
        <div className="border-t border-slate-100 p-3 space-y-3 bg-slate-50/50">
          {/* Coordinating atoms - table-like layout */}
          <div>
            <h4 className="text-xs font-semibold text-slate-600 mb-2 uppercase tracking-wide">
              Coordinating Atoms ({metal.coordinating.length})
            </h4>
            <div className="bg-white rounded border border-slate-200 divide-y divide-slate-100">
              {metal.coordinating.slice(0, 10).map((coord, i) => (
                <div key={i} className="flex items-center justify-between px-2.5 py-1.5 text-xs hover:bg-slate-50">
                  <div className="flex items-center gap-2">
                    <span className={`font-medium ${coord.isWater ? 'text-cyan-600' : 'text-slate-700'}`}>
                      {coord.residue}
                    </span>
                    <span className="text-slate-500">{coord.atom}</span>
                    {coord.isWater && (
                      <span className="px-1.5 py-0.5 rounded text-[10px] bg-cyan-50 text-cyan-600 border border-cyan-200">
                        water
                      </span>
                    )}
                  </div>
                  <span className="font-mono text-slate-600">{coord.distance.toFixed(2)} A</span>
                </div>
              ))}
              {metal.coordinating.length > 10 && (
                <div className="px-2.5 py-1.5 text-xs text-slate-400 text-center">
                  ... and {metal.coordinating.length - 10} more
                </div>
              )}
            </div>
          </div>

          {/* Hydration analysis */}
          <div className="bg-white rounded border border-slate-200 p-2.5">
            <h4 className="text-xs font-semibold text-slate-600 mb-1.5">Hydration State</h4>
            <p className="text-xs text-slate-500 mb-2">{metal.hydrationAnalysis.hydrationNote}</p>
            <div className="flex gap-3 text-xs">
              <span className="flex items-center gap-1">
                <span className="w-2 h-2 rounded-full bg-cyan-400"></span>
                <span className="text-slate-600">{metal.hydrationAnalysis.waterCount} water</span>
              </span>
              <span className="flex items-center gap-1">
                <span className="w-2 h-2 rounded-full bg-emerald-400"></span>
                <span className="text-slate-600">{metal.hydrationAnalysis.proteinLigandCount} protein</span>
              </span>
            </div>
          </div>

          {/* Geometry angles (if available) */}
          {metal.geometry && metal.geometry.angles.length > 0 && (
            <div className="bg-white rounded border border-slate-200 p-2.5">
              <h4 className="text-xs font-semibold text-slate-600 mb-1.5">
                L-M-L Angles <span className="font-normal text-slate-400">(avg: {metal.geometry.avgAngle.toFixed(1)}deg)</span>
              </h4>
              <div className="flex flex-wrap gap-1.5">
                {metal.geometry.angles.slice(0, 6).map((a, i) => (
                  <span key={i} className="px-2 py-0.5 rounded text-xs bg-slate-100 text-slate-600 font-mono">
                    {a.angle.toFixed(1)}°
                  </span>
                ))}
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  );
}

// Ligand card component - light theme with expanded details by default
function LigandCard({
  ligand,
  index,
  isFocused,
}: {
  ligand: LigandData;
  index: number;
  isFocused: boolean;
  onFocus: (index: number | null) => void;
}) {
  const [expanded, setExpanded] = useState(true); // Expanded by default

  return (
    <div
      className={`rounded-lg border shadow-sm transition-all ${
        isFocused
          ? 'border-blue-400 bg-blue-50 ring-2 ring-blue-200'
          : 'border-slate-200 bg-white'
      }`}
    >
      {/* Header */}
      <div className="p-3">
        <div className="flex items-center justify-between mb-2">
          <div className="flex items-center gap-2">
            <button
              onClick={() => setExpanded(!expanded)}
              className="p-0.5 hover:bg-slate-200 rounded text-slate-500"
            >
              {expanded ? <ChevronDown className="w-4 h-4" /> : <ChevronRight className="w-4 h-4" />}
            </button>
            <Pill className="w-4 h-4 text-emerald-600" />
            <span className="font-semibold text-slate-800">{ligand.name}</span>
            <span className="text-xs text-slate-500">
              Chain {ligand.chainId}, {ligand.atoms} atoms
            </span>
          </div>
          {isFocused && (
            <div className="p-1.5 rounded bg-blue-600 text-white">
              <Focus className="w-4 h-4" />
            </div>
          )}
        </div>

        {/* Quick info */}
        <div className="flex flex-wrap gap-1.5 mb-2">
          <BindingSiteTypeBadge type={ligand.bindingSiteType} />
          <span className="px-2 py-0.5 rounded text-xs bg-slate-100 text-slate-600 border border-slate-200">
            {ligand.contacts.length} contacts
          </span>
        </div>

        {/* Reason */}
        <p className="text-xs text-slate-500">{ligand.bindingSiteReason}</p>
      </div>

      {/* Expanded content */}
      {expanded && (
        <div className="border-t border-slate-100 p-3 space-y-3 bg-slate-50/50">
          {/* Protein contacts - table-like layout */}
          <div>
            <h4 className="text-xs font-semibold text-slate-600 mb-2 uppercase tracking-wide">
              Protein Contacts ({ligand.contacts.length})
            </h4>
            <div className="bg-white rounded border border-slate-200 divide-y divide-slate-100">
              {ligand.contacts.slice(0, 15).map((contact, i) => (
                <div key={i} className="flex items-center justify-between px-2.5 py-1.5 text-xs hover:bg-slate-50">
                  <div className="flex items-center gap-2">
                    <span
                      className="w-2.5 h-2.5 rounded-full flex-shrink-0"
                      style={{ backgroundColor: getInteractionColor(contact.interactionType) }}
                    />
                    <span className="font-medium text-slate-700">
                      {contact.residue}
                    </span>
                    <span className="text-slate-500">{contact.atom}</span>
                    <span className="px-1.5 py-0.5 rounded text-[10px] bg-slate-100 text-slate-500">
                      {getInteractionLabel(contact.interactionType)}
                    </span>
                  </div>
                  <span className="font-mono text-slate-600">{contact.distance.toFixed(2)} A</span>
                </div>
              ))}
              {ligand.contacts.length > 15 && (
                <div className="px-2.5 py-1.5 text-xs text-slate-400 text-center">
                  ... and {ligand.contacts.length - 15} more
                </div>
              )}
            </div>
          </div>

          {/* Interaction type legend */}
          <div className="bg-white rounded border border-slate-200 p-2.5">
            <h4 className="text-xs font-semibold text-slate-600 mb-2">Interaction Types</h4>
            <div className="flex flex-wrap gap-3 text-xs text-slate-600">
              <span className="flex items-center gap-1.5">
                <span className="w-2.5 h-2.5 rounded-full" style={{ backgroundColor: getInteractionColor('hydrogen_bond') }} />
                H-bond
              </span>
              <span className="flex items-center gap-1.5">
                <span className="w-2.5 h-2.5 rounded-full" style={{ backgroundColor: getInteractionColor('hydrophobic') }} />
                Hydrophobic
              </span>
              <span className="flex items-center gap-1.5">
                <span className="w-2.5 h-2.5 rounded-full" style={{ backgroundColor: getInteractionColor('salt_bridge') }} />
                Salt bridge
              </span>
              <span className="flex items-center gap-1.5">
                <span className="w-2.5 h-2.5 rounded-full" style={{ backgroundColor: getInteractionColor('pi_stacking') }} />
                Pi-stack
              </span>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}

export function BindingSitePanel({
  metalCoordination,
  ligandData,
  focusedMetalIndex,
  focusedLigandIndex,
  onFocusMetal,
  onFocusLigand,
  loading = false
}: BindingSitePanelProps) {
  const [metalsExpanded, setMetalsExpanded] = useState(true);
  const [ligandsExpanded, setLigandsExpanded] = useState(true);

  const hasMetals = metalCoordination && metalCoordination.length > 0;
  const hasLigands = ligandData && ligandData.ligandCount > 0;

  if (loading) {
    return (
      <div className="w-80 bg-slate-50 border-l border-slate-200 p-4">
        <div className="flex items-center gap-2 text-slate-500">
          <div className="w-4 h-4 border-2 border-blue-500 border-t-transparent rounded-full animate-spin" />
          <span className="text-sm font-medium">Analyzing structure...</span>
        </div>
      </div>
    );
  }

  if (!hasMetals && !hasLigands) {
    return (
      <div className="w-80 bg-slate-50 border-l border-slate-200 p-4">
        <div className="flex items-center gap-2 text-slate-500">
          <Info className="w-4 h-4" />
          <span className="text-sm font-medium">No binding sites detected</span>
        </div>
        <p className="text-xs text-slate-400 mt-2">
          Load a structure with metal ions or ligands to see binding site analysis.
        </p>
      </div>
    );
  }

  return (
    <div className="w-80 bg-slate-50 border-l border-slate-200 overflow-y-auto h-[564px]">
      {/* Header */}
      <div className="p-4 border-b border-slate-200 bg-white sticky top-0 z-10">
        <h2 className="font-semibold text-slate-800">Binding Site Analysis</h2>
        <p className="text-xs text-slate-500 mt-0.5">
          {hasMetals && `${metalCoordination!.length} metal${metalCoordination!.length !== 1 ? 's' : ''}`}
          {hasMetals && hasLigands && ' · '}
          {hasLigands && `${ligandData!.ligandCount} ligand${ligandData!.ligandCount !== 1 ? 's' : ''}`}
        </p>
      </div>

      <div className="p-4 space-y-4">
        {/* Metal ions section */}
        {hasMetals && (
          <div>
            <button
              onClick={() => setMetalsExpanded(!metalsExpanded)}
              className="flex items-center gap-2 w-full text-left mb-3 group"
            >
              <span className="text-slate-400 group-hover:text-slate-600 transition-colors">
                {metalsExpanded ? <ChevronDown className="w-4 h-4" /> : <ChevronRight className="w-4 h-4" />}
              </span>
              <Atom className="w-4 h-4 text-purple-600" />
              <span className="font-semibold text-slate-700">Metal Ions</span>
              <span className="text-xs text-slate-400 bg-slate-200 px-1.5 py-0.5 rounded">
                {metalCoordination!.length}
              </span>
            </button>

            {metalsExpanded && (
              <div className="space-y-3">
                {metalCoordination!.map((metal, i) => (
                  <MetalCard
                    key={`${metal.chainId}-${metal.resSeq}-${i}`}
                    metal={metal}
                    index={i}
                    isFocused={focusedMetalIndex === i}
                    onFocus={onFocusMetal}
                  />
                ))}
              </div>
            )}
          </div>
        )}

        {/* Ligands section */}
        {hasLigands && (
          <div>
            <button
              onClick={() => setLigandsExpanded(!ligandsExpanded)}
              className="flex items-center gap-2 w-full text-left mb-3 group"
            >
              <span className="text-slate-400 group-hover:text-slate-600 transition-colors">
                {ligandsExpanded ? <ChevronDown className="w-4 h-4" /> : <ChevronRight className="w-4 h-4" />}
              </span>
              <Pill className="w-4 h-4 text-emerald-600" />
              <span className="font-semibold text-slate-700">Ligands</span>
              <span className="text-xs text-slate-400 bg-slate-200 px-1.5 py-0.5 rounded">
                {ligandData!.ligandCount}
              </span>
            </button>

            {ligandsExpanded && (
              <div className="space-y-3">
                {ligandData!.ligandDetails.map((ligand, i) => (
                  <LigandCard
                    key={`${ligand.chainId}-${ligand.resSeq}-${i}`}
                    ligand={ligand}
                    index={i}
                    isFocused={focusedLigandIndex === i}
                    onFocus={onFocusLigand}
                  />
                ))}
              </div>
            )}
          </div>
        )}
      </div>
    </div>
  );
}

export default BindingSitePanel;
