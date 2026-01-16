'use client';

import { useState, useMemo } from 'react';
import {
  ChevronDown,
  ChevronRight,
  Download,
  FileJson,
  FileSpreadsheet,
  Target,
  Focus,
  CheckCircle,
  AlertCircle,
  HelpCircle,
  Lightbulb,
  Eye,
  EyeOff,
  Layers,
} from 'lucide-react';
import { useStore } from '@/lib/store';
import {
  extractPharmacophoreFeatures,
  classifyShellResidues,
  calculateBinderQualityScore,
  type LigandData,
  type QualityRating,
  type PharmacophoreType,
} from '@/lib/ligandAnalysis';
import {
  downloadJSON,
  downloadCSV,
  downloadRFdiffusionHotspots,
} from '@/lib/ligandExport';

interface LigandBindingPanelProps {
  onToggle3D?: (show: boolean) => void;
  show3D?: boolean;
}

// Quality rating color styles
function getRatingColorClasses(rating: QualityRating): string {
  switch (rating) {
    case 'EXCELLENT':
      return 'bg-emerald-100 text-emerald-700 border-emerald-200';
    case 'GOOD':
      return 'bg-blue-100 text-blue-700 border-blue-200';
    case 'MODERATE':
      return 'bg-amber-100 text-amber-700 border-amber-200';
    case 'POOR':
      return 'bg-red-100 text-red-700 border-red-200';
    default:
      return 'bg-slate-100 text-slate-700 border-slate-200';
  }
}

// Score bar color for quality score
function getScoreBarColor(score: number): string {
  if (score >= 70) return 'bg-emerald-500';
  if (score >= 50) return 'bg-blue-500';
  if (score >= 30) return 'bg-amber-500';
  return 'bg-red-500';
}

// Binding site type badge
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
        <span className="inline-flex items-center gap-1 px-2 py-0.5 rounded text-xs bg-slate-100 text-slate-600 border border-slate-200">
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

// Pharmacophore type color
function getPharmacophoreColor(type: PharmacophoreType): string {
  switch (type) {
    case 'donor':
      return 'bg-blue-100 text-blue-700';
    case 'acceptor':
      return 'bg-red-100 text-red-700';
    case 'aromatic':
      return 'bg-purple-100 text-purple-700';
    case 'hydrophobic':
      return 'bg-yellow-100 text-yellow-700';
    case 'positive':
      return 'bg-cyan-100 text-cyan-700';
    case 'negative':
      return 'bg-orange-100 text-orange-700';
    default:
      return 'bg-slate-100 text-slate-700';
  }
}

// Ligand card with expanded analysis
function LigandCard({
  ligand,
  index,
  isFocused,
  onFocus,
  onToggle3D,
  show3D,
}: {
  ligand: LigandData;
  index: number;
  isFocused: boolean;
  onFocus: (index: number | null) => void;
  onToggle3D?: (show: boolean) => void;
  show3D?: boolean;
}) {
  // Compute analysis for this ligand
  const quality = useMemo(() => calculateBinderQualityScore(ligand.contacts), [ligand.contacts]);
  const shell = useMemo(() => classifyShellResidues(ligand.contacts), [ligand.contacts]);
  const pharmacophore = useMemo(() => extractPharmacophoreFeatures(ligand.contacts), [ligand.contacts]);

  // Count pharmacophore features by type
  const pharmacophoreCounts = useMemo(() => {
    const counts: Record<PharmacophoreType, number> = {
      donor: 0,
      acceptor: 0,
      aromatic: 0,
      hydrophobic: 0,
      positive: 0,
      negative: 0,
    };
    for (const feature of pharmacophore) {
      counts[feature.type]++;
    }
    return counts;
  }, [pharmacophore]);

  return (
    <div
      className={`rounded-lg border shadow-sm transition-all cursor-pointer ${
        isFocused
          ? 'border-blue-400 bg-blue-50 ring-2 ring-blue-200'
          : 'border-slate-200 bg-white hover:border-slate-300'
      }`}
      onClick={() => onFocus(isFocused ? null : index)}
    >
      {/* Header - always visible */}
      <div className="p-3">
        <div className="flex items-center justify-between mb-2">
          <div className="flex items-center gap-2">
            <BindingSiteTypeBadge type={ligand.bindingSiteType} />
            <span className="font-semibold text-slate-800">{ligand.name}</span>
          </div>
          {isFocused && (
            <div className="p-1.5 rounded bg-blue-600 text-white">
              <Focus className="w-4 h-4" />
            </div>
          )}
        </div>

        <div className="flex items-center gap-3 text-xs text-slate-500">
          <span>{ligand.chainId}:{ligand.resSeq}</span>
          <span className="px-2 py-0.5 rounded bg-slate-100 text-slate-600 border border-slate-200">
            {ligand.contacts.length} contacts
          </span>
        </div>
      </div>

      {/* Expanded analysis - only when focused */}
      {isFocused && (
        <div className="border-t border-slate-100 p-3 space-y-4 bg-slate-50/50">
          {/* Binder Quality Score */}
          <div className="bg-white rounded-lg border border-slate-200 p-3">
            <div className="flex items-center justify-between mb-2">
              <h4 className="text-xs font-semibold text-slate-600 uppercase tracking-wide">
                Binder Quality Score
              </h4>
              <span className={`px-2 py-0.5 rounded text-xs font-medium border ${getRatingColorClasses(quality.rating)}`}>
                {quality.rating}
              </span>
            </div>

            {/* Score bar */}
            <div className="mb-3">
              <div className="flex items-center justify-between mb-1">
                <span className="text-2xl font-bold text-slate-800">{quality.score}</span>
                <span className="text-xs text-slate-500">/ 100</span>
              </div>
              <div className="h-2 bg-slate-100 rounded-full overflow-hidden">
                <div
                  className={`h-full ${getScoreBarColor(quality.score)} transition-all`}
                  style={{ width: `${quality.score}%` }}
                />
              </div>
            </div>

            {/* Score breakdown grid */}
            <div className="grid grid-cols-5 gap-2 text-center">
              <div className="p-2 bg-slate-50 rounded">
                <div className="text-sm font-bold text-blue-600">{quality.breakdown.hbondScore}</div>
                <div className="text-[10px] text-slate-500">H-bond</div>
              </div>
              <div className="p-2 bg-slate-50 rounded">
                <div className="text-sm font-bold text-green-600">{quality.breakdown.hydrophobicScore}</div>
                <div className="text-[10px] text-slate-500">Hydro</div>
              </div>
              <div className="p-2 bg-slate-50 rounded">
                <div className="text-sm font-bold text-purple-600">{quality.breakdown.piStackScore}</div>
                <div className="text-[10px] text-slate-500">Pi-stack</div>
              </div>
              <div className="p-2 bg-slate-50 rounded">
                <div className="text-sm font-bold text-red-600">{quality.breakdown.saltBridgeScore}</div>
                <div className="text-[10px] text-slate-500">Salt</div>
              </div>
              <div className="p-2 bg-slate-50 rounded">
                <div className="text-sm font-bold text-slate-600">{quality.breakdown.burialnessScore}</div>
                <div className="text-[10px] text-slate-500">Burial</div>
              </div>
            </div>
          </div>

          {/* Shell Analysis */}
          <div className="bg-white rounded-lg border border-slate-200 p-3">
            <h4 className="text-xs font-semibold text-slate-600 uppercase tracking-wide mb-2 flex items-center gap-1">
              <Layers className="w-3 h-3" />
              Shell Analysis
            </h4>
            <div className="grid grid-cols-2 gap-3">
              <div className="p-2 bg-emerald-50 rounded border border-emerald-100">
                <div className="text-lg font-bold text-emerald-700">{shell.stats.primaryCount}</div>
                <div className="text-xs text-emerald-600">Primary (&lt;4A)</div>
                {shell.stats.avgPrimaryDistance > 0 && (
                  <div className="text-[10px] text-emerald-500 mt-1">
                    Avg: {shell.stats.avgPrimaryDistance.toFixed(2)}A
                  </div>
                )}
              </div>
              <div className="p-2 bg-blue-50 rounded border border-blue-100">
                <div className="text-lg font-bold text-blue-700">{shell.stats.secondaryCount}</div>
                <div className="text-xs text-blue-600">Secondary (4-6A)</div>
                {shell.stats.avgSecondaryDistance > 0 && (
                  <div className="text-[10px] text-blue-500 mt-1">
                    Avg: {shell.stats.avgSecondaryDistance.toFixed(2)}A
                  </div>
                )}
              </div>
            </div>
          </div>

          {/* Pharmacophore Features */}
          <div className="bg-white rounded-lg border border-slate-200 p-3">
            <div className="flex items-center justify-between mb-2">
              <h4 className="text-xs font-semibold text-slate-600 uppercase tracking-wide">
                Pharmacophore Features
              </h4>
              {onToggle3D && (
                <button
                  onClick={(e) => {
                    e.stopPropagation();
                    onToggle3D(!show3D);
                  }}
                  className={`flex items-center gap-1 px-2 py-1 rounded text-xs font-medium transition-colors ${
                    show3D
                      ? 'bg-purple-100 text-purple-700 border border-purple-200'
                      : 'bg-slate-100 text-slate-600 hover:bg-slate-200 border border-slate-200'
                  }`}
                >
                  {show3D ? <EyeOff className="w-3 h-3" /> : <Eye className="w-3 h-3" />}
                  {show3D ? 'Hide 3D' : 'Show 3D'}
                </button>
              )}
            </div>
            <div className="flex flex-wrap gap-2">
              {(Object.entries(pharmacophoreCounts) as [PharmacophoreType, number][])
                .filter(([, count]) => count > 0)
                .map(([type, count]) => (
                  <span
                    key={type}
                    className={`px-2 py-1 rounded text-xs font-medium ${getPharmacophoreColor(type)}`}
                  >
                    {type}: {count}
                  </span>
                ))}
              {Object.values(pharmacophoreCounts).every((c) => c === 0) && (
                <span className="text-xs text-slate-500">No features detected</span>
              )}
            </div>
          </div>

          {/* Suggestions */}
          {quality.suggestions.length > 0 && (
            <div className="bg-amber-50 rounded-lg border border-amber-200 p-3">
              <h4 className="text-xs font-semibold text-amber-700 uppercase tracking-wide mb-2 flex items-center gap-1">
                <Lightbulb className="w-3 h-3" />
                Suggestions
              </h4>
              <ul className="space-y-1">
                {quality.suggestions.map((suggestion, i) => (
                  <li key={i} className="text-xs text-amber-800 flex items-start gap-1.5">
                    <ChevronRight className="w-3 h-3 mt-0.5 flex-shrink-0" />
                    {suggestion}
                  </li>
                ))}
              </ul>
            </div>
          )}
        </div>
      )}
    </div>
  );
}

export function LigandBindingPanel({ onToggle3D, show3D }: LigandBindingPanelProps) {
  const { ligandData, focusedLigandIndex, setFocusedLigandIndex } = useStore();
  const [exportMenuOpen, setExportMenuOpen] = useState(false);

  const hasLigands = ligandData && ligandData.ligandCount > 0;

  // Handle exports
  const handleExportJSON = () => {
    if (ligandData) {
      downloadJSON(ligandData);
      setExportMenuOpen(false);
    }
  };

  const handleExportCSV = () => {
    if (ligandData) {
      downloadCSV(ligandData);
      setExportMenuOpen(false);
    }
  };

  const handleExportRFdiffusion = () => {
    if (ligandData) {
      downloadRFdiffusionHotspots(ligandData);
      setExportMenuOpen(false);
    }
  };

  if (!hasLigands) {
    return (
      <div className="bg-white rounded-lg border border-slate-200 p-4">
        <div className="flex items-center gap-2 text-slate-500">
          <HelpCircle className="w-4 h-4" />
          <span className="text-sm font-medium">No ligands detected</span>
        </div>
        <p className="text-xs text-slate-400 mt-2">
          Load a structure with ligands to see binding analysis.
        </p>
      </div>
    );
  }

  return (
    <div className="bg-white rounded-lg border border-slate-200 overflow-hidden">
      {/* Header */}
      <div className="px-4 py-3 border-b border-slate-200 bg-slate-50">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-2">
            <h3 className="font-semibold text-slate-800">Ligand Binding Analysis</h3>
            <span className="px-2 py-0.5 rounded-full text-xs font-medium bg-blue-100 text-blue-700">
              {ligandData.ligandCount}
            </span>
          </div>

          {/* Export dropdown */}
          <div className="relative">
            <button
              onClick={() => setExportMenuOpen(!exportMenuOpen)}
              className="flex items-center gap-1 px-2 py-1.5 rounded text-xs font-medium bg-slate-100 text-slate-600 hover:bg-slate-200 transition-colors border border-slate-200"
            >
              <Download className="w-3.5 h-3.5" />
              Export
              <ChevronDown className={`w-3 h-3 transition-transform ${exportMenuOpen ? 'rotate-180' : ''}`} />
            </button>

            {exportMenuOpen && (
              <>
                {/* Backdrop to close menu */}
                <div
                  className="fixed inset-0 z-10"
                  onClick={() => setExportMenuOpen(false)}
                />
                {/* Menu */}
                <div className="absolute right-0 mt-1 w-48 bg-white rounded-lg border border-slate-200 shadow-lg z-20 py-1">
                  <button
                    onClick={handleExportJSON}
                    className="w-full px-3 py-2 text-left text-sm text-slate-700 hover:bg-slate-50 flex items-center gap-2"
                  >
                    <FileJson className="w-4 h-4 text-blue-500" />
                    Export as JSON
                  </button>
                  <button
                    onClick={handleExportCSV}
                    className="w-full px-3 py-2 text-left text-sm text-slate-700 hover:bg-slate-50 flex items-center gap-2"
                  >
                    <FileSpreadsheet className="w-4 h-4 text-green-500" />
                    Export as CSV
                  </button>
                  <div className="border-t border-slate-100 my-1" />
                  <button
                    onClick={handleExportRFdiffusion}
                    className="w-full px-3 py-2 text-left text-sm text-slate-700 hover:bg-slate-50 flex items-center gap-2"
                  >
                    <Target className="w-4 h-4 text-purple-500" />
                    RFdiffusion Hotspots
                  </button>
                </div>
              </>
            )}
          </div>
        </div>
      </div>

      {/* Ligand list */}
      <div className="p-4 space-y-3 max-h-[500px] overflow-y-auto">
        {ligandData.ligandDetails.map((ligand, index) => (
          <LigandCard
            key={`${ligand.chainId}-${ligand.resSeq}-${index}`}
            ligand={ligand}
            index={index}
            isFocused={focusedLigandIndex === index}
            onFocus={setFocusedLigandIndex}
            onToggle3D={onToggle3D}
            show3D={show3D}
          />
        ))}
      </div>
    </div>
  );
}

export default LigandBindingPanel;
