'use client';

import { useState, useRef } from 'react';
import { useStore } from '@/lib/store';
import api from '@/lib/api';
import { saveJob as saveJobToSupabase, updateJob as updateJobInSupabase } from '@/lib/supabase';
import { Sparkles, FlaskConical, Dna, Info, SlidersHorizontal, ChevronDown, AlertCircle, Loader2, Play, CheckCircle, Download, Upload } from 'lucide-react';

type ModelType = 'ligand_mpnn' | 'protein_mpnn';

export function MPNNPanel() {
  const {
    health,
    addJob,
    updateJob,
    addNotification,
    latestDesignPdb,
    latestRfd3Design,
    setLastCompletedJobType
  } = useStore();
  const [pdbContent, setPdbContent] = useState('');
  const [numSequences, setNumSequences] = useState(8);
  const [temperature, setTemperature] = useState(0.1);
  const [modelType, setModelType] = useState<ModelType>('ligand_mpnn');
  const [removeWaters, setRemoveWaters] = useState(true);
  const [submitting, setSubmitting] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<string | null>(null);
  const [resultModelType, setResultModelType] = useState<string | null>(null);
  const [isDragOver, setIsDragOver] = useState(false);
  const fileInputRef = useRef<HTMLInputElement>(null);

  // Advanced LigandMPNN parameters (Nature Methods 2025)
  const [showAdvanced, setShowAdvanced] = useState(false);
  const [packSideChains, setPackSideChains] = useState(false);
  const [packWithLigandContext, setPackWithLigandContext] = useState(true);
  const [numberOfPacksPerDesign, setNumberOfPacksPerDesign] = useState(4);
  const [biasAA, setBiasAA] = useState('');
  const [omitAA, setOmitAA] = useState('');
  const [modelNoiseLevel, setModelNoiseLevel] = useState<'005' | '010' | '020' | '030'>('010');
  const [ligandCutoffForScore, setLigandCutoffForScore] = useState(8.0);
  const [useSideChainContext, setUseSideChainContext] = useState(false);
  const [saveStats, setSaveStats] = useState(false);

  const handleFileUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (file) {
      const reader = new FileReader();
      reader.onload = (e) => {
        setPdbContent(e.target?.result as string);
      };
      reader.readAsText(file);
    }
  };

  const handleSubmit = async () => {
    if (!health) {
      setError('Backend not connected');
      return;
    }

    if (!pdbContent || pdbContent.length < 100) {
      setError('Please upload or paste a valid PDB file');
      return;
    }

    setError(null);
    setResult(null);
    setSubmitting(true);

    try {
      const request = {
        pdb_content: pdbContent,
        num_sequences: numSequences,
        temperature,
        model_type: modelType,
        remove_waters: removeWaters,
        // Advanced LigandMPNN parameters (only include if set)
        ...(packSideChains && { pack_side_chains: packSideChains }),
        ...(packSideChains && { pack_with_ligand_context: packWithLigandContext }),
        ...(packSideChains && { number_of_packs_per_design: numberOfPacksPerDesign }),
        ...(biasAA && { bias_AA: biasAA }),
        ...(omitAA && { omit_AA: omitAA }),
        ...(modelNoiseLevel !== '010' && { model_noise_level: modelNoiseLevel }),
        ...(ligandCutoffForScore !== 8.0 && { ligand_cutoff_for_score: ligandCutoffForScore }),
        ...(useSideChainContext && { use_side_chain_context: useSideChainContext }),
        ...(saveStats && { save_stats: saveStats }),
      };
      const response = await api.submitMPNNDesign(request);

      // Add to local store (UI state)
      addJob({
        id: response.job_id,
        type: 'mpnn',
        status: 'pending',
        createdAt: new Date().toISOString(),
      });

      // Persist to Supabase (async, non-blocking) - exclude pdb_content to save space
      saveJobToSupabase({
        runpod_id: response.job_id,
        type: 'mpnn',
        request: { num_sequences: numSequences, temperature, model_type: modelType },
      });

      const jobResult = await api.waitForJob(response.job_id, (status) => {
        // Update local store
        updateJob(response.job_id, {
          status: status.status,
          completedAt: status.completed_at,
          result: status.result,
          error: status.error,
        });
        // Update Supabase (async, non-blocking)
        updateJobInSupabase(response.job_id, {
          status: status.status,
          completed_at: status.completed_at || null,
          result: status.result || null,
        });
      });

      if (jobResult.status === 'completed' && jobResult.result?.sequences?.[0]) {
        setResult(jobResult.result.sequences[0].content);
        setResultModelType(jobResult.result.model_type || modelType);
        setLastCompletedJobType('mpnn');

        // Notify user with next step suggestion
        addNotification({
          type: 'success',
          title: 'Sequences designed!',
          message: `Generated ${numSequences} sequences using ${modelType}. Predict structure with RF3 to validate.`,
          action: {
            label: 'Predict Structure',
            tab: 'rf3',
          },
        });
      } else if (jobResult.status === 'failed') {
        addNotification({
          type: 'error',
          title: 'Sequence design failed',
          message: jobResult.error || 'Unknown error occurred',
        });
      }
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to submit job');
      addNotification({
        type: 'error',
        title: 'Submission failed',
        message: err instanceof Error ? err.message : 'Failed to submit job',
      });
    } finally {
      setSubmitting(false);
    }
  };

  const handleUseLatestDesign = () => {
    if (latestRfd3Design?.pdbContent) {
      setPdbContent(latestRfd3Design.pdbContent);
      addNotification({
        type: 'info',
        title: 'RFD3 design loaded',
        message: 'Latest RFD3 backbone has been loaded. Adjust parameters and submit.',
      });
    } else if (latestDesignPdb) {
      setPdbContent(latestDesignPdb);
      addNotification({
        type: 'info',
        title: 'Structure loaded',
        message: 'Latest design has been loaded. Adjust parameters and submit.',
      });
    }
  };

  const downloadFasta = () => {
    if (!result) return;
    const blob = new Blob([result], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'designed_sequences.fasta';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  const hasDesignAvailable = latestRfd3Design?.pdbContent || latestDesignPdb;

  const handleDrop = (e: React.DragEvent) => {
    e.preventDefault();
    setIsDragOver(false);
    const file = e.dataTransfer.files[0];
    if (file && (file.name.endsWith('.pdb') || file.name.endsWith('.cif'))) {
      const reader = new FileReader();
      reader.onload = (e) => setPdbContent(e.target?.result as string);
      reader.readAsText(file);
    }
  };

  return (
    <div className="space-y-8">
      {/* Header */}
      <div>
        <div className="flex items-center gap-3 mb-3">
          <span className="inline-flex items-center justify-center w-7 h-7 rounded-full bg-violet-100 text-violet-700 text-xs font-bold">2</span>
          <h2 className="text-lg font-bold text-slate-900">ProteinMPNN — Sequence Design</h2>
        </div>
        <p className="text-slate-600 text-sm leading-relaxed">
          Design amino acid sequences for your RFD3 backbone. After design, validate with RF3 structure prediction.
        </p>
      </div>

      {/* Use Latest Design button */}
      {hasDesignAvailable && (
        <button
          onClick={handleUseLatestDesign}
          className="w-full py-3 bg-gradient-to-r from-blue-600 to-violet-600 hover:from-blue-700 hover:to-violet-700 text-white rounded-xl font-semibold flex items-center justify-center gap-2 transition-all shadow-sm hover:shadow-md"
        >
          <Sparkles className="w-5 h-5" />
          Use Latest Design from RFD3/RF3
        </button>
      )}

      <div className="grid grid-cols-1 lg:grid-cols-12 gap-8">
        {/* Left Column - Parameters */}
        <div className="lg:col-span-7 space-y-6">
          {/* Model Type Selection */}
          <div className="space-y-3">
            <label className="text-xs font-bold text-slate-700 uppercase tracking-wider">Model Type</label>
            <div className="grid grid-cols-2 gap-3">
              <button
                onClick={() => setModelType('ligand_mpnn')}
                className={`p-4 rounded-xl border-2 transition-all text-left ${
                  modelType === 'ligand_mpnn'
                    ? 'border-violet-500 bg-violet-50 shadow-sm'
                    : 'border-slate-200 hover:border-slate-300 bg-white'
                }`}
              >
                <div className="flex items-center gap-2 mb-1">
                  <FlaskConical className="w-4 h-4 text-violet-600" />
                  <span className="font-semibold text-slate-900">LigandMPNN</span>
                </div>
                <p className="text-xs text-slate-500">Ligand-aware design (recommended)</p>
              </button>
              <button
                onClick={() => setModelType('protein_mpnn')}
                className={`p-4 rounded-xl border-2 transition-all text-left ${
                  modelType === 'protein_mpnn'
                    ? 'border-violet-500 bg-violet-50 shadow-sm'
                    : 'border-slate-200 hover:border-slate-300 bg-white'
                }`}
              >
                <div className="flex items-center gap-2 mb-1">
                  <Dna className="w-4 h-4 text-slate-600" />
                  <span className="font-semibold text-slate-900">ProteinMPNN</span>
                </div>
                <p className="text-xs text-slate-500">Original protein-only model</p>
              </button>
            </div>
            <div className="flex items-start gap-2 text-xs text-slate-500 bg-slate-50 p-3 rounded-lg">
              <Info className="w-4 h-4 mt-0.5 flex-shrink-0" />
              <span>LigandMPNN is recommended for most use cases. Use ProteinMPNN for pure protein design without ligands.</span>
            </div>
          </div>

          {/* Parameters */}
          <div className="grid grid-cols-2 gap-4">
            <div className="space-y-2">
              <label className="text-xs font-bold text-slate-700 uppercase tracking-wider">Number of Sequences</label>
              <input
                type="number"
                value={numSequences}
                onChange={(e) => setNumSequences(Math.max(1, Math.min(100, parseInt(e.target.value) || 8)))}
                min={1}
                max={100}
                className="w-full px-4 py-3 bg-slate-50 rounded-xl border border-slate-200 focus:border-violet-500 focus:ring-2 focus:ring-violet-500/20 focus:outline-none text-slate-900 transition-all"
              />
            </div>
            <div className="space-y-2">
              <label className="text-xs font-bold text-slate-700 uppercase tracking-wider">Temperature</label>
              <input
                type="number"
                value={temperature}
                onChange={(e) => setTemperature(Math.max(0.01, Math.min(2, parseFloat(e.target.value) || 0.1)))}
                min={0.01}
                max={2}
                step={0.01}
                className="w-full px-4 py-3 bg-slate-50 rounded-xl border border-slate-200 focus:border-violet-500 focus:ring-2 focus:ring-violet-500/20 focus:outline-none text-slate-900 transition-all"
              />
              <p className="text-xs text-slate-500">Lower = conservative, Higher = diverse</p>
            </div>
          </div>

          {/* Remove Waters Toggle */}
          <label className="flex items-center gap-3 p-4 bg-slate-50 rounded-xl cursor-pointer hover:bg-slate-100 transition-colors">
            <input
              type="checkbox"
              checked={removeWaters}
              onChange={(e) => setRemoveWaters(e.target.checked)}
              className="w-4 h-4 rounded border-slate-300 text-violet-600 focus:ring-violet-500"
            />
            <div>
              <span className="text-sm font-medium text-slate-900">Remove water molecules</span>
              <p className="text-xs text-slate-500">Exclude water from structural context</p>
            </div>
          </label>

          {/* Advanced Options Toggle */}
          <button
            onClick={() => setShowAdvanced(!showAdvanced)}
            className="w-full flex items-center justify-between p-4 bg-slate-50 rounded-xl hover:bg-slate-100 transition-colors"
          >
            <div className="flex items-center gap-2">
              <SlidersHorizontal className="w-4 h-4 text-slate-600" />
              <span className="text-sm font-medium text-slate-900">Advanced LigandMPNN Options</span>
            </div>
            <ChevronDown className={`w-4 h-4 text-slate-400 transition-transform ${showAdvanced ? 'rotate-180' : ''}`} />
          </button>

          {/* Advanced Options Panel */}
          {showAdvanced && (
            <div className="space-y-4 p-4 bg-gradient-to-b from-violet-50 to-slate-50 rounded-xl border border-violet-100">
              <div className="flex items-center gap-2 mb-2">
                <FlaskConical className="w-4 h-4 text-violet-600" />
                <span className="text-xs font-bold text-violet-700 uppercase tracking-wider">LigandMPNN Settings (Nature Methods 2025)</span>
              </div>

              {/* Model Noise Level */}
              <div className="space-y-2">
                <label className="text-xs font-medium text-slate-700">Model Noise Level</label>
                <select
                  value={modelNoiseLevel}
                  onChange={(e) => setModelNoiseLevel(e.target.value as '005' | '010' | '020' | '030')}
                  className="w-full px-3 py-2 bg-white rounded-lg border border-slate-200 focus:border-violet-500 focus:ring-2 focus:ring-violet-500/20 focus:outline-none text-sm"
                >
                  <option value="005">0.05Å - High accuracy (may overfit)</option>
                  <option value="010">0.10Å - Default (balanced)</option>
                  <option value="020">0.20Å - Robust to backbone errors</option>
                  <option value="030">0.30Å - Very robust (lower recovery)</option>
                </select>
              </div>

              {/* Ligand Cutoff */}
              <div className="space-y-2">
                <label className="text-xs font-medium text-slate-700">Ligand Cutoff Distance: {ligandCutoffForScore}Å</label>
                <input
                  type="range"
                  min={4}
                  max={12}
                  step={0.5}
                  value={ligandCutoffForScore}
                  onChange={(e) => setLigandCutoffForScore(parseFloat(e.target.value))}
                  className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-violet-600"
                />
                <p className="text-xs text-slate-500">Distance for ligand-adjacent residue scoring (6Å recommended for small molecules)</p>
              </div>

              {/* AA Bias */}
              <div className="space-y-2">
                <label className="text-xs font-medium text-slate-700">Amino Acid Bias</label>
                <input
                  type="text"
                  value={biasAA}
                  onChange={(e) => setBiasAA(e.target.value)}
                  placeholder="W:3.0,Y:2.0,F:1.5,C:-5.0"
                  className="w-full px-3 py-2 bg-white rounded-lg border border-slate-200 focus:border-violet-500 focus:ring-2 focus:ring-violet-500/20 focus:outline-none text-sm font-mono"
                />
                <p className="text-xs text-slate-500">Positive = favor, Negative = penalize (e.g., favor aromatics for hydrophobic ligands)</p>
              </div>

              {/* Omit AA */}
              <div className="space-y-2">
                <label className="text-xs font-medium text-slate-700">Omit Amino Acids</label>
                <input
                  type="text"
                  value={omitAA}
                  onChange={(e) => setOmitAA(e.target.value)}
                  placeholder="C"
                  className="w-full px-3 py-2 bg-white rounded-lg border border-slate-200 focus:border-violet-500 focus:ring-2 focus:ring-violet-500/20 focus:outline-none text-sm font-mono"
                />
                <p className="text-xs text-slate-500">Completely exclude these residues (e.g., C to avoid disulfides)</p>
              </div>

              {/* Sidechain Packing */}
              <div className="space-y-3 pt-2 border-t border-violet-200">
                <label className="flex items-center gap-3 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={packSideChains}
                    onChange={(e) => setPackSideChains(e.target.checked)}
                    className="w-4 h-4 rounded border-slate-300 text-violet-600 focus:ring-violet-500"
                  />
                  <div>
                    <span className="text-sm font-medium text-slate-900">Enable Sidechain Packing</span>
                    <p className="text-xs text-slate-500">Generate sidechain conformations (chi angles)</p>
                  </div>
                </label>

                {packSideChains && (
                  <div className="ml-7 space-y-3">
                    <label className="flex items-center gap-3 cursor-pointer">
                      <input
                        type="checkbox"
                        checked={packWithLigandContext}
                        onChange={(e) => setPackWithLigandContext(e.target.checked)}
                        className="w-4 h-4 rounded border-slate-300 text-violet-600 focus:ring-violet-500"
                      />
                      <span className="text-sm text-slate-700">Pack with ligand context</span>
                    </label>

                    <div className="space-y-1">
                      <label className="text-xs font-medium text-slate-700">Packs per design: {numberOfPacksPerDesign}</label>
                      <input
                        type="range"
                        min={1}
                        max={10}
                        value={numberOfPacksPerDesign}
                        onChange={(e) => setNumberOfPacksPerDesign(parseInt(e.target.value))}
                        className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-violet-600"
                      />
                    </div>
                  </div>
                )}
              </div>

              {/* Additional Options */}
              <div className="space-y-2 pt-2 border-t border-violet-200">
                <label className="flex items-center gap-3 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={useSideChainContext}
                    onChange={(e) => setUseSideChainContext(e.target.checked)}
                    className="w-4 h-4 rounded border-slate-300 text-violet-600 focus:ring-violet-500"
                  />
                  <div>
                    <span className="text-sm font-medium text-slate-900">Use sidechain context</span>
                    <p className="text-xs text-slate-500">Include fixed residue sidechains as design context</p>
                  </div>
                </label>

                <label className="flex items-center gap-3 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={saveStats}
                    onChange={(e) => setSaveStats(e.target.checked)}
                    className="w-4 h-4 rounded border-slate-300 text-violet-600 focus:ring-violet-500"
                  />
                  <div>
                    <span className="text-sm font-medium text-slate-900">Save confidence metrics</span>
                    <p className="text-xs text-slate-500">Return ligand_confidence and overall_confidence scores</p>
                  </div>
                </label>
              </div>
            </div>
          )}
        </div>

        {/* Right Column - File Upload */}
        <div className="lg:col-span-5 flex flex-col">
          <label className="text-xs font-bold text-slate-700 uppercase tracking-wider mb-3">Input Structure</label>
          <div
            onDragOver={(e) => { e.preventDefault(); setIsDragOver(true); }}
            onDragLeave={() => setIsDragOver(false)}
            onDrop={handleDrop}
            onClick={() => fileInputRef.current?.click()}
            className={`flex-1 min-h-[200px] border-2 border-dashed rounded-2xl cursor-pointer transition-all flex flex-col items-center justify-center gap-3 ${
              isDragOver
                ? 'border-violet-500 bg-violet-50'
                : pdbContent
                ? 'border-emerald-300 bg-emerald-50'
                : 'border-slate-200 bg-slate-50 hover:border-violet-400 hover:bg-violet-50/50'
            }`}
          >
            <input
              ref={fileInputRef}
              type="file"
              accept=".pdb,.cif"
              onChange={handleFileUpload}
              className="hidden"
            />
            <div className={`w-14 h-14 rounded-2xl flex items-center justify-center ${
              pdbContent ? 'bg-emerald-100' : 'bg-white shadow-sm'
            }`}>
              {pdbContent ? (
                <CheckCircle className="w-8 h-8 text-emerald-600" />
              ) : (
                <Upload className="w-8 h-8 text-slate-400" />
              )}
            </div>
            <div className="text-center">
              <p className={`font-semibold ${pdbContent ? 'text-emerald-700' : 'text-slate-700'}`}>
                {pdbContent ? 'Structure loaded' : 'Drop PDB/CIF here'}
              </p>
              <p className="text-xs text-slate-500 mt-1">or click to browse</p>
            </div>
          </div>

          {/* Paste PDB */}
          <div className="mt-4 space-y-2">
            <label className="text-xs font-medium text-slate-600">Or paste PDB content</label>
            <textarea
              value={pdbContent}
              onChange={(e) => setPdbContent(e.target.value)}
              placeholder="ATOM      1  N   ALA A   1      ..."
              rows={4}
              className="w-full px-4 py-3 bg-slate-50 rounded-xl border border-slate-200 focus:border-violet-500 focus:ring-2 focus:ring-violet-500/20 focus:outline-none font-mono text-xs text-slate-900 transition-all"
            />
          </div>
        </div>
      </div>

      {/* Error Display */}
      {error && (
        <div className="p-4 bg-red-50 border border-red-200 rounded-xl text-sm text-red-700 flex items-start gap-3">
          <AlertCircle className="w-5 h-5 text-red-500 flex-shrink-0" />
          {error}
        </div>
      )}

      {/* Submit Button */}
      <button
        onClick={handleSubmit}
        disabled={!health || submitting || !pdbContent}
        className="w-full py-4 bg-violet-600 hover:bg-violet-700 disabled:bg-slate-300 disabled:cursor-not-allowed text-white rounded-xl font-semibold flex items-center justify-center gap-2 transition-all shadow-sm hover:shadow-md"
      >
        {submitting ? (
          <>
            <Loader2 className="w-5 h-5 animate-spin" />
            Designing Sequences...
          </>
        ) : (
          <>
            <Play className="w-5 h-5" />
            Design Sequences
          </>
        )}
      </button>

      {/* Results */}
      {result && (
        <div className="space-y-3 pt-4 border-t border-slate-200">
          <div className="flex items-center justify-between">
            <div className="flex items-center gap-2">
              <CheckCircle className="w-5 h-5 text-emerald-600" />
              <label className="text-sm font-semibold text-slate-900">
                Designed Sequences
                {resultModelType && (
                  <span className="ml-2 text-xs font-normal text-slate-500">via {resultModelType}</span>
                )}
              </label>
            </div>
            <button
              onClick={downloadFasta}
              className="px-3 py-1.5 text-xs bg-slate-100 hover:bg-slate-200 text-slate-700 rounded-lg flex items-center gap-1.5 transition font-medium"
            >
              <Download className="w-4 h-4" />
              Download FASTA
            </button>
          </div>
          <pre className="p-4 bg-slate-900 rounded-xl text-xs font-mono overflow-auto max-h-64 text-emerald-400">
            {result}
          </pre>
        </div>
      )}
    </div>
  );
}
