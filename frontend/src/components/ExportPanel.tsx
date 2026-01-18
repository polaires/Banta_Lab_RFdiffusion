'use client';

import { useState } from 'react';
import { useStore } from '@/lib/store';
import api from '@/lib/api';
import { Download, FileText, Database, Loader2 } from 'lucide-react';

interface ExportPanelProps {
  pdbContent: string | null;
  cifContent?: string | null;
  source: 'rfd3' | 'rf3' | 'mpnn';
  filename?: string;
}

export function ExportPanel({ pdbContent, cifContent, source, filename = 'structure' }: ExportPanelProps) {
  const { latestConfidences } = useStore();
  const [exporting, setExporting] = useState(false);
  const [error, setError] = useState<string | null>(null);

  if (!pdbContent) return null;

  const downloadBlob = (content: string, name: string, type: string) => {
    const blob = new Blob([content], { type });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = name;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  const handleExportPdb = () => {
    downloadBlob(pdbContent, `${filename}.pdb`, 'text/plain');
  };

  const handleExportCif = async () => {
    // If we already have CIF content, use it
    if (cifContent) {
      downloadBlob(cifContent, `${filename}.cif`, 'text/plain');
      return;
    }

    // Otherwise, convert via API
    setExporting(true);
    setError(null);

    try {
      const result = await api.exportStructure({
        pdb_content: pdbContent,
        format: 'cif',
        include_confidences: source === 'rf3' && !!latestConfidences,
        ...(source === 'rf3' && latestConfidences && { confidences: latestConfidences }),
      });

      downloadBlob(result.content, result.filename, 'text/plain');
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Export failed');
    } finally {
      setExporting(false);
    }
  };

  const handleExportConfidences = () => {
    if (!latestConfidences) return;

    const jsonContent = JSON.stringify(latestConfidences, null, 2);
    downloadBlob(jsonContent, `${filename}_confidences.json`, 'application/json');
  };

  return (
    <div className="p-4 bg-muted/50 rounded-xl border border-border space-y-3">
      <h4 className="text-sm font-semibold text-foreground flex items-center gap-2">
        <Download className="w-4 h-4 text-muted-foreground" />
        Export Structure
      </h4>

      <div className="flex flex-wrap gap-2">
        {/* PDB Export */}
        <button
          onClick={handleExportPdb}
          className="px-3 py-1.5 bg-card border border-border hover:bg-muted rounded-lg text-sm text-foreground flex items-center gap-1.5 transition-colors"
        >
          <FileText className="w-3.5 h-3.5 text-muted-foreground" />
          PDB
        </button>

        {/* CIF Export */}
        <button
          onClick={handleExportCif}
          disabled={exporting}
          className="px-3 py-1.5 bg-card border border-border hover:bg-muted disabled:bg-muted disabled:text-muted-foreground rounded-lg text-sm text-foreground flex items-center gap-1.5 transition-colors"
        >
          {exporting ? (
            <Loader2 className="w-3.5 h-3.5 animate-spin" />
          ) : (
            <Database className="w-3.5 h-3.5 text-muted-foreground" />
          )}
          CIF
        </button>

        {/* Confidences JSON (only for RF3) */}
        {source === 'rf3' && latestConfidences && (
          <button
            onClick={handleExportConfidences}
            className="px-3 py-1.5 bg-violet-100 border border-violet-200 hover:bg-violet-200 rounded-lg text-sm text-violet-700 flex items-center gap-1.5 transition-colors"
          >
            <FileText className="w-3.5 h-3.5 text-violet-500" />
            Confidences (JSON)
          </button>
        )}
      </div>

      {error && (
        <p className="text-xs text-red-600">{error}</p>
      )}
    </div>
  );
}

export default ExportPanel;
