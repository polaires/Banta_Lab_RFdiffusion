'use client';

import { useRef, useState } from 'react';
import { Check, Upload, Download, Loader2 } from 'lucide-react';
import { useStore } from '@/lib/store';
import { fetchCatalyticSuggestions, extractPdbId, filterPdbToSingleAssembly } from '@/lib/catalyticDetection';

interface PdbUploaderProps {
  label: string;
  description?: string;
  required?: boolean;
  value: string | null;
  fileName: string | null;
  onChange: (content: string | null, fileName: string | null) => void;
  accept?: string;
  className?: string;
}

export function PdbUploader({
  label,
  description,
  required = false,
  value,
  fileName,
  onChange,
  accept = '.pdb,.cif',
  className = '',
}: PdbUploaderProps) {
  const [isDragOver, setIsDragOver] = useState(false);
  const [pdbIdInput, setPdbIdInput] = useState('');
  const [isLoadingFromRcsb, setIsLoadingFromRcsb] = useState(false);
  const [loadError, setLoadError] = useState<string | null>(null);
  const fileInputRef = useRef<HTMLInputElement>(null);

  const {
    backendUrl,
    setCatalyticSuggestions,
    setSuggestionsLoading,
    setSuggestionsError,
    clearSuggestions,
    setSelectedPdb,
  } = useStore();

  // Fetch catalytic suggestions when structure is loaded
  const fetchSuggestions = async (content: string) => {
    setSuggestionsLoading(true);
    try {
      const pdbId = extractPdbId(content);
      const result = await fetchCatalyticSuggestions(content, pdbId, backendUrl);
      setCatalyticSuggestions(result.residues, result.source);
    } catch (error) {
      console.error('Failed to fetch catalytic suggestions:', error);
      setSuggestionsError(error instanceof Error ? error.message : 'Unknown error');
    } finally {
      setSuggestionsLoading(false);
    }
  };

  const handleFileUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (file) {
      const reader = new FileReader();
      reader.onload = (e) => {
        const rawContent = e.target?.result as string;
        // Filter to single biological assembly to avoid duplicate ligands
        const content = filterPdbToSingleAssembly(rawContent);
        onChange(content, file.name);
        setSelectedPdb(content);
        // Auto-fetch catalytic suggestions
        fetchSuggestions(content);
      };
      reader.readAsText(file);
    }
  };

  const handleDrop = (e: React.DragEvent) => {
    e.preventDefault();
    setIsDragOver(false);
    const file = e.dataTransfer.files[0];
    if (file) {
      const reader = new FileReader();
      reader.onload = (ev) => {
        const rawContent = ev.target?.result as string;
        // Filter to single biological assembly to avoid duplicate ligands
        const content = filterPdbToSingleAssembly(rawContent);
        onChange(content, file.name);
        setSelectedPdb(content);
        // Auto-fetch catalytic suggestions
        fetchSuggestions(content);
      };
      reader.readAsText(file);
    }
  };

  const handleClear = (e: React.MouseEvent) => {
    e.stopPropagation();
    onChange(null, null);
    clearSuggestions();
    setSelectedPdb(null);
    if (fileInputRef.current) {
      fileInputRef.current.value = '';
    }
  };

  // Fetch PDB from RCSB via our proxy (to avoid CORS issues)
  const fetchFromRcsb = async (pdbId: string): Promise<string> => {
    const proxyUrl = `/api/rcsb/${pdbId}`;
    const response = await fetch(proxyUrl);

    if (!response.ok) {
      const errorData = await response.json().catch(() => ({}));
      throw new Error(errorData.error || `PDB ID "${pdbId}" not found`);
    }

    return response.text();
  };

  // Load PDB from RCSB by ID
  const handleLoadFromRcsb = async (e: React.MouseEvent) => {
    e.stopPropagation();
    const pdbId = pdbIdInput.trim().toUpperCase();

    if (!pdbId || pdbId.length !== 4) {
      setLoadError('Please enter a valid 4-character PDB ID');
      return;
    }

    setIsLoadingFromRcsb(true);
    setLoadError(null);

    try {
      const rawContent = await fetchFromRcsb(pdbId);
      // Filter to single biological assembly to avoid duplicate ligands
      const content = filterPdbToSingleAssembly(rawContent);
      const generatedFileName = `${pdbId}.pdb`;

      onChange(content, generatedFileName);
      setSelectedPdb(content);
      setPdbIdInput('');

      // Fetch catalytic suggestions with the known PDB ID
      setSuggestionsLoading(true);
      try {
        const result = await fetchCatalyticSuggestions(content, pdbId, backendUrl);
        setCatalyticSuggestions(result.residues, result.source);
      } catch (error) {
        console.error('Failed to fetch catalytic suggestions:', error);
        setSuggestionsError(error instanceof Error ? error.message : 'Unknown error');
      } finally {
        setSuggestionsLoading(false);
      }
    } catch (error) {
      console.error('Failed to load PDB from RCSB:', error);
      setLoadError(error instanceof Error ? error.message : 'Failed to load PDB');
    } finally {
      setIsLoadingFromRcsb(false);
    }
  };

  return (
    <div className={`space-y-2 ${className}`}>
      <label className="block text-sm font-medium text-foreground">
        {label}
        {required && <span className="text-red-500 ml-1">*</span>}
      </label>
      {description && (
        <p className="text-xs text-muted-foreground">{description}</p>
      )}
      <div
        onDragOver={(e) => { e.preventDefault(); setIsDragOver(true); }}
        onDragLeave={() => setIsDragOver(false)}
        onDrop={handleDrop}
        onClick={() => fileInputRef.current?.click()}
        className={`relative border-2 border-dashed rounded-lg transition-all cursor-pointer
                   flex flex-col items-center justify-center p-5 text-center
                   ${isDragOver
                     ? 'border-primary bg-primary/10'
                     : value
                       ? 'border-border bg-muted/50'
                       : 'border-border hover:border-primary hover:bg-muted'
                   }`}
      >
        <input
          ref={fileInputRef}
          type="file"
          accept={accept}
          onChange={handleFileUpload}
          className="hidden"
        />

        {value ? (
          <>
            <div className="w-10 h-10 bg-muted rounded-full flex items-center justify-center mb-2">
              <Check className="h-5 w-5 text-foreground" />
            </div>
            <div className="text-sm text-foreground font-medium">{fileName}</div>
            <button
              onClick={handleClear}
              className="mt-2 text-xs text-muted-foreground hover:text-red-600 transition-colors"
            >
              Remove file
            </button>
          </>
        ) : (
          <>
            <div className="w-10 h-10 bg-muted rounded-full flex items-center justify-center mb-2 group-hover:scale-110 transition-transform">
              <Upload className="h-5 w-5 text-muted-foreground" />
            </div>
            <div className="text-sm text-muted-foreground">
              Drop file here or <span className="text-primary font-medium">browse</span>
            </div>
            <p className="text-xs text-muted-foreground mt-1">
              {accept.replace(/\./g, '').toUpperCase()} files
            </p>
          </>
        )}
      </div>

      {/* PDB ID Input Section */}
      <div className="flex items-center gap-2">
        <div className="relative flex-1">
          <input
            type="text"
            value={pdbIdInput}
            onChange={(e) => {
              setPdbIdInput(e.target.value.toUpperCase());
              setLoadError(null);
            }}
            placeholder="Enter PDB ID (e.g., 1TIM)"
            maxLength={4}
            className="w-full px-3 py-2 text-sm border border-border rounded-md bg-background
                       focus:outline-none focus:ring-2 focus:ring-primary/50 focus:border-primary
                       placeholder:text-muted-foreground"
            disabled={isLoadingFromRcsb}
          />
        </div>
        <button
          onClick={handleLoadFromRcsb}
          disabled={isLoadingFromRcsb || !pdbIdInput.trim()}
          className="flex items-center gap-2 px-3 py-2 text-sm font-medium rounded-md
                     bg-primary text-primary-foreground hover:bg-primary/90
                     disabled:opacity-50 disabled:cursor-not-allowed
                     transition-colors"
        >
          {isLoadingFromRcsb ? (
            <>
              <Loader2 className="h-4 w-4 animate-spin" />
              Loading...
            </>
          ) : (
            <>
              <Download className="h-4 w-4" />
              Load from RCSB
            </>
          )}
        </button>
      </div>

      {/* Error message */}
      {loadError && (
        <p className="text-xs text-red-500">{loadError}</p>
      )}
    </div>
  );
}
