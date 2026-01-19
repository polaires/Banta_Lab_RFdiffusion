'use client';

import { useRef, useState } from 'react';
import { Check, Upload } from 'lucide-react';
import { useStore } from '@/lib/store';
import { fetchCatalyticSuggestions, extractPdbId } from '@/lib/catalyticDetection';

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
        const content = e.target?.result as string;
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
        const content = ev.target?.result as string;
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
    if (fileInputRef.current) {
      fileInputRef.current.value = '';
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
    </div>
  );
}
