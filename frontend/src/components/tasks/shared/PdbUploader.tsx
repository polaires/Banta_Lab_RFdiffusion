'use client';

import { useRef, useState } from 'react';
import { Check, Upload } from 'lucide-react';

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

  const handleFileUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (file) {
      const reader = new FileReader();
      reader.onload = (e) => {
        onChange(e.target?.result as string, file.name);
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
        onChange(ev.target?.result as string, file.name);
      };
      reader.readAsText(file);
    }
  };

  const handleClear = (e: React.MouseEvent) => {
    e.stopPropagation();
    onChange(null, null);
    if (fileInputRef.current) {
      fileInputRef.current.value = '';
    }
  };

  return (
    <div className={`space-y-2 ${className}`}>
      <label className="block text-sm font-medium text-slate-700">
        {label}
        {required && <span className="text-red-500 ml-1">*</span>}
      </label>
      {description && (
        <p className="text-xs text-slate-500">{description}</p>
      )}
      <div
        onDragOver={(e) => { e.preventDefault(); setIsDragOver(true); }}
        onDragLeave={() => setIsDragOver(false)}
        onDrop={handleDrop}
        onClick={() => fileInputRef.current?.click()}
        className={`relative border-2 border-dashed rounded-lg transition-all cursor-pointer
                   flex flex-col items-center justify-center p-5 text-center
                   ${isDragOver
                     ? 'border-blue-400 bg-blue-50/30'
                     : value
                       ? 'border-slate-300 bg-slate-50/50'
                       : 'border-slate-200 hover:border-blue-400 hover:bg-slate-50'
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
            <div className="w-10 h-10 bg-slate-200 rounded-full flex items-center justify-center mb-2">
              <Check className="h-5 w-5 text-slate-600" />
            </div>
            <div className="text-sm text-slate-700 font-medium">{fileName}</div>
            <button
              onClick={handleClear}
              className="mt-2 text-xs text-slate-500 hover:text-red-600 transition-colors"
            >
              Remove file
            </button>
          </>
        ) : (
          <>
            <div className="w-10 h-10 bg-slate-100 rounded-full flex items-center justify-center mb-2 group-hover:scale-110 transition-transform">
              <Upload className="h-5 w-5 text-slate-400" />
            </div>
            <div className="text-sm text-slate-600">
              Drop file here or <span className="text-blue-600 font-medium">browse</span>
            </div>
            <p className="text-xs text-slate-400 mt-1">
              {accept.replace(/\./g, '').toUpperCase()} files
            </p>
          </>
        )}
      </div>
    </div>
  );
}
