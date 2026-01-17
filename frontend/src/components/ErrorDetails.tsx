'use client';

import { useState } from 'react';
import { AlertCircle, ChevronUp, ChevronDown, Check, Copy, Cpu, HardDrive, AlertTriangle } from 'lucide-react';

interface ErrorContext {
  task?: string;
  input_keys?: string[];
  gpu_info?: {
    available?: boolean;
    name?: string;
    memory_gb?: number;
  };
  gpu_memory_used_mb?: number;
  gpu_memory_total_mb?: number;
  foundry_available?: boolean;
  checkpoint_dir?: string;
}

interface ErrorDetailsProps {
  error: string;
  errorType?: string;
  traceback?: string;
  context?: ErrorContext;
  className?: string;
  compact?: boolean;
}

export function ErrorDetails({
  error,
  errorType,
  traceback,
  context,
  className = '',
  compact = false,
}: ErrorDetailsProps) {
  const [expanded, setExpanded] = useState(false);
  const [copied, setCopied] = useState(false);

  const copyErrorDetails = () => {
    const details = [
      `Error: ${error}`,
      errorType && `Type: ${errorType}`,
      context?.task && `Task: ${context.task}`,
      context?.gpu_info?.name && `GPU: ${context.gpu_info.name}`,
      context?.gpu_memory_used_mb && context?.gpu_memory_total_mb &&
        `GPU Memory: ${context.gpu_memory_used_mb}/${context.gpu_memory_total_mb} MB`,
      traceback && `\nTraceback:\n${traceback}`,
    ].filter(Boolean).join('\n');

    navigator.clipboard.writeText(details);
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
  };

  if (compact) {
    return (
      <div className={`bg-red-50 border border-red-200 rounded-lg p-3 ${className}`}>
        <div className="flex items-start gap-2">
          <AlertCircle className="w-4 h-4 text-red-500 flex-shrink-0 mt-0.5" />
          <div className="flex-1 min-w-0">
            <p className="text-sm text-red-700 font-medium truncate" title={error}>
              {errorType ? `${errorType}: ` : ''}{error}
            </p>
            {traceback && (
              <button
                onClick={() => setExpanded(!expanded)}
                className="text-xs text-red-600 hover:text-red-800 mt-1 flex items-center gap-1"
              >
                {expanded ? <ChevronUp className="w-3 h-3" /> : <ChevronDown className="w-3 h-3" />}
                {expanded ? 'Hide details' : 'Show details'}
              </button>
            )}
          </div>
        </div>
        {expanded && traceback && (
          <pre className="mt-2 text-xs text-red-800 bg-red-100 rounded p-2 overflow-x-auto max-h-48 overflow-y-auto font-mono">
            {traceback}
          </pre>
        )}
      </div>
    );
  }

  return (
    <div className={`bg-red-50 border border-red-200 rounded-xl overflow-hidden ${className}`}>
      {/* Header */}
      <div className="p-4 flex items-start gap-3">
        <div className="w-10 h-10 bg-red-100 rounded-lg flex items-center justify-center flex-shrink-0">
          <AlertCircle className="w-5 h-5 text-red-600" />
        </div>
        <div className="flex-1 min-w-0">
          <div className="flex items-center gap-2 flex-wrap">
            <h3 className="text-sm font-semibold text-red-800">
              {errorType || 'Error'}
            </h3>
            {context?.task && (
              <span className="text-xs bg-red-100 text-red-700 px-2 py-0.5 rounded-full font-medium">
                {context.task.toUpperCase()}
              </span>
            )}
          </div>
          <p className="text-sm text-red-700 mt-1">{error}</p>
        </div>
        <button
          onClick={copyErrorDetails}
          className="p-2 hover:bg-red-100 rounded-lg transition-colors text-red-600 flex-shrink-0"
          title="Copy error details"
        >
          {copied ? <Check className="w-4 h-4" /> : <Copy className="w-4 h-4" />}
        </button>
      </div>

      {/* Context Info */}
      {context && (context.gpu_info || context.gpu_memory_used_mb) && (
        <div className="px-4 pb-3 flex flex-wrap gap-3 text-xs">
          {context.gpu_info?.name && (
            <div className="flex items-center gap-1.5 text-red-600">
              <Cpu className="w-3 h-3" />
              <span>{context.gpu_info.name}</span>
            </div>
          )}
          {context.gpu_memory_used_mb && context.gpu_memory_total_mb && (
            <div className="flex items-center gap-1.5 text-red-600">
              <HardDrive className="w-3 h-3" />
              <span>
                {Math.round(context.gpu_memory_used_mb / 1024 * 10) / 10}/
                {Math.round(context.gpu_memory_total_mb / 1024 * 10) / 10} GB VRAM
              </span>
              {context.gpu_memory_used_mb / context.gpu_memory_total_mb > 0.9 && (
                <span className="text-red-700 font-semibold">(OOM likely)</span>
              )}
            </div>
          )}
          {context.foundry_available === false && (
            <div className="flex items-center gap-1.5 text-red-700 font-medium">
              <AlertTriangle className="w-3 h-3" />
              <span>Foundry not available</span>
            </div>
          )}
        </div>
      )}

      {/* Expandable Traceback */}
      {traceback && (
        <div className="border-t border-red-200">
          <button
            onClick={() => setExpanded(!expanded)}
            className="w-full px-4 py-2 flex items-center justify-between text-xs text-red-700 hover:bg-red-100 transition-colors"
          >
            <span className="font-medium">Stack Trace</span>
            {expanded ? <ChevronUp className="w-4 h-4" /> : <ChevronDown className="w-4 h-4" />}
          </button>
          {expanded && (
            <div className="px-4 pb-4">
              <pre className="text-xs text-red-800 bg-red-100 rounded-lg p-3 overflow-x-auto max-h-96 overflow-y-auto font-mono whitespace-pre-wrap break-words">
                {traceback}
              </pre>
            </div>
          )}
        </div>
      )}
    </div>
  );
}
