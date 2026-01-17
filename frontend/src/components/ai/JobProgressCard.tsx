'use client';

import { Clock, RefreshCw, CheckCircle, AlertCircle, Loader2, Info } from 'lucide-react';

interface JobProgressCardProps {
  jobId: string;
  status: 'pending' | 'running' | 'completed' | 'failed';
  progress?: number;
  message?: string;
}

export function JobProgressCard({ jobId, status, progress, message }: JobProgressCardProps) {
  const getStatusConfig = () => {
    switch (status) {
      case 'pending':
        return {
          color: 'text-blue-500',
          bgColor: 'bg-blue-50 border-blue-200',
          label: 'Queued',
        };
      case 'running':
        return {
          color: 'text-violet-500',
          bgColor: 'bg-violet-50 border-violet-200',
          label: 'Running',
        };
      case 'completed':
        return {
          color: 'text-emerald-500',
          bgColor: 'bg-emerald-50 border-emerald-200',
          label: 'Completed',
        };
      case 'failed':
        return {
          color: 'text-red-500',
          bgColor: 'bg-red-50 border-red-200',
          label: 'Failed',
        };
    }
  };

  const getStatusIcon = () => {
    switch (status) {
      case 'pending':
        return <Clock className="h-6 w-6" />;
      case 'running':
        return <Loader2 className="h-6 w-6 animate-spin" />;
      case 'completed':
        return <CheckCircle className="h-6 w-6" />;
      case 'failed':
        return <AlertCircle className="h-6 w-6" />;
    }
  };

  const config = getStatusConfig();

  return (
    <div className={`rounded-xl p-4 border ${config.bgColor}`}>
      <div className="flex items-center gap-3">
        {/* Status icon with animation for running */}
        <div className={`${config.color}`}>
          {getStatusIcon()}
        </div>

        {/* Job info */}
        <div className="flex-1">
          <div className="flex items-center gap-2">
            <span className="font-medium text-slate-800">RFD3 Design</span>
            <span className={`text-xs px-2 py-0.5 rounded-full ${
              status === 'running' ? 'bg-violet-100 text-violet-700' :
              status === 'completed' ? 'bg-emerald-100 text-emerald-700' :
              status === 'failed' ? 'bg-red-100 text-red-700' :
              'bg-blue-100 text-blue-700'
            }`}>
              {config.label}
            </span>
          </div>
          <div className="text-xs text-slate-500 font-mono mt-0.5">
            Job ID: {jobId.slice(0, 12)}...
          </div>
        </div>

        {/* Progress percentage */}
        {status === 'running' && progress !== undefined && (
          <div className="text-right">
            <div className="text-lg font-bold text-violet-600">{progress}%</div>
          </div>
        )}
      </div>

      {/* Progress bar */}
      {status === 'running' && (
        <div className="mt-3">
          <div className="h-2 bg-violet-100 rounded-full overflow-hidden">
            <div
              className="h-full bg-gradient-to-r from-violet-500 to-purple-600 rounded-full transition-all duration-500"
              style={{ width: `${progress || 0}%` }}
            />
          </div>
          {message && (
            <p className="text-xs text-violet-600 mt-2">{message}</p>
          )}
        </div>
      )}

      {/* Completed message */}
      {status === 'completed' && (
        <div className="mt-2 text-sm text-emerald-700 flex items-center gap-1">
          <Info className="h-4 w-4" />
          Design complete! View results in the structure viewer.
        </div>
      )}

      {/* Error message */}
      {status === 'failed' && message && (
        <div className="mt-2 text-sm text-red-600 bg-red-100/50 p-2 rounded-lg">
          {message}
        </div>
      )}
    </div>
  );
}
