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
          color: 'text-info',
          bgColor: 'bg-info/10 border-info/20',
          label: 'Queued',
        };
      case 'running':
        return {
          color: 'text-primary',
          bgColor: 'bg-primary/5 border-primary/20',
          label: 'Running',
        };
      case 'completed':
        return {
          color: 'text-success',
          bgColor: 'bg-success/10 border-success/20',
          label: 'Completed',
        };
      case 'failed':
        return {
          color: 'text-destructive',
          bgColor: 'bg-destructive/10 border-destructive/20',
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
            <span className="font-medium text-foreground">RFD3 Design</span>
            <span className={`text-xs px-2 py-0.5 rounded-full ${
              status === 'running' ? 'bg-primary/10 text-primary' :
              status === 'completed' ? 'bg-success/10 text-success' :
              status === 'failed' ? 'bg-destructive/10 text-destructive' :
              'bg-info/10 text-info'
            }`}>
              {config.label}
            </span>
          </div>
          <div className="text-xs text-muted-foreground font-mono mt-0.5">
            Job ID: {jobId.slice(0, 12)}...
          </div>
        </div>

        {/* Progress percentage */}
        {status === 'running' && progress !== undefined && (
          <div className="text-right">
            <div className="text-lg font-bold text-primary">{progress}%</div>
          </div>
        )}
      </div>

      {/* Progress bar */}
      {status === 'running' && (
        <div className="mt-3">
          <div className="h-2 bg-primary/20 rounded-full overflow-hidden">
            <div
              className="h-full bg-primary rounded-full transition-all duration-500"
              style={{ width: `${progress || 0}%` }}
            />
          </div>
          {message && (
            <p className="text-xs text-primary mt-2">{message}</p>
          )}
        </div>
      )}

      {/* Completed message */}
      {status === 'completed' && (
        <div className="mt-2 text-sm text-success flex items-center gap-1">
          <Info className="h-4 w-4" />
          Design complete! View results in the structure viewer.
        </div>
      )}

      {/* Error message */}
      {status === 'failed' && message && (
        <div className="mt-2 text-sm text-destructive bg-destructive/10 p-2 rounded-lg">
          {message}
        </div>
      )}
    </div>
  );
}
