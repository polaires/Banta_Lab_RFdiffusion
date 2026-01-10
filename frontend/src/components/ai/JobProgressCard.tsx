'use client';

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
          icon: 'schedule',
          color: 'text-blue-500',
          bgColor: 'bg-blue-50 border-blue-200',
          label: 'Queued',
        };
      case 'running':
        return {
          icon: 'sync',
          color: 'text-violet-500',
          bgColor: 'bg-violet-50 border-violet-200',
          label: 'Running',
        };
      case 'completed':
        return {
          icon: 'check_circle',
          color: 'text-emerald-500',
          bgColor: 'bg-emerald-50 border-emerald-200',
          label: 'Completed',
        };
      case 'failed':
        return {
          icon: 'error',
          color: 'text-red-500',
          bgColor: 'bg-red-50 border-red-200',
          label: 'Failed',
        };
    }
  };

  const config = getStatusConfig();

  return (
    <div className={`rounded-xl p-4 border ${config.bgColor}`}>
      <div className="flex items-center gap-3">
        {/* Status icon with animation for running */}
        <div className={`${config.color}`}>
          {status === 'running' ? (
            <span className="material-symbols-outlined text-xl animate-spin">progress_activity</span>
          ) : (
            <span className="material-symbols-outlined text-xl">{config.icon}</span>
          )}
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
          <span className="material-symbols-outlined text-sm">info</span>
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
