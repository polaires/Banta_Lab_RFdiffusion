'use client';

import { useStore } from '@/lib/store';
import api from '@/lib/api';

const statusConfig: Record<string, { icon: string; color: string; bg: string; label: string; animate?: boolean }> = {
  pending: { icon: 'schedule', color: 'text-amber-500', bg: 'bg-amber-50', label: 'Pending' },
  running: { icon: 'progress_activity', color: 'text-blue-500', bg: 'bg-blue-50', label: 'Running', animate: true },
  completed: { icon: 'check_circle', color: 'text-emerald-500', bg: 'bg-emerald-50', label: 'Done' },
  failed: { icon: 'error', color: 'text-red-500', bg: 'bg-red-50', label: 'Failed' },
};

const typeConfig = {
  rfd3: { color: 'bg-blue-600', label: 'RFD3' },
  rf3: { color: 'bg-emerald-600', label: 'RF3' },
  mpnn: { color: 'bg-violet-600', label: 'MPNN' },
};

export function JobsPanel() {
  const { jobs, removeJob, setSelectedPdb } = useStore();

  const handleDelete = async (jobId: string) => {
    try {
      await api.deleteJob(jobId);
    } catch {
      // Job might not exist on backend, that's ok
    }
    removeJob(jobId);
  };

  const handleView = (job: typeof jobs[0]) => {
    if (job.result?.designs?.[0]) {
      setSelectedPdb(job.result.designs[0].content);
    } else if (job.result?.predictions?.[0]) {
      setSelectedPdb(job.result.predictions[0].content);
    }
  };

  if (jobs.length === 0) {
    return (
      <div className="text-center py-16">
        <div className="w-16 h-16 rounded-2xl bg-slate-100 flex items-center justify-center mx-auto mb-4">
          <span className="material-symbols-outlined text-3xl text-slate-400">inbox</span>
        </div>
        <p className="text-slate-500 font-medium">No jobs yet</p>
        <p className="text-slate-400 text-sm mt-1">Submit a design or prediction to get started</p>
      </div>
    );
  }

  return (
    <div className="space-y-6">
      <div className="flex items-center justify-between">
        <div className="flex items-center gap-3">
          <span className="material-symbols-outlined text-slate-600">history</span>
          <h2 className="text-lg font-bold text-slate-900">Job History</h2>
        </div>
        <span className="text-xs font-medium text-slate-500 bg-slate-100 px-2 py-1 rounded-lg">
          {jobs.length} job{jobs.length !== 1 ? 's' : ''}
        </span>
      </div>

      <div className="space-y-3">
        {jobs.map((job) => {
          const status = statusConfig[job.status];
          const type = typeConfig[job.type];

          return (
            <div
              key={job.id}
              className="flex items-center gap-4 p-4 bg-slate-50 hover:bg-slate-100 rounded-xl transition-colors"
            >
              {/* Status Icon */}
              <div className={`w-10 h-10 rounded-xl ${status.bg} flex items-center justify-center flex-shrink-0`}>
                <span className={`material-symbols-outlined ${status.color} ${status.animate ? 'animate-spin' : ''}`}>
                  {status.icon}
                </span>
              </div>

              {/* Type Badge */}
              <span className={`px-2.5 py-1 text-xs font-semibold text-white rounded-lg ${type.color}`}>
                {type.label}
              </span>

              {/* ID and Time */}
              <div className="flex-1 min-w-0">
                <p className="text-sm font-mono text-slate-700 truncate" title={job.id}>
                  {job.id.slice(0, 8)}...
                </p>
                <p className="text-xs text-slate-500">
                  {new Date(job.createdAt).toLocaleString()}
                </p>
              </div>

              {/* Status Badge */}
              <div className={`px-3 py-1 rounded-lg text-xs font-medium ${status.bg} ${status.color}`}>
                {job.status === 'failed' && job.error ? (
                  <span title={job.error}>{status.label}</span>
                ) : (
                  status.label
                )}
              </div>

              {/* Actions */}
              <div className="flex gap-1">
                {job.status === 'completed' && job.result && (
                  <button
                    onClick={() => handleView(job)}
                    className="p-2 hover:bg-white rounded-lg transition-colors text-slate-500 hover:text-blue-600"
                    title="View result"
                  >
                    <span className="material-symbols-outlined text-xl">visibility</span>
                  </button>
                )}
                <button
                  onClick={() => handleDelete(job.id)}
                  className="p-2 hover:bg-white rounded-lg transition-colors text-slate-400 hover:text-red-500"
                  title="Delete job"
                >
                  <span className="material-symbols-outlined text-xl">delete</span>
                </button>
              </div>
            </div>
          );
        })}
      </div>
    </div>
  );
}
