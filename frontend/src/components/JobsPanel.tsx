'use client';

import { useStore } from '@/lib/store';
import { Clock, CheckCircle, XCircle, Loader2, Trash2, Eye } from 'lucide-react';
import api from '@/lib/api';

const statusIcons = {
  pending: <Clock className="w-4 h-4 text-yellow-400" />,
  running: <Loader2 className="w-4 h-4 text-blue-400 animate-spin" />,
  completed: <CheckCircle className="w-4 h-4 text-green-400" />,
  failed: <XCircle className="w-4 h-4 text-red-400" />,
};

const typeColors = {
  rfd3: 'bg-blue-600',
  rf3: 'bg-green-600',
  mpnn: 'bg-purple-600',
};

const typeLabels = {
  rfd3: 'RFD3',
  rf3: 'RF3',
  mpnn: 'MPNN',
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
      <div className="text-center py-12 text-gray-400">
        <p>No jobs yet. Submit a design or prediction to get started.</p>
      </div>
    );
  }

  return (
    <div className="space-y-4">
      <h2 className="text-xl font-bold">Job History</h2>

      <div className="space-y-2">
        {jobs.map((job) => (
          <div
            key={job.id}
            className="flex items-center gap-4 p-4 bg-gray-800 rounded-lg"
          >
            {/* Status */}
            <div className="flex-shrink-0">
              {statusIcons[job.status]}
            </div>

            {/* Type Badge */}
            <span className={`px-2 py-1 text-xs font-medium rounded ${typeColors[job.type]}`}>
              {typeLabels[job.type]}
            </span>

            {/* ID and Time */}
            <div className="flex-1 min-w-0">
              <p className="text-sm font-mono truncate" title={job.id}>
                {job.id.slice(0, 8)}...
              </p>
              <p className="text-xs text-gray-500">
                {new Date(job.createdAt).toLocaleString()}
              </p>
            </div>

            {/* Status Text */}
            <div className="text-sm">
              {job.status === 'failed' && job.error && (
                <span className="text-red-400" title={job.error}>
                  Failed
                </span>
              )}
              {job.status === 'completed' && (
                <span className="text-green-400">Done</span>
              )}
              {job.status === 'running' && (
                <span className="text-blue-400">Running</span>
              )}
              {job.status === 'pending' && (
                <span className="text-yellow-400">Pending</span>
              )}
            </div>

            {/* Actions */}
            <div className="flex gap-2">
              {job.status === 'completed' && job.result && (
                <button
                  onClick={() => handleView(job)}
                  className="p-2 hover:bg-gray-700 rounded transition"
                  title="View result"
                >
                  <Eye className="w-4 h-4" />
                </button>
              )}
              <button
                onClick={() => handleDelete(job.id)}
                className="p-2 hover:bg-gray-700 rounded transition text-red-400"
                title="Delete job"
              >
                <Trash2 className="w-4 h-4" />
              </button>
            </div>
          </div>
        ))}
      </div>
    </div>
  );
}
