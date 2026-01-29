'use client';

import { useEffect, useState, useMemo, useCallback } from 'react';
import { useStore } from '@/lib/store';
import api from '@/lib/api';
import { getJobs as getJobsFromSupabase, deleteJob as deleteJobFromSupabase, isSupabaseConfigured } from '@/lib/supabase';
import { ErrorDetails } from './ErrorDetails';
import { useAuth } from '@/contexts/AuthContext';
import { Clock, Loader2, CheckCircle, AlertCircle, Building2, Brain, Type, History, Search, Inbox, ChevronDown, Timer, Eye, Trash2, type LucideIcon } from 'lucide-react';

// Status display configuration
const statusConfig: Record<string, { Icon: LucideIcon; color: string; bg: string; label: string; animate?: boolean }> = {
  pending: { Icon: Clock, color: 'text-muted-foreground', bg: 'bg-muted', label: 'Pending' },
  running: { Icon: Loader2, color: 'text-primary', bg: 'bg-primary/10', label: 'Running', animate: true },
  completed: { Icon: CheckCircle, color: 'text-primary', bg: 'bg-primary/10', label: 'Completed' },
  failed: { Icon: AlertCircle, color: 'text-destructive', bg: 'bg-destructive/10', label: 'Failed' },
};

// Type display configuration
const typeConfig: Record<string, { color: string; label: string; Icon: LucideIcon }> = {
  rfd3: { color: 'bg-primary', label: 'RFD3', Icon: Building2 },
  rf3: { color: 'bg-primary', label: 'RF3', Icon: Brain },
  mpnn: { color: 'bg-primary', label: 'MPNN', Icon: Type },
};

type SortOption = 'date' | 'status' | 'type';

// Helper to group jobs by date
function groupJobsByDate<T extends { createdAt: string }>(jobs: T[]): { label: string; jobs: T[] }[] {
  const now = new Date();
  const today = new Date(now.getFullYear(), now.getMonth(), now.getDate());
  const yesterday = new Date(today.getTime() - 24 * 60 * 60 * 1000);
  const thisWeek = new Date(today.getTime() - 7 * 24 * 60 * 60 * 1000);

  const groups: { label: string; jobs: T[] }[] = [
    { label: 'Today', jobs: [] },
    { label: 'Yesterday', jobs: [] },
    { label: 'This Week', jobs: [] },
    { label: 'Older', jobs: [] },
  ];

  jobs.forEach((job) => {
    const jobDate = new Date(job.createdAt);
    if (jobDate >= today) {
      groups[0].jobs.push(job);
    } else if (jobDate >= yesterday) {
      groups[1].jobs.push(job);
    } else if (jobDate >= thisWeek) {
      groups[2].jobs.push(job);
    } else {
      groups[3].jobs.push(job);
    }
  });

  return groups.filter((g) => g.jobs.length > 0);
}

// Helper to calculate expiration time for failed jobs
function getExpirationInfo(createdAt: string): { remaining: string; isExpiringSoon: boolean } | null {
  const created = new Date(createdAt).getTime();
  const expiresAt = created + 24 * 60 * 60 * 1000; // 24 hours
  const now = Date.now();
  const remaining = expiresAt - now;

  if (remaining <= 0) {
    return { remaining: 'Expired', isExpiringSoon: true };
  }

  const hours = Math.floor(remaining / (60 * 60 * 1000));
  const minutes = Math.floor((remaining % (60 * 60 * 1000)) / (60 * 1000));

  if (hours > 0) {
    return { remaining: `${hours}h ${minutes}m`, isExpiringSoon: hours < 2 };
  }
  return { remaining: `${minutes}m`, isExpiringSoon: true };
}

export function DesignHistoryPanel() {
  const { jobs, removeJob, setSelectedPdb, addJob, cleanupExpiredJobs, clearJobs } = useStore();
  const { user, isConfigured: authConfigured } = useAuth();
  const [loadingHistory, setLoadingHistory] = useState(false);
  const [historyLoaded, setHistoryLoaded] = useState(false);
  const [sortBy, setSortBy] = useState<SortOption>('date');
  const [filterStatus, setFilterStatus] = useState<string | null>(null);
  const [searchQuery, setSearchQuery] = useState('');
  const [collapsedGroups, setCollapsedGroups] = useState<Set<string>>(new Set());

  // Run cleanup on mount and every hour
  useEffect(() => {
    cleanupExpiredJobs();
    const interval = setInterval(cleanupExpiredJobs, 60 * 60 * 1000); // Every hour
    return () => clearInterval(interval);
  }, [cleanupExpiredJobs]);

  // Reset history when user changes (sign in/out)
  const [lastUserId, setLastUserId] = useState<string | null | undefined>(undefined);
  useEffect(() => {
    const currentUserId = user?.id ?? null;
    if (lastUserId !== undefined && lastUserId !== currentUserId) {
      // User changed — clear local jobs and re-fetch from Supabase
      clearJobs();
      setHistoryLoaded(false);
    }
    setLastUserId(currentUserId);
  }, [user?.id, lastUserId, clearJobs]);

  // Load job history from Supabase on mount and when user changes
  useEffect(() => {
    if (historyLoaded || !isSupabaseConfigured()) return;

    const loadHistory = async () => {
      setLoadingHistory(true);
      try {
        // RLS will automatically filter by user when authenticated
        const supabaseJobs = await getJobsFromSupabase({ limit: 100 });
        const existingIds = new Set(jobs.map(j => j.id));

        for (const sj of supabaseJobs) {
          if (!existingIds.has(sj.runpod_id)) {
            addJob({
              id: sj.runpod_id,
              type: sj.type,
              status: sj.status,
              createdAt: sj.created_at,
              completedAt: sj.completed_at || undefined,
              result: sj.result || undefined,
            });
          }
        }
      } catch (err) {
        console.error('[DesignHistoryPanel] Failed to load history:', err);
      } finally {
        setLoadingHistory(false);
        setHistoryLoaded(true);
      }
    };

    loadHistory();
  }, [historyLoaded, jobs, addJob]);

  // Filter and sort jobs
  const processedJobs = useMemo(() => {
    let filtered = [...jobs];

    // Apply status filter
    if (filterStatus) {
      filtered = filtered.filter((j) => j.status === filterStatus);
    }

    // Apply search filter
    if (searchQuery) {
      const query = searchQuery.toLowerCase();
      filtered = filtered.filter((j) =>
        j.id.toLowerCase().includes(query) ||
        j.type.toLowerCase().includes(query)
      );
    }

    // Sort
    switch (sortBy) {
      case 'date':
        filtered.sort((a, b) => new Date(b.createdAt).getTime() - new Date(a.createdAt).getTime());
        break;
      case 'status':
        const statusOrder = { running: 0, pending: 1, completed: 2, failed: 3 };
        filtered.sort((a, b) => (statusOrder[a.status] ?? 4) - (statusOrder[b.status] ?? 4));
        break;
      case 'type':
        filtered.sort((a, b) => a.type.localeCompare(b.type));
        break;
    }

    return filtered;
  }, [jobs, filterStatus, searchQuery, sortBy]);

  // Group jobs by date
  const groupedJobs = useMemo(() => groupJobsByDate(processedJobs), [processedJobs]);

  const handleDelete = useCallback(async (jobId: string) => {
    try {
      await api.deleteJob(jobId);
    } catch {
      // Job might not exist on backend
    }
    await deleteJobFromSupabase(jobId);
    removeJob(jobId);
  }, [removeJob]);

  const handleView = useCallback((job: typeof jobs[0]) => {
    if (job.result?.designs?.[0]) {
      const design = job.result.designs[0];
      // Handle both standard RFD3 (content) and interface_ligand (pdb_content via any)
      const designAny = design as unknown as Record<string, string>;
      const pdbContent = design.content || designAny.pdb_content;
      if (pdbContent) setSelectedPdb(pdbContent);
    } else if (job.result?.predictions?.[0]) {
      setSelectedPdb(job.result.predictions[0].content);
    } else {
      // Handle interface_ligand full dimer results
      const resultAny = job.result as unknown as Record<string, Record<string, string>>;
      if (resultAny?.dimer?.pdb_content) {
        setSelectedPdb(resultAny.dimer.pdb_content);
      }
    }
  }, [setSelectedPdb]);

  const toggleGroup = (label: string) => {
    setCollapsedGroups((prev) => {
      const next = new Set(prev);
      if (next.has(label)) {
        next.delete(label);
      } else {
        next.add(label);
      }
      return next;
    });
  };

  // Stats
  const stats = useMemo(() => ({
    total: jobs.length,
    completed: jobs.filter((j) => j.status === 'completed').length,
    failed: jobs.filter((j) => j.status === 'failed').length,
    running: jobs.filter((j) => j.status === 'running' || j.status === 'pending').length,
  }), [jobs]);

  if (loadingHistory) {
    return (
      <div className="p-8">
        <div className="text-center py-16">
          <div className="w-16 h-16 rounded-2xl bg-muted flex items-center justify-center mx-auto mb-4">
            <Loader2 className="w-8 h-8 text-muted-foreground animate-spin" />
          </div>
          <p className="text-muted-foreground font-medium">Loading design history...</p>
        </div>
      </div>
    );
  }

  return (
    <div className="p-8 space-y-6">
      {/* Header */}
      <div className="border-b border-border pb-6">
        <div className="flex items-center gap-3 mb-2">
          <History className="w-6 h-6 text-primary" />
          <h2 className="text-xl font-bold text-foreground">Design History</h2>
        </div>
        <p className="text-muted-foreground text-sm">
          {authConfigured && user
            ? `Designs for ${user.email}. Failed designs are automatically removed after 24 hours.`
            : 'All your protein designs in one place. Failed designs are automatically removed after 24 hours.'}
        </p>
      </div>

      {/* Stats Bar */}
      <div className="grid grid-cols-4 gap-3">
        <div className="bg-muted rounded-xl p-3 text-center">
          <p className="text-2xl font-bold text-foreground">{stats.total}</p>
          <p className="text-xs text-muted-foreground">Total</p>
        </div>
        <div className="bg-primary/10 rounded-xl p-3 text-center">
          <p className="text-2xl font-bold text-primary">{stats.completed}</p>
          <p className="text-xs text-primary">Completed</p>
        </div>
        <div className="bg-muted rounded-xl p-3 text-center">
          <p className="text-2xl font-bold text-foreground">{stats.running}</p>
          <p className="text-xs text-muted-foreground">In Progress</p>
        </div>
        <div className="bg-destructive/10 rounded-xl p-3 text-center">
          <p className="text-2xl font-bold text-destructive">{stats.failed}</p>
          <p className="text-xs text-destructive">Failed</p>
        </div>
      </div>

      {/* Filters and Search */}
      <div className="flex flex-col sm:flex-row gap-3">
        {/* Search */}
        <div className="flex-1 relative">
          <Search className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-muted-foreground" />
          <input
            type="text"
            value={searchQuery}
            onChange={(e) => setSearchQuery(e.target.value)}
            placeholder="Search by ID or type..."
            className="w-full pl-10 pr-4 py-2 border border-border rounded-lg text-sm focus:ring-2 focus:ring-ring focus:border-transparent bg-card"
          />
        </div>

        {/* Sort dropdown */}
        <select
          value={sortBy}
          onChange={(e) => setSortBy(e.target.value as SortOption)}
          className="px-3 py-2 border border-border rounded-lg text-sm bg-card focus:ring-2 focus:ring-ring"
        >
          <option value="date">Sort by Date</option>
          <option value="status">Sort by Status</option>
          <option value="type">Sort by Type</option>
        </select>

        {/* Status filter */}
        <select
          value={filterStatus || ''}
          onChange={(e) => setFilterStatus(e.target.value || null)}
          className="px-3 py-2 border border-border rounded-lg text-sm bg-card focus:ring-2 focus:ring-ring"
        >
          <option value="">All Statuses</option>
          <option value="completed">Completed</option>
          <option value="running">Running</option>
          <option value="pending">Pending</option>
          <option value="failed">Failed</option>
        </select>
      </div>

      {/* Empty State */}
      {jobs.length === 0 && (
        <div className="text-center py-16">
          <div className="w-16 h-16 rounded-2xl bg-muted flex items-center justify-center mx-auto mb-4">
            <Inbox className="w-8 h-8 text-muted-foreground" />
          </div>
          <p className="text-foreground font-medium">No designs yet</p>
          <p className="text-muted-foreground text-sm mt-1">Submit a design task to get started</p>
        </div>
      )}

      {/* Grouped Job List */}
      {groupedJobs.map((group) => (
        <div key={group.label} className="space-y-2">
          {/* Group Header */}
          <button
            onClick={() => toggleGroup(group.label)}
            className="w-full flex items-center justify-between py-2 text-left"
          >
            <div className="flex items-center gap-2">
              <ChevronDown
                className={`w-4 h-4 text-muted-foreground transition-transform ${
                  collapsedGroups.has(group.label) ? '-rotate-90' : ''
                }`}
              />
              <span className="text-sm font-semibold text-foreground">{group.label}</span>
              <span className="text-xs text-muted-foreground bg-muted px-2 py-0.5 rounded-full">
                {group.jobs.length}
              </span>
            </div>
          </button>

          {/* Group Content */}
          {!collapsedGroups.has(group.label) && (
            <div className="space-y-2 ml-6">
              {group.jobs.map((job) => {
                const status = statusConfig[job.status];
                const type = typeConfig[job.type];
                const hasFailed = job.status === 'failed';
                const expiration = hasFailed ? getExpirationInfo(job.createdAt) : null;
                const StatusIcon = status.Icon;
                const TypeIcon = type.Icon;

                return (
                  <div
                    key={job.id}
                    className={`rounded-xl transition-colors overflow-hidden ${
                      hasFailed ? 'bg-destructive/5 border border-destructive/20' : 'bg-muted hover:bg-muted/80'
                    }`}
                  >
                    <div className="flex items-center gap-3 p-3">
                      {/* Status Icon */}
                      <div className={`w-9 h-9 rounded-lg ${status.bg} flex items-center justify-center flex-shrink-0`}>
                        <StatusIcon className={`w-4 h-4 ${status.color} ${status.animate ? 'animate-spin' : ''}`} />
                      </div>

                      {/* Type Badge */}
                      <div className={`flex items-center gap-1.5 px-2 py-1 text-xs font-semibold text-white rounded-lg ${type.color}`}>
                        <TypeIcon className="w-3 h-3" />
                        {type.label}
                      </div>

                      {/* ID and Time */}
                      <div className="flex-1 min-w-0">
                        <p className="text-sm font-mono text-foreground truncate" title={job.id}>
                          {job.id.slice(0, 12)}...
                        </p>
                        <p className="text-xs text-muted-foreground">
                          {new Date(job.createdAt).toLocaleTimeString()}
                        </p>
                      </div>

                      {/* Expiration for failed jobs */}
                      {expiration && (
                        <div className={`text-xs px-2 py-1 rounded-lg flex items-center gap-1 ${
                          expiration.isExpiringSoon ? 'bg-destructive/10 text-destructive' : 'bg-muted text-muted-foreground'
                        }`}>
                          <Timer className="w-3 h-3" />
                          {expiration.remaining}
                        </div>
                      )}

                      {/* Status Badge */}
                      <div className={`px-2.5 py-1 rounded-lg text-xs font-medium ${status.bg} ${status.color}`}>
                        {status.label}
                      </div>

                      {/* Actions */}
                      <div className="flex gap-1">
                        {job.status === 'completed' && job.result && (
                          <button
                            onClick={() => handleView(job)}
                            className="p-1.5 hover:bg-card rounded-lg transition-colors text-muted-foreground hover:text-primary"
                            title="View result"
                          >
                            <Eye className="w-4 h-4" />
                          </button>
                        )}
                        <button
                          onClick={() => handleDelete(job.id)}
                          className="p-1.5 hover:bg-card rounded-lg transition-colors text-muted-foreground hover:text-red-500"
                          title="Delete"
                        >
                          <Trash2 className="w-4 h-4" />
                        </button>
                      </div>
                    </div>

                    {/* Error Details for Failed Jobs */}
                    {hasFailed && job.error && (
                      <div className="px-3 pb-3">
                        <ErrorDetails
                          error={job.error}
                          errorType={job.errorType}
                          traceback={job.traceback}
                          context={job.errorContext}
                          compact
                        />
                      </div>
                    )}
                  </div>
                );
              })}
            </div>
          )}
        </div>
      ))}
    </div>
  );
}

export default DesignHistoryPanel;
