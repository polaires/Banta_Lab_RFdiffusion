/**
 * Supabase Client for Job Persistence
 *
 * Stores job history for 24+ hours, surviving serverless worker restarts.
 *
 * Setup:
 * 1. Create a Supabase project at https://supabase.com
 * 2. Run the SQL schema below in the SQL editor
 * 3. Add NEXT_PUBLIC_SUPABASE_URL and NEXT_PUBLIC_SUPABASE_ANON_KEY to .env.local
 *
 * SQL Schema:
 * ```sql
 * CREATE TABLE jobs (
 *   id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
 *   runpod_id TEXT UNIQUE NOT NULL,
 *   type TEXT NOT NULL,
 *   status TEXT DEFAULT 'pending',
 *   request JSONB,
 *   result JSONB,
 *   created_at TIMESTAMPTZ DEFAULT NOW(),
 *   completed_at TIMESTAMPTZ
 * );
 *
 * CREATE INDEX idx_jobs_runpod_id ON jobs(runpod_id);
 * CREATE INDEX idx_jobs_created_at ON jobs(created_at DESC);
 *
 * -- Enable Row Level Security (optional, for public access)
 * ALTER TABLE jobs ENABLE ROW LEVEL SECURITY;
 * CREATE POLICY "Allow all" ON jobs FOR ALL USING (true);
 * ```
 */

import { createClient, SupabaseClient } from '@supabase/supabase-js';

// Types
export interface JobRecord {
  id: string;
  runpod_id: string;
  type: 'rfd3' | 'rf3' | 'mpnn';
  status: 'pending' | 'running' | 'completed' | 'failed';
  request: Record<string, any> | null;
  result: Record<string, any> | null;
  created_at: string;
  completed_at: string | null;
}

export type JobInsert = Omit<JobRecord, 'id' | 'created_at'>;
export type JobUpdate = Partial<Omit<JobRecord, 'id' | 'runpod_id' | 'created_at'>>;

// Supabase client singleton
let supabase: SupabaseClient | null = null;

function getSupabase(): SupabaseClient | null {
  if (supabase) return supabase;

  const url = process.env.NEXT_PUBLIC_SUPABASE_URL;
  const anonKey = process.env.NEXT_PUBLIC_SUPABASE_ANON_KEY;

  if (!url || !anonKey) {
    console.warn('[Supabase] Not configured. Job history will not persist.');
    return null;
  }

  supabase = createClient(url, anonKey);
  return supabase;
}

/**
 * Check if Supabase is configured
 */
export function isSupabaseConfigured(): boolean {
  return !!(
    process.env.NEXT_PUBLIC_SUPABASE_URL &&
    process.env.NEXT_PUBLIC_SUPABASE_ANON_KEY
  );
}

/**
 * Save a new job to the database
 */
export async function saveJob(job: {
  runpod_id: string;
  type: 'rfd3' | 'rf3' | 'mpnn';
  request?: Record<string, any>;
}): Promise<JobRecord | null> {
  const client = getSupabase();
  if (!client) return null;

  try {
    const { data, error } = await client
      .from('jobs')
      .insert({
        runpod_id: job.runpod_id,
        type: job.type,
        status: 'pending',
        request: job.request || null,
      })
      .select()
      .single();

    if (error) {
      console.error('[Supabase] Error saving job:', error);
      return null;
    }

    return data as JobRecord;
  } catch (err) {
    console.error('[Supabase] Error saving job:', err);
    return null;
  }
}

/**
 * Update a job's status and result
 */
export async function updateJob(
  runpodId: string,
  updates: JobUpdate
): Promise<JobRecord | null> {
  const client = getSupabase();
  if (!client) return null;

  try {
    const { data, error } = await client
      .from('jobs')
      .update(updates)
      .eq('runpod_id', runpodId)
      .select()
      .single();

    if (error) {
      console.error('[Supabase] Error updating job:', error);
      return null;
    }

    return data as JobRecord;
  } catch (err) {
    console.error('[Supabase] Error updating job:', err);
    return null;
  }
}

/**
 * Get a job by RunPod ID
 */
export async function getJob(runpodId: string): Promise<JobRecord | null> {
  const client = getSupabase();
  if (!client) return null;

  try {
    const { data, error } = await client
      .from('jobs')
      .select('*')
      .eq('runpod_id', runpodId)
      .single();

    if (error) {
      if (error.code === 'PGRST116') return null; // Not found
      console.error('[Supabase] Error getting job:', error);
      return null;
    }

    return data as JobRecord;
  } catch (err) {
    console.error('[Supabase] Error getting job:', err);
    return null;
  }
}

/**
 * Get recent jobs, optionally filtered by type
 */
export async function getJobs(options?: {
  limit?: number;
  type?: 'rfd3' | 'rf3' | 'mpnn';
  status?: string;
}): Promise<JobRecord[]> {
  const client = getSupabase();
  if (!client) return [];

  try {
    let query = client
      .from('jobs')
      .select('*')
      .order('created_at', { ascending: false })
      .limit(options?.limit || 50);

    if (options?.type) {
      query = query.eq('type', options.type);
    }

    if (options?.status) {
      query = query.eq('status', options.status);
    }

    const { data, error } = await query;

    if (error) {
      console.error('[Supabase] Error getting jobs:', error);
      return [];
    }

    return (data || []) as JobRecord[];
  } catch (err) {
    console.error('[Supabase] Error getting jobs:', err);
    return [];
  }
}

/**
 * Delete a job by RunPod ID
 */
export async function deleteJob(runpodId: string): Promise<boolean> {
  const client = getSupabase();
  if (!client) return false;

  try {
    const { error } = await client
      .from('jobs')
      .delete()
      .eq('runpod_id', runpodId);

    if (error) {
      console.error('[Supabase] Error deleting job:', error);
      return false;
    }

    return true;
  } catch (err) {
    console.error('[Supabase] Error deleting job:', err);
    return false;
  }
}

/**
 * Delete old jobs (cleanup)
 */
export async function deleteOldJobs(olderThanDays: number = 7): Promise<number> {
  const client = getSupabase();
  if (!client) return 0;

  try {
    const cutoffDate = new Date();
    cutoffDate.setDate(cutoffDate.getDate() - olderThanDays);

    const { data, error } = await client
      .from('jobs')
      .delete()
      .lt('created_at', cutoffDate.toISOString())
      .select('id');

    if (error) {
      console.error('[Supabase] Error deleting old jobs:', error);
      return 0;
    }

    return data?.length || 0;
  } catch (err) {
    console.error('[Supabase] Error deleting old jobs:', err);
    return 0;
  }
}
