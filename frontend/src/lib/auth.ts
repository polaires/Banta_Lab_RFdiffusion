/**
 * Authentication Service using Supabase Auth
 *
 * Setup in Supabase Dashboard:
 * 1. Go to Authentication > Providers > Google
 * 2. Enable Google provider
 * 3. Add your Google OAuth credentials (from Google Cloud Console)
 * 4. Add authorized redirect URL to Google OAuth consent screen
 *
 * SQL Schema (run in Supabase SQL Editor):
 * ```sql
 * -- Create profiles table linked to auth.users
 * CREATE TABLE profiles (
 *   id UUID PRIMARY KEY REFERENCES auth.users(id) ON DELETE CASCADE,
 *   email TEXT,
 *   full_name TEXT,
 *   avatar_url TEXT,
 *   created_at TIMESTAMPTZ DEFAULT NOW(),
 *   updated_at TIMESTAMPTZ DEFAULT NOW()
 * );
 *
 * -- Enable RLS
 * ALTER TABLE profiles ENABLE ROW LEVEL SECURITY;
 *
 * -- Users can read their own profile
 * CREATE POLICY "Users can view own profile" ON profiles
 *   FOR SELECT USING (auth.uid() = id);
 *
 * -- Users can update their own profile
 * CREATE POLICY "Users can update own profile" ON profiles
 *   FOR UPDATE USING (auth.uid() = id);
 *
 * -- Auto-create profile on signup
 * CREATE OR REPLACE FUNCTION handle_new_user()
 * RETURNS TRIGGER AS $$
 * BEGIN
 *   INSERT INTO public.profiles (id, email, full_name, avatar_url)
 *   VALUES (
 *     NEW.id,
 *     NEW.email,
 *     NEW.raw_user_meta_data->>'full_name',
 *     NEW.raw_user_meta_data->>'avatar_url'
 *   );
 *   RETURN NEW;
 * END;
 * $$ LANGUAGE plpgsql SECURITY DEFINER;
 *
 * CREATE TRIGGER on_auth_user_created
 *   AFTER INSERT ON auth.users
 *   FOR EACH ROW EXECUTE FUNCTION handle_new_user();
 *
 * -- Update jobs table to link to users (optional)
 * ALTER TABLE jobs ADD COLUMN user_id UUID REFERENCES auth.users(id);
 * CREATE INDEX idx_jobs_user_id ON jobs(user_id);
 *
 * -- Update jobs RLS to filter by user
 * DROP POLICY IF EXISTS "Allow all" ON jobs;
 * CREATE POLICY "Users can view own jobs" ON jobs
 *   FOR SELECT USING (auth.uid() = user_id OR user_id IS NULL);
 * CREATE POLICY "Users can insert own jobs" ON jobs
 *   FOR INSERT WITH CHECK (auth.uid() = user_id OR user_id IS NULL);
 * CREATE POLICY "Users can update own jobs" ON jobs
 *   FOR UPDATE USING (auth.uid() = user_id OR user_id IS NULL);
 * CREATE POLICY "Users can delete own jobs" ON jobs
 *   FOR DELETE USING (auth.uid() = user_id OR user_id IS NULL);
 * ```
 */

import { createClient, SupabaseClient, User, Session } from '@supabase/supabase-js';

// Types
export interface UserProfile {
  id: string;
  email: string | null;
  full_name: string | null;
  avatar_url: string | null;
  created_at: string;
}

// Get or create Supabase client
let supabaseAuth: SupabaseClient | null = null;

function getSupabaseClient(): SupabaseClient | null {
  if (supabaseAuth) return supabaseAuth;

  const url = process.env.NEXT_PUBLIC_SUPABASE_URL;
  const anonKey = process.env.NEXT_PUBLIC_SUPABASE_ANON_KEY;

  if (!url || !anonKey) {
    console.warn('[Auth] Supabase not configured');
    return null;
  }

  supabaseAuth = createClient(url, anonKey, {
    auth: {
      autoRefreshToken: true,
      persistSession: true,
      detectSessionInUrl: true,
    },
  });

  return supabaseAuth;
}

/**
 * Sign in with Google OAuth
 */
export async function signInWithGoogle(): Promise<{ error: Error | null }> {
  const client = getSupabaseClient();
  if (!client) {
    return { error: new Error('Supabase not configured') };
  }

  const { error } = await client.auth.signInWithOAuth({
    provider: 'google',
    options: {
      redirectTo: typeof window !== 'undefined'
        ? `${window.location.origin}/`
        : undefined,
    },
  });

  return { error: error ? new Error(error.message) : null };
}

/**
 * Sign out current user
 */
export async function signOut(): Promise<{ error: Error | null }> {
  const client = getSupabaseClient();
  if (!client) {
    return { error: new Error('Supabase not configured') };
  }

  const { error } = await client.auth.signOut();
  return { error: error ? new Error(error.message) : null };
}

/**
 * Get current session
 */
export async function getSession(): Promise<Session | null> {
  const client = getSupabaseClient();
  if (!client) return null;

  const { data: { session } } = await client.auth.getSession();
  return session;
}

/**
 * Get current user
 */
export async function getCurrentUser(): Promise<User | null> {
  const client = getSupabaseClient();
  if (!client) return null;

  const { data: { user } } = await client.auth.getUser();
  return user;
}

/**
 * Get user profile from profiles table
 */
export async function getUserProfile(userId: string): Promise<UserProfile | null> {
  const client = getSupabaseClient();
  if (!client) return null;

  const { data, error } = await client
    .from('profiles')
    .select('*')
    .eq('id', userId)
    .single();

  if (error) {
    console.error('[Auth] Error fetching profile:', error);
    return null;
  }

  return data as UserProfile;
}

/**
 * Subscribe to auth state changes
 */
export function onAuthStateChange(
  callback: (event: string, session: Session | null) => void
): (() => void) | null {
  const client = getSupabaseClient();
  if (!client) return null;

  const { data: { subscription } } = client.auth.onAuthStateChange(callback);
  return () => subscription.unsubscribe();
}

/**
 * Check if auth is configured
 */
export function isAuthConfigured(): boolean {
  return !!(
    process.env.NEXT_PUBLIC_SUPABASE_URL &&
    process.env.NEXT_PUBLIC_SUPABASE_ANON_KEY
  );
}
