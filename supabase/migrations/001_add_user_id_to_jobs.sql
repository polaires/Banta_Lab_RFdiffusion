-- Add user_id to jobs table for per-user job tracking
ALTER TABLE jobs ADD COLUMN user_id UUID REFERENCES auth.users(id);
CREATE INDEX idx_jobs_user_id ON jobs(user_id);

-- Replace open RLS with per-user policies (allow NULL for backward compat)
DROP POLICY IF EXISTS "Allow all" ON jobs;
CREATE POLICY "Users can view own jobs" ON jobs
  FOR SELECT USING (auth.uid() = user_id OR user_id IS NULL);
CREATE POLICY "Users can insert own jobs" ON jobs
  FOR INSERT WITH CHECK (auth.uid() = user_id OR user_id IS NULL);
CREATE POLICY "Users can update own jobs" ON jobs
  FOR UPDATE USING (auth.uid() = user_id OR user_id IS NULL);
CREATE POLICY "Users can delete own jobs" ON jobs
  FOR DELETE USING (auth.uid() = user_id OR user_id IS NULL);
