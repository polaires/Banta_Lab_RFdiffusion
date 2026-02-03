/**
 * Vercel Edge Function Proxy for RunPod Serverless
 *
 * This route proxies requests to RunPod Serverless endpoints,
 * hiding the API key from the frontend and maintaining API compatibility.
 *
 * Endpoints:
 *   POST /api/runpod/rfd3/design -> RunPod task: "rfd3"
 *   POST /api/runpod/rf3/predict -> RunPod task: "rf3"
 *   POST /api/runpod/mpnn/design -> RunPod task: "mpnn"
 *   POST /api/runpod/validate/rmsd -> RunPod task: "rmsd"
 *   GET  /api/runpod/jobs/{jobId} -> RunPod status
 *   GET  /api/runpod/health -> RunPod task: "health"
 */

import { NextRequest, NextResponse } from 'next/server';

// Environment variables (set in Vercel dashboard)
const RUNPOD_API_KEY = process.env.RUNPOD_API_KEY;
const RUNPOD_ENDPOINT_ID = process.env.RUNPOD_ENDPOINT_ID;

// Task mapping from API paths to RunPod task names
const TASK_MAP: Record<string, string> = {
  'rfd3/design': 'rfd3',
  'rf3/predict': 'rf3',
  'mpnn/design': 'mpnn',
  'validate/rmsd': 'rmsd',
  'health': 'health',
  'runsync': 'runsync', // Special: task comes from request body
  'async': 'async',     // Special: async submit, frontend polls for result
};

// RunPod status to internal status mapping
function mapRunPodStatus(runpodStatus: string): string {
  const statusMap: Record<string, string> = {
    'IN_QUEUE': 'pending',
    'IN_PROGRESS': 'running',
    'COMPLETED': 'completed',
    'FAILED': 'failed',
    'CANCELLED': 'failed',
    'TIMED_OUT': 'failed',
  };
  return statusMap[runpodStatus] || 'pending';
}

export async function POST(request: NextRequest) {
  // Check configuration
  if (!RUNPOD_API_KEY || !RUNPOD_ENDPOINT_ID) {
    return NextResponse.json(
      { error: 'Serverless not configured. Set RUNPOD_API_KEY and RUNPOD_ENDPOINT_ID.' },
      { status: 503 }
    );
  }

  try {
    const body = await request.json();

    // Extract path after /api/runpod/
    const pathname = request.nextUrl.pathname;
    const pathMatch = pathname.match(/\/api\/runpod\/(.+)/);
    const path = pathMatch ? pathMatch[1] : '';

    // Map path to task
    const task = TASK_MAP[path];
    if (!task) {
      return NextResponse.json(
        { error: `Unknown endpoint: ${path}` },
        { status: 404 }
      );
    }

    // Special handling for runsync - synchronous execution
    if (path === 'runsync') {
      const runpodResponse = await fetch(
        `https://api.runpod.ai/v2/${RUNPOD_ENDPOINT_ID}/runsync`,
        {
          method: 'POST',
          headers: {
            'Authorization': `Bearer ${RUNPOD_API_KEY}`,
            'Content-Type': 'application/json',
          },
          body: JSON.stringify(body), // Pass body as-is, it contains {input: {task, ...}}
        }
      );

      if (!runpodResponse.ok) {
        const errorText = await runpodResponse.text();
        console.error('[RunPod Proxy] runsync error:', runpodResponse.status, errorText);
        return NextResponse.json(
          { error: `RunPod error: ${errorText}` },
          { status: runpodResponse.status }
        );
      }

      const data = await runpodResponse.json();

      // Return full response for sync endpoint
      return NextResponse.json(data);
    }

    // Async submit — same body as runsync but uses /run (returns job ID immediately)
    // Frontend polls /api/runpod/jobs/{id} for the result.
    if (path === 'async') {
      const runpodResponse = await fetch(
        `https://api.runpod.ai/v2/${RUNPOD_ENDPOINT_ID}/run`,
        {
          method: 'POST',
          headers: {
            'Authorization': `Bearer ${RUNPOD_API_KEY}`,
            'Content-Type': 'application/json',
          },
          body: JSON.stringify(body),
        }
      );

      if (!runpodResponse.ok) {
        const errorText = await runpodResponse.text();
        console.error('[RunPod Proxy] async submit error:', runpodResponse.status, errorText);
        return NextResponse.json(
          { error: `RunPod error: ${errorText}` },
          { status: runpodResponse.status }
        );
      }

      const data = await runpodResponse.json();
      return NextResponse.json({ job_id: data.id, status: 'pending' });
    }

    // Submit to RunPod serverless (async)
    const runpodResponse = await fetch(
      `https://api.runpod.ai/v2/${RUNPOD_ENDPOINT_ID}/run`,
      {
        method: 'POST',
        headers: {
          'Authorization': `Bearer ${RUNPOD_API_KEY}`,
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          input: {
            task,
            ...body,
          },
        }),
      }
    );

    if (!runpodResponse.ok) {
      const errorText = await runpodResponse.text();
      console.error('[RunPod Proxy] Error:', runpodResponse.status, errorText);
      return NextResponse.json(
        { error: `RunPod error: ${errorText}` },
        { status: runpodResponse.status }
      );
    }

    const data = await runpodResponse.json();

    // Return job ID in compatible format
    return NextResponse.json({
      job_id: data.id,
      status: 'pending',
      message: `${task} job submitted successfully`,
    });

  } catch (error) {
    console.error('[RunPod Proxy] Error:', error);
    return NextResponse.json(
      { error: error instanceof Error ? error.message : 'Unknown error' },
      { status: 500 }
    );
  }
}

export async function GET(request: NextRequest) {
  // Check configuration
  if (!RUNPOD_API_KEY || !RUNPOD_ENDPOINT_ID) {
    return NextResponse.json(
      { error: 'Serverless not configured' },
      { status: 503 }
    );
  }

  try {
    const pathname = request.nextUrl.pathname;

    // Health check: /api/runpod/health
    if (pathname.endsWith('/health')) {
      const runpodResponse = await fetch(
        `https://api.runpod.ai/v2/${RUNPOD_ENDPOINT_ID}/run`,
        {
          method: 'POST',
          headers: {
            'Authorization': `Bearer ${RUNPOD_API_KEY}`,
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({
            input: { task: 'health' },
          }),
        }
      );

      if (!runpodResponse.ok) {
        return NextResponse.json({
          status: 'unhealthy',
          mode: 'serverless',
          error: 'Failed to reach RunPod endpoint',
        });
      }

      const data = await runpodResponse.json();

      // For health, we want sync response, so poll immediately
      // (health check is fast)
      const statusResponse = await pollForResult(data.id, 30000);
      if (statusResponse.output?.result) {
        return NextResponse.json({
          status: 'healthy',
          mode: 'serverless',
          ...statusResponse.output.result,
        });
      }

      return NextResponse.json({
        status: 'healthy',
        mode: 'serverless',
        gpu_available: true,
      });
    }

    // Job status: /api/runpod/jobs/{jobId}
    const jobMatch = pathname.match(/\/api\/runpod\/jobs\/([^/]+)/);
    if (jobMatch) {
      const jobId = jobMatch[1];

      const runpodResponse = await fetch(
        `https://api.runpod.ai/v2/${RUNPOD_ENDPOINT_ID}/status/${jobId}`,
        {
          headers: {
            'Authorization': `Bearer ${RUNPOD_API_KEY}`,
          },
        }
      );

      if (!runpodResponse.ok) {
        if (runpodResponse.status === 404) {
          return NextResponse.json(
            { error: 'Job not found' },
            { status: 404 }
          );
        }
        return NextResponse.json(
          { error: 'Failed to get job status' },
          { status: runpodResponse.status }
        );
      }

      const data = await runpodResponse.json();

      // Transform RunPod response to internal format
      const result = data.output?.result || data.output;
      const error = data.output?.error || data.error;

      return NextResponse.json({
        job_id: jobId,
        status: mapRunPodStatus(data.status),
        created_at: data.createdAt,
        completed_at: data.completedAt,
        result: result,
        error: error,
      });
    }

    return NextResponse.json(
      { error: 'Unknown endpoint' },
      { status: 404 }
    );

  } catch (error) {
    console.error('[RunPod Proxy] Error:', error);
    return NextResponse.json(
      { error: error instanceof Error ? error.message : 'Unknown error' },
      { status: 500 }
    );
  }
}

// Helper to poll for result (used for sync endpoints like health)
async function pollForResult(jobId: string, timeoutMs: number = 30000): Promise<any> {
  const startTime = Date.now();
  const pollInterval = 1000;

  while (Date.now() - startTime < timeoutMs) {
    const response = await fetch(
      `https://api.runpod.ai/v2/${RUNPOD_ENDPOINT_ID}/status/${jobId}`,
      {
        headers: {
          'Authorization': `Bearer ${RUNPOD_API_KEY}`,
        },
      }
    );

    if (response.ok) {
      const data = await response.json();
      if (data.status === 'COMPLETED' || data.status === 'FAILED') {
        return data;
      }
    }

    await new Promise(resolve => setTimeout(resolve, pollInterval));
  }

  return { status: 'TIMED_OUT' };
}

// Edge runtime for better performance
export const runtime = 'edge';
