/**
 * Next.js API Proxy for Traditional Backend Mode
 *
 * This route proxies requests to a traditional RunPod backend (always-on GPU pod),
 * avoiding CORS issues by routing through Next.js.
 *
 * The backend URL is passed via query parameter or stored in localStorage on frontend.
 *
 * Endpoints:
 *   POST /api/traditional/runsync?url=http://... -> Proxy to backend /runsync
 *   GET  /api/traditional/health?url=http://... -> Proxy health check
 */

import { NextRequest, NextResponse } from 'next/server';
import { Agent } from 'undici';

// Timeout for design requests (15 minutes) - full binder pipeline can take 10+ minutes
// for proteins like GFP (8 designs × structure gen + MPNN + ESM + FastRelax + analysis)
const DESIGN_TIMEOUT_MS = 15 * 60 * 1000;
// Timeout for health checks (2 minutes) - backend is single-threaded so health checks
// can be delayed when processing long-running design jobs
const HEALTH_TIMEOUT_MS = 2 * 60 * 1000;

// Custom undici agent with extended timeouts for long-running design jobs
// Default headersTimeout is 5 minutes, bodyTimeout is 5 minutes
const longTimeoutAgent = new Agent({
  headersTimeout: DESIGN_TIMEOUT_MS,
  bodyTimeout: DESIGN_TIMEOUT_MS,
  keepAliveTimeout: DESIGN_TIMEOUT_MS,
});

export async function POST(request: NextRequest) {
  try {
    const body = await request.json();
    const backendUrl = request.nextUrl.searchParams.get('url');

    if (!backendUrl) {
      return NextResponse.json(
        { error: 'Missing backend URL parameter' },
        { status: 400 }
      );
    }

    // Extract path after /api/traditional/
    const pathname = request.nextUrl.pathname;
    const pathMatch = pathname.match(/\/api\/traditional\/(.+)/);
    const path = pathMatch ? pathMatch[1] : 'runsync';

    // Use AbortController for timeout - designs can take several minutes
    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), DESIGN_TIMEOUT_MS);

    // Forward to backend with extended timeout dispatcher
    const response = await fetch(`${backendUrl}/${path}`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify(body),
      signal: controller.signal,
      // @ts-expect-error - dispatcher is a valid undici option in Node.js 18+
      dispatcher: longTimeoutAgent,
    });

    clearTimeout(timeoutId);

    if (!response.ok) {
      const errorText = await response.text();
      console.error('[Traditional Proxy] Error:', response.status, errorText);
      return NextResponse.json(
        { error: `Backend error: ${errorText}` },
        { status: response.status }
      );
    }

    const data = await response.json();
    return NextResponse.json(data);

  } catch (error) {
    console.error('[Traditional Proxy] Error:', error);
    // Check for timeout/abort error
    if (error instanceof Error && error.name === 'AbortError') {
      return NextResponse.json(
        { error: 'Request timed out - design may still be running on backend' },
        { status: 504 }
      );
    }
    return NextResponse.json(
      { error: error instanceof Error ? error.message : 'Failed to connect to backend' },
      { status: 500 }
    );
  }
}

export async function GET(request: NextRequest) {
  try {
    const backendUrl = request.nextUrl.searchParams.get('url');

    if (!backendUrl) {
      return NextResponse.json(
        { error: 'Missing backend URL parameter' },
        { status: 400 }
      );
    }

    // Extract path after /api/traditional/
    const pathname = request.nextUrl.pathname;
    const pathMatch = pathname.match(/\/api\/traditional\/(.+)/);
    const path = pathMatch ? pathMatch[1] : 'health';

    // For health check, use POST /runsync with health task
    if (path === 'health') {
      const controller = new AbortController();
      const timeoutId = setTimeout(() => controller.abort(), HEALTH_TIMEOUT_MS);

      const response = await fetch(`${backendUrl}/runsync`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ input: { task: 'health' } }),
        signal: controller.signal,
      });

      clearTimeout(timeoutId);

      if (!response.ok) {
        return NextResponse.json({
          status: 'unhealthy',
          error: 'Failed to reach backend',
        });
      }

      const data = await response.json();

      // Extract health info from RunPod response format
      if (data.output?.result) {
        return NextResponse.json({
          status: 'healthy',
          mode: data.output.result.mode,
          gpu_available: data.output.result.gpu_available,
          gpu_name: data.output.result.gpu_name,
          gpu_memory_gb: data.output.result.gpu_memory_gb,
        });
      }

      return NextResponse.json({
        status: 'healthy',
        mode: 'traditional',
      });
    }

    // For other GET requests, forward directly
    const response = await fetch(`${backendUrl}/${path}`, {
      method: 'GET',
      headers: {
        'Content-Type': 'application/json',
      },
    });

    if (!response.ok) {
      const errorText = await response.text();
      return NextResponse.json(
        { error: `Backend error: ${errorText}` },
        { status: response.status }
      );
    }

    const data = await response.json();
    return NextResponse.json(data);

  } catch (error) {
    console.error('[Traditional Proxy] Error:', error);
    return NextResponse.json(
      { error: error instanceof Error ? error.message : 'Failed to connect to backend' },
      { status: 500 }
    );
  }
}

// Use Node.js runtime for better fetch support
export const runtime = 'nodejs';

// Maximum execution time for Vercel serverless functions
// Hobby plan limit is 300 seconds (5 minutes). Long-running jobs use RunPod backend instead.
export const maxDuration = 300;
