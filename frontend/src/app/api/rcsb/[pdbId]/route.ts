/**
 * RCSB PDB Proxy API Route
 *
 * This route proxies requests to RCSB to avoid CORS issues.
 * The RCSB server doesn't include CORS headers, so we need to
 * fetch from the server side and return the content.
 *
 * Endpoint:
 *   GET /api/rcsb/[pdbId] -> Fetches PDB file from RCSB
 */

import { NextRequest, NextResponse } from 'next/server';

export async function GET(
  request: NextRequest,
  { params }: { params: Promise<{ pdbId: string }> }
) {
  const { pdbId } = await params;

  // Validate PDB ID format (4 characters, alphanumeric)
  if (!pdbId || !/^[A-Za-z0-9]{4}$/.test(pdbId)) {
    return NextResponse.json(
      { error: 'Invalid PDB ID. Must be exactly 4 alphanumeric characters.' },
      { status: 400 }
    );
  }

  const normalizedPdbId = pdbId.toUpperCase();

  try {
    const rcsbUrl = `https://files.rcsb.org/download/${normalizedPdbId}.pdb`;

    const response = await fetch(rcsbUrl, {
      headers: {
        'Accept': 'text/plain, chemical/x-pdb',
      },
    });

    if (!response.ok) {
      if (response.status === 404) {
        return NextResponse.json(
          { error: `PDB ID "${normalizedPdbId}" not found in RCSB database` },
          { status: 404 }
        );
      }
      return NextResponse.json(
        { error: `Failed to fetch PDB: ${response.statusText}` },
        { status: response.status }
      );
    }

    const pdbContent = await response.text();

    // Return the PDB content as plain text
    return new NextResponse(pdbContent, {
      status: 200,
      headers: {
        'Content-Type': 'text/plain',
        'Content-Disposition': `attachment; filename="${normalizedPdbId}.pdb"`,
        'Cache-Control': 'public, max-age=86400', // Cache for 24 hours
      },
    });
  } catch (error) {
    console.error('[RCSB Proxy] Error fetching PDB:', error);
    return NextResponse.json(
      { error: error instanceof Error ? error.message : 'Failed to fetch PDB from RCSB' },
      { status: 500 }
    );
  }
}

// Edge runtime for better performance
export const runtime = 'edge';
