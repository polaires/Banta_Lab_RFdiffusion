import type { NextConfig } from "next";

const nextConfig: NextConfig = {
  // Transpile molstar to handle its internal circular references
  transpilePackages: ['molstar'],

  // Turbopack config for molstar compatibility
  turbopack: {
    resolveExtensions: ['.tsx', '.ts', '.jsx', '.js', '.mjs', '.json'],
  },
};

export default nextConfig;
