import type { NextConfig } from "next";

const nextConfig: NextConfig = {
  // Transpile molstar to handle its internal circular references
  transpilePackages: ['molstar'],

  // Turbopack config for molstar compatibility (dev mode)
  turbopack: {
    resolveExtensions: ['.tsx', '.ts', '.jsx', '.js', '.mjs', '.json'],
  },

  // Exclude molstar from automatic import optimization (prevents tree-shaking issues)
  experimental: {
    optimizePackageImports: [], // Don't optimize any packages (molstar has issues with this)
  },

  // Webpack config for production builds - handle Molstar's complex module structure
  webpack: (config, { isServer, dev }) => {
    if (!isServer) {
      // Mark molstar as having side effects to prevent incorrect tree-shaking
      // This is critical - Molstar modules have side effects that must be preserved
      config.module.rules.push({
        test: /[\\/]node_modules[\\/]molstar[\\/]/,
        sideEffects: true,
      });

      // In production, prevent aggressive optimizations that break Molstar
      if (!dev && config.optimization) {
        // Disable module concatenation for molstar (breaks circular dependencies)
        config.optimization.concatenateModules = false;

        // Ensure used exports analysis doesn't remove Molstar internals
        config.optimization.usedExports = false;

        // Disable inner graph analysis which can break Molstar's dynamic requires
        config.optimization.innerGraph = false;

        // Ensure Molstar isn't put in a chunk that gets optimized away
        if (config.optimization.splitChunks) {
          config.optimization.splitChunks = {
            ...config.optimization.splitChunks,
            cacheGroups: {
              ...(config.optimization.splitChunks as any).cacheGroups,
              molstar: {
                test: /[\\/]node_modules[\\/]molstar[\\/]/,
                name: 'molstar',
                chunks: 'all',
                priority: 20,
                // Don't minimize this chunk
                enforce: true,
              },
            },
          };
        }
      }
    }

    return config;
  },
};

export default nextConfig;
