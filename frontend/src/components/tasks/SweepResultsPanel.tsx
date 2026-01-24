'use client';

import {
  Trophy,
  TrendingUp,
  BarChart3,
  ChevronDown,
  ChevronUp,
  CheckCircle,
  AlertCircle,
  XCircle,
} from 'lucide-react';
import { useState } from 'react';
import { cn } from '@/lib/utils';

import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { Button } from '@/components/ui/button';
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from '@/components/ui/table';
import { Collapsible, CollapsibleContent, CollapsibleTrigger } from '@/components/ui/collapsible';

// Quality tier colors
const TIER_COLORS: Record<string, string> = {
  S: 'bg-yellow-500 text-black',
  A: 'bg-green-500 text-white',
  B: 'bg-blue-500 text-white',
  C: 'bg-orange-500 text-white',
  F: 'bg-red-500 text-white',
};

interface ConfigRanking {
  rank: number;
  config_name: string;
  contig_range: string;
  cfg_scale: number;
  total_designs: number;
  passing: number;
  review: number;
  failed: number;
  pass_rate: number;
  avg_plddt: number;
  avg_ptm: number;
  tier_distribution: Record<string, number>;
}

interface DesignResult {
  name: string;
  config_name: string;
  sequence: string;
  plddt: number;
  ptm: number;
  pae: number;
  tier: string;
  status: 'pass' | 'review' | 'fail';
  seed?: number;
}

interface SweepResultsPanelProps {
  configRankings: ConfigRanking[];
  results: DesignResult[];
  bestDesign?: { name: string; plddt: number; ptm: number; tier: string } | null;
  onRunProduction?: (configName: string, config: { contig_range: string; cfg_scale: number }) => void;
}

export function SweepResultsPanel({
  configRankings,
  results,
  bestDesign,
  onRunProduction,
}: SweepResultsPanelProps) {
  const [expandedConfig, setExpandedConfig] = useState<string | null>(null);
  const [showAllConfigs, setShowAllConfigs] = useState(false);

  // Get results for a specific config
  const getConfigResults = (configName: string) => {
    return results.filter((r) => r.config_name === configName);
  };

  // Display limited configs unless expanded
  const displayedConfigs = showAllConfigs ? configRankings : configRankings.slice(0, 5);

  return (
    <div className="space-y-4">
      {/* Best config highlight */}
      {configRankings.length > 0 && (
        <Card className="border-primary/50 bg-primary/5">
          <CardHeader className="pb-3">
            <CardTitle className="text-sm flex items-center gap-2">
              <Trophy className="w-4 h-4 text-primary" />
              Best Configuration
            </CardTitle>
          </CardHeader>
          <CardContent>
            <div className="flex items-center justify-between">
              <div>
                <p className="font-mono font-medium">{configRankings[0].config_name}</p>
                <p className="text-sm text-muted-foreground">
                  {configRankings[0].contig_range} residues, CFG {configRankings[0].cfg_scale.toFixed(1)}
                </p>
              </div>
              <div className="text-right">
                <p className="text-2xl font-bold text-primary">
                  {(configRankings[0].pass_rate * 100).toFixed(0)}%
                </p>
                <p className="text-xs text-muted-foreground">pass rate</p>
              </div>
              {onRunProduction && (
                <Button
                  size="sm"
                  onClick={() => onRunProduction(configRankings[0].config_name, {
                    contig_range: configRankings[0].contig_range,
                    cfg_scale: configRankings[0].cfg_scale,
                  })}
                >
                  Run Production
                </Button>
              )}
            </div>

            {/* Tier distribution */}
            <div className="mt-4 flex items-center gap-2">
              <span className="text-xs text-muted-foreground">Quality tiers:</span>
              {Object.entries(configRankings[0].tier_distribution).map(([tier, count]) => (
                count > 0 && (
                  <Badge key={tier} variant="secondary" className={cn('text-xs', TIER_COLORS[tier])}>
                    {tier}: {count}
                  </Badge>
                )
              ))}
            </div>
          </CardContent>
        </Card>
      )}

      {/* Config rankings table */}
      <Card>
        <CardHeader>
          <CardTitle className="text-sm flex items-center gap-2">
            <BarChart3 className="w-4 h-4" />
            Configuration Rankings
          </CardTitle>
          <CardDescription>
            Ranked by pass rate and average pLDDT
          </CardDescription>
        </CardHeader>
        <CardContent>
          <Table>
            <TableHeader>
              <TableRow>
                <TableHead className="w-12">Rank</TableHead>
                <TableHead>Config</TableHead>
                <TableHead>Size</TableHead>
                <TableHead>CFG</TableHead>
                <TableHead className="text-right">Pass Rate</TableHead>
                <TableHead className="text-right">Avg pLDDT</TableHead>
                <TableHead className="text-right">Avg pTM</TableHead>
                <TableHead className="w-24"></TableHead>
              </TableRow>
            </TableHeader>
            <TableBody>
              {displayedConfigs.map((config, idx) => (
                <Collapsible
                  key={config.config_name}
                  open={expandedConfig === config.config_name}
                  onOpenChange={(open) => setExpandedConfig(open ? config.config_name : null)}
                >
                  <TableRow
                    className={cn(
                      'cursor-pointer hover:bg-muted/50',
                      idx === 0 && 'bg-primary/5'
                    )}
                  >
                    <TableCell>
                      <span className={cn(
                        'inline-flex items-center justify-center w-6 h-6 rounded-full text-xs font-medium',
                        idx === 0 ? 'bg-primary text-primary-foreground' : 'bg-muted'
                      )}>
                        {config.rank}
                      </span>
                    </TableCell>
                    <TableCell className="font-mono text-sm">{config.config_name}</TableCell>
                    <TableCell className="text-sm">{config.contig_range}</TableCell>
                    <TableCell className="text-sm">{config.cfg_scale.toFixed(1)}</TableCell>
                    <TableCell className="text-right">
                      <span className={cn(
                        'font-medium',
                        config.pass_rate >= 0.5 ? 'text-green-600' :
                        config.pass_rate >= 0.3 ? 'text-yellow-600' : 'text-red-600'
                      )}>
                        {(config.pass_rate * 100).toFixed(0)}%
                      </span>
                    </TableCell>
                    <TableCell className="text-right font-mono text-sm">
                      {config.avg_plddt.toFixed(2)}
                    </TableCell>
                    <TableCell className="text-right font-mono text-sm">
                      {config.avg_ptm.toFixed(2)}
                    </TableCell>
                    <TableCell>
                      <CollapsibleTrigger asChild>
                        <Button variant="ghost" size="sm">
                          {expandedConfig === config.config_name ? (
                            <ChevronUp className="w-4 h-4" />
                          ) : (
                            <ChevronDown className="w-4 h-4" />
                          )}
                        </Button>
                      </CollapsibleTrigger>
                    </TableCell>
                  </TableRow>
                  <CollapsibleContent asChild>
                    <TableRow className="bg-muted/30">
                      <TableCell colSpan={8} className="p-4">
                        <div className="space-y-3">
                          {/* Tier distribution */}
                          <div className="flex items-center gap-4">
                            <span className="text-sm text-muted-foreground">Tier distribution:</span>
                            {Object.entries(config.tier_distribution).map(([tier, count]) => (
                              <div key={tier} className="flex items-center gap-1">
                                <Badge className={cn('text-xs', TIER_COLORS[tier])}>{tier}</Badge>
                                <span className="text-sm">{count}</span>
                              </div>
                            ))}
                          </div>

                          {/* Design results for this config */}
                          <div className="text-sm text-muted-foreground">
                            {config.passing} passing, {config.review} review, {config.failed} failed
                          </div>

                          {/* Top designs from this config */}
                          <div className="grid grid-cols-3 gap-2">
                            {getConfigResults(config.config_name)
                              .filter((r) => r.status === 'pass')
                              .slice(0, 3)
                              .map((design) => (
                                <div
                                  key={design.name}
                                  className="p-2 rounded bg-background border text-xs"
                                >
                                  <div className="flex items-center justify-between">
                                    <span className="font-mono truncate">{design.name}</span>
                                    <Badge className={cn('text-xs ml-1', TIER_COLORS[design.tier])}>
                                      {design.tier}
                                    </Badge>
                                  </div>
                                  <div className="text-muted-foreground mt-1">
                                    pLDDT: {design.plddt.toFixed(2)} | pTM: {design.ptm.toFixed(2)}
                                  </div>
                                </div>
                              ))}
                          </div>
                        </div>
                      </TableCell>
                    </TableRow>
                  </CollapsibleContent>
                </Collapsible>
              ))}
            </TableBody>
          </Table>

          {configRankings.length > 5 && (
            <Button
              variant="ghost"
              className="w-full mt-2"
              onClick={() => setShowAllConfigs(!showAllConfigs)}
            >
              {showAllConfigs ? 'Show less' : `Show all ${configRankings.length} configs`}
            </Button>
          )}
        </CardContent>
      </Card>

      {/* Best design highlight */}
      {bestDesign && (
        <Card className="border-green-500/50 bg-green-500/5">
          <CardHeader className="pb-3">
            <CardTitle className="text-sm flex items-center gap-2">
              <CheckCircle className="w-4 h-4 text-green-500" />
              Best Overall Design
            </CardTitle>
          </CardHeader>
          <CardContent>
            <div className="flex items-center gap-4">
              <Badge className={cn(TIER_COLORS[bestDesign.tier])}>
                Tier {bestDesign.tier}
              </Badge>
              <span className="font-mono">{bestDesign.name}</span>
              <span className="text-sm text-muted-foreground">
                pLDDT: {bestDesign.plddt.toFixed(2)} | pTM: {bestDesign.ptm.toFixed(2)}
              </span>
            </div>
          </CardContent>
        </Card>
      )}
    </div>
  );
}

export default SweepResultsPanel;
