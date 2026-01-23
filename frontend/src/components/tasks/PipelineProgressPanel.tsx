'use client';

import { cn } from '@/lib/utils';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { Progress } from '@/components/ui/progress';
import {
  Activity,
  CheckCircle,
  XCircle,
  AlertTriangle,
  Pause,
  X,
  Trophy,
  TrendingUp,
} from 'lucide-react';
import type { PipelineProgress, PipelineMode } from '@/lib/store';

interface PipelineProgressPanelProps {
  progress: PipelineProgress;
  mode: PipelineMode;
  isRunning: boolean;
  onCancel?: () => void;
  onPause?: () => void;
}

export function PipelineProgressPanel({
  progress,
  mode,
  isRunning,
  onCancel,
  onPause,
}: PipelineProgressPanelProps) {
  // Calculate overall progress percentage
  const totalExpected = progress.totalConfigs * progress.designsPerConfig;
  const overallProgress = totalExpected > 0
    ? (progress.totalGenerated / totalExpected) * 100
    : 0;

  // Current config progress
  const configProgress = progress.designsPerConfig > 0
    ? (progress.currentDesign / progress.designsPerConfig) * 100
    : 0;

  // Get config name from current index
  const getConfigName = () => {
    if (mode === 'production') return 'Production';

    const sizes = ['small', 'medium', 'large'];
    const cfgs = ['low_cfg', 'mid_cfg', 'high_cfg'];

    const sizeIndex = Math.floor((progress.currentConfig - 1) / 3);
    const cfgIndex = (progress.currentConfig - 1) % 3;

    if (sizeIndex < 0 || sizeIndex >= sizes.length) return `Config ${progress.currentConfig}`;
    if (cfgIndex < 0 || cfgIndex >= cfgs.length) return `Config ${progress.currentConfig}`;

    return `${sizes[sizeIndex]}_${cfgs[cfgIndex]}`;
  };

  return (
    <Card className={cn(
      "border-2",
      isRunning ? "border-primary/50" : "border-border"
    )}>
      <CardHeader className="pb-3">
        <CardTitle className="text-sm flex items-center justify-between">
          <div className="flex items-center gap-2">
            <Activity className={cn(
              "w-4 h-4",
              isRunning && "animate-pulse text-primary"
            )} />
            Pipeline Progress
          </div>
          <Badge variant={isRunning ? "default" : "secondary"}>
            {isRunning ? 'Running' : 'Idle'}
          </Badge>
        </CardTitle>
      </CardHeader>
      <CardContent className="space-y-4">
        {/* Current Configuration */}
        {isRunning && (
          <div className="space-y-2">
            <div className="flex justify-between items-center text-sm">
              <span className="text-muted-foreground">Current Config:</span>
              <span className="font-medium">
                {getConfigName()} ({progress.currentConfig}/{progress.totalConfigs})
              </span>
            </div>
            <Progress value={configProgress} className="h-2" />
            <div className="text-xs text-muted-foreground text-right">
              Design {progress.currentDesign}/{progress.designsPerConfig}
            </div>
          </div>
        )}

        {/* Overall Progress */}
        <div className="space-y-2">
          <div className="flex justify-between items-center text-sm">
            <span className="text-muted-foreground">Overall Progress:</span>
            <span className="font-medium">{overallProgress.toFixed(1)}%</span>
          </div>
          <Progress value={overallProgress} className="h-3" />
        </div>

        {/* Statistics */}
        <div className="grid grid-cols-2 gap-4">
          <StatBox
            icon={CheckCircle}
            iconColor="text-success"
            label="Passing"
            value={progress.totalPassing}
          />
          <StatBox
            icon={AlertTriangle}
            iconColor="text-warning"
            label="Review"
            value={progress.totalReview}
          />
          <StatBox
            icon={XCircle}
            iconColor="text-destructive"
            label="Failed"
            value={progress.totalFailed}
          />
          <StatBox
            icon={TrendingUp}
            iconColor="text-info"
            label="Pass Rate"
            value={`${(progress.passRate * 100).toFixed(1)}%`}
          />
        </div>

        {/* Best Design */}
        {progress.bestDesign && (
          <div className="p-3 rounded-lg bg-muted/50 border border-primary/20">
            <div className="flex items-center gap-2 mb-2">
              <Trophy className="w-4 h-4 text-yellow-500" />
              <span className="text-sm font-medium">Best Design So Far</span>
            </div>
            <div className="text-sm font-mono text-muted-foreground">
              {progress.bestDesign.name}
            </div>
            <div className="flex gap-4 mt-2 text-xs">
              <span>
                pLDDT: <span className="font-medium text-foreground">
                  {(progress.bestDesign.plddt * 100).toFixed(1)}%
                </span>
              </span>
              <span>
                pTM: <span className="font-medium text-foreground">
                  {(progress.bestDesign.ptm * 100).toFixed(1)}%
                </span>
              </span>
              <span>
                PAE: <span className="font-medium text-foreground">
                  {progress.bestDesign.pae.toFixed(1)} A
                </span>
              </span>
            </div>
          </div>
        )}

        {/* Control Buttons */}
        {isRunning && (
          <div className="flex gap-2 pt-2">
            {onPause && (
              <Button
                variant="outline"
                size="sm"
                onClick={onPause}
                className="flex-1"
              >
                <Pause className="w-4 h-4 mr-2" />
                Pause
              </Button>
            )}
            {onCancel && (
              <Button
                variant="destructive"
                size="sm"
                onClick={onCancel}
                className="flex-1"
              >
                <X className="w-4 h-4 mr-2" />
                Cancel
              </Button>
            )}
          </div>
        )}
      </CardContent>
    </Card>
  );
}

function StatBox({
  icon: Icon,
  iconColor,
  label,
  value,
}: {
  icon: React.ComponentType<{ className?: string }>;
  iconColor: string;
  label: string;
  value: number | string;
}) {
  return (
    <div className="p-3 rounded-lg bg-muted/50">
      <div className="flex items-center gap-2 text-xs text-muted-foreground mb-1">
        <Icon className={cn("w-3 h-3", iconColor)} />
        {label}
      </div>
      <div className="text-lg font-semibold">{value}</div>
    </div>
  );
}
