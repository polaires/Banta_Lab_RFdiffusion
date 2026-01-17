'use client';

import { RotateCcw, Maximize2, Play, Pause } from 'lucide-react';
import { Button } from '@/components/ui/button';
import { Separator } from '@/components/ui/separator';
import { cn } from '@/lib/utils';

interface StructureInfo {
  residues?: number;
  plddt?: number;
  rmsd?: number;
}

interface ViewerPanelProps {
  children: React.ReactNode;
  structureInfo?: StructureInfo;
  isSpinning?: boolean;
  onToggleSpin?: () => void;
  onReset?: () => void;
  onExpand?: () => void;
  showControls?: boolean;
}

export function ViewerPanel({
  children,
  structureInfo,
  isSpinning = false,
  onToggleSpin,
  onReset,
  onExpand,
  showControls = true,
}: ViewerPanelProps) {
  return (
    <div className="h-full flex flex-col">
      <div className="flex-1 relative bg-muted/30">
        {children}
      </div>

      {showControls && (
        <>
          <Separator />
          <div className="p-2 flex items-center gap-1">
            <Button variant="ghost" size="sm" onClick={onToggleSpin} className="h-8 px-2">
              {isSpinning ? <Pause className="h-4 w-4" /> : <Play className="h-4 w-4" />}
              <span className="ml-1 text-xs">Spin</span>
            </Button>
            <Button variant="ghost" size="sm" onClick={onReset} className="h-8 px-2">
              <RotateCcw className="h-4 w-4" />
              <span className="ml-1 text-xs">Reset</span>
            </Button>
            <div className="flex-1" />
            <Button variant="ghost" size="sm" onClick={onExpand} className="h-8 px-2">
              <Maximize2 className="h-4 w-4" />
            </Button>
          </div>
        </>
      )}

      {structureInfo && (
        <>
          <Separator />
          <div className="p-3 space-y-2">
            <div className="text-xs font-medium text-muted-foreground uppercase tracking-wider">
              Structure Info
            </div>
            <div className="grid grid-cols-3 gap-2 text-sm">
              {structureInfo.residues !== undefined && (
                <div>
                  <div className="text-muted-foreground text-xs">Residues</div>
                  <div className="font-mono">{structureInfo.residues}</div>
                </div>
              )}
              {structureInfo.plddt !== undefined && (
                <div>
                  <div className="text-muted-foreground text-xs">pLDDT</div>
                  <div className={cn(
                    'font-mono',
                    structureInfo.plddt >= 90 ? 'text-green-600' :
                    structureInfo.plddt >= 70 ? 'text-amber-600' : 'text-red-600'
                  )}>
                    {structureInfo.plddt.toFixed(1)}
                  </div>
                </div>
              )}
              {structureInfo.rmsd !== undefined && (
                <div>
                  <div className="text-muted-foreground text-xs">RMSD</div>
                  <div className="font-mono">{structureInfo.rmsd.toFixed(2)} Å</div>
                </div>
              )}
            </div>
          </div>
        </>
      )}

      <Separator />
      <div className="p-3">
        <div className="text-xs font-medium text-muted-foreground uppercase tracking-wider mb-2">
          Confidence Coloring
        </div>
        <div className="grid grid-cols-2 gap-1 text-xs">
          <div className="flex items-center gap-1.5">
            <div className="w-3 h-3 rounded bg-blue-600" />
            <span>Very High (&gt;90)</span>
          </div>
          <div className="flex items-center gap-1.5">
            <div className="w-3 h-3 rounded bg-cyan-500" />
            <span>High (70-90)</span>
          </div>
          <div className="flex items-center gap-1.5">
            <div className="w-3 h-3 rounded bg-yellow-500" />
            <span>Medium (50-70)</span>
          </div>
          <div className="flex items-center gap-1.5">
            <div className="w-3 h-3 rounded bg-orange-500" />
            <span>Low (&lt;50)</span>
          </div>
        </div>
      </div>
    </div>
  );
}

export default ViewerPanel;
