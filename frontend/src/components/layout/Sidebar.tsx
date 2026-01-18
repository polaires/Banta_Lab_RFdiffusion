'use client';

import {
  Plus,
  Settings,
  Command,
  Check,
  Circle,
  RefreshCw,
  FlaskConical,
  Wifi,
  WifiOff,
  User,
  ChevronDown,
} from 'lucide-react';
import { cn } from '@/lib/utils';
import { Button } from '@/components/ui/button';
import { Switch } from '@/components/ui/switch';
import { ScrollArea } from '@/components/ui/scroll-area';
import { Separator } from '@/components/ui/separator';
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from '@/components/ui/tooltip';

type WorkflowStage = 'task' | 'rfd3' | 'mpnn' | 'rf3';

interface WorkflowStep {
  id: WorkflowStage;
  label: string;
  status: 'pending' | 'active' | 'completed';
  iterations?: number;
}

interface DesignHistoryItem {
  id: string;
  name: string;
  timestamp: string;
}

interface SidebarProps {
  currentStage: WorkflowStage;
  workflowSteps: WorkflowStep[];
  history: DesignHistoryItem[];
  manualMode: boolean;
  onManualModeChange: (enabled: boolean) => void;
  onStageClick: (stage: WorkflowStage) => void;
  onNewDesign: () => void;
  onHistoryClick: (id: string) => void;
  onSettingsClick: () => void;
  // Connection & user props (previously in header)
  connectionStatus?: 'connected' | 'disconnected' | 'connecting';
  onConnectionClick?: () => void;
  userName?: string;
  onUserClick?: () => void;
}

interface WorkflowStepItemProps {
  step: WorkflowStep;
  isActive: boolean;
  manualMode: boolean;
  onClick: () => void;
}

export function Sidebar({
  currentStage,
  workflowSteps,
  history,
  manualMode,
  onManualModeChange,
  onStageClick,
  onNewDesign,
  onHistoryClick,
  onSettingsClick,
  connectionStatus = 'disconnected',
  onConnectionClick,
  userName,
  onUserClick,
}: SidebarProps) {
  return (
    <TooltipProvider>
      <div className="h-full flex flex-col">
        {/* Logo + User */}
        <div className="p-3 border-b border-sidebar-border">
          <div className="flex items-center justify-between">
            <div className="flex items-center gap-2">
              <FlaskConical className="h-5 w-5 text-primary" />
              <div>
                <div className="font-semibold text-sm leading-tight">Banta Lab</div>
                <div className="text-[10px] text-muted-foreground">RFdiffusion</div>
              </div>
            </div>
            {/* User menu */}
            <button
              onClick={onUserClick}
              className="p-1 rounded-md hover:bg-sidebar-accent transition-colors"
            >
              <div className="h-7 w-7 rounded-full bg-primary flex items-center justify-center">
                <User className="h-3.5 w-3.5 text-primary-foreground" />
              </div>
            </button>
          </div>
        </div>

        {/* New Design Button */}
        <div className="p-3">
          <Button onClick={onNewDesign} className="w-full" size="sm">
            <Plus className="h-4 w-4 mr-2" />
            New Design
          </Button>
        </div>

        <Separator />

        {/* Workflow Steps */}
        <div className="p-3">
          <div className="text-xs font-medium text-muted-foreground uppercase tracking-wider mb-2">
            Current Workflow
          </div>
          <div className="space-y-1">
            {workflowSteps.map((step) => (
              <WorkflowStepItem
                key={step.id}
                step={step}
                isActive={currentStage === step.id}
                manualMode={manualMode}
                onClick={() => onStageClick(step.id)}
              />
            ))}
          </div>
        </div>

        <Separator />

        {/* History */}
        <div className="flex-1 p-3 overflow-hidden">
          <div className="text-xs font-medium text-muted-foreground uppercase tracking-wider mb-2">
            History
          </div>
          <ScrollArea className="h-[calc(100%-24px)]">
            <div className="space-y-1">
              {history.map((item) => (
                <button
                  key={item.id}
                  onClick={() => onHistoryClick(item.id)}
                  className="w-full flex items-center gap-2 px-2 py-1.5 rounded text-sm hover:bg-sidebar-accent text-left"
                >
                  <Check className="h-3.5 w-3.5 text-muted-foreground" />
                  <span className="flex-1 truncate">{item.name}</span>
                  <span className="text-xs text-muted-foreground">{item.timestamp}</span>
                </button>
              ))}
            </div>
          </ScrollArea>
        </div>

        <Separator />

        {/* Bottom Actions */}
        <div className="p-3 space-y-1">
          {/* Connection Status */}
          <button
            onClick={onConnectionClick}
            className={cn(
              'w-full flex items-center gap-2 px-2 py-1.5 rounded text-sm transition-colors',
              'hover:bg-sidebar-accent'
            )}
          >
            {connectionStatus === 'connected' ? (
              <Wifi className="h-4 w-4 text-green-600" />
            ) : connectionStatus === 'connecting' ? (
              <Wifi className="h-4 w-4 text-amber-500 animate-pulse" />
            ) : (
              <WifiOff className="h-4 w-4 text-destructive" />
            )}
            <span className={cn(
              'flex-1 text-left',
              connectionStatus === 'connected' ? 'text-foreground' : 'text-muted-foreground'
            )}>
              {connectionStatus === 'connected' ? 'Connected' :
               connectionStatus === 'connecting' ? 'Connecting...' : 'Disconnected'}
            </span>
          </button>

          <button
            onClick={onSettingsClick}
            className="w-full flex items-center gap-2 px-2 py-1.5 rounded text-sm hover:bg-sidebar-accent"
          >
            <Settings className="h-4 w-4" />
            Settings
          </button>

          <Separator className="my-2" />

          <div className="flex items-center gap-2 px-2 py-1.5 text-sm text-muted-foreground">
            <Command className="h-3.5 w-3.5" />
            <span className="text-xs">K</span>
            <span className="flex-1">AI Assistant</span>
          </div>

          <div className="flex items-center justify-between px-2 py-1.5">
            <span className="text-sm">Manual Mode</span>
            <Switch
              checked={manualMode}
              onCheckedChange={onManualModeChange}
            />
          </div>
        </div>
      </div>
    </TooltipProvider>
  );
}

function WorkflowStepItem({
  step,
  isActive,
  manualMode,
  onClick,
}: WorkflowStepItemProps) {
  // In manual mode, all steps are clickable. Otherwise, only completed/active steps are clickable.
  const isClickable = manualMode || step.status !== 'pending';

  const content = (
    <button
      onClick={isClickable ? onClick : undefined}
      disabled={!isClickable}
      className={cn(
        'w-full flex items-center gap-2 px-2 py-1.5 rounded text-sm transition-colors',
        isActive && 'bg-sidebar-accent',
        isClickable ? 'hover:bg-sidebar-accent cursor-pointer' : 'cursor-not-allowed opacity-60'
      )}
    >
      {step.status === 'completed' ? (
        <Check className="h-3.5 w-3.5 text-primary" />
      ) : step.status === 'active' ? (
        <Circle className="h-3.5 w-3.5 fill-primary text-primary" />
      ) : (
        <Circle className="h-3.5 w-3.5 text-muted-foreground" />
      )}
      <span className="flex-1 text-left">{step.label}</span>
      {step.iterations && step.iterations > 1 && (
        <span className="flex items-center gap-0.5 text-xs text-muted-foreground">
          <RefreshCw className="h-3 w-3" />
          {step.iterations}
        </span>
      )}
    </button>
  );

  if (!isClickable) {
    return (
      <Tooltip>
        <TooltipTrigger asChild>{content}</TooltipTrigger>
        <TooltipContent side="right">
          <p>Complete previous stage first, or enable Manual Mode</p>
        </TooltipContent>
      </Tooltip>
    );
  }

  return content;
}

export default Sidebar;
