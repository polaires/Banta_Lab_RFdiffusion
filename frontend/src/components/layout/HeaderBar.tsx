'use client';

import { Wifi, WifiOff, User, ChevronDown } from 'lucide-react';
import { cn } from '@/lib/utils';

interface HeaderBarProps {
  connectionStatus: 'connected' | 'disconnected' | 'connecting';
  userName?: string;
  onConnectionClick: () => void;
  onUserClick: () => void;
}

export function HeaderBar({
  connectionStatus,
  userName,
  onConnectionClick,
  onUserClick,
}: HeaderBarProps) {
  return (
    <div className="flex-1 flex items-center justify-between">
      <div />
      <div className="flex items-center gap-3">
        <button
          onClick={onConnectionClick}
          className={cn(
            'flex items-center gap-2 px-3 py-1.5 rounded-md text-sm transition-colors',
            'hover:bg-accent'
          )}
        >
          {connectionStatus === 'connected' ? (
            <>
              <Wifi className="h-4 w-4 text-success" />
              <span className="text-muted-foreground">Connected</span>
            </>
          ) : connectionStatus === 'connecting' ? (
            <>
              <Wifi className="h-4 w-4 text-warning animate-pulse" />
              <span className="text-muted-foreground">Connecting...</span>
            </>
          ) : (
            <>
              <WifiOff className="h-4 w-4 text-destructive" />
              <span className="text-muted-foreground">Disconnected</span>
            </>
          )}
        </button>

        <button
          onClick={onUserClick}
          className="flex items-center gap-2 px-3 py-1.5 rounded-md text-sm hover:bg-accent transition-colors"
        >
          <div className="h-6 w-6 rounded-full bg-primary flex items-center justify-center">
            <User className="h-3.5 w-3.5 text-primary-foreground" />
          </div>
          {userName && <span>{userName}</span>}
          <ChevronDown className="h-4 w-4 text-muted-foreground" />
        </button>
      </div>
    </div>
  );
}

export default HeaderBar;
