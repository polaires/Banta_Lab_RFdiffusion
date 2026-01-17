'use client';

import { useState } from 'react';
import { Cloud, Server, Wrench, Loader2, CheckCircle, XCircle } from 'lucide-react';
import {
  Sheet,
  SheetContent,
  SheetDescription,
  SheetHeader,
  SheetTitle,
} from '@/components/ui/sheet';
import { Button } from '@/components/ui/button';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Alert, AlertDescription } from '@/components/ui/alert';
import { cn } from '@/lib/utils';

type ConnectionMode = 'runpod' | 'traditional' | 'local';
type ConnectionStatus = 'idle' | 'connecting' | 'connected' | 'failed';

interface ConnectionSheetProps {
  open: boolean;
  onOpenChange: (open: boolean) => void;
  onConnect: (mode: ConnectionMode, url: string) => Promise<boolean>;
  defaultUrl?: string;
  defaultMode?: ConnectionMode;
}

const modes = [
  { id: 'runpod' as const, label: 'RunPod Serverless', description: 'Recommended for most users', icon: Cloud },
  { id: 'traditional' as const, label: 'Traditional Server', description: 'Self-hosted backend', icon: Server },
  { id: 'local' as const, label: 'Local Development', description: 'For contributors', icon: Wrench },
];

export function ConnectionSheet({
  open,
  onOpenChange,
  onConnect,
  defaultUrl = '',
  defaultMode = 'runpod',
}: ConnectionSheetProps) {
  const [mode, setMode] = useState<ConnectionMode>(defaultMode);
  const [url, setUrl] = useState(defaultUrl);
  const [status, setStatus] = useState<ConnectionStatus>('idle');
  const [error, setError] = useState<string | null>(null);

  const handleConnect = async () => {
    setStatus('connecting');
    setError(null);
    try {
      const success = await onConnect(mode, url);
      if (success) {
        setStatus('connected');
        setTimeout(() => onOpenChange(false), 1000);
      } else {
        setStatus('failed');
        setError('Connection failed. Please check the URL and try again.');
      }
    } catch (err) {
      setStatus('failed');
      setError(err instanceof Error ? err.message : 'Connection failed');
    }
  };

  return (
    <Sheet open={open} onOpenChange={onOpenChange}>
      <SheetContent className="w-[400px] sm:w-[540px]">
        <SheetHeader>
          <SheetTitle>Backend Connection</SheetTitle>
          <SheetDescription>Configure how to connect to the RFdiffusion backend</SheetDescription>
        </SheetHeader>

        <div className="mt-6 space-y-6">
          <div className="space-y-3">
            <Label>Choose your setup</Label>
            <div className="space-y-2">
              {modes.map((m) => (
                <button
                  key={m.id}
                  onClick={() => setMode(m.id)}
                  className={cn(
                    'w-full flex items-center gap-3 p-3 rounded-lg border transition-colors text-left',
                    mode === m.id ? 'border-primary bg-accent' : 'border-border hover:bg-accent/50'
                  )}
                >
                  <m.icon className="h-5 w-5 text-muted-foreground" />
                  <div>
                    <div className="font-medium text-sm">{m.label}</div>
                    <div className="text-xs text-muted-foreground">{m.description}</div>
                  </div>
                </button>
              ))}
            </div>
          </div>

          <div className="space-y-2">
            <Label htmlFor="server-url">Server URL</Label>
            <Input
              id="server-url"
              placeholder={mode === 'runpod' ? 'https://api.runpod.ai/v2/...' : mode === 'local' ? 'http://localhost:8000' : 'https://your-server.com/api'}
              value={url}
              onChange={(e) => setUrl(e.target.value)}
            />
          </div>

          {error && (
            <Alert variant="destructive">
              <XCircle className="h-4 w-4" />
              <AlertDescription>{error}</AlertDescription>
            </Alert>
          )}

          {status === 'connected' && (
            <Alert>
              <CheckCircle className="h-4 w-4 text-green-600" />
              <AlertDescription>Connected successfully!</AlertDescription>
            </Alert>
          )}

          <Button onClick={handleConnect} disabled={!url || status === 'connecting'} className="w-full">
            {status === 'connecting' ? (
              <><Loader2 className="h-4 w-4 mr-2 animate-spin" />Connecting...</>
            ) : 'Test Connection'}
          </Button>

          <div className="flex items-center gap-2 text-sm">
            <span className="text-muted-foreground">Status:</span>
            {status === 'connected' ? (
              <span className="flex items-center gap-1 text-green-600"><span className="w-2 h-2 rounded-full bg-green-600" />Connected</span>
            ) : status === 'connecting' ? (
              <span className="flex items-center gap-1 text-amber-600"><span className="w-2 h-2 rounded-full bg-amber-600 animate-pulse" />Connecting</span>
            ) : status === 'failed' ? (
              <span className="flex items-center gap-1 text-destructive"><span className="w-2 h-2 rounded-full bg-destructive" />Failed</span>
            ) : (
              <span className="flex items-center gap-1 text-muted-foreground"><span className="w-2 h-2 rounded-full bg-muted-foreground" />Not connected</span>
            )}
          </div>
        </div>
      </SheetContent>
    </Sheet>
  );
}

export default ConnectionSheet;
