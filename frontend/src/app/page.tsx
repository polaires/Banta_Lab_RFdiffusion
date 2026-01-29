'use client';

import { useEffect, useState } from 'react';
import { useStore } from '@/lib/store';
import api from '@/lib/api';

// New layout components
import { MainLayout } from '@/components/layout/MainLayout';
import { Sidebar } from '@/components/layout/Sidebar';
import { ViewerPanel } from '@/components/layout/ViewerPanel';
import { ConnectionSheet } from '@/components/connection/ConnectionSheet';

// Existing panels
import { AIDesignAssistantPanel } from '@/components/AIDesignAssistantPanel';
import { TaskPanel } from '@/components/TaskPanel';
import { RFD3Panel } from '@/components/RFD3Panel';
import { RF3Panel } from '@/components/RF3Panel';
import { MPNNPanel } from '@/components/MPNNPanel';
import { DesignHistoryPanel } from '@/components/DesignHistoryPanel';

// Existing viewer - use dynamic import wrapper to prevent bundling issues
import { ProteinViewer } from '@/components/ProteinViewer';

// Map TabId to workflow step display
const WORKFLOW_STEPS = [
  { id: 'task' as const, label: 'Task Selection', status: 'pending' as const },
  { id: 'rfd3' as const, label: 'RFdiffusion', status: 'pending' as const },
  { id: 'mpnn' as const, label: 'ProteinMPNN', status: 'pending' as const },
  { id: 'rf3' as const, label: 'Validation', status: 'pending' as const },
];

export default function Home() {
  const {
    activeTab,
    setActiveTab,
    backendUrl,
    setBackendUrl,
    health,
    setHealth,
    jobs,
    manualMode,
    setManualMode,
    connectionModalOpen,
    setConnectionModalOpen,
    selectedPdb,
    // Analysis state for viewer
    focusedMetalIndex,
    focusedLigandIndex,
    metalCoordination,
    ligandData,
    pharmacophoreFeatures,
    showPharmacophores3D,
  } = useStore();

  const [isSpinning, setIsSpinning] = useState(false);

  // Check connection on mount and periodically
  useEffect(() => {
    const checkConnection = async () => {
      try {
        api.setBaseUrl(backendUrl);
        const healthResponse = await api.checkHealth();
        setHealth(healthResponse);
      } catch {
        setHealth(null);
      }
    };

    checkConnection();
    const interval = setInterval(checkConnection, 30000);
    return () => clearInterval(interval);
  }, [backendUrl, setHealth]);

  // Derive connection status from health
  const connectionStatus = health ? 'connected' : 'disconnected';

  // Build workflow steps with dynamic status based on completed jobs
  const workflowSteps = WORKFLOW_STEPS.map((step) => {
    const hasCompletedJob = jobs.some(
      (job) => job.type === step.id && job.status === 'completed'
    );
    return {
      ...step,
      status: activeTab === step.id ? 'active' as const : hasCompletedJob ? 'completed' as const : 'pending' as const,
    };
  });

  // Build history from jobs
  const history = jobs
    .filter((job) => job.status === 'completed')
    .slice(0, 10)
    .map((job) => ({
      id: job.id,
      name: `Design ${job.id.slice(0, 6)}`,
      timestamp: formatTimestamp(job.completedAt || job.createdAt),
    }));

  // Show viewer for rfd3, rf3, jobs tabs, and AI mode when a structure is loaded
  const showViewer = activeTab === 'rfd3' || activeTab === 'rf3' || activeTab === 'jobs' || (!manualMode && !!selectedPdb);

  // Handle connection
  const handleConnect = async (mode: 'runpod' | 'traditional' | 'local', url: string) => {
    try {
      api.setBaseUrl(url);
      const healthResponse = await api.checkHealth();
      if (healthResponse) {
        setBackendUrl(url);
        setHealth(healthResponse);
        return true;
      }
      return false;
    } catch {
      return false;
    }
  };

  // Handle manual mode toggle - switch view accordingly
  const handleManualModeChange = (enabled: boolean) => {
    setManualMode(enabled);
    // When toggling manual mode, switch to appropriate default view
    if (enabled && activeTab === 'ai') {
      // Switching to manual mode from AI view -> go to task selection
      setActiveTab('task');
    } else if (!enabled && activeTab !== 'jobs') {
      // Switching to AI mode (except from jobs) -> go to AI panel
      setActiveTab('ai');
    }
  };

  // Render main content based on mode and active tab
  const renderMainContent = () => {
    // AI mode (manualMode OFF) - show AI assistant for most tabs
    if (!manualMode) {
      // Jobs panel is always accessible
      if (activeTab === 'jobs') {
        return <DesignHistoryPanel />;
      }
      // Everything else shows AI assistant
      return <AIDesignAssistantPanel />;
    }

    // Manual mode (manualMode ON) - show form panels
    switch (activeTab) {
      case 'ai':
        // In manual mode, 'ai' tab shows task selection
        return <TaskPanel />;
      case 'task':
        return <TaskPanel />;
      case 'rfd3':
        return <RFD3Panel />;
      case 'mpnn':
        return <MPNNPanel />;
      case 'rf3':
        return <RF3Panel />;
      case 'jobs':
        return <DesignHistoryPanel />;
      default:
        return <TaskPanel />;
    }
  };

  return (
    <>
      <MainLayout
        sidebar={
          <Sidebar
            currentStage={activeTab === 'ai' || activeTab === 'jobs' ? 'task' : activeTab}
            workflowSteps={workflowSteps}
            history={history}
            manualMode={manualMode}
            onManualModeChange={handleManualModeChange}
            onStageClick={(stage) => setActiveTab(stage)}
            onNewDesign={() => {
              setManualMode(false); // Switch to AI mode
              setActiveTab('ai');
              // Could also clear state here
            }}
            onHistoryClick={(id) => {
              // Navigate to jobs and select the job
              setActiveTab('jobs');
            }}
            onSettingsClick={() => setConnectionModalOpen(true)}
            // Connection props (moved from header)
            connectionStatus={connectionStatus}
            onConnectionClick={() => setConnectionModalOpen(true)}
          />
        }
        main={
          <div className="h-full p-6 overflow-auto">
            <div className="bg-card rounded-xl border border-border p-6">
              {renderMainContent()}
            </div>
          </div>
        }
        viewer={showViewer ? (
          <ViewerPanel
            isSpinning={isSpinning}
            onToggleSpin={() => setIsSpinning(!isSpinning)}
            onReset={() => {/* TODO: reset viewer */}}
            onExpand={() => {/* TODO: expand viewer */}}
          >
            {selectedPdb ? (
              <ProteinViewer
                pdbContent={selectedPdb}
                className="h-full w-full"
                focusedMetalIndex={focusedMetalIndex}
                focusedLigandIndex={focusedLigandIndex}
                metalCoordination={metalCoordination}
                ligandData={ligandData}
                pharmacophoreFeatures={pharmacophoreFeatures ?? undefined}
                showPharmacophores={showPharmacophores3D}
              />
            ) : (
              <div className="h-full flex items-center justify-center text-muted-foreground text-sm">
                No structure loaded
              </div>
            )}
          </ViewerPanel>
        ) : undefined}
      />

      {/* Connection Sheet */}
      <ConnectionSheet
        open={connectionModalOpen}
        onOpenChange={setConnectionModalOpen}
        onConnect={handleConnect}
        defaultUrl={backendUrl}
      />
    </>
  );
}

// Helper to format timestamps
function formatTimestamp(dateStr: string): string {
  const date = new Date(dateStr);
  const now = new Date();
  const diffMs = now.getTime() - date.getTime();
  const diffMins = Math.floor(diffMs / 60000);
  const diffHours = Math.floor(diffMins / 60);
  const diffDays = Math.floor(diffHours / 24);

  if (diffMins < 1) return 'just now';
  if (diffMins < 60) return `${diffMins}m ago`;
  if (diffHours < 24) return `${diffHours}h ago`;
  return `${diffDays}d ago`;
}
