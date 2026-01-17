'use client';

import { useEffect, useState } from 'react';
import { useStore } from '@/lib/store';
import api from '@/lib/api';

// New layout components
import { MainLayout } from '@/components/layout/MainLayout';
import { Sidebar } from '@/components/layout/Sidebar';
import { HeaderBar } from '@/components/layout/HeaderBar';
import { ViewerPanel } from '@/components/layout/ViewerPanel';
import { ConnectionSheet } from '@/components/connection/ConnectionSheet';

// Existing panels
import { AIDesignAssistantPanel } from '@/components/AIDesignAssistantPanel';
import { TaskPanel } from '@/components/TaskPanel';
import { RFD3Panel } from '@/components/RFD3Panel';
import { RF3Panel } from '@/components/RF3Panel';
import { MPNNPanel } from '@/components/MPNNPanel';
import { DesignHistoryPanel } from '@/components/DesignHistoryPanel';

// Existing viewer
import { ProteinViewerClient } from '@/components/ProteinViewerClient';

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
    latestConfidences,
    latestRmsdResult,
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

  // Show viewer for rfd3, rf3, and jobs tabs
  const showViewer = activeTab === 'rfd3' || activeTab === 'rf3' || activeTab === 'jobs';

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

  // Render main content based on mode and active tab
  const renderMainContent = () => {
    if (!manualMode && activeTab === 'ai') {
      return <AIDesignAssistantPanel />;
    }

    // Manual mode or non-AI tabs
    switch (activeTab) {
      case 'ai':
        return <AIDesignAssistantPanel />;
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
        return <AIDesignAssistantPanel />;
    }
  };

  return (
    <>
      <MainLayout
        header={
          <HeaderBar
            connectionStatus={connectionStatus}
            onConnectionClick={() => setConnectionModalOpen(true)}
            onUserClick={() => {/* TODO: user menu */}}
          />
        }
        sidebar={
          <Sidebar
            currentStage={activeTab === 'ai' || activeTab === 'jobs' ? 'task' : activeTab}
            workflowSteps={workflowSteps}
            history={history}
            manualMode={manualMode}
            onManualModeChange={setManualMode}
            onStageClick={(stage) => setActiveTab(stage)}
            onNewDesign={() => {
              setActiveTab('ai');
              // Could also clear state here
            }}
            onHistoryClick={(id) => {
              // Navigate to jobs and select the job
              setActiveTab('jobs');
            }}
            onSettingsClick={() => setConnectionModalOpen(true)}
          />
        }
        main={
          <div className="h-full p-6">
            <div className="bg-card rounded-xl border border-border p-6 h-full overflow-auto">
              {renderMainContent()}
            </div>
          </div>
        }
        viewer={showViewer ? (
          <ViewerPanel
            structureInfo={selectedPdb ? {
              residues: 150, // TODO: extract from PDB
              plddt: latestConfidences?.summary_confidences?.overall_plddt,
              rmsd: latestRmsdResult?.rmsd,
            } : undefined}
            isSpinning={isSpinning}
            onToggleSpin={() => setIsSpinning(!isSpinning)}
            onReset={() => {/* TODO: reset viewer */}}
            onExpand={() => {/* TODO: expand viewer */}}
          >
            {selectedPdb ? (
              <ProteinViewerClient pdbContent={selectedPdb} />
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
