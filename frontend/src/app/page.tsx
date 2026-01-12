'use client';

import { useEffect } from 'react';
import { useStore } from '@/lib/store';
import api from '@/lib/api';
import { Header } from '@/components/Header';
import { ConnectionModal } from '@/components/ConnectionModal';
import { CollapsibleViewer } from '@/components/CollapsibleViewer';
import { AIDesignAssistantPanel } from '@/components/AIDesignAssistantPanel';
import { TaskPanel } from '@/components/TaskPanel';
import { RFD3Panel } from '@/components/RFD3Panel';
import { RF3Panel } from '@/components/RF3Panel';
import { MPNNPanel } from '@/components/MPNNPanel';
import { JobsPanel } from '@/components/JobsPanel';
import { NotificationToast } from '@/components/NotificationToast';

export default function Home() {
  const { activeTab, backendUrl, setHealth } = useStore();

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

  // Show viewer for rfd3 and rf3 phases, hide for task selection and mpnn
  const showViewer = activeTab === 'rfd3' || activeTab === 'rf3' || activeTab === 'jobs';

  return (
    <div className="min-h-screen bg-slate-50 text-slate-900 flex flex-col">
      {/* Header with workflow stepper */}
      <Header />

      {/* Main Content */}
      <main className="flex-1 max-w-6xl mx-auto w-full px-6 py-8">
        {/* Active Panel */}
        <div className="bg-white rounded-2xl shadow-card p-8">
          {activeTab === 'ai' && <AIDesignAssistantPanel />}
          {activeTab === 'task' && <TaskPanel />}
          {activeTab === 'rfd3' && <RFD3Panel />}
          {activeTab === 'mpnn' && <MPNNPanel />}
          {activeTab === 'rf3' && <RF3Panel />}
          {activeTab === 'jobs' && <JobsPanel />}
        </div>

        {/* Collapsible Structure Viewer - shown for RFD3, RF3, and Jobs */}
        {showViewer && <CollapsibleViewer />}
      </main>

      {/* Connection Modal */}
      <ConnectionModal />

      {/* Notifications */}
      <NotificationToast />
    </div>
  );
}
