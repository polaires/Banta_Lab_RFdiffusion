'use client';

import { useStore } from '@/lib/store';
import { ConnectionStatus } from '@/components/ConnectionStatus';
import { ProteinViewer } from '@/components/ProteinViewer';
import { RFD3Panel } from '@/components/RFD3Panel';
import { RF3Panel } from '@/components/RF3Panel';
import { MPNNPanel } from '@/components/MPNNPanel';
import { JobsPanel } from '@/components/JobsPanel';
import { Dna, Atom, FlaskConical, History } from 'lucide-react';

const tabs = [
  { id: 'rfd3' as const, label: 'RFD3 Design', icon: Dna, color: 'text-blue-400' },
  { id: 'rf3' as const, label: 'RF3 Predict', icon: Atom, color: 'text-green-400' },
  { id: 'mpnn' as const, label: 'ProteinMPNN', icon: FlaskConical, color: 'text-purple-400' },
  { id: 'jobs' as const, label: 'Jobs', icon: History, color: 'text-gray-400' },
];

export default function Home() {
  const { activeTab, setActiveTab, selectedPdb } = useStore();

  return (
    <div className="min-h-screen bg-gray-900 text-white">
      {/* Header */}
      <header className="border-b border-gray-800 bg-gray-950">
        <div className="max-w-7xl mx-auto px-4 py-4 flex items-center justify-between">
          <div className="flex items-center gap-3">
            <Dna className="w-8 h-8 text-blue-400" />
            <div>
              <h1 className="text-xl font-bold">Foundry Protein Design</h1>
              <p className="text-xs text-gray-500">RFdiffusion3 | RosettaFold3 | ProteinMPNN</p>
            </div>
          </div>
          <div className="text-sm text-gray-400">
            Banta Lab
          </div>
        </div>
      </header>

      <div className="max-w-7xl mx-auto px-4 py-6">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {/* Left Sidebar - Connection & Settings */}
          <div className="lg:col-span-1 space-y-6">
            <ConnectionStatus />

            {/* Tab Navigation */}
            <div className="bg-gray-800 rounded-lg p-2">
              <nav className="space-y-1">
                {tabs.map((tab) => (
                  <button
                    key={tab.id}
                    onClick={() => setActiveTab(tab.id)}
                    className={`w-full flex items-center gap-3 px-4 py-3 rounded-lg transition ${
                      activeTab === tab.id
                        ? 'bg-gray-700 text-white'
                        : 'text-gray-400 hover:bg-gray-700/50 hover:text-gray-200'
                    }`}
                  >
                    <tab.icon className={`w-5 h-5 ${activeTab === tab.id ? tab.color : ''}`} />
                    <span className="font-medium">{tab.label}</span>
                  </button>
                ))}
              </nav>
            </div>
          </div>

          {/* Main Content Area */}
          <div className="lg:col-span-2 space-y-6">
            {/* Tool Panel */}
            <div className="bg-gray-800 rounded-lg p-6">
              {activeTab === 'rfd3' && <RFD3Panel />}
              {activeTab === 'rf3' && <RF3Panel />}
              {activeTab === 'mpnn' && <MPNNPanel />}
              {activeTab === 'jobs' && <JobsPanel />}
            </div>

            {/* Structure Viewer */}
            <div className="bg-gray-800 rounded-lg p-4">
              <h3 className="text-lg font-semibold mb-4">Structure Viewer</h3>
              <ProteinViewer pdbContent={selectedPdb} className="h-96" />
            </div>
          </div>
        </div>
      </div>

      {/* Footer */}
      <footer className="border-t border-gray-800 mt-12 py-6">
        <div className="max-w-7xl mx-auto px-4 text-center text-sm text-gray-500">
          <p>Powered by RosettaCommons Foundry on RunPod GPU</p>
          <p className="mt-1">
            <a href="https://github.com/RosettaCommons/foundry" className="hover:text-gray-300 transition">
              Foundry GitHub
            </a>
            {' | '}
            <a href="https://www.runpod.io" className="hover:text-gray-300 transition">
              RunPod
            </a>
          </p>
        </div>
      </footer>
    </div>
  );
}
