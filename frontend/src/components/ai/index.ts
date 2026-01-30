export { InterviewMode } from './InterviewMode';
export { PreferenceSummaryCard } from './PreferenceSummaryCard';
export { EvaluationCard } from './EvaluationCard';
export { JobProgressCard } from './JobProgressCard';

// Metal scaffold design (Round 7b monomer workflow)
export { MetalScaffoldInterviewMode } from './MetalScaffoldInterviewMode';
export type { MetalScaffoldPreferences } from './MetalScaffoldInterviewMode';
export { MetalScaffoldDesignGallery } from './MetalScaffoldDesignGallery';

// Metal scaffold design result type
export interface MetalScaffoldDesignResult {
  id: string;
  name: string;
  configName: string;
  sequence: string;
  tier: 'S' | 'A' | 'B' | 'C' | 'F';
  plddt: number;
  ptm: number;
  pae: number;
  status: 'passed' | 'failed';
  pdbContent?: string;
}

// Ligand interface design components
export { LigandInterviewMode } from './LigandInterviewMode';
export type { LigandPreferences } from './LigandInterviewMode';
export { LigandAnalysisCard } from './LigandAnalysisCard';
export type { LigandAnalysis } from './LigandAnalysisCard';
export { LigandEvaluationCard } from './LigandEvaluationCard';
export type { LigandEvaluation } from './LigandEvaluationCard';

// Protein binder design components
export { BinderInterviewMode } from './BinderInterviewMode';
export type { BinderPreferences } from './BinderInterviewMode';
export { BinderEvaluationCard } from './BinderEvaluationCard';
export type { BinderEvaluation } from './BinderEvaluationCard';

// Design gallery for viewing multiple designs
export { DesignGallery } from './DesignGallery';
export type { DesignResult } from './DesignGallery';

// AI Design Pipeline workflow visualization
export { AIDesignPipelineWorkflow, AIDesignPipelineStatus } from './AIDesignPipelineWorkflow';
export type { AIDesignStage } from './AIDesignPipelineWorkflow';

// AIDesignPanel removed — replaced by AIDesignAssistantPanel + PipelineRunner
