'use client';

import { useState, useCallback } from 'react';
import { AIDesignStage } from '@/components/ai/AIDesignPipelineWorkflow';

interface AIDesignResult {
  success: boolean;
  error?: string;
  design_intent?: {
    metal_type?: string;
    ligand_name?: string;
    design_goal: string;
    target_topology: string;
    chain_length_range: string;
    confidence: number;
    warnings: string[];
    suggestions: string[];
  };
  num_backbones: number;
  num_sequences: number;
  pass_rate: number;
  best_rmsd?: number;
  best_plddt?: number;
  best_sequence?: {
    sequence: string;
    confidence: number;
  };
  best_sequence_pdb?: string;
  backbone_pdbs?: string[];
  sequences?: Array<{ sequence: string; confidence: number }>;
  analysis?: Record<string, any>;
  report: string;
  recommendations: string[];
  timings: Record<string, number>;
  total_time: number;
}

interface UseAIDesignOptions {
  apiUrl?: string;
  claudeApiKey?: string;
  onStageChange?: (stage: AIDesignStage) => void;
  onProgress?: (message: string) => void;
}

interface AIDesignState {
  stage: AIDesignStage;
  stageInfo: {
    metal?: string;
    ligand?: string;
    confidence?: number;
    numDesigns?: number;
    numSequences?: number;
    passRate?: number;
    bestRmsd?: number;
    bestPlddt?: number;
    message?: string;
  };
  result?: AIDesignResult;
  error?: string;
  isRunning: boolean;
}

export function useAIDesign(options: UseAIDesignOptions = {}) {
  const {
    apiUrl = '/api/runpod/runsync',
    claudeApiKey,
    onStageChange,
    onProgress
  } = options;

  const [state, setState] = useState<AIDesignState>({
    stage: 'idle',
    stageInfo: {},
    isRunning: false
  });

  const updateStage = useCallback((stage: AIDesignStage, message?: string, extraInfo?: Partial<AIDesignState['stageInfo']>) => {
    setState(prev => ({
      ...prev,
      stage,
      stageInfo: {
        ...prev.stageInfo,
        message,
        ...extraInfo
      }
    }));
    onStageChange?.(stage);
    if (message) onProgress?.(message);
  }, [onStageChange, onProgress]);

  const runAIDesign = useCallback(async (
    query: string,
    options?: {
      numDesigns?: number;
      numSequences?: number;
      validate?: boolean;
      sessionName?: string;
    }
  ): Promise<AIDesignResult | null> => {
    if (state.isRunning) return null;

    setState(prev => ({
      ...prev,
      isRunning: true,
      error: undefined,
      result: undefined,
      stageInfo: {}
    }));

    try {
      // Stage 1: Parsing
      updateStage('parsing', 'Understanding your design request...');

      // Make API request
      const response = await fetch(apiUrl, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          input: {
            task: 'ai_design',
            query,
            claude_api_key: claudeApiKey,
            num_designs: options?.numDesigns ?? 4,
            num_sequences: options?.numSequences ?? 8,
            validate: options?.validate ?? true,
            session_name: options?.sessionName
          }
        })
      });

      if (!response.ok) {
        throw new Error(`API error: ${response.status}`);
      }

      const data = await response.json();

      // Check for errors
      if (data.status === 'failed' || data.error) {
        throw new Error(data.error || data.output?.error || 'Design failed');
      }

      // Extract result
      const result = data.output?.result ?? data.result;
      if (!result) {
        throw new Error('No result in response');
      }

      // Update state with result
      const finalResult: AIDesignResult = {
        success: result.success ?? true,
        error: result.error,
        design_intent: result.design_intent,
        num_backbones: result.num_backbones ?? 0,
        num_sequences: result.num_sequences ?? 0,
        pass_rate: result.pass_rate ?? 0,
        best_rmsd: result.best_rmsd,
        best_plddt: result.best_plddt,
        best_sequence: result.best_sequence,
        best_sequence_pdb: result.best_sequence_pdb,
        backbone_pdbs: result.backbone_pdbs,
        sequences: result.sequences,
        analysis: result.analysis,
        report: result.report ?? '',
        recommendations: result.recommendations ?? [],
        timings: result.timings ?? {},
        total_time: result.total_time ?? 0
      };

      // Update final state
      setState(prev => ({
        ...prev,
        stage: finalResult.success ? 'complete' : 'error',
        stageInfo: {
          metal: finalResult.design_intent?.metal_type,
          ligand: finalResult.design_intent?.ligand_name,
          confidence: finalResult.design_intent?.confidence,
          numDesigns: finalResult.num_backbones,
          numSequences: finalResult.num_sequences,
          passRate: finalResult.pass_rate,
          bestRmsd: finalResult.best_rmsd,
          bestPlddt: finalResult.best_plddt
        },
        result: finalResult,
        error: finalResult.error,
        isRunning: false
      }));

      return finalResult;

    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : 'Unknown error';

      setState(prev => ({
        ...prev,
        stage: 'error',
        stageInfo: {
          ...prev.stageInfo,
          message: errorMessage
        },
        error: errorMessage,
        isRunning: false
      }));

      return null;
    }
  }, [state.isRunning, apiUrl, claudeApiKey, updateStage]);

  const reset = useCallback(() => {
    setState({
      stage: 'idle',
      stageInfo: {},
      isRunning: false
    });
  }, []);

  return {
    ...state,
    runAIDesign,
    reset
  };
}

export default useAIDesign;
