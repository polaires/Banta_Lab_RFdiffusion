'use client';

import { useState, useEffect } from 'react';
import { Sparkles, ChevronUp, ChevronDown, Send, Loader2, FlaskConical } from 'lucide-react';
import type {
  TaskType,
  FormFieldConfig,
  AIAnalyzeResponse,
} from '@/lib/ai-types';
import { METAL_OPTIONS, SYMMETRY_OPTIONS, getConfidenceColor, formatConfidence } from '@/lib/ai-types';

interface SmartFormProps {
  taskType: TaskType | null;
  formConfig: FormFieldConfig[];
  suggestedParams: Record<string, unknown>;
  reasoning: string;
  confidence: number;
  onSubmit: (params: Record<string, unknown>) => void;
  onChat: (message: string) => void;
  isLoading?: boolean;
}

/**
 * Dynamic form that adapts based on AI-detected task type.
 * Shows only relevant fields with AI-suggested values.
 */
export function SmartForm({
  taskType,
  formConfig,
  suggestedParams,
  reasoning,
  confidence,
  onSubmit,
  onChat,
  isLoading = false,
}: SmartFormProps) {
  const [values, setValues] = useState<Record<string, unknown>>(suggestedParams);
  const [chatInput, setChatInput] = useState('');
  const [expandedReasoning, setExpandedReasoning] = useState(false);

  // Update values when suggested params change
  useEffect(() => {
    setValues(suggestedParams);
  }, [suggestedParams]);

  const handleChange = (fieldId: string, value: unknown) => {
    setValues(prev => ({ ...prev, [fieldId]: value }));
  };

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    onSubmit(values);
  };

  const handleChatSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    if (chatInput.trim()) {
      onChat(chatInput.trim());
      setChatInput('');
    }
  };

  const renderField = (field: FormFieldConfig) => {
    const value = values[field.id] ?? field.suggested_value ?? field.default_value;
    const hasAiSuggestion = field.ai_reasoning || field.suggested_value !== undefined;

    return (
      <div key={field.id} className="mb-4">
        <label className="block text-sm font-medium text-foreground mb-1.5">
          {field.label}
          {field.required && <span className="text-red-500 ml-1">*</span>}
          {hasAiSuggestion && (
            <span className="ml-2 text-xs text-primary bg-primary/10 px-2 py-0.5 rounded-full">
              AI suggested
            </span>
          )}
        </label>

        {/* Field input based on type */}
        {field.type === 'text' && (
          <input
            type="text"
            value={String(value || '')}
            onChange={(e) => handleChange(field.id, e.target.value)}
            className="w-full px-3 py-2 rounded-lg border border-border focus:border-primary focus:ring-1 focus:ring-primary text-sm"
            required={field.required}
          />
        )}

        {field.type === 'number' && (
          <input
            type="number"
            value={Number(value) || ''}
            onChange={(e) => handleChange(field.id, parseFloat(e.target.value))}
            className="w-full px-3 py-2 rounded-lg border border-border focus:border-primary focus:ring-1 focus:ring-primary text-sm"
            required={field.required}
            min={field.range_config?.min}
            max={field.range_config?.max}
            step={field.range_config?.step}
          />
        )}

        {field.type === 'select' && (
          <select
            value={String(value || '')}
            onChange={(e) => handleChange(field.id, e.target.value)}
            className="w-full px-3 py-2 rounded-lg border border-border focus:border-primary focus:ring-1 focus:ring-primary text-sm"
            required={field.required}
          >
            <option value="">Select...</option>
            {(field.options || []).map(opt => (
              <option key={opt.value} value={opt.value}>{opt.label}</option>
            ))}
          </select>
        )}

        {field.type === 'range' && field.range_config && (
          <div className="flex items-center gap-3">
            <input
              type="range"
              value={Number(value) || field.range_config.min}
              onChange={(e) => handleChange(field.id, parseFloat(e.target.value))}
              min={field.range_config.min}
              max={field.range_config.max}
              step={field.range_config.step}
              className="flex-1 accent-primary"
            />
            <span className="text-sm font-medium text-foreground w-16 text-right">
              {Number(value)?.toFixed(1) || field.range_config.min}
            </span>
          </div>
        )}

        {/* Help text */}
        <p className="mt-1 text-xs text-muted-foreground">{field.help_text}</p>

        {/* AI reasoning for this field */}
        {field.ai_reasoning && (
          <p className="mt-1 text-xs text-primary bg-primary/10 px-2 py-1 rounded">
            {field.ai_reasoning}
          </p>
        )}
      </div>
    );
  };

  return (
    <div className="bg-card rounded-xl border border-border shadow-sm overflow-hidden">
      {/* Header with task type and confidence */}
      <div className="bg-muted px-5 py-4 border-b border-border">
        <div className="flex items-center justify-between">
          <div>
            <div className="flex items-center gap-2">
              <Sparkles className="h-5 w-5 text-primary" />
              <h3 className="font-semibold text-foreground">
                {taskType ? taskType.replace('_', ' ').replace(/\b\w/g, l => l.toUpperCase()) : 'Design Parameters'}
              </h3>
            </div>
            {taskType && (
              <p className="text-sm text-muted-foreground mt-1">
                AI-suggested configuration for your design
              </p>
            )}
          </div>
          <div className={`text-sm font-medium ${getConfidenceColor(confidence)}`}>
            {formatConfidence(confidence)} confidence
          </div>
        </div>

        {/* AI Reasoning */}
        {reasoning && (
          <div className="mt-3">
            <button
              type="button"
              onClick={() => setExpandedReasoning(!expandedReasoning)}
              className="text-sm text-primary hover:text-primary/80 flex items-center gap-1"
            >
              {expandedReasoning ? <ChevronUp className="h-4 w-4" /> : <ChevronDown className="h-4 w-4" />}
              {expandedReasoning ? 'Hide' : 'Show'} AI reasoning
            </button>
            {expandedReasoning && (
              <p className="mt-2 text-sm text-foreground bg-card/70 p-3 rounded-lg">
                {reasoning}
              </p>
            )}
          </div>
        )}
      </div>

      {/* Form fields */}
      <form onSubmit={handleSubmit} className="p-5">
        {formConfig.length > 0 ? (
          formConfig.map(renderField)
        ) : (
          <p className="text-sm text-muted-foreground text-center py-4">
            No parameters configured. Describe your goal to get AI suggestions.
          </p>
        )}

        {/* Chat input for refinement */}
        <div className="mt-6 pt-4 border-t border-border">
          <label className="block text-sm font-medium text-foreground mb-2">
            Refine with AI
          </label>
          <div className="flex gap-2">
            <input
              type="text"
              value={chatInput}
              onChange={(e) => setChatInput(e.target.value)}
              placeholder="e.g., make it more aggressive, increase designs to 10..."
              className="flex-1 px-3 py-2 rounded-lg border border-border focus:border-primary focus:ring-1 focus:ring-primary text-sm"
              onKeyDown={(e) => {
                if (e.key === 'Enter' && !e.shiftKey) {
                  e.preventDefault();
                  handleChatSubmit(e);
                }
              }}
            />
            <button
              type="button"
              onClick={handleChatSubmit}
              disabled={!chatInput.trim()}
              className="px-4 py-2 rounded-lg bg-primary/10 text-primary hover:bg-primary/20 disabled:opacity-50 disabled:cursor-not-allowed"
            >
              <Send className="h-4 w-4" />
            </button>
          </div>
        </div>

        {/* Submit button */}
        <div className="mt-6 flex justify-end">
          <button
            type="submit"
            disabled={isLoading || formConfig.length === 0}
            className="px-6 py-2.5 rounded-xl font-medium bg-primary text-primary-foreground hover:bg-primary/90 shadow-sm disabled:opacity-50 disabled:cursor-not-allowed flex items-center gap-2"
          >
            {isLoading ? (
              <>
                <Loader2 className="h-4 w-4 animate-spin" />
                Processing...
              </>
            ) : (
              <>
                <FlaskConical className="h-4 w-4" />
                Run Design
              </>
            )}
          </button>
        </div>
      </form>
    </div>
  );
}
