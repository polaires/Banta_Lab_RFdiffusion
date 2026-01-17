'use client';

import { useState, useRef, useEffect } from 'react';
import { MessageSquare, FlaskConical, Sparkles, Loader2, Send } from 'lucide-react';
import type {
  TaskType,
  AIChatMessage,
  FormFieldConfig,
} from '@/lib/ai-types';
import { getConfidenceColor, formatConfidence } from '@/lib/ai-types';

interface AIChatProps {
  messages: AIChatMessage[];
  onSendMessage: (content: string) => Promise<void>;
  onTaskDetected?: (taskType: TaskType, params: Record<string, unknown>, formConfig: FormFieldConfig[]) => void;
  isLoading?: boolean;
  structureInfo?: {
    pdbId?: string;
    chains?: string[];
    numResidues?: number;
  } | null;
}

/**
 * Conversational interface that works alongside SmartForm.
 *
 * Modes:
 * 1. Initial: User describes goal -> AI suggests task + params
 * 2. Refinement: User asks to adjust -> AI updates suggestions
 * 3. Explanation: User asks why -> AI explains reasoning
 * 4. Alternative: User wants different approach -> AI regenerates
 */
export function AIChat({
  messages,
  onSendMessage,
  onTaskDetected,
  isLoading = false,
  structureInfo,
}: AIChatProps) {
  const [input, setInput] = useState('');
  const messagesEndRef = useRef<HTMLDivElement>(null);
  const inputRef = useRef<HTMLInputElement>(null);

  // Scroll to bottom on new messages
  useEffect(() => {
    messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  }, [messages]);

  // Focus input on mount
  useEffect(() => {
    inputRef.current?.focus();
  }, []);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    if (!input.trim() || isLoading) return;

    const message = input.trim();
    setInput('');
    await onSendMessage(message);
  };

  const renderMessage = (message: AIChatMessage) => {
    const isUser = message.role === 'user';
    const isSystem = message.role === 'system';

    return (
      <div
        key={message.id}
        className={`flex ${isUser ? 'justify-end' : 'justify-start'} mb-4`}
      >
        <div
          className={`max-w-[85%] rounded-2xl px-4 py-3 ${
            isUser
              ? 'bg-violet-600 text-white rounded-br-md'
              : isSystem
              ? 'bg-amber-50 text-amber-900 border border-amber-200'
              : 'bg-slate-100 text-slate-900 rounded-bl-md'
          }`}
        >
          {/* Message content */}
          <div className={`text-sm ${isUser ? '' : 'prose prose-sm prose-slate max-w-none'}`}>
            {message.content}
          </div>

          {/* Data payload (if analysis response) */}
          {message.data?.type === 'analysis' && message.data.payload && (
            <div className="mt-3 pt-3 border-t border-slate-200/50">
              <div className="flex items-center gap-2 text-xs">
                <span className={`font-medium ${getConfidenceColor((message.data.payload as any).confidence || 0)}`}>
                  {formatConfidence((message.data.payload as any).confidence || 0)}
                </span>
                <span className="text-slate-500">confidence</span>
              </div>
            </div>
          )}

          {/* Timestamp */}
          <div
            className={`text-xs mt-2 ${
              isUser ? 'text-violet-200' : 'text-slate-400'
            }`}
          >
            {new Date(message.timestamp).toLocaleTimeString([], {
              hour: '2-digit',
              minute: '2-digit',
            })}
          </div>
        </div>
      </div>
    );
  };

  const suggestedPrompts = [
    'Design a protein binder for my target',
    'Create a 60-residue binder',
    'Design a zinc-binding enzyme',
    'Create a C2 symmetric dimer',
  ];

  return (
    <div className="flex flex-col h-full bg-white rounded-xl border border-slate-200 shadow-sm overflow-hidden">
      {/* Header */}
      <div className="bg-gradient-to-r from-violet-50 to-purple-50 px-5 py-4 border-b border-violet-100">
        <div className="flex items-center gap-2">
          <MessageSquare className="h-5 w-5 text-violet-600" />
          <h3 className="font-semibold text-slate-900">AI Assistant</h3>
        </div>
        <p className="text-sm text-slate-600 mt-1">
          Describe your protein design goal in natural language
        </p>
        {structureInfo?.pdbId && (
          <div className="mt-2 flex items-center gap-2 text-xs text-violet-600 bg-violet-100/50 rounded-full px-3 py-1 w-fit">
            <FlaskConical className="h-4 w-4" />
            Working with {structureInfo.pdbId}
            {structureInfo.numResidues && ` (${structureInfo.numResidues} residues)`}
          </div>
        )}
      </div>

      {/* Messages */}
      <div className="flex-1 overflow-y-auto p-4 min-h-[300px]">
        {messages.length === 0 ? (
          <div className="text-center py-8">
            <div className="inline-flex items-center justify-center w-16 h-16 rounded-full bg-violet-100 mb-4">
              <Sparkles className="h-8 w-8 text-violet-600" />
            </div>
            <h4 className="font-medium text-slate-900 mb-2">
              What would you like to design?
            </h4>
            <p className="text-sm text-slate-500 mb-4">
              Describe your goal and I&apos;ll suggest the best approach
            </p>

            {/* Suggested prompts */}
            <div className="flex flex-wrap justify-center gap-2">
              {suggestedPrompts.map((prompt, i) => (
                <button
                  key={i}
                  onClick={() => {
                    setInput(prompt);
                    inputRef.current?.focus();
                  }}
                  className="text-xs px-3 py-1.5 rounded-full bg-slate-100 text-slate-600 hover:bg-violet-100 hover:text-violet-700 transition-colors"
                >
                  {prompt}
                </button>
              ))}
            </div>
          </div>
        ) : (
          <>
            {messages.map(renderMessage)}
            <div ref={messagesEndRef} />
          </>
        )}

        {/* Loading indicator */}
        {isLoading && (
          <div className="flex justify-start mb-4">
            <div className="bg-slate-100 rounded-2xl rounded-bl-md px-4 py-3">
              <div className="flex items-center gap-2 text-sm text-slate-500">
                <Loader2 className="h-4 w-4 animate-spin" />
                Thinking...
              </div>
            </div>
          </div>
        )}
      </div>

      {/* Input */}
      <form onSubmit={handleSubmit} className="p-4 border-t border-slate-200 bg-slate-50">
        <div className="flex gap-2">
          <input
            ref={inputRef}
            type="text"
            value={input}
            onChange={(e) => setInput(e.target.value)}
            placeholder="Describe your design goal..."
            disabled={isLoading}
            className="flex-1 px-4 py-2.5 rounded-xl border border-slate-300 focus:border-violet-500 focus:ring-2 focus:ring-violet-200 text-sm disabled:opacity-50 disabled:bg-slate-100"
          />
          <button
            type="submit"
            disabled={!input.trim() || isLoading}
            className="px-4 py-2.5 rounded-xl bg-violet-600 text-white hover:bg-violet-700 disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
          >
            <Send className="h-5 w-5" />
          </button>
        </div>
        <p className="text-xs text-slate-400 mt-2 text-center">
          Press Enter to send, or click suggested prompts above
        </p>
      </form>
    </div>
  );
}
