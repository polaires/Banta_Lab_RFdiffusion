'use client';

import { useRef, useEffect } from 'react';
import { ScrollArea } from '@/components/ui/scroll-area';
import type { ChatMessage as ChatMessageType } from '@/lib/chat-types';
import { ChatMessage } from './ChatMessage';

interface ChatThreadProps {
  messages: ChatMessageType[];
  onSelectDesign?: (pdbContent: string) => void;
  /** Index of the last message that existed before this render (for typing animation) */
  previousMessageCount?: number;
  children?: React.ReactNode;
}

export function ChatThread({ messages, onSelectDesign, previousMessageCount = 0, children }: ChatThreadProps) {
  const bottomRef = useRef<HTMLDivElement>(null);

  // Auto-scroll on new messages
  useEffect(() => {
    bottomRef.current?.scrollIntoView({ behavior: 'smooth' });
  }, [messages.length]);

  return (
    <div className="flex-1 overflow-y-auto min-h-0">
      <div className="space-y-4 p-4">
        {messages.map((msg, i) => (
          <ChatMessage
            key={msg.id}
            message={msg}
            isNew={i >= previousMessageCount && msg.role === 'assistant'}
            onSelectDesign={onSelectDesign}
          />
        ))}

        {/* Inline children (pipeline card, interview, quick-start cards, etc.) */}
        {children}

        <div ref={bottomRef} />
      </div>
    </div>
  );
}
