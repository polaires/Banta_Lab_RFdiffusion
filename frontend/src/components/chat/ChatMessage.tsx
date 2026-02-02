'use client';

import { useState, useEffect } from 'react';
import { Sparkles, User, Info } from 'lucide-react';
import { cn } from '@/lib/utils';
import type { ChatMessage as ChatMessageType } from '@/lib/chat-types';
import { AttachmentRenderer } from './AttachmentRenderer';

interface ChatMessageProps {
  message: ChatMessageType;
  isNew?: boolean;
  onSelectDesign?: (pdbContent: string) => void;
}

// Typing animation for new assistant messages
function useTypingAnimation(text: string, speed: number = 15, enabled: boolean = true) {
  const [displayedText, setDisplayedText] = useState('');
  const [isComplete, setIsComplete] = useState(false);

  useEffect(() => {
    if (!enabled) {
      setDisplayedText(text);
      setIsComplete(true);
      return;
    }
    setDisplayedText('');
    setIsComplete(false);
    let index = 0;
    const timer = setInterval(() => {
      if (index < text.length) {
        setDisplayedText(text.slice(0, index + 1));
        index++;
      } else {
        setIsComplete(true);
        clearInterval(timer);
      }
    }, speed);
    return () => clearInterval(timer);
  }, [text, speed, enabled]);

  return { displayedText, isComplete };
}

export function ChatMessage({ message, isNew = false, onSelectDesign }: ChatMessageProps) {
  if (message.role === 'system') {
    return (
      <div className="flex justify-center">
        <div className="bg-muted text-muted-foreground rounded-full px-4 py-1.5 text-xs flex items-center gap-2">
          <Info className="h-3.5 w-3.5" />
          {message.content}
        </div>
      </div>
    );
  }

  if (message.role === 'user') {
    return (
      <div className="flex gap-3 justify-end">
        <div className="bg-primary text-primary-foreground rounded-2xl rounded-tr-sm px-4 py-3 text-sm max-w-[80%] whitespace-pre-wrap">
          {message.content}
        </div>
        <div className="flex-shrink-0 w-8 h-8 rounded-full bg-muted flex items-center justify-center">
          <User className="h-4 w-4 text-muted-foreground" />
        </div>
      </div>
    );
  }

  // Assistant message
  return <AssistantMessage message={message} isNew={isNew} onSelectDesign={onSelectDesign} />;
}

function AssistantMessage({
  message,
  isNew,
  onSelectDesign,
}: {
  message: ChatMessageType;
  isNew: boolean;
  onSelectDesign?: (pdbContent: string) => void;
}) {
  const { displayedText, isComplete } = useTypingAnimation(message.content, 15, isNew);

  return (
    <div className="flex gap-3">
      <div className="flex-shrink-0 w-8 h-8 rounded-full bg-gradient-to-br from-primary to-primary/80 flex items-center justify-center">
        <Sparkles className="h-4 w-4 text-primary-foreground" />
      </div>
      <div className="flex-1 space-y-3 min-w-0">
        {message.content && (
          <div className="bg-muted rounded-2xl rounded-tl-sm px-4 py-3 text-sm text-foreground leading-relaxed whitespace-pre-wrap">
            {displayedText}
            {!isComplete && <span className="animate-pulse">|</span>}
          </div>
        )}
        {isComplete && message.attachment && (
          <AttachmentRenderer attachment={message.attachment} onSelectDesign={onSelectDesign} />
        )}
      </div>
    </div>
  );
}
