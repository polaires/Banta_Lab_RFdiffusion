'use client';

import { useState, useRef, useCallback, useEffect } from 'react';
import { Send } from 'lucide-react';
import { Button } from '@/components/ui/button';
import { cn } from '@/lib/utils';

interface ChatInputProps {
  onSend: (text: string) => void;
  placeholder?: string;
  className?: string;
}

export function ChatInput({ onSend, placeholder = 'What protein would you like to design?', className }: ChatInputProps) {
  const [value, setValue] = useState('');
  const [isSending, setIsSending] = useState(false);
  const textareaRef = useRef<HTMLTextAreaElement>(null);

  // Auto-resize textarea
  const adjustHeight = useCallback(() => {
    const el = textareaRef.current;
    if (!el) return;
    el.style.height = 'auto';
    el.style.height = Math.min(el.scrollHeight, 150) + 'px';
  }, []);

  useEffect(() => {
    adjustHeight();
  }, [value, adjustHeight]);

  const handleSend = useCallback(async () => {
    const trimmed = value.trim();
    if (!trimmed || isSending) return;

    setIsSending(true);
    setValue('');
    try {
      onSend(trimmed);
    } finally {
      setIsSending(false);
    }

    // Refocus
    requestAnimationFrame(() => textareaRef.current?.focus());
  }, [value, isSending, onSend]);

  const handleKeyDown = useCallback(
    (e: React.KeyboardEvent) => {
      if (e.key === 'Enter' && !e.shiftKey) {
        e.preventDefault();
        handleSend();
      }
    },
    [handleSend]
  );

  return (
    <div className={cn('flex gap-3 items-end', className)}>
      <div className="flex-1 relative">
        <textarea
          ref={textareaRef}
          value={value}
          onChange={(e) => setValue(e.target.value)}
          onKeyDown={handleKeyDown}
          placeholder={placeholder}
          className="w-full px-4 py-3 bg-muted rounded-xl border border-border focus:border-primary focus:ring-2 focus:ring-ring/20 focus:outline-none text-foreground text-sm transition-all resize-none min-h-[48px] max-h-[150px]"
          rows={1}
        />
      </div>
      <Button
        onClick={handleSend}
        disabled={!value.trim() || isSending}
        size="icon"
        className="h-12 w-12 rounded-xl shrink-0"
      >
        <Send className="h-5 w-5" />
      </Button>
    </div>
  );
}
