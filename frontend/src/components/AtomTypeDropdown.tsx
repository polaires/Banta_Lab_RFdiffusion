'use client';

import { useState } from 'react';
import {
  Popover,
  PopoverContent,
  PopoverTrigger,
} from '@/components/ui/popover';

interface AtomTypeDropdownProps {
  onSelect: (atomType: string) => void;
  children: React.ReactNode;
  align?: 'start' | 'center' | 'end';
}

const ATOM_TYPES = [
  { value: 'BKBN', label: 'Backbone only', desc: 'N, CA, C, O' },
  { value: 'ALL', label: 'All atoms', desc: 'Entire residue' },
  { value: 'TIP', label: 'Functional tip', desc: 'Side chain tip' },
  { value: '', label: 'Flexible', desc: 'No constraints' },
];

export function AtomTypeDropdown({ onSelect, children, align = 'end' }: AtomTypeDropdownProps) {
  const [open, setOpen] = useState(false);

  const handleSelect = (atomType: string) => {
    onSelect(atomType);
    setOpen(false);
  };

  return (
    <Popover open={open} onOpenChange={setOpen}>
      <PopoverTrigger asChild>
        {children}
      </PopoverTrigger>
      <PopoverContent align={align} className="w-48 p-1">
        <div className="space-y-0.5">
          {ATOM_TYPES.map((type) => (
            <button
              key={type.value}
              onClick={() => handleSelect(type.value)}
              className="w-full text-left px-3 py-2 rounded-md hover:bg-accent transition-colors"
            >
              <div className="text-sm font-medium">{type.label}</div>
              <div className="text-xs text-muted-foreground">{type.desc}</div>
            </button>
          ))}
        </div>
      </PopoverContent>
    </Popover>
  );
}
