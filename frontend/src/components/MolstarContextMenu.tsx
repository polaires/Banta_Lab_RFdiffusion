'use client';

import { useState, useEffect, useCallback } from 'react';

interface ResidueInfo {
  chain: string;
  residue: number;
  name: string;
}

interface ContextMenuState {
  visible: boolean;
  x: number;
  y: number;
  residue: ResidueInfo | null;
}

interface MolstarContextMenuProps {
  onAddResidue: (residue: ResidueInfo, atomType: string) => void;
  onRemoveResidue?: (residue: ResidueInfo) => void;
  isResidueSelected?: (residue: ResidueInfo) => boolean;
}

export function useMolstarContextMenu({
  onAddResidue,
  onRemoveResidue,
  isResidueSelected,
}: MolstarContextMenuProps) {
  const [menuState, setMenuState] = useState<ContextMenuState>({
    visible: false,
    x: 0,
    y: 0,
    residue: null,
  });

  const showMenu = useCallback((x: number, y: number, residue: ResidueInfo) => {
    setMenuState({ visible: true, x, y, residue });
  }, []);

  const hideMenu = useCallback(() => {
    setMenuState((prev) => ({ ...prev, visible: false }));
  }, []);

  // Close menu on click outside
  useEffect(() => {
    if (!menuState.visible) return;

    const handleClick = () => hideMenu();
    const handleEscape = (e: KeyboardEvent) => {
      if (e.key === 'Escape') hideMenu();
    };

    document.addEventListener('click', handleClick);
    document.addEventListener('keydown', handleEscape);

    return () => {
      document.removeEventListener('click', handleClick);
      document.removeEventListener('keydown', handleEscape);
    };
  }, [menuState.visible, hideMenu]);

  const MenuComponent = menuState.visible && menuState.residue && (
    <div
      className="fixed z-50 bg-popover border rounded-md shadow-lg py-1 min-w-[180px]"
      style={{ left: menuState.x, top: menuState.y }}
      onClick={(e) => e.stopPropagation()}
    >
      {/* Header */}
      <div className="px-3 py-1.5 text-sm font-medium border-b">
        {menuState.residue.name} {menuState.residue.chain}{menuState.residue.residue}
      </div>

      {/* Add as catalytic */}
      {!isResidueSelected?.(menuState.residue) && (
        <div className="relative group">
          <button className="w-full text-left px-3 py-2 text-sm hover:bg-accent flex items-center justify-between">
            Add as catalytic residue
            <span className="text-muted-foreground">→</span>
          </button>
          <div className="absolute left-full top-0 ml-1 hidden group-hover:block">
            <div className="bg-popover border rounded-md shadow-lg py-1 min-w-[140px]">
              {[
                { value: 'BKBN', label: 'Backbone only' },
                { value: 'ALL', label: 'All atoms' },
                { value: 'TIP', label: 'Functional tip' },
                { value: '', label: 'Flexible' },
              ].map((opt) => (
                <button
                  key={opt.value}
                  className="w-full text-left px-3 py-1.5 text-sm hover:bg-accent"
                  onClick={() => {
                    if (menuState.residue) {
                      onAddResidue(menuState.residue, opt.value);
                    }
                    hideMenu();
                  }}
                >
                  {opt.label}
                </button>
              ))}
            </div>
          </div>
        </div>
      )}

      {/* Remove from catalytic */}
      {isResidueSelected?.(menuState.residue) && onRemoveResidue && (
        <button
          className="w-full text-left px-3 py-2 text-sm hover:bg-accent text-destructive"
          onClick={() => {
            if (menuState.residue) {
              onRemoveResidue(menuState.residue);
            }
            hideMenu();
          }}
        >
          Remove from catalytic residues
        </button>
      )}
    </div>
  );

  return { showMenu, hideMenu, MenuComponent };
}
