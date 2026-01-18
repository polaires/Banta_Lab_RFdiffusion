/**
 * Molstar Focus View Manager
 * Manages focus views without reloading the structure
 */

import type { PluginUIContext } from 'molstar/lib/mol-plugin-ui/context';
import type { StateObjectRef } from 'molstar/lib/mol-state';
import { Color } from 'molstar/lib/mol-util/color';
import { Script } from 'molstar/lib/mol-script/script';
import { StructureSelection } from 'molstar/lib/mol-model/structure';
import {
  WATER_EXPRESSION,
  POLYMER_EXPRESSION,
  residueExpression,
  surroundingsExpression,
  exceptExpression,
  intersectExpression,
} from './molstar-expressions';

export interface FocusOptions {
  radius: number;
  showWaters: boolean;
  showLabels: boolean;
  backgroundAlpha: number;
  carbonColor?: number;
}

interface ComponentRefs {
  primary: StateObjectRef | null;
  surroundings: StateObjectRef | null;
  water: StateObjectRef | null;
  background: StateObjectRef | null;
}

/**
 * Manages focus view state and transitions
 */
export class FocusViewManager {
  private plugin: PluginUIContext;
  private structureRef: StateObjectRef;
  private componentRefs: ComponentRefs = {
    primary: null,
    surroundings: null,
    water: null,
    background: null,
  };
  private currentFocus: { type: 'metal' | 'ligand'; index: number } | null = null;

  constructor(plugin: PluginUIContext, structureRef: StateObjectRef) {
    this.plugin = plugin;
    this.structureRef = structureRef;
  }

  /**
   * Focus on a metal coordination site
   */
  async focusOnMetal(
    element: string,
    chainId: string,
    resSeq: number,
    resName: string,
    options: FocusOptions
  ): Promise<void> {
    await this.clearFocusComponents();

    const metalExpr = residueExpression(resName, chainId, resSeq);
    const coordSphereExpr = surroundingsExpression(metalExpr, options.radius + 2.0);
    const coordWithoutMetal = exceptExpression(coordSphereExpr, metalExpr);
    const coordProtein = exceptExpression(coordWithoutMetal, WATER_EXPRESSION);
    const backgroundExpr = exceptExpression(POLYMER_EXPRESSION, coordSphereExpr);

    // 1. Metal ion - large purple spacefill
    this.componentRefs.primary = await this.createComponent(
      metalExpr,
      `metal-${element}`,
      { type: 'spacefill', color: 0x7C3AED, sizeFactor: 1.5 }
    );

    // 2. Coordination sphere - ball-and-stick
    this.componentRefs.surroundings = await this.createComponent(
      coordProtein,
      'coordination-sphere',
      { type: 'ball-and-stick', color: 'element-symbol', sizeFactor: 0.3 }
    );

    // Add labels if enabled
    if (options.showLabels && this.componentRefs.surroundings) {
      await this.addLabels(this.componentRefs.surroundings);
    }

    // 3. Waters - conditionally
    if (options.showWaters) {
      const waterExpr = intersectExpression(coordWithoutMetal, WATER_EXPRESSION);
      this.componentRefs.water = await this.createComponent(
        waterExpr,
        'coord-water',
        { type: 'ball-and-stick', color: 'element-symbol', sizeFactor: 0.2, alpha: 0.7 }
      );
    }

    // 4. Background protein
    this.componentRefs.background = await this.createComponent(
      backgroundExpr,
      'protein-background',
      { type: 'cartoon', color: 'chain-id', alpha: options.backgroundAlpha }
    );

    // Focus camera
    await this.focusCamera(coordSphereExpr);

    this.currentFocus = { type: 'metal', index: resSeq };
  }

  /**
   * Focus on a ligand binding site
   */
  async focusOnLigand(
    name: string,
    chainId: string,
    resSeq: number,
    options: FocusOptions
  ): Promise<void> {
    await this.clearFocusComponents();

    const ligandExpr = residueExpression(name, chainId, resSeq);
    const bindingSiteExpr = surroundingsExpression(ligandExpr, options.radius);
    const siteWithoutLigand = exceptExpression(bindingSiteExpr, ligandExpr);
    const siteProtein = exceptExpression(siteWithoutLigand, WATER_EXPRESSION);
    const backgroundExpr = exceptExpression(POLYMER_EXPRESSION, bindingSiteExpr);

    // 1. Ligand - ball-and-stick with custom carbon color
    const ligandColor = options.carbonColor || 0x50C878;
    this.componentRefs.primary = await this.createComponent(
      ligandExpr,
      `ligand-${name}`,
      { type: 'ball-and-stick', color: ligandColor, sizeFactor: 0.4 }
    );

    // 2. Binding site residues
    this.componentRefs.surroundings = await this.createComponent(
      siteProtein,
      'binding-site',
      { type: 'ball-and-stick', color: 'element-symbol', sizeFactor: 0.25 }
    );

    if (options.showLabels && this.componentRefs.surroundings) {
      await this.addLabels(this.componentRefs.surroundings);
    }

    // 3. Waters - conditionally
    if (options.showWaters) {
      const waterExpr = intersectExpression(siteWithoutLigand, WATER_EXPRESSION);
      this.componentRefs.water = await this.createComponent(
        waterExpr,
        'binding-water',
        { type: 'ball-and-stick', color: 'element-symbol', sizeFactor: 0.15, alpha: 0.7 }
      );
    }

    // 4. Background protein
    this.componentRefs.background = await this.createComponent(
      backgroundExpr,
      'protein-background',
      { type: 'cartoon', color: 'chain-id', alpha: options.backgroundAlpha }
    );

    // Focus camera
    await this.focusCamera(bindingSiteExpr);

    this.currentFocus = { type: 'ligand', index: resSeq };
  }

  /**
   * Clear focus and return to default view
   */
  async resetFocus(): Promise<void> {
    await this.clearFocusComponents();
    this.currentFocus = null;
    this.plugin.canvas3d?.requestCameraReset();
  }

  /**
   * Update water visibility without full reload
   */
  async toggleWaters(show: boolean, expression: any): Promise<void> {
    if (show && !this.componentRefs.water) {
      this.componentRefs.water = await this.createComponent(
        expression,
        'focus-water',
        { type: 'ball-and-stick', color: 'element-symbol', sizeFactor: 0.15, alpha: 0.7 }
      );
    } else if (!show && this.componentRefs.water) {
      await this.deleteComponent(this.componentRefs.water);
      this.componentRefs.water = null;
    }
  }

  private async createComponent(
    expression: any,
    label: string,
    repr: {
      type: string;
      color: string | number;
      sizeFactor?: number;
      alpha?: number;
    }
  ): Promise<StateObjectRef | null> {
    try {
      const comp = await this.plugin.builders.structure.tryCreateComponentFromExpression(
        this.structureRef,
        expression,
        label
      );

      if (!comp) return null;

      const colorParams = typeof repr.color === 'number'
        ? { value: Color(repr.color) }
        : undefined;

      await this.plugin.builders.structure.representation.addRepresentation(comp, {
        type: repr.type as any,
        color: (typeof repr.color === 'string' ? repr.color : 'uniform') as any,
        colorParams,
        typeParams: {
          sizeFactor: repr.sizeFactor,
          alpha: repr.alpha,
        },
      });

      return comp.ref;
    } catch (err) {
      console.warn(`[FocusView] Failed to create component ${label}:`, err);
      return null;
    }
  }

  private async addLabels(componentRef: StateObjectRef): Promise<void> {
    try {
      const refString = typeof componentRef === 'string' ? componentRef : (componentRef as any).ref;
      const cell = this.plugin.state.data.cells.get(refString);
      if (cell) {
        await this.plugin.builders.structure.representation.addRepresentation(
          { ref: componentRef } as any,
          {
            type: 'label',
            typeParams: { level: 'residue' },
            color: 'uniform',
            colorParams: { value: Color(0x333333) },
          }
        );
      }
    } catch (err) {
      console.warn('[FocusView] Failed to add labels:', err);
    }
  }

  private async deleteComponent(ref: StateObjectRef): Promise<void> {
    try {
      await this.plugin.build().delete(ref).commit();
    } catch (err) {
      console.warn('[FocusView] Failed to delete component:', err);
    }
  }

  private async clearFocusComponents(): Promise<void> {
    const refs = Object.values(this.componentRefs).filter(Boolean) as StateObjectRef[];
    for (const ref of refs) {
      await this.deleteComponent(ref);
    }
    this.componentRefs = {
      primary: null,
      surroundings: null,
      water: null,
      background: null,
    };
  }

  private async focusCamera(expression: any): Promise<void> {
    try {
      const refString = typeof this.structureRef === 'string' ? this.structureRef : (this.structureRef as any).ref;
      const cell = this.plugin.state.data.cells.get(refString);
      if (!cell?.obj?.data) return;

      const selection = Script.getStructureSelection(expression, cell.obj.data);
      if (!StructureSelection.isEmpty(selection)) {
        const loci = StructureSelection.toLociWithSourceUnits(selection);
        this.plugin.managers.camera.focusLoci(loci, { durationMs: 400 });
      }
    } catch (err) {
      console.warn('[FocusView] Failed to focus camera:', err);
    }
  }
}
