/**
 * Molstar Type Definitions
 * Centralizes Molstar type imports for the application
 */

// Core plugin types
export type { PluginUIContext } from 'molstar/lib/mol-plugin-ui/context';
export type { PluginContext } from 'molstar/lib/mol-plugin/context';

// State types
export type { StateObjectRef, StateObjectCell } from 'molstar/lib/mol-state';
export type { StateTransform } from 'molstar/lib/mol-state/transform';

// Structure types
export type { Structure, StructureElement } from 'molstar/lib/mol-model/structure';
export type { StructureSelection } from 'molstar/lib/mol-model/structure/query/selection';

// Script/Query types
export type { MolScriptBuilder } from 'molstar/lib/mol-script/language/builder';
export type { Expression } from 'molstar/lib/mol-script/language/expression';

// Geometry types
export type { Mesh } from 'molstar/lib/mol-geo/geometry/mesh/mesh';
export type { Shape } from 'molstar/lib/mol-model/shape';

// Color types
export { Color } from 'molstar/lib/mol-util/color';

// Re-export commonly used builders
export { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder';
export { Script } from 'molstar/lib/mol-script/script';
export { StructureSelection as StructSel } from 'molstar/lib/mol-model/structure';
