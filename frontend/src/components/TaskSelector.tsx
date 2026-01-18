'use client';

import { LucideIcon, PlusCircle, Network, FlaskConical, Dna, Beaker, Hexagon, SlidersHorizontal, Link, Circle, Sparkles } from 'lucide-react';

export type DesignTask =
  | 'denovo'
  | 'protein_binder'
  | 'small_molecule'
  | 'nucleic_acid'
  | 'enzyme'
  | 'symmetric'
  | 'refinement'
  | 'interface_ligand'
  | 'interface_metal'
  | 'interface_metal_ligand';

interface TaskConfig {
  id: DesignTask;
  name: string;
  description: string;
  Icon: LucideIcon;
  requirements: string[];
}

const TASKS: TaskConfig[] = [
  {
    id: 'denovo',
    name: 'De Novo Protein',
    description: 'Generate a new protein structure from scratch',
    Icon: PlusCircle,
    requirements: ['Length specification'],
  },
  {
    id: 'protein_binder',
    name: 'Protein Binder',
    description: 'Design a protein that binds to a target protein',
    Icon: Network,
    requirements: ['Target PDB', 'Hotspot residues'],
  },
  {
    id: 'small_molecule',
    name: 'Small Molecule Binder',
    description: 'Design a binding pocket for a ligand (ATP, NAD, etc.)',
    Icon: FlaskConical,
    requirements: ['PDB with ligand', 'Ligand code'],
  },
  {
    id: 'nucleic_acid',
    name: 'Nucleic Acid Binder',
    description: 'Design a protein that binds DNA or RNA',
    Icon: Dna,
    requirements: ['PDB with DNA/RNA', 'NA chain selection'],
  },
  {
    id: 'enzyme',
    name: 'Enzyme Scaffold',
    description: 'Build a protein scaffold around an active site',
    Icon: Beaker,
    requirements: ['Theozyme PDB', 'Catalytic residues'],
  },
  {
    id: 'symmetric',
    name: 'Symmetric Oligomer',
    description: 'Design homo-oligomers (dimers, trimers, etc.)',
    Icon: Hexagon,
    requirements: ['Symmetry type', 'Subunit length'],
  },
  {
    id: 'refinement',
    name: 'Structure Refinement',
    description: 'Refine an existing structure with partial diffusion',
    Icon: SlidersHorizontal,
    requirements: ['Input PDB', 'Noise level'],
  },
  {
    id: 'interface_ligand',
    name: 'Interface Ligand Dimer',
    description: 'Design separable protein dimers with ligand at the interface',
    Icon: Link,
    requirements: ['Ligand SMILES', 'Approach'],
  },
  {
    id: 'interface_metal',
    name: 'Interface Metal Dimer',
    description: 'Design protein heterodimers with metal coordination at the interface',
    Icon: Circle,
    requirements: ['Metal ion', 'Coordination split'],
  },
  {
    id: 'interface_metal_ligand',
    name: 'Metal-Ligand Complex Dimer',
    description: 'Design homodimers binding metal-ligand complexes (e.g., citrate-Tb, PQQ-Ca)',
    Icon: Sparkles,
    requirements: ['Complex template', 'Ligand SMILES'],
  },
];

interface TaskSelectorProps {
  onTaskSelect: (task: DesignTask) => void;
}

export function TaskSelector({ onTaskSelect }: TaskSelectorProps) {
  return (
    <div className="space-y-6">
      {/* Header */}
      <div className="text-center pb-4">
        <h2 className="text-lg font-bold text-foreground mb-2">
          What would you like to design?
        </h2>
        <p className="text-muted-foreground text-sm">
          Select a design task to see relevant options and recommended settings
        </p>
      </div>

      {/* Task Grid */}
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-3">
        {TASKS.map((task) => (
          <button
            key={task.id}
            onClick={() => onTaskSelect(task.id)}
            className="group p-4 rounded-lg border border-border bg-card
                       text-left transition-all hover:border-primary hover:shadow-md
                       focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-ring"
          >
            {/* Icon and Title */}
            <div className="flex items-start gap-3 mb-2">
              <div className="w-9 h-9 rounded-lg bg-muted flex items-center justify-center
                             group-hover:bg-primary transition-colors flex-shrink-0">
                <task.Icon className="h-5 w-5 text-muted-foreground group-hover:text-primary-foreground transition-colors" />
              </div>
              <div className="flex-1 min-w-0">
                <h3 className="font-semibold text-foreground text-sm group-hover:text-primary transition-colors">
                  {task.name}
                </h3>
                <p className="text-xs text-muted-foreground mt-0.5 line-clamp-2">
                  {task.description}
                </p>
              </div>
            </div>

            {/* Requirements */}
            <div className="flex flex-wrap gap-1 mt-2 pl-12">
              {task.requirements.map((req) => (
                <span
                  key={req}
                  className="text-[10px] px-1.5 py-0.5 rounded bg-muted text-muted-foreground font-medium"
                >
                  {req}
                </span>
              ))}
            </div>
          </button>
        ))}
      </div>

      {/* Help Text */}
      <div className="text-center pt-2">
        <p className="text-xs text-muted-foreground">
          Not sure which to choose? Start with{' '}
          <button
            onClick={() => onTaskSelect('denovo')}
            className="text-primary hover:underline font-medium"
          >
            De Novo Protein
          </button>{' '}
          for the simplest workflow.
        </p>
      </div>
    </div>
  );
}
