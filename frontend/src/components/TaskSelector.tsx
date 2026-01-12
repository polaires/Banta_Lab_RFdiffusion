'use client';

export type DesignTask =
  | 'denovo'
  | 'protein_binder'
  | 'small_molecule'
  | 'nucleic_acid'
  | 'enzyme'
  | 'symmetric'
  | 'refinement'
  | 'interface_ligand';

interface TaskConfig {
  id: DesignTask;
  name: string;
  description: string;
  icon: string;
  requirements: string[];
}

const TASKS: TaskConfig[] = [
  {
    id: 'denovo',
    name: 'De Novo Protein',
    description: 'Generate a new protein structure from scratch',
    icon: 'add_circle',
    requirements: ['Length specification'],
  },
  {
    id: 'protein_binder',
    name: 'Protein Binder',
    description: 'Design a protein that binds to a target protein',
    icon: 'hub',
    requirements: ['Target PDB', 'Hotspot residues'],
  },
  {
    id: 'small_molecule',
    name: 'Small Molecule Binder',
    description: 'Design a binding pocket for a ligand (ATP, NAD, etc.)',
    icon: 'science',
    requirements: ['PDB with ligand', 'Ligand code'],
  },
  {
    id: 'nucleic_acid',
    name: 'Nucleic Acid Binder',
    description: 'Design a protein that binds DNA or RNA',
    icon: 'gesture',
    requirements: ['PDB with DNA/RNA', 'NA chain selection'],
  },
  {
    id: 'enzyme',
    name: 'Enzyme Scaffold',
    description: 'Build a protein scaffold around an active site',
    icon: 'biotech',
    requirements: ['Theozyme PDB', 'Catalytic residues'],
  },
  {
    id: 'symmetric',
    name: 'Symmetric Oligomer',
    description: 'Design homo-oligomers (dimers, trimers, etc.)',
    icon: 'hexagon',
    requirements: ['Symmetry type', 'Subunit length'],
  },
  {
    id: 'refinement',
    name: 'Structure Refinement',
    description: 'Refine an existing structure with partial diffusion',
    icon: 'tune',
    requirements: ['Input PDB', 'Noise level'],
  },
  {
    id: 'interface_ligand',
    name: 'Interface Ligand Dimer',
    description: 'Design separable protein dimers with ligand at the interface',
    icon: 'link',
    requirements: ['Ligand SMILES', 'Approach'],
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
        <h2 className="text-lg font-bold text-slate-900 mb-2">
          What would you like to design?
        </h2>
        <p className="text-slate-500 text-sm">
          Select a design task to see relevant options and recommended settings
        </p>
      </div>

      {/* Task Grid */}
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-3">
        {TASKS.map((task) => (
          <button
            key={task.id}
            onClick={() => onTaskSelect(task.id)}
            className="group p-4 rounded-lg border border-slate-200 bg-white
                       text-left transition-all hover:border-blue-400 hover:shadow-md
                       focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-blue-500"
          >
            {/* Icon and Title */}
            <div className="flex items-start gap-3 mb-2">
              <div className="w-9 h-9 rounded-lg bg-slate-100 flex items-center justify-center
                             group-hover:bg-blue-600 transition-colors flex-shrink-0">
                <span className="material-symbols-outlined text-slate-600 group-hover:text-white transition-colors text-lg">
                  {task.icon}
                </span>
              </div>
              <div className="flex-1 min-w-0">
                <h3 className="font-semibold text-slate-900 text-sm group-hover:text-blue-700 transition-colors">
                  {task.name}
                </h3>
                <p className="text-xs text-slate-500 mt-0.5 line-clamp-2">
                  {task.description}
                </p>
              </div>
            </div>

            {/* Requirements */}
            <div className="flex flex-wrap gap-1 mt-2 pl-12">
              {task.requirements.map((req) => (
                <span
                  key={req}
                  className="text-[10px] px-1.5 py-0.5 rounded bg-slate-100 text-slate-500 font-medium"
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
        <p className="text-xs text-slate-400">
          Not sure which to choose? Start with{' '}
          <button
            onClick={() => onTaskSelect('denovo')}
            className="text-blue-600 hover:underline font-medium"
          >
            De Novo Protein
          </button>{' '}
          for the simplest workflow.
        </p>
      </div>
    </div>
  );
}
