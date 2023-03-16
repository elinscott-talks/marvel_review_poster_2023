from ase.build import bulk
from koopmans.kpoints import Kpoints
from koopmans.workflows import SinglepointWorkflow

# Use ASE to create bulk silicon
atoms = bulk('Si')

# Create the workflow
workflow = SinglepointWorkflow(atoms=atoms,
        parameters={'pseudo_library': 'pseudo_dojo_standard', 'base_functional': 'lda', 'init_orbitals': 'mlwfs'}
        kpoints = Kpoints(grid=[4, 4, 4], path='LGXKG', cell=atoms.cell),
        calculator_parameters = {'pw': {'nbnd': 20}, 'w90': {'auto_projections': True}}

# Run the workflow
workflow.run()
