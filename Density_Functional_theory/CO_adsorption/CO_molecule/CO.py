from ase.build import molecule
from ase.calculators.vasp import Vasp

atoms = molecule('CO')
atoms.center(vacuum=10)
atoms.pbc=True


calc = Vasp(xc='PBE', #Level of theory (density functional)
            encut=400, #plane-wave kinetic energy cutoff (convergence test; I usually use about 520)
            kpts=[1,1,1], #number of k-points, e.g. [3,3,3] (convergence test; fewer kpoints are needed if lattice dimension is large)
            isif=2, #relaxation method (2 = positions, 3 = positions+volume)
            prec='Accurate', #accuracy (don't change this)
            sigma=0.01, #don't need to change (smearing width)
            ismear=0, #don't need to change (smearing method)
            isym=0, #disables symmetry
            ediffg=-0.03, #force tolerance, |ediffg| (eV/A) (I generally use 0.01 --> 0.03; lower is better but will take longer)
            nsw=1000, #maximum number of steps
            ivdw=12, #enable van der Waals interactions
            setups='recommended',
            ibrion=2) #geom opt algorithm (dont need change; default is conjugate gradient -- works fine)

atoms.set_calculator(calc)

atoms.get_potential_energy()