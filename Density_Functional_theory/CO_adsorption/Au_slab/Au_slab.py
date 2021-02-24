from ase.build import fcc111
from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms

Au = fcc111(symbol='Au',a=4.15402,size=(2,2,3),vacuum=10)

# Fixing 1-bottom layer of Au in Z-direction
constraint=FixAtoms(mask=[atom.z < 11 for atom in Au])
Au.set_constraint(constraint)

Au.pbc = True
calc4 = Vasp(xc='PBE', #Level of theory (density functional)
            encut=400, #plane-wave kinetic energy cutoff (convergence test; I usually use about 520)
            kpts=[4,4,1], #number of k-points, e.g. [3,3,3] (convergence test; fewer kpoints are needed if lattice dimension is large)
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

Au.set_calculator(calc4)

Au.get_potential_energy()