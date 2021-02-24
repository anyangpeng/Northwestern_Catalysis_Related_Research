from ase.build import fcc111
from ase.calculators.vasp import Vasp
Au = fcc111(symbol='Au',a=4.15402,size=(2,2,6),vacuum=10)# Using optimized lattice constant to build Au_fcc111 surface with 6 layers of Au 


Au.pbc = True
calc4 = Vasp(xc='PBE', #Level of theory (density functional)
            encut=520, #plane-wave kinetic energy cutoff (convergence test; I usually use about 520)
            kpts=[8,8,1], #number of k-points, e.g. [3,3,3] (convergence test; fewer kpoints are needed if lattice dimension is large)
            prec='Accurate', #accuracy (don't change this)
            sigma=0.01, #don't need to change (smearing width)
            ismear=0, #don't need to change (smearing method)
            isym=0, #disables symmetry
            ivdw=12, #enable van der Waals interactions
            setups='recommended')

Au.set_calculator(calc4)

Au.get_potential_energy()