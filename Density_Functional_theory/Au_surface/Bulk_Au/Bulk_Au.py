
from ase.lattice.cubic import FaceCenteredCubic

Bulk_Au = FaceCenteredCubic(symbol='Au',size=(1,1,1),pbc=True)

# Attach VASP calculation to Optimize Nuclei position
''' 
    From the calculation, the lattice constants can be obtained,
    and the constants are used to creat Surface_Au in the next
    step.
'''
from ase.calculators.vasp import Vasp

calc = Vasp(xc='PBE', #Level of theory (density functional)
            encut=520, #plane-wave kinetic energy cutoff (convergence test; I usually use about 520)
            kpts=[10,10,10], #number of k-points, e.g. [3,3,3] (convergence test; fewer kpoints are needed if lattice dimension is large)
            isif=3, #relaxation method (2 = positions, 3 = positions+volume)
            prec='Accurate', #accuracy (don't change this)
            sigma=0.01, #don't need to change (smearing width)
            ismear=0, #don't need to change (smearing method)
            isym=0, #disables symmetry
            ediffg=-0.03, #force tolerance, |ediffg| (eV/A) (I generally use 0.01 --> 0.03; lower is better but will take longer)
            nsw=1000, #maximum number of steps
            #ivdw=12, #enable van der Waals interactions
            ibrion=2) #geom opt algorithm (dont need change; default is conjugate gradient -- works fine)

Bulk_Au.set_calculator(calc)

Bulk_Au.get_potential_energy()