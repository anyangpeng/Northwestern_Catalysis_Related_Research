#!/bin/bash
#SBATCH --job-name="ASE"
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH -t 4:00:00
#SBATCH -p short
#SBATCH -A p31374
#SBATCH --mail-user=l1y9y3@u.northwestern.edu
#SBATCH --mail-type=END
#SBATCH --constraint=[quest5|quest6|quest8]

module load vasp/5.4.4-vtst-openmpi-4.0.5-intel-19.0.5.281
export VASP_COMMAND='mpirun -n 24 vasp_std > vasp.out'
python encut580.py > encut580.out