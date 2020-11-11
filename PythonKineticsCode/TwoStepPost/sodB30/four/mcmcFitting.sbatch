#!/bin/bash
#SBATCH --job-name=sodB30_seven
#SBATCH --partition=broadwl
#SBATCH --time=33:59:00
#SBATCH --ntasks=50
#SBATCH --output=output.out
#SBATCH --error=err.out

module load python/3.6.1+intel-16.0
module load intelmpi/5.1+intel-16.0
mpirun python3 sodB130_plus_delta.py

