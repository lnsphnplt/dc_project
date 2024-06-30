#!/bin/bash

#SBATCH --time=700
#SBATCH --mem-per-cpu=150G
#SBATCH --job-name=CD_testRun
#SBATCH --open-mode=truncate
#SBATCH --output=output.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emathier@ethz.ch

module purge
module load gcc/11.4.0
module load python/3.11.6
echo "Loaded modules"
python /cluster/scratch/emathier/CD2.py
