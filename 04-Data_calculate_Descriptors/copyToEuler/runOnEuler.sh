#!/bin/bash

#SBATCH --time=25
#SBATCH --mem-per-cpu=300M
#SBATCH --job-name=mp_CD
#SBATCH --open-mode=append
#SBATCH --output=output.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emathier@ethz.ch
#SBATCH --ntasks=48

echo "Loading modules"
module purge
module load gcc/11.4.0
module load python/3.11.6
echo "Loaded modules"
echo "Installing dependencies"
pip install pandas numpy rdkit
python $SCRATCH/CD2.py
