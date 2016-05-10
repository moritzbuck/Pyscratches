#!/bin/bash
#SBATCH -D /home/moritz/repos/Pyscratches/20151130_read_shuffler/
#SBATCH -J gambit
#SBATCH -o /home/moritz/repos/Pyscratches/20151130_read_shuffler/gambit.out
#SBATCH -e /home/moritz/repos/Pyscratches/20151130_read_shuffler/gambit.err
#SBATCH -A b2011105
#SBATCH -t 7-00:00:00
#SBATCH -n 1
#SBATCH -p node
#SBATCH --mail-user murumbii@gmail.com
#SBATCH --mail-type=ALL

python read_shuffler.py
