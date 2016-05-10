#!/bin/bash
#SBATCH -D /home/moritz/repos/Pyscratches/20150610_orthomcl_tools/
#SBATCH -J raxml
#SBATCH -o /home/moritz/repos/Pyscratches/20150610_orthomcl_tools/slurm_scratch.out
#SBATCH -e /home/moritz/repos/Pyscratches/20150610_orthomcl_tools/slurm_scratch.err
#SBATCH -A b2014318
#SBATCH -t 6-13:30:00
#SBATCH -n 16
#SBATCH -p node
#SBATCH --mail-user murumbii@gmail.com
#SBATCH --mail-type=ALL

module load raxml

python orthmcl_tools/__init__.py
