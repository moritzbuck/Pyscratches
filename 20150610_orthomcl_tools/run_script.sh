#!/bin/bash
#SBATCH -D /home/moritz/repos/Pyscratches/20150610_orthomcl_tools/
#SBATCH -J stuff
#SBATCH -o /home/moritz/repos/Pyscratches/20150610_orthomcl_tools/stuff.out
#SBATCH -e /home/moritz/repos/Pyscratches/20150610_orthomcl_tools/stuff.err
#SBATCH -A b2011032
#SBATCH -t 7-00:00:00
#SBATCH -n 16
#SBATCH -p node
#SBATCH --mail-user murumbii@gmail.com
#SBATCH --mail-type=ALL


module load raxml
module load blast/2.2.26
python /home/moritz/repos/Pyscratches/20150610_orthomcl_tools/orthmcl_tools/__init__.py
