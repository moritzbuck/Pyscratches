#!/home/moritz/.pyenv/shims/python -u
#SBATCH -D /home/moritz/repos/Pyscratches/20151130_read_shuffler/
#SBATCH -J gambit
#SBATCH -o /home/moritz/repos/Pyscratches/20151130_read_shuffler/gambit_2.out
#SBATCH -e /home/moritz/repos/Pyscratches/20151130_read_shuffler/gambit_2.err
#SBATCH -A b2013127
#SBATCH -t 7-00:00:00
#SBATCH -n 1
#SBATCH -p core
#SBATCH --mail-user murumbii@gmail.com
#SBATCH --mail-type=ALL

""" splits the table into taxonomical levels

Usage:
  metaphlantablesplit.py -i <infile> -o <outheader>
  read_shuffler.py -h
  
Options:
    -o <outheader>   output header
    -i <infile>    input file
"""

from docopt import docopt
import pandas
import sys
import os


if __name__ == '__main__':
    arguments = docopt(__doc__)
    input_file = arguments['-i']
    output_head = arguments['-o']

raw_data = pandas.read_table(input_file,sep="\t", comment="#", index_col=0)
raw_data.columns = pandas.Index([c.replace(".metaphlan","") for c in raw_data.columns])
levels = [ "Kingdom" , "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"  ]

for i in range(len(levels)):
    rows = [len(f.split("|")) == i+1 for f in raw_data.index]
    subtable = raw_data.iloc[rows]
    subtable.to_csv(output_head + "_" + levels[i] + "_" + str(i+1) + ".csv")
