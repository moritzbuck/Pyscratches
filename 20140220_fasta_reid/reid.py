""" small pipeline generating an fasta file with long orfs using glimmer from a fasta files containing multiple contigs

Usage:
  reid.py  -i <infile>

Options:
   -h --help       Show this screen.
   -i <infile>     input fasta file
"""
from docopt import docopt
import os
import sh
from gefes.fasta.single import FASTA

updist=25-1

if __name__ == '__main__':
    arguments = docopt(__doc__)
    print(arguments)

in_file=FASTA(arguments['-i'])        
with FASTA('./temp.file') as out_file:
    for seq in in_file:
        seq.description=seq.description.replace("  ",";",1)
        seq.name=seq.description.split("  ")[0]
        seq.id=seq.name
        out_file.add_seq(seq)
sh.mv('temp.file',arguments['-i'])
