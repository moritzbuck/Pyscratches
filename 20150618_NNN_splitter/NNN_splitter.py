""" to split records in a fasta file if stretches of N are in it


Usage:
  NNN_splitter.py [-N <int>  -o <outfile>] -i <infile> 
  NNN_splitter.py -h
  
Options:
    -o <outfile>    outfile (by default adds a ".XN.split" before the .fasta)
    -i <infile>     infile
    -N <int>        number of Ns
    -h              this help
"""

from docopt import docopt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re

if __name__ == '__main__':
    arguments = docopt(__doc__)
    infile = arguments['-i']
    Ns = 10
    outfile = ".".join(infile.split(".")[:-1] + [ str(Ns) + "Ns.split"] + [ infile.split(".")[-1]])
    if arguments['-o']:
        outfile = arguments['-o']
    if arguments['-N']:
        Ns = int(arguments['-N'])
        
else:
    Ns = 10
    infile = "bin_62.fasta"
    outfile = ".".join(infile.split(".")[:-1] + [ str(Ns) + "Ns.split"] + [ infile.split(".")[-1]])

with open(infile) as handle:
    seqs = [s for s in SeqIO.parse(handle, "fasta")]


out_seqs = []
pattern = 'N{%i,}' % Ns
for s in seqs:
    chars = str(s.seq)
    if "N"*Ns in chars:
        sub_seqs = re.split(pattern, str(s.seq))
        out_seqs +=[SeqRecord(Seq(t), id=s.id+"_"+str(i), description ="") for i,t in enumerate(sub_seqs)]
    else :
        out_seqs += [s]

with open(outfile, "w") as handle:
    SeqIO.write(out_seqs, handle, "fasta")        
    
