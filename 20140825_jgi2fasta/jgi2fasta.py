""" to split a fasta file exported from a IMG scaffold cart into individual genomes


Usage:
  jgi2fasta.py [-o <outfolder>] -f <infasta> -c <scaffoldCart.xls> 
  jgi2fasta.py -h
  
Options:
    -o <outfolder>  output folder (default current working dir)
    -f <infasta>    input fasta file obtained from IMG
    -c <scaff.xls>  input table obtained by exporting the scaffold cart from IMG 
    -h              this help
"""


from Bio import SeqIO
from pandas import DataFrame
import os
from docopt import docopt

if __name__ == '__main__':
    arguments = docopt(__doc__)
    input_fasta = arguments['-f']
    input_metadata = arguments['-c']
    path = "."
    if arguments['-o']:
        path = arguments['-o']
        
else:
    input_fasta = "od1s.fasta"
    input_metadata = "od1s.xls"
    path = "/home/moritz/glob/data/genomes/OD1s/from_IMG/"

dir = path if path[-1] == "/" else path + "/"
if not os.path.exists(dir) :
    os.makedirs(dir)
    
with open(input_fasta,"r") as file:
    seqs = [seq for seq in SeqIO.parse(file,"fasta") ]
    
md = DataFrame.from_csv(input_metadata,sep="\t",header=0,index_col=0)

seq_sets = {}
for g in set(md['Genome'].values):
    seq_sets[g] = []
    
for i,s in enumerate(seqs):
    seq_sets[md.iloc[i]['Genome']].append(s)

for g in seq_sets:
    clean_g = g.replace(" ","_") 
    clean_g = clean_g.replace(",","_")
    clean_g = clean_g.replace("-","_")
    clean_g = clean_g.replace("(","_")
    clean_g = clean_g.replace(")","_")
    clean_g = clean_g.replace("/","_")


    
    with open(dir +  clean_g + ".fasta","w") as file:
        SeqIO.write(seq_sets[g],file,"fasta")
