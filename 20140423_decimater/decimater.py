from Bio import SeqIO
from random import uniform
from tqdm import tqdm
import sys
import os

#simple script to keep only x % of a fasta file
# usage : python decimater.py file_to_decimate.fasta probabilty_to_keep_a_seq
# [propabilty_to_keep_a_seq is a float between 0 and 1, obiviously... basicall 1/desired_size_reduction]

ifile=sys.argv[1]
input = open(ifile, "r")
fraction = float(sys.argv[2])
ofile="deci_"+str(fraction) + "_" + os.path.basename(ifile)
russians = []

#pick which ones too keep
for record in tqdm(SeqIO.parse(ifile, "fasta")) :
    if uniform(0,1) < fraction: russians.append(record)
    
#write them to a file         
with open(ofile, "w") as output:
    SeqIO.write(russians,output,"fasta")   
    
