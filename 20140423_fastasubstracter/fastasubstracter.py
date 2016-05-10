""" small script substracting fasta files from each other using the ID line (e.g. removes from file1 all the sequences that are in file2 and saves the resulting fasta into outpu


Usage:
  fastasubstracter.py <file1> <file2> <output>
  fastasubstracter.py -h
  
Options:
   -h --help       Show this screen.
"""
from docopt import docopt
from Bio import SeqIO
from tqdm import tqdm
import sys
import os


if __name__ == '__main__':
    arguments = docopt(__doc__)
    print(arguments)

big_one = arguments['<file1>']
small_one = arguments['<file2>']

print "parse the small file"
small_rec_iter = SeqIO.parse(open(small_one,"rU"),"fasta")
to_remove = [r.id for r in tqdm(small_rec_iter)]

big_rec_iter = SeqIO.parse(open(big_one,"rU"),"fasta")

print "parse the large file"
to_output = (r for r in tqdm(big_rec_iter) if r.id not in to_remove)

print "write the out-file"
#write them to a file         
with open(arguments['<output>'], "w") as output:
    SeqIO.write(to_output, output, "fasta")   
    
