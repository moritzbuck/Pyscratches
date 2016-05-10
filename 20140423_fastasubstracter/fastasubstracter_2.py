""" small script lines from a a mothur count file <file1> according to a list file <file2>


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

    
count_file = arguments['<file1>']
list_file = arguments['<file2>']

print "parse the list file"
with open(list_file) as ofile:
    id_list = ofile.readlines()

id_list = id_list[0]
id_list = id_list.split("\t")[2]
id_list = id_list.split(",")

with open(count_file) as ifile:
    lines = ifile.readlines()

with open("reduced.count_table","w") as ofile:
    for line in lines :
        if line.split("\t")[0] in id_list:
            ofile.write(line)

