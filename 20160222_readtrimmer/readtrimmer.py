""" trim n-first and last bases from a bunch of reads

Usage:
  readtrimmer.py [-f -g -n <bps>] -i <infastq> -o <outfastq> 
  readtrimmer.py -h
  
Options:
    -o <outfastq>   output fastq
    -i <infastq>    input fastq
    -g              input is gzipped
    -n <bps>        number of basepairs to trim (default 10)
    -h              this help
    -f              the infastq is actually a folder
"""


from Bio import SeqIO
from docopt import docopt
import gzip
from tqdm import tqdm
import os
from os.path import join

if __name__ == '__main__':
    arguments = docopt(__doc__)
    input_fastq = arguments['-i']
    output_fastq = arguments['-o']
    bps = int(arguments['-n']) if arguments['-n'] else 10
    file_opener = gzip.open if arguments['-g'] else open
    if arguments['-f']:
        input_fastq = [join(input_fastq,f) for f in os.listdir(input_fastq) if "fastq" in f]
    else :
        input_fastq = [input_fastq]
        
    print "parsing the reads"
    reads = []
    for file in input_fastq:
        with file_opener(file) as handle:
            reads += [s[bps:-bps] for s in tqdm(SeqIO.parse(handle,"fastq"))]
    print "writing the reads"
    with file_opener(output_fastq,"w") as handle: SeqIO.write(reads,handle,"fastq")



