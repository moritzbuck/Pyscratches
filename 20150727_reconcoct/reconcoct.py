""" rerun concoct on a bin 
Usage:
  reconcoct.py [-m <min_size> -p <prefix>] <mapping> <fasta>
  reconcoct.py -h
  
Options:
   -h --help       Show this screen.
   -m <min_size>   min length of accepted subbin [default 0]
   -p <prefix>     prefix of subbin name [default bin_]
"""
from docopt import docopt
from Bio import SeqIO
import os
import sh
from pandas import DataFrame
from generic.pre_and_post import make_bins

    
if __name__ == '__main__':
    arguments = docopt(__doc__)
    print(arguments)
    fasta = arguments['<fasta>']
    print("opening mappin")
    mapping = DataFrame.from_csv(arguments['<mapping>'])
    min_s = int(arguments['-m']) if arguments['-m'] else 0
    pref = arguments['-p'] if arguments['-p'] else "bin_"
    with open(fasta) as handle:
        ids = [i.id for i in SeqIO.parse(handle, "fasta")]
    print "filter mapping"
    mapping = mapping.loc[ids][[l for l in mapping if l not in ["length", "GC"]]]

    ### removing samples where less than 95% of the contig have reads
    mapping = mapping[[c for c in mapping if float(sum(mapping[c] > 0 ))/len(mapping[c])>0.95]]
    mapping.to_csv("submapping.csv", sep="\t")

    print "actual binning"
   sh.concoct( "--coverage_file", "submapping.csv", "--composition_file", fasta,  "-k",  "4", "-b", "reconcoct/", "-l", "999")
    print "making fastas"
    DataFrame.from_dict(make_bins("./", "reconcoct/clustering_gt999.csv", fasta, "submapping.csv", prefix= pref, min_size=min_s)[1], orient = "index").to_csv("reconc_namemap.csv")
    

