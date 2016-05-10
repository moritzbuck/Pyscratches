""" rename the IDs in a phylophlan tree


Usage:
  treenamer.py [-p <phylophlan_dir> -t <taxa_level>] -i <intree> -o <outree> 
  treenamer.py -h
  
Options:
    -o <outree>          output tree-file
    -i <intree>          input tree outputed by phylophlan
    -p <phylophlan_dir>  directory where phylphlan is
    -t <taxa_level>      taxonomic level of annotation (default: -1)
    -h                   this help
    
"""

from pandas import DataFrame
import re
from docopt import docopt

if __name__ == '__main__':
    arguments = docopt(__doc__)
    itree = arguments['-i']
    otree = arguments['-o']
    if arguments['-p']:
        taxa_file = arguments['-p'] + "/data/ppafull.tax.txt"
    else :
        taxa_file = "/home/moritz/repos/phylophlan/data/ppafull.tax.txt"

    taxa_level = int(arguments['-t']) if arguments['-t'] else -1
print arguments

taxas = DataFrame.from_csv(taxa_file, sep="\t", header=None)
taxas = {l[0] : l[1].split(".")[taxa_level].split("_")[-1] for l in taxas.itertuples()}

with open(itree) as handle:
	tree = [l for l in handle.readlines()][0]

g_keys = [k for k in taxas.keys() if k in tree]
for k in g_keys:
	tree=tree.replace(k,taxas[k])

with open(otree, "w") as handle:
	handle.write(tree)

