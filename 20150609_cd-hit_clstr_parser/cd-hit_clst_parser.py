""" small script doing some stats and plots from a spades generated fasta


Usage:
  cd-hit_clst_parser.py -o <outhead> -i <in_clstr> -m <mappin>
  hmmparser.py -h
  
Options:
    -o <outhead>    path + name_header of output files
    -i <in_clstr>   cd-hit clstr file
    -m <mappin>     mapping file (tsv format)
"""

from docopt import docopt
from tqdm import tqdm
from pandas import DataFrame
import json
from pandas import Series

if __name__ == '__main__':
    arguments = docopt(__doc__)
    print arguments
    infile = arguments['-i']
    outhead = arguments['-o']
    map_file = arguments['-m']

else:
    infile = "/home/moritz/b2011138/TroutBogHypo/assemblies/assemblies/megahit/bins/cdhit_0.95_clusters.txt.clstr"
    outhead = "troubog_clstrs"
    map_file = "/home/moritz/b2011138/TroutBogHypo/assemblies/assemblies/megahit/bins/mapping_all_cdss.ffn/mapping.tsv"

with open(infile) as handle:
    lines = handle.readlines()

clusters = {}
name = lines[0][1:-1]
cdss = {}

for i in tqdm(lines[1:]):
    if i[0]==">":
        name = i[1:-1]
        clusters[name] = cdss
        cdss = {}
    else:
        temp = i.split()
        length = int(temp[1][:-3])
        gene = temp[2][1:-3]
        identity = None if len(temp) == 4 else float(temp[-1][2:-1])/100
        cdss[gene] = (length, identity) 


non_singletons = {k:v for k,v in clusters.iteritems() if len(v.values()) > 1}
bins = [b for b in bins if len(b) > 1]
bins = [list(b) for b in bins]
bin_list = set(sum(bins,[]))
bin_cooc = {b: sum([c for c in bins if b in c],[]) for b in tqdm(bin_list)}
cutoff = 0
bin_cooc_filter = {k: {b: v.count(b) for b in set(v) if v.count(b) > cutoff and b != k} for k,v in tqdm(bin_cooc.iteritems())}
df = DataFrame.from_dict(bin_cooc_filter).fillna(0)


mapping = DataFrame.from_csv(map_file, sep="\t")
cds2cov = mapping[[m for m in mapping.columns if "IH" in m]].apply(sum,1).to_dict()
