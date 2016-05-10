""" small script doing some stats and plots from a spades generated fasta


Usage:
  hmmparser.py [-b -C <coverages>] -o <outhead> -H <hmmerout> -G <gff> 
  hmmparser.py -h
  
Options:
    -h --help       Show this screen.
    -o <outhead>    path + name_header of output files
    -H <hmmerout>   output from hmmersearch
    -G <gff>        GFF file (from prokka) where I can find which contig each predicted AA-sequence is from
    -C <coverages>  CSV-file with contigs as lines and samples as columns (with head)
    -b              do a bin_wise table
"""

from docopt import docopt
from tqdm import tqdm
from pandas import DataFrame
import json
from pandas import Series
import sys

if __name__ == '__main__':
    arguments = docopt(__doc__)
    print arguments
    infile = arguments['-H']
    gff = arguments['-G']
    coverage_file = arguments['-C']
    outhead = arguments['-o']
    binwise = arguments['-b']

else:
    infile = "/home/moritz/b2011138/TroutBog/assemblies/assemblies/diginormalised_megahit/bins/all_bins_Pfam-A.hmm.raw"
    gff = "/home/moritz/b2011138/TroutBog/assemblies/assemblies/diginormalised_megahit/bins/all_proteins.gff"
    coverage_file = "/home/moritz/b2011138/TroutBog/assemblies/assemblies/diginormalised_megahit/bins/yes_bin_mapping.csv"
#    coverage_file = None
    outhead = "troutbog_pfam"
    binwise = True


#parse gff for gene => contig map
with open(gff) as handle:
	r_gff = [l[:-1] for l in handle.readlines() if "CDS" in l]

r_gff = [l.split("\t") for l in r_gff]
contig_dict = {[v.replace("ID=","") for v in l[8].split(";") if "ID" in v][0] : l[0] for l in r_gff}


# parse file with coverages to get a contig => coverages/sample map
if coverage_file != None:
    coverages = DataFrame.from_csv(coverage_file, sep=",")
    coverages = coverages[[c for c in coverages.columns if c not in ['length', 'GC']]]
    
if len(coverages.columns) < 1:
    print "your coverage file is not the right kind probably. it needs to be a csv, e.g. with commas not tabs"
    sys.exit()


# parse single copy pfams file 
with open("/home/moritz/repos/Pyscratches/20150211_hmmparser/sc_pfams.txt") as handle:
    sc_pfams = [l[:-1] for l in handle.readlines()]
    

# organise and parse hmmsearch output
entries = []
t = []
post = False
with open(infile) as handle:
    print "Parsing hmmsearch output!"
    for l in tqdm(handle.readlines()):
        if l[0:17] == "Domain annotation" and not post :
            post = True
        if l[0] != "#" and len(l[:-1]) > 0 and not post :
            t.append(l[:-1])
        if l == "//\n":
            entries.append(t)
            t = []
            post = False 


# make a BMFD with all entries and periferal info
hmm_dict = {}
print "Making a nice usable data structure!"
for e in tqdm(entries):
    acc = e[1].split()[1].split(".")[0]
    hmm_dict[acc] = {}
    hmm_dict[acc]['version'] = "NA"
    if "KO:" in e[1].split()[1]:
        if len(e[2].split()) > 1:
            hmm_dict[acc]['version'] = " ".join(e[2].split()[2:])
    if "PF" in e[1].split()[1]:
        hmm_dict[acc]['version'] = e[1].split()[1].split(".")[1]
        
    hmm_dict[acc]['short_name'] = e[0].split()[1]
    hmm_dict[acc]['description'] = " ".join(e[2].split()[1:])
    hmm_dict[acc]['cdss'] = []
    for l in e[7:]:
        if "inclusion" in l or "[No hits detected that satisfy reporting thresholds]" in l :
            break 
        hmm_dict[acc]['cdss'].append(l.split()[8])
    if coverage_file != None:
        coverages.loc[[contig_dict[cds] for cds in hmm_dict[acc]['cdss']]].apply(sum)
        hmm_dict[acc]['per_sample'] = dict(coverages.loc[[contig_dict[cds] for cds in hmm_dict[acc]['cdss']]].apply(sum))
        temp = coverages.loc[[contig_dict[cds] for cds in hmm_dict[acc]['cdss']]].apply(sum,1)
        temp.index =  hmm_dict[acc]['cdss']
        hmm_dict[acc]['per_cds'] = dict(temp)

# normalisation by single-copy PFAMS to get genome equivalents

if coverage_file != None:
    print "Normalisation"
    genome_equivalents = DataFrame.from_dict({p:hmm_dict[p]['per_sample'] for p in sc_pfams}).transpose().median()
#    norm_factor = DataFrame.from_dict({p:hmm_dict[p]['per_sample'] for p in sc_pfams}).transpose().mean()

    norm_factor = 1.0/genome_equivalents

    genome_equivalents.to_csv(outhead +"_genome_equivalents.csv")
    for k in tqdm(hmm_dict):
        hmm_dict[k]['normed_per_sample'] = dict(Series(hmm_dict[k]['per_sample'])*norm_factor)
        temp = (coverages.loc[[contig_dict[cds] for cds in hmm_dict[k]['cdss']]]*norm_factor).apply(sum,1)
        temp.index =  hmm_dict[k]['cdss']
        hmm_dict[k]['normed_per_cds'] = dict(temp)

        
# making pretty table and writing outputs
    normed_data = DataFrame.from_dict({k : v['normed_per_sample'] for k,v in hmm_dict.iteritems()}).transpose()
    normed_data[normed_data != normed_data]=0
    non_nulls = normed_data.sum(1) != 0
    periferal_data = DataFrame.from_dict({k : {'description' : v['description'], 'short_name' : v['short_name'] } for k,v in hmm_dict.iteritems()}).transpose()
    gene_lists = DataFrame.from_dict({k : {'hits' : ";".join(v['cdss']) } for k,v in hmm_dict.iteritems()}).transpose()
    data = periferal_data.join(normed_data).join(gene_lists).loc[non_nulls]

    if binwise:
        print "Bin-wise table making"
        bins = list(set(["_".join(k.split("_")[:2]) for k in contig_dict.keys()]))
        binwise_dict = {}
        for pfam in tqdm(hmm_dict):
            binwise_dict[pfam] = {b : sum([v for k,v in hmm_dict[pfam]['normed_per_cds'].iteritems() if b+"_" in k]) for b in bins}

        binwise = DataFrame.from_dict(binwise_dict)
        binwise = binwise.transpose().loc[ binwise.sum() != 0 ]
        binwise = binwise.fillna(0)
        periferal_data.loc[binwise.index].join(binwise).to_csv(outhead + "_binwise.csv")
        for bin in bins:
            rows = [sum([l[8].split(";"), ["contig=" + l[0]]],[]) for l in r_gff if bin +"_" in l[8]]
            dicts = [{v.split("=")[0] : v.split("=")[1] for v in r} for r in rows]
    print "Writing output..."
    data.to_csv(outhead + ".csv")

else :
    print "Writing output..."
    with open(outhead + ".csv","w") as handle:
        handle.writelines("PFAM\tshort_name\tversion\tname\tgenes\n")
        lines = [ "\t".join([k, v["short_name"], v["version"], v["description"], ";".join(sorted(v['cdss']))])+ "\n" for k,v in hmm_dict.iteritems() if len(v['cdss']) != 0]
        lines.sort()
        handle.writelines(lines)


with open(outhead + ".json","w") as handle:
    handle.write(json.dumps({k:v for k,v in hmm_dict.iteritems() if len(v['cdss']) != 0}, sort_keys=True, indent=4, separators=(',', ': ')))


    
    
