from Bio import SeqIO
from pandas import DataFrame
from pandas import Series
import pandas
from tqdm import tqdm
import sh
import os
from Bio import SeqIO
from pandas import DataFrame
import json
from collections import Counter
from numpy import mean
from math import isnan

class GeneCluster(object):

    def __repr__(self): return '<%s object %s with %i genes from %i genomes>' % (self.__class__.__name__, self.annotation, len(self.genes), len(self.genomes))

    def __init__(self, id2name_map, genes ):
        self.id2name_map = id2name_map
        
        if type(genes) == dict :
            self.from_dict(genes)
        else:
            self.from_list(genes)
           
    def to_dict(self):
        return {u'annot_fraction': self.annot_fraction, u'annotation': self.annotation,  u'genes': self.genome_2_gene_map,  u'mapping': self.mapping}

    def from_dict(self,imp):
        self.genes = imp['mapping'].keys()
        self.genomes = imp['genes'].keys()
        self.genome_2_gene_map =  imp['genes']
        
        self.annotation = imp['annotation']
        self.annot_fraction = imp['annot_fraction']
        self.mapping = imp['mapping']
        
    def from_list(self, genes):
        self.genes = genes
        self.genomes = list(set(["_".join(g.split('_')[:-1]) for g in genes]))
        self.genome_2_gene_map =  {go : [ge for ge in genes if go in ge] for go in self.genomes}
        
        sub_dict = {g : self.id2name_map[g] for g in self.genes}
        name_counts = Counter(sub_dict.values())
        total = sum([name_counts[z] for z in name_counts])
        annot_frac = float(name_counts.most_common()[0][1])/float(total)

        self.annotation = name_counts.most_common(1)[0][0]
        self.annot_fraction = annot_frac
        self.mapping = sub_dict



blast_head = ["query","subject","identity","length","mismatches", "gapopenings", "querystart","queryend","subjectstart","subjectend","Evalue","bitscore"]

def make_db(path , seq_type = "proteins" ):
    #seq_type can be "genes", "genomes", or "proteins"         
    sh.makeblastdb("-in", path, "-dbtype" , "prot" if seq_type == "proteins" else "nucl")

########################################################################################################3
    
def blast(query,output_path, db, alg = "blastp", outfmt = 6, evalue = 0.1):
    itype = "gene" if alg == "blastn" or alg == "tblastx" else "protein"
    otype = "protein" if alg == "blastp" or alg == "tblastn" else "gene"
    print "Blasting", query, " to the", otype, "db"
    blasty = sh.Command(alg)

    blasty("-evalue", evalue, "-db", db, "-query",query, "-outfmt", outfmt, "-num_threads", 1, _out =  output_path + ".raw.blast")

    raw_blast = DataFrame.from_csv(output_path+".raw.blast",header=-1,index_col=[0,1],sep="\t")
    raw_blast.columns = blast_head[2:12]
    raw_blast.index.names = blast_head[0:2]
    processed_blast = raw_blast
        
    return processed_blast

                
########################################################################################################3
genes = "/home/moritz/people/valerie/OilSands/syntrophacae/all_syntrophos.faa"
pat = "/home/moritz/people/valerie/OilSands/assemblies/merged/temp/"
maps = [ "map_run003-pool01-smds.coverage.orfs" , "map_run003-pool03-smds.coverage.orfs" , "map_run003-pool05-smds.coverage.orfs" , "map_run003-pool07-smds.coverage.orfs" , "map_run003-pool09-smds.coverage.orfs" , "map_run003-pool02-smds.coverage.orfs" , "map_run003-pool04-smds.coverage.orfs" , "map_run003-pool06-smds.coverage.orfs" , "map_run003-pool08-smds.coverage.orfs" , "map_run003-pool10-smds.coverage.orfs"]

maps = { m.split("-smds")[0] : pat + m for m in maps}


covs_dict = {}
feat_lens = {}
for k,v in maps.iteritems():
    with open(v) as handle:
        print "Sample:", k
        lines = handle.readlines()
        covs_dict[k] = {}
        for l in tqdm(lines):
            temp = l.split("\t")
            if temp[0] != "all":
                orf = temp[8].split(";")[0].split("=")[1]
                cov = int(temp[9])
                width = int(temp[10])
                orf_len = int(temp[11])
                if not covs_dict[k].has_key(orf):
                    covs_dict[k][orf]=0
                covs_dict[k][orf] += cov*width
                feat_lens[orf]=orf_len
                
covs = DataFrame.from_dict(covs_dict)

## print "blasting"
all_syntrophos = genes 
make_db(all_syntrophos)

blast(all_syntrophos, all_syntrophos, db = all_syntrophos , alg = "blastp", outfmt=6)

raw_blast = DataFrame.from_csv(all_syntrophos + ".raw.blast", header=-1,index_col=[0,1],sep="\t")
raw_blast.columns = blast_head[2:12]
raw_blast.index.names = blast_head[0:2]

with open(all_syntrophos) as file:
    dicti = {s.id : len(s) for s in SeqIO.parse(file,"fasta")}

qcov = 0.40
identity = 30


blast_covs = (raw_blast['subjectend']-raw_blast['subjectstart']+1)/ [float(dicti[n]) for n in raw_blast.index.get_level_values(0)]
# -raw_blast['mismatches']- raw_blast['gapopenings']
filtered_blast = raw_blast[blast_covs > qcov]
filtered_blast = filtered_blast[filtered_blast['identity'] > identity]
print "MCL ing"

filtered_blast.to_csv(genes +".blast",sep="\t", header=False)
resource = 4
inflation = 1.1
sh.cut("-f", "1,2,11", genes +".blast", _out = genes +".abc")
sh.mcxload("-abc", genes +".abc", "--stream-mirror", "--stream-neg-log10", "-stream-tf", "ceil(200)", "-o",  genes +".mci", "-write-tab"  , genes +".tab")
sh.mcl(genes +".mci",  "-use-tab",  genes +".tab", "-o",  genes +".raw.clust")

with open(genes +".raw.clust") as c_file:
    clusters=[l[:-1].split("\t") for l in c_file.readlines()]
    clusters = [[g.split("(")[0] for g in v ] for v in clusters]
clusters_p = []

with open(all_syntrophos,"r") as fas:
    id2name_map = {s.description.split(" ")[0] : " ".join(s.description.split(" ")[1:]) for s in SeqIO.parse(fas, "fasta")}
                    
for i,c in enumerate(clusters):
    clusters_p.append(GeneCluster(id2name_map , c))
            
with open(genes + ".proc.clust", 'w') as outfile:
    json.dump([c.to_dict() for c in clusters_p], outfile,  indent=4, sort_keys=True)


c_to_g = {i : c.genes for i, c in enumerate(clusters_p)}
g_to_c = {}
for k,v in c_to_g.iteritems():
    for g in v:
        if g in covs.index:
            g_to_c[g] = k

covs['cluster'] = Series(g_to_c)      
covs['length'] = Series(feat_lens)

SCCs = sum([c.genes for c in clusters_p if len([g for g in c.genomes if "bin_" not in g]) == 7 and len([g for g in c.genes if "bin_" not in g]) == len([g for g in c.genomes if "bin_" not in g])],[])

 covs['single'] = Series({g : True for g in SCCs })

covs.to_csv("/home/moritz/people/valerie/OilSands/gene_wise_coverages.csv")
