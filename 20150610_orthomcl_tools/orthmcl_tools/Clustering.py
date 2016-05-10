
import os
import sh
import json
from Bio import SeqIO
import Bio
from collections import Counter
from tqdm import tqdm
import sys
import shutil
from pandas import DataFrame
from numpy import nan_to_num, prod
from orthmcl_tools.orthoMCL import orthoMCL;
import dendropy
import re
import operator
import pandas
from numpy import power

def json2nexus(filename, out_name, trait_file, updated_json = None):
    with open(filename) as handle:
        data=json.loads(handle.read())

    genomes = list(set(sum([d['genes'].keys() for d in data],[])))

    trait_vectors = {g: [len(c['genes'][g]) if( c['genes'].has_key(g) ) else 0 for c in data] for g in genomes }
    
    formated_trait_vector = {g : "".join([hex(15 if t >15 else t)[-1].capitalize() for t in ts]) for g,ts in trait_vectors.iteritems()}

    with open(out_name, "w") as handle:
        for g, ts in formated_trait_vector.iteritems():
            handle.write(">" + g + "\n")
            handle.write(ts + "\n")

    zf = len(str(len(data)))
    with open(trait_file,"w") as handle:
        handle.write(
"""#NEXUS
BEGIN CHARACTERS;
DIMENSIONS NCHAR=%d;
CHARSTATELABELS""" %(len(data)))
        for i,d in enumerate(data):
            handle.write("\n%d %s," %(i, "orthoMCL_"+str(i).zfill(zf)))
            d['name'] = "orthoMCL_"+str(i+1).zfill(zf)
        handle.write("\b;\nEND;\n")

    if updated_json:
        with open(updated_json,"w") as handle:
            json.dump(data, handle,  indent=4, sort_keys=True)

def mergeJSONwPAUP(json_file, paup_file):
    with open(json_file) as handle:
        data=json.loads(handle.read())
    changes = DataFrame.from_csv(paup_file, sep="\t", index_col=None, header=0)
    for d in tqdm(data):
        if len(d['genes']) > 1:
            sub = changes[changes['orthoMCL'] == (int(d['name'].split("_")[1])+1)]
            d['ancestral'] = {r[1] : {'from' : r[3], 'to': r[4]} for r in sub.itertuples()}
    return data

def make_bmft(data, genomes_of_interest, ancestry, clustering):
    data_sub = [d for d in data if any([g in genomes_of_interest for g in d['genes'].keys()] ) ]
    head = ["ID", "hypothetical.name", "name.confidence"]+  genomes_of_interest+ [ "ancestor" , ancestry , "genomes.of.interest", "other.genomes", "genes"]
    lines = []
    genomes = set(clustering.gene2genome.values())
    for d in data_sub:
        name = d['name'] 
        hname = d['annotation']
        c_name = d['annot_fraction']
        if d.has_key('ancestral') and d['ancestral'].has_key(ancestry):
            ancestor = d['ancestral'][ancestry]
            frmo = ancestor['from']
            diff = ancestor['to']-ancestor['from']
        else :
            frmo = "NA"
            diff = 0
        counts = [len(d['genes'][g]) if d['genes'].has_key(g) else 0 for g in genomes_of_interest]
        icounts = len([g for g in genomes if g in genomes_of_interest and d['genes'].has_key(g)])
        ocounts = len([g for g in genomes if g not in genomes_of_interest and d['genes'].has_key(g)])
        genes = ";".join(d['mapping'].keys())
        line = [ name, hname, c_name ] + counts + [ frmo, diff , icounts, ocounts, genes ] 
        lines +=  [[str(l) for l in line]]
    bmft = DataFrame.from_records(lines, columns=head, index="ID")
    return bmft
   
        
class GeneCluster(object):

    def __repr__(self): return '<%s object %s, annotated as  %s with %i genes from %i genomes>' % (self.__class__.__name__, self.name, self.annotation, len(self.genes), len(self.genomes))

    def __init__(self, clustering , genes, name = None ):

        self.clustering = clustering

        if type(genes) == dict :
            self.from_dict(genes)
        else:
            self.name = name
            self.from_list(genes)
           
    def to_dict(self):
        return {u'name': self.name, u'annot_fraction': self.annot_fraction, u'annotation': self.annotation,  u'genes': self.genome_2_gene_map,  u'mapping': self.mapping}

    def from_dict(self,imp):
        self.name = imp['name']
        self.genes = imp['mapping'].keys()
        self.genomes = imp['genes'].keys()
        self.genome_2_gene_map =  imp['genes']
        
        self.annotation = imp['annotation']
        self.annot_fraction = imp['annot_fraction']
        self.mapping = imp['mapping']
        
    def from_list(self, genes):
        self.genomes = list(set([g.split("|")[0] for g in genes]))
        self.genes = [g.split("|")[1] for g in genes]
        self.genome_2_gene_map = {go : [ge.split("|")[1] for ge in genes if go == ge.split("|")[0]] for go in self.genomes}
       
        sub_dict = {g : self.clustering.id2name_map[g] for g in self.genes}
        name_counts = Counter(sub_dict.values())
        total = sum([name_counts[z] for z in name_counts])
        annot_frac = float(name_counts.most_common()[0][1])/float(total)

        self.annotation = name_counts.most_common(1)[0][0]
        self.annot_fraction = annot_frac
        self.mapping = sub_dict

    def to_sequences(self, short=False):
        seqs = []
        for g in self.genomes:
            genome = [f for f in self.clustering.genomes if f.split("/")[-1].split(".")[0] == g ][0]
            with open(genome, "r") as handle:
                seqs += [s for s in SeqIO.parse(handle, "fasta") if s.id in self.genes]
                if short:
                    for s in seqs:
                        s.description = ""
        return seqs

    def calc_checksum(self, s):
        return str(sum(ord(c)*i for i,c in enumerate(s)))

    def align(self,output_file):
        with open("temp.faa","w") as unalign:
            temp_seqs = self.to_sequences(short=True)
            SeqIO.write(temp_seqs, unalign, "fasta")
        sh.muscle("-in", "temp.faa","-out", "temp_aligned.faa")
        os.remove("temp.faa")
        try:
            if self.clustering.seq_type == "genes":
                sh.Gblocks("temp_aligned.faa", "-tD")
            else:
                sh.Gblocks("temp_aligned.faa", "-t=p", "-b5=h", "-b4=2", "-b2=0", "-b3=2000")
        except:
            pass
        
        os.remove("temp_aligned.faa")
        if os.path.exists("temp_aligned.faa-gb.htm"):
            os.remove("temp_aligned.faa-gb.htm")
            shutil.move("temp_aligned.faa-gb", output_file)
            return 1
        else : 
            return 0

    def tree_construction(self,alignment, outputtree, sccs = False):
        if sccs:
            with open(alignment) as handle:
                seqs = [s for s in SeqIO.parse(handle, "fasta")]
            for s in seqs:
                s.id = "_".join(s.id.split("_")[:-1])
                s.description = ""
            with open("tempblabla.faa","w") as handle:
                SeqIO.write(seqs,handle,"fasta")

        sh.FastTree("-out", outputtree, "tempblabla.faa" if sccs else alignment)
        if sccs:
            sh.rm("tempblabla.faa")
            
        
    def core_probability(self):
        present = prod([self.clustering.completnesses[g] for g in self.genomes])
        abscent = prod([1-v for k,v in self.clustering.completnesses.iteritems() if k not in self.genomes])
        return present*abscent

    def non_core_probability(self):
        prob_of_random_pres = lambda g: 1.0 - power(float(len(self.clustering)-1)/len(self.clustering) , self.clustering.genome2len[g])
        prob_of_random_absc = lambda g: power(float(len(self.clustering)-1)/len(self.clustering) , self.clustering.genome2len[g])
        present = prod([prob_of_random_pres(g) for g in self.genomes])
        abscent = prod([prob_of_random_absc(k) for k in self.clustering.completnesses.keys() if k not in self.genomes])
        return present*abscent

    def non_core_probability_plural(self):
        prob_of_random_pres = lambda g: 1.0 - power(float(sum(self.clustering.genome2len.values())-len(self.genes))/sum(self.clustering.genome2len.values()) , self.clustering.genome2len[g])
        prob_of_random_absc = lambda g: power(float(sum(self.clustering.genome2len.values())-len(self.genes))/sum(self.clustering.genome2len.values()) , self.clustering.genome2len[g])
        present = prod([prob_of_random_pres(g) for g in self.genomes])
        abscent = prod([prob_of_random_absc(k) for k in self.clustering.completnesses.keys() if k not in self.genomes])
        return present*abscent
        
    
class Clustering(object):

    def __repr__(self): return '<%s object "%s with %i clusters">' % (self.__class__.__name__, self.name, len(self))
    def __iter__(self): return iter(self.clusters)
    def __len__(self): return len(self.clusters)
    def __getitem__(self, key): return self.clusters[key]

        
    def __init__(self,proteoms,  out_path, name, mcl, gff = None, seq_type="proteins", checkm = None):

        self.genomes = proteoms
        self.seq_type = seq_type
        self.path = out_path
        self.oMCL_path = mcl.out_dir
        if not os.path.exists(out_path):
            os.makedirs(out_path)

        if checkm:
            self.checkm = pandas.read_table(checkm, skiprows = 3, index_col = 0, sep = r"\s*", names = ["genome", "lineage","nb_gen", "nb_markers", "nb_sets", "0", "1","2","3","4","5+", "completness", "contamination", "heterogeneity"], comment = '-')
            self.completnesses = (self.checkm['completness']/100.0).to_dict()
            self.completnesses = {unicode(k): v if v < 0.95 else 0.95 for k,v in self.completnesses.iteritems()}
        self.seed = 23
        self.name = name
        self.base = self.path + self.name
        self.scg_tree =  self.base +"_nodes_labeled.tree"
        
        self.processed_clusters = self.base + ".json"
        self.align_path = self.base + "_align/"
        self.scc_align_path = self.base + "_scc_align/"
        self.db = self.oMCL_path + "goodProteins.fasta"
        if os.path.exists(self.processed_clusters):
            with open(self.processed_clusters) as file:
                self.clusters= [GeneCluster(self,c) for c in json.load(file)]
        else :
            self.clusters = []

        if gff:
            with open(gff,"r") as handle:
                temp = [ {v.split("=")[0] : v.split("=")[1] for v in l.split("\t")[-1].split(";")} for l in handle.readlines() ]
                self.id2name_map = {cds['locus_tag'] : cds['product'] for cds in temp if cds.has_key('product')}
        else :
            self.id2name_map = {}
            self.gene2genome = {}
            self.genome2len = {}
            for g in self.genomes:
                with open(g) as handle:
                    temp = {s.id : " ".join(s.description.split()[1:]) for s in SeqIO.parse(handle,"fasta")}
                    self.genome2len[g.split("/")[-1].split(".")[0]] = len(temp.keys())
                    self.gene2genome.update({gene : ".".join(g.split(".")[:-1]).split("/")[-1] for gene in temp.keys()})
                    self.id2name_map.update(temp)

                    
        self.orthoMCL = mcl #orthoMCL(self.oMCL_path, self.genomes, self.name)
        self.raw_clusters = self.orthoMCL.out_mcl
        self.anc_rec_path = out_path + "AncRec/"
    
    def single_copy_clusters(self):
        return [c for c in self.clusters if len(c.genes) == len(self.genomes) and len(c.genomes)==len(self.genomes)]

    def almost_single_copy_clusters(self):
        return [c for c in self.clusters if len(c.genes) == (len(self.assemblies)-1) and (len(c.genomes)==len(self.assemblies)-1)]
    

    def post_process(self):
        print "Post processing cluster:"
        with open(self.raw_clusters) as c_file:
            self.clusters=[GeneCluster(self, name=l[:-1].split(": ")[0], genes = l[:-1].split(": ")[1].split()) for l in tqdm(c_file.readlines())]
            
        print "Post processing single genes:"
        non_singletons = set(sum([c.genes for c in self.clusters],[]))
        for i in tqdm(self.id2name_map.keys()):
            if i not in non_singletons:
                self.clusters += [GeneCluster(self, name = i,  genes =  [self.gene2genome[i] + "|" + i ])]
            
        with open(self.processed_clusters, 'w') as outfile:
            json.dump([c.to_dict() for c in self.clusters], outfile,  indent=4, sort_keys=True)

    def run(self):
#        self.orthoMCL.full_pipe()
        print "post-process"
        self.post_process()
        
    def align_all(self):
        print "aligning the hell out of it"
        if os.path.exists(self.align_path):
             shutil.rmtree(self.align_path)
        os.makedirs(self.align_path)

        for i,clust in tqdm(enumerate(self.single_copy_clusters())):
            clust.align(self.align_path + str(i) + ".fasta")

        
    def cat_align(self):
        print "CONCATENATE!!!11"
        cat_seqs = None
        for f in os.listdir(self.align_path):
            if "fasta" in f:
                with open(self.align_path + f,"r") as file:
                    seqs = [s for s in SeqIO.parse(file,"fasta")]
                    if not cat_seqs:
                        cat_seqs = seqs
                        order = [s.id.split("|")[0] for s in cat_seqs]
                        for i,s in enumerate(cat_seqs):
                            s.id = order[i]
                    else :
                        seqs_sorted = [[z for z in seqs if o in z.id][0] for o in order]
                        cat_seqs = [ s + seqs_sorted[i] for i,s in enumerate(cat_seqs) ]

        for i,s in enumerate(cat_seqs):
            s.id = order[i]
            s.description = "composite of " + str(len([f for f in os.listdir(self.align_path) if "fasta" in f])) + " single copy genes-clusters"
        with open(self.base + "_cat_align.fasta","w") as outp:
            SeqIO.write(cat_seqs,outp,"fasta")
        return cat_seqs
        
    def tree_construction(self,root = None, sccs = False):
        threads = 16 
        print "build a tree"
        if os.path.exists(self.base + "RAxML/" ):
            sh.rm("-r", self.base + "RAxML/")
        os.makedirs(self.base + "RAxML/")

        if self.seq_type == "proteins" :
            model = "PROTGAMMALG"
        else:
            model = "GTRGAMMA"

        alignment = self.base + "_scc_cat_align.fasta" if sccs else self.base + "_cat_align.fasta"
        
        sh.raxmlHPC_PTHREADS_AVX("-w", self.base + "RAxML/", "-T", threads-2, "-m", model, "-p", self.seed, "-#", 20, "-s", alignment, "-n", "T13", "-o", root) 
        print "boostrap dat tree"
        sh.raxmlHPC_PTHREADS_AVX("-w", self.base + "RAxML/", "-T", threads-2, "-m", model, "-p", self.seed, "-b", self.seed, "-#", 100, "-s", alignment, "-n", "T14", "-o", root)
        print "combine"
        sh.raxmlHPC_AVX("-m", "GTRCAT", "-w", self.base + "RAxML/", "-p", self.seed, "-f", "b", "-t", self.base + "RAxML/"+"RAxML_bestTree.T13", "-z",self.base + "RAxML/"+ "RAxML_bootstrap.T14", "-n", "T15", "-o", root)
        print "clean up"
        if os.path.exists(self.base + "_branches_labeled.tree"):
            os.remove(self.base + "_branches_labeled.tree")
            os.remove(self.base + "_nodes_labeled.tree")
        sh.ln("-s",  self.base + "RAxML/RAxML_bipartitionsBranchLabels.T15", self.base +"_branches_labeled.tree")
        sh.ln("-s",  self.base + "RAxML/RAxML_bipartitions.T15", self.scg_tree)

        
    def rm_genome(self, name):
        self.assemblies = [a for a in self.assemblies if name not in a.name]
        for c in self:                                                             
            c.genomes = [g for g in  c.genomes if name not in g]
        for c in self:                                                             
            c.genes = [g for g in  c.genes if name not in g]


    def keep_genomes(self, genomes):
        self.assemblies = [a for a in self.assemblies if a.name in genomes]
        for c in self:                                                             
            c.genomes = [g for g in  c.genomes if g in genomes ]
        for c in self:                                                             
            c.genes = [g for g in  c.genes if len([o for o in genomes if g.count(o) ==1]) > 0 ]


    def cooccurence_matrix(self):
        matrix = DataFrame(data=0, index = [a.name for a in self.assemblies], columns=[a.name for a in self.assemblies])
        for a in self.assemblies:
            for b in self.assemblies:
                count = 0
                for c in self:
                    if a.name in c.genomes and b.name in c.genomes:
                        count += 1
                matrix[a.name][b.name] = count
                
        return matrix

    def make_cluster_bmft(self):
        cluster_table = DataFrame.from_dict({i : {k : len(v) for k,v in c.to_dict()['genes'].iteritems()} for i,c in enumerate(self)}, orient='index')
        cluster_table = cluster_table.apply(nan_to_num)
        cluster_table['annotations'] = [c.annotation for c in self]
        cluster_table['qual_annot'] = [c.annot_fraction for c in self]
        cluster_table['genes'] = [";".join(c.genes) for c in self]
        return cluster_table
        

    def ancestral_reconstr(self, outgroup = None):
        print "Make trait File"
        if not os.path.exists(self.anc_rec_path):
            os.makedirs(self.anc_rec_path)
        json2nexus(self.processed_clusters, self.anc_rec_path + "Trait_genome.fasta", self.anc_rec_path + "Traits.nex"  )
        with open(self.anc_rec_path + "Trait_genome.fasta") as handle:
            seqs=[s for s in SeqIO.parse(handle,"fasta") ]
        for s in seqs:
            s.seq.alphabet = Bio.Alphabet.DNAAlphabet()
        with open(self.anc_rec_path + "Matrix.nex","w") as handle:
            SeqIO.write(seqs,handle,"nexus")
        os.system("sed -i 's/format datatype=dna missing=? gap=- interleave;/format interleave datatype=standard   gap=- symbols=\"0123456789ABCDEF\";/' " + self.anc_rec_path + "Matrix.nex")

        print "Make tree file"
        new = dendropy.Tree.get(path=self.scg_tree, schema="newick")
#        new.write(path = self.anc_rec_path + "Tree.nex", schema="nexus")
        tree_string = re.sub("[0-9]{2,}","", re.sub(r':[0123456789.]+','', str(new)))
        taxons = {x : i for i,x in enumerate([t.taxon.label for t in new.nodes() if t.is_leaf()])}
        for k in taxons.keys():
            tree_string=re.sub(k+"([),])",str(taxons[k]+1)+r"\1", tree_string)

        taxons_t = sorted(taxons.items(), key=operator.itemgetter(1))

        t_table = "\n".join([("\t\t" + str(p[1]+1) + "\t" + p[0].replace(" ","_") +",")  for p in taxons_t])
        paup_tree = """#NEXUS
Set AllowPunct=Yes;
BEGIN TREES;
\tTRANSLATE
%s\b
;
TREE Bioperl_1 = [&U] %s ;
END; 
""" % (t_table , tree_string)
            
        paup_script = """
#NEXUS
BEGIN PAUP;
      log file=%sPaup.log replace=yes start=yes;
      exe %sMatrix.nex;
      exe %sAssumptions.nex;
      exe %sTree.nex;
      ctype GeneCopy:All;
      pset Opt=AccTran;

      taxset out= %d;
      outgroup %s;
      set root=outgroup outroot=monophyl;
      root root=outgroup;      

      describetrees 1/Xout=both ApoList=Yes ChgList=yes; [brLens=yes;]
      savetrees from=1 to=1 format=nexus root=yes file=%sOutTree.rtree;
      q;
END;
""" % (self.anc_rec_path, self.anc_rec_path, self.anc_rec_path, self.anc_rec_path, len(seqs), outgroup if outgroup else "out",  self.anc_rec_path)

        assumptions="""
#NEXUS

BEGIN ASSUMPTIONS;
USERTYPE GeneCopy (STEPMATRIX) = 16
      0   1   2   3   4   5   6   7   8   9   A   B   C   D   E   F
[0]   .   10.0 11.0 11.2 11.4 11.6 11.8 12.0 12.2 12.4 12.6 12.8 13.0 13.2 13.4 13.6
[1]   5.0 .   1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4 3.6
[2]   5.2 0.2 .   0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6
[3]   5.4 0.4 0.2 .   0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4
[4]   5.6 0.6 0.4 0.2 .   0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2
[5]   5.8 0.8 0.6 0.4 0.2 .   0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0
[6]   6.0 1.0 0.8 0.6 0.4 0.2 .   0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8
[7]   6.2 1.2 1.0 0.8 0.6 0.4 0.2 .   0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6
[8]   6.4 1.4 1.2 1.0 0.8 0.6 0.4 0.2 .   0.2 0.4 0.6 0.8 1.0 1.2 1.4
[9]   6.6 1.6 1.4 1.2 1.0 0.8 0.6 0.4 0.2 .   0.2 0.4 0.6 0.8 1.0 1.2
[A]   6.8 1.8 1.6 1.4 1.2 1.0 0.8 0.6 0.4 0.2 .   0.2 0.4 0.6 0.8 1.0
[B]   7.0 2.0 1.8 1.6 1.4 1.2 1.0 0.8 0.6 0.4 0.2 .   0.2 0.4 0.6 0.8
[C]   7.2 2.2 2.0 1.8 1.6 1.4 1.2 1.0 0.8 0.6 0.4 0.2 .   0.2 0.4 0.6
[D]   7.4 2.4 2.2 2.0 1.8 1.6 1.4 1.2 1.0 0.8 0.6 0.4 0.2 .   0.2 0.4
[E]   7.6 2.6 2.4 2.2 2.0 1.8 1.6 1.4 1.2 1.0 0.8 0.6 0.4 0.2 .   0.2
[F]   7.8 2.8 2.6 2.4 2.2 2.0 1.8 1.6 1.4 1.2 1.0 0.8 0.6 0.4 0.2 .
;
charset All =1-.;
END;
"""
        print "Write nexus files"
        with open(self.anc_rec_path + "Script.nex", "w") as handle:
            handle.writelines(paup_script)

        with open(self.anc_rec_path + "Assumptions.nex", "w") as handle:
            handle.writelines(assumptions)

        with open(self.anc_rec_path + "Tree.nex", "w") as handle:
            handle.writelines(paup_tree)

        print "Run and parse paup"
        os.system("paup4a146_centos64 -f -n " + self.anc_rec_path + "Script.nex > " +  self.anc_rec_path + "OutPut.nex")
        sh.perl("/pica/h1/moritz/repos/Pyscratches/20150610_orthomcl_tools/orthmcl_tools/parsePaupLog.pl", self.anc_rec_path + "OutPut.nex")
        os.remove(self.anc_rec_path + "OutPut.nex.nodes.johan")
        print "Merge data with clusters"
        merged_json = mergeJSONwPAUP(self.processed_clusters, self.anc_rec_path + "OutPut.nex.changes")
        with open(self.anc_rec_path + "clusters.json","w") as handle:
            json.dump(merged_json, handle,  indent=4, sort_keys=True)

        return merged_json
    
        
    def sccs_cat_align(self):
        all_sccs = {}
        for a in self.assemblies:
            with open(a.sccs_genes) as fasta:
                all_sccs[a.name] = [s for s in SeqIO.parse(fasta,"fasta")]
        cogs = list(set(sum([[s.id for s in seqs] for seqs in all_sccs.values()],[])))
        cogs.sort()

        a2id = {a.name : self.assembly_id(a.name + a.name)  for a in self.assemblies}
        id2a = {a2id[a] : a for a in a2id}
        
        for c in cogs:
            c_seqs = []
            for a,seqs in all_sccs.iteritems():
                seq = [s for s in seqs if s.id == c]
                if len(seq) == 1:
                    seq = seq[0]
                else:
                    continue
                seq.id = a2id[a]
                seq.description = ""
                c_seqs.append(seq)
            if len(c_seqs) == len(self.assemblies):
                with open("temp.faa","w") as fasta:
                    SeqIO.write(c_seqs,fasta,"fasta")
                print "aligning genes for",c
                sh.muscle("-in", "temp.faa","-out", "temp_aligned.faa")
                os.remove("temp.faa")
                try:
                    sh.Gblocks("temp_aligned.faa", "-tD")
                except:
                    pass
        
                os.remove("temp_aligned.faa")
                if os.path.exists("temp_aligned.faa-gb.htm"):
                    os.remove("temp_aligned.faa-gb.htm")
                    if not os.path.exists(self.scc_align_path):
                        os.makedirs(self.scc_align_path)
                    shutil.move("temp_aligned.faa-gb", self.scc_align_path + c + ".ffn")


        print "CONCATENATE!!!11"
        cat_seqs = None
        for f in os.listdir(self.scc_align_path):
            if "ffn" in f:
                with open(self.scc_align_path + f,"r") as file:
                    seqs = [s for s in SeqIO.parse(file,"fasta")]
                    if not cat_seqs:
                        cat_seqs = seqs
                        order =  [s.id for s in cat_seqs]
                        for i,s in enumerate(cat_seqs):
                            s.id = order[i]
                    else :
                        seqs_sorted = [[z for z in seqs if o in z.id][0] for o in order]
                        cat_seqs = [ s + seqs_sorted[i] for i,s in enumerate(cat_seqs) ]
            for i,s in enumerate(cat_seqs):
                s.id = id2a[order[i]]
                s.description = "composite of " + str(len([f for f in os.listdir(self.scc_align_path) if "ffn" in f])) + " single copy cogs"
            with open(self.base + "_scc_cat_align.fasta","w") as outp:
                SeqIO.write(cat_seqs,outp,"fasta")
                
        return all_sccs



def trash(data):
    lines = []
    data_sub = [d for d in data if any([g in genomes_of_interest for g in d['genes'].keys()] ) ]
    head = ["ID", "hypothetical.name", "name.confidence"]+  genomes_of_interest+ [ "ancestor" , ancestry , "genomes.of.interest", "other.genomes", "genes"]
    for d in data_sub:
        name = d['name'] 
        hname = d['annotation']
        c_name = d['annot_fraction']
        if d.has_key('ancestral') and d['ancestral'].has_key("node_55_to_node_54"):
            ancestor = d['ancestral']["node_55_to_node_54"]
            frmo = ancestor['from']
            diff = ancestor['to']-ancestor['from']
        else :
            frmo = "NA"
            diff = 0
        if d.has_key('ancestral') and d['ancestral'].has_key("node_54_to_node_53"):
            ancestor = d['ancestral']["node_54_to_node_53"]
            frmo2 = ancestor['from']
            diff2 = ancestor['to']-ancestor['from']
        else :
            frmo2 = "NA"
            diff2 = 0
    
        counts = [len(d['genes'][g]) if d['genes'].has_key(g) else 0 for g in genomes_of_interest]
        icounts = len([g for g in genomes if g in genomes_of_interest and d['genes'].has_key(g)])
        ocounts = len([g for g in genomes if g not in genomes_of_interest and d['genes'].has_key(g)])
        genes = ";".join(d['mapping'].keys())
        line = [ name, hname, c_name ] + counts + [ frmo, diff, frmo2, diff2, icounts, ocounts, genes ] 
        lines +=  [[str(l) for l in line]]
