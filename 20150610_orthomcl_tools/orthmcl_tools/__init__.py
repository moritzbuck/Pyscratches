import os
import sh
from orthmcl_tools.Clustering import Clustering
from orthmcl_tools.orthoMCL import orthoMCL
from tqdm import tqdm
from numpy import log

if __name__ == '__main__':
    prots = [ l[:-1] for l in sh.find("/home/moritz/people/moritz/bacteriodetes/") if ".faa" in l  and not "checkm" in l]
    mcl = orthoMCL("/home/moritz/people/moritz/bacteriodetes/clustering/", prots,"bacteroidetes")
    dat = Clustering(prots, mcl.out_dir, mcl.name, mcl, checkm = mcl.out_dir + "../checkm.txt")
    bayes_facts = { c : c.core_probability()/c.non_core_probability_plural()  for c in  tqdm(dat) }
    temp = {k:2*log(l) for k,l in bayes_facts.iteritems() if 2*log(l) > 0 }
    for t in tqdm(tt):
        if not os.path.exists(dat.path + "trees/" + t.name) : os.makedirs(dat.path + "trees/" + t.name)
        t.align(dat.path + "trees/" + t.name + "/aligned.faa" )
        t.tree_construction(dat.path + "trees/" + t.name + "/aligned.faa", dat.path + "trees/" + t.name + "/tree.nwk", sccs=True )

    #mcl.post_blast()
    #tho = Clustering([ "/home/moritz/people/olle/thorsellia_reconstruction/proteoms/" + g for g in os.listdir("/home/moritz/people/olle/thorsellia_reconstruction/proteoms/") if ".faa" in g], "/home/moritz/people/olle/thorsellia_reconstruction/", "ThoAnRec")
    #tho.tree_construction(root="PseAer")
