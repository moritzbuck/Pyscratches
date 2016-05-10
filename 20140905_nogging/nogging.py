from pandas import DataFrame
from tqdm import tqdm
import json
from Bio import SeqIO

all_nogg_prots_path = "/home/moritz/glob/data/NOG/eggnogv4.proteins.all.fa"
all_nogg_member = "/home/moritz/glob/data/NOG/NOGs/NOG.members.txt"

nogg_fasta_db = "/home/moritz/glob/data/NOG/NOGs/fastas/"

def generate_nogg_dict():
    all_dat = DataFrame.from_csv(all_nogg_member,sep="\t",header=0,index_col=None)
    nogs_idx = list(set(all_dat["#nog name"].values))
    nogg_dict = { n : [] for n in nogs_idx}
    for k,l in tqdm(all_dat.iterrows()):
        nogg_dict[l['#nog name']].append(l['protein name'])

    with open('nogg_as_dict.json', 'w') as outfile:
        json.dump(nogg_dict, outfile)

def generate_id2nogg_dict():
    all_dat = DataFrame.from_csv(all_nogg_member,sep="\t",header=0,index_col=None)
    prot_idx = list(set(all_dat["protein name"].values))
    id2nogg_dict = { n : [] for n in prot_idx}
    print len(all_dat), "number of lines to be parsed"
    for k,l in tqdm(all_dat.iterrows()):
        id2nogg_dict[l['protein name']].append(l['#nog name'])

    with open('id2nogg_as_dict.json', 'w') as outfile:
        json.dump(nogg_dict, outfile)


         

def make_fastas(max_len = 10000, min_len = 200):
    noggs = [k for k,n in nogg_dict.iteritems() if len(n) < max_len and len(n) > min_len]
    with open(all_nogg_prots_path) as handle:
        print "(harcoded:) 14 875 530 proteins in db"
        all_prots = SeqIO.parse(handle, "fasta")
        prot_dict = { n : [] for n in nogs_idx}
        for i,prot in tqdm(enumerate(all_prots)):
            if(prot.id in prot_idx):
                for nogg in id2nogg_dict[prot.id]:
                    prot_dict[nogg].append(prot)
                if i >5000:
                    return 1

        for nogg in tqdm(nogs_idx):
            with open(nogg_fasta_db + nogg + ".faa","w") as handle:
                SeqIO.write(prot_dict[nogg],handle, "fasta")
