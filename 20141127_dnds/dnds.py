from Bio import SeqIO
from pandas import DataFrame
import pandas
from tqdm import tqdm
import sh
import os
from Bio import SeqIO
from pandas import DataFrame
import json
from collections import Counter
from numpy import mean
from numpy import argmin
from numpy import argmax
from Bio.Seq import Seq
from Bio.Data.CodonTable import standard_dna_table
import re


def get_seq(g, cutof=0.95, rare_filter = 0.05):
    vec = ['fA','fT','fG','fC']
    comp= {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    codon_table = standard_dna_table.forward_table
    for c in standard_dna_table.stop_codons:
        codon_table[c] = "*"   
    start = pos_dict[g]['start']-1
    end = pos_dict[g]['end']
    c = pos_dict[g]['contig']
    sense = pos_dict[g]['sense']
    dat = nucl_compo[c][start:end]
    if sense == "-":
         dat = dat[::-1]
    str_seq = "".join([v['base'] for v in dat])
    max_freq = [ max([v['fA'],v['fT'],v['fG'],v['fC']]) for v in dat if v['base'] != 'X']
    rate = mean(max_freq)
    seq = Seq(str_seq)
    if sense == "-":
        seq = seq.complement()
    variants = [(i,dat[i]['base'],[v.split("f")[1] for v in vec if dat[i]['base']!='X' and dat[i][v] > rare_filter and dat[i]['base'] not in v]) for i,f in enumerate(max_freq) if f<cutof if dat[i]['coverage'] != 0]
    if sense == "-":
        variants = [(a, comp[b], [comp[d] for d in c]) for  a,b,c in variants]
    codons = re.findall('...',str(seq))

    syn = 0
    non = 0
    stop = False

    for i,b,vv in variants:
        off = i-(i/3)*3
        codon = codons[i/3]
        if codon_table.has_key(codon):
            aa = codon_table[codon]
            for v in vv:
                c = list(codon)
                c[off]=v
                if codon_table["".join(c)] == aa:
                    syn += 1
                else:
                    non += 1
                    if codon_table["".join(c)] == "*":
                        stop = True
                
    return {"rate": rate,"syn":syn, "non":non , "stop": stop, "snp_freq" : float(syn+non)/float(end-start), "len" : end-start  }


with open("bin_7.gff") as handle:
    gff = handle.readlines()

gff = [g.split("\t") for g in gff if g[0] != "#"]
gff = [g for g in gff if len(g) == 9]
gff = [g for g in gff if g[2] != "gene"]

pos_dict = {g[8].split(";")[0].split("=")[1] : {'contig': g[0], 'start': int(g[3]), 'end': int(g[4]) , 'sense' : g[6], 'rna' : "RNA" in g[8]}  for g in gff}

genes =genes = "bin_7.fna"
with open(genes) as handle:
    seq_dict = {s.name : s for s in SeqIO.parse(handle,"fasta")}

with open("readcount.txt") as handle:
    data=handle.readlines()

prots = "bin_7.faa"
with open(prots) as handle:
    aa_dict = {s.name : s for s in SeqIO.parse(handle,"fasta")}

    

data_dict = {(l.split("\t")[0], int(l.split("\t")[1])) : l[:-1].split("\t")  for l in data}


nucl_compo = {}
for c in set([v['contig'] for v in pos_dict.values()]):
    nucl_compo[c] = [{} for i in range(len(seq_dict[c]))]

for c in tqdm(set([v['contig'] for v in pos_dict.values()])):
    for i in range(len(seq_dict[c])):
        if not data_dict.has_key((c,i+1)):
            nucl_compo[c][i]['base'] = 'X'
            nucl_compo[c][i]['coverage'] = 0
        else :
            lise = data_dict[(c,i+1)]
            nucl_compo[c][i]['base'] = lise[2]
            nucl_compo[c][i]['coverage'] = lise[3]
            for bi in lise[4:10]:
                b=bi.split(":")
                nucl_compo[c][i][b[0]] = int(b[1])
                nucl_compo[c][i]["f" + b[0]] = float(b[1])/float(lise[3])


mut_rates = {ga : get_seq(ga) for ga in tqdm(pos_dict) if not pos_dict[ga]["rna"]}

