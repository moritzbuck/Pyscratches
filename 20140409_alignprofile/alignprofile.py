from Bio import SeqIO
from numpy import median
from numpy import mean
from tqdm import tqdm 
import os.path
from pandas import Series


file ="test.fasta"
file ="/home/moritz/dropbox/stability.trim.contigs.good.unique.align"

if not os.path.isfile("coverages.csv") :
    print "compute coverages"
    if 'coverage' in locals() : del coverage
    handle = open(file, "rU")
    for record in tqdm(SeqIO.parse(handle, "fasta")) :
        seq = str(record.seq)
        l = len(seq)
        if 'coverage' not in locals():
            coverage = [0]*l

        for (i,c) in enumerate(seq):
            if c not in ['.','-']:
                coverage[i] = coverage[i] +1
    coverage=Series(coverage)
    coverage.to_csv("coverages.csv",index=False)
    handle.close()
else :
    print "import coverages"
    coverage = Series.from_csv("coverages.csv",header=-1, index_col=False)

print "compute median-ish things"
medians = []
means = []
maxs = []
mins = []
lens = []
left = []
right = []
unsure = []
handle = open(file, "rU")
positions=list(coverage[coverage > 500000].index)
l = len(positions)
for record in tqdm(SeqIO.parse(handle, "fasta")) :
    seq = str(record.seq)
    poss =[]
    for (i,c) in enumerate(positions):
        if seq[c] not in ['.','-']:
            poss.append(i)
    if len(poss) >0 :
        medians.append(median(poss))
        means.append(mean(poss))
        mins.append(min(poss))
        maxs.append(max(poss))
        lens.append(len(poss))
        if mean(poss) < 300:
            left.append(record)
        else :
            right.append(record)
    else:
        unsure.append(record)
handle.close()

Series(medians).to_csv("low_cov_removed_median.csv",index=False)
Series(means).to_csv("low_cov_removed_means.csv",index=False)
Series(maxs).to_csv("low_cov_removed_maxs.csv",index=False)
Series(mins).to_csv("low_cov_removed_mins.csv",index=False)
Series(lens).to_csv("low_cov_removed_lens.csv",index=False)

print "Write fastas"

with open("left_side.fasta","w") as lefty:
    SeqIO.write(left,lefty,"fasta")

with open("right_side.fasta","w") as righty:
    SeqIO.write(right,righty,"fasta")

with open("unsure.fasta","w") as unsurey:
    SeqIO.write(unsure,unsurey,"fasta")

