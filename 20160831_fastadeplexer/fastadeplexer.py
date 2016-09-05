import shutil
from tqdm import tqdm

pairs = open("/scratch/pairs.fa","w")
singles = open("/scratch/singles.fa","w")
with open("/scratch/canada_qual_paired.keep.noabund.keep.fa") as handle:
    odds = True
    first = True
    for l in tqdm(handle):
        if first :
            if odds:
                id1 = l
                odds = False
            else :
                seq1 = l
                odds = True
                first = False
        else : 
            if odds :
                id2 = l
                odds = False
            else :
                seq2 = l
                odds = True
                if id1.split()[0] == id2.split()[0]:
                    pairs.writelines([id1,seq1,id2, seq2])
                    first = True
                else :
                    singles.writelines([id2,seq2])
                    seq1 = seq2
                    id1 = id2
                    first = False
#shutil.move("/scratch/pairs.fa", "/home/moritz/pairs.fa")
#shutil.move("/scratch/singles.fa", "/home/moritz/singles.fa")
pairs.close()
singles.close()
