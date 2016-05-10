# Futures #
from __future__ import division

# Built-in modules #
from itertools import product

# External modules #
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from tqdm import tqdm
from numpy.random import permutation
from numpy.random import seed
import sh
import os, glob


seed(1234)

illum_bcs = ["TAAGGCGA",
             "CGTACTAG",
             "AGGCAGAA",
             "TCCTGAGC",
             "GGACTCCT",
             "TAGGCATG",
             "CTCTCTAC",
             "CAGAGAGG",
             "GCTACGCT",
             "CGAGGCTG",
             "AAGAGGCA",
             "GTAGAGGA",
             "TAGATCGC",
             "CTCTCTAT",
             "TATCCTCT",
             "AGAGTAGA",
             "GTAAGGAG",
             "ACTGCATA",
             "AAGGAGTA",
             "CTAAGCCT"]


illum_names = {"TAAGGCGA" : "N701",
               "CGTACTAG" : "N702",
               "AGGCAGAA" : "N703",
               "TCCTGAGC" : "N704",
               "GGACTCCT" : "N705",
               "TAGGCATG" : "N706",
               "CTCTCTAC" : "N707",
               "CAGAGAGG" : "N708",
               "GCTACGCT" : "N709",
               "CGAGGCTG" : "N710",
               "AAGAGGCA" : "N711",
               "GTAGAGGA" : "N712",
               "TAGATCGC" : "N501",
               "CTCTCTAT" : "N502",
               "TATCCTCT" : "N503",
               "AGAGTAGA" : "N504",
               "GTAAGGAG" : "N505",
               "ACTGCATA" : "N506",
               "AAGGAGTA" : "N507",
               "CTAAGCCT" : "N508"}




tetramers = ["".join(tetramer) for tetramer in product('ACGT', repeat = 8)]
tetramers = list(permutation(tetramers))

print "putting the already used barcodes to the end"

for i in tqdm(range(len(tetramers))):
    seq = tetramers.pop(0)
    if seq not in illum_bcs:
        tetramers.append(seq)

tetramers.extend(illum_bcs)

# first path, removing bad gc and symetricality
print "total len: " + str(len(tetramers))

for i in tqdm(range(len(tetramers))) :
    candidate = Seq(tetramers.pop(0))
    gc_content = candidate.count('G') + candidate.count('C')
    gc_content = gc_content / len(candidate)
    if gc_content > 0.6 or gc_content < 0.3:
        continue

    if str(candidate.reverse_complement()) in tetramers:
        continue

    if str(candidate.complement().reverse_complement()) in tetramers :
        continue

    if str(candidate.reverse_complement()) == str(candidate):
        continue

    if str(candidate.complement().reverse_complement()) == str(candidate) :
        continue


    tetramers.append(str(candidate))

print ""
print "First weeding: "
print "list len: " + str(len(tetramers))



def hdist(c1,c2):
    dist = 0
    for i in range(len(c1)):
        if c1[i] != c2[i]:
            dist = dist + 1
    return dist

tetramers2 = list(permutation(tetramers))

for i in tqdm(range(len(tetramers2))):
    seq = tetramers2.pop(0)
    if seq not in illum_bcs:
        tetramers2.append(seq)

tetramers2.extend(illum_bcs)



def weed2(tetras):
    tetras2 = list(tetras)
    for i in tqdm(range(len(tetras2))) :
        candidate = tetras2.pop(0)
        for other in tetras2:
            dist = hdist(candidate,other)
            if dist < 2:
                break
        if dist < 2 :
            continue
        tetras2.append(candidate)
    return tetras2

tetramers2=weed2(tetramers2)

for i in tqdm(range(len(tetramers2))):
    seq = tetramers2.pop(0)
    if seq not in illum_bcs:
        tetramers2.append(seq)

tetramers2.extend(illum_bcs)


print "calculating dG"

dGs_fwd={}
for seq in tqdm(tetramers2):
        with open("temp.seq","w") as file:
                head="AATGATACGGCGACCACCGAGATCTACAC"
                tail="ACACTCTTTCCCTACACGACG"
                file.write(head + seq + tail)
        try:
                sh.mfold("SEQ=temp.seq","NA=DNA")
                result=str(sh.grep("dG","temp.out"))
                for x in glob.glob("temp*"): sh.rm(x)
                sh.rm("date.test")
                dGs_fwd[seq] = float(result.split("=")[1].split("dH")[0])
        except:
                for x in glob.glob("*.37"): sh.rm(x)
                for x in glob.glob("temp.*"): sh.rm(x)
                sh.rm("date.test")
                dGs_fwd[seq] = float("inf")

dGs_rew={}
for seq in tqdm(tetramers2):
        with open("temp.seq","w") as file:
            
                head="CAAGCAGAAGACGGCATACGAGAT"
                tail="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
                file.write(head + seq + tail)
        try:
                sh.mfold("SEQ=temp.seq","NA=DNA")
                result=str(sh.grep("dG","temp.out"))
                for x in glob.glob("temp.*"): sh.rm(x)
                sh.rm("date.test")
                dGs_rew[seq] = float(result.split("=")[1].split("dH")[0])
        except:
                for x in glob.glob("*.37"): sh.rm(x)
                for x in glob.glob("temp*"): sh.rm(x)
                sh.rm("date.test")
                dGs_rew[seq] = float("inf")

                
        
print ""
print "Second weeding: "
print "list len: " + str(len(tetramers2))
dGs_tot = {k : dGs_fwd[k] + dGs_rew[k] for k in dGs_fwd}

temp = dGs_tot.values()
temp.sort()
temp = temp[-189]
selecteds = [x for (x,z) in dGs_tot.iteritems() if z > temp and x not in illum_bcs]
selecteds.extend(illum_bcs)

<<<<<<< HEAD
selecteds = [x for (x,z) in dGs.iteritems() if z == float("inf") or x in illum_bcs]

for i in tqdm(range(len(selecteds))):
    seq = selecteds.pop(0)
    if seq not in illum_bcs:
        selecteds.append(seq)

selecteds.extend(illum_bcs)


=======
>>>>>>> c682bc734bbd5df733589ebb3fefe32aae3ec7d3
selecteds.reverse()
names =[illum_names[seq] if seq in illum_names.keys() else str(i+1) for (i,seq) in enumerate(selecteds)]

with open("barcodes.txt","w") as file:
    file.write("".join([ ">barcode_" + names[i] + "\n" + seq +"\n" for (i,seq) in enumerate(selecteds) ]))
