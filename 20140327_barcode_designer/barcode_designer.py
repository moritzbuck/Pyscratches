from __future__ import division
from Bio.Seq import Seq
from itertools import product


def gc_content(kmer):
   return (kmer.count("G")+kmer.count("C"))/len(kmer)


octamers = [Seq("".join(octo)) for octo in product('ACGT', repeat=8)]
octamers = [octs for octs in octamers if gc_content(octs) < 0.60 and gc_content(octs) > 0.30]

weeded=[octamers[0]]
for octo in octamers[1:]:
    if sum([str(octa) == str(octo.reverse_complement()) for octa in weeded]) == 0 :
        weeded.append(octo)


octamers = weeded
weeded=[octamers[0]]
for octo in octamers[1:]:
    if sum([str(octa) == str(octo.complement().reverse_complement()) for octa in weeded]) == 0 :
        weeded.append(octo)
