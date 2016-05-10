""" converts ifile from clustalw to phylip


Usage:
  fastasubstracter.py <ifile> <ofile>
  fastasubstracter.py -h
  
Options:
   -h --help       Show this screen.
"""
from docopt import docopt
from Bio import SeqIO

if __name__ == '__main__':
    arguments = docopt(__doc__)
    print(arguments)

    
i_file = arguments['<ifile>']
o_file = arguments['<ofile>']

with open(i_file,"r") as ifile:
    alis = [s for s in SeqIO.parse(ifile,"clustal")]

dups = [s for s in alis if len([s2 for s2 in alis if s2.id[0:10] == s.id[0:10] ]) > 1  ]
sets = set([s.id[0:10] for s in alis if len([s2 for s2 in alis if s2.id[0:10] == s.id[0:10] ]) > 1  ])
alis = [a for a in alis if a not in dups]

for s in sets:
    print "Have to rename the ones starting with " + s
    sub_set = [a for a in dups if a.id[0:10] == s]
    for i in range(1,len(sub_set)+1):
        print("renaming " + sub_set[i-1].id + "to"),
        mod = list(sub_set[i-1].id)
        mod[9] =  str(i)
        sub_set[i-1].id = "".join(mod[0:10])
        print("to " + sub_set[i-1].id)
    alis.extend(sub_set)
        
        
with open(o_file,"w") as ofile:
        SeqIO.write(alis,ofile,"phylip")


