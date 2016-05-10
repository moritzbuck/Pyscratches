""" extracts sequences related to a taxa 
Usage:
  reconcoct.py <taxa> <assigns> <centers> <map> <all_reads>
  reconcoct.py -h
  
Options:
   -h --help       Show this screen.
"""
from docopt import docopt
from Bio import SeqIO
    
if __name__ == '__main__':
    arguments = docopt(__doc__)
    print(arguments)
    taxa = arguments['<taxa>']
    centers = arguments['<centers>']
    assign = arguments['<assigns>']
    mappy = arguments['<map>']
    all_reads = arguments['<all_reads>']

    print "looking for assigneds"    
    with open(assign) as handle:
        otus = [l.split()[0] for l in handle.readlines() if taxa in l]

    print "pull out centers"    
    with open(centers) as handle:
        centers_seqs = [s for s in SeqIO.parse(handle,"fasta") if s.id in otus]

    with open(taxa + "_centers.fasta", "w") as handle:
        SeqIO.write(centers_seqs, handle, "fasta")

    print "pull out read names"
    reads = {o : [] for o in otus}
    all_map_reads = []
    with open(mappy) as handle:
        for l in handle.readlines():
            if l.split()[-1] in otus:
                reads[l.split()[-1]] += [l.split()[-2]]
                all_map_reads += [l.split()[-2]]
                
    print "open all reads"
    with open(all_reads) as handle:
        all_reads = [s for s in SeqIO.parse(handle,"fasta") if s.id in all_map_reads]

    print "write reads of:"
    for k,v in reads.iteritems():
        with open(taxa +"_"+k+ "_reads.fasta","w") as handle:
            print k
            SeqIO.parse([s for s in all_reads if s.id in v],handle, "fasta")
            
